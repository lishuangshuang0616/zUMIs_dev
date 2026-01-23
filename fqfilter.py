import sys
import yaml
import subprocess
import re
import os

def get_config(yaml_file):
    with open(yaml_file, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)

def parse_definition(definition):
    if isinstance(definition, list):
        parts = definition
    else:
        parts = str(definition).split(';')
    def_dict = {}
    for part in parts:
        part = str(part).strip()
        match = re.match(r'(\w+)\((.*)\)', part)
        if match:
            key, val = match.groups()
            ranges = []
            for r in val.split(','):
                start, end = map(int, r.split('-'))
                ranges.append((start - 1, end)) # 0-indexed
            def_dict[key] = ranges
    return def_dict

def extract_seq(seq, qual, definition, ss3_no_pattern=False):
    bc_seq = b""
    bc_qual = b""
    umi_seq = b""
    umi_qual = b""
    cdna_seq = b""
    cdna_qual = b""
    
    # Handle BC
    if 'BC' in definition:
        for start, end in definition['BC']:
            bc_seq += seq[start:end]
            bc_qual += qual[start:end]
            
    # Handle UMI
    if 'UMI' in definition:
        if ss3_no_pattern:
            umi_seq = b""
            umi_qual = b""
        else:
            for start, end in definition['UMI']:
                umi_seq += seq[start:end]
                umi_qual += qual[start:end]
                
    # Handle cDNA
    if 'cDNA' in definition:
        # Assuming only one cDNA range as per original perl script logic
        start, end = definition['cDNA'][0]
        if ss3_no_pattern:
            start = 0 # If smart-seq3 pattern not found, take full read
        cdna_seq = seq[start:end]
        cdna_qual = qual[start:end]
        
    return bc_seq, bc_qual, umi_seq, umi_qual, cdna_seq, cdna_qual

def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        return len(s1)
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))

def fastq_iter(handle):
    while True:
        header = handle.readline()
        if not header: break
        seq = handle.readline().rstrip(b'\n\r')
        _plus = handle.readline()
        qual = handle.readline().rstrip(b'\n\r')
        
        if len(seq) != len(qual):
            sys.stderr.write(f"Warning: SEQ/QUAL length mismatch in {header.decode().strip()}: {len(seq)} vs {len(qual)}\n")
            # Pad or truncate qual to match seq length to prevent crashes/offsets
            if len(qual) < len(seq):
                qual += b'#' * (len(seq) - len(qual))
            else:
                qual = qual[:len(seq)]
                
        yield header.rstrip(b'\n\r'), seq, qual

def main():
    if len(sys.argv) < 7:
        print("Usage: python3 fqfilter.py <yaml> <samtools> <rscript> <pigz> <zumis_dir> <tmp_prefix> [--limit N]")
        sys.exit(1)

    yaml_file = sys.argv[1]
    samtools = sys.argv[2]
    pigz = sys.argv[4]
    tmp_prefix = sys.argv[6]
    
    # Parse optional --limit
    read_limit = 0
    if len(sys.argv) > 7 and sys.argv[7] == '--limit':
        try:
            read_limit = int(sys.argv[8])
        except (IndexError, ValueError):
            pass
            
    config = get_config(yaml_file)
    project = config['project']
    out_dir = config['out_dir']
    
    def sorted_sequence_file_entries(sequence_files):
        if isinstance(sequence_files, dict):
            def key_fn(k):
                m = re.search(r'(\d+)', str(k))
                return int(m.group(1)) if m else str(k)

            for k in sorted(sequence_files.keys(), key=key_fn):
                yield sequence_files[k]
            return
        if isinstance(sequence_files, list):
            for entry in sequence_files:
                yield entry

    seq_entries = list(sorted_sequence_file_entries(config.get('sequence_files', {})))
    filenames = [e['name'] for e in seq_entries if isinstance(e, dict) and e.get('name')]
    base_definitions = [parse_definition(e.get('base_definition', [])) for e in seq_entries if isinstance(e, dict)]
    patterns = [e.get('find_pattern', 'character(0)') for e in seq_entries if isinstance(e, dict)]
    pattern_bytes = []
    for p in patterns:
        if p is None or p == 'character(0)':
            pattern_bytes.append(None)
            continue
        p = str(p)
        mm = None
        if ';' in p:
            p, mm = p.split(';', 1)
            mm = int(mm)
        else:
            mm = 1
        pattern_bytes.append((p.encode('ascii'), mm))

    if isinstance(config.get('filter_cutoffs', {}).get('BC_filter'), dict):
        bc_filter = [
            int(config['filter_cutoffs']['BC_filter']['num_bases']),
            int(config['filter_cutoffs']['BC_filter']['phred']),
        ]
    else:
        bc_filter = list(map(int, str(config['filter_cutoffs']['BC_filter']).split()))

    if isinstance(config.get('filter_cutoffs', {}).get('UMI_filter'), dict):
        umi_filter = [
            int(config['filter_cutoffs']['UMI_filter']['num_bases']),
            int(config['filter_cutoffs']['UMI_filter']['phred']),
        ]
    else:
        umi_filter = list(map(int, str(config['filter_cutoffs']['UMI_filter']).split()))
    
    out_bam = os.path.join(out_dir, "zUMIs_output/.tmpMerge", f"{project}.{tmp_prefix}.raw.tagged.bam")
    out_bc_stats = os.path.join(out_dir, "zUMIs_output/.tmpMerge", f"{project}.{tmp_prefix}.BCstats.txt")

    pigz_procs = []
    handles = []
    for f in filenames:
        base_name = os.path.basename(f)
        if base_name.endswith('.gz'):
            base_name = base_name[:-3]
        chunk_path = os.path.join(out_dir, "zUMIs_output/.tmpMerge", f"{base_name}{tmp_prefix}.gz")

        p = subprocess.Popen([pigz, '-p', '2', '-dc', chunk_path], stdout=subprocess.PIPE, text=False, bufsize=1024*1024)
        pigz_procs.append(p)
        handles.append(p.stdout)

    iters = [fastq_iter(h) for h in handles]
    
    bc_stats = {}
    
    out_bam_fh = open(out_bam, 'wb')
    samtools_proc = subprocess.Popen([samtools, 'view', '-Sb', '-'], stdin=subprocess.PIPE, stdout=out_bam_fh)
    bam_out = samtools_proc.stdin

    # PG Header
    pg_line = f"@PG\tID:zUMIs-fqfilter\tPN:zUMIs-fqfilter\tVN:3.0\tCL:python3 fqfilter.py {' '.join(sys.argv[1:])}\n"
    bam_out.write(pg_line.encode("utf-8"))

    total = 0
    filtered = 0
    
    try:
        for records in zip(*iters):
            if read_limit > 0 and total >= read_limit:
                break
                
            total += 1
            
            # Record 0 is usually the one with BC/UMI in this pipeline's convention
            # But we need to aggregate across all files as per perl script
            
            final_bc = b""
            final_bc_q = b""
            final_umi = b""
            final_umi_q = b""
            final_cdna1 = b""
            final_cdna1_q = b""
            final_cdna2 = b""
            final_cdna2_q = b""
            
            go_ahead = True
            layout = "SE"
            
            for i, (_header, seq, qual) in enumerate(records):
                ss3_status = "yespattern"
                pat = pattern_bytes[i] if i < len(pattern_bytes) else None
                if pat is not None:
                    pat_seq, mm = pat
                    if pat_seq == b"ATTGCGCAATG":
                        if hamming_distance(seq[:len(pat_seq)], pat_seq) <= mm:
                            ss3_status = "yespattern"
                        else:
                            ss3_status = "nopattern"
                            go_ahead = False
                    else:
                        if not seq.startswith(pat_seq):
                            go_ahead = False
                
                # Extract parts
                bc, bc_q, umi, umi_q, c1, c1_q = extract_seq(seq, qual, base_definitions[i], ss3_status == "nopattern")
                
                final_bc += bc
                final_bc_q += bc_q
                final_umi += umi
                final_umi_q += umi_q
                
                if i == 0:
                    final_cdna1, final_cdna1_q = c1, c1_q
                else:
                    final_cdna2, final_cdna2_q = c1, c1_q
                    layout = "PE"

            if not go_ahead:
                continue

            # Quality filtering
            def check_qual(q_str, threshold_count, threshold_val):
                low_quals = sum(1 for q in q_str if (q - 33) < threshold_val)
                return low_quals < threshold_count

            if not check_qual(final_bc_q, bc_filter[0], bc_filter[1]):
                continue
            if not check_qual(final_umi_q, umi_filter[0], umi_filter[1]):
                continue
            
            # Success
            filtered += 1
            bc_stats[final_bc] = bc_stats.get(final_bc, 0) + 1
            
            rid = records[0][0].split()[0]
            if rid.startswith(b'@'):
                rid = rid[1:]
            
            tags = (
                b"\tCR:Z:" + final_bc +
                b"\tUR:Z:" + final_umi +
                b"\tCY:Z:" + final_bc_q +
                b"\tUY:Z:" + final_umi_q +
                b"\n"
            )
            
            # Ensure valid SAM fields for SEQ and QUAL
            seq1_out = final_cdna1 if final_cdna1 else b"*"
            qual1_out = final_cdna1_q if final_cdna1_q else b"*"
            
            if seq1_out != b"*" and qual1_out == b"*":
                sys.stderr.write(f"DEBUG: R1 SEQ present but QUAL missing! ID: {rid}\n")
            
            seq2_out = final_cdna2 if final_cdna2 else b"*"
            qual2_out = final_cdna2_q if final_cdna2_q else b"*"

            if layout == "SE":
                line = b"\t".join([
                    rid, b"4", b"*", b"0", b"0", b"*", b"*", b"0", b"0", seq1_out, qual1_out
                ]) + tags
                bam_out.write(line)
            else:
                line1 = b"\t".join([
                    rid, b"77", b"*", b"0", b"0", b"*", b"*", b"0", b"0", seq1_out, qual1_out
                ]) + tags
                line2 = b"\t".join([
                    rid, b"141", b"*", b"0", b"0", b"*", b"*", b"0", b"0", seq2_out, qual2_out
                ]) + tags
                bam_out.write(line1)
                bam_out.write(line2)

    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
    finally:
        try:
            bam_out.close()
        except Exception:
            pass

        samtools_proc.wait()
        out_bam_fh.close()

        for h in handles:
            try:
                h.close()
            except Exception:
                pass

        # We close the pipes above, so pigz might get SIGPIPE.
        # This is expected behavior when we limit reads.
        for p in pigz_procs:
            p.wait()

        if samtools_proc.returncode != 0:
            raise RuntimeError(f"samtools failed (rc={samtools_proc.returncode}) writing {out_bam}")

        for p in pigz_procs:
            # Allow SIGPIPE (141 or -13)
            if p.returncode not in (0, -13, 141):
                raise RuntimeError(f"pigz failed (rc={p.returncode}) while reading chunk(s) for prefix {tmp_prefix}")

    # Write BC stats
    with open(out_bc_stats, 'w') as f:
        for bc, count in bc_stats.items():
            f.write(f"{bc.decode('ascii')}\t{count}\n")

if __name__ == "__main__":
    main()
