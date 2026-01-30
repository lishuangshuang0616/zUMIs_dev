#!/usr/bin/env python3
import sys
import os
import yaml
import subprocess
import csv
import shutil
import json
import bisect
import math
import collections
import gzip

def load_config(yaml_file):
    with open(yaml_file, 'r') as f:
        return yaml.safe_load(f)

def check_dependencies(samtools_exec, featurecounts_exec):
    def check_one(tool, name):
        if os.path.isabs(tool) or os.path.sep in str(tool):
            if not (os.path.exists(tool) and os.access(tool, os.X_OK)):
                print(f"Error: {name} is not executable: {tool}")
                sys.exit(1)
        else:
            if shutil.which(tool) is None:
                print(f"Error: {name} is not found in PATH: {tool}")
                sys.exit(1)

    check_one(samtools_exec, "samtools")
    check_one(featurecounts_exec, "featureCounts")

def get_bam_chromosomes(bam_file, samtools_exec='samtools'):
    """Reads chromosome names from BAM header."""
    cmd = [samtools_exec, 'view', '-H', bam_file]
    try:
        output = subprocess.check_output(cmd, universal_newlines=True)
        chroms = set()
        for line in output.splitlines():
            if line.startswith('@SQ'):
                parts = line.split('\t')
                for part in parts:
                    if part.startswith('SN:'):
                        chroms.add(part[3:])
        return chroms
    except subprocess.CalledProcessError:
        print(f"Warning: Could not read header from {bam_file}")
        return set()

def load_gene_models(gtf_file):
    """
    Loads exon models from GTF, merges them, and calculates 100 genomic percentile points for each gene.
    Returns: dict {gene_id: {'chrom': str, 'strand': str, 'percentiles': list of 100 ints}}
    """
    print(f"Loading gene models for stats from {gtf_file}...")
    gene_exons = collections.defaultdict(list)
    gene_strand = {}
    gene_chrom = {}
    
    if not os.path.exists(gtf_file):
        print("Warning: GTF file not found. Coverage stats will be skipped.")
        return {}

    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            if parts[2] != 'exon': continue

            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attributes = parts[8]

            gene_id = None
            if 'gene_id "' in attributes:
                gene_id = attributes.split('gene_id "')[1].split('"')[0]
            elif 'gene_id' in attributes:
                # Fallback for some GTF formats
                try:
                    gene_id = attributes.split('gene_id')[1].strip().split(';')[0].strip().strip('"')
                except:
                    pass
            
            if gene_id:
                gene_exons[gene_id].append((start, end))
                gene_strand[gene_id] = strand
                gene_chrom[gene_id] = chrom
    
    models = {}
    for gene_id, exons in gene_exons.items():
        exons.sort()
        merged = []
        if not exons: continue
        
        curr_s, curr_e = exons[0]
        for s, e in exons[1:]:
            if s <= curr_e + 1:
                curr_e = max(curr_e, e)
            else:
                merged.append((curr_s, curr_e))
                curr_s, curr_e = s, e
        merged.append((curr_s, curr_e))
        
        gene_all_base = []
        for s, e in merged:
            gene_all_base.extend(range(s, e + 1))
        
        # Need enough bases for percentile calculation
        if len(gene_all_base) < 100:
            continue
            
        gene_all_base.sort()
        strand = gene_strand.get(gene_id, '+')
        if strand == '-':
            gene_all_base.reverse()
            
        points = []
        size = len(gene_all_base)
        for i in range(1, 101):
            idx = int(math.ceil(size * i / 100.0)) - 1
            points.append(gene_all_base[idx])
            
        models[gene_id] = {
            "chrom": gene_chrom.get(gene_id),
            "strand": strand,
            "percentiles": sorted(points) # Sorted for bisect
        }

    print(f"Loaded {len(models)} gene models.")
    return models

def parse_gtf_and_create_saf(gtf_file, out_prefix, valid_chroms=None):
    """
    Parses GTF, merges exons per gene, and creates a Combined SAF file (Exon + Intron).
    Exons use GeneID. Introns use GeneID__INTRON__.
    Returns path to combined_saf, and a dictionary mapping gene_id to gene_name.
    """
    print(f"Parsing GTF: {gtf_file}...")
    
    genes = {}
    gene_id_to_name = {}
    
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            
            feature_type = parts[2]
            if feature_type != 'exon': continue
            
            chrom = parts[0]
            if valid_chroms and chrom not in valid_chroms:
                continue
                
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attributes = parts[8]
            
            gene_id = None
            if 'gene_id "' in attributes:
                gene_id = attributes.split('gene_id "')[1].split('"')[0]
            elif 'gene_id' in attributes:
                pass
            
            if not gene_id: continue

            # Extract gene_name if available
            gene_name = gene_id # Default to gene_id
            if 'gene_name "' in attributes:
                gene_name = attributes.split('gene_name "')[1].split('"')[0]
            
            # Store mapping
            if gene_id not in gene_id_to_name:
                gene_id_to_name[gene_id] = gene_name
            
            if gene_id not in genes:
                genes[gene_id] = {'chrom': chrom, 'strand': strand, 'intervals': []}
            
            genes[gene_id]['intervals'].append((start, end))
            
    print(f"Loaded {len(genes)} genes. Generating Combined SAF file...")
    
    combined_saf_path = f"{out_prefix}.combined.saf"
    
    with open(combined_saf_path, 'w') as f_out:
        header = "GeneID\tChr\tStart\tEnd\tStrand\n"
        f_out.write(header)
        
        for gene_id, data in genes.items():
            chrom = data['chrom']
            strand = data['strand']
            intervals = sorted(data['intervals'])
            
            merged = []
            if intervals:
                curr_start, curr_end = intervals[0]
                for next_start, next_end in intervals[1:]:
                    if next_start <= curr_end + 1:
                        curr_end = max(curr_end, next_end)
                    else:
                        merged.append((curr_start, curr_end))
                        curr_start, curr_end = next_start, next_end
                merged.append((curr_start, curr_end))
            
            # Write Exons
            for start, end in merged:
                f_out.write(f"{gene_id}\t{chrom}\t{start}\t{end}\t{strand}\n")
            
            # Write Introns
            if len(merged) > 1:
                for i in range(len(merged) - 1):
                    intron_start = merged[i][1] + 1
                    intron_end = merged[i+1][0] - 1
                    
                    if intron_end >= intron_start:
                        if (intron_end - intron_start + 1) > 10:
                            # Use suffix to distinguish
                            f_out.write(f"{gene_id}__INTRON__\t{chrom}\t{intron_start}\t{intron_end}\t{strand}\n")

    return combined_saf_path, gene_id_to_name

def run_featurecounts_cmd(featurecounts_exec, input_bam, saf_file, out_prefix, threads, strand_mode, feature_type):
    """
    Runs featureCounts.
    """
    print(f"Running featureCounts for {feature_type} (Strand: {strand_mode})...")
    output_counts = f"{out_prefix}.counts.txt"
    
    cmd = [
        featurecounts_exec,
        '-M',
        '-a', saf_file,
        '-F', 'SAF',
        '-o', output_counts,
        '-T', str(threads),
        '-R', 'BAM',
        '-s', str(strand_mode),
        '-p',                # Enable Paired-End mode
        '--primary',
        '-Q', '0',
        input_bam
    ]
    cmd.append('--largestOverlap')
    
    subprocess.check_call(cmd)
    
    # Cleanup counts.txt and .summary files (intermediate outputs not needed)
    if os.path.exists(output_counts): os.remove(output_counts)
    if os.path.exists(output_counts + ".summary"): os.remove(output_counts + ".summary")
    
    generated_bam = f"{input_bam}.featureCounts.bam"
    if not os.path.exists(generated_bam):
        raise FileNotFoundError(f"featureCounts did not generate {generated_bam}")
        
    target_bam = f"{out_prefix}.bam"
    os.rename(generated_bam, target_bam)
    
    return target_bam

def process_bam_and_calculate_stats(input_bam, out_bam, samtools_exec, threads=4, gene_map=None, source_label=None, gene_models=None):
    """
    Processes a single BAM from Combined SAF featureCounts.
    Parses XT tag to distinguish Exon (GeneID) vs Intron (GeneID__INTRON__).
    Adds RE:Z:E/N/I tag.
    Adds GN:Z:GeneName tag.
    Calculates Stats on the fly.
    """
    print(f"Processing BAM {input_bam} -> {out_bam} (Source: {source_label}, Threads: {threads})...")
    
    read_stats = collections.defaultdict(lambda: collections.defaultdict(int))
    cov_arr = [0] * 100
    cov_count = 0
    MAX_COV_READS = 500000

    def update_stats(read_obj, category, source_lbl):
        bc = None
        if isinstance(read_obj, str):
            cb_idx = read_obj.find("CB:Z:")
            if cb_idx != -1:
                end_cb = read_obj.find('\t', cb_idx)
                if end_cb == -1: end_cb = len(read_obj)
                bc = read_obj[cb_idx+5:end_cb].strip()
            
            if bc:
                parts = read_obj.split('\t')
                flag = int(parts[1])
                is_r1 = (flag & 0x40) != 0
                is_r2 = (flag & 0x80) != 0
                if not (is_r1 or is_r2): is_r1 = True
                
                if is_r1:
                    read_stats[bc][category] += 1
                    if source_lbl == 'UMI': read_stats[bc]['UMI_Reads'] += 1
                    elif source_lbl == 'Internal': read_stats[bc]['Internal_Reads'] += 1
                    
        else: # pysam
            if read_obj.has_tag("CB"):
                bc = read_obj.get_tag("CB")
                if read_obj.is_read1:
                    read_stats[bc][category] += 1
                    if source_lbl == 'UMI': read_stats[bc]['UMI_Reads'] += 1
                    elif source_lbl == 'Internal': read_stats[bc]['Internal_Reads'] += 1
            else:
                 if read_obj.is_read1:
                     read_stats["__NO_CB__"]["Unused BC"] += 1

    def update_coverage(read_obj, gene_models, cov_arr):
        if not gene_models: return False
        
        gene_id = None
        if isinstance(read_obj, str):
            gx_idx = read_obj.find("GX:Z:")
            if gx_idx != -1:
                end_gx = read_obj.find('\t', gx_idx)
                if end_gx == -1: end_gx = len(read_obj)
                gene_id = read_obj[gx_idx+5:end_gx].strip()
        else:
            if read_obj.has_tag('GX'):
                gene_id = read_obj.get_tag('GX')
        
        if not gene_id: return False
        
        model = gene_models.get(gene_id)
        if not model: return False
        
        blocks = []
        if isinstance(read_obj, str):
            parts = read_obj.split('\t')
            pos = int(parts[3])
            cigar = parts[5]
            curr_pos = pos
            num = 0
            for ch in cigar:
                if ch.isdigit():
                    num = num * 10 + int(ch)
                else:
                    if ch in "M=X":
                        blocks.append((curr_pos, curr_pos + num))
                        curr_pos += num
                    elif ch in "DN":
                        curr_pos += num
                    elif ch in "SH":
                        pass
                    num = 0
        else:
            blocks = read_obj.get_blocks()
            
        if not blocks: return False
        
        pct_points = model["percentiles"]
        strand_minus = (model["strand"] == '-')
        
        hit = False
        for b_start, b_end in blocks:
            idx_start = bisect.bisect_left(pct_points, b_start)
            idx_end = bisect.bisect_right(pct_points, b_end - 1)
            
            if idx_end > idx_start:
                hit = True
                indices = range(idx_start, idx_end)
                if strand_minus:
                    for i in indices:
                        bin_idx = 99 - i
                        cov_arr[bin_idx] += 1
                else:
                    for i in indices:
                        bin_idx = i
                        cov_arr[bin_idx] += 1
        return hit

    # Try pysam
    try:
        import pysam
        print("Using pysam for BAM processing...")
        
        # Input is likely name sorted from featureCounts
        with pysam.AlignmentFile(input_bam, "rb", threads=int(threads)) as f_in, \
             pysam.AlignmentFile(out_bam, "wb", template=f_in, threads=int(threads)) as f_out:
            
            count = 0
            for read in f_in:
                count += 1
                if count % 1000000 == 0: print(f"Processed {count} reads...", end='\r')
                
                category = "Intergenic"
                final_read = read
                
                if read.has_tag('XT'):
                    xt_val = read.get_tag('XT')
                    
                    if "__INTRON__" in xt_val:
                        # Intron Assignment
                        real_gene = xt_val.split("__INTRON__")[0]
                        final_read.set_tag('GX', real_gene)
                        final_read.set_tag('RE', 'N') # Intron
                        category = "Intron"
                    else:
                        # Exon Assignment
                        final_read.set_tag('GX', xt_val)
                        final_read.set_tag('RE', 'E') # Exon
                        category = "Exon"
                        
                else:
                    # Unassigned
                    status = "Unassigned"
                    if read.has_tag('XS'):
                        xs_val = read.get_tag('XS')
                        if "Unassigned_" in xs_val:
                            status = xs_val.replace("Unassigned_", "")
                        elif xs_val == "Assigned": 
                             pass
                    
                    if status == "NoFeatures":
                         final_read.set_tag('RE', 'I')
                         category = "Intergenic"
                    elif read.is_unmapped:
                         category = "Unmapped"
                    else:
                         category = status

                if source_label:
                    final_read.set_tag('SR', source_label)
                
                if final_read.has_tag('GX'):
                    g_id = final_read.get_tag('GX')
                    if gene_map:
                        g_name = gene_map.get(g_id, g_id)
                        final_read.set_tag('GN', g_name)
                
                # Stats
                update_stats(final_read, category, source_label)
                if cov_count < MAX_COV_READS:
                    if update_coverage(final_read, gene_models, cov_arr):
                        cov_count += 1
                
                f_out.write(final_read)
        
        print("\nProcessing complete (via pysam).")
        return read_stats, cov_arr

    except ImportError:
        print("pysam not found. Falling back to samtools pipe method...")
    except Exception as e:
        print(f"pysam processing failed: {e}. Falling back to samtools pipe...")

    # Fallback to samtools pipe
    buf_size = 64 * 1024 * 1024
    
    cmd_in = [samtools_exec, 'view', '-h', '-@', str(threads), input_bam]
    proc_in = subprocess.Popen(cmd_in, stdout=subprocess.PIPE, text=True, bufsize=buf_size)
    
    cmd_out = [samtools_exec, 'view', '-b', '-@', str(threads), '-o', out_bam, '-']
    proc_out = subprocess.Popen(cmd_out, stdin=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=buf_size)
    
    try:
        count = 0
        for line in proc_in.stdout:
            if line.startswith('@'):
                proc_out.stdin.write(line)
                continue
            
            count += 1
            if count % 1000000 == 0: print(f"Processed {count} reads...", end='\r')
            
            final_line = None
            category = "Intergenic"
            
            has_xt = 'XT:Z:' in line
            
            if has_xt:
                # Find XT value
                xt_start = line.find('XT:Z:') + 5
                xt_end = line.find('\t', xt_start)
                if xt_end == -1: xt_val = line[xt_start:].strip()
                else: xt_val = line[xt_start:xt_end]
                
                if "__INTRON__" in xt_val:
                    real_gene = xt_val.split("__INTRON__")[0]
                    final_line = line.rstrip().replace(f'XT:Z:{xt_val}', f'GX:Z:{real_gene}') + "\tRE:Z:N"
                    category = "Intron"
                else:
                    final_line = line.rstrip().replace(f'XT:Z:{xt_val}', f'GX:Z:{xt_val}') + "\tRE:Z:E"
                    category = "Exon"
            else:
                # Unassigned
                reason = None
                xs_pos = line.find('XS:Z:')
                if xs_pos != -1:
                    xs_end = line.find('\t', xs_pos)
                    reason = line[xs_pos+5 : xs_end] if xs_end != -1 else line[xs_pos+5:].rstrip()
                
                if reason and "Unassigned_" in reason:
                    status = reason.replace("Unassigned_", "")
                    if status == "NoFeatures":
                        final_line = line.rstrip() + "\tRE:Z:I"
                        category = "Intergenic"
                    else:
                        final_line = line.rstrip()
                        category = status
                else:
                    # Check unmapped flag
                    parts = line.split('\t')
                    flag = int(parts[1])
                    if not (flag & 0x4):
                         final_line = line.rstrip() + "\tRE:Z:I"
                         category = "Intergenic"
                    else:
                         final_line = line.rstrip()
                         category = "Unmapped"

            if source_label:
                final_line += f"\tSR:Z:{source_label}"
            
            # GN Tag
            if gene_map and ('GX:Z:' in final_line):
                start_gx = final_line.find('GX:Z:') + 5
                end_gx = final_line.find('\t', start_gx)
                gene_id = final_line[start_gx:] if end_gx == -1 else final_line[start_gx:end_gx]
                
                gene_name = gene_map.get(gene_id, gene_id)
                final_line += f"\tGN:Z:{gene_name}"

            # Stats
            update_stats(final_line, category, source_label)
            if cov_count < MAX_COV_READS:
                 if update_coverage(final_line, gene_models, cov_arr):
                     cov_count += 1

            proc_out.stdin.write(final_line + "\n")
            
    except BrokenPipeError:
        outs, errs = proc_out.communicate()
        if errs: print(errs)
        raise
    finally:
        if proc_in.stdout: proc_in.stdout.close()
        if proc_out.stdin: proc_out.stdin.close()
        proc_in.wait()
        proc_out.wait()
    
    print("\nProcessing complete.")
    return read_stats, cov_arr

def split_bam_smartseq3(bam_file, threads, samtools_exec):
    print("Splitting BAM for Smart-seq3 processing (One-pass Optimized)...")
    prefix = bam_file.replace('.bam', '')
    umi_bam = f"{prefix}.UMI.bam"
    internal_bam = f"{prefix}.internal.bam"
    
    # Method 1: Try pysam (Fastest/Cleanest if installed)
    try:
        import pysam
        print("Using pysam for splitting...")
        read_threads = max(1, int(threads) // 2)
        write_threads = max(1, int(threads) // 4)
        
        with pysam.AlignmentFile(bam_file, "rb", threads=read_threads) as infile:
            with pysam.AlignmentFile(umi_bam, "wb", template=infile, threads=write_threads) as out_umi, \
                 pysam.AlignmentFile(internal_bam, "wb", template=infile, threads=write_threads) as out_int:
                
                for read in infile:
                    if read.has_tag('UR'):
                        val = read.get_tag('UR')
                        if val: 
                            out_umi.write(read)
                        else:
                            out_int.write(read)
                    else:
                        out_int.write(read)
        return internal_bam, umi_bam
    except ImportError:
        pass # Fallback to subprocess

    # Method 2: Subprocess Pipe
    print("pysam not found, using samtools pipe...")
    buf_size = 64 * 1024 * 1024
    
    cmd_in = [samtools_exec, 'view', '-h', '-@', str(max(1, int(threads)//2)), bam_file]
    proc_in = subprocess.Popen(cmd_in, stdout=subprocess.PIPE, text=True, bufsize=buf_size)
    
    cmd_umi = [samtools_exec, 'view', '-b', '-@', str(max(1, int(threads)//4)), '-o', umi_bam, '-']
    proc_umi = subprocess.Popen(cmd_umi, stdin=subprocess.PIPE, text=True, bufsize=buf_size)
    
    cmd_int = [samtools_exec, 'view', '-b', '-@', str(max(1, int(threads)//4)), '-o', internal_bam, '-']
    proc_int = subprocess.Popen(cmd_int, stdin=subprocess.PIPE, text=True, bufsize=buf_size)
    
    try:
        for line in proc_in.stdout:
            if line.startswith('@'):
                proc_umi.stdin.write(line)
                proc_int.stdin.write(line)
                continue
            
            if 'UR:Z:' in line and 'UR:Z:\t' not in line:
                proc_umi.stdin.write(line)
            else:
                proc_int.stdin.write(line)
                
    except BrokenPipeError:
        print("Error: Broken Pipe during BAM splitting.")
        raise
    finally:
        if proc_in.stdout: proc_in.stdout.close()
        if proc_umi.stdin: proc_umi.stdin.close()
        if proc_int.stdin: proc_int.stdin.close()
        proc_in.wait()
        proc_umi.wait()
        proc_int.wait()
        
    return internal_bam, umi_bam

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('yaml_file')
    parser.add_argument('--umi_bam', required=False, help="Aligned UMI BAM")
    parser.add_argument('--internal_bam', required=False, help="Aligned Internal BAM")
    args = parser.parse_args()
        
    yaml_file = args.yaml_file
    config = load_config(yaml_file)
    
    project = config['project']
    out_dir = config['out_dir']
    num_threads = int(config.get('num_threads', 4))
    samtools_exec = config.get('samtools_exec', 'samtools')
    featurecounts_exec = config.get('featureCounts_exec', 'featureCounts')

    check_dependencies(samtools_exec, featurecounts_exec)
    
    gtf_file = os.path.join(out_dir, f"{project}.final_annot.gtf")
    final_bam = os.path.join(out_dir, f"{project}.filtered.Aligned.GeneTagged.bam")
    
    counting_opts = config.get('counting_opts', {})
    
    print(f"Processing Project: {project} with {num_threads} threads.")
    
    umi_bam = args.umi_bam
    internal_bam = args.internal_bam
    
    if not umi_bam or not internal_bam:
        umi_bam = os.path.join(out_dir, f"{project}.filtered.tagged.umi.Aligned.out.bam")
        internal_bam = os.path.join(out_dir, f"{project}.filtered.tagged.internal.Aligned.out.bam")
        
    # Check existence
    if not os.path.exists(umi_bam) and not os.path.exists(internal_bam):
         print("Error: Input BAMs not found.")
         sys.exit(1)

    header_bam = umi_bam if os.path.exists(umi_bam) else internal_bam
    valid_chroms = get_bam_chromosomes(header_bam, samtools_exec)
    
    saf_dir = os.path.join(out_dir, "zUMIs_output", "expression")
    if not os.path.exists(saf_dir): os.makedirs(saf_dir)
        
    saf_prefix = os.path.join(saf_dir, f"{project}")
    combined_saf, gene_map = parse_gtf_and_create_saf(gtf_file, saf_prefix, valid_chroms)
    
    gene_models = load_gene_models(gtf_file)
    
    bams_to_merge = []
    total_read_stats = collections.defaultdict(lambda: collections.defaultdict(int))
    total_cov_umi = [0] * 100
    total_cov_int = [0] * 100
    
    def merge_stats(dest_stats, src_stats):
        for bc, counts in src_stats.items():
            for cat, val in counts.items():
                dest_stats[bc][cat] += val

    def merge_coverage(dest_cov, src_cov):
        for i in range(100):
            dest_cov[i] += src_cov[i]
    
    # Process Internal (Strand 0)
    if os.path.exists(internal_bam):
        print("Running featureCounts for Internal Reads...")
        fc_prefix_int = internal_bam + ".fc"
        fc_out_int = run_featurecounts_cmd(featurecounts_exec, internal_bam, combined_saf, fc_prefix_int, num_threads, 0, "Combined_Internal")
        
        processed_int = internal_bam + ".processed.bam"
        r_stats, cov = process_bam_and_calculate_stats(
            fc_out_int, processed_int, samtools_exec, num_threads, 
            gene_map, source_label="Internal", gene_models=gene_models
        )
        
        merge_stats(total_read_stats, r_stats)
        merge_coverage(total_cov_int, cov)
        os.remove(fc_out_int)
        bams_to_merge.append(processed_int)

    # Process UMI (Strand 1)
    if os.path.exists(umi_bam):
        print("Running featureCounts for UMI Reads...")
        fc_prefix_umi = umi_bam + ".fc"
        fc_out_umi = run_featurecounts_cmd(featurecounts_exec, umi_bam, combined_saf, fc_prefix_umi, num_threads, 1, "Combined_UMI")
        
        processed_umi = umi_bam + ".processed.bam"
        r_stats, cov = process_bam_and_calculate_stats(
            fc_out_umi, processed_umi, samtools_exec, num_threads,
            gene_map, source_label="UMI", gene_models=gene_models
        )
        
        merge_stats(total_read_stats, r_stats)
        merge_coverage(total_cov_umi, cov)
        os.remove(fc_out_umi)
        bams_to_merge.append(processed_umi)

    # Save Stats
    stats_out = os.path.join(out_dir, "zUMIs_output", "stats", f"{project}.read_stats.json")
    if not os.path.exists(os.path.dirname(stats_out)):
        os.makedirs(os.path.dirname(stats_out))
        
    stats_data = {
        "read_stats": total_read_stats,
        "coverage_umi": total_cov_umi,
        "coverage_int": total_cov_int
    }
    with open(stats_out, 'w') as f:
        json.dump(stats_data, f)

    # Final Merge
    if len(bams_to_merge) == 1:
        os.rename(bams_to_merge[0], final_bam)
    elif len(bams_to_merge) > 1:
        print(f"Merging BAMs with {num_threads} threads...")
        cmd = [samtools_exec, 'cat', '-@', str(num_threads), '-o', final_bam] + bams_to_merge
        subprocess.check_call(cmd)
        for b in bams_to_merge:
            if os.path.exists(b): os.remove(b)
    else:
        print("Error: No BAMs processed.")
        sys.exit(1)

    print("FeatureCounts pipeline finished successfully.")

if __name__ == "__main__":
    main()
