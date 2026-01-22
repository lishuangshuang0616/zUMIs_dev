import sys
import os
try:
    import pysam
except ImportError:
    sys.stderr.write("Error: pysam module is required for this script. Please install it (pip install pysam).\n")
    sys.exit(1)

def load_bc_map(binmap_file):
    bc_map = {}
    try:
        with open(binmap_file, 'r') as f:
            for line in f:
                if 'falseBC' in line:
                    continue
                parts = line.replace(',', '\t').split('\t')
                if len(parts) >= 3:
                    raw = parts[0].strip().upper()
                    fixed = parts[2].strip().upper()
                    bc_map[raw] = fixed
    except Exception as e:
        sys.stderr.write(f"Error loading BC map: {e}\n")
    return bc_map

def load_id_map(id_map_file):
    id_map = {}
    internal_bcs = set() # Set to store internal barcodes
    try:
        with open(id_map_file, 'r') as f:
            for line in f:
                if line.startswith('wellID'): continue
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    well_id = parts[0]
                    umi_seqs = parts[1].split(',')
                    int_seqs = parts[2].split(',')
                    
                    for u in umi_seqs:
                        u = u.strip().upper()
                        if not u: continue
                        id_map[u] = well_id
                        
                    for i in int_seqs:
                        i = i.strip().upper()
                        if not i: continue
                        id_map[i] = well_id
                        internal_bcs.add(i)
                            
    except Exception as e:
        sys.stderr.write(f"Error loading ID map: {e}\n")
    return id_map, internal_bcs

def main():
    if len(sys.argv) < 4:
        # Args: 1:in_bam, 2:out_bam, 3:binmap1, 4+:id_map
        print("Usage: python3 correct_BCtag.py <inbam> <outbam> <BCbinmap> [ID_map_file]")
        sys.exit(1)

    in_bam = sys.argv[1]
    out_bam = sys.argv[2]
    binmap = sys.argv[3]  
    id_map_file = sys.argv[4]

    print(f"Loading maps... (BCbinmap: {bool(binmap)}, ID_Map: {bool(id_map_file)})")
    bc_map = load_bc_map(binmap)
    id_map = {}
    internal_bcs = set()
    if id_map_file:
        id_map, internal_bcs = load_id_map(id_map_file)
        print(f"Loaded {len(id_map)} ID mappings and {len(internal_bcs)} internal barcodes.")

    # Open BAM files
    try:
        infile = pysam.AlignmentFile(in_bam, "rb", check_sq=False)
    except ValueError:
        infile = pysam.AlignmentFile(in_bam, "rb", check_sq=False)

    # Prepare header
    header = infile.header.to_dict()
    pg_entry = {
        'ID': 'zUMIs-correct_BCtag',
        'PN': 'zUMIs-correct_BCtag',
        'VN': '3.0-pysam-mgi-custom',
        'CL': 'python3 ' + ' '.join(sys.argv)
    }
    if 'PG' in header:
        header['PG'].append(pg_entry)
    else:
        header['PG'] = [pg_entry]
        
    outfile = pysam.AlignmentFile(out_bam, "wb", header=header)

    processed = 0
    
    for read in infile:
        processed += 1
        if processed % 1000000 == 0:
            print(f"Processed {processed} reads...", flush=True)

        try:
            try:
                raw_bc = read.get_tag("BC")
                if isinstance(raw_bc, str):
                    raw_bc = raw_bc.upper()
                else:
                    raw_bc = None
            except KeyError:
                raw_bc = None
            
            if not raw_bc:
                outfile.write(read)
                continue

            correct_bc = bc_map.get(raw_bc, raw_bc)
            is_internal = correct_bc in internal_bcs
            flag = read.flag

            seq = read.query_sequence
            qual = read.query_qualities  # None or array('B')

            # -------------------------
            # R1: flag == 77
            # -------------------------
            if flag == 77 and seq:
                # ---------- UMI read ----------
                if not is_internal:
                    if len(seq) > 3:
                        new_seq = seq[3:]
                        if qual is not None and len(qual) == len(seq):
                            new_qual = qual[3:]
                        else:
                            new_qual = None
                        read.query_sequence = new_seq
                        read.query_qualities = new_qual

                # ---------- Internal read ----------
                else:
                    try:
                        ub = read.get_tag("UB")
                        qu = read.get_tag("QU")
                    except KeyError:
                        ub, qu = None, None

                    if ub:
                        new_seq = ub + seq
                        if (
                            qual is not None
                            and qu is not None
                            and len(qu) == len(ub)
                        ):
                            qu_ints = [ord(c) - 33 for c in qu]
                            new_qual = qu_ints + list(qual)
                        else:
                            new_qual = None
                        read.query_sequence = new_seq
                        read.query_qualities = new_qual
                        # clear tags
                        read.set_tag("UB", None)
                        read.set_tag("QU", None)

            # -------------------------
            # R2: flag == 141
            # -------------------------
            elif flag == 141:
                if is_internal:
                    read.set_tag("UB", None)
                    read.set_tag("QU", None)

            # -------------------------
            # Set tags
            # -------------------------
            # BX: Original raw barcode (Always kept)
            read.set_tag("BX", raw_bc)
            
            final_bc_out = None
            if id_map:
                final_bc_out = id_map.get(correct_bc)

            if final_bc_out:
                # Valid cell found in ID map
                read.set_tag("BZ", correct_bc)   # Corrected Sequence
                read.set_tag("BC", final_bc_out) # Well ID
            else:
                # Invalid/Junk cell
                # Remove tags if they exist to keep BAM clean
                read.set_tag("BZ", None)
                read.set_tag("BC", None)

            outfile.write(read)

        except Exception as e:
            # sys.stderr.write(f"Warning: error processing read {read.query_name}: {e}\n")
            outfile.write(read)
            
    infile.close()
    outfile.close()

if __name__ == "__main__":
    main()
