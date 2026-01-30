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
        raise e # Critical error: ID map is required for correction
    return id_map, internal_bcs

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--binning', required=True, help="Barcode binning file")
    parser.add_argument('--idmap', required=True, help="ID map file")
    parser.add_argument('--type', choices=['umi', 'internal'], required=True, help="Output type filter")
    parser.add_argument('bam_files', nargs='*', default=['-'], help="Input BAM files")
    args = parser.parse_args()

    # Load Maps
    bc_map = load_bc_map(args.binning)
    id_map, internal_bcs = load_id_map(args.idmap)
    target_type = args.type

    # Open Output (Standard Output) - initialized on first file
    outfile = None

    try:
        for bam_path in args.bam_files:
            # Handle '-' for stdin
            if bam_path == '-':
                f_obj = sys.stdin.buffer
            else:
                f_obj = bam_path # pysam accepts path string

            try:
                infile = pysam.AlignmentFile(f_obj, "rb", check_sq=False)
            except ValueError:
                continue

            # Initialize output using header from first file
            # Output SAM (Uncompressed Text) to stdout for STAR compatibility
            if outfile is None:
                try:
                    # Mode "w" = SAM text. File "-" = stdout.
                    outfile = pysam.AlignmentFile("-", "w", template=infile)
                except (BrokenPipeError, IOError):
                    infile.close()
                    return

            try:
                for read in infile:
                    try:
                        raw_bc = None
                        if read.has_tag("CR"):
                            raw_bc = read.get_tag("CR")
                            if isinstance(raw_bc, str):
                                raw_bc = raw_bc.upper()
                            else:
                                raw_bc = None
                        
                        # If no barcode, default to UMI stream?
                        if not raw_bc:
                            if target_type == 'umi':
                                outfile.write(read)
                            continue

                        correct_bc = bc_map.get(raw_bc, raw_bc)
                        is_internal = correct_bc in internal_bcs
                        
                        # Filter Logic: Only output if matches target type
                        if target_type == 'umi' and is_internal:
                            continue
                        if target_type == 'internal' and not is_internal:
                            continue

                        flag = read.flag
                        seq = read.query_sequence
                        qual = read.query_qualities

                        # Logic copied from correct_BCtag.py
                        if flag == 77 and seq:
                            if not is_internal:
                                if len(seq) > 3:
                                    new_seq = seq[3:]
                                    new_qual = qual[3:] if qual else None
                                    read.query_sequence = new_seq
                                    read.query_qualities = new_qual
                            else:
                                try:
                                    if read.has_tag("UR"):
                                        ub = read.get_tag("UR")
                                    else:
                                        ub = None
                                    
                                    if read.has_tag("UY"):
                                        qu = read.get_tag("UY")
                                    else:
                                        qu = None
                                except KeyError:
                                    ub, qu = None, None
                                if ub:
                                    new_seq = ub + seq
                                    if qual and qu and len(qu) == len(ub):
                                        qu_ints = [ord(c) - 33 for c in qu]
                                        new_qual = qu_ints + list(qual)
                                    else:
                                        new_qual = None
                                    read.query_sequence = new_seq
                                    read.query_qualities = new_qual
                                    read.set_tag("UR", None)
                                    read.set_tag("UY", None)

                        elif flag == 141:
                            if is_internal:
                                read.set_tag("UR", None)
                                read.set_tag("UY", None)

                        read.set_tag("CR", raw_bc)
                        final_bc_out = id_map.get(correct_bc)

                        if final_bc_out:
                            read.set_tag("CC", correct_bc)
                            read.set_tag("CB", final_bc_out)
                        else:
                            read.set_tag("CC", None)
                            read.set_tag("CB", None)

                        outfile.write(read)

                    except (BrokenPipeError, IOError) as e:
                        if e.errno == 32: # EPIPE
                            break
                        else:
                            raise
                    except Exception:
                        if target_type == 'umi':
                            outfile.write(read)
            
            finally:
                infile.close()

    except (BrokenPipeError, KeyboardInterrupt):
        pass
    finally:
        if outfile:
            try:
                outfile.close()
            except (BrokenPipeError, IOError):
                pass

if __name__ == "__main__":
    main()
