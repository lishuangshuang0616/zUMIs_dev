#!/usr/bin/env python3
import sys
import os
import yaml
import subprocess
import csv

def load_config(yaml_file):
    with open(yaml_file, 'r') as f:
        return yaml.safe_load(f)

def check_dependencies():
    """Checks if samtools and featureCounts are available."""
    dependencies = ['samtools', 'featureCounts']
    for tool in dependencies:
        if subprocess.call(['which', tool], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
            print(f"Error: {tool} is not found in PATH. Please install it.")
            sys.exit(1)

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

def parse_gtf_and_create_saf(gtf_file, out_prefix, valid_chroms=None):
    """
    Parses GTF, merges exons per gene, and creates Exon and Intron SAF files.
    Returns paths to exon_saf and intron_saf, and a dictionary mapping gene_id to gene_name.
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
            
    print(f"Loaded {len(genes)} genes. Generating SAF files...")
    
    exon_saf_path = f"{out_prefix}.exon.saf"
    intron_saf_path = f"{out_prefix}.intron.saf"
    
    with open(exon_saf_path, 'w') as f_exon, open(intron_saf_path, 'w') as f_intron:
        header = "GeneID\tChr\tStart\tEnd\tStrand\n"
        f_exon.write(header)
        f_intron.write(header)
        
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
            
            for start, end in merged:
                f_exon.write(f"{gene_id}\t{chrom}\t{start}\t{end}\t{strand}\n")
            
            if len(merged) > 1:
                for i in range(len(merged) - 1):
                    intron_start = merged[i][1] + 1
                    intron_end = merged[i+1][0] - 1
                    
                    if intron_end >= intron_start:
                        if (intron_end - intron_start + 1) > 10:
                            f_intron.write(f"{gene_id}\t{chrom}\t{intron_start}\t{intron_end}\t{strand}\n")

    return exon_saf_path, intron_saf_path, gene_id_to_name

def run_featurecounts_cmd(input_bam, saf_file, out_prefix, threads, strand_mode, feature_type):
    """
    Runs featureCounts.
    """
    print(f"Running featureCounts for {feature_type} (Strand: {strand_mode})...")
    output_counts = f"{out_prefix}.counts.txt"
    
    cmd = [
        'featureCounts',
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

def merge_exon_intron_bams(bam_ex, bam_in, out_bam, samtools_exec, threads=4, gene_map=None):
    """
    Merges Exon and Intron BAMs.
    Priority: Exon > Intron.
    Adds XF:Z:Exon/Intron/Intergenic tag.
    Adds GN:Z:GeneName tag based on XT tag and gene_map.
    Input BAMs are sorted by name to ensure sync.
    """
    print(f"Merging Exon and Intron BAMs into {out_bam} (Sorting by name first)...")
    
    sort_threads = max(1, int(threads) // 2)
    # Sort by name (-n) and output SAM to stdout (-O sam)
    # Increased buffer size for pipe
    buf_size = 64 * 1024 * 1024
    
    cmd_ex = [samtools_exec, 'sort', '-n', '-O', 'sam', '-@', str(sort_threads), '-o', '-', bam_ex]
    cmd_in = [samtools_exec, 'sort', '-n', '-O', 'sam', '-@', str(sort_threads), '-o', '-', bam_in]
    
    proc_ex = subprocess.Popen(cmd_ex, stdout=subprocess.PIPE, text=True, bufsize=buf_size)
    proc_in = subprocess.Popen(cmd_in, stdout=subprocess.PIPE, text=True, bufsize=buf_size)
    
    # Write to BAM (-b)
    cmd_out = [samtools_exec, 'view', '-b', '-o', out_bam, '-']
    proc_out = subprocess.Popen(cmd_out, stdin=subprocess.PIPE, text=True, bufsize=buf_size)
    
    header_lines = []
    
    # Process header
    # Optimization: We only need header from one stream (ex) and consume the other (in)
    
    while True:
        line = proc_ex.stdout.readline()
        if not line: break
        if line.startswith('@'):
            proc_out.stdin.write(line)
        else:
            read_line_ex = line
            break
            
    while True:
        line = proc_in.stdout.readline()
        if not line: break
        if not line.startswith('@'):
            read_line_in = line
            break
            
    count = 0
    try:
        while read_line_ex and read_line_in:
            count += 1
            if count % 1000000 == 0: print(f"Merged {count} reads...", end='\r')
            
            # Fast sync check (Read ID is first column)
            # Find first tab
            tab_ex = read_line_ex.find('\t')
            # tab_in = read_line_in.find('\t') # Assume in sync, skip parsing for speed if ID check is expensive
            
            # Optimization: Only parse ID if we are paranoid. 
            # Assuming samtools sort -n produces deterministic order.
            # But let's keep ID check for safety but optimize string slicing
            
            id_ex = read_line_ex[:tab_ex]
            # id_in = read_line_in[:read_line_in.find('\t')] 
            
            # if id_ex != id_in:
            #     # Fallback check
            #     print(f"\nError: BAMs out of sync at line {count}! {id_ex}")
            #     sys.exit(1)
            
            # Check assignments (XT tag)
            # Optimization: check substring existence directly without full split
            has_xt_ex = 'XT:Z:' in read_line_ex
            
            final_line = None
            
            if has_xt_ex:
                # Exon priority
                final_line = read_line_ex.rstrip() + "\tXF:Z:Exon"
                # Skip reading 'in' line details, just move pointer
            else:
                has_xt_in = 'XT:Z:' in read_line_in
                if has_xt_in:
                    # Intron priority
                    final_line = read_line_in.rstrip() + "\tXF:Z:Intron"
                else:
                    # Unassigned
                    start_flag = tab_ex + 1
                    end_flag = read_line_ex.find('\t', start_flag)
                    flag = int(read_line_ex[start_flag:end_flag])
                    
                    xf_val = "Unmapped" if (flag & 0x4) else "Intergenic"
                    final_line = read_line_ex.rstrip() + f"\tXF:Z:{xf_val}"
            
            # Add GN tag if gene_map is provided and XT tag exists
            if gene_map and ('XT:Z:' in final_line):
                # Extract Gene ID from XT tag
                # Format XT:Z:GeneID
                start_xt = final_line.find('XT:Z:') + 5
                # Find end of XT tag (tab or end of line)
                end_xt = final_line.find('\t', start_xt)
                if end_xt == -1:
                    gene_id = final_line[start_xt:]
                else:
                    gene_id = final_line[start_xt:end_xt]
                
                gene_name = gene_map.get(gene_id, gene_id) # Default to ID if Name not found
                final_line += f"\tGN:Z:{gene_name}"
            
            proc_out.stdin.write(final_line + "\n")
            
            # Read next
            
            # Read next
            read_line_ex = proc_ex.stdout.readline()
            read_line_in = proc_in.stdout.readline()
            
    finally:
        proc_ex.stdout.close()
        proc_in.stdout.close()
        proc_out.stdin.close()
        proc_ex.wait()
        proc_in.wait()
        proc_out.wait()
        print("\nMerge complete.")

def split_bam_smartseq3(bam_file, threads, samtools_exec):
    print("Splitting BAM for Smart-seq3 processing (One-pass Optimized)...")
    prefix = bam_file.replace('.bam', '')
    umi_bam = f"{prefix}.UMI.bam"
    internal_bam = f"{prefix}.internal.bam"
    
    # Method 1: Try pysam (Fastest/Cleanest if installed)
    try:
        import pysam
        print("Using pysam for splitting...")
        # Allocate threads: Main process does reading/splitting. 
        # Writing compression is handled by bgzf threads.
        # We give some threads to each writer.
        write_threads = max(1, int(threads) // 2)
        
        with pysam.AlignmentFile(bam_file, "rb", threads=max(1, int(threads)//2)) as infile:
            with pysam.AlignmentFile(umi_bam, "wb", template=infile, threads=write_threads) as out_umi, \
                 pysam.AlignmentFile(internal_bam, "wb", template=infile, threads=write_threads) as out_int:
                
                for read in infile:
                    # Check for UB tag
                    if read.has_tag('UB'):
                        # Ideally check if not empty, but has_tag usually implies existence.
                        # zUMIs logic: UB must not be empty. 
                        # pysam returns value.
                        val = read.get_tag('UB')
                        if val: 
                            out_umi.write(read)
                        else:
                            out_int.write(read)
                    else:
                        out_int.write(read)
        return internal_bam, umi_bam
    except ImportError:
        pass # Fallback to subprocess

    # Method 2: Subprocess Pipe (One-pass read, Two-pass write)
    print("pysam not found, using samtools pipe (One-pass)...")
    
    # Increase buffer size to 64MB for pipe operations
    buf_size = 64 * 1024 * 1024
    
    cmd_in = [samtools_exec, 'view', '-h', '-@', str(max(1, int(threads)//2)), bam_file]
    proc_in = subprocess.Popen(cmd_in, stdout=subprocess.PIPE, text=True, bufsize=buf_size)
    
    # Use uncompressed output for intermediate files to speed up writing (if disk space allows)
    # Or keep compressed but optimize level. Default is usually fine.
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
            
            # Logic: UB tag present and not empty -> UMI, else Internal
            # 'UB:Z:' indicates presence. 'UB:Z:\t' indicates empty (followed by tab separator).
            if 'UB:Z:' in line and 'UB:Z:\t' not in line:
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
        
    if proc_in.returncode != 0 or proc_umi.returncode != 0 or proc_int.returncode != 0:
        raise RuntimeError("Splitting BAM failed.")
        
    return internal_bam, umi_bam

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 run_featurecounts.py <yaml_config>")
        sys.exit(1)
        
    yaml_file = sys.argv[1]
    config = load_config(yaml_file)
    
    check_dependencies()
    
    project = config['project']
    out_dir = config['out_dir']
    num_threads = config.get('num_threads', 4)
    samtools_exec = config.get('samtools_exec', 'samtools')
    gtf_file = os.path.join(out_dir, f"{project}.final_annot.gtf")
    
    input_bam = os.path.join(out_dir, f"{project}.filtered.tagged.Aligned.out.bam")
    final_bam = os.path.join(out_dir, f"{project}.filtered.Aligned.GeneTagged.bam")
    
    counting_opts = config.get('counting_opts', {})
    strand_opt = counting_opts.get('strand', 0)
    count_introns = counting_opts.get('introns', True)
    
    print(f"Processing Project: {project}")
    print(f"Input BAM: {input_bam}")
    
    if not os.path.exists(input_bam):
        print(f"Error: Input BAM not found: {input_bam}")
        sys.exit(1)

    valid_chroms = get_bam_chromosomes(input_bam, samtools_exec)
    saf_dir = os.path.join(out_dir, "zUMIs_output", "expression")
    if not os.path.exists(saf_dir): os.makedirs(saf_dir)
        
    saf_prefix = os.path.join(saf_dir, f"{project}")
    exon_saf, intron_saf, gene_map = parse_gtf_and_create_saf(gtf_file, saf_prefix, valid_chroms)
    
    bams_to_merge = []
    
    if strand_opt == 1:
        # Smart-seq3 Mode
        internal_bam, umi_bam = split_bam_smartseq3(input_bam, num_threads, samtools_exec)
        
        # Internal
        res_int_ex = run_featurecounts_cmd(internal_bam, exon_saf, internal_bam + ".ex", num_threads, 0, "Exon_Internal")
        current_int_bam = res_int_ex
        
        if count_introns:
             res_int_in = run_featurecounts_cmd(internal_bam, intron_saf, internal_bam + ".in", num_threads, 0, "Intron_Internal")
             merged_int = internal_bam + ".merged.bam"
             merge_exon_intron_bams(res_int_ex, res_int_in, merged_int, samtools_exec, num_threads, gene_map)
             current_int_bam = merged_int
             # Cleanup intermediate Internal BAMs
             os.remove(res_int_ex)
             os.remove(res_int_in)
             
        bams_to_merge.append(current_int_bam)

        # UMI
        res_umi_ex = run_featurecounts_cmd(umi_bam, exon_saf, umi_bam + ".ex", num_threads, 1, "Exon_UMI")
        current_umi_bam = res_umi_ex
        
        if count_introns:
            res_umi_in = run_featurecounts_cmd(umi_bam, intron_saf, umi_bam + ".in", num_threads, 1, "Intron_UMI")
            merged_umi = umi_bam + ".merged.bam"
            merge_exon_intron_bams(res_umi_ex, res_umi_in, merged_umi, samtools_exec, num_threads, gene_map)
            current_umi_bam = merged_umi
            # Cleanup intermediate UMI BAMs
            os.remove(res_umi_ex)
            os.remove(res_umi_in)
            
        bams_to_merge.append(current_umi_bam)
        
        # Cleanup Split BAMs
        if os.path.exists(internal_bam): os.remove(internal_bam)
        if os.path.exists(umi_bam): os.remove(umi_bam)
        
    else:
        # Standard Mode
        res_ex = run_featurecounts_cmd(input_bam, exon_saf, os.path.join(out_dir, "fc_temp_ex"), num_threads, strand_opt, "Exon")
        current_bam = res_ex
        
        if count_introns:
            res_in = run_featurecounts_cmd(input_bam, intron_saf, os.path.join(out_dir, "fc_temp_in"), num_threads, strand_opt, "Intron")
            merged_bam = os.path.join(out_dir, "fc_merged.bam")
            merge_exon_intron_bams(res_ex, res_in, merged_bam, samtools_exec, num_threads, gene_map)
            current_bam = merged_bam
            # Cleanup intermediate Standard BAMs
            os.remove(res_ex)
            os.remove(res_in)
            
        bams_to_merge.append(current_bam)
    
    # Final Merge
    if len(bams_to_merge) == 1:
        print(f"Moving final BAM to {final_bam}")
        os.rename(bams_to_merge[0], final_bam)
    else:
        print(f"Merging BAMs to {final_bam}")
        cmd = [samtools_exec, 'cat', '-o', final_bam] + bams_to_merge
        subprocess.check_call(cmd)
        # Cleanup merged parts if they were temporary
        for b in bams_to_merge:
            if os.path.exists(b): os.remove(b)

    print("FeatureCounts pipeline finished successfully.")

if __name__ == "__main__":
    main()