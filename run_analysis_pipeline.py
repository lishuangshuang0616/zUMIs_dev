#!/usr/bin/env python3
#-*-coding:utf-8-*- 

import os
import argparse
import yaml
import copy
import sys
import subprocess
import shutil
from datetime import datetime
import multiprocessing
import gzip
import math
import glob

# Import constants
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'src'))
import constant

def get_read_length(fastq_path):
    """
    Reads the first sequence from a FASTQ file to determine read length.
    Handles gzip if needed.
    """
    try:
        if fastq_path.endswith('.gz'):
            opener = gzip.open
        else:
            opener = open
            
        with opener(fastq_path, 'rt') as f:
            # Skip header
            f.readline()
            # Read sequence
            seq = f.readline().strip()
            if not seq:
                raise ValueError(f"Empty sequence in {fastq_path}")
            return len(seq)
    except Exception as e:
        raise ValueError(f"Failed to determine read length for {fastq_path}: {e}")

def make_dir(data):
    out_path = data['out_dir'] # This is XPRESS_PROCESSING
    root_path = os.path.dirname(out_path) # Project root
    outs_path = os.path.join(root_path, 'outs')
    
    if os.path.exists(out_path):
        print(f"Warning: Processing directory '{out_path}' already exists. Resuming/Overwriting analysis.")
    
    os.makedirs(out_path, exist_ok=True)
    os.makedirs(outs_path, exist_ok=True)
    os.makedirs(os.path.join(out_path, 'config'), exist_ok=True)
    
    print(f"Directory 'XPRESS_PROCESSING' (out_dir) created/verified at: {out_path}")
    print(f"Directory 'outs' created/verified at: {outs_path}")
    
    # Create zUMIs output dirs inside XPRESS_PROCESSING
    for dirs in ['zUMIs_output', 'zUMIs_output/expression/downsampling', 'zUMIs_output/stats', 'zUMIs_output/.tmpMerge']:
        os.makedirs(f'{out_path}/{dirs}', exist_ok=True)

def create_barcode(data):
    sample_type = data['sample']['sample_type'].lower()
    out_path = data['out_dir'] # XPRESS_PROCESSING
    script_path = data['zUMIs_directory']
    
    # Custom/External Mode
    if sample_type == 'custom' or sample_type == 'external':
        provided_bc = data['barcodes']['barcode_file']
        print(f"Using custom barcode file: {provided_bc}")
        dest_summary = os.path.join(out_path, 'config', 'expect_id_barcode.tsv')
        dest_pipe = os.path.join(out_path, 'config', 'expect_barcode.tsv')
        
        shutil.copy(provided_bc, dest_summary)
        
        # Extract just barcodes for pipe file
        with open(provided_bc, 'r') as infile, open(dest_pipe, 'w') as outfile:
            for line in infile:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                     if parts[0].lower() == 'wellid': continue
                     outfile.write('\t'.join(parts[1:]) + '\n')
        
        # Update barcode file path in config
        data['barcodes']['barcode_file'] = dest_pipe
        return

    sample_id_str = str(data['sample']['sample_id'])
    sample_ids = [s.strip() for s in sample_id_str.split(',')]

    with open(os.path.join(out_path, 'config', 'expect_barcode.tsv'), 'w') as pipe_file, \
         open(os.path.join(out_path, 'config', 'expect_id_barcode.tsv'), 'w') as summary_file:
        
        print('\t'.join(['wellID','umi_barcodes','internal_barcodes']), file=summary_file)

        if sample_type == 'manual':
            # Manual mode: sample_id corresponds to manual_barcode_list.yaml keys (e.g. 20, 21...)
            with open(f'{script_path}/yaml/manual_barcode_list.yaml', 'r', encoding='utf-8') as y:
                barcode_set = yaml.safe_load(y)
            
            for c in sample_ids:
                if c not in barcode_set: continue
                bc_t_i5 = barcode_set[c][0]
                bc_n_i5 = barcode_set[c][1]
                bc_n_i7 = barcode_set[c][2]
                
                # Generate all combinations
                umi_bcs = []
                int_bcs = []
                for k, v in zip(bc_t_i5, bc_n_i5):
                    for j in bc_n_i7:
                        umi_bcs.append(k + j)
                        int_bcs.append(v + j)
                
                well_id = f"MANUAL{c}"
                umi_str = ",".join(umi_bcs)
                int_str = ",".join(int_bcs)
                
                # Write one line per manual sample
                print(f"{well_id}\t{umi_str}\t{int_str}", file=summary_file)
                print(f"{umi_str}\t{int_str}", file=pipe_file)
        else:
            # Auto Mode (Multi-plate support): sample_id corresponds to plate number (e.g. 1, 2...)
            with open(f'{script_path}/yaml/auto_barcode_list.yaml', 'r', encoding='utf-8') as y:
                barcode_set = yaml.safe_load(y)
            
            for sid in sample_ids:
                plate_key = f'plate{sid}'
                if plate_key not in barcode_set:
                    raise ValueError(f"Plate ID {sid} not found in auto_barcode_list.yaml")
                
                for well in barcode_set[plate_key]:
                    well_id = f"P{sid}{well}"
                    # Assuming auto_barcode_list values are already lists or strings
                    val = barcode_set[plate_key][well]
                    # Format depends on structure: [umi_seq, int_seq] or [[umi1, umi2], [int1, int2]]
                    umi_part = val[0]
                    int_part = val[1]
                    
                    umi_str = ",".join(umi_part) if isinstance(umi_part, list) else str(umi_part)
                    int_str = ",".join(int_part) if isinstance(int_part, list) else str(int_part)
                    
                    print(f"{well_id}\t{umi_str}\t{int_str}", file=summary_file)
                    print(f"{umi_str}\t{int_str}", file=pipe_file)    
    
    # Update barcode file path in config
    data['barcodes']['barcode_file'] = os.path.join(out_path, 'config', 'expect_barcode.tsv')

def check_file_exists(data):
    files_to_check=[data['reference']['STAR_index'], data['reference']['GTF_file']]
    for f in files_to_check:
        if not f or not os.path.exists(f):
            raise FileNotFoundError(f'Reference file not found: {f}')
    
    # Check sequence files
    fq1 = data['sequence_files']['file1']['name']
    fq2 = data['sequence_files']['file2']['name']
    
    for f in [fq1, fq2]:
        if f:
             for subf in f.split(','):
                 if not os.path.exists(subf.strip()):
                     raise FileNotFoundError(f'Fastq file not found: {subf}')

def process_fq(data):
    # User requested to use original data location, skipping symlinks.
    pass

def run_pipeline_stages(yaml_file):
    """
    Orchestrates the pipeline stages.
    """
    print(f"Loading config from {yaml_file}...")
    with open(yaml_file, 'r') as f:
        config = yaml.safe_load(f)

    project = config['project']
    out_dir = config['out_dir']
    num_threads = int(config['num_threads'])
    which_stage = config['which_Stage']
    
    # Executables
    samtools = config.get('samtools_exec', 'samtools')
    pigz = config.get('pigz_exec', 'pigz')
    seqkit = config.get('seqkit_exec', 'seqkit')
    zumis_dir = config.get('zUMIs_directory', '.')

    exec_env = os.environ.copy()
    software_dir = os.path.join(zumis_dir, 'software')
    if sys.platform.startswith('linux') and os.path.isdir(software_dir):
        exec_env['PATH'] = software_dir + os.pathsep + exec_env.get('PATH', '')
    zumis_src_dir = os.path.join(zumis_dir, 'src')
    if os.path.isdir(zumis_src_dir) and zumis_src_dir not in sys.path:
        sys.path.insert(0, zumis_src_dir)
    import pipeline_modules

    def resolve_script(script_name):
        direct = os.path.join(zumis_dir, script_name)
        if os.path.exists(direct):
            return direct
        in_src = os.path.join(zumis_dir, 'src', script_name)
        if os.path.exists(in_src):
            return in_src
        raise FileNotFoundError(f"Script not found: {script_name}. Tried: {direct}, {in_src}")
    
    class Tee:
        def __init__(self, *streams):
            self._streams = streams

        def write(self, s):
            for stream in self._streams:
                stream.write(s)

        def flush(self):
            for stream in self._streams:
                stream.flush()

    def estimate_avg_line_len(path, sample_lines=1000):
        opener = gzip.open if path.endswith('.gz') else open
        total = 0
        n = 0
        with opener(path, 'rb') as fh:
            for _ in range(sample_lines):
                line = fh.readline()
                if not line:
                    break
                total += len(line)
                n += 1
        return (total / n) if n else 0.0

    
    # config['out_dir'] is already set to XPRESS_PROCESSING in main()
    analysis_dir = out_dir
    log_path = os.path.join(analysis_dir, 'zUMIs_run.log')
    tmp_merge_dir = os.path.join(analysis_dir, 'zUMIs_output', '.tmpMerge')
    
    # Ensure log directory exists
    os.makedirs(os.path.dirname(log_path), exist_ok=True)

    original_stdout = sys.stdout
    with open(log_path, 'a') as run_log:
        sys.stdout = Tee(original_stdout, run_log)
        try:
            print(f"Starting Pipeline for project: {project}")
            print(f"Stage: {which_stage}")

            chunk_suffixes = []

            def run_stage_cmd(cmd, stage_name, shell=False):
                if isinstance(cmd, list) and not shell:
                    cmd_str = " ".join(cmd)
                else:
                    cmd_str = str(cmd)
                res = subprocess.run(cmd, stdout=run_log, stderr=subprocess.STDOUT, shell=shell, env=exec_env)
                if res.returncode != 0:
                    run_log.flush()
                    try:
                        with open(log_path, 'r') as lr:
                            print(f"\n[ERROR] {stage_name} failed (rc={res.returncode}). Last 30 lines of log ({log_path}):\n", file=sys.stderr)
                            print("".join(lr.readlines()[-30:]), file=sys.stderr)
                    except Exception:
                        pass
                    raise RuntimeError(f"{stage_name} failed with exit code {res.returncode}.")

            def remove_path(path):
                if not path:
                    return
                if os.path.exists(path):
                    os.remove(path)
                bai = path + ".bai"
                if os.path.exists(bai):
                    os.remove(bai)

            if which_stage == "Filtering":
                print(">>> Starting Filtering Stage")
                
                f1_str = config.get('sequence_files', {}).get('file1', {}).get('name', '')
                f2_str = config.get('sequence_files', {}).get('file2', {}).get('name', '')
                
                fq1_files = [f.strip() for f in f1_str.split(',')] if f1_str else []
                fq2_files = [f.strip() for f in f2_str.split(',')] if f2_str else []

                if not fq1_files:
                    raise ValueError("No file1 found in YAML configuration.")

                total_size_bytes = 0
                first_fq = fq1_files[0]
                
                for f in fq1_files:
                    if f.endswith('.gz'):
                        total_size_bytes += os.path.getsize(f) * 3
                    else:
                        total_size_bytes += os.path.getsize(f)

                avg_line_len = estimate_avg_line_len(first_fq, sample_lines=1000)
                if avg_line_len <= 0:
                    raise ValueError(f"Failed to estimate average line length for {first_fq}")

                total_lines_est = total_size_bytes / avg_line_len
                
                lines_per_chunk = int(math.ceil(total_lines_est / num_threads))
                rem = lines_per_chunk % 4
                if rem != 0: lines_per_chunk += (4 - rem)
                if lines_per_chunk < 4000: lines_per_chunk = 4000
                
                print(f"Total input estimation: {int(total_lines_est/4)} reads.")
                print(f"Split config: {lines_per_chunk} lines per chunk.")

                pool = multiprocessing.Pool(processes=min(2, num_threads)) 
                results = []
                
                if fq2_files:
                    res = pool.apply_async(
                        pipeline_modules.split_fastq,
                        (fq1_files, num_threads, lines_per_chunk, tmp_merge_dir, project, pigz, seqkit, fq2_files),
                    )
                    results.append(res)
                else:
                    res = pool.apply_async(
                        pipeline_modules.split_fastq,
                        (fq1_files, num_threads, lines_per_chunk, tmp_merge_dir, project, pigz, seqkit),
                    )
                    results.append(res)

                pool.close()
                pool.join()

                chunk_suffixes = results[0].get() 

                print(">>> Running fqfilter.py on chunks")
                max_reads = config.get('counting_opts', {}).get('max_reads', 0)
                
                processes = []
                for suffix in chunk_suffixes:
                    cmd = ['python3', resolve_script('fqfilter.py'), yaml_file, samtools, pigz, zumis_dir, suffix]
                    
                    if max_reads and int(max_reads) > 0:
                        chunk_limit = int(int(max_reads) / len(chunk_suffixes))
                        if chunk_limit < 1: chunk_limit = 1
                        cmd.extend(['--limit', str(chunk_limit)])
                        
                    processes.append(subprocess.Popen(cmd, stdout=run_log, stderr=subprocess.STDOUT, env=exec_env))

                for p in processes:
                    p.wait()
                    if p.returncode != 0:
                        raise RuntimeError(f"fqfilter failed (rc={p.returncode}). Check {log_path} for details.")

                print(">>> Cleaning up temporary FASTQ chunks...")
                import glob
                cleanup_candidates = (
                    glob.glob(os.path.join(tmp_merge_dir, "*.part_*"))
                    + glob.glob(os.path.join(tmp_merge_dir, "*.part_*.gz"))
                    + glob.glob(os.path.join(tmp_merge_dir, "*.fq.part_*"))
                    + glob.glob(os.path.join(tmp_merge_dir, "*.fq.part_*.gz"))
                    + glob.glob(os.path.join(tmp_merge_dir, "*.fastq.part_*"))
                    + glob.glob(os.path.join(tmp_merge_dir, "*.fastq.part_*.gz"))
                )
                for f in cleanup_candidates:
                    if not os.path.exists(f): continue
                    base = os.path.basename(f)
                    if base.endswith(".bam") or base.endswith(".bai") or base.endswith(".txt"): continue
                    if ".raw.tagged." in base or ".filtered.tagged." in base: continue
                    os.remove(f)

                print(">>> Merging BAM Stats")
                pipeline_modules.merge_bam_stats(tmp_merge_dir, project, analysis_dir, yaml_file, samtools)

                print(">>> Running Barcode Detection")
                run_stage_cmd(["python3", resolve_script("zUMIs_BCdetection.py"), yaml_file], "BCdetection")

                bc_bin_table = os.path.join(analysis_dir, 'zUMIs_output', f"{project}.BCbinning.txt")
                expect_id_barcode_file = os.path.join(out_dir, 'config', 'expect_id_barcode.tsv')
                
                if os.path.exists(bc_bin_table):
                    print(">>> Correcting BC Tags")
                    correct_processes = []
                    umi_chunks = []
                    int_chunks = []

                    for suffix in chunk_suffixes:
                        raw_bam = os.path.join(tmp_merge_dir, f"{project}{suffix}.raw.tagged.bam")
                        fixed_bam_umi = os.path.join(tmp_merge_dir, f"{project}{suffix}.filtered.tagged.umi.bam")
                        fixed_bam_int = os.path.join(tmp_merge_dir, f"{project}{suffix}.filtered.tagged.internal.bam")
                        
                        umi_chunks.append(fixed_bam_umi)
                        int_chunks.append(fixed_bam_int)

                        if os.path.exists(fixed_bam_umi): os.remove(fixed_bam_umi) 
                        if os.path.exists(fixed_bam_int): os.remove(fixed_bam_int)
                        
                        cmd_args = ['python3', resolve_script('correct_BCtag.py'), raw_bam, fixed_bam_umi, fixed_bam_int, bc_bin_table, expect_id_barcode_file]
                        correct_processes.append(subprocess.Popen(cmd_args, stdout=run_log, stderr=subprocess.STDOUT, env=exec_env))

                    for p in correct_processes:
                        p.wait()
                        if p.returncode != 0:
                            raise RuntimeError(f"correct_BCtag failed (rc={p.returncode}). Check {log_path} for details.")

                    for suffix in chunk_suffixes:
                        raw_bam = os.path.join(tmp_merge_dir, f"{project}{suffix}.raw.tagged.bam")
                        if os.path.exists(raw_bam): os.remove(raw_bam)
                            
                    print(">>> Skipping physical merge of chunks (will stream to STAR)...")
                    
            if which_stage in ["Filtering", "Mapping"]:
                print(">>> Starting Mapping Stage")
                
                umi_arg = ""
                int_arg = ""
                
                # Helper to find chunks if not in memory (e.g. restarting from Mapping)
                def find_chunks(suffix_pattern):
                    found = glob.glob(os.path.join(tmp_merge_dir, suffix_pattern))
                    return sorted(found)

                # Resolve UMI inputs
                if 'umi_chunks' in locals() and umi_chunks:
                    umi_arg = ",".join(umi_chunks)
                else:
                    # Try to find chunks on disk
                    disk_umi_chunks = find_chunks(f"{project}*.filtered.tagged.umi.bam")
                    if disk_umi_chunks:
                        print(f"Found {len(disk_umi_chunks)} UMI chunks on disk.")
                        umi_chunks = disk_umi_chunks # Update local var for cleanup later
                        umi_arg = ",".join(disk_umi_chunks)
                    else:
                        legacy_umi = os.path.join(analysis_dir, f"{project}.filtered.tagged.umi.unmapped.bam")
                        if os.path.exists(legacy_umi): 
                            umi_arg = legacy_umi
                        else:
                            # If we are starting at Mapping, we expect inputs.
                            if which_stage == "Mapping":
                                raise FileNotFoundError(f"Could not find input BAMs for Mapping stage. Checked for chunks in {tmp_merge_dir} and merged file {legacy_umi}")

                # Resolve Internal inputs
                if 'int_chunks' in locals() and int_chunks:
                    int_arg = ",".join(int_chunks)
                else:
                    disk_int_chunks = find_chunks(f"{project}*.filtered.tagged.internal.bam")
                    if disk_int_chunks:
                        print(f"Found {len(disk_int_chunks)} Internal chunks on disk.")
                        int_chunks = disk_int_chunks # Update local var for cleanup later
                        int_arg = ",".join(disk_int_chunks)
                    else:
                        legacy_int = os.path.join(analysis_dir, f"{project}.filtered.tagged.internal.unmapped.bam")
                        if os.path.exists(legacy_int): 
                            int_arg = legacy_int

                map_cmd = ['python3', resolve_script('mapping_analysis.py'), yaml_file, '--umi_bam', umi_arg, '--internal_bam', int_arg]
                expect_id_file = os.path.join(out_dir, 'config', 'expect_id_barcode.tsv')
                map_cmd.extend(['--expect_id_file', expect_id_file])
                run_stage_cmd(map_cmd, "mapping_analysis.py")
                                
                if 'umi_chunks' in locals() and umi_chunks:
                     for f in umi_chunks:
                         if os.path.exists(f): os.remove(f)
                
                if 'int_chunks' in locals() and int_chunks:
                     for f in int_chunks:
                         if os.path.exists(f): os.remove(f)

            if which_stage in ["Filtering", "Mapping", "Counting"]:
                print(">>> Starting Counting Stage")
                
                umi_aligned = os.path.join(analysis_dir, f"{project}.filtered.tagged.umi.Aligned.out.bam")
                int_aligned = os.path.join(analysis_dir, f"{project}.filtered.tagged.internal.Aligned.out.bam")
                umi_to_tx = os.path.join(analysis_dir, f"{project}.filtered.tagged.umi.Aligned.toTranscriptome.out.bam")
                int_to_tx = os.path.join(analysis_dir, f"{project}.filtered.tagged.internal.Aligned.toTranscriptome.out.bam")
                
                featurecounts_cmd = ['python3', resolve_script('run_featurecounts.py'), yaml_file, '--umi_bam', umi_aligned, '--internal_bam', int_aligned]
                run_stage_cmd(featurecounts_cmd, "FeatureCounts (Python)")

                remove_path(umi_aligned)
                remove_path(int_aligned)
                remove_path(umi_to_tx)
                remove_path(int_to_tx)

                print(">>> Starting DGE Analysis (Python)")
                dge_cmd = ['python3', resolve_script('dge_analysis.py'), yaml_file, samtools]
                run_stage_cmd(dge_cmd, "dge_analysis.py")

                gene_tagged_bam = os.path.join(analysis_dir, f"{project}.filtered.Aligned.GeneTagged.bam")
                stats_enabled = str(config.get('make_stats', 'yes')).lower() in ['yes', 'true']
                if not stats_enabled:
                    remove_path(gene_tagged_bam)

            if which_stage in ["Filtering", "Mapping", "Counting", "Summarising"]:
                if str(config.get('make_stats', 'yes')).lower() in ['yes', 'true']:
                    print(">>> Starting Statistics Stage")
                    stats_cmd = ['python3', resolve_script('generate_stats.py'), yaml_file]
                    run_stage_cmd(stats_cmd, "Stats (Python)")
                    gene_tagged_bam = os.path.join(analysis_dir, f"{project}.filtered.Aligned.GeneTagged.bam")
                    remove_path(gene_tagged_bam)

            print("Pipeline Finished Successfully.")
        finally:
            sys.stdout = original_stdout

def main():
    parser = argparse.ArgumentParser(description='Mhsflt Data Analysis Pipeline')
    parser.add_argument('--fastq', nargs='+', required=True, help='Input FASTQ files (R1 R2)')
    parser.add_argument('--genomeDir', required=True, help='Directory containing STAR index and GTF file')
    parser.add_argument('--sample', required=True, help='Sample name')
    parser.add_argument('--outdir', help='Output directory (default: ./<sample_name>)')
    parser.add_argument('--threads', type=int, default=20, help='Number of threads')
    parser.add_argument('--stage', choices=['Filtering', 'Mapping', 'Counting', 'Summarising'], default='Filtering', help='Analysis stage to start from')
    
    # Mutually exclusive group for sample mode
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument('--manual', help='Manual sample IDs (comma separated, e.g. "20,21"). Sets sample_type=manual.')
    mode_group.add_argument('--plate', help='Plate ID (e.g. "1"). Sets sample_type=auto.')
    mode_group.add_argument('--expectBarcode', help='Path to custom barcode file. Sets sample_type=custom.')

    args = parser.parse_args()

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Start analysis for {args.sample}.', flush=True)

    # Load default config
    config = copy.deepcopy(constant.DEFAULT_CONFIG)

    # Populate Config
    config['project'] = args.sample
    config['num_threads'] = args.threads
    config['zUMIs_directory'] = os.path.dirname(os.path.abspath(__file__))
    config['which_Stage'] = args.stage

    # Set Output Directory
    if args.outdir:
        root_out = os.path.abspath(args.outdir)
    else:
        root_out = os.path.join(os.getcwd(), args.sample)
        
    # Standardize structure: 
    # root_out/XPRESS_PROCESSING (Pipeline work dir)
    # root_out/outs (Final outputs)
    config['out_dir'] = os.path.join(root_out, 'XPRESS_PROCESSING')

    # Sample Type & ID
    if args.manual:
        config['sample']['sample_type'] = 'manual'
        config['sample']['sample_id'] = args.manual
    elif args.plate:
        config['sample']['sample_type'] = 'auto'
        config['sample']['sample_id'] = args.plate
    elif args.expectBarcode:
        config['sample']['sample_type'] = 'custom'
        config['sample']['sample_id'] = '1' # Dummy ID for custom
        if not os.path.exists(args.expectBarcode):
             raise FileNotFoundError(f"Custom barcode file not found: {args.expectBarcode}")
        config['barcodes']['barcode_file'] = os.path.abspath(args.expectBarcode)

    # Sequence Files & Read Length Detection
    if len(args.fastq) < 2:
        raise ValueError("At least 2 FASTQ files (R1 and R2) are required for PE analysis.")
    
    r1_file = os.path.abspath(args.fastq[0])
    r2_file = os.path.abspath(args.fastq[1])
    
    config['sequence_files']['file1']['name'] = r1_file
    config['sequence_files']['file2']['name'] = r2_file

    print(f"Detecting read lengths...")
    len_r1 = get_read_length(r1_file)
    len_r2 = get_read_length(r2_file)
    print(f"Detected R1 Length: {len_r1}, R2 Length: {len_r2}")
    
    # 1-based indexing for YAML
    umi_len = 10
    bc_len = 20
    r2_cdna_end = len_r2 - bc_len
    r2_bc_start = len_r2 - bc_len + 1
    
    if len_r1 <= umi_len:
        raise ValueError(f"R1 length ({len_r1}) must be > UMI length ({umi_len}).")
    config['sequence_files']['file1']['base_definition'] = [
        f"cDNA({umi_len + 1}-{len_r1})",
        f"UMI(1-{umi_len})"
    ]
    config['sequence_files']['file2']['base_definition'] = [
        f"cDNA(1-{r2_cdna_end})",
        f"BC({r2_bc_start}-{len_r2})"
    ]

    # Reference Files
    config['reference']['STAR_index'] = os.path.join(os.path.abspath(args.genomeDir), "star")
    
    # Auto-detect GTF in genomeDir
    gtf_files = glob.glob(os.path.join(args.genomeDir, "genes", "genes.gtf"))
    if not gtf_files:
        # Fallback to search in root if not in genes/
        gtf_files = glob.glob(os.path.join(args.genomeDir, "*.gtf"))
        
    if not gtf_files:
        raise FileNotFoundError(f"No GTF file found in {args.genomeDir} (checked 'genes/genes.gtf' and '*.gtf')")
    if len(gtf_files) > 1:
        print(f"Warning: Multiple GTF files found. Using: {gtf_files[0]}")
    config['reference']['GTF_file'] = os.path.abspath(gtf_files[0])

    # Setup Directories
    make_dir(config)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Directories created.', flush=True)

    # Validate Files
    check_file_exists(config)

    # Process Fastq (Link to data dir)
    process_fq(config)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Fastq processed.', flush=True)

    # Create Barcode File
    create_barcode(config)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Barcode files created.', flush=True)
    
    # We will replicate this behavior for the yaml passed to the pipeline stages
    run_config = copy.deepcopy(config)
    # run_config['out_dir'] was set to root_out which is correct.
    # make_dir created subdirectories inside root_out.
    # The pipeline scripts expect 'out_dir' to be the base where zUMIs_output etc live.
    
    # But we can check for software dir.
    software_dir = os.path.join(config['zUMIs_directory'], 'software')
    def resolve_tool(name):
        candidate = os.path.join(software_dir, name)
        if os.path.exists(candidate) and os.access(candidate, os.X_OK):
            return os.path.abspath(candidate)
        return name

    run_config['samtools_exec'] = resolve_tool('samtools')
    run_config['pigz_exec'] = resolve_tool('pigz')
    run_config['seqkit_exec'] = resolve_tool('seqkit')
    run_config['STAR_exec'] = resolve_tool('STAR')
    run_config['featureCounts_exec'] = resolve_tool('featureCounts')

    # Yaml Dumping Customizations
    class ForceStr:
        def __init__(self, value):
            self.value = value
    def force_str_representer(dumper, data):
        return dumper.represent_scalar('tag:yaml.org,2002:str', data.value, style='"')
    class ZumisDumper(yaml.SafeDumper): pass
    def bool_representer(dumper, value):
        return dumper.represent_scalar('tag:yaml.org,2002:bool', 'yes' if value else 'no')
    def none_representer(dumper, _value):
        return dumper.represent_scalar('tag:yaml.org,2002:null', '~')

    yaml.add_representer(ForceStr, force_str_representer, Dumper=ZumisDumper)
    yaml.add_representer(bool, bool_representer, Dumper=ZumisDumper)
    yaml.add_representer(type(None), none_representer, Dumper=ZumisDumper)

    # Ensure downsampling is string
    ds = run_config['counting_opts'].get('downsampling', '0')
    if isinstance(ds, list): ds = ",".join(map(str, ds))
    else: ds = str(ds)
    run_config['counting_opts']['downsampling'] = ForceStr(ds)

    final_yaml_path = os.path.join(config['out_dir'], 'config', 'run_config.yaml')
    os.makedirs(os.path.dirname(final_yaml_path), exist_ok=True)
    with open(final_yaml_path, 'w') as f:
        yaml.dump(run_config, f, Dumper=ZumisDumper, default_flow_style=False, sort_keys=False)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Config generated: {final_yaml_path}', flush=True)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Starting Pipeline...', flush=True)

    run_pipeline_stages(final_yaml_path)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} All analysis finished.', flush=True)

if __name__ == '__main__':
    main()
