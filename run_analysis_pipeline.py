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

def load_yaml(yaml_file):
    if not os.path.exists(yaml_file):
        raise FileNotFoundError(f'YAML file not exist: {yaml_file}')
    
    with open(yaml_file, 'r', encoding='utf-8') as f:
        data = yaml.safe_load(f)

    return data


def check_nested_empty(data, path=''):
    if isinstance(data, dict):
        for k, v in data.items():
            check_nested_empty(v, f"{path}.{k}" if path else k)
    elif isinstance(data, (list, tuple)):
        for i, v in enumerate(data):
            check_nested_empty(v, f"{path}[{i}]")
    elif data is None:
        if path not in ['reference.additional_files', 'barcodes.barcode_num', 'reference.additional_STAR_params', 'barcodes.barcode_file']:
            raise ValueError(f'Find empty value: {path}')


def check_species(data):
    species=data['sample']['sample_species'].lower()
    if species not in ['human', 'mouse']:
        raise ValueError(f'Species now only support "human" and "mouse".')
    

def make_dir(data):
    out_path=data['out_dir']
    if os.path.exists(out_path):
        print(f"Warning: Output directory '{out_path}' already exists. Resuming/Overwriting analysis.")
    else:
        os.makedirs(out_path)
    
    for dirs in ['data','config','analysis','results']:
        os.makedirs(f'{out_path}/{dirs}', exist_ok=True) 
    
    # Create zUMIs output dirs here as well
    for dirs in ['zUMIs_output', 'zUMIs_output/expression/downsampling', 'zUMIs_output/stats', 'zUMIs_output/.tmpMerge']:
        os.makedirs(f'{out_path}/analysis/{dirs}', exist_ok=True)


def validate_type_id(data):
    sample_type = data['sample']['sample_type'].lower()
    
    if sample_type not in ['manual', 'auto', 'custom', 'external']:
        raise ValueError('Wrong sample type in yaml file, must be "manual", "auto" or "custom".')
    
    if sample_type == 'custom' or sample_type == 'external':
        if 'barcode_file' not in data['barcodes'] or not data['barcodes']['barcode_file']:
             raise ValueError('For "custom" sample type, "barcode_file" must be specified in "barcodes" section.')
        if not os.path.exists(data['barcodes']['barcode_file']):
             raise FileNotFoundError(f"Custom barcode file not found: {data['barcodes']['barcode_file']}")
        return # No ID validation needed for custom

    if data['sample']['sample_id'] == '':      
        raise ValueError('sample_id must be set in yaml file for manual/auto modes')

    if sample_type == 'manual':
        s_id = str(data['sample']['sample_id']).split(',')
        for ind in s_id:
            ind = ind.strip() # Remove whitespace
            if not ind.isdigit() or int(ind) < 1 or int(ind) > 24:
                raise ValueError('Manual version sample id must be set within the range of 1-24, with multiple ids separated by commas.')                
    else:
        # Auto mode
        s_id = str(data['sample']['sample_id']).split(',')
        for ind in s_id:
            ind = ind.strip()
            try: 
                i_val = int(ind)
                if i_val < 1 or i_val > 12:
                    raise ValueError('Automatic version sample ID must be set within the range 1-12.')
            except ValueError:
                raise TypeError(f'Automatic version sample ID must be integers (1-12), got "{ind}".')


def create_barcode(data):
    sample_type = data['sample']['sample_type'].lower()
    out_path = data['out_dir']
    script_path = data['zUMIs_directory']
    
    # Custom/External Mode
    if sample_type == 'custom' or sample_type == 'external':
        provided_bc = data['barcodes']['barcode_file']
        print(f"Using custom barcode file: {provided_bc}")
        dest_summary = f'{out_path}/config/expect_id_barcode.tsv'
        dest_pipe = f'{out_path}/config/expect_barcode.tsv'
        
        shutil.copy(provided_bc, dest_summary)
        
        # Extract just barcodes for pipe file
        with open(provided_bc, 'r') as infile, open(dest_pipe, 'w') as outfile:
            for line in infile:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                     if parts[0].lower() == 'wellid': continue
                     outfile.write('\t'.join(parts[1:]) + '\n')
        return

    sample_id_str = str(data['sample']['sample_id'])
    sample_ids = [s.strip() for s in sample_id_str.split(',')]

    with open(f'{out_path}/config/expect_barcode.tsv','w') as pipe_file, \
         open(f'{out_path}/config/expect_id_barcode.tsv','w') as summary_file:
        
        print('\t'.join(['wellID','umi_barcodes','internal_barcodes']), file=summary_file)

        if sample_type == 'manual':
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
            # Auto Mode (Multi-plate support)
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
    
    
def check_file_exists(data):
    files_to_check=[data['reference']['STAR_index'], data['reference']['GTF_file']]
    
    # Check sequence files (comma separated)
    def check_seq_entry(entry):
        if not entry: return
        for f in entry.split(','):
            if not os.path.exists(f.strip()):
                raise FileNotFoundError(f'{f} not exist, please check!')

    check_seq_entry(data['sequence_files']['file1']['name'])
    check_seq_entry(data['sequence_files']['file2']['name'])
           

def process_fq(data):
    sample_type = data['sample']['sample_type'].lower()
    out_path = data['out_dir']
    
    def handle_file_entry(entry_name):
        if not entry_name: return ""
        # Handle comma-separated list
        src_files = [f.strip() for f in entry_name.split(',')]
        dest_names = []
        
        for src in src_files:
            base = os.path.basename(src)
            dest = f'{out_path}/data/{base}'
            
            if os.path.lexists(dest):
                try:
                    if os.path.islink(dest) and os.readlink(dest) == src:
                        pass
                    else:
                        os.remove(dest)
                        os.symlink(src, dest)
                except OSError:
                    pass # Race condition or permission
            else:
                os.symlink(src, dest)
            dest_names.append(base)
            
        return ",".join(dest_names)

    fq1_in = data['sequence_files']['file1']['name']
    fq2_in = data['sequence_files']['file2']['name']
    
    fq1_out = handle_file_entry(fq1_in)
    fq2_out = handle_file_entry(fq2_in)

    return fq1_out, fq2_out


def modify_yaml(data, fq1_names, fq2_names):
    out_path=data['out_dir']
    
    # fq1_names is comma-separated basenames. Prepend path.
    def prepend_path(names):
        return ",".join([f'{out_path}/data/{n}' for n in names.split(',')])
        
    data['sequence_files']['file1']['name'] = prepend_path(fq1_names)
    data['sequence_files']['file2']['name'] = prepend_path(fq2_names)

    # Unified barcode file
    data['barcodes']['barcode_file'] = f'{out_path}/config/expect_barcode.tsv'
    
    new_data=copy.deepcopy(data)
    del new_data['sample']
    new_data['out_dir'] = f'{out_path}/analysis'

    # Detect and enforce software paths from 'software' directory if present
    software_dir = os.path.join(data['zUMIs_directory'], 'software')
    
    def resolve_tool(tool_name, config_key):
        # 1. Check software dir
        candidate = os.path.join(software_dir, tool_name)
        if os.path.exists(candidate) and os.access(candidate, os.X_OK):
            return os.path.abspath(candidate)
        # 2. Keep existing config if set
        if config_key in new_data:
            return new_data[config_key]
        # 3. Default to name
        return tool_name

    new_data['samtools_exec'] = resolve_tool('samtools', 'samtools_exec')
    new_data['pigz_exec'] = resolve_tool('pigz', 'pigz_exec')
    new_data['seqkit_exec'] = resolve_tool('seqkit', 'seqkit_exec')
    new_data['STAR_exec'] = resolve_tool('STAR', 'STAR_exec')
    new_data['featureCounts_exec'] = resolve_tool('featureCounts', 'featureCounts_exec')
    
    # Rscript usually not in local software dir, but check if needed. Keeping as is for now.
    if 'Rscript_exec' in new_data:
        # Optional: could also verify Rscript
        pass

    if 'zUMIs_directory' not in new_data: new_data['zUMIs_directory'] = data['zUMIs_directory']

    class ForceStr:
        def __init__(self, value):
            self.value = value

    def force_str_representer(dumper, data):
        return dumper.represent_scalar('tag:yaml.org,2002:str', data.value, style='"')

    class ZumisDumper(yaml.SafeDumper):
        pass

    def bool_representer(dumper, value):
        return dumper.represent_scalar('tag:yaml.org,2002:bool', 'yes' if value else 'no')

    def none_representer(dumper, _value):
        return dumper.represent_scalar('tag:yaml.org,2002:null', '~')

    yaml.add_representer(ForceStr, force_str_representer, Dumper=ZumisDumper)
    yaml.add_representer(bool, bool_representer, Dumper=ZumisDumper)
    yaml.add_representer(type(None), none_representer, Dumper=ZumisDumper)

    # Ensure downsampling is a string
    ds_val = new_data['counting_opts'].get('downsampling', '0')
    if isinstance(ds_val, list):
        ds_val = ",".join(map(str, ds_val))
    else:
        ds_val = str(ds_val)
    
    new_data['counting_opts']['downsampling'] = ForceStr(ds_val)

    config_path = f'{out_path}/config/final_config.yaml'
    with open(config_path, 'w') as out_yaml:
        yaml.dump(new_data, out_yaml, Dumper=ZumisDumper, default_flow_style=False, sort_keys=False)
    
    return config_path


def run_pipeline_stages(yaml_file):
    """
    Orchestrates the pipeline stages previously handled by zUMIs.sh.
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

    log_path = os.path.join(out_dir, 'zUMIs_run.log')
    tmp_merge_dir = os.path.join(out_dir, 'zUMIs_output', '.tmpMerge')
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
                
                # Get File Lists (Handle comma-separated)
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
                        # Estimate uncompressed size (approx 3x compressed) for fast calculation
                        total_size_bytes += os.path.getsize(f) * 3
                    else:
                        total_size_bytes += os.path.getsize(f)

                # Estimate avg line len from first file
                avg_line_len = estimate_avg_line_len(first_fq, sample_lines=1000)
                if avg_line_len <= 0:
                    raise ValueError(f"Failed to estimate average line length for {first_fq}")

                total_lines_est = total_size_bytes / avg_line_len
                
                # Calculate strict lines_per_chunk for R1/R2 sync
                lines_per_chunk = int(math.ceil(total_lines_est / num_threads))
                # Multiple of 4
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
                    if not os.path.exists(f):
                        continue
                    base = os.path.basename(f)
                    if base.endswith(".bam") or base.endswith(".bai") or base.endswith(".txt"):
                        continue
                    if ".raw.tagged." in base or ".filtered.tagged." in base:
                        continue
                    os.remove(f)

                print(">>> Merging BAM Stats")
                pipeline_modules.merge_bam_stats(tmp_merge_dir, project, out_dir, yaml_file, samtools)

                print(">>> Running Barcode Detection")
                run_stage_cmd(["python3", resolve_script("zUMIs_BCdetection.py"), yaml_file], "BCdetection")

                bc_bin_table = os.path.join(out_dir, 'zUMIs_output', f"{project}.BCbinning.txt")
                expect_id_barcode_file = os.path.join(out_dir, '../config', 'expect_id_barcode.tsv')
                
                if os.path.exists(bc_bin_table):
                    print(">>> Correcting BC Tags")
                    correct_processes = []
                    
                    umi_chunks = []
                    int_chunks = []

                    for suffix in chunk_suffixes:
                        raw_bam = os.path.join(tmp_merge_dir, f"{project}{suffix}.raw.tagged.bam")
                        # Output files
                        fixed_bam_umi = os.path.join(tmp_merge_dir, f"{project}{suffix}.filtered.tagged.umi.bam")
                        fixed_bam_int = os.path.join(tmp_merge_dir, f"{project}{suffix}.filtered.tagged.internal.bam")
                        
                        umi_chunks.append(fixed_bam_umi)
                        int_chunks.append(fixed_bam_int)

                        # If fixed exists, remove it to ensure clean run
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
                        if os.path.exists(raw_bam):
                            os.remove(raw_bam)
                            
                    print(">>> Skipping physical merge of chunks (will stream to STAR)...")
                    
            if which_stage in ["Filtering", "Mapping"]:
                print(">>> Starting Mapping Stage")
                
                umi_arg = ""
                int_arg = ""
                
                # Check if we have chunks from the Filtering step
                if 'umi_chunks' in locals() and umi_chunks:
                    umi_arg = ",".join(umi_chunks)
                else:
                    # Fallback: check for legacy merged file
                    legacy_umi = os.path.join(out_dir, f"{project}.filtered.tagged.umi.unmapped.bam")
                    if os.path.exists(legacy_umi):
                        umi_arg = legacy_umi
                    else:
                        pass # Rely on Mapping script handling (or fail there)

                if 'int_chunks' in locals() and int_chunks:
                    int_arg = ",".join(int_chunks)
                else:
                    legacy_int = os.path.join(out_dir, f"{project}.filtered.tagged.internal.unmapped.bam")
                    if os.path.exists(legacy_int):
                        int_arg = legacy_int

                map_cmd = ['python3', resolve_script('mapping_analysis.py'), yaml_file, '--umi_bam', umi_arg, '--internal_bam', int_arg]
                run_stage_cmd(map_cmd, "mapping_analysis.py")
                                
                # Cleanup chunks after successful mapping
                if 'umi_chunks' in locals() and umi_chunks:
                     print(">>> Cleaning up UMI chunks...")
                     for f in umi_chunks:
                         if os.path.exists(f): os.remove(f)
                
                if 'int_chunks' in locals() and int_chunks:
                     print(">>> Cleaning up Internal chunks...")
                     for f in int_chunks:
                         if os.path.exists(f): os.remove(f)

            if which_stage in ["Filtering", "Mapping", "Counting"]:
                print(">>> Starting Counting Stage")
                
                # Updated to pass aligned BAMs (names determined by mapping_analysis output)
                umi_aligned = os.path.join(out_dir, f"{project}.filtered.tagged.umi.Aligned.out.bam")
                int_aligned = os.path.join(out_dir, f"{project}.filtered.tagged.internal.Aligned.out.bam")
                umi_to_tx = os.path.join(out_dir, f"{project}.filtered.tagged.umi.Aligned.toTranscriptome.out.bam")
                int_to_tx = os.path.join(out_dir, f"{project}.filtered.tagged.internal.Aligned.toTranscriptome.out.bam")
                
                featurecounts_cmd = ['python3', resolve_script('run_featurecounts.py'), yaml_file, '--umi_bam', umi_aligned, '--internal_bam', int_aligned]
                run_stage_cmd(featurecounts_cmd, "FeatureCounts (Python)")

                remove_path(umi_aligned)
                remove_path(int_aligned)
                remove_path(umi_to_tx)
                remove_path(int_to_tx)

                print(">>> Starting DGE Analysis (Python)")
                dge_cmd = ['python3', resolve_script('dge_analysis.py'), yaml_file, samtools]
                run_stage_cmd(dge_cmd, "dge_analysis.py")

                gene_tagged_bam = os.path.join(out_dir, f"{project}.filtered.Aligned.GeneTagged.bam")
                stats_enabled = str(config.get('make_stats', 'yes')).lower() in ['yes', 'true']
                if not stats_enabled:
                    remove_path(gene_tagged_bam)

            if which_stage in ["Filtering", "Mapping", "Counting", "Summarising"]:
                # Force stats to run by default if not explicitly disabled
                if str(config.get('make_stats', 'yes')).lower() in ['yes', 'true']:
                    print(">>> Starting Statistics Stage")
                    stats_cmd = ['python3', resolve_script('generate_stats.py'), yaml_file]
                    run_stage_cmd(stats_cmd, "Stats (Python)")
                    gene_tagged_bam = os.path.join(out_dir, f"{project}.filtered.Aligned.GeneTagged.bam")
                    remove_path(gene_tagged_bam)

            print("Pipeline Finished Successfully.")
        finally:
            sys.stdout = original_stdout

def get_args():
    parser = argparse.ArgumentParser(description='MGIEasy high sensitive full length transcriptome data analysis pipeline.')
    parser.add_argument('-y', '--yaml', required=True, type=str, help='Path of yaml file.')
    args = parser.parse_args()
    return args


def main():

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Start analysis.', flush=True)
    args = get_args()
    
    data=load_yaml(args.yaml)    
    check_nested_empty(data)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} YAML file load and check complete, all required fields are not empty.', flush=True)

    check_species(data)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Species is support.', flush=True)

    make_dir(data)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Create output directory finish.', flush=True)

    check_file_exists(data)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Fastq, reference and annotation files all exist.', flush=True)

    validate_type_id(data)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Check sample_type and sample_id are matched successfully.', flush=True)

    create_barcode(data)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Create barcode file finsh.', flush=True)

    fq1_names, fq2_names = process_fq(data)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Process fastq finish.', flush=True)

    # This now returns the path to the final config
    final_config_path = modify_yaml(data, fq1_names, fq2_names)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Create new yaml finish.', flush=True)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Start run zUMIs Pipeline.', flush=True)
    
    # Run the new Python pipeline logic
    run_pipeline_stages(final_config_path)
    
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Finish run zUMIs Pipeline.', flush=True)

    print('All analysis finish, bye.')

    
if __name__ == '__main__':
    main()
