#!/usr/bin/env python3
#-*-coding:utf-8-*- 

import os
import argparse
import yaml
import re, copy
import sys
import subprocess
import shutil
import pipeline_modules
from datetime import datetime
import multiprocessing
import gzip
import struct

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
    files_to_check=[data['sequence_files']['file1']['name'],
                   data['sequence_files']['file2']['name'],
                   data['reference']['STAR_index'],
                   data['reference']['GTF_file']]
    for filepath in files_to_check:
        if not os.path.exists(filepath):
            raise FileNotFoundError(f'{filepath} not exist, please check!')
           

def process_fq(data):
    sample_type = data['sample']['sample_type'].lower()
    fq1 = data['sequence_files']['file1']['name']
    fq1_name = os.path.basename(fq1)
    fq2 = data['sequence_files']['file2']['name']
    fq2_name = os.path.basename(fq2) 
    out_path=data['out_dir']

    def ensure_symlink(src, dst):
        if os.path.lexists(dst):
            try:
                if os.path.islink(dst) and os.readlink(dst) == src:
                    return
            except OSError:
                pass
            os.remove(dst)
        os.symlink(src, dst)

    ensure_symlink(fq1, f'{out_path}/data/{fq1_name}')
    ensure_symlink(fq2, f'{out_path}/data/{fq2_name}')
    
    # Unified logic: always symlink, no preprocessing with transfer_barcode
    out_fq2_name = fq2_name

    return fq1_name, out_fq2_name


def modify_yaml(data, fq1_names, fq2_names):
    out_path=data['out_dir']
    
    data['sequence_files']['file1']['name'] = f'{out_path}/data/{fq1_names}'
    data['sequence_files']['file2']['name'] = f'{out_path}/data/{fq2_names}'

    # Unified barcode file
    data['barcodes']['barcode_file'] = f'{out_path}/config/expect_barcode.tsv'
    
    new_data=copy.deepcopy(data)
    del new_data['sample']
    new_data['out_dir'] = f'{out_path}/analysis'

    # Add default executables if not present
    if 'samtools_exec' not in new_data: new_data['samtools_exec'] = 'samtools'
    if 'pigz_exec' not in new_data: new_data['pigz_exec'] = 'pigz'
    if 'STAR_exec' not in new_data: new_data['STAR_exec'] = 'STAR'
    if 'Rscript_exec' not in new_data: new_data['Rscript_exec'] = 'Rscript'
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

    # Ensure downsampling is a string, even if parsed as list or number
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
    rscript = config.get('Rscript_exec', 'Rscript')
    zumis_dir = config.get('zUMIs_directory', '.')
    
    class Tee:
        def __init__(self, *streams):
            self._streams = streams

        def write(self, s):
            for stream in self._streams:
                stream.write(s)

        def flush(self):
            for stream in self._streams:
                stream.flush()

    def gzip_uncompressed_size(path):
        with open(path, 'rb') as f:
            f.seek(-4, os.SEEK_END)
            return struct.unpack('<I', f.read(4))[0]

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
                
                # print(f">>> Running {stage_name}: {cmd_str}") # Already printed by stage headers mostly
                
                res = subprocess.run(cmd, stdout=run_log, stderr=subprocess.STDOUT, shell=shell)
                if res.returncode != 0:
                    run_log.flush()
                    try:
                        with open(log_path, 'r') as lr:
                            print(f"\n[ERROR] {stage_name} failed (rc={res.returncode}). Last 30 lines of log ({log_path}):\n", file=sys.stderr)
                            print("".join(lr.readlines()[-30:]), file=sys.stderr)
                    except Exception:
                        pass
                    raise RuntimeError(f"{stage_name} failed with exit code {res.returncode}.")

            if which_stage == "Filtering":
                print(">>> Starting Filtering Stage")
                # ... (Filtering code remains mostly same, but fqfilter is manual Popen) ...
                
                file_paths = []
                for entry in sorted_sequence_file_entries(config.get('sequence_files', {})):
                    if isinstance(entry, dict) and entry.get('name'):
                        file_paths.append(entry['name'])

                if not file_paths:
                    raise ValueError("No sequence file paths found in YAML configuration.")

                first_fq = file_paths[0]
                print(f"Estimating read count from {first_fq}...")

                if first_fq.endswith('.gz'):
                    file_size = gzip_uncompressed_size(first_fq) or os.path.getsize(first_fq)
                else:
                    file_size = os.path.getsize(first_fq)

                avg_line_len = estimate_avg_line_len(first_fq, sample_lines=1000)
                if avg_line_len <= 0:
                    raise ValueError(f"Failed to estimate average line length for {first_fq}")

                total_lines_est = file_size / avg_line_len
                total_reads_est = max(1, total_lines_est / 4)
                print(f"Estimated reads: {int(total_reads_est)}")

                pool = multiprocessing.Pool(processes=min(len(file_paths), num_threads))
                results = []
                for fq in file_paths:
                    res = pool.apply_async(
                        pipeline_modules.split_fastq,
                        (fq, num_threads, total_reads_est, tmp_merge_dir, project, pigz),
                    )
                    results.append(res)

                pool.close()
                pool.join()

                chunk_suffixes = results[0].get()

                print(">>> Running fqfilter.py on chunks")
                
                # Check for max_reads limit in config
                max_reads = config.get('counting_opts', {}).get('max_reads', 0)
                if not max_reads:
                    max_reads = config.get('max_reads', 0) # Support top-level as well
                
                processes = []
                for suffix in chunk_suffixes:
                    cmd = ['python3', f'{zumis_dir}/fqfilter.py', yaml_file, samtools, rscript, pigz, zumis_dir, suffix]
                    
                    if max_reads and int(max_reads) > 0:
                        # If we have multiple chunks, we should probably divide the limit?
                        # Or just apply the limit to each chunk?
                        # If we split the file, each chunk has a fraction of reads.
                        # If we apply total limit to EACH chunk, we get N * limit total.
                        # But since we use round-robin, reads are distributed evenly.
                        # So if we want TOTAL 100k reads, and we have 10 chunks, we should limit each to 10k.
                        
                        chunk_limit = int(int(max_reads) / len(chunk_suffixes))
                        if chunk_limit < 1: chunk_limit = 1
                        cmd.extend(['--limit', str(chunk_limit)])
                        
                    processes.append(subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True))

                for p in processes:
                    _, stderr = p.communicate()
                    if p.returncode != 0:
                        raise RuntimeError(f"fqfilter failed (rc={p.returncode}): {stderr.strip()}")

                print(">>> Merging BAM Stats")
                pipeline_modules.merge_bam_stats(tmp_merge_dir, project, out_dir, yaml_file, samtools)

                print(">>> Running Barcode Detection")
                # pipeline_modules.run_shell_cmd([rscript, f"{zumis_dir}/zUMIs-BCdetection.R", yaml_file], "BCdetection", log_path)
                # Updated to use the faster Python implementation for BC detection
                run_stage_cmd(["python3", f"{zumis_dir}/zUMIs_BCdetection.py", yaml_file], "BCdetection")

                bc_bin_table = os.path.join(out_dir, 'zUMIs_output', f"{project}.BCbinning.txt")
                expect_id_barcode_file = os.path.join(out_dir, '../config', 'expect_id_barcode.tsv')
                
                if os.path.exists(bc_bin_table):
                    print(">>> Correcting BC Tags")
                    correct_processes = []
                    
                    umi_chunks = []
                    int_chunks = []

                    for suffix in chunk_suffixes:
                        raw_bam = os.path.join(tmp_merge_dir, f"{project}.{suffix}.raw.tagged.bam")
                        # Fixed BAMs now split
                        fixed_bam_umi = os.path.join(tmp_merge_dir, f"{project}.{suffix}.filtered.tagged.umi.bam")
                        fixed_bam_int = os.path.join(tmp_merge_dir, f"{project}.{suffix}.filtered.tagged.internal.bam")
                        
                        umi_chunks.append(fixed_bam_umi)
                        int_chunks.append(fixed_bam_int)

                        if os.path.exists(fixed_bam_umi): os.remove(fixed_bam_umi) # Cleanup/Rename logic from before was weird, just process raw
                        # Note: The previous logic renamed fixed->raw if fixed existed (rerun?).
                        # Assuming raw exists or was renamed from fixed.
                        # For safety, let's assume we proceed from raw.
                        
                        # Check if old single filtered exists and we are re-running? 
                        # Ignore complex re-run logic for now, stick to standard flow.
                        
                        cmd_args = ['python3', f'{zumis_dir}/correct_BCtag.py', raw_bam, fixed_bam_umi, fixed_bam_int, bc_bin_table, expect_id_barcode_file]

                        correct_processes.append(subprocess.Popen(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True))

                    for p in correct_processes:
                        _, stderr = p.communicate()
                        if p.returncode != 0:
                            raise RuntimeError(f"correct_BCtag failed (rc={p.returncode}): {stderr.strip()}")
                            
                    # Merge UMI Chunks
                    print(">>> Merging UMI and Internal BAM chunks...")
                    umi_unmapped = os.path.join(out_dir, f"{project}.filtered.tagged.umi.unmapped.bam")
                    int_unmapped = os.path.join(out_dir, f"{project}.filtered.tagged.internal.unmapped.bam")
                    
                    # Using samtools cat
                    subprocess.check_call([samtools, 'cat', '-o', umi_unmapped] + umi_chunks)
                    subprocess.check_call([samtools, 'cat', '-o', int_unmapped] + int_chunks)
                    
                    # Cleanup chunks
                    for f in umi_chunks + int_chunks:
                         if os.path.exists(f): os.remove(f)

            if which_stage in ["Filtering", "Mapping"]:
                print(">>> Starting Mapping Stage")
                # Updated to pass the two split BAMs
                umi_unmapped = os.path.join(out_dir, f"{project}.filtered.tagged.umi.unmapped.bam")
                int_unmapped = os.path.join(out_dir, f"{project}.filtered.tagged.internal.unmapped.bam")
                
                map_cmd = ['python3', f'{zumis_dir}/mapping_analysis.py', yaml_file, '--umi_bam', umi_unmapped, '--internal_bam', int_unmapped]
                run_stage_cmd(map_cmd, "mapping_analysis.py")

            if which_stage in ["Filtering", "Mapping", "Counting"]:
                print(">>> Starting Counting Stage")
                # Replaced R script with Python implementation
                
                # Updated to pass aligned BAMs (names determined by mapping_analysis output)
                # Mapping analysis will produce:
                umi_aligned = os.path.join(out_dir, f"{project}.filtered.tagged.umi.Aligned.out.bam")
                int_aligned = os.path.join(out_dir, f"{project}.filtered.tagged.internal.Aligned.out.bam")
                
                featurecounts_cmd = ['python3', f'{zumis_dir}/run_featurecounts.py', yaml_file, '--umi_bam', umi_aligned, '--internal_bam', int_aligned]
                run_stage_cmd(featurecounts_cmd, "FeatureCounts (Python)")

                print(">>> Starting DGE Analysis (Python)")
                dge_cmd = ['python3', f'{zumis_dir}/dge_analysis.py', yaml_file, samtools]
                run_stage_cmd(dge_cmd, "dge_analysis.py")

            if which_stage in ["Filtering", "Mapping", "Counting", "Summarising"]:
                # Force stats to run by default if not explicitly disabled
                if str(config.get('make_stats', 'yes')).lower() in ['yes', 'true']:
                    print(">>> Starting Statistics Stage")
                    stats_cmd = ['python3', f'{zumis_dir}/generate_stats.py', yaml_file]
                    run_stage_cmd(stats_cmd, "Stats (Python)")

            print("Pipeline Finished Successfully.")
        finally:
            sys.stdout = original_stdout


def run_create_reports(data):
    sample_name = data['project']
    sample_species = data['sample']['sample_species'].upper()
    script_path = data['zUMIs_directory']
    result_dir = data['out_dir']

    summary_script = f'{script_path}/report/create_summary'
    if not os.path.exists(summary_script):
        print(f"Warning: Summary script not found at {summary_script}. Skipping report generation.")
        # Try to copy stats PDF if available even if report generation is skipped
        pdf_source = f'{result_dir}/analysis/zUMIs_output/stats/{sample_name}.features.pdf'
        if os.path.exists(pdf_source):
             shutil.copy(pdf_source, f'{result_dir}/results/{sample_name}.features.pdf')
        return

    if sample_species == 'HUMAN':
        sample_species = 'Human'
    else:
        sample_species = 'Mouse'
    
    cmd = [
        f'{script_path}/report/create_summary',
        '--sample', str(sample_name),
        '--indir', f'{result_dir}/analysis',
        '--species', str(sample_species),
        '--well', f'{result_dir}/config/expect_id_barcode.tsv',
    ]
    with open(f'{result_dir}/analysis/report_run.log','w') as run_log, open(f'{result_dir}/analysis/report_run_error.log','w') as error_log:
        process_status = subprocess.run(cmd, stdout=run_log, stderr=error_log, text=True)

    returncode = process_status.returncode
    if int(returncode) != 0:
        raise Exception(f'Error happen when cerate reports, see {result_dir}/analysis/report_run_error.log for detail.')
    
    shutil.copytree(f'{result_dir}/analysis/summary/div', f'{result_dir}/results/html')
    shutil.copy(f'{result_dir}/analysis/summary/{sample_name}_stat.xls', f'{result_dir}/results/{sample_name}_stat.xls')
    
    pdf_source = f'{result_dir}/analysis/zUMIs_output/stats/{sample_name}.features.pdf'
    if os.path.exists(pdf_source):
        shutil.copy(pdf_source, f'{result_dir}/results/{sample_name}.features.pdf')
    else:
        print(f"Warning: Features PDF not found at {pdf_source}, skipping copy.")


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

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Start create results.', flush=True)
    run_create_reports(data)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Finish create results.', flush=True)

    print('All analysis finish, bye.')

    
if __name__ == '__main__':
    main()
