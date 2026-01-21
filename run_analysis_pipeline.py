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
        raise OSError('Out dirctory in yaml is exist, please delete it')
    
    for dirs in ['data','config','analysis','results']:
        os.makedirs(f'{out_path}/{dirs}') 
    
    # Create zUMIs output dirs here as well
    for dirs in ['zUMIs_output', 'zUMIs_output/expression/downsampling', 'zUMIs_output/stats', 'zUMIs_output/.tmpMerge']:
        os.makedirs(f'{out_path}/analysis/{dirs}', exist_ok=True)


def validate_type_id(data):
    if data['sample']['sample_type'].lower() != 'manual' and data['sample']['sample_type'].lower() != 'auto':
        raise ValueError('Wrong sample type in yaml file, must be "manual" or "auto".')
    
    if data['sample']['sample_type'] == '' or data['sample']['sample_id'] == '':      
        raise ValueError('sample type and sample id must set in yaml file')

    if data['sample']['sample_type'].lower() == 'manual':
        s_id=str(data['sample']['sample_id']).split(',')
        for ind in s_id:
            if int(ind) < 1 or int(ind) > 24:
                raise ValueError('Manual version sample id must be set within the range of 1-24, with multiple ids separated by commas.')                
    else:
        try: 
            ind=int(data['sample']['sample_id'])
            if ind < 1 or ind > 12:
                raise ValueError('Automatic version sample ID must be set within the range 1-12.')
        except:
            raise TypeError(f'Automatic version sample ID only support single sample id.')


def create_barcode(data):
    sample_type=data['sample']['sample_type'].lower()
    sample_id=str(data['sample']['sample_id'])
    out_path=data['out_dir']
    script_path=data['zUMIs_directory']

    if sample_type=='manual':
        with open(f'{script_path}/manual_barcode_list.yaml', 'r', encoding='utf-8') as y:
            barcode_set = yaml.safe_load(y)
        sample_id=sample_id.split(',')
        with open(f'{out_path}/config/expect_barcode.tsv','w') as pipe_file, open(f'{out_path}/config/expect_bin_barcode.tsv','w') as bin_file, open(f'{out_path}/config/expect_id_barcode.tsv','w') as summary_file:
            print('\t'.join(['wellID','umi_barcodes','internal_barcodes']),file=summary_file)
            for c in sample_id:
                bc_t_i5=barcode_set[c][0];bc_n_i5=barcode_set[c][1];bc_n_i7=barcode_set[c][2]
                for k,v in zip(bc_t_i5,bc_n_i5):
                    for j in bc_n_i7:
                        print('\t'.join([k+j,v+j]),file=pipe_file)
                for kk,vv in zip(bc_t_i5,bc_n_i5):
                    for jj in bc_n_i7:
                        print('\t'.join([kk+jj,vv+jj]),file=bin_file)
                        print('\t'.join(['A'+c,kk+jj,vv+jj]),file=summary_file)
                        break
                    break
    else:
        with open(f'{script_path}/auto_barcode_list.yaml', 'r', encoding='utf-8') as y:
            barcode_set = yaml.safe_load(y)
        with open(f'{out_path}/config/expect_barcode.tsv','w') as pipe_file, open(f'{out_path}/config/expect_id_barcode.tsv','w') as summary_file:
            print('\t'.join(['wellID','umi_barcodes','internal_barcodes']),file=summary_file)
            for k in barcode_set[f'plate{sample_id}']:
                print('\t'.join(barcode_set[f'plate{sample_id}'][k]),file=pipe_file)
                print(k+'\t'+'\t'.join(barcode_set[f'plate{sample_id}'][k]),file=summary_file)    
    
    
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
    script_path = data['zUMIs_directory']
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

    if sample_type == 'manual':
        out_fq2_name = f'{fq2_name}.process.fq.gz'
        ham_dist=data['barcodes']['BarcodeBinning']
        threads=data['num_threads']
        print('Process fastq binning......', flush=True)
        cmd = [
            os.path.join(script_path, 'transfer_barcode'),
            '-i', fq2,
            '-l', f'{out_path}/config/expect_barcode.tsv',
            '-o', out_fq2_name,
            '-d', str(ham_dist),
            '-p', f'{out_path}/data',
            '-t', str(threads),
            '-e', script_path,
        ]
        process_status = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if process_status.returncode != 0:
            raise RuntimeError(f'Bin barcode error: {process_status.stderr}')
    else:
        ensure_symlink(fq2, f'{out_path}/data/{fq2_name}')
        out_fq2_name = fq2_name

    return fq1_name, out_fq2_name


def modify_yaml(data, fq1_names, fq2_names):
    sample_type = data['sample']['sample_type'].lower()
    out_path=data['out_dir']
    
    data['sequence_files']['file1']['name'] = f'{out_path}/data/{fq1_names}'
    data['sequence_files']['file2']['name'] = f'{out_path}/data/{fq2_names}'

    if sample_type == 'manual':
        data['barcodes']['barcode_file'] = f'{out_path}/config/expect_bin_barcode.tsv'
    else:
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
        return dumper.represent_scalar('tag:yaml.org,2002:str', data.value, style='\'')

    class ZumisDumper(yaml.SafeDumper):
        pass

    def bool_representer(dumper, value):
        return dumper.represent_scalar('tag:yaml.org,2002:bool', 'yes' if value else 'no')

    def none_representer(dumper, _value):
        return dumper.represent_scalar('tag:yaml.org,2002:null', '~')

    yaml.add_representer(ForceStr, force_str_representer, Dumper=ZumisDumper)
    yaml.add_representer(bool, bool_representer, Dumper=ZumisDumper)
    yaml.add_representer(type(None), none_representer, Dumper=ZumisDumper)

    new_data['counting_opts']['downsampling'] = ForceStr(new_data['counting_opts']['downsampling'])

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

            if which_stage == "Filtering":
                print(">>> Starting Filtering Stage")

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
                processes = []
                for suffix in chunk_suffixes:
                    cmd = ['python3', f'{zumis_dir}/fqfilter.py', yaml_file, samtools, rscript, pigz, zumis_dir, suffix]
                    processes.append(subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True))

                for p in processes:
                    _, stderr = p.communicate()
                    if p.returncode != 0:
                        raise RuntimeError(f"fqfilter failed (rc={p.returncode}): {stderr.strip()}")

                print(">>> Merging BAM Stats")
                pipeline_modules.merge_bam_stats(tmp_merge_dir, project, out_dir, yaml_file, samtools)

                print(">>> Running Barcode Detection")
                pipeline_modules.run_shell_cmd([rscript, f"{zumis_dir}/zUMIs-BCdetection.R", yaml_file], "BCdetection", log_path)

                bc_bin_table = os.path.join(out_dir, 'zUMIs_output', f"{project}.BCbinning.txt")
                bc_bin_raw_table = os.path.join(out_dir, 'zUMIs_output', f"{project}.BCbinning.raw.txt")
                if os.path.exists(bc_bin_table):
                    print(">>> Correcting BC Tags")
                    chemistry = config.get('chemistry', '')
                    correct_processes = []

                    for suffix in chunk_suffixes:
                        raw_bam = os.path.join(tmp_merge_dir, f"{project}.{suffix}.raw.tagged.bam")
                        fixed_bam = os.path.join(tmp_merge_dir, f"{project}.{suffix}.filtered.tagged.bam")

                        if os.path.exists(fixed_bam):
                            os.rename(fixed_bam, raw_bam)

                        cmd_args = ['python3', f'{zumis_dir}/correct_BCtag.py', raw_bam, fixed_bam, bc_bin_table, samtools]
                        if chemistry == 'MGI':
                            cmd_args = ['python3', f'{zumis_dir}/correct_BCtag.py', raw_bam, fixed_bam, bc_bin_table, samtools, bc_bin_raw_table]

                        correct_processes.append(subprocess.Popen(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True))

                    for p in correct_processes:
                        _, stderr = p.communicate()
                        if p.returncode != 0:
                            raise RuntimeError(f"correct_BCtag failed (rc={p.returncode}): {stderr.strip()}")

            if which_stage in ["Filtering", "Mapping"]:
                print(">>> Starting Mapping Stage")
                map_cmd = ['python3', f'{zumis_dir}/mapping_analysis.py', yaml_file]
                subprocess.run(map_cmd, stdout=run_log, stderr=subprocess.STDOUT, check=True)

            if which_stage in ["Filtering", "Mapping", "Counting"]:
                print(">>> Starting Counting Stage")
                pipeline_modules.run_shell_cmd([rscript, f"{zumis_dir}/zUMIs-dge2.R", yaml_file], "FeatureCounts (R)", log_path)

                print(">>> Starting DGE Analysis (Python)")
                dge_cmd = ['python3', f'{zumis_dir}/dge_analysis.py', yaml_file, samtools]
                subprocess.run(dge_cmd, stdout=run_log, stderr=subprocess.STDOUT, check=True)

                if config.get('velocyto', 'no') == 'yes':
                    pipeline_modules.run_shell_cmd([rscript, f"{zumis_dir}/runVelocyto.R", yaml_file], "Velocyto", log_path)

            if which_stage in ["Filtering", "Mapping", "Counting", "Summarising"]:
                if config.get('make_stats', 'no') == 'yes':
                    print(">>> Starting Statistics Stage")
                    pipeline_modules.run_shell_cmd([rscript, f"{zumis_dir}/zUMIs-stats2.R", yaml_file], "Stats", log_path)

            print("Pipeline Finished Successfully.")
        finally:
            sys.stdout = original_stdout


def run_create_reports(data):
    sample_name = data['project']
    sample_species = data['sample']['sample_species'].upper()
    script_path = data['zUMIs_directory']
    result_dir = data['out_dir']

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
    shutil.copy(f'{result_dir}/analysis/zUMIs_output/stats/{sample_name}.features.pdf', f'{result_dir}/results/{sample_name}.features.pdf')


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
