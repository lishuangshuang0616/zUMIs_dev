#!/usr/bin/env python3
import sys
import os
import glob
import subprocess
import yaml
import math
import shutil
import collections
import itertools

def load_config(yaml_file):
    with open(yaml_file, 'r') as f:
        return yaml.safe_load(f)

def run_cmd(cmd, shell=True):
    print(f"Running: {cmd}")
    subprocess.check_call(cmd, shell=shell)

def get_bam_read_length(bam_file, samtools, n_reads=1000):
    proc = subprocess.Popen(
        [samtools, 'view', bam_file],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
    )
    lengths = collections.Counter()
    try:
        for line in itertools.islice(proc.stdout, n_reads):
            parts = line.rstrip('\n').split('\t')
            if len(parts) > 9:
                lengths[len(parts[9])] += 1
    finally:
        try:
            proc.stdout.close()
        except Exception:
            pass
        proc.terminate()
        proc.wait()

    # rc=0: success
    # rc=-15: SIGTERM (we terminated it)
    # rc=-13: SIGPIPE (we closed the pipe while it was writing, expected)
    if proc.returncode not in (0, -15, -13):
        stderr = proc.stderr.read().strip() if proc.stderr else ''
        raise RuntimeError(f"samtools view failed (rc={proc.returncode}) for {bam_file}: {stderr}")

    if not lengths:
        return 0

    return lengths.most_common(1)[0][0]


def get_dir_size_gb(path):
    try:
        out = subprocess.check_output(['du', '-sk', path], text=True).split('\t', 1)[0]
        kb = float(out.strip())
        return kb / (1024 * 1024)
    except Exception:
        return 25.0

def setup_gtf(config, project, out_dir, samtools):
    # Handle additional files
    gtf = config['reference']['GTF_file']
    additional_files = config['reference'].get('additional_files', [])
    if not additional_files:
        final_gtf = os.path.join(out_dir, f"{project}.final_annot.gtf")
        shutil.copyfile(gtf, final_gtf)
        return final_gtf, ""
    
    # Process additional fasta
    add_gtf_path = os.path.join(out_dir, "additional_sequence_annot.gtf")
    with open(add_gtf_path, 'w') as out_f:
        for fa in additional_files:
            # Get lengths using samtools faidx
            # Assuming faidx exists? If not generate it.
            if not os.path.exists(fa + ".fai"):
                subprocess.run([samtools, "faidx", fa], check=True)
            
            with open(fa + ".fai") as fai:
                for line in fai:
                    parts = line.split('\t')
                    name = parts[0]
                    length = parts[1]
                    # Write custom GTF line
                    # R code: gene_id "name"; transcript_id "name"; ...
                    attr = f'gene_id "{name}"; transcript_id "{name}"; exon_number "1"; gene_name "{name}"; gene_biotype "User"; transcript_name "{name}"; exon_id "{name}"'
                    out_f.write(f"{name}\tUser\texon\t1\t{length}\t.\t+\t.\t{attr}\n")

    final_gtf = os.path.join(out_dir, f"{project}.final_annot.gtf")
    with open(final_gtf, 'w') as outfile:
        # Cat original GTF
        with open(gtf, 'r') as infile:
            shutil.copyfileobj(infile, outfile)
        # Cat additional
        with open(add_gtf_path, 'r') as infile:
            shutil.copyfileobj(infile, outfile)
            
    param_add = f"--genomeFastaFiles {' '.join(additional_files)}"
    return final_gtf, param_add

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('yaml_file')
    parser.add_argument('--umi_bam', required=False, help="Merged UMI Unmapped BAM")
    parser.add_argument('--internal_bam', required=False, help="Merged Internal Unmapped BAM")
    args = parser.parse_args()
        
    yaml_file = args.yaml_file
    config = load_config(yaml_file)
    
    project = config['project']
    out_dir = config['out_dir']
    num_threads = int(config['num_threads'])
    samtools = config.get('samtools_exec', 'samtools')
    star_exec = config.get('STAR_exec', 'STAR')
    star_index = config['reference']['STAR_index']
    
    # Inputs from args or fallback (though pipeline should provide them)
    umi_bam = args.umi_bam
    internal_bam = args.internal_bam
    
    if not umi_bam or not internal_bam:
        # Fallback to single file logic or error?
        # Given the refactor, let's assume pipeline is updated.
        # But if running manually, maybe check defaults.
        umi_bam = os.path.join(out_dir, f"{project}.filtered.tagged.umi.unmapped.bam")
        internal_bam = os.path.join(out_dir, f"{project}.filtered.tagged.internal.unmapped.bam")
    
    if not os.path.exists(umi_bam) and not os.path.exists(internal_bam):
        raise FileNotFoundError("Input unmapped BAMs not found.")

    # 2. Setup GTF
    final_gtf, param_add_fa = setup_gtf(config, project, out_dir, samtools)
    
    # 3. Determine Read Length (Use UMI bam as representative)
    read_len = 0
    if os.path.exists(umi_bam):
        read_len = get_bam_read_length(umi_bam, samtools)
    elif os.path.exists(internal_bam):
        read_len = get_bam_read_length(internal_bam, samtools)
        
    print(f"Detected Read Length: {read_len}")
    
    # 4. Resource Allocation
    # User requested: 3/4 threads for UMI, 1/4 for Internal
    # Reserve threads for SAMtools inside STAR? Usually --runThreadN is mapping threads.
    # STAR also spawns sorting/BAM writing threads.
    
    # Let's allocate raw threads to --runThreadN
    
    t_umi = max(1, int(num_threads * 0.75))
    t_int = max(1, num_threads - t_umi)
    
    print(f"Allocating threads: UMI={t_umi}, Internal={t_int} (Total {num_threads})")

    # 5. Build STAR Commands
    read_layout = config.get('read_layout', 'SE')
    
    defaults_umi = f"--readFilesCommand {samtools} view -@ {max(1, int(t_umi/4))} --outSAMmultNmax 1 --outFilterMultimapNmax 50 --outSAMunmapped Within --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM --limitOutSJcollapsed 5000000"
    defaults_int = f"--readFilesCommand {samtools} view -@ {max(1, int(t_int/4))} --outSAMmultNmax 1 --outFilterMultimapNmax 50 --outSAMunmapped Within --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM --limitOutSJcollapsed 5000000"
    
    # Note: Using separate genomeDir loading might cause memory issues if RAM is tight.
    # But parallel execution requires it unless using LoadAndKeep.
    # We will proceed with standard parallel execution.
    
    misc_umi = f"--genomeDir {star_index} --sjdbGTFfile {final_gtf} --runThreadN {t_umi} --sjdbOverhang {read_len-1} --readFilesType SAM {read_layout}"
    misc_int = f"--genomeDir {star_index} --sjdbGTFfile {final_gtf} --runThreadN {t_int} --sjdbOverhang {read_len-1} --readFilesType SAM {read_layout}"
    
    extra_params = config['reference'].get('additional_STAR_params', '') or ""
    
    star_cmd_umi_base = f"{star_exec} {defaults_umi} {misc_umi} {extra_params} {param_add_fa}"
    star_cmd_int_base = f"{star_exec} {defaults_int} {misc_int} {extra_params} {param_add_fa}"
    
    if config['counting_opts'].get('twoPass', False):
        star_cmd_umi_base += " --twopassMode Basic"
        star_cmd_int_base += " --twopassMode Basic"

    procs = []

    # Run UMI
    if os.path.exists(umi_bam):
        prefix_umi = os.path.join(out_dir, f"{project}.filtered.tagged.umi.")
        cmd_umi = f"{star_cmd_umi_base} --readFilesIn {umi_bam} --outFileNamePrefix {prefix_umi}"
        print("Starting STAR for UMI...")
        p_umi = subprocess.Popen(cmd_umi, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        procs.append(("UMI", p_umi))

    # Run Internal
    if os.path.exists(internal_bam):
        prefix_int = os.path.join(out_dir, f"{project}.filtered.tagged.internal.")
        cmd_int = f"{star_cmd_int_base} --readFilesIn {internal_bam} --outFileNamePrefix {prefix_int}"
        print("Starting STAR for Internal...")
        p_int = subprocess.Popen(cmd_int, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        procs.append(("Internal", p_int))
        
    # Wait
    failed = False
    for name, p in procs:
        _, stderr = p.communicate()
        if p.returncode != 0:
            print(f"STAR {name} failed (rc={p.returncode}): {stderr.strip()}")
            failed = True
        else:
            print(f"STAR {name} finished successfully.")
            
    if failed:
        raise RuntimeError("One or more STAR instances failed.")

if __name__ == "__main__":
    main()