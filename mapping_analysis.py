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
    if len(sys.argv) < 2:
        print("Usage: python3 mapping_analysis.py <yaml_config>")
        sys.exit(1)
        
    yaml_file = sys.argv[1]
    config = load_config(yaml_file)
    
    project = config['project']
    out_dir = config['out_dir']
    num_threads = int(config['num_threads'])
    samtools = config.get('samtools_exec', 'samtools')
    star_exec = config.get('STAR_exec', 'STAR')
    star_index = config['reference']['STAR_index']
    mem_limit = int(config.get('mem_limit', 0))
    if mem_limit == 0: mem_limit = 100 # default 100G in R script?
    
    which_stage = config['which_Stage']
    
    tmp_merge_dir = os.path.join(out_dir, "zUMIs_output/.tmpMerge")
    
    # 1. Collect BAMs
    filtered_bams = []
    if which_stage == "Filtering":
        filtered_bams = glob.glob(os.path.join(tmp_merge_dir, f"{project}.*.filtered.tagged.bam"))
        # Merge unmapped in parallel?
        # R script runs STAR and samtools cat in parallel.
    else:
        filtered_bams = [os.path.join(out_dir, f"{project}.filtered.tagged.unmapped.bam")]

    if not filtered_bams:
        msg = f"No BAM files found for mapping in {tmp_merge_dir} matching pattern {project}.*.filtered.tagged.bam"
        print(msg)
        raise FileNotFoundError(msg)

    # 2. Setup GTF
    final_gtf, param_add_fa = setup_gtf(config, project, out_dir, samtools)
    
    # 3. Determine Read Length
    read_len = get_bam_read_length(filtered_bams[0], samtools)
    print(f"Detected Read Length: {read_len}")
    
    if read_len <= 0:
        raise ValueError(f"Detected read length is {read_len}. Please check if the input BAM file {filtered_bams[0]} is empty or corrupted.")
    
    # 4. Calc STAR instances
    # R script logic: genome size from du -sh. 
    # Let's approximate.
    genome_size = get_dir_size_gb(star_index)
        
    num_instances = math.floor(mem_limit / genome_size)
    if num_instances < 1: num_instances = 1
    if num_instances > 5: num_instances = 5 # Limit to max 5 chunks per user request
    if num_instances > num_threads: num_instances = num_threads

    print(f"DEBUG: Genome Size: {genome_size:.2f} GB")
    print(f"DEBUG: Memory Limit: {mem_limit} GB")
    print(f"DEBUG: Calculated STAR instances: {num_instances} (Max allowed: 5)")
    print(f"DEBUG: Threads per instance: {num_threads} (total) -> {max(1, math.floor((num_threads - (2 if num_threads > 8 else 1)) / num_instances))} (per STAR)")
    
    samtools_cores = 2 if num_threads > 8 else 1
    avail_cores = num_threads - samtools_cores
    if which_stage == "Filtering":
        avail_cores = max(1, math.floor(avail_cores / num_instances))
    
    # 5. Build STAR Command
    read_layout = config.get('read_layout', 'SE')
    
    # Note: R script used --readFilesType SAM <read_layout> 
    # Usually valid values are 'SAM SE' or 'SAM PE' 
    
    defaults = f"--readFilesCommand {samtools} view -@{samtools_cores} --outSAMmultNmax 1 --outFilterMultimapNmax 50 --outSAMunmapped Within --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM --limitOutSJcollapsed 5000000"
    misc = f"--genomeDir {star_index} --sjdbGTFfile {final_gtf} --runThreadN {avail_cores} --sjdbOverhang {read_len-1} --readFilesType SAM {read_layout}"
    
    extra_params = config['reference'].get('additional_STAR_params', '')
    if extra_params is None: extra_params = ""
    
    star_cmd_base = f"{star_exec} {defaults} {misc} {extra_params} {param_add_fa}"
    
    if config['counting_opts'].get('twoPass', False):
        star_cmd_base += " --twopassMode Basic"

    # 6. Run
    # Prepare unmapped merge command if Filtering stage
    unmapped_bam = os.path.join(out_dir, f"{project}.filtered.tagged.unmapped.bam")
    sam_merge_cmd = [samtools, "cat", "-o", unmapped_bam] + filtered_bams
    
    if num_instances > 1 and which_stage == "Filtering":
        map_tmp_dir = os.path.join(out_dir, "zUMIs_output/.tmpMap")
        if not os.path.exists(map_tmp_dir): os.makedirs(map_tmp_dir)
        
        # Split inputs
        chunk_size = math.ceil(len(filtered_bams) / num_instances)
        # Split filtered_bams list into chunks
        chunks = [filtered_bams[i:i + chunk_size] for i in range(0, len(filtered_bams), chunk_size)]
        
        star_procs = []
        for i, chunk in enumerate(chunks):
            chunk_files = ",".join(chunk)
            prefix = os.path.join(map_tmp_dir, f"tmp.{project}.{i+1}.")
            
            cmd = f"{star_cmd_base} --readFilesIn {chunk_files} --outFileNamePrefix {prefix}"
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            star_procs.append(p)
            
        # Run sam merge in parallel
        sam_proc = subprocess.Popen(sam_merge_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        # Wait all
        for p in star_procs:
            _, stderr = p.communicate()
            if p.returncode != 0:
                raise RuntimeError(f"STAR failed (rc={p.returncode}): {stderr.strip()}")
        _, sam_stderr = sam_proc.communicate()
        if sam_proc.returncode != 0:
            raise RuntimeError(f"samtools cat failed (rc={sam_proc.returncode}): {sam_stderr.strip()}")
        
        # Merge outputs
        # SJs
        sjs = glob.glob(os.path.join(map_tmp_dir, "*.SJ.out.tab"))
        if sjs:
            for p in sjs:
                shutil.copy2(p, out_dir)
            
        # Logs
        logs = glob.glob(os.path.join(map_tmp_dir, "*.Log.final.out"))
        if logs:
            out_log = os.path.join(out_dir, f"{project}.filtered.tagged.Log.final.out")
            with open(out_log, "w") as out_f:
                for p in logs:
                    with open(p, "r") as in_f:
                        shutil.copyfileobj(in_f, out_f)
            
        # BAMs
        bams = glob.glob(os.path.join(map_tmp_dir, "*.Aligned.out.bam"))
        tx_bams = glob.glob(os.path.join(map_tmp_dir, "*.Aligned.toTranscriptome.out.bam"))
        
        final_bam = os.path.join(out_dir, f"{project}.filtered.tagged.Aligned.out.bam")
        final_tx_bam = os.path.join(out_dir, f"{project}.filtered.tagged.Aligned.toTranscriptome.out.bam")
        
        p1 = subprocess.Popen([samtools, "cat", "-o", final_bam] + bams, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        p2 = subprocess.Popen([samtools, "cat", "-o", final_tx_bam] + tx_bams, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        _, p1_stderr = p1.communicate()
        _, p2_stderr = p2.communicate()
        if p1.returncode != 0:
            raise RuntimeError(f"samtools cat failed (rc={p1.returncode}) for {final_bam}: {p1_stderr.strip()}")
        if p2.returncode != 0:
            raise RuntimeError(f"samtools cat failed (rc={p2.returncode}) for {final_tx_bam}: {p2_stderr.strip()}")
        
        shutil.rmtree(map_tmp_dir)
        
    else:
        # Single instance
        input_files = ",".join(filtered_bams)
        prefix = os.path.join(out_dir, f"{project}.filtered.tagged.")
        cmd = f"{star_cmd_base} --readFilesIn {input_files} --outFileNamePrefix {prefix}"
        
        star_proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        if which_stage == "Filtering":
            sam_proc = subprocess.Popen(sam_merge_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            _, sam_stderr = sam_proc.communicate()
            if sam_proc.returncode != 0:
                raise RuntimeError(f"samtools cat failed (rc={sam_proc.returncode}): {sam_stderr.strip()}")
            
        _, star_stderr = star_proc.communicate()
        if star_proc.returncode != 0:
            raise RuntimeError(f"STAR failed (rc={star_proc.returncode}): {star_stderr.strip()}")

    # Clean up chunks
    if which_stage == "Filtering":
        for f in filtered_bams:
            try: os.remove(f)
            except: pass

if __name__ == "__main__":
    main()