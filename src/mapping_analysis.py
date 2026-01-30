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

def run_star_pipe(corrector_args, star_cmd_str):
    """
    Runs a pipeline: Corrector (Python) -> STAR
    """
    print(f"Starting Pipeline: {' '.join(corrector_args[:3])}... -> STAR")
    
    # Start Producer (Corrector)
    # Use list args to avoid shell quoting issues with many files
    p1 = subprocess.Popen(corrector_args, stdout=subprocess.PIPE)
    
    # Start Consumer (STAR)
    # star_cmd_str should use --readFilesIn /dev/stdin
    p2 = subprocess.Popen(star_cmd_str, shell=True, stdin=p1.stdout)
    
    # Close p1's stdout in this parent process so only p2 holds it
    p1.stdout.close()
    
    # Wait for completion
    p2.wait()
    p1.wait()
    
    if p2.returncode != 0:
        raise RuntimeError(f"STAR failed with return code {p2.returncode}")
    
    # If STAR succeeds, p1 should also succeed (0) or SIGPIPE (141/-13)
    if p1.returncode not in (0, -13, 141):
         print(f"Warning: Corrector script exited with code {p1.returncode}")

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('yaml_file')
    parser.add_argument('--umi_bam', required=False, help="Merged UMI Unmapped BAM")
    parser.add_argument('--internal_bam', required=False, help="Merged Internal Unmapped BAM")
    parser.add_argument('--expect_id_file', required=False, help="Path to expect_id_barcode.tsv")
    args = parser.parse_args()
        
    yaml_file = args.yaml_file
    config = load_config(yaml_file)
    
    project = config['project']
    out_dir = config['out_dir']
    num_threads = int(config.get('num_threads', 1))
    
    samtools = config.get('samtools_exec', 'samtools')
    star_exec = config.get('STAR_exec', 'STAR')
    star_index = config['reference']['STAR_index']
    
    # Executables Check
    if not shutil.which(star_exec) and not os.path.exists(star_exec):
         print(f"Error: STAR executable not found: {star_exec}")
         sys.exit(1)

    # 1. Parse Inputs
    # Support both command line args (comma separated list) and legacy/yaml lookup
    umi_bams = []
    if args.umi_bam:
        umi_bams = [x.strip() for x in args.umi_bam.split(',') if x.strip()]
    
    internal_bams = []
    if args.internal_bam:
        internal_bams = [x.strip() for x in args.internal_bam.split(',') if x.strip()]
        
    # Explicit expect_id_file argument overrides everything
    if args.expect_id_file:
        expect_id_file = args.expect_id_file
    else:
        # Fallback (legacy logic)
        barcode_config_path = config['barcodes'].get('barcode_file')
        if barcode_config_path:
            config_dir = os.path.dirname(barcode_config_path)
            expect_id_file = os.path.join(config_dir, "expect_id_barcode.tsv")
        else:
            root_dir = os.path.dirname(out_dir.rstrip(os.sep))
            expect_id_file = os.path.join(root_dir, "config", "expect_id_barcode.tsv")

    if not os.path.exists(expect_id_file):
        raise FileNotFoundError(f"ID Map file not found: {expect_id_file}")
    
    # Check existence
    for f in umi_bams:
        if not os.path.exists(f):
            raise FileNotFoundError(f"Input UMI BAM not found: {f}")
    for f in internal_bams:
        if not os.path.exists(f):
            raise FileNotFoundError(f"Input Internal BAM not found: {f}")

    if not umi_bams and not internal_bams:
        raise FileNotFoundError("No input unmapped BAMs found.")

    # 2. Setup GTF
    final_gtf, param_add_fa = setup_gtf(config, project, out_dir, samtools)
    
    # 3. Determine Read Length (Use first available UMI bam)
    read_len = 0
    if umi_bams:
        read_len = get_bam_read_length(umi_bams[0], samtools)
    elif internal_bams:
        read_len = get_bam_read_length(internal_bams[0], samtools)
        
    print(f"Detected Read Length: {read_len}")
    
    # 4. Resource Allocation
    print(f"Allocating {num_threads} threads for sequential execution.")

    # 5. Build STAR Commands
    read_layout = config.get('read_layout', 'SE')
    
    # Define paths for corrector
    corrector_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "stream_corrector.py")
    bc_bin_file = os.path.join(out_dir, "zUMIs_output", f"{project}.BCbinning.txt")
    
    # Base params (Common)
    # Important: readFilesType SAM because our corrector outputs uncompressed BAM (which STAR treats as SAM/BAM stream)
    # STAR auto-detects BAM vs SAM if we say SAM usually, or we can use BAM Unsorted
    # Fix: Set sjdbOverhang to 100 to match typical index generation, avoiding mismatch errors with shorter reads
    misc_base = f"--genomeDir {star_index} --sjdbGTFfile {final_gtf} --runThreadN {num_threads} --sjdbOverhang {read_len - 1} --readFilesType SAM {read_layout} --outSAMmultNmax 1 --outFilterMultimapNmax 50 --outSAMunmapped Within --outSAMtype BAM Unsorted --limitOutSJcollapsed 5000000"
    
    extra_params = config['reference'].get('additional_STAR_params', '') or ""
    
    # Two-pass mode
    twopass = ""
    if config['counting_opts'].get('twoPass', False):
        twopass = "--twopassMode Basic"

    # Run UMI
    if umi_bams:
        prefix_umi = os.path.join(out_dir, f"{project}.filtered.tagged.umi.")
        
        # Corrector Args (Producer)
        corrector_args = ['python3', corrector_script, '--binning', bc_bin_file, '--idmap', expect_id_file, '--type', 'umi'] + umi_bams
        
        # STAR Command (Consumer)
        # --readFilesIn /dev/stdin
        cmd_umi = f"{star_exec} {misc_base} {extra_params} {param_add_fa} {twopass} --readFilesIn /dev/stdin --outFileNamePrefix {prefix_umi}"
        
        run_star_pipe(corrector_args, cmd_umi)
        print("STAR UMI finished.")

    # Run Internal
    if internal_bams:
        prefix_int = os.path.join(out_dir, f"{project}.filtered.tagged.internal.")
        
        # Corrector Args
        corrector_args = ['python3', corrector_script, '--binning', bc_bin_file, '--idmap', expect_id_file, '--type', 'internal'] + internal_bams
        
        # STAR Command
        cmd_int = f"{star_exec} {misc_base} {extra_params} {param_add_fa} {twopass} --readFilesIn /dev/stdin --outFileNamePrefix {prefix_int}"
        
        run_star_pipe(corrector_args, cmd_int)
        print("STAR Internal finished.")
        


if __name__ == "__main__":
    main()
