import os
import subprocess
import yaml
import glob
import math
import shlex

def run_shell_cmd(cmd, step_name, log_file=None):
    is_shell = isinstance(cmd, str)
    cmd_str = cmd if is_shell else " ".join(shlex.quote(str(x)) for x in cmd)

    print(f"[{step_name}] Running: {cmd_str}")
    if log_file:
        with open(log_file, 'a') as f:
            f.write(f"Running: {cmd_str}\n")
            process = subprocess.run(cmd, shell=is_shell, stdout=f, stderr=subprocess.STDOUT)
    else:
        process = subprocess.run(cmd, shell=is_shell)
    
    if process.returncode != 0:
        raise Exception(f"Error in step [{step_name}]. Command failed: {cmd_str}")

def split_fastq(fq_file, n_threads, n_reads, out_dir, project, pigz_exec="pigz"):
    """
    Splits a FastQ file into exactly n_threads chunks using Batch Round Robin.
    This ensures perfectly balanced chunks without needing accurate read count estimation.
    """
    base_name = os.path.basename(fq_file)
    if base_name.endswith('.gz'):
        base_name = base_name[:-3]
    
    batch_size = 200000 # lines per batch (50,000 reads)
    print(f"Splitting {fq_file} into {n_threads} balanced chunks (batch size: {batch_size} lines)...")
    
    # Open input stream
    if fq_file.endswith('.gz'):
        p_in = subprocess.Popen([pigz_exec, '-dc', fq_file], stdout=subprocess.PIPE, bufsize=1024*1024)
        input_stream = p_in.stdout
    else:
        input_stream = open(fq_file, 'rb')

    # Initialize all output processes and file handles
    out_procs = []
    out_fhs = []
    prefixes = []
    buffers = [bytearray() for _ in range(n_threads)]
    flush_threshold = 1024 * 1024 # 1MB internal buffer before pipe write
    
    try:
        for i in range(n_threads):
            suffix = f"{chr(ord('a') + i // 26)}{chr(ord('a') + i % 26)}"
            prefix_suffix = f"{project}{suffix}"
            full_prefix = f"{base_name}{prefix_suffix}"
            out_path = os.path.join(out_dir, f"{full_prefix}.gz")
            
            prefixes.append(prefix_suffix)
            fh = open(out_path, 'wb')
            # Use pigz for fast parallel compression
            p = subprocess.Popen([pigz_exec, '-c'], stdin=subprocess.PIPE, stdout=fh, bufsize=1024*1024)
            
            out_fhs.append(fh)
            out_procs.append(p)

        current_chunk = 0
        line_count = 0
        
        for line in input_stream:
            buffers[current_chunk].extend(line)
            line_count += 1
            
            # Flush buffer to pipe if threshold reached
            if len(buffers[current_chunk]) >= flush_threshold:
                out_procs[current_chunk].stdin.write(buffers[current_chunk])
                buffers[current_chunk].clear()
            
            # Switch chunk after reaching batch size (must be multiple of 4 lines)
            if line_count >= batch_size:
                # Flush remaining buffer for current chunk before switching
                if buffers[current_chunk]:
                    out_procs[current_chunk].stdin.write(buffers[current_chunk])
                    buffers[current_chunk].clear()
                
                current_chunk = (current_chunk + 1) % n_threads
                line_count = 0

    finally:
        # Final flush and close all processes
        for i in range(n_threads):
            if buffers[i]:
                try: out_procs[i].stdin.write(buffers[i])
                except: pass
            
            if out_procs[i].stdin:
                out_procs[i].stdin.close()
            
        for p in out_procs:
            p.wait()
            
        for fh in out_fhs:
            fh.close()
            
        if fq_file.endswith('.gz'):
            p_in.wait()
        else:
            input_stream.close()

    return prefixes


def merge_bam_stats(tmp_dir, project, out_dir, yaml_file, samtools_exec):
    """
    Replaces mergeBAM.sh
    1. Concatenates BCstats.txt
    2. Checks read layout (SE/PE) from first BAM and updates YAML.
    """
    
    # 1. Cat stats
    stats_files = glob.glob(os.path.join(tmp_dir, f"{project}.*.BCstats.txt"))
    out_stats = os.path.join(out_dir, f"{project}.BCstats.txt")
    
    bc_counts = {}
    for fname in stats_files:
        with open(fname, 'r') as infile:
            for line in infile:
                parts = line.strip().split()
                if len(parts) >= 2:
                    bc = parts[0]
                    try:
                        count = int(parts[1])
                        bc_counts[bc] = bc_counts.get(bc, 0) + count
                    except ValueError:
                        continue

    with open(out_stats, 'w') as outfile:
        for bc in sorted(bc_counts.keys()):
            outfile.write(f"{bc}\t{bc_counts[bc]}\n")
                
    # 2. Check Layout
    bam_files = glob.glob(os.path.join(tmp_dir, f"{project}.*.filtered.tagged.bam"))
    if not bam_files:
        print("No BAM files found to check layout.")
        return

    first_bam = bam_files[0]
    
    # Check flag of first read
    try:
        proc = subprocess.Popen(
            [samtools_exec, 'view', first_bam],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )
        line = proc.stdout.readline() if proc.stdout else ""
        proc.terminate()
        _, stderr = proc.communicate()
        
        if line:
            parts = line.split('\t')
            if len(parts) > 1:
                flag = int(parts[1])
                layout = "SE" if flag == 4 else "PE"
                
                # Update YAML
                with open(yaml_file, 'r') as f:
                    ydata = yaml.safe_load(f)
                
                ydata['read_layout'] = layout

                class ZumisDumper(yaml.SafeDumper):
                    pass

                def bool_representer(dumper, value):
                    return dumper.represent_scalar('tag:yaml.org,2002:bool', 'yes' if value else 'no')

                def none_representer(dumper, _value):
                    return dumper.represent_scalar('tag:yaml.org,2002:null', '~')

                yaml.add_representer(bool, bool_representer, Dumper=ZumisDumper)
                yaml.add_representer(type(None), none_representer, Dumper=ZumisDumper)

                with open(yaml_file, 'w') as f:
                    yaml.dump(ydata, f, Dumper=ZumisDumper, default_flow_style=False, sort_keys=False)
                print(f"Detected Read Layout: {layout}")
        elif proc.returncode not in (0, -15):
            raise RuntimeError(f"samtools view produced no output (rc={proc.returncode}): {stderr.strip()}")

    except Exception as e:
        print(f"Error checking BAM layout: {e}")