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
    Splits a FastQ file into chunks based on number of threads.
    Replaces splitfq.sh logic but without depending on GNU split --filter.
    """
    base_name = os.path.basename(fq_file)
    if base_name.endswith('.gz'):
        base_name = base_name[:-3]
    
    # Calculate lines per chunk
    # n_reads is total reads. n_chunks = n_threads.
    # Each read is 4 lines.
    reads_per_chunk = math.ceil(n_reads / n_threads)
    lines_per_chunk = reads_per_chunk * 4
    
    print(f"Splitting {fq_file} into {n_threads} chunks (~{reads_per_chunk} reads each)...")
    
    # Open input stream
    if fq_file.endswith('.gz'):
        # Use pigz for fast decompression
        p_in = subprocess.Popen([pigz_exec, '-dc', fq_file], stdout=subprocess.PIPE, bufsize=1024*1024)
        input_stream = p_in.stdout
    else:
        input_stream = open(fq_file, 'rb')

    chunk_idx = 0
    line_count = 0
    current_p_out = None
    current_out_fh = None
    buffer = bytearray()
    flush_bytes = 1024 * 1024
    prefixes = []

    try:
        # Create first chunk
        prefix_suffix = f"{base_name}{project}{chr(ord('a') + chunk_idx // 26)}{chr(ord('a') + chunk_idx % 26)}"
        out_path = os.path.join(out_dir, f"{prefix_suffix}.gz")
        prefixes.append(prefix_suffix)
        
        # Use pigz for fast compression
        current_out_fh = open(out_path, 'wb')
        current_p_out = subprocess.Popen([pigz_exec, '-c'], stdin=subprocess.PIPE, stdout=current_out_fh)
        
        for line in input_stream:
            buffer.extend(line)
            line_count += 1
            if len(buffer) >= flush_bytes:
                current_p_out.stdin.write(buffer)
                buffer.clear()
            
            if line_count >= lines_per_chunk and (line_count % 4 == 0):
                if buffer:
                    current_p_out.stdin.write(buffer)
                    buffer.clear()
                # Close current chunk
                current_p_out.stdin.close()
                current_p_out.wait()
                current_out_fh.close()
                if current_p_out.returncode != 0:
                    raise RuntimeError(f"{pigz_exec} failed (rc={current_p_out.returncode}) while writing {out_path}")
                
                # Start next chunk if not end of stream (handled by next loop iteration usually, but we check logic)
                chunk_idx += 1
                if chunk_idx >= n_threads: 
                     # If we exceeded expected chunks (due to integer math), just keep writing to the last one?
                     # Or create a new one. The original script creates up to n_threads.
                     # But split command creates as many as needed. Let's just create new one.
                     pass
                
                line_count = 0
                prefix_suffix = f"{base_name}{project}{chr(ord('a') + chunk_idx // 26)}{chr(ord('a') + chunk_idx % 26)}"
                out_path = os.path.join(out_dir, f"{prefix_suffix}.gz")
                prefixes.append(prefix_suffix)
                
                current_out_fh = open(out_path, 'wb')
                current_p_out = subprocess.Popen([pigz_exec, '-c'], stdin=subprocess.PIPE, stdout=current_out_fh)

    finally:
        if current_p_out and current_p_out.stdin:
            try:
                if buffer:
                    current_p_out.stdin.write(buffer)
                    buffer.clear()
                current_p_out.stdin.close()
            except Exception:
                pass
            current_p_out.wait()
            if current_p_out.returncode != 0:
                raise RuntimeError(f"{pigz_exec} failed (rc={current_p_out.returncode}) while writing {out_path}")
        if current_out_fh:
            try:
                current_out_fh.close()
            except Exception:
                pass
        if fq_file.endswith('.gz'):
            p_in.wait()
            if p_in.returncode != 0:
                raise RuntimeError(f"{pigz_exec} failed (rc={p_in.returncode}) while reading {fq_file}")
        else:
            input_stream.close()

    # Generate listPrefix file expected by fqfilter
    # The suffix logic above needs to match what fqfilter expects or we pass the list explicitly.
    # The original script generated suffixes like 'aa', 'ab'.
    # My logic: base + project + suffix.
    # Original: $t$pref$project + suffix.
    # Returns the list of SUFFIXES (or full prefixes) for fqfilter to iterate.
    
    # Actually, fqfilter_v2.pl iterated over the output of `ls`.
    # Our new fqfilter.py iterates over inputs passed to it?
    # Wait, the new fqfilter.py (python version) I wrote calculates path: 
    # chunk_path = os.path.join(out_dir, "zUMIs_output/.tmpMerge", f"{base_name}{tmp_prefix}.gz")
    # So I just need to return the list of 'suffix' strings (e.g., 'projectaa', 'projectab').
    
    return [f"{project}{chr(ord('a') + i // 26)}{chr(ord('a') + i % 26)}" for i in range(len(prefixes))]


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
                # Flag 4 means unmapped. But assuming fqfilter outputs standard flags?
                # fqfilter.py writes: SE -> 4 (unmapped? wait), PE -> 77/141.
                # In fqfilter.py:
                # SE: flag 4 (segment unmapped). Wait, if it's unmapped, STAR will map it later.
                # PE: 77 (paired, unmapped, etc), 141 (paired, unmapped, second in pair).
                
                # Logic from mergeBAM.sh:
                # if [[ $flag == 4 ]]; then SE else PE
                
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
