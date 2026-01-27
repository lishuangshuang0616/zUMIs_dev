import os
import subprocess
import yaml
import glob
import math
import shlex
import shutil

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

def split_fastq(fq_files, n_threads, lines_per_chunk, out_dir, project, pigz_exec="pigz", seqkit_exec="seqkit", fq2_files=None, compress_chunks=False):
    """
    Splits FastQ files.
    Optimizations:
    1. If multiple files, split them in parallel (using subprocess).
    2. If seqkit is available, use it (faster than split).
    3. Fallback to GNU split.
    4. Handles PE data via fq2_files if provided (using seqkit -1 -2 or fallback).
    """
    if isinstance(fq_files, str):
        fq_files = [fq_files]
    if fq2_files and isinstance(fq2_files, str):
        fq2_files = [fq2_files]
        
    if not fq_files:
        raise ValueError("No input files provided for splitting.")
    
    # Check for seqkit
    has_seqkit = False
    if os.path.isfile(seqkit_exec) and os.access(seqkit_exec, os.X_OK):
        has_seqkit = True
    elif shutil.which(seqkit_exec):
        has_seqkit = True
    
    # Check if files are gzipped
    is_gzipped = fq_files[0].endswith('.gz')
    
    # Determine concurrency
    # We want to run multiple split jobs in parallel if we have multiple files
    # Each job consumes some threads (pigz -dc + split + pigz).
    # Approx 3 threads per job if gzipped.
    max_jobs = max(1, n_threads // 3)
    
    # Adjust threads per job based on ACTUAL number of files we process
    num_concurrent_jobs = min(len(fq_files), max_jobs)
    
    # Process files in chunks of max_jobs
    import time
    
    jobs = []
    
    mode = "SE"
    if fq2_files and len(fq2_files) > 0:
        mode = "PE"
        if len(fq_files) != len(fq2_files):
             raise ValueError(f"Mismatch in R1/R2 file counts: {len(fq_files)} vs {len(fq2_files)}")
    
    print(f"Splitting {len(fq_files)} files ({mode}). Method: {'SeqKit' if has_seqkit else 'GNU Split'}. Parallel Jobs: {num_concurrent_jobs}. Compress Output: {compress_chunks}")

    # Prepare suffixes list to return
    # We need to collect ALL suffixes generated.
    # Since we process files in parallel, we need to gather them carefully.
    # Actually, the original code returns a list of suffixes from the LAST file? 
    # Or merges them? 
    # The return value is used by fqfilter to iterate.
    # fqfilter expects chunks to exist for all files for a given suffix?
    # No, fqfilter iterates over suffixes, and for each suffix, it iterates over input files.
    # So we must ensure that for every file, the same set of suffixes exists.
    # If we split independently, we must ensure line counts match exactly so chunks match.
    # Seqkit -1 -2 guarantees this for PE.
    
    # We will collect suffixes from the first file/pair and assume consistency.
    first_file_suffixes = []

    for i, fpath1 in enumerate(fq_files):
        fpath2 = fq2_files[i] if mode == "PE" else None
        
        # Base names
        def get_base(p):
            b = os.path.basename(p)
            if b.endswith('.gz'): b = b[:-3]
            if b.endswith('.fastq'): b = b[:-6]
            if b.endswith('.fq'): b = b[:-3]
            return b

        base1 = get_base(fpath1)
        
        # File prefix logic
        # We need a unique prefix for this file (or pair).
        # We use {project}.F{i}. as the "unique" part.
        # But we must ensure the FINAL filenames match what fqfilter expects.
        # fqfilter expects: {original_base_name}{suffix}.gz
        # So we will rename the outputs to match this.
        
        # Intermediate prefix for splitting tools
        file_prefix = f"{project}.F{i:02d}."
        split_prefix_path = os.path.join(out_dir, file_prefix)
        split_prefix = split_prefix_path
        
        # Cleanup existing for this prefix
        existing = glob.glob(f"{split_prefix_path}*")
        for e in existing:
            try: os.remove(e)
            except: pass

        cmd = ""
        if has_seqkit:
            # seqkit split2 -p N (by parts)
            # The user requested splitting by thread count directly.
            # We split each file into n_threads parts.
            # This ensures maximum parallelism for downstream steps even if files differ in size.
            # Note: -p triggers 2-pass reading (counting then splitting), but it's robust.
            
            # Distribute threads among concurrent jobs
            seqkit_threads = max(1, n_threads // num_concurrent_jobs)
            
            seqkit_cmd = shlex.quote(seqkit_exec) if os.path.exists(seqkit_exec) else seqkit_exec
            
            ext_flag = "-e .fastq" if not compress_chunks else ""
            
            if mode == "PE":
                # PE Split
                # -p n_threads
                cmd = f"{seqkit_cmd} split2 -p {n_threads} -1 {shlex.quote(fpath1)} -2 {shlex.quote(fpath2)} -O {shlex.quote(out_dir)} -f -j {seqkit_threads} {ext_flag}"
            else:
                # SE Split
                cmd = f"{seqkit_cmd} split2 -p {n_threads} -O {shlex.quote(out_dir)} -f {shlex.quote(fpath1)} -j {seqkit_threads} {ext_flag}"
        
        else:
            # GNU split Fallback
            pigz_q = shlex.quote(str(pigz_exec))
            pigz_threads = max(1, n_threads // num_concurrent_jobs // 2) 
            
            output_ext = ".gz" if compress_chunks else ""
            
            if mode == "PE":
                # PE Split with GNU split is hard to sync perfectly if we run separate processes.
                # But since we used calculated line counts, it SHOULD match.
                # We will run two commands in parallel? Or sequential?
                # To avoid desync issues if one fails, sequential is safer but slower.
                # We can chain them.
                
                # R1
                dc1 = f"{pigz_q} -p {pigz_threads} -dc {shlex.quote(fpath1)}" if is_gzipped else f"cat {shlex.quote(fpath1)}"
                
                if compress_chunks:
                    cmp1 = f"{pigz_q} -p {pigz_threads} > $FILE.gz"
                else:
                    cmp1 = "cat > $FILE"
                    
                # Use a specific suffix for R1 temp
                prefix1 = f"{split_prefix}R1."
                cmd1 = f"{dc1} | split -l {lines_per_chunk} --filter='{cmp1}' - {shlex.quote(prefix1)}"
                
                # R2
                dc2 = f"{pigz_q} -p {pigz_threads} -dc {shlex.quote(fpath2)}" if is_gzipped else f"cat {shlex.quote(fpath2)}"
                
                if compress_chunks:
                    cmp2 = f"{pigz_q} -p {pigz_threads} > $FILE.gz"
                else:
                    cmp2 = "cat > $FILE"

                prefix2 = f"{split_prefix}R2."
                cmd2 = f"{dc2} | split -l {lines_per_chunk} --filter='{cmp2}' - {shlex.quote(prefix2)}"
                
                cmd = f"{cmd1} && {cmd2}"
            else:
                decompress_cmd = f"{pigz_q} -p {pigz_threads} -dc {shlex.quote(fpath1)}" if is_gzipped else f"cat {shlex.quote(fpath1)}"
                
                if compress_chunks:
                     compress_cmd = f"{pigz_q} -p {pigz_threads} > $FILE.gz"
                else:
                     compress_cmd = "cat > $FILE"
                     
                # Standard prefix
                cmd = f"{decompress_cmd} | split -l {lines_per_chunk} --filter='{compress_cmd}' - {shlex.quote(split_prefix)}"

        jobs.append({
            'cmd': cmd, 
            'file1': fpath1, 
            'file2': fpath2,
            'base1': base1,
            'base2': get_base(fpath2) if fpath2 else None,
            'prefix': file_prefix, # For GNU split renaming
            'mode': mode,
            'has_seqkit': has_seqkit,
            'compress': compress_chunks
        })

    # Execute Jobs
    active_procs = []
    
    def start_job(job_idx):
        job = jobs[job_idx]
        name = os.path.basename(job['file1'])
        if job['mode'] == "PE":
            name += f" & {os.path.basename(job['file2'])}"
        print(f"  [Split {job_idx+1}/{len(jobs)}] {name}")
        p = subprocess.Popen(job['cmd'], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return p

    next_job_idx = 0
    
    while next_job_idx < len(jobs) and len(active_procs) < num_concurrent_jobs:
        active_procs.append(start_job(next_job_idx))
        next_job_idx += 1
        
    while active_procs:
        still_active = []
        for p in active_procs:
            if p.poll() is None:
                still_active.append(p)
            else:
                if p.returncode != 0:
                    raise RuntimeError(f"Split command failed with rc={p.returncode}")
                if next_job_idx < len(jobs):
                    still_active.append(start_job(next_job_idx))
                    next_job_idx += 1
        active_procs = still_active
        if active_procs:
            time.sleep(0.5)

    # Post-process filenames to match fqfilter expectations
    # fqfilter expects: {original_base}{suffix}.gz
    # We need to scan output dir and rename files.
    
    collected_suffixes = set()

    for job in jobs:
        # Identify generated files
        # If SeqKit:
        #   Default naming: {input_basename}.part_{part_id}.{ext}
        #   We need to rename to: {input_basename}.part_{part_id}.{ext} -> {input_basename}.part_{part_id}.gz
        #   (Wait, fqfilter logic: base_name does NOT include .gz)
        #   Input: A.fq.gz -> Base: A.fq
        #   Seqkit out: A.fq.part_001.gz (if it inserts part)
        #   fqfilter expects: A.fq + suffix + .gz
        #   If suffix is .part_001, then A.fq.part_001.gz
        #   So SeqKit default naming MIGHT JUST WORK if it inserts before extension?
        #   Let's check seqkit docs pattern:
        #   "R1.fq.gz" -> "R1.part_001.fq.gz"
        #   Base (fqfilter) = "R1.fq"
        #   Base + suffix + ".gz" = "R1.fq" + ".part_001" + ".gz" = "R1.fq.part_001.gz"
        #   This is CLOSE to "R1.part_001.fq.gz" but position of extension differs.
        #   "R1.fq.part_001.gz" != "R1.part_001.fq.gz"
        
        #   So we MUST rename SeqKit output to append suffix at the end (before .gz).
        
        targets = []
        if job['mode'] == "PE":
            targets.append((job['file1'], job['base1']))
            targets.append((job['file2'], job['base2']))
        else:
            targets.append((job['file1'], job['base1']))
            
        for fpath, base in targets:
            # Find files generated from this input
            # Seqkit default: starts with input basename?
            # Or strict filename?
            input_fname = os.path.basename(fpath)
            # Pattern: {input_fname without extension?}.part_*.gz ?
            # User example: R1.fq.gz -> R1.part_001.fq.gz
            # It seems to preserve the secondary extension (.fq)?
            
            # We search for files starting with the input filename part
            # Or just glob *part_*.gz and match?
            
            # More robust: glob the out_dir
            candidates = glob.glob(os.path.join(out_dir, "*"))
            
            for c in candidates:
                c_name = os.path.basename(c)
                if job['compress']:
                     if not c_name.endswith('.gz'): continue
                else:
                     # If NOT compressed, we expect NO .gz
                     if c_name.endswith('.gz'): continue
                
                # Check if this file belongs to current input
                # Heuristic: starts with base name?
                # base = R1.fq (from R1.fq.gz)
                # Seqkit out: R1.part_001.fq.gz
                
                # Check if it matches SeqKit pattern or GNU split pattern
                
                new_suffix = None
                
                if job['has_seqkit']:
                    # Seqkit Pattern: R1.part_001.fq.gz
                    # We want: R1.fq.part_001.gz
                    # Check if c_name looks like seqkit output for this file
                    # input_fname = R1.fq.gz
                    # Check if c_name contains input_fname parts
                    
                    # Hack: split by .part_
                    if ".part_" in c_name:
                        parts = c_name.split(".part_")
                        prefix_part = parts[0] # R1
                        suffix_part = parts[1] # 001.fq.gz
                        
                        # Is prefix_part related to our base?
                        # base = R1.fq
                        # prefix_part = R1
                        
                        if base.startswith(prefix_part):
                            # It's likely ours.
                            # Extract ID: 001
                            # suffix_part might be "001.fq.gz"
                            
                            # Extract numeric ID
                            m_id = suffix_part.split('.')[0] # 001
                            
                            # Construct new name: {base}.part_{id}.gz
                            # R1.fq.part_001.gz
                            
                            ext_str = ".gz" if job['compress'] else ""
                            new_name = f"{base}.part_{m_id}{ext_str}"
                            new_path = os.path.join(out_dir, new_name)
                            
                            if c != new_path:
                                os.rename(c, new_path)
                            
                            new_suffix = f".part_{m_id}"
                            collected_suffixes.add(new_suffix)

                else:
                    # GNU Split Pattern
                    # We used prefix: {project}.F{i}.R1. OR {project}.F{i}.
                    # Output: {prefix}xaa.gz
                    
                    if job['mode'] == "PE":
                         # Prefix was {split_prefix}R1. -> ...F00.R1.
                         if f"R1." in c_name and c_name.startswith(f"{project}.F"):
                             # It's R1
                             # c_name: Project.F00.R1.xaa.gz
                             # Suffix: xaa
                             # We want: {base}.xaa.gz
                             
                             # Extract suffix from end
                             # remove .gz -> xaa
                             # remove prefix -> ...
                             
                             # Identify if this file belongs to THIS job
                             # Job prefix: job['prefix'] -> Project.F00.
                             if job['prefix'] in c_name:
                                 # ...F00.R1.xaa.gz
                                 # suffix is xaa
                                 parts = c_name.split('.')
                                 # parts: [Project, F00, R1, xaa, gz]
                                 
                                 # If compressed: [Project, F00, R1, xaa, gz] -> suffix at -2
                                 # If uncompressed: [Project, F00, R1, xaa] -> suffix at -1
                                 
                                 if job['compress']:
                                     suffix_code = parts[-2]
                                 else:
                                     suffix_code = parts[-1]
                                 
                                 ext_str = ".gz" if job['compress'] else ""
                                 new_name = f"{base}.{suffix_code}{ext_str}"
                                 new_path = os.path.join(out_dir, new_name)
                                 os.rename(c, new_path)
                                 
                                 new_suffix = f".{suffix_code}"
                                 collected_suffixes.add(new_suffix)
                                 
                         elif f"R2." in c_name and c_name.startswith(f"{project}.F"):
                             # R2
                             if job['prefix'] in c_name:
                                 parts = c_name.split('.')
                                 if job['compress']:
                                     suffix_code = parts[-2]
                                 else:
                                     suffix_code = parts[-1]
                                 
                                 ext_str = ".gz" if job['compress'] else ""
                                 new_name = f"{base}.{suffix_code}{ext_str}"
                                 new_path = os.path.join(out_dir, new_name)
                                 os.rename(c, new_path)
                                 # Don't add to collected_suffixes (R1 adds it)
                    else:
                        # SE GNU Split
                        # Prefix: Project.F00.
                        # Output: Project.F00.xaa.gz
                        if c_name.startswith(job['prefix']):
                            # Suffix: xaa
                            parts = c_name.split('.')
                            if job['compress']:
                                suffix_code = parts[-2]
                            else:
                                suffix_code = parts[-1]
                            
                            ext_str = ".gz" if job['compress'] else ""
                            new_name = f"{base}.{suffix_code}{ext_str}"
                            new_path = os.path.join(out_dir, new_name)
                            os.rename(c, new_path)
                            
                            new_suffix = f".{suffix_code}"
                            collected_suffixes.add(new_suffix)

    return sorted(list(collected_suffixes))

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
    bam_files = glob.glob(os.path.join(tmp_dir, f"{project}.*.raw.tagged.bam"))
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
