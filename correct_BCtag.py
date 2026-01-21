import sys
import subprocess
import os

def load_bc_map(binmap_file):
    bc_map = {}
    # Handles both comma and tab separated, and filters out 'falseBC'
    try:
        with open(binmap_file, 'r') as f:
            for line in f:
                if 'falseBC' in line:
                    continue
                parts = line.replace(',', '\t').split('\t')
                if len(parts) >= 3:
                    raw = parts[0].strip()
                    fixed = parts[2].strip()
                    bc_map[raw] = fixed
    except Exception as e:
        sys.stderr.write(f"Error loading BC map: {e}\n")
    return bc_map

def load_mgi_map(binmap2_file):
    mgi_map = {}
    try:
        with open(binmap2_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    mgi_map[parts[0]] = parts[1]
    except Exception as e:
        sys.stderr.write(f"Error loading MGI map: {e}\n")
    return mgi_map

def main():
    if len(sys.argv) < 5:
        print("Usage: python3 correct_BCtag.py <inbam> <outbam> <BCbinmap> <samtools> [BCbinmap2_for_MGI]")
        sys.exit(1)

    in_bam = sys.argv[1]
    out_bam = sys.argv[2]
    binmap1 = sys.argv[3]
    samtools = sys.argv[4]
    binmap2 = sys.argv[5] if len(sys.argv) > 5 else None

    bc_map1 = {k.encode(): v.encode() for k, v in load_bc_map(binmap1).items()}
    bc_map2 = {k.encode(): v.encode() for k, v in (load_mgi_map(binmap2) if binmap2 else {}).items()}

    # Open readers and writers in binary mode
    sam_view_cmd = [samtools, 'view', '-h', in_bam]
    sam_proc = subprocess.Popen(sam_view_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=1024*1024)
    
    bam_out_cmd = [samtools, 'view', '-b', '-o', out_bam, '-']
    bam_out_proc = subprocess.Popen(bam_out_cmd, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    bam_out = bam_out_proc.stdin

    header_finished = False

    try:
        for line in sam_proc.stdout:
            if line.startswith(b'@'):
                bam_out.write(line)
                continue
        
            if not header_finished:
                pg_line = f"@PG\tID:zUMIs-fqfilter\tPN:zUMIs-correct_BCtag\tVN:3.0\tCL:python3 {' '.join(sys.argv)}\n"
                bam_out.write(pg_line.encode())
                header_finished = True

            fields = line.strip().split(b'\t')
            if len(fields) < 12:
                continue

            tags = fields[11:]
            bc_tag_idx = -1
            ub_tag_idx = -1
            qu_tag_idx = -1
        
            for i, t in enumerate(tags):
                if t.startswith(b'BC:Z:'):
                    bc_tag_idx = i
                elif t.startswith(b'UB:Z:'):
                    ub_tag_idx = i
                elif t.startswith(b'QU:Z:'):
                    qu_tag_idx = i

            if bc_tag_idx == -1:
                raw_bc = b""
                for t in tags:
                    if t.startswith(b'BC:Z:'):
                        raw_bc = t[5:]
                        break
            else:
                raw_bc = tags[bc_tag_idx][5:]

            correct_bc = bc_map1.get(raw_bc, raw_bc)
        
            final_bc = bc_map2.get(correct_bc, correct_bc)

            if binmap2:
                flag = int(fields[1])
                if flag == 77:
                    if final_bc == correct_bc:
                        fields[9] = fields[9][3:]
                        fields[10] = fields[10][3:]
                    else:
                        ub_val = tags[ub_tag_idx][5:] if ub_tag_idx != -1 else b""
                        qu_val = tags[qu_tag_idx][5:] if qu_tag_idx != -1 else b""
                        if qu_val:
                            fields[9] = ub_val + fields[9]
                            fields[10] = qu_val + fields[10]
                            if ub_tag_idx != -1:
                                tags[ub_tag_idx] = b"UB:Z:"
                            if qu_tag_idx != -1:
                                tags[qu_tag_idx] = b"QU:Z:"
                elif flag == 141:
                    if final_bc != correct_bc:
                        if ub_tag_idx != -1:
                            tags[ub_tag_idx] = b"UB:Z:"
                        if qu_tag_idx != -1:
                            tags[qu_tag_idx] = b"QU:Z:"

            new_tags = [b"BX:Z:" + raw_bc, b"BC:Z:" + final_bc]
            for i, t in enumerate(tags):
                if i != bc_tag_idx:
                    new_tags.append(t)
        
            new_line = b"\t".join(fields[:11] + new_tags) + b"\n"
            bam_out.write(new_line)
    except BrokenPipeError:
        stderr = bam_out_proc.stderr.read().decode().strip() if bam_out_proc.stderr else ""
        raise RuntimeError(f"samtools output process terminated early while writing {out_bam}: {stderr}")

    bam_out.close()
    bam_out_proc.wait()
    sam_proc.wait()
    if sam_proc.returncode != 0:
        stderr = sam_proc.stderr.read().decode().strip() if sam_proc.stderr else ""
        raise RuntimeError(f"samtools view failed (rc={sam_proc.returncode}) for {in_bam}: {stderr}")
    if bam_out_proc.returncode != 0:
        stderr = bam_out_proc.stderr.read().decode().strip() if bam_out_proc.stderr else ""
        raise RuntimeError(f"samtools view -b failed (rc={bam_out_proc.returncode}) writing {out_bam}: {stderr}")
    except BrokenPipeError:
        stderr = bam_out_proc.stderr.read().strip() if bam_out_proc.stderr else ""
        raise RuntimeError(f"samtools output process terminated early while writing {out_bam}: {stderr}")

    bam_out.close()
    bam_out_proc.wait()
    sam_proc.wait()
    if sam_proc.returncode != 0:
        stderr = sam_proc.stderr.read().strip() if sam_proc.stderr else ""
        raise RuntimeError(f"samtools view failed (rc={sam_proc.returncode}) for {in_bam}: {stderr}")
    if bam_out_proc.returncode != 0:
        stderr = bam_out_proc.stderr.read().strip() if bam_out_proc.stderr else ""
        raise RuntimeError(f"samtools view -b failed (rc={bam_out_proc.returncode}) writing {out_bam}: {stderr}")

if __name__ == "__main__":
    main()
