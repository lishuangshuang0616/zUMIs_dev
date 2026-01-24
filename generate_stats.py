#!/usr/bin/env python3
import sys
import os
import yaml
import subprocess
import csv
import collections
import statistics
import math
import gzip

import numpy as np
import bisect

# Try to import plotting libraries
try:
    import matplotlib
    matplotlib.use('Agg') # Non-interactive backend
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("Warning: matplotlib not found. Plots will not be generated.")

def load_config(yaml_file):
    with open(yaml_file, 'r') as f:
        return yaml.safe_load(f)

def load_gene_coords(gtf_file):
    """
    Loads gene coordinates (start, end, strand) from GTF for coverage calc.
    Returns: dict {gene_id: (chrom, start, end, strand)}
    """
    print(f"Loading gene coordinates from {gtf_file}...")
    gene_coords = {}
    
    if not os.path.exists(gtf_file):
        return {}

    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            if parts[2] != 'gene': continue 
            
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attributes = parts[8]
            
            gene_id = None
            if 'gene_id "' in attributes:
                gene_id = attributes.split('gene_id "')[1].split('"')[0]
            
            if gene_id:
                gene_coords[gene_id] = (start, end, strand)
    return gene_coords

def cigar_ref_len(cigar):
    """
    Calculates the reference genome length corresponding to the CIGAR string, 
    handling S and H operators.
    """
    if not cigar or cigar == '*':
        return 0
    num = 0
    ref_len = 0
    for ch in cigar:
        o = ord(ch)
        if 48 <= o <= 57:
            num = num * 10 + (o - 48)  # Digit
            continue
        if ch in ('M', 'D', 'N', '=', 'X'):
            ref_len += num  # Add to reference length
        num = 0
    return ref_len

def load_exon_models(gtf_file):
    """
    Loads exon models, merges adjacent exons, and generates gene models.
    Returns: dict {gene_id: {'chrom': str, 'strand': str, 'starts': list, 'ends': list, 'cum': list, 'total_len': int}}
    """
    print(f"Loading exon models from {gtf_file}...")
    gene_exons = collections.defaultdict(list)
    gene_strand = {}
    gene_chrom = {}
    
    if not os.path.exists(gtf_file):
        return {}

    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            if parts[2] != 'exon':
                continue

            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attributes = parts[8]

            gene_id = None
            if 'gene_id "' in attributes:
                gene_id = attributes.split('gene_id "')[1].split('"')[0]
            
            if gene_id:
                gene_exons[gene_id].append((start, end))
                gene_strand[gene_id] = strand
                gene_chrom[gene_id] = chrom
    
    # Merge adjacent exons
    models = {}
    for gene_id, exons in gene_exons.items():
        exons.sort()  # Sort exons
        merged = []
        curr_s, curr_e = exons[0]
        for s, e in exons[1:]:
            if s <= curr_e + 1:  # Adjacent exons
                curr_e = max(curr_e, e)
            else:
                merged.append((curr_s, curr_e))
                curr_s, curr_e = s, e
        merged.append((curr_s, curr_e))  # Last exon

        starts = [s for s, _ in merged]
        ends = [e for _, e in merged]
        cum = []
        total = 0
        for s, e in merged:
            cum.append(total)
            total += (e - s + 1)  # Cumulative length

        models[gene_id] = {
            "chrom": gene_chrom.get(gene_id),
            "strand": gene_strand.get(gene_id, '+'),
            "starts": starts,
            "ends": ends,
            "cum": cum,
            "total_len": total,
        }

    return models

def load_sparse_stats(matrix_dir):
    """
    Reads matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz from a directory
    and calculates stats directly without loading full dense matrix.
    """
    matrix_path = os.path.join(matrix_dir, "matrix.mtx.gz")
    barcodes_path = os.path.join(matrix_dir, "barcodes.tsv.gz")
    features_path = os.path.join(matrix_dir, "features.tsv.gz")
    
    if not (os.path.exists(matrix_path) and os.path.exists(barcodes_path) and os.path.exists(features_path)):
        # Fallback to unzipped
        matrix_path = matrix_path[:-3]
        barcodes_path = barcodes_path[:-3]
        features_path = features_path[:-3]
        if not (os.path.exists(matrix_path) and os.path.exists(barcodes_path) and os.path.exists(features_path)):
            return {}, {}

    # Load Barcodes
    barcodes = []
    opener = gzip.open if barcodes_path.endswith('.gz') else open
    with opener(barcodes_path, 'rt') as f:
        barcodes = [line.strip() for line in f]
        
    # Load Features (Genes)
    genes = []
    opener = gzip.open if features_path.endswith('.gz') else open
    with opener(features_path, 'rt') as f:
        genes = [line.strip().split('\t')[0] for line in f]
        
    # Stats containers
    genes_per_cell = collections.defaultdict(int)
    counts_per_gene = collections.defaultdict(int)
    
    # Parse MTX
    opener = gzip.open if matrix_path.endswith('.gz') else open
    with opener(matrix_path, 'rt') as f:
        header_passed = False
        for line in f:
            if line.startswith('%'): continue
            if not header_passed:
                header_passed = True
                continue
            
            parts = line.split()
            if len(parts) < 3: continue
            
            # MTX is 1-based
            row_idx = int(parts[0]) - 1
            col_idx = int(parts[1]) - 1
            val = float(parts[2]) # Can be int or float
            
            if val > 0:
                bc = barcodes[col_idx]
                gene = genes[row_idx]
                
                genes_per_cell[bc] += 1
                counts_per_gene[gene] += int(val)
                
    return genes_per_cell, counts_per_gene

def parse_bam_stats(bam_file, samtools_exec, kept_barcodes, gene_models):
    """
    Calculates Mapping Stats & Coverage from BAM file.
    Returns: dict {Barcode: {Category: Count}} and coverage arrays.
    """
    print(f"Calculating Mapping Stats & Coverage from {bam_file}...")
    stats = collections.defaultdict(lambda: collections.defaultdict(int))
    cov_umi = np.zeros(100, dtype=np.int64)
    cov_int = np.zeros(100, dtype=np.int64)
    
    cmd = [samtools_exec, 'view', bam_file]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1)
    
    count = 0
    try:
        for line in proc.stdout:
            count += 1
            if count % 1000000 == 0:
                print(f"Processed {count} reads...", end='\r')

            parts = line.strip().split('\t')
            if len(parts) < 12: continue

            flag = int(parts[1])
            if flag & 0x100 or flag & 0x800: continue  # Exclude secondary/supplementary
            
            # Count only Read 1 to avoid double counting fragments in PE data
            if flag & 0x1: # Paired
                if not (flag & 0x40): # If not Read 1, skip
                    continue
            
            tags = {}
            for t in parts[11:]:
                if ':' not in t:
                    continue
                tag_parts = t.split(':', 2)
                if len(tag_parts) != 3:
                    continue
                tags[tag_parts[0]] = tag_parts[2]
            
            bc = tags.get('CB')
            if not bc:
                stats["__NO_CB__"]["Unused BC"] += 1
                continue
            bc = bc.strip()
            
            # --- 1. Mapping Stats ---
            category = "Intergenic"
            source = tags.get('SR')
            if bc in kept_barcodes:
                if source == 'UMI':
                    stats[bc]['UMI_Reads'] += 1
                elif source == 'Internal':
                    stats[bc]['Internal_Reads'] += 1
            
            if 'RE' in tags:
                re_val = tags['RE']
                if re_val == 'E':
                    category = 'Exon'
                elif re_val == 'N':
                    category = 'Intron'
                elif re_val == 'I':
                    category = 'Intergenic'
                else:
                    category = re_val
            else:
                if 'XS' in tags:
                    status = tags['XS']
                    if status.startswith("Unassigned_"):
                        status = status.replace("Unassigned_", "")
                    if status == "NoFeatures": category = "Intergenic"
                    else: category = status
                else:
                    flag = int(parts[1])
                    if flag & 0x4: category = "Unmapped"
            
            stats[bc][category] += 1
            
            # --- 2. Coverage Calculation ---
            if bc in kept_barcodes and 'GX' in tags and source in ['UMI', 'Internal']:
                gene_id = tags['GX']
                model = gene_models.get(gene_id)
                if not model:
                    continue
                if model["total_len"] < 100:
                    continue

                pos_start = int(parts[3])
                ref_len = cigar_ref_len(parts[5])
                if ref_len <= 0:
                    continue
                pos = pos_start + (ref_len - 1) // 2

                idx = bisect.bisect_right(model["starts"], pos) - 1
                if idx < 0:
                    continue
                if pos > model["ends"][idx]:
                    continue

                forward_pos = model["cum"][idx] + (pos - model["starts"][idx])
                if forward_pos < 0 or forward_pos >= model["total_len"]:
                    continue

                if model["strand"] == '-':
                    transcript_pos = model["total_len"] - 1 - forward_pos
                else:
                    transcript_pos = forward_pos

                denom = model["total_len"] - 1
                if denom <= 0:
                    continue
                rel_pos = transcript_pos / denom
                if rel_pos < 0.0:
                    rel_pos = 0.0
                elif rel_pos > 1.0:
                    rel_pos = 1.0
                bin_idx = int(rel_pos * 99)
                if source == 'UMI':
                    cov_umi[bin_idx] += 1
                else:
                    cov_int[bin_idx] += 1

    finally:
        proc.stdout.close()
        proc.wait()
        print(f"\nFinished parsing {count} reads.")
        
    return stats, cov_umi, cov_int

def plot_coverage(cov_umi, cov_int, out_prefix):
    if not HAS_MATPLOTLIB: return
    
    def smooth(y, window=7):
        y = np.asarray(y, dtype=np.float64)
        if window <= 1:
            return y
        if window % 2 == 0:
            window += 1
        if y.size < window:
            return y
        kernel = np.ones(window, dtype=np.float64) / window
        pad = window // 2
        y_pad = np.pad(y, (pad, pad), mode='edge')
        return np.convolve(y_pad, kernel, mode='valid')

    sum_umi = float(np.sum(cov_umi))
    sum_int = float(np.sum(cov_int))
    frac_umi = (cov_umi / sum_umi) if sum_umi > 0 else np.zeros_like(cov_umi, dtype=np.float64)
    frac_int = (cov_int / sum_int) if sum_int > 0 else np.zeros_like(cov_int, dtype=np.float64)

    max_umi = float(np.max(cov_umi)) if cov_umi.size else 0.0
    max_int = float(np.max(cov_int)) if cov_int.size else 0.0
    maxnorm_umi = (cov_umi / max_umi) if max_umi > 0 else np.zeros_like(cov_umi, dtype=np.float64)
    maxnorm_int = (cov_int / max_int) if max_int > 0 else np.zeros_like(cov_int, dtype=np.float64)

    maxnorm_umi_s = smooth(maxnorm_umi, window=7)
    maxnorm_int_s = smooth(maxnorm_int, window=7)

    x = np.arange(1, 101)

    plt.figure(figsize=(10, 8))
    plt.plot(x, maxnorm_int_s, label='Inter Coverage', color='blue', linewidth=2)
    plt.plot(x, maxnorm_umi_s, label='UMI Coverage', color='red', linewidth=2)

    plt.ylim(0, 1.05)
    plt.xlabel("Gene Body Percentile (5'->3')")
    plt.ylabel("Coverage")
    plt.title("Gene Body Coverage Comparison")
    plt.grid(True, linestyle='-', alpha=0.2)
    plt.legend(frameon=False, loc="upper right")

    plt.tight_layout()
    plt.savefig(f"{out_prefix}.geneBodyCoverage.pdf")
    plt.close()

    with open(f"{out_prefix}.geneBodyCoverage.txt", 'w') as f:
        f.write("Percentile\tUMI_Reads\tInternal_Reads\tUMI_Frac\tInternal_Frac\tUMI_MaxNorm\tInternal_MaxNorm\n")
        for i in range(100):
            f.write(f"{i+1}\t{int(cov_umi[i])}\t{int(cov_int[i])}\t{frac_umi[i]:.6f}\t{frac_int[i]:.6f}\t{maxnorm_umi[i]:.6f}\t{maxnorm_int[i]:.6f}\n")

def plot_gene_counts(genes_per_cell_umi, genes_per_cell_read, out_pdf):
    if not HAS_MATPLOTLIB: return
    
    fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    
    # UMI Data
    data_umi = list(genes_per_cell_umi.values())
    if data_umi:
        ax[0].boxplot(data_umi)
        ax[0].set_title("Genes per Cell (UMI)")
        ax[0].set_ylabel("Number of Genes")
        med = statistics.median(data_umi)
        ax[0].text(1.1, med, f"Median: {med}", verticalalignment='center')
    
    # Read Data
    data_read = list(genes_per_cell_read.values())
    if data_read:
        ax[1].boxplot(data_read)
        ax[1].set_title("Genes per Cell (Read)")
        ax[1].set_ylabel("Number of Genes")
        med = statistics.median(data_read)
        ax[1].text(1.1, med, f"Median: {med}", verticalalignment='center')
        
    plt.tight_layout()
    plt.savefig(out_pdf)
    plt.close()

def plot_gene_umi_counts_by_type(stats_exon, stats_intron, stats_inex, wells, out_pdf):
    if not HAS_MATPLOTLIB:
        return
    types = ["Exon", "Intron+Exon", "Intron"]
    # zUMIs colors
    colors = ["#1A5084", "#914614", "#118730"]
    
    gene_data = [
        [stats_exon.get(w, {}).get("genes", 0) for w in wells],
        [stats_inex.get(w, {}).get("genes", 0) for w in wells],
        [stats_intron.get(w, {}).get("genes", 0) for w in wells],
    ]
    umi_data = [
        [stats_exon.get(w, {}).get("umis", 0) for w in wells],
        [stats_inex.get(w, {}).get("umis", 0) for w in wells],
        [stats_intron.get(w, {}).get("umis", 0) for w in wells],
    ]

    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    
    # Gene Counts Plot
    bp0 = axes[0].boxplot(gene_data, labels=types, notch=True, patch_artist=True)
    for patch, color in zip(bp0['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.8)
        
    axes[0].set_title("Number of Genes")
    axes[0].set_ylabel("Count")
    
    # UMI Counts Plot
    bp1 = axes[1].boxplot(umi_data, labels=types, notch=True, patch_artist=True)
    for patch, color in zip(bp1['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.8)
        
    axes[1].set_title("Number of UMIs")
    axes[1].set_ylabel("Count")

    for ax, data in zip(axes, [gene_data, umi_data]):
        for i, arr in enumerate(data, start=1):
            if not arr:
                continue
            med = statistics.median(arr)
            ax.text(i, med, f"{int(med)}", ha="center", va="bottom", color="orange", fontweight='bold')

    plt.tight_layout()
    plt.savefig(out_pdf)
    plt.close()

def plot_features(read_stats, wells, out_pdf):
    if not HAS_MATPLOTLIB:
        return
    feat_colors = {
        "Exon": "#1A5084",
        "Intron": "#118730",
        "Unmapped": "#545454",
        "Ambiguity": "#FFA54F",
        "MultiMapping": "#631879FF",
        "Intergenic": "#FFD700",
        "Unused BC": "#BABABA",
    }

    categories_order = ["Exon", "Intron", "Intergenic", "Ambiguity", "MultiMapping", "Unmapped"]
    totals = {c: 0 for c in categories_order}
    per_cell_frac = {c: [] for c in categories_order}
    unused_total = int(read_stats.get("__NO_CB__", {}).get("Unused BC", 0))

    for bc in wells:
        counts = read_stats.get(bc, {})
        cell_total = 0
        for c in categories_order:
            v = int(counts.get(c, 0))
            totals[c] += v
            cell_total += v
        if cell_total > 0:
            for c in categories_order:
                per_cell_frac[c].append(int(counts.get(c, 0)) / cell_total)
        else:
            for c in categories_order:
                per_cell_frac[c].append(0.0)

    total_sum = sum(totals.values()) + unused_total
    total_fracs = {c: (totals[c] / total_sum if total_sum > 0 else 0.0) for c in categories_order}
    unused_frac = unused_total / total_sum if total_sum > 0 else 0.0

    bar_categories = [c for c in categories_order if totals.get(c, 0) > 0]
    box_categories = [c for c in categories_order if any(v > 0 for v in per_cell_frac[c])]
    if not box_categories:
        box_categories = ["Exon", "Intron", "Intergenic", "Ambiguity", "Unmapped"]

    fig, axes = plt.subplots(2, 1, figsize=(12, 7), gridspec_kw={"height_ratios": [0.6, 1.0]})

    left = 0.0
    for c in bar_categories:
        w = total_fracs.get(c, 0.0)
        if w > 0:
            axes[0].barh([0], [w * 100.0], left=[left * 100.0], color=feat_colors[c], edgecolor="white", linewidth=0.8, label=c, height=0.5)
            left += w
    if unused_frac > 0:
        axes[0].barh([0], [unused_frac * 100.0], left=[left * 100.0], color=feat_colors["Unused BC"], edgecolor="white", linewidth=0.8, label="Unused BC", height=0.5)

    axes[0].set_xlim(0, 100)
    axes[0].set_yticks([])
    axes[0].set_xlabel("% of total reads")
    axes[0].set_title("Total Read Distribution")
    axes[0].grid(axis="x", linestyle="-", alpha=0.2)
    legend_n = len(bar_categories) + (1 if unused_frac > 0 else 0)
    axes[0].legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.30),
        ncol=max(1, legend_n),
        frameon=False,
        fontsize=9,
        handlelength=1.1,
        handletextpad=0.4,
        columnspacing=0.9,
        borderaxespad=0.0,
    )

    box_data = [[v * 100.0 for v in per_cell_frac[c]] for c in box_categories]
    bp = axes[1].boxplot(box_data, labels=box_categories, notch=True, patch_artist=True, widths=0.6, showfliers=True, flierprops={"marker": "o", "markersize": 2, "markerfacecolor": "black", "markeredgecolor": "black", "alpha": 0.5})
    for patch, c in zip(bp["boxes"], box_categories):
        patch.set_facecolor(feat_colors[c])
        patch.set_alpha(0.8)
        patch.set_edgecolor("black")
        patch.set_linewidth(0.8)

    if len(wells) <= 30:
        for i, c in enumerate(box_categories, start=1):
            vals = box_data[i - 1]
            if not vals:
                continue
            x = np.random.normal(loc=i, scale=0.04, size=len(vals))
            axes[1].scatter(x, vals, s=8, color="black", alpha=0.25, linewidths=0)

    axes[1].set_ylim(0, 100)
    axes[1].set_title("Reads per Cell")
    axes[1].set_ylabel("% reads/cell")
    axes[1].tick_params(axis="x", rotation=0)
    axes[1].grid(axis="y", linestyle="-", alpha=0.25)

    fig.tight_layout()
    plt.savefig(out_pdf)
    plt.close()

def load_barcode_mapping(out_dir, project):
    """
    Loads mapping from WellID to UMI/Internal sequences.
    Returns: dict {wellID: {'umi': seq, 'int': seq}}
    """
    mapping = {}
    
    # Priority: expect_id_barcode.tsv
    config_dir = os.path.join(os.path.dirname(out_dir.rstrip('/')), "config")
    expect_file = os.path.join(config_dir, "expect_id_barcode.tsv")
    
    if os.path.exists(expect_file):
        with open(expect_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    if parts[0] in ['wellID', 'WellID']: continue
                    well = parts[0].strip() # Critical strip
                    umi_bc = parts[1].strip()
                    int_bc = parts[2].strip()
                    mapping[well] = {'umi': umi_bc, 'int': int_bc}
    else:
        # Fallback: kept_barcodes.txt
        kept_file = os.path.join(out_dir, "zUMIs_output", f"{project}kept_barcodes.txt")
        if os.path.exists(kept_file):
            with open(kept_file, 'r') as f:
                f.readline() # header
                for line in f:
                    parts = line.strip().replace(',','\t').split('\t')
                    if parts[0]:
                        well = parts[0].strip()
                        mapping[well] = {'umi': well, 'int': "NA"}
    return mapping

def calculate_matrix_stats(matrix_dir):
    """
    Calculates UMI count and Gene count per cell from a sparse matrix directory.
    Returns: dict {bc: {'umis': int, 'genes': int, 'gene_set': set}}
    """
    stats = collections.defaultdict(lambda: {'umis': 0, 'genes': 0, 'gene_set': set()})
    
    matrix_path = os.path.join(matrix_dir, "matrix.mtx.gz")
    barcodes_path = os.path.join(matrix_dir, "barcodes.tsv.gz")
    features_path = os.path.join(matrix_dir, "features.tsv.gz")
    
    if not (os.path.exists(matrix_path) and os.path.exists(barcodes_path)):
         matrix_path = matrix_path[:-3]
         barcodes_path = barcodes_path[:-3]
         if not (os.path.exists(matrix_path) and os.path.exists(barcodes_path)):
             return stats

    # Load Barcodes
    opener = gzip.open if barcodes_path.endswith('.gz') else open
    with opener(barcodes_path, 'rt') as f:
        barcodes = [line.strip() for line in f]
        
    # Parse MTX
    opener = gzip.open if matrix_path.endswith('.gz') else open
    with opener(matrix_path, 'rt') as f:
        header_passed = False
        for line in f:
            if line.startswith('%'): continue
            if not header_passed:
                header_passed = True
                continue
            
            parts = line.split()
            if len(parts) < 3: continue
            
            # col_idx is 1-based index for barcode
            col_idx = int(parts[1]) - 1
            row_idx = int(parts[0]) - 1 # gene index
            val = float(parts[2])
            
            if val > 0:
                bc = barcodes[col_idx]
                stats[bc]['umis'] += int(val)
                stats[bc]['genes'] += 1
                stats[bc]['gene_set'].add(row_idx)
                
    return stats

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 generate_stats.py <yaml_config>")
        sys.exit(1)
        
    yaml_file = sys.argv[1]
    config = load_config(yaml_file)
    
    project = config['project']
    out_dir = config['out_dir']
    samtools = config.get('samtools_exec', 'samtools')
    
    stats_dir = os.path.join(out_dir, "zUMIs_output", "stats")
    if not os.path.exists(stats_dir): os.makedirs(stats_dir)
    
    # --- 1. Load Barcode Mapping ---
    bc_mapping = load_barcode_mapping(out_dir, project)
    all_wells = sorted(list(bc_mapping.keys()))
    
    # Use expected wells as the valid set for stats
    kept_barcodes = set(all_wells)
    print(f"Validating against {len(kept_barcodes)} expected barcodes.")
    
    # --- 2. Calculate Matrix Stats (UMI/Gene counts) ---
    print("Calculating Matrix Statistics...")
    
    umi_exon_dir = os.path.join(out_dir, "zUMIs_output", "expression", f"{project}.exon.umi")
    umi_intron_dir = os.path.join(out_dir, "zUMIs_output", "expression", f"{project}.intron.umi")
    
    stats_exon = calculate_matrix_stats(umi_exon_dir)
    stats_intron = calculate_matrix_stats(umi_intron_dir)
    
    # Calculate Intron+Exon Stats (Union)
    stats_inex = collections.defaultdict(lambda: {'umis': 0, 'genes': 0})
    for bc in set(stats_exon.keys()) | set(stats_intron.keys()):
        ex_st = stats_exon.get(bc, {'umis': 0, 'genes': 0, 'gene_set': set()})
        in_st = stats_intron.get(bc, {'umis': 0, 'genes': 0, 'gene_set': set()})
        
        stats_inex[bc]['umis'] = ex_st['umis'] + in_st['umis']
        stats_inex[bc]['genes'] = len(ex_st['gene_set'] | in_st['gene_set'])

    # Write gene counts (legacy support)
    if stats_exon:
         # Load features for legacy format
         features_path = os.path.join(umi_exon_dir, "features.tsv.gz")
         if not os.path.exists(features_path): features_path = features_path[:-3]
         genes = []
         opener = gzip.open if features_path.endswith('.gz') else open
         if os.path.exists(features_path):
             with opener(features_path, 'rt') as f:
                 genes = [line.strip().split('\t')[0] for line in f]
         
         # Count per gene
         counts_per_gene = collections.defaultdict(int)
         # Re-parsing just for gene counts is inefficient but robust
         # Or just skip legacy format? Let's skip for now or rely on matrix stats if needed.
         # The requested table is per-well, so we focus on that.
    
    # --- 3. Calculate Mapping Stats (Reads) & Coverage ---
    gtf_file = os.path.join(out_dir, f"{project}.final_annot.gtf")
    gene_models = load_exon_models(gtf_file)
    
    bam_file = os.path.join(out_dir, f"{project}.filtered.Aligned.GeneTagged.bam")
    read_stats = collections.defaultdict(lambda: collections.defaultdict(int))
    cov_umi = np.zeros(100, dtype=np.int64)
    cov_int = np.zeros(100, dtype=np.int64)
    
    if os.path.exists(bam_file):
        read_stats, cov_umi, cov_int = parse_bam_stats(bam_file, samtools, kept_barcodes, gene_models)
        
        # Plot Coverage
        out_prefix = os.path.join(stats_dir, project)
        plot_coverage(cov_umi, cov_int, out_prefix)
        
        # Write Reads per Cell (Legacy)
        with open(os.path.join(stats_dir, f"{project}.readspercell.txt"), 'w') as f:
            f.write("Barcode\tExon\tIntron\tIntergenic\tUnmapped\tAmbiguity\n")
            for bc, counts in read_stats.items():
                if bc in kept_barcodes:
                    f.write(f"{bc}\t{counts['Exon']}\t{counts['Intron']}\t{counts['Intergenic']}\t{counts['Unmapped']}\t{counts['Ambiguity']}\n")

    # --- 4. Compile Comprehensive Stats Table ---
    print(f"Compiling stats for {len(all_wells)} wells...")
    
    output_table = os.path.join(stats_dir, f"{project}.stats.tsv")
    
    with open(output_table, 'w') as f:
        # Header
        headers = [
            "wellID", "internal_barcodes", "umi_barcodes", 
            "internal_reads", "umi_reads", 
            "Ambiguity_reads", "Exon_reads", "Intergenic_reads", "intron_reads", 
            "Unmapped_reads", "Multimapped_reads",
            "Exon_umis", "Intron_umis", "Intron_Exon_umis",
            "Exon_genes", "Intron_genes", "Intron_Exon_genes"
        ]
        f.write("\t".join(headers) + "\n")
        
        for well in all_wells:
            # Metadata
            int_bc = bc_mapping[well]['int']
            umi_bc = bc_mapping[well]['umi']
            
            # Read Stats
            r_stats = read_stats.get(well, {})
            ambig = r_stats.get('Ambiguity', 0)
            exon_r = r_stats.get('Exon', 0)
            inter_r = r_stats.get('Intergenic', 0)
            intron_r = r_stats.get('Intron', 0)
            unmapped = r_stats.get('Unmapped', 0)
            multi = r_stats.get('MultiMapping', 0)
            
            internal_r = r_stats.get('Internal_Reads', 0)
            umi_r = r_stats.get('UMI_Reads', 0)
            
            # UMI Stats
            ex_st = stats_exon.get(well, {'umis': 0, 'genes': 0})
            in_st = stats_intron.get(well, {'umis': 0, 'genes': 0})
            inex_st = stats_inex.get(well, {'umis': 0, 'genes': 0})
            
            ex_umis = ex_st['umis']
            in_umis = in_st['umis']
            inex_umis = inex_st['umis']
            
            ex_genes = ex_st['genes']
            in_genes = in_st['genes']
            inex_genes = inex_st['genes']
            
            row = [
                well, int_bc, umi_bc,
                internal_r, umi_r,
                ambig, exon_r, inter_r, intron_r,
                unmapped, multi,
                ex_umis, in_umis, inex_umis,
                ex_genes, in_genes, inex_genes
            ]
            
            f.write("\t".join(map(str, row)) + "\n")

    if HAS_MATPLOTLIB:
        print("Generating Plots...")
        kept_set = set(kept_barcodes)
        plot_gene_umi_counts_by_type(
            stats_exon=stats_exon,
            stats_intron=stats_intron,
            stats_inex=stats_inex,
            wells=all_wells,
            out_pdf=os.path.join(stats_dir, f"{project}.geneUMIcounts.pdf"),
        )
        if os.path.exists(bam_file):
            plot_features(
                read_stats=read_stats,
                wells=kept_set,
                out_pdf=os.path.join(stats_dir, f"{project}.features.pdf"),
            )
            
    print("Statistics generation finished.")

if __name__ == "__main__":
    main()
