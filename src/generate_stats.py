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
import json

import numpy as np
import bisect

try:
    import pysam
except ImportError:
    # Not required if loading JSON stats
    pass

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

def load_config(yaml_file):
    with open(yaml_file, 'r') as f:
        return yaml.safe_load(f)

def load_barcode_mapping(out_dir, project):
    mapping = {}
    expect_candidates = [
        os.path.join(out_dir, "config", "expect_id_barcode.tsv"),
        os.path.join(out_dir, "expect_id_barcode.tsv"),
        os.path.join(os.path.dirname(out_dir.rstrip('/')), "config", "expect_id_barcode.tsv"),
    ]
    expect_file = next((p for p in expect_candidates if os.path.exists(p)), None)
    
    if not expect_file:
        return mapping

    with open(expect_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                if parts[0] in ['wellID', 'WellID']: continue
                well = parts[0].strip()
                umi_bc = parts[1].strip()
                int_bc = parts[2].strip()
                mapping[well] = {'umi': umi_bc, 'int': int_bc}
    return mapping

def calculate_matrix_stats(matrix_dir):
    """
    Calculates UMI count and Gene count per cell from a sparse matrix directory.
    Returns: dict {bc: {'umis': int, 'genes': int, 'gene_set': set}}
    """
    stats = collections.defaultdict(lambda: {'umis': 0, 'genes': 0, 'gene_set': set()})
    
    matrix_path = os.path.join(matrix_dir, "matrix.mtx.gz")
    barcodes_path = os.path.join(matrix_dir, "barcodes.tsv.gz")
    
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
            
            col_idx = int(parts[1]) - 1
            row_idx = int(parts[0]) - 1 
            val = float(parts[2])
            
            if val > 0:
                bc = barcodes[col_idx]
                stats[bc]['umis'] += int(val)
                stats[bc]['genes'] += 1
                stats[bc]['gene_set'].add(row_idx)
                
    return stats

def plot_coverage(cov_umi, cov_int, out_prefix):
    if not HAS_MATPLOTLIB: return
    
    def smooth(y, window=7):
        y = np.asarray(y, dtype=np.float64)
        if window <= 1: return y
        if window % 2 == 0: window += 1
        if y.size < window: return y
        kernel = np.ones(window, dtype=np.float64) / window
        pad = window // 2
        y_pad = np.pad(y, (pad, pad), mode='edge')
        return np.convolve(y_pad, kernel, mode='valid')

    sum_umi = float(np.sum(cov_umi))
    sum_int = float(np.sum(cov_int))
    # Avoid div by zero
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

def plot_gene_umi_counts_by_type(stats_exon, stats_intron, stats_inex, wells, out_pdf):
    if not HAS_MATPLOTLIB: return
    types = ["Exon", "Intron+Exon", "Intron"]
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
    
    bp0 = axes[0].boxplot(gene_data, labels=types, notch=True, patch_artist=True)
    for patch, color in zip(bp0['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.8)
    axes[0].set_title("Number of Genes")
    axes[0].set_ylabel("Count")
    
    bp1 = axes[1].boxplot(umi_data, labels=types, notch=True, patch_artist=True)
    for patch, color in zip(bp1['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.8)
    axes[1].set_title("Number of UMIs")
    axes[1].set_ylabel("Count")

    plt.tight_layout()
    plt.savefig(out_pdf)
    plt.close()

def plot_features(read_stats, wells, out_pdf):
    if not HAS_MATPLOTLIB: return
    feat_colors = {
        "Exon": "#1A5084", "Intron": "#118730", "Unmapped": "#545454",
        "Ambiguity": "#FFA54F", "MultiMapping": "#631879FF",
        "Intergenic": "#FFD700", "Unused BC": "#BABABA",
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
            for c in categories_order: per_cell_frac[c].append(0.0)

    total_sum = sum(totals.values()) + unused_total
    total_fracs = {c: (totals[c] / total_sum if total_sum > 0 else 0.0) for c in categories_order}
    unused_frac = unused_total / total_sum if total_sum > 0 else 0.0

    bar_categories = [c for c in categories_order if totals.get(c, 0) > 0]
    box_categories = [c for c in categories_order if len(per_cell_frac[c]) > 0]

    fig, axes = plt.subplots(2, 1, figsize=(12, 7), gridspec_kw={"height_ratios": [0.6, 1.0]})

    left = 0.0
    for c in bar_categories:
        w = total_fracs.get(c, 0.0)
        if w > 0:
            axes[0].barh([0], [w * 100.0], left=[left * 100.0], color=feat_colors[c], label=c, height=0.5)
            left += w
    if unused_frac > 0:
        axes[0].barh([0], [unused_frac * 100.0], left=[left * 100.0], color=feat_colors["Unused BC"], label="Unused BC", height=0.5)

    axes[0].set_xlim(0, 100)
    axes[0].set_yticks([])
    axes[0].set_xlabel("% of total reads")
    axes[0].set_title("Total Read Distribution")
    axes[0].legend(ncol=max(1, len(bar_categories)+1), loc="upper center", bbox_to_anchor=(0.5, -0.30))

    box_pairs = []
    for c in box_categories:
        vals = [v * 100.0 for v in per_cell_frac[c]]
        if len(vals) == 0:
            continue
        box_pairs.append((c, vals))

    if not box_pairs or total_sum <= 0:
        axes[1].axis("off")
        axes[1].text(0.5, 0.5, "No reads available for per-cell distribution", ha="center", va="center")
    else:
        plot_labels = [c for c, _vals in box_pairs]
        box_data = [_vals for _c, _vals in box_pairs]
        bp = axes[1].boxplot(box_data, labels=plot_labels, notch=True, patch_artist=True)
        for patch, c in zip(bp["boxes"], plot_labels):
            patch.set_facecolor(feat_colors[c])

    axes[1].set_title("Reads per Cell")
    axes[1].set_ylabel("% reads/cell")
    
    fig.tight_layout()
    plt.savefig(out_pdf)
    plt.close()

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 generate_stats.py <yaml_config>")
        sys.exit(1)
        
    yaml_file = sys.argv[1]
    config = load_config(yaml_file)
    project = config['project']
    out_dir = config['out_dir']
    
    stats_dir = os.path.join(out_dir, "zUMIs_output", "stats")
    if not os.path.exists(stats_dir): os.makedirs(stats_dir)
    
    # 1. Load Barcode Mapping
    bc_mapping = load_barcode_mapping(out_dir, project)
    all_wells = sorted(list(bc_mapping.keys()))
    kept_barcodes = set(all_wells)
    
    # 2. Calculate Matrix Stats (Fast, independent)
    print("Calculating Matrix Statistics...")
    umi_exon_dir = os.path.join(out_dir, "zUMIs_output", "expression", f"{project}.exon.umi")
    umi_intron_dir = os.path.join(out_dir, "zUMIs_output", "expression", f"{project}.intron.umi")
    
    stats_exon = calculate_matrix_stats(umi_exon_dir)
    stats_intron = calculate_matrix_stats(umi_intron_dir)
    
    stats_inex = collections.defaultdict(lambda: {'umis': 0, 'genes': 0})
    for bc in set(stats_exon.keys()) | set(stats_intron.keys()):
        ex_st = stats_exon.get(bc, {'umis': 0, 'genes': 0, 'gene_set': set()})
        in_st = stats_intron.get(bc, {'umis': 0, 'genes': 0, 'gene_set': set()})
        stats_inex[bc]['umis'] = ex_st['umis'] + in_st['umis']
        stats_inex[bc]['genes'] = len(ex_st['gene_set'] | in_st['gene_set'])

    # 3. Load Pre-calculated Read Stats & Coverage
    stats_json_path = os.path.join(out_dir, "zUMIs_output", "stats", f"{project}.read_stats.json")
    
    read_stats = collections.defaultdict(lambda: collections.defaultdict(int))
    cov_umi = np.zeros(100, dtype=np.int64)
    cov_int = np.zeros(100, dtype=np.int64)
    
    if os.path.exists(stats_json_path):
        print(f"Loading pre-calculated stats from {stats_json_path}")
        with open(stats_json_path, 'r') as f:
            data = json.load(f)
            # Reconstruct defaultdict from JSON dict
            raw_stats = data.get("read_stats", {})
            for k, v in raw_stats.items():
                read_stats[k] = v
                
            cov_umi_list = data.get("coverage_umi", [])
            cov_int_list = data.get("coverage_int", [])
            
            if cov_umi_list:
                cov_umi = np.array(cov_umi_list, dtype=np.int64)
            if cov_int_list:
                cov_int = np.array(cov_int_list, dtype=np.int64)
    else:
        print("Warning: Pre-calculated stats not found. Plots will be missing.")

    # 4. Plots and Tables
    out_prefix = os.path.join(stats_dir, project)
    plot_coverage(cov_umi, cov_int, out_prefix)
    
    # Write Table
    output_table = os.path.join(stats_dir, f"{project}.stats.tsv")
    print(f"Writing stats table to {output_table}...")
    
    with open(output_table, 'w') as f:
        headers = [
            "wellID", "internal_barcodes", "umi_barcodes", 
            "internal_reads", "umi_reads", 
            "Ambiguity_reads", "Exon_reads", "Intergenic_reads", "intron_reads", 
            "Unmapped_reads", 
            "Exon_umis", "Intron_umis", "Intron_Exon_umis",
            "Exon_genes", "Intron_genes", "Intron_Exon_genes"
        ]
        f.write("\t".join(headers) + "\n")
        
        for well in all_wells:
            int_bc = bc_mapping[well]['int']
            umi_bc = bc_mapping[well]['umi']
            
            r_stats = read_stats.get(well, {})
            ambig = r_stats.get('Ambiguity', 0)
            exon_r = r_stats.get('Exon', 0)
            inter_r = r_stats.get('Intergenic', 0)
            intron_r = r_stats.get('Intron', 0)
            unmapped = r_stats.get('Unmapped', 0)

            internal_r = r_stats.get('Internal_Reads', 0)
            umi_r = r_stats.get('UMI_Reads', 0)
            
            ex_st = stats_exon.get(well, {'umis': 0, 'genes': 0})
            in_st = stats_intron.get(well, {'umis': 0, 'genes': 0})
            inex_st = stats_inex.get(well, {'umis': 0, 'genes': 0})
            
            row = [
                well, int_bc, umi_bc,
                internal_r, umi_r,
                ambig, exon_r, inter_r, intron_r,
                unmapped, 
                ex_st['umis'], in_st['umis'], inex_st['umis'],
                ex_st['genes'], in_st['genes'], inex_st['genes']
            ]
            f.write("\t".join(map(str, row)) + "\n")

    if HAS_MATPLOTLIB:
        plot_gene_umi_counts_by_type(stats_exon, stats_intron, stats_inex, all_wells, os.path.join(stats_dir, f"{project}.geneUMIcounts.pdf"))
        plot_features(read_stats, kept_barcodes, os.path.join(stats_dir, f"{project}.features.pdf"))
            
    print("Statistics generation finished.")

if __name__ == "__main__":
    main()
