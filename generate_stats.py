#!/usr/bin/env python3
import sys
import os
import yaml
import subprocess
import csv
import collections
import statistics
import math

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

def calculate_gene_stats(csv_file):
    """
    Reads dgecounts.csv/readcounts.csv and calculates:
    - Number of genes detected per cell.
    - Total reads/UMIs per gene.
    """
    # Format: Gene,Barcode,Count
    
    genes_per_cell = collections.defaultdict(int)
    counts_per_gene = collections.defaultdict(int)
    
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            bc = row['Barcode']
            gene = row['Gene']
            count = int(row['Count'])
            
            if count > 0:
                genes_per_cell[bc] += 1
                counts_per_gene[gene] += count
                
    return genes_per_cell, counts_per_gene

def parse_bam_stats(bam_file, samtools_exec, kept_barcodes):
    """
    Iterates BAM and counts read assignments (Exon, Intron, etc.) per cell.
    Uses XF tag if available, else infers from XS/XT.
    """
    print(f"Calculating Mapping Stats from {bam_file}...")
    
    # Stats: {Barcode: {Category: Count}}
    stats = collections.defaultdict(lambda: collections.defaultdict(int))
    
    # Categories: Exon, Intron, Intergenic, Unmapped, Ambiguity
    
    cmd = [samtools_exec, 'view', bam_file]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1)
    
    count = 0
    try:
        for line in proc.stdout:
            count += 1
            if count % 1000000 == 0: print(f"Processed {count} reads...", end='\r')
            
            parts = line.strip().split('\t')
            if len(parts) < 12: continue
            
            # Extract Tags
            tags = {}
            for t in parts[11:]:
                if ':' in t:
                    k = t.split(':')[0]
                    v = t.split(':')[2]
                    tags[k] = v
                    
            bc = tags.get('BC')
            if not bc:
                # Try to extract from read name if not in tag? zUMIs usually has BC tag.
                continue
            
            # Determine Category
            category = "Intergenic"
            
            # Use XF tag if I added it
            if 'XF' in tags:
                category = tags['XF']
            else:
                # Fallback Logic
                if 'XT' in tags:
                    category = "Exon" # Assume exon if assigned and no XF
                elif 'XS' in tags:
                    status = tags['XS']
                    if "Unassigned_NoFeatures" in status:
                        category = "Intergenic"
                    elif "Unassigned_Ambiguity" in status:
                        category = "Ambiguity"
                    elif "Unassigned_Unmapped" in status:
                        category = "Unmapped"
                    else:
                        category = status
                else:
                    # Check flag for unmapped
                    flag = int(parts[1])
                    if flag & 0x4:
                        category = "Unmapped"
            
            # Check if BC is valid
            if bc not in kept_barcodes:
                # Aggregate to "bad" or keep specific? zUMIs aggregates bad BCs.
                # But stats usually show Top BCs?
                # zUMIs stats shows "Unused BC" category in total barplot.
                # Here we just mark the BC as "Unused" for the key?
                # Actually, zUMIs accumulates stats for ALL reads.
                # If BC is not in kept list, it counts towards "Unused BC" totals?
                # Let's count it under the specific BC first, then aggregate later.
                pass

            stats[bc][category] += 1
            
    finally:
        proc.stdout.close()
        proc.wait()
        print(f"\nFinished parsing {count} reads.")
        
    return stats

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
    
    # Inputs
    expr_dir = os.path.join(out_dir, "zUMIs_output", "expression", project)
    dge_file = f"{expr_dir}.dgecounts.csv"
    read_file = f"{expr_dir}.readcounts.csv"
    
    # Valid Barcodes
    bc_file = os.path.join(out_dir, "zUMIs_output", f"{project}kept_barcodes.txt")
    kept_barcodes = set()
    if os.path.exists(bc_file):
        with open(bc_file, 'r') as f:
            # Skip header if present? zUMIs kept_barcodes usually has header XC,n...
            # Check first line
            header = f.readline()
            if 'XC' in header:
                pass
            else:
                kept_barcodes.add(header.strip().split(',')[0]) # Assuming CSV? Or TSV?
                # zUMIs output is usually CSV or TSV.
            
            for line in f:
                parts = line.strip().split(',')
                if parts[0]: kept_barcodes.add(parts[0])
    
    # 1. Gene Stats
    print("Calculating Gene Statistics...")
    genes_per_cell_umi = {}
    genes_per_cell_read = {}
    
    if os.path.exists(dge_file):
        genes_per_cell_umi, counts_per_gene_umi = calculate_gene_stats(dge_file)
        # Write gene counts
        with open(os.path.join(stats_dir, f"{project}.genecounts.txt"), 'w') as f:
            f.write("GeneID\tCount\n")
            for gene, count in counts_per_gene_umi.items():
                f.write(f"{gene}\t{count}\n")
    
    if os.path.exists(read_file):
        genes_per_cell_read, counts_per_gene_read = calculate_gene_stats(read_file)
        
    # 2. Mapping Stats
    bam_file = os.path.join(out_dir, f"{project}.filtered.Aligned.GeneTagged.bam")
    if os.path.exists(bam_file):
        stats = parse_bam_stats(bam_file, samtools, kept_barcodes)
        
        # Write Reads per Cell
        with open(os.path.join(stats_dir, f"{project}.readspercell.txt"), 'w') as f:
            f.write("Barcode\tExon\tIntron\tIntergenic\tUnmapped\tAmbiguity\n")
            for bc, counts in stats.items():
                if bc in kept_barcodes:
                    f.write(f"{bc}\t{counts['Exon']}\t{counts['Intron']}\t{counts['Intergenic']}\t{counts['Unmapped']}\t{counts['Ambiguity']}\n")
                    
        # 3. Plots
        if HAS_MATPLOTLIB:
            print("Generating Plots...")
            # Gene Counts Plot
            plot_gene_counts(genes_per_cell_umi, genes_per_cell_read, os.path.join(stats_dir, f"{project}.geneUMIcounts.pdf"))
            
            # Features Plot (Summary)
            # Create a stacked bar chart of mapping stats for top BCs?
            # zUMIs makes "features.pdf"
            
            fig, ax = plt.subplots(figsize=(10, 6))
            
            # Aggregate totals for valid BCs
            totals = collections.defaultdict(int)
            for bc, counts in stats.items():
                if bc in kept_barcodes:
                    for cat, val in counts.items():
                        totals[cat] += val
                else:
                    # Count as Unused
                    totals["Unused BC"] += sum(counts.values())
                    
            categories = list(totals.keys())
            values = list(totals.values())
            
            ax.bar(categories, values)
            ax.set_title("Total Read Distribution")
            ax.set_ylabel("Number of Reads")
            plt.xticks(rotation=45)
            plt.tight_layout()
            
            plt.savefig(os.path.join(stats_dir, f"{project}.features.pdf"))
            plt.close()
            
    print("Statistics generation finished.")

if __name__ == "__main__":
    main()
