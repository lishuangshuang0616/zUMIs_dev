#!/usr/bin/env python3
import sys
import os
import yaml
import subprocess
import csv
import gzip
import shutil
from collections import defaultdict, Counter
import scipy.io
import scipy.sparse
import numpy as np
import pandas as pd

try:
    import pysam
except ImportError:
    print("Error: 'pysam' module is required. Please install it via pip install pysam")
    sys.exit(1)

def load_config(yaml_file):
    with open(yaml_file, 'r') as f:
        return yaml.safe_load(f)

def hamming_distance(s1, s2):
    if len(s1) != len(s2): return len(s1)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def cluster_umis(umis, threshold=1):
    if not umis: return {}
    if threshold == 0:
        return {u: u for u in umis}

    counts = Counter(umis)
    unique_umis = sorted(counts.keys(), key=lambda x: (-counts[x], x))
    
    parent_map = {} 
    visited = set()
    
    for parent in unique_umis:
        if parent in visited: continue
        parent_map[parent] = parent
        visited.add(parent)
        
        children = []
        for candidate in unique_umis:
            if candidate in visited: continue
            if hamming_distance(parent, candidate) <= threshold:
                children.append(candidate)
        
        for child in children:
            parent_map[child] = parent
            visited.add(child)
            
    return parent_map

def write_sparse_matrix(counts_dict, gene_names_map, out_dir, subdir_name):
    full_out_dir = os.path.join(out_dir, "zUMIs_output", "expression", subdir_name)
    print(f"Generating sparse matrix in {full_out_dir}...")
    
    if not os.path.exists(full_out_dir):
        os.makedirs(full_out_dir)
    
    barcodes = sorted(counts_dict.keys())
    genes_set = set()
    for bc in counts_dict:
        genes_set.update(counts_dict[bc].keys())
    genes = sorted(list(genes_set))
    
    gene_to_idx = {g: i for i, g in enumerate(genes)}
    bc_to_idx = {b: i for i, b in enumerate(barcodes)}
    
    rows = []
    cols = []
    data = []
    
    for bc, gene_counts in counts_dict.items():
        bc_idx = bc_to_idx[bc]
        for gene, count in gene_counts.items():
            if count > 0:
                gene_idx = gene_to_idx[gene]
                rows.append(gene_idx)
                cols.append(bc_idx)
                data.append(count)
                
    mat = scipy.sparse.coo_matrix((data, (rows, cols)), shape=(len(genes), len(barcodes)))
    
    matrix_file = os.path.join(full_out_dir, "matrix.mtx")
    scipy.io.mmwrite(matrix_file, mat)
    with open(matrix_file, 'rb') as f_in:
        with gzip.open(matrix_file + ".gz", 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(matrix_file)
    
    with gzip.open(os.path.join(full_out_dir, "barcodes.tsv.gz"), "wt") as f:
        for bc in barcodes:
            f.write(f"{bc}\n")
            
    with gzip.open(os.path.join(full_out_dir, "features.tsv.gz"), "wt") as f:
        for g_id in genes:
            g_name = gene_names_map.get(g_id, g_id)
            f.write(f"{g_id}\t{g_name}\tGene Expression\n")

def process_bam_and_matrix(bam_file, out_bam, config, threads):
    project = config['project']
    out_dir = config['out_dir']
    ham_dist = int(config['counting_opts'].get('Ham_Dist', 0))
    count_introns = config.get('counting_opts', {}).get('introns', True)
    
    print(f"Processing {bam_file}")
    
    # Structure: data[category][barcode][gene] = list_of_umis
    umi_data = {
        'exon': defaultdict(lambda: defaultdict(list)),
        'intron': defaultdict(lambda: defaultdict(list))
    }
    
    # Structure: read_counts_raw[category][barcode][gene] = int
    # Counts raw reads BEFORE UMI collapsing
    read_counts_raw = {
        'exon': defaultdict(lambda: defaultdict(int)),
        'intron': defaultdict(lambda: defaultdict(int))
    }
    
    gene_names_map = {} 
    correction_map = defaultdict(lambda: defaultdict(dict))
    
    print("Pass 1: Reading BAM to aggregate UMIs and Reads...")
    with pysam.AlignmentFile(bam_file, "rb", threads=threads) as bam:
        for read in bam:
            if read.is_unmapped: continue
            try:
                bc = read.get_tag("BC")
                gene_id = read.get_tag("XT")
            except KeyError:
                continue 
            
            if read.has_tag("GN"):
                gene_names_map[gene_id] = read.get_tag("GN")
            elif gene_id not in gene_names_map:
                gene_names_map[gene_id] = gene_id
            
            ftype = "exon"
            if read.has_tag("XF"):
                xf = read.get_tag("XF")
                if xf == "Intron": ftype = "intron"
                elif xf == "Exon": ftype = "exon"
                else: continue
            
            if not count_introns and ftype == 'intron': continue
            
            # Increment raw read count
            read_counts_raw[ftype][bc][gene_id] += 1
            
            # Store UMI
            if read.has_tag("UB"):
                umi_data[ftype][bc][gene_id].append(read.get_tag("UB"))

    print("Pass 1 Complete. Calculating Statistics...")

    # Final count containers
    final_umi_counts = {
        'exon': defaultdict(lambda: defaultdict(int)),
        'intron': defaultdict(lambda: defaultdict(int)),
        'inex': defaultdict(lambda: defaultdict(int))
    }
    
    final_read_counts = {
        'exon': read_counts_raw['exon'],
        'intron': read_counts_raw['intron'],
        'inex': defaultdict(lambda: defaultdict(int))
    }

    # Calculate Inex Read Counts (Sum of Exon + Intron)
    # Reads are mutually exclusive in BAM (a read is either Exon OR Intron), so simple sum works.
    all_bcs_reads = set(read_counts_raw['exon'].keys()) | set(read_counts_raw['intron'].keys())
    for bc in all_bcs_reads:
        genes_exon = set(read_counts_raw['exon'].get(bc, {}).keys())
        genes_intron = set(read_counts_raw['intron'].get(bc, {}).keys())
        all_genes = genes_exon | genes_intron
        for gene in all_genes:
            c_ex = read_counts_raw['exon'].get(bc, {}).get(gene, 0)
            c_in = read_counts_raw['intron'].get(bc, {}).get(gene, 0)
            final_read_counts['inex'][bc][gene] = c_ex + c_in

    # Calculate UMI Counts (with Clustering)
    all_bcs_umis = set(umi_data['exon'].keys()) | set(umi_data['intron'].keys())
    total_bcs = len(all_bcs_umis)
    
    for i, bc in enumerate(all_bcs_umis):
        if i % 100 == 0: print(f"Clustering Barcode {i}/{total_bcs}...", end='\r')
        
        genes_exon = set(umi_data['exon'].get(bc, {}).keys())
        genes_intron = set(umi_data['intron'].get(bc, {}).keys())
        all_genes = genes_exon | genes_intron
        
        for gene in all_genes:
            umis_ex = umi_data['exon'].get(bc, {}).get(gene, [])
            umis_in = umi_data['intron'].get(bc, {}).get(gene, [])
            umis_total = umis_ex + umis_in
            if not umis_total: continue
            
            mapping = cluster_umis(umis_total, threshold=ham_dist)
            if ham_dist > 0:
                correction_map[bc][gene] = mapping
            
            # Count Unique Corrected UMIs
            unique_total = set([mapping[u] for u in umis_total])
            final_umi_counts['inex'][bc][gene] = len(unique_total)
            
            if umis_ex:
                unique_ex = set([mapping[u] for u in umis_ex])
                final_umi_counts['exon'][bc][gene] = len(unique_ex)
            
            if umis_in:
                unique_in = set([mapping[u] for u in umis_in])
                final_umi_counts['intron'][bc][gene] = len(unique_in)

    print("\nWriting Matrices...")
    
    # Write UMI Counts
    write_sparse_matrix(final_umi_counts['exon'], gene_names_map, out_dir, f"{project}.exon.umi")
    if count_introns:
        write_sparse_matrix(final_umi_counts['intron'], gene_names_map, out_dir, f"{project}.intron.umi")
        write_sparse_matrix(final_umi_counts['inex'], gene_names_map, out_dir, f"{project}.inex.umi")
        
    # Write Read Counts
    write_sparse_matrix(final_read_counts['exon'], gene_names_map, out_dir, f"{project}.exon.read")
    if count_introns:
        write_sparse_matrix(final_read_counts['intron'], gene_names_map, out_dir, f"{project}.intron.read")
        write_sparse_matrix(final_read_counts['inex'], gene_names_map, out_dir, f"{project}.inex.read")
    
    # --- PASS 2: BAM Correction ---
    print("Pass 2: Writing corrected BAM...")
    tmp_bam_name = out_bam + ".tmp.bam"
    
    with pysam.AlignmentFile(bam_file, "rb", threads=threads) as infile:
        with pysam.AlignmentFile(tmp_bam_name, "wb", template=infile, threads=threads) as outfile:
            for read in infile:
                # Process UMI tags if UB exists (Apply to ALL runs, Ham_Dist 0 or >0)
                if read.has_tag("UB"):
                    raw_umi = read.get_tag("UB")
                    # Always set UX to raw UMI
                    read.set_tag("UX", raw_umi)
                    
                    # Check if assigned to gene
                    if read.has_tag("BC") and read.has_tag("XT"):
                        bc = read.get_tag("BC")
                        gene = read.get_tag("XT")
                        
                        # Apply correction to UB if available (only relevant if Ham_Dist > 0)
                        if ham_dist > 0:
                            if bc in correction_map and gene in correction_map[bc]:
                                if raw_umi in correction_map[bc][gene]:
                                    corrected_umi = correction_map[bc][gene][raw_umi]
                                    read.set_tag("UB", corrected_umi)
                        # else: UB remains raw_umi (correct behavior for Ham_Dist=0)
                    else:
                        # NOT assigned to gene: Remove UB tag
                        read.set_tag("UB", None)
                
                outfile.write(read)
    
    print(f"Sorting corrected BAM to {out_bam}...")
    pysam.sort("-@", str(threads), "-o", out_bam, tmp_bam_name)
    if os.path.exists(tmp_bam_name): os.remove(tmp_bam_name)
    print("Indexing...")
    pysam.index(out_bam)

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 dge_analysis.py <yaml_config> <samtools_exec>")
        sys.exit(1)
    yaml_file = sys.argv[1]
    config = load_config(yaml_file)
    project = config['project']
    out_dir = config['out_dir']
    num_threads = int(config.get('num_threads', 1))
    input_bam = os.path.join(out_dir, f"{project}.filtered.Aligned.GeneTagged.bam")
    if not os.path.exists(input_bam):
        print(f"Error: Input BAM {input_bam} not found.")
        sys.exit(1)
    final_bam = os.path.join(out_dir, f"{project}.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam")
    process_bam_and_matrix(input_bam, final_bam, config, threads=num_threads)
    print("DGE Analysis pipeline finished.")

if __name__ == "__main__":
    main()
