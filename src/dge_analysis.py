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
import multiprocessing

try:
    import pysam
except ImportError:
    print("Error: 'pysam' module is required. Please install it via pip install pysam")
    sys.exit(1)

def load_config(yaml_file):
    with open(yaml_file, 'r') as f:
        return yaml.safe_load(f)

def cluster_umis(umis, threshold=1):
    """
    Optimized UMI clustering using hash-based neighbor search.
    Complexity: O(N * L) instead of O(N^2).
    """
    if not umis: return {}
    
    # 1. Frequency counting
    counts = Counter(umis)
    
    # 2. Sort unique UMIs: descending count, then lexicographical
    # This ensures we always map lower-count/variant UMIs to higher-count/canonical ones
    unique_umis = sorted(counts.keys(), key=lambda x: (-counts[x], x))
    
    parent_map = {u: u for u in unique_umis}
    
    # Optimization: If threshold is 0, no clustering needed
    if threshold == 0:
        return parent_map
        
    # Set for fast lookup of valid UMIs
    umi_set = set(unique_umis)
    
    # We only mark children as 'visited' so they don't become parents
    # But in the sorted loop, once a child is merged, we can skip processing it as a parent?
    # Yes. 'visited' tracks UMIs that have been assigned to a parent.
    visited = set()

    bases = {'A', 'C', 'G', 'T', 'N'}

    for parent in unique_umis:
        if parent in visited:
            continue
            
        # Parent is valid. Now look for children.
        # Instead of iterating all other UMIs, generate neighbors.
        
        parent_len = len(parent)
        parent_count = counts[parent]
        
        # Generate all 1-mismatch neighbors
        # For L=10, there are 30 neighbors (excluding Ns logic, say max 40)
        neighbors = set()
        p_chars = list(parent)
        
        for i in range(parent_len):
            orig = p_chars[i]
            for b in bases:
                if b == orig: continue
                p_chars[i] = b
                neighbors.add("".join(p_chars))
            p_chars[i] = orig # Restore
            
        # Check which neighbors are in our dataset
        for child in neighbors:
            if child in umi_set and child not in visited:
                child_count = counts[child]
                
                # Logic: directional adjacency
                # Usually we merge if parent count >= (2 * child count) - 1
                # But here we use zUMIs-like simple logic: merge if parent count >= child count
                # Since we sorted by count, parent is guaranteed to be >= child (mostly), 
                # or equal but lexicographically smaller.
                
                # To follow strict "directional" method from UMI-tools:
                # if parent_count >= (2 * child_count) - 1:
                
                # Using simple logic consistent with previous implementation:
                # Merge if valid neighbor found. The sorting priority handles the "best parent" selection.
                
                parent_map[child] = parent
                visited.add(child)

    return parent_map

def process_barcode_worker(args):
    """
    Worker function for parallel UMI clustering.
    Args:
        args: tuple (bc, umis_exon_map, umis_intron_map, ham_dist)
    Returns:
        bc, res_umi_counts_exon, res_umi_counts_intron, res_umi_counts_inex, res_correction
    """
    bc, umis_exon_map, umis_intron_map, ham_dist = args
    
    res_umi_counts_exon = {}
    res_umi_counts_intron = {}
    res_umi_counts_inex = {}
    res_correction = {}
    
    genes_exon = set(umis_exon_map.keys())
    genes_intron = set(umis_intron_map.keys())
    all_genes = genes_exon | genes_intron
    
    for gene in all_genes:
        umis_ex = umis_exon_map.get(gene, [])
        umis_in = umis_intron_map.get(gene, [])
        umis_total = umis_ex + umis_in
        
        if not umis_total: continue
        
        # Cluster
        mapping = cluster_umis(umis_total, threshold=ham_dist)
        
        if ham_dist > 0:
            res_correction[gene] = mapping
            
        # Count
        unique_total = set([mapping[u] for u in umis_total])
        res_umi_counts_inex[gene] = len(unique_total)
        
        if umis_ex:
            unique_ex = set([mapping[u] for u in umis_ex])
            res_umi_counts_exon[gene] = len(unique_ex)
            
        if umis_in:
            unique_in = set([mapping[u] for u in umis_in])
            res_umi_counts_intron[gene] = len(unique_in)
            
    return bc, res_umi_counts_exon, res_umi_counts_intron, res_umi_counts_inex, res_correction

def natural_sort_key(s):
    """
    Key for natural sorting (e.g., A2 < A10).
    Splits string into mixed list of strings and integers.
    """
    import re
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split('([0-9]+)', s)]

def load_barcodes(out_dir, project):
    """
    Loads barcodes defining the matrix columns.
    Priority 1: expect_id_barcode.tsv (Well IDs) - Ensures all wells are present & consistent.
    Priority 2: kept_barcodes.txt (Detected Seqs) - Fallback for droplet/unstructured data.
    """
    # 1. Try Expect ID File (Config Dir)
    # out_dir is usually ".../analysis"
    # config is ".../config"
    config_dir = os.path.join(os.path.dirname(out_dir.rstrip('/')), "config")
    expect_file = os.path.join(config_dir, "expect_id_barcode.tsv")
    
    barcodes = []
    source = ""
    
    if os.path.exists(expect_file):
        print(f"Loading reference barcodes (Well IDs) from {expect_file}...")
        source = "expect"
        with open(expect_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                # Skip header if present (check commonly used header names)
                if parts[0] in ['wellID', 'WellID', 'CellID', 'Barcode']:
                    continue
                if parts[0]:
                    barcodes.append(parts[0])
    else:
        # 2. Fallback to Kept Barcodes (Analysis Dir)
        kept_file = os.path.join(out_dir, "zUMIs_output", f"{project}kept_barcodes.txt")
        print(f"Loading reference barcodes (Detected) from {kept_file}...")
        source = "kept"
        if os.path.exists(kept_file):
            with open(kept_file, 'r') as f:
                # Check header
                first = f.readline()
                if not ('XC' in first or 'n' in first):
                    p = first.replace(',', '\t').split('\t')
                    if p[0]: barcodes.append(p[0])
                
                for line in f:
                    p = line.replace(',', '\t').split('\t')
                    if p[0]: barcodes.append(p[0])
        else:
            print("Warning: No barcode file found. Matrix will be empty.")
            return [], set()

    # Remove duplicates just in case
    barcodes = list(set(barcodes))
    
    # Sort
    # Use natural sort for consistency (e.g. P1A2 before P1A10)
    barcodes.sort(key=natural_sort_key)
    
    print(f"Loaded {len(barcodes)} barcodes from {source}.")
    return barcodes, set(barcodes)

def load_genes_from_gtf(gtf_file):
    print(f"Loading reference genes from {gtf_file}...")
    gene_order = []
    gene_map = {}
    seen_genes = set()
    
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            if parts[2] != 'exon': continue # Use exons to find gene entries
            
            attributes = parts[8]
            gene_id = None
            if 'gene_id "' in attributes:
                gene_id = attributes.split('gene_id "')[1].split('"')[0]
            elif 'gene_id' in attributes: # Fallback
                try: gene_id = attributes.split('gene_id')[1].strip().split(';')[0].strip('"')
                except: pass
            
            if not gene_id or gene_id in seen_genes: continue
            
            gene_name = gene_id
            if 'gene_name "' in attributes:
                gene_name = attributes.split('gene_name "')[1].split('"')[0]
            
            seen_genes.add(gene_id)
            gene_order.append(gene_id)
            gene_map[gene_id] = gene_name
            
    return gene_order, gene_map

def write_sparse_matrix(counts_dict, gene_list, gene_names_map, barcode_list, out_dir, subdir_name):
    full_out_dir = os.path.join(out_dir, "zUMIs_output", "expression", subdir_name)
    print(f"Generating deterministic, sorted sparse matrix in {full_out_dir}...")
    
    if not os.path.exists(full_out_dir):
        os.makedirs(full_out_dir)
    
    # Use fixed lists for indices
    gene_to_idx = {g: i for i, g in enumerate(gene_list)}
    bc_to_idx = {b: i for i, b in enumerate(barcode_list)}
    
    triplets = []
    
    # Collect data
    for bc, gene_counts in counts_dict.items():
        if bc not in bc_to_idx: continue
        bc_idx = bc_to_idx[bc]
        
        for gene, count in gene_counts.items():
            if count > 0 and gene in gene_to_idx:
                gene_idx = gene_to_idx[gene]
                # Store as (bc_idx, gene_idx, count) for sorting
                triplets.append((bc_idx, gene_idx, count))
                
    # Sort: Primary = Barcode index (Col), Secondary = Gene index (Row)
    # This produces the requested "Barcode then Feature" order in the file
    triplets.sort(key=lambda x: (x[0], x[1]))
    
    # Re-extract for matrix creation
    # Note: coo_matrix expects (data, (row_indices, col_indices))
    rows = [t[1] for t in triplets]
    cols = [t[0] for t in triplets]
    data = [t[2] for t in triplets]
                
    # Shape is fixed to the full reference lists
    mat = scipy.sparse.coo_matrix((data, (rows, cols)), shape=(len(gene_list), len(barcode_list)))
    
    # Optimized: Write directly to gzip stream to avoid temp file I/O
    matrix_file_gz = os.path.join(full_out_dir, "matrix.mtx.gz")
    
    # Manually write MatrixMarket format to gzip stream
    # Header: %%MatrixMarket matrix coordinate integer general
    # Size: Rows Cols Entries
    with gzip.open(matrix_file_gz, 'wt') as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n")
        f.write("%\n")
        f.write(f"{mat.shape[0]} {mat.shape[1]} {mat.nnz}\n")
        
        # Data: row col val (1-based)
        # Zip iterators are faster than index access
        for r, c, d in zip(mat.row, mat.col, mat.data):
            f.write(f"{r+1} {c+1} {d}\n")

    with gzip.open(os.path.join(full_out_dir, "barcodes.tsv.gz"), "wt") as f:
        for bc in barcode_list:
            f.write(f"{bc}\n")
            
    with gzip.open(os.path.join(full_out_dir, "features.tsv.gz"), "wt") as f:
        for g_id in gene_list:
            g_name = gene_names_map.get(g_id, g_id)
            f.write(f"{g_id}\t{g_name}\tGene Expression\n")

def process_bam_and_matrix(bam_file, out_bam, config, threads):
    project = config['project']
    out_dir = config['out_dir']
    ham_dist = int(config['counting_opts'].get('Ham_Dist', 0))
    count_introns = config.get('counting_opts', {}).get('introns', True)
    
    # Load Reference Lists
    gtf_file = os.path.join(out_dir, f"{project}.final_annot.gtf")
    
    if not os.path.exists(gtf_file):
        raise FileNotFoundError(f"GTF file not found: {gtf_file}")
    
    # Updated barcode loading logic
    barcode_list, barcode_set = load_barcodes(out_dir, project)
        
    gene_list, gene_names_ref = load_genes_from_gtf(gtf_file)
    gene_set = set(gene_list)
    
    print(f"Reference: {len(barcode_list)} Barcodes (Cols), {len(gene_list)} Genes (Rows)")
    
    print(f"Processing {bam_file}")
    
    # Structure: data[category][barcode][gene] = list_of_umis
    umi_data = {
        'exon': defaultdict(lambda: defaultdict(list)),
        'intron': defaultdict(lambda: defaultdict(list))
    }
    
    # Structure: read_counts_raw[category][barcode][gene] = int
    read_counts_raw = {
        'exon': defaultdict(lambda: defaultdict(int)),
        'intron': defaultdict(lambda: defaultdict(int))
    }
    
    correction_map = defaultdict(lambda: defaultdict(dict))
    
    print("Pass 1: Reading BAM to aggregate UMIs and Reads...")
    with pysam.AlignmentFile(bam_file, "rb", threads=threads) as bam:
        for read in bam:
            if read.is_unmapped: continue
            try:
                bc = read.get_tag("CB")
                gene_id = read.get_tag("GX")
            except KeyError:
                continue 
            
            # Filter strictly by Reference Lists
            if bc not in barcode_set: continue
            if gene_id not in gene_set: continue
            
            ftype = "exon"
            if read.has_tag("RE"):
                xf = read.get_tag("RE")
                if xf == "N": ftype = "intron"
                elif xf == "E": ftype = "exon"
                else: continue
            
            if not count_introns and ftype == 'intron': continue
            
            # Increment raw read count
            read_counts_raw[ftype][bc][gene_id] += 1
            
            # Store UMI
            if read.has_tag("UR"):
                umi_data[ftype][bc][gene_id].append(read.get_tag("UR"))

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

    # Calculate Inex Read Counts
    for bc in read_counts_raw['exon']:
        for gene in read_counts_raw['exon'][bc]:
            final_read_counts['inex'][bc][gene] += read_counts_raw['exon'][bc][gene]
            
    for bc in read_counts_raw['intron']:
        for gene in read_counts_raw['intron'][bc]:
            final_read_counts['inex'][bc][gene] += read_counts_raw['intron'][bc][gene]

    # Calculate UMI Counts (with Clustering)
    all_bcs_umis = list(set(umi_data['exon'].keys()) | set(umi_data['intron'].keys()))
    total_bcs = len(all_bcs_umis)
    
    print(f"Clustering UMIs for {total_bcs} barcodes using {threads} threads...")
    
    pool_args = []
    for bc in all_bcs_umis:
        pool_args.append((bc, umi_data['exon'][bc], umi_data['intron'][bc], ham_dist))
        
    if threads > 1:
        # Parallel Processing
        with multiprocessing.Pool(threads) as pool:
            # imap_unordered yields results as they complete
            for i, res in enumerate(pool.imap_unordered(process_barcode_worker, pool_args, chunksize=20)):
                if i % 100 == 0: 
                    print(f"Clustering Barcode {i}/{total_bcs}...", end='\r')
                
                bc, c_ex, c_in, c_inex, corr = res
                
                if c_ex: final_umi_counts['exon'][bc].update(c_ex)
                if c_in: final_umi_counts['intron'][bc].update(c_in)
                if c_inex: final_umi_counts['inex'][bc].update(c_inex)
                if corr:
                    correction_map[bc].update(corr)
    else:
        # Serial Processing
        for i, args in enumerate(pool_args):
            if i % 100 == 0: 
                print(f"Clustering Barcode {i}/{total_bcs}...", end='\r')
            
            bc, c_ex, c_in, c_inex, corr = process_barcode_worker(args)
            
            if c_ex: final_umi_counts['exon'][bc].update(c_ex)
            if c_in: final_umi_counts['intron'][bc].update(c_in)
            if c_inex: final_umi_counts['inex'][bc].update(c_inex)
            if corr:
                correction_map[bc].update(corr)

    print("\nWriting Matrices...")
    
    # Write UMI Counts
    write_sparse_matrix(final_umi_counts['exon'], gene_list, gene_names_ref, barcode_list, out_dir, f"{project}.exon.umi")
    if count_introns:
        write_sparse_matrix(final_umi_counts['intron'], gene_list, gene_names_ref, barcode_list, out_dir, f"{project}.intron.umi")
        write_sparse_matrix(final_umi_counts['inex'], gene_list, gene_names_ref, barcode_list, out_dir, f"{project}.inex.umi")
        
    # Write Read Counts
    write_sparse_matrix(final_read_counts['exon'], gene_list, gene_names_ref, barcode_list, out_dir, f"{project}.exon.read")
    if count_introns:
        write_sparse_matrix(final_read_counts['intron'], gene_list, gene_names_ref, barcode_list, out_dir, f"{project}.intron.read")
        write_sparse_matrix(final_read_counts['inex'], gene_list, gene_names_ref, barcode_list, out_dir, f"{project}.inex.read")
    
    # --- PASS 2: BAM Correction & Streaming Sort ---
    print("Pass 2: Correcting BAM and streaming to samtools sort...")
    
    # Construct samtools sort command
    # Sorts from stdin (-) to out_bam
    sort_cmd = [
        sys.argv[2], # samtools_exec passed from main
        "sort",
        "-@", str(threads),
        "-T", out_bam + ".sort_tmp",
        "-o", out_bam,
        "-"
    ]
    
    sort_process = subprocess.Popen(sort_cmd, stdin=subprocess.PIPE, bufsize=10*1024*1024)
    
    try:
        with pysam.AlignmentFile(bam_file, "rb", threads=threads) as infile:
            # Open a pysam writer pointing to the sort process's stdin
            # We use "wb" mode and pass the process's stdin as the file
            # template=infile ensures header is copied
            with pysam.AlignmentFile(sort_process.stdin, "wb", template=infile, threads=threads) as outfile:
                for read in infile:
                    # Process UMI tags if UR exists
                    if read.has_tag("UR"):
                        raw_umi = read.get_tag("UR")
                        
                        # Check if assigned to gene
                        if read.has_tag("CB") and read.has_tag("GX"):
                            bc = read.get_tag("CB")
                            gene = read.get_tag("GX")
                            
                            final_umi = raw_umi
                            # Apply correction
                            if ham_dist > 0:
                                # Safe dictionary lookup without explicit existence check
                                # correction_map is defaultdict(lambda: defaultdict(dict))
                                # Using .get() on standard dicts is safer, but for defaultdict we must be careful
                                # Accessing correction_map[bc] creates it if missing, which is bad for memory if bc is junk.
                                # Use get() on the defaultdict itself to avoid creation.
                                bc_map = correction_map.get(bc)
                                if bc_map:
                                    gene_map = bc_map.get(gene)
                                    if gene_map:
                                        # Use get() for final lookup, default to raw_umi
                                        final_umi = gene_map.get(raw_umi, raw_umi)
                            
                            read.set_tag("UB", final_umi)
                        else:
                            read.set_tag("UB", None)
                    
                    outfile.write(read)
                    
    except Exception as e:
        print(f"Error during BAM streaming: {e}")
        sort_process.kill()
        raise
    finally:
        # Close stdin to signal EOF to samtools sort
        if sort_process.stdin:
            sort_process.stdin.close()
        
        # Wait for sort to finish
        ret = sort_process.wait()
        if ret != 0:
            raise RuntimeError(f"samtools sort failed with return code {ret}")
    
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
