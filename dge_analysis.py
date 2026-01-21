#!/usr/bin/env python3
import sys
import os
import yaml
import subprocess
import csv
import gzip
from collections import defaultdict
import glob
import math

def load_config(yaml_file):
    with open(yaml_file, 'r') as f:
        return yaml.safe_load(f)

def hamming_distance(s1, s2):
    if len(s1) != len(s2): return len(s1)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def parse_tags(tag_str):
    tags = {}
    for t in tag_str:
        parts = t.split(':')
        if len(parts) >= 3:
            tags[parts[0]] = parts[2]
    return tags

def umi_collapse(umis, threshold=1):
    """
    Simple adjacency-based UMI collapsing.
    Returns count of unique molecules.
    """
    if not umis: return 0
    if threshold == 0: return len(set(umis))
    
    # Sort by abundance (descending) to prioritize common UMIs as 'parents'
    umi_counts = collections.Counter(umis)
    sorted_umis = sorted(umi_counts.keys(), key=lambda x: umi_counts[x], reverse=True)
    
    kept_umis = []
    
    for umi in sorted_umis:
        if not kept_umis:
            kept_umis.append(umi)
            continue
            
        # Check distance to already kept UMIs
        # If close to any kept UMI, merge into it (discard current)
        is_merged = False
        for parent in kept_umis:
            if hamming_distance(umi, parent) <= threshold:
                is_merged = True
                break
        
        if not is_merged:
            kept_umis.append(umi)
            
    return len(kept_umis)

# Optimized version of collapse using graph clustering is better but slower to implement in pure python without networkx.
# Using a simple greedy approach:
# 1. Take most abundant UMI.
# 2. Remove all UMIs within distance d.
# 3. Repeat.
import collections

def greedy_umi_collapse(umis, threshold=1):
    if not umis: return 0
    if threshold == 0: return len(set(umis))
    
    counts = collections.Counter(umis)
    # Sort UMIs by count (desc) then by sequence (for stability)
    sorted_umis = sorted(counts.keys(), key=lambda x: (-counts[x], x))
    
    unique_molecules = 0
    
    # We mark UMIs as visited
    visited = set()
    
    for i, umi in enumerate(sorted_umis):
        if umi in visited:
            continue
            
        unique_molecules += 1
        visited.add(umi)
        
        # Find network of errors
        # To optimize, we only look at other UMIs if we haven't visited them
        for other_umi in sorted_umis[i+1:]:
            if other_umi in visited:
                continue
            
            if hamming_distance(umi, other_umi) <= threshold:
                visited.add(other_umi)
                
    return unique_molecules

def process_bam(bam_file, samtools_exec, ham_dist, threads, out_prefix, is_smart3=False):
    """
    Reads BAM, aggregates counts by (BC, Gene).
    Outputs CSVs.
    """
    print(f"Processing {bam_file} with Hamming Distance {ham_dist}...")
    
    # Data structures:
    # counts[barcode][gene] = list_of_umis
    # For smart-seq3 internal reads, we might just count reads?
    # Based on zUMIs logic: internal reads are also counted.
    
    umi_data = defaultdict(lambda: defaultdict(list))
    read_data = defaultdict(lambda: defaultdict(int))
    internal_data = defaultdict(lambda: defaultdict(int)) # For Smart-seq3 internal
    
    # Open BAM stream
    # Using -F 0x4 (mapped only) is handled by featureCounts usually, but good to be safe.
    # We need tags: BC, UB, XT (gene)
    # If featureCounts used, gene tag might be 'XT'
    
    cmd = [samtools_exec, 'view', '-@', str(threads), bam_file]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1)
    
    line_count = 0
    for line in proc.stdout:
        line_count += 1
        if line_count % 1000000 == 0:
            print(f"Processed {line_count} reads...", end='\r')
            
        parts = line.strip().split('\t')
        if len(parts) < 12: continue
        
        # Parse Tags
        tags = {}
        # Tags start from index 11
        for t in parts[11:]:
            if t.startswith('BC:Z:'): tags['BC'] = t[5:]
            elif t.startswith('UB:Z:'): tags['UB'] = t[5:]
            elif t.startswith('XT:Z:'): tags['XT'] = t[5:] # FeatureCounts gene
            
        bc = tags.get('BC')
        gene = tags.get('XT')
        umi = tags.get('UB')
        
        if not bc or not gene:
            continue
            
        if gene == 'NA': continue # Unassigned
        
        # Check for Smart-seq3 internal reads
        # zUMIs logic: if UMI is missing or pattern not matched?
        # Actually in fqfilter.py we might have set UB to empty if it's internal.
        
        is_internal = False
        if not umi:
             is_internal = True
        
        # If specific Smart-seq3 logic is needed, we check here.
        # Assuming UB present = UMI read. UB absent = Internal read.
        
        read_data[bc][gene] += 1
        
        if not is_internal:
            umi_data[bc][gene].append(umi)
        else:
            if is_smart3:
                internal_data[bc][gene] += 1
                
    print(f"\nFinished reading BAM. Total lines: {line_count}")
    proc.wait()
    if proc.returncode != 0:
        stderr = proc.stderr.read().strip() if proc.stderr else ""
        raise RuntimeError(f"samtools view failed (rc={proc.returncode}) for {bam_file}: {stderr}")
    try:
        if proc.stdout:
            proc.stdout.close()
    except Exception:
        pass
    try:
        if proc.stderr:
            proc.stderr.close()
    except Exception:
        pass
    
    # Collapse and Write
    print("Collapsing UMIs and writing outputs...")
    
    # Write UMI counts
    with open(f"{out_prefix}.dgecounts.csv", 'w') as f:
        f.write("Gene,Barcode,Count\n")
        for bc, genes in umi_data.items():
            for gene, umis in genes.items():
                count = greedy_umi_collapse(umis, ham_dist)
                f.write(f"{gene},{bc},{count}\n")
                
    # Write Read counts
    with open(f"{out_prefix}.readcounts.csv", 'w') as f:
        f.write("Gene,Barcode,Count\n")
        for bc, genes in read_data.items():
            for gene, count in genes.items():
                f.write(f"{gene},{bc},{count}\n")
                
    if is_smart3:
        with open(f"{out_prefix}.internal_readcounts.csv", 'w') as f:
            f.write("Gene,Barcode,Count\n")
            for bc, genes in internal_data.items():
                for gene, count in genes.items():
                    f.write(f"{gene},{bc},{count}\n")

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 dge_analysis.py <yaml_config> <samtools_exec>")
        sys.exit(1)
        
    yaml_file = sys.argv[1]
    samtools = sys.argv[2]
    
    config = load_config(yaml_file)
    project = config['project']
    out_dir = config['out_dir']
    num_threads = config['num_threads']
    ham_dist = config['counting_opts'].get('Ham_Dist', 0)
    
    # Determine input BAM (output of featureCounts)
    # The R script names it: project.filtered.Aligned.GeneTagged.bam
    # Wait, R script sorts it to: project.filtered.Aligned.GeneTagged.sorted.bam
    
    # We need to coordinate with the modified R script.
    # Let's assume R script generates: out_dir/project.filtered.Aligned.GeneTagged.bam
    
    bam_file = os.path.join(out_dir, f"{project}.filtered.Aligned.GeneTagged.bam")
    
    if not os.path.exists(bam_file):
        # Maybe it is sorted?
        bam_file = os.path.join(out_dir, f"{project}.filtered.Aligned.GeneTagged.sorted.bam")
        
    if not os.path.exists(bam_file):
        print(f"Error: Input BAM not found: {bam_file}")
        sys.exit(1)
        
    # Check for Smart-seq3
    # Original logic: check if sequence files contain pattern "ATTGCGCAATG"
    # We can pass this as arg or check config
    is_smart3 = False
    # Simplified check
    seq_files = config.get('sequence_files', {})
    for k, v in seq_files.items():
        if 'find_pattern' in v and "ATTGCGCAATG" in str(v['find_pattern']):
            is_smart3 = True
            break
            
    out_prefix = os.path.join(out_dir, "zUMIs_output", "expression", project)
    if not os.path.exists(os.path.dirname(out_prefix)):
        os.makedirs(os.path.dirname(out_prefix))
        
    process_bam(bam_file, samtools, int(ham_dist), num_threads, out_prefix, is_smart3)

if __name__ == "__main__":
    main()
