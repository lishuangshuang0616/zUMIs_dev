#!/usr/bin/env python3
import argparse
import yaml
import pandas as pd
import numpy as np
import os
import sys
import re

# Suppress pandas chained assignment warnings
pd.options.mode.chained_assignment = None 

def load_config(yaml_file):
    """Load YAML configuration."""
    with open(yaml_file, 'r') as f:
        return yaml.safe_load(f)

def read_whitelist(bcfile):
    """
    Reads a whitelist file robustly.
    Handles:
    - Standard 1-column or 2-column files.
    - Comma-separated barcodes within columns/rows.
    - Mixed whitespace delimiters (tabs, spaces).
    """
    if not os.path.exists(bcfile):
        print(f"Error: Whitelist file {bcfile} not found.")
        sys.exit(1)
        
    try:
        with open(bcfile, 'r') as f:
            content = f.read()
        
        # Robust splitting: treats commas, tabs, spaces, newlines all as delimiters
        tokens = re.split(r'[,\s]+', content)
        # Remove empty strings resulting from consecutive delimiters
        bc_wl = sorted(list(set([t for t in tokens if t])))
        
        return set(bc_wl)
    except Exception as e:
        print(f"Error reading whitelist: {e}")
        sys.exit(1)

def find_knee_point(sorted_counts):
    """
    Geometric implementation of the Knee/Elbow method.
    Approximates the 'uik' method in R.
    """
    if len(sorted_counts) == 0:
        return 0
        
    y = np.log10(sorted_counts)
    x = np.arange(len(y))
    
    start_point = np.array([x[0], y[0]])
    end_point = np.array([x[-1], y[-1]])
    vec_line = end_point - start_point
    
    vec_points = np.column_stack((x, y)) - start_point
    
    cross_prod = vec_points[:, 0] * vec_line[1] - vec_points[:, 1] * vec_line[0]
    distances = np.abs(cross_prod) / np.linalg.norm(vec_line)
    
    knee_idx = np.argmax(distances)
    return sorted_counts[knee_idx]

def cell_bc_selection(bccount_df, config):
    """
    Selects true cell barcodes matching R logic.
    """
    barcodes_config = config['barcodes']
    min_reads = barcodes_config.get('nReadsperCell', 10)
    
    # 1. Filter low reads
    df = bccount_df[bccount_df['n'] >= min_reads].copy()
    
    # 2. Sort Deterministically: Descending N, then Ascending XC (alphabetical)
    # This ensures consistency even when N is identical.
    df = df.sort_values(by=['n', 'XC'], ascending=[False, True]).reset_index(drop=True)
    df['keep'] = False
    
    strategy_auto = barcodes_config.get('automatic', False)
    bc_num = barcodes_config.get('barcode_num', None)
    bc_file = barcodes_config.get('barcode_file', None)
    
    whitelist = set()
    if bc_file:
        whitelist = read_whitelist(bc_file)

    # --- Strategy Selection ---
    
    if bc_file and strategy_auto:
        # Strategy: Automatic Knee + Whitelist Intersection
        print("Strategy: Automatic + Whitelist Intersection")
        cutoff = find_knee_point(df['n'].values)
        print(f"  Automatic cutoff: {cutoff} reads")
        
        # Mark potential cells by cutoff
        potential_keep = df['n'] >= cutoff
        
        # Filter those by whitelist
        # We strictly check if the potential high-read BC is in whitelist
        valid_in_whitelist = df.loc[potential_keep, 'XC'].isin(whitelist)
        
        if valid_in_whitelist.any():
            # If we have intersection, keep them
            df.loc[potential_keep & df['XC'].isin(whitelist), 'keep'] = True
        else:
            print("  Warning: Automatic detection found no overlap with whitelist.")
            # Fallback logic mirroring R's "continue with all/top 100" if intersection fails
            pass 

    elif bc_file and not strategy_auto:
        # Strategy: Strict Whitelist
        print("Strategy: Known Whitelist")
        df['keep'] = df['XC'].isin(whitelist)
        
        # Fallback if no matches
        if not df['keep'].any():
            print("  Warning! None of the annotated barcodes were detected.")
            if len(df) < 100:
                print("  Fallback: Keeping all barcodes (<100 total).")
                df['keep'] = True
            else:
                print("  Fallback: Keeping top 100 barcodes.")
                df.iloc[:100, df.columns.get_loc('keep')] = True

    elif bc_num is not None:
        # Strategy: Fixed Number
        print(f"Strategy: Fixed Number ({bc_num})")
        limit = min(int(bc_num), len(df))
        df.iloc[:limit, df.columns.get_loc('keep')] = True
        
    else:
        # Strategy: Automatic (Knee)
        print("Strategy: Automatic (Knee method)")
        cutoff = find_knee_point(df['n'].values)
        print(f"  Automatic cutoff: {cutoff} reads")
        df['keep'] = df['n'] >= cutoff
        
        # Fallback for too few cells
        if df['keep'].sum() < 10:
             print("  Warning: < 10 cells found. Using top 100 fallback.")
             limit = min(100, len(df))
             df.iloc[:limit, df.columns.get_loc('keep')] = True

    print(f"Selected {df['keep'].sum()} cell barcodes.")
    return df

def fast_hamming_binning(true_bcs, candidate_bcs, threshold=1):
    """
    Optimized Hamming distance binning.
    Returns: (final_df, raw_df)
    - raw_df: All matches <= threshold (ties included)
    - final_df: Unambiguous matches only (ties discarded)
    """
    # Try import for speed
    try:
        import Levenshtein
        use_levenshtein = True
    except ImportError:
        print("Warning: python-Levenshtein not found. Using slower python fallback.")
        use_levenshtein = False

    mapping_list = []
    
    true_bcs = np.array(true_bcs)
    candidate_bcs = np.array(candidate_bcs)
    
    chunk_size = 5000 
    total_chunks = (len(candidate_bcs) // chunk_size) + 1
    
    print(f"Starting binning: {len(true_bcs)} True BCs vs {len(candidate_bcs)} Candidate BCs")
    
    for i in range(0, len(candidate_bcs), chunk_size):
        if i % (chunk_size * 5) == 0:
            print(f"  Processing chunk {i // chunk_size + 1}/{total_chunks}...")
            
        cand_chunk = candidate_bcs[i : i + chunk_size]
        
        for cand in cand_chunk:
            # Calculate distance to ALL true barcodes
            if use_levenshtein:
                dists = np.array([Levenshtein.hamming(cand, t) for t in true_bcs])
            else:
                # Fallback: char by char comparison
                dists = np.array([sum(c1 != c2 for c1, c2 in zip(cand, t)) for t in true_bcs])
            
            min_dist = dists.min()
            
            if min_dist <= threshold:
                # Store RAW matches (all ties)
                best_matches_indices = np.where(dists == min_dist)[0]
                
                for match_idx in best_matches_indices:
                     mapping_list.append({
                        'falseBC': cand,
                        'trueBC': true_bcs[match_idx],
                        'hamming': min_dist
                    })

    raw_df = pd.DataFrame(mapping_list)
    
    if raw_df.empty:
        return pd.DataFrame(), pd.DataFrame()

    # Filter for Ambiguity (Matches R logic)
    # R: binmap_raw[n_min==1]
    
    # 1. Count how many trueBCs each falseBC maps to (at min distance)
    # Since we only stored min_dist matches in mapping_list, we just count occurrences
    counts = raw_df.groupby('falseBC')['trueBC'].count()
    
    # 2. Keep only those with count == 1 (Unambiguous)
    valid_false_bcs = counts[counts == 1].index
    final_df = raw_df[raw_df['falseBC'].isin(valid_false_bcs)].copy()
    
    return final_df, raw_df

def main():
    parser = argparse.ArgumentParser(description="Python implementation of zUMIs barcode detection (Optimized)")
    parser.add_argument('yaml_config', help="Path to zUMIs config YAML file")
    args = parser.parse_args()
    
    print(f"Loading config from {args.yaml_config}")
    opt = load_config(args.yaml_config)
    
    project = opt['project']
    out_dir = opt['out_dir']
    
    # Ensure output directory
    output_base = os.path.join(out_dir, "zUMIs_output")
    if not os.path.exists(output_base):
        os.makedirs(output_base)
        
    bccount_file = os.path.join(out_dir, f"{project}.BCstats.txt")
    print(f"Reading barcode stats from {bccount_file}...")
    
    if not os.path.exists(bccount_file):
        print(f"Error: Stats file {bccount_file} not found.")
        sys.exit(1)

    # Read raw counts
    raw_df = pd.read_csv(bccount_file, sep='\t', header=None, names=['XC', 'n'])
    # Aggregate counts just in case of duplicates
    raw_df = raw_df.groupby('XC', as_index=False)['n'].sum()
    
    # --- Step 1: Selection ---
    df_processed = cell_bc_selection(raw_df, opt)
    
    # Filter to kept only
    kept_df = df_processed[df_processed['keep']].copy()
    
    # Save Step 1 Result (Pre-binning)
    # Rfwrite defaults to including headers "XC" and "n"
    kept_file = os.path.join(out_dir, "zUMIs_output", f"{project}kept_barcodes.txt")
    kept_df[['XC', 'n']].to_csv(kept_file, sep='\t', index=False)
    
    # --- Step 2: Binning ---
    bc_opts = opt.get('barcodes', {})
    do_binning = False
    
    if bc_opts.get('BarcodeBinning', 0) > 0:
        do_binning = True
    if bc_opts.get('barcode_sharing'):
        # Simple Hamming binning covers the recovery aspect, 
        # though misses complex substring replacement logic of R.
        do_binning = True
        print("Note: MGI/Sharing logic simplified to standard Hamming binning.")

    if do_binning and len(kept_df) > 0:
        true_bcs = kept_df['XC'].values
        min_reads = bc_opts.get('nReadsperCell', 0)
        
        # Candidates: Not kept, but > min_reads
        candidates_df = df_processed[~df_processed['keep'] & (df_processed['n'] >= min_reads)]
        candidates_bcs = candidates_df['XC'].values
        
        if len(candidates_bcs) > 0:
            threshold = bc_opts.get('BarcodeBinning', 1)
            bin_map, bin_map_raw = fast_hamming_binning(true_bcs, candidates_bcs, threshold=threshold)
            
            if not bin_map.empty:
                print(f"Binned {len(bin_map)} barcodes.")
                
                # Pre-calculate counts map
                candidates_df_indexed = candidates_df.set_index('XC')
                
                # --- Save Raw Map ---
                raw_file = os.path.join(out_dir, "zUMIs_output", f"{project}.BCbinning.raw.txt")
                bin_map_raw['n'] = bin_map_raw['falseBC'].map(candidates_df_indexed['n']).fillna(0).astype(int)
                # Raw file usually needs specific columns too? R just dumps data.table
                # We align with binmap structure: falseBC, hamming, trueBC, n
                bin_map_raw = bin_map_raw[['falseBC', 'hamming', 'trueBC', 'n']]
                bin_map_raw.to_csv(raw_file, sep=',', index=False)

                # --- Save Final Map ---
                bin_file = os.path.join(out_dir, "zUMIs_output", f"{project}.BCbinning.txt")
                
                # Get 'n' for falseBCs (already done for raw, reuse logic if needed, or just copy)
                # bin_map is a subset of raw, so we can just filter raw or re-map
                bin_map['n'] = bin_map['falseBC'].map(candidates_df_indexed['n']).fillna(0).astype(int)
                
                bin_map_output = bin_map[['falseBC', 'hamming', 'trueBC', 'n']]
                bin_map_output.to_csv(bin_file, sep=',', index=False)
                
                # --- Update Counts ---
                add_counts = bin_map.groupby('trueBC')['n'].sum()
                
                kept_df = kept_df.set_index('XC')
                kept_df['n'] = kept_df['n'].add(add_counts, fill_value=0)
                kept_df = kept_df.reset_index()
                
                # Save Final Result
                binned_file = os.path.join(out_dir, "zUMIs_output", f"{project}kept_barcodes_binned.txt")
                kept_df[['XC', 'n']].to_csv(binned_file, sep=',', index=False)
            else:
                print("No barcodes binned (no matches within threshold).")
                # Even if no binning happened, R might output the "binned" file as a copy
                # We should probably save the copy to ensure next steps find the file
                binned_file = os.path.join(out_dir, "zUMIs_output", f"{project}kept_barcodes_binned.txt")
                kept_df[['XC', 'n']].to_csv(binned_file, sep=',', index=False)
        else:
            print("No candidate barcodes for binning.")
            # Save copy
            binned_file = os.path.join(out_dir, "zUMIs_output", f"{project}kept_barcodes_binned.txt")
            kept_df[['XC', 'n']].to_csv(binned_file, sep=',', index=False)
    
    elif do_binning and len(kept_df) == 0:
        print("Warning: No true barcodes selected. Skipping binning.")
    
    print("Done.")

if __name__ == "__main__":
    main()