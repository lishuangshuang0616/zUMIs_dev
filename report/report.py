import pandas as pd
import json
import numpy as np
import math
import collections
from pathlib import Path
from string import Template
from dnbc4tools.tools.plotly_draw import downsample_scatterplot_by_density, segment_log_plot
from dnbc4tools.rna.src.generate_report import get_args_from_file as get_rna_args
from dnbc4tools.vdj.src.generate_report import get_args_from_file as get_vdj_args
from dnbc4tools.atac.src.generate_report import get_args_from_file as get_atac_args
from dnbc4tools.__init__ import __version__
from dnbc4tools.multi.src.config import infer_enabled_sections

def _as_bool(v):
    if isinstance(v, bool):
        return v
    if isinstance(v, str):
        return v.strip().lower() in ('true', '1', 'yes', 'y', 't')
    if isinstance(v, (int, float)):
        return bool(v)
    return False



def plot_cmap(density):
    # Gray-green to deep green palette (40 stops), deepest is #0D9488
    plot_colors =  [
        "#DDDDDD","#D6D9DC","#CFD6DB","#C8D3DA","#C1D0D9",
        "#BACDD9","#B3C9D8","#ACC6D7","#A5C3D6","#9EC0D6",
        "#97BDD5","#90BAD4","#89B6D3","#82B3D3","#7BB0D2",
        "#74ADD1","#6DAAD0","#66A6CF","#5FA3CF","#58A0CE",
        "#539DCC","#4F99CA","#4C95C8","#4992C6","#458EC3",
        "#428AC1","#3F87BF","#3B83BD","#3880BA","#347CB8",
        "#3178B6","#2E75B4","#2A71B1","#276DAF","#236AAD",
        "#2066AB","#1D62A8","#195FA6","#165BA4","#1358A2"
    ]

    levels = len(plot_colors)
    ind = min(levels - 1, int(math.floor(levels * density)))
    return plot_colors[ind]

def create_report_directories(sample_outdir, config):
    """
    Create necessary REPORT directory structure for multi-omics analysis
    based on the configuration
    """
    sample_path = Path(sample_outdir)
    
    # Map config keys to directory names
    omics_mapping = {
        'rna': 'RNA',
        'vdj-t': 'VDJ-T', 
        'vdj-b': 'VDJ-B',
        'atac': 'ATAC'
    }
    
    enabled = set(infer_enabled_sections(config))

    for config_key, dir_name in omics_mapping.items():
        if config_key in enabled:
            report_dir = sample_path / f'{dir_name}_ANALYSIS_WORKFLOW_PROCESSING' / 'REPORT'
            div_dir = report_dir / 'div'
            table_dir = report_dir / 'table'
            
            # Create directories if they don't exist
            div_dir.mkdir(parents=True, exist_ok=True)
            table_dir.mkdir(parents=True, exist_ok=True)
            
            print(f"Created REPORT directories for {dir_name}")
    
    # Create final output directory 'outs'
    outs_dir = sample_path / 'outs'
    outs_dir.mkdir(parents=True, exist_ok=True)
    print(f"Created output directory: {outs_dir}")
    
    return outs_dir

def get_omics_data(outdir, sample, config):
    omics_data = {}
    outdir = Path(outdir)
    sample_outdir = outdir / sample

    enabled = set(infer_enabled_sections(config))

    # RNA
    if 'rna' in enabled:
        try:
            rna_section = config.get('rna', {})
            genome_dir_val = rna_section.get('genomeDir')
            species = Path(genome_dir_val).name if genome_dir_val else 'undefined'
            intron = not config['rna'].get('no_introns', False)
            end5 = config['rna'].get('end5', False)
            rna_proc_dir = sample_outdir / 'RNA_ANALYSIS_WORKFLOW_PROCESSING'
            if rna_proc_dir.exists():
                stat, plot_dict, table = get_rna_args(sample_outdir, intron, species, sample, end5, prefix='rna_')
                omics_data['rna'] = {'stat': stat, 'plot_dict': plot_dict, 'table': table}
            else:
                print("Warning: RNA processing directory missing; skipping RNA in report.")
        except Exception as e:
            print(f"Warning: Failed to load RNA report data: {e}")

    # VDJ-T
    if 'vdj-t' in enabled:
        try:
            species = config.get('vdj-t', {}).get('ref', 'undefined')
            vdj_t_proc_dir = sample_outdir / 'VDJ-T_ANALYSIS_WORKFLOW_PROCESSING'
            if vdj_t_proc_dir.exists():
                stat, plot_dict, table = get_vdj_args(sample_outdir, species, sample, 'TR', prefix='vdj-t_')
                omics_data['vdj-t'] = {'stat': stat, 'plot_dict': plot_dict, 'table': table}
            else:
                print("Warning: VDJ-T processing directory missing; skipping VDJ-T in report.")
        except Exception as e:
            print(f"Warning: Failed to load VDJ-T report data: {e}")

    # VDJ-B
    if 'vdj-b' in enabled:
        try:
            species = config.get('vdj-b', {}).get('ref', 'undefined')
            vdj_b_proc_dir = sample_outdir / 'VDJ-B_ANALYSIS_WORKFLOW_PROCESSING'
            if vdj_b_proc_dir.exists():
                stat, plot_dict, table = get_vdj_args(sample_outdir, species, sample, 'IG', prefix='vdj-b_')
                omics_data['vdj-b'] = {'stat': stat, 'plot_dict': plot_dict, 'table': table}
            else:
                print("Warning: VDJ-B processing directory missing; skipping VDJ-B in report.")
        except Exception as e:
            print(f"Warning: Failed to load VDJ-B report data: {e}")

    # ATAC
    if 'atac' in enabled:
        try:
            atac_section = config.get('atac', {})
            genome_dir_val = atac_section.get('genomeDir')
            species = Path(genome_dir_val).name if genome_dir_val else 'undefined'
            atac_proc_dir = sample_outdir / 'ATAC_ANALYSIS_WORKFLOW_PROCESSING'
            if atac_proc_dir.exists():
                stat, plot_dict = get_atac_args(atac_proc_dir, sample, species)
                omics_data['atac'] = {'stat': stat, 'plot_dict': plot_dict}
            else:
                print("Warning: ATAC processing directory missing; skipping ATAC in report.")
        except Exception as e:
            print(f"Warning: Failed to load ATAC report data: {e}")

    return omics_data

def load_clonotype_data_as_json(sample_outdir, immuetype):
    """Load clonotype data from CSV and convert to JSON format for template"""
    try:
        if immuetype.upper() in ('TR', 'TCR'):
            processing_dir_name = 'VDJ-T_ANALYSIS_WORKFLOW_PROCESSING'
        elif immuetype.upper() in ('IG', 'BCR'):
            processing_dir_name = 'VDJ-B_ANALYSIS_WORKFLOW_PROCESSING'
        else:
            return '[]'
            
        clonotype_file = sample_outdir / processing_dir_name / 'CLONOTYPE' / 'clonotypes.csv'
        
        if not clonotype_file.exists():
            return '[]'
            
        clonotype_df = pd.read_csv(clonotype_file, encoding='utf8')
        
        # Take top 10 clonotypes and format for JavaScript
        top_clonotypes = clonotype_df.head(10)
        clonotype_list = []
        
        for _, row in top_clonotypes.iterrows():
            clonotype_id = str(row['clonotype_id']).replace('clonotype', '')
            cdr3s_aa = str(row['cdr3s_aa']) if pd.notna(row['cdr3s_aa']) else ''
            frequency = int(row['frequency']) if pd.notna(row['frequency']) else 0
            proportion = float(row['proportion']) if pd.notna(row['proportion']) else 0.0
            
            clonotype_list.append({
                'id': clonotype_id,
                'cdr3s': cdr3s_aa,
                'frequency': frequency,
                'proportion': proportion
            })
            
        return json.dumps(clonotype_list)
    except Exception as e:
        print(f"Warning: Failed to load clonotype data: {e}")
        return '[]'

def load_marker_table_data_as_json(sample_outdir):
    """Load marker table data from CSV and convert to JSON format for template"""
    try:
        marker_file = sample_outdir / 'RNA_ANALYSIS_WORKFLOW_PROCESSING' / 'ANALYSIS' / 'marker_table.csv'
        
        if not marker_file.exists():
            return '[]'
            
        marker_df = pd.read_csv(marker_file, encoding='utf8')
        
        # Convert DataFrame to list of dictionaries for JavaScript
        marker_list = []
        
        for _, row in marker_df.iterrows():
            marker_dict = {
                'id': str(row.get('gene_ids', '')),
                'name': str(row.get('gene_name', '')),
                'gene_symbols': str(row.get('gene_symbols', ''))
            }
            
            # Add cluster-specific data (log2FC and p-values)
            for col in marker_df.columns:
                if col.startswith('cluster') and ('_lg2fc' in col or '_pval' in col):
                    marker_dict[col] = str(row[col]) if pd.notna(row[col]) else ''
            
            marker_list.append(marker_dict)
            
        return json.dumps(marker_list)
    except Exception as e:
        print(f"Warning: Failed to load marker table data: {e}")
        return '[]'

# Helper functions for processing omics data
def _process_rna_data(omics_data, sample_outdir, combined_context):
    """Process RNA data and add to combined context"""
    if 'rna' not in omics_data:
        return
        
    # Add basic RNA statistics and plots
    for key, value in omics_data['rna']['stat'].items():
        combined_context[f'rna_{key}'] = value
    combined_context.update(omics_data['rna']['plot_dict'])
    combined_context['rna_table'] = omics_data['rna']['table']
    
    # Process RNA-specific data
    _process_rna_saturation_data(sample_outdir, combined_context)
    _process_rna_display_flags(combined_context)
    _process_rna_qc_data(sample_outdir, combined_context)
    _process_rna_marker_data(sample_outdir, combined_context)
    _process_rna_cell_annotation_data(sample_outdir, combined_context)
    _process_rna_barcode_rank_plot(sample_outdir, combined_context)
    _process_rna_bead_count_data(sample_outdir, combined_context)
    _process_rna_cluster_assignment(sample_outdir, combined_context)

def _process_rna_saturation_data(sample_outdir, combined_context):
    """Prepare RNA saturation arrays for plotting"""
    try:
        sat_path = sample_outdir / 'logs/.temp' / '.rna_saturation_metrics.tsv'
        colnames = [
            'sampling_fraction',
            'Median Genes per Cell',
            'Mean Gene per Cell',
            'Total Gene',
            'Sequencing Saturation',
            'UMI Saturation'
        ]
        if sat_path.exists():
            sat_df = pd.read_table(sat_path, encoding='utf8', names=colnames, sep='\t')
            if not sat_df.empty:
                # Coerce to numeric and scale percent columns
                sat_df['sampling_fraction'] = pd.to_numeric(sat_df['sampling_fraction'], errors='coerce').fillna(0)
                sat_df['Median Genes per Cell'] = pd.to_numeric(sat_df['Median Genes per Cell'], errors='coerce').fillna(0)
                sat_df['Sequencing Saturation'] = pd.to_numeric(sat_df['Sequencing Saturation'], errors='coerce').fillna(0) * 100.0
                combined_context['rna_saturation_sampling_fraction'] = json.dumps(sat_df['sampling_fraction'].tolist())
                combined_context['rna_saturation_seq_pct'] = json.dumps(sat_df['Sequencing Saturation'].tolist())
                combined_context['rna_saturation_median_genes'] = json.dumps(sat_df['Median Genes per Cell'].tolist())
                combined_context['rna_saturation_available'] = True
            else:
                combined_context['rna_saturation_sampling_fraction'] = '[]'
                combined_context['rna_saturation_seq_pct'] = '[]'
                combined_context['rna_saturation_median_genes'] = '[]'
                combined_context['rna_saturation_available'] = False
        else:
            combined_context['rna_saturation_sampling_fraction'] = '[]'
            combined_context['rna_saturation_seq_pct'] = '[]'
            combined_context['rna_saturation_median_genes'] = '[]'
            combined_context['rna_saturation_available'] = False
    except Exception as e:
        print(f"Warning: Failed to prepare RNA saturation arrays: {e}")

def _process_rna_display_flags(combined_context):
    """Set RNA display flags for template"""
    try:
        intron_flag = _as_bool(combined_context.get('rna_intron_boolean', False))
        combined_context['rna_intron_display'] = '✓ True' if intron_flag else '✗ False'
    except Exception as e:
        print(f"Warning: failed to set rna_intron_display: {e}")
        combined_context['rna_intron_display'] = '✗ False'
        
    try:
        end5_flag = _as_bool(combined_context.get('rna_end5', False))
        combined_context['rna_end_display'] = "Single Cell 5'" if end5_flag else "Single Cell 3'"
    except Exception as e:
        print(f"Warning: failed to set rna_end_display: {e}")
        combined_context['rna_end_display'] = "Single Cell 3'"

def _process_rna_qc_data(sample_outdir, combined_context):
    """Process RNA QC data for violin plots"""
    raw_qc_path = sample_outdir / 'RNA_ANALYSIS_WORKFLOW_PROCESSING/ANALYSIS/raw_qc.xls'
    if raw_qc_path.exists():
        raw_qc_df = pd.read_table(raw_qc_path, encoding='utf8')
        # Provide arrays for three violin plots sourced from raw_qc columns
        try:
            def col_to_list(df, col, fallback=None):
                if col in df.columns:
                    return df[col].dropna().astype(float).tolist()
                if fallback and fallback in df.columns:
                    return df[fallback].dropna().astype(float).tolist()
                return []

            # Map to plots: Genes -> n_genes_by_counts; UMIs -> total_counts_mt; Mito percent -> pct_counts_mt
            genes_vals = col_to_list(raw_qc_df, 'n_genes_by_counts')
            umis_vals = col_to_list(raw_qc_df, 'total_counts')
            mito_vals = col_to_list(raw_qc_df, 'pct_counts_mt')

            combined_context['rna_violin_mid'] = json.dumps(genes_vals)
            combined_context['rna_violin_gene'] = json.dumps(umis_vals)
            combined_context['rna_violin_mito'] = json.dumps(mito_vals)
        except Exception as e:
            print(f"Warning: Failed to prepare raw_qc arrays: {e}")
            combined_context['rna_violin_mid'] = '[]'
            combined_context['rna_violin_gene'] = '[]'
            combined_context['rna_violin_mito'] = '[]'
    else:
        combined_context['rna_violin_mid'] = '[]'
        combined_context['rna_violin_gene'] = '[]'
        combined_context['rna_violin_mito'] = '[]'

def _process_rna_marker_data(sample_outdir, combined_context):
    """Add RNA marker table data as JSON"""
    combined_context['rna_marker_data'] = load_marker_table_data_as_json(sample_outdir)

def _process_rna_cell_annotation_data(sample_outdir, combined_context):
    """Check if cell annotation is available"""
    cluster_csv_path = sample_outdir / 'RNA_ANALYSIS_WORKFLOW_PROCESSING/ANALYSIS/cluster.csv'
    cell_annotation_available = False
    if cluster_csv_path.exists():
        try:
            cluster_df = pd.read_csv(cluster_csv_path, encoding='utf8')
            cell_annotation_available = 'Predicted cell type' in cluster_df.columns
        except Exception as e:
            print(f"Warning: Failed to check cell annotation availability: {e}")
    combined_context['cell_annotation_available'] = cell_annotation_available

def _process_rna_barcode_rank_plot(sample_outdir, combined_context):
    """Prepare RNA Barcode Rank Plot data arrays"""
    try:
        singlecell_csv_path = sample_outdir / 'RNA_ANALYSIS_WORKFLOW_PROCESSING/WRITE_MATRIX/singlecell.csv'
        if singlecell_csv_path.exists():
            try:
                sc_df = pd.read_csv(
                    singlecell_csv_path,
                    encoding='utf8',
                    usecols=['is_cell_barcode', 'umi'],
                    low_memory=False
                )
            except Exception:
                # Fallback to full read if selective columns fail
                sc_df = pd.read_csv(singlecell_csv_path, encoding='utf8', low_memory=False)
        else:
            sc_df = pd.DataFrame(columns=['is_cell_barcode', 'umi'])

        # Normalize and filter
        if 'umi' not in sc_df.columns:
            sc_df['umi'] = 0
        if 'is_cell_barcode' not in sc_df.columns:
            sc_df['is_cell_barcode'] = 0
        sc_df['umi'] = pd.to_numeric(sc_df['umi'], errors='coerce').fillna(0)
        sc_df['is_cell_barcode'] = pd.to_numeric(sc_df['is_cell_barcode'], errors='coerce').fillna(0).astype(int)
        sc_df = sc_df[sc_df['umi'] > 0].copy()

        # Sort by UMI desc and build rank
        sc_df = sc_df.sort_values(by='umi', ascending=False).reset_index(drop=True)
        sc_df['rank'] = sc_df.index + 1
        total_bc = len(sc_df)
        # Compute boundary indices ix1 and ix2 consistent with plotly_draw
        cell_bc = np.array(sc_df.index[sc_df['is_cell_barcode'] == 1])
        if len(cell_bc) == 0:
            dup_first = sc_df.drop_duplicates('is_cell_barcode', keep='first').index
            dup_last = sc_df.drop_duplicates('is_cell_barcode', keep='last').index
            ix1 = int(dup_first[0]) if len(dup_first) > 0 else 0
            ix2 = int(dup_last[0]) if len(dup_last) > 0 else ix1
        else:
            dup_first = sc_df.drop_duplicates('is_cell_barcode', keep='first').index
            dup_last = sc_df.drop_duplicates('is_cell_barcode', keep='last').index
            # Last cell index approximated as the position before the first 0
            if len(dup_first) >= 2:
                ix1 = int(dup_first[1] - 1)
            else:
                ix1 = int(cell_bc.max())
            # Index where noise starts (fallback to last cell index)
            ix2 = int(dup_last[0]) if len(dup_last) > 0 else ix1

        # Build segments: start TRUE, mixed in [ix1..ix2], and final NOISE
        BarcodeSegment = collections.namedtuple('BarcodeSegment', ['start', 'end', 'density', 'legend'])
        segments = []
        if total_bc > 0:
            # TRUE leading segment: density forced to 1.0 for color/hover consistency
            segments.append(BarcodeSegment(start=0, end=max(ix1, 0), density=1.0, legend=True))
            # Mixed segments based on log-length partition
            if ix2 > ix1:
                sorted_counts = sc_df['umi'].to_numpy()
                seg_bounds = segment_log_plot(sorted_counts, ix1, ix2)
                for i in range(len(seg_bounds) - 1):
                    s = int(seg_bounds[i])
                    e = int(seg_bounds[i + 1])
                    # Count cells within segment using original cell indices
                    num_cells = sum(1 for j in range(s, e) if j in cell_bc)
                    density = float(num_cells) / float(max(e - s, 1))
                    segments.append(BarcodeSegment(start=s, end=e, density=density, legend=False))
            # NOISE trailing segment
            if (ix2 + 1) < total_bc:
                segments.append(BarcodeSegment(start=ix2 + 1, end=total_bc, density=0.0, legend=True))

        # Convert segments to traces, pruning duplicate UMIs like plotly_draw
        traces = []
        for seg in segments:
            start = max(0, seg.start - 1)
            end = seg.end
            sel = sc_df.iloc[start:end]
            # Drop indices that are duplicated for UMI both at first and last occurrence
            dp_first = set(sel[sel[["umi"]].duplicated(keep="first")].index)
            dp_last = set(sel[sel[["umi"]].duplicated(keep="last")].index)
            dp_inter = dp_first & dp_last
            if dp_inter:
                sel = sel.drop(list(dp_inter), axis=0)

            x_vals = sel['rank'].tolist()
            y_vals = sel['umi'].tolist()
            trace_name = "TRUE" if seg.density > 0 else "NOISE"
            if seg.density > 0:
                n_barcodes = max(seg.end - seg.start, 1)
                n_cells = int(round(seg.density * n_barcodes))
                hover = "{:.0f}% Cell<br>({}/{})".format(100 * seg.density, n_cells, n_barcodes)
            else:
                hover = "NOISE"

            hover_prefix = [f"Rank: {rv} | UMI: {uv}" for rv, uv in zip(x_vals, y_vals)]
            hover_texts = [f"{p}<br>{hover}" for p in hover_prefix]

            traces.append({
                "x": x_vals,
                "y": y_vals,
                "name": trace_name,
                "hovertemplate": "%{text}<extra></extra>",
                "text": hover_texts,
                "type": "scattergl",
                "mode": "lines",
                "line": {"width": 3, "color": plot_cmap(seg.density)},
                "showlegend": seg.legend
            })

        combined_context['rna_rankplot_data'] = json.dumps(traces)
    except Exception as e:
        print(f"Warning: Failed to prepare RNA barcode rank plot data: {e}")
        combined_context['rna_rankplot_data'] = '[]'

def _process_rna_bead_count_data(sample_outdir, combined_context):
    """Provide Bead Count Distribution data"""
    try:
        filtered_cellid_path = sample_outdir / 'RNA_ANALYSIS_WORKFLOW_PROCESSING/WRITE_MATRIX/filtered_cellid.csv'
        x_vals = list(range(1, 10))
        y_vals = [0] * 9

        if filtered_cellid_path.exists():
            df = pd.read_csv(filtered_cellid_path, sep='\t', header=None)
            # Extract bead count after '_N' from cell id string
            # If split fails (no '_N'), coerce to NaN and drop
            counts_series = df[0].astype(str).str.split('_N', expand=True)
            if counts_series.shape[1] >= 2:
                df['Count'] = pd.to_numeric(counts_series[1], errors='coerce')
                # Build frequency table for 1..9, fill missing with 0
                freq = df['Count'].value_counts().to_dict()
                y_vals = [int(freq.get(i, 0)) for i in x_vals]
            else:
                # Fallback: no '_N' found, keep zeros
                y_vals = [0] * 9
        # Inject as JSON strings for template parsing
        combined_context['bead_count_x'] = json.dumps(x_vals)
        combined_context['bead_count_y'] = json.dumps(y_vals)
    except Exception as e:
        print(f"Warning: Failed to prepare bead count distribution: {e}")
        combined_context['bead_count_x'] = json.dumps(list(range(1, 10)))
        combined_context['bead_count_y'] = json.dumps([0] * 9)

def _process_rna_cluster_assignment(sample_outdir, combined_context):
    """Provide Cluster Assignment and RNA overlay data"""
    try:
        umap_x = []
        umap_y = []
        leiden_labels = []
        numi_values = []
        unique_clusters = []
        pred_celltypes = []
        pred_types = []
        pred_counts = {}
        vdj_t_flags = None
        vdj_b_flags = None

        cluster_csv_path = sample_outdir / 'RNA_ANALYSIS_WORKFLOW_PROCESSING/ANALYSIS/cluster.csv'
        if cluster_csv_path.exists():
            cluster_df = pd.read_csv(cluster_csv_path, encoding='utf8')
            required_cols = {'UMAP_1', 'UMAP_2', 'leiden', 'nUMI'}
            if required_cols.issubset(cluster_df.columns):
                # Drop rows with missing required values
                cluster_df = cluster_df.dropna(subset=list(required_cols))
                umap_x = cluster_df['UMAP_1'].astype(float).tolist()
                umap_y = cluster_df['UMAP_2'].astype(float).tolist()
                # Convert leiden to string labels (robust to numeric)
                leiden_labels = cluster_df['leiden'].astype(str).tolist()
                # UMI values as floats
                numi_values = cluster_df['nUMI'].astype(float).tolist()
                # Unique cluster labels for legend ordering
                try:
                    # Sort numerically when possible
                    def try_float(s):
                        try:
                            return float(s)
                        except Exception:
                            return s
                    unique_clusters = sorted(cluster_df['leiden'].astype(str).unique(), key=lambda v: try_float(v))
                except Exception:
                    unique_clusters = cluster_df['leiden'].astype(str).unique().tolist()

                # Predicted cell type (optional)
                if 'Predicted cell type' in cluster_df.columns:
                    ct_series = cluster_df['Predicted cell type'].fillna('Unknown').astype(str)
                    # Clean trailing count patterns like ": 498" or "：498" from labels
                    ct_series_clean = ct_series.str.replace(r"\s*[:：]\s*\d+\s*$", '', regex=True)
                    pred_celltypes = ct_series_clean.tolist()
                    # Counts per cleaned type for legend labels
                    pred_counts = ct_series_clean.value_counts().to_dict()
                    # Order types by descending count
                    pred_types = [t for t, _ in sorted(pred_counts.items(), key=lambda kv: (-kv[1], kv[0]))]

                # VDJ flags derived from filtered_contig_annotations.csv (preferred)
                try:
                    def read_vdj_filtered_barcodes(base_outdir: Path):
                        """Read VDJ cell IDs from CELL_CONVERT only and return TR/IG sets.
                        Chain information is not required; CLONOTYPE is ignored.
                        """
                        tr_set, ig_set = set(), set()
                        tr_path = base_outdir / 'VDJ-T_ANALYSIS_WORKFLOW_PROCESSING' / 'CELL_CONVERT' / 'filtered_contig_annotations.csv'
                        ig_path = base_outdir / 'VDJ-B_ANALYSIS_WORKFLOW_PROCESSING' / 'CELL_CONVERT' / 'filtered_contig_annotations.csv'

                        def read_one(path: Path):
                            if not path.exists():
                                return pd.DataFrame()
                            try:
                                return pd.read_csv(path, encoding='utf8')
                            except Exception:
                                return pd.DataFrame()

                        tr_df = read_one(tr_path)
                        ig_df = read_one(ig_path)

                        # Normalize expected columns (barcode and is_cell)
                        def normalize(df: pd.DataFrame) -> pd.DataFrame:
                            if df.empty:
                                return df
                            cols_lower = {c.lower(): c for c in df.columns}
                            # Map possible variants
                            bc_col = cols_lower.get('barcode') or cols_lower.get('#barcode') or cols_lower.get('cell')
                            is_cell_col = cols_lower.get('is_cell') or cols_lower.get('iscell')
                            result = pd.DataFrame()
                            if bc_col:
                                # barcode in filtered file is already the cell id
                                result['barcode'] = df[bc_col].astype(str)
                            # Default to True when missing, otherwise respect provided flag
                            result['is_cell'] = df[is_cell_col].astype(bool) if is_cell_col else True
                            return result

                        tr_df = normalize(tr_df)
                        ig_df = normalize(ig_df)

                        # Keep only cell barcodes and deduplicate by barcode (no chain filtering)
                        if not tr_df.empty and {'barcode','is_cell'}.issubset(tr_df.columns):
                            valid_tr = tr_df[tr_df['is_cell']]
                            tr_set = set(valid_tr['barcode'].drop_duplicates())
                        if not ig_df.empty and {'barcode','is_cell'}.issubset(ig_df.columns):
                            valid_ig = ig_df[ig_df['is_cell']]
                            ig_set = set(valid_ig['barcode'].drop_duplicates())
                        return tr_set, ig_set

                    tr_barcodes, ig_barcodes = read_vdj_filtered_barcodes(sample_outdir)
                    # Intersect with RNA filtered_cellid.csv cells (exact match, no suffix stripping)
                    rna_cells_list = []
                    rna_cells_set = set()
                    filtered_cellid_path = sample_outdir / 'RNA_ANALYSIS_WORKFLOW_PROCESSING/WRITE_MATRIX/filtered_cellid.csv'
                    if filtered_cellid_path.exists():
                        try:
                            cellid_df = pd.read_csv(filtered_cellid_path, sep='\t', header=None)
                            rna_cells_list = cellid_df[0].astype(str).tolist()
                            rna_cells_set = set(rna_cells_list)
                        except Exception:
                            rna_cells_list = []
                            rna_cells_set = set()

                    if rna_cells_set:
                        tr_barcodes = tr_barcodes & rna_cells_set
                        ig_barcodes = ig_barcodes & rna_cells_set

                    # Determine cluster barcodes to align flags in UMAP order
                    barcode_list = None
                    # Prefer 'Unnamed: 0' (case-insensitive); avoid looping over candidates
                    cols_lower_map = {col.lower(): col for col in cluster_df.columns}
                    if 'unnamed: 0' in cols_lower_map:
                        col_name = cols_lower_map['unnamed: 0']
                        barcode_list = cluster_df[col_name].astype(str).str.strip().tolist()

                    # Fallbacks: use index if not default RangeIndex; otherwise first column
                    if barcode_list is None:
                        try:
                            if not isinstance(cluster_df.index, pd.RangeIndex):
                                barcode_list = cluster_df.index.astype(str).str.strip().tolist()
                            else:
                                barcode_list = cluster_df.iloc[:, 0].astype(str).str.strip().tolist()
                        except Exception:
                            barcode_list = None

                    # Build flags aligned to UMAP points; pad zeros if length mismatched
                    if barcode_list is not None:
                        limit = min(len(barcode_list), len(umap_x))
                        if len(tr_barcodes) > 0:
                            vdj_t_flags = [1 if barcode_list[i] in tr_barcodes else 0 for i in range(limit)]
                            if limit < len(umap_x):
                                vdj_t_flags += [0] * (len(umap_x) - limit)
                        if len(ig_barcodes) > 0:
                            vdj_b_flags = [1 if barcode_list[i] in ig_barcodes else 0 for i in range(limit)]
                            if limit < len(umap_x):
                                vdj_b_flags += [0] * (len(umap_x) - limit)
                except Exception as e:
                    print(f"Warning: Failed to derive VDJ flags from filtered_contig_annotations.csv: {e}")
            else:
                print("Warning: cluster.csv is missing required columns: UMAP_1, UMAP_2, leiden, nUMI")
        # Inject JSON strings for template parsing
        combined_context['rna_umap_x'] = json.dumps(umap_x)
        combined_context['rna_umap_y'] = json.dumps(umap_y)
        combined_context['rna_leiden'] = json.dumps(leiden_labels)
        combined_context['rna_numi'] = json.dumps(numi_values)
        combined_context['rna_clusters'] = json.dumps(unique_clusters)
        combined_context['rna_pred_celltype'] = json.dumps(pred_celltypes)
        combined_context['rna_pred_types'] = json.dumps(pred_types)
        combined_context['rna_pred_counts'] = json.dumps(pred_counts)
        if vdj_t_flags is not None:
            combined_context['vdj_t_flags'] = json.dumps(vdj_t_flags)
        if vdj_b_flags is not None:
            combined_context['vdj_b_flags'] = json.dumps(vdj_b_flags)
    except Exception as e:
        print(f"Warning: Failed to prepare cluster assignment and RNA overlay: {e}")
        combined_context['rna_umap_x'] = '[]'
        combined_context['rna_umap_y'] = '[]'
        combined_context['rna_leiden'] = '[]'
        combined_context['rna_numi'] = '[]'
        combined_context['rna_clusters'] = '[]'
        combined_context['rna_pred_celltype'] = '[]'
        combined_context['rna_pred_types'] = '[]'
        combined_context['rna_pred_counts'] = '{}'
        combined_context['vdj_t_flags'] = '[]'
        combined_context['vdj_b_flags'] = '[]'

def _process_vdj_t_data(omics_data, sample_outdir, combined_context):
    """Process VDJ-T data and add to combined context"""
    if 'vdj-t' not in omics_data:
        return
        
    # Add basic VDJ-T statistics and plots
    for key, value in omics_data['vdj-t']['stat'].items():
        combined_context[f'vdj_t_{key}'] = value
    combined_context.update(omics_data['vdj-t']['plot_dict'])
    combined_context['vdj_t_table'] = omics_data['vdj-t']['table']
    
    # Add clonotype data as JSON
    combined_context['vdj_t_clonotype_data'] = load_clonotype_data_as_json(sample_outdir, 'TR')
    
    # Prepare VDJ-T Barcode Rank Plot data arrays
    _process_vdj_barcode_rank_plot(
        sample_outdir, 
        'VDJ-T_ANALYSIS_WORKFLOW_PROCESSING', 
        'vdj_t_rankplot_data', 
        combined_context
    )

def _process_vdj_b_data(omics_data, sample_outdir, combined_context):
    """Process VDJ-B data and add to combined context"""
    if 'vdj-b' not in omics_data:
        return
        
    # Add basic VDJ-B statistics and plots
    for key, value in omics_data['vdj-b']['stat'].items():
        combined_context[f'vdj_b_{key}'] = value
    combined_context.update(omics_data['vdj-b']['plot_dict'])
    combined_context['vdj_b_table'] = omics_data['vdj-b']['table']
    
    # Add clonotype data as JSON
    combined_context['vdj_b_clonotype_data'] = load_clonotype_data_as_json(sample_outdir, 'IG')
    
    # Prepare VDJ-B Barcode Rank Plot data arrays
    _process_vdj_barcode_rank_plot(
        sample_outdir, 
        'VDJ-B_ANALYSIS_WORKFLOW_PROCESSING', 
        'vdj_b_rankplot_data', 
        combined_context
    )

def _process_vdj_barcode_rank_plot(sample_outdir, processing_dir, context_key, combined_context):
    """Generic function to process VDJ barcode rank plot data"""
    try:
        barcode_info_path = sample_outdir / processing_dir / 'CLONOTYPE' / 'barcode_info.csv'
        if barcode_info_path.exists():
            try:
                vdj_df = pd.read_csv(
                    barcode_info_path,
                    encoding='utf8',
                    usecols=['is_cell', 'num_umis'],
                    low_memory=False
                ).rename(columns={'is_cell': 'is_cell_barcode', 'num_umis': 'umi'})
            except Exception:
                vdj_df = pd.read_csv(barcode_info_path, encoding='utf8', low_memory=False)
                if 'is_cell_barcode' not in vdj_df.columns and 'is_cell' in vdj_df.columns:
                    vdj_df['is_cell_barcode'] = vdj_df['is_cell']
                if 'umi' not in vdj_df.columns:
                    for cand in ['num_umis', 'umis', 'UMI']:
                        if cand in vdj_df.columns:
                            vdj_df['umi'] = vdj_df[cand]
                            break
        else:
            vdj_df = pd.DataFrame(columns=['is_cell_barcode', 'umi'])

        if 'umi' not in vdj_df.columns:
            vdj_df['umi'] = 0
        if 'is_cell_barcode' not in vdj_df.columns:
            vdj_df['is_cell_barcode'] = 0

        vdj_df['umi'] = pd.to_numeric(vdj_df['umi'], errors='coerce').fillna(0)
        vdj_df['is_cell_barcode'] = pd.to_numeric(vdj_df['is_cell_barcode'], errors='coerce').fillna(0).astype(int)
        vdj_df = vdj_df[vdj_df['umi'] > 0].copy()

        vdj_df = vdj_df.sort_values(by='umi', ascending=False).reset_index(drop=True)
        vdj_df['rank'] = vdj_df.index + 1
        total_bc = len(vdj_df)

        cell_bc = np.array(vdj_df.index[vdj_df['is_cell_barcode'] == 1])
        if len(cell_bc) == 0:
            dup_first = vdj_df.drop_duplicates('is_cell_barcode', keep='first').index
            dup_last = vdj_df.drop_duplicates('is_cell_barcode', keep='last').index
            ix1 = int(dup_first[0]) if len(dup_first) > 0 else 0
            ix2 = int(dup_last[0]) if len(dup_last) > 0 else ix1
        else:
            dup_first = vdj_df.drop_duplicates('is_cell_barcode', keep='first').index
            dup_last = vdj_df.drop_duplicates('is_cell_barcode', keep='last').index
            if len(dup_first) >= 2:
                ix1 = int(dup_first[1] - 1)
            else:
                ix1 = int(cell_bc.max())
            ix2 = int(dup_last[0]) if len(dup_last) > 0 else ix1

        # Follow plotly_draw VDJ behavior: force ix1 to 1 regardless of computed value
        ix1 = 1

        BarcodeSegment = collections.namedtuple('BarcodeSegment', ['start', 'end', 'density', 'legend'])
        segments = []
        if total_bc > 0:
            segments.append(BarcodeSegment(start=0, end=max(ix1 + 1, 1), density=1.0, legend=True))
            if ix2 > ix1:
                sorted_counts = vdj_df['umi'].to_numpy()
                seg_bounds = segment_log_plot(sorted_counts, ix1, ix2)
                for i in range(len(seg_bounds) - 1):
                    s = int(seg_bounds[i])
                    e = int(seg_bounds[i + 1])
                    num_cells = sum(1 for j in range(s, e) if j in cell_bc)
                    density = float(num_cells) / float(max(e - s, 1))
                    segments.append(BarcodeSegment(start=s, end=e, density=density, legend=False))
            if (ix2) < total_bc:
                segments.append(BarcodeSegment(start=ix2, end=total_bc, density=0.0, legend=True))

        traces = []
        for seg in segments:
            start = max(0, seg.start - 1)
            end = seg.end
            sel = vdj_df.iloc[start:end]
            dp_first = set(sel[sel[["umi"]].duplicated(keep="first")].index)
            dp_last = set(sel[sel[["umi"]].duplicated(keep="last")].index)
            dp_inter = dp_first & dp_last
            if dp_inter:
                sel = sel.drop(list(dp_inter), axis=0)

            x_vals = sel['rank'].tolist()
            y_vals = sel['umi'].tolist()
            trace_name = "TRUE" if seg.density > 0 else "NOISE"
            if seg.density > 0:
                n_barcodes = max(seg.end - seg.start, 1)
                n_cells = int(round(seg.density * n_barcodes))
                hover = "{:.0f}% Cell<br>({}/{})".format(100 * seg.density, n_cells, n_barcodes)
            else:
                hover = "NOISE"

            hover_prefix = [f"Rank: {rv} | UMI: {uv}" for rv, uv in zip(x_vals, y_vals)]
            hover_texts = [f"{p}<br>{hover}" for p in hover_prefix]

            traces.append({
                "x": x_vals,
                "y": y_vals,
                "name": trace_name,
                "hovertemplate": "%{text}<extra></extra>",
                "text": hover_texts,
                "type": "scattergl",
                "mode": "lines",
                "line": {"width": 3, "color": plot_cmap(seg.density)},
                "showlegend": seg.legend
            })

        combined_context[context_key] = json.dumps(traces)
    except Exception as e:
        print(f"Warning: Failed to prepare {processing_dir} barcode rank plot data: {e}")
        combined_context[context_key] = '[]'

def _process_atac_data(omics_data, sample_outdir, combined_context):
    """Process ATAC data and add to combined context"""
    if 'atac' not in omics_data:
        return
        
    # Add basic ATAC statistics and plots
    for key, value in omics_data['atac']['stat'].items():
        combined_context[f'atac_{key}'] = value
    combined_context.update(omics_data['atac']['plot_dict'])
    
    # Process ATAC-specific data
    atac_processing_dir = sample_outdir / 'ATAC_ANALYSIS_WORKFLOW_PROCESSING'
    _process_atac_violin_plots(atac_processing_dir, combined_context)
    _process_atac_fragment_rank_plot(atac_processing_dir, combined_context)
    _process_atac_tss_enrichment(atac_processing_dir, combined_context)
    _process_atac_tss_targeting(atac_processing_dir, combined_context)
    _process_atac_bead_count(atac_processing_dir, combined_context)
    _process_atac_cluster_assignment(atac_processing_dir, combined_context)
    _process_atac_insert_size(atac_processing_dir, combined_context)
    _process_atac_saturation(atac_processing_dir, combined_context, combined_context.get('atac_mean_rawreads_percell', '0'))
    _process_atac_jaccard_index(atac_processing_dir, combined_context, combined_context.get('atac_jaccard_thres', '0.02'))

def _process_atac_violin_plots(atac_processing_dir, combined_context):
    """Provide arrays for ATAC mini violin plots in template"""
    try:
        singlecell_csv_path = atac_processing_dir / 'SUMMARY_METRIC' / 'singlecell.csv'
        if singlecell_csv_path.exists():
            sc_df = pd.read_csv(singlecell_csv_path, encoding='utf8')

            # Filter to valid cells if available
            if 'is_cell_barcode' in sc_df.columns:
                try:
                    sc_df = sc_df[sc_df['is_cell_barcode'] == 1]
                except Exception:
                    pass

            # Determine column names robustly
            frag_col = 'fragments' if 'fragments' in sc_df.columns else None
            tss_col = 'TSS_region_fragments' if 'TSS_region_fragments' in sc_df.columns else None
            peak_col = 'peak_region_fragments' if 'peak_region_fragments' in sc_df.columns else None
            mt_col = 'mt_region_fragments' if 'mt_region_fragments' in sc_df.columns else None
            cp_col = 'chloroplast_region_fragments' if 'chloroplast_region_fragments' in sc_df.columns else None

            frag_vals: list = []
            tss_vals: list = []
            peak_vals: list = []

            if frag_col:
                frags_series = pd.to_numeric(sc_df[frag_col], errors='coerce').fillna(0)
                frags_series = frags_series.clip(lower=0)
                # Use log10(frags + 1) for fragments violin
                frag_vals = np.log10(frags_series + 1).tolist()

            if tss_col and frag_col:
                tss_series = pd.to_numeric(sc_df[tss_col], errors='coerce').fillna(0)
                frags_series = pd.to_numeric(sc_df[frag_col], errors='coerce').fillna(0)
                mt_series = pd.to_numeric(sc_df[mt_col], errors='coerce').fillna(0) if mt_col else 0
                cp_series = pd.to_numeric(sc_df[cp_col], errors='coerce').fillna(0) if cp_col else 0
                eff_frags = (frags_series - mt_series - cp_series)
                eff_frags = eff_frags.replace([np.inf, -np.inf], np.nan).fillna(0).clip(lower=0)
                ratio = (tss_series / eff_frags).replace([np.inf, -np.inf], np.nan).dropna()
                tss_vals = ratio.tolist()

            if peak_col and frag_col:
                peak_series = pd.to_numeric(sc_df[peak_col], errors='coerce').fillna(0)
                frags_series = pd.to_numeric(sc_df[frag_col], errors='coerce').fillna(0)
                mt_series = pd.to_numeric(sc_df[mt_col], errors='coerce').fillna(0) if mt_col else 0
                cp_series = pd.to_numeric(sc_df[cp_col], errors='coerce').fillna(0) if cp_col else 0
                eff_frags = (frags_series - mt_series - cp_series)
                eff_frags = eff_frags.replace([np.inf, -np.inf], np.nan).fillna(0).clip(lower=0)
                ratio = (peak_series / eff_frags).replace([np.inf, -np.inf], np.nan).dropna()
                peak_vals = ratio.tolist()

            combined_context['atac_violin_frag'] = json.dumps(frag_vals)
            combined_context['atac_violin_tss'] = json.dumps(tss_vals)
            combined_context['atac_violin_peak'] = json.dumps(peak_vals)
        else:
            combined_context['atac_violin_frag'] = '[]'
            combined_context['atac_violin_tss'] = '[]'
            combined_context['atac_violin_peak'] = '[]'
    except Exception as e:
        print(f"Warning: Failed to prepare ATAC violin arrays: {e}")
        combined_context['atac_violin_frag'] = '[]'
        combined_context['atac_violin_tss'] = '[]'
        combined_context['atac_violin_peak'] = '[]'

def _process_atac_fragment_rank_plot(atac_processing_dir, combined_context):
    """Prepare ATAC Fragment Rank Plot data arrays"""
    try:
        sc_path = atac_processing_dir / 'SUMMARY_METRIC' / 'singlecell.csv'
        if sc_path.exists():
            atac_sc_df = pd.read_csv(sc_path, encoding='utf8', low_memory=False)
        else:
            atac_sc_df = pd.DataFrame(columns=['is_cell_barcode', 'peak_region_fragments'])

        peak_col = 'peak_region_fragments' if 'peak_region_fragments' in atac_sc_df.columns else None
        if peak_col is None:
            combined_context['atac_rankplot_data'] = '[]'
        else:
            if 'is_cell_barcode' not in atac_sc_df.columns:
                atac_sc_df['is_cell_barcode'] = 0
            atac_sc_df[peak_col] = pd.to_numeric(atac_sc_df[peak_col], errors='coerce').fillna(0)
            atac_sc_df['is_cell_barcode'] = pd.to_numeric(atac_sc_df['is_cell_barcode'], errors='coerce').fillna(0).astype(int)
            # Keep only positive fragment counts
            atac_sc_df = atac_sc_df[atac_sc_df[peak_col] > 0].copy()

            # Sort by fragments in peaks (desc) and build rank starting at 1 (for display)
            atac_sc_df = atac_sc_df.sort_values(by=peak_col, ascending=False).reset_index(drop=True)
            atac_sc_df['rank'] = atac_sc_df.index + 1

            total_bc = len(atac_sc_df)
            traces_atac = []
            if total_bc > 0:
                cell_bc_idx = np.array(atac_sc_df.index[atac_sc_df['is_cell_barcode'] == 1])
                has_cells = len(cell_bc_idx) > 0
                has_noise = np.any(atac_sc_df['is_cell_barcode'] == 0)

                dup_first = atac_sc_df.drop_duplicates('is_cell_barcode', keep='first').index
                dup_last = atac_sc_df.drop_duplicates('is_cell_barcode', keep='last').index

                # Determine boundaries like tools.plotly_draw._plot_barcoderanks_atac
                if has_cells and has_noise and len(dup_first) >= 2 and len(dup_last) >= 1:
                    ix1 = int(dup_first[1] - 1)
                    ix2 = int(dup_last[0])
                elif has_cells and not has_noise:
                    ix1 = total_bc - 1
                    ix2 = total_bc - 1
                else:
                    ix1 = 0
                    ix2 = -1

                BarcodeSegment = collections.namedtuple('BarcodeSegment', ['start', 'end', 'density', 'legend'])
                segments_atac = []

                # TRUE leading segment
                if has_cells:
                    segments_atac.append(BarcodeSegment(start=0, end=max(ix1, 0), density=1.0, legend=True))

                # Mixed segments
                if has_cells and has_noise and ix2 > ix1:
                    sorted_counts = atac_sc_df[peak_col].to_numpy()
                    seg_bounds = segment_log_plot(sorted_counts, ix1, ix2)
                    for i in range(len(seg_bounds) - 1):
                        s = int(seg_bounds[i])
                        e = int(seg_bounds[i + 1])
                        n_cells = sum(1 for j in range(s, e) if j in cell_bc_idx)
                        density = float(n_cells) / float(max(e - s, 1))
                        segments_atac.append(BarcodeSegment(start=s, end=e, density=density, legend=False))

                # NOISE trailing segment
                if (ix2 + 1) < total_bc:
                    segments_atac.append(BarcodeSegment(start=ix2 + 1, end=total_bc, density=0.0, legend=True))

                # Convert segments into traces (pruning duplicates like tools.plotly_draw)
                for seg in segments_atac:
                    start = max(0, seg.start - 1)
                    end = seg.end
                    sel = atac_sc_df.iloc[start:end]

                    dp_first = set(sel[sel[[peak_col]].duplicated(keep="first")].index)
                    dp_last = set(sel[sel[[peak_col]].duplicated(keep="last")].index)
                    dp_inter = dp_first & dp_last
                    if dp_inter:
                        sel = sel.drop(list(dp_inter), axis=0)

                    x_vals = sel['rank'].tolist()
                    y_vals = sel[peak_col].tolist()
                    trace_name = "TRUE" if seg.density > 0 else "NOISE"
                    if seg.density > 0:
                        n_barcodes = max(seg.end - seg.start, 1)
                        n_cells = int(round(seg.density * n_barcodes))
                        hover = "{:.0f}% Cell<br>({}/{})".format(100 * seg.density, n_cells, n_barcodes)
                    else:
                        hover = "NOISE"

                    hover_prefix = [f"Rank: {rv} | Fragments: {yv}" for rv, yv in zip(x_vals, y_vals)]
                    hover_texts = [f"{p}<br>{hover}" for p in hover_prefix]

                    traces_atac.append({
                        "x": x_vals,
                        "y": y_vals,
                        "name": trace_name,
                        "hoverinfo": "text",
                        "text": hover_texts,
                        "type": "scattergl",
                        "mode": "lines",
                        "line": {"width": 3, "color": plot_cmap(seg.density)},
                        "showlegend": seg.legend
                    })

            combined_context['atac_rankplot_data'] = json.dumps(traces_atac)
    except Exception as e:
        print(f"Warning: Failed to prepare ATAC fragment rank plot data: {e}")
        combined_context['atac_rankplot_data'] = '[]'

def _process_atac_tss_enrichment(atac_processing_dir, combined_context):
    """Provide ATAC TSS enrichment profile arrays"""
    try:
        tss_flankwindow_csv_path = atac_processing_dir / 'SUMMARY_METRIC' / 'tss.flankwindow.csv'
        default_x = list(range(-1000, 1001))
        if tss_flankwindow_csv_path.exists():
            tss_df = pd.read_csv(tss_flankwindow_csv_path, encoding='utf8', header=None)
            # Use first column as y-values, coerce to float and drop NaNs
            y_series = pd.to_numeric(tss_df.iloc[:, 0], errors='coerce').dropna()
            y_vals = y_series.tolist()
            # Align x length to y length if it differs from 2001
            if len(y_vals) != len(default_x):
                x_vals = list(range(-1000, -1000 + len(y_vals)))
            else:
                x_vals = default_x
            combined_context['atac_tss_x'] = json.dumps(x_vals)
            combined_context['atac_tss_y'] = json.dumps(y_vals)
        else:
            combined_context['atac_tss_x'] = json.dumps(default_x)
            combined_context['atac_tss_y'] = '[]'
    except Exception as e:
        print(f"Warning: Failed to prepare ATAC TSS enrichment arrays: {e}")
        combined_context['atac_tss_x'] = json.dumps(list(range(-1000, 1001)))
        combined_context['atac_tss_y'] = '[]'

def _process_atac_tss_targeting(atac_processing_dir, combined_context):
    """Provide ATAC TSS targeting arrays for scatter plot"""
    try:
        singlecell_csv_path_target = atac_processing_dir / 'SUMMARY_METRIC' / 'singlecell.csv'
        if singlecell_csv_path_target.exists():
            target_df = pd.read_csv(singlecell_csv_path_target, encoding='utf8')
            # Robustly extract columns and coerce types
            if 'fragments' in target_df.columns and 'TSS_region_fragments' in target_df.columns:
                totals_series = pd.to_numeric(target_df['fragments'], errors='coerce').fillna(0).clip(lower=0)
                tss_series = pd.to_numeric(target_df['TSS_region_fragments'], errors='coerce').fillna(0).clip(lower=0)
                if 'is_cell_barcode' in target_df.columns:
                    is_cell_series = pd.to_numeric(target_df['is_cell_barcode'], errors='coerce').fillna(0).astype(int)
                else:
                    is_cell_series = pd.Series([0] * len(target_df))
                if 'Cell' in target_df.columns:
                    cell_label_series = target_df['Cell'].astype(str).fillna('NO_BARCODE')
                else:
                    cell_label_series = pd.Series(['NO_BARCODE'] * len(target_df))

                # Ensure equal lengths and build base dataframe
                n = int(min(len(totals_series), len(tss_series), len(is_cell_series), len(cell_label_series)))
                base_df = pd.DataFrame({
                    'total': totals_series.iloc[:n],
                    'subtype': tss_series.iloc[:n],
                    'is_cell': is_cell_series.iloc[:n],
                    'cell_label': cell_label_series.iloc[:n],
                })
                # Filter invalid totals (avoid division by zero on log axis)
                base_df = base_df[base_df['total'] > 0]

                # Masks per plotly_draw: Non-cells = not cell & Cell != 'NO_BARCODE'; Cells = is_cell == 1
                cell_mask = base_df['is_cell'] == 1
                non_cell_mask = (base_df['is_cell'] == 0) & (base_df['cell_label'] != 'NO_BARCODE')

                # Build per-group dataframes and drop duplicate points
                cells_df = base_df.loc[cell_mask, ['total', 'subtype']].copy()
                non_cells_df = base_df.loc[non_cell_mask, ['total', 'subtype']].copy()
                cells_df.drop_duplicates(inplace=True)
                non_cells_df.drop_duplicates(inplace=True)

                # Downsample by density per group if combined points exceed threshold
                if (len(cells_df) + len(non_cells_df)) > 2000:
                    if len(cells_df) > 0:
                        cells_df = downsample_scatterplot_by_density(cells_df, 2000, 'total', 'subtype')
                    if len(non_cells_df) > 0:
                        non_cells_df = downsample_scatterplot_by_density(non_cells_df, 2000, 'total', 'subtype')

                # Compute x and y for each group (y = subtype / total)
                cells_x = cells_df['total'].astype(float).tolist()
                cells_y = (cells_df['subtype'] / cells_df['total']).astype(float).tolist()
                non_cells_x = non_cells_df['total'].astype(float).tolist()
                non_cells_y = (non_cells_df['subtype'] / non_cells_df['total']).astype(float).tolist()

                combined_context['atac_target_non_cells_x'] = json.dumps(non_cells_x)
                combined_context['atac_target_non_cells_y'] = json.dumps(non_cells_y)
                combined_context['atac_target_cells_x'] = json.dumps(cells_x)
                combined_context['atac_target_cells_y'] = json.dumps(cells_y)
            else:
                combined_context['atac_target_non_cells_x'] = '[]'
                combined_context['atac_target_non_cells_y'] = '[]'
                combined_context['atac_target_cells_x'] = '[]'
                combined_context['atac_target_cells_y'] = '[]'
        else:
            combined_context['atac_target_non_cells_x'] = '[]'
            combined_context['atac_target_non_cells_y'] = '[]'
            combined_context['atac_target_cells_x'] = '[]'
            combined_context['atac_target_cells_y'] = '[]'
    except Exception as e:
        print(f"Warning: Failed to prepare ATAC targeting arrays: {e}")
        combined_context['atac_target_non_cells_x'] = '[]'
        combined_context['atac_target_non_cells_y'] = '[]'
        combined_context['atac_target_cells_x'] = '[]'
        combined_context['atac_target_cells_y'] = '[]'

def _process_atac_bead_count(atac_processing_dir, combined_context):
    """Prepare ATAC beads per droplet distribution"""
    try:
        beads_barcodes_path = atac_processing_dir / 'PEAK_MATRIX' / 'beads_barcodes.txt'
        x_vals = list(range(1, 10))
        y_vals = [0] * 9
        if beads_barcodes_path.exists():
            df = pd.read_csv(beads_barcodes_path, sep='\t', header=None)
            counts_series = df[0].astype(str).str.split('_N', expand=True)
            if counts_series.shape[1] >= 2:
                df['Count'] = pd.to_numeric(counts_series[1], errors='coerce')
                freq = df['Count'].value_counts().to_dict()
                y_vals = [int(freq.get(i, 0)) for i in x_vals]
        combined_context['atac_bead_count_x'] = json.dumps(x_vals)
        combined_context['atac_bead_count_y'] = json.dumps(y_vals)
    except Exception as e:
        print(f"Warning: Failed to prepare ATAC bead count distribution: {e}")
        combined_context['atac_bead_count_x'] = json.dumps(list(range(1, 10)))
        combined_context['atac_bead_count_y'] = json.dumps([0] * 9)

def _process_atac_cluster_assignment(atac_processing_dir, combined_context):
    """Provide ATAC Cluster Assignment and overlay data"""
    try:
        atac_umap_x = []
        atac_umap_y = []
        atac_leiden_labels = []
        atac_unique_clusters = []
        atac_frags_vals = []

        atac_cluster_csv_path = atac_processing_dir / 'ANALYSIS' / 'cluster.csv'
        if atac_cluster_csv_path.exists():
            atac_cluster_df = pd.read_csv(atac_cluster_csv_path, encoding='utf8')
            required_cols = {'UMAP_1', 'UMAP_2', 'leiden', 'log10_uniqueFrags'}
            if required_cols.issubset(atac_cluster_df.columns):
                atac_cluster_df = atac_cluster_df.dropna(subset=list(required_cols))
                atac_umap_x = atac_cluster_df['UMAP_1'].astype(float).tolist()
                atac_umap_y = atac_cluster_df['UMAP_2'].astype(float).tolist()
                atac_leiden_labels = atac_cluster_df['leiden'].astype(str).tolist()
                atac_frags_vals = atac_cluster_df['log10_uniqueFrags'].astype(float).tolist()
                # Unique cluster labels for legend ordering
                try:
                    def try_float(v):
                        try:
                            return float(v)
                        except Exception:
                            return v
                    atac_unique_clusters = sorted(atac_cluster_df['leiden'].astype(str).unique(), key=lambda v: try_float(v))
                except Exception:
                    atac_unique_clusters = atac_cluster_df['leiden'].astype(str).unique().tolist()
            else:
                print("Warning: ATAC cluster.csv is missing required columns: UMAP_1, UMAP_2, leiden, log10_uniqueFrags")

        combined_context['atac_umap_x'] = json.dumps(atac_umap_x)
        combined_context['atac_umap_y'] = json.dumps(atac_umap_y)
        combined_context['atac_leiden'] = json.dumps(atac_leiden_labels)
        combined_context['atac_clusters'] = json.dumps(atac_unique_clusters)
        combined_context['atac_frags'] = json.dumps(atac_frags_vals)
    except Exception as e:
        print(f"Warning: Failed to prepare ATAC cluster assignment and fragment overlay: {e}")
        combined_context['atac_umap_x'] = '[]'
        combined_context['atac_umap_y'] = '[]'
        combined_context['atac_leiden'] = '[]'
        combined_context['atac_clusters'] = '[]'
        combined_context['atac_frags'] = '[]'

def _process_atac_insert_size(atac_processing_dir, combined_context):
    """Provide ATAC insert size distribution arrays"""
    try:
        fraglength_csv_path = atac_processing_dir / 'SUMMARY_METRIC' / 'fraglength.ratio.csv'
        x_vals, y_vals = [], []
        if fraglength_csv_path.exists():
            df = pd.read_csv(fraglength_csv_path, encoding='utf8')
            if 'Width' in df.columns and 'Nr_frag' in df.columns:
                x_series = pd.to_numeric(df['Width'], errors='coerce').dropna()
                y_series = pd.to_numeric(df['Nr_frag'], errors='coerce').dropna()
                # Align lengths safely
                min_len = min(len(x_series), len(y_series))
                x_vals = x_series.tolist()[:min_len]
                y_vals = y_series.tolist()[:min_len]
        combined_context['atac_insert_x'] = json.dumps(x_vals)
        combined_context['atac_insert_y'] = json.dumps(y_vals)
    except Exception as e:
        print(f"Warning: Failed to prepare ATAC insert size arrays: {e}")
        combined_context['atac_insert_x'] = '[]'
        combined_context['atac_insert_y'] = '[]'

def _process_atac_saturation(atac_processing_dir, combined_context, mean_reads_str):
    """Provide ATAC saturation arrays"""
    try:
        sampling_stats_path = atac_processing_dir / 'SUMMARY_METRIC' / 'sampling.stats.xls'
        atac_sat_x_vals, atac_sat_y_vals = [], []
        if sampling_stats_path.exists():
            sat_df = pd.read_table(sampling_stats_path, encoding='utf8', sep='\t')
            # Expect columns: sampling_fraction, median_uniq_frag_per_bc
            if 'sampling_fraction' in sat_df.columns and 'median_uniq_frag_per_bc' in sat_df.columns:
                sf_series = pd.to_numeric(sat_df['sampling_fraction'], errors='coerce').fillna(0).clip(lower=0)
                med_uniq_series = pd.to_numeric(sat_df['median_uniq_frag_per_bc'], errors='coerce').fillna(0).clip(lower=0)
                # Mean read pairs per cell comes from ATAC stat
                try:
                    mean_reads_num = int(float(mean_reads_str.replace(',', '').replace('%', '')))
                except Exception:
                    mean_reads_num = 0
                atac_sat_x_vals = (sf_series * mean_reads_num).tolist()
                atac_sat_y_vals = med_uniq_series.tolist()
        combined_context['atac_saturation_x'] = json.dumps(atac_sat_x_vals)
        combined_context['atac_saturation_y'] = json.dumps(atac_sat_y_vals)
        combined_context['atac_saturation_available'] = bool(atac_sat_x_vals and atac_sat_y_vals and len(atac_sat_x_vals) == len(atac_sat_y_vals))
    except Exception as e:
        print(f"Warning: Failed to prepare ATAC saturation arrays: {e}")
        combined_context['atac_saturation_x'] = '[]'
        combined_context['atac_saturation_y'] = '[]'
        combined_context['atac_saturation_available'] = False

def _process_atac_jaccard_index(atac_processing_dir, combined_context, thres_str):
    """Provide ATAC jaccard index distribution arrays"""
    try:
        jaccard_tsv_path = atac_processing_dir / 'JACCARD_MERGE' / 'cb.jaccard.index.tsv'
        true_x, true_y, false_x, false_y = [], [], [], []
        if jaccard_tsv_path.exists():
            jacc_df = pd.read_csv(jaccard_tsv_path, encoding='utf8', sep='\t')
            # Robust rounding and duplicate removal on jaccard values
            if 'jaccard' in jacc_df.columns:
                jacc_df['jaccard'] = pd.to_numeric(jacc_df['jaccard'], errors='coerce').round(4)
                jacc_df = jacc_df.dropna(subset=['jaccard'])
                dp_first = set(jacc_df[jacc_df[["jaccard"]].duplicated(keep="first")].index)
                dp_last = set(jacc_df[jacc_df[["jaccard"]].duplicated(keep="last")].index)
                dp_inter = dp_first & dp_last
                jacc_df = jacc_df.drop(list(dp_inter), axis=0)
            # Determine threshold from existing context or default
            try:
                thres_val = float(thres_str)
            except Exception:
                thres_val = 0.02
            # Column names for rank/x
            rank_col = 'jaccard_rank' if 'jaccard_rank' in jacc_df.columns else None
            if rank_col is None:
                # Fallback: if no rank column, create a 1-based rank by order
                jacc_df = jacc_df.reset_index(drop=True)
                jacc_df['jaccard_rank'] = jacc_df.index + 1
                rank_col = 'jaccard_rank'
            jacc_df[rank_col] = pd.to_numeric(jacc_df[rank_col], errors='coerce').fillna(0).astype(int)
            # Split into TRUE/FALSE traces
            if 'jaccard' in jacc_df.columns:
                df_true = jacc_df[jacc_df['jaccard'] >= thres_val]
                df_false = jacc_df[jacc_df['jaccard'] < thres_val]
                true_x = df_true[rank_col].tolist()
                true_y = df_true['jaccard'].tolist()
                false_x = df_false[rank_col].tolist()
                false_y = df_false['jaccard'].tolist()
        combined_context['atac_jaccard_true_x'] = json.dumps(true_x)
        combined_context['atac_jaccard_true_y'] = json.dumps(true_y)
        combined_context['atac_jaccard_false_x'] = json.dumps(false_x)
        combined_context['atac_jaccard_false_y'] = json.dumps(false_y)
    except Exception as e:
        print(f"Warning: Failed to prepare ATAC jaccard distribution arrays: {e}")
        combined_context['atac_jaccard_true_x'] = '[]'
        combined_context['atac_jaccard_true_y'] = '[]'
        combined_context['atac_jaccard_false_x'] = '[]'
        combined_context['atac_jaccard_false_y'] = '[]'

def _generate_fastq_display_html(sample_outdir, available_libraries):
    logs_dir = Path(sample_outdir) / 'logs'
    
    def read_log_and_extract_reads(log_name, prefixes):
        """Helper to read a log and extract specified read paths."""
        log_file = logs_dir / f'{log_name}_cmd.log'
        extracted_reads = {prefix: '' for prefix in prefixes}
        
        if not log_file.exists():
            return extracted_reads
            
        try:
            with open(log_file, 'r', encoding='utf-8') as f: # Use utf-8 for safety
                for line in f:
                    for prefix in prefixes:
                        if line.startswith(prefix):
                            extracted_reads[prefix] = line.split(prefix, 1)[1].strip()
        except Exception:
            pass
            
        return extracted_reads

    def format_text_for_reads(title, reads_dict):
        """Helper to format a plain text block for a set of reads."""
        valid_reads = {k: v for k, v in reads_dict.items() if v}
        if not valid_reads:
            return ""
        
        # Create a plain text block with indentation
        lines = [f"{title}:"]
        for read_path in valid_reads.values():
            lines.append(f"  {read_path}") # Indent for readability
        return '\n'.join(lines)

    text_parts = []
    
    # Process RNA
    if available_libraries.get('gene-expression', False):
        rna_reads = read_log_and_extract_reads('rna', ['Parsed cDNA Read1: ', 'Parsed cDNA Read2: '])
        oligo_reads = read_log_and_extract_reads('rna', ['Parsed Oligo Read1: ', 'Parsed Oligo Read2: '])
        
        cnda_text = format_text_for_reads('cDNA', rna_reads)
        if cnda_text:
            text_parts.append(cnda_text)
            
        oligo_text = format_text_for_reads('Oligo', oligo_reads)
        if oligo_text:
            text_parts.append(oligo_text)

    # Process other libraries
    for lib_type in ['vdj-t', 'vdj-b', 'atac']:
        if available_libraries.get(lib_type, False):
            reads = read_log_and_extract_reads(lib_type, ['Parsed Read1: ', 'Parsed Read2: '])
            lib_text = format_text_for_reads(lib_type.upper(), reads)
            if lib_text:
                text_parts.append(lib_text)
            
    # Join the different sections with a blank line
    return '\n\n'.join(filter(None, text_parts))


def generate_multi_report(name, outdir, config):
    """
    Generates a multi-omics report.
    """
    print(f"Generating multi-omics report for sample {name} in {outdir}")
    
    sample_outdir = Path(outdir) / name
    
    # Create necessary REPORT directories and outs directory based on config
    outs_dir = create_report_directories(sample_outdir, config)
    
    omics_data = get_omics_data(outdir, name, config)
    
    template_path = Path(__file__).parent.parent.parent / 'config' / 'template' / 'template_multi.html'
    with open(template_path, 'r') as f:
        template_str = f.read()
    
    template = Template(template_str)
    
    combined_context = {}
    sample_outdir = Path(outdir) / name
    
    # Process each omics data type
    _process_rna_data(omics_data, sample_outdir, combined_context)
    _process_vdj_t_data(omics_data, sample_outdir, combined_context)
    _process_vdj_b_data(omics_data, sample_outdir, combined_context)
    _process_atac_data(omics_data, sample_outdir, combined_context)
    
    # Add available libraries information
    available_libraries = {
        'gene-expression': 'rna' in omics_data,
        'vdj-t': 'vdj-t' in omics_data,
        'vdj-b': 'vdj-b' in omics_data,
        'atac': 'atac' in omics_data
    }
    combined_context['available_libraries'] = json.dumps(available_libraries)

    end5_ok = _as_bool(combined_context.get('rna_end5', False)) if available_libraries['gene-expression'] else False
    vdj_t_has_beadstrans = bool(config.get('vdj-t', {}).get('beadstrans'))
    vdj_b_has_beadstrans = bool(config.get('vdj-b', {}).get('beadstrans'))
    vdj_t_target_enabled = available_libraries['gene-expression'] and end5_ok and available_libraries['vdj-t'] and not vdj_t_has_beadstrans
    vdj_b_target_enabled = available_libraries['gene-expression'] and end5_ok and available_libraries['vdj-b'] and not vdj_b_has_beadstrans
    combined_context['vdj_t_target_enabled'] = vdj_t_target_enabled
    combined_context['vdj_b_target_enabled'] = vdj_b_target_enabled
    combined_context['vdj_t_target_display'] = '' if vdj_t_target_enabled else 'display:none;'
    combined_context['vdj_b_target_display'] = '' if vdj_b_target_enabled else 'display:none;'

    if vdj_t_target_enabled:
        combined_context['vdj_t_target_umap_x'] = combined_context.get('rna_umap_x', '[]')
        combined_context['vdj_t_target_umap_y'] = combined_context.get('rna_umap_y', '[]')
        combined_context['vdj_t_target_leiden'] = combined_context.get('rna_leiden', '[]')
        combined_context['vdj_t_target_pred_celltype'] = combined_context.get('rna_pred_celltype', '[]')
        combined_context['vdj_t_target_flags'] = combined_context.get('vdj_t_flags', '[]')
    if vdj_b_target_enabled:
        combined_context['vdj_b_target_umap_x'] = combined_context.get('rna_umap_x', '[]')
        combined_context['vdj_b_target_umap_y'] = combined_context.get('rna_umap_y', '[]')
        combined_context['vdj_b_target_leiden'] = combined_context.get('rna_leiden', '[]')
        combined_context['vdj_b_target_pred_celltype'] = combined_context.get('rna_pred_celltype', '[]')
        combined_context['vdj_b_target_flags'] = combined_context.get('vdj_b_flags', '[]')

    combined_context['fastq_display_html'] = _generate_fastq_display_html(sample_outdir, available_libraries)

    # Add sample name and version
    combined_context['samplename'] = name
    combined_context['version'] = __version__

    # Add input CSV and config parameters
    combined_context['input_csv_data'] = config.get('csv_content', '')
    config_for_display = config.copy()
    config_for_display.pop('csv_content', None)
    combined_context['config_parameters'] = json.dumps(config_for_display, indent=4)

    template_dir = template_path.parent
    plotly_candidates = [
        template_dir / 'plotly.js',
        template_dir / 'plotly-2.26.0.min.js',
        template_dir / 'plotly.min.js'
    ]
    raw_plotly_js = ''
    for js_path in plotly_candidates:
        if js_path.exists():
            try:
                with open(js_path, 'r', encoding='utf-8') as f:
                    raw_plotly_js = f.read()
                break
            except Exception:
                pass
    if raw_plotly_js:
        combined_context['plotly_loader_tag'] = '<script>' + raw_plotly_js + '</script>'
    else:
        combined_context['plotly_loader_tag'] = '<script src="https://cdn.plot.ly/plotly-2.26.0.min.js"></script>'

    report_html = template.safe_substitute(combined_context)
    
    with open(outs_dir/ f'{name}_multi_report.html', 'w') as f:
        f.write(report_html)
    
    print(f"Multi-omics report saved to: {outs_dir/ f'{name}_multi_report.html'}")
    print("Multi-omics report generation complete.")
