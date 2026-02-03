import json
import csv
import statistics
import gzip
from pathlib import Path
from string import Template
import re

# Hardcoded version since src/__init__.py might not exist
__version__ = "1.0.0"

_JS_TEMPLATE_PLACEHOLDERS = {
    "id",
    "plotId",
    "clusterNum",
    "libraryId",
    "coef",
    "exp",
    "cluster",
    "emptyMessage",
}

_JSON_ARRAY_PLACEHOLDERS = {
    "bead_count_x",
    "bead_count_y",
    "rna_umap_x",
    "rna_umap_y",
    "rna_leiden",
    "rna_numi",
    "rna_clusters",
    "rna_gene_body_percentile",
    "rna_gene_body_umi_maxnorm",
    "rna_gene_body_internal_maxnorm",
    "rna_saturation_fraction",
    "rna_saturation_lib_pct",
    "rna_saturation_gene_pct",
    "rna_saturation_median_genes_umi",
    "rna_saturation_median_genes_read",
    "rna_saturation_sampling_fraction",
    "rna_saturation_seq_pct",
    "rna_saturation_median_genes",
    "rna_rankplot_data",
    "rna_stats_table_data",
    "rna_marker_data",
}

_JSON_OBJECT_PLACEHOLDERS = {
    "rna_read_distribution_bar_data",
    "rna_read_distribution_box_data",
}

_CANONICAL_STATS_KEYS = {
    "exon_reads": "Exon_reads",
    "intron_reads": "intron_reads",
    "intergenic_reads": "Intergenic_reads",
    "ambiguity_reads": "Ambiguity_reads",
    "unmapped_reads": "Unmapped_reads",
    "antisense_reads": "Antisense_reads",
    "exon_genes": "Exon_genes",
    "intron_genes": "Intron_genes",
    "intron_exon_genes": "Intron_Exon_genes",
    "exon_umis": "Exon_umis",
    "intron_umis": "Intron_umis",
    "intron_exon_umis": "Intron_Exon_umis",
    "umi_reads": "umi_reads",
    "internal_reads": "internal_reads",
    "wellid": "wellID",
}


def _default_placeholder_value(name):
    if name in _JSON_ARRAY_PLACEHOLDERS:
        return "[]"
    if name in _JSON_OBJECT_PLACEHOLDERS:
        return "{}"
    if name.endswith("_enabled") or name.endswith("_available"):
        return "false"
    if name.endswith("_html"):
        return ""
    if name.endswith("_data"):
        return "[]"
    if name.endswith("_pct") or name.endswith("_percent") or name.endswith("_percentage"):
        return ""
    if name == "rna_species":
        return "Unknown"
    if name == "rna_fraction_transcriptome_in_cells":
        return ""
    if name == "rna_saturation":
        return ""
    return ""


def _normalize_stats_record(record):
    if not isinstance(record, dict):
        return {}
    out = {}
    for k, v in record.items():
        if k is None:
            continue
        key = str(k).strip()
        out[key] = v
        lk = key.lower()
        canon = _CANONICAL_STATS_KEYS.get(lk)
        if canon:
            if canon not in out or out.get(canon) in (None, ""):
                out[canon] = v
    return out


def _count_lines_in_tsv_gz(path):
    try:
        p = Path(path)
        if not p.exists():
            return None
        with gzip.open(p, "rt", encoding="utf-8", errors="ignore") as f:
            n = 0
            for _ in f:
                n += 1
            return n
    except Exception:
        return None


def _infer_total_genes_from_expression(sample_outdir, project):
    if not project:
        return None
    base = Path(sample_outdir) / "zUMIs_output" / "expression"
    candidates = [
        base / f"{project}.inex.umi" / "genes.tsv.gz",
        base / f"{project}.inex.umi" / "features.tsv.gz",
        base / f"{project}.exon.umi" / "genes.tsv.gz",
        base / f"{project}.exon.umi" / "features.tsv.gz",
    ]
    for p in candidates:
        n = _count_lines_in_tsv_gz(p)
        if n is not None and n > 0:
            return n
    return None


def _process_rna_cluster_assignment(_sample_outdir, combined_context):
    combined_context["rna_umap_x"] = "[]"
    combined_context["rna_umap_y"] = "[]"
    combined_context["rna_leiden"] = "[]"
    combined_context["rna_numi"] = "[]"
    combined_context["rna_clusters"] = "[]"
    combined_context["cell_annotation_available"] = "false"
    return


def _process_rna_bead_count_data(_sample_outdir, combined_context):
    combined_context["bead_count_x"] = "[]"
    combined_context["bead_count_y"] = "[]"
    return


def _process_rna_rankplot_data(_sample_outdir, combined_context):
    combined_context["rna_rankplot_data"] = "[]"
    return


def infer_enabled_sections(config):
    """
    Infer which sections (rna, vdj, etc.) are enabled based on config
    or directory structure.
    """
    enabled = []
    # For now, we mainly support RNA in this pipeline
    if config.get('sequence_files'):
        enabled.append('rna')
    
    # Check if VDJ/ATAC outputs exist (future proofing)
    if 'out_dir' in config:
        out_dir = Path(config['out_dir'])
        if (out_dir / 'VDJ-T_ANALYSIS_WORKFLOW_PROCESSING').exists():
            enabled.append('vdj-t')
        if (out_dir / 'VDJ-B_ANALYSIS_WORKFLOW_PROCESSING').exists():
            enabled.append('vdj-b')
        if (out_dir / 'ATAC_ANALYSIS_WORKFLOW_PROCESSING').exists():
            enabled.append('atac')
            
    return enabled

def create_report_directories(sample_outdir, _config=None):
    """
    Create necessary REPORT directory structure for multi-omics analysis
    based on the configuration
    """
    sample_path = Path(sample_outdir)

    # For Mhsflt_toolkit, we might not need strict directory structure for REPORT 
    # if we are just generating one HTML, but we'll keep it for consistency.
    # The pipeline outputs are in zUMIs_output mostly.
    
    # Create final output directory 'outs'
    # In the pipeline, out_dir is XPRESS_PROCESSING. 
    # sample_outdir passed here is likely XPRESS_PROCESSING.
    # The final report should go to 'outs' in the project root.
    
    project_root = sample_path.parent # Assuming sample_outdir is project/XPRESS_PROCESSING
    outs_dir = project_root / 'outs'
    outs_dir.mkdir(parents=True, exist_ok=True)
    
    return outs_dir

def calculate_summary_metrics(sample_outdir, project=""):
    """
    Calculate summary metrics from stats.tsv
    """
    stats_tsv = list((sample_outdir / 'zUMIs_output' / 'stats').glob('*.stats.tsv'))
    if not stats_tsv:
        stats_tsv = list(sample_outdir.rglob('*.stats.tsv'))
        
    metrics = {
        "rna_estm_Num_cell": "",
        "rna_median_umis_percell": "",
        "rna_mean_umis_percell": "",
        "rna_median_genes_percell": "",
        "rna_mean_genes_percell": "",
        "rna_mean_reads_percell": "",
        "rna_total_gene": "",
        "rna_fraction_transcriptome_in_cells": "",
        "rna_species": "Unknown",
    }
    
    if not stats_tsv:
        print("Warning: No stats.tsv found for summary metrics.")
        return metrics
        
    def to_float(v):
        if v is None:
            return None
        s = str(v).strip()
        if not s:
            return None
        try:
            return float(s.replace(",", ""))
        except Exception:
            return None

    def fmt_int(n):
        try:
            return f"{int(n):,}"
        except Exception:
            return ""

    try:
        stats_path = stats_tsv[0]
        with open(stats_path, "r", encoding="utf-8", errors="ignore", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = [_normalize_stats_record(r) for r in reader if r]
        if not rows:
            return metrics

        metrics["rna_estm_Num_cell"] = fmt_int(len(rows))

        umis = []
        genes = []
        reads = []
        exon_reads = []
        intron_reads = []
        intergenic_reads = []
        ambiguity_reads = []
        unmapped_reads = []

        for r in rows:
            u = to_float(r.get("Intron_Exon_umis")) if r.get("Intron_Exon_umis") is not None else None
            if u is None:
                u = to_float(r.get("Exon_umis")) if r.get("Exon_umis") is not None else None
            if u is not None:
                umis.append(u)

            g = to_float(r.get("Intron_Exon_genes")) if r.get("Intron_Exon_genes") is not None else None
            if g is None:
                g = to_float(r.get("Exon_genes")) if r.get("Exon_genes") is not None else None
            if g is not None:
                genes.append(g)

            ru = to_float(r.get("umi_reads")) if r.get("umi_reads") is not None else 0.0
            ri = to_float(r.get("internal_reads")) if r.get("internal_reads") is not None else 0.0
            if ru is not None or ri is not None:
                reads.append(float(ru or 0.0) + float(ri or 0.0))

            v = to_float(r.get("Exon_reads"))
            if v is not None:
                exon_reads.append(v)
            v = to_float(r.get("intron_reads"))
            if v is not None:
                intron_reads.append(v)
            v = to_float(r.get("Intergenic_reads"))
            if v is not None:
                intergenic_reads.append(v)
            v = to_float(r.get("Ambiguity_reads"))
            if v is not None:
                ambiguity_reads.append(v)
            v = to_float(r.get("Unmapped_reads"))
            if v is not None:
                unmapped_reads.append(v)

        if umis:
            metrics["rna_median_umis_percell"] = fmt_int(statistics.median(umis))
            metrics["rna_mean_umis_percell"] = fmt_int(sum(umis) / len(umis))
        if genes:
            metrics["rna_median_genes_percell"] = fmt_int(statistics.median(genes))
            metrics["rna_mean_genes_percell"] = fmt_int(sum(genes) / len(genes))
        if reads:
            metrics["rna_mean_reads_percell"] = fmt_int(sum(reads) / len(reads))

        if exon_reads or intron_reads or intergenic_reads or ambiguity_reads or unmapped_reads:
            exon_sum = sum(exon_reads) if exon_reads else 0.0
            intron_sum = sum(intron_reads) if intron_reads else 0.0
            intergenic_sum = sum(intergenic_reads) if intergenic_reads else 0.0
            ambig_sum = sum(ambiguity_reads) if ambiguity_reads else 0.0
            unmapped_sum = sum(unmapped_reads) if unmapped_reads else 0.0
            denom = exon_sum + intron_sum + intergenic_sum + ambig_sum + unmapped_sum
            transcriptome = exon_sum + intron_sum
            if denom > 0:
                metrics["rna_fraction_transcriptome_in_cells"] = f"{(transcriptome / denom) * 100.0:.2f}%"

        total_genes = _infer_total_genes_from_expression(sample_outdir, str(project))
        if total_genes is not None:
            metrics["rna_total_gene"] = fmt_int(total_genes)
    except Exception as e:
        print(f"Error calculating summary metrics: {e}")
        
    return metrics

def get_omics_data(outdir, _sample, config):
    omics_data = {}
    outdir = Path(outdir)
    # In run_analysis_pipeline, outdir is passed as the XPRESS_PROCESSING directory.
    
    enabled = set(infer_enabled_sections(config))

    # RNA
    if 'rna' in enabled:
        try:
            # We calculate stats manually now
            project = config.get("project") or _sample or ""
            stat = calculate_summary_metrics(outdir, project=str(project))
            
            # Metadata
            rna_section = config.get('reference', {}) # In new pipeline config, genomeDir is in reference
            # Try to guess species from genome dir path
            genome_dir = rna_section.get('STAR_index', '')
            species = 'Unknown'
            if 'human' in str(genome_dir).lower() or 'hg' in str(genome_dir).lower():
                species = 'Homo sapiens'
            elif 'mouse' in str(genome_dir).lower() or 'mm' in str(genome_dir).lower():
                species = 'Mus musculus'
            stat['rna_species'] = species
            
            omics_data['rna'] = {'stat': stat, 'plot_dict': {}, 'table': None}
            
        except Exception as e:
            print(f"Warning: Failed to load RNA report data: {e}")

    return omics_data

# Helper functions for processing omics data
def _process_rna_data(omics_data, sample_outdir, combined_context):
    """Process RNA data and add to combined context"""
    if 'rna' not in omics_data:
        return
        
    # Add basic RNA statistics
    for key, value in omics_data['rna']['stat'].items():
        combined_context[key] = value
        
    # Process RNA-specific data
    _process_rna_stats_table_data(sample_outdir, combined_context)
    _process_rna_saturation_data(sample_outdir, combined_context)
    _process_rna_gene_body_coverage_data(sample_outdir, combined_context)
    _process_rna_read_distribution_data(sample_outdir, combined_context)
    _process_rna_bead_count_data(sample_outdir, combined_context)
    _process_rna_cluster_assignment(sample_outdir, combined_context)
    _process_rna_rankplot_data(sample_outdir, combined_context)

def _process_rna_stats_table_data(sample_outdir, combined_context):
    combined_context['rna_stats_table_data'] = '[]'
    combined_context['rna_stats_table_available'] = False
    try:
        candidates = list((sample_outdir / 'zUMIs_output' / 'stats').glob('*.stats.tsv'))
        if not candidates:
            candidates = list(sample_outdir.rglob('*.stats.tsv'))
        if not candidates:
            return

        stats_path = candidates[0]
        with open(stats_path, "r", encoding="utf-8", errors="ignore", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            records = [_normalize_stats_record(r) for r in reader if r]
        combined_context["rna_stats_table_data"] = json.dumps(records)
        combined_context["rna_stats_table_available"] = bool(len(records) > 0)
    except Exception as e:
        print(f"Warning: Failed to prepare RNA stats table data: {e}")

def _process_rna_gene_body_coverage_data(sample_outdir, combined_context):
    combined_context['rna_gene_body_percentile'] = '[]'
    combined_context['rna_gene_body_umi_maxnorm'] = '[]'
    combined_context['rna_gene_body_internal_maxnorm'] = '[]'
    combined_context['rna_gene_body_coverage_available'] = False
    try:
        candidates = list((sample_outdir / 'zUMIs_output' / 'stats').glob('*.geneBodyCoverage.txt'))
        if not candidates:
            candidates = list(sample_outdir.rglob('*.geneBodyCoverage.txt'))
        if not candidates:
            return

        gb_path = candidates[0]
        x_vals = []
        umi_vals = []
        internal_vals = []
        with open(gb_path, "r", encoding="utf-8", errors="ignore", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                if not r:
                    continue
                try:
                    x = float(str(r.get("Percentile", "")).strip() or 0)
                    umi = float(str(r.get("UMI_MaxNorm", "")).strip() or 0)
                    internal = float(str(r.get("Internal_MaxNorm", "")).strip() or 0)
                except Exception:
                    continue
                x_vals.append(x)
                umi_vals.append(umi)
                internal_vals.append(internal)

        combined_context["rna_gene_body_percentile"] = json.dumps(x_vals)
        combined_context["rna_gene_body_umi_maxnorm"] = json.dumps(umi_vals)
        combined_context["rna_gene_body_internal_maxnorm"] = json.dumps(internal_vals)
        combined_context["rna_gene_body_coverage_available"] = bool(len(x_vals) > 0 and len(x_vals) == len(umi_vals) == len(internal_vals))
    except Exception as e:
        print(f"Warning: Failed to prepare RNA gene body coverage data: {e}")

def _process_rna_read_distribution_data(sample_outdir, combined_context):
    combined_context['rna_read_distribution_bar_data'] = '{}'
    combined_context['rna_read_distribution_box_data'] = '{}'
    try:
        candidates = list((sample_outdir / 'zUMIs_output' / 'stats').glob('*.read_stats.json'))
        if not candidates:
            candidates = list(sample_outdir.rglob('*.read_stats.json'))
        read_stats = {}
        if candidates:
            stats_path = candidates[0]
            with open(stats_path, 'r', encoding='utf-8') as f:
                stats_json = json.load(f)
            read_stats = stats_json.get('read_stats', {}) or {}

        categories_order = ["Exon", "Intron", "Intergenic", "Ambiguity", "Unmapped"]
        feat_colors = {
            "Exon": "#1A5084",
            "Intron": "#118730",
            "Intergenic": "#FFD700",
            "Ambiguity": "#FFA54F",
            "Unmapped": "#545454",
            "Unused BC": "#BABABA",
        }

        total_fracs_pct = {c: 0.0 for c in categories_order}
        unused_frac_pct = 0.0
        
        per_cell_pct = {c: [] for c in categories_order}
        
        # Calculate totals from read_stats.json
        if isinstance(read_stats, dict) and read_stats:
            total_barcodes = [k for k in read_stats.keys() if isinstance(k, str) and not k.startswith('__')]
            unused_total = 0
            try:
                unused_total = int((read_stats.get("__NO_CB__", {}) or {}).get("Unused BC", 0))
            except Exception:
                unused_total = 0

            totals = {c: 0 for c in categories_order}
            
            # Calculate per-cell percentages and total sums
            for bc in total_barcodes:
                counts = read_stats.get(bc, {}) or {}
                cell_total = 0
                for c in categories_order:
                    val = int(counts.get(c, 0) or 0)
                    totals[c] += val
                    cell_total += val
                
                if cell_total > 0:
                    for c in categories_order:
                        val = int(counts.get(c, 0) or 0)
                        per_cell_pct[c].append((val / cell_total) * 100.0)

            total_sum = sum(totals.values()) + unused_total
            total_fracs_pct = {c: (totals[c] / total_sum * 100.0) if total_sum > 0 else 0.0 for c in categories_order}
            unused_frac_pct = (unused_total / total_sum * 100.0) if total_sum > 0 else 0.0

        bar_payload = {
            "order": categories_order + ["Unused BC"],
            "values": {**total_fracs_pct, "Unused BC": unused_frac_pct},
            "colors": feat_colors,
        }
        box_payload = {
            "order": categories_order,
            "values": per_cell_pct,
            "colors": feat_colors,
        }

        combined_context['rna_read_distribution_bar_data'] = json.dumps(bar_payload)
        combined_context['rna_read_distribution_box_data'] = json.dumps(box_payload)
    except Exception as e:
        print(f"Warning: Failed to prepare RNA read distribution data: {e}")

def _process_rna_saturation_data(sample_outdir, combined_context):
    """Prepare RNA saturation arrays for plotting"""
    try:
        combined_context['rna_saturation_fraction'] = '[]'
        combined_context['rna_saturation_lib_pct'] = '[]'
        combined_context['rna_saturation_gene_pct'] = '[]'
        combined_context['rna_saturation_median_genes_umi'] = '[]'
        combined_context['rna_saturation_median_genes_read'] = '[]'
        combined_context['rna_saturation_available'] = False

        candidates = list((sample_outdir / 'zUMIs_output' / 'stats').glob('*.saturation.tsv'))
        if not candidates:
            candidates = list(sample_outdir.rglob('*.saturation.tsv'))
        if candidates:
            sat_path = candidates[0]
            frac_vals = []
            sat_lib_pct = []
            sat_gene_pct = []
            med_umi_vals = []
            med_read_vals = []
            with open(sat_path, "r", encoding="utf-8", errors="ignore", newline="") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for r in reader:
                    if not r:
                        continue
                    try:
                        frac = float(str(r.get("Fraction", "")).strip() or 0)
                        lib = float(str(r.get("Seq_Saturation_Library", "")).strip() or 0) * 100.0
                        gene = float(str(r.get("Seq_Saturation_Gene", "")).strip() or 0) * 100.0
                        med_umi = float(str(r.get("Median_Genes_UMI", "")).strip() or 0)
                        med_read = float(str(r.get("Median_Genes_Read", "")).strip() or 0)
                    except Exception:
                        continue
                    frac_vals.append(frac)
                    sat_lib_pct.append(lib)
                    sat_gene_pct.append(gene)
                    med_umi_vals.append(med_umi)
                    med_read_vals.append(med_read)

            if frac_vals and len(frac_vals) == len(sat_lib_pct) == len(sat_gene_pct):
                combined_context["rna_saturation_fraction"] = json.dumps(frac_vals)
                combined_context["rna_saturation_lib_pct"] = json.dumps(sat_lib_pct)
                combined_context["rna_saturation_gene_pct"] = json.dumps(sat_gene_pct)
                combined_context["rna_saturation_median_genes_umi"] = json.dumps(med_umi_vals)
                combined_context["rna_saturation_median_genes_read"] = json.dumps(med_read_vals)
                combined_context["rna_saturation_available"] = True

        combined_context['rna_saturation_sampling_fraction'] = combined_context.get('rna_saturation_fraction', '[]')
        combined_context['rna_saturation_seq_pct'] = combined_context.get('rna_saturation_lib_pct', '[]')
        combined_context['rna_saturation_median_genes'] = combined_context.get('rna_saturation_median_genes_umi', '[]')
        try:
            lib_pct = json.loads(combined_context.get("rna_saturation_lib_pct", "[]") or "[]")
            if isinstance(lib_pct, list) and lib_pct:
                combined_context["rna_saturation"] = f"{float(lib_pct[-1]):.1f}%"
        except Exception:
            pass
    except Exception as e:
        print(f"Warning: Failed to prepare RNA saturation arrays: {e}")


def generate_multi_report(name, outdir, config):
    """
    Generates a multi-omics report.
    """
    print(f"Generating multi-omics report for sample {name} in {outdir}")
    
    sample_outdir = Path(outdir)
    
    # Create necessary REPORT directories and outs directory based on config
    outs_dir = create_report_directories(sample_outdir, config)
    
    omics_data = get_omics_data(outdir, name, config)
    
    # Locate template
    # Assuming report.py is in src/, and template is in report/
    template_path = Path(__file__).parent.parent / 'report' / 'template_multi.html'
    
    if not template_path.exists():
        print(f"Error: Template not found at {template_path}")
        return

    with open(template_path, 'r', encoding='utf-8') as f:
        template_str = f.read()
    
    template = Template(template_str)
    
    combined_context = {}
    
    # Process each omics data type
    # Add sample name and version
    combined_context['samplename'] = name
    combined_context['version'] = __version__
    combined_context["sample_type"] = str((config.get("sample") or {}).get("sample_type") or "").strip().lower()
    combined_context["vdj_t_target_enabled"] = "false"
    combined_context["vdj_b_target_enabled"] = "false"
    combined_context["fastq_display_html"] = ""

    # Add input CSV and config parameters
    combined_context['input_csv_data'] = config.get('csv_content', '')
    config_for_display = config.copy()
    config_for_display.pop('csv_content', None)
    combined_context['config_parameters'] = json.dumps(config_for_display, indent=4, default=str)

    # Handle Plotly JS loading
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
        # Fallback to CDN if local file not found
        combined_context['plotly_loader_tag'] = '<script src="https://cdn.plot.ly/plotly-2.26.0.min.js"></script>'

    _process_rna_data(omics_data, sample_outdir, combined_context)

    identifiers = set(re.findall(r"\$\{([a-zA-Z0-9_]+)\}", template_str))
    for key in identifiers:
        if key in _JS_TEMPLATE_PLACEHOLDERS:
            continue
        if key not in combined_context:
            combined_context[key] = _default_placeholder_value(key)

    try:
        report_html = template.safe_substitute(combined_context)
        
        out_file = outs_dir / f'{name}_multi_report.html'
        with open(out_file, 'w', encoding='utf-8') as f:
            f.write(report_html)
        
        print(f"Multi-omics report saved to: {out_file}")
        print("Multi-omics report generation complete.")
    except Exception as e:
        print(f"Error substituting template: {e}")
