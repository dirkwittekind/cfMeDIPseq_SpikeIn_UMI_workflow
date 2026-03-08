#!/usr/bin/env python3
"""
MEDIPIPE Sample QC Report Generator

Generates per-sample HTML/TSV reports with:
- Readcount (deduplicated reads)
- Coverage statistics
- Overamplification QC
- UMI metrics (unique fragments)
- Spike-In QC values
- Fragmentomics summary
- Tissue of origin deconvolution

Author: MEDIPIPE
"""

import argparse
import json
import os
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
import numpy as np


def load_bam_stats(stats_file: str) -> Dict:
    """Load BAM stats from samtools stats output."""
    stats = {}
    if not os.path.exists(stats_file):
        return stats
    
    with open(stats_file) as f:
        for line in f:
            if line.startswith('SN'):
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    key = parts[1].rstrip(':')
                    value = parts[2]
                    try:
                        stats[key] = int(value)
                    except ValueError:
                        try:
                            stats[key] = float(value)
                        except ValueError:
                            stats[key] = value
    return stats


def load_coverage_stats(coverage_file: str) -> Dict:
    """Load coverage statistics from mosdepth output."""
    stats = {}
    if not os.path.exists(coverage_file):
        return stats
    
    # Try mosdepth summary format
    try:
        df = pd.read_csv(coverage_file, sep='\t')
        if 'mean' in df.columns and 'chrom' in df.columns:
            genome_row = df[df['chrom'] == 'total']
            if len(genome_row) > 0:
                stats['mean_coverage'] = genome_row['mean'].values[0]
                stats['min_coverage'] = genome_row.get('min', [0]).values[0] if 'min' in genome_row else 0
                stats['max_coverage'] = genome_row.get('max', [0]).values[0] if 'max' in genome_row else 0
    except:
        pass
    
    return stats


def load_overamp_metrics(overamp_file: str) -> Dict:
    """Load overamplification/duplication metrics."""
    stats = {}
    if not os.path.exists(overamp_file):
        return stats
    
    try:
        with open(overamp_file) as f:
            content = f.read()
            
        # Parse MEDIPIPE overamp format
        for line in content.strip().split('\n'):
            if ':' in line:
                key, value = line.split(':', 1)
                key = key.strip()
                value = value.strip()
                try:
                    stats[key] = float(value.replace('%', ''))
                except:
                    stats[key] = value
    except:
        pass
    
    return stats


def load_umi_metrics(umi_file: str) -> Dict:
    """Load UMI metrics from deduplication output."""
    stats = {}
    if not os.path.exists(umi_file):
        return stats
    
    try:
        # Try JSON format (umi_tools output)
        if umi_file.endswith('.json'):
            with open(umi_file) as f:
                stats = json.load(f)
        else:
            # Try TSV format
            df = pd.read_csv(umi_file, sep='\t')
            if len(df) > 0:
                stats = df.iloc[0].to_dict()
    except:
        pass
    
    return stats


def load_spikein_metrics(spikein_file: str) -> Dict:
    """Load spike-in normalization metrics."""
    stats = {}
    if not os.path.exists(spikein_file):
        return stats
    
    try:
        df = pd.read_csv(spikein_file, sep='\t')
        if len(df) > 0:
            stats = df.iloc[0].to_dict()
    except:
        pass
    
    return stats


def load_fragmentomics(frag_dir: str, sample_id: str) -> Dict:
    """Load fragmentomics metrics from fragment profile output."""
    stats = {}
    
    # Check for insert size metrics
    insert_file = os.path.join(frag_dir, f"{sample_id}_insert_size_metrics.txt")
    if os.path.exists(insert_file):
        try:
            with open(insert_file) as f:
                lines = f.readlines()
            
            # Parse Picard-style metrics
            for i, line in enumerate(lines):
                if line.startswith('MEDIAN_INSERT_SIZE'):
                    header = line.strip().split('\t')
                    values = lines[i+1].strip().split('\t') if i+1 < len(lines) else []
                    for h, v in zip(header, values):
                        try:
                            stats[h.lower()] = float(v)
                        except:
                            stats[h.lower()] = v
                    break
        except:
            pass
    
    # Check for fragment length distribution
    frag_dist_file = os.path.join(frag_dir, f"{sample_id}_fragment_length_hist.tsv")
    if os.path.exists(frag_dist_file):
        try:
            df = pd.read_csv(frag_dist_file, sep='\t')
            if 'count' in df.columns:
                stats['total_fragments'] = df['count'].sum()
                # Calculate mono/di-nucleosomal ratios
                if 'length' in df.columns:
                    mono = df[(df['length'] >= 120) & (df['length'] <= 180)]['count'].sum()
                    di = df[(df['length'] >= 250) & (df['length'] <= 400)]['count'].sum()
                    total = df['count'].sum()
                    if total > 0:
                        stats['mono_nucleosomal_pct'] = mono / total * 100
                        stats['di_nucleosomal_pct'] = di / total * 100
        except:
            pass
    
    return stats


def load_tissue_of_origin(tof_file: str) -> Dict:
    """Load tissue of origin deconvolution results."""
    stats = {}
    if not os.path.exists(tof_file):
        return stats
    
    try:
        df = pd.read_csv(tof_file, sep='\t')
        if len(df) > 0:
            # Get top tissues
            if 'tissue' in df.columns and 'fraction' in df.columns:
                df = df.sort_values('fraction', ascending=False)
                for i, row in df.head(5).iterrows():
                    stats[f"tissue_{i+1}"] = row['tissue']
                    stats[f"fraction_{i+1}"] = row['fraction']
    except:
        pass
    
    return stats


def collect_sample_metrics(sample_id: str, output_dir: str, config: dict) -> Dict:
    """Collect all QC metrics for a sample."""
    metrics = {
        'sample_id': sample_id,
        'report_date': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    }
    
    # BAM stats (readcount)
    bam_stats_file = os.path.join(output_dir, 'dedup_bam_pe', f'{sample_id}.bam.stats.txt')
    if not os.path.exists(bam_stats_file):
        bam_stats_file = os.path.join(output_dir, 'dedup_bam_pe', f'{sample_id}_dedup.bam.stats.txt')
    bam_stats = load_bam_stats(bam_stats_file)
    if bam_stats:
        metrics['total_reads'] = bam_stats.get('raw total sequences', 0)
        metrics['mapped_reads'] = bam_stats.get('reads mapped', 0)
        metrics['duplicate_reads'] = bam_stats.get('reads duplicated', 0)
        metrics['deduplicated_reads'] = metrics['mapped_reads'] - metrics.get('duplicate_reads', 0)
        if metrics['mapped_reads'] > 0:
            metrics['duplication_rate_pct'] = metrics['duplicate_reads'] / metrics['mapped_reads'] * 100
    
    # Coverage
    coverage_file = os.path.join(output_dir, 'coverage', f'{sample_id}.mosdepth.summary.txt')
    coverage_stats = load_coverage_stats(coverage_file)
    metrics.update({f"coverage_{k}": v for k, v in coverage_stats.items()})
    
    # Overamplification
    overamp_file = os.path.join(output_dir, 'meth_qc_quant', f'{sample_id}_meth_qc.txt')
    overamp_stats = load_overamp_metrics(overamp_file)
    metrics.update({f"overamp_{k}": v for k, v in overamp_stats.items()})
    
    # UMI metrics
    umi_file = os.path.join(output_dir, 'umi_dedup', f'{sample_id}_umi_stats.json')
    umi_stats = load_umi_metrics(umi_file)
    if umi_stats:
        metrics['unique_umi_count'] = umi_stats.get('unique_umis', 0)
        metrics['umi_duplication_rate'] = umi_stats.get('duplication_rate', 0)
    
    # Spike-in QC
    spikein_file = os.path.join(output_dir, 'spikein', f'{sample_id}_spikein_stats.tsv')
    spikein_stats = load_spikein_metrics(spikein_file)
    metrics.update({f"spikein_{k}": v for k, v in spikein_stats.items()})
    
    # Fragmentomics
    frag_dir = os.path.join(output_dir, 'fragment_size')
    frag_stats = load_fragmentomics(frag_dir, sample_id)
    metrics.update({f"frag_{k}": v for k, v in frag_stats.items()})
    
    # Tissue of origin
    tof_file = os.path.join(output_dir, 'tissue_of_origin', f'{sample_id}_tissue_fractions.tsv')
    tof_stats = load_tissue_of_origin(tof_file)
    metrics.update({f"tof_{k}": v for k, v in tof_stats.items()})
    
    return metrics


def generate_html_report(metrics: Dict, output_file: str):
    """Generate HTML report for a sample."""
    
    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>MEDIPIPE Sample Report - {metrics['sample_id']}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; background-color: #f5f5f5; }}
        .container {{ max-width: 1200px; margin: 0 auto; background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        h1 {{ color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }}
        h2 {{ color: #34495e; margin-top: 30px; }}
        .metric-section {{ background: #ecf0f1; padding: 15px; border-radius: 5px; margin: 10px 0; }}
        .metric {{ display: inline-block; width: 200px; margin: 5px; padding: 10px; background: white; border-radius: 3px; }}
        .metric-label {{ font-size: 12px; color: #7f8c8d; }}
        .metric-value {{ font-size: 18px; font-weight: bold; color: #2c3e50; }}
        .status-good {{ color: #27ae60; }}
        .status-warn {{ color: #f39c12; }}
        .status-bad {{ color: #e74c3c; }}
        table {{ border-collapse: collapse; width: 100%; margin: 15px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #3498db; color: white; }}
        tr:nth-child(even) {{ background-color: #f2f2f2; }}
        .footer {{ margin-top: 30px; padding-top: 15px; border-top: 1px solid #ddd; color: #7f8c8d; font-size: 12px; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 MEDIPIPE Sample QC Report</h1>
        <p><strong>Sample ID:</strong> {metrics['sample_id']}</p>
        <p><strong>Report Generated:</strong> {metrics['report_date']}</p>
        
        <h2>📊 Read Statistics</h2>
        <div class="metric-section">
            <div class="metric">
                <div class="metric-label">Total Reads</div>
                <div class="metric-value">{metrics.get('total_reads', 'N/A'):,}</div>
            </div>
            <div class="metric">
                <div class="metric-label">Mapped Reads</div>
                <div class="metric-value">{metrics.get('mapped_reads', 'N/A'):,}</div>
            </div>
            <div class="metric">
                <div class="metric-label">Deduplicated Reads</div>
                <div class="metric-value">{metrics.get('deduplicated_reads', 'N/A'):,}</div>
            </div>
            <div class="metric">
                <div class="metric-label">Duplication Rate</div>
                <div class="metric-value">{metrics.get('duplication_rate_pct', 0):.1f}%</div>
            </div>
        </div>
        
        <h2>📈 Coverage</h2>
        <div class="metric-section">
            <div class="metric">
                <div class="metric-label">Mean Coverage</div>
                <div class="metric-value">{metrics.get('coverage_mean_coverage', 'N/A')}</div>
            </div>
        </div>
        
        <h2>🔬 Overamplification QC</h2>
        <div class="metric-section">
"""
    
    # Add overamp metrics
    for key, value in metrics.items():
        if key.startswith('overamp_'):
            label = key.replace('overamp_', '').replace('_', ' ').title()
            html += f"""            <div class="metric">
                <div class="metric-label">{label}</div>
                <div class="metric-value">{value}</div>
            </div>
"""
    
    html += """        </div>
        
        <h2>🔖 UMI Metrics</h2>
        <div class="metric-section">
"""
    
    # Add UMI metrics
    if metrics.get('unique_umi_count'):
        html += f"""            <div class="metric">
                <div class="metric-label">Unique UMI Count</div>
                <div class="metric-value">{metrics.get('unique_umi_count', 'N/A'):,}</div>
            </div>
"""
    
    html += """        </div>
        
        <h2>📐 Fragmentomics</h2>
        <div class="metric-section">
"""
    
    # Add fragmentomics metrics
    for key, value in metrics.items():
        if key.startswith('frag_'):
            label = key.replace('frag_', '').replace('_', ' ').title()
            if isinstance(value, float):
                html += f"""            <div class="metric">
                <div class="metric-label">{label}</div>
                <div class="metric-value">{value:.2f}</div>
            </div>
"""
            else:
                html += f"""            <div class="metric">
                <div class="metric-label">{label}</div>
                <div class="metric-value">{value}</div>
            </div>
"""
    
    html += """        </div>
        
        <h2>🧫 Tissue of Origin</h2>
        <div class="metric-section">
"""
    
    # Add tissue of origin
    for key, value in metrics.items():
        if key.startswith('tof_tissue_'):
            idx = key.replace('tof_tissue_', '')
            fraction = metrics.get(f'tof_fraction_{idx}', 0)
            html += f"""            <div class="metric">
                <div class="metric-label">Top Tissue {idx}</div>
                <div class="metric-value">{value} ({fraction:.1%})</div>
            </div>
"""
    
    html += f"""        </div>
        
        <div class="footer">
            <p>Generated by MEDIPIPE cfMeDIP-seq Analysis Pipeline</p>
            <p>Report version: 1.0 | {datetime.now().strftime('%Y-%m-%d')}</p>
        </div>
    </div>
</body>
</html>
"""
    
    with open(output_file, 'w') as f:
        f.write(html)


def main():
    parser = argparse.ArgumentParser(
        description="Generate MEDIPIPE sample QC report"
    )
    parser.add_argument("--sample-id", required=True, help="Sample ID")
    parser.add_argument("--output-dir", required=True, help="MEDIPIPE output directory")
    parser.add_argument("--report-dir", required=True, help="Output directory for reports")
    parser.add_argument("--config", help="Path to config.yaml")
    parser.add_argument("--format", choices=['html', 'tsv', 'both'], default='both',
                        help="Output format")
    args = parser.parse_args()
    
    # Load config if provided
    config = {}
    if args.config and os.path.exists(args.config):
        import yaml
        with open(args.config) as f:
            config = yaml.safe_load(f)
    
    # Collect metrics
    print(f"Collecting QC metrics for sample: {args.sample_id}")
    metrics = collect_sample_metrics(args.sample_id, args.output_dir, config)
    
    # Create output directory
    os.makedirs(args.report_dir, exist_ok=True)
    
    # Generate reports
    if args.format in ['tsv', 'both']:
        tsv_file = os.path.join(args.report_dir, f"{args.sample_id}_qc_report.tsv")
        pd.DataFrame([metrics]).to_csv(tsv_file, sep='\t', index=False)
        print(f"  TSV report: {tsv_file}")
    
    if args.format in ['html', 'both']:
        html_file = os.path.join(args.report_dir, f"{args.sample_id}_qc_report.html")
        generate_html_report(metrics, html_file)
        print(f"  HTML report: {html_file}")
    
    print("Done!")


if __name__ == "__main__":
    main()
