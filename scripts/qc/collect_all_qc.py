#!/usr/bin/env python3
"""
Collect all QC data for MEDIPIPE samples:
- UMI fragment counts and deduplication stats
- Spike-in analysis 
- Tissue of Origin (ToO) deconvolution
"""

import os
import sys
import json
import glob
import pandas as pd
from pathlib import Path

def collect_umi_dedup_stats(output_dir):
    """Collect UMI deduplication statistics from UMI-tools logs"""
    stats = {}
    
    # Check for UMI dedup logs
    log_pattern = os.path.join(output_dir, "../work/logs/umi_dedup/*.log")
    # Also check in work directory
    log_pattern2 = os.path.join(output_dir, "../../medipipe_warp/work/logs/umi_dedup/*.log")
    
    for pattern in [log_pattern, log_pattern2]:
        for log_file in glob.glob(pattern):
            sample = os.path.basename(log_file).replace(".log", "").replace("_dedup", "")
            try:
                with open(log_file) as f:
                    content = f.read()
                    # Parse UMI-tools output
                    for line in content.split('\n'):
                        if 'Input Reads:' in line:
                            stats.setdefault(sample, {})['input_reads'] = int(line.split(':')[1].strip())
                        elif 'Reads Output:' in line or 'Number of reads out:' in line:
                            stats.setdefault(sample, {})['output_reads'] = int(line.split(':')[1].strip())
                        elif 'unique UMIs' in line.lower():
                            # Try to parse unique UMI count
                            parts = line.split()
                            for i, p in enumerate(parts):
                                if p.isdigit():
                                    stats.setdefault(sample, {})['unique_umis'] = int(p)
                                    break
            except Exception as e:
                print(f"Warning: Could not parse {log_file}: {e}")
    
    # Also check BAM stats for dedup counts
    stats_pattern = os.path.join(output_dir, "dedup_bam_umi/*_dedup.bam.stats.txt")
    stats_pattern2 = os.path.join(output_dir, "dedup_bam_pe/*_dedup.bam.stats.txt")
    
    for pattern in [stats_pattern, stats_pattern2]:
        for stats_file in glob.glob(pattern):
            sample = os.path.basename(stats_file).replace("_dedup.bam.stats.txt", "")
            try:
                with open(stats_file) as f:
                    for line in f:
                        if line.startswith('SN\tsequences:'):
                            stats.setdefault(sample, {})['dedup_reads'] = int(line.split('\t')[2].strip())
                        elif line.startswith('SN\treads mapped:'):
                            stats.setdefault(sample, {})['mapped_reads'] = int(line.split('\t')[2].strip())
            except Exception as e:
                print(f"Warning: Could not parse {stats_file}: {e}")
    
    return stats

def collect_spikein_stats(output_dir):
    """Collect spike-in analysis statistics"""
    stats = {}
    
    # Check for spike-in MEDIPS QC output
    spikein_dir = os.path.join(output_dir, "spikein_quant")
    
    if os.path.exists(spikein_dir):
        for qc_file in glob.glob(os.path.join(spikein_dir, "*_spikein_qc.txt")):
            sample = os.path.basename(qc_file).replace("_spikein_qc.txt", "")
            try:
                with open(qc_file) as f:
                    for line in f:
                        line = line.strip()
                        if '=' in line:
                            key, value = line.split('=', 1)
                            key = key.strip().lower().replace(' ', '_')
                            try:
                                value = float(value.strip())
                            except:
                                value = value.strip()
                            stats.setdefault(sample, {})[f'spikein_{key}'] = value
            except Exception as e:
                print(f"Warning: Could not parse {qc_file}: {e}")
    
    # Also check for spike-in counts in meth_quant
    for count_file in glob.glob(os.path.join(output_dir, "meth_quant/*_spikein_count.txt")):
        sample = os.path.basename(count_file).replace("_spikein_count.txt", "")
        try:
            df = pd.read_csv(count_file, sep='\t')
            if 'meth_reads' in df.columns:
                stats.setdefault(sample, {})['spikein_meth_reads'] = df['meth_reads'].sum()
            if 'unmeth_reads' in df.columns:
                stats.setdefault(sample, {})['spikein_unmeth_reads'] = df['unmeth_reads'].sum()
        except Exception as e:
            pass
    
    return stats

def collect_too_data(output_dir):
    """Collect Tissue of Origin deconvolution data"""
    stats = {}
    
    # Check for NNLS deconvolution output
    too_file = os.path.join(output_dir, "tissue_deconv/nnls_proportions.tsv")
    
    if os.path.exists(too_file):
        try:
            df = pd.read_csv(too_file, sep='\t', index_col=0)
            
            # Define tissue categories
            brain_tissues = ['brain', 'neuron', 'astrocyte', 'oligodendrocyte', 'cerebellum', 'cortex']
            blood_tissues = ['blood', 'monocyte', 'neutrophil', 'lymphocyte', 'cd4', 'cd8', 'b_cell', 'nk_cell']
            liver_tissues = ['liver', 'hepatocyte']
            
            for sample in df.columns:
                sample_clean = sample.replace('.bw', '').replace('_S1_L003', '').replace('_S', '_')
                
                stats[sample_clean] = {}
                
                # Calculate aggregates
                brain_total = 0
                blood_total = 0
                liver_total = 0
                
                for tissue, proportion in df[sample].items():
                    tissue_lower = tissue.lower()
                    if any(bt in tissue_lower for bt in brain_tissues):
                        brain_total += proportion
                    elif any(bt in tissue_lower for bt in blood_tissues):
                        blood_total += proportion
                    elif any(lt in tissue_lower for lt in liver_tissues):
                        liver_total += proportion
                
                stats[sample_clean]['too_brain'] = round(brain_total * 100, 1)
                stats[sample_clean]['too_blood'] = round(blood_total * 100, 1)
                stats[sample_clean]['too_liver'] = round(liver_total * 100, 1)
                stats[sample_clean]['too_other'] = round((1 - brain_total - blood_total - liver_total) * 100, 1)
                
                # Get top tissues
                top_tissues = df[sample].nlargest(5)
                stats[sample_clean]['too_top_tissues'] = ', '.join([f"{t}:{v:.1%}" for t, v in top_tissues.items()])
                
        except Exception as e:
            print(f"Warning: Could not parse {too_file}: {e}")
    
    return stats

def collect_meth_qc_stats(output_dir):
    """Collect methylation QC statistics from MEDIPS"""
    stats = {}
    
    for summary_file in glob.glob(os.path.join(output_dir, "meth_quant/*.summary.txt")):
        sample = os.path.basename(summary_file).replace(".summary.txt", "")
        try:
            with open(summary_file) as f:
                for line in f:
                    if ':' in line:
                        key, value = line.strip().split(':', 1)
                        key = key.strip().lower().replace(' ', '_')
                        try:
                            value = float(value.strip())
                        except:
                            value = value.strip()
                        stats.setdefault(sample, {})[key] = value
        except Exception as e:
            print(f"Warning: Could not parse {summary_file}: {e}")
    
    return stats

def main():
    output_dir = "/data/medipipe_data/output"
    
    print("=" * 60)
    print("MEDIPIPE QC Data Collection")
    print("=" * 60)
    
    # Collect all stats
    umi_stats = collect_umi_dedup_stats(output_dir)
    spikein_stats = collect_spikein_stats(output_dir)
    too_stats = collect_too_data(output_dir)
    meth_stats = collect_meth_qc_stats(output_dir)
    
    # Merge all stats
    all_samples = set(umi_stats.keys()) | set(spikein_stats.keys()) | set(too_stats.keys()) | set(meth_stats.keys())
    
    print(f"\nFound {len(all_samples)} samples with QC data\n")
    
    # Create summary table
    rows = []
    for sample in sorted(all_samples):
        row = {'sample': sample}
        row.update(umi_stats.get(sample, {}))
        row.update(spikein_stats.get(sample, {}))
        row.update(too_stats.get(sample, {}))
        row.update(meth_stats.get(sample, {}))
        rows.append(row)
    
    df = pd.DataFrame(rows)
    
    # Print summary
    print("\n=== UMI Deduplication Statistics ===")
    dedup_cols = ['sample', 'input_reads', 'output_reads', 'dedup_reads', 'unique_umis']
    available_dedup = [c for c in dedup_cols if c in df.columns]
    if len(available_dedup) > 1:
        print(df[available_dedup].to_string(index=False))
    else:
        print("No UMI deduplication data found")
    
    print("\n=== Spike-in Statistics ===")
    spikein_cols = [c for c in df.columns if 'spikein' in c.lower()]
    if spikein_cols:
        print(df[['sample'] + spikein_cols].to_string(index=False))
    else:
        print("No spike-in data found")
    
    print("\n=== Tissue of Origin (ToO) ===")
    too_cols = [c for c in df.columns if c.startswith('too_')]
    if too_cols:
        print(df[['sample'] + too_cols].to_string(index=False))
    else:
        print("No ToO data found")
    
    print("\n=== Methylation QC ===")
    meth_cols = ['sample', 'total_reads', 'total_windows', 'genome-wide_coverage']
    available_meth = [c for c in meth_cols if c in df.columns]
    if len(available_meth) > 1:
        print(df[available_meth].to_string(index=False))
    else:
        print("No methylation QC data found")
    
    # Save full report
    report_path = os.path.join(output_dir, "qc_reports/all_qc_summary.tsv")
    os.makedirs(os.path.dirname(report_path), exist_ok=True)
    df.to_csv(report_path, sep='\t', index=False)
    print(f"\n\nFull report saved to: {report_path}")
    
    return df

if __name__ == "__main__":
    main()
