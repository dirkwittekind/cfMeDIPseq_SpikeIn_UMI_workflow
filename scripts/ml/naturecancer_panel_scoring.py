#!/usr/bin/env python3
"""
NatureCancer Panel Scoring Analysis

Calculates z-scores for AEG samples based on NatureCancer reference DMR panels.
This computes the mean methylation within each panel's regions and compares
AEG samples to controls using z-score normalization.

Based on Zeng et al. Nature Cancer methodology:
- |z| > 2 is treated as a descriptive deviation flag
"""

import os
import subprocess
import pandas as pd
import numpy as np
import tempfile
from scipy import stats

# Paths
NATURE_CANCER_DIR = "/home/dirk/medipipe_warp/resources/Code/15191455/pan_cancer_cfMeDIP_Figures_source_data/pan_cancer_cfMDIP_Figures_source_data/DMRs_Source_Data/BED_files"
MATRIX_FILE = "/data/medipipe_data/output/ml_discrimination/matrices/windows/hg38_w2000/all/matrix_autosomal.tsv"
ANNOTATION_FILE = "/data/medipipe_data/output/_upload_package_meta/sample_annotation.active.with_covariates.clean.tsv"

def load_matrix_with_coords():
    """Load the methylation matrix with coordinates."""
    print("Loading methylation matrix...")
    df = pd.read_csv(MATRIX_FILE, sep="\t")
    print(f"  Loaded {len(df)} regions, {len(df.columns)-3} samples")
    return df

def load_annotation():
    """Load sample annotation."""
    ann = pd.read_csv(ANNOTATION_FILE, sep="\t")
    return ann

def get_panel_regions(panel_bed_path):
    """Load regions from a NatureCancer panel BED file."""
    regions = []
    with open(panel_bed_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                regions.append({
                    'chr': parts[0],
                    'start': int(parts[1]),
                    'end': int(parts[2])
                })
    return pd.DataFrame(regions)

def find_overlapping_matrix_regions(matrix_df, panel_regions):
    """
    Find matrix regions that overlap with panel regions.
    Uses bedtools intersect for efficiency.
    """
    # Create temp BED files
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        matrix_bed = f.name
        for _, row in matrix_df.iterrows():
            f.write(f"{row['chr']}\t{row['start']}\t{row['end']}\n")
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        panel_bed = f.name
        for _, row in panel_regions.iterrows():
            f.write(f"{row['chr']}\t{row['start']}\t{row['end']}\n")
    
    try:
        # Find overlapping regions
        result = subprocess.run(
            ["bedtools", "intersect", "-a", matrix_bed, "-b", panel_bed, "-wa"],
            capture_output=True, text=True, check=True
        )
        
        # Parse overlapping regions
        overlapping = set()
        for line in result.stdout.strip().split('\n'):
            if line:
                parts = line.split('\t')
                overlapping.add(f"{parts[0]}:{parts[1]}-{parts[2]}")
        
        return overlapping
    finally:
        os.unlink(matrix_bed)
        os.unlink(panel_bed)

def calculate_panel_scores(matrix_df, sample_cols, panel_regions):
    """
    Calculate mean methylation score for each sample within panel regions.
    """
    # Create region ID for matching
    matrix_df['region_id'] = matrix_df['chr'] + ':' + matrix_df['start'].astype(str) + '-' + matrix_df['end'].astype(str)
    
    # Find overlapping regions
    overlapping = find_overlapping_matrix_regions(matrix_df, panel_regions)
    
    if len(overlapping) == 0:
        return None, 0
    
    # Filter to overlapping regions
    mask = matrix_df['region_id'].isin(overlapping)
    panel_data = matrix_df.loc[mask, sample_cols]
    
    # Calculate mean methylation per sample
    scores = panel_data.mean(axis=0)
    
    return scores, len(overlapping)

def main():
    print("=" * 80)
    print("NATURECANCER PANEL SCORING ANALYSIS")
    print("Z-score calculation for AEG vs CTRL")
    print("=" * 80)
    
    # Load data
    matrix_df = load_matrix_with_coords()
    ann = load_annotation()
    
    # Get sample columns (match matrix columns to annotation)
    sample_cols = []
    sample_groups = {}
    
    matrix_samples = [c for c in matrix_df.columns if c not in ['chr', 'start', 'end']]
    
    for m_sample in matrix_samples:
        for _, row in ann.iterrows():
            a_sample = row['sample_id']
            if a_sample.startswith(m_sample + '_') or a_sample == m_sample:
                sample_cols.append(m_sample)
                sample_groups[m_sample] = row['group']
                break
    
    print(f"\nMatched {len(sample_cols)} samples")
    
    aeg_samples = [s for s in sample_cols if sample_groups.get(s) == 'AEG']
    ctrl_samples = [s for s in sample_cols if sample_groups.get(s) == 'CTRL']
    
    print(f"  AEG: {len(aeg_samples)} samples")
    print(f"  CTRL: {len(ctrl_samples)} samples")
    
    # Get all panel BED files
    panel_files = {}
    for f in sorted(os.listdir(NATURE_CANCER_DIR)):
        if f.endswith('.bed') and not f.startswith('.'):
            name = f.replace('.bed', '')
            # Skip very large files (PBL has 2M regions)
            if 'PBL' not in name:
                panel_files[name] = os.path.join(NATURE_CANCER_DIR, f)
    
    print(f"\nAnalyzing {len(panel_files)} NatureCancer panels...")
    
    # Calculate scores for each panel
    results = []
    
    for panel_name, panel_path in sorted(panel_files.items()):
        panel_regions = get_panel_regions(panel_path)
        
        if len(panel_regions) == 0:
            continue
        
        scores, n_overlap = calculate_panel_scores(matrix_df, sample_cols, panel_regions)
        
        if scores is None or n_overlap < 5:
            continue
        
        # Calculate z-score: (AEG_mean - CTRL_mean) / CTRL_std
        aeg_scores = scores[aeg_samples].values
        ctrl_scores = scores[ctrl_samples].values
        
        aeg_mean = np.mean(aeg_scores)
        ctrl_mean = np.mean(ctrl_scores)
        ctrl_std = np.std(ctrl_scores, ddof=1)
        
        if ctrl_std > 0:
            z_score = (aeg_mean - ctrl_mean) / ctrl_std
        else:
            z_score = 0
        
        # Also calculate t-test p-value
        if len(aeg_scores) > 1 and len(ctrl_scores) > 1:
            t_stat, p_value = stats.ttest_ind(aeg_scores, ctrl_scores)
        else:
            t_stat, p_value = 0, 1
        
        results.append({
            'Panel': panel_name,
            'N_Regions': len(panel_regions),
            'N_Overlap': n_overlap,
            'AEG_Mean': round(aeg_mean, 6),
            'CTRL_Mean': round(ctrl_mean, 6),
            'CTRL_Std': round(ctrl_std, 6),
            'Z_Score': round(z_score, 2),
            'Direction': 'Higher in AEG' if z_score > 0 else 'Lower in AEG',
            'Significant': '**' if abs(z_score) >= 2 else ('*' if abs(z_score) >= 1.5 else ''),
            'P_Value': round(p_value, 4)
        })
    
    # Create results DataFrame
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('Z_Score', key=abs, ascending=False)
    
    print("\n" + "=" * 80)
    print("RESULTS: Panel Z-Scores (AEG vs CTRL)")
    print("=" * 80)
    print("Threshold: |z| > 2 = significant deviation (per Zeng et al.)")
    print("           |z| > 1.5 = trending\n")
    
    # Show key cancer panels
    key_panels = ['Colorectal', 'panCancer', 'Pancreatic', 'Brain', 'Head and Neck', 
                  'Lung', 'Breast', 'Prostate', 'AML', 'Bladder', 'Renal', 'Uveal']
    
    print("Cancer-Specific Panels:")
    print("-" * 80)
    for _, row in results_df.iterrows():
        if any(k in row['Panel'] for k in key_panels):
            sig = row['Significant']
            print(f"  {row['Panel'][:45]:<45} z={row['Z_Score']:>6.2f} {sig:<2} (n={row['N_Overlap']})")
    
    print("\n" + "-" * 80)
    print("Other Panels:")
    for _, row in results_df.iterrows():
        if not any(k in row['Panel'] for k in key_panels):
            sig = row['Significant']
            print(f"  {row['Panel'][:45]:<45} z={row['Z_Score']:>6.2f} {sig:<2} (n={row['N_Overlap']})")
    
    # Save results
    output_dir = "/data/medipipe_data/output/ml_discrimination/Outputs/AEG_vs_CTRL/naturecancer_overlap"
    os.makedirs(output_dir, exist_ok=True)
    
    output_file = os.path.join(output_dir, "panel_zscore_results.tsv")
    results_df.to_csv(output_file, sep="\t", index=False)
    print(f"\nSaved: {output_file}")
    
    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    
    sig_panels = results_df[results_df['Z_Score'].abs() >= 2]
    trend_panels = results_df[(results_df['Z_Score'].abs() >= 1.5) & (results_df['Z_Score'].abs() < 2)]
    
    print(f"Panels with |z| >= 2 (significant): {len(sig_panels)}")
    for _, row in sig_panels.iterrows():
        print(f"  - {row['Panel']}: z = {row['Z_Score']}")
    
    print(f"\nPanels with 1.5 <= |z| < 2 (trending): {len(trend_panels)}")
    for _, row in trend_panels.iterrows():
        print(f"  - {row['Panel']}: z = {row['Z_Score']}")
    
    print("\n" + "=" * 80)

if __name__ == "__main__":
    main()
