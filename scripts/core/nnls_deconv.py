#!/usr/bin/env python3
"""
NNLS Deconvolution for Tissue of Origin Analysis
Estimates tissue contributions to cfDNA methylation signal.

Based on Loyfer et al. 2023 (Nature) methodology for cfDNA tissue deconvolution.
Reference: https://doi.org/10.1038/s41586-022-05580-6
"""

import argparse
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from sklearn.preprocessing import MinMaxScaler, StandardScaler
import warnings
import os
warnings.filterwarnings('ignore')


def load_marker_bed(bed_path: str) -> pd.DataFrame:
    """
    Load marker regions BED file to create coordinate-to-marker_id mapping.
    
    BED format: chr, start, end, marker_id
    """
    df = pd.read_csv(bed_path, sep='\t', header=None, 
                     names=['chr', 'start', 'end', 'marker_id'])
    # Create coordinate key for matching
    df['coord_key'] = df['chr'] + ':' + df['start'].astype(str) + '-' + df['end'].astype(str)
    print(f"Loaded {len(df)} marker regions from BED file")
    return df


def load_reference_matrix(filepath: str) -> pd.DataFrame:
    """Load reference tissue methylation matrix (Loyfer et al. 2023 format)."""
    df = pd.read_csv(filepath, sep='\t', index_col=0)
    print(f"Loaded reference matrix: {df.shape[0]} markers x {df.shape[1]} tissues")
    return df


def load_sample_matrix_with_coords(filepath: str, marker_bed: pd.DataFrame = None) -> pd.DataFrame:
    """
    Load sample methylation matrix, handling coordinate-based format from multiBigwigSummary.
    
    If marker_bed is provided, converts coordinates to marker IDs.
    """
    df = pd.read_csv(filepath, sep='\t', comment=None)
    
    # Clean column names (remove quotes from multiBigwigSummary output)
    df.columns = [c.strip("'") for c in df.columns]
    
    # Check if first column is a comment line indicator
    if df.columns[0].startswith('#'):
        df.columns = [c.lstrip('#').strip("'") for c in df.columns]
    
    # Detect coordinate-based format (chr, start, end as first 3 columns)
    first_cols = [c.lower() for c in df.columns[:3]]
    is_coord_format = ('chr' in first_cols and 'start' in first_cols and 'end' in first_cols)
    
    if is_coord_format:
        print("Detected coordinate-based sample matrix format")
        
        # Standardize column names
        col_map = {}
        for c in df.columns[:3]:
            if c.lower() == 'chr':
                col_map[c] = 'chr'
            elif c.lower() == 'start':
                col_map[c] = 'start'
            elif c.lower() == 'end':
                col_map[c] = 'end'
        df = df.rename(columns=col_map)
        
        # Create coordinate key
        df['coord_key'] = df['chr'] + ':' + df['start'].astype(str) + '-' + df['end'].astype(str)
        
        if marker_bed is not None:
            # Merge with marker BED to get marker IDs
            # Also try matching without exact coordinates (within tolerance)
            merged = df.merge(marker_bed[['coord_key', 'marker_id']], on='coord_key', how='left')
            
            # For unmatched, try matching by chr and overlapping coordinates
            unmatched_mask = merged['marker_id'].isna()
            n_unmatched = unmatched_mask.sum()
            
            if n_unmatched > 0:
                print(f"  {n_unmatched} regions without exact coordinate match, trying overlap matching...")
                # Try to match by overlap
                for idx in merged[unmatched_mask].index:
                    row = merged.loc[idx]
                    chrom, start, end = row['chr'], row['start'], row['end']
                    # Find markers on same chromosome with overlapping coordinates
                    candidates = marker_bed[
                        (marker_bed['chr'] == chrom) & 
                        (marker_bed['start'] <= end) & 
                        (marker_bed['end'] >= start)
                    ]
                    if len(candidates) == 1:
                        merged.loc[idx, 'marker_id'] = candidates.iloc[0]['marker_id']
                    elif len(candidates) > 1:
                        # Pick best overlap
                        overlaps = candidates.apply(
                            lambda r: min(end, r['end']) - max(start, r['start']), axis=1
                        )
                        best_idx = overlaps.idxmax()
                        merged.loc[idx, 'marker_id'] = candidates.loc[best_idx, 'marker_id']
            
            # Set marker_id as index and drop coordinate columns
            matched = merged.dropna(subset=['marker_id'])
            print(f"  Matched {len(matched)} of {len(df)} regions to marker IDs")
            
            # Handle duplicates by averaging values for the same marker
            # Drop coordinate columns first
            for col in ['chr', 'start', 'end', 'coord_key']:
                if col in matched.columns:
                    matched = matched.drop(columns=[col])
            
            # Check for duplicates
            n_duplicates = matched.duplicated(subset=['marker_id']).sum()
            if n_duplicates > 0:
                print(f"  Found {n_duplicates} duplicate markers, averaging values...")
                # Group by marker_id and take mean of numeric columns
                numeric_cols = matched.select_dtypes(include=[np.number]).columns.tolist()
                matched = matched.groupby('marker_id')[numeric_cols].mean()
            else:
                matched = matched.set_index('marker_id')
            
            return matched
        else:
            # No marker BED provided, use coordinates as index
            df = df.set_index('coord_key')
            df = df.drop(columns=['chr', 'start', 'end'], errors='ignore')
            return df
    else:
        # Standard format with marker IDs
        return df


def load_sample_methylation(filepath: str, marker_col: str = None, marker_bed: pd.DataFrame = None) -> pd.Series:
    """Load sample methylation values (single sample)."""
    df = pd.read_csv(filepath, sep='\t')
    
    # Clean column names
    df.columns = [c.strip("'").lstrip('#') for c in df.columns]
    
    # Try to identify the marker and value columns
    if 'region' in df.columns:
        df = df.set_index('region')
    elif 'marker' in df.columns:
        df = df.set_index('marker')
    elif 'marker_id' in df.columns:
        df = df.set_index('marker_id')
    elif marker_col:
        df = df.set_index(marker_col)
    else:
        # Assume first column is marker, last is value
        df = df.set_index(df.columns[0])
    
    # Get the methylation values (try common column names)
    for col in ['methylation', 'meth', 'beta', 'rms', 'value']:
        if col in df.columns:
            return df[col]
    
    # Use last numeric column
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    if len(numeric_cols) > 0:
        return df[numeric_cols[-1]]
    
    raise ValueError("Could not identify methylation values in sample file")


def scale_data(reference: pd.DataFrame, sample: pd.Series, method: str = 'minmax_by_marker') -> tuple:
    """Scale reference and sample data."""
    
    # Align markers
    common_markers = reference.index.intersection(sample.index)
    print(f"Common markers: {len(common_markers)}")
    
    if len(common_markers) < 10:
        raise ValueError(f"Too few common markers: {len(common_markers)}")
    
    ref_aligned = reference.loc[common_markers]
    sample_aligned = sample.loc[common_markers]
    
    if method == 'minmax_by_marker':
        # Scale each marker independently to [0, 1]
        ref_min = ref_aligned.min(axis=1)
        ref_max = ref_aligned.max(axis=1)
        ref_range = ref_max - ref_min
        ref_range[ref_range == 0] = 1  # Avoid division by zero
        
        ref_scaled = ref_aligned.sub(ref_min, axis=0).div(ref_range, axis=0)
        sample_scaled = (sample_aligned - ref_min) / ref_range
        
    elif method == 'zscore':
        # Z-score normalization
        ref_mean = ref_aligned.mean(axis=1)
        ref_std = ref_aligned.std(axis=1)
        ref_std[ref_std == 0] = 1
        
        ref_scaled = ref_aligned.sub(ref_mean, axis=0).div(ref_std, axis=0)
        sample_scaled = (sample_aligned - ref_mean) / ref_std
        
    elif method == 'none':
        ref_scaled = ref_aligned
        sample_scaled = sample_aligned
        
    else:
        raise ValueError(f"Unknown scaling method: {method}")
    
    return ref_scaled, sample_scaled


def run_nnls(reference: pd.DataFrame, sample: pd.Series) -> dict:
    """
    Run NNLS deconvolution.
    
    Args:
        reference: Reference matrix (markers x tissues)
        sample: Sample methylation values
    
    Returns:
        Dictionary with tissue fractions and residual
    """
    # Convert to numpy
    A = reference.values
    b = sample.values
    
    # Handle any NaN values
    valid_mask = ~(np.isnan(b) | np.any(np.isnan(A), axis=1))
    A = A[valid_mask]
    b = b[valid_mask]
    
    print(f"Running NNLS on {A.shape[0]} markers, {A.shape[1]} tissues")
    
    # Run NNLS
    x, residual = nnls(A, b)
    
    # Normalize to sum to 1
    x_sum = x.sum()
    if x_sum > 0:
        x_normalized = x / x_sum
    else:
        x_normalized = x
    
    # Calculate fit statistics
    predicted = A @ x
    ss_res = np.sum((b - predicted) ** 2)
    ss_tot = np.sum((b - np.mean(b)) ** 2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    
    result = {
        'fractions': dict(zip(reference.columns, x_normalized)),
        'raw_coefficients': dict(zip(reference.columns, x)),
        'residual': residual,
        'r_squared': r_squared,
        'n_markers': len(b)
    }
    
    return result


def run_bootstrap(reference: pd.DataFrame, sample: pd.Series, n_iterations: int = 100) -> dict:
    """Run bootstrap analysis to estimate confidence intervals."""
    
    n_markers = len(sample)
    bootstrap_results = {tissue: [] for tissue in reference.columns}
    
    for i in range(n_iterations):
        # Sample with replacement
        idx = np.random.choice(n_markers, size=n_markers, replace=True)
        ref_boot = reference.iloc[idx]
        sample_boot = sample.iloc[idx]
        
        # Run NNLS
        result = run_nnls(ref_boot, sample_boot)
        
        for tissue, frac in result['fractions'].items():
            bootstrap_results[tissue].append(frac)
    
    # Calculate statistics
    stats = {}
    for tissue, fracs in bootstrap_results.items():
        fracs = np.array(fracs)
        stats[tissue] = {
            'mean': np.mean(fracs),
            'std': np.std(fracs),
            'ci_lower': np.percentile(fracs, 2.5),
            'ci_upper': np.percentile(fracs, 97.5)
        }
    
    return stats


def scale_cfmedip_coverage(sample_matrix: pd.DataFrame) -> pd.DataFrame:
    """
    Scale cfMeDIP-seq coverage data using per-marker min-max scaling.
    
    This converts coverage values to [0,1] range per marker, making them
    comparable to reference methylation beta values.
    
    Based on Loyfer et al. 2023 methodology.
    """
    # Per-marker min-max scaling
    mins = sample_matrix.min(axis=1)
    maxs = sample_matrix.max(axis=1)
    rng = maxs - mins
    rng[rng < 1e-6] = 1e-6  # Avoid division by zero
    
    scaled = sample_matrix.sub(mins, axis=0).div(rng, axis=0)
    
    # Clip to [0,1] and fill NaN with 0 (neutral value)
    scaled = scaled.clip(lower=0, upper=1)
    scaled = scaled.fillna(0)
    
    return scaled


def scale_data_multi(reference: pd.DataFrame, sample_matrix: pd.DataFrame, method: str = 'minmax_by_marker') -> tuple:
    """
    Scale reference and multi-sample data.
    
    Args:
        reference: Reference matrix (markers x tissues)
        sample_matrix: Sample matrix (markers x samples)
        method: Scaling method
    
    Returns:
        Tuple of (scaled_reference, scaled_samples)
    """
    # Align markers
    common_markers = reference.index.intersection(sample_matrix.index)
    print(f"Common markers: {len(common_markers)}")
    
    if len(common_markers) < 10:
        raise ValueError(f"Too few common markers: {len(common_markers)}")
    
    ref_aligned = reference.loc[common_markers]
    samples_aligned = sample_matrix.loc[common_markers]
    
    if method == 'cfmedip':
        # cfMeDIP-seq specific: scale coverage to [0,1] per marker
        # Reference is already beta values, so no scaling needed
        print("Using cfMeDIP-specific scaling (per-marker min-max on samples)")
        ref_scaled = ref_aligned
        samples_scaled = scale_cfmedip_coverage(samples_aligned)
        
    elif method == 'minmax_by_marker':
        # Scale each marker independently to [0, 1] based on reference range
        ref_min = ref_aligned.min(axis=1)
        ref_max = ref_aligned.max(axis=1)
        ref_range = ref_max - ref_min
        ref_range[ref_range == 0] = 1  # Avoid division by zero
        
        ref_scaled = ref_aligned.sub(ref_min, axis=0).div(ref_range, axis=0)
        samples_scaled = samples_aligned.sub(ref_min, axis=0).div(ref_range, axis=0)
        
    elif method == 'zscore':
        # Z-score normalization
        ref_mean = ref_aligned.mean(axis=1)
        ref_std = ref_aligned.std(axis=1)
        ref_std[ref_std == 0] = 1
        
        ref_scaled = ref_aligned.sub(ref_mean, axis=0).div(ref_std, axis=0)
        samples_scaled = samples_aligned.sub(ref_mean, axis=0).div(ref_std, axis=0)
        
    elif method == 'none':
        ref_scaled = ref_aligned
        samples_scaled = samples_aligned
        
    else:
        raise ValueError(f"Unknown scaling method: {method}")
    
    return ref_scaled, samples_scaled


def run_multi_sample_deconv(reference: pd.DataFrame, sample_matrix: pd.DataFrame, 
                            scale_method: str = 'minmax_by_marker') -> pd.DataFrame:
    """
    Run NNLS deconvolution on multiple samples.
    
    Args:
        reference: Reference matrix (markers x tissues)
        sample_matrix: Sample matrix (markers x samples)  
        scale_method: Scaling method
    
    Returns:
        DataFrame with tissue fractions for each sample
    """
    # Scale data
    ref_scaled, samples_scaled = scale_data_multi(reference, sample_matrix, scale_method)
    
    results = []
    for sample_name in samples_scaled.columns:
        sample_values = samples_scaled[sample_name]
        
        try:
            result = run_nnls(ref_scaled, sample_values)
            row = {'sample': sample_name}
            row.update(result['fractions'])
            row['r_squared'] = result['r_squared']
            row['n_markers'] = result['n_markers']
            results.append(row)
        except Exception as e:
            print(f"  Warning: Failed for {sample_name}: {e}")
    
    return pd.DataFrame(results)


def main():
    parser = argparse.ArgumentParser(description='NNLS Deconvolution for Tissue of Origin (Loyfer et al. 2023)')
    parser.add_argument('-s', '--sample', required=True, 
                       help='Sample methylation file (TSV) - can be single sample or multi-sample matrix')
    parser.add_argument('-r', '--reference', required=True, 
                       help='Reference matrix file (TSV) with marker_id index')
    parser.add_argument('-o', '--output', required=True, help='Output file path')
    parser.add_argument('-m', '--markers_bed', default=None,
                       help='Marker regions BED file (chr, start, end, marker_id) for coordinate-to-ID conversion')
    parser.add_argument('--scale', default='cfmedip', 
                       choices=['cfmedip', 'minmax_by_marker', 'zscore', 'none'],
                       help='Scaling method: cfmedip (per-marker min-max on samples for cfMeDIP-seq), '
                            'minmax_by_marker (scale both by reference range), zscore, none (default: cfmedip)')
    parser.add_argument('--bootstrap', type=int, default=0,
                       help='Number of bootstrap iterations (0 to disable)')
    parser.add_argument('--marker_col', default=None, help='Column name for markers in sample file')
    
    args = parser.parse_args()
    
    print("=== NNLS Tissue of Origin Deconvolution ===")
    print(f"Based on Loyfer et al. 2023 (Nature) methodology")
    print(f"Sample: {args.sample}")
    print(f"Reference: {args.reference}")
    print(f"Markers BED: {args.markers_bed}")
    print(f"Scaling: {args.scale}")
    
    # Load marker BED if provided (for coordinate-to-marker_id conversion)
    marker_bed = None
    if args.markers_bed:
        marker_bed = load_marker_bed(args.markers_bed)
    
    # Load reference matrix
    reference = load_reference_matrix(args.reference)
    
    # Try loading as multi-sample matrix first
    try:
        sample_matrix = load_sample_matrix_with_coords(args.sample, marker_bed)
        
        # Check if it's a multi-sample matrix (more than 1 column after dropping coords)
        sample_cols = [c for c in sample_matrix.columns if c not in ['chr', 'start', 'end', 'coord_key']]
        
        if len(sample_cols) > 1:
            print(f"\nProcessing {len(sample_cols)} samples...")
            
            # Run multi-sample deconvolution
            results_df = run_multi_sample_deconv(reference, sample_matrix, args.scale)
            
            # Save results
            results_df.to_csv(args.output, sep='\t', index=False)
            print(f"\nResults saved to: {args.output}")
            
            # Print summary
            print(f"\nProcessed {len(results_df)} samples successfully")
            print(f"Tissue columns: {[c for c in results_df.columns if c not in ['sample', 'r_squared', 'n_markers']]}")
            
            # Also save wide format for backward compatibility
            wide_output = args.output.replace('.tsv', '_wide.tsv')
            if wide_output == args.output:
                wide_output = args.output + '.wide.tsv'
            results_df.to_csv(wide_output, sep='\t', index=False)
            print(f"Wide format saved to: {wide_output}")
            
            print("\nDeconvolution complete!")
            return
            
    except Exception as e:
        print(f"Note: Could not load as multi-sample matrix ({e}), trying single-sample format...")
    
    # Fall back to single-sample mode
    sample = load_sample_methylation(args.sample, args.marker_col, marker_bed)
    
    # Scale data
    ref_scaled, sample_scaled = scale_data(reference, sample, args.scale)
    
    # Run NNLS
    result = run_nnls(ref_scaled, sample_scaled)
    
    print(f"\nR-squared: {result['r_squared']:.4f}")
    print(f"Residual: {result['residual']:.4f}")
    print(f"\nTissue Fractions:")
    
    # Sort by fraction
    fractions_sorted = sorted(result['fractions'].items(), key=lambda x: x[1], reverse=True)
    for tissue, frac in fractions_sorted:
        print(f"  {tissue}: {frac:.4f} ({frac*100:.2f}%)")
    
    # Save results
    results_df = pd.DataFrame([
        {'tissue': tissue, 'fraction': frac, 'percentage': frac * 100}
        for tissue, frac in fractions_sorted
    ])
    
    results_df.to_csv(args.output, sep='\t', index=False)
    print(f"\nResults saved to: {args.output}")
    
    # Save summary
    summary = {
        'sample': args.sample,
        'n_markers': result['n_markers'],
        'r_squared': result['r_squared'],
        'residual': result['residual'],
        'scaling': args.scale
    }
    summary.update({f"frac_{k}": v for k, v in result['fractions'].items()})
    
    summary_df = pd.DataFrame([summary])
    summary_file = args.output.replace('.tsv', '.summary.tsv')
    if summary_file == args.output:
        summary_file = args.output + '.summary.tsv'
    summary_df.to_csv(summary_file, sep='\t', index=False)
    print(f"Summary saved to: {summary_file}")
    
    # Run bootstrap if requested
    if args.bootstrap > 0:
        print(f"\nRunning {args.bootstrap} bootstrap iterations...")
        boot_stats = run_bootstrap(ref_scaled, sample_scaled, args.bootstrap)
        
        boot_df = pd.DataFrame([
            {
                'tissue': tissue,
                'mean': stats['mean'],
                'std': stats['std'],
                'ci_lower': stats['ci_lower'],
                'ci_upper': stats['ci_upper']
            }
            for tissue, stats in boot_stats.items()
        ])
        
        boot_file = args.output.replace('.tsv', '.bootstrap.tsv')
        if boot_file == args.output:
            boot_file = args.output + '.bootstrap.tsv'
        boot_df.to_csv(boot_file, sep='\t', index=False)
        print(f"Bootstrap statistics saved to: {boot_file}")
    
    print("\nDeconvolution complete!")


if __name__ == '__main__':
    main()
