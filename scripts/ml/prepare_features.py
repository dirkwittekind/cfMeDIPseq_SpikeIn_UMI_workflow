#!/usr/bin/env python3
"""
Feature Preparation for ML Classification
Prepares methylation features from multiple samples for ML training.
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.feature_selection import VarianceThreshold, SelectKBest, f_classif
from sklearn.preprocessing import StandardScaler
import pickle


def load_methylation_matrix(input_files: list, sample_ids: list = None) -> pd.DataFrame:
    """
    Load and combine methylation data from multiple samples.
    
    Args:
        input_files: List of methylation TSV files
        sample_ids: Optional list of sample IDs
    
    Returns:
        DataFrame with samples as rows and regions as columns
    """
    all_samples = []
    
    for i, filepath in enumerate(input_files):
        df = pd.read_csv(filepath, sep='\t')
        
        # Identify region column
        region_cols = ['region', 'chr', 'window']
        region_col = None
        for col in region_cols:
            if col in df.columns:
                region_col = col
                break
        
        if region_col is None:
            # Create region from chr, start, end
            if all(c in df.columns for c in ['chr', 'start', 'stop']):
                df['region'] = df['chr'] + ':' + df['start'].astype(str) + '-' + df['stop'].astype(str)
                region_col = 'region'
            elif all(c in df.columns for c in ['chr', 'start', 'end']):
                df['region'] = df['chr'] + ':' + df['start'].astype(str) + '-' + df['end'].astype(str)
                region_col = 'region'
            else:
                raise ValueError(f"Cannot identify region column in {filepath}")
        
        # Identify value column
        value_cols = ['rms', 'methylation', 'meth', 'beta', 'value', 'MSets1.rms']
        value_col = None
        for col in value_cols:
            if col in df.columns:
                value_col = col
                break
        
        if value_col is None:
            numeric_cols = df.select_dtypes(include=[np.number]).columns
            if len(numeric_cols) > 0:
                value_col = numeric_cols[-1]
            else:
                raise ValueError(f"Cannot identify value column in {filepath}")
        
        # Create sample series
        sample_id = sample_ids[i] if sample_ids else Path(filepath).stem
        sample_data = df.set_index(region_col)[value_col]
        sample_data.name = sample_id
        all_samples.append(sample_data)
    
    # Combine all samples
    matrix = pd.concat(all_samples, axis=1).T
    print(f"Combined matrix: {matrix.shape[0]} samples x {matrix.shape[1]} regions")
    
    return matrix


def filter_missing(matrix: pd.DataFrame, max_missing: float = 0.1) -> pd.DataFrame:
    """Filter regions with too many missing values."""
    missing_frac = matrix.isna().mean()
    keep_cols = missing_frac[missing_frac <= max_missing].index
    matrix_filtered = matrix[keep_cols]
    print(f"After missing filter: {matrix_filtered.shape[1]} regions")
    return matrix_filtered


def impute_missing(matrix: pd.DataFrame, method: str = 'median') -> pd.DataFrame:
    """Impute remaining missing values."""
    if method == 'median':
        return matrix.fillna(matrix.median())
    elif method == 'mean':
        return matrix.fillna(matrix.mean())
    elif method == 'zero':
        return matrix.fillna(0)
    else:
        raise ValueError(f"Unknown imputation method: {method}")


def select_features_variance(matrix: pd.DataFrame, threshold: float = 0.01) -> tuple:
    """Select features by variance threshold."""
    selector = VarianceThreshold(threshold=threshold)
    matrix_selected = selector.fit_transform(matrix)
    
    selected_cols = matrix.columns[selector.get_support()]
    matrix_df = pd.DataFrame(matrix_selected, index=matrix.index, columns=selected_cols)
    
    print(f"After variance filter: {matrix_df.shape[1]} features")
    return matrix_df, selector


def select_features_kbest(matrix: pd.DataFrame, labels: pd.Series, k: int = 1000) -> tuple:
    """Select top k features by ANOVA F-value."""
    k = min(k, matrix.shape[1])
    
    selector = SelectKBest(score_func=f_classif, k=k)
    matrix_selected = selector.fit_transform(matrix, labels)
    
    selected_cols = matrix.columns[selector.get_support()]
    matrix_df = pd.DataFrame(matrix_selected, index=matrix.index, columns=selected_cols)
    
    print(f"After SelectKBest: {matrix_df.shape[1]} features")
    return matrix_df, selector


def scale_features(matrix: pd.DataFrame) -> tuple:
    """Standardize features to zero mean and unit variance."""
    scaler = StandardScaler()
    matrix_scaled = scaler.fit_transform(matrix)
    matrix_df = pd.DataFrame(matrix_scaled, index=matrix.index, columns=matrix.columns)
    return matrix_df, scaler


def main():
    parser = argparse.ArgumentParser(description='Prepare features for ML classification')
    parser.add_argument('-i', '--input', nargs='+', required=True, help='Input methylation files')
    parser.add_argument('-l', '--labels', required=True, help='Labels file (TSV with sample_id and label columns)')
    parser.add_argument('-o', '--output', required=True, help='Output prefix')
    parser.add_argument('--max_missing', type=float, default=0.1, help='Max missing fraction per feature')
    parser.add_argument('--impute', default='median', choices=['median', 'mean', 'zero'], help='Imputation method')
    parser.add_argument('--var_threshold', type=float, default=0.01, help='Variance threshold')
    parser.add_argument('--n_features', type=int, default=1000, help='Number of features to select')
    parser.add_argument('--feature_selection', default='kbest', choices=['variance', 'kbest', 'both'], 
                       help='Feature selection method')
    
    args = parser.parse_args()
    
    print("=== Feature Preparation for ML ===")
    print(f"Input files: {len(args.input)}")
    
    # Load labels
    labels_df = pd.read_csv(args.labels, sep='\t')
    if 'sample_id' not in labels_df.columns:
        labels_df.columns = ['sample_id', 'label'] + list(labels_df.columns[2:])
    
    sample_ids = labels_df['sample_id'].tolist()
    labels = labels_df.set_index('sample_id')['label']
    
    print(f"Labels: {labels.value_counts().to_dict()}")
    
    # Load methylation matrix
    matrix = load_methylation_matrix(args.input, sample_ids)
    
    # Align samples with labels
    common_samples = matrix.index.intersection(labels.index)
    matrix = matrix.loc[common_samples]
    labels = labels.loc[common_samples]
    
    print(f"Common samples: {len(common_samples)}")
    
    # Filter missing values
    matrix = filter_missing(matrix, args.max_missing)
    
    # Impute remaining missing
    matrix = impute_missing(matrix, args.impute)
    
    # Feature selection
    if args.feature_selection in ['variance', 'both']:
        matrix, var_selector = select_features_variance(matrix, args.var_threshold)
    
    if args.feature_selection in ['kbest', 'both']:
        matrix, kbest_selector = select_features_kbest(matrix, labels, args.n_features)
    
    # Scale features
    matrix_scaled, scaler = scale_features(matrix)
    
    # Save outputs
    # Feature matrix
    matrix_file = f"{args.output}.features.tsv"
    matrix_scaled.to_csv(matrix_file, sep='\t')
    print(f"Features saved to: {matrix_file}")
    
    # Labels (aligned)
    labels_file = f"{args.output}.labels.tsv"
    labels.to_frame().to_csv(labels_file, sep='\t')
    print(f"Labels saved to: {labels_file}")
    
    # Feature names
    features_file = f"{args.output}.feature_names.txt"
    with open(features_file, 'w') as f:
        f.write('\n'.join(matrix_scaled.columns))
    print(f"Feature names saved to: {features_file}")
    
    # Preprocessing objects
    preproc_file = f"{args.output}.preprocessing.pkl"
    preproc = {
        'scaler': scaler,
        'feature_names': list(matrix_scaled.columns),
        'n_features': matrix_scaled.shape[1],
        'n_samples': matrix_scaled.shape[0]
    }
    
    if args.feature_selection in ['variance', 'both']:
        preproc['var_selector'] = var_selector
    if args.feature_selection in ['kbest', 'both']:
        preproc['kbest_selector'] = kbest_selector
    
    with open(preproc_file, 'wb') as f:
        pickle.dump(preproc, f)
    print(f"Preprocessing objects saved to: {preproc_file}")
    
    print(f"\nFinal feature matrix: {matrix_scaled.shape}")
    print("Feature preparation complete!")


if __name__ == '__main__':
    main()
