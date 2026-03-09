#!/usr/bin/env python3
"""
Validate File Compatibility

Checks that all MEDIPIPE output files have consistent:
1. Sample IDs between matrices and annotation
2. Header format (chr, start, end)
3. Column naming conventions
"""

import argparse
import pandas as pd
from pathlib import Path
import sys


def load_annotation(ann_file: str) -> set:
    """Load sample IDs from annotation file."""
    df = pd.read_csv(ann_file, sep='\t')
    return set(df['sample_id'].tolist()), df


def check_matrix(path: str, ann_samples: set) -> dict:
    """Check a single matrix file for compatibility."""
    result = {
        'path': path,
        'exists': False,
        'header_ok': False,
        'n_samples': 0,
        'samples_in_ann': 0,
        'missing_from_ann': [],
        'header': [],
    }
    
    if not Path(path).exists():
        return result
    
    result['exists'] = True
    
    df = pd.read_csv(path, sep='\t', nrows=1)
    result['header'] = list(df.columns[:3])
    
    # Check header format
    result['header_ok'] = (
        df.columns[0] == 'chr' and 
        df.columns[1] == 'start' and 
        df.columns[2] == 'end'
    )
    
    # Get sample columns
    samples = [c for c in df.columns if c not in ['chr', 'start', 'end']]
    result['n_samples'] = len(samples)
    
    # Check against annotation
    missing = set(samples) - ann_samples
    result['samples_in_ann'] = len(samples) - len(missing)
    result['missing_from_ann'] = sorted(missing)
    
    return result


def main():
    parser = argparse.ArgumentParser(description='Validate MEDIPIPE file compatibility')
    parser.add_argument('--output-dir', '-o', default='/mnt/data1/medipipe_data/output',
                       help='MEDIPIPE output directory')
    parser.add_argument('--annotation', '-a', 
                       default='/mnt/data1/medipipe_data/output/_upload_package_meta/sample_annotation.active.with_covariates.clean.tsv',
                       help='Sample annotation file')
    parser.add_argument('--verbose', '-v', action='store_true', help='Show detailed output')
    
    args = parser.parse_args()
    
    print("=== MEDIPIPE File Compatibility Validation ===\n")
    
    # Load annotation
    try:
        ann_samples, ann_df = load_annotation(args.annotation)
        print(f"✅ Annotation: {len(ann_samples)} samples")
        if args.verbose:
            print(f"   Groups: {ann_df['group'].value_counts().to_dict()}")
    except Exception as e:
        print(f"❌ Error loading annotation: {e}")
        sys.exit(1)
    
    # Define matrices to check
    output_dir = args.output_dir
    matrices = {
        "windows/matrix.tsv": f"{output_dir}/ml_discrimination/matrices/windows/hg38_w2000/all/matrix.tsv",
        "windows/matrix_autosomal.tsv": f"{output_dir}/ml_discrimination/matrices/windows/hg38_w2000/all/matrix_autosomal.tsv",
        "promoters/matrix.tsv": f"{output_dir}/ml_discrimination/matrices/aggregated/promoters/matrix.tsv",
        "gene_bodies/matrix.tsv": f"{output_dir}/ml_discrimination/matrices/aggregated/gene_bodies/matrix.tsv",
        "enhancers/matrix.tsv": f"{output_dir}/ml_discrimination/matrices/aggregated/enhancers/matrix.tsv",
        "cpg_islands/matrix.tsv": f"{output_dir}/ml_discrimination/matrices/aggregated/cpg_islands/matrix.tsv",
        "tissue_of_origin/marker_matrix.tsv": f"{output_dir}/tissue_of_origin/marker_matrix.tsv",
    }
    
    # Check each matrix
    print("\n--- Matrix Files ---")
    all_ok = True
    all_missing = set()
    
    for name, path in matrices.items():
        result = check_matrix(path, ann_samples)
        
        if not result['exists']:
            print(f"⏭️  {name}: Not found")
            continue
        
        issues = []
        if not result['header_ok']:
            issues.append(f"bad header {result['header']}")
            all_ok = False
        
        if result['missing_from_ann']:
            issues.append(f"{len(result['missing_from_ann'])} samples not in annotation")
            all_missing.update(result['missing_from_ann'])
            all_ok = False
        
        if issues:
            print(f"⚠️  {name}: {', '.join(issues)}")
            if args.verbose and result['missing_from_ann']:
                for s in result['missing_from_ann']:
                    print(f"       - {s}")
        else:
            print(f"✅ {name}: {result['n_samples']} samples")
    
    # Summary
    print("\n--- Summary ---")
    if all_ok:
        print("✅ All files are compatible!")
    else:
        print("⚠️  Issues found:")
        if all_missing:
            print(f"   - {len(all_missing)} samples missing from annotation")
            if args.verbose:
                for s in sorted(all_missing):
                    print(f"       {s}")
    
    # Return exit code
    sys.exit(0 if all_ok else 1)


if __name__ == '__main__':
    main()
