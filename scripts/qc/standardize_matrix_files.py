#!/usr/bin/env python3
"""
Standardize Matrix File Naming Conventions

This script ensures consistent naming across all coverage/annotation matrix files:
1. Column headers: chr, start, end (no quotes, no #)
2. Sample columns: Short sample IDs without .bw suffix (e.g., "148" not "'148.bw'")
3. Consistent with sample annotation file
"""

import argparse
import pandas as pd
import re
import os
from pathlib import Path


def extract_sample_id(col_name: str) -> str:
    """
    Extract clean sample ID from various formats.
    
    Handles:
    - '148.bw' -> 148
    - '014-CSF_S1_L003.bw' -> 014-CSF
    - 148_S6_L004 -> 148
    - 2001CSF.bw -> 2001CSF
    """
    # Remove quotes
    col_name = col_name.strip("'\"")
    
    # Remove .bw suffix
    col_name = re.sub(r'\.bw$', '', col_name)
    
    # Handle Illumina sample sheet naming (_S*_L00*)
    # Pattern: SAMPLE_S#_L### -> SAMPLE
    match = re.match(r'^(.+?)_S\d+_L\d+$', col_name)
    if match:
        return match.group(1)
    
    # Handle CSF samples with lane info (014-CSF_S1_L003 -> 014-CSF)
    match = re.match(r'^(\d{3}-CSF)_S\d+_L\d+$', col_name)
    if match:
        return match.group(1)
    
    return col_name


def standardize_header(col_name: str) -> str:
    """Standardize coordinate column headers."""
    col_lower = col_name.strip("'\"#").lower()
    
    if col_lower == 'chr' or col_lower == 'chrom':
        return 'chr'
    elif col_lower == 'start':
        return 'start'
    elif col_lower in ['end', 'stop']:
        return 'end'
    else:
        # Sample column - extract clean ID
        return extract_sample_id(col_name)


def standardize_matrix(input_file: str, output_file: str = None, 
                       sample_mapping: dict = None, dry_run: bool = False) -> dict:
    """
    Standardize a matrix file.
    
    Args:
        input_file: Path to input matrix TSV
        output_file: Path to output (default: overwrite input)
        sample_mapping: Optional dict mapping old names to new names
        dry_run: If True, only report changes without writing
    
    Returns:
        Dict with standardization report
    """
    if output_file is None:
        output_file = input_file
    
    # Detect comment character
    with open(input_file, 'r') as f:
        first_line = f.readline()
        has_comment = first_line.startswith('#')
    
    # Read matrix
    df = pd.read_csv(input_file, sep='\t', comment=None if not has_comment else None)
    
    # Get original columns
    orig_cols = list(df.columns)
    
    # Standardize column names
    new_cols = []
    changes = []
    
    for col in orig_cols:
        new_col = standardize_header(col)
        
        # Apply custom mapping if provided
        if sample_mapping and new_col in sample_mapping:
            new_col = sample_mapping[new_col]
        
        if col != new_col:
            changes.append((col, new_col))
        
        new_cols.append(new_col)
    
    # Check for duplicates after standardization
    if len(new_cols) != len(set(new_cols)):
        # Find duplicates
        seen = {}
        for i, col in enumerate(new_cols):
            if col in seen:
                print(f"WARNING: Duplicate column after standardization: {col}")
                print(f"  Original: {orig_cols[seen[col]]} and {orig_cols[i]}")
            else:
                seen[col] = i
    
    # Apply new column names
    df.columns = new_cols
    
    # Report
    report = {
        'file': input_file,
        'n_cols': len(orig_cols),
        'n_changes': len(changes),
        'changes': changes,
        'coord_cols': [c for c in new_cols[:3]],
        'sample_cols': new_cols[3:],
    }
    
    if dry_run:
        print(f"\n=== {input_file} ===")
        print(f"Columns: {len(orig_cols)}")
        print(f"Changes needed: {len(changes)}")
        for old, new in changes[:10]:
            print(f"  {old} -> {new}")
        if len(changes) > 10:
            print(f"  ... and {len(changes) - 10} more")
    else:
        # Write output
        df.to_csv(output_file, sep='\t', index=False)
        print(f"Standardized: {output_file} ({len(changes)} changes)")
    
    return report


def create_sample_mapping(annotation_file: str) -> dict:
    """
    Create sample ID mapping from annotation file.
    
    Maps various formats to the canonical sample_id in the annotation.
    """
    df = pd.read_csv(annotation_file, sep='\t')
    
    mapping = {}
    for sample_id in df['sample_id']:
        # Base ID (without _S*_L*)
        base_id = extract_sample_id(sample_id)
        
        # Map various forms to canonical
        mapping[base_id] = base_id  # e.g., 148 -> 148
        mapping[f"{base_id}.bw"] = base_id  # e.g., 148.bw -> 148
        mapping[f"'{base_id}.bw'"] = base_id  # e.g., '148.bw' -> 148
        mapping[sample_id] = base_id  # e.g., 148_S6_L004 -> 148
    
    return mapping


def main():
    parser = argparse.ArgumentParser(description='Standardize matrix file naming conventions')
    parser.add_argument('files', nargs='+', help='Matrix files to standardize')
    parser.add_argument('--annotation', '-a', help='Sample annotation file for ID mapping')
    parser.add_argument('--dry-run', '-n', action='store_true', help='Only report, do not modify')
    parser.add_argument('--backup', '-b', action='store_true', help='Create .bak backup before modifying')
    
    args = parser.parse_args()
    
    # Load sample mapping if annotation provided
    sample_mapping = None
    if args.annotation:
        sample_mapping = create_sample_mapping(args.annotation)
        print(f"Loaded {len(sample_mapping)} sample ID mappings from annotation")
    
    # Process each file
    all_reports = []
    for filepath in args.files:
        if not os.path.exists(filepath):
            print(f"WARNING: File not found: {filepath}")
            continue
        
        if args.backup and not args.dry_run:
            backup_path = filepath + '.bak'
            if not os.path.exists(backup_path):
                import shutil
                shutil.copy(filepath, backup_path)
                print(f"Created backup: {backup_path}")
        
        report = standardize_matrix(
            filepath, 
            sample_mapping=sample_mapping,
            dry_run=args.dry_run
        )
        all_reports.append(report)
    
    # Summary
    print(f"\n=== Summary ===")
    print(f"Files processed: {len(all_reports)}")
    total_changes = sum(r['n_changes'] for r in all_reports)
    print(f"Total changes: {total_changes}")
    
    if args.dry_run:
        print("\n(Dry run - no files modified)")


if __name__ == '__main__':
    main()
