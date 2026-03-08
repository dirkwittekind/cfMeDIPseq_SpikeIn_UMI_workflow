#!/usr/bin/env python3
"""
Univariate Feature Driver Analysis
Ranks windows/features by effect size (Cohen's d, delta mean) for cfMeDIP-seq data.
Based on standard workflow for AEG vs CTRL discrimination.
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
import os


def load_annotation(ann_path: str, group_col: str = "group") -> pd.DataFrame:
    """Load and clean sample annotation file."""
    ann = pd.read_csv(ann_path, sep="\t", dtype=str)
    ann.columns = [c.strip() for c in ann.columns]
    ann["sample_id"] = ann["sample_id"].astype(str).str.strip()
    ann[group_col] = ann[group_col].astype(str).str.strip()
    # Map common control names to CTRL
    ann[group_col] = ann[group_col].replace({"Kontrolle": "CTRL", "Control": "CTRL", "control": "CTRL"})
    return ann


def load_matrix(matrix_path: str, sample_ids: list) -> pd.DataFrame:
    """Load feature matrix and orient it correctly (features x samples)."""
    m = pd.read_csv(matrix_path, sep="\t")
    
    # Clean column names: remove quotes, #, and .bw suffix from BigWig output
    def clean_col(c):
        c = str(c).strip("'").strip('"').replace("#", "").strip()
        if c.endswith(".bw"):
            c = c[:-3]
        return c
    
    m.columns = [clean_col(c) for c in m.columns]
    feat_col = m.columns[0]
    m[feat_col] = m[feat_col].astype(str).str.strip()
    
    # Check if samples are in columns
    sample_cols = [c for c in m.columns if c in sample_ids]
    
    if len(sample_cols) >= 2:
        # Samples in columns - standard orientation
        X = m[[feat_col] + sample_cols].set_index(feat_col)
    else:
        # Samples might be in rows - transpose
        X = m.set_index(feat_col).T
        X = X.loc[[s for s in sample_ids if s in X.index]]
        X = X.T
    
    X = X.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    
    # Align columns to sample_ids order
    available_samples = [s for s in sample_ids if s in X.columns]
    X = X[available_samples]
    
    return X


def calculate_statistics(X: pd.DataFrame, group_labels: pd.Series, 
                         case_group: str, ctrl_group: str) -> pd.DataFrame:
    """
    Calculate univariate statistics for each feature.
    
    Returns DataFrame with:
    - mean_case, mean_ctrl
    - delta (case - ctrl)
    - cohen_d
    - t_statistic, p_value
    """
    # Get indices for each group
    idx_case = (group_labels == case_group).values
    idx_ctrl = (group_labels == ctrl_group).values
    
    # Values arrays
    case_vals = X.values[:, idx_case]
    ctrl_vals = X.values[:, idx_ctrl]
    
    n_case = case_vals.shape[1]
    n_ctrl = ctrl_vals.shape[1]
    
    # Means
    mean_case = case_vals.mean(axis=1)
    mean_ctrl = ctrl_vals.mean(axis=1)
    delta = mean_case - mean_ctrl
    
    # Standard deviations
    sd_case = case_vals.std(axis=1, ddof=1)
    sd_ctrl = ctrl_vals.std(axis=1, ddof=1)
    
    # Pooled SD for Cohen's d
    pooled_sd = np.sqrt(((n_case - 1) * sd_case**2 + (n_ctrl - 1) * sd_ctrl**2) / (n_case + n_ctrl - 2))
    cohen_d = np.divide(delta, pooled_sd, out=np.zeros_like(delta), where=pooled_sd > 0)
    
    # T-test
    t_stats = np.zeros(X.shape[0])
    p_vals = np.ones(X.shape[0])
    
    for i in range(X.shape[0]):
        if sd_case[i] > 0 or sd_ctrl[i] > 0:
            t_stat, p_val = stats.ttest_ind(case_vals[i, :], ctrl_vals[i, :], equal_var=False)
            t_stats[i] = t_stat if not np.isnan(t_stat) else 0
            p_vals[i] = p_val if not np.isnan(p_val) else 1
    
    # Build output DataFrame
    result = pd.DataFrame({
        "feature": X.index,
        f"mean_{case_group}": mean_case,
        f"mean_{ctrl_group}": mean_ctrl,
        f"delta_{case_group}_minus_{ctrl_group}": delta,
        "cohen_d": cohen_d,
        "abs_cohen_d": np.abs(cohen_d),
        "t_statistic": t_stats,
        "p_value": p_vals,
        f"sd_{case_group}": sd_case,
        f"sd_{ctrl_group}": sd_ctrl,
    })
    
    # Multiple testing correction (Benjamini-Hochberg)
    result = result.sort_values("p_value")
    n_tests = len(result)
    result["rank"] = range(1, n_tests + 1)
    result["p_adjusted"] = result["p_value"] * n_tests / result["rank"]
    result["p_adjusted"] = result["p_adjusted"].clip(upper=1.0)
    # Make monotonic
    result["p_adjusted"] = result["p_adjusted"][::-1].cummin()[::-1]
    
    # Sort by absolute effect size
    result = result.sort_values("abs_cohen_d", ascending=False)
    result = result.drop(columns=["rank"])
    
    # Add direction interpretation
    result["direction"] = result["cohen_d"].apply(
        lambda x: f"hyper_in_{case_group}" if x > 0 else f"hypo_in_{case_group}"
    )
    
    return result


def annotate_with_genes(features_df: pd.DataFrame, gene_bed: str = None) -> pd.DataFrame:
    """
    Optionally annotate features with nearest genes.
    Feature IDs expected in format 'chr:start-end' or 'chr_start_end'.
    """
    if gene_bed is None or not os.path.exists(gene_bed):
        return features_df
    
    try:
        import subprocess
        import tempfile
        
        # Parse feature IDs to BED format
        bed_lines = []
        for feat in features_df["feature"]:
            # Try chr:start-end format
            if ":" in feat and "-" in feat:
                chrom, coords = feat.split(":")
                start, end = coords.split("-")
            elif "_" in feat:
                parts = feat.split("_")
                if len(parts) >= 3:
                    chrom = parts[0]
                    start = parts[1]
                    end = parts[2]
                else:
                    continue
            else:
                continue
            bed_lines.append(f"{chrom}\t{start}\t{end}\t{feat}")
        
        if not bed_lines:
            return features_df
        
        # Write temp BED file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
            f.write("\n".join(bed_lines))
            temp_bed = f.name
        
        # Run bedtools closest
        result = subprocess.run(
            ["bedtools", "closest", "-a", temp_bed, "-b", gene_bed, "-d"],
            capture_output=True, text=True
        )
        
        os.unlink(temp_bed)
        
        if result.returncode == 0:
            # Parse output
            gene_map = {}
            for line in result.stdout.strip().split("\n"):
                if line:
                    parts = line.split("\t")
                    if len(parts) >= 8:
                        feat_id = parts[3]
                        gene_name = parts[7] if len(parts) > 7 else "."
                        distance = parts[-1] if parts[-1].isdigit() else "0"
                        gene_map[feat_id] = {"nearest_gene": gene_name, "gene_distance": int(distance)}
            
            features_df["nearest_gene"] = features_df["feature"].map(lambda x: gene_map.get(x, {}).get("nearest_gene", "."))
            features_df["gene_distance"] = features_df["feature"].map(lambda x: gene_map.get(x, {}).get("gene_distance", -1))
    
    except Exception as e:
        print(f"Warning: Could not annotate genes: {e}")
    
    return features_df


def main():
    parser = argparse.ArgumentParser(
        description="Rank windows/features by univariate effect size for cfMeDIP-seq discrimination"
    )
    parser.add_argument("--ann", "-a", required=True, help="Sample annotation TSV (sample_id, group columns required)")
    parser.add_argument("--matrix", "-m", required=True, help="Feature matrix TSV (features x samples)")
    parser.add_argument("--out", "-o", required=True, help="Output file path")
    parser.add_argument("--case-group", default="AEG", help="Case group name (default: AEG)")
    parser.add_argument("--ctrl-group", default="CTRL", help="Control group name (default: CTRL)")
    parser.add_argument("--top-n", type=int, default=2000, help="Number of top features to output (default: 2000)")
    parser.add_argument("--gene-bed", help="Optional gene BED file for annotation")
    parser.add_argument("--all-features", action="store_true", help="Output all features, not just top N")
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("Univariate Feature Driver Analysis")
    print("=" * 60)
    
    # Load annotation
    print(f"\nLoading annotation: {args.ann}")
    ann = load_annotation(args.ann)
    print(f"  Samples: {len(ann)}")
    print(f"  Groups: {ann['group'].value_counts().to_dict()}")
    
    keep = ann["sample_id"].tolist()
    
    # Load matrix
    print(f"\nLoading matrix: {args.matrix}")
    X = load_matrix(args.matrix, keep)
    print(f"  Features: {X.shape[0]}")
    print(f"  Samples loaded: {X.shape[1]}")
    
    # Get group labels aligned to matrix columns
    ann_aligned = ann.set_index("sample_id").loc[X.columns]
    group_labels = ann_aligned["group"]
    
    # Check groups
    unique_groups = group_labels.unique()
    print(f"\nGroups in data: {list(unique_groups)}")
    
    if args.case_group not in unique_groups:
        print(f"Warning: Case group '{args.case_group}' not found. Available: {list(unique_groups)}")
        args.case_group = [g for g in unique_groups if g != args.ctrl_group][0]
        print(f"Using case group: {args.case_group}")
    
    if args.ctrl_group not in unique_groups:
        print(f"Warning: Control group '{args.ctrl_group}' not found. Available: {list(unique_groups)}")
        args.ctrl_group = [g for g in unique_groups if g != args.case_group][0]
        print(f"Using control group: {args.ctrl_group}")
    
    # Calculate statistics
    print(f"\nCalculating statistics: {args.case_group} vs {args.ctrl_group}")
    result = calculate_statistics(X, group_labels, args.case_group, args.ctrl_group)
    
    # Annotate with genes if requested
    if args.gene_bed:
        print(f"\nAnnotating with genes: {args.gene_bed}")
        result = annotate_with_genes(result, args.gene_bed)
    
    # Output
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    
    if args.all_features:
        output_df = result
    else:
        output_df = result.head(args.top_n)
    
    output_df.to_csv(args.out, sep="\t", index=False)
    print(f"\nWrote {len(output_df)} features to: {args.out}")
    
    # Summary statistics
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    sig_features = result[result["p_adjusted"] < 0.05]
    print(f"Significant features (padj < 0.05): {len(sig_features)}")
    
    hyper = result[result["cohen_d"] > 0.5]
    hypo = result[result["cohen_d"] < -0.5]
    print(f"Strong effect (|d| > 0.5): {len(hyper)} hyper, {len(hypo)} hypo in {args.case_group}")
    
    print("\nTop 10 by effect size:")
    top10 = result.head(10)[["feature", "cohen_d", "p_adjusted", "direction"]]
    print(top10.to_string(index=False))


if __name__ == "__main__":
    main()
