#!/usr/bin/env python3
"""
ML Robustness Checks for cfMeDIP-seq Classification

Performs four key validation checks:
A) Autosomal-only analysis (exclude chrX/chrY/chrM)
B) Depth/complexity confounder detection
C) Permutation test for statistical significance
D) Feature driver chromosome distribution analysis

Author: MEDIPIPE
"""

import argparse
import pandas as pd
import numpy as np
from scipy import stats
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
import warnings
import json
from pathlib import Path

warnings.filterwarnings('ignore')


def load_data(matrix_path, annotation_path, case_group, control_group):
    """Load and prepare matrix and annotation data."""
    m = pd.read_csv(matrix_path, sep="\t")
    
    # Clean column names
    def clean_col(c):
        c = str(c).replace("#", "").replace("'", "").replace('"', "").strip()
        if c.endswith(".bw"):
            c = c[:-3]
        return c
    
    m.columns = [clean_col(c) for c in m.columns]
    
    ann = pd.read_csv(annotation_path, sep="\t")
    ann["sample_id"] = ann["sample_id"].astype(str)
    sample_groups = ann.set_index("sample_id")["group"]
    
    sample_cols = [c for c in m.columns if c in ann["sample_id"].values]
    case_samples = [s for s in sample_cols if sample_groups.get(s) == case_group]
    ctrl_samples = [s for s in sample_cols if sample_groups.get(s) == control_group]
    
    return m, sample_cols, case_samples, ctrl_samples, sample_groups


def run_loo_cv(X, y, k_features=5000):
    """Run leave-one-out cross-validation."""
    loo = LeaveOneOut()
    probs = []
    labels = []
    
    for train_idx, test_idx in loo.split(X):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        
        scaler = StandardScaler()
        X_train_s = scaler.fit_transform(X_train)
        X_test_s = scaler.transform(X_test)
        
        k = min(k_features, X_train_s.shape[1])
        selector = SelectKBest(f_classif, k=k)
        X_train_sel = selector.fit_transform(X_train_s, y_train)
        X_test_sel = selector.transform(X_test_s)
        
        model = LogisticRegression(C=1.0, solver="liblinear", class_weight="balanced", random_state=1)
        model.fit(X_train_sel, y_train)
        
        probs.append(model.predict_proba(X_test_sel)[0, 1])
        labels.append(y_test[0])
    
    return probs, labels


def check_autosomal(m, sample_cols, y):
    """Check A: Autosomal-only analysis."""
    autosomal_chrs = [f"chr{i}" for i in range(1, 23)]
    m_auto = m[m["chr"].isin(autosomal_chrs)].copy()
    
    n_excluded_x = len(m[m["chr"] == "chrX"])
    n_excluded_y = len(m[m["chr"] == "chrY"])
    
    X_auto = m_auto[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0).T.values
    probs, labels = run_loo_cv(X_auto, y, k_features=5000)
    auroc = roc_auc_score(labels, probs)
    
    return {
        "auroc_autosomal": auroc,
        "n_windows_autosomal": len(m_auto),
        "n_excluded_chrX": n_excluded_x,
        "n_excluded_chrY": n_excluded_y,
        "pct_autosomal": len(m_auto) / len(m) * 100
    }, probs


def check_depth_confounder(m, sample_cols, sample_groups, probs):
    """Check B: Depth/complexity as confounder."""
    X_full = m[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    
    depth_metrics = pd.DataFrame({
        "sample_id": sample_cols,
        "total_signal": X_full.sum(axis=0).values,
        "mean_signal": X_full.mean(axis=0).values,
        "nonzero_rate": (X_full > 0).sum(axis=0).values / len(X_full),
        "zero_rate": (X_full == 0).sum(axis=0).values / len(X_full),
        "group": [sample_groups.get(s) for s in sample_cols],
        "predicted_prob": probs
    })
    
    corr_signal, p_signal = stats.pearsonr(depth_metrics["total_signal"], depth_metrics["predicted_prob"])
    corr_nonzero, p_nonzero = stats.pearsonr(depth_metrics["nonzero_rate"], depth_metrics["predicted_prob"])
    
    groups = depth_metrics.groupby("group")["total_signal"]
    group_means = groups.mean().to_dict()
    
    return {
        "depth_pred_correlation_r": corr_signal,
        "depth_pred_correlation_p": p_signal,
        "nonzero_pred_correlation_r": corr_nonzero,
        "nonzero_pred_correlation_p": p_nonzero,
        "depth_confounder_detected": p_signal < 0.05 or p_nonzero < 0.05,
        "group_mean_depths": group_means
    }, depth_metrics


def check_permutation(X, y, observed_auroc, n_permutations=1000):
    """Check C: Permutation test."""
    np.random.seed(42)
    null_aurocs = []
    loo = LeaveOneOut()
    
    for i in range(n_permutations):
        y_perm = np.random.permutation(y)
        probs_perm = []
        
        for train_idx, test_idx in loo.split(X):
            X_train, X_test = X[train_idx], X[test_idx]
            y_train = y_perm[train_idx]
            
            scaler = StandardScaler()
            X_train_s = scaler.fit_transform(X_train)
            X_test_s = scaler.transform(X_test)
            
            k = min(1000, X_train_s.shape[1])
            selector = SelectKBest(f_classif, k=k)
            X_train_sel = selector.fit_transform(X_train_s, y_train)
            X_test_sel = selector.transform(X_test_s)
            
            model = LogisticRegression(C=1.0, solver="liblinear", class_weight="balanced", random_state=1)
            model.fit(X_train_sel, y_train)
            probs_perm.append(model.predict_proba(X_test_sel)[0, 1])
        
        try:
            null_aurocs.append(roc_auc_score(y_perm, probs_perm))
        except:
            null_aurocs.append(0.5)
    
    null_aurocs = np.array(null_aurocs)
    p_value = (null_aurocs >= observed_auroc).sum() / n_permutations
    
    return {
        "permutation_p_value": p_value,
        "null_auroc_mean": null_aurocs.mean(),
        "null_auroc_std": null_aurocs.std(),
        "null_auroc_95pct": np.percentile(null_aurocs, 95),
        "null_auroc_99pct": np.percentile(null_aurocs, 99),
        "significant_at_0.05": p_value < 0.05
    }


def check_feature_drivers(drivers_path):
    """Check D: Feature driver chromosome distribution."""
    drivers = pd.read_csv(drivers_path, sep="\t")
    
    # Extract chr from feature column if chr column doesn't exist
    if "chr" not in drivers.columns and "feature" in drivers.columns:
        drivers["chr"] = drivers["feature"].str.split(":").str[0]
    
    results = {}
    for n in [100, 500, 2000]:
        n_actual = min(n, len(drivers))
        top_n = drivers.head(n_actual)
        chr_counts = top_n["chr"].value_counts()
        sex_chr = chr_counts.get("chrX", 0) + chr_counts.get("chrY", 0)
        results[f"sex_chr_pct_top{n}"] = sex_chr / n_actual * 100
        results[f"chrX_top{n}"] = chr_counts.get("chrX", 0)
        results[f"chrY_top{n}"] = chr_counts.get("chrY", 0)
    
    top_100 = drivers.head(min(100, len(drivers)))
    results["hyper_in_top100"] = (top_100["cohen_d"] > 0).sum()
    results["hypo_in_top100"] = (top_100["cohen_d"] < 0).sum()
    results["sex_chr_enrichment_concern"] = results["sex_chr_pct_top100"] > 20
    
    return results


def main():
    parser = argparse.ArgumentParser(description="ML Robustness Checks")
    parser.add_argument("--matrix", required=True, help="Path to feature matrix")
    parser.add_argument("--annotation", required=True, help="Path to sample annotation")
    parser.add_argument("--drivers", required=True, help="Path to feature drivers file")
    parser.add_argument("--case-group", required=True, help="Case group name")
    parser.add_argument("--control-group", required=True, help="Control group name")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--n-permutations", type=int, default=1000, help="Number of permutations")
    args = parser.parse_args()
    
    Path(args.output).mkdir(parents=True, exist_ok=True)
    
    print("="*80)
    print("ML ROBUSTNESS CHECKS")
    print("="*80)
    
    # Load data
    print("\n[1] Loading data...")
    m, sample_cols, case_samples, ctrl_samples, sample_groups = load_data(
        args.matrix, args.annotation, args.case_group, args.control_group
    )
    y = np.array([1 if sample_groups.get(s) == args.case_group else 0 for s in sample_cols])
    
    print(f"    Samples: {len(case_samples)} {args.case_group}, {len(ctrl_samples)} {args.control_group}")
    print(f"    Windows: {len(m)}")
    
    # Check A: Autosomal-only
    print("\n[A] Autosomal-only analysis...")
    autosomal_results, probs_auto = check_autosomal(m, sample_cols, y)
    print(f"    AUROC (autosomal): {autosomal_results['auroc_autosomal']:.4f}")
    
    # Check B: Depth confounder
    print("\n[B] Depth confounder check...")
    depth_results, depth_df = check_depth_confounder(m, sample_cols, sample_groups, probs_auto)
    print(f"    Correlation (depth vs pred): r={depth_results['depth_pred_correlation_r']:.3f}, p={depth_results['depth_pred_correlation_p']:.4f}")
    
    # Check C: Permutation test
    print(f"\n[C] Permutation test ({args.n_permutations} iterations)...")
    autosomal_chrs = [f"chr{i}" for i in range(1, 23)]
    m_auto = m[m["chr"].isin(autosomal_chrs)]
    X_auto = m_auto[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0).T.values
    perm_results = check_permutation(X_auto, y, autosomal_results["auroc_autosomal"], args.n_permutations)
    print(f"    Permutation p-value: {perm_results['permutation_p_value']:.4f}")
    
    # Check D: Feature drivers
    print("\n[D] Feature driver analysis...")
    driver_results = check_feature_drivers(args.drivers)
    print(f"    Sex chr in top 100: {driver_results['sex_chr_pct_top100']:.1f}%")
    
    # Combine results
    all_results = {
        **autosomal_results,
        **depth_results,
        **perm_results,
        **driver_results
    }
    del all_results["group_mean_depths"]  # Not serializable simply
    
    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    checks_passed = sum([
        autosomal_results["auroc_autosomal"] > 0.7,
        not depth_results["depth_confounder_detected"],
        perm_results["significant_at_0.05"],
        not driver_results["sex_chr_enrichment_concern"]
    ])
    print(f"Checks passed: {checks_passed}/4")
    
    # Save results
    results_df = pd.DataFrame([all_results])
    results_df.to_csv(f"{args.output}/robustness_checks.tsv", sep="\t", index=False)
    
    depth_df.to_csv(f"{args.output}/sample_depth_metrics.tsv", sep="\t", index=False)
    
    with open(f"{args.output}/robustness_checks.json", "w") as f:
        json.dump(all_results, f, indent=2)
    
    print(f"\nResults saved to: {args.output}/")


if __name__ == "__main__":
    main()
