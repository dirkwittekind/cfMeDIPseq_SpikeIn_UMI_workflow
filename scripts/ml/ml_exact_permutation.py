#!/usr/bin/env python3
"""
Exact Permutation Test for cfMeDIP-seq ML Classification

Enumerates ALL possible labelings preserving class sizes (n=6 vs n=5 → 462 labelings)
and runs full LOO-CV pipeline for each to compute exact p-values.

Key features:
- Strict leakage control: ALL preprocessing inside each fold
- Dual model support: L2 Logistic Regression + Elastic-Net (parallel)
- Parallelized labeling enumeration using joblib (scalable to 48+ cores)
- Computes exact p-values for AUROC and AUPRC
- Detects if signal is driven by single sample (leave-one-out stability)

Author: MEDIPIPE
"""

import argparse
import itertools
import json
import math
import os
from dataclasses import dataclass, asdict
from typing import List, Tuple, Dict, Optional
from functools import partial

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.feature_selection import SelectKBest, VarianceThreshold, f_classif
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.metrics import average_precision_score, roc_auc_score, roc_curve
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler
import warnings

warnings.filterwarnings('ignore')

# Model configurations
# NOTE: SGDClassifier alpha is NOT equivalent to 1/C from LogisticRegression.
# For cfMeDIP-seq data with 1.5M features, alpha=0.01 works well (tested empirically).
MODEL_CONFIGS = {
    'L2': {
        'class': LogisticRegression,
        'params': {
            'penalty': 'l2',
            'solver': 'liblinear',
            'class_weight': 'balanced',
            'max_iter': 5000,
            'random_state': 42
        }
    },
    'ElasticNet': {
        'class': SGDClassifier,
        'params': {
            'loss': 'log_loss',
            'penalty': 'elasticnet',
            'l1_ratio': 0.5,  # Balance between L1 and L2
            'alpha': 0.01,    # Regularization strength (empirically tuned for cfMeDIP-seq)
            'class_weight': 'balanced',
            'max_iter': 5000,
            'random_state': 42,
            'early_stopping': False
        }
    }
}


def load_matrix(matrix_path: str) -> pd.DataFrame:
    """Load feature matrix and return samples x features DataFrame."""
    m = pd.read_csv(matrix_path, sep="\t")
    
    # Clean column names (handle BigWig naming quirks)
    def clean_col(c):
        c = str(c).replace("#", "").replace("'", "").replace('"', "").strip()
        if c.endswith(".bw"):
            c = c[:-3]
        return c
    
    m.columns = [clean_col(c) for c in m.columns]
    
    # Identify coordinate columns vs sample columns
    coord_cols = [c for c in m.columns if c in ["chr", "start", "end"]]
    
    if coord_cols:
        # Create feature IDs from coordinates
        m["feature_id"] = m["chr"].astype(str) + ":" + m["start"].astype(str) + "-" + m["end"].astype(str)
        sample_cols = [c for c in m.columns if c not in coord_cols + ["feature_id"]]
        X = m.set_index("feature_id")[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    else:
        # Assume first column is feature ID
        feature_col = m.columns[0]
        X = m.set_index(feature_col).apply(pd.to_numeric, errors="coerce").fillna(0.0)
    
    return X.T  # samples x features


def load_annotation(annotation_path: str, case_group: str, ctrl_group: str) -> pd.DataFrame:
    """Load sample annotation and filter to case/control groups."""
    ann = pd.read_csv(annotation_path, sep="\t", dtype=str).fillna("")
    ann.columns = [c.strip() for c in ann.columns]
    ann["sample_id"] = ann["sample_id"].astype(str).str.strip()
    ann["group"] = ann["group"].astype(str).str.strip()
    
    # Filter to case/control only
    ann = ann[ann["group"].isin([case_group, ctrl_group])].copy()
    return ann


def fit_fold_pipeline(Xtr: np.ndarray, ytr: np.ndarray, top_k: int, C: float, 
                      model_type: str = 'L2', alpha: float = 0.01):
    """
    Fit preprocessing + classifier pipeline on training fold.
    
    ALL steps are fit on training data only to prevent leakage:
    1. Z-standardization (FIRST - critical for stable F-scores)
    2. Univariate k-best selection (ANOVA F-score)
    3. Classifier (L2 LogReg or Elastic-Net)
    
    Args:
        Xtr: Training features
        ytr: Training labels
        top_k: Number of top features to select
        C: Regularization parameter (for L2) or inverse alpha (for ElasticNet)
        model_type: 'L2' or 'ElasticNet'
    
    Note: Scaling before selection is important for numerical stability
    with high-dimensional cfMeDIP-seq data.
    """
    # Step 1: Z-standardization FIRST (improves F-score stability)
    scaler = StandardScaler()
    Xtr_s = scaler.fit_transform(Xtr)
    
    # Step 2: Univariate k-best selection
    k = min(top_k, Xtr_s.shape[1])
    sel = SelectKBest(score_func=f_classif, k=k)
    Xtr_sel = sel.fit_transform(Xtr_s, ytr)
    
    if Xtr_sel.shape[1] == 0:
        raise ValueError("No features selected")
    
    # Step 3: Fit classifier based on model type
    config = MODEL_CONFIGS[model_type]
    clf = config['class'](**config['params'])
    
    # Adjust regularization parameter
    if model_type == 'L2':
        clf.C = C
    else:  # ElasticNet
        clf.alpha = alpha  # Use alpha from command line argument
    
    clf.fit(Xtr_sel, ytr)
    
    return scaler, sel, clf


def transform_with_pipeline(X: np.ndarray, scaler, sel) -> np.ndarray:
    """Apply fitted preprocessing pipeline to new data."""
    Xs = scaler.transform(X)
    Xsel = sel.transform(Xs)
    return Xsel


def loo_oof_predictions(X: pd.DataFrame, y: np.ndarray, top_k: int, C: float,
                        model_type: str = 'L2', alpha: float = 0.01) -> Tuple[np.ndarray, np.ndarray]:
    """
    Run leave-one-out CV and return out-of-fold predictions.
    
    Args:
        X: Feature matrix (samples x features)
        y: Labels
        top_k: Number of top features to select
        C: Regularization parameter
        model_type: 'L2' or 'ElasticNet'
    
    Returns:
        scores: predicted probability of positive class for each sample
        labels: true labels
    """
    loo = LeaveOneOut()
    n = len(y)
    scores = np.zeros(n, dtype=float)
    
    for train_idx, test_idx in loo.split(X):
        Xtr = X.iloc[train_idx].values
        Xte = X.iloc[test_idx].values
        ytr = y[train_idx]
        
        try:
            scaler, sel, clf = fit_fold_pipeline(Xtr, ytr, top_k=top_k, C=C, model_type=model_type, alpha=alpha)
            Xte_trans = transform_with_pipeline(Xte, scaler, sel)
            # Handle both predict_proba (LogisticRegression) and decision_function (SGDClassifier)
            if hasattr(clf, 'predict_proba'):
                scores[test_idx[0]] = clf.predict_proba(Xte_trans)[0, 1]
            else:
                # Convert decision function to probability-like score
                decision = clf.decision_function(Xte_trans)[0]
                scores[test_idx[0]] = 1 / (1 + np.exp(-decision))  # Sigmoid
        except Exception as e:
            # If pipeline fails (e.g., no features selected), assign 0.5
            scores[test_idx[0]] = 0.5
    
    return scores, y


def _evaluate_single_labeling(labeling_tuple: Tuple[int, np.ndarray], X: pd.DataFrame, 
                               y_true: np.ndarray, top_k: int, C: float, 
                               model_type: str = 'L2', alpha: float = 0.01) -> dict:
    """
    Evaluate a single labeling permutation (for parallel execution).
    
    Args:
        labeling_tuple: (index, permuted_labels)
        X: Feature matrix
        y_true: True labels (for checking if this is the observed labeling)
        top_k: Number of features to select
        C: Regularization parameter
        model_type: Model type to use
    
    Returns:
        Dict with labeling_id, auroc, auprc, is_observed_labeling
    """
    i, yp = labeling_tuple
    
    scores_p, _ = loo_oof_predictions(X, yp, top_k=top_k, C=C, model_type=model_type, alpha=alpha)
    
    try:
        auc_p = roc_auc_score(yp, scores_p)
        auprc_p = average_precision_score(yp, scores_p)
    except:
        n_pos = int(yp.sum())
        auc_p = 0.5
        auprc_p = n_pos / len(yp)  # baseline AUPRC
    
    is_observed = np.array_equal(yp, y_true)
    
    return {
        "labeling_id": i,
        "auroc": auc_p,
        "auprc": auprc_p,
        "is_observed_labeling": is_observed
    }


def all_labelings_same_class_size(n: int, n_pos: int):
    """
    Generate all unique labelings preserving class sizes.
    
    For n=11 and n_pos=6: C(11,6) = 462 unique labelings
    """
    for pos_idx in itertools.combinations(range(n), n_pos):
        y = np.zeros(n, dtype=int)
        y[list(pos_idx)] = 1
        yield y


def exact_permutation_test(X: pd.DataFrame, y: np.ndarray, top_k: int, C: float, 
                           model_type: str = 'L2', alpha: float = 0.01, n_jobs: int = -1,
                           verbose: bool = True) -> dict:
    """
    Run exact permutation test by enumerating all possible labelings.
    
    PARALLELIZED version using joblib for ~10x speedup on multi-core systems.
    
    Args:
        X: Feature matrix (samples x features)
        y: Labels (binary)
        top_k: Number of top features to select
        C: Regularization parameter
        model_type: 'L2' or 'ElasticNet'
        n_jobs: Number of parallel jobs (-1 = all cores)
        verbose: Print progress
    
    Returns dict with:
    - observed metrics (AUROC, AUPRC)
    - null distribution
    - exact p-values
    - per-sample stability metrics
    """
    n = len(y)
    n_pos = int(y.sum())
    n_labelings = math.comb(n, n_pos)
    
    if verbose:
        print(f"Running exact permutation test ({model_type}):")
        print(f"  n={n}, n_pos={n_pos}, n_neg={n-n_pos}")
        print(f"  Total labelings: C({n},{n_pos}) = {n_labelings}")
        print(f"  Parallel jobs: {n_jobs}")
    
    # Observed performance
    obs_scores, obs_labels = loo_oof_predictions(X, y, top_k=top_k, C=C, model_type=model_type, alpha=alpha)
    obs_auroc = roc_auc_score(obs_labels, obs_scores)
    obs_auprc = average_precision_score(obs_labels, obs_scores)
    
    if verbose:
        print(f"\n  Observed LOO-CV: AUROC={obs_auroc:.4f}, AUPRC={obs_auprc:.4f}")
        print(f"\n  Enumerating all {n_labelings} labelings (PARALLEL)...")
    
    # Generate all labelings as list of tuples (index, labels)
    all_labelings = [(i, yp.copy()) for i, yp in enumerate(all_labelings_same_class_size(n, n_pos), start=1)]
    
    # PARALLEL execution of permutation test
    null_rows = Parallel(n_jobs=n_jobs, verbose=10 if verbose else 0)(
        delayed(_evaluate_single_labeling)(labeling, X, y, top_k, C, model_type, alpha)
        for labeling in all_labelings
    )
    
    null_df = pd.DataFrame(null_rows)
    
    # Compute exact p-values (proportion of labelings with metric >= observed)
    p_auroc = (null_df["auroc"] >= obs_auroc).mean()
    p_auprc = (null_df["auprc"] >= obs_auprc).mean()
    
    if verbose:
        print(f"\n  Computing per-sample stability (PARALLEL)...")
    
    # Per-sample leave-one-out stability (also parallelized)
    sample_ids = X.index.tolist()
    
    def _compute_stability_for_sample(i):
        X_loo = X.drop(X.index[i])
        y_loo = np.delete(y, i)
        try:
            scores_loo, _ = loo_oof_predictions(X_loo, y_loo, top_k=top_k, C=C, model_type=model_type, alpha=alpha)
            auc_loo = roc_auc_score(y_loo, scores_loo)
        except:
            auc_loo = 0.5
        return {
            "sample_id": sample_ids[i],
            "label": int(y[i]),
            "auroc_without": auc_loo,
            "auroc_drop": obs_auroc - auc_loo
        }
    
    stability_results = Parallel(n_jobs=n_jobs, verbose=0)(
        delayed(_compute_stability_for_sample)(i) for i in range(n)
    )
    stability_df = pd.DataFrame(stability_results)
    
    if verbose:
        print(f"\n  Results:")
        print(f"    Observed AUROC: {obs_auroc:.4f}")
        print(f"    Exact p-value (AUROC): {p_auroc:.4f}")
        print(f"    Observed AUPRC: {obs_auprc:.4f}")
        print(f"    Exact p-value (AUPRC): {p_auprc:.4f}")
        print(f"    Null distribution: mean={null_df['auroc'].mean():.3f}, std={null_df['auroc'].std():.3f}")
        
        # Check for single-sample driver
        max_drop = stability_df["auroc_drop"].max()
        if max_drop > 0.15:
            driver_sample = stability_df.loc[stability_df["auroc_drop"].idxmax(), "sample_id"]
            print(f"\n  ⚠️ WARNING: Sample {driver_sample} may be driving classification (AUROC drop={max_drop:.3f})")
    
    return {
        "model_type": model_type,
        "observed_auroc": obs_auroc,
        "observed_auprc": obs_auprc,
        "exact_p_auroc": p_auroc,
        "exact_p_auprc": p_auprc,
        "n_labelings": n_labelings,
        "null_auroc_mean": null_df["auroc"].mean(),
        "null_auroc_std": null_df["auroc"].std(),
        "null_auroc_95pct": null_df["auroc"].quantile(0.95),
        "null_df": null_df,
        "stability_df": stability_df,
        "observed_scores": obs_scores,
        "observed_labels": obs_labels,
        "sample_ids": X.index.tolist()
    }


def exact_permutation_test_dual_model(X: pd.DataFrame, y: np.ndarray, top_k: int, C: float,
                                       alpha: float = 0.01, n_jobs: int = -1, verbose: bool = True) -> Dict[str, dict]:
    """
    Run exact permutation test with BOTH L2 and Elastic-Net models in parallel.
    
    Returns dict with results for each model type.
    """
    results = {}
    
    for model_type in ['L2', 'ElasticNet']:
        if verbose:
            print(f"\n{'='*70}")
            print(f"  MODEL: {model_type}")
            print(f"{'='*70}")
        
        results[model_type] = exact_permutation_test(
            X, y, top_k=top_k, C=C, model_type=model_type, alpha=alpha, n_jobs=n_jobs, verbose=verbose
        )
    
    # Compare models
    if verbose:
        print(f"\n{'='*70}")
        print(f"  MODEL COMPARISON")
        print(f"{'='*70}")
        print(f"  {'Model':<12} {'AUROC':<8} {'p-value':<10} {'Significant?'}")
        print(f"  {'-'*50}")
        for model_type, res in results.items():
            sig = "✓ YES" if res['exact_p_auroc'] < 0.05 else "✗ No"
            print(f"  {model_type:<12} {res['observed_auroc']:<8.4f} {res['exact_p_auroc']:<10.4f} {sig}")
    
    return results


def make_plots(results: dict, output_dir: str, feature_space: str):
    """Generate diagnostic plots."""
    os.makedirs(output_dir, exist_ok=True)
    
    null_df = results["null_df"]
    obs_auroc = results["observed_auroc"]
    obs_auprc = results["observed_auprc"]
    
    # Plot 1: AUROC null distribution
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(null_df["auroc"], bins=30, alpha=0.7, color="steelblue", edgecolor="black")
    ax.axvline(obs_auroc, color="red", linestyle="--", linewidth=2, 
               label=f"Observed={obs_auroc:.3f}")
    ax.axvline(null_df["auroc"].quantile(0.95), color="orange", linestyle=":", linewidth=2,
               label=f"95th pct={null_df['auroc'].quantile(0.95):.3f}")
    ax.set_xlabel("AUROC", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.set_title(f"{feature_space}: Exact Permutation Null Distribution (AUROC)\n"
                 f"n={results['n_labelings']} labelings, p={results['exact_p_auroc']:.4f}", fontsize=11)
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{feature_space}.permutation_auroc.pdf"))
    plt.savefig(os.path.join(output_dir, f"{feature_space}.permutation_auroc.png"), dpi=150)
    plt.close()
    
    # Plot 2: AUPRC null distribution
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(null_df["auprc"], bins=30, alpha=0.7, color="seagreen", edgecolor="black")
    ax.axvline(obs_auprc, color="red", linestyle="--", linewidth=2,
               label=f"Observed={obs_auprc:.3f}")
    ax.set_xlabel("AUPRC", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.set_title(f"{feature_space}: Exact Permutation Null Distribution (AUPRC)\n"
                 f"p={results['exact_p_auprc']:.4f}", fontsize=11)
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{feature_space}.permutation_auprc.pdf"))
    plt.close()
    
    # Plot 3: Sample stability (leave-one-out impact)
    stability_df = results["stability_df"]
    fig, ax = plt.subplots(figsize=(10, 5))
    colors = ["red" if l == 1 else "blue" for l in stability_df["label"]]
    bars = ax.bar(range(len(stability_df)), stability_df["auroc_drop"], color=colors, alpha=0.7)
    ax.axhline(0, color="black", linewidth=0.5)
    ax.axhline(0.15, color="orange", linestyle="--", linewidth=1, label="Warning threshold")
    ax.set_xticks(range(len(stability_df)))
    ax.set_xticklabels(stability_df["sample_id"], rotation=45, ha="right")
    ax.set_xlabel("Sample", fontsize=12)
    ax.set_ylabel("AUROC drop when removed", fontsize=12)
    ax.set_title(f"{feature_space}: Leave-One-Out Sample Stability\n"
                 f"Red=Case, Blue=Control", fontsize=11)
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{feature_space}.sample_stability.pdf"))
    plt.savefig(os.path.join(output_dir, f"{feature_space}.sample_stability.png"), dpi=150)
    plt.close()
    
    # Plot 4: ROC curve with observed predictions
    fig, ax = plt.subplots(figsize=(6, 6))
    fpr, tpr, _ = roc_curve(results["observed_labels"], results["observed_scores"])
    ax.plot(fpr, tpr, color="darkblue", linewidth=2, 
            label=f"LOO-CV ROC (AUROC={obs_auroc:.3f})")
    ax.plot([0, 1], [0, 1], color="gray", linestyle="--", linewidth=1)
    ax.set_xlabel("False Positive Rate", fontsize=12)
    ax.set_ylabel("True Positive Rate", fontsize=12)
    ax.set_title(f"{feature_space}: ROC Curve (LOO-CV)", fontsize=11)
    ax.legend(loc="lower right")
    ax.set_aspect("equal")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{feature_space}.roc_curve.pdf"))
    plt.close()


def _save_results(results: dict, args, model_type: str = None):
    """
    Save results to files. Helper for both single and dual-model modes.
    """
    # Determine file prefix
    if model_type:
        prefix = f"{args.feature_space}.{model_type}"
    else:
        prefix = args.feature_space
    
    print(f"\nSaving results to {args.output}/{prefix}.*")
    
    # Null distribution
    results["null_df"].to_csv(
        os.path.join(args.output, f"{prefix}.exact_null.tsv"),
        sep="\t", index=False
    )
    
    # Sample stability
    results["stability_df"].to_csv(
        os.path.join(args.output, f"{prefix}.sample_stability.tsv"),
        sep="\t", index=False
    )
    
    # Summary metrics
    summary = {
        "feature_space": args.feature_space,
        "model_type": model_type or results.get('model_type', 'L2'),
        "n_samples": len(results["observed_labels"]),
        "n_case": int(results["observed_labels"].sum()),
        "n_ctrl": len(results["observed_labels"]) - int(results["observed_labels"].sum()),
        "top_k": args.top_k,
        "C": args.C,
        "observed_auroc": results["observed_auroc"],
        "observed_auprc": results["observed_auprc"],
        "exact_p_auroc": results["exact_p_auroc"],
        "exact_p_auprc": results["exact_p_auprc"],
        "n_labelings": results["n_labelings"],
        "null_auroc_mean": results["null_auroc_mean"],
        "null_auroc_std": results["null_auroc_std"],
        "null_auroc_95pct": results["null_auroc_95pct"],
        "max_sample_auroc_drop": results["stability_df"]["auroc_drop"].max(),
        "driver_sample_warning": results["stability_df"]["auroc_drop"].max() > 0.15
    }
    
    pd.DataFrame([summary]).to_csv(
        os.path.join(args.output, f"{prefix}.summary.tsv"),
        sep="\t", index=False
    )
    
    with open(os.path.join(args.output, f"{prefix}.exact_p_values.json"), "w") as f:
        json.dump({
            "feature_space": args.feature_space,
            "model_type": model_type or results.get('model_type', 'L2'),
            "observed_auroc": float(results["observed_auroc"]),
            "observed_auprc": float(results["observed_auprc"]),
            "exact_p_auroc": float(results["exact_p_auroc"]),
            "exact_p_auprc": float(results["exact_p_auprc"]),
            "n_labelings": int(results["n_labelings"]),
            "significant_at_0.05": bool(results["exact_p_auroc"] < 0.05)
        }, f, indent=2)
    
    # Generate plots
    make_plots(results, args.output, prefix)
    
    print(f"  Observed AUROC: {results['observed_auroc']:.4f}")
    print(f"  Exact p-value: {results['exact_p_auroc']:.4f}")
    print(f"  Significant at α=0.05: {'YES ✓' if results['exact_p_auroc'] < 0.05 else 'NO'}")


def main():
    parser = argparse.ArgumentParser(
        description="Exact permutation test for cfMeDIP-seq ML classification (PARALLELIZED)"
    )
    parser.add_argument("--matrix", required=True, help="Path to feature matrix TSV")
    parser.add_argument("--annotation", required=True, help="Path to sample annotation TSV")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--feature-space", required=True, help="Name of feature space (for output files)")
    parser.add_argument("--case-group", required=True, help="Case group name (e.g., AEG)")
    parser.add_argument("--ctrl-group", required=True, help="Control group name (e.g., CTRL)")
    parser.add_argument("--top-k", type=int, default=5000, help="Number of top features to select (default: 5000)")
    parser.add_argument("--C", type=float, default=1.0, help="Regularization parameter for L2 (default: 1.0)")
    parser.add_argument("--alpha", type=float, default=0.01, help="Regularization parameter for ElasticNet (default: 0.01)")
    parser.add_argument("--model-type", choices=['L2', 'ElasticNet', 'both'], default='L2',
                        help="Model type: L2 (default), ElasticNet, or both")
    parser.add_argument("--n-jobs", type=int, default=-1, 
                        help="Number of parallel jobs (-1 = all cores, default: -1)")
    parser.add_argument("--quiet", action="store_true", help="Suppress verbose output")
    args = parser.parse_args()
    
    os.makedirs(args.output, exist_ok=True)
    
    print("="*70)
    print("EXACT PERMUTATION TEST FOR cfMeDIP-seq CLASSIFICATION")
    print(f"  PARALLELIZED with {args.n_jobs} jobs (use -1 for all {os.cpu_count()} cores)")
    print("="*70)
    
    # Load data
    print(f"\nLoading data...")
    X_all = load_matrix(args.matrix)
    ann = load_annotation(args.annotation, args.case_group, args.ctrl_group)
    
    # Align samples
    common = [s for s in ann["sample_id"] if s in X_all.index]
    if len(common) == 0:
        raise ValueError("No matching samples between matrix and annotation!")
    
    ann = ann[ann["sample_id"].isin(common)].copy()
    X = X_all.loc[common].copy()
    y = (ann.set_index("sample_id").loc[common, "group"] == args.case_group).astype(int).values
    
    print(f"  Samples: {len(common)} ({int(y.sum())} {args.case_group}, {len(y)-int(y.sum())} {args.ctrl_group})")
    print(f"  Features: {X.shape[1]}")
    print(f"  Parameters: top_k={args.top_k}, C={args.C}, model={args.model_type}")
    
    # Run exact permutation test
    if args.model_type == 'both':
        # Run both L2 and Elastic-Net
        all_results = exact_permutation_test_dual_model(
            X, y, top_k=args.top_k, C=args.C, alpha=args.alpha, n_jobs=args.n_jobs, verbose=not args.quiet
        )
        # Save results for both models
        for model_type, results in all_results.items():
            _save_results(results, args, model_type)
        return
    else:
        results = exact_permutation_test(
            X, y, top_k=args.top_k, C=args.C, model_type=args.model_type,
            alpha=args.alpha, n_jobs=args.n_jobs, verbose=not args.quiet
        )
    
    # Save results using helper
    _save_results(results, args)
    
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"  Feature space: {args.feature_space}")
    print(f"  Model: {args.model_type}")
    print(f"  Observed AUROC: {results['observed_auroc']:.4f}")
    print(f"  Exact p-value: {results['exact_p_auroc']:.4f}")
    print(f"  Significant at α=0.05: {'YES ✓' if results['exact_p_auroc'] < 0.05 else 'NO'}")
    driver_warn = results["stability_df"]["auroc_drop"].max() > 0.15
    print(f"  Single-sample driver warning: {'YES ⚠️' if driver_warn else 'NO ✓'}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
