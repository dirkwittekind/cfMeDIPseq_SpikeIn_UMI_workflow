#!/usr/bin/env python3
"""
ML Cross-Validation for cfMeDIP-seq Group Discrimination
Implements LOO + RepeatedStratifiedKFold with per-sample mean aggregation.
Based on standard workflow for AEG vs CTRL discrimination.

Features:
- Hypo/hyper encoding (split features by direction)
- K-best feature selection
- Orientation QC (check if predictions are flipped)
- Multiple CV strategies: LOO, RepeatedStratifiedKFold
- Per-sample mean aggregation for robust metrics
"""

import argparse
import os
import json
import pickle
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional

from sklearn.model_selection import (
    LeaveOneOut, RepeatedStratifiedKFold, StratifiedKFold
)
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve, precision_recall_curve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')


def load_annotation(ann_path: str, group_col: str = "group") -> pd.DataFrame:
    """Load and clean sample annotation file."""
    ann = pd.read_csv(ann_path, sep="\t", dtype=str)
    ann.columns = [c.strip() for c in ann.columns]
    ann["sample_id"] = ann["sample_id"].astype(str).str.strip()
    ann[group_col] = ann[group_col].astype(str).str.strip()
    ann[group_col] = ann[group_col].replace({"Kontrolle": "CTRL", "Control": "CTRL"})
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
    
    sample_cols = [c for c in m.columns if c in sample_ids]
    
    if len(sample_cols) >= 2:
        X = m[[feat_col] + sample_cols].set_index(feat_col)
    else:
        X = m.set_index(feat_col).T
        X = X.loc[[s for s in sample_ids if s in X.index]]
        X = X.T
    
    X = X.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    available_samples = [s for s in sample_ids if s in X.columns]
    X = X[available_samples]
    
    return X


def hypo_hyper_encode(X: pd.DataFrame, y: pd.Series, case_label: int = 1) -> pd.DataFrame:
    """
    Split features into hypo and hyper based on mean difference.
    
    For each feature:
    - If mean(case) > mean(ctrl): hyper feature
    - If mean(case) < mean(ctrl): hypo feature
    
    Creates separate columns: feature_hyper and feature_hypo
    """
    case_mask = (y == case_label).values
    ctrl_mask = ~case_mask
    
    case_mean = X.values[:, case_mask].mean(axis=1)
    ctrl_mean = X.values[:, ctrl_mask].mean(axis=1)
    
    hyper_mask = case_mean > ctrl_mean
    
    X_encoded = pd.DataFrame(index=X.columns)
    
    for i, feat in enumerate(X.index):
        if hyper_mask[i]:
            X_encoded[f"{feat}_hyper"] = X.loc[feat].values
        else:
            X_encoded[f"{feat}_hypo"] = X.loc[feat].values
    
    return X_encoded


def run_cv_pipeline(
    X: pd.DataFrame,
    y: np.ndarray,
    sample_ids: np.ndarray,
    kbest: int = 5000,
    C: float = 1.0,
    hypo_hyper: bool = True,
    seed: int = 1
) -> Dict:
    """
    Run the full CV pipeline:
    1. LOO CV
    2. RepeatedStratifiedKFold CV
    3. Per-sample mean aggregation
    
    Returns predictions and metrics for each strategy.
    """
    np.random.seed(seed)
    results = {
        "loo": {"preds": [], "probs": [], "labels": [], "samples": []},
        "rskfold": {"preds": [], "probs": [], "labels": [], "samples": [], "fold": [], "repeat": []},
    }
    
    # Feature names for coefficient tracking
    feature_coefs = []
    selected_features_all = []
    
    # Prepare data
    X_arr = X.values.astype(np.float64)
    
    # =========================================================================
    # LOO Cross-Validation
    # =========================================================================
    print("\n--- Leave-One-Out CV ---")
    loo = LeaveOneOut()
    
    for train_idx, test_idx in loo.split(X_arr):
        X_train, X_test = X_arr[train_idx], X_arr[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        
        # Feature selection on training set only
        if hypo_hyper:
            # Encode based on training set
            X_train_df = pd.DataFrame(X_train.T, index=X.columns)
            y_train_series = pd.Series(y_train)
            X_train_enc = hypo_hyper_encode(X_train_df, y_train_series, case_label=1)
            X_train = X_train_enc.values
            
            # Apply same encoding to test
            X_test_df = pd.DataFrame(X_test.T, index=X.columns)
            X_test_enc = pd.DataFrame(index=X_test_df.columns)
            for col in X_train_enc.columns:
                base_feat = col.replace("_hyper", "").replace("_hypo", "")
                if base_feat in X_test_df.index:
                    X_test_enc[col] = X_test_df.loc[base_feat].values
            X_test = X_test_enc.values
        
        # Standardize
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        # Feature selection
        k = min(kbest, X_train_scaled.shape[1])
        selector = SelectKBest(f_classif, k=k)
        X_train_sel = selector.fit_transform(X_train_scaled, y_train)
        X_test_sel = selector.transform(X_test_scaled)
        
        # Train
        model = LogisticRegression(C=C, solver="liblinear", class_weight="balanced", random_state=seed)
        model.fit(X_train_sel, y_train)
        
        # Predict
        prob = model.predict_proba(X_test_sel)[:, 1]
        pred = model.predict(X_test_sel)
        
        results["loo"]["probs"].append(prob[0])
        results["loo"]["preds"].append(pred[0])
        results["loo"]["labels"].append(y_test[0])
        results["loo"]["samples"].append(sample_ids[test_idx[0]])
    
    # =========================================================================
    # Repeated Stratified K-Fold CV
    # =========================================================================
    print("--- Repeated Stratified K-Fold CV ---")
    rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=50, random_state=seed)
    
    for repeat_fold, (train_idx, test_idx) in enumerate(rskf.split(X_arr, y)):
        repeat = repeat_fold // 5
        fold = repeat_fold % 5
        
        X_train, X_test = X_arr[train_idx], X_arr[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        
        # Feature encoding
        if hypo_hyper:
            X_train_df = pd.DataFrame(X_train.T, index=X.columns)
            y_train_series = pd.Series(y_train)
            X_train_enc = hypo_hyper_encode(X_train_df, y_train_series, case_label=1)
            X_train = X_train_enc.values
            
            X_test_df = pd.DataFrame(X_test.T, index=X.columns)
            X_test_enc = pd.DataFrame(index=X_test_df.columns)
            for col in X_train_enc.columns:
                base_feat = col.replace("_hyper", "").replace("_hypo", "")
                if base_feat in X_test_df.index:
                    X_test_enc[col] = X_test_df.loc[base_feat].values
            X_test = X_test_enc.values
        
        # Standardize
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        # Feature selection
        k = min(kbest, X_train_scaled.shape[1])
        selector = SelectKBest(f_classif, k=k)
        X_train_sel = selector.fit_transform(X_train_scaled, y_train)
        X_test_sel = selector.transform(X_test_scaled)
        
        # Train
        model = LogisticRegression(C=C, solver="liblinear", class_weight="balanced", random_state=seed)
        model.fit(X_train_sel, y_train)
        
        # Store coefficients
        if hypo_hyper:
            selected_mask = selector.get_support()
            selected_features = np.array(X_train_enc.columns)[selected_mask]
            coefs = model.coef_[0]
            feature_coefs.append(dict(zip(selected_features, coefs)))
            selected_features_all.extend(selected_features)
        
        # Predict
        probs = model.predict_proba(X_test_sel)[:, 1]
        preds = model.predict(X_test_sel)
        
        for i, idx in enumerate(test_idx):
            results["rskfold"]["probs"].append(probs[i])
            results["rskfold"]["preds"].append(preds[i])
            results["rskfold"]["labels"].append(y_test[i])
            results["rskfold"]["samples"].append(sample_ids[idx])
            results["rskfold"]["fold"].append(fold)
            results["rskfold"]["repeat"].append(repeat)
    
    # Store feature importance info
    results["feature_coefs"] = feature_coefs
    results["selected_features"] = list(set(selected_features_all))
    
    return results


def compute_metrics(results: Dict) -> Dict:
    """Compute AUROC, AUPRC for each CV strategy."""
    metrics = {}
    
    for cv_type in ["loo", "rskfold"]:
        probs = np.array(results[cv_type]["probs"])
        labels = np.array(results[cv_type]["labels"])
        
        auroc = roc_auc_score(labels, probs)
        auprc = average_precision_score(labels, probs)
        
        # Also compute flipped metrics for orientation QC
        auroc_flipped = roc_auc_score(labels, 1 - probs)
        
        metrics[cv_type] = {
            "auroc": auroc,
            "auroc_flipped": auroc_flipped,
            "auprc": auprc,
            "n_samples": len(labels),
            "n_pos": int(labels.sum()),
            "n_neg": int((1 - labels).sum()),
        }
    
    # Per-sample mean aggregation for RSKFold
    rskf_df = pd.DataFrame({
        "sample": results["rskfold"]["samples"],
        "prob": results["rskfold"]["probs"],
        "label": results["rskfold"]["labels"],
    })
    sample_means = rskf_df.groupby("sample").agg({
        "prob": "mean",
        "label": "first"
    }).reset_index()
    
    auroc_persample = roc_auc_score(sample_means["label"], sample_means["prob"])
    metrics["rskfold_persample"] = {
        "auroc": auroc_persample,
        "n_samples": len(sample_means),
    }
    
    results["sample_means"] = sample_means
    
    return metrics


def plot_curves(results: Dict, metrics: Dict, output_dir: str, prefix: str):
    """Generate ROC and PR curve plots."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # LOO curves
    loo_probs = np.array(results["loo"]["probs"])
    loo_labels = np.array(results["loo"]["labels"])
    
    fpr, tpr, _ = roc_curve(loo_labels, loo_probs)
    axes[0].plot(fpr, tpr, label=f'LOO (AUROC={metrics["loo"]["auroc"]:.3f})', color='blue')
    
    precision, recall, _ = precision_recall_curve(loo_labels, loo_probs)
    axes[1].plot(recall, precision, label=f'LOO (AUPRC={metrics["loo"]["auprc"]:.3f})', color='blue')
    
    # RSKFold per-sample mean
    sample_means = results["sample_means"]
    fpr, tpr, _ = roc_curve(sample_means["label"], sample_means["prob"])
    axes[0].plot(fpr, tpr, label=f'RSKFold PerSample (AUROC={metrics["rskfold_persample"]["auroc"]:.3f})', color='green')
    
    precision, recall, _ = precision_recall_curve(sample_means["label"], sample_means["prob"])
    axes[1].plot(recall, precision, label='RSKFold PerSample', color='green')
    
    # Formatting
    axes[0].plot([0, 1], [0, 1], 'k--', alpha=0.5)
    axes[0].set_xlabel('False Positive Rate')
    axes[0].set_ylabel('True Positive Rate')
    axes[0].set_title('ROC Curve')
    axes[0].legend(loc='lower right')
    axes[0].grid(True, alpha=0.3)
    
    axes[1].set_xlabel('Recall')
    axes[1].set_ylabel('Precision')
    axes[1].set_title('Precision-Recall Curve')
    axes[1].legend(loc='upper right')
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{prefix}.curves.png"), dpi=150)
    plt.close()


def save_predictions(results: Dict, output_dir: str, prefix: str):
    """Save prediction tables."""
    # LOO predictions
    loo_df = pd.DataFrame({
        "sample_id": results["loo"]["samples"],
        "true_label": results["loo"]["labels"],
        "predicted_prob": results["loo"]["probs"],
        "predicted_class": results["loo"]["preds"],
    })
    loo_df.to_csv(os.path.join(output_dir, f"{prefix}.preds.LOO.tsv"), sep="\t", index=False)
    
    # RSKFold predictions
    rskf_df = pd.DataFrame({
        "sample_id": results["rskfold"]["samples"],
        "true_label": results["rskfold"]["labels"],
        "predicted_prob": results["rskfold"]["probs"],
        "predicted_class": results["rskfold"]["preds"],
        "repeat": results["rskfold"]["repeat"],
        "fold": results["rskfold"]["fold"],
    })
    rskf_df.to_csv(os.path.join(output_dir, f"{prefix}.preds.RSKFold.tsv"), sep="\t", index=False)
    
    # Per-sample mean predictions
    results["sample_means"].to_csv(
        os.path.join(output_dir, f"{prefix}.preds.RSKFold_perSampleMean.tsv"),
        sep="\t", index=False
    )


def aggregate_feature_importance(results: Dict, output_dir: str, prefix: str):
    """Aggregate feature coefficients across folds."""
    if not results.get("feature_coefs"):
        return
    
    # Collect all coefficients
    all_coefs = {}
    for fold_coefs in results["feature_coefs"]:
        for feat, coef in fold_coefs.items():
            if feat not in all_coefs:
                all_coefs[feat] = []
            all_coefs[feat].append(coef)
    
    # Aggregate
    importance_data = []
    for feat, coefs in all_coefs.items():
        importance_data.append({
            "feature": feat,
            "mean_coef": np.mean(coefs),
            "mean_abs_coef": np.mean(np.abs(coefs)),
            "std_coef": np.std(coefs),
            "selection_freq": len(coefs) / len(results["feature_coefs"]),
            "sign_consistency": np.mean(np.sign(coefs) == np.sign(np.mean(coefs))),
        })
    
    imp_df = pd.DataFrame(importance_data)
    imp_df = imp_df.sort_values("mean_abs_coef", ascending=False)
    imp_df.to_csv(os.path.join(output_dir, f"{prefix}.feature_importance.tsv"), sep="\t", index=False)
    
    print(f"\nTop 10 features by importance:")
    print(imp_df.head(10)[["feature", "mean_abs_coef", "sign_consistency"]].to_string(index=False))


def main():
    parser = argparse.ArgumentParser(
        description="ML CV for cfMeDIP-seq group discrimination"
    )
    parser.add_argument("--ann", required=True, help="Sample annotation TSV")
    parser.add_argument("--matrix", required=True, help="Feature matrix TSV")
    parser.add_argument("--out-dir", required=True, help="Output directory")
    parser.add_argument("--prefix", default="ml_cv", help="Output file prefix")
    parser.add_argument("--case-group", default="AEG", help="Case group name")
    parser.add_argument("--ctrl-group", default="CTRL", help="Control group name")
    parser.add_argument("--kbest", type=int, default=5000, help="K-best features")
    parser.add_argument("--C", type=float, default=1.0, help="Regularization parameter")
    parser.add_argument("--seed", type=int, default=1, help="Random seed")
    parser.add_argument("--mode", choices=["methylOnly", "methylPlusCov"], default="methylOnly",
                       help="Feature mode")
    parser.add_argument("--covariates", nargs="+", help="Covariate columns for methylPlusCov mode")
    parser.add_argument("--no-hypo-hyper", action="store_true", help="Disable hypo/hyper encoding")
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("ML Cross-Validation for cfMeDIP-seq Discrimination")
    print("=" * 70)
    
    # Load data
    print(f"\nLoading annotation: {args.ann}")
    ann = load_annotation(args.ann)
    print(f"  Samples: {len(ann)}")
    print(f"  Groups: {ann['group'].value_counts().to_dict()}")
    
    sample_ids = ann["sample_id"].tolist()
    
    print(f"\nLoading matrix: {args.matrix}")
    X = load_matrix(args.matrix, sample_ids)
    print(f"  Features: {X.shape[0]}")
    print(f"  Samples: {X.shape[1]}")
    
    # Align annotation to matrix
    ann_aligned = ann.set_index("sample_id").loc[X.columns].reset_index()
    # Ensure sample_id column name is correct after reset
    if ann_aligned.columns[0] != "sample_id":
        ann_aligned = ann_aligned.rename(columns={ann_aligned.columns[0]: "sample_id"})
    
    # Encode labels (CTRL=0, case=1)
    y = (ann_aligned["group"] == args.case_group).astype(int).values
    sample_ids_arr = ann_aligned["sample_id"].values
    
    print(f"\nLabels: {args.ctrl_group}=0 ({(y==0).sum()}), {args.case_group}=1 ({(y==1).sum()})")
    
    # Transpose X for sklearn (samples x features)
    X_samples = X.T
    
    # Add covariates if requested
    if args.mode == "methylPlusCov" and args.covariates:
        print(f"\nAdding covariates: {args.covariates}")
        for cov in args.covariates:
            if cov in ann_aligned.columns:
                cov_vals = pd.to_numeric(ann_aligned[cov], errors="coerce").fillna(0).values
                X_samples[f"cov_{cov}"] = cov_vals
    
    # Create output directory
    os.makedirs(args.out_dir, exist_ok=True)
    
    # Run CV pipeline
    print("\n" + "=" * 70)
    print("Running CV Pipeline")
    print("=" * 70)
    results = run_cv_pipeline(
        X_samples, y, sample_ids_arr,
        kbest=args.kbest,
        C=args.C,
        hypo_hyper=not args.no_hypo_hyper,
        seed=args.seed
    )
    
    # Compute metrics
    metrics = compute_metrics(results)
    
    # Print summary
    print("\n" + "=" * 70)
    print("Results Summary")
    print("=" * 70)
    print(f"\nLOO CV:")
    print(f"  AUROC: {metrics['loo']['auroc']:.4f}")
    print(f"  AUROC (flipped): {metrics['loo']['auroc_flipped']:.4f}")
    print(f"  AUPRC: {metrics['loo']['auprc']:.4f}")
    
    print(f"\nRepeatedStratifiedKFold (5x50):")
    print(f"  AUROC (pooled): {metrics['rskfold']['auroc']:.4f}")
    print(f"  AUROC (per-sample mean): {metrics['rskfold_persample']['auroc']:.4f}")
    
    # Orientation QC
    print(f"\nOrientation QC:")
    if metrics['loo']['auroc'] > metrics['loo']['auroc_flipped']:
        print(f"  ✓ Predictions correctly oriented (raw > flipped)")
    else:
        print(f"  ⚠ Predictions may be flipped (raw < flipped)")
    
    # Save outputs
    print(f"\nSaving outputs to: {args.out_dir}")
    
    # Predictions
    save_predictions(results, args.out_dir, args.prefix)
    
    # Metrics summary
    summary_df = pd.DataFrame([
        {"metric": "auroc_loo", "value": metrics["loo"]["auroc"]},
        {"metric": "auroc_loo_flipped", "value": metrics["loo"]["auroc_flipped"]},
        {"metric": "auprc_loo", "value": metrics["loo"]["auprc"]},
        {"metric": "auroc_rskfold_pooled", "value": metrics["rskfold"]["auroc"]},
        {"metric": "auroc_rskfold_persample", "value": metrics["rskfold_persample"]["auroc"]},
        {"metric": "n_case", "value": metrics["loo"]["n_pos"]},
        {"metric": "n_ctrl", "value": metrics["loo"]["n_neg"]},
    ])
    summary_df.to_csv(os.path.join(args.out_dir, f"{args.prefix}.summary.tsv"), sep="\t", index=False)
    
    # Orientation QC
    qc_df = pd.DataFrame([{
        "auroc_raw": metrics["loo"]["auroc"],
        "auroc_flipped": metrics["loo"]["auroc_flipped"],
        "orientation_ok": metrics["loo"]["auroc"] > metrics["loo"]["auroc_flipped"],
    }])
    qc_df.to_csv(os.path.join(args.out_dir, f"{args.prefix}.orientation_qc.tsv"), sep="\t", index=False)
    
    # Plots
    plot_curves(results, metrics, args.out_dir, args.prefix)
    
    # Feature importance
    aggregate_feature_importance(results, args.out_dir, args.prefix)
    
    print("\nML CV complete!")


if __name__ == "__main__":
    main()
