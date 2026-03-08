#!/usr/bin/env python3
"""
ML Visualization for cfMeDIP-seq Discrimination Analysis

Generates:
1. Volcano plots: Feature importance/coefficient vs statistical significance
2. Heatmaps: Top-K driver regions with z-scores per sample
3. Feature importance bar plots

Works with both L2 and Elastic-Net models, extracting feature coefficients
and combining with univariate statistics for comprehensive visualization.

Author: MEDIPIPE
"""

import argparse
import os
from typing import List, Tuple, Dict, Optional

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.preprocessing import StandardScaler
from collections import defaultdict
import warnings

warnings.filterwarnings('ignore')

# Model configurations (same as ml_exact_permutation.py)
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
            'l1_ratio': 0.5,
            'alpha': 0.01,
            'class_weight': 'balanced',
            'max_iter': 5000,
            'random_state': 42,
            'early_stopping': False
        }
    }
}


def load_matrix(matrix_path: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load feature matrix and return samples x features DataFrame + feature info.
    
    Returns:
        X: DataFrame with samples as rows, features as columns
        feature_info: DataFrame with feature coordinates (chr, start, end)
    """
    m = pd.read_csv(matrix_path, sep="\t")
    
    # Clean column names
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
        
        # Feature info for annotation
        feature_info = m[["feature_id", "chr", "start", "end"]].copy()
        feature_info.set_index("feature_id", inplace=True)
        
        X = m.set_index("feature_id")[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    else:
        # Assume first column is feature ID
        feature_col = m.columns[0]
        X = m.set_index(feature_col).apply(pd.to_numeric, errors="coerce").fillna(0.0)
        feature_info = pd.DataFrame(index=X.index)
        feature_info["chr"] = "unknown"
        feature_info["start"] = 0
        feature_info["end"] = 0
    
    return X.T, feature_info  # samples x features, feature_info


def load_gene_annotation(gene_bed_path: str, gene_map_path: Optional[str] = None) -> Dict[str, str]:
    """
    Load gene annotation BED file and optional gene name mapping.
    
    Returns a dictionary mapping genomic regions to gene names.
    Genomic regions are stored as "chr:start-end" keys.
    """
    # Load ENSEMBL to HGNC mapping if provided
    ensembl_to_hgnc = {}
    if gene_map_path and os.path.exists(gene_map_path):
        gene_map = pd.read_csv(gene_map_path, sep="\t")
        for _, row in gene_map.iterrows():
            ensembl_id = str(row.iloc[0]).strip()
            gene_name = str(row.iloc[1]).strip()
            if ensembl_id and gene_name and gene_name != 'nan':
                ensembl_to_hgnc[ensembl_id] = gene_name
        print(f"  Loaded {len(ensembl_to_hgnc)} gene name mappings")
    
    # Load BED file and build interval lookup
    region_to_gene = {}
    gene_intervals = defaultdict(list)  # chr -> list of (start, end, gene_name)
    
    if os.path.exists(gene_bed_path):
        with open(gene_bed_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    gene_id = parts[3]
                    
                    # Convert ENSEMBL ID to HGNC if possible
                    if gene_id.startswith('ENSG') and gene_id in ensembl_to_hgnc:
                        gene_name = ensembl_to_hgnc[gene_id]
                    else:
                        gene_name = gene_id
                    
                    gene_intervals[chrom].append((start, end, gene_name))
        
        # Sort intervals for efficient lookup
        for chrom in gene_intervals:
            gene_intervals[chrom].sort(key=lambda x: x[0])
        
        print(f"  Loaded {sum(len(v) for v in gene_intervals.values())} gene regions")
    
    return gene_intervals


def find_gene_for_region(chrom: str, start: int, end: int, 
                          gene_intervals: Dict[str, list]) -> Optional[str]:
    """
    Find the gene that overlaps with a given genomic region.
    Returns the gene name or None if no overlap found.
    """
    if chrom not in gene_intervals:
        return None
    
    # Find overlapping genes
    overlapping = []
    for g_start, g_end, g_name in gene_intervals[chrom]:
        # Check for overlap
        if g_start < end and g_end > start:
            overlap_size = min(end, g_end) - max(start, g_start)
            overlapping.append((g_name, overlap_size))
    
    if overlapping:
        # Return the gene with the largest overlap
        overlapping.sort(key=lambda x: -x[1])
        return overlapping[0][0]
    
    return None


def get_feature_gene_name(feature_id: str, feature_info: pd.DataFrame, 
                          gene_intervals: Optional[Dict[str, list]] = None) -> str:
    """
    Get gene name for a feature, falling back to coordinates if no gene found.
    """
    if feature_id in feature_info.index and gene_intervals:
        chrom = str(feature_info.loc[feature_id, "chr"])
        start = int(feature_info.loc[feature_id, "start"])
        end = int(feature_info.loc[feature_id, "end"])
        
        gene_name = find_gene_for_region(chrom, start, end, gene_intervals)
        if gene_name:
            return gene_name
    
    # Fallback to abbreviated coordinate format
    if feature_id in feature_info.index:
        chrom = str(feature_info.loc[feature_id, "chr"])
        start = int(feature_info.loc[feature_id, "start"])
        return f"{chrom}:{start//1000}k"
    
    return feature_id[:15]


def load_annotation(annotation_path: str, case_group: str, ctrl_group: str) -> pd.DataFrame:
    """Load sample annotation and filter to case/control groups."""
    ann = pd.read_csv(annotation_path, sep="\t", dtype=str).fillna("")
    ann.columns = [c.strip() for c in ann.columns]
    ann["sample_id"] = ann["sample_id"].astype(str).str.strip()
    ann["group"] = ann["group"].astype(str).str.strip()
    ann = ann[ann["group"].isin([case_group, ctrl_group])].copy()
    return ann


def compute_univariate_stats(X: pd.DataFrame, y: np.ndarray) -> pd.DataFrame:
    """
    Compute univariate statistics for each feature.
    
    Returns DataFrame with:
    - F-score (ANOVA)
    - p-value
    - -log10(p-value)
    - Mean difference (case - control)
    - Fold change
    """
    n_features = X.shape[1]
    feature_names = X.columns.tolist()
    
    # ANOVA F-test
    f_scores, p_values = f_classif(X, y)
    
    # Mean difference
    case_means = X.iloc[y == 1].mean(axis=0)
    ctrl_means = X.iloc[y == 0].mean(axis=0)
    mean_diff = case_means - ctrl_means
    
    # Log fold change (add small constant for stability)
    epsilon = 1e-6
    fold_change = np.log2((case_means + epsilon) / (ctrl_means + epsilon))
    
    stats_df = pd.DataFrame({
        "feature_id": feature_names,
        "f_score": f_scores,
        "p_value": p_values,
        "neg_log10_p": -np.log10(p_values + 1e-300),
        "mean_diff": mean_diff.values,
        "log2_fc": fold_change.values,
        "case_mean": case_means.values,
        "ctrl_mean": ctrl_means.values
    })
    
    # Multiple testing correction (Benjamini-Hochberg)
    stats_df = stats_df.sort_values("p_value")
    n = len(stats_df)
    stats_df["rank"] = range(1, n + 1)
    stats_df["fdr"] = stats_df["p_value"] * n / stats_df["rank"]
    stats_df["fdr"] = stats_df["fdr"].clip(upper=1.0)
    # Make FDR monotonic
    stats_df["fdr"] = stats_df["fdr"][::-1].cummin()[::-1]
    
    return stats_df.set_index("feature_id")


def fit_model_and_get_coefficients(X: pd.DataFrame, y: np.ndarray, 
                                    top_k: int, C: float, alpha: float,
                                    model_type: str = 'L2') -> pd.DataFrame:
    """
    Fit model on full data and extract feature coefficients.
    
    Returns DataFrame with:
    - coefficient (from model)
    - abs_coefficient
    - selected (was feature in top-k?)
    """
    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Select top-k features
    k = min(top_k, X_scaled.shape[1])
    sel = SelectKBest(score_func=f_classif, k=k)
    X_sel = sel.fit_transform(X_scaled, y)
    
    # Get selected feature mask
    selected_mask = sel.get_support()
    selected_features = X.columns[selected_mask].tolist()
    
    # Fit classifier
    config = MODEL_CONFIGS[model_type]
    clf = config['class'](**config['params'])
    
    if model_type == 'L2':
        clf.C = C
    else:  # ElasticNet
        clf.alpha = alpha
    
    clf.fit(X_sel, y)
    
    # Extract coefficients
    coef = clf.coef_.flatten()
    
    # Map back to all features
    coef_df = pd.DataFrame({
        "feature_id": X.columns,
        "coefficient": 0.0,
        "abs_coefficient": 0.0,
        "selected": selected_mask
    })
    coef_df.set_index("feature_id", inplace=True)
    
    # Fill in coefficients for selected features
    for i, feat in enumerate(selected_features):
        if i < len(coef):
            coef_df.loc[feat, "coefficient"] = coef[i]
            coef_df.loc[feat, "abs_coefficient"] = abs(coef[i])
    
    return coef_df


def create_volcano_plot(stats_df: pd.DataFrame, coef_df: pd.DataFrame,
                        feature_info: pd.DataFrame, gene_intervals: Optional[Dict[str, list]],
                        output_path: str, feature_space: str, model_type: str,
                        case_group: str, ctrl_group: str,
                        top_n_label: int = 10, fdr_threshold: float = 0.05):
    """
    Create volcano plot: log2(fold change) vs -log10(p-value).
    Points colored by model coefficient direction, labeled with gene names.
    """
    # Merge statistics with coefficients
    merged = stats_df.join(coef_df, how="left")
    merged = merged[merged["selected"] == True].copy()
    
    if len(merged) == 0:
        print(f"  Warning: No selected features for volcano plot")
        return
    
    # Categorize features
    merged["category"] = "Not significant"
    merged.loc[(merged["coefficient"] > 0) & (merged["p_value"] < 0.05), "category"] = f"{case_group}-enriched"
    merged.loc[(merged["coefficient"] < 0) & (merged["p_value"] < 0.05), "category"] = f"{ctrl_group}-enriched"
    
    # Color palette
    color_map = {
        f"{case_group}-enriched": "#D62728",  # Red
        f"{ctrl_group}-enriched": "#1F77B4",  # Blue  
        "Not significant": "#999999"  # Gray
    }
    
    fig, ax = plt.subplots(figsize=(12, 9))
    
    # Plot each category separately for legend
    for category, color in color_map.items():
        mask = merged["category"] == category
        if mask.sum() == 0:
            continue
        subset = merged[mask]
        
        # Size by coefficient magnitude
        max_coef = merged["abs_coefficient"].max()
        if max_coef > 0:
            sizes = 30 + 120 * (subset["abs_coefficient"] / max_coef).fillna(0)
        else:
            sizes = 50
        
        ax.scatter(
            subset["log2_fc"],
            subset["neg_log10_p"],
            c=color,
            s=sizes,
            alpha=0.7,
            edgecolors="white",
            linewidths=0.5,
            label=f"{category} (n={mask.sum()})"
        )
    
    # Add significance threshold lines
    ax.axhline(-np.log10(0.05), color="#666666", linestyle="--", linewidth=1.5, alpha=0.7)
    ax.text(ax.get_xlim()[1] * 0.95, -np.log10(0.05) + 0.1, "p = 0.05", 
            fontsize=10, ha="right", va="bottom", color="#666666")
    
    if fdr_threshold < 1.0:
        fdr_line_y = merged.loc[merged["fdr"] < fdr_threshold, "neg_log10_p"].min() if (merged["fdr"] < fdr_threshold).any() else None
        if fdr_line_y:
            ax.axhline(fdr_line_y, color="#FFA500", linestyle=":", linewidth=1.5, alpha=0.7)
            ax.text(ax.get_xlim()[1] * 0.95, fdr_line_y + 0.1, f"FDR = {fdr_threshold}", 
                    fontsize=10, ha="right", va="bottom", color="#FFA500")
    
    ax.axvline(0, color="#CCCCCC", linestyle="-", linewidth=1, alpha=0.5)
    
    # Label top features with gene names
    top_features = merged.nlargest(top_n_label, "abs_coefficient")
    texts = []
    for idx, row in top_features.iterrows():
        gene_name = get_feature_gene_name(idx, feature_info, gene_intervals)
        
        # Color label by direction
        label_color = "#D62728" if row["coefficient"] > 0 else "#1F77B4"
        
        text = ax.annotate(
            gene_name,
            xy=(row["log2_fc"], row["neg_log10_p"]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=10,
            fontweight="bold",
            color=label_color,
            alpha=0.9,
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.7)
        )
        texts.append(text)
    
    # Styling
    ax.set_xlabel(f"Log₂ Fold Change ({case_group} vs {ctrl_group})", fontsize=14, fontweight="bold")
    ax.set_ylabel("-Log₁₀(p-value)", fontsize=14, fontweight="bold")
    ax.set_title(f"{feature_space} - {model_type} Model\nVolcano Plot of Discriminative Regions", 
                 fontsize=16, fontweight="bold", pad=15)
    
    # Legend
    legend = ax.legend(loc="upper right", fontsize=11, framealpha=0.9, 
                       title="Classification", title_fontsize=12)
    
    # Add annotation for point size
    ax.text(0.02, 0.98, "Point size = |coefficient|", 
            transform=ax.transAxes, fontsize=9, va="top", 
            style="italic", color="#666666")
    
    ax.tick_params(axis="both", labelsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(output_path.replace(".pdf", ".pdf"), dpi=150, bbox_inches="tight")
    plt.savefig(output_path.replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    plt.close()
    
    print(f"  Saved volcano plot: {output_path}")


def create_driver_heatmap(X: pd.DataFrame, y: np.ndarray, stats_df: pd.DataFrame,
                          coef_df: pd.DataFrame, feature_info: pd.DataFrame,
                          gene_intervals: Optional[Dict[str, list]],
                          output_path: str, feature_space: str, model_type: str,
                          case_group: str, ctrl_group: str,
                          sample_ids: List[str], top_n: int = 20):
    """
    Create heatmap of top-N driver regions with z-scores per sample.
    Uses gene names for labels when available.
    """
    # Merge statistics with coefficients
    merged = stats_df.join(coef_df, how="left")
    merged = merged[merged["selected"] == True].copy()
    
    if len(merged) == 0:
        print(f"  Warning: No selected features for heatmap")
        return
    
    # Get top-N features by absolute coefficient
    top_features = merged.nlargest(top_n, "abs_coefficient").index.tolist()
    
    # Extract data for these features
    X_top = X[top_features].copy()
    
    # Compute z-scores per feature
    X_zscore = (X_top - X_top.mean()) / X_top.std()
    X_zscore = X_zscore.fillna(0)
    
    # Create feature labels with gene names
    feature_labels = []
    gene_names_for_table = []
    for feat in top_features:
        coef = coef_df.loc[feat, "coefficient"]
        direction = "↑" if coef > 0 else "↓"
        
        # Get gene name
        gene_name = get_feature_gene_name(feat, feature_info, gene_intervals)
        gene_names_for_table.append(gene_name)
        
        # Create label with direction indicator
        label = f"{gene_name} {direction}"
        feature_labels.append(label)
    
    # Create sample labels with group info
    sample_labels = []
    for i, sid in enumerate(sample_ids):
        group = case_group if y[i] == 1 else ctrl_group
        sample_labels.append(f"{sid} ({group})")
    
    # Sort samples by group
    sort_idx = np.argsort(y)[::-1]  # Cases first
    X_zscore_sorted = X_zscore.iloc[sort_idx]
    sample_labels_sorted = [sample_labels[i] for i in sort_idx]
    y_sorted = y[sort_idx]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(max(14, len(sample_labels) * 0.9), max(10, top_n * 0.55)))
    
    # Create heatmap with improved styling
    heatmap = sns.heatmap(
        X_zscore_sorted.T.values,
        ax=ax,
        cmap="RdBu_r",
        center=0,
        vmin=-3, vmax=3,
        xticklabels=sample_labels_sorted,
        yticklabels=feature_labels,
        cbar_kws={"label": "Z-score", "shrink": 0.7},
        linewidths=0.5,
        linecolor="white"
    )
    
    # Style the colorbar
    cbar = heatmap.collections[0].colorbar
    cbar.ax.tick_params(labelsize=10)
    cbar.ax.set_ylabel("Z-score", fontsize=12, fontweight="bold")
    
    # Add group separator
    n_case = int(y.sum())
    ax.axvline(n_case, color="black", linewidth=3)
    
    # Color y-axis labels by direction
    for i, label in enumerate(ax.get_yticklabels()):
        if "↑" in label.get_text():
            label.set_color("#D62728")  # Red for case-enriched
        else:
            label.set_color("#1F77B4")  # Blue for control-enriched
        label.set_fontweight("bold")
    
    ax.set_xlabel("Samples", fontsize=14, fontweight="bold")
    ax.set_ylabel(f"Top {top_n} Driver Genes", fontsize=14, fontweight="bold")
    ax.set_title(f"{feature_space} - {model_type} Model\n"
                 f"Top {top_n} Discriminative Regions (Z-score normalized)", 
                 fontsize=16, fontweight="bold", pad=15)
    
    # Add group labels at top
    ax.text(n_case / 2, -0.5, case_group, fontsize=12, fontweight="bold", 
            ha="center", va="bottom", color="#D62728")
    ax.text(n_case + (len(y) - n_case) / 2, -0.5, ctrl_group, fontsize=12, fontweight="bold", 
            ha="center", va="bottom", color="#1F77B4")
    
    plt.xticks(rotation=45, ha="right", fontsize=10)
    plt.yticks(fontsize=11)
    
    plt.tight_layout()
    plt.savefig(output_path.replace(".pdf", ".pdf"), dpi=150, bbox_inches="tight")
    plt.savefig(output_path.replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    plt.close()
    
    print(f"  Saved driver heatmap: {output_path}")
    
    # Save the top features table with gene names
    top_table = merged.loc[top_features].copy()
    top_table["feature_id"] = top_features
    top_table["gene_name"] = gene_names_for_table
    if len(feature_info) > 0:
        top_table = top_table.join(feature_info[["chr", "start", "end"]], how="left")
    
    # Reorder columns for clarity
    first_cols = ["gene_name", "chr", "start", "end", "coefficient", "p_value", "fdr", "log2_fc"]
    other_cols = [c for c in top_table.columns if c not in first_cols]
    top_table = top_table[[c for c in first_cols if c in top_table.columns] + other_cols]
    
    table_path = output_path.replace(".pdf", ".tsv")
    top_table.to_csv(table_path, sep="\t")
    print(f"  Saved top features table: {table_path}")


def create_feature_importance_barplot(stats_df: pd.DataFrame, coef_df: pd.DataFrame,
                                       feature_info: pd.DataFrame,
                                       gene_intervals: Optional[Dict[str, list]],
                                       output_path: str, feature_space: str, 
                                       model_type: str, case_group: str, ctrl_group: str,
                                       top_n: int = 30):
    """
    Create horizontal bar plot of feature importance (coefficient magnitude).
    Uses gene names for labels when available.
    """
    # Merge statistics with coefficients
    merged = stats_df.join(coef_df, how="left")
    merged = merged[merged["selected"] == True].copy()
    
    if len(merged) == 0:
        print(f"  Warning: No selected features for bar plot")
        return
    
    # Get top-N features by absolute coefficient
    top_features = merged.nlargest(top_n, "abs_coefficient")
    
    # Create labels with gene names
    labels = []
    for feat in top_features.index:
        gene_name = get_feature_gene_name(feat, feature_info, gene_intervals)
        labels.append(gene_name)
    
    fig, ax = plt.subplots(figsize=(12, max(9, top_n * 0.4)))
    
    # Colors based on coefficient sign
    colors = ["#D62728" if c > 0 else "#1F77B4" for c in top_features["coefficient"]]
    
    y_pos = range(len(top_features))
    bars = ax.barh(y_pos, top_features["abs_coefficient"], color=colors, alpha=0.85, 
                   edgecolor="white", linewidth=0.5)
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=11, fontweight="bold")
    
    # Color y-axis labels by direction
    for i, (label, coef) in enumerate(zip(ax.get_yticklabels(), top_features["coefficient"])):
        label.set_color("#D62728" if coef > 0 else "#1F77B4")
    
    ax.invert_yaxis()  # Highest at top
    
    ax.set_xlabel("Absolute Model Coefficient", fontsize=14, fontweight="bold")
    ax.set_title(f"{feature_space} - {model_type} Model\nTop {top_n} Discriminative Features", 
                 fontsize=16, fontweight="bold", pad=15)
    
    # Add legend
    case_patch = mpatches.Patch(color="#D62728", label=f"{case_group}-enriched")
    ctrl_patch = mpatches.Patch(color="#1F77B4", label=f"{ctrl_group}-enriched")
    ax.legend(handles=[case_patch, ctrl_patch], loc="lower right", fontsize=11)
    
    ax.tick_params(axis="x", labelsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(output_path.replace(".pdf", ".pdf"), dpi=150, bbox_inches="tight")
    plt.savefig(output_path.replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    plt.close()
    
    print(f"  Saved feature importance plot: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate ML visualizations for cfMeDIP-seq discrimination analysis"
    )
    parser.add_argument("--matrix", required=True, help="Path to feature matrix TSV")
    parser.add_argument("--annotation", required=True, help="Path to sample annotation TSV")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--feature-space", required=True, help="Name of feature space")
    parser.add_argument("--case-group", required=True, help="Case group name (e.g., AEG)")
    parser.add_argument("--ctrl-group", required=True, help="Control group name (e.g., CTRL)")
    parser.add_argument("--top-k", type=int, default=5000, help="Number of features to select (default: 5000)")
    parser.add_argument("--C", type=float, default=1.0, help="L2 regularization parameter (default: 1.0)")
    parser.add_argument("--alpha", type=float, default=0.01, help="ElasticNet alpha parameter (default: 0.01)")
    parser.add_argument("--model-type", choices=['L2', 'ElasticNet', 'both'], default='both',
                        help="Model type for visualization (default: both)")
    parser.add_argument("--top-n-heatmap", type=int, default=20, help="Top features for heatmap (default: 20)")
    parser.add_argument("--top-n-barplot", type=int, default=30, help="Top features for bar plot (default: 30)")
    parser.add_argument("--gene-annotation", default=None, help="Path to gene annotation BED file (for gene names)")
    parser.add_argument("--gene-map", default=None, help="Path to ENSEMBL->HGNC gene name mapping TSV")
    args = parser.parse_args()
    
    os.makedirs(args.output, exist_ok=True)
    
    print("=" * 70)
    print("ML VISUALIZATION FOR cfMeDIP-seq DISCRIMINATION")
    print("=" * 70)
    
    # Load data
    print("\nLoading data...")
    X_all, feature_info = load_matrix(args.matrix)
    ann = load_annotation(args.annotation, args.case_group, args.ctrl_group)
    
    # Align samples - handle partial matching (e.g., "148" vs "148_S6_L004")
    matrix_samples = X_all.index.tolist()
    ann_samples = ann["sample_id"].tolist()
    
    # Try exact match first
    common = [s for s in ann_samples if s in matrix_samples]
    
    # If no exact matches, try partial matching
    if len(common) == 0:
        sample_mapping = {}  # matrix_sample -> annotation_sample
        for m_sample in matrix_samples:
            for a_sample in ann_samples:
                # Check if matrix sample is prefix of annotation sample
                if a_sample.startswith(m_sample + "_") or a_sample == m_sample:
                    sample_mapping[m_sample] = a_sample
                    break
                # Check if annotation sample is prefix of matrix sample
                elif m_sample.startswith(a_sample + "_") or m_sample.startswith(a_sample + "."):
                    sample_mapping[m_sample] = a_sample
                    break
        
        if sample_mapping:
            # Rename matrix index to match annotation
            X_all = X_all.rename(index={m: a for m, a in sample_mapping.items()})
            common = [s for s in ann_samples if s in X_all.index]
            print(f"  Matched {len(sample_mapping)} samples using partial matching")
    
    if len(common) == 0:
        print(f"  Matrix samples: {matrix_samples[:5]}...")
        print(f"  Annotation samples: {ann_samples[:5]}...")
        raise ValueError("No matching samples between matrix and annotation!")
    
    ann = ann[ann["sample_id"].isin(common)].copy()
    X = X_all.loc[common].copy()
    y = (ann.set_index("sample_id").loc[common, "group"] == args.case_group).astype(int).values
    sample_ids = common
    
    print(f"  Samples: {len(common)} ({int(y.sum())} {args.case_group}, {len(y)-int(y.sum())} {args.ctrl_group})")
    print(f"  Features: {X.shape[1]}")
    
    # Load gene annotation if provided
    gene_intervals = None
    if args.gene_annotation:
        print(f"\nLoading gene annotation...")
        gene_intervals = load_gene_annotation(args.gene_annotation, args.gene_map)
    else:
        # Try default paths - load both promoters and gene bodies for better coverage
        default_gene_map = "/home/dirk/medipipe_warp/resources/annotations/gene_map/gene_map.gencode_v46_basic.tsv"
        gene_beds = [
            "/home/dirk/medipipe_warp/resources/annotations/hg38_promoters_2kb.bed",
            "/home/dirk/medipipe_warp/resources/annotations/hg38_gencode_v49_genebody.bed"
        ]
        
        gene_intervals = defaultdict(list)
        for bed_path in gene_beds:
            if os.path.exists(bed_path):
                print(f"\nLoading gene annotation from {os.path.basename(bed_path)}...")
                intervals = load_gene_annotation(bed_path, default_gene_map if os.path.exists(default_gene_map) else None)
                # Merge intervals
                for chrom in intervals:
                    gene_intervals[chrom].extend(intervals[chrom])
        
        # Sort merged intervals
        for chrom in gene_intervals:
            gene_intervals[chrom].sort(key=lambda x: x[0])
        
        if not gene_intervals:
            gene_intervals = None
    
    # Compute univariate statistics
    print("\nComputing univariate statistics...")
    stats_df = compute_univariate_stats(X, y)
    
    # Save univariate stats
    stats_path = os.path.join(args.output, f"{args.feature_space}.univariate_stats.tsv")
    stats_df.to_csv(stats_path, sep="\t")
    print(f"  Saved: {stats_path}")
    
    # Process each model type
    model_types = ['L2', 'ElasticNet'] if args.model_type == 'both' else [args.model_type]
    
    for model_type in model_types:
        print(f"\n{'='*50}")
        print(f"Model: {model_type}")
        print(f"{'='*50}")
        
        # Fit model and get coefficients
        print(f"\nFitting {model_type} model and extracting coefficients...")
        coef_df = fit_model_and_get_coefficients(
            X, y, top_k=args.top_k, C=args.C, alpha=args.alpha, model_type=model_type
        )
        
        # Save coefficients
        coef_path = os.path.join(args.output, f"{args.feature_space}.{model_type}.coefficients.tsv")
        coef_df.to_csv(coef_path, sep="\t")
        print(f"  Saved: {coef_path}")
        
        # Create volcano plot
        print(f"\nCreating volcano plot...")
        volcano_path = os.path.join(args.output, f"{args.feature_space}.{model_type}.volcano.pdf")
        create_volcano_plot(
            stats_df, coef_df, feature_info, gene_intervals,
            volcano_path, args.feature_space, model_type,
            args.case_group, args.ctrl_group
        )
        
        # Create driver heatmap
        print(f"\nCreating driver heatmap (top {args.top_n_heatmap})...")
        heatmap_path = os.path.join(args.output, f"{args.feature_space}.{model_type}.driver_heatmap.pdf")
        create_driver_heatmap(
            X, y, stats_df, coef_df, feature_info, gene_intervals,
            heatmap_path, args.feature_space, model_type,
            args.case_group, args.ctrl_group, sample_ids,
            top_n=args.top_n_heatmap
        )
        
        # Create feature importance bar plot
        print(f"\nCreating feature importance plot (top {args.top_n_barplot})...")
        barplot_path = os.path.join(args.output, f"{args.feature_space}.{model_type}.feature_importance.pdf")
        create_feature_importance_barplot(
            stats_df, coef_df, feature_info, gene_intervals,
            barplot_path, args.feature_space, model_type,
            args.case_group, args.ctrl_group,
            top_n=args.top_n_barplot
        )
    
    print(f"\n{'='*70}")
    print("VISUALIZATION COMPLETE")
    print(f"Output directory: {args.output}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
