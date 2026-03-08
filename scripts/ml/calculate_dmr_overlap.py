#!/usr/bin/env python3
"""
Calculate overlap between AEG DMRs and NatureCancer reference DMRs.

This uses the actual DMR results (limma moderated t-test) rather than ML coefficients.
DMRs are selected by p-value threshold and direction (hyper/hypo).
"""

import os
import subprocess
import pandas as pd
import numpy as np
import gzip
import tempfile

# Paths
NATURE_CANCER_DIR = "/home/dirk/medipipe_warp/resources/Code/15191455/pan_cancer_cfMeDIP_Figures_source_data/pan_cancer_cfMDIP_Figures_source_data/DMRs_Source_Data/BED_files"
DMR_FILE = "/data/medipipe_data/output/ml_discrimination/Outputs/AEG_vs_CTRL/DMR/windows2000_all/dmr.tsv.gz"
MATRIX_FILE = "/data/medipipe_data/output/ml_discrimination/matrices/windows/hg38_w2000/all/matrix_autosomal.tsv"

def load_coordinate_mapping():
    """Load matrix to get feature index -> coordinates mapping."""
    print("Loading coordinate mapping from matrix...")
    
    # Read just the coordinate columns
    df = pd.read_csv(MATRIX_FILE, sep="\t", usecols=["chr", "start", "end"])
    df["feature_idx"] = range(len(df))
    
    print(f"  Loaded {len(df)} regions")
    return df

def load_dmrs(coord_df, p_threshold=0.05, top_n=None):
    """Load DMR results and map to coordinates."""
    print(f"Loading DMRs (p < {p_threshold})...")
    
    # Read DMR file
    with gzip.open(DMR_FILE, 'rt') as f:
        dmr_df = pd.read_csv(f, sep="\t")
    
    print(f"  Total DMRs tested: {len(dmr_df)}")
    
    # Filter by p-value
    sig_dmrs = dmr_df[dmr_df["P.Value"] < p_threshold].copy()
    print(f"  Significant (p < {p_threshold}): {len(sig_dmrs)}")
    
    # Take top N if specified
    if top_n and len(sig_dmrs) > top_n:
        sig_dmrs = sig_dmrs.nsmallest(top_n, "P.Value")
        print(f"  Using top {top_n} by p-value")
    
    # Map feature IDs to coordinates
    sig_dmrs["feature_idx"] = sig_dmrs["feature"].astype(int)
    sig_dmrs = sig_dmrs.merge(coord_df, on="feature_idx", how="left")
    
    # Remove unmapped
    sig_dmrs = sig_dmrs.dropna(subset=["chr", "start", "end"])
    
    print(f"  With coordinates: {len(sig_dmrs)}")
    print(f"  Hyper in AEG: {(sig_dmrs['direction'] == 'hyper_in_AEG').sum()}")
    print(f"  Hypo in AEG: {(sig_dmrs['direction'] == 'hypo_in_AEG').sum()}")
    
    return sig_dmrs

def count_overlaps(query_bed, reference_bed):
    """Count overlapping regions using bedtools."""
    try:
        result = subprocess.run(
            ["bedtools", "intersect", "-a", query_bed, "-b", reference_bed, "-u"],
            capture_output=True, text=True, check=True
        )
        return len([l for l in result.stdout.strip().split("\n") if l])
    except:
        return 0

def get_reference_beds():
    """Get NatureCancer reference BED files."""
    beds = {}
    for f in os.listdir(NATURE_CANCER_DIR):
        if f.endswith(".bed") and not f.startswith("."):
            path = os.path.join(NATURE_CANCER_DIR, f)
            with open(path) as fp:
                n = sum(1 for line in fp if line.strip() and not line.startswith('#'))
            beds[f.replace(".bed", "")] = {"path": path, "n": n}
    return beds

def calculate_overlaps(dmr_df, label="All DMRs"):
    """Calculate overlaps with all NatureCancer references."""
    
    if len(dmr_df) == 0:
        return pd.DataFrame()
    
    # Create temp BED
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        bed_path = f.name
        for _, row in dmr_df.iterrows():
            f.write(f"{row['chr']}\t{int(row['start'])}\t{int(row['end'])}\n")
    
    try:
        refs = get_reference_beds()
        results = []
        
        for name, info in sorted(refs.items()):
            n_overlap = count_overlaps(bed_path, info["path"])
            pct = 100 * n_overlap / len(dmr_df) if len(dmr_df) > 0 else 0
            results.append({
                "Query": label,
                "Reference": name,
                "Ref_N": info["n"],
                "Query_N": len(dmr_df),
                "Overlap": n_overlap,
                "Pct": round(pct, 2)
            })
        
        return pd.DataFrame(results)
    finally:
        os.unlink(bed_path)

def main():
    print("=" * 80)
    print("DMR-BASED OVERLAP: AEG DMRs vs NatureCancer Reference")
    print("=" * 80)
    
    # Load coordinate mapping
    coord_df = load_coordinate_mapping()
    
    # Analyze at different p-value thresholds
    for p_thresh in [0.001, 0.01, 0.05]:
        print(f"\n{'='*60}")
        print(f"P-value threshold: {p_thresh}")
        print(f"{'='*60}")
        
        dmrs = load_dmrs(coord_df, p_threshold=p_thresh)
        
        if len(dmrs) == 0:
            print("  No significant DMRs at this threshold")
            continue
        
        # All DMRs
        print(f"\n--- All DMRs (n={len(dmrs)}) ---")
        all_results = calculate_overlaps(dmrs, f"All (n={len(dmrs)})")
        sig = all_results[all_results["Overlap"] > 0].sort_values("Pct", ascending=False)
        if len(sig) > 0:
            print(sig.head(15).to_string(index=False))
        
        # Hyper DMRs (higher methylation in AEG = tumor signal)
        hyper = dmrs[dmrs["direction"] == "hyper_in_AEG"]
        if len(hyper) > 5:
            print(f"\n--- Hypermethylated in AEG (n={len(hyper)}) ---")
            hyper_results = calculate_overlaps(hyper, f"Hyper (n={len(hyper)})")
            sig = hyper_results[hyper_results["Overlap"] > 0].sort_values("Pct", ascending=False)
            if len(sig) > 0:
                print(sig.head(10).to_string(index=False))
        
        # Hypo DMRs (lower methylation in AEG)
        hypo = dmrs[dmrs["direction"] == "hypo_in_AEG"]
        if len(hypo) > 5:
            print(f"\n--- Hypomethylated in AEG (n={len(hypo)}) ---")
            hypo_results = calculate_overlaps(hypo, f"Hypo (n={len(hypo)})")
            sig = hypo_results[hypo_results["Overlap"] > 0].sort_values("Pct", ascending=False)
            if len(sig) > 0:
                print(sig.head(10).to_string(index=False))
    
    # Save results for top 1000 DMRs
    print("\n" + "=" * 80)
    print("Saving detailed results for top 1000 DMRs...")
    
    output_dir = "/data/medipipe_data/output/ml_discrimination/Outputs/AEG_vs_CTRL/naturecancer_overlap"
    os.makedirs(output_dir, exist_ok=True)
    
    dmrs = load_dmrs(coord_df, p_threshold=1.0, top_n=1000)  # Top 1000 by p-value
    
    all_results = calculate_overlaps(dmrs, "All_Top1000")
    all_results.to_csv(os.path.join(output_dir, "DMR_top1000_all_overlap.tsv"), sep="\t", index=False)
    
    hyper = dmrs[dmrs["direction"] == "hyper_in_AEG"]
    if len(hyper) > 0:
        hyper_results = calculate_overlaps(hyper, "Hyper_Top1000")
        hyper_results.to_csv(os.path.join(output_dir, "DMR_top1000_hyper_overlap.tsv"), sep="\t", index=False)
    
    hypo = dmrs[dmrs["direction"] == "hypo_in_AEG"]
    if len(hypo) > 0:
        hypo_results = calculate_overlaps(hypo, "Hypo_Top1000")
        hypo_results.to_csv(os.path.join(output_dir, "DMR_top1000_hypo_overlap.tsv"), sep="\t", index=False)
    
    print(f"  Saved to: {output_dir}")
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)

if __name__ == "__main__":
    main()
