#!/usr/bin/env python3
"""
Calculate overlap between AEG driver regions and NatureCancer reference DMRs.

Uses bedtools intersect to find overlapping regions and reports:
1. Number of AEG regions overlapping each cancer type
2. Percentage overlap
3. Enrichment statistics
"""

import os
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
import tempfile

# Paths
NATURE_CANCER_DIR = "/home/dirk/medipipe_warp/resources/Code/15191455/pan_cancer_cfMeDIP_Figures_source_data/pan_cancer_cfMDIP_Figures_source_data/DMRs_Source_Data/BED_files"
AEG_RESULTS_DIR = "/data/medipipe_data/output/ml_discrimination/Outputs/AEG_vs_CTRL/visualizations"

def load_aeg_drivers(model_type="L2", top_n=100):
    """Load top AEG driver regions from ML analysis."""
    coef_file = os.path.join(AEG_RESULTS_DIR, f"autosomal_windows.{model_type}.coefficients.tsv")
    
    if not os.path.exists(coef_file):
        raise FileNotFoundError(f"Coefficients file not found: {coef_file}")
    
    df = pd.read_csv(coef_file, sep="\t")
    
    # Filter to selected features and sort by absolute coefficient
    df = df[df["selected"] == True].copy()
    df = df.nlargest(top_n, "abs_coefficient")
    
    # Parse feature_id to get chr, start, end
    regions = []
    for feat_id in df["feature_id"]:
        try:
            parts = feat_id.split(":")
            chrom = parts[0]
            coords = parts[1].split("-")
            start = int(coords[0])
            end = int(coords[1])
            regions.append({
                "chr": chrom,
                "start": start,
                "end": end,
                "feature_id": feat_id,
                "coefficient": df[df["feature_id"] == feat_id]["coefficient"].values[0]
            })
        except:
            continue
    
    return pd.DataFrame(regions)

def create_bed_from_df(df, output_path):
    """Create BED file from DataFrame."""
    bed_df = df[["chr", "start", "end"]].copy()
    bed_df.to_csv(output_path, sep="\t", header=False, index=False)
    return output_path

def count_overlaps(query_bed, reference_bed):
    """Count number of query regions overlapping with reference using bedtools."""
    try:
        result = subprocess.run(
            ["bedtools", "intersect", "-a", query_bed, "-b", reference_bed, "-u"],
            capture_output=True,
            text=True,
            check=True
        )
        # Count lines in output
        overlapping = len([l for l in result.stdout.strip().split("\n") if l])
        return overlapping
    except subprocess.CalledProcessError as e:
        print(f"  Warning: bedtools error - {e.stderr}")
        return 0
    except FileNotFoundError:
        print("  Error: bedtools not found. Please install bedtools.")
        return 0

def get_reference_bed_files():
    """Get all cancer-specific BED files from NatureCancer."""
    bed_files = {}
    
    for f in os.listdir(NATURE_CANCER_DIR):
        if f.endswith(".bed") and not f.startswith("."):
            # Parse cancer type from filename
            name = f.replace(".bed", "")
            full_path = os.path.join(NATURE_CANCER_DIR, f)
            
            # Count regions in file
            with open(full_path, 'r') as fp:
                n_regions = sum(1 for line in fp if line.strip() and not line.startswith('#'))
            
            bed_files[name] = {
                "path": full_path,
                "n_regions": n_regions
            }
    
    return bed_files

def calculate_all_overlaps(aeg_df, top_n):
    """Calculate overlaps with all NatureCancer reference sets."""
    
    # Create temporary BED file for AEG regions
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        aeg_bed = f.name
        for _, row in aeg_df.iterrows():
            f.write(f"{row['chr']}\t{row['start']}\t{row['end']}\n")
    
    try:
        reference_beds = get_reference_bed_files()
        
        results = []
        for name, info in sorted(reference_beds.items()):
            n_overlap = count_overlaps(aeg_bed, info["path"])
            pct_overlap = 100 * n_overlap / len(aeg_df) if len(aeg_df) > 0 else 0
            
            results.append({
                "Reference": name,
                "Reference_Regions": info["n_regions"],
                "AEG_Overlapping": n_overlap,
                "AEG_Total": len(aeg_df),
                "Overlap_Pct": round(pct_overlap, 2)
            })
        
        return pd.DataFrame(results)
    
    finally:
        os.unlink(aeg_bed)

def main():
    print("=" * 80)
    print("OVERLAP ANALYSIS: AEG Driver Regions vs NatureCancer Reference DMRs")
    print("=" * 80)
    
    # Analyze both models
    for model in ["L2", "ElasticNet"]:
        print(f"\n{'='*60}")
        print(f"Model: {model}")
        print(f"{'='*60}")
        
        for top_n in [50, 100, 500, 1000, 5000]:
            print(f"\n--- Top {top_n} Driver Regions ---")
            
            try:
                aeg_df = load_aeg_drivers(model, top_n)
                print(f"Loaded {len(aeg_df)} AEG driver regions")
                
                # Separate hyper (positive coef) and hypo (negative coef)
                hyper_df = aeg_df[aeg_df["coefficient"] > 0]
                hypo_df = aeg_df[aeg_df["coefficient"] < 0]
                
                print(f"  Hypermethylated (AEG-enriched): {len(hyper_df)}")
                print(f"  Hypomethylated (CTRL-enriched): {len(hypo_df)}")
                
                # Calculate overlaps for all regions
                print("\nAll AEG Driver Regions:")
                results = calculate_all_overlaps(aeg_df, top_n)
                
                # Filter to show meaningful overlaps
                sig_results = results[results["AEG_Overlapping"] > 0].sort_values("Overlap_Pct", ascending=False)
                
                if len(sig_results) > 0:
                    print(sig_results.to_string(index=False))
                else:
                    print("  No overlaps found")
                
                # Show hyper-specific overlaps
                if len(hyper_df) > 5:
                    print(f"\nHypermethylated Regions Only (n={len(hyper_df)}):")
                    hyper_results = calculate_all_overlaps(hyper_df, top_n)
                    hyper_sig = hyper_results[hyper_results["AEG_Overlapping"] > 0].sort_values("Overlap_Pct", ascending=False).head(10)
                    if len(hyper_sig) > 0:
                        print(hyper_sig.to_string(index=False))
                
                # Show hypo-specific overlaps
                if len(hypo_df) > 5:
                    print(f"\nHypomethylated Regions Only (n={len(hypo_df)}):")
                    hypo_results = calculate_all_overlaps(hypo_df, top_n)
                    hypo_sig = hypo_results[hypo_results["AEG_Overlapping"] > 0].sort_values("Overlap_Pct", ascending=False).head(10)
                    if len(hypo_sig) > 0:
                        print(hypo_sig.to_string(index=False))
                        
            except Exception as e:
                print(f"  Error: {e}")
    
    # Save detailed results
    print("\n" + "=" * 80)
    print("Saving detailed results...")
    
    output_dir = "/data/medipipe_data/output/ml_discrimination/Outputs/AEG_vs_CTRL/naturecancer_overlap"
    os.makedirs(output_dir, exist_ok=True)
    
    for model in ["L2", "ElasticNet"]:
        for top_n in [100, 500, 5000]:
            try:
                aeg_df = load_aeg_drivers(model, top_n)
                results = calculate_all_overlaps(aeg_df, top_n)
                
                output_file = os.path.join(output_dir, f"{model}_top{top_n}_overlap.tsv")
                results.to_csv(output_file, sep="\t", index=False)
                print(f"  Saved: {output_file}")
            except Exception as e:
                print(f"  Error for {model} top{top_n}: {e}")
    
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)

if __name__ == "__main__":
    main()
