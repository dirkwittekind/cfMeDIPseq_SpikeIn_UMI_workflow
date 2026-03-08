#!/usr/bin/env python3
"""
Run AlphaGenome analysis on intergenic/unannotated driver regions.
This script prepares the input BED file and runs the AlphaGenome workflow.

Usage:
    python3 run_alphagenome_drivers.py --drivers DRIVERS_TSV --output-dir OUTPUT_DIR
    
Example:
    python3 run_alphagenome_drivers.py \
        --drivers /data/medipipe_data/output/ml_discrimination/Outputs/AEG_vs_CTRL/visualizations/top10_L2_drivers_table.tsv \
        --output-dir /data/medipipe_data/output/ml_discrimination/Outputs/AEG_vs_CTRL/alphagenome

Author: MEDIPIPE Team
"""

import argparse
import os
import sys
import subprocess
import pandas as pd
import yaml
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="Run AlphaGenome on driver regions")
    parser.add_argument("--drivers", required=True, help="TSV file with driver regions")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--filter-intergenic", action="store_true", 
                        help="Only analyze regions without gene annotation")
    parser.add_argument("--slop", type=int, default=2000,
                        help="Extend regions by this many bp on each side (default: 2000)")
    parser.add_argument("--top-n", type=int, default=10,
                        help="Maximum number of regions to analyze (default: 10)")
    parser.add_argument("--api-key-env", default="ALPHAGENOME_API_KEY",
                        help="Environment variable containing AlphaGenome API key")
    parser.add_argument("--dry-run", action="store_true",
                        help="Only prepare files, don't run Snakemake")
    return parser.parse_args()


def prepare_driver_bed(drivers_path: str, output_dir: str, 
                       filter_intergenic: bool, slop: int, top_n: int) -> str:
    """Prepare BED file from driver regions."""
    df = pd.read_csv(drivers_path, sep="\t")
    
    # Expected columns: Gene, Region, Chr, Start, End, ...
    required_cols = ["Region", "Chr", "Start", "End"]
    for col in required_cols:
        if col not in df.columns:
            # Try lowercase
            col_lower = col.lower()
            if col_lower in df.columns:
                df = df.rename(columns={col_lower: col})
            else:
                raise ValueError(f"Missing required column: {col}")
    
    # Filter intergenic regions if requested
    if filter_intergenic:
        if "Gene" in df.columns:
            # Keep rows where gene contains coordinates (chr:pos) or ENSEMBL IDs
            intergenic_mask = (
                df["Gene"].str.contains(r"^chr\d+:", regex=True, na=False) |
                df["Gene"].str.startswith("ENSG", na=False) |
                df["Gene"].str.lower().str.contains("intergenic", na=False)
            )
            df = df[intergenic_mask]
            print(f"Filtered to {len(df)} intergenic/unannotated regions")
    
    # Limit to top N
    df = df.head(top_n)
    
    # Add slop (extend regions)
    df["Start_slop"] = (df["Start"] - slop).clip(lower=0)
    df["End_slop"] = df["End"] + slop
    
    # Create BED file
    bed_path = os.path.join(output_dir, "driver_regions.bed")
    os.makedirs(output_dir, exist_ok=True)
    
    bed_df = pd.DataFrame({
        "chrom": df["Chr"],
        "start": df["Start_slop"],
        "end": df["End_slop"],
        "name": df["Region"] if "Region" in df.columns else df.index.astype(str),
        "score": df.get("L2_Coef", 0).abs() if "L2_Coef" in df.columns else 0,
        "strand": "."
    })
    
    bed_df.to_csv(bed_path, sep="\t", header=False, index=False)
    print(f"Created BED file: {bed_path} ({len(bed_df)} regions)")
    
    return bed_path


def create_alphagenome_config(output_dir: str, bed_path: str, 
                              api_key_env: str, top_n: int) -> str:
    """Create AlphaGenome workflow config file."""
    config = {
        "paths": {
            "base": output_dir,
            "outputs": os.path.join(output_dir, "Outputs")
        },
        "input": {
            "regions_bed": bed_path
        },
        "alphagenome": {
            "enable": True,
            "api_key_env": api_key_env,
            "ontology_terms": [
                "UBERON:0007650",  # esophagus muscularis
                "UBERON:0001043",  # esophagus
                "UBERON:0000945",  # stomach
                "UBERON:0000178"   # blood
            ],
            "requested_outputs": ["DNASE", "CAGE", "CHIP_TF"],
            "sequence_length_key": "SEQUENCE_LENGTH_16KB",
            "top_n_regions": top_n,
            "ism_max_len_bp": 512,
            "ism_quantile": 0.995,
            "min_hotspot_bp": 8,
            "ism_output": "DNASE"
        }
    }
    
    config_path = os.path.join(output_dir, "config.yaml")
    with open(config_path, "w") as f:
        yaml.dump(config, f, default_flow_style=False)
    
    print(f"Created config file: {config_path}")
    return config_path


def setup_workflow(output_dir: str) -> str:
    """Setup Snakemake workflow files."""
    # Source workflow directory
    src_workflow = "/home/dirk/MEDIPIPE_PROJECT/work/post_analysis/DMR_ML_Discrimination/AlphaGenome_Motif/workflow"
    
    # Create workflow directory structure
    workflow_dir = os.path.join(output_dir, "workflow")
    os.makedirs(os.path.join(workflow_dir, "rules"), exist_ok=True)
    os.makedirs(os.path.join(workflow_dir, "scripts"), exist_ok=True)
    os.makedirs(os.path.join(workflow_dir, "envs"), exist_ok=True)
    
    # Copy/link workflow files
    import shutil
    for subdir in ["rules", "scripts", "envs"]:
        src = os.path.join(src_workflow, subdir)
        dst = os.path.join(workflow_dir, subdir)
        if os.path.exists(src):
            for f in os.listdir(src):
                src_file = os.path.join(src, f)
                dst_file = os.path.join(dst, f)
                if os.path.isfile(src_file):
                    shutil.copy2(src_file, dst_file)
    
    # Create Snakefile
    snakefile_content = '''configfile: "config.yaml"

include: "workflow/rules/alphagenome_motif.smk"

rule all:
    input:
        "Outputs/tables/region_tracks.tsv",
        "Outputs/tables/ism_hotspots.tsv",
        "Outputs/tables/ism_bedgraph.tsv",
        "Outputs/tables/final_report.tsv"
'''
    
    snakefile_path = os.path.join(output_dir, "snakefile")
    with open(snakefile_path, "w") as f:
        f.write(snakefile_content)
    
    print(f"Created workflow in: {output_dir}")
    return snakefile_path


def run_workflow(output_dir: str, dry_run: bool = False):
    """Run the Snakemake workflow."""
    cmd = ["snakemake", "--use-conda", "--cores", "4"]
    if dry_run:
        cmd.append("-n")
    
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=output_dir)
    
    if result.returncode != 0:
        print(f"ERROR: Workflow failed with return code {result.returncode}")
        sys.exit(1)
    
    print("Workflow completed successfully!")


def main():
    args = parse_args()
    
    # Check API key
    api_key = os.environ.get(args.api_key_env, "")
    if not api_key and not args.dry_run:
        print(f"WARNING: Environment variable {args.api_key_env} not set!")
        print("Set it with: export ALPHAGENOME_API_KEY='your-api-key'")
    
    # Prepare BED file
    bed_path = prepare_driver_bed(
        args.drivers, 
        args.output_dir,
        args.filter_intergenic,
        args.slop,
        args.top_n
    )
    
    # Create config
    config_path = create_alphagenome_config(
        args.output_dir,
        bed_path,
        args.api_key_env,
        args.top_n
    )
    
    # Setup workflow
    snakefile_path = setup_workflow(args.output_dir)
    
    print(f"\n{'='*60}")
    print("AlphaGenome workflow prepared!")
    print(f"{'='*60}")
    print(f"Output directory: {args.output_dir}")
    print(f"BED file: {bed_path}")
    print(f"Config: {config_path}")
    print(f"\nTo run manually:")
    print(f"  cd {args.output_dir}")
    print(f"  export {args.api_key_env}='your-api-key'")
    print(f"  snakemake --use-conda --cores 4")
    
    if not args.dry_run:
        print("\nRunning workflow...")
        run_workflow(args.output_dir, dry_run=False)
    else:
        print("\n(Dry run - workflow not executed)")


if __name__ == "__main__":
    main()
