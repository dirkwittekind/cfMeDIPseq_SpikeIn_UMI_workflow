#!/usr/bin/env python3
"""
Generate QC Summary for cfMeDIP-seq Reports

Extracts and formats:
1. Unique fragments (UMI deduplication stats)
2. IP efficiency (spike-in metrics: relH, GoGe)

Usage:
    python3 generate_qc_summary.py --samples 148,149,150 --output report_qc.tsv
    python3 generate_qc_summary.py --annotation annotation.tsv --output report_qc.tsv
"""

import os
import argparse
import pandas as pd
import glob

# Default paths
UMI_STATS_DIR = "/data/medipipe_data/output/stats"
SPIKEIN_STATS_DIR = "/data/medipipe_data/output_spikein/meth_qc_quant_spikein"

# Literature benchmarks
BENCHMARKS = {
    "unique_fragments": {
        "minimum": 5_000_000,
        "optimal": 10_000_000,
        "typical_range": "5-30 million"
    },
    "relH": {
        "threshold": 2.5,
        "description": "Relative CpG enrichment"
    },
    "GoGe": {
        "threshold": 1.5,
        "description": "Observed/expected CpG ratio"
    }
}


def get_unique_fragments(sample_id):
    """Extract unique fragment count from UMI dedup stats."""
    # Try different naming patterns
    patterns = [
        f"{UMI_STATS_DIR}/{sample_id}_dedup_umi.stats.txt",
        f"{UMI_STATS_DIR}/{sample_id}_S*_L*_dedup_umi.stats.txt",
    ]
    
    for pattern in patterns:
        files = glob.glob(pattern)
        if files:
            stats_file = files[0]
            break
    else:
        return None, None
    
    try:
        with open(stats_file, 'r') as f:
            for line in f:
                if line.startswith("SN\traw total sequences:"):
                    total_reads = int(line.strip().split('\t')[2])
                elif line.startswith("SN\t1st fragments:"):
                    unique_fragments = int(line.strip().split('\t')[2])
                    return unique_fragments, total_reads
    except:
        pass
    
    return None, None


def get_spikein_metrics(sample_id):
    """Extract spike-in IP efficiency metrics."""
    # Try different naming patterns
    patterns = [
        f"{SPIKEIN_STATS_DIR}/{sample_id}_meth_qc.txt",
        f"{SPIKEIN_STATS_DIR}/{sample_id}_S*_L*_meth_qc.txt",
    ]
    
    for pattern in patterns:
        files = glob.glob(pattern)
        if files:
            qc_file = files[0]
            break
    else:
        return None, None
    
    try:
        df = pd.read_csv(qc_file, sep='\t')
        if len(df) > 0:
            relH = df['enrichment.relH'].values[0]
            GoGe = df['enrichment.GoGe'].values[0]
            return relH, GoGe
    except:
        pass
    
    return None, None


def assess_quality(unique_frags, relH, GoGe):
    """Assess sample quality based on benchmarks."""
    issues = []
    
    if unique_frags is not None:
        if unique_frags < BENCHMARKS["unique_fragments"]["minimum"]:
            issues.append(f"Low fragments ({unique_frags/1e6:.1f}M < 5M)")
        elif unique_frags < BENCHMARKS["unique_fragments"]["optimal"]:
            issues.append(f"Moderate fragments ({unique_frags/1e6:.1f}M)")
    
    if relH is not None and relH < BENCHMARKS["relH"]["threshold"]:
        issues.append(f"Low relH ({relH:.2f} < 2.5)")
    
    if GoGe is not None and GoGe < BENCHMARKS["GoGe"]["threshold"]:
        issues.append(f"Low GoGe ({GoGe:.2f} < 1.5)")
    
    if not issues:
        return "PASS"
    else:
        return "; ".join(issues)


def generate_qc_summary(sample_ids, output_file=None):
    """Generate QC summary for given samples."""
    
    results = []
    
    for sample_id in sample_ids:
        unique_frags, total_reads = get_unique_fragments(sample_id)
        relH, GoGe = get_spikein_metrics(sample_id)
        
        assessment = assess_quality(unique_frags, relH, GoGe)
        
        results.append({
            "Sample": sample_id,
            "Unique_Fragments": unique_frags,
            "Unique_Fragments_M": f"{unique_frags/1e6:.1f}M" if unique_frags else "N/A",
            "Total_Reads": total_reads,
            "Dedup_Rate": f"{100*unique_frags/total_reads:.1f}%" if unique_frags and total_reads else "N/A",
            "relH": round(relH, 2) if relH else None,
            "GoGe": round(GoGe, 2) if GoGe else None,
            "Assessment": assessment
        })
    
    df = pd.DataFrame(results)
    
    # Print summary
    print("=" * 80)
    print("QC SUMMARY")
    print("=" * 80)
    print(f"\nBenchmarks:")
    print(f"  Unique fragments: minimum {BENCHMARKS['unique_fragments']['minimum']/1e6:.0f}M, optimal >{BENCHMARKS['unique_fragments']['optimal']/1e6:.0f}M")
    print(f"  relH (CpG enrichment): >{BENCHMARKS['relH']['threshold']}")
    print(f"  GoGe (observed/expected CpG): >{BENCHMARKS['GoGe']['threshold']}")
    print()
    
    # Print table
    print(df[["Sample", "Unique_Fragments_M", "relH", "GoGe", "Assessment"]].to_string(index=False))
    
    # Summary stats
    frags = [r["Unique_Fragments"] for r in results if r["Unique_Fragments"]]
    if frags:
        print(f"\nFragment count: {min(frags)/1e6:.1f}M - {max(frags)/1e6:.1f}M (mean: {sum(frags)/len(frags)/1e6:.1f}M)")
    
    pass_count = sum(1 for r in results if r["Assessment"] == "PASS")
    print(f"Quality: {pass_count}/{len(results)} samples PASS all thresholds")
    
    # Save if output specified
    if output_file:
        df.to_csv(output_file, sep='\t', index=False)
        print(f"\nSaved: {output_file}")
    
    return df


def main():
    parser = argparse.ArgumentParser(description="Generate QC summary for cfMeDIP-seq reports")
    parser.add_argument("--samples", help="Comma-separated sample IDs")
    parser.add_argument("--annotation", help="Annotation TSV with sample_id column")
    parser.add_argument("--output", help="Output TSV file")
    parser.add_argument("--aeg-ctrl", action="store_true", help="Use standard AEG vs CTRL samples")
    args = parser.parse_args()
    
    # Determine samples
    if args.aeg_ctrl:
        sample_ids = ["148", "149", "150", "152", "154", "155", "C1", "C3", "C4", "C5", "C6"]
    elif args.samples:
        sample_ids = [s.strip() for s in args.samples.split(",")]
    elif args.annotation:
        ann = pd.read_csv(args.annotation, sep='\t')
        sample_ids = ann["sample_id"].tolist()
        # Try to extract short IDs
        sample_ids = [s.split("_")[0] for s in sample_ids]
    else:
        # Default to AEG vs CTRL
        sample_ids = ["148", "149", "150", "152", "154", "155", "C1", "C3", "C4", "C5", "C6"]
    
    generate_qc_summary(sample_ids, args.output)


if __name__ == "__main__":
    main()
