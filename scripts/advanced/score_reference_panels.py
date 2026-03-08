#!/usr/bin/env python3
"""
Score samples against Nature Cancer cfMeDIP-seq reference panels.
Computes within-cohort z-scores vs controls for each DMR panel.

Based on: AEG_vs_CTRL_ReferenceDB_NatureCancer_SignatureMatch_HOWTO.txt
"""

import argparse
import os
import glob
import numpy as np
import pandas as pd
from pathlib import Path


def read_bw_map(path):
    """Read BigWig mapping file (sample_id, group, bw_path)"""
    m = pd.read_csv(path, sep="\t", header=None, names=["sample_id", "group", "bw"])
    # Label is the sample_id (since we passed --labels to multiBigwigSummary)
    m["label"] = m["sample_id"].astype(str)
    return m


def load_panel_raw(path):
    """Load deepTools multiBigwigSummary output"""
    df = pd.read_csv(path, sep="\t")
    # Clean column names (remove quotes added by deepTools)
    df.columns = [c.strip("'").strip('"').replace("#", "").strip() for c in df.columns]
    coord_cols = {"chr", "start", "end", "gene", "name", "score", "strand"}
    sample_cols = [c for c in df.columns if c not in coord_cols]
    return df, sample_cols


def panel_score(matrix, log1p=True):
    """Calculate panel score as mean(log1p(signal)) across regions"""
    x = matrix.copy()
    if log1p:
        x = np.log1p(x)
    return x.mean(axis=0)


def main():
    ap = argparse.ArgumentParser(description="Score samples against reference DMR panels")
    ap.add_argument("--bw-map", required=True, help="BigWig mapping file (sample_id, group, bw_path)")
    ap.add_argument("--panel-raw-glob", required=True, help="Glob pattern for panel raw.tsv files")
    ap.add_argument("--control-label", default="CTRL", help="Label for control group [default: CTRL]")
    ap.add_argument("--out-prefix", required=True, help="Output file prefix")
    args = ap.parse_args()

    # Load BigWig mapping
    bw = read_bw_map(args.bw_map)
    label2sid = dict(zip(bw["label"], bw["sample_id"]))
    sid2group = dict(zip(bw["sample_id"], bw["group"]))

    print(f"Loaded {len(bw)} samples from BigWig map")
    print(f"Groups: {bw['group'].value_counts().to_dict()}")
    print(f"Control label: {args.control_label}")

    rows = []
    raw_files = sorted(glob.glob(args.panel_raw_glob))
    if not raw_files:
        raise SystemExit(f"No files matched: {args.panel_raw_glob}")

    print(f"\nProcessing {len(raw_files)} panels...")
    
    for raw_path in raw_files:
        panel_id = os.path.basename(raw_path).replace(".raw.tsv", "")
        print(f"  - {panel_id}")
        
        df, sample_cols = load_panel_raw(raw_path)

        keep, rename = [], {}
        for c in sample_cols:
            if c in label2sid:
                sid = label2sid[c]
                keep.append(c)
                rename[c] = sid
        
        if not keep:
            print(f"    WARNING: No sample columns matched bw_map labels")
            continue

        mat = df[keep].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        mat = mat.rename(columns=rename)

        score = panel_score(mat.values, log1p=True)
        s = pd.Series(score, index=mat.columns, name="score")

        tmp = pd.DataFrame({
            "panel_id": panel_id,
            "sample_id": s.index,
            "score": s.values,
            "group": [sid2group.get(i, "NA") for i in s.index]
        })

        # Calculate z-scores vs controls
        ctrl = tmp[tmp["group"] == args.control_label]["score"].astype(float)
        mu = float(ctrl.mean()) if len(ctrl) else float("nan")
        sd = float(ctrl.std(ddof=1)) if len(ctrl) > 1 else float("nan")
        tmp["z_vs_ctrl"] = (tmp["score"] - mu) / sd if (sd and sd > 0) else float("nan")
        tmp["ctrl_mean"] = mu
        tmp["ctrl_sd"] = sd
        
        rows.append(tmp)
        print(f"    {len(tmp)} samples, ctrl_mean={mu:.4f}, ctrl_sd={sd:.4f}")

    # Combine all panels
    long = pd.concat(rows, ignore_index=True)
    long_output = args.out_prefix + ".scores.long.tsv"
    long.to_csv(long_output, sep="\t", index=False)
    print(f"\nLong-format scores saved to: {long_output}")

    # Create summary table
    summ = []
    for pid, g in long.groupby("panel_id", sort=True):
        a = g[g["group"] == args.control_label]["z_vs_ctrl"].dropna()
        b = g[g["group"] != args.control_label]["z_vs_ctrl"].dropna()
        summ.append({
            "panel_id": pid,
            "n_ctrl": int((g["group"] == args.control_label).sum()),
            "n_case": int((g["group"] != args.control_label).sum()),
            "mean_z_ctrl": float(a.mean()) if len(a) else float("nan"),
            "mean_z_case": float(b.mean()) if len(b) else float("nan"),
            "median_z_ctrl": float(a.median()) if len(a) else float("nan"),
            "median_z_case": float(b.median()) if len(b) else float("nan"),
        })
    
    summ_df = pd.DataFrame(summ)
    summ_output = args.out_prefix + ".scores.summary.tsv"
    summ_df.to_csv(summ_output, sep="\t", index=False)
    print(f"Summary saved to: {summ_output}")

    # Print interpretation
    print("\n=== Interpretation Guide ===")
    print("pancancer_hyper: Case samples should show POSITIVE z-scores (hyper-methylated)")
    print("pancancer_hypo: Case samples should show NEGATIVE z-scores (hypo-methylated)")
    print("Cancer-type panels: Treat as similarity signals, not definitive type calls")


if __name__ == "__main__":
    main()
