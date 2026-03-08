import os
import re
import pandas as pd
import numpy as np
from tqdm import tqdm

from alphagenome.data import genome
from alphagenome.models import dna_client


def read_bed_any(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    if df.shape[1] < 3:
        raise ValueError("BED requires >=3 columns: chrom, start, end")
    df = df.rename(columns={0: "chrom", 1: "start", 2: "end"})
    for i in range(3, df.shape[1]):
        df = df.rename(columns={i: f"c{i}"})
    return df


def get_seq_len(seq_key: str) -> int:
    if not hasattr(dna_client, seq_key):
        raise ValueError(f"Unknown sequence_length_key={seq_key}. "
                         f"Expected one of SUPPORTED_SEQUENCE_LENGTHS keys like SEQUENCE_LENGTH_16KB.")
    return int(getattr(dna_client, seq_key))


def to_output_types(names):
    outs = []
    for n in names:
        if not hasattr(dna_client.OutputType, n):
            raise ValueError(f"Unknown OutputType '{n}'. Check dna_client.OutputType list.")
        outs.append(getattr(dna_client.OutputType, n))
    return outs


def slice_track_to_interval(tdata, interval: genome.Interval) -> np.ndarray | None:
    """
    TrackData stores:
      - tdata.interval (Interval)
      - tdata.resolution (bp per value)
      - tdata.values shape (sequence_length, n_tracks)
    We slice by converting genomic coords -> indices.
    """
    if tdata is None or tdata.interval is None or tdata.resolution is None:
        return None

    iv = tdata.interval
    res = int(tdata.resolution)

    ov_start = max(interval.start, iv.start)
    ov_end = min(interval.end, iv.end)
    if ov_end <= ov_start:
        return None

    i0 = (ov_start - iv.start) // res
    i1 = int(np.ceil((ov_end - iv.start) / res))
    i0 = max(i0, 0)
    i1 = min(i1, tdata.values.shape[0])
    if i1 <= i0:
        return None
    return tdata.values[i0:i1, :]


def main():
    in_bed = snakemake.input["regions"]
    out_tsv = snakemake.output["tsv"]
    logf = snakemake.log[0]

    api_key_env = snakemake.params["api_key_env"]
    ontology_terms = list(snakemake.params["ontology_terms"])
    requested_outputs = list(snakemake.params["requested_outputs"])
    seq_len_key = snakemake.params["seq_len_key"]
    top_n = int(snakemake.params["top_n"])

    api_key = os.environ.get(api_key_env, "").strip()
    if not api_key:
        raise RuntimeError(f"Missing API key. Set env var {api_key_env} قبل running.")

    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    os.makedirs(os.path.dirname(logf), exist_ok=True)

    df = read_bed_any(in_bed)

    # Optional sorting by 4th BED column if numeric
    sort_col = "c3" if "c3" in df.columns else None
    if sort_col is not None:
        vals = pd.to_numeric(df[sort_col], errors="coerce")
        if vals.notna().any():
            df["_score"] = vals
            df = df.sort_values("_score", ascending=False).drop(columns=["_score"])

    df = df.head(top_n).reset_index(drop=True)
    df["region_id"] = df.index.astype(int)

    seq_len = get_seq_len(seq_len_key)

    model = dna_client.create(api_key)  # official client creation :contentReference[oaicite:2]{index=2}
    out_types = to_output_types(requested_outputs)

    rows = []
    with open(logf, "w") as log:
        log.write(f"regions={len(df)} seq_len={seq_len} ontology_terms={ontology_terms} outputs={requested_outputs}\n")

        for _, r in tqdm(df.iterrows(), total=len(df)):
            rid = int(r["region_id"])
            chrom = str(r["chrom"])
            start = int(r["start"])
            end = int(r["end"])

            mid = (start + end) // 2
            half = seq_len // 2
            ctx_start = max(mid - half, 0)
            ctx_end = ctx_start + seq_len

            ctx_iv = genome.Interval(chromosome=chrom, start=ctx_start, end=ctx_end)
            reg_iv = genome.Interval(chromosome=chrom, start=start, end=end)

            try:
                output = model.predict_interval(
                    interval=ctx_iv,
                    requested_outputs=out_types,
                    ontology_terms=ontology_terms,
                )  # :contentReference[oaicite:3]{index=3}
            except Exception as e:
                log.write(f"[WARN] predict_interval failed region_id={rid} {chrom}:{start}-{end}: {e}\n")
                continue

            for ot in out_types:
                attr = ot.name.lower()
                if not hasattr(output, attr):
                    continue

                tdata = getattr(output, attr)
                arr = slice_track_to_interval(tdata, reg_iv)
                if arr is None:
                    continue

                means = np.nanmean(arr, axis=0)
                maxs = np.nanmax(arr, axis=0)
                meta = tdata.metadata.copy() if tdata.metadata is not None else pd.DataFrame()

                # Ensure meta rows = tracks
                for j in range(arr.shape[1]):
                    m = meta.iloc[j].to_dict() if (not meta.empty and j < len(meta)) else {}
                    row = {
                        "region_id": rid,
                        "chrom": chrom, "start": start, "end": end,
                        "context_start": ctx_start, "context_end": ctx_end,
                        "output_type": ot.name,
                        "track_idx": j,
                        "track_name": m.get("name", None),
                        "strand": m.get("strand", None),
                        "ontology_curie": m.get("ontology_curie", None),
                        "biosample_name": m.get("biosample_name", None),
                        "data_source": m.get("data_source", None),
                        "region_mean": float(means[j]),
                        "region_max": float(maxs[j]),
                    }
                    # carry over optional BED extra cols
                    for c in [c for c in df.columns if c.startswith("c")]:
                        row[c] = r.get(c, None)
                    rows.append(row)

    pd.DataFrame(rows).to_csv(out_tsv, sep="\t", index=False)


def cli_main():
    """CLI entry point for standalone execution."""
    import argparse
    parser = argparse.ArgumentParser(description="AlphaGenome region track prediction")
    parser.add_argument("--regions", required=True, help="Input BED file with regions")
    parser.add_argument("--output", required=True, help="Output TSV file")
    parser.add_argument("--api-key-env", default="ALPHAGENOME_API_KEY", help="Env var for API key")
    parser.add_argument("--ontology-terms", default="UBERON:0007650,UBERON:0001043", help="Comma-separated ontology terms")
    parser.add_argument("--requested-outputs", default="DNASE,CAGE,CHIP_TF", help="Comma-separated output types")
    parser.add_argument("--sequence-length-key", default="SEQUENCE_LENGTH_16KB", help="Sequence length key")
    parser.add_argument("--top-n", type=int, default=100, help="Top N regions to analyze")
    args = parser.parse_args()
    
    # Create mock snakemake object for compatibility
    class MockSnakemake:
        def __init__(self, args):
            self.input = {"regions": args.regions}
            self.output = {"tsv": args.output}
            self.log = [args.output.replace(".tsv", ".log")]
            self.params = {
                "api_key_env": args.api_key_env,
                "ontology_terms": args.ontology_terms.replace("[", "").replace("]", "").replace("'", "").split(","),
                "requested_outputs": args.requested_outputs.replace("[", "").replace("]", "").replace("'", "").split(","),
                "seq_len_key": args.sequence_length_key,
                "top_n": args.top_n
            }
    
    global snakemake
    snakemake = MockSnakemake(args)
    main()


if __name__ == "__main__":
    # Check if running via Snakemake or CLI
    try:
        snakemake
        main()
    except NameError:
        cli_main()

