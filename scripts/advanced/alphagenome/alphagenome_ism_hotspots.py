import os
import pandas as pd
import numpy as np
from tqdm import tqdm

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers
from alphagenome.interpretation import ism


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
        raise ValueError(f"Unknown sequence_length_key={seq_key}.")
    return int(getattr(dna_client, seq_key))


def clamp_ism_interval(chrom: str, start: int, end: int, max_len: int) -> tuple[int, int]:
    """Return an ISM sub-interval within [start,end) with length <= max_len, centered."""
    length = end - start
    if length <= max_len:
        return start, end
    mid = (start + end) // 2
    half = max_len // 2
    s = mid - half
    e = s + max_len
    if s < start:
        s = start
        e = s + max_len
    if e > end:
        e = end
        s = e - max_len
    return s, e


def main():
    in_bed = snakemake.input["regions"]
    out_hotspots = snakemake.output["hotspots"]
    out_bedgraph = snakemake.output["bedgraph"]
    logf = snakemake.log[0]

    api_key_env = snakemake.params["api_key_env"]
    ontology_terms = list(snakemake.params["ontology_terms"])
    seq_len_key = snakemake.params["seq_len_key"]
    ism_output = snakemake.params["ism_output"]
    top_n = int(snakemake.params["top_n"])
    ism_max_len = int(snakemake.params["ism_max_len"])
    q = float(snakemake.params["q"])
    min_hotspot = int(snakemake.params["min_hotspot"])

    api_key = os.environ.get(api_key_env, "").strip()
    if not api_key:
        raise RuntimeError(f"Missing API key. Set env var {api_key_env} before running.")

    os.makedirs(os.path.dirname(out_hotspots), exist_ok=True)
    os.makedirs(os.path.dirname(logf), exist_ok=True)

    df = read_bed_any(in_bed).head(top_n).reset_index(drop=True)
    df["region_id"] = df.index.astype(int)

    seq_len = get_seq_len(seq_len_key)

    model = dna_client.create(api_key)

    # Score on one modality (DNASE/ATAC etc.)
    if not hasattr(dna_client.OutputType, ism_output):
        raise ValueError(f"Unknown ism_output={ism_output}")
    requested_output = getattr(dna_client.OutputType, ism_output)

    # CenterMaskScorer is a standard, lightweight scorer for regulatory activity (mask-based aggregation). :contentReference[oaicite:4]{index=4}
    scorer = variant_scorers.CenterMaskScorer(requested_output=requested_output)

    bedgraph_rows = []
    hotspot_rows = []

    with open(logf, "w") as log:
        log.write(f"regions={len(df)} seq_len={seq_len} ontology_terms={ontology_terms} ism_output={ism_output}\n")

        for _, r in tqdm(df.iterrows(), total=len(df)):
            rid = int(r["region_id"])
            chrom = str(r["chrom"])
            start = int(r["start"])
            end = int(r["end"])

            mid = (start + end) // 2
            half = seq_len // 2
            ctx_start = max(mid - half, 0)
            ctx_end = ctx_start + seq_len
            sequence_interval = genome.Interval(chromosome=chrom, start=ctx_start, end=ctx_end)

            ism_s, ism_e = clamp_ism_interval(chrom, start, end, ism_max_len)
            if ism_e <= ism_s:
                continue
            ism_interval = genome.Interval(chromosome=chrom, start=ism_s, end=ism_e)

            # ISM variants list (all SNVs) – sequence is fetched internally by server for interval-based scoring
            # and can also be provided in other modes; for score_ism_variants the key objects are the intervals. :contentReference[oaicite:5]{index=5}
            try:
                scores_nested = model.score_ism_variants(
                    interval=sequence_interval,
                    ism_interval=ism_interval,
                    variant_scorers=[scorer],
                    ontology_terms=ontology_terms,
                    progress_bar=False,
                )  # :contentReference[oaicite:6]{index=6}
            except TypeError:
                # Some versions may not accept ontology_terms here; fallback without it.
                scores_nested = model.score_ism_variants(
                    interval=sequence_interval,
                    ism_interval=ism_interval,
                    variant_scorers=[scorer],
                    progress_bar=False,
                )
            except Exception as e:
                log.write(f"[WARN] score_ism_variants failed region_id={rid} {chrom}:{start}-{end}: {e}\n")
                continue

            # scores_nested: list[list[AnnData]]: one per variant, one per scorer. :contentReference[oaicite:7]{index=7}
            scalar_scores = []
            for per_variant in scores_nested:
                ad = per_variant[0]  # single scorer
                tidy = variant_scorers.tidy_scores([ad], match_gene_strand=False)  # :contentReference[oaicite:8]{index=8}
                if "raw_score" not in tidy.columns or tidy.empty:
                    scalar_scores.append(0.0)
                else:
                    scalar_scores.append(float(np.nanmax(np.abs(tidy["raw_score"].values))))

            # Build the list of all possible SNVs for the ISM interval, then map scalar scores -> (pos, base) matrix. :contentReference[oaicite:9]{index=9}
            variants = ism.ism_variants(ism_interval, skip_n=True)
            M = ism.ism_matrix(
                variant_scores=scalar_scores,
                variants=variants,
                interval=ism_interval,
                multiply_by_sequence=False,
                require_fully_filled=False,
            )

            # Importance per position: max |score| over bases (A/C/G/T)
            imp = np.nanmax(np.abs(M), axis=1)

            # Bedgraph (1bp)
            for i, val in enumerate(imp):
                gstart = ism_s + i
                bedgraph_rows.append({
                    "region_id": rid,
                    "chrom": chrom,
                    "start": gstart,
                    "end": gstart + 1,
                    "importance": float(val),
                })

            # Hotspots: contiguous segments above per-region quantile
            thr = float(np.quantile(imp, q)) if len(imp) else 0.0
            above = imp >= thr

            i = 0
            while i < len(above):
                if not above[i]:
                    i += 1
                    continue
                j = i
                while j < len(above) and above[j]:
                    j += 1
                if (j - i) >= min_hotspot:
                    hs_start = ism_s + i
                    hs_end = ism_s + j
                    hotspot_rows.append({
                        "region_id": rid,
                        "chrom": chrom,
                        "start": hs_start,
                        "end": hs_end,
                        "hotspot_len": int(hs_end - hs_start),
                        "hotspot_max_importance": float(np.nanmax(imp[i:j])),
                        "threshold": thr,
                        "parent_region_start": start,
                        "parent_region_end": end,
                        "ism_interval_start": ism_s,
                        "ism_interval_end": ism_e,
                    })
                i = j

    pd.DataFrame(bedgraph_rows).to_csv(out_bedgraph, sep="\t", index=False)
    pd.DataFrame(hotspot_rows).to_csv(out_hotspots, sep="\t", index=False)


def cli_main():
    """CLI entry point for standalone execution."""
    import argparse
    parser = argparse.ArgumentParser(description="AlphaGenome ISM hotspot detection")
    parser.add_argument("--regions", required=True, help="Input BED file with regions")
    parser.add_argument("--output-hotspots", required=True, help="Output hotspots TSV file")
    parser.add_argument("--output-bedgraph", required=True, help="Output bedgraph TSV file")
    parser.add_argument("--api-key-env", default="ALPHAGENOME_API_KEY", help="Env var for API key")
    parser.add_argument("--ontology-terms", default="UBERON:0007650,UBERON:0001043", help="Comma-separated ontology terms")
    parser.add_argument("--sequence-length-key", default="SEQUENCE_LENGTH_16KB", help="Sequence length key")
    parser.add_argument("--ism-output", default="DNASE", help="ISM output modality")
    parser.add_argument("--top-n", type=int, default=100, help="Top N regions to analyze")
    parser.add_argument("--ism-max-len", type=int, default=512, help="Max ISM interval length")
    parser.add_argument("--quantile", type=float, default=0.995, help="Quantile threshold for hotspots")
    parser.add_argument("--min-hotspot", type=int, default=8, help="Minimum hotspot length (bp)")
    args = parser.parse_args()
    
    # Create mock snakemake object for compatibility
    class MockSnakemake:
        def __init__(self, args):
            self.input = {"regions": args.regions}
            self.output = {
                "hotspots": args.output_hotspots,
                "bedgraph": args.output_bedgraph
            }
            self.log = [args.output_hotspots.replace(".tsv", ".log")]
            self.params = {
                "api_key_env": args.api_key_env,
                "ontology_terms": args.ontology_terms.replace("[", "").replace("]", "").replace("'", "").split(","),
                "seq_len_key": args.sequence_length_key,
                "ism_output": args.ism_output,
                "top_n": args.top_n,
                "ism_max_len": args.ism_max_len,
                "q": args.quantile,
                "min_hotspot": args.min_hotspot
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

