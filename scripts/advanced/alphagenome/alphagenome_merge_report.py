import os
import pandas as pd


def read_bed_any(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    df = df.rename(columns={0: "chrom", 1: "start", 2: "end"})
    for i in range(3, df.shape[1]):
        df = df.rename(columns={i: f"c{i}"})
    df["region_id"] = range(len(df))
    return df


def main():
    tracks = pd.read_csv(snakemake.input["tracks"], sep="\t")
    hotspots = pd.read_csv(snakemake.input["hotspots"], sep="\t")
    regions = read_bed_any(snakemake.input["regions"])

    out_report = snakemake.output["report"]
    os.makedirs(os.path.dirname(out_report), exist_ok=True)

    # 1) DNASE (oder gewählte modality) summary: take max region_max across tracks
    dnase = tracks[tracks["output_type"] == "DNASE"].copy()
    dnase_sum = (dnase.groupby("region_id")["region_max"]
                      .max()
                      .reset_index()
                      .rename(columns={"region_max": "dnase_region_max"}))

    # 2) Top TF tracks: CHIP_TF sorted by region_max
    tf = tracks[tracks["output_type"] == "CHIP_TF"].copy()
    if not tf.empty:
        tf = tf.sort_values(["region_id", "region_max"], ascending=[True, False])
        top_tf = (tf.groupby("region_id")
                    .head(10)
                    .groupby("region_id")["track_name"]
                    .apply(lambda x: ",".join([str(v) for v in x.tolist()]))
                    .reset_index()
                    .rename(columns={"track_name": "top_TF_tracks"}))
    else:
        top_tf = pd.DataFrame({"region_id": [], "top_TF_tracks": []})

    # 3) Hotspot summary
    if not hotspots.empty:
        hs_sum = (hotspots.groupby("region_id")
                          .agg(hotspot_count=("hotspot_len", "size"),
                               hotspot_total_bp=("hotspot_len", "sum"),
                               hotspot_max_importance=("hotspot_max_importance", "max"))
                          .reset_index())
    else:
        hs_sum = pd.DataFrame({"region_id": [], "hotspot_count": [], "hotspot_total_bp": [], "hotspot_max_importance": []})

    rep = regions.merge(dnase_sum, on="region_id", how="left")
    rep = rep.merge(top_tf, on="region_id", how="left")
    rep = rep.merge(hs_sum, on="region_id", how="left")

    rep["dnase_region_max"] = rep["dnase_region_max"].fillna(0.0)
    rep["top_TF_tracks"] = rep["top_TF_tracks"].fillna("")
    rep["hotspot_count"] = rep["hotspot_count"].fillna(0).astype(int)
    rep["hotspot_total_bp"] = rep["hotspot_total_bp"].fillna(0).astype(int)
    rep["hotspot_max_importance"] = rep["hotspot_max_importance"].fillna(0.0)

    # Optional: sort by ISM + DNASE (pragmatisch)
    rep = rep.sort_values(["hotspot_max_importance", "dnase_region_max"], ascending=False)

    rep.to_csv(out_report, sep="\t", index=False)


def cli_main():
    """CLI entry point for standalone execution."""
    import argparse
    parser = argparse.ArgumentParser(description="AlphaGenome merge report")
    parser.add_argument("--tracks", required=True, help="Input region tracks TSV")
    parser.add_argument("--hotspots", required=True, help="Input ISM hotspots TSV")
    parser.add_argument("--regions", required=True, help="Input regions BED file")
    parser.add_argument("--output", required=True, help="Output report TSV")
    args = parser.parse_args()
    
    # Create mock snakemake object for compatibility
    class MockSnakemake:
        def __init__(self, args):
            self.input = {
                "tracks": args.tracks,
                "hotspots": args.hotspots,
                "regions": args.regions
            }
            self.output = {"report": args.output}
    
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

