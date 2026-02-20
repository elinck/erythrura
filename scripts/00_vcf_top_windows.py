#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--quantile", type=float, default=0.999)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--tsv", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.inp, sep=r"\s+")
    need = ["CHROM","BIN_START","BIN_END","N_VARIANTS","WEIGHTED_FST","MEAN_FST"]
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise SystemExit(f"Missing columns {missing}. Found: {list(df.columns)}")

    # Coerce and clean MEAN_FST
    df = df.replace([np.inf, -np.inf], np.nan)
    df["MEAN_FST"] = pd.to_numeric(df["MEAN_FST"], errors="coerce")
    df = df.dropna(subset=["MEAN_FST"]).copy()

    # Keep plausible MEAN_FST (vcftools can yield negative or >1 in some edge cases;
    # if you want to keep those, remove this filter)
    df = df[(df["MEAN_FST"] >= 0)].copy()

    q = float(args.quantile)
    if not (0 < q < 1):
        raise SystemExit("--quantile must be between 0 and 1")

    cutoff = df["MEAN_FST"].quantile(q)
    out = df[df["MEAN_FST"] >= cutoff].copy()

    # Rank
    out = out.sort_values("MEAN_FST", ascending=False).reset_index(drop=True)
    out["rank"] = np.arange(1, len(out) + 1)
    out["quantile"] = q
    out["cutoff_mean_fst"] = cutoff

    # Make BED for intersections later
    # vcftools BIN_START/BIN_END are 1-based coordinates;
    # BED expects 0-based start, half-open end
    bed = pd.DataFrame({
        "chrom": out["CHROM"].astype(str),
        "start": out["BIN_START"].astype(int) - 1,
        "end": out["BIN_END"].astype(int),
        "win_id": out.apply(lambda r: f"fst_win_{r['CHROM']}:{int(r['BIN_START'])}-{int(r['BIN_END'])}", axis=1),
        "mean_fst": out["MEAN_FST"].astype(float),
        "rank": out["rank"].astype(int),
    })

    bed.to_csv(args.bed, sep="\t", header=False, index=False)
    out.to_csv(args.tsv, sep="\t", index=False)

if __name__ == "__main__":
    main()

