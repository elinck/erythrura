#!/usr/bin/env python3
import argparse
import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--intersections", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--candidate", default="BMP4")
    args = ap.parse_args()

    cols = [
        "chrom","start","end","win_id","mean_fst","rank",
        "g_chrom","g_start","g_end","gene","g_score","g_strand"
    ]
    df = pd.read_csv(args.intersections, sep="\t", header=None, names=cols)

    # Summarize to one row per window
    df["gene"] = df["gene"].fillna(".").astype(str)
    summ = (
        df.groupby(["chrom","start","end","win_id","mean_fst","rank"], as_index=False)
          .agg(
              gene_hits=("gene", lambda x: len({g for g in x if g != "."})),
              genes=("gene", lambda x: ",".join(sorted({g for g in x if g != "."})) if any(g != "." for g in x) else ".")
          )
    )
    summ["candidate_hit"] = summ["genes"].apply(lambda s: args.candidate in s.split(",") if s != "." else False)
    summ = summ.sort_values(["mean_fst","rank"], ascending=[False, True])
    summ.to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    main()

