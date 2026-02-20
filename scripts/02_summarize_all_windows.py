#!/usr/bin/env python3
import argparse
import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fst", required=True)          # original vcftools tsv
    ap.add_argument("--bed", required=True)          # all_windows.bed
    ap.add_argument("--intersections", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    # Read original vcftools output (keeps all numeric columns)
    fst = pd.read_csv(args.fst, sep=r"\s+")
    # Make helper coords matching BED conversion:
    fst["chrom"] = fst["CHROM"].astype(str)
    fst["start0"] = fst["BIN_START"].astype(int) - 1
    fst["end"] = fst["BIN_END"].astype(int)

    # Read intersections:
    # win bed: chrom start end win_id mean_fst
    # genes bed: g_chrom g_start g_end gene . strand
    cols = ["chrom","start0","end","win_id","mean_fst",
            "g_chrom","g_start","g_end","gene","g_score","g_strand"]
    inter = pd.read_csv(args.intersections, sep="\t", header=None, names=cols)

    # Summarize genes per window
    inter["gene"] = inter["gene"].fillna(".").astype(str)
    summ = (
        inter.groupby(["chrom","start0","end"], as_index=False)
             .agg(
                 gene_hits=("gene", lambda x: len({g for g in x if g != "."})),
                 genes=("gene", lambda x: ",".join(sorted({g for g in x if g != "."}))
                                 if any(g != "." for g in x) else ".")
             )
    )

    merged = fst.merge(summ, on=["chrom","start0","end"], how="left")
    merged["gene_hits"] = merged["gene_hits"].fillna(0).astype(int)
    merged["genes"] = merged["genes"].fillna(".")

    # Drop helpers
    merged = merged.drop(columns=["chrom","start0"])

    merged.to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    main()

