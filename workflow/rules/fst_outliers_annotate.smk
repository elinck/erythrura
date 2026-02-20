ACC = config["reference"]["accession"]
FST = config["fst_scan"]["vcftools_windows"]
Q   = config["fst_scan"]["quantile"]
CAND = config["fst_scan"]["candidate_gene"]

DATASET_ZIP = f"resources/ncbi/{ACC}.zip"
DATASET_DIR = f"resources/ncbi/{ACC}"

rule fst_outliers_with_genes:
    input:
        "results/outliers/top0.1pct_fst_windows_with_genes.tsv.gz"

rule download_ncbi_dataset:
    conda:
        "envs/fst_outliers.yaml"
    output:
        DATASET_ZIP
    params:
        acc=ACC
    shell:
        r"""
        set -euo pipefail
        mkdir -p resources/ncbi
        datasets download genome accession {params.acc} \
          --include genome,gff3,gtf,protein,cds,rna \
          --filename {output}
        """

rule unpack_ncbi_dataset:
    input:
        DATASET_ZIP
    output:
        touch(f"{DATASET_DIR}/.ready")
    shell:
        r"""
        set -euo pipefail
        mkdir -p {DATASET_DIR}
        unzip -o {input} -d {DATASET_DIR}
        touch {output}
        """

rule genes_bed_from_gff:
    conda:
        "envs/fst_outliers.yaml"
    input:
        ready=f"{DATASET_DIR}/.ready"
    output:
        bed="resources/annotation/genes.bed"
    shell:
        r"""
        set -euo pipefail
        mkdir -p resources/annotation

        GFF=$(find {DATASET_DIR} -type f \( -name "*.gff" -o -name "*.gff3" \) | head -n 1)
        if [ -z "$GFF" ]; then
          echo "ERROR: No GFF/GFF3 found under {DATASET_DIR}" >&2
          exit 1
        fi

        awk -F'\t' 'BEGIN{{OFS="\t"}}
          $0 !~ /^#/ && $3=="gene" {{
            chr=$1; start=$4-1; end=$5; strand=$7; attr=$9;

            name=".";
            if (match(attr, /Name=([^;]+)/, m)) name=m[1];
            else if (match(attr, /gene=([^;]+)/, m)) name=m[1];
            else if (match(attr, /gene_name=([^;]+)/, m)) name=m[1];
            else if (match(attr, /ID=gene-?([^;]+)/, m)) name=m[1];
            else if (match(attr, /ID=([^;]+)/, m)) name=m[1];

            print chr, start, end, name, ".", strand;
          }}' "$GFF" \
        | sort -k1,1 -k2,2n > {output.bed}
        """

rule vcftools_windows_to_bed:
    conda:
        "envs/fst_outliers.yaml"
    input:
        FST
    output:
        bed="results/fst/all_windows.bed"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/fst
        # BED: chrom start0 end name mean_fst
        awk 'BEGIN{{OFS="\t"}} NR==1{{next}}
             $6!="NA" && $6!="-nan" {{
               print $1, $2-1, $3, "win_"$1":"$2"-"$3, $6
             }}' {input} \
        | sort -k1,1 -k2,2n > {output.bed}
        """

rule top_windows_from_vcftools:
    conda:
        "envs/fst_outliers.yaml"
    input:
        FST
    output:
        bed="results/outliers/top0.1pct_fst_windows.bed",
        tsv="results/outliers/top0.1pct_fst_windows.tsv"
    params:
        q=Q
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/outliers
        python scripts/00_vcf_top_windows.py \
          --in {input} \
          --quantile {params.q} \
          --bed {output.bed} \
          --tsv {output.tsv}
        """

rule intersect_all_windows_genes:
    conda:
        "envs/fst_outliers.yaml"
    input:
        win="results/fst/all_windows.bed",
        genes="resources/annotation/genes.bed"
    output:
        "results/fst/all_window_gene_intersections.tsv"
    shell:
        r"""
        set -euo pipefail
        bedtools intersect -wa -wb -a {input.win} -b {input.genes} > {output}
        """

rule intersect_outliers_genes:
    conda:
        "envs/fst_outliers.yaml"
    input:
        outliers="results/outliers/top0.1pct_fst_windows.bed",
        genes="resources/annotation/genes.bed"
    output:
        "results/outliers/outlier_gene_intersections.tsv"
    shell:
        r"""
        set -euo pipefail
        bedtools intersect -wa -wb -a {input.outliers} -b {input.genes} > {output}
        """
rule all_windows_with_genes:
    conda:
        "envs/fst_outliers.yaml"
    input:
        fst=FST,
        bed="results/fst/all_windows.bed",
        inter="results/fst/all_window_gene_intersections.tsv"
    output:
        "results/fst/all_windows_with_genes.tsv.gz"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/fst
        python scripts/02_summarize_all_windows.py \
          --fst {input.fst} \
          --bed {input.bed} \
          --intersections {input.inter} \
          --out results/fst/all_windows_with_genes.tsv

        gzip -f results/fst/all_windows_with_genes.tsv
        """

rule summarize_outliers_with_genes:
    conda:
        "envs/fst_outliers.yaml"
    input:
        out="results/outliers/top0.1pct_fst_windows.tsv",
        inter="results/outliers/outlier_gene_intersections.tsv"
    output:
        "results/outliers/top0.1pct_fst_windows_with_genes.tsv.gz"
    params:
        cand=CAND,
        acc=ACC
    shell:
        r"""
        set -euo pipefail

        python scripts/01_summarize_intersections.py \
          --intersections {input.inter} \
          --candidate {params.cand} \
          --out results/outliers/_tmp_with_genes.tsv

        python - <<'PY'
import pandas as pd
acc = "{params.acc}"

outliers = pd.read_csv("{input.out}", sep="\t")

# create BED coords from vcftools coords
outliers["chrom"] = outliers["CHROM"].astype(str)
outliers["start"] = outliers["BIN_START"].astype(int) - 1
outliers["end"]   = outliers["BIN_END"].astype(int)

summ = pd.read_csv("results/outliers/_tmp_with_genes.tsv", sep="\t")

merged = outliers.merge(
    summ[["chrom","start","end","gene_hits","genes","candidate_hit"]],
    on=["chrom","start","end"],
    how="left"
)

merged["gene_hits"] = merged["gene_hits"].fillna(0).astype(int)
merged["genes"] = merged["genes"].fillna(".")
merged["candidate_hit"] = merged["candidate_hit"].fillna(False)
merged["ref_accession"] = acc
merged = merged.drop(columns=["chrom","start","end"])

# nice ordering
cols_first = ["ref_accession","rank","CHROM","BIN_START","BIN_END","N_VARIANTS","WEIGHTED_FST","MEAN_FST",
              "gene_hits","genes","candidate_hit","cutoff_mean_fst","quantile"]
cols = [c for c in cols_first if c in merged.columns] + [c for c in merged.columns if c not in cols_first]
merged = merged[cols].sort_values(["MEAN_FST","rank"], ascending=[False, True])

merged.to_csv("results/outliers/top0.1pct_fst_windows_with_genes.tsv", sep="\t", index=False)
PY

        bgzip -f results/outliers/top0.1pct_fst_windows_with_genes.tsv
        rm -f results/outliers/_tmp_with_genes.tsv
        """

