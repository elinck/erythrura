configfile: "config/config.yaml"

include: "workflow/rules/all_sites_vcf.smk"
include: "workflow/rules/fst_outliers_annotate.smk"
include: "workflow/rules/admixture.smk"

K_LIST = list(range(1, 6)) 
REPS = 10

# Default target (pick what you want "snakemake" to do)
rule all:
    input:
        "results/outliers/top0.1pct_fst_windows_with_genes.tsv.gz"
        # If/when you want all-sites as default too, add:
        # "results/all_sites/all_sites.vcf.gz.tbi"

rule outliers:
    input:
        "results/outliers/top0.1pct_fst_windows_with_genes.tsv.gz"

rule allsites:
    input:
        "results/all_sites/all_sites.vcf.gz.tbi"

rule windows_annotated:
    input:
        "results/fst/all_windows_with_genes.tsv.gz"

rule admixture:
    input:
        expand("results/admixture/runs/K{K}/rep{rep}/erythrura.{K}.Q",
               K=K_LIST, rep=range(1, REPS+1)),
        expand("results/admixture/runs/K{K}/rep{rep}/erythrura.{K}.P",
               K=K_LIST, rep=range(1, REPS+1)),
        expand("results/admixture/runs/K{K}/rep{rep}/run.log",
               K=K_LIST, rep=range(1, REPS+1))

