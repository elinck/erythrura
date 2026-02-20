VCF = "/home/k14m234/erythrura_assembly/results/GCF_005870125.1/QC/erythrura.pruned.vcf.gz"

# Ks + replicates
K_LIST = list(range(1, 6))   # change as desired
REPS = 10                    # replicates per K
SEED_BASE = 1000

# Filtering knobs (tune if too strict)
GENO = 0.10   # drop variants with >10% missing
MAF  = 0.05   # drop rare variants

OUTLIERS = "config/admixture_outliers.txt"

rule vcf_to_plink:
    input:
        vcf=VCF
    output:
        bed="results/admixture/erythrura.bed",
        bim="results/admixture/erythrura.bim",
        fam="results/admixture/erythrura.fam"
    conda:
        "envs/admixture.yaml"
    shell:
        r"""
        mkdir -p results/admixture
        plink2 --vcf {input.vcf} \
               --double-id \
               --allow-extra-chr \
               --make-bed \
               --out results/admixture/erythrura
        """

rule make_numeric_chr:
    input:
        bed="results/admixture/erythrura.bed",
        bim="results/admixture/erythrura.bim",
        fam="results/admixture/erythrura.fam"
    output:
        bed="results/admixture/erythrura.numchr.bed",
        bim="results/admixture/erythrura.numchr.bim",
        fam="results/admixture/erythrura.numchr.fam"
    conda:
        "envs/admixture.yaml"
    shell:
        r"""
        cp -f {input.bed} {output.bed}
        cp -f {input.fam} {output.fam}
        awk 'BEGIN{{OFS="\t"}} {{$1=1; print}}' {input.bim} > {output.bim}
        """

rule drop_outliers:
    input:
        bed="results/admixture/erythrura.numchr.bed",
        bim="results/admixture/erythrura.numchr.bim",
        fam="results/admixture/erythrura.numchr.fam",
        rm=OUTLIERS
    output:
        bed="results/admixture/erythrura.no_outliers.bed",
        bim="results/admixture/erythrura.no_outliers.bim",
        fam="results/admixture/erythrura.no_outliers.fam"
    conda:
        "envs/admixture.yaml"
    shell:
        r"""
        plink2 --bfile results/admixture/erythrura.numchr \
               --remove {input.rm} \
               --make-bed \
               --out results/admixture/erythrura.no_outliers
        """

rule admixture_friendly:
    input:
        bed="results/admixture/erythrura.no_outliers.bed",
        bim="results/admixture/erythrura.no_outliers.bim",
        fam="results/admixture/erythrura.no_outliers.fam"
    output:
        bed="results/admixture/erythrura.admix.bed",
        bim="results/admixture/erythrura.admix.bim",
        fam="results/admixture/erythrura.admix.fam"
    conda:
        "envs/admixture.yaml"
    shell:
        r"""
        plink2 --bfile results/admixture/erythrura.no_outliers \
               --max-alleles 2 \
               --snps-only just-acgt \
               --geno {GENO} \
               --maf {MAF} \
               --make-bed \
               --out results/admixture/erythrura.admix
        """

rule run_admixture:
    input:
        bed="results/admixture/erythrura.admix.bed",
        bim="results/admixture/erythrura.admix.bim",
        fam="results/admixture/erythrura.admix.fam"
    output:
        log="results/admixture/runs/K{K}/rep{rep}/run.log",
        Q="results/admixture/runs/K{K}/rep{rep}/erythrura.{K}.Q",
        P="results/admixture/runs/K{K}/rep{rep}/erythrura.{K}.P"
    threads: 8
    conda:
        "envs/admixture.yaml"
    params:
        seed=lambda wc: SEED_BASE + int(wc.K) * 1000 + int(wc.rep)
    shell:
        r"""
        set -euo pipefail
        d=results/admixture/runs/K{wildcards.K}/rep{wildcards.rep}
        mkdir -p "$d"
        cd "$d"

        ln -sf ../../../erythrura.admix.bed erythrura.bed
        ln -sf ../../../erythrura.admix.bim erythrura.bim
        ln -sf ../../../erythrura.admix.fam erythrura.fam

        admixture -j{threads} --seed={params.seed} --cv erythrura.bed {wildcards.K} \
          2>&1 | tee run.log
        """

# Optional convenience target (DON'T name it "all" if your main Snakefile already has one!)
rule admixture_done:
    input:
        expand("results/admixture/runs/K{K}/rep{rep}/run.log", K=K_LIST, rep=range(1, REPS+1)),
        "results/admixture/erythrura.admix.bed",
        "results/admixture/erythrura.admix.bim",
        "results/admixture/erythrura.admix.fam"

