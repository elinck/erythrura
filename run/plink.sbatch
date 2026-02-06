#!/bin/bash
#SBATCH --account=priority-ethanlinck
#SBATCH --partition=priority
#SBATCH --job-name=plink
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=0-12:00:00
#SBATCH --output=plink-%j.out
#SBATCH --error=plink-%j.err
#SBATCH --mail-user=ethan.linck@montana.edu
#SBATCH --mail-type=ALL

set -euo pipefail
source ~/.bashrc
mamba activate plink

plink2 \
  --vcf /home/k14m234/erythrura_assembly/results/GCF_005870125.1/QC/erythrura.pruned.vcf.gz \
  --max-alleles 2 \
  --snps-only \
  --mac 2 \
  --geno 0.2 \
  --thin-count 100000 \
  --recode vcf bgz \
  --threads 8 \
  --out /home/k14m234/erythrura/results/plink/erythrura_pca

tabix -p vcf /home/k14m234/erythrura/results/plink/erythrura_pca.vcf.gz

zgrep -vc '^#' /home/k14m234/erythrura/results/plink/erythrura_pca.vcf.gz > snp_count.txt




