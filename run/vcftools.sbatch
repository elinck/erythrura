#!/bin/bash
#SBATCH --account=priority-ethanlinck
#SBATCH --partition=priority
#SBATCH --job-name=vcftools
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0-12:00:00
#SBATCH --output=vcftools-%j.out
#SBATCH --error=vcftools-%j.err
#SBATCH --mail-user=ethan.linck@montana.edu
#SBATCH --mail-type=ALL

set -euo pipefail
source ~/.bashrc
mamba activate vcftools

vcftools --gzvcf /home/k14m234/erythrura_assembly/results/GCF_005870125.1/QC/erythrura.pruned.vcf.gz \
  --weir-fst-pop /home/k14m234/erythrura/config/trichroa.txt --weir-fst-pop /home/k14m234/erythrura/config/papuana.txt \
  --fst-window-size 50000 --fst-window-step 10000 \
  --out /home/k14m234/erythrura/results/vcftools/species_fst_w50k_s10k
