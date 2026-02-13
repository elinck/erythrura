#!/bin/bash
#SBATCH --account=priority-ethanlinck
#SBATCH --partition=priority
#SBATCH --job-name=vcftools_in
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0-12:00:00
#SBATCH --output=log/vcftools_in-%j.out
#SBATCH --error=log/vcftools_in-%j.err
#SBATCH --mail-user=ethan.linck@montana.edu
#SBATCH --mail-type=ALL

set -euo pipefail
source ~/.bashrc
mamba activate vcftools

vcftools --gzvcf /home/k14m234/erythrura_assembly/results/GCF_005870125.1/QC/erythrura.pruned.vcf.gz \
  --keep /home/k14m234/erythrura/config/keep_samples.txt \
  --weir-fst-pop /home/k14m234/erythrura/config/trichroa_in.txt --weir-fst-pop /home/k14m234/erythrura/config/papuana_in.txt \
  --fst-window-size 50000 --fst-window-step 10000 \
  --out /home/k14m234/erythrura/results/vcftools/ingroup_fst_w50k_s10k
