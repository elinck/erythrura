#!/bin/bash
#SBATCH --account=priority-ethanlinck
#SBATCH --partition=priority
#SBATCH --job-name=admixture
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01-00:00:00
#SBATCH --output=log/admixture-%j.out
#SBATCH --error=log/admixture-%j.err
#SBATCH --mail-user=ethan.linck@montana.edu
#SBATCH --mail-type=ALL

set -euo pipefail
mkdir -p log

source ~/.bashrc
mamba activate snakemake

snakemake admixture \
  -s /home/k14m234/erythrura/Snakefile \
  -d /home/k14m234/erythrura \
  --workflow-profile /home/k14m234/erythrura/workflow-profiles/ \
  --cores 4 \
  -j 50 \
  --rerun-incomplete \
  --latency-wait 60 \
  --use-conda \
  --rerun-incomplete
