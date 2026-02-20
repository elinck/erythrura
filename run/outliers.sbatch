#!/bin/bash
#SBATCH --account=priority-ethanlinck
#SBATCH --partition=priority
#SBATCH --job-name=outliers
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=06-00:00:00
#SBATCH --output=log/outliers-%j.out
#SBATCH --error=log/outliers-%j.err
#SBATCH --mail-user=ethan.linck@montana.edu
#SBATCH --mail-type=ALL

set -euo pipefail
mkdir -p log

source ~/.bashrc
mamba activate snakemake

snakemake windows_annotated \
  -s /home/k14m234/erythrura/Snakefile \
  -d /home/k14m234/erythrura \
  --workflow-profile /home/k14m234/erythrura/workflow-profiles/ \
  --cores 4 \
  -j 8 \
  --rerun-incomplete \
  --latency-wait 60 \
  --use-conda \
  --rerun-incomplete
