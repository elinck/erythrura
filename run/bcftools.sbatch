#!/bin/bash
#SBATCH --account=priority-ethanlinck
#SBATCH --partition=priority
#SBATCH --job-name=bcftools
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=6-00:00:00
#SBATCH --output=bcftools-%j.out
#SBATCH --error=bcftools-%j.err
#SBATCH --mail-user=ethan.linck@montana.edu
#SBATCH --mail-type=ALL

set -euo pipefail
source ~/.bashrc
mamba activate snakemake

# Run from your existing project directory (in /home)
snakemake --executor slurm \
  -s /home/k14m234/erythrura/Snakefile \
  -d /home/k14m234/erythrura \
  --workflow-profile /home/k14m234/erythrura/workflow-profiles/


