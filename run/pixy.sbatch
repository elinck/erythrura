#!/bin/bash
#SBATCH --account=priority-ethanlinck
#SBATCH --partition=priority
#SBATCH --job-name=pixy
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0-06:00:00
#SBATCH --output=pixy-%j.out
#SBATCH --error=pixy-%j.err
#SBATCH --mail-user=ethan.linck@montana.edu
#SBATCH --mail-type=ALL

set -euo pipefail
source ~/.bashrc
mamba activate pixy
tabix -f -p vcf snparcher/variants/erythrura.pruned.vcf.gz

pixy --stats fst \
--vcf snparcher/variants/erythrura.pruned.vcf.gz \
--populations config/pixy_pop.txt \
--window_size 10000 \
--n_cores 4 \
--output_folder results \
--output_prefix pixy_output\
--bypass_invariant_check 'yes'

