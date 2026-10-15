#!/usr/bin/env bash
# Refresh one pair or one TF without launching the full EXP008 array.
# Usage: sbatch run_exp008_pair.sbatch OPC_MTC MAF
#SBATCH --job-name=EXP008_pair
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/project/zduren/durenlab/naqing_palmetto/Alzheimer_project_explore/experiments/EXP008_results/logs/exp008_pair_%j.out
#SBATCH --error=/project/zduren/durenlab/naqing_palmetto/Alzheimer_project_explore/experiments/EXP008_results/logs/exp008_pair_%j.err
set -euo pipefail
cd /project/zduren/durenlab/naqing_palmetto/Alzheimer_project_explore
PAIR=${1:?Pass a celltype_region pair, e.g. OPC_MTC}
TF=${2:-ALL}
mkdir -p experiments/EXP008_results/logs
ml anaconda3
source activate bench_r
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript --vanilla experiments/exp008/scripts/plot_tf_tg_expression.R \
  "$PAIR" experiments/EXP008_results celltype_specific "$TF"
