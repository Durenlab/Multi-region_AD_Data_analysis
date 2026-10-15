#!/usr/bin/env bash
# Submit with an afterok dependency on run_exp008_cross_celltype.sbatch.
#SBATCH --job-name=EXP008_cross_summary
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00
#SBATCH --output=/project/zduren/durenlab/naqing_palmetto/Alzheimer_project_explore/experiments/EXP008_results/cross_celltype/logs/finalize_%j.out
#SBATCH --error=/project/zduren/durenlab/naqing_palmetto/Alzheimer_project_explore/experiments/EXP008_results/cross_celltype/logs/finalize_%j.err
set -euo pipefail
cd /project/zduren/durenlab/naqing_palmetto/Alzheimer_project_explore
ml anaconda3
source activate graph_clean
python experiments/exp008/scripts/finalize_results.py \
  experiments/EXP008_results/cross_celltype
