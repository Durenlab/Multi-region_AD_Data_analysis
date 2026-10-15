#!/usr/bin/env bash
# Run the separate region-shared cross-celltype universe after manifest review.
#SBATCH --job-name=EXP008_cross
#SBATCH --array=0-35%2
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/project/zduren/durenlab/naqing_palmetto/Alzheimer_project_explore/experiments/EXP008_results/cross_celltype/logs/cross_%A_%a.out
#SBATCH --error=/project/zduren/durenlab/naqing_palmetto/Alzheimer_project_explore/experiments/EXP008_results/cross_celltype/logs/cross_%A_%a.err
set -euo pipefail
cd /project/zduren/durenlab/naqing_palmetto/Alzheimer_project_explore
MANIFEST=experiments/EXP008_results/cross_celltype/audit/cross_celltype_pairs.tsv
[[ -s "$MANIFEST" ]] || { echo "Run prepare_cross_celltype_manifest.py and review its audit first" >&2; exit 1; }
N=$(( $(wc -l < "$MANIFEST") - 1 ))
[[ "$N" -eq 36 ]] || { echo "Expected 36 non-Vasc region-consistent pairs; found $N" >&2; exit 1; }
awk -F '\t' 'NR>1 && $8 != "ready" {bad=1} END {exit bad}' "$MANIFEST" || { echo "Cross-celltype input preflight failed" >&2; exit 1; }
PAIR=$(sed -n "$((SLURM_ARRAY_TASK_ID + 2))p" "$MANIFEST" | cut -f1)
[[ -n "$PAIR" ]] || { echo "No pair for array index $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
ml anaconda3
source activate bench_r
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript --vanilla experiments/exp008/scripts/plot_tf_tg_expression.R \
  "$PAIR" experiments/EXP008_results/cross_celltype cross_celltype ALL
