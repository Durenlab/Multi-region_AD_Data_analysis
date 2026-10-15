#!/usr/bin/env bash
# Run prepare_manifest.py once and review its audit tables before submitting.
#SBATCH --job-name=EXP008_TGexp
#SBATCH --array=0-41%2
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/project/zduren/durenlab/naqing_palmetto/Alzheimer_project_explore/experiments/EXP008_results/logs/exp008_%A_%a.out
#SBATCH --error=/project/zduren/durenlab/naqing_palmetto/Alzheimer_project_explore/experiments/EXP008_results/logs/exp008_%A_%a.err
set -euo pipefail
cd /project/zduren/durenlab/naqing_palmetto/Alzheimer_project_explore
MANIFEST=experiments/EXP008_results/audit/pairs.tsv
[[ -s "$MANIFEST" ]] || { echo "Run experiments/exp008/scripts/prepare_manifest.py first" >&2; exit 1; }
N=$(( $(wc -l < "$MANIFEST") - 1 ))
[[ "$N" -eq 42 ]] || { echo "Expected 42 expression pairs; found $N" >&2; exit 1; }
awk -F '\t' 'NR>1 && $8 != "ready" {bad=1} END {exit bad}' "$MANIFEST" || { echo "Pair input preflight failed" >&2; exit 1; }
PAIR=$(sed -n "$((SLURM_ARRAY_TASK_ID + 2))p" "$MANIFEST" | cut -f1)
[[ -n "$PAIR" ]] || { echo "No pair for array index $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
ml anaconda3
source activate bench_r
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript --vanilla experiments/exp008/scripts/plot_tf_tg_expression.R \
  "$PAIR" experiments/EXP008_results
