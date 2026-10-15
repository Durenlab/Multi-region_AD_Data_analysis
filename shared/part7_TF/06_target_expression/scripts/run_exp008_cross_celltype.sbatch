#!/usr/bin/env bash
# Run the separate region-shared cross-celltype universe after manifest review.
#SBATCH --job-name=EXP008_cross
#SBATCH --array=0-35%2
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/cross_%A_%a.out
#SBATCH --error=logs/cross_%A_%a.err
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REGA_SHARED_ROOT="${REGA_SHARED_ROOT:-$(cd "${SCRIPT_DIR}/../../.." && pwd)}"
OUT_ROOT="${REGA_SHARED_ROOT}/part7_TF/06_target_expression/results"
cd "${OUT_ROOT}"
MANIFEST="${OUT_ROOT}/audit/cross_celltype_pairs.tsv"
[[ -s "$MANIFEST" ]] || { echo "Run prepare_cross_celltype_manifest.py and review its audit first" >&2; exit 1; }
N=$(( $(wc -l < "$MANIFEST") - 1 ))
[[ "$N" -eq 36 ]] || { echo "Expected 36 non-Vasc region-consistent pairs; found $N" >&2; exit 1; }
awk -F '\t' 'NR>1 && $8 != "ready" {bad=1} END {exit bad}' "$MANIFEST" || { echo "Cross-celltype input preflight failed" >&2; exit 1; }
PAIR=$(sed -n "$((SLURM_ARRAY_TASK_ID + 2))p" "$MANIFEST" | cut -f1)
[[ -n "$PAIR" ]] || { echo "No pair for array index $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
ml anaconda3
source activate bench_r
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript --vanilla "${SCRIPT_DIR}/plot_tf_tg_expression.R" \
  "$PAIR" "${OUT_ROOT}" cross_celltype ALL
