#!/usr/bin/env bash
#SBATCH --job-name=TF_TG_shared
#SBATCH --array=0-35%2
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export AD_MANUSCRIPT_WORK_ROOT="${AD_MANUSCRIPT_WORK_ROOT:-$(realpath "${SCRIPT_DIR}/../../../../..")}"
export REGA_DATA_ROOT="${REGA_DATA_ROOT:-$(realpath "${SCRIPT_DIR}/../../../../../..")}"
RESULTS_ROOT="${TF_TARGET_RESULTS_ROOT:-${SCRIPT_DIR}/../results}"
RESULTS_ROOT="$(realpath -m "$RESULTS_ROOT")"
export TF_TARGET_RESULTS_ROOT="$RESULTS_ROOT"
MANIFEST="$RESULTS_ROOT/cross_celltype/audit/cross_celltype_pairs.tsv"
[[ -s "$MANIFEST" ]] || { echo "Run prepare_cross_celltype_manifest.py and review its audit first" >&2; exit 1; }
N=$(( $(wc -l < "$MANIFEST") - 1 ))
[[ "$N" -eq 36 ]] || { echo "Expected 36 nonvascular pairs; found $N" >&2; exit 1; }
awk -F '\t' 'NR>1 && $8 != "ready" {bad=1} END {exit bad}' "$MANIFEST" || { echo "Shared input preflight failed" >&2; exit 1; }
PAIR=$(sed -n "$((SLURM_ARRAY_TASK_ID + 2))p" "$MANIFEST" | cut -f1)
[[ -n "$PAIR" ]] || { echo "No pair for array index $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
ml anaconda3
source activate bench_r
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript --vanilla "$SCRIPT_DIR/plot_tf_tg_expression.R" "$PAIR" "$RESULTS_ROOT/cross_celltype" cross_celltype ALL
