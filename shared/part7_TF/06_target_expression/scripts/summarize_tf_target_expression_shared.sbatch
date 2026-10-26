#!/usr/bin/env bash
# Submit after the shared pair array completes successfully.
#SBATCH --job-name=TF_TG_shared_summary
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export AD_MANUSCRIPT_WORK_ROOT="${AD_MANUSCRIPT_WORK_ROOT:-$(realpath "${SCRIPT_DIR}/../../../../..")}"
export REGA_DATA_ROOT="${REGA_DATA_ROOT:-$(realpath "${SCRIPT_DIR}/../../../../../..")}"
RESULTS_ROOT="${TF_TARGET_RESULTS_ROOT:-${SCRIPT_DIR}/../results}"
RESULTS_ROOT="$(realpath -m "$RESULTS_ROOT")"
export TF_TARGET_RESULTS_ROOT="$RESULTS_ROOT"
ml anaconda3
source activate graph_clean
python "$SCRIPT_DIR/finalize_results.py" "$RESULTS_ROOT/cross_celltype"
