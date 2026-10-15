#!/usr/bin/env bash
# Submit with an afterok dependency on run_exp008_cross_celltype.sbatch.
#SBATCH --job-name=EXP008_cross_summary
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00
#SBATCH --output=logs/finalize_%j.out
#SBATCH --error=logs/finalize_%j.err
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REGA_SHARED_ROOT="${REGA_SHARED_ROOT:-$(cd "${SCRIPT_DIR}/../../.." && pwd)}"
OUT_ROOT="${REGA_SHARED_ROOT}/part7_TF/06_target_expression/results"
cd "${OUT_ROOT}"
ml anaconda3
source activate graph_clean
python "${SCRIPT_DIR}/finalize_results.py" \
  "${OUT_ROOT}"
