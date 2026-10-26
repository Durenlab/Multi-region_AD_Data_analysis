#!/usr/bin/env bash
#SBATCH --job-name=rega_region
#SBATCH --array=1-6
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=h100:1
#SBATCH --time=72:00:00
#SBATCH --mem=600G
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REGA_SHARED_ROOT="${REGA_SHARED_ROOT:-$(dirname "${SCRIPT_DIR}")}"
REGA_SEPARATE_ROOT="${REGA_SEPARATE_ROOT:?Set REGA_SEPARATE_ROOT}"

set -euo pipefail

mkdir -p "${SCRIPT_DIR}/logs_rega"

ml anaconda3/2023.09-0
source activate graph_clean
cd "${SCRIPT_DIR}"

export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

REGIONS=(PFC TH AG EC HC MTC)
REGION=${REGIONS[$((SLURM_ARRAY_TASK_ID - 1))]}

CELL_KEYS="Ast_${REGION},Exc_${REGION},OPC_${REGION},Oli_${REGION},Inh_${REGION},Mic_Immune_${REGION}"

echo "=========================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-NA}"
echo "Region: ${REGION}"
echo "Cell keys: ${CELL_KEYS}"
echo "=========================================="

python "${SCRIPT_DIR}/07run_REGA_blocks_multicell_shareH.py" \
  --blocks-root "${REGA_SHARED_ROOT}/part1_prepare_REGA/05REGA_blocks_region_consistent" \
  --out-root "${REGA_SHARED_ROOT}/part2_run_REGA/07REGA_region_shareH_output" \
  --run-name AD_region_shareH_${REGION}_K20 \
  --cells "${CELL_KEYS}" \
  --K 20 \
  --inner1 2000 \
  --outer1 30 \
  --inner2 500 \
  --outer2 3000

echo "=========================================="
echo "Finished region: ${REGION}"
echo "End: $(date)"
echo "=========================================="
