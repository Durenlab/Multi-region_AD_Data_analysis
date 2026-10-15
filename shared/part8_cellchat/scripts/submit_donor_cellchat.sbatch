#!/bin/bash
# Draft only. Review this resource request and script before submission.
# One region per array task; at most two regions run concurrently.
# The earlier pooled CellChat tasks used up to 6 CPUs and 100 GB; donor-level
# objects are much smaller, so this draft requests 4 CPUs and 32 GB per task.
# Runtime is uncertain because each region contains 45-100 donors.

#SBATCH --job-name=CC_donor
#SBATCH --array=1-6%2
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/donor_%A_%a.out
#SBATCH --error=logs/donor_%A_%a.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REGA_SHARED_ROOT="${REGA_SHARED_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
CELLCHAT_UPSTREAM_ROOT="${CELLCHAT_UPSTREAM_ROOT:-/project/zduren/durenlab/naqing_palmetto/AD_REGA_region_wise/part8_cellchat}"
BASE_DIR="${CELLCHAT_UPSTREAM_ROOT}"
OUT_DIR="${REGA_SHARED_ROOT}/part8_cellchat"
SCRIPT="${SCRIPT_DIR}/run_donor_cellchat_by_region.R"

REGIONS=(AG EC HC MTC PFC TH)
REGION="${REGIONS[$((SLURM_ARRAY_TASK_ID - 1))]}"
INPUT_H5AD="${BASE_DIR}/01_region_inputs/${REGION}_cellchat_input.h5ad"
POOLED_RDS="${BASE_DIR}/02_cellchat_results/${REGION}/${REGION}_CellChat_final.rds"

cd "${OUT_DIR}"
ml anaconda3
source activate cellchat

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

test -s "${INPUT_H5AD}"
test -s "${POOLED_RDS}"
test -s "${SCRIPT}"

echo "Job ${SLURM_JOB_ID}; task ${SLURM_ARRAY_TASK_ID}; region ${REGION}"
echo "Input: ${INPUT_H5AD}"
echo "Interaction universe: ${POOLED_RDS}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}; memory request: 32G"
Rscript --vanilla "${SCRIPT}" \
  "${REGION}" "${INPUT_H5AD}" "${POOLED_RDS}" "${OUT_DIR}" 20 100
