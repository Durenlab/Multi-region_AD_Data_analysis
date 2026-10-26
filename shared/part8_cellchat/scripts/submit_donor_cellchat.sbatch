#!/bin/bash
# Draft only. Review this resource request and script before submission.
# One region per array task; at most two regions run concurrently.
# Donor-level network fits; review resource requests before submission.

#SBATCH --job-name=CC_donor
#SBATCH --array=1-6%2
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
: "${CELLCHAT_WORK_DIR:?Set CELLCHAT_WORK_DIR to a results directory outside the code release}"
OUT_DIR="${CELLCHAT_WORK_DIR}"
SCRIPT="${SCRIPT_DIR}/run_donor_cellchat_by_region.R"

REGIONS=(AG EC HC MTC PFC TH)
REGION="${REGIONS[$((SLURM_ARRAY_TASK_ID - 1))]}"
INPUT_H5AD="${OUT_DIR}/01_region_inputs/${REGION}_cellchat_input.h5ad"
CANDIDATE_RDS="${OUT_DIR}/candidate_universe/${REGION}_candidate_lr.rds"

cd "${OUT_DIR}"
ml anaconda3
source activate cellchat

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

test -s "${INPUT_H5AD}"
test -s "${CANDIDATE_RDS}"
test -s "${SCRIPT}"

echo "Job ${SLURM_JOB_ID}; task ${SLURM_ARRAY_TASK_ID}; region ${REGION}"
echo "Input: ${INPUT_H5AD}"
echo "Interaction universe: ${CANDIDATE_RDS}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}; memory request: 32G"
Rscript --vanilla "${SCRIPT}" \
  "${REGION}" "${INPUT_H5AD}" "${CANDIDATE_RDS}" "${OUT_DIR}" 20 100
