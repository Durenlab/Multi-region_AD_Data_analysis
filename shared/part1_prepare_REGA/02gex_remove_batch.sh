#!/bin/bash
#SBATCH --job-name=limma_array
#SBATCH --array=1-36
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=limma_logs/limma_%A_%a.out
#SBATCH --error=limma_logs/limma_%A_%a.err

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"
export REGA_SHARED_ROOT="${REGA_SHARED_ROOT:-$(dirname "${SCRIPT_DIR}")}"

mkdir -p limma_logs

ml anaconda3
source activate bench_r

# Pass the SLURM array task ID to the R script
# Each task processes one line from rna_files_list.txt
Rscript 02gex_remove_batch.r ${SLURM_ARRAY_TASK_ID}