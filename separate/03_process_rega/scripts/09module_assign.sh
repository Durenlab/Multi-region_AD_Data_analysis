#!/bin/bash
#SBATCH --job-name=assign_module
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --time=48:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
: "${REGA_SEPARATE_ROOT:?Set REGA_SEPARATE_ROOT}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

ml anaconda3
source activate bench_r

# Pass the array ID to the R script
Rscript "${SCRIPT_DIR}/09module_assign.r"
