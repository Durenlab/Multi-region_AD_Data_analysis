#!/bin/bash
#SBATCH --job-name=assign_module
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --time=48:00:00
#SBATCH --output=logs/assign_M.out
#SBATCH --error=logs/assign_M.err

ml anaconda3
source activate bench_r

# Pass the array ID to the R script
Rscript /project/zduren/durenlab/naqing_palmetto/REGA_separate/03_process_rega/scripts/09module_assign.r