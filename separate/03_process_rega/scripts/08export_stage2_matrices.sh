#!/bin/bash
#SBATCH --job-name=extract_matrix
#SBATCH --cpus-per-task=4
#SBATCH --mem=800G
#SBATCH --time=48:00:00
#SBATCH --output=logs/extract.out
#SBATCH --error=logs/extract.err

mkdir -p logs

ml anaconda3 
source activate graph_clean
python /project/zduren/durenlab/naqing_palmetto/REGA_separate/03_process_rega/scripts/08export_stage2_matrices.py
