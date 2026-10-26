#!/bin/bash
#SBATCH --job-name=extract_matrix
#SBATCH --cpus-per-task=4
#SBATCH --mem=800G
#SBATCH --time=48:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
: "${REGA_SEPARATE_ROOT:?Set REGA_SEPARATE_ROOT}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p logs

ml anaconda3 
source activate graph_clean
python "${SCRIPT_DIR}/08export_stage2_matrices.py"
