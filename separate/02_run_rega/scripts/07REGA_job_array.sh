#!/usr/bin/env bash
#SBATCH --job-name=rega_array
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=h100:1
#SBATCH --time=24:00:00
#SBATCH --mem=500G
#SBATCH --array=0-35
#SBATCH --output=logs_rega/%A_%a.out
#SBATCH --error=logs_rega/%A_%a.err

set -euo pipefail

ml anaconda3

source activate graph_clean

export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

regions=("AG" "EC" "HC" "MTC" "PFC" "TH")
celltypes=("Ast" "Exc" "OPC" "Oli" "Inh" "Mic_Immune")
ct_index=$(( SLURM_ARRAY_TASK_ID % 6 ))
rg_index=$(( SLURM_ARRAY_TASK_ID / 6 ))

CELLTYPE=${celltypes[$ct_index]}
REGION=${regions[$rg_index]}

echo "Running job ${SLURM_ARRAY_TASK_ID}: ${REGION} - ${CELLTYPE}"

python /project/zduren/durenlab/naqing_palmetto/REGA_separate/02_run_rega/scripts/07REGA_full_gpu_run_array_save_mem_version.py --region ${REGION} --celltype ${CELLTYPE}
