#!/usr/bin/env bash
#SBATCH --job-name=rega_array
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=h100:1
#SBATCH --time=24:00:00
#SBATCH --mem=500G
#SBATCH --array=0-35
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
: "${REGA_SEPARATE_ROOT:?Set REGA_SEPARATE_ROOT}"

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

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

python "${SCRIPT_DIR}/07REGA_full_gpu_run_array_save_mem_version.py" --region "${REGION}" --celltype "${CELLTYPE}"
