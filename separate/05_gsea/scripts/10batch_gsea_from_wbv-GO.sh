#!/bin/bash
#SBATCH --job-name=wbv_gsea
#SBATCH --output=logs/wbv_gsea_%A_%a.out
#SBATCH --error=logs/wbv_gsea_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --array=1-18

set -euo pipefail

# =========================
# Paths
# =========================
BASE_DIR="/project/zduren/durenlab/naqing_palmetto/REGA_separate/05_gsea/"
SCRIPT="${BASE_DIR}/scripts/10batch_gsea_from_wbv.py"
LIST_FILE="${BASE_DIR}/inputs/selected_context.txt"
GMT="${BASE_DIR}/data/GO_Biological_Process_2026_Human.gmt"

WBV_ROOT="/project/zduren/durenlab/naqing_palmetto/REGA_separate/03_process_rega/matrices"
OUT_ROOT="${BASE_DIR}/results/go_bp"

# =========================
# Environment
# =========================

# Modify this to your actual conda setup if needed

ml anaconda3
source activate graph_clean

mkdir -p logs
mkdir -p "${OUT_ROOT}"

# =========================
# Get context from list
# =========================

CONTEXT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${LIST_FILE}")

if [[ -z "${CONTEXT}" ]]; then
    echo "[ERROR] No context found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "[INFO] SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}"
echo "[INFO] CONTEXT = ${CONTEXT}"

# =========================
# Input / output paths
# =========================

WBV="${WBV_ROOT}/${CONTEXT}/${CONTEXT}_WBV_genes_by_module.tsv"
CONTEXT_OUT="${OUT_ROOT}/${CONTEXT}"

if [[ ! -f "${WBV}" ]]; then
    echo "[ERROR] WBV file not found: ${WBV}"
    exit 1
fi

if [[ ! -f "${GMT}" ]]; then
    echo "[ERROR] GMT file not found: ${GMT}"
    exit 1
fi

mkdir -p "${CONTEXT_OUT}"

echo "[INFO] WBV = ${WBV}"
echo "[INFO] OUT = ${CONTEXT_OUT}"
echo "[INFO] GMT = ${GMT}"

# =========================
# Run GSEA
# =========================

python "${SCRIPT}" \
  --wbv "${WBV}" \
  --out-root "${CONTEXT_OUT}" \
  --gmt "${GMT}" \
  --perm 1000 \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --upper \
  --jitter

echo "[DONE] Finished GSEA for ${CONTEXT}"
