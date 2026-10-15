#!/bin/bash
#SBATCH --job-name=wbv_gsea
#SBATCH --output=logs/wbv_gsea_%A_%a.out
#SBATCH --error=logs/wbv_gsea_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --array=1-36
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REGA_SHARED_ROOT="${REGA_SHARED_ROOT:-$(dirname "${SCRIPT_DIR}")}"
REGA_SEPARATE_ROOT="${REGA_SEPARATE_ROOT:-/project/zduren/durenlab/naqing_palmetto/REGA_separate}"

set -euo pipefail

# =========================
# Paths
# =========================
BASE_DIR="${REGA_SHARED_ROOT}/part5_gsea"
SCRIPT="${BASE_DIR}/11batch_gsea_from_wbv.py"
GMT="${BASE_DIR}/KEGG_2026.Human.gmt"

WBV_ROOT="${REGA_SHARED_ROOT}/part3_process_REGA_output/09extracted_matrices_region_shareH"

OUT_ROOT="${BASE_DIR}/GSEA_from_WBV_KEGG2026_all_celltype_region"

# =========================
# Environment
# =========================
source activate graph_clean

mkdir -p logs
mkdir -p "${OUT_ROOT}"

# =========================
# Define all celltype-region contexts
# =========================

regions=(
  "PFC"
  "TH"
  "AG"
  "EC"
  "HC"
  "MTC"
)

celltypes=(
  "Ast"
  "Exc"
  "Inh"
  "Mic_Immune"
  "OPC"
  "Oli"
)

# 6 regions x 6 celltypes = 36 contexts

N_REGIONS=${#regions[@]}
N_CELLTYPES=${#celltypes[@]}
N_TOTAL=$((N_REGIONS * N_CELLTYPES))

if [[ "${SLURM_ARRAY_TASK_ID}" -lt 1 || "${SLURM_ARRAY_TASK_ID}" -gt "${N_TOTAL}" ]]; then
    echo "[ERROR] SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} outside valid range 1-${N_TOTAL}"
    exit 1
fi

IDX=$((SLURM_ARRAY_TASK_ID - 1))

REGION_IDX=$((IDX / N_CELLTYPES))
CELLTYPE_IDX=$((IDX % N_CELLTYPES))

REGION="${regions[$REGION_IDX]}"
CELLTYPE="${celltypes[$CELLTYPE_IDX]}"
CONTEXT="${CELLTYPE}_${REGION}"

# =========================
# Input / output paths
# =========================

WBV="${WBV_ROOT}/${REGION}/${CONTEXT}_WBV_genes_by_module.tsv"
CONTEXT_OUT="${OUT_ROOT}/${CONTEXT}"

if [[ ! -f "${WBV}" ]]; then
    echo "[ERROR] WBV file not found:"
    echo "  REGION=${REGION}"
    echo "  CELLTYPE=${CELLTYPE}"
    echo "  CONTEXT=${CONTEXT}"
    echo "  WBV=${WBV}"
    exit 1
fi

if [[ ! -f "${GMT}" ]]; then
    echo "[ERROR] GMT file not found: ${GMT}"
    exit 1
fi

mkdir -p "${CONTEXT_OUT}"

echo "[INFO] SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}"
echo "[INFO] REGION = ${REGION}"
echo "[INFO] CELLTYPE = ${CELLTYPE}"
echo "[INFO] CONTEXT = ${CONTEXT}"
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
