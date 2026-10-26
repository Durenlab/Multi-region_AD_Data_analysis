#!/usr/bin/env bash
#SBATCH --job-name=aging_gsea
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --array=1-36

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REGA_SHARED_ROOT="${REGA_SHARED_ROOT:-$(dirname "${SCRIPT_DIR}")}"
regions=(PFC TH AG EC HC MTC)
celltypes=(Ast Exc Inh Mic_Immune OPC Oli)

if [[ $# -eq 2 ]]; then
  REGION="$1"
  CELLTYPE="$2"
else
  task_id="${SLURM_ARRAY_TASK_ID:?Pass REGION CELLTYPE when running outside a Slurm array}"
  (( task_id >= 1 && task_id <= 36 )) || { echo "Invalid array task: ${task_id}" >&2; exit 2; }
  idx=$((task_id - 1))
  REGION="${regions[$((idx / 6))]}"
  CELLTYPE="${celltypes[$((idx % 6))]}"
fi

CONTEXT="${CELLTYPE}_${REGION}"
WBV="${REGA_SHARED_ROOT}/part3_process_REGA_output/09extracted_matrices_region_shareH/${REGION}/${CONTEXT}_WBV_genes_by_module.tsv"
GMT="${SCRIPT_DIR}/aging_up_gene_sets.gmt"
OUT="${SCRIPT_DIR}/GSEA_all_celltype_region/${CONTEXT}"
[[ -s "${WBV}" && -s "${GMT}" ]] || { echo "Missing WBV or GMT: ${WBV} ${GMT}" >&2; exit 1; }
mkdir -p "${OUT}"
ml anaconda3
source activate graph_clean
python "${SCRIPT_DIR}/03aging_gsea_from_wbv.py" \
  --wbv "${WBV}" --gmt "${GMT}" --celltype "${CELLTYPE}" \
  --context "${CONTEXT}" --out-root "${OUT}" \
  --perm 1000 --threads "${SLURM_CPUS_PER_TASK:-8}" --jitter
