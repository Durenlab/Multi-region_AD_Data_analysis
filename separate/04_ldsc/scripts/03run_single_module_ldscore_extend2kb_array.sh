#!/bin/bash
set -euo pipefail

ldsc_dir="/project/zduren/durenlab/naqing_palmetto/REGA_separate/04_ldsc/"
ctm_file="${REGA_LDSC_MODULE_FILE:-${ldsc_dir}/inputs/selected_modules.txt}"

PHASE3_DIR="/project/zduren/durenlab/palmetto/cham/GWAS/1000G_Phase3"
LDSC_PATH="/project/zduren/durenlab/palmetto/cham/GWAS/ldsc"
PLINK_DIR="${PHASE3_DIR}/1000G_EUR_Phase3_plink"
PRINT_SNP="${PHASE3_DIR}/list.txt"

task_id="${SLURM_ARRAY_TASK_ID}"

raw=$(sed -n "${task_id}p" "${ctm_file}")
raw=$(echo "${raw}" | tr -d '\r' | xargs)
module_id=$(echo "${raw}" | awk '{print $NF}')

if [ -z "${module_id}" ]; then
    echo "ERROR: empty module ID for task ${task_id}"
    exit 1
fi
if [[ "${module_id}" == Vasc_* ]]; then
    echo "ERROR: Vasc is excluded: ${module_id}" >&2
    exit 1
fi

if [[ ! "${module_id}" =~ _M[0-9]+$ ]]; then
    echo "ERROR: cannot parse module ID: ${module_id}"
    echo "Expected format: Ast_EC_M4"
    exit 1
fi

annot_dir="${REGA_LDSC_ANNOT_ROOT:-${ldsc_dir}/results/anno_phase3_single_module_extend2kb}/${module_id}"
out_dir="${REGA_LDSC_LDSCORE_ROOT:-${ldsc_dir}/results/ldscore_phase3_single_module_extend2kb}/${module_id}"

mkdir -p "${out_dir}"

echo "=================================================="
echo "Task ID: ${task_id}"
echo "Module ID: ${module_id}"
echo "Annotation dir: ${annot_dir}"
echo "Output dir: ${out_dir}"
echo "=================================================="

for chr in {1..22}
do
    annot_file="${annot_dir}/Base_CelltypePeakLD.${chr}.annot"
    bfile="${PLINK_DIR}/1000G.EUR.QC.${chr}"
    out_prefix="${out_dir}/Base_CelltypePeakLD.${chr}"

    if [ ! -f "${annot_file}" ]; then
        echo "ERROR: missing annotation file: ${annot_file}"
        exit 1
    fi

    if [ -f "${out_prefix}.l2.ldscore.gz" ] && \
       [ -f "${out_prefix}.l2.M" ] && \
       [ -f "${out_prefix}.l2.M_5_50" ]; then
        echo "Output exists for ${module_id} chr${chr}; skipping."
        continue
    fi

    python "${LDSC_PATH}/ldsc.py" \
      --print-snps "${PRINT_SNP}" \
      --l2 \
      --bfile "${bfile}" \
      --ld-wind-cm 1 \
      --annot "${annot_file}" \
      --out "${out_prefix}"
done

echo "Done LD scores for ${module_id}"
