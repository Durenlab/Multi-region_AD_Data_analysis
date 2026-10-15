#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REGA_SHARED_ROOT="${REGA_SHARED_ROOT:-$(dirname "$(dirname "${SCRIPT_DIR}")")}"
REGA_SEPARATE_ROOT="${REGA_SEPARATE_ROOT:-/project/zduren/durenlab/naqing_palmetto/REGA_separate}"
set -euo pipefail

ldsc_dir="${REGA_SHARED_ROOT}/part4_ldsc/"
ctm_file="${ldsc_dir}/selected_modules_FC.txt"

task_id="${SLURM_ARRAY_TASK_ID}"

raw=$(sed -n "${task_id}p" "${ctm_file}")
raw=$(echo "${raw}" | tr -d '\r' | xargs)
module_id=$(echo "${raw}" | awk '{print $NF}')

if [ -z "${module_id}" ]; then
    echo "ERROR: empty module ID for task ${task_id}"
    exit 1
fi

if [[ ! "${module_id}" =~ _M[0-9]+$ ]]; then
    echo "ERROR: cannot parse module ID: ${module_id}"
    echo "Expected format: Ast_EC_M4"
    exit 1
fi

module=$(echo "${module_id}" | grep -oE 'M[0-9]+$')
ct=$(echo "${module_id}" | sed -E 's/_M[0-9]+$//')

echo "=================================================="
echo "Task ID: ${task_id}"
echo "Module ID: ${module_id}"
echo "Celltype-region: ${ct}"
echo "Module: ${module}"
echo "=================================================="

Phase3_dir="${LDSC_PHASE3_DIR:-/project/zduren/durenlab/palmetto/cham/GWAS/1000G_Phase3}"
ldsc_path="${LDSC_CODE_DIR:-/project/zduren/durenlab/palmetto/cham/GWAS/ldsc}"

GWAS_SUMSTATS="${REGA_SEPARATE_ROOT}/external_inputs/AD_Bellenguez.sumstats.gz"

base_L2="${Phase3_dir}/baselineLD/baselineLD."
custom_L2="${ldsc_dir}/ld_scores_phase3/${module_id}/Base_CelltypePeakLD."
Weight_dir="${Phase3_dir}/1000G_Phase3_weights_hm3_no_MHC"
Frq_dir="${Phase3_dir}/1000G_Phase3_frq"

annot_dir="${ldsc_dir}/annotations_phase3/${module_id}"
ldscore_dir="${ldsc_dir}/ld_scores_phase3/${module_id}"

out_dir="${ldsc_dir}/heritability_results/${module_id}"
out_prefix="${out_dir}/AD_Bellenguez_extend2kb_${module_id}"

mkdir -p "${out_dir}"

if [ ! -f "${GWAS_SUMSTATS}" ]; then
    echo "ERROR: missing GWAS sumstats: ${GWAS_SUMSTATS}"
    exit 1
fi

if [ ! -f "${ldscore_dir}/Base_CelltypePeakLD.1.l2.ldscore.gz" ]; then
    echo "ERROR: missing LD score file: ${ldscore_dir}/Base_CelltypePeakLD.1.l2.ldscore.gz"
    exit 1
fi

# LDSC --overlap-annot expects .annot files next to the LD score prefix.
for chr in {1..22}
do
    if [ ! -f "${annot_dir}/Base_CelltypePeakLD.${chr}.annot" ]; then
        echo "ERROR: missing annotation file: ${annot_dir}/Base_CelltypePeakLD.${chr}.annot"
        exit 1
    fi

    ln -sf "${annot_dir}/Base_CelltypePeakLD.${chr}.annot" \
      "${ldscore_dir}/Base_CelltypePeakLD.${chr}.annot"
done

python "${ldsc_path}/ldsc.py" \
  --h2 "${GWAS_SUMSTATS}" \
  --ref-ld-chr "${base_L2},${custom_L2}" \
  --w-ld-chr "${Weight_dir}/weights.hm3_noMHC." \
  --overlap-annot \
  --frqfile-chr "${Frq_dir}/1000G.EUR.QC." \
  --print-coefficients \
  --out "${out_prefix}"

echo "Done h2 for ${module_id}"