#!/bin/bash
set -euo pipefail

ldsc_dir="/project/zduren/durenlab/naqing_palmetto/REGA_separate/04_ldsc"

ctm_file="${REGA_LDSC_MODULE_FILE:-${ldsc_dir}/inputs/selected_modules.txt}"
peak_dir="${REGA_LDSC_EXTENDED_ROOT:-${ldsc_dir}/data/extended_peaks/extend_bed}"
out_root="${REGA_LDSC_ANNOT_ROOT:-${ldsc_dir}/results/anno_phase3_single_module_extend2kb}"

task_id="${SLURM_ARRAY_TASK_ID}"
#task_id=24

raw=$(sed -n "${task_id}p" "${ctm_file}")
raw=$(echo "${raw}" | tr -d '\r' | xargs)

if [ -z "${raw}" ]; then
    echo "ERROR: empty module ID for task ${task_id}"
    exit 1
fi
if [[ "${raw}" == Vasc_* ]]; then
    echo "ERROR: Vasc is excluded: ${raw}" >&2
    exit 1
fi

if [[ "${raw}" =~ _M[0-9]+$ ]]; then
    module=$(echo "${raw}" | grep -oE 'M[0-9]+$')
    ct=$(echo "${raw}" | sed -E 's/_M[0-9]+$//')
    module_id="${raw}"
else
    echo "ERROR: cannot parse module ID: ${raw}"
    echo "Expected format: Ast_EC_M4"
    exit 1
fi

peak_bed="${peak_dir}/${module_id}.bed"
out_dir="${out_root}/${module_id}"

echo "=================================================="
echo "Task ID: ${task_id}"
echo "Raw ID: ${raw}"
echo "Celltype-region: ${ct}"
echo "Module: ${module}"
echo "Module ID: ${module_id}"
echo "Peak BED: ${peak_bed}"
echo "=================================================="

if [ ! -f "${peak_bed}" ]; then
    echo "ERROR: missing peak BED: ${peak_bed}"
    exit 1
fi

bash "${ldsc_dir}/scripts/make_single_module_annot_phase3.sh" \
  "${module_id}" \
  "${peak_bed}" \
  "${out_dir}" \
  "${module_id}"
