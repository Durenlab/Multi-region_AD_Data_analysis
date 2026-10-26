#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REGA_SHARED_ROOT="${REGA_SHARED_ROOT:-$(dirname "$(dirname "${SCRIPT_DIR}")")}"
REGA_SEPARATE_ROOT="${REGA_SEPARATE_ROOT:?Set REGA_SEPARATE_ROOT}"
set -euo pipefail

ldsc_dir="${REGA_SHARED_ROOT}/part4_ldsc"

ctm_file="${ldsc_dir}/selected_modules_FC.txt"
peak_dir="${ldsc_dir}/peaks_extended_2kb"
out_root="${ldsc_dir}/annotations_phase3"

task_id="${SLURM_ARRAY_TASK_ID}"
#task_id=24

raw=$(sed -n "${task_id}p" "${ctm_file}")
raw=$(echo "${raw}" | tr -d '\r' | xargs)

if [ -z "${raw}" ]; then
    echo "ERROR: empty module ID for task ${task_id}"
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