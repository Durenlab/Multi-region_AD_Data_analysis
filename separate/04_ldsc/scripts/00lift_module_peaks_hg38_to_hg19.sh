#!/bin/bash
set -euo pipefail

# Activate conda env
source activate crossmap_env
base_dir="/project/zduren/durenlab/naqing_palmetto/REGA_separate/"
list_file="${REGA_LDSC_CONTEXT_FILE:-${base_dir}/04_ldsc/inputs/selected_context.txt}"

chain_file="/project/zduren/durenlab/palmetto/cham/util_data/position_mapping_hg/hg38ToHg19.over.chain.gz"

assignment_dir="${REGA_ASSIGNMENTS_ROOT:-${base_dir}/03_process_rega/module_assignments}"
out_dir="${REGA_LDSC_HG19_ROOT:-${base_dir}/04_ldsc/data/module_peaks_hg19}"
log_dir="${out_dir}/logs"

mkdir -p "${out_dir}" "${log_dir}"

echo "Input list: ${list_file}"
echo "Output dir: ${out_dir}"
echo "Chain file: ${chain_file}"
echo

while read -r celltype_region; do
    # skip empty lines
    [ -z "${celltype_region}" ] && continue
    [[ "${celltype_region}" == Vasc_* ]] && continue

    echo "Processing ${celltype_region}"

    for m in $(seq 1 20); do
        if [[ -n "${REGA_LDSC_MODULE_FILE:-}" ]] && ! grep -Fxq "${celltype_region}_M${m}" "$REGA_LDSC_MODULE_FILE"; then
            continue
        fi
        in_bed="${assignment_dir}/${celltype_region}/module_peaks/${celltype_region}_M${m}.bed"
        out_bed="${out_dir}/${celltype_region}_M${m}.bed"
        unmapped_bed="${out_dir}/${celltype_region}_M${m}.unmap.bed"
        log_file="${log_dir}/${celltype_region}_M${m}.log"

        if [ ! -s "${in_bed}" ]; then
            echo "  Missing or empty: ${in_bed}" | tee -a "${log_dir}/missing_files.log"
            continue
        fi

        echo "  Lifting M${m}"

        CrossMap bed \
            "${chain_file}" \
            "${in_bed}" \
            "${out_bed}" \
            > "${log_file}" 2>&1

        # CrossMap often creates unmapped file as output.bed.unmap
        if [ -f "${out_bed}.unmap" ]; then
            mv "${out_bed}.unmap" "${unmapped_bed}"
        fi
    done

    echo
done < "${list_file}"

echo "Done."
echo "Lifted hg19 BED files are in:"
echo "${out_dir}"
