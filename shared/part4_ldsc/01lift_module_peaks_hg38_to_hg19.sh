#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REGA_SHARED_ROOT="${REGA_SHARED_ROOT:-$(dirname "${SCRIPT_DIR}")}"
REGA_SEPARATE_ROOT="${REGA_SEPARATE_ROOT:-/project/zduren/durenlab/naqing_palmetto/REGA_separate}"
set -euo pipefail

# Activate conda env
ml anaconda3
source activate crossmap_env

base_dir="${REGA_SHARED_ROOT}"

chain_file="${LDSC_CHAIN_FILE:-/project/zduren/durenlab/palmetto/cham/util_data/position_mapping_hg/hg38ToHg19.over.chain.gz}"

input_root="${base_dir}/part4_ldsc/module_peaks_hg38"

out_dir="${base_dir}/part4_ldsc/module_peaks_hg19"
log_dir="${out_dir}/logs"

mkdir -p "${out_dir}" "${log_dir}"

echo "Input root: ${input_root}"
echo "Output dir: ${out_dir}"
echo "Chain file: ${chain_file}"
echo

# Find all celltype-region module_peaks directories:
# example:
# ${input_root}/AG/Ast_AG/module_peaks
find "${input_root}" -mindepth 3 -maxdepth 3 -type d -name "module_peaks" | sort | while read -r module_peak_dir; do

    celltype_region="$(basename "$(dirname "${module_peak_dir}")")"
    region="$(basename "$(dirname "$(dirname "${module_peak_dir}")")")"

    echo "Processing region=${region}, celltype_region=${celltype_region}"

    region_out_dir="${out_dir}"
    region_log_dir="${log_dir}/${region}/${celltype_region}"

    mkdir -p "${region_out_dir}" "${region_log_dir}"

    for m in $(seq 1 20); do
        in_bed="${module_peak_dir}/${celltype_region}_M${m}.bed"

        out_bed="${region_out_dir}/${celltype_region}_M${m}.bed"
        unmapped_bed="${region_out_dir}/${celltype_region}_M${m}.unmap.bed"
        log_file="${region_log_dir}/${celltype_region}_M${m}.log"

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

done

echo "Done."
echo "Lifted hg19 BED files are in:"
echo "${out_dir}"