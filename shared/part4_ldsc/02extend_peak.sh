#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REGA_SHARED_ROOT="${REGA_SHARED_ROOT:-$(dirname "${SCRIPT_DIR}")}"
REGA_SEPARATE_ROOT="${REGA_SEPARATE_ROOT:?Set REGA_SEPARATE_ROOT}"
set -euo pipefail

cd "${REGA_SHARED_ROOT}/part4_ldsc/"
mkdir -p peaks_extended_2kb/temp

BEDTOOLS="${BEDTOOLS:?Set BEDTOOLS}"
GENOME="${REGA_SEPARATE_ROOT}/04_ldsc/data/extended_peaks/hg19.autosomes.chrom.sizes"

IN_DIR="${REGA_SHARED_ROOT}/part4_ldsc/module_peaks_hg19"
OUT_DIR="${REGA_SHARED_ROOT}/part4_ldsc/peaks_extended_2kb"
TMP_DIR="${OUT_DIR}/temp"

echo -e "module_id\tinput\tautosome\tremoved\textended_merged" > "${OUT_DIR}/autosome_filter_counts.tsv"

for in_bed in "${IN_DIR}"/*.bed
do
    [ -e "$in_bed" ] || {
        echo "WARNING no bed files found in ${IN_DIR}"
        break
    }

    file_name=$(basename "$in_bed")
    module_id="${file_name%.bed}"

    out_bed="${OUT_DIR}/${module_id}.bed"
    tmp_bed="${TMP_DIR}/${module_id}.autosome.tmp.bed"

    echo "Extending ${module_id}"

    awk 'BEGIN{OFS="\t"} $1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ {print $0}' "$in_bed" > "$tmp_bed"

    n_in=$(wc -l < "$in_bed")
    n_auto=$(wc -l < "$tmp_bed")

    if [ "$n_auto" -eq 0 ]; then
        echo "WARNING ${module_id} has no autosome peaks after filtering"
        : > "$out_bed"
        n_ext=0
    else
        "$BEDTOOLS" slop \
          -i "$tmp_bed" \
          -g "$GENOME" \
          -b 2000 | \
        sort -k1,1 -k2,2n | \
        "$BEDTOOLS" merge -i - \
          > "$out_bed"

        n_ext=$(wc -l < "$out_bed")
    fi

    echo -e "${module_id}\t${n_in}\t${n_auto}\t$((n_in-n_auto))\t${n_ext}" \
      >> "${OUT_DIR}/autosome_filter_counts.tsv"

done
