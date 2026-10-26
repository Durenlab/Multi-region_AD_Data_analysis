#!/bin/bash
set -euo pipefail

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <module_id> <peak_bed> <out_dir> <annotation_name>"
    echo "Example: $0 Ast_EC_M4 module_peaks_hg19/Ast_EC_M4.bed anno_phase3_single_module/Ast_EC_M4 Ast_EC_M4"
    exit 1
fi

module_id=$1
peak_bed=$2
out_dir=$3
annot_name=$4

bim_dir="${LDSC_BIM_DIR:?Set LDSC_BIM_DIR}"

mkdir -p "${out_dir}" "${out_dir}/tmp"

if [ ! -f "${peak_bed}" ]; then
    echo "ERROR: missing peak bed: ${peak_bed}"
    exit 1
fi

if [[ "${peak_bed}" == *".unmap.bed" ]]; then
    echo "ERROR: input peak bed is an unmapped file: ${peak_bed}"
    exit 1
fi

echo "Module ID: ${module_id}"
echo "Peak BED: ${peak_bed}"
echo "Output dir: ${out_dir}"
echo "Annotation name: ${annot_name}"

for chr in $(seq 1 22)
do
    echo "Processing chr${chr}"

    snp_bed="${out_dir}/tmp/snps.chr${chr}.bed"
    base_annot="${out_dir}/tmp/base.chr${chr}.annot"
    module_tsv="${out_dir}/tmp/${module_id}.chr${chr}.tsv"
    out_annot="${out_dir}/Base_CelltypePeakLD.${chr}.annot"

    sed 's/^/chr/' "${bim_dir}/1000G.EUR.QC.${chr}.bim" | \
    awk 'BEGIN{OFS="\t"} {print $1,$4-1,$4,$2,$3}' \
    > "${snp_bed}"

    awk 'BEGIN{OFS="\t"; print "CHR","BP","SNP","CM"} {print $1,$4,$2,$3}' \
      "${bim_dir}/1000G.EUR.QC.${chr}.bim" \
      > "${base_annot}"

    "${BEDTOOLS:?Set BEDTOOLS}" intersect \
      -a "${snp_bed}" \
      -b "${peak_bed}" \
      -c | \
    cut -f 6 | \
    awk -v name="${annot_name}" 'BEGIN{print name}{print $0}' \
      > "${module_tsv}"

    paste "${base_annot}" "${module_tsv}" \
      > "${out_annot}"
done

echo "Done annotation for ${module_id}"