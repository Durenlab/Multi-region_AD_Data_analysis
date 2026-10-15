#!/usr/bin/env bash
# Intersect each cell-type cCRE set with coding-gene TSS windows.
set -euo pipefail
if (( $# != 2 )); then
  echo "Usage: bash 05TF_RE_prepare.sh TSS_WINDOWS_BED PREP_DIR" >&2
  exit 2
fi

tss_bed=$1
prep_dir=$2
out_dir="$prep_dir/re_tg_matrices"
mkdir -p "$out_dir"
for celltype in Ast Exc Inh Mic_Immune Oli OPC; do
  peak_bed="$prep_dir/peak_resources/merged_${celltype}.bed"
  [[ -s $peak_bed ]] || { echo "Missing cCRE BED: $peak_bed" >&2; exit 1; }
  "${BEDTOOLS_BIN:-bedtools}" intersect -a "$tss_bed" -b "$peak_bed" -wa -wb \
    | uniq > "$out_dir/${celltype}_RE_TG_matrix.bed"
  echo "$celltype: $out_dir/${celltype}_RE_TG_matrix.bed"
done
