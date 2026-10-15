#!/usr/bin/env bash
# Merge subcell-type cCREs into six nonvascular cell-type BED files.
set -euo pipefail
if (( $# != 2 )); then
  echo "Usage: bash merge_bedfiles.sh CCRE_DIR PREP_DIR" >&2
  exit 2
fi

ccre_dir=$1
prep_dir=$2
out_dir="$prep_dir/peak_resources"
mkdir -p "$out_dir"
shopt -s nullglob
for celltype in Ast Exc Inh Mic_Immune Oli OPC; do
  prefix=$celltype
  if [[ $celltype == Mic_Immune ]]; then prefix=Mic; fi
  files=("$ccre_dir/${prefix}_"*.bed)
  if (( ${#files[@]} == 0 )); then
    echo "No cCRE BED files for $celltype" >&2
    exit 1
  fi
  cat "${files[@]}" \
    | sort -k1,1 -k2,2n \
    | "${BEDTOOLS_BIN:-bedtools}" merge -i stdin \
    | awk -v OFS='\t' '{print $1,$2,$3,$1":"$2"-"$3}' \
    > "$out_dir/merged_${celltype}.bed"
  echo "$celltype: $out_dir/merged_${celltype}.bed"
done
