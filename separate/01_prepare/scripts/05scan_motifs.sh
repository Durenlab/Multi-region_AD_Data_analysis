#!/usr/bin/env bash
# Scan a nonvascular cell-type cCRE BED with HOMER.
set -euo pipefail
if (( $# != 3 )); then
  echo "Usage: bash 05scan_motifs.sh PREP_DIR MOTIF_FILE CELLTYPE" >&2
  exit 2
fi

prep_dir=$1
motif_file=$2
celltype=$3
case "$celltype" in
  Ast|Exc|Inh|Mic_Immune|Oli|OPC) ;;
  *) echo "Unknown or excluded cell type: $celltype" >&2; exit 2 ;;
esac
peak_bed="$prep_dir/peak_resources/merged_${celltype}.bed"
[[ -s $peak_bed ]] || { echo "Missing cCRE BED: $peak_bed" >&2; exit 1; }
out_dir="$prep_dir/motif_scan_results"
mkdir -p "$out_dir/$celltype"
findMotifsGenome.pl "$peak_bed" hg38 "$out_dir/$celltype" \
  -size given -mask -find "$motif_file" \
  > "$out_dir/${celltype}_scan_results.tsv"
