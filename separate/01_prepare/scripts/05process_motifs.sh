#!/usr/bin/env bash
# Convert HOMER motif hits to peak–motif pairs for all six cell types.
set -euo pipefail
if (( $# != 1 )); then
  echo "Usage: bash 05process_motifs.sh PREP_DIR" >&2
  exit 2
fi

out_dir="$1/motif_scan_results"
for celltype in Ast Exc Inh Mic_Immune Oli OPC; do
  input="$out_dir/${celltype}_scan_results.tsv"
  [[ -s $input ]] || { echo "Missing motif scan: $input" >&2; exit 1; }
  awk -F'\t' 'NR > 1 {print $1 "\t" $4}' "$input" \
    | sort -V -k1,1 -k2,2 \
    | uniq > "$out_dir/${celltype}_peak_motif.txt"
  echo "$celltype: $out_dir/${celltype}_peak_motif.txt"
done
