#!/bin/bash
set -euo pipefail

base_dir="${REGA_ASSIGNMENTS_ROOT:-/project/zduren/durenlab/naqing_palmetto/REGA_separate/03_process_rega/module_assignments}"

for ct_dir in "${base_dir}"/*
do
    [ -d "$ct_dir" ] || continue

    ct=$(basename "$ct_dir")
    [[ "$ct" == Vasc_* ]] && continue
    assignment_file="${ct_dir}/RE_module_assignments.tsv"
    outdir="${ct_dir}/module_peaks"

    if [ ! -f "$assignment_file" ]; then
        echo "WARNING: missing ${assignment_file}"
        continue
    fi

    echo "Processing ${ct}"

    mkdir -p "$outdir"
    if compgen -G "${outdir}/${ct}_*.bed" > /dev/null; then
        echo "ERROR: module peak BEDs already exist for ${ct}; use a fresh version directory" >&2
        exit 1
    fi

    awk -v ct="$ct" -v outdir="$outdir" '
    BEGIN {
        FS=OFS="\t"
    }
    NR == 1 {
        next
    }
    $2 != "NA" && $2 != "" {
        peak=$1
        module=$2

        split(peak, a, ":")
        chr=a[1]

        split(a[2], b, "-")
        start=b[1]
        end=b[2]

        out = outdir "/" ct "_" module ".bed"
        print chr, start, end, peak >> out
    }
    ' "$assignment_file"

done
