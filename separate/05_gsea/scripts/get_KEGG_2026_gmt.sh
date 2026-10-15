#!/usr/bin/env bash
# GMT source: Ma'ayan Lab Enrichr human gene-set libraries.
# KEGG: KEGG_2026, downloaded directly from
#   https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2026
# GO Biological Process: GO_Biological_Process_2026, retrieved from Enrichr
#   through GSEApy get_library(name="GO_Biological_Process_2026", organism="Human").
#   Enrichr library URL:
#   https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2026
# Run from the directory where the two GMT files should be saved.

set -euo pipefail

wget -O KEGG_2026.Human.gmt \
    'https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2026'

python - <<'PY'
import gseapy as gp

library = "GO_Biological_Process_2026"
organism = "Human"
out_gmt = "GO_Biological_Process_2026_Human.gmt"

go_bp = gp.get_library(name=library, organism=organism)

with open(out_gmt, "w") as f:
    for term, genes in go_bp.items():
        f.write(term + "\t" + "na" + "\t" + "\t".join(genes) + "\n")

print(f"Saved {len(go_bp)} gene sets to {out_gmt}")
PY
