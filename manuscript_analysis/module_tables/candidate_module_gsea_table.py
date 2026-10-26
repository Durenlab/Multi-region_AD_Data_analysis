#!/usr/bin/env python3
"""Filter positive GSEA hit tables to curated candidate modules."""
from pathlib import Path
import os

import pandas as pd


ROOT = Path(os.environ.get("AD_MANUSCRIPT_WORK_ROOT", Path(__file__).resolve().parents[3])).expanduser().resolve()
BASE = Path(os.environ.get("REGA_DATA_ROOT", ROOT.parent)).expanduser().resolve()

SEPARATE = BASE / "REGA_separate"
SHARED = BASE / "REGA_shared"
CELLTYPE_CANDIDATES = set((SEPARATE / "02_run_rega/pathology_association/selected_module.txt").read_text().split())
CROSS_CANDIDATES = set((SHARED / "part2_run_REGA/08pathology_tests/selected_module.txt").read_text().split())

SOURCES = [
    {
        "name": "celltype_specific_GOBP",
        "path": SEPARATE / "05_gsea/results/go_bp/GO_positive_p005_hits.tsv",
        "scope": "celltype",
        "output": ROOT / "manuscript/supplementary_table2-celltype_specific_GOBP.tsv",
        "mirror": ROOT / "manuscript/supp/Supplementary_Table_celltype_specific_GOBP_GSEA.tsv",
    },
    {
        "name": "celltype_specific_KEGG",
        "path": SEPARATE / "05_gsea/results/kegg/KEGG_positive_p005_hits.tsv",
        "scope": "celltype",
        "output": ROOT / "manuscript/supplementary_table2-celltype_specific_KEGG.tsv",
        "mirror": ROOT / "manuscript/supp/Supplementary_Table_celltype_specific_KEGG_GSEA.tsv",
    },
    {
        "name": "cross_celltype_GOBP",
        "path": SHARED / "part5_gsea/visualize/region_wise_GOBP_positive_p005_hits.tsv",
        "scope": "cross",
        "output": ROOT / "manuscript/supplementary_table2-cross_celltype_GOBP.tsv",
        "mirror": ROOT / "manuscript/supp/Supplementary_Table_cross_celltype_GOBP_GSEA.tsv",
    },
    {
        "name": "cross_celltype_KEGG",
        "path": SHARED / "part5_gsea/visualize/region_wise_KEGG_positive_p005_hits.tsv",
        "scope": "cross",
        "output": ROOT / "manuscript/supplementary_table2-cross_celltype_KEGG.tsv",
        "mirror": ROOT / "manuscript/supp/Supplementary_Table_cross_celltype_KEGG_GSEA.tsv",
    },
]


def filter_table(path: Path, scope: str) -> pd.DataFrame:
    data = pd.read_csv(path, sep="\t")
    required = {"module", "celltype_region"}
    if not required.issubset(data.columns):
        raise ValueError(f"{path} lacks required columns: {sorted(required - set(data.columns))}")
    data["module"] = data["module"].astype(str).str.strip()
    data["celltype_region"] = data["celltype_region"].astype(str).str.strip()
    data = data.loc[~data["celltype_region"].str.startswith("Vasc_")].copy()

    if scope == "celltype":
        candidate_pairs = set(CELLTYPE_CANDIDATES)
        data["candidate_module"] = data["celltype_region"] + "_" + data["module"]
        selected = data["candidate_module"].isin(candidate_pairs)
    else:
        data["_region"] = data["celltype_region"].str.rsplit("_", n=1).str[-1]
        data["candidate_module"] = data["_region"] + "_" + data["module"]
        selected = data["candidate_module"].isin(CROSS_CANDIDATES)
    result = data.loc[selected].copy()
    if "_region" in result:
        result = result.drop(columns="_region")
    columns = ["candidate_module"] + [c for c in result.columns if c != "candidate_module"]
    return result[columns].sort_values(["candidate_module", "celltype_region", "term"], kind="stable").reset_index(drop=True)


if len(CELLTYPE_CANDIDATES) != 23 or len(CROSS_CANDIDATES) != 6:
    raise ValueError("Unexpected current selected-module inventory")

for config in SOURCES:
    result = filter_table(config["path"], config["scope"])
    config["output"].parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(config["output"], sep="\t", index=False, na_rep="NA")
    result.to_csv(config["mirror"], sep="\t", index=False, na_rep="NA")
    print(f"{config['name']}: {len(result)} rows across {result['candidate_module'].nunique()} candidate modules -> {config['output']}")
