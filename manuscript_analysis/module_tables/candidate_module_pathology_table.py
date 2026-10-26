#!/usr/bin/env python3
"""Prepare candidate-module linear and nonAD-vs-lateAD pathology tables."""
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
        "scope": "celltype",
        "path": SEPARATE / "02_run_rega/pathology_association/selected_module_significant_test_details.tsv",
        "output": ROOT / "manuscript/supplementary_table1_celltype_candidate_module_pathology_tests.tsv",
        "mirror": ROOT / "manuscript/supp/Supplementary_Table_celltype_candidate_module_pathology_tests.tsv",
    },
    {
        "scope": "cross",
        "path": SHARED / "part2_run_REGA/08pathology_tests/ALL_regions_sharedH_All_four_tests_H_pathology_results.tsv",
        "output": ROOT / "manuscript/supplementary_table1_cross_celltype_candidate_module_pathology_tests.tsv",
        "mirror": ROOT / "manuscript/supp/Supplementary_Table_cross_celltype_candidate_module_pathology_tests.tsv",
    },
]


def select(path: Path, scope: str) -> pd.DataFrame:
    data = pd.read_csv(path, sep="\t")
    required = {"test", "Module", "Comparison"}
    if scope == "celltype":
        required |= {"celltype_region", "celltype_region_module"}
    else:
        required.add("Region")
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f"{path} lacks required columns: {sorted(missing)}")

    linear = data["test"].astype(str).eq("Linear_Regression")
    pairwise = (
        data["test"].astype(str).eq("Pairwise_MannWhitneyU")
        & data["Comparison"].astype(str).str.strip().str.casefold().eq("nonad vs latead")
    )
    data = data.loc[linear | pairwise].copy()

    if scope == "celltype":
        data["candidate_module"] = data["celltype_region_module"].astype(str)
        data = data[data["candidate_module"].isin(CELLTYPE_CANDIDATES)]
    else:
        data["candidate_module"] = data["Region"].astype(str).str.strip() + "_" + data["Module"].astype(str).str.strip()
        data = data[data["candidate_module"].isin(CROSS_CANDIDATES)]

    columns = ["candidate_module"] + [column for column in data.columns if column != "candidate_module"]
    return data[columns].sort_values(["candidate_module", "test", "Comparison"], kind="stable").reset_index(drop=True)


if len(CELLTYPE_CANDIDATES) != 23 or len(CROSS_CANDIDATES) != 6:
    raise ValueError("Unexpected current selected-module inventory")

for config in SOURCES:
    result = select(config["path"], config["scope"])
    expected = CELLTYPE_CANDIDATES if config["scope"] == "celltype" else CROSS_CANDIDATES
    if set(result["candidate_module"]) != expected:
        raise ValueError(f"Missing pathology tests for {expected - set(result['candidate_module'])}")
    config["output"].parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(config["output"], sep="\t", index=False, na_rep="NA")
    result.to_csv(config["mirror"], sep="\t", index=False, na_rep="NA")
    print(f"{config['output']}: {len(result)} rows, {result['candidate_module'].nunique()} candidate modules")
