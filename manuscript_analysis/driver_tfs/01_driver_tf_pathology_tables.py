#!/usr/bin/env python3
"""Build compact supplementary tables of significant pathology-associated TFs."""
from pathlib import Path
import os

import pandas as pd


ROOT = Path(os.environ.get("AD_MANUSCRIPT_WORK_ROOT", Path(__file__).resolve().parents[3])).expanduser().resolve()
BASE = Path(os.environ.get("REGA_DATA_ROOT", ROOT.parent)).expanduser().resolve()
SEPARATE = BASE / "REGA_separate"
SHARED = BASE / "REGA_shared"
SOURCES = {
    "celltype_specific": {
        "association": SEPARATE / "07_tf/results/activity_pathology/all_celltype_region_TF_pathology_association.tsv",
        "assignments": SEPARATE / "07_tf/results/tf_assignment",
        "output": ROOT / "manuscript/supplementary_table5-Cell_type_specific_Driver_TF.tsv",
    },
    "cross_celltype": {
        "association": SHARED / "part7_TF/04DriverTF_lm/all_celltype_region_TF_pathology_association.tsv",
        "assignments": SHARED / "part7_TF/02tf_assignment_FC",
        "output": ROOT / "manuscript/supplementary_table5-Cross_celltype_Driver_TF.tsv",
    },
}


def read_assignments(scope: str, root: Path) -> tuple[dict[tuple[str, str], str], set[str]]:
    """Map the correct celltype(-region) and TF universe to its module."""
    mapping: dict[tuple[str, str], str] = {}
    available_contexts: set[str] = set()
    paths = sorted(root.glob("*_TF_module_assignments.tsv")) if scope == "celltype_specific" else sorted(root.glob("*/*_ALL_cell_TF_FC_assignment_with_unassigned.tsv"))
    if scope == "celltype_specific":
        key_col, tf_col, module_col = "celltype_region", "entity", "assigned_module"
    else:
        key_col, tf_col, module_col = "celltype_region", "tf", "assigned_module"

    for path in paths:
        frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        required = {tf_col, module_col} | ({key_col} if scope == "cross_celltype" else set())
        if not required.issubset(frame.columns):
            raise ValueError(f"{path} lacks required assignment columns: {sorted(required - set(frame.columns))}")
        for _, row in frame.iterrows():
            tf = row[tf_col].strip().upper()
            module = row[module_col].strip()
            if scope == "celltype_specific":
                # The per-pair celltype-specific assignment file is authoritative.
                pair = path.name.removesuffix("_TF_module_assignments.tsv")
            else:
                pair = row[key_col].strip()
            if pair:
                available_contexts.add(pair)
            if not pair or not tf or not module:
                continue
            key = (pair, tf)
            if key in mapping and mapping[key] != module:
                raise ValueError(f"Multiple modules found for {scope} assignment {key}: {mapping[key]}, {module}")
            mapping[key] = module
    return mapping, available_contexts


def build(scope: str, config: dict[str, Path]) -> pd.DataFrame:
    data = pd.read_csv(config["association"], sep="\t")
    required = {
        "celltype_region", "celltype", "region", "TF", "n_samples",
        "lm_beta_pathology", "lm_pval_pathology", "lm_fdr_pathology",
        "mean_nonAD", "mean_earlyAD", "mean_lateAD", "lateAD_minus_nonAD",
        "mannwhitney_lateAD_vs_nonAD_pval", "mannwhitney_lateAD_vs_nonAD_fdr",
        "n_lateAD", "n_nonAD",
    }
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f"{config['association']} lacks required columns: {sorted(missing)}")

    data = data.loc[~data["celltype_region"].astype(str).str.startswith("Vasc_")].copy()
    data["TF"] = data["TF"].astype(str).str.strip().str.upper()
    lm = (data["lm_beta_pathology"] > 0) & (data["lm_fdr_pathology"] < 0.05)
    pairwise = (data["lateAD_minus_nonAD"] > 0) & (data["mannwhitney_lateAD_vs_nonAD_fdr"] < 0.05)
    data = data.loc[lm | pairwise].copy()
    data["significant_test"] = [
        "Linear model and pairwise" if a and b else "Linear model" if a else "Pairwise"
        for a, b in zip(lm.loc[data.index], pairwise.loc[data.index])
    ]

    assignments, available_contexts = read_assignments(scope, config["assignments"])
    data["TF_assignment"] = [
        assignments.get(
            (str(pair), tf),
            "Unassigned" if str(pair) in available_contexts else "Not assessed (no candidate module)",
        )
        for pair, tf in zip(data["celltype_region"], data["TF"])
    ]
    columns = [
        "celltype_region", "celltype", "region", "TF", "TF_assignment", "significant_test",
        "n_samples", "n_nonAD", "n_lateAD", "mean_nonAD", "mean_earlyAD", "mean_lateAD",
        "lm_beta_pathology", "lm_pval_pathology", "lm_fdr_pathology",
        "lateAD_minus_nonAD", "mannwhitney_lateAD_vs_nonAD_pval",
        "mannwhitney_lateAD_vs_nonAD_fdr",
    ]
    return data[columns].sort_values(["region", "celltype", "TF"], kind="stable").reset_index(drop=True)


for scope, config in SOURCES.items():
    table = build(scope, config)
    config["output"].parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(config["output"], sep="\t", index=False, na_rep="NA")
    print(f"{config['output']}: {len(table)} rows")
