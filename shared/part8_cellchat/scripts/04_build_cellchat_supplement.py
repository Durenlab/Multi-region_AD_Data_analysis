#!/usr/bin/env python3
"""Write the donor-level linear CellChat supplementary table."""

from pathlib import Path
import os

import pandas as pd


BASE = Path(os.environ["CELLCHAT_WORK_DIR"])
REGIONS = ("AG", "EC", "HC", "MTC", "PFC", "TH")
CELLTYPES = {"Exc", "Inh", "Ast", "Oli", "OPC", "Mic_Immune", "Vasc"}
RESULTS = BASE / "finalized" / "tables"


def main() -> None:
    if (RESULTS / "Supplementary_Table_6_CellChat_significant.tsv").exists():
        raise FileExistsError("Supplementary Table 6 already exists; use a new CELLCHAT_WORK_DIR")
    tests = pd.read_csv(RESULTS / "pathology_progression_tests.tsv", sep="\t")
    if not set(tests.source).issubset(CELLTYPES) or not set(tests.target).issubset(CELLTYPES):
        raise ValueError("Progression table includes an unknown cell class")
    metadata = []
    for region in REGIONS:
        path = BASE / "tables" / f"{region}_candidate_interactions.tsv"
        table = pd.read_csv(path, sep="\t")
        table.insert(0, "region", region)
        required = ["region", "interaction_id", "pathway_name", "ligand", "receptor"]
        metadata.append(table[required].drop_duplicates(["region", "interaction_id"]))
    metadata = pd.concat(metadata, ignore_index=True)
    merged = tests.merge(metadata, on=["region", "interaction_id"], how="left", validate="many_to_one")
    if merged.pathway_name.isna().any():
        raise ValueError("Missing candidate interaction metadata")
    columns = [
        "region", "source", "target", "pathway_name", "ligand", "receptor",
        "interaction_id", "n_nonAD", "n_earlyAD", "n_lateAD",
        "mean_nonAD", "mean_earlyAD", "mean_lateAD",
        "beta_pathology_stage", "se_HC3", "ci_low", "ci_high",
        "p_value", "q_value_within_family", "FDR_family", "status",
    ]
    merged = merged[columns].sort_values(
        ["region", "source", "target", "q_value_within_family", "interaction_id"],
        na_position="last",
    )
    RESULTS.mkdir(parents=True, exist_ok=True)
    merged.to_csv(RESULTS / "CellChat_all_linear_tests.tsv", sep="\t", index=False, na_rep="NA")
    significant = merged.loc[
        merged.status.eq("tested") & merged.q_value_within_family.lt(0.05)
    ].copy()
    significant.to_csv(
        RESULTS / "Supplementary_Table_6_CellChat_significant.tsv",
        sep="\t", index=False,
    )
    print(f"Saved {len(merged):,} candidate rows and {len(significant):,} significant rows")


if __name__ == "__main__":
    main()
