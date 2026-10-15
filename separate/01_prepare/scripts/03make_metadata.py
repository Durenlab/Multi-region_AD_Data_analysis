#!/usr/bin/env python3
"""Write donor metadata aligned to each normalized RNA matrix."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True, help="Source snRNA AnnData")
    parser.add_argument("--matrix-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.source, backed="r")
    columns = [
        "ROSMAP_IndividualID", "RNA.Class", "BrainRegion", "Study",
        "age_death", "Sex", "pmi", "Pathology",
    ]
    obs = adata.obs[columns].copy()
    adata.file.close()
    pfc_study = "Mathys, et al., Cell, 2023 and/or Xiong et al., Cell, 2023."
    other_study = "Mathys et al., Nature, 2024."
    study_mask = (
        ((obs["BrainRegion"] == "PFC") & (obs["Study"] == pfc_study))
        | ((obs["BrainRegion"] != "PFC") & (obs["Study"] == other_study))
    )
    obs = obs.loc[study_mask & (obs["RNA.Class"] != "Vasc")].copy()
    obs["age_death"] = obs["age_death"].replace("90+", 90).astype(float)
    obs["Study"] = obs["Study"].astype(str).replace({
        other_study: "Mathys_2024", pfc_study: "Combined_2023",
    })
    obs["Sex"] = obs["Sex"].astype(str)
    obs["Pathology"] = obs["Pathology"].astype(str)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    matrices = sorted(p for p in args.matrix_dir.glob("*.tsv") if p.name != "summary_table.tsv")
    if len(matrices) != 36:
        raise RuntimeError(f"Expected 36 normalized matrices, found {len(matrices)}")
    for matrix_path in matrices:
        columns = pd.read_csv(matrix_path, sep="\t", index_col=0, nrows=0).columns
        label = matrix_path.stem
        celltype, region = label.rsplit("_", 1)
        current = obs.loc[
            (obs["RNA.Class"] == celltype) & (obs["BrainRegion"] == region)
        ]
        donor_metadata = current.groupby("ROSMAP_IndividualID", observed=True).first()
        missing = set(columns) - set(donor_metadata.index)
        if missing:
            raise RuntimeError(f"Missing metadata for {label}: {sorted(missing)}")
        selected = donor_metadata.loc[columns].copy()
        selected.index.name = "ROSMAP_IndividualID"
        selected.to_csv(args.output_dir / f"{label}_metadata.tsv", sep="\t")
        print(f"{label}: {len(selected)} donor metadata rows")


if __name__ == "__main__":
    main()
