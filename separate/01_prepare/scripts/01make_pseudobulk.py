#!/usr/bin/env python3
"""Create mean RNA pseudobulk for nonvascular donor groups with >=10 nuclei."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True, help="Source snRNA AnnData")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    rna = sc.read_h5ad(args.source)
    pfc_study = "Mathys, et al., Cell, 2023 and/or Xiong et al., Cell, 2023."
    other_study = "Mathys et al., Nature, 2024."
    study_mask = (
        ((rna.obs["BrainRegion"] == "PFC") & (rna.obs["Study"] == pfc_study))
        | ((rna.obs["BrainRegion"] != "PFC") & (rna.obs["Study"] == other_study))
    )
    rna = rna[study_mask & (rna.obs["RNA.Class"] != "Vasc")].copy()
    rna.obs["pb_group"] = (
        rna.obs["RNA.Class"].astype(str)
        + "__" + rna.obs["ROSMAP_IndividualID"].astype(str)
        + "__" + rna.obs["BrainRegion"].astype(str)
    )
    group_counts = rna.obs.groupby("pb_group", observed=True).size().rename("n_nuclei")
    retained = group_counts[group_counts >= 10]
    rna = rna[rna.obs["pb_group"].isin(retained.index)].copy()

    pseudobulk = sc.get.aggregate(rna, by="pb_group", func="mean")
    metadata = (
        rna.obs.drop_duplicates("pb_group").set_index("pb_group")
        [["RNA.Class", "ROSMAP_IndividualID", "BrainRegion", "Pathology", "pmi", "age_death", "Sex"]]
        .rename(columns={"RNA.Class": "Celltype", "ROSMAP_IndividualID": "Individual"})
    )
    if metadata.groupby(["Celltype", "BrainRegion"], observed=True).ngroups != 36:
        raise RuntimeError("Expected 36 nonvascular cell type–region groups")
    pseudobulk.obs = metadata.loc[pseudobulk.obs_names].copy()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pseudobulk.write_h5ad(
        args.output_dir / "rna_pseudobulk_celltype_individual_region_mean.h5ad"
    )
    retained.rename_axis("pb_group").to_csv(
        args.output_dir / "retained_group_nuclei.tsv", sep="\t"
    )
    print(f"Retained {pseudobulk.n_obs} donor groups across six nonvascular cell types")


if __name__ == "__main__":
    main()
