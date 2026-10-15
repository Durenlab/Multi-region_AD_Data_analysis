#!/usr/bin/env python3
"""Normalize each retained cell type–region RNA pseudobulk matrix."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pseudobulk", type=Path, required=True)
    parser.add_argument("--protein-coding-genes", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    pb = sc.read_h5ad(args.pseudobulk)
    if "mean" not in pb.layers:
        raise RuntimeError("Pseudobulk mean layer is missing")
    if (pb.obs["Celltype"] == "Vasc").any():
        raise RuntimeError("Vasc is excluded from this analysis")
    pb.X = pb.layers["mean"].copy().astype(float)
    protein_coding = set(
        pd.read_csv(args.protein_coding_genes, sep="\t", header=None)[0].astype(str)
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary = []

    for (celltype, region), pair_obs in pb.obs.groupby(
        ["Celltype", "BrainRegion"], observed=True
    ):
        label = f"{celltype}_{region}"
        subset = pb[pair_obs.index].copy()
        if subset.obs["Individual"].duplicated().any():
            raise RuntimeError(f"Duplicate donor in {label}")
        sc.pp.normalize_total(subset, target_sum=1e6)
        values = subset.X.toarray() if hasattr(subset.X, "toarray") else np.asarray(subset.X)
        fraction = np.asarray((values > 1).sum(axis=0)).ravel() / subset.n_obs
        subset = subset[:, fraction >= 0.20].copy()
        sc.pp.log1p(subset)
        subset = subset[:, subset.var_names.isin(protein_coding)].copy()
        values = subset.X.toarray() if hasattr(subset.X, "toarray") else np.asarray(subset.X)
        matrix = pd.DataFrame(
            values.T, index=subset.var_names, columns=subset.obs["Individual"].astype(str)
        )
        matrix.to_csv(args.output_dir / f"{label}.tsv", sep="\t")
        summary.append({
            "Celltype": celltype,
            "Region": region,
            "Individual_Count": subset.n_obs,
            "Gene_Count": subset.n_vars,
        })
        print(f"{label}: {subset.n_obs} donors, {subset.n_vars} genes")

    if len(summary) != 36:
        raise RuntimeError(f"Expected 36 nonvascular pairs, found {len(summary)}")
    pd.DataFrame(summary).to_csv(
        args.output_dir / "summary_table.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
