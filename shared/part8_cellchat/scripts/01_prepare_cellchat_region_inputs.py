#!/usr/bin/env python3
"""
Prepare region-specific, donor-balanced, normalized AnnData files for CellChat.

INPUT
-----
A single AnnData/H5AD object containing all nuclei/cells, regions, and cell types.

Expected expression matrix:
    adata.X = raw UMI counts

Required adata.obs columns:
    BrainRegion
    RNA.Class
    ROSMAP_IndividualID
    Sample

Optional adata.obs columns retained when present:
    RNA.Subclass, RNA.Subtype, Pathology, Sex, age_death, pmi

OUTPUT
------
For each region, this script writes:
    <REGION>_cellchat_input.h5ad
        X: normalized to 10,000 counts/cell and log1p transformed
        obs: CellChat metadata, including:
             celltype = RNA.Class
             donor_id = ROSMAP_IndividualID
             samples = ROSMAP_IndividualID
             sample_id = Sample
        var: original gene metadata after removing all-zero genes

    <REGION>_celltype_counts_before.tsv
    <REGION>_celltype_counts_after.tsv
    <REGION>_donor_celltype_counts_before.tsv
    <REGION>_donor_celltype_counts_after.tsv
    <REGION>_sample_donor_celltype_counts_after.tsv

It also writes combined tables across all regions and a manifest.

Sampling design
---------------
Cells are balanced by biological donor, not by library/sample:
    region x RNA.Class x ROSMAP_IndividualID

For each donor-cell-type combination, at most
MAX_CELLS_PER_DONOR_CELLTYPE cells are retained. If a cell type still exceeds
MAX_CELLS_PER_CELLTYPE in a region, a donor-balanced round-robin cap is used.

The large source H5AD is read in backed mode, and only selected cells from one
region are loaded into memory at a time.
"""

from __future__ import annotations

import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy import sparse


# =============================================================================
# STEP 0. USER CONFIGURATION
# =============================================================================

# INPUT: study-harmonized whole-RNA AnnData with raw counts in X.
INPUT_H5AD = Path(os.environ["CELLCHAT_SOURCE_H5AD"])

# OUTPUT: one normalized H5AD per brain region plus QC/count tables.
OUTPUT_DIR = Path(os.environ["CELLCHAT_REGION_INPUT_DIR"])

REGION_ORDER = ["AG", "EC", "HC", "MTC", "PFC", "TH"]

CELLTYPE_ORDER = [
    "Exc",
    "Inh",
    "Ast",
    "Oli",
    "OPC",
    "Mic_Immune",
    "Vasc",
]

# Required metadata columns in the source object.
REGION_COL = "BrainRegion"
CELLTYPE_COL = "RNA.Class"
DONOR_COL = "ROSMAP_IndividualID"
SAMPLE_COL = "Sample"
TOTAL_COUNTS_COL = "total_counts"

OPTIONAL_OBS_COLUMNS = [
    "RNA.Subclass",
    "RNA.Subtype",
    "Pathology",
    "Sex",
    "age_death",
    "pmi",
]

# Donor-balanced downsampling parameters.
MAX_CELLS_PER_DONOR_CELLTYPE = 100
MAX_CELLS_PER_CELLTYPE = 5000
MIN_CELLS_PER_CELLTYPE = 50

NORMALIZATION_TARGET_SUM = 10_000
RANDOM_SEED = 20260727

# CellChat only needs normalized/log1p X. Keep False unless raw counts are also
# needed in each output H5AD.
SAVE_RAW_COUNTS_LAYER = False


# =============================================================================
# STEP 1. VALIDATION AND HELPER FUNCTIONS
# =============================================================================


def ensure_required_metadata(adata: AnnData) -> None:
    """Stop early if required source metadata are missing."""

    required = [REGION_COL, CELLTYPE_COL, DONOR_COL, SAMPLE_COL]
    missing = [column for column in required if column not in adata.obs.columns]

    if missing:
        raise KeyError(
            "The input AnnData is missing required obs columns: "
            f"{missing}"
        )

    for column in required:
        values = adata.obs[column]
        missing_count = int(values.isna().sum())
        blank_count = int(values.astype(str).str.strip().eq("").sum())

        if missing_count > 0 or blank_count > 0:
            raise ValueError(
                f"obs[{column!r}] contains {missing_count} missing and "
                f"{blank_count} blank values."
            )


def check_celltype_labels(obs: pd.DataFrame) -> None:
    """Confirm that RNA.Class uses the seven expected source labels."""

    observed = set(obs[CELLTYPE_COL].dropna().astype(str).str.strip().unique())
    expected = set(CELLTYPE_ORDER)

    unexpected = sorted(observed - expected)
    missing = sorted(expected - observed)

    if unexpected:
        raise ValueError(
            f"Unexpected values in obs[{CELLTYPE_COL!r}]: {unexpected}. "
            "Update CELLTYPE_ORDER or correct the source metadata."
        )

    if missing:
        print(
            "Warning: the full dataset does not contain these expected "
            f"cell types: {missing}"
        )


def check_sample_to_donor_mapping(obs: pd.DataFrame) -> None:
    """Require each original Sample label to map to exactly one donor."""

    mapping = (
        obs[[SAMPLE_COL, DONOR_COL]]
        .drop_duplicates()
        .groupby(SAMPLE_COL, observed=True)[DONOR_COL]
        .nunique()
    )

    bad = mapping[mapping > 1]
    if not bad.empty:
        raise ValueError(
            "Some Sample labels map to more than one donor. Examples:\n"
            f"{bad.head(20)}"
        )

    print(
        "Sample-to-donor mapping check passed: "
        f"{mapping.shape[0]} Sample labels, each mapped to one donor."
    )


def check_gene_names(adata: AnnData) -> None:
    """Confirm that var_names are unique human gene symbols."""

    genes = pd.Index(adata.var_names.astype(str))

    if not genes.is_unique:
        duplicated = genes[genes.duplicated()].unique()[:20].tolist()
        raise ValueError(
            "adata.var_names are duplicated. Resolve or aggregate duplicated "
            f"gene symbols before CellChat. Examples: {duplicated}"
        )

    ensembl_fraction = float(
        genes.str.match(r"^ENSG[0-9]+(?:\.[0-9]+)?$").mean()
    )
    print(
        "Fraction of var_names resembling Ensembl IDs: "
        f"{ensembl_fraction:.4f}"
    )

    if ensembl_fraction > 0.50:
        raise ValueError(
            "Most adata.var_names look like Ensembl IDs. CellChatDB.human "
            "expects human gene symbols."
        )


def validate_raw_counts(adata: AnnData, rng: np.random.Generator) -> None:
    """Check that adata.X behaves like raw, nonnegative integer UMI counts."""

    n_check = min(500, adata.n_obs)
    positions = np.sort(rng.choice(adata.n_obs, size=n_check, replace=False))
    x_test = adata[positions, :].X

    if sparse.issparse(x_test):
        x_test = x_test.tocsr()
        values = x_test.data
    else:
        x_test = np.asarray(x_test)
        values = x_test.ravel()
        values = values[values != 0]

    if values.size == 0:
        raise ValueError("The sampled expression matrix contains no nonzero values.")

    minimum = float(values.min())
    maximum = float(values.max())
    integer_fraction = float(
        np.mean(np.isclose(values, np.round(values), rtol=0, atol=1e-8))
    )

    print("Raw-count validation:")
    print(f"  sampled cells: {n_check}")
    print(f"  sampled nonzero values: {values.size:,}")
    print(f"  minimum nonzero value: {minimum}")
    print(f"  maximum value: {maximum}")
    print(f"  integer-valued fraction: {integer_fraction:.6f}")

    if minimum < 0:
        raise ValueError("Negative values detected in adata.X.")

    if integer_fraction < 0.999999:
        raise ValueError(
            "adata.X is not fully integer-valued and does not look like raw "
            "UMI counts."
        )

    if TOTAL_COUNTS_COL in adata.obs.columns:
        x_sums = np.asarray(x_test.sum(axis=1)).ravel()
        obs_sums = (
            adata.obs.iloc[positions][TOTAL_COUNTS_COL]
            .astype(float)
            .to_numpy()
        )

        relative_difference = np.abs(x_sums - obs_sums) / np.maximum(obs_sums, 1)
        median_relative_difference = float(np.median(relative_difference))
        correlation = float(np.corrcoef(x_sums, obs_sums)[0, 1])

        print(
            "  median relative difference versus "
            f"obs[{TOTAL_COUNTS_COL!r}]: {median_relative_difference:.6g}"
        )
        print(
            f"  correlation versus obs[{TOTAL_COUNTS_COL!r}]: "
            f"{correlation:.6f}"
        )

        if median_relative_difference > 1e-8 or correlation < 0.999999:
            raise ValueError(
                "adata.X row sums do not match obs['total_counts']; verify "
                "which matrix contains raw counts."
            )
    else:
        print(
            f"  {TOTAL_COUNTS_COL!r} is absent; row-sum comparison skipped."
        )


def prepare_metadata(adata: AnnData) -> pd.DataFrame:
    """Create standardized metadata for donor-balanced sampling and CellChat."""

    obs = adata.obs.copy()

    # Original row position is required for backed-mode matrix slicing later.
    obs["source_position"] = np.arange(adata.n_obs, dtype=np.int64)

    obs["region"] = obs[REGION_COL].astype(str).str.strip()
    obs["celltype"] = obs[CELLTYPE_COL].astype(str).str.strip()

    # Biological individual used for balanced sampling and as the CellChat
    # biological replicate.
    obs["donor_id"] = obs[DONOR_COL].astype(str).str.strip()
    obs["samples"] = obs["donor_id"]

    # Original tissue/library sample identifier, retained for traceability only.
    obs["sample_id"] = obs[SAMPLE_COL].astype(str).str.strip()
    obs["sample"] = obs["sample_id"]

    obs["celltype_region"] = obs["celltype"] + "_" + obs["region"]

    obs["celltype"] = pd.Categorical(
        obs["celltype"],
        categories=CELLTYPE_ORDER,
        ordered=True,
    )

    return obs


def donor_balanced_cap(
    donor_chunks: dict[str, np.ndarray],
    max_total: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Cap pooled cells while keeping donor contributions as even as possible."""

    total = sum(len(chunk) for chunk in donor_chunks.values())
    if total <= max_total:
        return np.concatenate(list(donor_chunks.values()))

    shuffled = {
        donor: rng.permutation(chunk)
        for donor, chunk in donor_chunks.items()
    }
    pointers = {donor: 0 for donor in shuffled}
    donors = list(shuffled)
    selected: list[int] = []

    while len(selected) < max_total:
        active = False
        donor_order = donors.copy()
        rng.shuffle(donor_order)

        for donor in donor_order:
            pointer = pointers[donor]
            donor_cells = shuffled[donor]

            if pointer < len(donor_cells):
                selected.append(int(donor_cells[pointer]))
                pointers[donor] += 1
                active = True

                if len(selected) == max_total:
                    break

        if not active:
            break

    if len(selected) != max_total:
        raise RuntimeError(
            f"Donor-balanced cap selected {len(selected)} cells instead of "
            f"the requested {max_total}."
        )

    return np.asarray(selected, dtype=np.int64)


def balanced_sample_region(
    region_obs: pd.DataFrame,
    rng: np.random.Generator,
) -> np.ndarray:
    """Select donor-balanced cells for one region."""

    selected_by_celltype: list[np.ndarray] = []

    for celltype in CELLTYPE_ORDER:
        celltype_obs = region_obs.loc[region_obs["celltype"] == celltype]

        if celltype_obs.empty:
            print(f"  {celltype}: absent in this region")
            continue

        donor_chunks: dict[str, np.ndarray] = {}

        for donor_id, donor_obs in celltype_obs.groupby(
            "donor_id",
            observed=True,
        ):
            positions = donor_obs["source_position"].to_numpy(dtype=np.int64)
            n_take = min(len(positions), MAX_CELLS_PER_DONOR_CELLTYPE)

            donor_chunks[str(donor_id)] = rng.choice(
                positions,
                size=n_take,
                replace=False,
            )

        pooled = donor_balanced_cap(
            donor_chunks=donor_chunks,
            max_total=MAX_CELLS_PER_CELLTYPE,
            rng=rng,
        )

        if len(pooled) < MIN_CELLS_PER_CELLTYPE:
            print(
                f"  {celltype}: excluded because only {len(pooled)} cells "
                f"remain (< {MIN_CELLS_PER_CELLTYPE})"
            )
            continue

        print(
            f"  {celltype}: {len(celltype_obs):,} original cells from "
            f"{len(donor_chunks)} donors -> {len(pooled):,} selected"
        )
        selected_by_celltype.append(pooled)

    if not selected_by_celltype:
        return np.array([], dtype=np.int64)

    return np.sort(np.concatenate(selected_by_celltype))


def build_output_obs(selected_obs: pd.DataFrame) -> pd.DataFrame:
    """Build metadata saved in each region-specific CellChat H5AD."""

    output_obs = pd.DataFrame(index=selected_obs.index.copy())

    output_obs["region"] = selected_obs["region"].astype(str)
    output_obs["celltype"] = pd.Categorical(
        selected_obs["celltype"].astype(str),
        categories=CELLTYPE_ORDER,
        ordered=True,
    )
    output_obs["celltype_region"] = selected_obs["celltype_region"].astype(str)

    output_obs["donor_id"] = selected_obs["donor_id"].astype(str)

    # Exact column name recognized by CellChat v2.
    output_obs["samples"] = selected_obs["donor_id"].astype(str)

    output_obs["sample_id"] = selected_obs["sample_id"].astype(str)
    output_obs["sample"] = selected_obs["sample"].astype(str)

    # Original source labels retained for traceability.
    output_obs[REGION_COL] = selected_obs[REGION_COL].astype(str)
    output_obs[CELLTYPE_COL] = selected_obs[CELLTYPE_COL].astype(str)
    output_obs[DONOR_COL] = selected_obs[DONOR_COL].astype(str)
    output_obs[SAMPLE_COL] = selected_obs[SAMPLE_COL].astype(str)

    for column in OPTIONAL_OBS_COLUMNS:
        if column in selected_obs.columns:
            output_obs[column] = selected_obs[column]

    if not np.array_equal(
        output_obs["samples"].to_numpy(),
        output_obs["donor_id"].to_numpy(),
    ):
        raise RuntimeError("'samples' and 'donor_id' do not match.")

    return output_obs


def extract_expression_matrix(
    adata: AnnData,
    selected_positions: np.ndarray,
):
    """Load only selected rows of X from the backed source H5AD."""

    x_selected = adata[selected_positions, :].X

    if sparse.issparse(x_selected):
        return x_selected.tocsr().astype(np.float32)

    return sparse.csr_matrix(np.asarray(x_selected, dtype=np.float32))


def summarize_celltypes(obs: pd.DataFrame, count_column: str) -> pd.DataFrame:
    """Count cells by region and cell type."""

    return (
        obs.groupby(["region", "celltype"], observed=True)
        .size()
        .reset_index(name=count_column)
    )


def summarize_donors(obs: pd.DataFrame, count_column: str) -> pd.DataFrame:
    """Count cells by region, cell type, and biological donor."""

    return (
        obs.groupby(
            ["region", "celltype", "donor_id"],
            observed=True,
        )
        .size()
        .reset_index(name=count_column)
    )


def summarize_samples(obs: pd.DataFrame, count_column: str) -> pd.DataFrame:
    """Count cells by region, cell type, donor, and original Sample."""

    return (
        obs.groupby(
            ["region", "celltype", "donor_id", "sample_id"],
            observed=True,
        )
        .size()
        .reset_index(name=count_column)
    )


def assert_no_null_encoding(h5ad_path: Path) -> None:
    """Fail if the written H5AD contains a newer null encoding."""

    null_objects: list[str] = []

    with h5py.File(h5ad_path, "r") as handle:
        def inspect(name: str, obj: h5py.Group | h5py.Dataset) -> None:
            encoding = obj.attrs.get("encoding-type", "")
            if isinstance(encoding, bytes):
                encoding = encoding.decode("utf-8")
            if encoding == "null":
                null_objects.append("/" + name)

        handle.visititems(inspect)

    if null_objects:
        raise RuntimeError(
            f"Null-encoded objects found in {h5ad_path}: {null_objects}"
        )


# =============================================================================
# STEP 2. READ SOURCE OBJECT AND PROCESS EACH REGION
# =============================================================================


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    existing = sorted(OUTPUT_DIR.glob("*_cellchat_input.h5ad"))
    if existing:
        raise FileExistsError(
            f"Regional inputs already exist in {OUTPUT_DIR}; use a new output directory"
        )

    print("=" * 80)
    print("STEP 2: Read source H5AD in backed mode")
    print(f"INPUT : {INPUT_H5AD}")
    print(f"OUTPUT: {OUTPUT_DIR}")

    if not INPUT_H5AD.exists():
        raise FileNotFoundError(f"Input H5AD does not exist: {INPUT_H5AD}")

    adata = sc.read_h5ad(INPUT_H5AD, backed="r")
    print(f"Source object: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    try:
        ensure_required_metadata(adata)
        check_celltype_labels(adata.obs)
        check_sample_to_donor_mapping(adata.obs)
        check_gene_names(adata)

        rng = np.random.default_rng(RANDOM_SEED)
        validate_raw_counts(adata, rng)

        if not adata.obs_names.is_unique:
            raise ValueError("adata.obs_names are not unique.")

        print("\n" + "=" * 80)
        print("STEP 3: Standardize metadata and add CellChat labels")

        meta = prepare_metadata(adata)

        print("Cells by BrainRegion:")
        print(meta["region"].value_counts().sort_index())

        print("\nCells by RNA.Class:")
        print(meta["celltype"].value_counts(sort=False))

        before_all_celltype = summarize_celltypes(meta, "n_cells_before")
        before_all_donor = summarize_donors(meta, "n_cells_before")

        before_all_celltype.to_csv(
            OUTPUT_DIR / "ALL_regions_celltype_counts_before.tsv",
            sep="\t",
            index=False,
        )
        before_all_donor.to_csv(
            OUTPUT_DIR / "ALL_regions_donor_celltype_counts_before.tsv",
            sep="\t",
            index=False,
        )

        all_after_celltype: list[pd.DataFrame] = []
        all_after_donor: list[pd.DataFrame] = []
        all_after_sample: list[pd.DataFrame] = []
        manifest_rows: list[dict[str, object]] = []

        for region in REGION_ORDER:
            print("\n" + "=" * 80)
            print(f"REGION: {region}")

            print("STEP 4: Select cells from this brain region")
            region_meta = meta.loc[meta["region"] == region].copy()

            if region_meta.empty:
                print(f"Warning: no cells found for {region}; skipping.")
                continue

            before_region_celltype = summarize_celltypes(
                region_meta,
                "n_cells_before",
            )
            before_region_donor = summarize_donors(
                region_meta,
                "n_cells_before",
            )

            before_region_celltype.to_csv(
                OUTPUT_DIR / f"{region}_celltype_counts_before.tsv",
                sep="\t",
                index=False,
            )
            before_region_donor.to_csv(
                OUTPUT_DIR / f"{region}_donor_celltype_counts_before.tsv",
                sep="\t",
                index=False,
            )

            print("STEP 5: Donor-balanced downsampling by donor x RNA.Class")
            print(
                "  maximum per donor x cell type: "
                f"{MAX_CELLS_PER_DONOR_CELLTYPE}"
            )
            print(f"  maximum per cell type: {MAX_CELLS_PER_CELLTYPE}")

            selected_positions = balanced_sample_region(region_meta, rng)

            if selected_positions.size == 0:
                print(f"Warning: no cells retained for {region}; skipping.")
                continue

            selected_meta = meta.iloc[selected_positions].copy()
            output_obs = build_output_obs(selected_meta)

            print("STEP 6: Load selected raw counts into memory")
            print(f"  selected cells: {len(selected_positions):,}")

            x_selected = extract_expression_matrix(adata, selected_positions)
            var_selected = adata.var.copy()

            region_adata = AnnData(
                X=x_selected,
                obs=output_obs,
                var=var_selected,
            )

            if SAVE_RAW_COUNTS_LAYER:
                region_adata.layers["counts"] = region_adata.X.copy()

            genes_before = region_adata.n_vars
            sc.pp.filter_genes(region_adata, min_cells=1)
            genes_after = region_adata.n_vars

            print(
                "  genes retained after removing all-zero genes: "
                f"{genes_after:,}/{genes_before:,}"
            )

            print("STEP 7: Normalize to 10,000 counts/cell and apply log1p")
            sc.pp.normalize_total(
                region_adata,
                target_sum=NORMALIZATION_TARGET_SUM,
            )
            sc.pp.log1p(region_adata)

            # Older zellkonverter/anndata combinations cannot read the newer
            # null encoding that may be written for uns['log1p']['base'].
            # CellChat does not need this metadata, so remove it before writing.
            region_adata.uns.pop("log1p", None)

            if sparse.issparse(region_adata.X):
                normalized_values = region_adata.X.data
            else:
                normalized_values = np.asarray(region_adata.X).ravel()

            if normalized_values.size > 0 and normalized_values.min() < 0:
                raise ValueError(
                    "Negative values detected after normalization for region "
                    f"{region}."
                )

            region_adata.uns["cellchat_preparation"] = {
                "source_h5ad": str(INPUT_H5AD),
                "source_matrix": "raw UMI counts from X",
                "region": region,
                "celltype_column": CELLTYPE_COL,
                "region_column": REGION_COL,
                "donor_id_definition": DONOR_COL,
                "cellchat_samples_definition": DONOR_COL,
                "sample_id_definition": SAMPLE_COL,
                "max_cells_per_donor_celltype": int(
                    MAX_CELLS_PER_DONOR_CELLTYPE
                ),
                "max_cells_per_celltype": int(MAX_CELLS_PER_CELLTYPE),
                "min_cells_per_celltype": int(MIN_CELLS_PER_CELLTYPE),
                "normalization_target_sum": int(NORMALIZATION_TARGET_SUM),
                "log1p_applied": True,
                "random_seed": int(RANDOM_SEED),
            }

            print("STEP 8: Save normalized CellChat input and QC tables")
            output_h5ad = OUTPUT_DIR / f"{region}_cellchat_input.h5ad"
            region_adata.write_h5ad(output_h5ad, compression="gzip")
            assert_no_null_encoding(output_h5ad)

            after_celltype = summarize_celltypes(
                region_adata.obs,
                "n_cells_after",
            )
            after_donor = summarize_donors(
                region_adata.obs,
                "n_cells_after",
            )
            after_sample = summarize_samples(
                region_adata.obs,
                "n_cells_after",
            )

            after_celltype.to_csv(
                OUTPUT_DIR / f"{region}_celltype_counts_after.tsv",
                sep="\t",
                index=False,
            )
            after_donor.to_csv(
                OUTPUT_DIR / f"{region}_donor_celltype_counts_after.tsv",
                sep="\t",
                index=False,
            )
            after_sample.to_csv(
                OUTPUT_DIR / f"{region}_sample_donor_celltype_counts_after.tsv",
                sep="\t",
                index=False,
            )

            all_after_celltype.append(after_celltype)
            all_after_donor.append(after_donor)
            all_after_sample.append(after_sample)

            manifest_rows.append(
                {
                    "region": region,
                    "output_h5ad": str(output_h5ad),
                    "n_cells": int(region_adata.n_obs),
                    "n_genes": int(region_adata.n_vars),
                    "n_celltypes": int(region_adata.obs["celltype"].nunique()),
                    "n_donors": int(region_adata.obs["donor_id"].nunique()),
                    "normalization": (
                        f"normalize_total(target_sum={NORMALIZATION_TARGET_SUM}) "
                        "+ log1p"
                    ),
                    "sampling_unit": "ROSMAP_IndividualID x RNA.Class",
                }
            )

            print(f"  OUTPUT H5AD: {output_h5ad}")
            print("  Final cells by cell type:")
            print(region_adata.obs["celltype"].value_counts(sort=False))

            del region_adata, x_selected, selected_meta, output_obs

        print("\n" + "=" * 80)
        print("STEP 9: Write combined summary tables")

        if all_after_celltype:
            pd.concat(all_after_celltype, ignore_index=True).to_csv(
                OUTPUT_DIR / "ALL_regions_celltype_counts_after.tsv",
                sep="\t",
                index=False,
            )

        if all_after_donor:
            pd.concat(all_after_donor, ignore_index=True).to_csv(
                OUTPUT_DIR / "ALL_regions_donor_celltype_counts_after.tsv",
                sep="\t",
                index=False,
            )

        if all_after_sample:
            pd.concat(all_after_sample, ignore_index=True).to_csv(
                OUTPUT_DIR / "ALL_regions_sample_donor_celltype_counts_after.tsv",
                sep="\t",
                index=False,
            )

        pd.DataFrame(manifest_rows).to_csv(
            OUTPUT_DIR / "cellchat_input_manifest.tsv",
            sep="\t",
            index=False,
        )

        print("\nFinished preparing region-specific CellChat inputs.")
        print(f"All outputs are under: {OUTPUT_DIR}")

    finally:
        if getattr(adata, "file", None) is not None:
            adata.file.close()


if __name__ == "__main__":
    main()
