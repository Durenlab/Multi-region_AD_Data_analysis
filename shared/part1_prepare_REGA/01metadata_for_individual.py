from pathlib import Path
import os
SHARED_ROOT = Path(os.environ["REGA_SHARED_ROOT"]) if "REGA_SHARED_ROOT" in os.environ else next(
    p for p in Path(__file__).resolve().parents if p.name in {"REGA_shared", "shared"}
)
SEPARATE_ROOT = Path(os.environ["REGA_SEPARATE_ROOT"])
import pandas as pd
import scanpy as ad
import os

# ============================================================
# Paths
# ============================================================

MASTER_H5AD = str(SEPARATE_ROOT / '01_prepare/data/source_snRNA_matrix.h5ad')

input_dir = str(SHARED_ROOT / "part1_prepare_REGA/00region_consistent_individuals")
output_dir = str(SHARED_ROOT / "part1_prepare_REGA/01region_consistent_metadata")

os.makedirs(output_dir, exist_ok=True)

# ============================================================
# Load master AnnData
# ============================================================

print("Loading master AnnData...")

adata = ad.read_h5ad(MASTER_H5AD, backed="r")

# ============================================================
# Study name cleanup
# ============================================================

study_map = {
    "Mathys et al., Nature, 2024.": "Mathys_2024",
    "Mathys, et al., Cell, 2023 and/or Xiong et al., Cell, 2023.": "Combined_2023",
    "Xiong et al., Cell, 2023.": "Xiong_2023",
    "This study": "Internal"
}

# ============================================================
# Extract useful metadata
# ============================================================

obs_subset = adata.obs[
    [
        "ROSMAP_IndividualID",
        "RNA.Class",
        "BrainRegion",
        "Study",
        "age_death",
        "Sex",
        "pmi",
        "Pathology"
    ]
].copy()
adata.file.close()

obs_subset["age_death"] = (
    obs_subset["age_death"]
    .replace("90+", 90)
    .astype(float)
)

obs_subset["Study"] = (
    obs_subset["Study"]
    .astype(str)
    .map(study_map)
    .fillna("Other")
)

obs_subset["Sex"] = obs_subset["Sex"].astype(str)
obs_subset["Pathology"] = obs_subset["Pathology"].astype(str)
obs_subset["BrainRegion"] = obs_subset["BrainRegion"].astype(str)
obs_subset["RNA.Class"] = obs_subset["RNA.Class"].astype(str)
obs_subset["ROSMAP_IndividualID"] = obs_subset["ROSMAP_IndividualID"].astype(str)

# ============================================================
# Cell types and regions
# ============================================================

cell_types = ["Ast", "Exc", "Inh", "Mic_Immune", "Oli", "OPC"]
regions = ["AG", "EC", "HC", "MTC", "PFC", "TH"]

# ============================================================
# Helper function
# ============================================================

def make_metadata_for_matrix(ct, reg, rna_file, output_dir):
    """
    Make metadata file matching one region-consistent RNA matrix.
    Rows are individuals.
    Row order exactly matches RNA matrix columns.
    """

    # Read only header to get individual IDs from RNA matrix
    rna_cols = pd.read_csv(
        rna_file,
        sep="\t",
        index_col=0,
        nrows=0
    ).columns.astype(str).tolist()

    # Filter cells for this cell type and brain region
    current_obs = obs_subset[
        (obs_subset["RNA.Class"] == ct) &
        (obs_subset["BrainRegion"] == reg)
    ].copy()

    if current_obs.shape[0] == 0:
        print(f"WARNING: no obs metadata found for {ct}_{reg}")
        return None

    # Collapse cell-level metadata to individual-level metadata
    stratified_meta = (
        current_obs
        .groupby("ROSMAP_IndividualID")
        .first()
    )

    # Check whether all RNA individuals exist in metadata
    missing_ids = sorted(list(set(rna_cols) - set(stratified_meta.index)))

    if len(missing_ids) > 0:
        print(f"WARNING: {ct}_{reg} has {len(missing_ids)} RNA individuals missing from metadata")
        print(missing_ids[:10])
        return None

    # Reorder metadata to exactly match RNA matrix columns
    final_meta = stratified_meta.loc[rna_cols].copy()

    # Optional sanity check
    assert list(final_meta.index) == rna_cols

    # Save metadata
    out_name = os.path.join(
        output_dir,
        f"{ct}_{reg}_region_consistent_metadata.tsv"
    )

    final_meta.to_csv(out_name, sep="\t")

    return {
        "Celltype": ct,
        "Region": reg,
        "RNA_File": rna_file,
        "Metadata_File": out_name,
        "Individual_Count": len(rna_cols),
        "Metadata_Row_Count": final_meta.shape[0]
    }

# ============================================================
# Generate metadata for each region-consistent RNA matrix
# ============================================================

summary_rows = []

for ct in cell_types:
    for reg in regions:

        rna_file = os.path.join(
            input_dir,
            f"{ct}_{reg}_region_consistent.tsv"
        )

        if not os.path.exists(rna_file):
            print(f"Skipping missing RNA file: {rna_file}")
            continue

        print(f"Making metadata for: {ct}_{reg}")

        result = make_metadata_for_matrix(
            ct=ct,
            reg=reg,
            rna_file=rna_file,
            output_dir=output_dir
        )

        if result is not None:
            summary_rows.append(result)

# ============================================================
# Save summary table
# ============================================================

summary_df = pd.DataFrame(summary_rows)

summary_df.to_csv(
    os.path.join(output_dir, "region_consistent_metadata_summary.tsv"),
    sep="\t",
    index=False
)

print("\nDone.")
print(f"Metadata saved to: {output_dir}")
print(f"Metadata files generated: {len(summary_rows)}")
