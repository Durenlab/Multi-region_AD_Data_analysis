from pathlib import Path
import os
SHARED_ROOT = Path(os.environ["REGA_SHARED_ROOT"]) if "REGA_SHARED_ROOT" in os.environ else next(
    p for p in Path(__file__).resolve().parents if p.name in {"REGA_shared", "shared"}
)
SEPARATE_ROOT = Path(os.environ.get("REGA_SEPARATE_ROOT", "/project/zduren/durenlab/naqing_palmetto/REGA_separate"))
import os
import pandas as pd
from collections import defaultdict

# ============================================================
# Paths
# ============================================================

INPUT_DIR = str(SEPARATE_ROOT / '01_prepare/data/original_individuals')
OUT_DIR = str(SHARED_ROOT / "part1_prepare_REGA/00region_consistent_individuals")

os.makedirs(OUT_DIR, exist_ok=True)

# Skip summary table
matrix_files = [
    f for f in os.listdir(INPUT_DIR)
    if f.endswith(".tsv") and f != "summary_table.tsv"
]
if len(matrix_files) != 36 or any(f.startswith("Vasc_") for f in matrix_files):
    raise ValueError("Expected exactly 36 RNA matrices across six cell types and regions")

# ============================================================
# Helper: parse celltype and region from filename
# Example:
#   Ast_EC.tsv
#   Mic_Immune_MTC.tsv
# ============================================================

def parse_celltype_region(filename):
    name = filename.replace(".tsv", "")

    # Assumes region is the last underscore-separated field
    # and celltype may itself contain underscores, e.g. Mic_Immune
    parts = name.split("_")
    region = parts[-1]
    celltype = "_".join(parts[:-1])

    return celltype, region

# ============================================================
# 1. Read matrices and group by region
# ============================================================

region_to_matrices = defaultdict(dict)

for f in matrix_files:
    celltype, region = parse_celltype_region(f)

    path = os.path.join(INPUT_DIR, f)

    df = pd.read_csv(path, sep="\t", index_col=0)

    # Make sure individual IDs are strings
    df.columns = df.columns.astype(str)

    key = f"{celltype}_{region}"

    region_to_matrices[region][key] = df

print("Regions found:")
for region, matrices in region_to_matrices.items():
    print(region, len(matrices), "celltype matrices")

# ============================================================
# 2. For each region, find shared individuals and trim matrices
# ============================================================

summary_rows = []

for region, matrices in region_to_matrices.items():

    print(f"\nProcessing region: {region}")

    individual_sets = [
        set(df.columns)
        for df in matrices.values()
    ]

    shared_individuals = sorted(list(set.intersection(*individual_sets)))

    print(f"Shared individuals in {region}: {len(shared_individuals)}")

    # Save shared individual list
    pd.DataFrame({"Individual": shared_individuals}).to_csv(
        os.path.join(OUT_DIR, f"{region}_shared_individuals.tsv"),
        sep="\t",
        index=False
    )

    # Trim and save each matrix
    for key, df in matrices.items():

        trimmed_df = df.loc[:, shared_individuals]

        out_file = os.path.join(
            OUT_DIR,
            f"{key}_region_consistent.tsv"
        )

        trimmed_df.to_csv(out_file, sep="\t")

        celltype = "_".join(key.split("_")[:-1])

        summary_rows.append({
            "Celltype": celltype,
            "Region": region,
            "Original_Individual_Count": df.shape[1],
            "Region_Shared_Individual_Count": len(shared_individuals),
            "Gene_Count": df.shape[0],
            "Output_File": out_file
        })

# ============================================================
# 3. Save summary
# ============================================================

summary_df = pd.DataFrame(summary_rows)

summary_df.to_csv(
    os.path.join(OUT_DIR, "region_wise_summary_table.tsv"),
    sep="\t",
    index=False
)

print("\nDone.")
print(f"Region-wise consistent matrices saved to: {OUT_DIR}")
