#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# ============================================================
# Paths
# ============================================================
magma_file = "/project/zduren/durenlab/naqing_palmetto/REGA_separate/06_magma/results/merged/AD_Bellenguez_hg19_merged.genes.out"

motif_tf_file = "/project/zduren/durenlab/naqing_palmetto/REGA_separate/07_tf/data/motif_by_TF_matrix.tsv"

driver_file = "/project/zduren/durenlab/naqing_palmetto/REGA_separate/07_tf/results/activity_pathology/all_celltype_region_TF_pathology_association.tsv"

outdir = "/project/zduren/durenlab/naqing_palmetto/REGA_separate/07_tf/results/tf_snp_overlap"

os.makedirs(outdir, exist_ok=True)

gene_with_snp_file = os.path.join(outdir, "gene_with_snp.tsv")
tf_with_snp_file = os.path.join(outdir, "TF_with_SNP.tsv")
summary_file = os.path.join(outdir, "driverTF_SNP_overlap_summary.tsv")
detail_file = os.path.join(outdir, "significant_driverTF_with_SNP_annotation.tsv")
top20_file = os.path.join(outdir, "top20_driverTF_with_SNP_annotation.tsv")

fig_png = os.path.join(outdir, "top20_driverTF_SNP_overview.png")
fig_pdf = os.path.join(outdir, "top20_driverTF_SNP_overview.pdf")


# ============================================================
# 1. MAGMA genes: P < 1e-6
# ============================================================
magma = pd.read_csv(magma_file, sep=r"\s+", engine="python")

required_magma_cols = {"GENE", "P"}
missing = required_magma_cols - set(magma.columns)
if missing:
    raise ValueError(f"MAGMA file missing columns: {missing}")

magma["P"] = pd.to_numeric(magma["P"], errors="coerce")

gene_with_snp = (
    magma.loc[magma["P"] < 1e-6, ["GENE", "P"]]
    .dropna(subset=["GENE"])
    .drop_duplicates(subset=["GENE"])
    .sort_values("P")
    .reset_index(drop=True)
)

# User requested GENE name as gene_with_snp.tsv.
# Keep only one column to make it a clean gene list.
gene_with_snp[["GENE"]].to_csv(
    gene_with_snp_file,
    sep="\t",
    index=False
)

gene_snp_set = set(gene_with_snp["GENE"].astype(str))

print(f"MAGMA genes with P < 1e-6: {len(gene_snp_set):,}")
print(f"Saved: {gene_with_snp_file}")


# ============================================================
# 2. TF list from motif_by_TF_matrix.tsv
#    TFs = column names starting from the second column
# ============================================================
# Only need header, so nrows=0 avoids loading a potentially large matrix.
motif_header = pd.read_csv(
    motif_tf_file,
    sep="\t",
    nrows=0
)

if motif_header.shape[1] < 2:
    raise ValueError(
        "motif_by_TF_matrix.tsv has fewer than 2 columns; "
        "cannot extract TF names from columns 2 onward."
    )

tf_list = pd.Index(motif_header.columns[1:].astype(str))
tf_list = tf_list[tf_list.notna()]
tf_list = pd.Index(pd.unique(tf_list))

tf_with_snp = sorted(set(tf_list) & gene_snp_set)

pd.DataFrame({"TF": tf_with_snp}).to_csv(
    tf_with_snp_file,
    sep="\t",
    index=False
)

tf_snp_set = set(tf_with_snp)

print(f"TFs in motif matrix: {len(tf_list):,}")
print(f"TFs overlapping MAGMA P < 1e-6 genes: {len(tf_snp_set):,}")
print(f"Saved: {tf_with_snp_file}")


# ============================================================
# 3. Significant driver TFs
#    criterion: mannwhitney_lateAD_vs_nonAD_fdr < 0.05
# ============================================================
driver = pd.read_csv(driver_file, sep="\t")

required_driver_cols = {
    "celltype_region",
    "celltype",
    "region",
    "TF",
    "mannwhitney_lateAD_vs_nonAD_fdr",
}

missing = required_driver_cols - set(driver.columns)
if missing:
    raise ValueError(f"Driver-TF file missing columns: {missing}")

driver["mannwhitney_lateAD_vs_nonAD_fdr"] = pd.to_numeric(
    driver["mannwhitney_lateAD_vs_nonAD_fdr"],
    errors="coerce"
)

# Helpful ranking/effect columns if present
for col in [
    "mannwhitney_lateAD_vs_nonAD_pval",
    "lateAD_minus_nonAD",
    "spearman_rho_pathology",
    "spearman_fdr_pathology",
]:
    if col in driver.columns:
        driver[col] = pd.to_numeric(driver[col], errors="coerce")

sig = driver.loc[
    driver["mannwhitney_lateAD_vs_nonAD_fdr"] < 0.05
].copy()

sig["TF"] = sig["TF"].astype(str)
sig["has_SNP"] = sig["TF"].isin(tf_snp_set)

# In case the source file contains duplicate rows for the same TF in a group,
# keep the row with the smallest Mann-Whitney FDR.
sort_cols = [
    "celltype_region",
    "mannwhitney_lateAD_vs_nonAD_fdr",
]

if "mannwhitney_lateAD_vs_nonAD_pval" in sig.columns:
    sort_cols.append("mannwhitney_lateAD_vs_nonAD_pval")

sig = (
    sig.sort_values(sort_cols, na_position="last")
    .drop_duplicates(subset=["celltype_region", "TF"], keep="first")
    .reset_index(drop=True)
)

# Save all significant driver TFs with SNP annotation
sig.to_csv(detail_file, sep="\t", index=False)

print(f"Significant driver TF rows: {len(sig):,}")
print(f"Celltype-region groups with >=1 significant driver TF: "
      f"{sig['celltype_region'].nunique():,}")


# ============================================================
# 4. Summary table per celltype_region
# ============================================================
summary = (
    sig.groupby(
        ["celltype_region", "celltype", "region"],
        as_index=False
    )
    .agg(
        n_driver_TF=("TF", "nunique"),
        n_driver_TF_with_SNP=("has_SNP", "sum"),
    )
)

summary["fraction_driver_TF_with_SNP"] = (
    summary["n_driver_TF_with_SNP"] /
    summary["n_driver_TF"]
)

summary = summary.sort_values(
    ["region", "celltype", "celltype_region"]
).reset_index(drop=True)

summary.to_csv(summary_file, sep="\t", index=False)

print(f"Saved summary table: {summary_file}")


# ============================================================
# 5. Top 20 driver TFs per celltype_region
#
# Ranking:
#   1) smaller Mann-Whitney FDR
#   2) smaller raw Mann-Whitney p-value, if available
#   3) larger absolute lateAD-vs-nonAD difference, if available
# ============================================================
rank_cols = ["mannwhitney_lateAD_vs_nonAD_fdr"]
ascending = [True]

if "mannwhitney_lateAD_vs_nonAD_pval" in sig.columns:
    rank_cols.append("mannwhitney_lateAD_vs_nonAD_pval")
    ascending.append(True)

if "lateAD_minus_nonAD" in sig.columns:
    sig["abs_lateAD_minus_nonAD"] = sig["lateAD_minus_nonAD"].abs()
    rank_cols.append("abs_lateAD_minus_nonAD")
    ascending.append(False)

top20 = (
    sig.sort_values(
        ["celltype_region"] + rank_cols,
        ascending=[True] + ascending,
        na_position="last"
    )
    .groupby("celltype_region", group_keys=False)
    .head(20)
    .copy()
)

top20["rank"] = (
    top20.groupby("celltype_region")
    .cumcount() + 1
)

top20.to_csv(top20_file, sep="\t", index=False)

print(f"Saved top-20 table: {top20_file}")


# ============================================================
# 6. Plot
#
# One large overview:
#   y = celltype_region
#   x = rank among significant driver TFs (1..20)
#   text = TF name
#   yellow box = TF overlaps TF_with_SNP.tsv
#
# Only celltype_region groups with >=1 significant driver TF are included.
# ============================================================

if top20.empty:
    print("No significant driver TFs found at FDR < 0.05; no plot generated.")
else:
    # Order groups by region, then celltype.
    order_df = (
        top20[["celltype_region", "celltype", "region"]]
        .drop_duplicates()
        .sort_values(["region", "celltype", "celltype_region"])
        .reset_index(drop=True)
    )

    group_order = order_df["celltype_region"].tolist()

    # Give enough vertical room to each row.
    n_groups = len(group_order)
    fig_height = max(5.0, 0.48 * n_groups + 2.0)
    fig_width = 24

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    y_map = {
        group: n_groups - 1 - i
        for i, group in enumerate(group_order)
    }

    # light horizontal guide lines
    for group, y in y_map.items():
        ax.axhline(
            y=y,
            linewidth=0.5,
            alpha=0.15,
            zorder=0
        )

    # text labels
    for _, row in top20.iterrows():
        x = row["rank"]
        y = y_map[row["celltype_region"]]

        bbox = None
        fontweight = "normal"

        if bool(row["has_SNP"]):
            bbox = dict(
                boxstyle="round,pad=0.18",
                facecolor="#FFD54F",
                edgecolor="none",
                alpha=0.95
            )
            fontweight = "bold"

        ax.text(
            x,
            y,
            str(row["TF"]),
            ha="center",
            va="center",
            fontsize=7.5,
            rotation=55,
            rotation_mode="anchor",
            fontweight=fontweight,
            bbox=bbox,
            clip_on=False
        )

    ax.set_xlim(0.4, 20.6)
    ax.set_xticks(range(1, 21))
    ax.set_xlabel(
        "Rank among significant driver TFs "
        "(Mann–Whitney lateAD vs nonAD FDR < 0.05)",
        fontsize=11
    )

    ax.set_yticks([y_map[g] for g in group_order])
    ax.set_yticklabels(group_order, fontsize=9)
    ax.set_ylabel("Cell type–region", fontsize=11)

    ax.set_title(
        "Top 20 significant driver TFs by cell type–region\n"
        "Yellow highlight: TF overlaps MAGMA gene P < 1e-6",
        fontsize=13,
        pad=16
    )

    ax.grid(
        axis="x",
        linestyle=":",
        linewidth=0.5,
        alpha=0.25
    )

    # cleaner frame
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    legend_handles = [
        Patch(
            facecolor="#FFD54F",
            edgecolor="none",
            label="TF with MAGMA P < 1e-6"
        )
    ]
    ax.legend(
        handles=legend_handles,
        loc="upper right",
        frameon=False,
        fontsize=9
    )

    plt.tight_layout()

    fig.savefig(
        fig_png,
        dpi=300,
        bbox_inches="tight"
    )
    fig.savefig(
        fig_pdf,
        bbox_inches="tight"
    )
    plt.close(fig)

    print(f"Saved figure: {fig_png}")
    print(f"Saved figure: {fig_pdf}")


# ============================================================
# 7. Console summary
# ============================================================
print("\n=== Driver TF / SNP overlap summary ===")
print(
    summary[
        [
            "celltype_region",
            "n_driver_TF",
            "n_driver_TF_with_SNP",
            "fraction_driver_TF_with_SNP"
        ]
    ].to_string(index=False)
)
