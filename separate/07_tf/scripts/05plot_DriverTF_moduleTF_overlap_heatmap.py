#!/usr/bin/env python3

from pathlib import Path
import argparse
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch


# ============================================================
# Command-line options
# ============================================================
def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Identify positively AD-associated TFs in selected celltype-region contexts, "
            "overlap them with celltype-region-specific module assignments, and plot "
            "a top bar plot plus heatmap."
        )
    )
    parser.add_argument(
        "--significance",
        choices=["fdr", "nominal"],
        default="fdr",
        help=(
            "Significance measure for the Mann-Whitney lateAD vs nonAD test. "
            "Use 'fdr' for mannwhitney_lateAD_vs_nonAD_fdr or 'nominal' for "
            "mannwhitney_lateAD_vs_nonAD_pval. Default: fdr."
        ),
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=0.05,
        help="Significance cutoff. Default: 0.05.",
    )
    return parser.parse_args()


ARGS = parse_args()
SIGNIFICANCE_MODE = ARGS.significance
SIGNIFICANCE_CUTOFF = ARGS.cutoff


# ============================================================
# User settings
# ============================================================
BASE_DIR = Path(
    "/project/zduren/durenlab/naqing_palmetto/REGA_separate"
)

ASSOC_FILE = (
    BASE_DIR
    / "07_tf"
    / "results/activity_pathology"
    / "all_celltype_region_TF_pathology_association.tsv"
)

# This count table is used directly for the top bars only in FDR mode.
# In nominal mode, bar counts are derived from ASSOC_FILE because the supplied
# count table contains only FDR-based counts.
COUNT_FILE = (
    BASE_DIR
    / "07_tf"
    / "results/activity_pathology"
    / "mannwhitney_TF_count_plots"
    / "mannwhitney_significant_TF_counts.tsv"
)

ASSIGN_ROOT = (
    BASE_DIR
    / "07_tf"
    / "results/tf_assignment"
)

OUTDIR = (
    BASE_DIR
    / "07_tf"
    / "results/driver_tf_module_overlap"
    / f"mannwhitney_{SIGNIFICANCE_MODE}_{SIGNIFICANCE_CUTOFF:g}"
)

REGION_ORDER = ["AG", "EC", "HC", "MTC", "PFC", "TH"]
CELLTYPE_ORDER = ["Exc", "Inh", "Ast", "Oli", "OPC", "Mic_Immune"]
MODULE_ORDER = [f"M{i}" for i in range(1, 21)]
REGION_GAP = 0.65

# Only these celltype-region contexts were included in this REGA run.
SELECTED_CONTEXTS = [
    "Ast_EC",
    "Exc_MTC",
    "Exc_PFC",
    "Inh_AG",
    "Mic_Immune_AG",
    "Mic_Immune_EC",
    "Mic_Immune_MTC",
    "OPC_PFC",
    "Ast_AG",
    "Ast_HC",
    "Ast_MTC",
    "Exc_AG",
    "Exc_EC",
    "OPC_HC",
    "OPC_MTC",
    "Oli_EC",
    "Oli_MTC",
    "Oli_AG",
]
SELECTED_CONTEXT_SET = set(SELECTED_CONTEXTS)

# Plot selected contexts grouped by region, then by CELLTYPE_ORDER.
CONTEXTS = [
    f"{celltype}_{region}"
    for region in REGION_ORDER
    for celltype in CELLTYPE_ORDER
    if f"{celltype}_{region}" in SELECTED_CONTEXT_SET
]

if set(CONTEXTS) != SELECTED_CONTEXT_SET:
    missing = sorted(SELECTED_CONTEXT_SET - set(CONTEXTS))
    raise ValueError(
        "Some SELECTED_CONTEXTS do not match REGION_ORDER/CELLTYPE_ORDER: "
        f"{missing}"
    )

CELLTYPE_COLORS = {
    "Exc": "#3182bd",
    "Inh": "#ff7f0e",
    "Ast": "#31a354",
    "Oli": "#e41a1c",
    "OPC": "#9467bd",
    "Mic_Immune": "#9c755f",
}

if SIGNIFICANCE_MODE == "fdr":
    SIGNIFICANCE_COLUMN = "mannwhitney_lateAD_vs_nonAD_fdr"
    SIGNIFICANCE_LABEL = "FDR"
    COUNT_COLUMN = "n_positive_TFs_fdr_lt_0.05"
else:
    SIGNIFICANCE_COLUMN = "mannwhitney_lateAD_vs_nonAD_pval"
    SIGNIFICANCE_LABEL = "nominal p-value"
    COUNT_COLUMN = None


# ============================================================
# Helpers
# ============================================================
def module_number(value):
    match = re.search(r"(\d+)$", str(value))
    return int(match.group(1)) if match else 10**9


def require_columns(df, required, label):
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(f"{label} is missing columns: {missing}")


def get_celltype(celltype_region):
    """Extract cell type safely, including Mic_Immune."""
    for celltype in sorted(CELLTYPE_COLORS, key=len, reverse=True):
        if str(celltype_region).startswith(celltype + "_"):
            return celltype
    return str(celltype_region).rsplit("_", 1)[0]


def read_assignments():
    """Read one TF-assignment file for every celltype_region.

    Expected file example:
      ASSIGN_ROOT/Ast_AG_TF_module_assignments.tsv

    Expected columns:
      entity, assigned_module
    """
    frames = []

    for context in CONTEXTS:
        region = context.rsplit("_", 1)[-1]
        path = ASSIGN_ROOT / f"{context}_TF_module_assignments.tsv"

        if not path.exists():
            raise FileNotFoundError(f"Missing assignment file: {path}")

        table = pd.read_csv(path, sep="\t")
        require_columns(
            table,
            ["entity", "assigned_module"],
            str(path),
        )

        table = table[["entity", "assigned_module"]].copy()
        table = table.rename(
            columns={
                "entity": "tf",
                "assigned_module": "module",
            }
        )
        table["region"] = region
        table["celltype_region"] = context
        frames.append(
            table[["region", "celltype_region", "module", "tf"]]
        )

    assignments = pd.concat(frames, ignore_index=True)
    assignments["module"] = assignments["module"].astype(str).str.strip()
    assignments["tf"] = assignments["tf"].astype(str).str.strip()
    assignments = assignments[
        assignments["module"].isin(MODULE_ORDER)
        & assignments["tf"].ne("")
    ].copy()
    assignments["TF_key"] = assignments["tf"].str.upper()
    return assignments


def make_x_layout(contexts):
    x_positions = []
    current = 0.0
    previous_region = None
    region_centers = {}
    region_bounds = {}

    for context in contexts:
        region = context.rsplit("_", 1)[-1]

        if previous_region is not None and region != previous_region:
            current += REGION_GAP

        x_positions.append(current)
        current += 1.0
        previous_region = region

    for region in REGION_ORDER:
        indices = [
            i
            for i, context in enumerate(contexts)
            if context.endswith(f"_{region}")
        ]
        if indices:
            values = [x_positions[i] for i in indices]
            region_centers[region] = (min(values) + max(values)) / 2
            region_bounds[region] = (min(values) - 0.5, max(values) + 0.5)

    return np.asarray(x_positions), region_centers, region_bounds


# ============================================================
# Read association table and select positively associated TFs
# ============================================================
OUTDIR.mkdir(parents=True, exist_ok=True)

assoc = pd.read_csv(ASSOC_FILE, sep="\t")
require_columns(
    assoc,
    [
        "region",
        "celltype_region",
        "TF",
        "lateAD_minus_nonAD",
        SIGNIFICANCE_COLUMN,
    ],
    str(ASSOC_FILE),
)

assoc[SIGNIFICANCE_COLUMN] = pd.to_numeric(
    assoc[SIGNIFICANCE_COLUMN], errors="coerce"
)
assoc["lateAD_minus_nonAD"] = pd.to_numeric(
    assoc["lateAD_minus_nonAD"], errors="coerce"
)

positive = assoc.loc[
    assoc["celltype_region"].isin(CONTEXTS)
    & assoc[SIGNIFICANCE_COLUMN].notna()
    & (assoc[SIGNIFICANCE_COLUMN] < SIGNIFICANCE_CUTOFF)
    & (assoc["lateAD_minus_nonAD"] > 0)
].copy()

positive["TF_key"] = positive["TF"].astype(str).str.upper()
positive = positive.sort_values(
    ["region", "celltype_region", SIGNIFICANCE_COLUMN, "TF"]
)

mode_tag = f"{SIGNIFICANCE_MODE}_lt_{SIGNIFICANCE_CUTOFF:g}"

positive.to_csv(
    OUTDIR / f"positive_AD_associated_TFs_by_celltype_region_{mode_tag}.tsv",
    sep="\t",
    index=False,
)

(
    positive[["TF"]]
    .drop_duplicates()
    .sort_values("TF")
    .to_csv(
        OUTDIR / f"positive_AD_associated_TF_names_{mode_tag}.txt",
        index=False,
        header=False,
    )
)


# ============================================================
# Read module assignments and calculate TF overlaps
# ============================================================
assignments = read_assignments()

# Match each positive TF only within its own celltype_region.
positive_keys = positive[["celltype_region", "TF_key"]].drop_duplicates()

overlap_rows = assignments.merge(
    positive_keys,
    on=["celltype_region", "TF_key"],
    how="inner",
)

overlap_rows = overlap_rows.drop_duplicates(
    ["region", "celltype_region", "module", "TF_key"]
)

overlap_rows["module_number"] = overlap_rows["module"].map(module_number)
overlap_rows = overlap_rows.sort_values(
    ["region", "celltype_region", "module_number", "TF_key"]
).drop(columns="module_number")

overlap_rows.to_csv(
    OUTDIR / f"positive_AD_TF_module_overlap_details_{mode_tag}.tsv",
    sep="\t",
    index=False,
)

heat_long_observed = (
    overlap_rows.groupby(
        ["region", "celltype_region", "module"], observed=True
    )["TF_key"]
    .nunique()
    .rename("n_positive_AD_TF_overlap")
    .reset_index()
)

contexts = CONTEXTS.copy()

full_grid = pd.MultiIndex.from_product(
    [contexts, MODULE_ORDER],
    names=["celltype_region", "module"],
).to_frame(index=False)

full_grid["region"] = full_grid["celltype_region"].str.rsplit(
    "_", n=1
).str[-1]

heat_long = full_grid.merge(
    heat_long_observed,
    on=["region", "celltype_region", "module"],
    how="left",
)

heat_long["n_positive_AD_TF_overlap"] = (
    heat_long["n_positive_AD_TF_overlap"].fillna(0).astype(int)
)

heat_long.to_csv(
    OUTDIR / f"positive_AD_TF_module_overlap_counts_{mode_tag}.tsv",
    sep="\t",
    index=False,
)

heat = (
    heat_long.pivot(
        index="module",
        columns="celltype_region",
        values="n_positive_AD_TF_overlap",
    )
    .reindex(index=MODULE_ORDER, columns=contexts, fill_value=0)
)


# ============================================================
# Get counts for the top bar plot
# ============================================================
derived_bar_values = (
    positive.groupby("celltype_region")["TF_key"]
    .nunique()
    .reindex(contexts, fill_value=0)
    .astype(int)
)

if SIGNIFICANCE_MODE == "fdr" and COUNT_FILE.exists():
    counts = pd.read_csv(COUNT_FILE, sep="\t")
    require_columns(
        counts,
        ["region", "celltype_region", COUNT_COLUMN],
        str(COUNT_FILE),
    )

    supplied_bar_values = (
        counts.drop_duplicates("celltype_region")
        .set_index("celltype_region")[COUNT_COLUMN]
        .reindex(contexts)
    )

    if supplied_bar_values.isna().any():
        missing = list(supplied_bar_values.index[supplied_bar_values.isna()])
        raise ValueError(
            f"Count table lacks these celltype_region values: {missing}"
        )

    supplied_bar_values = supplied_bar_values.astype(int)

    check = pd.DataFrame(
        {
            "celltype_region": contexts,
            "positive_count_from_association_table": derived_bar_values.values,
            "positive_count_from_supplied_count_table": supplied_bar_values.values,
        }
    )
    check["counts_match"] = (
        check["positive_count_from_association_table"]
        == check["positive_count_from_supplied_count_table"]
    )
    check.to_csv(
        OUTDIR / f"positive_TF_barplot_count_check_{mode_tag}.tsv",
        sep="\t",
        index=False,
    )

    if check["counts_match"].all():
        bar_values = supplied_bar_values
        print("FDR count table matches association-derived counts.")
    else:
        print(
            "WARNING: supplied FDR count table disagrees with the current "
            "association table. Association-derived counts will be plotted."
        )
        print(check.loc[~check["counts_match"]].to_string(index=False))
        bar_values = derived_bar_values
else:
    # Nominal mode always derives counts from the association table.
    bar_values = derived_bar_values

bar_count_table = pd.DataFrame(
    {
        "region": [context.rsplit("_", 1)[-1] for context in contexts],
        "celltype": [get_celltype(context) for context in contexts],
        "celltype_region": contexts,
        "n_positive_TFs": bar_values.values,
        "significance_mode": SIGNIFICANCE_MODE,
        "cutoff": SIGNIFICANCE_CUTOFF,
    }
)
bar_count_table.to_csv(
    OUTDIR / f"positive_TF_barplot_counts_{mode_tag}.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# Combined top bar plot and heatmap
# ============================================================
x, region_centers, region_bounds = make_x_layout(contexts)
heat_array = heat.to_numpy(dtype=float)
bar_colors = [
    CELLTYPE_COLORS.get(get_celltype(context), "#808080")
    for context in contexts
]

fig = plt.figure(figsize=(18, 10.5))
gs = gridspec.GridSpec(
    2,
    1,
    height_ratios=[1.25, 4.8],
    hspace=0.035,
)

ax_bar = fig.add_subplot(gs[0])
ax_heat = fig.add_subplot(gs[1], sharex=ax_bar)

ax_bar.bar(
    x,
    bar_values.values,
    width=0.82,
    color=bar_colors,
    edgecolor="black",
    linewidth=0.45,
)

label_offset = max(float(bar_values.max()) * 0.018, 0.25)
for x_value, value in zip(x, bar_values.values):
    if value > 0:
        ax_bar.text(
            x_value,
            value + label_offset,
            str(value),
            ha="center",
            va="bottom",
            fontsize=7.5,
        )

ax_bar.set_ylabel(
    f"Positive AD-associated TFs\n({SIGNIFICANCE_LABEL} < {SIGNIFICANCE_CUTOFF:g})",
    fontsize=10,
)
ax_bar.set_title(
    "Positive AD-associated TFs and overlap with selected cell type-region regulatory modules",
    fontsize=15,
    pad=34,
)
ax_bar.spines[["top", "right", "bottom"]].set_visible(False)
ax_bar.tick_params(axis="x", bottom=False, labelbottom=False)
ax_bar.grid(axis="y", linestyle=":", alpha=0.35)
ax_bar.set_ylim(0, max(float(bar_values.max()) * 1.25, 1))

present_celltypes = [
    celltype
    for celltype in CELLTYPE_ORDER
    if any(get_celltype(context) == celltype for context in contexts)
]
legend_handles = [
    Patch(
        facecolor=CELLTYPE_COLORS[celltype],
        edgecolor="black",
        linewidth=0.45,
        label=celltype,
    )
    for celltype in present_celltypes
]
ax_bar.legend(
    handles=legend_handles,
    title="Cell type",
    frameon=False,
    ncol=len(legend_handles),
    loc="upper center",
    bbox_to_anchor=(0.5, 1.15),
    fontsize=8.5,
    title_fontsize=9,
)

base_cmap = LinearSegmentedColormap.from_list(
    "white_to_blue",
    ["#FFFFFF", "#DCEAF7", "#8AB6D9", "#2F6DA3"],
)

vmax = max(int(np.nanmax(heat_array)), 1)
image = None

for column_index, x_value in enumerate(x):
    column = heat_array[:, column_index : column_index + 1]

    image = ax_heat.imshow(
        column,
        aspect="auto",
        interpolation="nearest",
        cmap=base_cmap,
        vmin=0,
        vmax=vmax,
        extent=(
            x_value - 0.41,
            x_value + 0.41,
            len(MODULE_ORDER) - 0.5,
            -0.5,
        ),
    )

    for row_index, value in enumerate(column[:, 0]):
        if value > 0:
            text_color = "white" if value >= 0.62 * vmax else "black"
            ax_heat.text(
                x_value,
                row_index,
                int(value),
                ha="center",
                va="center",
                fontsize=7,
                color=text_color,
            )

ax_heat.set_yticks(np.arange(len(MODULE_ORDER)))
ax_heat.set_yticklabels(MODULE_ORDER, fontsize=9)
ax_heat.set_ylabel("Cell type-region module", fontsize=11)
ax_heat.set_xticks(x)
ax_heat.set_xticklabels(contexts, rotation=90, fontsize=8)
ax_heat.tick_params(length=0)
ax_heat.spines[["top", "right", "left", "bottom"]].set_visible(False)

# Draw separators only between regions that are present in the selected run.
present_regions = [region for region in REGION_ORDER if region in region_bounds]
for region in present_regions[:-1]:
    right = region_bounds[region][1]
    ax_bar.axvline(right, color="0.78", linewidth=0.9)
    ax_heat.axvline(right, color="0.78", linewidth=0.9)

for region, center in region_centers.items():
    ax_bar.text(
        center,
        1.03,
        region,
        transform=ax_bar.get_xaxis_transform(),
        ha="center",
        va="bottom",
        fontsize=11,
        fontweight="bold",
    )

ax_heat.set_xlim(x.min() - 0.65, x.max() + 0.65)

cbar = fig.colorbar(
    image,
    ax=[ax_bar, ax_heat],
    fraction=0.018,
    pad=0.012,
)
cbar.set_label(
    "Number of overlapping positive AD-associated TFs",
    fontsize=10,
)

output_stem = f"positive_AD_TF_celltype_region_module_overlap_heatmap_with_top_bars_{mode_tag}"
fig.savefig(OUTDIR / f"{output_stem}.pdf", bbox_inches="tight")
fig.savefig(
    OUTDIR / f"{output_stem}.png",
    dpi=300,
    bbox_inches="tight",
)
plt.close(fig)

print(f"Significance mode: {SIGNIFICANCE_MODE}")
print(f"Significance column: {SIGNIFICANCE_COLUMN}")
print(f"Cutoff: {SIGNIFICANCE_CUTOFF}")
print(f"Selected celltype-region contexts: {len(contexts)}")
print(f"Positive association rows in selected contexts: {len(positive):,}")
print(f"Unique positive TF names: {positive['TF'].nunique():,}")
print(f"Output directory: {OUTDIR}")
