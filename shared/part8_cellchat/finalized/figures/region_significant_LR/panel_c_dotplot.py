"""Create the significant-region CellChat Panel C dot-and-whisker plot."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REGIONS = ["AG", "EC", "MTC", "TH"]


def find_project_root(start=None):
    """Return the curated CellChat root."""
    return Path(__file__).resolve().parents[3]


def load_significant_results(project_root=None, regions=REGIONS):
    """Load tested results, attach ligand/receptor labels, and retain pair-FDR hits."""
    root = Path(project_root) if project_root else find_project_root()
    base = root
    tables = base / "finalized/tables"
    tests = pd.read_csv(tables / "pathology_progression_tests.tsv", sep="\t")
    tests = tests[
        tests.status.eq("tested")
        & tests.region.isin(regions)
        & tests.source.notna()
        & tests.target.notna()
        & tests.p_value.notna()
        & tests.beta_pathology_stage.notna()
        & tests.q_value_within_family.notna()
        & tests.ci_low.notna()
        & tests.ci_high.notna()
        & tests.source.ne("Vasc")
        & tests.target.ne("Vasc")
        & tests.q_value_within_family.lt(0.05)
    ].copy()

    metadata = []
    for region in regions:
        path = base / "tables" / f"{region}_candidate_interactions.tsv"
        frame = pd.read_csv(path, sep="\t", usecols=["interaction_id", "pathway_name", "ligand", "receptor"])
        frame["region"] = region
        metadata.append(frame)
    metadata = pd.concat(metadata, ignore_index=True).drop_duplicates(["region", "interaction_id"])
    results = tests.merge(metadata, on=["region", "interaction_id"], how="left", validate="many_to_one")
    if results.pathway_name.isna().any():
        missing = results.loc[results.pathway_name.isna(), ["region", "interaction_id"]].drop_duplicates()
        raise ValueError(f"Missing candidate metadata for significant results:\n{missing.to_string(index=False)}")
    results["ligand_receptor"] = (
        results.ligand.fillna(results.interaction_id).astype(str)
        + "–"
        + results.receptor.fillna("").astype(str).str.replace("_", "/", regex=False)
    )
    results["feature_label"] = (
        results.region.astype(str) + "  |  " + results.source.astype(str) + "→" + results.target.astype(str)
        + "  |  " + results.ligand_receptor
    )
    region_order = {region: i for i, region in enumerate(regions)}
    results["_region_order"] = results.region.map(region_order)
    return results.sort_values(["_region_order", "q_value_within_family", "p_value"]).drop(columns="_region_order").reset_index(drop=True)


def plot_dot_whisker(results, output_dir, *, region, width=10, row_height=0.38,
                     label_fontsize=11, tick_fontsize=11, axis_fontsize=14,
                     title_fontsize=18, stats_fontsize=10, dpi=300):
    """Save one regional slope dot plot with HC3 95% confidence intervals."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    results = results[results.region.eq(region)].copy().reset_index(drop=True)
    if results.empty:
        raise ValueError(f"No FDR-significant features found for {region}")

    n = len(results)
    height = max(5.0, row_height * n + 2.4)
    fig, ax = plt.subplots(figsize=(width, height))
    y = np.arange(n)
    xerr = np.vstack([
        results.beta_pathology_stage.to_numpy() - results.ci_low.to_numpy(),
        results.ci_high.to_numpy() - results.beta_pathology_stage.to_numpy(),
    ])
    colors = np.where(results.beta_pathology_stage.to_numpy() >= 0, "#D62728", "#0072B2")
    for i, row in enumerate(results.itertuples(index=False)):
        ax.errorbar(
            row.beta_pathology_stage, i,
            xerr=[[row.beta_pathology_stage - row.ci_low], [row.ci_high - row.beta_pathology_stage]],
            fmt="o", color=colors[i], ecolor=colors[i], markerfacecolor=colors[i],
            markeredgecolor=colors[i], markersize=6, elinewidth=1.2,
            capsize=3, capthick=1.0, zorder=3,
        )
    ax.axvline(0, color="#333333", linewidth=0.9, zorder=1)
    # Shorten this cell-type name for the figure only; the results table keeps
    # its original source/target labels.
    display_celltype = lambda value: "Mic" if value == "Mic_Immune" else value
    labels = results.apply(
        lambda row: f"{display_celltype(row.source)}→{display_celltype(row.target)}  |  {row.ligand_receptor}",
        axis=1,
    )
    ax.set_yticks(y, labels, fontsize=label_fontsize)
    ax.invert_yaxis()
    ax.tick_params(axis="x", labelsize=tick_fontsize)
    ax.set_xlabel("Pathology slope β", fontsize=axis_fontsize)
    ax.set_title(f"{region}", fontsize=title_fontsize, weight="bold", pad=10)
    ax.grid(axis="x", color="#e6e6e6", linewidth=0.7)
    ax.set_axisbelow(True)

    ax.text(1.01, -0.8, "BH FDR", transform=ax.get_yaxis_transform(),
            ha="left", va="bottom", fontsize=stats_fontsize, weight="bold")
    for i, row in enumerate(results.itertuples(index=False)):
        ax.text(1.01, i, f"{row.q_value_within_family:.3g}",
                transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=stats_fontsize,
                clip_on=False)
    ax.set_xlim(
        min(float(results.ci_low.min()), 0) - 0.08 * max(abs(results.ci_low.min()), abs(results.ci_high.max()), 1e-8),
        max(float(results.ci_high.max()), 0) + 0.08 * max(abs(results.ci_low.min()), abs(results.ci_high.max()), 1e-8),
    )
    fig.subplots_adjust(left=0.38, right=0.86, top=0.96, bottom=0.06)
    for suffix, save_kwargs in [("pdf", {}), ("png", {"dpi": dpi})]:
        fig.savefig(output_dir / f"figure2_style_panel_C_{region}_dotplot.{suffix}", bbox_inches="tight", **save_kwargs)
    plt.close(fig)


def save_plot_table(results, output_dir):
    """Save the plotted rows and statistics as a tab-separated table."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    columns = [
        "region", "source", "target", "interaction_id", "pathway_name", "ligand", "receptor",
        "beta_pathology_stage", "ci_low", "ci_high", "p_value", "q_value_within_family",
        "n_donors", "feature_label",
    ]
    available = [column for column in columns if column in results.columns]
    path = output_dir / "figure2_style_all_regions_panel_C_dotplot_table.tsv"
    results[available].to_csv(path, sep="\t", index=False)
    return path
