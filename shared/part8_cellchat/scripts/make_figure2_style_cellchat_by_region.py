#!/usr/bin/env python3
"""Create a Figure 2-style CellChat summary separately for each region."""
from pathlib import Path
import os
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Circle, FancyArrowPatch, Patch
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


BASE = Path(os.environ["CELLCHAT_WORK_DIR"])
TABLES = BASE / "finalized/tables"
FIGURES = BASE / "finalized/figures"
REGIONS = ["AG", "EC", "HC", "MTC", "PFC", "TH"]
CELLTYPES = ["Exc", "Inh", "Ast", "Oli", "OPC", "Mic_Immune"]
CELL_COLORS = {
    "Exc": "#0072B2", "Inh": "#E69F00", "Ast": "#009E73",
    "Oli": "#D55E00", "OPC": "#CC79A7", "Mic_Immune": "#6B5B53",
}
UP_COLOR = "#C44E52"
DOWN_COLOR = "#4C72B0"

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                     "pdf.fonttype": 42, "ps.fonttype": 42})


def read_results():
    tests = pd.read_csv(TABLES / "pathology_progression_tests.tsv", sep="\t")
    tests = tests[(tests.status == "tested") & tests.p_value.notna()
                  & tests.beta_pathology_stage.notna()
                  & tests.source.isin(CELLTYPES) & tests.target.isin(CELLTYPES)].copy()
    meta = []
    for region in REGIONS:
        candidate_path = BASE / "tables" / f"{region}_candidate_interactions.tsv"
        c = pd.read_csv(candidate_path, sep="\t", usecols=[
            "interaction_id", "pathway_name", "ligand", "receptor"
        ])
        c["region"] = region
        meta.append(c)
    meta = pd.concat(meta, ignore_index=True).drop_duplicates(["region", "interaction_id"])
    d = tests.merge(meta, on=["region", "interaction_id"], how="left", validate="many_to_one")
    if d.pathway_name.isna().any():
        raise ValueError("Missing pathway/ligand-receptor metadata for tested results")
    d["lr"] = d.ligand.fillna(d.interaction_id).astype(str) + "–" + d.receptor.fillna("").astype(str).str.replace("_", "/", regex=False)
    d["neglog10p"] = -np.log10(d.p_value.clip(lower=np.nextafter(0.0, 1.0)))
    d["abs_slope"] = d.beta_pathology_stage.abs()
    # Rank both ingredients within region, then multiply the percentiles.
    # This prevents the p-value and raw slope units from arbitrarily setting scale.
    d["p_percentile"] = d.groupby("region").neglog10p.rank(pct=True, method="average")
    d["slope_percentile"] = d.groupby("region").abs_slope.rank(pct=True, method="average")
    d["combined_score"] = d.p_percentile * d.slope_percentile
    d["rank_combined"] = d.groupby("region").combined_score.rank(
        method="first", ascending=False).astype(int)
    return d.sort_values(["region", "rank_combined"])


def panel_network(ax, d, region):
    theta = np.linspace(np.pi / 2, np.pi / 2 - 2 * np.pi, len(CELLTYPES), endpoint=False)
    pos = {cell: np.array([np.cos(t), np.sin(t)]) for cell, t in zip(CELLTYPES, theta)}
    edges = d.groupby(["source", "target"], as_index=False).agg(
        n_features=("interaction_id", "size"))
    if edges.empty:
        ax.text(.5, .5, "No pair-specific FDR-significant interactions",
                transform=ax.transAxes, ha="center", va="center", color="#666666")
    else:
        max_n = int(edges.n_features.max())
        for edge in edges.itertuples(index=False):
            start_cell, end_cell = pos[edge.source], pos[edge.target]
            if edge.source == edge.target:
                start = start_cell + [-.09, .06]
                end = end_cell + [.09, .06]
                rad = -2.0
            else:
                start = start_cell * 1.08
                end = end_cell * 1.08
                rad = .13 if CELLTYPES.index(edge.source) < CELLTYPES.index(edge.target) else -.13
            ax.add_patch(FancyArrowPatch(
                start, end, arrowstyle="-|>", mutation_scale=9,
                connectionstyle=f"arc3,rad={rad}",
                linewidth=1.0 + 1.5 * np.sqrt(edge.n_features / max_n),
                color=CELL_COLORS[edge.source], alpha=.8, zorder=1))
            if edge.source != edge.target:
                mid = (start_cell + end_cell) / 2
                ax.text(mid[0] * .9, mid[1] * .9, str(edge.n_features),
                        ha="center", va="center", fontsize=7, color="#222222",
                        bbox={"facecolor": "white", "edgecolor": "none", "alpha": .85, "pad": .5}, zorder=3)
    for cell, xy in pos.items():
        ax.add_patch(Circle(xy, .12, facecolor=CELL_COLORS[cell], edgecolor="white", linewidth=.8, zorder=4))
        label_xy = xy * 1.28
        ha = "left" if label_xy[0] > .1 else "right" if label_xy[0] < -.1 else "center"
        ax.text(label_xy[0], label_xy[1], cell, fontsize=7, ha=ha, va="center", zorder=5)
    ax.set_xlim(-1.55, 1.55)
    ax.set_ylim(-1.45, 1.42)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(f"A  Significant directed paths · {len(d)} features", loc="left", weight="bold", fontsize=10)


def panel_matrix(ax, d):
    # Using one shared order on both axes places each cell type on the diagonal.
    mat = np.zeros((len(CELLTYPES), len(CELLTYPES)), dtype=int)
    for row in d.itertuples(index=False):
        mat[CELLTYPES.index(row.source), CELLTYPES.index(row.target)] += 1
    im = ax.imshow(mat, cmap="YlOrRd", vmin=0, vmax=max(1, int(mat.max())), aspect="equal")
    ax.set_xticks(range(len(CELLTYPES)), CELLTYPES, rotation=38, ha="right", fontsize=7)
    ax.set_yticks(range(len(CELLTYPES)), CELLTYPES, fontsize=7)
    ax.set_xlabel("Receiver")
    ax.set_ylabel("Sender")
    ax.set_title("B  Significant interaction-feature counts", loc="left", weight="bold", fontsize=10)
    for i in range(len(CELLTYPES)):
        for j in range(len(CELLTYPES)):
            val = mat[i, j]
            ax.text(j, i, str(val), ha="center", va="center", fontsize=7,
                    color="white" if val >= max(2, mat.max() * .55) else "#333333")
    ax.set_xticks(np.arange(-.5, len(CELLTYPES), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(CELLTYPES), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=.8)
    ax.tick_params(which="minor", bottom=False, left=False)


def panel_trend(ax, d):
    if d.empty:
        ax.text(.5, .5, "No significant ligand–receptor features",
                transform=ax.transAxes, ha="center", va="center", color="#666666")
        ax.set_title("C  Donor-level linear trend", loc="left", weight="bold", fontsize=10)
        ax.set_axis_off()
        return
    d = d.sort_values("combined_score", ascending=True).reset_index(drop=True)
    y = np.arange(len(d))
    labels = [textwrap.fill(f"{r.source}→{r.target} | {r.lr}", width=25)
              for r in d.itertuples(index=False)]
    for i, row in enumerate(d.itertuples(index=False)):
        color = UP_COLOR if row.beta_pathology_stage > 0 else DOWN_COLOR
        ax.barh(i, row.beta_pathology_stage, height=.56, color=color, zorder=2)
        xerr = np.array([[row.beta_pathology_stage - row.ci_low],
                         [row.ci_high - row.beta_pathology_stage]])
        ax.errorbar(row.beta_pathology_stage, i, xerr=xerr, fmt="none", ecolor="#242424",
                    elinewidth=.65, capsize=1.7, zorder=3)
    ax.axvline(0, color="#444444", linewidth=.7, zorder=1)
    ax.set_yticks(y, labels, fontsize=5.6)
    ax.set_xlabel("Slope β per pathology step", fontsize=7)
    ax.tick_params(axis="x", labelsize=6)
    ax.grid(axis="x", color="#e8e8e8", linewidth=.5)
    ax.set_axisbelow(True)
    ax.set_title("C  Donor-level linear trend", loc="left", weight="bold", fontsize=9)


def plot_region_panels_ab(all_results, region):
    d = all_results[(all_results.region == region) &
                    (all_results.q_value_within_family < .05)].copy()
    fig = plt.figure(figsize=(15, 7.2))
    gs = GridSpec(1, 5, figure=fig, wspace=.55)
    ax_a = fig.add_subplot(gs[0, :2])
    panel_network(ax_a, d, region)
    ax_b = fig.add_subplot(gs[0, 2:])
    panel_matrix(ax_b, d)
    fig.suptitle(f"{region}: significant donor-level CellChat features", fontsize=14, weight="bold", y=.98)
    fig.text(.5, .92,
             "Simple donor-level OLS (nonAD=0, earlyAD=1, lateAD=2; HC3 SE); pair-specific BH FDR < 0.05. "
             "Six nonvascular cell classes.", ha="center", fontsize=8)
    fig.text(.5, .025,
             "A: edge width/labels count significant interaction features; edge color identifies sender. "
             "B: identical cell-type order on sender and receiver axes aligns matching types diagonally. "
             "Inferred communication associations do not establish causality.",
             ha="center", fontsize=8)
    for fmt, kwargs in [("pdf", {}), ("png", {"dpi": 300})]:
        fig.savefig(FIGURES / f"figure2_style_{region}_pathology.{fmt}", bbox_inches="tight", **kwargs)
    plt.close(fig)


def plot_all_region_trends(all_results):
    top = all_results.groupby("region", sort=False).head(10).copy()
    top = top.sort_values(["region", "rank_combined"], ascending=[True, True]).reset_index(drop=True)
    top.to_csv(TABLES / "figure2_style_all_regions_top10_features.tsv", sep="\t", index=False)
    for region in REGIONS:
        top[top.region.eq(region)].to_csv(
            TABLES / f"figure2_style_{region}_top10_features.tsv", sep="\t", index=False)

    # Ten rows per region, one combined panel. The signed coefficient determines
    # direction; p/FDR are printed in a separate aligned text column.
    n = len(top)
    y = np.arange(n)
    row_height = .38
    fig_height = max(14, row_height * n + 3)
    fig = plt.figure(figsize=(15, fig_height))
    gs = GridSpec(1, 3, figure=fig, width_ratios=[3.5, 5.0, 1.8], wspace=.08,
                  left=.02, right=.98, top=.97, bottom=.05)
    label_ax = fig.add_subplot(gs[0, 0])
    bar_ax = fig.add_subplot(gs[0, 1], sharey=label_ax)
    stats_ax = fig.add_subplot(gs[0, 2], sharey=label_ax)
    label_ax.set_xlim(0, 1)
    label_ax.set_ylim(n - .5, -.5)
    label_ax.axis("off")
    stats_ax.set_xlim(0, 1)
    stats_ax.set_ylim(n - .5, -.5)
    stats_ax.axis("off")
    max_ci = float(np.nanmax(np.abs(top[["ci_low", "ci_high"]].to_numpy(float))))
    max_ci = max(max_ci, 1e-8)
    bar_ax.set_xlim(-max_ci * 1.12, max_ci * 1.12)
    bar_ax.axvline(0, color="#333333", linewidth=.8, zorder=1)
    bar_ax.grid(axis="x", color="#e5e5e5", linewidth=.6)
    bar_ax.set_axisbelow(True)
    for i, row in enumerate(top.itertuples(index=False)):
        label = f"{row.region}  |  {row.source}→{row.target}  |  {row.lr}"
        label_ax.text(.99, i, label, ha="right", va="center", fontsize=11.2)
        bar_ax.barh(i, row.beta_pathology_stage, height=.28, color="#9E9E9E",
                    edgecolor="#686868", linewidth=.35, zorder=2)
        xerr = np.array([[row.beta_pathology_stage - row.ci_low],
                         [row.ci_high - row.beta_pathology_stage]])
        bar_ax.errorbar(row.beta_pathology_stage, i, xerr=xerr, fmt="none",
                        ecolor="#333333", elinewidth=.8, capsize=2.0, zorder=3)
        stats_ax.text(.02, i, f"p={row.p_value:.2g}   FDR={row.q_value_within_family:.2g}",
                      ha="left", va="center", fontsize=11.2)
    bar_ax.set_yticks(y)
    bar_ax.tick_params(axis="y", left=False, labelleft=False)
    bar_ax.tick_params(axis="x", labelsize=11.2)
    bar_ax.set_xlabel("Pathology-stage slope β (CellChat probability per step)", fontsize=14)
    for boundary in range(10, n, 10):
        for ax in [label_ax, bar_ax, stats_ax]:
            ax.axhline(boundary - .5, color="#bdbdbd", linewidth=.8, zorder=4)
    fig.suptitle("Donor-level linear trend", fontsize=18, weight="bold", y=.995)
    fig.text(.5, .012,
             "Top 10 tested interaction × sender/receiver features per region, ranked by the product of "
             "within-region percentile ranks for −log10(p) and |slope|. All bars are gray; direction retains the signed slope. "
             "Whiskers show 95% HC3 confidence intervals; FDR is BH within region × sender/receiver pair.",
             ha="center", va="bottom", fontsize=11.2, wrap=True)
    for fmt, kwargs in [("pdf", {}), ("png", {"dpi": 300})]:
        fig.savefig(FIGURES / f"figure2_style_all_regions_panel_C.{fmt}", bbox_inches="tight", **kwargs)
    plt.close(fig)


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)
    all_results = read_results()
    all_results = all_results.sort_values(["region", "rank_combined"])
    for region in REGIONS:
        plot_region_panels_ab(all_results, region)
    plot_all_region_trends(all_results)
    all_sig = all_results[all_results.q_value_within_family < .05]
    summary = (all_sig.groupby("region", as_index=False)
               .agg(n_significant_features=("interaction_id", "size"),
                    n_unique_interactions=("interaction_id", "nunique")))
    summary = summary.set_index("region").reindex(REGIONS, fill_value=0).reset_index()
    summary["n_displayed_panel_c"] = 10
    summary.to_csv(TABLES / "figure2_style_region_summary.tsv", sep="\t", index=False)
    legend_path = BASE / "finalized/methods/figure_legends.md"
    legend_path.parent.mkdir(parents=True, exist_ok=True)
    legend = legend_path.read_text() if legend_path.exists() else ""
    marker = "**Regional Figure 2-style CellChat panels.**"
    if marker in legend:
        legend = legend[:legend.index(marker)].rstrip()
    legend += ("\n\n**Regional Figure 2-style CellChat panels.** Panels A and B are saved separately for each region. "
               "Panel A displays all pair-specific FDR-significant directed paths among six nonvascular classes; edge color "
               "identifies sender and edge width/labels count significant interaction features. Panel B counts the "
               "same features in a sender-by-receiver matrix with the same six-cell-type order on both axes, so "
               "matching cell types lie on the diagonal. One combined 15-inch-wide Panel C, titled ‘Donor-level "
               "linear trend’, shows the top ten tested interaction × sender/receiver features per region. Ranking "
               "uses the product of within-region percentile ranks for −log10(p) and |slope|, computed across all "
               "tested features. Bars are gray, extend according to signed slope, and are half the prior thickness; "
               "font sizes are doubled. Nominal p-values and pair-specific BH FDR values are printed for each row; "
               "whiskers show 95% HC3 confidence intervals. Donors are the replicates; pathology is coded "
               "nonAD=0, earlyAD=1, lateAD=2.\n")
    legend_path.write_text(legend)
    print("Wrote regional A/B panels and combined all-region panel C")
    print(summary.to_string(index=False))
    ec_top = all_results[all_results.region.eq("EC")].head(10)
    ec_spp1 = ec_top[ec_top.pathway_name.eq("SPP1")]
    print("EC SPP1 features in combined-score top ten:", len(ec_spp1))
    if not ec_spp1.empty:
        print(ec_spp1[["source", "target", "interaction_id", "rank_combined", "beta_pathology_stage",
                       "p_value", "q_value_within_family"]].to_string(index=False))


if __name__ == "__main__":
    main()
