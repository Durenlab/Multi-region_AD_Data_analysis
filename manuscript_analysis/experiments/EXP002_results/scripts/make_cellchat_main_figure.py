#!/usr/bin/env python
"""Create a compact all-region CellChat overview and assembled main figure."""
from __future__ import annotations

import base64
import hashlib
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "experiments" / "EXP002_results" / "finalized"
TABLES = OUT / "tables"
FIGURES = OUT / "figures"
REGIONS = ["AG", "EC", "HC", "MTC", "PFC", "TH"]
CELLTYPES = ["Exc", "Inh", "Ast", "Oli", "OPC", "Mic_Immune", "Vasc"]
PATHS = ["nonAD", "earlyAD", "lateAD"]
COLORS = {"nonAD": "#4477AA", "earlyAD": "#EE9944", "lateAD": "#CC6677"}
SEED = 20260929


def draw_overview(fig, spec, summary):
    grid = GridSpecFromSubplotSpec(2, 4, subplot_spec=spec,
                                   width_ratios=[1, 1, 1, .055],
                                   wspace=.24, hspace=.38)
    axes = [fig.add_subplot(grid[i // 3, i % 3]) for i in range(6)]
    cax = fig.add_subplot(grid[:, 3])
    pair_order = [f"{s} → {t}" for s in CELLTYPES for t in CELLTYPES]
    pivot = summary.pivot(index="pair", columns="region", values="median_beta").reindex(index=pair_order, columns=REGIONS)
    max_abs = float(np.nanquantile(np.abs(pivot.to_numpy(float)), .98))
    max_abs = max(max_abs, 1e-10)
    norm = Normalize(vmin=-max_abs, vmax=max_abs)
    cmap = plt.get_cmap("coolwarm")
    for panel_index, (ax, region) in enumerate(zip(axes, REGIONS)):
        region_data = summary[summary.region == region].set_index(["source", "target"])
        mat = np.full((len(CELLTYPES), len(CELLTYPES)), np.nan)
        hits = np.zeros_like(mat, dtype=int)
        for i, source in enumerate(CELLTYPES):
            for j, target in enumerate(CELLTYPES):
                if (source, target) in region_data.index:
                    row = region_data.loc[(source, target)]
                    mat[i, j] = row.median_beta
                    hits[i, j] = int(row.n_q_lt_0_05)
        ax.imshow(mat, cmap=cmap, norm=norm, aspect="equal")
        ax.set_title(region, fontsize=10, weight="bold", pad=4)
        ax.set_xticks(range(len(CELLTYPES)), CELLTYPES, rotation=50, ha="right", fontsize=6.5)
        ax.set_yticks(range(len(CELLTYPES)), CELLTYPES, fontsize=6.5)
        if panel_index >= 3:
            ax.set_xlabel("Receiver", fontsize=7)
        if panel_index in (0, 3):
            ax.set_ylabel("Sender", fontsize=7)
        ax.set_xticks(np.arange(-.5, len(CELLTYPES), 1), minor=True)
        ax.set_yticks(np.arange(-.5, len(CELLTYPES), 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=.8)
        ax.tick_params(which="minor", bottom=False, left=False)
        for i in range(len(CELLTYPES)):
            for j in range(len(CELLTYPES)):
                n_hit = hits[i, j]
                if n_hit:
                    ax.scatter(j, i, s=25 + 30 * n_hit, color="#202020", alpha=.88,
                               edgecolor="white", linewidth=.35, zorder=4)
                    ax.text(j, i, str(n_hit), ha="center", va="center", color="white",
                            fontsize=5.8, weight="bold", zorder=5)
        ax.set_xlim(-.5, len(CELLTYPES)-.5)
        ax.set_ylim(len(CELLTYPES)-.5, -.5)
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cb.set_label("Median pathology slope\nper candidate interaction", fontsize=7.5)
    cb.ax.tick_params(labelsize=6.5)
    return axes


def draw_significant(fig, spec, selected, scores):
    sig = selected[selected.q_value_within_family < .05].copy()
    ncols = 2
    nrows = max(1, int(np.ceil(len(sig) / ncols)))
    grid = GridSpecFromSubplotSpec(nrows, ncols, subplot_spec=spec, wspace=.25, hspace=.42)
    axes = [fig.add_subplot(grid[i // ncols, i % ncols]) for i in range(nrows * ncols)]
    if sig.empty:
        axes[0].text(.5, .5, "No fixed selected features pass pair-specific BH q < 0.05",
                     transform=axes[0].transAxes, ha="center", va="center", wrap=True, color="#555555")
        axes[0].axis("off")
        for ax in axes[1:]:
            ax.set_visible(False)
        return axes, 0
    rng = np.random.default_rng(SEED)
    for ax, row in zip(axes, sig.itertuples(index=False)):
        obs = scores[(scores.region == row.region) & (scores.source == row.source) &
                     (scores.target == row.target) & (scores.interaction_id == row.interaction_id)]
        for x, path in enumerate(PATHS):
            values = obs.loc[obs.pathology == path, "score"].to_numpy(float)
            ax.scatter(x + rng.normal(0, .055, len(values)), values,
                       color=COLORS[path], alpha=.8, s=22, edgecolor="white",
                       linewidth=.3, zorder=3)
            ax.text(x, .98, f"n={len(values)}", transform=ax.get_xaxis_transform(),
                    ha="center", va="top", fontsize=8)
        fitted = [row.prediction_nonAD_mean_covars,
                  row.prediction_earlyAD_mean_covars,
                  row.prediction_lateAD_mean_covars]
        if np.isfinite(fitted).all():
            ax.plot(range(3), fitted, color="#777777", linestyle="--", linewidth=1.8)
        ax.set_title(f"{row.region}: {row.source} → {row.target} | {row.interaction_id}\n"
                     f"β={row.beta_pathology_stage:.3g} [{row.ci_low:.3g}, {row.ci_high:.3g}]; "
                     f"pair-specific BH q={row.q_value_within_family:.3g}", fontsize=8.5)
        ax.set_xticks(range(3), PATHS, rotation=12, ha="right")
        ax.set_ylabel("CellChat probability", fontsize=8)
        ax.grid(axis="y", color="#dddddd", linewidth=.55)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes[len(sig):]:
        ax.set_visible(False)
    return axes, len(sig)


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)
    summary = pd.read_csv(TABLES / "all_pair_overview_summary.tsv", sep="\t")
    selected = pd.read_csv(TABLES / "fixed_selected_features_refreshed.tsv", sep="\t")
    scores = pd.read_csv(TABLES / "fixed_selected_feature_donor_scores.tsv", sep="\t")
    if len(summary) != len(REGIONS) * len(CELLTYPES) ** 2:
        raise ValueError(f"Expected {len(REGIONS) * len(CELLTYPES) ** 2} region × directed-pair rows; got {len(summary)}")
    sig_n = int((selected.q_value_within_family < .05).sum())
    overview_fig = plt.figure(figsize=(16.5, 8.1))
    overview_grid = overview_fig.add_gridspec(1, 1, left=.075, right=.94, bottom=.09, top=.88)
    draw_overview(overview_fig, overview_grid[0], summary)
    overview_fig.suptitle("CellChat pathology overview across six regions", fontsize=15, y=.975)
    overview_fig.text(.5, .925,
                      "Each matrix cell is a directed sender → receiver pair. Color: median donor-level slope across ligand–receptor tests. "
                      "Dark dot and number: count of interactions with pair-specific BH q < 0.05.",
                      ha="center", va="center", fontsize=8.5)
    overview_fig.savefig(FIGURES / "all_region_compact_overview.pdf", bbox_inches="tight")
    overview_fig.savefig(FIGURES / "all_region_compact_overview.png", dpi=300, bbox_inches="tight")
    plt.close(overview_fig)

    panel_rows = max(1, int(np.ceil(sig_n / 2)))
    fig = plt.figure(figsize=(17, 8 + 3.1 * panel_rows))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.04, 0.55 + 0.5 * panel_rows], hspace=.30,
                             left=.075, right=.94, bottom=.055, top=.94)
    draw_overview(fig, outer[0], summary)
    _, n_sig = draw_significant(fig, outer[1], selected, scores)
    fig.text(.015, .935, "A", fontsize=16, weight="bold")
    fig.text(.015, .475, "B", fontsize=16, weight="bold")
    fig.suptitle("Donor-level CellChat communication changes with AD pathology", fontsize=16, y=.99)
    fig.savefig(FIGURES / "cellchat_main_figure.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / "cellchat_main_figure.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    legend_path = OUT / "methods" / "figure_legends.md"
    legend = legend_path.read_text()
    if "**Figure 1. CellChat pathology overview" in legend:
        legend = legend.replace(legend[legend.index("**Figure 1. CellChat pathology overview"):], "")
    legend = ("**Figure 1. CellChat pathology overview across six regions.** Panel A displays all 49 directed sender → receiver cell-class pairs in each region. Tile color represents the median simple OLS pathology-stage slope across estimable ligand–receptor interactions in that region and pair; the dark-dot number gives the count with pair-specific BH q < 0.05. Pathology is coded nonAD=0, earlyAD=1, lateAD=2. The donor is the replicate, cell groups require >20 cells, and donor networks must have completed. P values use HC3 robust standard errors without donor covariates; BH adjustment is done separately for each region × sender/receiver pair across the estimable ligand–receptor interactions in that pair. Panel B shows all fixed selected features that pass their pair-specific BH threshold, with donor scores and simple linear trends.\n\n" + legend)
    legend_path.write_text(legend)

    report_path = OUT / "report" / "cellchat_final_report.html"
    report = report_path.read_text()
    image = base64.b64encode((FIGURES / "cellchat_main_figure.png").read_bytes()).decode("ascii")
    heading = "<h2>Figures</h2>"
    block = (heading + "<h2>Main figure: regional overview and significant selected features</h2>"
             f"<img alt='CellChat main figure with six-region overview and significant features' src='data:image/png;base64,{image}'>")
    if heading not in report:
        raise ValueError("Could not find Figures section in finalized report")
    existing_block = (r"<h2>Figures</h2>(?:<h2>Main figure: regional overview and significant selected features</h2>"
                      r"<img alt='CellChat main figure with six-region overview and significant features' "
                      r"src='data:image/png;base64,[^']+'>)+")
    report, n_existing = re.subn(existing_block, block, report, count=1)
    if n_existing == 0:
        report = report.replace(heading, block, 1)
    report_path.write_text(report)

    script = Path(__file__).resolve()
    revision = {"updated_utc": datetime.now(timezone.utc).isoformat(),
                "script": str(script),
                "script_sha256": hashlib.sha256(script.read_bytes()).hexdigest(),
                "overview_layout": "six small-multiple 7x7 directed sender-by-receiver matrices",
                "overview_rows": int(len(summary)),
                "significant_selected_features_in_panel_b": n_sig,
                "main_figure": "cellchat_main_figure.pdf/png",
                "statistical_results_recomputed": False}
    (OUT / "config" / "main_figure_revision.json").write_text(json.dumps(revision, indent=2) + "\n")
    cfg_path = OUT / "config" / "analysis_config.json"
    cfg = json.loads(cfg_path.read_text())
    cfg["main_figure_revision"] = revision
    cfg_path.write_text(json.dumps(cfg, indent=2) + "\n")
    print("MAIN_FIGURE_COMPLETE")
    print("overview region-pair rows:", len(summary))
    print("significant selected features in panel B:", n_sig)
    print("output:", FIGURES / "cellchat_main_figure.pdf")
    print("report refreshed:", report_path)


if __name__ == "__main__":
    main()
