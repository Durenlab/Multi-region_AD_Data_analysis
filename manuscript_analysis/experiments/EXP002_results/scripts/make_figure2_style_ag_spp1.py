#!/usr/bin/env python3
"""Figure 2-inspired display for the strongest finalized pathway example."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Circle, FancyArrowPatch
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "experiments/EXP002_results"
TABLES = BASE / "finalized/tables"
FIGURES = BASE / "finalized/figures"
CELL_COLORS = {"Oli": "#AA3377", "OPC": "#BBBBBB", "Mic_Immune": "#228833"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                     "pdf.fonttype": 42, "ps.fonttype": 42})


def get_data():
    p = TABLES / "pathology_progression_tests.tsv"
    d = pd.read_csv(p, sep="\t")
    d = d[(d.region == "AG") & (d.status == "tested") &
          (d.q_value_within_family < 0.05)].copy()
    md = pd.read_csv(BASE / "tables/AG_candidate_interactions.tsv", sep="\t",
                     usecols=["interaction_id", "pathway_name", "ligand", "receptor"])
    d = d.merge(md, on="interaction_id", how="left", validate="many_to_one")
    d = d[d.pathway_name == "SPP1"].copy()
    if d.empty:
        raise ValueError("No AG SPP1 features pass pair-specific BH q < 0.05")
    d["lr"] = d["ligand"].astype(str) + " – " + d["receptor"].astype(str).str.replace("_", "/", regex=False)
    d["pair"] = d.source + " → " + d.target
    return d.sort_values(["source", "target", "beta_pathology_stage"], ascending=[True, True, False])


def panel_network(ax, d):
    cells = ["Oli", "OPC", "Mic_Immune"]
    theta = np.linspace(np.pi / 2, np.pi / 2 - 2 * np.pi, len(cells), endpoint=False)
    pos = {cell: np.array([np.cos(t), np.sin(t)]) for cell, t in zip(cells, theta)}
    edges = d.groupby(["source", "target"], as_index=False).agg(n_lr=("interaction_id", "size"))
    for row in edges.itertuples(index=False):
        a, b = pos[row.source], pos[row.target]
        if row.source == row.target:
            start, end, rad = a + [-0.08, 0.06], a + [0.08, 0.06], -2.0
        else:
            start, end = a + 0.13 * a, b + 0.13 * b
            rad = 0.18 if cells.index(row.source) < cells.index(row.target) else -0.18
        ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=11,
                       connectionstyle=f"arc3,rad={rad}", linewidth=1.3 + 0.9 * row.n_lr,
                       color=CELL_COLORS[row.source], alpha=.8, zorder=1))
        if row.source != row.target:
            mid = (a + b) / 2
            ax.text(mid[0] * .88, mid[1] * .88, str(row.n_lr), fontsize=8,
                    color=CELL_COLORS[row.source], weight="bold", ha="center", va="center",
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": .9, "pad": 1})
    for cell, xy in pos.items():
        ax.add_patch(Circle(xy, .13, facecolor=CELL_COLORS[cell], edgecolor="white", linewidth=1, zorder=3))
        lab = xy * 1.27
        ha = "left" if lab[0] > .1 else "right" if lab[0] < -.1 else "center"
        ax.text(lab[0], lab[1], cell, fontsize=8, weight="bold", ha=ha, va="center")
    ax.set_xlim(-1.52, 1.52); ax.set_ylim(-1.42, 1.4); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("A  Significant AG SPP1 paths", loc="left", weight="bold", fontsize=11)
    ax.text(.5, -.05, "Edge color identifies sender; edge width/label gives the number of significant LR features.",
            transform=ax.transAxes, ha="center", va="top", fontsize=7.5, wrap=True)


def panel_matrix(ax, d):
    cells = ["Oli", "Mic_Immune", "OPC"]
    mat = np.zeros((len(cells), len(cells)), dtype=int)
    for r in d.itertuples(index=False):
        mat[cells.index(r.source), cells.index(r.target)] += 1
    im = ax.imshow(mat, cmap="YlOrRd", vmin=0, vmax=max(1, mat.max()), aspect="equal")
    ax.set_xticks(range(len(cells)), cells, rotation=25, ha="right", fontsize=8)
    ax.set_yticks(range(len(cells)), cells, fontsize=8)
    ax.set_xlabel("Receiver"); ax.set_ylabel("Sender")
    ax.set_title("B  Pathology-associated LR feature counts", loc="left", weight="bold", fontsize=11)
    for i in range(len(cells)):
        for j in range(len(cells)):
            if mat[i, j]:
                ax.text(j, i, str(mat[i, j]), ha="center", va="center", weight="bold",
                        color="white" if mat[i, j] >= 2 else "#333333", fontsize=9)
    ax.set_xticks(np.arange(-.5, len(cells), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(cells), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1)
    ax.tick_params(which="minor", bottom=False, left=False)


def panel_effects(ax, d):
    d = d.sort_values("beta_pathology_stage", ascending=True).reset_index(drop=True)
    y = np.arange(len(d))
    labels = [f"{r.pair}  |  {r.lr}" for r in d.itertuples(index=False)]
    xerr = np.vstack([d.beta_pathology_stage - d.ci_low,
                      d.ci_high - d.beta_pathology_stage])
    colors = [CELL_COLORS[s] for s in d.source]
    ax.axvline(0, color="#555555", linewidth=.8)
    for i, (x, lo, hi, color) in enumerate(zip(d.beta_pathology_stage, d.ci_low, d.ci_high, colors)):
        ax.errorbar(x, i, xerr=[[x - lo], [hi - x]], fmt="o", color=color,
                    ecolor=color, elinewidth=1.4, capsize=2.4, markersize=5.5)
    ax.set_yticks(y, labels, fontsize=7.5)
    ax.set_xlabel("Adjusted pathology-stage slope β (CellChat probability units)")
    ax.grid(axis="x", alpha=.2)
    ax.set_title("C  Donor-level simple linear trends", loc="left", weight="bold", fontsize=11)
    right = d.ci_high.max() * 1.08
    ax.set_xlim(min(0, d.ci_low.min() * 1.2), right * 1.3)
    for i, r in enumerate(d.itertuples(index=False)):
        ax.text(right, i, f"q={r.q_value_within_family:.3g}", va="center", fontsize=7.5, color="#444444")


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)
    d = get_data()
    fig = plt.figure(figsize=(14, 9))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.2], height_ratios=[1, 1.25],
                          left=.07, right=.97, top=.88, bottom=.14, hspace=.47, wspace=.30)
    ax_a = fig.add_subplot(gs[0, 0]); panel_network(ax_a, d)
    ax_b = fig.add_subplot(gs[0, 1]); panel_matrix(ax_b, d)
    ax_c = fig.add_subplot(gs[1, :]); panel_effects(ax_c, d)
    fig.suptitle("AG SPP1 communication features associated with AD pathology progression",
                 fontsize=14, weight="bold", y=.97)
    fig.text(.5, .925,
             "Figure 2-inspired summary of finalized donor-level tests. Pathology is coded nonAD=0, earlyAD=1, lateAD=2; "
             "models use pathology stage only (HC3 robust standard errors; no donor covariates).", ha="center", fontsize=8.5)
    fig.text(.5, .025,
             "Only interaction-level tests with pair-specific BH q < 0.05 are displayed. Network counts are not probability-weighted; "
             "results are associations, not proof of ligand–receptor signaling or causality.", ha="center", fontsize=8)
    fig.savefig(FIGURES / "figure2_style_AG_SPP1_pathology.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / "figure2_style_AG_SPP1_pathology.png", dpi=250, bbox_inches="tight")
    plt.close(fig)
    print("Saved Figure 2-inspired AG SPP1 figure; significant donor-level tests:", len(d))


if __name__ == "__main__":
    main()
