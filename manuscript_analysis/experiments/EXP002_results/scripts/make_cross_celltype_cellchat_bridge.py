#!/usr/bin/env python3
"""Create a cautious visual bridge from cross-celltype modules to CellChat."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "experiments/EXP002_results/finalized/figures"
TABLES = ROOT / "experiments/EXP002_results/finalized/tables"
EXP004 = ROOT / "experiments/EXP004_results"

regions = ["AG", "EC", "HC", "MTC", "PFC", "TH"]
modules = pd.read_csv(EXP004 / "tables/module_pathology.tsv", sep="\t")
cellchat = pd.read_csv(TABLES / "all_pair_overview_summary.tsv", sep="\t")
donors = pd.read_csv(
    EXP004 / "ec_comparison/tables/cellchat_figure_donor_scores.tsv", sep="\t"
)
donors = donors[
    (donors["source"] == "OPC")
    & (donors["target"] == "Ast")
    & (donors["interaction_id"] == "COL6A1_CD44")
].copy()
stats = pd.read_csv(TABLES / "pathology_progression_tests.tsv", sep="\t")
stats = stats[
    (stats["region"] == "EC")
    & (stats["source"] == "OPC")
    & (stats["target"] == "Ast")
    & (stats["interaction_id"] == "COL6A1_CD44")
].iloc[0]
pairwise = pd.read_csv(TABLES / "pathology_pairwise_contrasts.tsv", sep="\t")
pairwise = pairwise[
    (pairwise["region"] == "EC")
    & (pairwise["source"] == "OPC")
    & (pairwise["target"] == "Ast")
    & (pairwise["interaction_id"] == "COL6A1_CD44")
    & (pairwise["comparison"] == "lateAD_minus_nonAD")
].iloc[0]

# Earlier cross-celltype screen defined these region/module candidates. Display
# the module associations as a separate evidence layer, not as signaling proof.
candidate_modules = {
    "AG": ["M13"],
    "EC": ["M14", "M17"],
    "HC": [],
    "MTC": ["M6", "M13"],
    "PFC": ["M14", "M18"],
    "TH": [],
}

plt.rcParams.update(
    {"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False,
     "axes.spines.right": False, "pdf.fonttype": 42, "ps.fonttype": 42}
)
fig = plt.figure(figsize=(12.6, 7.0), constrained_layout=True)
gs = fig.add_gridspec(2, 2, height_ratios=[0.95, 1.25], width_ratios=[1, 1.15])
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, :])

# Panel A: previous cross-celltype module screen, retained by region.
ax1.set_title("A  Cross-celltype regulatory modules", loc="left", weight="bold")
ax1.set_xlim(-0.5, 2.7)
ax1.set_ylim(-0.65, len(regions) - 0.35)
ax1.invert_yaxis()
ax1.set_yticks(range(len(regions)), regions)
ax1.set_xticks([])
ax1.set_xlabel("Pathology-associated shared-H module candidates")
for i, region in enumerate(regions):
    names = candidate_modules[region]
    if not names:
        ax1.text(0.05, i, "—", va="center", color="#9a9a9a")
        continue
    for j, module in enumerate(names):
        ax1.scatter(j * 1.05, i, s=130, color="#4477AA", edgecolor="white", zorder=2)
        ax1.text(j * 1.05, i, module, ha="center", va="center", color="white",
                 weight="bold", fontsize=8)
ax1.text(
    0, 1.36,
    "Shared module activity is a coordinated regulatory axis; it does not specify a sender or receiver.",
    transform=ax1.transAxes, fontsize=8, color="#444444", va="top", wrap=True,
)

# Panel B: CellChat regional overview distilled to candidate hits per region.
counts = cellchat.groupby("region", observed=True)["n_q_lt_0_05"].sum().reindex(regions)
bars = ax2.bar(regions, counts.values, color="#CC6677", width=0.68)
ax2.set_title("B  CellChat communication screen", loc="left", weight="bold")
ax2.set_ylabel("Pathology-associated LR features\n(pair-specific q < 0.05)")
ax2.set_ylim(0, max(counts.max() * 1.25, 2))
ax2.grid(axis="y", alpha=0.2)
for b, val in zip(bars, counts.values):
    ax2.text(b.get_x() + b.get_width() / 2, val + 0.18, str(int(val)),
             ha="center", va="bottom", fontsize=9)
ax2.text(
    0, 1.36,
    "Counts are significant ligand–receptor features, not counts of significant cell pairs.",
    transform=ax2.transAxes, fontsize=8, color="#444444", va="top", wrap=True,
)

# Panel C: detailed example, with donors shown individually and the gene-level
# mapping caveat made explicit.
ax3.set_title("C  EC example: pathway-aligned candidate, with donor-level evidence", loc="left", weight="bold")
path_order = ["nonAD", "earlyAD", "lateAD"]
colors = {"nonAD": "#4477AA", "earlyAD": "#EEAA33", "lateAD": "#BB5566"}
rng = np.random.default_rng(31)
for x, group in enumerate(path_order):
    vals = donors.loc[donors.pathology == group, "communication_probability"].to_numpy()
    jitter = rng.uniform(-0.14, 0.14, len(vals))
    ax3.scatter(vals, x + jitter, s=34, alpha=0.82, color=colors[group],
                edgecolor="white", linewidth=0.4, zorder=3)
    ax3.plot(np.median(vals), x, marker="|", markersize=18, markeredgewidth=2,
             color="#222222", zorder=4)
ax3.set_yticks(range(3), [f"{p}  (n={int((donors.pathology == p).sum())})" for p in path_order])
ax3.invert_yaxis()
ax3.set_xlabel("Donor-level CellChat communication probability (OPC → Ast, COL6A1–CD44)")
ax3.set_ylabel("")
ax3.grid(axis="x", alpha=0.2)
ax3.text(
    0.99, 0.95,
    f"EC_M14 shared module: pathology-associated\n"
    f"Final CellChat model: β={stats['beta_pathology_stage']:.4f}, "
    f"progression pair-specific q={stats['q_value_within_family']:.3f}; "
    f"lateAD vs nonAD q={pairwise['q_value_within_family']:.3f}\n"
    "Not significant after pair-specific BH; pathway-level match only: "
    "COL6A1 and CD44 were not leading-edge genes assigned to EC_M14.",
    transform=ax3.transAxes, ha="right", va="top", fontsize=8.5,
    bbox={"boxstyle": "round,pad=0.45", "facecolor": "white", "edgecolor": "#aaaaaa", "alpha": 0.95},
)
ax3.text(
    0.01, -0.27,
    "Each dot is a donor; vertical ticks mark group medians. Absent sparse scores are zero only for eligible completed donor networks. "
    "The plotted association is observational and does not establish signaling or regulatory causality.",
    transform=ax3.transAxes, fontsize=8, color="#444444", va="top", wrap=True,
)

fig.suptitle("From shared regulatory programs to inferred cell–cell communication", fontsize=14, weight="bold")
OUT.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT / "cross_celltype_cellchat_bridge.pdf", bbox_inches="tight")
fig.savefig(OUT / "cross_celltype_cellchat_bridge.png", dpi=220, bbox_inches="tight")
plt.close(fig)

print(f"Saved {OUT / 'cross_celltype_cellchat_bridge.pdf'}")
print(f"Saved {OUT / 'cross_celltype_cellchat_bridge.png'}")
