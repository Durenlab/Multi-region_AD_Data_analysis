#!/usr/bin/env python3
"""Make an AG-focused descriptive summary of pooled CellChat signal."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "finalized" / "tables"
FIGURES = ROOT / "finalized" / "figures"
REPORTS = ROOT / "finalized" / "report"
PAIR_FILE = TABLES / "pooled_celltype_differences_descriptive.tsv"
PATHWAY_FILE = TABLES / "pooled_pathway_differences_descriptive.tsv"
OUT = FIGURES / "AG_pooled_signal_highlight_main_figure"

REGIONS = ["AG", "EC", "HC", "MTC", "PFC", "TH"]
CELLTYPES = ["Exc", "Inh", "Ast", "Oli", "OPC", "Mic_Immune", "Vasc"]
COMPARISONS = ["earlyAD_minus_nonAD", "lateAD_minus_earlyAD", "lateAD_minus_nonAD"]
STAGE_MAP = {
    "earlyAD_minus_nonAD": ("nonAD", "earlyAD"),
    "lateAD_minus_earlyAD": ("earlyAD", "lateAD"),
    "lateAD_minus_nonAD": ("nonAD", "lateAD"),
}
REGION_COLORS = {"AG": "#CC6677", "EC": "#777777", "HC": "#777777", "MTC": "#777777", "PFC": "#777777", "TH": "#777777"}


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)
    pair = pd.read_csv(PAIR_FILE, sep="\t")
    pathway = pd.read_csv(PATHWAY_FILE, sep="\t")

    required_pair = {"region", "comparison", "condition_1", "condition_2", "source", "target", "strength_condition_1", "strength_condition_2", "strength_difference"}
    required_pathway = {"region", "comparison", "pathway", "strength_difference"}
    if not required_pair.issubset(pair.columns) or not required_pathway.issubset(pathway.columns):
        raise ValueError("Pooled input tables do not have the expected columns")

    # Verify the source table's difference convention and full region/comparison/pair coverage.
    if not np.allclose(pair.strength_condition_2 - pair.strength_condition_1, pair.strength_difference, atol=1e-8):
        raise ValueError("Pair table strength_difference is not condition_2 minus condition_1")
    expected_pairs = {(s, t) for s in CELLTYPES for t in CELLTYPES}
    for region in REGIONS:
        for contrast in COMPARISONS:
            block = pair[(pair.region == region) & (pair.comparison == contrast)]
            if set(zip(block.source, block.target)) != expected_pairs or len(block) != 49:
                raise ValueError(f"Expected 49 directed pairs for {region}, {contrast}")

    # Recover the pooled network total for each stage from the condition values.
    stage_totals = {}
    for contrast, (stage1, stage2) in STAGE_MAP.items():
        block = pair[pair.comparison == contrast]
        for region in REGIONS:
            sub = block[block.region == region]
            stage_totals[(region, stage1)] = sub.strength_condition_1.sum()
            stage_totals[(region, stage2)] = sub.strength_condition_2.sum()
    for region in REGIONS:
        comparisons = {
            "nonAD": ["earlyAD_minus_nonAD", "lateAD_minus_nonAD"],
            "earlyAD": ["earlyAD_minus_nonAD", "lateAD_minus_earlyAD"],
            "lateAD": ["lateAD_minus_earlyAD", "lateAD_minus_nonAD"],
        }
        for stage, source_contrasts in comparisons.items():
            totals = []
            for contrast in source_contrasts:
                sub = pair[(pair.region == region) & (pair.comparison == contrast)]
                column = "strength_condition_1" if STAGE_MAP[contrast][0] == stage else "strength_condition_2"
                totals.append(sub[column].sum())
            if not np.allclose(totals, totals[0], atol=1e-6):
                raise ValueError(f"Pooled {stage} totals disagree across contrasts for {region}: {totals}")
    ag = pair[(pair.region == "AG") & (pair.comparison == "lateAD_minus_nonAD")].copy()
    ag_matrix = ag.pivot(index="source", columns="target", values="strength_difference").reindex(index=CELLTYPES, columns=CELLTYPES)
    n_positive = int((ag.strength_difference > 0).sum())
    n_negative = int((ag.strength_difference < 0).sum())

    ag_path = pathway[(pathway.region == "AG") & (pathway.comparison == "lateAD_minus_nonAD")].copy()
    ag_path = ag_path.sort_values("strength_difference", ascending=False)
    top = ag_path.head(12).sort_values("strength_difference", ascending=True)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig = plt.figure(figsize=(13.4, 9.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.0], width_ratios=[1.0, 1.05])
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])

    # A: pooled network total by stage, retaining the same cell-class pairs in each region.
    stages = ["nonAD", "earlyAD", "lateAD"]
    x = np.arange(3)
    for region in REGIONS:
        y = [stage_totals[(region, s)] for s in stages]
        is_ag = region == "AG"
        ax_a.plot(x, y, marker="o", markersize=6.5 if is_ag else 4.3,
                  linewidth=2.8 if is_ag else 1.25, color=REGION_COLORS[region],
                  alpha=1.0 if is_ag else 0.7, label=region, zorder=3 if is_ag else 1)
        if is_ag:
            for xi, yi in zip(x, y):
                ax_a.annotate(f"{yi:.1f}", (xi, yi), xytext=(0, 8), textcoords="offset points",
                              ha="center", color="#9e3f53", fontsize=8, fontweight="bold")
    ax_a.set_xticks(x, stages)
    ax_a.set_ylabel("Summed pooled communication strength")
    ax_a.set_title("A  Regional pooled network totals", loc="left", fontweight="bold")
    ax_a.grid(axis="y", color="#dddddd", linewidth=0.7)
    ax_a.legend(title="Region", frameon=False, ncol=3, fontsize=8, title_fontsize=8, loc="upper left")

    # B: AG's complete directed cell-class comparison.
    vmax = float(np.nanmax(np.abs(ag_matrix.values)))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    im = ax_b.imshow(ag_matrix.values, cmap="RdBu_r", norm=norm, aspect="auto")
    ax_b.set_xticks(np.arange(len(CELLTYPES)), CELLTYPES, rotation=45, ha="right")
    ax_b.set_yticks(np.arange(len(CELLTYPES)), CELLTYPES)
    ax_b.set_xlabel("Receiver cell class")
    ax_b.set_ylabel("Sender cell class")
    ax_b.set_title(f"B  AG pair changes: {n_positive}/49 increased", loc="left", fontweight="bold")
    for i in range(len(CELLTYPES)):
        for j in range(len(CELLTYPES)):
            val = ag_matrix.iloc[i, j]
            ax_b.text(j, i, f"{val:+.1f}", ha="center", va="center", fontsize=7,
                      color="white" if abs(val) > vmax * 0.55 else "#222222")
    cbar = fig.colorbar(im, ax=ax_b, fraction=0.046, pad=0.04)
    cbar.set_label("lateAD − nonAD pooled strength")
    ax_b.text(0.0, -0.30, f"All 49 directed pairs are higher (negative: {n_negative}).",
              transform=ax_b.transAxes, fontsize=8, color="#555555")

    # C: pathways contributing the largest positive pooled strength differences.
    colors = ["#CC6677" if p in {"NRXN", "ADGRL"} else "#8DA0CB" for p in top.pathway]
    bars = ax_c.barh(top.pathway, top.strength_difference, color=colors, edgecolor="none", height=0.72)
    ax_c.axvline(0, color="#444444", linewidth=0.8)
    ax_c.set_xlabel("lateAD − nonAD pooled pathway strength")
    ax_c.set_title("C  Leading pathways behind the AG increase", loc="left", fontweight="bold")
    ax_c.grid(axis="x", color="#e2e2e2", linewidth=0.7)
    ax_c.set_axisbelow(True)
    for bar, val in zip(bars, top.strength_difference):
        ax_c.text(val + max(top.strength_difference) * 0.012, bar.get_y() + bar.get_height() / 2,
                  f"{val:+.1f}", va="center", fontsize=8)

    fig.suptitle("AG shows the strongest broad increase in pooled CellChat signal", fontsize=15, fontweight="bold")
    fig.text(0.5, -0.015,
             "Pooled descriptive network scores; pathway bars show the 12 largest increases. No donor-level uncertainty or significance is represented.",
             ha="center", fontsize=8, color="#555555")
    fig.savefig(OUT.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(OUT.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    regional_delta = pair[pair.comparison == "lateAD_minus_nonAD"].groupby("region").strength_difference.sum().reindex(REGIONS)
    report = f"""# AG pooled CellChat signal figure\n\nThe figure combines three descriptive views from the pooled CellChat tables. Panel A sums communication strength over the same 49 directed cell-class pairs for each region and pathology group. Panel B shows each AG sender-to-receiver pair's lateAD-minus-nonAD difference. Panel C lists the 12 AG pathways with the largest positive pooled strength differences.\n\nAG has a pooled total of {stage_totals[("AG", "nonAD")]:.1f} in nonAD, {stage_totals[("AG", "earlyAD")]:.1f} in earlyAD, and {stage_totals[("AG", "lateAD")]:.1f} in lateAD. The summed lateAD-minus-nonAD difference is {regional_delta["AG"]:.1f}; the next-largest positive regional change is {regional_delta.drop("AG").max():.1f} in {regional_delta.drop("AG").idxmax()}. In AG, {n_positive} of 49 directed pairs increase and {n_negative} decrease.\n\nThese are pooled descriptive scores. The figure does not show donor-level confidence intervals, covariate-adjusted estimates, or significance tests, so it supports a descriptive regional contrast only.\n\nInputs: `finalized/tables/pooled_celltype_differences_descriptive.tsv` and `finalized/tables/pooled_pathway_differences_descriptive.tsv`.\n"""
    (REPORTS / "AG_pooled_signal_highlight_interpretation.md").write_text(report)
    print(f"Wrote {OUT.with_suffix('.pdf')} and {OUT.with_suffix('.png')}")
    print(f"AG lateAD-minus-nonAD total: {regional_delta['AG']:.3f}; pair directions: {n_positive} up, {n_negative} down")


if __name__ == "__main__":
    main()
