#!/usr/bin/env python3
"""Plot descriptive pooled CellChat strength changes by region."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "experiments/EXP002_results"
TABLE = BASE / "finalized/tables/pooled_celltype_differences_descriptive.tsv"
OUT = BASE / "finalized/figures"
REGIONS = ["AG", "EC", "HC", "MTC", "PFC", "TH"]
CELLTYPES = ["Exc", "Inh", "Ast", "Oli", "OPC", "Mic_Immune", "Vasc"]
CONTRASTS = ["earlyAD_minus_nonAD", "lateAD_minus_earlyAD", "lateAD_minus_nonAD"]
CONTRAST_LABELS = ["earlyAD − nonAD", "lateAD − earlyAD", "lateAD − nonAD"]
CONTRAST_COLORS = ["#4477AA", "#EEAA33", "#BB5566"]

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "ps.fonttype": 42})


def load_summary():
    d = pd.read_csv(TABLE, sep="\t")
    expected = len(REGIONS) * len(CONTRASTS) * len(CELLTYPES) ** 2
    if len(d) != expected:
        raise ValueError(f"Expected {expected} pooled directed-pair rows, found {len(d)}")
    if not np.allclose(d.strength_difference,
                       d.strength_condition_2 - d.strength_condition_1,
                       atol=1e-8, equal_nan=True):
        raise ValueError("Pooled strength_difference does not match condition_2 - condition_1")
    return d


def make_figure(d):
    late = d[d.comparison == "lateAD_minus_nonAD"]
    piv = late.pivot(index=["source", "target"], columns="region", values="strength_difference")
    pair_index = pd.MultiIndex.from_product([CELLTYPES, CELLTYPES], names=["source", "target"])
    piv = piv.reindex(index=pair_index, columns=REGIONS)
    limit = float(np.nanmax(np.abs(piv.to_numpy(float))))
    norm = TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)

    fig = plt.figure(figsize=(14.4, 11.4))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.45, .85], hspace=.40,
                             left=.075, right=.94, top=.93, bottom=.16)
    grid = outer[0].subgridspec(2, 3, wspace=.30, hspace=.43)
    cmap = plt.get_cmap("RdBu_r")
    ims = []
    for i, region in enumerate(REGIONS):
        ax = fig.add_subplot(grid[i // 3, i % 3])
        mat = piv[region].to_numpy(float).reshape(len(CELLTYPES), len(CELLTYPES))
        im = ax.imshow(mat, cmap=cmap, norm=norm, aspect="equal")
        ims.append(im)
        ax.set_title(region, loc="left", fontsize=10, weight="bold", pad=3)
        ax.set_xticks(range(len(CELLTYPES)), CELLTYPES, rotation=48, ha="right", fontsize=6.5)
        ax.set_yticks(range(len(CELLTYPES)), CELLTYPES, fontsize=6.5)
        ax.set_xlabel("Receiver", fontsize=7)
        ax.set_ylabel("Sender", fontsize=7)
        ax.set_xticks(np.arange(-.5, len(CELLTYPES), 1), minor=True)
        ax.set_yticks(np.arange(-.5, len(CELLTYPES), 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=.65)
        ax.tick_params(which="minor", bottom=False, left=False)
    cb = fig.colorbar(ims[0], ax=[fig.axes[j] for j in range(6)], fraction=.025, pad=.025)
    cb.set_label("Pooled communication strength difference\n(lateAD − nonAD)", fontsize=8)

    axb = fig.add_subplot(outer[1])
    totals = d.groupby(["region", "comparison"], observed=True)["strength_difference"].sum()
    x = np.arange(len(REGIONS))
    width = .22
    for j, (comparison, label, color) in enumerate(zip(CONTRASTS, CONTRAST_LABELS, CONTRAST_COLORS)):
        vals = np.array([totals.loc[(r, comparison)] for r in REGIONS])
        xpos = x + (j - 1) * width
        bars = axb.bar(xpos, vals, width=width, color=color, label=label, edgecolor="white", linewidth=.4)
        for b, val in zip(bars, vals):
            y = val + (1.0 if val >= 0 else -1.0)
            axb.text(b.get_x() + b.get_width()/2, y, f"{val:+.1f}", ha="center",
                     va="bottom" if val >= 0 else "top", fontsize=7, color="#333333")
    axb.axhline(0, color="#333333", linewidth=.8)
    axb.set_xticks(x, REGIONS)
    axb.set_ylabel("Sum of pooled edge-strength differences\nacross 49 directed pairs")
    axb.set_title("B  Overall pooled network strength change", loc="left", fontsize=11, weight="bold")
    axb.grid(axis="y", color="#dddddd", linewidth=.55)
    axb.set_axisbelow(True)
    axb.legend(ncol=3, frameon=False, loc="upper right", fontsize=8)

    fig.suptitle("Regional CellChat signal-strength changes in pooled networks", fontsize=15, weight="bold", y=.985)
    fig.text(.5, .045,
             "A: Pair-level strength difference (lateAD − nonAD); red = higher pooled strength, blue = lower. "
             "B: Pairwise contrast sums across all 49 directed sender–receiver pairs.\n"
             "Descriptive pooled summaries only: no donor-level confidence intervals, covariate adjustment, or significance tests.",
             ha="center", va="center", fontsize=8)
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "pooled_strength_change_by_region_main_figure.pdf", bbox_inches="tight")
    fig.savefig(OUT / "pooled_strength_change_by_region_main_figure.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return totals


def main():
    d = load_summary()
    totals = make_figure(d)
    print("Saved pooled regional strength-change heatmap and bar plot.")
    print(totals.unstack().reindex(index=REGIONS, columns=CONTRASTS).round(4).to_string())
    print("Source:", TABLE)


if __name__ == "__main__":
    main()
