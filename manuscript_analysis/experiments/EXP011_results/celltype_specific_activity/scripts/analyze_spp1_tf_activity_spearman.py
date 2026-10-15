#!/usr/bin/env python3
"""EXP011 variant using cell-type-specific TF activity for every pair.

This version reports Spearman's rank correlation only. It keeps the original
EXP011 expression and metadata sources and substitutes cell-type-specific
activity matrices for all three TFs.
"""
from __future__ import annotations

import platform
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy

ROOT = Path(__file__).resolve().parents[4]
BASE_SCRIPT_DIR = ROOT / "experiments/EXP011_results/scripts"
sys.path.insert(0, str(BASE_SCRIPT_DIR))
import analyze_spp1_tf_activity_correlations as base  # noqa: E402

OUT = Path(__file__).resolve().parents[1]
TABLES = OUT / "tables"
FIGURES = OUT / "figures"
CELLTYPE_ACTIVITY_DIR = base.PALMETTO / "AD_REGA/part7_TF/01TF_activity"

PAIRS = []
for original in base.PAIRS:
    spec = original.copy()
    spec["activity"] = CELLTYPE_ACTIVITY_DIR / f"{spec['context']}_TF_activity_pcc_cutoff.csv"
    PAIRS.append(spec)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    frames, results = [], []
    for spec in PAIRS:
        frame, diagnostics = base.load_pair(spec, require_estimable=False)
        if frame["SPP1_expression"].nunique() < 2:
            status = "not estimable: constant SPP1 expression"
        elif frame["TF_activity"].nunique() < 2:
            status = "not estimable: constant TF activity"
        else:
            status = "estimable"
        if status == "estimable":
            result = base.analyze(frame, diagnostics)
        else:
            result = {
                **diagnostics,
                "n_used": int(len(frame)),
                "spearman_rho": np.nan,
                "permutation_p_two_sided": np.nan,
                "permutation_p_greater": np.nan,
                "permutations": base.PERMUTATIONS,
                "permutation_seed": base.PERMUTATION_SEED,
            }
        result["status"] = status
        frames.append(frame)
        results.append(result)

    joined = pd.concat(frames, ignore_index=True)
    summary = pd.DataFrame(results)
    joined.to_csv(TABLES / "spp1_tf_activity_donor_pairs.tsv", sep="\t", index=False)
    summary.to_csv(TABLES / "spp1_tf_activity_spearman_summary.tsv", sep="\t", index=False)

    for frame, result in zip(frames, results):
        if result["status"] == "estimable":
            base.plot_pair(
                frame, result,
                FIGURES / f"{result['context']}_SPP1_{result['TF']}_activity_spearman",
                figsize=(6.0, 4.7),
            )
        else:
            fig, ax = plt.subplots(figsize=(6.0, 4.7))
            for pathology in base.PATHOLOGY_ORDER:
                sub = frame[frame["pathology"] == pathology]
                ax.scatter(sub["SPP1_expression"], sub["TF_activity"],
                           color=base.PATHOLOGY_COLORS[pathology], alpha=0.85, s=48,
                           edgecolor="white", linewidth=0.45, label=f"{pathology} (n={len(sub)})")
            ax.set_title(f"{result['context']}: {result['TF']} activity is constant; Spearman rho undefined")
            ax.set_xlabel("Cleaned SPP1 gene expression")
            ax.set_ylabel("TF activity (cell-type-specific universe)")
            ax.legend(frameon=False, fontsize=8)
            fig.tight_layout()
            outbase = FIGURES / f"{result['context']}_SPP1_{result['TF']}_activity_spearman"
            fig.savefig(outbase.with_suffix(".pdf"), bbox_inches="tight")
            fig.savefig(outbase.with_suffix(".png"), dpi=300, bbox_inches="tight")
            plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(17.5, 5.0))
    for ax, frame, result in zip(axes, frames, results):
        for pathology in base.PATHOLOGY_ORDER:
            sub = frame[frame["pathology"] == pathology]
            ax.scatter(sub["SPP1_expression"], sub["TF_activity"],
                       color=base.PATHOLOGY_COLORS[pathology], alpha=0.82, s=32,
                       edgecolor="white", linewidth=0.35, label=pathology)
        x = frame["SPP1_expression"].to_numpy(float)
        y = frame["TF_activity"].to_numpy(float)
        if result["status"] == "estimable":
            slope, intercept, _, _ = base.theilslopes(y, x)
            grid = np.linspace(float(np.min(x)), float(np.max(x)), 100)
            ax.plot(grid, intercept + slope * grid, color="#333333", linestyle="--", linewidth=1.4)
            title = (f"Spearman rho={result['spearman_rho']:.2f}, "
                     f"two-sided p={result['permutation_p_two_sided']:.2g}, "
                     f"positive one-sided p={result['permutation_p_greater']:.2g}; n={result['n_used']}")
        else:
            title = f"rho undefined: constant TF activity; n={result['n_used']}"
        ax.set_title(f"{result['context']} ({result['TF']})\n{title}", fontsize=9)
        ax.set_xlabel("Cleaned SPP1 expression")
        ax.set_ylabel("TF activity (cell-type-specific universe)")
        ax.grid(color="#e5e5e5", linewidth=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 1.04))
    fig.suptitle("SPP1 expression vs cell-type-specific TF activity", y=1.12, fontsize=15)
    fig.tight_layout()
    fig.savefig(FIGURES / "spp1_tf_activity_spearman_all_pairs.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / "spp1_tf_activity_spearman_all_pairs.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    report = [
        "# EXP011 variant: cell-type-specific TF activity",
        "",
        "## Analysis",
        "",
        "This parallel version uses the cell-type-specific TF activity matrix from `AD_REGA/part7_TF/01TF_activity` for every pair. The expression and metadata inputs remain those specified in the original EXP011 plan. Only Spearman's rank correlation is performed; no pathology regression or interaction model is fit. If either variable is constant, the correlation is reported as not estimable.",
        "",
        "Two-sided permutation p-values are primary. Positive-direction one-sided permutation p-values are included as exploratory values because that direction was requested after observing the associations. Each p-value uses 100,000 permutations of donor ranks with tied ranks preserved and a fixed seed. No pathology adjustment or multiple-testing correction is applied. All complete matched donors are retained.",
        "",
        "## Results",
        "",
        "| Scope | Context | TF | n | Status | Spearman rho | two-sided p (primary) | positive one-sided p (exploratory) |",
        "|---|---|---:|---|---:|---:|---:|---:|",
    ]
    for result in results:
        report.append(
            f"| {result['scope']} | {result['context']} | {result['TF']} | {result['n_used']} | "
            f"{result['status']} | "
            f"{result['spearman_rho']:.3f} | {result['permutation_p_two_sided']:.3g} | "
            f"{result['permutation_p_greater']:.3g} |" if result["status"] == "estimable" else
            f"| {result['scope']} | {result['context']} | {result['TF']} | {result['n_used']} | "
            f"{result['status']} | NA | NA | NA |"
        )
    report += ["", "## Inputs", ""]
    for result in results:
        report.append(
            f"- {result['context']}: activity `{result['activity_input']}`; "
            f"expression `{result['expression_input']}`; metadata `{result['metadata_input']}`. "
            f"Matched complete donors: {result['n_used']} (shared IDs before complete-case filtering: {result['n_shared_donors']})."
        )
    report += [
        "",
        "## Outputs",
        "",
        "- `tables/spp1_tf_activity_donor_pairs.tsv`: matched observations.",
        "- `tables/spp1_tf_activity_spearman_summary.tsv`: Spearman rho, both permutation p-values, and estimability status.",
        "- `figures/*_activity_spearman.{pdf,png}`: one plot per pair.",
        "- `figures/spp1_tf_activity_spearman_all_pairs.{pdf,png}`: combined figure.",
        "- `scripts/analyze_spp1_tf_activity_spearman.py`: reproducible analysis script.",
        "",
        "## Software",
        "",
        f"Python {platform.python_version()}; pandas {pd.__version__}; NumPy {np.__version__}; SciPy {scipy.__version__}; Matplotlib {matplotlib.__version__}.",
    ]
    (OUT / "EXP011_celltype_specific_activity.md").write_text("\n".join(report) + "\n")
    (OUT / "session_info.txt").write_text(
        f"Python {platform.python_version()}\npandas {pd.__version__}\nNumPy {np.__version__}\n"
        f"SciPy {scipy.__version__}\nMatplotlib {matplotlib.__version__}\n"
    )
    print(summary[["context", "TF", "n_used", "spearman_rho",
                   "permutation_p_two_sided", "permutation_p_greater"]].to_string(index=False))
    print(f"\nWrote alternate EXP011 outputs under {OUT}")


if __name__ == "__main__":
    main()
