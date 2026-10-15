#!/usr/bin/env python3
"""EXP011 variant: Theil-Sen robust regression and Kendall tau-b."""
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
from scipy.stats import kendalltau, theilslopes

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "experiments/EXP011_results/scripts"))
import analyze_spp1_tf_activity_correlations as base  # noqa: E402

OUT = Path(__file__).resolve().parents[1]
TABLES = OUT / "tables"
FIGURES = OUT / "figures"


def analyze_pair(frame: pd.DataFrame, diagnostics: dict) -> dict:
    activity = frame["TF_activity"].to_numpy(float)
    expression = frame["SPP1_expression"].to_numpy(float)
    tau_two = kendalltau(activity, expression, variant="b", alternative="two-sided")
    tau_greater = kendalltau(activity, expression, variant="b", alternative="greater")
    slope, intercept, slope_low, slope_high = theilslopes(
        expression, activity, alpha=0.95, method="separate"
    )
    return {
        **diagnostics,
        "n_used": len(frame),
        "kendall_tau_b": float(tau_two.statistic),
        "kendall_p_two_sided": float(tau_two.pvalue),
        "kendall_p_greater": float(tau_greater.pvalue),
        "theil_sen_slope_expression_per_activity": float(slope),
        "theil_sen_intercept": float(intercept),
        "theil_sen_slope_ci95_low": float(slope_low),
        "theil_sen_slope_ci95_high": float(slope_high),
    }


def plot_pair(frame: pd.DataFrame, result: dict, outbase: Path, figsize=(6.0, 4.7)) -> None:
    fig, ax = plt.subplots(figsize=figsize)
    for pathology in base.PATHOLOGY_ORDER:
        sub = frame[frame["pathology"] == pathology]
        ax.scatter(sub["TF_activity"], sub["SPP1_expression"],
                   color=base.PATHOLOGY_COLORS[pathology], alpha=0.85, s=48,
                   edgecolor="white", linewidth=0.45,
                   label=f"{pathology} (n={len(sub)})")
    activity = frame["TF_activity"].to_numpy(float)
    grid = np.linspace(float(activity.min()), float(activity.max()), 100)
    ax.plot(grid, result["theil_sen_intercept"]
            + result["theil_sen_slope_expression_per_activity"] * grid,
            color="#333333", linestyle="--", linewidth=1.6, label="Theil–Sen fit")
    ax.set_title(
        f"{result['context']}: {result['TF']} activity vs SPP1 expression\n"
        f"Kendall tau-b={result['kendall_tau_b']:.3f}; "
        f"two-sided p={result['kendall_p_two_sided']:.3g}; "
        f"positive one-sided p={result['kendall_p_greater']:.3g}; n={result['n_used']}",
        pad=12,
    )
    ax.set_xlabel(f"{result['TF']} TF activity")
    ax.set_ylabel("Cleaned SPP1 gene expression")
    ax.grid(color="#e5e5e5", linewidth=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(outbase.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    frames, results = [], []
    for spec in base.PAIRS:
        # These are the original EXP011 activity sources: Oli_AG is
        # cell-type-specific, while Oli_EC/SPI1 and Mic_Immune_AG/USF1 are
        # cross-celltype-restricted.
        frame, diagnostics = base.load_pair(spec)
        result = analyze_pair(frame, diagnostics)
        frames.append(frame)
        results.append(result)

    pd.concat(frames, ignore_index=True).to_csv(
        TABLES / "spp1_tf_activity_donor_pairs.tsv", sep="\t", index=False)
    summary = pd.DataFrame(results)
    summary.to_csv(TABLES / "spp1_tf_activity_robust_summary.tsv", sep="\t", index=False)

    for frame, result in zip(frames, results):
        plot_pair(frame, result,
                  FIGURES / f"{result['context']}_SPP1_{result['TF']}_robust_association")

    fig, axes = plt.subplots(1, 3, figsize=(17.5, 5.0))
    for ax, frame, result in zip(axes, frames, results):
        for pathology in base.PATHOLOGY_ORDER:
            sub = frame[frame["pathology"] == pathology]
            ax.scatter(sub["TF_activity"], sub["SPP1_expression"],
                       color=base.PATHOLOGY_COLORS[pathology], alpha=0.82, s=32,
                       edgecolor="white", linewidth=0.35, label=pathology)
        x = frame["TF_activity"].to_numpy(float)
        grid = np.linspace(float(x.min()), float(x.max()), 100)
        ax.plot(grid, result["theil_sen_intercept"]
                + result["theil_sen_slope_expression_per_activity"] * grid,
                color="#333333", linestyle="--", linewidth=1.4)
        ax.set_title(
            f"{result['context']} ({result['TF']})\n"
            f"tau-b={result['kendall_tau_b']:.2f}, "
            f"two-sided p={result['kendall_p_two_sided']:.2g}, "
            f"positive p={result['kendall_p_greater']:.2g}; n={result['n_used']}",
            fontsize=9,
        )
        ax.set_xlabel("TF activity")
        ax.set_ylabel("Cleaned SPP1 expression")
        ax.grid(color="#e5e5e5", linewidth=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 1.04))
    fig.suptitle("Robust SPP1–TF associations: Theil–Sen and Kendall tau-b", y=1.12, fontsize=15)
    fig.tight_layout()
    fig.savefig(FIGURES / "spp1_tf_activity_robust_all_pairs.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / "spp1_tf_activity_robust_all_pairs.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    report = [
        "# EXP011 variant: robust regression and Kendall tau-b",
        "",
        "## Analysis",
        "",
        "For each donor-level pair, this version reports Kendall's tau-b (which accounts for tied ranks) and fits a Theil–Sen robust regression with SPP1 expression as the outcome and TF activity as the predictor. No pathology covariate or interaction is included; pathology is only shown by plot color. No donors are removed and no multiple-testing correction is applied.",
        "",
        "Two-sided Kendall p-values are primary. Positive-direction one-sided p-values are included as exploratory results, requested after the associations had been inspected. Kendall p-values use SciPy's tie-adjusted asymptotic test. Theil–Sen slope confidence limits are 95% intervals. These analyses describe association, not direct TF regulation.",
        "",
        "Activity sources follow the original EXP011 plan: the cell-type-specific matrix for Oli_AG/ETV5 and cross-celltype-restricted matrices for Oli_EC/SPI1 and Mic_Immune_AG/USF1.",
        "",
        "## Results",
        "",
        "| Scope | Context | TF | n | Kendall tau-b | two-sided p (primary) | positive one-sided p (exploratory) | Theil–Sen slope | slope 95% CI |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for result in results:
        report.append(
            f"| {result['scope']} | {result['context']} | {result['TF']} | {result['n_used']} | "
            f"{result['kendall_tau_b']:.3f} | {result['kendall_p_two_sided']:.3g} | "
            f"{result['kendall_p_greater']:.3g} | "
            f"{result['theil_sen_slope_expression_per_activity']:.4g} | "
            f"({result['theil_sen_slope_ci95_low']:.4g}, {result['theil_sen_slope_ci95_high']:.4g}) |"
        )
    report += ["", "## Inputs", ""]
    for result in results:
        report.append(
            f"- {result['context']}: activity `{result['activity_input']}`; "
            f"expression `{result['expression_input']}`; metadata `{result['metadata_input']}`; "
            f"{result['n_used']} matched complete donors."
        )
    report += [
        "",
        "## Outputs",
        "",
        "- `tables/spp1_tf_activity_donor_pairs.tsv`: matched donor observations.",
        "- `tables/spp1_tf_activity_robust_summary.tsv`: Kendall tau-b p-values and Theil–Sen estimates.",
        "- `figures/*_robust_association.{pdf,png}`: one plot per pair.",
        "- `figures/spp1_tf_activity_robust_all_pairs.{pdf,png}`: combined figure.",
        "- `scripts/analyze_spp1_tf_activity_robust.py`: reproducible analysis script.",
        "",
        "## Software",
        "",
        f"Python {platform.python_version()}; pandas {pd.__version__}; NumPy {np.__version__}; SciPy {scipy.__version__}; Matplotlib {matplotlib.__version__}.",
    ]
    (OUT / "EXP011_robust_regression_kendall_tau_b.md").write_text("\n".join(report) + "\n")
    (OUT / "session_info.txt").write_text(
        f"Python {platform.python_version()}\npandas {pd.__version__}\nNumPy {np.__version__}\n"
        f"SciPy {scipy.__version__}\nMatplotlib {matplotlib.__version__}\n"
    )
    print(summary[["context", "TF", "n_used", "kendall_tau_b", "kendall_p_two_sided",
                   "kendall_p_greater", "theil_sen_slope_expression_per_activity"]].to_string(index=False))
    print(f"\nWrote robust EXP011 outputs under {OUT}")


if __name__ == "__main__":
    main()
