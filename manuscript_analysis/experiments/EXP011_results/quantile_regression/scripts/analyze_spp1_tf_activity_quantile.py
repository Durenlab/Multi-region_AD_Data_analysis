#!/usr/bin/env python3
"""EXP011 quantile regression variant with pathology-specific TF slopes."""
from __future__ import annotations

import platform
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import statsmodels.api as sm
from scipy.stats import chi2, kendalltau

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "experiments/EXP011_results/scripts"))
import analyze_spp1_tf_activity_correlations as base  # noqa: E402

OUT = Path(__file__).resolve().parents[1]
TABLES = OUT / "tables"
FIGURES = OUT / "figures"
QUANTILES = (0.50, 0.75)
PATHOLOGY_CONTRASTS = {
    "nonAD": (np.array([0, 1, 0, 0, 0, 0.]), np.array([1, 0, 0, 0, 0, 0.])),
    "earlyAD": (np.array([0, 1, 0, 0, 1, 0.]), np.array([1, 0, 1, 0, 0, 0.])),
    "lateAD": (np.array([0, 1, 0, 0, 0, 1.]), np.array([1, 0, 0, 1, 0, 0.])),
}


def fit_quantile(frame: pd.DataFrame, quantile: float) -> tuple[dict, list[dict], object]:
    y = frame["SPP1_expression"].to_numpy(float)
    activity = frame["TF_activity"].to_numpy(float)
    groups = frame["pathology"].astype(str).to_numpy()
    mean_activity = float(activity.mean())
    sd_activity = float(activity.std(ddof=1))
    group_ranges = {
        group: (float(activity[groups == group].min()), float(activity[groups == group].max()))
        for group in base.PATHOLOGY_ORDER
    }
    common_min = max(bounds[0] for bounds in group_ranges.values())
    common_max = min(bounds[1] for bounds in group_ranges.values())
    if common_min > common_max:
        raise ValueError("No TF activity overlap across pathology groups for a common-support comparison")
    common_reference = (common_min + common_max) / 2
    x = (activity - mean_activity) / sd_activity
    early = (groups == "earlyAD").astype(float)
    late = (groups == "lateAD").astype(float)
    design = np.column_stack([np.ones(len(y)), x, early, late, x * early, x * late])
    if np.linalg.matrix_rank(design) < design.shape[1]:
        raise ValueError("Quantile-regression design matrix is rank deficient")
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        fit = sm.QuantReg(y, design).fit(q=quantile, max_iter=10_000, p_tol=1e-8)
    covariance = np.asarray(fit.cov_params())
    beta = np.asarray(fit.params)
    interaction = beta[[4, 5]]
    interaction_cov = covariance[np.ix_([4, 5], [4, 5])]
    wald = float(interaction @ np.linalg.solve(interaction_cov, interaction))
    model_result = {
        "quantile": quantile,
        "n": int(len(frame)),
        "activity_mean": mean_activity,
        "activity_sd": sd_activity,
        "common_activity_min": common_min,
        "common_activity_max": common_max,
        "common_activity_reference": common_reference,
        "interaction_wald_chi2": wald,
        "interaction_df": 2,
        "interaction_p_asymptotic": float(chi2.sf(wald, 2)),
        "iterations": int(fit.iterations),
        "warnings": "; ".join(str(item.message) for item in captured),
    }
    group_rows = []
    z95 = 1.959963984540054
    for group, (slope_contrast, location_contrast) in PATHOLOGY_CONTRASTS.items():
        group_n = int(np.sum(groups == group))
        slope = float(slope_contrast @ beta)
        slope_se = float(np.sqrt(slope_contrast @ covariance @ slope_contrast))
        prediction_contrast = location_contrast + ((common_reference - mean_activity) / sd_activity) * slope_contrast
        location = float(prediction_contrast @ beta)
        location_se = float(np.sqrt(prediction_contrast @ covariance @ prediction_contrast))
        group_rows.append({
            "pathology": group,
            "n_group": group_n,
            "common_activity_reference": common_reference,
            "slope_SPP1_per_1SD_TF_activity": slope,
            "slope_se_asymptotic": slope_se,
            "slope_ci95_low": slope - z95 * slope_se,
            "slope_ci95_high": slope + z95 * slope_se,
            "conditional_quantile_SPP1_at_common_activity": location,
            "location_se_asymptotic": location_se,
            "location_ci95_low": location - z95 * location_se,
            "location_ci95_high": location + z95 * location_se,
        })
    return model_result, group_rows, fit


def plot_quantile(frames: list[pd.DataFrame], results: list[dict], group_results: list[dict],
                  quantile: float) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(17.5, 5.0))
    for ax, frame, result, estimates in zip(axes, frames, results, group_results):
        for row in estimates:
            group = row["pathology"]
            sub = frame[frame["pathology"].astype(str) == group]
            ax.scatter(sub["TF_activity"], sub["SPP1_expression"],
                       color=base.PATHOLOGY_COLORS[group], alpha=0.8, s=33,
                       edgecolor="white", linewidth=0.35, label=f"{group} (n={len(sub)})")
            grid = np.linspace(float(sub["TF_activity"].min()),
                               float(sub["TF_activity"].max()), 100)
            center = row["common_activity_reference"]
            scale = result["activity_sd"]
            predicted = (row["conditional_quantile_SPP1_at_common_activity"]
                         + row["slope_SPP1_per_1SD_TF_activity"] * (grid - center) / scale)
            ax.plot(grid, predicted, color=base.PATHOLOGY_COLORS[group], linewidth=1.6)
        ax.set_title(f"{result['context']} ({result['TF']}); tau={quantile:.2f}\n"
                     f"interaction p={result['interaction_p_asymptotic']:.3g}; n={result['n']}")
        ax.set_xlabel("TF activity")
        ax.set_ylabel("Cleaned SPP1 expression")
        ax.grid(color="#e5e5e5", linewidth=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 1.04))
    fig.suptitle(f"Pathology-specific SPP1 quantile regression (tau={quantile:.2f})",
                 y=1.12, fontsize=15)
    fig.tight_layout()
    suffix = f"tau_{int(quantile * 100):02d}"
    fig.savefig(FIGURES / f"spp1_tf_activity_quantile_{suffix}.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / f"spp1_tf_activity_quantile_{suffix}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    frames, diagnostics = [], []
    for spec in base.PAIRS:
        frame, check = base.load_pair(spec)
        frames.append(frame)
        diagnostics.append(check)
    pd.concat(frames, ignore_index=True).to_csv(
        TABLES / "spp1_tf_activity_donor_pairs.tsv", sep="\t", index=False)

    kendall_rows = []
    for frame, check in zip(frames, diagnostics):
        x = frame["TF_activity"].to_numpy(float)
        y = frame["SPP1_expression"].to_numpy(float)
        tau2 = kendalltau(x, y, variant="b", alternative="two-sided")
        tau1 = kendalltau(x, y, variant="b", alternative="greater")
        kendall_rows.append({**check, "n_used": len(frame), "kendall_tau_b": float(tau2.statistic),
                             "kendall_p_two_sided": float(tau2.pvalue),
                             "kendall_p_greater": float(tau1.pvalue)})
    pd.DataFrame(kendall_rows).to_csv(
        TABLES / "spp1_tf_activity_kendall_tau_b.tsv", sep="\t", index=False)

    model_rows, slope_rows = [], []
    for quantile in QUANTILES:
        quantile_results, quantile_slopes, slopes_by_context = [], [], []
        for frame, check in zip(frames, diagnostics):
            model, slopes, _ = fit_quantile(frame, quantile)
            quantile_results.append({"scope": check["scope"], "context": check["context"],
                                     "TF": check["TF"], **model})
            slopes_by_context.append(slopes)
            quantile_slopes.extend({"scope": check["scope"], "context": check["context"],
                                    "TF": check["TF"], "quantile": quantile, **row}
                                   for row in slopes)
        model_rows.extend(quantile_results)
        slope_rows.extend(quantile_slopes)
        plot_quantile(frames, quantile_results, slopes_by_context, quantile)
    pd.DataFrame(model_rows).to_csv(
        TABLES / "spp1_tf_activity_quantile_model_summary.tsv", sep="\t", index=False)
    pd.DataFrame(slope_rows).to_csv(
        TABLES / "spp1_tf_activity_quantile_group_estimates.tsv", sep="\t", index=False)

    report = [
        "# EXP011 variant: pathology-specific quantile regression",
        "",
        "## Analysis",
        "",
        "This version uses the original activity sources (cell-type-specific for Oli_AG/ETV5 and cross-celltype-restricted for Oli_EC/SPI1 and Mic_Immune_AG/USF1). It fits quantile regressions at tau=0.50 (conditional median) and tau=0.75 (upper quartile): cleaned SPP1 expression ~ standardized TF activity × categorical pathology, with nonAD as reference. It estimates pathology-specific activity slopes and conditional quantiles. For group-level quantile comparisons, estimates are evaluated at the midpoint of the TF activity range shared across pathology groups; plotted lines are limited to each group's observed activity range. All complete donors are retained.",
        "",
        "The interaction test is an asymptotic two-degree-of-freedom Wald chi-square test using the QuantReg covariance estimate. Quantile-regression standard errors and p-values are approximate, especially for the small lateAD groups. Treat results as exploratory. The Kendall table provides an unadjusted donor-level tau-b comparison; two-sided p is primary and positive one-sided p is exploratory. No multiple-testing correction is applied.",
        "",
        "## Kendall tau-b",
        "",
        "| Context | TF | n | tau-b | two-sided p | positive one-sided p |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for row in kendall_rows:
        report.append(f"| {row['context']} | {row['TF']} | {row['n_used']} | {row['kendall_tau_b']:.3f} | {row['kendall_p_two_sided']:.3g} | {row['kendall_p_greater']:.3g} |")
    report += ["", "## Quantile-regression interaction tests", "",
               "| Context | TF | tau | N | Interaction Wald chi-square (df=2) | asymptotic p |",
               "|---|---|---:|---:|---:|---:|"]
    for row in model_rows:
        report.append(f"| {row['context']} | {row['TF']} | {row['quantile']:.2f} | {row['n']} | {row['interaction_wald_chi2']:.3f} | {row['interaction_p_asymptotic']:.3g} |")
    report += ["", "Pathology-specific estimates (slope per 1 SD TF activity; conditional SPP1 quantile at the midpoint of the TF activity range shared by all pathology groups):", "",
               "| Context | tau | Pathology | n | Common-support TF reference | Slope | slope 95% CI | conditional SPP1 quantile | quantile 95% CI |",
               "|---|---:|---|---:|---:|---:|---:|---:|---:|"]
    for row in slope_rows:
        report.append(f"| {row['context']} | {row['quantile']:.2f} | {row['pathology']} | {row['n_group']} | {row['common_activity_reference']:.4g} | {row['slope_SPP1_per_1SD_TF_activity']:.4g} | ({row['slope_ci95_low']:.4g}, {row['slope_ci95_high']:.4g}) | {row['conditional_quantile_SPP1_at_common_activity']:.4g} | ({row['location_ci95_low']:.4g}, {row['location_ci95_high']:.4g}) |")
    report += ["", "LateAD has only seven donors per context, so its conditional quantile estimates and intervals are particularly uncertain. Quantile regression reduces sensitivity to extreme outcome residuals; it does not eliminate high-leverage predictor observations. The model estimates association, not TF causation.",
               "", "## Inputs", ""]
    for row in diagnostics:
        report.append(f"- {row['context']}: activity `{row['activity_input']}`; expression `{row['expression_input']}`; metadata `{row['metadata_input']}`; {row['n_complete_cases']} complete donors.")
    report += ["", "## Outputs", "",
               "- `tables/spp1_tf_activity_donor_pairs.tsv`: matched donors.",
               "- `tables/spp1_tf_activity_kendall_tau_b.tsv`: unadjusted Kendall tau-b statistics.",
               "- `tables/spp1_tf_activity_quantile_model_summary.tsv`: joint interaction tests.",
               "- `tables/spp1_tf_activity_quantile_group_estimates.tsv`: group-specific slopes and conditional quantiles.",
               "- `figures/spp1_tf_activity_quantile_tau_50.{pdf,png}` and `...tau_75...`: pathology-colored fitted quantile lines.",
               "- `scripts/analyze_spp1_tf_activity_quantile.py`: reproduction script.",
               "", "## Software", "",
               f"Python {platform.python_version()}; pandas {pd.__version__}; NumPy {np.__version__}; SciPy {scipy.__version__}; Statsmodels {sm.__version__}; Matplotlib {matplotlib.__version__}."]
    (OUT / "EXP011_quantile_regression.md").write_text("\n".join(report) + "\n")
    (OUT / "session_info.txt").write_text(
        f"Python {platform.python_version()}\npandas {pd.__version__}\nNumPy {np.__version__}\n"
        f"SciPy {scipy.__version__}\nStatsmodels {sm.__version__}\nMatplotlib {matplotlib.__version__}\n")
    print(pd.DataFrame(model_rows)[["context", "TF", "quantile", "interaction_wald_chi2", "interaction_p_asymptotic"]].to_string(index=False))
    print(f"\nWrote quantile-regression outputs under {OUT}")


if __name__ == "__main__":
    main()
