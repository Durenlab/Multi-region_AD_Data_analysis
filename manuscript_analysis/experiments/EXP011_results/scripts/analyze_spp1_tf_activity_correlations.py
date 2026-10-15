#!/usr/bin/env python3
"""Correlate donor-level SPP1 expression and selected TF activity (EXP011).

Reads only existing derived donor matrices and metadata. Writes joined data,
correlation results, figures, and a compact report under EXP011_results.
"""
from __future__ import annotations

from pathlib import Path
import platform

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.stats import f as f_distribution, rankdata, t as t_distribution, theilslopes


ROOT = Path(__file__).resolve().parents[3]
PALMETTO = Path("/project/zduren/durenlab/naqing_palmetto")
OUT = ROOT / "experiments/EXP011_results"
TABLES = OUT / "tables"
FIGURES = OUT / "figures"
PATHOLOGY_ORDER = ["nonAD", "earlyAD", "lateAD"]
PATHOLOGY_COLORS = {"nonAD": "#4477AA", "earlyAD": "#EE9944", "lateAD": "#CC6677"}
PAIRS = [
    {
        "scope": "cell-type-specific", "context": "Oli_AG", "tf": "ETV5",
        "expression": PALMETTO / "AD_REGA/part1_prepare_REGA/04cleaned_gex/cleaned_Oli_AG.tsv",
        "activity": PALMETTO / "AD_REGA/part7_TF/01TF_activity/Oli_AG_TF_activity_pcc_cutoff.csv",
        "metadata": PALMETTO / "AD_REGA/part1_prepare_REGA/03metadata_by_group/Oli_AG_metadata.tsv",
    },
    {
        "scope": "cross-celltype", "context": "Oli_EC", "tf": "SPI1",
        "expression": PALMETTO / "AD_REGA_region_wise/part1_prepare_REGA/02cleaned_region_consistent_gex/cleaned_Oli_EC_region_consistent.tsv",
        "activity": PALMETTO / "AD_REGA_region_wise/part4_TF/03TF_activity/Oli_EC/TF_activity_pcc_pos_only/Oli_EC_TF_activity_pcc_cutoff.csv",
        "metadata": PALMETTO / "AD_REGA_region_wise/part1_prepare_REGA/01region_consistent_metadata/Oli_EC_region_consistent_metadata.tsv",
    },
    {
        "scope": "cross-celltype", "context": "Mic_Immune_AG", "tf": "USF1",
        "expression": PALMETTO / "AD_REGA_region_wise/part1_prepare_REGA/02cleaned_region_consistent_gex/cleaned_Mic_Immune_AG_region_consistent.tsv",
        "activity": PALMETTO / "AD_REGA_region_wise/part4_TF/03TF_activity/Mic_Immune_AG/TF_activity_pcc_pos_only/Mic_Immune_AG_TF_activity_pcc_cutoff.csv",
        "metadata": PALMETTO / "AD_REGA_region_wise/part1_prepare_REGA/01region_consistent_metadata/Mic_Immune_AG_region_consistent_metadata.tsv",
    },
]
PERMUTATIONS = 100_000
PERMUTATION_SEED = 1108


def load_pair(spec: dict, require_estimable: bool = True) -> tuple[pd.DataFrame, dict]:
    for key in ["expression", "activity", "metadata"]:
        if not spec[key].exists():
            raise FileNotFoundError(f"Missing {key} input: {spec[key]}")
    expression = pd.read_csv(spec["expression"], sep="\t", index_col=0)
    activity = pd.read_csv(spec["activity"], index_col=0)
    expression.index = expression.index.astype(str).str.strip().str.upper()
    activity.index = activity.index.astype(str).str.strip().str.upper()
    expression.columns = expression.columns.astype(str).str.strip()
    activity.columns = activity.columns.astype(str).str.strip()
    if expression.columns.duplicated().any() or activity.columns.duplicated().any():
        raise ValueError(f"Duplicate donor IDs in one of the matrices for {spec['context']}")
    if (expression.index == "SPP1").sum() != 1:
        raise ValueError(f"Expected exactly one SPP1 row in {spec['expression']}")
    if (activity.index == spec["tf"].upper()).sum() != 1:
        raise ValueError(f"Expected exactly one {spec['tf']} row in {spec['activity']}")

    metadata = pd.read_csv(spec["metadata"], sep="\t", dtype={"ROSMAP_IndividualID": str})
    needed = {"ROSMAP_IndividualID", "Pathology"}
    if not needed.issubset(metadata.columns):
        raise ValueError(f"Metadata lacks required columns: {sorted(needed - set(metadata.columns))}")
    metadata["ROSMAP_IndividualID"] = metadata["ROSMAP_IndividualID"].astype(str).str.strip()
    metadata["Pathology"] = metadata["Pathology"].astype(str).str.strip()
    if metadata["ROSMAP_IndividualID"].duplicated().any():
        raise ValueError(f"Duplicate metadata donor IDs for {spec['context']}")

    common = sorted(set(expression.columns) & set(activity.columns) &
                    set(metadata["ROSMAP_IndividualID"]))
    if not common:
        raise ValueError(f"No donor ID overlap for {spec['context']}")
    metadata = metadata.set_index("ROSMAP_IndividualID")
    frame = pd.DataFrame({
        "scope": spec["scope"],
        "context": spec["context"],
        "TF": spec["tf"],
        "donor_id": common,
        "pathology": metadata.loc[common, "Pathology"].to_numpy(),
        "SPP1_expression": pd.to_numeric(expression.loc["SPP1", common], errors="coerce").to_numpy(),
        "TF_activity": pd.to_numeric(activity.loc[spec["tf"].upper(), common], errors="coerce").to_numpy(),
    })
    frame["pathology"] = pd.Categorical(frame["pathology"], categories=PATHOLOGY_ORDER, ordered=True)
    n_before = len(frame)
    frame = frame[
        frame["pathology"].notna()
        & np.isfinite(frame["SPP1_expression"])
        & np.isfinite(frame["TF_activity"])
    ].copy()
    if frame["donor_id"].duplicated().any():
        raise ValueError(f"Duplicate matched donors in joined data for {spec['context']}")
    if len(frame) < 3 or (require_estimable and
                          (frame["SPP1_expression"].nunique() < 2
                           or frame["TF_activity"].nunique() < 2)):
        raise ValueError(f"Correlation is unestimable for {spec['context']}: insufficient N or zero variance")

    diagnostics = {
        "scope": spec["scope"],
        "context": spec["context"],
        "TF": spec["tf"],
        "n_expression_matrix_donors": int(expression.shape[1]),
        "n_activity_matrix_donors": int(activity.shape[1]),
        "n_metadata_donors": int(len(metadata)),
        "n_shared_donors": int(len(common)),
        "n_complete_cases": int(len(frame)),
        "n_excluded_after_matching": int(n_before - len(frame)),
        "expression_input": str(spec["expression"]),
        "activity_input": str(spec["activity"]),
        "metadata_input": str(spec["metadata"]),
    }
    return frame, diagnostics


def spearman_rho(x: np.ndarray, y: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    """Compute Spearman's rho from average ranks, including ties."""
    xr, yr = rankdata(x, method="average"), rankdata(y, method="average")
    xr = xr - xr.mean()
    yr = yr - yr.mean()
    denominator = np.sqrt(np.sum(xr ** 2) * np.sum(yr ** 2))
    if denominator == 0:
        raise ValueError("Spearman correlation is undefined for a constant rank vector")
    return float(np.sum(xr * yr) / denominator), xr, yr


def permutation_pvalues(xr: np.ndarray, yr: np.ndarray, observed: float) -> tuple[float, float]:
    """Permutation p-values for Spearman rho, preserving rank ties."""
    denominator = np.sqrt(np.sum(xr ** 2) * np.sum(yr ** 2))
    rng = np.random.default_rng(PERMUTATION_SEED)
    extreme_two_sided = 0
    extreme_greater = 0
    batch_size = 10_000
    for start in range(0, PERMUTATIONS, batch_size):
        size = min(batch_size, PERMUTATIONS - start)
        order = np.argsort(rng.random((size, len(yr))), axis=1)
        permuted = yr[order]
        permuted_r = (permuted @ xr) / denominator
        extreme_two_sided += int(np.count_nonzero(np.abs(permuted_r) >= abs(observed)))
        extreme_greater += int(np.count_nonzero(permuted_r >= observed))
    return (
        (extreme_two_sided + 1) / (PERMUTATIONS + 1),
        (extreme_greater + 1) / (PERMUTATIONS + 1),
    )


def analyze(frame: pd.DataFrame, diagnostics: dict) -> dict:
    x = frame["SPP1_expression"].to_numpy(float)
    y = frame["TF_activity"].to_numpy(float)
    rho, xr, yr = spearman_rho(x, y)
    p_two_sided, p_greater = permutation_pvalues(xr, yr, rho)
    return {
        **diagnostics,
        "n_used": int(len(frame)),
        "spearman_rho": rho,
        "permutation_p_two_sided": float(p_two_sided),
        "permutation_p_greater": float(p_greater),
        "permutations": PERMUTATIONS,
        "permutation_seed": PERMUTATION_SEED,
    }


def fit_pathology_interaction(frame: pd.DataFrame) -> tuple[dict, list[dict]]:
    """Fit SPP1 ~ standardized TF activity * pathology with HC3 covariance."""
    y = frame["SPP1_expression"].to_numpy(float)
    activity = frame["TF_activity"].to_numpy(float)
    pathology = frame["pathology"].astype(str).to_numpy()
    mean_activity = float(activity.mean())
    sd_activity = float(activity.std(ddof=1))
    if not np.isfinite(sd_activity) or sd_activity == 0:
        raise ValueError("Pathology interaction model is unestimable: TF activity has zero variance")
    x = (activity - mean_activity) / sd_activity
    early = (pathology == "earlyAD").astype(float)
    late = (pathology == "lateAD").astype(float)
    design = np.column_stack([np.ones(len(y)), x, early, late, x * early, x * late])
    n, n_parameters = design.shape
    residual_df = n - n_parameters
    if residual_df <= 0 or np.linalg.matrix_rank(design) < n_parameters:
        raise ValueError("Pathology interaction model is rank-deficient or has no residual degrees of freedom")

    xtx_inv = np.linalg.inv(design.T @ design)
    beta = xtx_inv @ design.T @ y
    residuals = y - design @ beta
    leverage = np.sum((design @ xtx_inv) * design, axis=1)
    hc3_weights = (residuals / (1 - leverage)) ** 2
    covariance = xtx_inv @ (design.T @ (hc3_weights[:, None] * design)) @ xtx_inv
    mse = float(np.sum(residuals ** 2) / residual_df)
    cooks_d = (residuals ** 2 / (n_parameters * mse)) * leverage / (1 - leverage) ** 2
    max_leverage_idx = int(np.argmax(leverage))
    max_cooks_idx = int(np.argmax(cooks_d))

    interaction_indices = np.array([4, 5])
    interaction_beta = beta[interaction_indices]
    interaction_cov = covariance[np.ix_(interaction_indices, interaction_indices)]
    wald_chi2 = float(interaction_beta @ np.linalg.solve(interaction_cov, interaction_beta))
    interaction_f = wald_chi2 / len(interaction_indices)
    interaction_p = float(f_distribution.sf(interaction_f, len(interaction_indices), residual_df))
    model_summary = {
        "n": n,
        "activity_mean": mean_activity,
        "activity_sd": sd_activity,
        "interaction_F_HC3": interaction_f,
        "interaction_df_num": len(interaction_indices),
        "interaction_df_den": residual_df,
        "interaction_p_HC3": interaction_p,
        "max_leverage": float(leverage[max_leverage_idx]),
        "max_leverage_donor": str(frame.iloc[max_leverage_idx]["donor_id"]),
        "max_cooks_d": float(cooks_d[max_cooks_idx]),
        "max_cooks_d_donor": str(frame.iloc[max_cooks_idx]["donor_id"]),
        "covariance": "HC3",
    }

    slope_rows = []
    contrasts = {
        "nonAD": np.array([0, 1, 0, 0, 0, 0], dtype=float),
        "earlyAD": np.array([0, 1, 0, 0, 1, 0], dtype=float),
        "lateAD": np.array([0, 1, 0, 0, 0, 1], dtype=float),
    }
    for group, contrast in contrasts.items():
        estimate = float(contrast @ beta)
        se = float(np.sqrt(contrast @ covariance @ contrast))
        t_stat = estimate / se
        p_value = float(2 * t_distribution.sf(abs(t_stat), residual_df))
        critical = float(t_distribution.ppf(0.975, residual_df))
        slope_rows.append({
            "pathology": group,
            "n_group": int(np.sum(pathology == group)),
            "slope_per_1sd_TF_activity": estimate,
            "HC3_se": se,
            "CI95_low": estimate - critical * se,
            "CI95_high": estimate + critical * se,
            "p_two_sided_HC3": p_value,
        })
    return model_summary, slope_rows


def plot_pair(frame: pd.DataFrame, result: dict, outbase: Path,
              figsize: tuple[float, float]) -> None:
    fig, ax = plt.subplots(figsize=figsize)
    for pathology in PATHOLOGY_ORDER:
        sub = frame[frame["pathology"] == pathology]
        ax.scatter(sub["SPP1_expression"], sub["TF_activity"],
                   color=PATHOLOGY_COLORS[pathology], alpha=0.85, s=48,
                   edgecolor="white", linewidth=0.45, label=f"{pathology} (n={len(sub)})")
    x = frame["SPP1_expression"].to_numpy(float)
    y = frame["TF_activity"].to_numpy(float)
    slope, intercept, _, _ = theilslopes(y, x)
    grid = np.linspace(float(np.min(x)), float(np.max(x)), 100)
    ax.plot(grid, intercept + slope * grid, color="#333333", linewidth=1.6,
            linestyle="--", label="Theil–Sen fit")
    annotation = (f"Spearman rho={result['spearman_rho']:.3f}, "
                  f"two-sided p={result['permutation_p_two_sided']:.3g}, "
                  f"positive-direction p={result['permutation_p_greater']:.3g}; n={result['n_used']}")
    ax.set_title(f"{result['context']}: SPP1 expression vs {result['TF']} activity\n{annotation}", pad=12)
    ax.set_xlabel("Cleaned SPP1 gene expression")
    ax.set_ylabel(f"{result['TF']} TF activity")
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
    frames, diagnostics, results = [], [], []
    interaction_results, pathology_slopes = [], []
    for spec in PAIRS:
        frame, check = load_pair(spec)
        if len(frame) < 3:
            raise ValueError(f"Correlation is unestimable for {spec['context']}: insufficient N")
        frames.append(frame)
        diagnostics.append(check)
        association = analyze(frame, check)
        interaction, slopes = fit_pathology_interaction(frame)
        loo_interaction_p = []
        for donor_index, donor in frame.iterrows():
            reduced = frame.drop(index=donor_index).reset_index(drop=True)
            p_loo = fit_pathology_interaction(reduced)[0]["interaction_p_HC3"]
            loo_interaction_p.append((str(donor["donor_id"]), p_loo))
        min_loo = min(loo_interaction_p, key=lambda value: value[1])
        max_loo = max(loo_interaction_p, key=lambda value: value[1])
        interaction["loo_interaction_p_min"] = float(min_loo[1])
        interaction["loo_interaction_p_min_donor"] = min_loo[0]
        interaction["loo_interaction_p_max"] = float(max_loo[1])
        interaction["loo_interaction_p_max_donor"] = max_loo[0]
        results.append({**association,
                        "interaction_F_HC3": interaction["interaction_F_HC3"],
                        "interaction_p_HC3": interaction["interaction_p_HC3"]})
        interaction_results.append({
            "scope": spec["scope"], "context": spec["context"], "TF": spec["tf"],
            **{key: value for key, value in interaction.items() if key != "covariance"},
            "covariance": interaction["covariance"],
        })
        pathology_slopes.extend({
            "scope": spec["scope"], "context": spec["context"], "TF": spec["tf"], **row,
        } for row in slopes)
        if len(frames) == 1:
            joined_frames = [frame]
        else:
            joined_frames.append(frame)
    joined = pd.concat(joined_frames, ignore_index=True)
    summary = pd.DataFrame(results)

    joined.to_csv(TABLES / "spp1_tf_activity_donor_pairs.tsv", sep="\t", index=False)
    summary.to_csv(TABLES / "spp1_tf_activity_correlation_summary.tsv", sep="\t", index=False)
    pd.DataFrame(interaction_results).to_csv(
        TABLES / "spp1_tf_activity_pathology_interaction_summary.tsv", sep="\t", index=False)
    pd.DataFrame(pathology_slopes).to_csv(
        TABLES / "spp1_tf_activity_pathology_group_slopes.tsv", sep="\t", index=False)

    for frame, row in zip(frames, summary.to_dict(orient="records")):
        plot_pair(frame, row, FIGURES / f"{row['context']}_SPP1_{row['TF']}_activity_correlation",
                  figsize=(6.0, 4.7))
    fig, axes = plt.subplots(1, 3, figsize=(17.5, 5.0))
    for ax, frame, row in zip(axes, frames, summary.to_dict(orient="records")):
        for pathology in PATHOLOGY_ORDER:
            sub = frame[frame["pathology"] == pathology]
            ax.scatter(sub["SPP1_expression"], sub["TF_activity"],
                       color=PATHOLOGY_COLORS[pathology], alpha=0.82, s=32,
                       edgecolor="white", linewidth=0.35, label=pathology)
        x = frame["SPP1_expression"].to_numpy(float)
        y = frame["TF_activity"].to_numpy(float)
        slope, intercept, _, _ = theilslopes(y, x)
        grid = np.linspace(float(np.min(x)), float(np.max(x)), 100)
        ax.plot(grid, intercept + slope * grid, color="#333333", linestyle="--", linewidth=1.4)
        ax.set_title(f"{row['context']} ({row['TF']})\nSpearman rho={row['spearman_rho']:.2f}, two-sided p={row['permutation_p_two_sided']:.2g}, positive p={row['permutation_p_greater']:.2g}; n={row['n_used']}",
                     fontsize=10)
        ax.set_xlabel("Cleaned SPP1 expression")
        ax.set_ylabel("TF activity")
        ax.grid(color="#e5e5e5", linewidth=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 1.04))
    fig.suptitle("Donor-level SPP1 expression and TF activity", y=1.12, fontsize=15)
    fig.tight_layout()
    fig.savefig(FIGURES / "spp1_tf_activity_correlations_all_pairs.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / "spp1_tf_activity_correlations_all_pairs.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    report = [
        "# EXP011: SPP1 expression and selected TF activity",
        "",
        "## Analysis",
        "",
        "The donor-level association analysis measured Spearman's rank correlation between cleaned SPP1 expression and the specified TF activity within each prespecified celltype-region. Two-sided and positive-direction one-sided nominal p-values were estimated by the same 100,000 permutations of donor ranks, preserving ties. The two-sided test is primary. The positive-direction test is an exploratory sensitivity analysis added after inspecting the positive observed correlations; it is not confirmatory. This marginal rank correlation is not adjusted for pathology. All complete donor pairs were retained. The association plot shows a Theil–Sen trend line.",
        "",
        "Donors were matched by ROSMAP_IndividualID. No missing values were imputed. All shared complete donors are included; pathology is shown by color for context and was not used in the correlation.",
        "",
        "## Results",
        "",
        "| Scope | Context | TF | N | Spearman rho | two-sided permutation p (primary) | positive-direction one-sided permutation p (exploratory) |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in summary.to_dict(orient="records"):
        report.append(
            f"| {row['scope']} | {row['context']} | {row['TF']} | {row['n_used']} | "
            f"{row['spearman_rho']:.3f} | {row['permutation_p_two_sided']:.3g} | "
            f"{row['permutation_p_greater']:.3g} |"
        )
    report += [
        "",
        "## Inputs and quality checks",
        "",
    ]
    for check in diagnostics:
        row = next(item for item in results if item["context"] == check["context"])
        report.append(
            f"- {check['context']}: {check['n_shared_donors']} shared donors; "
            f"{check['n_complete_cases']} complete cases; "
            f"{check['n_excluded_after_matching']} removed for missing/invalid data; "
            f"all remaining complete donor pairs retained. "
            f"Inputs: `{check['expression_input']}`, `{check['activity_input']}`, `{check['metadata_input']}`."
        )
    report += [
        "",
        "## Pathology interaction linear model",
        "",
        "Separately, a Gaussian linear model was fit in each context: cleaned SPP1 expression ~ standardized TF activity × categorical pathology. TF activity was centered and scaled to one standard deviation. nonAD is the reference group. HC3 heteroscedasticity-robust covariance was used. The joint two-sided F test evaluates whether the TF activity slope differs across pathology groups (both interaction terms equal zero). Group-specific slope estimates are in units of SPP1 expression per 1 SD TF activity, with 95% HC3 confidence intervals. HC3 addresses heteroscedasticity, not influential observations. The interaction p-value was also recalculated after omitting each donor in turn.",
        "",
        "| Context | N | Interaction F (2, residual df) | HC3 interaction p | Leave-one-out p range (donors at min/max) |",
        "|---|---:|---:|---:|---|",
    ]
    for row in interaction_results:
        report.append(
            f"| {row['context']} | {row['n']} | F({row['interaction_df_num']}, {row['interaction_df_den']})="
            f"{row['interaction_F_HC3']:.3f} | {row['interaction_p_HC3']:.3g} | "
            f"{row['loo_interaction_p_min']:.3g} ({row['loo_interaction_p_min_donor']}) to "
            f"{row['loo_interaction_p_max']:.3g} ({row['loo_interaction_p_max_donor']}) |"
        )
    report += [
        "",
        "Group-specific slopes (nominal, two-sided HC3 intervals):",
        "",
        "| Context | Pathology | n | Slope per 1 SD TF activity | 95% CI |",
        "|---|---|---:|---:|---:|",
    ]
    for row in pathology_slopes:
        report.append(
            f"| {row['context']} | {row['pathology']} | {row['n_group']} | "
            f"{row['slope_per_1sd_TF_activity']:.4g} | "
            f"({row['CI95_low']:.4g}, {row['CI95_high']:.4g}) |"
        )
    report += [
        "",
        "The interaction model is separate from the marginal Spearman analysis. With small lateAD groups, especially n=7 per context, slope estimates may be imprecise. Leave-one-out interaction p-values cross 0.05 for Oli_AG and Mic_Immune_AG, so the apparent pathology-specific slope patterns are sensitive to individual donors and should be treated as inconclusive. Oli_EC is especially sensitive to R9176076 (maximum leverage and Cook's distance; see diagnostics), despite its leave-one-out interaction p-values remaining high. HC3 protects standard errors against heteroscedasticity but does not eliminate sensitivity to influential observations.",
        "",
        "These correlations describe donor-level co-variation and do not demonstrate direct TF→SPP1 regulation.",
        "",
        "## Outputs",
        "",
        "- `tables/spp1_tf_activity_donor_pairs.tsv`: matched donor observations.",
        "- `tables/spp1_tf_activity_correlation_summary.tsv`: Spearman rho, two-sided and exploratory positive-direction permutation p-values, and sample sizes.",
        "- `tables/spp1_tf_activity_pathology_interaction_summary.tsv`: HC3 joint pathology-by-activity interaction tests.",
        "- `tables/spp1_tf_activity_pathology_group_slopes.tsv`: pathology-specific TF activity slopes and HC3 confidence intervals.",
        "- `figures/*_SPP1_*_activity_correlation.{pdf,png}`: one plot per pair.",
        "- `figures/spp1_tf_activity_correlations_all_pairs.{pdf,png}`: combined figure.",
        "- `scripts/analyze_spp1_tf_activity_correlations.py`: reproduction script.",
        "",
        "## Software",
        "",
        f"Python {platform.python_version()}; pandas {pd.__version__}; NumPy {np.__version__}; SciPy {scipy.__version__}; Matplotlib {matplotlib.__version__}.",
    ]
    (OUT / "EXP011_results.md").write_text("\n".join(report) + "\n")
    (OUT / "session_info.txt").write_text(
        f"Python {platform.python_version()}\npandas {pd.__version__}\nNumPy {np.__version__}\n"
        f"SciPy {scipy.__version__}\nMatplotlib {matplotlib.__version__}\n"
    )
    print(summary[["context", "TF", "n_used", "spearman_rho",
                   "permutation_p_two_sided", "permutation_p_greater"]].to_string(index=False))
    print(f"\nWrote EXP011 outputs under {OUT}")


if __name__ == "__main__":
    main()
