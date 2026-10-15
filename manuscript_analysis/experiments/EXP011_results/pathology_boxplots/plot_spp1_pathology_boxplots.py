#!/usr/bin/env python3
"""Plot cleaned SPP1 expression by pathology with a right-tailed linear trend test."""
from __future__ import annotations

from pathlib import Path
import platform

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.stats import t as student_t

ROOT = Path(__file__).resolve().parents[3]
PALMETTO = Path("/project/zduren/durenlab/naqing_palmetto")
OUT = Path(__file__).resolve().parent
FIGURES = OUT / "figures"
TABLES = OUT / "tables"
PATHOLOGY_ORDER = ["nonAD", "earlyAD", "lateAD"]
PATHOLOGY_SCORE = {"nonAD": 0, "earlyAD": 1, "lateAD": 2}
COLORS = {"nonAD": "#4477AA", "earlyAD": "#EE9944", "lateAD": "#CC6677"}

PAIRS = [
    {
        "context": "Oli_EC",
        "expression": PALMETTO / "AD_REGA_region_wise/part1_prepare_REGA/02cleaned_region_consistent_gex/cleaned_Oli_EC_region_consistent.tsv",
        "metadata": PALMETTO / "AD_REGA_region_wise/part1_prepare_REGA/01region_consistent_metadata/Oli_EC_region_consistent_metadata.tsv",
    },
    {
        "context": "Mic_Immune_AG",
        "expression": PALMETTO / "AD_REGA_region_wise/part1_prepare_REGA/02cleaned_region_consistent_gex/cleaned_Mic_Immune_AG_region_consistent.tsv",
        "metadata": PALMETTO / "AD_REGA_region_wise/part1_prepare_REGA/01region_consistent_metadata/Mic_Immune_AG_region_consistent_metadata.tsv",
    },
    {
        "context": "Oli_AG",
        "expression": PALMETTO / "AD_REGA/part1_prepare_REGA/04cleaned_gex/cleaned_Oli_AG.tsv",
        "metadata": PALMETTO / "AD_REGA/part1_prepare_REGA/03metadata_by_group/Oli_AG_metadata.tsv",
    },
]


def load_expression(spec: dict) -> pd.DataFrame:
    for key in ("expression", "metadata"):
        if not spec[key].exists():
            raise FileNotFoundError(f"Missing {key} input: {spec[key]}")
    expression = pd.read_csv(spec["expression"], sep="\t", index_col=0)
    expression.index = expression.index.astype(str).str.strip().str.upper()
    expression.columns = expression.columns.astype(str).str.strip()
    if (expression.index == "SPP1").sum() != 1:
        raise ValueError(f"Expected one SPP1 row in {spec['expression']}")

    metadata = pd.read_csv(spec["metadata"], sep="\t", dtype={"ROSMAP_IndividualID": str})
    metadata["ROSMAP_IndividualID"] = metadata["ROSMAP_IndividualID"].astype(str).str.strip()
    metadata["Pathology"] = metadata["Pathology"].astype(str).str.strip()
    if metadata["ROSMAP_IndividualID"].duplicated().any():
        raise ValueError(f"Duplicate donor IDs in {spec['metadata']}")
    meta = metadata.set_index("ROSMAP_IndividualID")
    donors = sorted(set(expression.columns) & set(meta.index))
    frame = pd.DataFrame({
        "context": spec["context"],
        "donor_id": donors,
        "pathology": meta.loc[donors, "Pathology"].to_numpy(),
        "SPP1_expression": pd.to_numeric(expression.loc["SPP1", donors], errors="coerce").to_numpy(),
    })
    frame = frame[frame["pathology"].isin(PATHOLOGY_ORDER)
                  & np.isfinite(frame["SPP1_expression"])].copy()
    if frame["donor_id"].duplicated().any():
        raise ValueError(f"Duplicate matched donor IDs in {spec['context']}")
    frame["pathology"] = pd.Categorical(frame["pathology"], categories=PATHOLOGY_ORDER, ordered=True)
    frame["pathology_score"] = frame["pathology"].astype(str).map(PATHOLOGY_SCORE).astype(int)
    return frame


def right_tailed_linear_trend(frame: pd.DataFrame) -> dict:
    """Test H1: SPP1 expression has a positive linear trend across pathology order."""
    y = frame["SPP1_expression"].to_numpy(float)
    x = frame["pathology_score"].to_numpy(float)
    design = np.column_stack([np.ones(len(x)), x])
    beta, _, rank, _ = np.linalg.lstsq(design, y, rcond=None)
    if rank < 2:
        raise ValueError(f"Ordered pathology trend is unestimable for {frame['context'].iloc[0]}")
    residual = y - design @ beta
    df = len(y) - 2
    sigma2 = float(residual @ residual / df)
    covariance = sigma2 * np.linalg.inv(design.T @ design)
    slope = float(beta[1])
    se = float(np.sqrt(covariance[1, 1]))
    t_stat = slope / se
    p_greater = float(student_t.sf(t_stat, df))
    p_two_sided = float(2 * student_t.sf(abs(t_stat), df))
    return {"n": int(len(y)), "slope_per_pathology_step": slope,
            "slope_se": se, "t_statistic": float(t_stat), "df": int(df),
            "p_right_tailed_greater": p_greater, "p_two_sided": p_two_sided}


def plot_panel(ax, frame: pd.DataFrame, result: dict) -> None:
    data = [frame.loc[frame["pathology"] == group, "SPP1_expression"].to_numpy(float)
            for group in PATHOLOGY_ORDER]
    artists = ax.boxplot(data, positions=np.arange(1, 4), widths=0.55,
                         patch_artist=True, showfliers=False,
                         medianprops={"color": "#222222", "linewidth": 1.5},
                         whiskerprops={"color": "#555555"},
                         capprops={"color": "#555555"})
    for patch, group in zip(artists["boxes"], PATHOLOGY_ORDER):
        patch.set_facecolor(COLORS[group])
        patch.set_alpha(0.35)
    rng = np.random.default_rng(1108)
    for xpos, group, values in zip(np.arange(1, 4), PATHOLOGY_ORDER, data):
        jitter = rng.uniform(-0.12, 0.12, len(values))
        ax.scatter(xpos + jitter, values, color=COLORS[group], s=25, alpha=0.82,
                   edgecolor="white", linewidth=0.35, zorder=3,
                   label=f"{group} (n={len(values)})")
    ax.set_xticks(np.arange(1, 4), PATHOLOGY_ORDER)
    ax.set_title(
        f"{frame['context'].iloc[0]}\n"
        f"p-value={result['p_two_sided']:.3g}; "
        f"right-tailed p-value={result['p_right_tailed_greater']:.3g}; "
        f"slope={result['slope_per_pathology_step']:.3g}; n={result['n']}",
        fontsize=10,
    )
    ax.set_xlabel("Pathology group (ordered)")
    ax.set_ylabel("Cleaned SPP1 expression")
    ax.grid(axis="y", color="#e5e5e5", linewidth=0.6)
    ax.spines[["top", "right"]].set_visible(False)


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    frames, results = [], []
    for spec in PAIRS:
        frame = load_expression(spec)
        result = {"context": spec["context"], **right_tailed_linear_trend(frame),
                  "expression_input": str(spec["expression"]),
                  "metadata_input": str(spec["metadata"])}
        frames.append(frame)
        results.append(result)

    donors = pd.concat(frames, ignore_index=True)
    donors.to_csv(TABLES / "spp1_expression_pathology_donors.tsv", sep="\t", index=False)
    pd.DataFrame(results).to_csv(TABLES / "spp1_expression_pathology_trend.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(14.5, 5.2), sharey=False)
    for ax, frame, result in zip(axes, frames, results):
        plot_panel(ax, frame, result)
    fig.suptitle("SPP1 expression by pathology", fontsize=15, y=1.02)
    fig.tight_layout()
    for ext, kwargs in [("pdf", {}), ("png", {"dpi": 300})]:
        fig.savefig(FIGURES / f"spp1_expression_by_pathology_three_contexts.{ext}",
                    bbox_inches="tight", **kwargs)
    plt.close(fig)

    for frame, result in zip(frames, results):
        fig, ax = plt.subplots(figsize=(5.5, 4.8))
        plot_panel(ax, frame, result)
        fig.tight_layout()
        stem = FIGURES / f"{result['context']}_SPP1_expression_by_pathology"
        fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
        fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
        plt.close(fig)

    report = [
        "# SPP1 expression by pathology",
        "",
        "Boxplots show cleaned SPP1 expression in nonAD, earlyAD, and lateAD donors. The linear trend regresses expression on ordered pathology scores (nonAD=0, earlyAD=1, lateAD=2). The figure reports the two-sided p-value as p-value and the positive-direction one-sided p-value as right-tailed p-value (H1: slope > 0). Donors with complete expression and pathology data were included; points are shown individually. No covariate adjustment or multiple-testing correction was applied.",
        "",
        "| Context | n | nonAD | earlyAD | lateAD | Slope per pathology step | p-value | Right-tailed p-value |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for frame, result in zip(frames, results):
        counts = frame["pathology"].value_counts().reindex(PATHOLOGY_ORDER, fill_value=0)
        report.append(
            f"| {result['context']} | {result['n']} | {counts['nonAD']} | {counts['earlyAD']} | "
            f"{counts['lateAD']} | {result['slope_per_pathology_step']:.4g} | "
            f"{result['p_two_sided']:.3g} | {result['p_right_tailed_greater']:.3g} |"
        )
    report += ["", "## Inputs", ""]
    for result in results:
        report.append(f"- {result['context']}: expression `{result['expression_input']}`; metadata `{result['metadata_input']}`.")
    report += ["", "## Software", "",
               f"Python {platform.python_version()}; pandas {pd.__version__}; NumPy {np.__version__}; SciPy {scipy.__version__}; Matplotlib {matplotlib.__version__}.",
               ""]
    (OUT / "SPP1_expression_by_pathology.md").write_text("\n".join(report))
    print(pd.DataFrame(results)[["context", "n", "slope_per_pathology_step", "p_right_tailed_greater"]].to_string(index=False))
    print(f"\nWrote figures and tables under {OUT}")


if __name__ == "__main__":
    main()
