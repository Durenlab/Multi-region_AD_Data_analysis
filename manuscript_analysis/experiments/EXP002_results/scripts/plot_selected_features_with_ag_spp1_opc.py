#!/usr/bin/env python3
"""Extend the significant selected-feature plot with AG SPP1→OPC features.

Existing selected-feature scores and model results are reused as finalized.
For the six additional AG SPP1–integrin features, donor scores are recovered
from the stored nonzero probabilities; eligible completed donor networks with
no recorded feature score are represented as zero, following the EXP002
missingness/zero convention. No statistical model is refit here.
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "experiments" / "EXP002_results" / "finalized"
TABLE = OUT / "tables"
INTERMEDIATE = ROOT / "experiments" / "EXP002_results" / "intermediate"
PATHOLOGY_ORDER = ["nonAD", "earlyAD", "lateAD"]
COLORS = {"nonAD": "#4477AA", "earlyAD": "#EE9944", "lateAD": "#CC6677"}
SEED = 20261001


def load_additional_ag_features() -> tuple[pd.DataFrame, pd.DataFrame]:
    tests = pd.read_csv(TABLE / "pathology_progression_tests.tsv", sep="\t")
    extra = tests.loc[
        tests["region"].eq("AG")
        & tests["source"].isin(["Mic_Immune", "Oli"])
        & tests["target"].eq("OPC")
        & tests["interaction_id"].str.startswith("SPP1_", na=False)
        & tests["q_value_within_family"].lt(0.05)
    ].copy()
    if len(extra) != 6:
        raise ValueError(f"Expected six significant AG SPP1→OPC features; found {len(extra)}")

    status = pd.read_csv(ROOT / "experiments/EXP002_results/tables/donor_network_status.tsv", sep="\t")
    coverage = pd.read_csv(ROOT / "experiments/EXP002_results/tables/AG_donor_celltype_coverage.tsv", sep="\t")
    coverage = coverage.loc[
        coverage["region"].eq("AG")
        & coverage["retained_for_network"].astype(str).str.lower().eq("true"),
        ["donor_id", "celltype"],
    ]
    completed = status.loc[
        status["region"].eq("AG") & status["donor_network_status"].eq("completed"),
        ["region", "donor_id", "pathology"],
    ].drop_duplicates()
    if completed["donor_id"].duplicated().any():
        raise ValueError("AG donor network status has duplicate donor records")

    nonzero = pd.read_csv(INTERMEDIATE / "AG_donor_nonzero_probabilities.tsv.gz", sep="\t")
    new_score_rows: list[pd.DataFrame] = []
    for row in extra.itertuples(index=False):
        eligible = completed.merge(
            coverage.loc[coverage["celltype"].eq(row.source), ["donor_id"]].assign(source_eligible=True),
            on="donor_id", how="inner",
        ).merge(
            coverage.loc[coverage["celltype"].eq(row.target), ["donor_id"]].assign(target_eligible=True),
            on="donor_id", how="inner",
        )
        observed = nonzero.loc[
            nonzero["region"].eq("AG")
            & nonzero["source"].eq(row.source)
            & nonzero["target"].eq(row.target)
            & nonzero["interaction_id"].eq(row.interaction_id),
            ["donor_id", "communication_probability"],
        ]
        if observed["donor_id"].duplicated().any():
            raise ValueError(f"Duplicate nonzero donor scores for {row.interaction_id}")
        scores = eligible.merge(observed, on="donor_id", how="left", validate="one_to_one")
        scores = scores.rename(columns={"communication_probability": "score"})
        scores["score_row_present"] = scores["score"].notna()
        scores["score_status"] = np.where(scores["score_row_present"], "nonzero_recorded", "eligible_zero")
        scores["score"] = scores["score"].fillna(0.0)
        scores["region"] = "AG"
        scores["source"] = row.source
        scores["target"] = row.target
        scores["interaction_id"] = row.interaction_id

        observed_counts = scores.groupby("pathology")["donor_id"].nunique().reindex(PATHOLOGY_ORDER, fill_value=0)
        expected_counts = [row.n_nonAD, row.n_earlyAD, row.n_lateAD]
        if observed_counts.tolist() != expected_counts:
            raise ValueError(
                f"Eligibility counts mismatch for {row.interaction_id}: "
                f"{observed_counts.tolist()} != {expected_counts}"
            )
        observed_means = scores.groupby("pathology")["score"].mean().reindex(PATHOLOGY_ORDER)
        expected_means = [row.mean_nonAD, row.mean_earlyAD, row.mean_lateAD]
        if not np.allclose(observed_means.to_numpy(float), expected_means, rtol=1e-7, atol=1e-9):
            raise ValueError(f"Reconstructed score means do not match finalized tests for {row.interaction_id}")
        new_score_rows.append(scores[[
            "region", "donor_id", "pathology", "source", "target", "interaction_id",
            "score", "score_row_present", "score_status",
        ]])

    return extra, pd.concat(new_score_rows, ignore_index=True)


def plot(features: pd.DataFrame, scores: pd.DataFrame, path: Path) -> None:
    ncols = 2
    nrows = math.ceil(len(features) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(8.5, 3.55 * nrows), squeeze=False)
    rng = np.random.default_rng(SEED)
    for ax, row in zip(axes.flat, features.itertuples(index=False)):
        obs = scores.loc[
            scores["region"].eq(row.region)
            & scores["source"].eq(row.source)
            & scores["target"].eq(row.target)
            & scores["interaction_id"].eq(row.interaction_id)
        ]
        for x, pathology in enumerate(PATHOLOGY_ORDER):
            values = obs.loc[obs["pathology"].eq(pathology), "score"].to_numpy(float)
            jitter = rng.normal(0, 0.055, len(values))
            ax.scatter(x + jitter, values, color=COLORS[pathology], alpha=0.8,
                       s=21, edgecolor="white", linewidth=0.3, zorder=3)
            ax.text(x, 0.98, f"n={len(values)}", transform=ax.get_xaxis_transform(),
                    ha="center", va="top", fontsize=8.5)
        fitted = [row.prediction_nonAD_mean_covars,
                  row.prediction_earlyAD_mean_covars,
                  row.prediction_lateAD_mean_covars]
        if np.isfinite(fitted).all():
            ax.plot(range(3), fitted, color="#777777", linestyle="--", linewidth=1.8, zorder=2)
        beta = "β=NA" if pd.isna(row.beta_pathology_stage) else (
            f"β={row.beta_pathology_stage:.3g} [{row.ci_low:.3g}, {row.ci_high:.3g}]"
        )
        qvalue = "BH q=NA" if pd.isna(row.q_value_within_family) else f"BH q={row.q_value_within_family:.3g}"
        ax.set_title(f"{row.region}: {row.source} → {row.target}\n{row.interaction_id}\n{beta}; {qvalue}", fontsize=9)
        ax.set_xticks(range(3), PATHOLOGY_ORDER, rotation=15, ha="right")
        ax.set_ylabel("CellChat probability")
        ax.grid(axis="y", color="#dddddd", linewidth=0.55)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes.flat[len(features):]:
        ax.set_visible(False)
    fig.suptitle("Significant selected features plus AG SPP1–integrin signals to OPC",
                 fontsize=15, y=1.003)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    fixed = pd.read_csv(TABLE / "fixed_selected_features_refreshed.tsv", sep="\t")
    fixed = fixed.loc[fixed["q_value_within_family"].lt(0.05)].copy()
    fixed_scores = pd.read_csv(TABLE / "fixed_selected_feature_donor_scores.tsv", sep="\t")
    extra, extra_scores = load_additional_ag_features()
    fixed["feature_set"] = "Original fixed selected significant feature"
    extra["feature_set"] = "Added AG SPP1–integrin to OPC"
    all_features = pd.concat([fixed, extra], ignore_index=True)
    all_features = all_features.sort_values(
        ["feature_set", "region", "source", "target", "interaction_id"],
        key=lambda s: s.map({"Mic_Immune": 0, "Oli": 1}) if s.name == "source" else s,
    )
    all_scores = pd.concat([fixed_scores, extra_scores], ignore_index=True)

    table_path = TABLE / "selected_features_significant_with_AG_SPP1_OPC.tsv"
    score_path = TABLE / "selected_feature_donor_scores_with_AG_SPP1_OPC.tsv"
    all_features.to_csv(table_path, sep="\t", index=False)
    all_scores.to_csv(score_path, sep="\t", index=False)
    figbase = OUT / "figures" / "selected_features_significant_with_AG_SPP1_OPC"
    plot(all_features, all_scores, figbase)
    print(f"features: {len(all_features)} total ({len(fixed)} existing, {len(extra)} added AG SPP1→OPC)")
    print(f"figure: {figbase}.pdf / .png")
    print(f"tables: {table_path}; {score_path}")


if __name__ == "__main__":
    main()
