#!/usr/bin/env python3
"""Plot significant AG SPP1-integrin features sent to Mic_Immune.

Reuses finalized donor-level CellChat tests and the selected-features plotting
style. Missing feature rows are zero for eligible completed donor networks,
matching the established EXP002 convention; no models are refit.
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
EXP = ROOT / "experiments" / "EXP002_results"
FINAL = EXP / "finalized"
TABLE = FINAL / "tables"
INTERMEDIATE = EXP / "intermediate"
FIGURES = FINAL / "figures"
PATHOLOGY = ["nonAD", "earlyAD", "lateAD"]
COLORS = {"nonAD": "#4477AA", "earlyAD": "#EE9944", "lateAD": "#CC6677"}
SEED = 20261001


def select_and_reconstruct() -> tuple[pd.DataFrame, pd.DataFrame]:
    tests = pd.read_csv(TABLE / "pathology_progression_tests.tsv", sep="\t")
    features = tests.loc[
        tests["region"].eq("AG")
        & tests["source"].isin(["Mic_Immune", "Oli"])
        & tests["target"].eq("Mic_Immune")
        & tests["interaction_id"].str.startswith("SPP1_ITG", na=False)
        & tests["q_value_within_family"].lt(0.05)
    ].copy()
    if len(features) != 2:
        raise ValueError(f"Expected two significant AG SPP1-integrin features to Mic_Immune; found {len(features)}")

    status = pd.read_csv(EXP / "tables/donor_network_status.tsv", sep="\t")
    coverage = pd.read_csv(EXP / "tables/AG_donor_celltype_coverage.tsv", sep="\t")
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
    reconstructed = []
    for feature in features.itertuples(index=False):
        eligible = completed.merge(
            coverage.loc[coverage["celltype"].eq(feature.source), ["donor_id"]]
                .assign(source_eligible=True), on="donor_id", how="inner"
        ).merge(
            coverage.loc[coverage["celltype"].eq(feature.target), ["donor_id"]]
                .assign(target_eligible=True), on="donor_id", how="inner"
        )
        observed = nonzero.loc[
            nonzero["region"].eq("AG")
            & nonzero["source"].eq(feature.source)
            & nonzero["target"].eq(feature.target)
            & nonzero["interaction_id"].eq(feature.interaction_id),
            ["donor_id", "communication_probability"],
        ]
        if observed["donor_id"].duplicated().any():
            raise ValueError(f"Duplicate donor scores for {feature.interaction_id}")
        scores = eligible.merge(observed, on="donor_id", how="left", validate="one_to_one")
        scores = scores.rename(columns={"communication_probability": "score"})
        scores["score_row_present"] = scores["score"].notna()
        scores["score"] = scores["score"].fillna(0.0)
        scores["region"] = "AG"
        scores["source"] = feature.source
        scores["target"] = feature.target
        scores["interaction_id"] = feature.interaction_id

        counts = scores.groupby("pathology")["donor_id"].nunique().reindex(PATHOLOGY, fill_value=0)
        expected_counts = [feature.n_nonAD, feature.n_earlyAD, feature.n_lateAD]
        if counts.tolist() != expected_counts:
            raise ValueError(f"Eligibility counts mismatch for {feature.interaction_id}: {counts.tolist()} != {expected_counts}")
        means = scores.groupby("pathology")["score"].mean().reindex(PATHOLOGY)
        expected_means = [feature.mean_nonAD, feature.mean_earlyAD, feature.mean_lateAD]
        if not np.allclose(means.to_numpy(float), expected_means, rtol=1e-7, atol=1e-9):
            raise ValueError(f"Reconstructed means do not match finalized test for {feature.interaction_id}")
        reconstructed.append(scores[[
            "region", "donor_id", "pathology", "source", "target", "interaction_id",
            "score", "score_row_present",
        ]])
    return features, pd.concat(reconstructed, ignore_index=True)


def plot(features: pd.DataFrame, scores: pd.DataFrame, path: Path) -> None:
    ncols = 2
    nrows = math.ceil(len(features) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(8.5, 3.55 * nrows), squeeze=False)
    rng = np.random.default_rng(SEED)
    for ax, row in zip(axes.flat, features.itertuples(index=False)):
        obs = scores.loc[
            scores["source"].eq(row.source)
            & scores["target"].eq(row.target)
            & scores["interaction_id"].eq(row.interaction_id)
        ]
        for x, pathology in enumerate(PATHOLOGY):
            values = obs.loc[obs["pathology"].eq(pathology), "score"].to_numpy(float)
            ax.scatter(x + rng.normal(0, 0.055, len(values)), values,
                       color=COLORS[pathology], alpha=0.8, s=21,
                       edgecolor="white", linewidth=0.3, zorder=3)
            ax.text(x, 0.98, f"n={len(values)}", transform=ax.get_xaxis_transform(),
                    ha="center", va="top", fontsize=8.5)
        fitted = [row.prediction_nonAD_mean_covars,
                  row.prediction_earlyAD_mean_covars,
                  row.prediction_lateAD_mean_covars]
        if np.isfinite(fitted).all():
            ax.plot(range(3), fitted, color="#777777", linestyle="--",
                    linewidth=1.8, zorder=2, label="Simple linear pathology trend")
        beta = "β=NA" if pd.isna(row.beta_pathology_stage) else (
            f"β={row.beta_pathology_stage:.3g} [{row.ci_low:.3g}, {row.ci_high:.3g}]"
        )
        qvalue = "BH q=NA" if pd.isna(row.q_value_within_family) else f"BH q={row.q_value_within_family:.3g}"
        ax.set_title(f"AG: {row.source} → {row.target}\n{row.interaction_id}\n{beta}; {qvalue}", fontsize=9)
        ax.set_xticks(range(3), PATHOLOGY, rotation=15, ha="right")
        ax.set_ylabel("CellChat probability")
        ax.grid(axis="y", color="#dddddd", linewidth=0.55)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes.flat[len(features):]:
        ax.set_visible(False)
    fig.suptitle("Significant AG SPP1–integrin features sent to Mic_Immune",
                 fontsize=15, y=1.003)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    features, scores = select_and_reconstruct()
    FIGURES.mkdir(parents=True, exist_ok=True)
    features.to_csv(TABLE / "AG_SPP1_integrin_to_Mic_Immune_significant.tsv", sep="\t", index=False)
    scores.to_csv(TABLE / "AG_SPP1_integrin_to_Mic_Immune_donor_scores.tsv", sep="\t", index=False)
    output = FIGURES / "AG_SPP1_integrin_to_Mic_Immune_significant"
    plot(features, scores, output)
    print(f"Significant features: {len(features)}")
    print(f"Figure: {output}.pdf and {output}.png")


if __name__ == "__main__":
    main()
