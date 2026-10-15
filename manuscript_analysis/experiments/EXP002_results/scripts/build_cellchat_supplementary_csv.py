#!/usr/bin/env python3
"""Combine donor-level CellChat tests into one supplementary CSV per region."""

from __future__ import annotations

import csv
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[3]
TABLE_DIR = PROJECT_ROOT / "experiments/EXP002_results/finalized/tables"
INPUT_DIR = PROJECT_ROOT / "experiments/EXP002_results/tables"
OUTPUT_DIR = TABLE_DIR / "supplementary_cellchat_csv"
REGIONS = ("AG", "EC", "HC", "MTC", "PFC", "TH")

OUTPUT_COLUMNS = (
    "analysis",
    "contrast",
    "condition_1",
    "condition_2",
    "region",
    "sender_cell_type",
    "receiver_cell_type",
    "interaction_id",
    "pathway",
    "ligand",
    "receptor",
    "n_nonAD",
    "n_earlyAD",
    "n_lateAD",
    "n_condition_1",
    "n_condition_2",
    "mean_score_nonAD",
    "mean_score_earlyAD",
    "mean_score_lateAD",
    "effect_measure",
    "effect_estimate",
    "HC3_standard_error",
    "CI95_low",
    "CI95_high",
    "p_value",
    "q_value_within_family",
    "FDR_family",
    "status",
    "covariates",
    "covariate_decisions",
    "model_formula",
    "design_rank",
    "residual_degrees_of_freedom",
)


def clean(value: object) -> str:
    text = "" if value is None else str(value)
    return "" if text in {"NA", "NaN", "nan", "NULL"} else text


def read_annotations(region: str) -> dict[str, dict[str, str]]:
    path = INPUT_DIR / f"{region}_candidate_interactions.tsv"
    with path.open(newline="", encoding="utf-8") as handle:
        return {
            row["interaction_id"]: {
                "pathway": clean(row.get("pathway_name")),
                "ligand": clean(row.get("ligand")),
                "receptor": clean(row.get("receptor")),
            }
            for row in csv.DictReader(handle, delimiter="\t")
        }


def progression_row(row: dict[str, str], annotation: dict[str, str]) -> dict[str, str]:
    result = {column: "" for column in OUTPUT_COLUMNS}
    result.update(
        analysis="Progression trend",
        contrast="nonAD (0) → earlyAD (1) → lateAD (2)",
        region=clean(row.get("region")),
        sender_cell_type=clean(row.get("source")),
        receiver_cell_type=clean(row.get("target")),
        interaction_id=clean(row.get("interaction_id")),
        n_nonAD=clean(row.get("n_nonAD")),
        n_earlyAD=clean(row.get("n_earlyAD")),
        n_lateAD=clean(row.get("n_lateAD")),
        mean_score_nonAD=clean(row.get("mean_nonAD")),
        mean_score_earlyAD=clean(row.get("mean_earlyAD")),
        mean_score_lateAD=clean(row.get("mean_lateAD")),
        effect_measure="OLS pathology-stage slope",
        effect_estimate=clean(row.get("beta_pathology_stage")),
        HC3_standard_error=clean(row.get("se_HC3")),
        CI95_low=clean(row.get("ci_low")),
        CI95_high=clean(row.get("ci_high")),
        p_value=clean(row.get("p_value")),
        q_value_within_family=clean(row.get("q_value_within_family")),
        FDR_family=clean(row.get("FDR_family")),
        status=clean(row.get("status")),
        covariates=clean(row.get("covariates")),
        covariate_decisions=clean(row.get("covariate_decisions")),
        model_formula=clean(row.get("model_formula")),
        design_rank=clean(row.get("design_rank")),
        residual_degrees_of_freedom=clean(row.get("residual_df")),
    )
    result.update(annotation)
    return result


def pairwise_row(row: dict[str, str], annotation: dict[str, str]) -> dict[str, str]:
    condition_1 = clean(row.get("condition_1"))
    condition_2 = clean(row.get("condition_2"))
    result = {column: "" for column in OUTPUT_COLUMNS}
    result.update(
        analysis="Categorical pathology contrast",
        contrast=f"{condition_2} − {condition_1}",
        condition_1=condition_1,
        condition_2=condition_2,
        region=clean(row.get("region")),
        sender_cell_type=clean(row.get("source")),
        receiver_cell_type=clean(row.get("target")),
        interaction_id=clean(row.get("interaction_id")),
        n_condition_1=clean(row.get("n_condition_1")),
        n_condition_2=clean(row.get("n_condition_2")),
        effect_measure="Covariate-adjusted score difference (condition 2 − condition 1)",
        effect_estimate=clean(row.get("effect_condition2_minus_condition1")),
        HC3_standard_error=clean(row.get("se_HC3")),
        CI95_low=clean(row.get("ci_low")),
        CI95_high=clean(row.get("ci_high")),
        p_value=clean(row.get("p_value")),
        q_value_within_family=clean(row.get("q_value_within_family")),
        status=clean(row.get("status")),
        covariates=clean(row.get("covariates")),
        covariate_decisions=clean(row.get("covariate_decisions")),
        design_rank=clean(row.get("design_rank")),
        residual_degrees_of_freedom=clean(row.get("residual_df")),
    )
    for group in ("nonAD", "earlyAD", "lateAD"):
        if group == condition_1:
            result[f"n_{group}"] = result["n_condition_1"]
        elif group == condition_2:
            result[f"n_{group}"] = result["n_condition_2"]
    result.update(annotation)
    return result


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    handles = {}
    writers = {}
    annotations = {}
    try:
        for region in REGIONS:
            annotations[region] = read_annotations(region)
            path = OUTPUT_DIR / f"Supplementary_Table_CellChat_{region}.csv"
            handle = path.open("w", newline="", encoding="utf-8")
            handles[region] = handle
            writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, extrasaction="ignore")
            writers[region] = writer
            writer.writeheader()

        progression_path = TABLE_DIR / "pathology_progression_tests.tsv"
        with progression_path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                region = row["region"]
                if region in writers:
                    annotation = annotations[region].get(row["interaction_id"], {})
                    writers[region].writerow(progression_row(row, annotation))

        pairwise_path = TABLE_DIR / "pathology_pairwise_contrasts.tsv"
        with pairwise_path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                region = row["region"]
                if region in writers:
                    annotation = annotations[region].get(row["interaction_id"], {})
                    writers[region].writerow(pairwise_row(row, annotation))
    finally:
        for handle in handles.values():
            handle.close()

    readme = """CellChat supplementary donor-level results

This directory contains one CSV per brain region. Each file combines the complete finalized donor-level progression tests and categorical pathology contrasts for that region. Pooled/descriptive network comparisons are excluded.

Progression rows report the donor-level OLS pathology-stage slope using nonAD=0, earlyAD=1, and lateAD=2, with HC3 standard errors and the finalized within-family BH q value. Categorical contrast rows report the previously finalized covariate-adjusted contrast (condition 2 minus condition 1) and its source p/q values. The two analysis types therefore retain their original models and correction values; they were not refit for this table.

Blank statistical fields indicate missing or non-estimable values in the source result. Status fields distinguish tested rows from no-score-variation and other untestable rows. Donor counts are given by pathology group for progression tests and by condition for categorical contrasts. Ligand, receptor, and pathway labels were joined using the corresponding region's candidate-interaction metadata.

Scores and interactions are CellChat-inferred communication hypotheses, not direct measurements of signaling or causal effects.
"""
    (OUTPUT_DIR / "README.txt").write_text(readme, encoding="utf-8")
    print(f"Wrote six regional supplementary CSV files to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
