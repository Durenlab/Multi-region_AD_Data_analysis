#!/usr/bin/env python3
"""Build concise, region-wise tables of FDR-significant CellChat trends."""

from __future__ import annotations

import csv
import math
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[3]
SOURCE_TABLES = PROJECT_ROOT / "experiments/EXP002_results/finalized/tables"
INPUT_TABLES = PROJECT_ROOT / "experiments/EXP002_results/tables"
OUTPUT_DIR = SOURCE_TABLES / "supplementary_cellchat_main_results"
REGIONS = ("AG", "EC", "HC", "MTC", "PFC", "TH")
OUTPUT_COLUMNS = (
    "region",
    "sender_cell_type",
    "receiver_cell_type",
    "pathway",
    "ligand",
    "receptor",
    "interaction_id",
    "n_nonAD",
    "n_earlyAD",
    "n_lateAD",
    "slope_per_pathology_stage",
    "CI95_low",
    "CI95_high",
    "p_value",
    "q_value_within_family",
)


def read_annotations(region: str) -> dict[str, dict[str, str]]:
    path = INPUT_TABLES / f"{region}_candidate_interactions.tsv"
    with path.open(newline="", encoding="utf-8") as handle:
        return {
            row["interaction_id"]: {
                "pathway": row.get("pathway_name", ""),
                "ligand": row.get("ligand", ""),
                "receptor": row.get("receptor", ""),
            }
            for row in csv.DictReader(handle, delimiter="\t")
        }


def is_significant(row: dict[str, str]) -> bool:
    if row.get("status") != "tested":
        return False
    try:
        q_value = float(row["q_value_within_family"])
    except (KeyError, TypeError, ValueError):
        return False
    return math.isfinite(q_value) and q_value < 0.05


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    annotations = {region: read_annotations(region) for region in REGIONS}
    writers = {}
    handles = {}
    for region in REGIONS:
        handle = (OUTPUT_DIR / f"CellChat_significant_progression_{region}.csv").open(
            "w", newline="", encoding="utf-8"
        )
        handles[region] = handle
        writers[region] = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS)
        writers[region].writeheader()

    source = SOURCE_TABLES / "pathology_progression_tests.tsv"
    try:
        with source.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                region = row.get("region", "")
                if region not in writers or not is_significant(row):
                    continue
                annotation = annotations[region].get(row.get("interaction_id", ""), {})
                writers[region].writerow(
                    {
                        "region": region,
                        "sender_cell_type": row.get("source", ""),
                        "receiver_cell_type": row.get("target", ""),
                        "pathway": annotation.get("pathway", ""),
                        "ligand": annotation.get("ligand", ""),
                        "receptor": annotation.get("receptor", ""),
                        "interaction_id": row.get("interaction_id", ""),
                        "n_nonAD": row.get("n_nonAD", ""),
                        "n_earlyAD": row.get("n_earlyAD", ""),
                        "n_lateAD": row.get("n_lateAD", ""),
                        "slope_per_pathology_stage": row.get("beta_pathology_stage", ""),
                        "CI95_low": row.get("ci_low", ""),
                        "CI95_high": row.get("ci_high", ""),
                        "p_value": row.get("p_value", ""),
                        "q_value_within_family": row.get("q_value_within_family", ""),
                    }
                )
    finally:
        for handle in handles.values():
            handle.close()

    readme = """Concise CellChat main-results tables

Each region has one CSV containing only donor-level pathology-progression tests with status=tested and within-family BH q<0.05. Regions without FDR-significant results contain only the header.

The displayed effect is the OLS slope in inferred CellChat communication probability per pathology-stage increment (nonAD=0, earlyAD=1, lateAD=2). Confidence limits are 95% HC3 intervals. The q value is BH-adjusted within the region × directed sender-to-receiver family; it is not a global correction across all regions and cell pairs. Donor counts are given for each pathology group.

These tables show only the primary progression findings. They omit non-significant and non-estimable tests, categorical pairwise contrasts, donor-level score rows, and pooled/descriptive network summaries. CellChat interactions are inferred communication hypotheses, not direct measurements or causal evidence.
"""
    (OUTPUT_DIR / "README.txt").write_text(readme, encoding="utf-8")
    print(f"Wrote concise regional tables to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
