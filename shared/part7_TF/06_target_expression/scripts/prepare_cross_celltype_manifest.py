#!/usr/bin/env python3
"""Preflight pairs for the region-shared target-expression analysis."""
import csv
from pathlib import Path

import os
DATA_ROOT = Path(os.environ.get("REGA_DATA_ROOT", Path(__file__).resolve().parents[5].parent)).expanduser().resolve()
SHARED_ROOT = Path(os.environ.get("REGA_SHARED_ROOT", DATA_ROOT / "REGA_shared")).expanduser().resolve()
ADREGA_CROSS = SHARED_ROOT
GRN_ROOT = SHARED_ROOT / "part7_TF/05GRN/results"
GEX_ROOT = ADREGA_CROSS / "part1_prepare_REGA/02cleaned_region_consistent_gex"
GENE_ASSIGNMENT_ROOT = ADREGA_CROSS / "part3_process_REGA_output/10-1gene_module_assignment_region_shareH_FC_no_Vasc"
TF_ASSIGNMENT_ROOT = ADREGA_CROSS / "part7_TF/02tf_assignment"
# Cross-celltype drivers must come from the region-wise TF activity analysis.
DRIVER_FILE = ADREGA_CROSS / "part7_TF/04TF_pathology/all_celltype_region_TF_pathology_association.tsv"
DRIVER_ROOT = Path(os.environ.get("AD_DRIVER_TF_ROOT", DATA_ROOT / "AD_Driver_TF")).expanduser().resolve()
PATHOLOGY = Path(os.environ.get("AD_PATHOLOGY_FILE", DRIVER_ROOT / "individual_pathology.tsv")).expanduser().resolve()
RESULTS_ROOT = Path(os.environ.get("TF_TARGET_RESULTS_ROOT", Path(__file__).resolve().parents[1] / "results")).expanduser().resolve()
OUT = RESULTS_ROOT / "audit/cross_celltype_pairs.tsv"


def main():
    if not DRIVER_FILE.is_file() or not PATHOLOGY.is_file():
        raise SystemExit("Missing shared pathology driver or individual pathology input")
    with DRIVER_FILE.open(newline="") as handle:
        drivers = list(csv.DictReader(handle, delimiter="\t"))
    driver_counts = {}
    for row in drivers:
        try:
            significant = float(row["lm_fdr_pathology"]) < 0.05
        except (KeyError, TypeError, ValueError):
            significant = False
        if significant:
            driver_counts[row["celltype_region"]] = driver_counts.get(row["celltype_region"], 0) + 1

    pairs = sorted(path.name.removeprefix("cleaned_").removesuffix("_region_consistent.tsv")
                   for path in GEX_ROOT.glob("cleaned_*_region_consistent.tsv"))
    if not pairs:
        raise SystemExit(f"No region-consistent expression matrices under {GEX_ROOT}")
    rows = []
    for pair in pairs:
        if pair.startswith("Vasc_"):
            continue
        region = pair.rsplit("_", 1)[1]
        gex = GEX_ROOT / f"cleaned_{pair}_region_consistent.tsv"
        grn = GRN_ROOT / f"{pair}_TG_TF.rds"
        gene_assignment = GENE_ASSIGNMENT_ROOT / region / f"{region}_ALL_cellgene_FC_assignment_with_unassigned.tsv"
        tf_assignment = TF_ASSIGNMENT_ROOT / region / f"{region}_ALL_cell_TF_FC_assignment_with_unassigned.tsv"
        checks = {
            "gex_exists": gex.is_file(),
            "grn_exists": grn.is_file(),
            "gene_assignment_exists": gene_assignment.is_file(),
            "tf_assignment_exists": tf_assignment.is_file(),
        }
        rows.append({
            "pair": pair, "region": region,
            "n_significant_driver_TF": driver_counts.get(pair, 0),
            **checks,
            "status": "ready" if all(checks.values()) else "missing_required_input",
        })
    OUT.parent.mkdir(parents=True, exist_ok=True)
    (SHARED_ROOT / "part7_TF/06_target_expression/results/logs").mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    ready = sum(row["status"] == "ready" for row in rows)
    print(f"Cross-celltype pairs: {len(rows)}; ready: {ready}; pairs with pathology-significant TF drivers: {sum(int(r['n_significant_driver_TF']) > 0 for r in rows)}")
    if ready != len(rows):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
