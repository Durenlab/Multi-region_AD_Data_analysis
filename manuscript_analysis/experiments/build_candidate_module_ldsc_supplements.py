#!/usr/bin/env python3
"""Refresh manuscript Table 3 from the promoted min10, non-Vasc LDSC results."""

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = Path("/project/zduren/durenlab/naqing_palmetto")
CONFIGS = (
    (
        BASE / "REGA_separate/02_run_rega/pathology_association/selected_module.txt",
        BASE / "REGA_separate/04_ldsc/results/h2/AD_Bellenguez_single_module_extend2kb_LDSC_results.tsv",
        ROOT / "manuscript/supplementary_table3-S_LDSC_Cell_type_specific.tsv",
        23,
    ),
    (
        BASE / "REGA_shared/current/part2_run_REGA/08pathology_tests/selected_module.txt",
        BASE / "REGA_shared/current/part4_ldsc/heritability_results/AD_Bellenguez_single_module_extend2kb_LDSC_results.tsv",
        ROOT / "manuscript/supplementary_table3-S_LDSC_Cross_celltype.tsv",
        36,
    ),
)


def read_tsv(path):
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        return reader.fieldnames, list(reader)


for selection, source, output, expected_rows in CONFIGS:
    selected = set(selection.read_text().split())
    columns, rows = read_tsv(source)
    if not {"Category", "module", "celltype_region", "FDR_BH_Enrichment_p"}.issubset(columns):
        raise ValueError(f"Missing LDSC columns in {source}")
    is_shared = expected_rows == 36
    def candidate(row):
        return f"{row['celltype_region'].rsplit('_', 1)[-1]}_{row['module']}" if is_shared else row["Category"]
    rows = [row for row in rows if candidate(row) in selected and not row["celltype_region"].startswith("Vasc_")]
    if len(rows) != expected_rows or len({row["Category"] for row in rows}) != expected_rows:
        raise ValueError(f"Expected {expected_rows} unique current LDSC rows, found {len(rows)} in {source}")
    # The source includes local result/log paths; keep the established manuscript schema.
    columns = [col for col in columns if col not in {"result_file", "log_file"}]
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"{output}: {len(rows)} rows from {source}")
