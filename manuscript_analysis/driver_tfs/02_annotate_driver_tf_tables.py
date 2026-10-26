#!/usr/bin/env python3
"""Add MAGMA-list and pathology-module flags to manuscript Supplementary Table 5."""
from __future__ import annotations

import csv
import os
from pathlib import Path


ROOT = Path(os.environ.get("AD_MANUSCRIPT_WORK_ROOT", Path(__file__).resolve().parents[3])).expanduser().resolve()
SNP_LIST = Path(os.environ.get("AD_SNP_TF_LIST", ROOT / "experiments/TF_prioritize/TF_with_SNP/TFs_with_SNP.txt")).expanduser().resolve()
BASE = Path(os.environ.get("REGA_DATA_ROOT", ROOT.parent)).expanduser().resolve()
TABLES = {
    ROOT / "manuscript/supplementary_table5-Cell_type_specific_Driver_TF.tsv": {
        "scope": "celltype-specific",
        "modules": BASE / "REGA_separate/02_run_rega/pathology_association/selected_module.txt",
        "mirror": ROOT / "manuscript/supp/Supplementary_Table_celltype_specific_driver_TFs.tsv",
    },
    ROOT / "manuscript/supplementary_table5-Cross_celltype_Driver_TF.tsv": {
        "scope": "cross-celltype",
        "modules": BASE / "REGA_shared/part2_run_REGA/08pathology_tests/selected_module.txt",
        "mirror": ROOT / "manuscript/supp/Supplementary_Table_cross_celltype_driver_TFs.tsv",
    },
}


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        rows = list(reader)
        if reader.fieldnames is None:
            raise ValueError(f"Missing header: {path}")
        return reader.fieldnames, rows


def main() -> None:
    snp_tfs = {
        line.strip().upper()
        for line in SNP_LIST.read_text().splitlines()
        if line.strip()
    }
    for path, config in TABLES.items():
        columns, rows = read_tsv(path)
        if not {"celltype_region", "TF", "TF_assignment"}.issubset(columns):
            raise ValueError(f"Unexpected Supplementary Table 5 schema: {path}")
        disease_modules = set(config["modules"].read_text().split())
        flags = ["is_TF_with_SNP", "is_assigned_to_disease_associated_module"]
        columns = [column for column in columns if column not in flags] + flags
        for row in rows:
            tf = row["TF"].strip().upper()
            context = (
                row["region"].strip()
                if config["scope"] == "cross-celltype"
                else row["celltype_region"].strip()
            )
            assigned = f'{context}_{row["TF_assignment"].strip()}'
            row["is_TF_with_SNP"] = "Yes" if tf in snp_tfs else "No"
            row["is_assigned_to_disease_associated_module"] = (
                "Yes" if assigned in disease_modules else "No"
            )
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        config["mirror"].write_bytes(path.read_bytes())
        print(f"{path}: {len(rows)} rows; flags added")


if __name__ == "__main__":
    main()
