#!/usr/bin/env python3
"""Select modules with region-level BH FDR < 0.05 in either prespecified test."""

import argparse
import csv
import math
from pathlib import Path


REGIONS = ("AG", "EC", "HC", "MTC", "PFC", "TH")
TABLE = "ALL_regions_sharedH_All_four_tests_H_pathology_results.tsv"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path(__file__).resolve().parent / "08pathology_tests" / TABLE)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent / "08pathology_tests" / "selected_module.txt")
    parser.add_argument("--fdr", type=float, default=0.05)
    args = parser.parse_args()
    selected = set()
    counts = {"Linear_Regression": 0, "Pairwise_MannWhitneyU": 0}
    with args.input.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"Region", "Module", "test", "Comparison", "pval_adj"}
        if not required <= set(reader.fieldnames or ()):
            raise ValueError(f"Missing columns: {sorted(required - set(reader.fieldnames or ()))}")
        for row in reader:
            test = row["test"]
            if test == "Linear_Regression":
                pass
            elif test == "Pairwise_MannWhitneyU" and row["Comparison"] == "nonAD vs lateAD":
                pass
            else:
                continue
            if row["Region"] not in REGIONS:
                raise ValueError(f"Unexpected region: {row['Region']}")
            fdr = float(row["pval_adj"]) if row["pval_adj"] not in ("", "NA", "NaN") else math.nan
            if math.isfinite(fdr) and fdr < args.fdr:
                selected.add((row["Region"], row["Module"]))
                counts[test] += 1
    ordered = sorted(selected, key=lambda item: (REGIONS.index(item[0]), int(item[1][1:])))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("".join(f"{region}_{module}\n" for region, module in ordered))
    print(f"Selected {len(ordered)} unique modules; contributing rows: {counts}")


if __name__ == "__main__":
    main()
