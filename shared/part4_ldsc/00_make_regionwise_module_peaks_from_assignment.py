#!/usr/bin/env python3
from pathlib import Path
import os
SHARED_ROOT = Path(os.environ["REGA_SHARED_ROOT"]) if "REGA_SHARED_ROOT" in os.environ else next(
    p for p in Path(__file__).resolve().parents if p.name in {"REGA_shared", "shared"}
)
SEPARATE_ROOT = Path(os.environ.get("REGA_SEPARATE_ROOT", "/project/zduren/durenlab/naqing_palmetto/REGA_separate"))

import os
import re
import argparse
import pandas as pd


def parse_re_to_bed(re_name):
    """
    Convert RE string like:
        chr19:43316629-43317129
    to:
        chr19, 43316629, 43317129, chr19:43316629-43317129
    """
    m = re.match(r"^(chr[0-9XYM]+):([0-9]+)-([0-9]+)$", str(re_name))
    if m is None:
        return None

    chrom = m.group(1)
    start = int(m.group(2))
    end = int(m.group(3))

    if end <= start:
        return None

    return chrom, start, end, f"{chrom}:{start}-{end}"


def chrom_sort_key(chrom):
    """
    Sort chromosomes in natural order:
    chr1, chr2, ..., chr22, chrX, chrY, chrM
    """
    chrom = str(chrom).replace("chr", "")

    if chrom.isdigit():
        return int(chrom)

    if chrom == "X":
        return 23
    if chrom == "Y":
        return 24
    if chrom in ["M", "MT"]:
        return 25

    return 99


def process_one_region(region, input_file, out_root):
    print(f"\n[INFO] Processing region: {region}")
    print(f"[INFO] Input: {input_file}")

    if not os.path.exists(input_file):
        print(f"[WARN] Missing input file, skip: {input_file}")
        return

    df = pd.read_csv(input_file, sep="\t")

    required_cols = {"region", "celltype_region", "module", "RE"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in {input_file}: {missing}")

    df = df[["region", "celltype_region", "module", "RE"]].copy()

    # Keep only this region, just in case
    df = df[df["region"] == region].copy()

    if df.empty:
        print(f"[WARN] No rows for region {region}, skip.")
        return

    # Parse RE coordinates
    parsed = df["RE"].apply(parse_re_to_bed)
    bad = parsed.isna().sum()

    if bad > 0:
        print(f"[WARN] {bad} rows have invalid RE format and will be dropped.")

    df = df[parsed.notna()].copy()
    parsed = parsed[parsed.notna()]

    df[["chr", "start", "end", "name"]] = pd.DataFrame(
        parsed.tolist(),
        index=df.index
    )

    # Remove exact duplicate peak assignments within the same celltype-region/module
    before = df.shape[0]
    df = df.drop_duplicates(
        subset=["celltype_region", "module", "chr", "start", "end", "name"]
    ).copy()
    after = df.shape[0]

    if before != after:
        print(f"[INFO] Removed duplicated rows: {before - after}")

    # Write one BED file per celltype_region x module
    summary_rows = []

    for (celltype_region, module), sub in df.groupby(["celltype_region", "module"]):
        out_dir = os.path.join(out_root, region, celltype_region, "module_peaks")
        os.makedirs(out_dir, exist_ok=True)

        out_bed = os.path.join(out_dir, f"{celltype_region}_{module}.bed")

        sub = sub[["chr", "start", "end", "name"]].copy()
        sub["chr_order"] = sub["chr"].map(chrom_sort_key)
        sub = sub.sort_values(["chr_order", "start", "end", "name"])
        sub = sub[["chr", "start", "end", "name"]]

        sub.to_csv(out_bed, sep="\t", header=False, index=False)

        summary_rows.append({
            "region": region,
            "celltype_region": celltype_region,
            "module": module,
            "n_peaks": sub.shape[0],
            "out_bed": out_bed
        })

    summary = pd.DataFrame(summary_rows)
    summary_file = os.path.join(out_root, region, f"{region}_module_peaks_summary.tsv")
    os.makedirs(os.path.dirname(summary_file), exist_ok=True)
    summary.to_csv(summary_file, sep="\t", index=False)

    print(f"[INFO] Wrote {summary.shape[0]} module peak BED files.")
    print(f"[INFO] Summary: {summary_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Make region-wise shared-H module peak BED files from P RE-module assignment table."
    )

    parser.add_argument(
        "--base-dir",
        default=str(SHARED_ROOT),
        help="Base directory for the region-shared REGA analysis."
    )

    parser.add_argument(
        "--assignment-subdir",
        default="part3_process_REGA_output/10-2RE_module_assignment_region_shareH_FC",
        help="Subdirectory containing region-wise FC RE assignment files."
    )

    parser.add_argument(
        "--out-root",
        default=None,
        help="Output root directory. If not provided, use BASE_DIR/part4_ldsc/module_peaks_hg38."
    )

    parser.add_argument(
        "--regions",
        nargs="+",
        default=None,
        help="Regions to process. If not provided, infer from folders under assignment-subdir."
    )

    args = parser.parse_args()

    base_dir = args.base_dir
    input_root = os.path.join(base_dir, args.assignment_subdir)

    if args.out_root is None:
        out_root = os.path.join(
            base_dir,
            "part4_ldsc",
            "module_peaks_hg38"
        )
    else:
        out_root = args.out_root

    if args.regions is None:
        regions = sorted([
            x for x in os.listdir(input_root)
            if os.path.isdir(os.path.join(input_root, x))
        ])
    else:
        regions = args.regions

    print("[INFO] Input root:", input_root)
    print("[INFO] Output root:", out_root)
    print("[INFO] Regions:", ", ".join(regions))

    os.makedirs(out_root, exist_ok=True)

    all_summary = []

    for region in regions:
        input_file = os.path.join(
            input_root,
            region,
            f"{region}_ALL_celltype_region_RE_module_assignments_FC.tsv"
        )

        process_one_region(region, input_file, out_root)

        summary_file = os.path.join(out_root, region, f"{region}_module_peaks_summary.tsv")
        if os.path.exists(summary_file):
            all_summary.append(pd.read_csv(summary_file, sep="\t"))

    if all_summary:
        all_summary = pd.concat(all_summary, axis=0, ignore_index=True)
        all_summary_file = os.path.join(out_root, "ALL_regions_module_peaks_summary.tsv")
        all_summary.to_csv(all_summary_file, sep="\t", index=False)
        print("\n[INFO] Wrote all-region summary:")
        print(all_summary_file)


if __name__ == "__main__":
    main()
