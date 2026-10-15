#!/usr/bin/env python
from pathlib import Path
import os
SHARED_ROOT = Path(os.environ["REGA_SHARED_ROOT"]) if "REGA_SHARED_ROOT" in os.environ else next(
    p for p in Path(__file__).resolve().parents if p.name in {"REGA_shared", "shared"}
)
SEPARATE_ROOT = Path(os.environ.get("REGA_SEPARATE_ROOT", "/project/zduren/durenlab/naqing_palmetto/REGA_separate"))

import os
import re
import pandas as pd
from statsmodels.stats.multitest import multipletests

ldsc_dir = str(SHARED_ROOT / 'part4_ldsc')
ctm_file = os.path.join(ldsc_dir, "selected_modules_FC.txt")

result_root = os.path.join(
    ldsc_dir,
    "heritability_results"
)

annot_root = os.path.join(
    ldsc_dir,
    "annotations_phase3"
)

out_all = os.path.join(
    result_root,
    "AD_Bellenguez_single_module_extend2kb_LDSC_results.tsv"
)

out_qc = os.path.join(
    result_root,
    "AD_Bellenguez_single_module_extend2kb_LDSC_QC.tsv"
)

def parse_module_id(raw):
    raw = str(raw).strip().split()[-1]

    if not re.search(r"_M\d+$", raw):
        raise ValueError(f"Cannot parse module ID: {raw}")

    module = re.search(r"M\d+$", raw).group(0)
    celltype_region = re.sub(r"_M\d+$", "", raw)
    return raw, celltype_region, module


def get_total_annotated_snps(module_id):
    total = 0
    missing = []

    for chrom in range(1, 23):
        f = os.path.join(
            annot_root,
            module_id,
            f"Base_CelltypePeakLD.{chrom}.annot"
        )

        if not os.path.exists(f):
            missing.append(f)
            continue

        with open(f) as handle:
            next(handle)
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 5:
                    try:
                        total += float(parts[4])
                    except ValueError:
                        pass

    return int(total), len(missing)


def get_m_5_50(module_id):
    total = 0
    missing = []

    for chrom in range(1, 23):
        f = os.path.join(
            ldsc_dir,
            "ld_scores_phase3",
            module_id,
            f"Base_CelltypePeakLD.{chrom}.l2.M_5_50"
        )

        if not os.path.exists(f):
            missing.append(f)
            continue

        with open(f) as handle:
            x = handle.read().strip().split()
            if len(x) > 0:
                total += float(x[0])

    return int(total), len(missing)


def find_module_row(df):
    if "Category" not in df.columns:
        return None

    cat = df["Category"].astype(str)

    # single custom annotation is usually L2_1
    sub = df.loc[cat.eq("L2_1")].copy()
    if sub.shape[0] > 0:
        return sub

    # fallback: any *_1 row
    sub = df.loc[cat.str.endswith("_1")].copy()
    if sub.shape[0] == 1:
        return sub

    return None


with open(ctm_file) as f:
    raw_modules = [x.strip() for x in f if x.strip() and not x.startswith("#")]

records = []
qc_records = []

for raw in raw_modules:
    module_id, celltype_region, module = parse_module_id(raw)

    result_file = os.path.join(
        result_root,
        module_id,
        f"AD_Bellenguez_extend2kb_{module_id}.results"
    )

    log_file = os.path.join(
        result_root,
        module_id,
        f"AD_Bellenguez_extend2kb_{module_id}.log"
    )

    total_snps, n_missing_annot = get_total_annotated_snps(module_id)
    m550, n_missing_m550 = get_m_5_50(module_id)

    qc = {
        "module_id": module_id,
        "celltype_region": celltype_region,
        "module": module,
        "result_file": result_file,
        "log_file": log_file,
        "result_exists": os.path.exists(result_file),
        "log_exists": os.path.exists(log_file),
        "total_annotated_SNPs": total_snps,
        "M_5_50": m550,
        "missing_annotation_chr_files": n_missing_annot,
        "missing_M_5_50_chr_files": n_missing_m550,
        "n_result_rows_matched": 0,
        "matched_category": "",
    }

    if not os.path.exists(result_file):
        qc_records.append(qc)
        print(f"WARNING missing result file: {result_file}")
        continue

    try:
        df = pd.read_csv(result_file, sep=r"\s+")
    except Exception as e:
        qc["read_error"] = str(e)
        qc_records.append(qc)
        print(f"WARNING could not read {result_file}: {e}")
        continue

    sub = find_module_row(df)

    if sub is None:
        qc_records.append(qc)
        print(f"WARNING no module row matched in {result_file}")
        print(df["Category"].astype(str).tolist())
        continue

    qc["n_result_rows_matched"] = sub.shape[0]
    qc["matched_category"] = ";".join(sub["Category"].astype(str).tolist())
    qc_records.append(qc)

    sub["original_Category_from_LDSC"] = sub["Category"]
    sub["Category"] = module_id
    sub["module_id"] = module_id
    sub["celltype_region"] = celltype_region
    sub["module"] = module
    sub["total_annotated_SNPs"] = total_snps
    sub["M_5_50"] = m550
    sub["result_file"] = result_file
    sub["log_file"] = log_file

    records.append(sub)

qc_df = pd.DataFrame(qc_records)
qc_df.to_csv(out_qc, sep="\t", index=False)

if len(records) == 0:
    raise SystemExit("No LDSC result rows collected. Check QC file.")

res = pd.concat(records, ignore_index=True)

metadata_cols = {
    "Category",
    "original_Category_from_LDSC",
    "module_id",
    "celltype_region",
    "module",
    "result_file",
    "log_file",
}

for col in res.columns:
    if col not in metadata_cols:
        res[col] = pd.to_numeric(res[col], errors="ignore")

if "Enrichment_p" in res.columns:
    p = pd.to_numeric(res["Enrichment_p"], errors="coerce")
    valid = p.notna()
    res["FDR_BH_Enrichment_p"] = pd.NA
    res.loc[valid, "FDR_BH_Enrichment_p"] = multipletests(p[valid], method="fdr_bh")[1]

if "Coefficient_P_value" in res.columns:
    p = pd.to_numeric(res["Coefficient_P_value"], errors="coerce")
    valid = p.notna()
    res["FDR_BH_Coefficient_P_value"] = pd.NA
    res.loc[valid, "FDR_BH_Coefficient_P_value"] = multipletests(p[valid], method="fdr_bh")[1]

if "Enrichment_p" in res.columns:
    res = res.sort_values("Enrichment_p")

res.to_csv(out_all, sep="\t", index=False)

print("==================================================")
print("Done collecting extend2kb single-module LDSC results")
print("==================================================")
print(f"Input modules: {len(raw_modules)}")
print(f"Collected rows: {res.shape[0]}")
print(f"Wrote result table: {out_all}")
print(f"Wrote QC table: {out_qc}")

cols_to_show = [
    "module_id",
    "original_Category_from_LDSC",
    "celltype_region",
    "module",
    "total_annotated_SNPs",
    "M_5_50",
]

for c in [
    "Prop._SNPs",
    "Prop._h2",
    "Enrichment",
    "Enrichment_p",
    "Coefficient",
    "Coefficient_z-score",
    "Coefficient_P_value",
    "FDR_BH_Enrichment_p",
    "FDR_BH_Coefficient_P_value",
]:
    if c in res.columns:
        cols_to_show.append(c)

print(res[cols_to_show].head(30).to_string(index=False))