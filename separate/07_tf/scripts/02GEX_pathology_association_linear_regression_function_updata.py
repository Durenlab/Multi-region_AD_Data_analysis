#!/usr/bin/env python3

import os
import glob
import numpy as np
import pandas as pd

from scipy.stats import mannwhitneyu
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests


# ============================================================
# Paths
# ============================================================

BASE_DIR = (os.environ["REGA_SEPARATE_ROOT"] + '')

TF_ACTIVITY_ROOT = os.environ.get(
    "REGA_TF_ACTIVITY_ROOT", os.path.join(BASE_DIR, "07_tf/results/activity")
)

META_DIR = os.environ.get(
    "REGA_METADATA_ROOT", os.path.join(BASE_DIR, "01_prepare/data/metadata_by_group")
)

OUT_DIR = os.environ.get(
    "REGA_TF_EXPRESSION_PATHOLOGY_ROOT", os.path.join(BASE_DIR, "07_tf/results/expression_pathology")
)

os.makedirs(OUT_DIR, exist_ok=True)

# ============================================================
# Settings
# ============================================================
TF_FILE_SUFFIX = "_TF_expression_pcc_cutoff.csv"

# Ordinal pathology score used as the predictor in linear regression.
PATHOLOGY_ORDER = {
    "nonAD": 0,
    "earlyAD": 1,
    "lateAD": 2,
}

MIN_NON_NA_SAMPLES = 10


# ============================================================
# Helper functions
# ============================================================

def parse_cell_region(label):
    """
    Parse label like:
        Ast_AG
        Mic_Immune_MTC

    Returns:
        celltype, region
    """
    regions = ["AG", "EC", "HC", "MTC", "PFC", "TH"]

    for region in regions:
        suffix = "_" + region
        if label.endswith(suffix):
            celltype = label[:-len(suffix)]
            return celltype, region

    raise ValueError(f"Cannot parse celltype/region from: {label}")


def safe_linear_regression(y, x):
    """
    Linear regression:
        y = intercept + beta * x + error

    Here:
        y = TF activity
        x = pathology score (nonAD=0, earlyAD=1, lateAD=2)

    Returns:
        beta, intercept, p-value for beta, R^2, n_used
    """
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)

    ok = np.isfinite(y) & np.isfinite(x)
    y = y[ok]
    x = x[ok]

    if len(y) < MIN_NON_NA_SAMPLES:
        return np.nan, np.nan, np.nan, np.nan, len(y)

    # Regression slope cannot be estimated if pathology score has no variation.
    # A constant TF activity vector is also not informative for association testing.
    if np.nanstd(x) == 0 or np.nanstd(y) == 0:
        return np.nan, np.nan, np.nan, np.nan, len(y)

    X = sm.add_constant(x)
    model = sm.OLS(y, X).fit()

    intercept = model.params[0]
    beta = model.params[1]
    pval = model.pvalues[1]
    r_squared = model.rsquared

    return beta, intercept, pval, r_squared, len(y)


def safe_mannwhitney(group1, group2):
    """
    Mann-Whitney U test.
    group1 and group2 are numeric arrays.
    """
    group1 = np.asarray(group1, dtype=float)
    group2 = np.asarray(group2, dtype=float)

    group1 = group1[np.isfinite(group1)]
    group2 = group2[np.isfinite(group2)]

    if len(group1) < 3 or len(group2) < 3:
        return np.nan, np.nan, len(group1), len(group2)

    if np.nanstd(np.concatenate([group1, group2])) == 0:
        return np.nan, np.nan, len(group1), len(group2)

    #stat, pval = mannwhitneyu(group1, group2, alternative="two-sided")
    stat, pval = mannwhitneyu(
    group1,
    group2,
    alternative="two-sided",
    method="asymptotic",
    use_continuity=True
)

    return stat, pval, len(group1), len(group2)


# ============================================================
# Main per celltype-region
# ============================================================

def run_one(label, tf_path, meta_path):
    print("\n" + "=" * 80)
    print(f"[RUN] {label}")
    print(f"[TF] {tf_path}")
    print(f"[META] {meta_path}")

    tf = pd.read_csv(tf_path, index_col=0)
    meta = pd.read_csv(meta_path, sep="\t")

    tf.index = tf.index.astype(str)
    tf.columns = tf.columns.astype(str)

    if "ROSMAP_IndividualID" not in meta.columns:
        raise ValueError(f"{meta_path} missing ROSMAP_IndividualID column")

    if "Pathology" not in meta.columns:
        raise ValueError(f"{meta_path} missing Pathology column")

    meta["ROSMAP_IndividualID"] = meta["ROSMAP_IndividualID"].astype(str)
    meta["Pathology"] = meta["Pathology"].astype(str)

    # Keep samples shared by TF activity and metadata.
    meta_sample_set = set(meta["ROSMAP_IndividualID"])
    common_samples = [
        s for s in tf.columns
        if s in meta_sample_set
    ]

    if len(common_samples) == 0:
        print(f"[WARN] {label}: no common samples, skip.")
        return None

    tf = tf.loc[:, common_samples]

    meta = (
        meta.set_index("ROSMAP_IndividualID")
            .loc[common_samples]
            .copy()
    )

    # Encode pathology as ordinal progression:
    # nonAD=0, earlyAD=1, lateAD=2.
    meta["Pathology_score"] = meta["Pathology"].map(PATHOLOGY_ORDER)

    valid = meta["Pathology_score"].notna()

    tf = tf.loc[:, valid.values]
    meta = meta.loc[valid].copy()

    print(f"[INFO] samples after alignment = {tf.shape[1]}")
    print("[INFO] pathology counts:")
    print(meta["Pathology"].value_counts().to_string())

    if tf.shape[1] < MIN_NON_NA_SAMPLES:
        print(f"[WARN] {label}: too few samples, skip.")
        return None

    pathology_score = meta["Pathology_score"].to_numpy(dtype=float)

    rows = []

    for tf_name in tf.index:
        values = pd.to_numeric(
            tf.loc[tf_name], errors="coerce"
        ).to_numpy(dtype=float)

        # Linear regression:
        # TF_activity ~ Pathology_score
        # Pathology_score: nonAD=0, earlyAD=1, lateAD=2
        beta, intercept, lm_pval, r_squared, n_used = safe_linear_regression(
            values,
            pathology_score
        )

        # Group means
        group_means = {}
        for group in ["nonAD", "earlyAD", "lateAD"]:
            idx = meta["Pathology"].values == group
            group_values = values[idx]
            finite_group_values = group_values[np.isfinite(group_values)]

            if len(finite_group_values) > 0:
                group_means[f"mean_{group}"] = np.mean(finite_group_values)
            else:
                group_means[f"mean_{group}"] = np.nan

        # lateAD vs nonAD Mann-Whitney U test
        late_vals = values[meta["Pathology"].values == "lateAD"]
        non_vals = values[meta["Pathology"].values == "nonAD"]

        mw_stat, mw_p, n_late, n_non = safe_mannwhitney(late_vals, non_vals)

        rows.append({
            "TF": tf_name,
            "n_samples": n_used,
            "lm_beta_pathology": beta,
            "lm_intercept_pathology": intercept,
            "lm_pval_pathology": lm_pval,
            "lm_r2_pathology": r_squared,
            "mean_nonAD": group_means["mean_nonAD"],
            "mean_earlyAD": group_means["mean_earlyAD"],
            "mean_lateAD": group_means["mean_lateAD"],
            "lateAD_minus_nonAD": (
                group_means["mean_lateAD"] - group_means["mean_nonAD"]
            ),
            "mannwhitney_lateAD_vs_nonAD_stat": mw_stat,
            "mannwhitney_lateAD_vs_nonAD_pval": mw_p,
            "n_lateAD": n_late,
            "n_nonAD": n_non,
        })

    res = pd.DataFrame(rows)

    # FDR correction for linear-regression pathology association
    mask = res["lm_pval_pathology"].notna()
    res["lm_fdr_pathology"] = np.nan
    if mask.sum() > 0:
        res.loc[mask, "lm_fdr_pathology"] = multipletests(
            res.loc[mask, "lm_pval_pathology"],
            method="fdr_bh"
        )[1]

    # FDR correction for lateAD vs nonAD Mann-Whitney test
    mask = res["mannwhitney_lateAD_vs_nonAD_pval"].notna()
    res["mannwhitney_lateAD_vs_nonAD_fdr"] = np.nan
    if mask.sum() > 0:
        res.loc[mask, "mannwhitney_lateAD_vs_nonAD_fdr"] = multipletests(
            res.loc[mask, "mannwhitney_lateAD_vs_nonAD_pval"],
            method="fdr_bh"
        )[1]

    # Add label columns
    celltype, region = parse_cell_region(label)

    res.insert(0, "celltype_region", label)
    res.insert(1, "celltype", celltype)
    res.insert(2, "region", region)

    # Sort: smallest regression p-value first;
    # among ties, strongest positive beta first.
    res = res.sort_values(
        ["lm_pval_pathology", "lm_beta_pathology"],
        ascending=[True, False],
        na_position="last"
    )

    out_cell_dir = os.path.join(OUT_DIR, label)
    os.makedirs(out_cell_dir, exist_ok=True)

    out_path = os.path.join(
        out_cell_dir,
        f"{label}_TF_pathology_association.tsv"
    )
    res.to_csv(out_path, sep="\t", index=False)

    print(f"[SAVE] {out_path}")
    print("[TOP positive beta]")
    print(
        res.sort_values("lm_beta_pathology", ascending=False)
           .head(10)[
               [
                   "TF",
                   "lm_beta_pathology",
                   "lm_pval_pathology",
                   "lm_fdr_pathology",
                   "lm_r2_pathology",
                   "mean_nonAD",
                   "mean_earlyAD",
                   "mean_lateAD",
                   "lateAD_minus_nonAD",
               ]
           ]
           .to_string(index=False)
    )

    return res


# ============================================================
# Main
# ============================================================

def main():
    all_results = []

    tf_paths = sorted(glob.glob(
        os.path.join(
        TF_ACTIVITY_ROOT,
        f"*{TF_FILE_SUFFIX}"
        )
    ))

    print(f"[INFO] TF activity files found = {len(tf_paths)}")

    if len(tf_paths) == 0:
        raise SystemExit("No TF activity files found. Check TF_ACTIVITY_ROOT.")

    for tf_path in tf_paths:

        filename = os.path.basename(tf_path)

        if not filename.endswith(TF_FILE_SUFFIX):
            print(f"[WARN] unexpected filename, skip: {filename}")
            continue

        # Example:
        # Ast_AG_TF_activity_pcc_cutoff.csv -> Ast_AG
        label = filename[:-len(TF_FILE_SUFFIX)]
        if label.startswith("Vasc_"):
            continue

        meta_path = os.path.join(
            META_DIR,
            f"{label}_metadata.tsv"
        )

        if not os.path.exists(meta_path):
            print(f"[WARN] missing metadata for {label}: {meta_path}")
            continue

        try:
            res = run_one(label, tf_path, meta_path)

            if res is not None:
                all_results.append(res)

        except Exception as e:
            print(f"[ERROR] {label}: {e}")

    if len(all_results) == 0:
        raise SystemExit("No successful results.")

    all_df = pd.concat(all_results, axis=0, ignore_index=True)

    out_all = os.path.join(
        OUT_DIR,
        "all_celltype_region_TF_pathology_association.tsv"
    )
    all_df.to_csv(out_all, sep="\t", index=False)

    print("\n" + "=" * 80)
    print(f"[SAVE] combined result: {out_all}")

    # Save results ranked toward positive pathology association.
    top_pos = all_df.sort_values(
        ["lm_fdr_pathology", "lm_beta_pathology"],
        ascending=[True, False],
        na_position="last"
    )

    out_top_pos = os.path.join(
        OUT_DIR,
        "top_positive_TF_pathology_association.tsv"
    )
    top_pos.to_csv(out_top_pos, sep="\t", index=False)

    # Save results ranked toward negative pathology association.
    top_neg = all_df.sort_values(
        ["lm_fdr_pathology", "lm_beta_pathology"],
        ascending=[True, True],
        na_position="last"
    )

    out_top_neg = os.path.join(
        OUT_DIR,
        "top_negative_TF_pathology_association.tsv"
    )
    top_neg.to_csv(out_top_neg, sep="\t", index=False)

    print(f"[SAVE] top positive: {out_top_pos}")
    print(f"[SAVE] top negative: {out_top_neg}")


if __name__ == "__main__":
    main()
