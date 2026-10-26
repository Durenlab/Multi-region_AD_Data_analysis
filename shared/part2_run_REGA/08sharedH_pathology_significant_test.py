#!/usr/bin/env python
from pathlib import Path
import os
import argparse
SHARED_ROOT = Path(os.environ["REGA_SHARED_ROOT"]) if "REGA_SHARED_ROOT" in os.environ else next(
    p for p in Path(__file__).resolve().parents if p.name in {"REGA_shared", "shared"}
)
SEPARATE_ROOT = Path(os.environ["REGA_SEPARATE_ROOT"])

import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import kruskal, mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm


# ============================================================
# Paths
# ============================================================

BASE_DIR = str(SHARED_ROOT)

# Region-wise shared-H REGA output directory
INPUT_ROOT = os.path.join(
    BASE_DIR,
    "part2_run_REGA/07REGA_region_shareH_output"
)

# Metadata generated from region-consistent RNA matrices
# Expected files:
#   Ast_AG_region_consistent_metadata.tsv
#   Exc_AG_region_consistent_metadata.tsv
#   ...
META_DIR = os.path.join(
    BASE_DIR,
    "part1_prepare_REGA/01region_consistent_metadata"
)

OUT_ROOT = os.path.join(
    BASE_DIR,
    "part2_run_REGA/08pathology_tests"
)

os.makedirs(OUT_ROOT, exist_ok=True)


# ============================================================
# Settings
# ============================================================

regions = ["PFC", "TH", "AG", "EC", "HC", "MTC"]

# Use one representative metadata file per region.
# Since H is shared across cell types within a region, sample list should be the same
# across celltype-region metadata files for that region.
representative_celltype = "Ast"

pathology_order = ["nonAD", "earlyAD", "lateAD"]

pathology_mapping = {
    "nonAD": 0,
    "earlyAD": 1,
    "lateAD": 2
}

pairwise_comparisons = [
    ("nonAD", "earlyAD"),
    ("nonAD", "lateAD"),
    ("earlyAD", "lateAD"),
]


# ============================================================
# Helper functions
# ============================================================

def add_fdr(df, p_col="pval", out_col="pval_adj"):
    """
    Add Benjamini-Hochberg FDR correction.
    NaN p-values are ignored.
    """
    df = df.copy()
    df[out_col] = np.nan

    valid = df[p_col].notna()

    if valid.sum() > 0:
        df.loc[valid, out_col] = multipletests(
            df.loc[valid, p_col],
            method="fdr_bh"
        )[1]

    return df


def p_to_star(p):
    """
    Convert p-value or adjusted p-value to significance label.
    """
    if pd.isna(p):
        return ""
    elif p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return ""


def find_sample_column(meta_df, h_sample_names):
    """
    Try to identify which metadata column matches H sample names.
    """
    candidate_cols = [
        "ROSMAP_IndividualID",
        "sample_id",
        "Sample",
        "SampleID",
        "sample",
        "donor_id",
        "Donor",
        "individual",
        "Individual",
        "IndividualID",
        "RNA_ID",
        "ID"
    ]

    h_set = set(map(str, h_sample_names))

    best_col = None
    best_overlap = 0

    for col in candidate_cols:
        if col not in meta_df.columns:
            continue

        overlap = len(h_set.intersection(set(meta_df[col].astype(str))))

        if overlap > best_overlap:
            best_col = col
            best_overlap = overlap

    return best_col, best_overlap


def load_region_metadata(region, sample_names):
    """
    Load one representative metadata file for this region and align to H samples.

    Because region-wise REGA shares H across cell types within the same region,
    all celltype-region matrices in the same region should have the same individuals.
    Therefore, one metadata file, for example Ast_REGION, is enough for sample-level
    pathology labels.
    """

    meta_path = os.path.join(
        META_DIR,
        f"{representative_celltype}_{region}_region_consistent_metadata.tsv"
    )

    if not os.path.exists(meta_path):
        raise FileNotFoundError(f"Missing metadata file: {meta_path}")

    meta = pd.read_csv(meta_path, sep="\t")

    print("\nLoaded metadata:")
    print(meta_path)
    print(meta.shape)
    print(meta.head())

    if "Pathology" not in meta.columns:
        raise KeyError("Metadata must contain a column named 'Pathology'.")

    sample_col, overlap = find_sample_column(meta, sample_names)

    if sample_col is not None and overlap == len(sample_names):
        print(f"Matching metadata using column: {sample_col}")

        meta[sample_col] = meta[sample_col].astype(str)

        meta = (
            meta
            .set_index(sample_col)
            .loc[sample_names]
            .reset_index()
        )

    elif meta.shape[0] == len(sample_names):
        print(
            "\nWARNING: Could not find a metadata sample column that fully matches H columns.\n"
            "Using metadata row order. Make sure metadata rows are in the exact same order as H columns."
        )

    else:
        raise ValueError(
            f"Cannot align metadata to H for region {region}.\n"
            f"H has {len(sample_names)} samples, but metadata has {meta.shape[0]} rows.\n"
            f"Best matching metadata column: {sample_col}, overlap: {overlap}."
        )

    return meta, meta_path


def analyze_one_region(region):
    """
    Run all H-pathology tests for one brain region.
    """

    print("\n============================================================")
    print(f"Region: {region}")
    print("============================================================")

    input_dir = os.path.join(
        INPUT_ROOT,
        f"AD_region_shareH_{region}_K20"
    )

    h_path = os.path.join(input_dir, "H_with_samples.csv")

    out_dir = os.path.join(OUT_ROOT, region)
    boxplot_dir = os.path.join(out_dir, "boxplot")

    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(boxplot_dir, exist_ok=True)

    if not os.path.exists(h_path):
        print(f"WARNING: missing H file, skipping region {region}: {h_path}")
        return None

    # ============================================================
    # Load H matrix
    # ============================================================

    H_df = pd.read_csv(h_path, index_col=0)

    print("Loaded H matrix:")
    print(H_df.shape)
    print(H_df.iloc[:5, :5])

    # Expected: modules x individuals.
    # If matrix appears samples x modules, transpose.
    n_rows, n_cols = H_df.shape

    if n_rows > n_cols:
        print(
            "\nH appears to have more rows than columns. "
            "Assuming rows are samples and columns are modules; transposing."
        )
        H_df = H_df.T

    H_df.index = [str(x) for x in H_df.index]

    old_module_names = H_df.index.tolist()
    module_order = [f"M{i + 1}" for i in range(H_df.shape[0])]

    module_name_map = dict(zip(old_module_names, module_order))
    H_df.index = module_order

    n_modules, n_samples = H_df.shape
    sample_names = H_df.columns.astype(str).tolist()

    print(f"\nFinal H orientation: {n_modules} modules x {n_samples} samples")
    print("Modules:", module_order)
    print("First 5 samples:", sample_names[:5])

    # ============================================================
    # Load and align metadata
    # ============================================================

    meta, meta_path = load_region_metadata(region, sample_names)

    # ============================================================
    # Prepare pathology labels
    # ============================================================

    pathology = meta["Pathology"].astype(str).str.strip()
    pathology_num = pathology.map(pathology_mapping)

    if pathology_num.isna().any():
        bad_values = sorted(pathology[pathology_num.isna()].unique())
        raise ValueError(
            f"Some Pathology values are not recognized in region {region}. "
            f"Expected {list(pathology_mapping.keys())}, got problematic values: {bad_values}"
        )

    n_counts = pathology.value_counts().to_dict()
    n_nonAD = n_counts.get("nonAD", 0)
    n_earlyAD = n_counts.get("earlyAD", 0)
    n_lateAD = n_counts.get("lateAD", 0)

    print("\nPathology counts:")
    print(f"nonAD   : {n_nonAD}")
    print(f"earlyAD : {n_earlyAD}")
    print(f"lateAD  : {n_lateAD}")

    # ============================================================
    # Reshape H to long format
    # ============================================================

    df_list = []

    for module in module_order:
        tmp = pd.DataFrame({
            "Region": region,
            "Module": module,
            "Value": H_df.loc[module, :].values,
            "sample_id": sample_names,
            "Pathology": pathology.values,
            "Pathology_num": pathology_num.values
        })

        for col in ["celltype", "region", "BrainRegion", "Sex", "age_death", "pmi", "Study"]:
            if col in meta.columns:
                tmp[col] = meta[col].values

        df_list.append(tmp)

    df = pd.concat(df_list, ignore_index=True)

    df["Pathology"] = pd.Categorical(
        df["Pathology"],
        categories=pathology_order,
        ordered=True
    )

    # ============================================================
    # Store results
    # ============================================================

    kw_all = []
    rs_all = []
    lm_all = []
    spearman_all = []

    # ============================================================
    # 1. Kruskal-Wallis test
    # ============================================================

    for m in module_order:

        sub = df[df["Module"] == m]

        groups_for_test = [
            sub[sub["Pathology"] == g]["Value"].dropna().values
            for g in pathology_order
        ]

        if all(len(x) > 0 for x in groups_for_test):
            stat, pval = kruskal(*groups_for_test)
        else:
            stat, pval = np.nan, np.nan

        kw_all.append({
            "Region": region,
            "test": "Kruskal_Wallis",
            "Module": m,
            "Comparison": "nonAD vs earlyAD vs lateAD",
            "statistic": stat,
            "effect_size": np.nan,
            "pval": pval,
            "n": sub["Value"].dropna().shape[0],
            "n_nonAD": n_nonAD,
            "n_earlyAD": n_earlyAD,
            "n_lateAD": n_lateAD
        })

    # ============================================================
    # 2. Pairwise Mann-Whitney U tests
    # ============================================================

    for m in module_order:

        sub = df[df["Module"] == m]

        for g1, g2 in pairwise_comparisons:

            x = sub[sub["Pathology"] == g1]["Value"].dropna().values
            y = sub[sub["Pathology"] == g2]["Value"].dropna().values

            if len(x) > 0 and len(y) > 0:
                stat, pval = mannwhitneyu(
                    x,
                    y,
                    alternative="two-sided"
                )
            else:
                stat, pval = np.nan, np.nan

            rs_all.append({
                "Region": region,
                "test": "Pairwise_MannWhitneyU",
                "Module": m,
                "Comparison": f"{g1} vs {g2}",
                "statistic": stat,
                "effect_size": np.nan,
                "pval": pval,
                "n1": len(x),
                "n2": len(y),
                "n_nonAD": n_nonAD,
                "n_earlyAD": n_earlyAD,
                "n_lateAD": n_lateAD
            })

    # ============================================================
    # 3. Linear regression: H value ~ pathology stage
    # ============================================================

    for m in module_order:

        sub = df[df["Module"] == m].copy()
        sub = sub.dropna(subset=["Value", "Pathology_num"])

        if sub.shape[0] >= 3 and sub["Pathology_num"].nunique() >= 2:

            X = sub["Pathology_num"].values
            y = sub["Value"].values

            X = sm.add_constant(X)

            model = sm.OLS(y, X).fit()

            slope = model.params[1]
            stat = model.tvalues[1]
            pval = model.pvalues[1]
            n = sub.shape[0]

        else:
            slope = np.nan
            stat = np.nan
            pval = np.nan
            n = sub.shape[0]

        lm_all.append({
            "Region": region,
            "test": "Linear_Regression",
            "Module": m,
            "Comparison": "H module value ~ pathology stage",
            "Pathology_encoding": "nonAD=0, earlyAD=1, lateAD=2",
            "statistic": stat,
            "effect_size": slope,
            "slope": slope,
            "pval": pval,
            "n": n,
            "n_nonAD": n_nonAD,
            "n_earlyAD": n_earlyAD,
            "n_lateAD": n_lateAD
        })

    # ============================================================
    # 4. Spearman correlation: H value vs pathology stage
    # ============================================================

    for m in module_order:

        sub = df[df["Module"] == m].copy()
        sub = sub.dropna(subset=["Value", "Pathology_num"])

        if sub.shape[0] >= 3 and sub["Pathology_num"].nunique() >= 2:

            rho, pval = spearmanr(
                sub["Pathology_num"].values,
                sub["Value"].values
            )

            n = sub.shape[0]
            n_pathology_groups = sub["Pathology_num"].nunique()

        else:
            rho = np.nan
            pval = np.nan
            n = sub.shape[0]
            n_pathology_groups = sub["Pathology_num"].nunique()

        spearman_all.append({
            "Region": region,
            "test": "Spearman",
            "Module": m,
            "Comparison": "Spearman: H module value vs pathology stage",
            "Pathology_encoding": "nonAD=0, earlyAD=1, lateAD=2",
            "statistic": rho,
            "effect_size": rho,
            "rho": rho,
            "pval": pval,
            "n": n,
            "n_pathology_groups": n_pathology_groups,
            "n_nonAD": n_nonAD,
            "n_earlyAD": n_earlyAD,
            "n_lateAD": n_lateAD
        })

    # ============================================================
    # Combine result tables
    # ============================================================

    kw_df = pd.DataFrame(kw_all)
    rs_df = pd.DataFrame(rs_all)
    lm_df = pd.DataFrame(lm_all)
    spearman_df = pd.DataFrame(spearman_all)

    # ============================================================
    # FDR correction within this region
    # ============================================================

    kw_df = add_fdr(kw_df, p_col="pval", out_col="pval_adj")
    lm_df = add_fdr(lm_df, p_col="pval", out_col="pval_adj")
    spearman_df = add_fdr(spearman_df, p_col="pval", out_col="pval_adj")

    rs_df = (
        rs_df
        .groupby("Comparison", group_keys=False)
        .apply(lambda x: add_fdr(x, p_col="pval", out_col="pval_adj"))
    )

    # Optional global FDR columns within this region/test table
    kw_df = add_fdr(kw_df, p_col="pval", out_col="pval_adj_global")
    rs_df = add_fdr(rs_df, p_col="pval", out_col="pval_adj_global")
    lm_df = add_fdr(lm_df, p_col="pval", out_col="pval_adj_global")
    spearman_df = add_fdr(spearman_df, p_col="pval", out_col="pval_adj_global")

    # ============================================================
    # Sort result tables
    # ============================================================

    kw_df = kw_df.sort_values(["pval_adj", "pval"], na_position="last")
    rs_df = rs_df.sort_values(["Comparison", "pval_adj", "pval"], na_position="last")
    lm_df = lm_df.sort_values(["pval_adj", "pval"], na_position="last")
    spearman_df = spearman_df.sort_values(["pval_adj", "pval"], na_position="last")

    # ============================================================
    # Save individual result tables
    # ============================================================

    kw_out = os.path.join(out_dir, f"{region}_sharedH_Kruskal_Wallis.tsv")
    rs_out = os.path.join(out_dir, f"{region}_sharedH_Pairwise_MannWhitneyU.tsv")
    lm_out = os.path.join(out_dir, f"{region}_sharedH_LinearRegression.tsv")
    spearman_out = os.path.join(out_dir, f"{region}_sharedH_Spearman_pathology_H_modules.tsv")

    kw_df.to_csv(kw_out, sep="\t", index=False)
    rs_df.to_csv(rs_out, sep="\t", index=False)
    lm_df.to_csv(lm_out, sep="\t", index=False)
    spearman_df.to_csv(spearman_out, sep="\t", index=False)

    # ============================================================
    # Save combined result table
    # ============================================================

    combined_df = pd.concat(
        [kw_df, rs_df, lm_df, spearman_df],
        ignore_index=True,
        sort=False
    )

    combined_df = combined_df.sort_values(
        ["Region", "Module", "test", "Comparison"],
        na_position="last"
    )

    combined_out = os.path.join(
        out_dir,
        f"{region}_sharedH_All_four_tests_H_pathology_results.tsv"
    )

    combined_df.to_csv(combined_out, sep="\t", index=False)

    # ============================================================
    # Select significant shared H modules
    # ============================================================

    sig_df = combined_df[
        combined_df["pval_adj"].notna() &
        (combined_df["pval_adj"] < 0.05)
    ].copy()

    selected_modules = (
        sig_df["Module"]
        .dropna()
        .drop_duplicates()
        .sort_values()
    )

    selected_module_out = os.path.join(
        out_dir,
        f"{region}_sharedH_selected_module.txt"
    )

    selected_modules.to_csv(
        selected_module_out,
        sep="\t",
        index=False,
        header=False
    )

    sig_detail_out = os.path.join(
        out_dir,
        f"{region}_sharedH_selected_module_significant_test_details.tsv"
    )

    sig_df = sig_df.sort_values(
        ["Region", "Module", "test", "Comparison", "pval_adj", "pval"],
        na_position="last"
    )

    sig_df.to_csv(sig_detail_out, sep="\t", index=False)

    # ============================================================
    # Boxplot
    # ============================================================

    plt.figure(figsize=(10, 4))

    ax = sns.boxplot(
        data=df,
        x="Module",
        y="Value",
        hue="Pathology",
        order=module_order,
        hue_order=pathology_order,
        showfliers=False
    )

    plt.xticks(rotation=45)
    plt.title(f"{region} shared H: Module activity vs Pathology")

    y_min = df["Value"].quantile(0.01)
    y_max = df["Value"].quantile(0.98)
    y_range = y_max - y_min

    if y_range == 0:
        y_range = 1

    label_y = y_max + 0.05 * y_range

    for i, module in enumerate(module_order):

        row = kw_df[kw_df["Module"] == module]

        if row.empty:
            continue

        fdr = row["pval_adj"].values[0]
        label = p_to_star(fdr)

        if label != "":
            ax.text(
                i,
                label_y,
                label,
                ha="center",
                va="bottom",
                fontsize=14,
                fontweight="bold"
            )

    ax.set_ylim(y_min, y_max + 0.20 * y_range)

    plt.tight_layout()

    boxplot_path = os.path.join(
        boxplot_dir,
        f"{region}_sharedH_module_activity_vs_pathology.pdf"
    )

    plt.savefig(boxplot_path)
    plt.close()

    # ============================================================
    # Summary for this region
    # ============================================================

    print("\nDone region:", region)
    print(f"Kruskal-Wallis results saved to: {kw_out}")
    print(f"Pairwise rank-sum results saved to: {rs_out}")
    print(f"Linear regression results saved to: {lm_out}")
    print(f"Spearman correlation results saved to: {spearman_out}")
    print(f"Combined result table saved to: {combined_out}")
    print(f"Selected significant shared H modules saved to: {selected_module_out}")
    print(f"Significant test details saved to: {sig_detail_out}")
    print(f"Boxplot saved to: {boxplot_path}")

    print("\nNumber of selected significant shared H modules:")
    print(selected_modules.shape[0])

    print("\nSelected significant shared H modules:")
    print(selected_modules.to_string(index=False))

    return {
        "Region": region,
        "H_path": h_path,
        "metadata_path": meta_path,
        "n_modules": n_modules,
        "n_samples": n_samples,
        "n_nonAD": n_nonAD,
        "n_earlyAD": n_earlyAD,
        "n_lateAD": n_lateAD,
        "n_selected_modules": selected_modules.shape[0],
        "selected_modules": ",".join(selected_modules.astype(str).tolist())
    }


# ============================================================
# Main
# ============================================================

parser = argparse.ArgumentParser(description="Test shared module pathology association")
parser.add_argument("--region", choices=regions, help="Run one region only")
args = parser.parse_args()
all_summary = []

for region in ([args.region] if args.region else regions):
    result = analyze_one_region(region)
    if result is not None:
        all_summary.append(result)

if len(all_summary) > 0:
    summary_df = pd.DataFrame(all_summary)
    summary_out = os.path.join(
        OUT_ROOT,
        f"{args.region}_pathology_test_summary.tsv" if args.region else "region_sharedH_pathology_test_summary.tsv"
    )
    summary_df.to_csv(summary_out, sep="\t", index=False)
    print("\nSaved region summary to:", summary_out)

print("\nAll done.")
