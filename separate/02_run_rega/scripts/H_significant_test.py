import os
import numpy as np
import pandas as pd
import scipy.io
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import kruskal, mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm

# ============================================================
# Paths
# ============================================================

BASE_DIR = "/project/zduren/durenlab/naqing_palmetto/REGA_separate/"
# Model outputs and matching donor metadata.
INPUT_DIR = os.environ.get("REGA_STAGE2_ROOT", os.path.join(BASE_DIR, "02_run_rega/rega_output"))
META_DIR = os.environ.get("REGA_METADATA_ROOT", os.path.join(BASE_DIR, "01_prepare/data/metadata_by_group"))

OUT_DIR = os.environ.get("REGA_PATHOLOGY_ROOT", os.path.join(BASE_DIR, "02_run_rega/pathology_association"))
BOXPLOT_DIR = os.path.join(OUT_DIR, "boxplot")

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(BOXPLOT_DIR, exist_ok=True)


# ============================================================
# Cell types and regions
# ============================================================

celltypes = ["Exc", "OPC", "Oli", "Ast", "Mic_Immune", "Inh"]
regions = ["AG", "EC", "HC", "MTC", "PFC", "TH"]

pathology_order = ["nonAD", "earlyAD", "lateAD"]

pathology_mapping = {
    "nonAD": 1,
    "earlyAD": 2,
    "lateAD": 3
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


# ============================================================
# Store all results
# ============================================================

kw_all = []
rs_all = []
lm_all = []
spearman_all = []


# ============================================================
# Main loop
# ============================================================

for celltype in celltypes:
    for region in regions:

        celltype_region = f"{celltype}_{region}"

        print(f"\nProcessing {celltype_region}")

        mat_path = os.path.join(INPUT_DIR, region, celltype, "stage2.mat")
        meta_path = os.path.join(
            META_DIR,
            f"{celltype_region}_metadata.tsv"
        )

        if not os.path.exists(mat_path):
            print(f"  Skipping: missing H matrix: {mat_path}")
            continue

        if not os.path.exists(meta_path):
            print(f"  Skipping: missing metadata: {meta_path}")
            continue

        # ------------------------------------------------------------
        # Load H matrix
        # ------------------------------------------------------------

        mat = scipy.io.loadmat(mat_path)

        if "H" not in mat:
            print(f"  Skipping: no key 'H' in {mat_path}")
            continue

        H = mat["H"]

        # H expected shape: modules x samples
        n_modules, n_samples = H.shape
        module_order = [f"M{i + 1}" for i in range(n_modules)]

        # ------------------------------------------------------------
        # Load metadata
        # ------------------------------------------------------------

        meta = pd.read_csv(meta_path, sep="\t")

        if H.shape[1] != meta.shape[0]:
            print(
                f"  Skipping: mismatch H columns != metadata rows. "
                f"H has {H.shape[1]} samples, metadata has {meta.shape[0]} rows."
            )
            continue

        if "Pathology" in meta.columns:
            pathology = meta["Pathology"].astype(str).str.strip()
        else:
            pathology = meta.iloc[:, 7].astype(str).str.strip()

        pathology_num = pathology.map(pathology_mapping)

        n_counts = pathology.value_counts().to_dict()
        n_nonAD = n_counts.get("nonAD", 0)
        n_earlyAD = n_counts.get("earlyAD", 0)
        n_lateAD = n_counts.get("lateAD", 0)

        # ------------------------------------------------------------
        # Reshape H to long format
        # ------------------------------------------------------------

        df_list = []

        for i in range(n_modules):
            tmp = pd.DataFrame({
                "Module": f"M{i + 1}",
                "Value": H[i, :],
                "Pathology": pathology.values,
                "Pathology_num": pathology_num.values,
                "celltype": celltype,
                "region": region,
                "celltype_region": celltype_region
            })
            df_list.append(tmp)

        df = pd.concat(df_list, ignore_index=True)

        df["Pathology"] = pd.Categorical(
            df["Pathology"],
            categories=pathology_order,
            ordered=True
        )

        # ============================================================
        # 1. Kruskal-Wallis test
        # ============================================================

        kw_results = []

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

            kw_results.append({
                "test": "Kruskal_Wallis",
                "celltype_region": celltype_region,
                "celltype": celltype,
                "region": region,
                "Module": m,
                "celltype_region_module": f"{celltype_region}_{m}",
                "Comparison": "nonAD vs earlyAD vs lateAD",
                "statistic": stat,
                "effect_size": np.nan,
                "pval": pval,
                "n": sub["Value"].dropna().shape[0],
                "n_nonAD": n_nonAD,
                "n_earlyAD": n_earlyAD,
                "n_lateAD": n_lateAD
            })

        kw_group_df = pd.DataFrame(kw_results)

        # FDR within this celltype-region across modules
        kw_group_df = add_fdr(
            kw_group_df,
            p_col="pval",
            out_col="pval_adj"
        )

        kw_all.extend(kw_group_df.to_dict("records"))

        # ============================================================
        # 2. Pairwise rank-sum tests
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
                    "test": "Pairwise_MannWhitneyU",
                    "celltype_region": celltype_region,
                    "celltype": celltype,
                    "region": region,
                    "Module": m,
                    "celltype_region_module": f"{celltype_region}_{m}",
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
                "test": "Linear_Regression",
                "celltype_region": celltype_region,
                "celltype": celltype,
                "region": region,
                "Module": m,
                "celltype_region_module": f"{celltype_region}_{m}",
                "Comparison": "H module value ~ pathology stage",
                "Pathology_encoding": "nonAD=1, earlyAD=2, lateAD=3",
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
                "test": "Spearman",
                "celltype_region": celltype_region,
                "celltype": celltype,
                "region": region,
                "Module": m,
                "celltype_region_module": f"{celltype_region}_{m}",
                "Comparison": "Spearman: H module value vs pathology stage",
                "Pathology_encoding": "nonAD=1, earlyAD=2, lateAD=3",
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
        # Boxplot
        # Use Kruskal-Wallis FDR labels
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
        plt.title(f"{celltype_region}: Module activity vs Pathology")

        y_min = df["Value"].quantile(0.01)
        y_max = df["Value"].quantile(0.98)
        y_range = y_max - y_min

        if y_range == 0:
            y_range = 1

        label_y = y_max + 0.05 * y_range

        for i, module in enumerate(module_order):

            row = kw_group_df[kw_group_df["Module"] == module]

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
            BOXPLOT_DIR,
            f"{celltype_region}.pdf"
        )

        plt.savefig(boxplot_path)
        plt.close()


# ============================================================
# Combine result tables
# ============================================================

kw_df = pd.DataFrame(kw_all)
rs_df = pd.DataFrame(rs_all)
lm_df = pd.DataFrame(lm_all)
spearman_df = pd.DataFrame(spearman_all)


# ============================================================
# FDR correction
# ============================================================

# Kruskal-Wallis:
# For each celltype-region, adjust across modules
kw_df = (
    kw_df
    .groupby("celltype_region", group_keys=False)
    .apply(lambda x: add_fdr(x, p_col="pval", out_col="pval_adj"))
)

# Pairwise rank-sum:
# For each celltype-region and each comparison, adjust across modules
rs_df = (
    rs_df
    .groupby(["celltype_region", "Comparison"], group_keys=False)
    .apply(lambda x: add_fdr(x, p_col="pval", out_col="pval_adj"))
)

# Linear regression:
# For each celltype-region, adjust across modules
lm_df = (
    lm_df
    .groupby("celltype_region", group_keys=False)
    .apply(lambda x: add_fdr(x, p_col="pval", out_col="pval_adj"))
)

# Spearman:
# For each celltype-region, adjust across modules
spearman_df = (
    spearman_df
    .groupby("celltype_region", group_keys=False)
    .apply(lambda x: add_fdr(x, p_col="pval", out_col="pval_adj"))
)


# Optional global FDR columns across all celltype-regions
kw_df = add_fdr(kw_df, p_col="pval", out_col="pval_adj_global")
rs_df = add_fdr(rs_df, p_col="pval", out_col="pval_adj_global")
lm_df = add_fdr(lm_df, p_col="pval", out_col="pval_adj_global")
spearman_df = add_fdr(spearman_df, p_col="pval", out_col="pval_adj_global")


# ============================================================
# Sort result tables
# ============================================================

kw_df = kw_df.sort_values(
    ["celltype_region", "pval_adj", "pval"],
    na_position="last"
)

rs_df = rs_df.sort_values(
    ["celltype_region", "Comparison", "pval_adj", "pval"],
    na_position="last"
)

lm_df = lm_df.sort_values(
    ["celltype_region", "pval_adj", "pval"],
    na_position="last"
)

spearman_df = spearman_df.sort_values(
    ["celltype_region", "pval_adj", "pval"],
    na_position="last"
)


# ============================================================
# Save individual result tables
# ============================================================

kw_out = os.path.join(
    OUT_DIR,
    "Kruskal_Wallis_all_celltype_region.tsv"
)

rs_out = os.path.join(
    OUT_DIR,
    "Pairwise_MannWhitneyU_all_celltype_region.tsv"
)

lm_out = os.path.join(
    OUT_DIR,
    "LinearRegression_all_celltype_region.tsv"
)

spearman_out = os.path.join(
    OUT_DIR,
    "Spearman_pathology_H_modules_all_celltype_region.tsv"
)

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
    ["celltype_region", "Module", "test", "Comparison"],
    na_position="last"
)

combined_out = os.path.join(
    OUT_DIR,
    "All_four_tests_H_pathology_results.tsv"
)

combined_df.to_csv(combined_out, sep="\t", index=False)


# ============================================================
# Diagnostic broad selection: significant in any of the four tests.
# Publication candidates are selected separately by select_pathology_modules.R
# using linear regression or nonAD-vs-lateAD Mann-Whitney only.
# ============================================================

sig_df = combined_df[
    combined_df["pval_adj"].notna() &
    (combined_df["pval_adj"] < 0.05)
].copy()

selected_modules = (
    sig_df["celltype_region_module"]
    .dropna()
    .drop_duplicates()
    .sort_values()
)

selected_module_out = os.path.join(
    OUT_DIR,
    "selected_module_any_test.txt"
)

selected_modules.to_csv(
    selected_module_out,
    sep="\t",
    index=False,
    header=False
)


# Also save detailed significant rows so you know why each module was selected
sig_detail_out = os.path.join(
    OUT_DIR,
    "selected_module_any_test_details.tsv"
)

sig_df = sig_df.sort_values(
    ["celltype_region_module", "test", "Comparison", "pval_adj", "pval"],
    na_position="last"
)

sig_df.to_csv(sig_detail_out, sep="\t", index=False)


# ============================================================
# Summary
# ============================================================

print("\nDone.")
print(f"Kruskal-Wallis results saved to: {kw_out}")
print(f"Pairwise rank-sum results saved to: {rs_out}")
print(f"Linear regression results saved to: {lm_out}")
print(f"Spearman correlation results saved to: {spearman_out}")
print(f"Combined result table saved to: {combined_out}")
print(f"Selected significant modules saved to: {selected_module_out}")
print(f"Significant test details saved to: {sig_detail_out}")
print(f"Boxplots saved to: {BOXPLOT_DIR}")

print("\nNumber of selected significant modules:")
print(selected_modules.shape[0])
