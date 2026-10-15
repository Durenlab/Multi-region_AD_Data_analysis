#!/usr/bin/env python3

import os
import argparse
import numpy as np
import pandas as pd
import gseapy as gp
from statsmodels.stats.multitest import multipletests


# =========================================================
# arguments
# =========================================================

ap = argparse.ArgumentParser()

ap.add_argument("--wbv", required=True)
ap.add_argument("--gmt", required=True)
ap.add_argument("--celltype", required=True)
ap.add_argument("--context", required=True)
ap.add_argument("--out-root", required=True)

ap.add_argument("--perm", type=int, default=1000)
ap.add_argument("--threads", type=int, default=4)
ap.add_argument("--jitter", action="store_true")

args = ap.parse_args()


# =========================================================
# setup
# =========================================================

os.makedirs(args.out_root, exist_ok=True)

GENESET = f"{args.celltype}_AGING_UP"

OUT_FILE = os.path.join(
    args.out_root,
    f"{args.context}_aging_GSEA.tsv"
)

print("Context:", args.context)
print("Gene set:", GENESET)


# =========================================================
# 1. read WBV
# gene x module
# =========================================================

WBV = pd.read_csv(
    args.wbv,
    sep="\t",
    index_col=0
)

WBV.index = WBV.index.astype(str).str.upper()


# =========================================================
# 2. normalize
# column L1 -> row L1
# =========================================================

WBV = WBV.div(
    WBV.abs().sum(axis=0).replace(0, np.nan),
    axis=1
)

WBV = WBV.div(
    WBV.abs().sum(axis=1).replace(0, np.nan),
    axis=0
)

WBV = WBV.fillna(0)


# =========================================================
# 3. read corresponding aging gene set
# =========================================================

aging_genes = None

with open(args.gmt) as f:

    for line in f:

        x = line.rstrip().split("\t")

        if x[0] == GENESET:

            aging_genes = [
                g.upper()
                for g in x[2:]
            ]

            break


if aging_genes is None:
    raise ValueError(
        f"{GENESET} not found in GMT"
    )

aging_genes = list(set(aging_genes))

overlap = set(aging_genes) & set(WBV.index)

print("Aging genes:", len(aging_genes))
print("Overlap with WBV:", len(overlap))


# =========================================================
# 4. GSEA for all 20 modules
# =========================================================

results = []
gsea_objects = {}

for mod in WBV.columns:

    print("Running:", mod)

    rnk = WBV[mod].sort_values(
        ascending=False
    )

    rnk = pd.DataFrame({
        "gene": rnk.index,
        "score": rnk.values
    })

    # optional tiny jitter for tied scores
    if args.jitter:

        rng = np.random.default_rng(123)

        rnk["score"] = (
            rnk["score"].to_numpy()
            + rng.normal(
                0,
                1e-12,
                size=len(rnk)
            )
        )

    gsea = gp.prerank(
        rnk=rnk,
        gene_sets={GENESET: aging_genes},
        permutation_num=args.perm,
        min_size=5,
        max_size=5000,
        threads=args.threads,
        seed=123,
        no_plot=True,
        verbose=False
    )

    res = gsea.res2d.iloc[0]

    results.append({
        "context": args.context,
        "celltype": args.celltype,
        "geneset": GENESET,
        "module": mod,
        "ES": res["ES"],
        "NES": res["NES"],
        "pval": res["NOM p-val"],
        "leading_genes": res["Lead_genes"]
    })

    gsea_objects[mod] = gsea


# =========================================================
# 5. FDR across 20 modules
# =========================================================

results = pd.DataFrame(results)

results["ES"] = pd.to_numeric(
    results["ES"]
)

results["NES"] = pd.to_numeric(
    results["NES"]
)

results["pval"] = pd.to_numeric(
    results["pval"]
)

results["FDR"] = multipletests(
    results["pval"],
    method="fdr_bh"
)[1]


results["significant_positive"] = (
    (results["NES"] > 0) &
    (results["FDR"] < 0.05)
)


results = results.sort_values(
    ["FDR", "NES"],
    ascending=[True, False]
)


# =========================================================
# 6. save table
# =========================================================

results.to_csv(
    OUT_FILE,
    sep="\t",
    index=False
)

print("\nSaved:", OUT_FILE)


# =========================================================
# 7. plot significant positive modules
# =========================================================

plot_dir = os.path.join(
    args.out_root,
    "GSEA_plots"
)

os.makedirs(
    plot_dir,
    exist_ok=True
)


hits = results[
    results["significant_positive"]
]


print(
    "Significant positive modules:",
    hits["module"].tolist()
)


from gseapy import gseaplot

for mod in hits["module"]:

    print("Plotting:", mod)

    gsea = gsea_objects[mod]

    # get our BH-FDR across 20 modules
    module_fdr = results.loc[
        results["module"] == mod,
        "FDR"
    ].iloc[0]

    # copy GSEA result
    plot_res = gsea.results[GENESET].copy()

    # replace GSEApy internal FDR with our BH-FDR
    plot_res["fdr"] = float(module_fdr)

    outfile = os.path.join(
        plot_dir,
        f"{args.context}_{mod}_aging_GSEA.pdf"
    )

    gseaplot(
        rank_metric=gsea.ranking,
        term=GENESET,
        ofname=outfile,
        **plot_res
    )