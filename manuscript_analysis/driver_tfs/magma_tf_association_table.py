#!/usr/bin/env python3
"""Create nominally significant MAGMA TF results and their TF-name list."""
from pathlib import Path
import os

import pandas as pd


PROJECT = Path(os.environ.get("AD_MANUSCRIPT_WORK_ROOT", Path(__file__).resolve().parents[3])).expanduser().resolve()
REGA_SEPARATE = Path(os.environ.get("REGA_SEPARATE_ROOT", PROJECT.parent / "REGA_separate")).expanduser().resolve()
MAGMA_RESULTS = REGA_SEPARATE / "06_magma/results/merged/AD_Bellenguez_hg19_merged.genes.out"
MOTIF_MATRIX = REGA_SEPARATE / "07_tf/data/motif_by_TF_matrix.tsv"
TABLE_OUT = PROJECT / "manuscript/Supplementary_Table_MAGMA_TF_SNP.tsv"
LIST_OUT = PROJECT / "manuscript/MAGMA_TFs_with_nominal_SNP_association.txt"


def main() -> None:
    magma = pd.read_csv(MAGMA_RESULTS, sep=r"\s+", engine="python")
    motif_header = pd.read_csv(MOTIF_MATRIX, sep="\t", nrows=0).columns
    tf_universe = {str(name).strip().upper() for name in motif_header if name != "motif"}
    magma["GENE"] = magma["GENE"].astype(str).str.strip().str.upper()
    magma = magma[magma["GENE"].isin(tf_universe)].copy()
    magma = magma.drop_duplicates("GENE", keep="first")
    magma["MAGMA_gene_FDR"] = magma["P"].rank(method="first", ascending=True)
    # Benjamini-Hochberg over every motif-derived TF represented in the MAGMA
    # output; correct before applying the requested nominal P < 0.05 filter.
    n_tests = len(magma)
    magma["MAGMA_gene_FDR"] = (magma["P"] * n_tests / magma["MAGMA_gene_FDR"]).clip(upper=1)
    magma = magma.sort_values("P", kind="stable").copy()
    # Enforce monotonicity of the BH adjusted values in ranked order.
    magma["MAGMA_gene_FDR"] = magma["MAGMA_gene_FDR"].iloc[::-1].cummin().iloc[::-1]
    magma["FDR_lt_0.05"] = magma["MAGMA_gene_FDR"] < 0.05
    significant = magma[magma["P"] < 0.05].copy()
    output_columns = [
        "GENE", "CHR", "START", "STOP", "NSNPS", "NPARAM", "N",
        "ZSTAT", "P", "MAGMA_gene_FDR", "FDR_lt_0.05",
    ]
    significant = significant[output_columns].rename(
        columns={"GENE": "TF", "P": "MAGMA_gene_P", "ZSTAT": "MAGMA_gene_ZSTAT"}
    )
    significant.to_csv(TABLE_OUT, sep="\t", index=False, na_rep="NA", float_format="%.8g")
    LIST_OUT.write_text("\n".join(significant["TF"].tolist()) + "\n")
    print(f"TF universe present in MAGMA: {n_tests}")
    print(f"Nominal MAGMA P < 0.05 TFs: {len(significant)}")
    print(f"BH FDR < 0.05 TFs among nominal hits: {int(significant['FDR_lt_0.05'].sum())}")
    print(f"Wrote {TABLE_OUT}")
    print(f"Wrote {LIST_OUT}")


if __name__ == "__main__":
    main()
