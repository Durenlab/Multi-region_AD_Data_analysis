#!/usr/bin/env python3
from pathlib import Path
import os
SHARED_ROOT = Path(os.environ["REGA_SHARED_ROOT"]) if "REGA_SHARED_ROOT" in os.environ else next(
    p for p in Path(__file__).resolve().parents if p.name in {"REGA_shared", "shared"}
)
SEPARATE_ROOT = Path(os.environ["REGA_SEPARATE_ROOT"])
# -*- coding: utf-8 -*-

import os
import re
import glob
import numpy as np
import pandas as pd


# ============================================================
# Paths
# ============================================================

REGA_ROOT = str(SHARED_ROOT / "part3_process_REGA_output/09extracted_matrices_region_shareH")

GEX_ROOT = str(SHARED_ROOT / "part1_prepare_REGA/03REGA_input_region_consistent")

MOTIF_TF_FILE = str(SEPARATE_ROOT / "07_tf/data/Motif_TF_name_correspondance_uppercase.txt")

OUT_SUBDIR = "TF_activity_pcc_pos_only"
CORR_THRESHOLD = 0.3

REGIONS = ["AG", "EC", "HC", "MTC", "PFC", "TH"]
if os.environ.get("REGA_REGION"):
    if os.environ["REGA_REGION"] not in REGIONS:
        raise ValueError(f"Unknown REGA_REGION: {os.environ['REGA_REGION']}")
    REGIONS = [os.environ["REGA_REGION"]]
CELLTYPES = ["Ast", "Exc", "OPC", "Oli", "Inh", "Mic_Immune"]



# ============================================================
# Helper functions
# ============================================================

def load_motif_tf_mapping(path):
    """
    Motif-TF mapping file:
        motif TF
    No header expected.
    """
    df = pd.read_csv(
        path,
        sep=r"\s+|\t",
        engine="python",
        header=None,
        names=["motif", "TF"]
    )
    df["motif"] = df["motif"].astype(str)
    df["TF"] = df["TF"].astype(str).str.upper()
    df = df.drop_duplicates()
    return df


def compute_corr_matrix(tf_expr, vh):
    """
    Pearson correlation between TF expression and motif activity.

    tf_expr: TF × sample
    vh:      motif × sample

    return:
        TF × motif correlation matrix
    """
    X = tf_expr.to_numpy(dtype=float)
    Y = vh.to_numpy(dtype=float)

    T, S = X.shape
    M, S2 = Y.shape

    if S != S2:
        raise ValueError(f"Sample mismatch: tf_expr has {S}, VH has {S2}")

    X_mean = X.mean(axis=1, keepdims=True)
    Y_mean = Y.mean(axis=1, keepdims=True)

    X_std = X.std(axis=1, ddof=1, keepdims=True)
    Y_std = Y.std(axis=1, ddof=1, keepdims=True)

    X_std[X_std == 0] = 1.0
    Y_std[Y_std == 0] = 1.0

    Xz = (X - X_mean) / X_std
    Yz = (Y - Y_mean) / Y_std

    corr = (Xz @ Yz.T) / (S - 1)
    corr = np.nan_to_num(corr, nan=0.0)

    return pd.DataFrame(corr, index=tf_expr.index, columns=vh.index)


def find_file(directory, filename_no_ext):
    """
    Try exact name, .csv, .tsv.
    """
    candidates = [
        os.path.join(directory, filename_no_ext),
        os.path.join(directory, filename_no_ext + ".csv"),
        os.path.join(directory, filename_no_ext + ".tsv"),
    ]
    for p in candidates:
        if os.path.exists(p):
            return p
    return None


def read_matrix(path):
    """
    Read csv/tsv automatically.
    First column is row name.
    """
    if path.endswith(".tsv"):
        return pd.read_csv(path, sep="\t", index_col=0)
    return pd.read_csv(path, index_col=0)


def parse_label_from_vfile(v_path, region):
    """
    From:
        Vasc_MTC_V_TF_by_module.tsv
    return:
        Vasc
    """
    base = os.path.basename(v_path)
    suffix = f"_{region}_V_TF_by_module"
    if base.endswith(".tsv"):
        base = base[:-4]
    elif base.endswith(".csv"):
        base = base[:-4]

    if not base.endswith(suffix):
        return None

    celltype = base[: -len(suffix)]
    return celltype if celltype else None


def orient_v(V, motif_tf):
    """
    Make sure V is motif × module.
    """
    motif_set = set(motif_tf["motif"].astype(str))

    row_match = len(set(V.index.astype(str)) & motif_set)
    col_match = len(set(V.columns.astype(str)) & motif_set)

    if row_match >= col_match:
        V.index = V.index.astype(str)
        V.columns = V.columns.astype(str)
        return V
    else:
        V = V.T
        V.index = V.index.astype(str)
        V.columns = V.columns.astype(str)
        return V


def orient_h(H, module_names):
    """
    Make sure H is module × sample.
    """
    module_set = set([str(x) for x in module_names])

    row_match = len(set(H.index.astype(str)) & module_set)
    col_match = len(set(H.columns.astype(str)) & module_set)

    if row_match >= col_match:
        H.index = H.index.astype(str)
        H.columns = H.columns.astype(str)
        return H
    else:
        H = H.T
        H.index = H.index.astype(str)
        H.columns = H.columns.astype(str)
        return H


def compute_vh(V, H, motif_tf):
    """
    Compute VH = V @ H.
    V: motif × module
    H: module × sample
    VH: motif × sample
    """
    V = orient_v(V, motif_tf)
    H = orient_h(H, V.columns)

    common_modules = [m for m in V.columns if m in set(H.index)]
    if len(common_modules) == 0:
        raise ValueError("No common modules between V columns and H rows.")

    V = V.loc[:, common_modules]
    H = H.loc[common_modules, :]

    vh_mat = V.to_numpy(dtype=float) @ H.to_numpy(dtype=float)

    return pd.DataFrame(
        vh_mat,
        index=V.index.astype(str),
        columns=H.columns.astype(str)
    )


def find_gex_path(region, celltype):
    """
    GEX file example:
        Ast_AG.GEX_sorted.csv
    """
    label = f"{celltype}_{region}"

    candidates = [
        os.path.join(GEX_ROOT, f"{label}.GEX_sorted.csv"),
        os.path.join(GEX_ROOT, f"{label}.GEX_sorted.tsv"),
        os.path.join(GEX_ROOT, f"{label}_GEX_sorted.csv"),
        os.path.join(GEX_ROOT, f"{label}_GEX_sorted.tsv"),
        os.path.join(GEX_ROOT, region, f"{label}.GEX_sorted.csv"),
        os.path.join(GEX_ROOT, region, f"{label}.GEX_sorted.tsv"),
        os.path.join(GEX_ROOT, region, celltype, f"{label}.GEX_sorted.csv"),
        os.path.join(GEX_ROOT, region, celltype, f"{label}.GEX_sorted.tsv"),
    ]

    for p in candidates:
        if os.path.exists(p):
            return p

    # fallback: search anywhere under GEX_ROOT
    hits = glob.glob(os.path.join(GEX_ROOT, "**", f"{label}*.csv"), recursive=True)
    if len(hits) > 0:
        return hits[0]

    hits = glob.glob(os.path.join(GEX_ROOT, "**", f"{label}*.tsv"), recursive=True)
    if len(hits) > 0:
        return hits[0]

    raise FileNotFoundError(f"Missing GEX file for {label} under {GEX_ROOT}")


def discover_jobs():
    """
    Scan region directories, then pair:
        {celltype}_{region}_V_TF_by_module(.csv/.tsv)
        {celltype}_{region}_H_module_by_sample(.csv/.tsv)
    """
    jobs = []

    for region in REGIONS:
        region_dir = os.path.join(REGA_ROOT, region)
        if not os.path.isdir(region_dir):
            print(f"[WARN] region dir not found: {region_dir}")
            continue

        v_candidates = sorted(glob.glob(os.path.join(region_dir, f"*_{region}_V_TF_by_module*")))
        if len(v_candidates) == 0:
            print(f"[WARN] no V files found in {region_dir}")
            continue

        for v_path in v_candidates:
            celltype = parse_label_from_vfile(v_path, region)
            if celltype is None:
                continue
            if celltype not in CELLTYPES:
                print(f"[WARN] unexpected celltype '{celltype}' from {v_path}, skip.")
                continue

            h_name = f"{region}_H_module_by_sample"
            h_path = find_file(region_dir, h_name)
            if h_path is None:
                print(f"[WARN] missing H file for {celltype}_{region}, skip.")
                continue

            jobs.append((celltype, region, region_dir, v_path, h_path))

    return jobs


# ============================================================
# Main computation for one celltype-region
# ============================================================

def run_one(celltype, region, region_dir, v_path, h_path, motif_tf):
    label = f"{celltype}_{region}"

    print("\n" + "=" * 80)
    print(f"[RUN] {label}")
    print(f"[V] {v_path}")
    print(f"[H] {h_path}")

    V = read_matrix(v_path)
    H = read_matrix(h_path)

    print(f"[V raw] {V.shape}")
    print(f"[H raw] {H.shape}")

    vh = compute_vh(V, H, motif_tf)
    print(f"[VH = V @ H] motif × sample = {vh.shape}")

    gex_path = find_gex_path(region, celltype)
    print(f"[GEX] {gex_path}")

    gex = pd.read_csv(gex_path, index_col=0)
    gex.index = gex.index.astype(str)
    gex.columns = gex.columns.astype(str)

    print(f"[GEX] gene × sample = {gex.shape}")

    # Align samples
    samples_vh = vh.columns.astype(str).tolist()
    samples_gex = set(gex.columns.astype(str))
    common_samples = [s for s in samples_vh if s in samples_gex]

    if len(common_samples) == 0:
        print(f"[WARN] {label}: no common samples between VH and GEX, skip.")
        return

    vh_sub = vh.loc[:, common_samples]
    gex_sub = gex.loc[:, common_samples]

    print(f"[INFO] common samples = {len(common_samples)}")

    # Match motif-TF mapping
    motifs_vh = set(vh_sub.index.astype(str))
    motif_tf_sub = motif_tf[motif_tf["motif"].isin(motifs_vh)].copy()

    if motif_tf_sub.empty:
        print(f"[WARN] {label}: no motif matched to mapping, skip.")
        return

    # Match TF to expression genes
    gene_upper_to_original = {}
    for gene in gex_sub.index.astype(str):
        gene_upper_to_original[gene.upper()] = gene

    tf_all = pd.unique(motif_tf_sub["TF"].astype(str).str.upper())
    tf_upper_in_expr = [tf for tf in tf_all if tf in gene_upper_to_original]

    if len(tf_upper_in_expr) == 0:
        print(f"[WARN] {label}: no TF found in GEX, skip.")
        return

    motif_tf_sub["TF"] = motif_tf_sub["TF"].astype(str).str.upper()
    motif_tf_sub = motif_tf_sub[motif_tf_sub["TF"].isin(tf_upper_in_expr)]

    tf_original_names = [gene_upper_to_original[tf] for tf in tf_upper_in_expr]
    tf_expr = gex_sub.loc[tf_original_names].copy()
    tf_expr.index = tf_upper_in_expr

    print(f"[INFO] TF in expr = {len(tf_upper_in_expr)}")
    print(f"[INFO] motif-TF pairs filtered = {motif_tf_sub.shape[0]}")
    print(f"[TF expression] TF × sample = {tf_expr.shape}")

    # M1: TF × motif binary matrix
    motifs_used = vh_sub.index.astype(str).tolist()

    motif_tf_sub["val"] = 1.0
    M1 = motif_tf_sub.pivot_table(
        index="TF",
        columns="motif",
        values="val",
        aggfunc="max",
        fill_value=0.0
    )

    M1 = M1.reindex(
        index=tf_upper_in_expr,
        columns=motifs_used,
        fill_value=0.0
    )

    print(f"[M1] TF × motif = {M1.shape}")

    # M2: TF expression vs motif activity PCC
    vh_sub = vh_sub.loc[M1.columns]
    M2 = compute_corr_matrix(tf_expr, vh_sub)

    print(f"[M2] TF × motif corr = {M2.shape}")

    M2_pos = M2.clip(lower=0.0)

    if CORR_THRESHOLD is not None and CORR_THRESHOLD > 0:
        before_nnz = (M2_pos.to_numpy() != 0).sum()
        M2_thr = M2_pos.mask(M2_pos < CORR_THRESHOLD, 0.0)
        after_nnz = (M2_thr.to_numpy() != 0).sum()

        print(
            f"[M2] negatives -> 0 and rho >= {CORR_THRESHOLD}, "
            f"nonzero {before_nnz} -> {after_nnz}"
        )
        M2 = M2_thr
    else:
        before_nnz = (M2_pos.to_numpy() != 0).sum()
        print(f"[M2] negatives -> 0, nonzero = {before_nnz}")
        M2 = M2_pos

    # TF activity
    W_tf_motif = M1.to_numpy(dtype=float) * M2.to_numpy(dtype=float)
    VH_mat = vh_sub.to_numpy(dtype=float)
    TF_activity = W_tf_motif @ VH_mat

    df_TF_activity = pd.DataFrame(
        TF_activity,
        index=M1.index,
        columns=common_samples
    )

    print(f"[TF activity] TF × sample = {df_TF_activity.shape}")

    # TF expression × activity
    TF_expr_times_activity = tf_expr.to_numpy(dtype=float) * TF_activity
    df_TF_expr_times_activity = pd.DataFrame(
        TF_expr_times_activity,
        index=tf_expr.index,
        columns=common_samples
    )

    # Save
    out_dir = os.path.join(SHARED_ROOT, "part7_TF", "03TF_activity", label, OUT_SUBDIR)
    os.makedirs(out_dir, exist_ok=True)

    out_vh = os.path.join(out_dir, f"{label}_VH_from_V_times_H.csv")
    out_act = os.path.join(out_dir, f"{label}_TF_activity_pcc_cutoff.csv")
    out_exp = os.path.join(out_dir, f"{label}_TF_expression_pcc_cutoff.csv")
    out_mul = os.path.join(out_dir, f"{label}_TF_expr_times_activity_pcc_cutoff.csv")
    out_m1 = os.path.join(out_dir, f"{label}_TF_motif_binary_M1.csv")
    out_m2 = os.path.join(out_dir, f"{label}_TF_motif_pcc_M2.csv")

    vh_sub.to_csv(out_vh)
    df_TF_activity.to_csv(out_act)
    tf_expr.to_csv(out_exp)
    df_TF_expr_times_activity.to_csv(out_mul)
    M1.to_csv(out_m1)
    M2.to_csv(out_m2)

    print(f"[SAVE] {label}:")
    print(f"    {out_vh}")
    print(f"    {out_act}")
    print(f"    {out_exp}")
    print(f"    {out_mul}")
    print(f"    {out_m1}")
    print(f"    {out_m2}")


# ============================================================
# Main
# ============================================================

def main():
    motif_tf = load_motif_tf_mapping(MOTIF_TF_FILE)

    print(f"[INFO] motif-TF mapping pairs = {motif_tf.shape[0]}")
    print(f"[INFO] REGA_ROOT = {REGA_ROOT}")
    print(f"[INFO] GEX_ROOT = {GEX_ROOT}")

    jobs = discover_jobs()
    print(f"[INFO] jobs found = {len(jobs)}")

    if len(jobs) == 0:
        raise SystemExit("No jobs found. Check REGA_ROOT and file naming.")

    for celltype, region, region_dir, v_path, h_path in jobs:
        run_one(celltype, region, region_dir, v_path, h_path, motif_tf)


if __name__ == "__main__":
    main()
