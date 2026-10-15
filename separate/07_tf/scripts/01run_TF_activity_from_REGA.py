#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import glob
import numpy as np
import pandas as pd


# ============================================================
# Paths
# ============================================================

REGA_ROOT = os.environ.get(
    "REGA_MATRICES_ROOT",
    "/project/zduren/durenlab/naqing_palmetto/REGA_separate/03_process_rega/matrices",
)

GEX_ROOT = os.environ.get(
    "REGA_INPUTS_ROOT",
    "/project/zduren/durenlab/naqing_palmetto/REGA_separate/01_prepare/data/rega_inputs",
)

MOTIF_TF_FILE = (
    "/project/zduren/durenlab/naqing_palmetto/REGA_separate/"
    "07_tf/data/Motif_TF_name_correspondance_uppercase.txt"
)

OUT_SUBDIR = "TF_activity_pcc_pos_only"

CORR_THRESHOLD = 0.3

REGIONS = ["AG", "EC", "HC", "MTC", "PFC", "TH"]
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
    Try both exact name and .csv.
    Example:
        Ast_AG_V_TF_by_module
        Ast_AG_V_TF_by_module.csv
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
    else:
        return pd.read_csv(path, index_col=0)


def parse_cell_region(dirname):
    """
    Parse directory name like:
        Ast_AG
        Mic_Immune_MTC
    """
    base = os.path.basename(dirname)

    for region in REGIONS:
        suffix = "_" + region
        if base.endswith(suffix):
            celltype = base[:-len(suffix)]
            return celltype, region

    return None, None


def orient_v(V, motif_tf):
    """
    Make sure V is motif × module.

    Your expected file:
        {celltype}_{region}_V_TF_by_module
    should already be motif × module.
    But this function checks and transposes if needed.
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

    Your expected file:
        {celltype}_{region}_H_module_by_sample
    should already be module × sample.
    But this function checks and transposes if needed.
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

    vh = pd.DataFrame(
        vh_mat,
        index=V.index.astype(str),
        columns=H.columns.astype(str)
    )

    return vh


def get_gex_path(region, celltype):
    path = os.path.join(GEX_ROOT, region, celltype, "input_GEX.csv")

    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing GEX file: {path}")

    return path


# ============================================================
# Main computation for one celltype-region
# ============================================================

def run_one(celltype, region, out_dir_cell, motif_tf):
    label = f"{celltype}_{region}"

    print("\n" + "=" * 80)
    print(f"[RUN] {label}")

    V_name = f"{label}_V_TF_by_module"
    H_name = f"{label}_H_module_by_sample"

    V_path = find_file(out_dir_cell, V_name)
    H_path = find_file(out_dir_cell, H_name)

    if V_path is None:
        print(f"[WARN] missing V file for {label}, skip.")
        return

    if H_path is None:
        print(f"[WARN] missing H file for {label}, skip.")
        return

    print(f"[V] {V_path}")
    print(f"[H] {H_path}")

    V = read_matrix(V_path)
    H = read_matrix(H_path)

    print(f"[V raw] {V.shape}")
    print(f"[H raw] {H.shape}")

    vh = compute_vh(V, H, motif_tf)

    print(f"[VH = V @ H] motif × sample = {vh.shape}")

    gex_path = get_gex_path(region, celltype)

    print(f"[GEX] {gex_path}")

    gex = pd.read_csv(gex_path, index_col=0)
    gex.index = gex.index.astype(str)
    gex.columns = gex.columns.astype(str)

    print(f"[GEX] gene × sample = {gex.shape}")

    # ------------------------------------------------------------
    # Align samples
    # ------------------------------------------------------------
    samples_vh = vh.columns.astype(str).tolist()
    samples_gex = set(gex.columns.astype(str))

    common_samples = [s for s in samples_vh if s in samples_gex]

    if len(common_samples) == 0:
        print(f"[WARN] {label}: no common samples between VH and GEX, skip.")
        return

    vh_sub = vh.loc[:, common_samples]
    gex_sub = gex.loc[:, common_samples]

    print(f"[INFO] common samples = {len(common_samples)}")

    # ------------------------------------------------------------
    # Match motif-TF mapping
    # ------------------------------------------------------------
    motifs_vh = set(vh_sub.index.astype(str))

    motif_tf_sub = motif_tf[motif_tf["motif"].isin(motifs_vh)].copy()

    if motif_tf_sub.empty:
        print(f"[WARN] {label}: no motif matched to mapping, skip.")
        return

    # ------------------------------------------------------------
    # Match TF to expression genes
    # ------------------------------------------------------------
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

    # ------------------------------------------------------------
    # M1: TF × motif binary matrix
    # ------------------------------------------------------------
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

    # ------------------------------------------------------------
    # M2: TF expression vs motif activity PCC
    # ------------------------------------------------------------
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

    # ------------------------------------------------------------
    # TF activity
    # ------------------------------------------------------------
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

    # ------------------------------------------------------------
    # Save
    # ------------------------------------------------------------
    out_dir = os.environ.get(
        "REGA_TF_ACTIVITY_ROOT",
        "/project/zduren/durenlab/naqing_palmetto/REGA_separate/07_tf/results/activity",
    )
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

    dirs = sorted(glob.glob(os.path.join(REGA_ROOT, "*_*")))

    jobs = []

    for d in dirs:
        if not os.path.isdir(d):
            continue

        celltype, region = parse_cell_region(d)

        if celltype is None or celltype not in CELLTYPES:
            continue

        jobs.append((celltype, region, d))

    print(f"[INFO] jobs found = {len(jobs)}")

    if len(jobs) == 0:
        raise SystemExit("No jobs found. Check REGA_ROOT.")

    for celltype, region, d in jobs:
        run_one(celltype, region, d, motif_tf)


if __name__ == "__main__":
    main()
