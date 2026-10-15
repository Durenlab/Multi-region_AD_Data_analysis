#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, re, h5py
from pathlib import Path
SHARED_ROOT = Path(os.environ.get("REGA_SHARED_ROOT", Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd
from scipy import sparse
# ---------------- Paths ----------------
BASE_IN  = str(SHARED_ROOT / "part1_prepare_REGA/03REGA_input_region_consistent")
GENE_DIR = os.path.join(BASE_IN, "coords")
PEAK_DIR = os.path.join(BASE_IN, "coords")
OUT_ROOT = str(SHARED_ROOT / "part1_prepare_REGA/05REGA_blocks_region_consistent")

os.makedirs(OUT_ROOT, exist_ok=True)


def decode_bytes(v):
    return v.decode() if isinstance(v, (bytes, bytearray)) else str(v)

def process_one_cell(cell):
    print(f"\n==== {cell} ====")
    # ---- input files ----
    H5_W  = os.path.join(BASE_IN, f"{cell}.re_tg_matrix.h5")
    CSV_B = os.path.join(BASE_IN, f"{cell}.peak_motif_dense.csv")
    CSV_X = os.path.join(BASE_IN, f"{cell}.GEX_sorted.csv")

    gene_bed = os.path.join(GENE_DIR, f"{cell}.expr_genes_300k.bed")
    re_bed_path = os.path.join(PEAK_DIR, f"{cell}.PeakCalls_named.bed")

    for p in (H5_W, CSV_B, CSV_X, gene_bed, re_bed_path):
        if not os.path.exists(p) or (p.endswith(".bed") and os.path.getsize(p) == 0):
            print(f"[SKIP] 缺少或空文件：{p}")
            return

    # ---- 读取 W（HDF5 三元组）----
    with h5py.File(H5_W, "r") as f:
        rows   = f["/i"][:].astype(np.int64) - 1
        cols   = f["/j"][:].astype(np.int64) - 1
        values = f["/x"][:].astype(np.float32)
        dims   = f["/dims"][:].astype(np.int64)
        W_gene_names = [decode_bytes(x) for x in f["/rownames"][:]]
        W_re_names   = [decode_bytes(x) for x in f["/colnames"][:]]

    W_sparse = sparse.csr_matrix((values, (rows, cols)),
                                 shape=(int(dims[0]), int(dims[1])))
    print(f"[W] shape={W_sparse.shape}, nnz={W_sparse.nnz}")

    # ---- 读取 B/X ----
    B = pd.read_csv(CSV_B, index_col=0)  # RE × Motif
    X = pd.read_csv(CSV_X, index_col=0)  # Gene × Sample
    print(f"[B] {B.shape}  [X] {X.shape}")

    # ---- 基因坐标 ----
    gene_loc = pd.read_csv(gene_bed, sep=r"\s+", header=None,
                           names=["chr","start","end","gene"])
    gene_loc["gene"] = gene_loc["gene"].astype(str)
    gene_loc["chr"]  = gene_loc["chr"].astype(str)
    gene_loc["_len"] = gene_loc["end"] - gene_loc["start"]
    gene_loc_unique = (gene_loc.sort_values(["gene","_len"], ascending=[True, False])
                                .drop_duplicates("gene").drop(columns="_len")
                                .set_index("gene"))

    # ---- RE 坐标（合并峰）----
    re_bed = pd.read_csv(re_bed_path, sep=r"\s+", header=None,
                         usecols=[0,1,2,3],
                         names=["chr","start","end","re_id"])
    re_bed["re_id"] = re_bed["re_id"].astype(str)
    re_bed["chr"]   = re_bed["chr"].astype(str)
    if re_bed["re_id"].duplicated().any():
        re_bed["_len"] = re_bed["end"] - re_bed["start"]
        re_bed = (re_bed.sort_values(["re_id","_len"], ascending=[True, False])
                         .drop_duplicates("re_id").drop(columns="_len"))
    re_bed_unique = re_bed.set_index("re_id")

    # ---- 对齐 ----
    genes_W   = pd.Index(W_gene_names)
    genes_set = X.index.intersection(genes_W).intersection(gene_loc_unique.index)
    RE_W      = pd.Index(W_re_names)
    RE_set    = B.index.intersection(RE_W).intersection(re_bed_unique.index)
    print(f"[对齐] genes={len(genes_set)}  RE={len(RE_set)}")
    if len(genes_set)==0 or len(RE_set)==0:
        print("[SKIP] 对齐后为空")
        return

    X2 = X.loc[genes_set]
    gene2row = pd.Index(W_gene_names).get_indexer(genes_set)
    re2col   = pd.Index(W_re_names).get_indexer(RE_set)
    W2 = W_sparse[gene2row, :][:, re2col]
    B2 = B.loc[RE_set]

    loc_X  = gene_loc_unique.loc[genes_set]
    loc_RE = re_bed_unique.loc[RE_set]

    # ---- 染色体公共集合 ----
    common_chrs = sorted(set(loc_X["chr"]).intersection(set(loc_RE["chr"])))
    if not common_chrs:
        print("[SKIP] 无公共染色体")
        return

    # ---- 输出目录 ----
    out_dir = os.path.join(OUT_ROOT, f"{cell}_blocks_aligned")
    os.makedirs(out_dir, exist_ok=True)

    # ---- 分块保存 ----
    total_g = total_r = 0
    for chr_name in common_chrs:
        gmask = (loc_X["chr"].values == chr_name)
        rmask = (loc_RE["chr"].values == chr_name)
        if not gmask.any() or not rmask.any():
            continue

        X_blk = X2.iloc[gmask, :].values.astype(np.float32)
        W_blk = W2[gmask, :][:, rmask].toarray().astype(np.float32)
        B_blk = B2.iloc[rmask, :].values.astype(np.float32)

        if X_blk.size == 0 or W_blk.size == 0 or B_blk.size == 0:
            continue

        gi, n = X_blk.shape
        ri, p = B_blk.shape
        total_g += gi; total_r += ri
        print(f"[{cell}|{chr_name}] X:{gi}×{n}  W:{gi}×{ri}  B:{ri}×{p}")

        np.savez_compressed(
            os.path.join(out_dir, f"{chr_name}.npz"),
            X=X_blk, W=W_blk, B=B_blk,
            genes=X2.index[gmask].to_numpy(),
            RE=B2.index[rmask].to_numpy(),
            samples=X2.columns.to_numpy(),
            motifs=B2.columns.to_numpy(),
            gene_chr=loc_X["chr"].values[gmask],
            re_chr=loc_RE["chr"].values[rmask]
        )

    print(f"[完成 {cell}] 染色体块={len(common_chrs)}  总基因={total_g}  总RE={total_r}")
    print(f"[输出] {out_dir}")

def main():
    # 以 *.re_tg_matrix.h5 列表作为待处理细胞类型
    cells = sorted({
        re.sub(r"\.re_tg_matrix\.h5$", "", fn)
        for fn in os.listdir(BASE_IN)
        if fn.endswith(".re_tg_matrix.h5") and not fn.startswith("Vasc_")
    })
    if not cells:
        raise SystemExit("未找到任何 *.re_tg_matrix.h5")

    print("待处理细胞类型：", ", ".join(cells))
    for cell in cells:
        process_one_cell(cell)

if __name__ == "__main__":
    main()

