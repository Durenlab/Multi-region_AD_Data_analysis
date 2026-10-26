#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import time
import contextlib
from typing import Optional, Tuple

import numpy as np
import pandas as pd
import h5py
import torch
from scipy import sparse
from scipy.io import savemat
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--region", type=str, required=True)
parser.add_argument("--celltype", type=str, required=True)
parser.add_argument(
    "--input-root", type=str,
    default=(os.environ["REGA_SEPARATE_ROOT"] + '/01_prepare/data/rega_inputs'),
    help="Directory containing REGION/CELLTYPE REGA input files",
)
parser.add_argument(
    "--output-root", type=str,
    default=(os.environ["REGA_SEPARATE_ROOT"] + '/02_run_rega/rega_output'),
    help="Directory for REGION/CELLTYPE stage1.mat and stage2.mat",
)
args = parser.parse_args()

REGION = args.region
CELLTYPE = args.celltype
if CELLTYPE == "Vasc":
    parser.error("Vasc is excluded from the active REGA analysis")

# =========================
# 路径（按需修改）
# =========================
BASE = args.input_root

PATH_X_CSV = f"{BASE}/{REGION}/{CELLTYPE}/input_GEX.csv"
PATH_W_H5  = f"{BASE}/{REGION}/{CELLTYPE}/re_tg_matrix.h5"
PATH_B_CSV = f"{BASE}/{REGION}/{CELLTYPE}/input_binding.csv"

OUTDIR = os.path.join(args.output_root, REGION, CELLTYPE)
os.makedirs(OUTDIR, exist_ok=True)

for p in [PATH_X_CSV, PATH_W_H5, PATH_B_CSV]:
    if not os.path.exists(p):
        raise FileNotFoundError(f"Missing input: {p}")

# =========================
# 训练与数值参数
# =========================
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
DTYPE  = torch.float32

K       = 20         # 主题/模块数
EPS     = 1e-10      # 数值稳定项（需要时可调大到 1e-8/1e-6）
USE_AMP = False      # 自动混合精度开关（需要 GPU/支持时可 True）

# 非负 L1 正则（非负约束下 sum == L1）；默认关闭
LAM_W = 0.0
LAM_B = 0.0
LAM_V = 0.0
LAM_H = 0.0

LOG_EVERY_INNER = 500

# 两阶段
INNER1, OUTER1 = 2000, 30
INNER2, OUTER2 = 500, 3000

# =========================
# 工具
# =========================
@contextlib.contextmanager
def timer(msg: str):
    t0 = time.time()
    yield
    print(f"[{msg}] {time.time()-t0:.3f}s")

def ensure_finite(a: np.ndarray, name: str) -> np.ndarray:
    bad = np.sum(~np.isfinite(a))
    if bad > 0:
        print(f"[sanitize] {name}: {bad} non-finite -> 0")
        a = np.nan_to_num(a, nan=0.0, posinf=0.0, neginf=0.0)
    return a

def load_X(csv_path: str) -> Tuple[np.ndarray, list, list]:
    df = pd.read_csv(csv_path, header=0, index_col=0)
    X = df.to_numpy(dtype=np.float32)
    X = ensure_finite(X, "X")
    return X, df.index.tolist(), df.columns.tolist()

def load_B(csv_path: str) -> np.ndarray:
    # 与 MATLAB readmatrix 后丢第一列一致：仅使用数值列
    df = pd.read_csv(csv_path, header=0, index_col=0)
    B = df.to_numpy(dtype=np.float32)
    return B

def load_W_h5_as_coo(h5_path: str) -> Tuple[sparse.coo_matrix, list, list]:
    with h5py.File(h5_path, "r") as f:
        rows = np.array(f["/i"]).ravel()
        cols = np.array(f["/j"]).ravel()
        vals = np.array(f["/x"]).ravel()
        dims = np.array(f["/dims"]).ravel()
        rnames = [s.decode("utf-8") if isinstance(s, (bytes, bytearray)) else str(s) for s in np.array(f["/rownames"])]
        cnames = [s.decode("utf-8") if isinstance(s, (bytes, bytearray)) else str(s) for s in np.array(f["/colnames"])]
    # 1-based -> 0-based 自动识别
    if rows.min() == 1 or cols.min() == 1 or rows.max() == dims[0] or cols.max() == dims[1]:
        rows -= 1; cols -= 1
    vals = ensure_finite(vals.astype(np.float32), "W.data")
    D = sparse.coo_matrix((vals, (rows, cols)), shape=(int(dims[0]), int(dims[1])))
    # 去除 0
    nz = D.data != 0
    if nz.sum() != D.nnz:
        D = sparse.coo_matrix((D.data[nz], (D.row[nz], D.col[nz])), shape=D.shape)
    return D, rnames, cnames

# 稀疏 W 的 torch 封装与乘法
@torch.no_grad()
def to_torch_sparse(D: sparse.coo_matrix, device, dtype):
    D = D.tocoo()
    idx = np.vstack([D.row, D.col]).astype(np.int64)
    idx = torch.from_numpy(idx).to(device)
    val = torch.from_numpy(D.data.astype(np.float32)).to(device=device, dtype=dtype)
    W_sp = torch.sparse_coo_tensor(idx, val, size=D.shape, device=device, dtype=dtype).coalesce()
    return W_sp

@torch.no_grad()
def spmm(W_sp: torch.Tensor, M_dn: torch.Tensor) -> torch.Tensor:
    return torch.sparse.mm(W_sp, M_dn)

@torch.no_grad()
def spmmT(W_sp: torch.Tensor, M_dn: torch.Tensor) -> torch.Tensor:
    idx = W_sp.indices()
    WT = torch.sparse_coo_tensor(idx.flip(0), W_sp.values(),
                                 (W_sp.shape[1], W_sp.shape[0]),
                                 device=W_sp.device, dtype=W_sp.dtype).coalesce()
    return torch.sparse.mm(WT, M_dn)

# 完整 Frobenius 损失（无 1/2）
@torch.no_grad()
def frob_loss_full(X: torch.Tensor, W_sp: torch.Tensor, B: torch.Tensor,
                   V: torch.Tensor, H: torch.Tensor,
                   lam_w=0.0, lam_b=0.0, lam_v=0.0, lam_h=0.0, eps=1e-10) -> torch.Tensor:
    WB = spmm(W_sp, B)          # (d×p)
    diff = X - WB @ V @ H       # (d×n)
    loss = (diff * diff).sum()
    if lam_w: loss = loss + lam_w * W_sp.values().sum()
    if lam_b: loss = loss + lam_b * B.sum()
    if lam_v: loss = loss + lam_v * V.sum()
    if lam_h: loss = loss + lam_h * H.sum()
    return loss

# =========================
# 乘法更新（稀疏 W 版本，支持 V/H 热启动）
# =========================
@torch.no_grad()
def mu_full_sparseW(
    X_np: np.ndarray,             # (d×n)
    W_coo: sparse.coo_matrix,     # (d×m) 稀疏
    B_np: np.ndarray,             # (m×p)
    K: int,
    inner_iters: int,
    outer_iters: int,
    device: str = DEVICE,
    dtype: torch.dtype = DTYPE,
    use_amp: bool = USE_AMP,
    eps: float = EPS,
    lam_w: float = LAM_W,
    lam_b: float = LAM_B,
    lam_v: float = LAM_V,
    lam_h: float = LAM_H,
    log_every_inner: int = LOG_EVERY_INNER,
    V0_np: np.ndarray | None = None,   # 热启动 V
    H0_np: np.ndarray | None = None,   # 热启动 H
):
    # to torch
    X = torch.from_numpy(X_np).to(device=device, dtype=dtype)   # (d×n)
    B = torch.from_numpy(B_np).to(device=device, dtype=dtype)   # (m×p)
    W_sp = to_torch_sparse(W_coo, device=device, dtype=dtype)   # (d×m)

    d, n = X.shape
    m, p = B.shape
    assert W_sp.shape == (d, m)
    
    #gen = torch.Generator(device=device).manual_seed(0)
    #V = torch.rand(p, K, generator=gen, device=device, dtype=dtype)
    #H = torch.rand(K, n, generator=gen, device=device, dtype=dtype)
    # ----- init V, H -----
    if V0_np is None:
        gen = torch.Generator(device=device).manual_seed(123)
        V = torch.rand(p, K, generator=gen, device=device, dtype=dtype)
    else:
        V = torch.from_numpy(V0_np).to(device=device, dtype=dtype)

    if H0_np is None:
        gen = torch.Generator(device=device).manual_seed(1234)
        H = torch.rand(K, n, generator=gen, device=device, dtype=dtype)
    else:
        H = torch.from_numpy(H0_np).to(device=device, dtype=dtype)

    # 初始 loss
    init_loss = float(frob_loss_full(X, W_sp, B, V, H, lam_w, lam_b, lam_v, lam_h, eps).cpu())
    print(f"[Init] loss = {init_loss:.6e}")

    scaler = torch.cuda.amp.GradScaler(enabled=use_amp)  # 预留

    for it in range(1, outer_iters + 1):
        with timer(f"Outer {it:04d}"):
            with torch.cuda.amp.autocast(enabled=use_amp):
                Z   = V @ H            # (p×n)

            # ---- W 更新（仅非零位置；避免生成 dense d x m 矩阵）----
            with torch.cuda.amp.autocast(enabled=use_amp):
                BVH = B @ Z                         # (m x n)
                W_BVH = spmm(W_sp, BVH)             # (d x n)

            idx = W_sp.indices()
            rows, cols = idx[0], idx[1]
            old_vals = W_sp.values()

            # number of nonzero W entries
            nnz = old_vals.numel()

            # chunk size controls GPU memory.
            # If still OOM, reduce to 200000 or 100000.
            edge_chunk_size = 500000

            new_vals_chunks = []

            for s in range(0, nnz, edge_chunk_size):
                e = min(s + edge_chunk_size, nnz)

                r = rows[s:e]
                c = cols[s:e]

                # numerator for nonzero edges:
                # numerW[r, c] = dot(X[r, :], BVH[c, :])
                numer_nz = (X[r, :] * BVH[c, :]).sum(dim=1)

                # denominator for nonzero edges:
                # denomW[r, c] = dot(W_BVH[r, :], BVH[c, :]) + lam_w
                denom_nz = (W_BVH[r, :] * BVH[c, :]).sum(dim=1) + lam_w + eps

                new_chunk = torch.clamp(
                    old_vals[s:e] * (numer_nz / denom_nz),
                    min=0.0
                )

                new_vals_chunks.append(new_chunk)

            new_vals = torch.cat(new_vals_chunks)

            W_sp = torch.sparse_coo_tensor(
                idx,
                new_vals,
                W_sp.shape,
                device=W_sp.device,
                dtype=W_sp.dtype
            ).coalesce()

            # ---- B 更新 ----
            with torch.cuda.amp.autocast(enabled=use_amp):
                WB      = spmm(W_sp, B)                # (d×p)
                VH      = V @ H                        # (p×n)
                numer_B = spmmT(W_sp, X @ VH.t())      # (m×p)
                WBVH    = WB @ VH                      # (d×n)
                denom_B = spmmT(W_sp, WBVH @ VH.t()) + lam_b   # (m×p)
                B.mul_(numer_B / (denom_B + eps)).clamp_(min=0.0)

            # ---- V/H 内层 ----
            for inner in range(1, inner_iters + 1):
                with torch.cuda.amp.autocast(enabled=use_amp):
                    WB = spmm(W_sp, B)                 # (d×p)
                    S1 = (WB.t() @ X)                  # (p×n)
                    S2 = (WB.t() @ WB)                 # (p×p)

                    HHt = H @ H.t()
                    V.mul_((S1 @ H.t()) / (S2 @ V @ HHt + lam_v + eps)).clamp_(min=0.0)
                    G = V.t() @ S2 @ V
                    H.mul_((V.t() @ S1) / (G @ H + lam_h + eps)).clamp_(min=0.0)

                if inner % log_every_inner == 0:
                    with torch.no_grad():
                        loss_val = float(frob_loss_full(X, W_sp, B, V, H, lam_w, lam_b, lam_v, lam_h, eps).cpu())
                    print(f"  [Inner {inner}/{inner_iters}] loss={loss_val:.6e}")

        with torch.no_grad():
            loss_val = float(frob_loss_full(X, W_sp, B, V, H, lam_w, lam_b, lam_v, lam_h, eps).cpu())
        print(f"[Outer {it:04d}] total loss = {loss_val:.6e}")

    # 导出 numpy/scipy
    W_idx = W_sp.indices().to("cpu").numpy()
    W_val = W_sp.values().to("cpu").numpy()
    W_out = sparse.coo_matrix((W_val, (W_idx[0], W_idx[1])), shape=W_sp.shape)

    return W_out, B.detach().cpu().numpy(), V.detach().cpu().numpy(), H.detach().cpu().numpy()

# =========================
# 主流程
# =========================
def main():
    print(f"[Device] {DEVICE}, dtype={DTYPE}, amp={USE_AMP}")

    with timer("Load X"):
        X_np, x_rows, x_cols = load_X(PATH_X_CSV)
    with timer("Load W (H5)"):
        W_coo, w_rows, w_cols = load_W_h5_as_coo(PATH_W_H5)
    with timer("Load B"):
        B_np = load_B(PATH_B_CSV)

    d, n = X_np.shape
    m, p = B_np.shape
    assert W_coo.shape == (d, m), f"W shape {W_coo.shape} != ({d},{m})"
    print(f"[Shapes] X={X_np.shape}  W={W_coo.shape}  B={B_np.shape}  -> V=(p×K)=({p},{K}), H=(K×n)=({K},{n})")

    print("\n========== Phase 1: (inner=2000, outer=30) ==========\n")
    W1, B1, V1, H1 = mu_full_sparseW(
        X_np=X_np, W_coo=W_coo, B_np=B_np,
        K=K, inner_iters=INNER1, outer_iters=OUTER1,
        device=DEVICE, dtype=DTYPE, use_amp=USE_AMP, eps=EPS,
        lam_w=LAM_W, lam_b=LAM_B, lam_v=LAM_V, lam_h=LAM_H,
    )
    savemat(os.path.join(OUTDIR, "stage1.mat"), {
        "W": W1.tocsr(), "B": B1.astype(np.float32), "V": V1.astype(np.float32), "H": H1.astype(np.float32),
        "X_rows": np.array(x_rows, dtype=object), "X_cols": np.array(x_cols, dtype=object),
        "W_rows": np.array(w_rows, dtype=object), "W_cols": np.array(w_cols, dtype=object),
    })
    print(f"[SAVE] stage1.mat -> {OUTDIR}")

    print("\n========== Phase 2: (inner=500, outer=3000) ==========\n")
    # 用阶段1的结果热启动 V/H；W 用阶段1的稀疏和值；B 从阶段1开始继续
    W1 = W1.tocoo()
    W2, B2, V2, H2 = mu_full_sparseW(
        X_np=X_np, W_coo=W1, B_np=B1,
        K=K, inner_iters=INNER2, outer_iters=OUTER2,
        device=DEVICE, dtype=DTYPE, use_amp=USE_AMP, eps=EPS,
        lam_w=LAM_W, lam_b=LAM_B, lam_v=LAM_V, lam_h=LAM_H,
        V0_np=V1, H0_np=H1,     # <<<< 关键：热启动
    )
    savemat(os.path.join(OUTDIR, "stage2.mat"), {
        "W": W2.tocsr(), "B": B2.astype(np.float32), "V": V2.astype(np.float32), "H": H2.astype(np.float32),
        "X_rows": np.array(x_rows, dtype=object), "X_cols": np.array(x_cols, dtype=object),
        "W_rows": np.array(w_rows, dtype=object), "W_cols": np.array(w_cols, dtype=object),
    })
    print(f"[SAVE] stage2.mat -> {OUTDIR}")

if __name__ == "__main__":
    main()
