#!/usr/bin/env python3
from pathlib import Path
import os
SHARED_ROOT = Path(os.environ["REGA_SHARED_ROOT"]) if "REGA_SHARED_ROOT" in os.environ else next(
    p for p in Path(__file__).resolve().parents if p.name in {"REGA_shared", "shared"}
)
SEPARATE_ROOT = Path(os.environ.get("REGA_SEPARATE_ROOT", "/project/zduren/durenlab/naqing_palmetto/REGA_separate"))
# -*- coding: utf-8 -*-

import os, glob, time, contextlib, argparse
from typing import List, Optional, Tuple, Dict

import numpy as np
import torch
import pandas as pd

# ---------------- Utils ----------------
@contextlib.contextmanager
def timer(msg: str):
    t0 = time.time(); yield; print(f"[{msg}] {time.time()-t0:.3f}s")

def to_str_list(x) -> List[str]:
    arr = list(x)
    out = []
    for v in arr:
        if isinstance(v, (bytes, np.bytes_)):
            out.append(v.decode())
        else:
            out.append(str(v))
    return out

def check_all_equal(name: str, items: List[List[str]]) -> List[str]:
    ref = items[0]
    for i, it in enumerate(items[1:], 1):
        if it != ref:
            raise ValueError(f"[Mismatch] {name} not identical across blocks (first!=#{i})")
    return ref

# ---------------- Loss ----------------
def frob_loss_blocks(
    X_blocks: List[torch.Tensor],
    W_blocks: List[torch.Tensor],
    B_blocks: List[torch.Tensor],
    V_cells: Dict[str, torch.Tensor],   # cell -> (p,K)
    H: torch.Tensor,                    # (K,n) shared
    block_cells: List[str],
    lam_w=0.0, lam_b=0.0, lam_v=0.0, lam_h=0.0, eps=1e-10
) -> torch.Tensor:
    dev, dt = H.device, H.dtype
    loss = torch.tensor(0.0, device=dev, dtype=dt)
    for Xc, Wc, Bc, cell in zip(X_blocks, W_blocks, B_blocks, block_cells):
        Zc = V_cells[cell] @ H
        diff = Xc - Wc @ (Bc @ Zc)
        loss = loss + (diff*diff).sum()
    if lam_w: loss = loss + lam_w * sum(w.sum() for w in W_blocks)
    if lam_b: loss = loss + lam_b * sum(b.sum() for b in B_blocks)
    if lam_v: loss = loss + lam_v * sum(V_cells[c].sum() for c in V_cells)
    if lam_h: loss = loss + lam_h * H.sum()
    return loss

# ---------------- Training ----------------
def mu_train_shareH(
    X_blocks: List[torch.Tensor],
    W_blocks: List[torch.Tensor],
    B_blocks: List[torch.Tensor],
    V_cells: Dict[str, torch.Tensor],   # per-cell V
    H: torch.Tensor,                    # shared H
    block_cells: List[str],
    *,
    inner_iters: int,
    outer_iters: int,
    device: str,
    use_amp: bool,
    eps: float,
    clamp_min: float,
    log_every_inner: int,
    ckpt_every_outer: int,
    outdir: str,
    block_chrs: List[str],
    samples: List[str],
    motifs: List[str],
    lam_w: float = 0.0, lam_b: float = 0.0, lam_v: float = 0.0, lam_h: float = 0.0,
):
    os.makedirs(outdir, exist_ok=True)
    cells = sorted(set(block_cells))
    p, K = next(iter(V_cells.values())).shape
    n = H.shape[1]

    with torch.no_grad():
        init_loss = float(frob_loss_blocks(
            X_blocks, W_blocks, B_blocks, V_cells, H, block_cells,
            lam_w, lam_b, lam_v, lam_h, eps
        ).cpu())
    print(f"[Init] loss={init_loss:.6e} | blocks={len(X_blocks)} | cells={len(cells)} | p={p} | K={K} | n={n}")

    for it in range(1, outer_iters+1):
        with timer(f"Outer {it:04d}"):
            # --- (A) 用当前 V_c,H 更新所有块的 W,B 一次 ---
            for Xc, Wc, Bc, cell in zip(X_blocks, W_blocks, B_blocks, block_cells):
                with torch.amp.autocast(device_type="cuda", enabled=(use_amp and device == "cuda")):
                    Zc  = V_cells[cell] @ H               # (p,n)
                    Yc  = Bc @ Zc                          # (r_i,n)
                    numW = Xc @ Yc.t()                     # (g_i,r_i)
                    denW = (Wc @ (Yc @ Yc.t())) + lam_w + eps
                    Wc.mul_(numW / denW).clamp_(min=clamp_min)

                    numB = (Wc.t() @ Xc) @ Zc.t()          # (r_i,p)
                    denB = ((Wc.t() @ Wc) @ Bc @ (Zc @ Zc.t())) + lam_b + eps
                    Bc.mul_(numB / denB).clamp_(min=clamp_min)

            # --- (B) 固定 W,B，聚合 S1_c,S2_c（只依赖 W,B,X；外层一次计算即可）---
            S1_c = {c: torch.zeros(p, n, device=device, dtype=H.dtype) for c in cells}
            S2_c = {c: torch.zeros(p, p, device=device, dtype=H.dtype) for c in cells}
            for Xc, Wc, Bc, cell in zip(X_blocks, W_blocks, B_blocks, block_cells):
                S1_c[cell] += Bc.t() @ (Wc.t() @ Xc)      # (p,n)
                S2_c[cell] += Bc.t() @ (Wc.t() @ Wc) @ Bc # (p,p)

            # --- (C) 内层：在固定 S1_c,S2_c 下，交替 inner_iters 次更新 V_c 与共享 H ---
            for inner in range(1, inner_iters + 1):
                with torch.amp.autocast(device_type="cuda", enabled=(use_amp and device == "cuda")):
                    # 更新各细胞的 V_c
                    HHt = H @ H.t()                        # (K,K)
                    for cell in cells:
                        Vc = V_cells[cell]
                        numV = S1_c[cell] @ H.t()          # (p,K)
                        denV = (S2_c[cell] @ Vc @ HHt) + lam_v + eps
                        Vc.mul_(numV / denV).clamp_(min=clamp_min)

                    # 更新共享 H
                    numH = torch.zeros(K, n, device=device, dtype=H.dtype)
                    denG = torch.zeros(K, K, device=device, dtype=H.dtype)
                    for cell in cells:
                        Vc = V_cells[cell]
                        numH += Vc.t() @ S1_c[cell]        # (K,n)
                        denG += Vc.t() @ S2_c[cell] @ Vc   # (K,K)
                    denH = (denG @ H) + lam_h + eps
                    H.mul_(numH / denH).clamp_(min=clamp_min)

                if (inner % log_every_inner) == 0:
                    with torch.no_grad():
                        l = float(frob_loss_blocks(
                            X_blocks, W_blocks, B_blocks, V_cells, H, block_cells,
                            lam_w, lam_b, lam_v, lam_h, eps
                        ).cpu())
                    print(f"  [Inner {inner}/{inner_iters}] loss={l:.6e}")

        with torch.no_grad():
            loss_val = float(frob_loss_blocks(
                X_blocks, W_blocks, B_blocks, V_cells, H, block_cells,
                lam_w, lam_b, lam_v, lam_h, eps
            ).cpu())
        print(f"[Outer {it:04d}] total loss = {loss_val:.6e}")

        if (it % ckpt_every_outer) == 0 or it == outer_iters:
            ck = os.path.join(outdir, f"mu_blocks_outer{it:04d}.pt")
            torch.save({
                "step": it,
                "W_blocks": [w.detach().cpu() for w in W_blocks],
                "B_blocks": [b.detach().cpu() for b in B_blocks],
                "V_cells": {c: V_cells[c].detach().cpu() for c in cells},
                "H": H.detach().cpu(),
                "loss": float(loss_val),
                "block_cells": block_cells,
                "block_chrs": block_chrs,
                "samples": samples,
                "motifs": motifs,
                "cells": cells,
                "inner_iters": inner_iters,
            }, ck)
            print(f"[CKPT] saved -> {ck}")

    return V_cells, H

# ---------------- Data loading ----------------
def load_all_blocks_aligned(
    blocks_root: str,
    cells_filter: Optional[List[str]] = None,
    include_chrs: Optional[List[str]] = None,
):
    dirs = sorted(d for d in glob.glob(os.path.join(blocks_root, "*_blocks_aligned")) if os.path.isdir(d))
    if cells_filter:
        keep = set(cells_filter)
        dirs = [d for d in dirs if os.path.basename(d).replace("_blocks_aligned","") in keep]
    if not dirs:
        raise SystemExit(f"No *_blocks_aligned in {blocks_root}")

    Xs, Ws, Bs = [], [], []
    block_cells, block_chrs = [], []
    samples_list, motifs_list = [], []

    for d in dirs:
        cell = os.path.basename(d).replace("_blocks_aligned","")
        npzs = sorted(glob.glob(os.path.join(d, "chr*.npz")))
        if not npzs:
            print(f"[WARN] {cell}: no chr*.npz"); continue
        for f in npzs:
            chrn = os.path.splitext(os.path.basename(f))[0]
            if include_chrs is not None and chrn not in include_chrs: continue
            z = np.load(f, allow_pickle=True)
            X = z["X"].astype(np.float32)
            W = z["W"].astype(np.float32)
            B = z["B"].astype(np.float32)
            if not np.isfinite(X).all(): X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
            if not np.isfinite(W).all(): W = np.nan_to_num(W, nan=0.0, posinf=0.0, neginf=0.0)
            if not np.isfinite(B).all(): B = np.nan_to_num(B, nan=0.0, posinf=0.0, neginf=0.0)
            print(
                f"[Block] cell={cell} chr={chrn} "
                f"X={X.shape} W={W.shape} B={B.shape}"
            )

            Xs.append(X); Ws.append(W); Bs.append(B)
            block_cells.append(cell); block_chrs.append(chrn)
            samples_list.append(to_str_list(z["samples"]))
            motifs_list.append(to_str_list(z["motifs"]))

    if not Xs: raise SystemExit("No valid blocks loaded")
    samples_ref = check_all_equal("samples", samples_list)
    motifs_ref  = check_all_equal("motifs",  motifs_list)

    n = Xs[0].shape[1]; p = Bs[0].shape[1]
    for i in range(len(Xs)):
        assert Xs[i].shape[1] == n
        assert Bs[i].shape[1] == p

    total_bytes = 0
    for X, W, B in zip(Xs, Ws, Bs):
        total_bytes += X.nbytes + W.nbytes + B.nbytes

    print(f"[Memory estimate] X/W/B numpy total = {total_bytes / 1024**3:.2f} GB")
    print("[Memory estimate] Actual GPU memory will be higher because training creates temporary matrices.")

    print(f"[Load] cells={sorted(set(block_cells))}; blocks={len(Xs)}; n={n}; p={p}")
    return Xs, Ws, Bs, block_cells, block_chrs, samples_ref, motifs_ref

# ---------------- Main ----------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blocks-root", default=str(SHARED_ROOT / "part1_prepare_REGA/05REGA_blocks_region_consistent"))
    ap.add_argument("--cells", default="", help="逗号分隔细胞子集；留空=全部")
    ap.add_argument("--include-chrs", default="", help="只用这些染色体，逗号分隔")
    ap.add_argument("--out-root", default=str(SHARED_ROOT / "part2_run_REGA/07REGA_region_shareH_output"))
    ap.add_argument("--run-name", default="REGA_shareH")
    ap.add_argument("--K", type=int, default=20)
    ap.add_argument("--inner1", type=int, default=2000)
    ap.add_argument("--outer1", type=int, default=30)
    ap.add_argument("--inner2", type=int, default=500)
    ap.add_argument("--outer2", type=int, default=3000)
    ap.add_argument("--eps", type=float, default=1e-10)
    ap.add_argument("--amp", action="store_true")
    ap.add_argument("--lam-w", type=float, default=0.0)
    ap.add_argument("--lam-b", type=float, default=0.0)
    ap.add_argument("--lam-v", type=float, default=0.0)
    ap.add_argument("--lam-h", type=float, default=0.0)
    args = ap.parse_args()

    include_chrs = [s.strip() for s in args.include_chrs.split(",") if s.strip()] or None
    cells_filter = [s.strip() for s in args.cells.split(",") if s.strip()] or None

    device = "cuda" if torch.cuda.is_available() else "cpu"
    dtype  = torch.float32
    
    print(f"[Device] {device}")
    if device == "cuda":
        print(f"[CUDA] device name: {torch.cuda.get_device_name(0)}")
        print(f"[CUDA] total memory GB: {torch.cuda.get_device_properties(0).total_memory / 1024**3:.2f}")    
    # Load
    X_np, W_np, B_np, block_cells, block_chrs, samples_ref, motifs_ref = load_all_blocks_aligned(
        args.blocks_root, cells_filter, include_chrs
    )
    X_blocks = [torch.tensor(x, device=device, dtype=dtype) for x in X_np]
    W_blocks = [torch.tensor(w, device=device, dtype=dtype) for w in W_np]
    B_blocks = [torch.tensor(b, device=device, dtype=dtype) for b in B_np]

    n = X_np[0].shape[1]
    p = B_np[0].shape[1]
    K = int(args.K)
    cells = sorted(set(block_cells))

    # Init V per cell; shared H
    gen = torch.Generator(device=device).manual_seed(123)
    V_cells = {c: torch.rand(p, K, generator=gen, device=device, dtype=dtype) for c in cells}
    gen2 = torch.Generator(device=device).manual_seed(1234)
    H = torch.rand(K, n, generator=gen2, device=device, dtype=dtype)

    out_base   = os.path.join(args.out_root, args.run_name)
    phase1_dir = os.path.join(out_base, "phase1")
    phase2_dir = os.path.join(out_base, "phase2")
    os.makedirs(phase1_dir, exist_ok=True)
    os.makedirs(phase2_dir, exist_ok=True)

    print("\n========== Phase 1 ==========\n")
    V_cells, H = mu_train_shareH(
        X_blocks, W_blocks, B_blocks, V_cells, H, block_cells,
        inner_iters=args.inner1, outer_iters=args.outer1,
        device=device, use_amp=args.amp, eps=args.eps, clamp_min=0.0,
        log_every_inner=500, ckpt_every_outer=10, outdir=phase1_dir,
        block_chrs=block_chrs, samples=samples_ref, motifs=motifs_ref,
        lam_w=args.lam_w, lam_b=args.lam_b, lam_v=args.lam_v, lam_h=args.lam_h
    )

    print("\n========== Phase 2 ==========\n")
    V_cells, H = mu_train_shareH(
        X_blocks, W_blocks, B_blocks, V_cells, H, block_cells,
        inner_iters=args.inner2, outer_iters=args.outer2,
        device=device, use_amp=args.amp, eps=args.eps, clamp_min=0.0,
        log_every_inner=500, ckpt_every_outer=1000, outdir=phase2_dir,
        block_chrs=block_chrs, samples=samples_ref, motifs=motifs_ref,
        lam_w=args.lam_w, lam_b=args.lam_b, lam_v=args.lam_v, lam_h=args.lam_h
    )

    with torch.no_grad():
        final_loss = float(frob_loss_blocks(
            X_blocks, W_blocks, B_blocks, V_cells, H, block_cells,
            lam_w=args.lam_w, lam_b=args.lam_b, lam_v=args.lam_v, lam_h=args.lam_h, eps=args.eps
        ).cpu())
    print(f"\n[Done] Final loss = {final_loss:.6e}")

    # Save final & labeled CSVs
    os.makedirs(out_base, exist_ok=True)
    torch.save({
        "V_cells": {c: V_cells[c].detach().cpu() for c in cells},
        "H": H.detach().cpu(),
        "samples": samples_ref,
        "motifs": motifs_ref,
        "cells": cells,
        "block_cells": block_cells,
        "block_chrs": block_chrs,
        "final_loss": final_loss,
        "inner1": args.inner1, "outer1": args.outer1,
        "inner2": args.inner2, "outer2": args.outer2
    }, os.path.join(out_base, "final.pt"))

    modules = [f"module{i}" for i in range(1, K+1)]
    pd.DataFrame(H.detach().cpu().numpy(), index=modules, columns=samples_ref)\
      .to_csv(os.path.join(out_base, "H_with_samples.csv"))

    for c in cells:
        Vc = V_cells[c].detach().cpu().numpy()
        pd.DataFrame(Vc, index=motifs_ref, columns=modules)\
          .to_csv(os.path.join(out_base, f"V_{c}_with_motifs.csv"))

    print("[SAVE] H_with_samples.csv & V_<cell>_with_motifs.csv @", out_base)

if __name__ == "__main__":
    main()
