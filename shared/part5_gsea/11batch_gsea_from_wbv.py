#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import argparse
import numpy as np
import pandas as pd

# try import gseapy lazily
try:
    import gseapy as gp
except Exception:
    gp = None

# ----------------- 默认 -----------------
DEF_OUT_ROOT = "./gsea_from_wbv"
DEF_GMT = None

# ----------------- FDR（BH） -----------------
def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    p = np.where(np.isfinite(p), p, 1.0)
    n = p.size
    order = np.argsort(p)
    p_sorted = p[order]
    q_sorted = p_sorted * n / (np.arange(1, n + 1))
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q = np.empty_like(q_sorted)
    q[order] = q_sorted
    return np.clip(q, 0.0, 1.0)

try:
    from statsmodels.stats.multitest import multipletests
    def fdr_vec(p):
        p = np.asarray(p, float)
        p = np.where(np.isfinite(p), p, 1.0)
        return multipletests(p, method="fdr_bh")[1]
except Exception:
    def fdr_vec(p):
        return bh_fdr(p)

# ----------------- 主逻辑 -----------------
def main():
    ap = argparse.ArgumentParser(description="Batch prerank GSEA for all modules in a WBV CSV (gene x module).")
    ap.add_argument("--wbv",        required=True, help="WBV CSV 文件：行=gene, 列=module")
    ap.add_argument("--out-root",   default=DEF_OUT_ROOT, help="GSEA 输出根目录")
    ap.add_argument("--gmt",        default=DEF_GMT, help="本地 GMT 文件路径；留空则使用 --library（需联网）")
    ap.add_argument("--library",    default="MSigDB_Hallmark_2020", help="当 --gmt 为空时使用的 gseapy 基因集名称（需要联网）")
    ap.add_argument("--min-size",   type=int, default=10,   help="geneset 最小大小")
    ap.add_argument("--max-size",   type=int, default=5000, help="geneset 最大大小")
    ap.add_argument("--perm",       type=int, default=1000, help="GSEA 置换次数（prerank）")
    ap.add_argument("--threads",    type=int, default=4,    help="gseapy 线程数")
    ap.add_argument("--upper",      action="store_true",    help="基因符号统一为大写并按大写去重（保留绝对值最大行）")
    ap.add_argument("--ranking",    choices=["raw","abs","zscore"], default="raw",
                    help="用于 ranking 的基因分数：raw=原始 WBV；abs=绝对值；zscore=按列 z-score")
    ap.add_argument("--cache-dir",  default="", help="重定向 gseapy 缓存目录（写入 XDG_CACHE_HOME）")
    ap.add_argument("--jitter",     action="store_true", help="对完全相同分数添加极小抖动以稳定排序")
    ap.add_argument("--alpha",      type=float, default=0.05, help="仅信息输出用，不再用作筛选")
    args = ap.parse_args()

    if args.cache_dir:
        os.environ["XDG_CACHE_HOME"] = args.cache_dir
        os.makedirs(args.cache_dir, exist_ok=True)

    os.makedirs(args.out_root, exist_ok=True)

    if gp is None:
        raise SystemExit("[ERR] gseapy 未安装或导入失败。请先安装 gseapy。")

    # 读 WBV：行 = gene，列 = module
    if not os.path.exists(args.wbv):
        raise SystemExit(f"[ERR] 找不到 WBV 文件: {args.wbv}")
    WBV = pd.read_csv(args.wbv, sep="\t", index_col=0)
    WBV.index = WBV.index.astype(str)
    WBV.columns = WBV.columns.astype(str)

    # 规范化：列 L1 normalize -> 行 L1 normalize（与原脚本一致）
    col_sum = WBV.abs().sum(axis=0)
    col_sum = col_sum.replace(0, np.nan)
    WBV = WBV.div(col_sum, axis=1)

    row_sum = WBV.abs().sum(axis=1)
    row_sum = row_sum.replace(0, np.nan)
    WBV = WBV.div(row_sum, axis=0)

    WBV = WBV.fillna(0)

    # 可选大写去重（保留绝对值最大的行）
    if args.upper:
        WBV = WBV.groupby(WBV.index.str.upper()).agg(lambda s: s.iloc[np.argmax(np.abs(s.values))])
        WBV.index.name = "gene"

    # 可选 zscore
    if args.ranking == "zscore":
        WBV = (WBV - WBV.mean(axis=0)) / (WBV.std(axis=0).replace(0, np.nan))
        WBV = WBV.replace([np.inf, -np.inf], np.nan)
        WBV = WBV.fillna(0)

    modules = WBV.columns.tolist()
    print(f"[INFO] 从 WBV 读取到 {len(modules)} 个模块，将对每个模块做 GSEA。")

    all_rows = []

    for mod in modules:
        print(f"[INFO] 处理模块: {mod}")
        scores = WBV[mod].astype(float)
        scores = scores[np.isfinite(scores)]
        if args.ranking == "abs":
            scores = scores.abs()

        # 跳过太少或全零
        if (scores.abs().sum() == 0) or (scores.shape[0] < max(args.min_size, 5)):
            print(f"[WARN] {mod} 可用基因过少或全零，跳过")
            continue

        # prerank 两列：gene, score（降序）
        rnk = scores.sort_values(ascending=False)
        rnk_df = pd.DataFrame({"gene": rnk.index, "score": rnk.values})

        # jitter（必要时）
        if args.jitter and rnk_df["score"].duplicated().any():
            rng = np.random.default_rng(123)
            rnk_df["score"] = rnk_df["score"].to_numpy() + rng.normal(0, 1e-12, size=len(rnk_df))

        mod_out = os.path.join(args.out_root, mod)
        os.makedirs(mod_out, exist_ok=True)
        rnk_path = os.path.join(mod_out, f"{mod}_prerank.rnk.tsv")
        rnk_df.to_csv(rnk_path, sep="\t", index=False)

        gene_sets_arg = args.gmt if args.gmt else args.library

        # 调用 gseapy.prerank
        try:
            pre_res = gp.prerank(
                rnk=rnk_df,
                gene_sets=gene_sets_arg,
                threads=args.threads,
                permutation_num=args.perm,
                min_size=args.min_size,
                max_size=args.max_size,
                outdir=mod_out,
                format="png",
                seed=1234,
                verbose=False,
            )
        except Exception as e:
            print(f"[ERROR] {mod} gseapy.prerank 运行失败: {e}")
            continue

        # 兼容不同版本的结果对象
        res = getattr(pre_res, "res2d", None)
        if res is None or not isinstance(res, pd.DataFrame) or res.empty:
            res = getattr(pre_res, "results", None)
            if res is None or not isinstance(res, pd.DataFrame) or res.empty:
                print(f"[WARN] {mod} GSEA 无结果")
                continue
        res = res.copy()

        # 列名标准化（与原脚本一致）
        lc2orig = {str(c).strip().lower(): c for c in res.columns}
        def pick(*cands):
            for k in cands:
                if k in lc2orig:
                    return lc2orig[k]
            return None

        col_term = pick("term","name","gs","geneset","geneset_name","pathway")
        col_nes  = pick("nes","normalized enrichment score","normalized enrichment score (nes)")
        col_p    = pick("pval","p-value","p_value","p value","p","pvalue")
        col_fdr  = pick("fdr","fdr q-val","fdr q value","q","qval","adj_pval","padj","fdr_bh","fdrq")
        col_es   = pick("es","enrichment score")
        col_size = pick("size","gs_size","setsize","matched_size")
        col_lead = pick("lead_genes","leading_edge","lead_genes (tags)")

        out = pd.DataFrame()
        if col_term is not None: out["term"] = res[col_term].astype(str)
        if col_nes  is not None: out["nes"]  = pd.to_numeric(res[col_nes], errors="coerce")
        if col_es   is not None: out["es"]   = pd.to_numeric(res[col_es],  errors="coerce")
        if col_size is not None: out["matched_size"] = pd.to_numeric(res[col_size], errors="coerce")
        if col_lead is not None: out["lead_genes"]   = res[col_lead].astype(str)
        if col_p    is not None: out["pval"] = pd.to_numeric(res[col_p],   errors="coerce")
        if col_fdr  is not None:
            out["fdr"] = pd.to_numeric(res[col_fdr], errors="coerce")
        elif "pval" in out.columns:
            out["fdr"] = fdr_vec(out["pval"].to_numpy())

        if "nes" in out.columns:
            out["_absNES"] = out["nes"].abs()

        # 保存完整结果
        res_path = os.path.join(mod_out, f"{mod}_gsea_results.csv")
        out.to_csv(res_path, index=False)

        # 汇总前若干条用于全局 summary
        sort_keys = [c for c in ["fdr","pval","_absNES"] if c in out.columns]
        if sort_keys:
            ascending = [(c != "_absNES") for c in sort_keys]
            res_sorted = out.sort_values(by=sort_keys, ascending=ascending)
        else:
            res_sorted = out

        head = res_sorted.head(min(50, len(res_sorted))).copy()
        head.insert(0, "module", mod)
        all_rows.append(head)

    # 写全局汇总
    if all_rows:
        summary = pd.concat(all_rows, axis=0, ignore_index=True)
        keep_cols = [c for c in ["module","term","nes","pval","fdr","es","matched_size","lead_genes","_absNES"] if c in summary.columns]
        summary = summary[keep_cols]
        out_all = os.path.join(args.out_root, "GSEA_summary_all_modules.csv")
        summary.to_csv(out_all, index=False)
        print(f"[SAVE] 汇总 -> {out_all}")
    else:
        print("[INFO] 无可写入的 GSEA 结果（可能 WBV/基因集不匹配或所有模块被跳过）")

if __name__ == "__main__":
    main()

