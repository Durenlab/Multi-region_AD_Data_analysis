#!/usr/bin/env python3
"""Combine completed pair-level EXP008 tables into review-friendly summaries."""
from pathlib import Path
import sys
import pandas as pd

DEFAULT_ROOT = Path(__file__).resolve().parents[1] / "results"
ROOT = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else DEFAULT_ROOT


def collect(relative_glob: str) -> pd.DataFrame:
    files = sorted(ROOT.glob(relative_glob))
    frames = [pd.read_csv(path, sep="\t") for path in files]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def skipped_pairs() -> set[str]:
    files = sorted((ROOT / "audit").glob("*_pair_status.tsv"))
    frames = [pd.read_csv(path, sep="\t", dtype=str) for path in files]
    if not frames:
        return set()
    status = pd.concat(frames, ignore_index=True)
    if not {"pair", "status"}.issubset(status.columns):
        raise SystemExit("Invalid pair-status audit schema; expected pair and status columns")
    return set(status.loc[status.status.eq("skipped_no_significant_driver"), "pair"].dropna())


def main():
    tf = collect("*/regression_results/*/*_TF_level_regression.tsv")
    tg = collect("*/regression_results/*/*_TG_level_regression.tsv.gz")
    targets = collect("*/target_genes/*/*_target_genes.tsv.gz")
    skipped = skipped_pairs()
    tf = tf.loc[~tf.pair.isin(skipped)] if not tf.empty and "pair" in tf else tf
    tg = tg.loc[~tg.pair.isin(skipped)] if not tg.empty and "pair" in tg else tg
    targets = targets.loc[~targets.pair.isin(skipped)] if not targets.empty and "pair" in targets else targets
    if tf.empty:
        raise SystemExit("No TF-level pair results found; ensure the pair array completed")
    summary = ROOT / "summary"
    summary.mkdir(parents=True, exist_ok=True)
    tf.to_csv(summary / "TF_level_regression_results.tsv", sep="\t", index=False)
    if not tg.empty:
        tg.to_csv(summary / "TG_level_regression_results.tsv.gz", sep="\t", index=False, compression="gzip")
    if not targets.empty:
        target_counts = targets.groupby(["pair", "TF", "mode"], as_index=False).agg(
            n_selected_genes=("target_gene", "nunique"),
            n_expression_matched=("expression_available", "sum"),
        )
        target_counts.to_csv(summary / "target_set_counts.tsv", sep="\t", index=False)
    key = ["pair", "TF"]
    wide = tf.pivot_table(index=key, columns="mode", values=["beta_pathology", "p_value", "n_selected_target_genes", "n_expression_matched_target_genes"], aggfunc="first")
    wide.columns = ["__".join(str(part) for part in col if str(part)) for col in wide.columns]
    wide.reset_index().to_csv(summary / "mode_comparison.tsv", sep="\t", index=False)
    print(f"Combined {len(tf)} TF-mode test rows and {len(tg)} target-gene test rows; excluded {len(skipped)} skipped pairs")


if __name__ == "__main__":
    main()
