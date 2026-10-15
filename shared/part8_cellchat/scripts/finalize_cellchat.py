#!/usr/bin/env python
"""Legacy full finalizer for EXP002 CellChat.

Reads existing donor CellChat probabilities and source-H5AD donor metadata;
does not rerun CellChat or modify source data. Outputs are written only below
EXP002_results/finalized/. Its original progression model/FDR settings are
superseded by retest_cellchat_linear_region_fdr.py. Do not use --overwrite to
recreate the current progression results without first updating this legacy
pipeline; it would restore the prior covariate-adjusted, sender/receiver-family
results.
"""
from __future__ import annotations

import base64
import hashlib
import html
import json
import os
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.stats import t as student_t


EXP = Path(__file__).resolve().parents[1]
OUT = EXP / "finalized"
SOURCE = Path(os.environ.get("CELLCHAT_SOURCE_H5AD", "/project/zduren/durenlab/naqing_palmetto/AD_REGA_region_wise/part8_cellchat/AD_scRNA_used.h5ad"))
REGIONS = ["AG", "EC", "HC", "MTC", "PFC", "TH"]
PATHS = ["nonAD", "earlyAD", "lateAD"]
STAGE = {"nonAD": 0.0, "earlyAD": 1.0, "lateAD": 2.0}
CELLTYPES = ["Exc", "Inh", "Ast", "Oli", "OPC", "Mic_Immune", "Vasc"]
MIN_CELLS = 20
MIN_GROUP_N = 3
MIN_RESID_DF = 5
BOOTSTRAPS = 1000
SEED = 20260923
CONTRASTS = {
    "earlyAD_minus_nonAD": ("nonAD", "earlyAD"),
    "lateAD_minus_nonAD": ("nonAD", "lateAD"),
    "lateAD_minus_earlyAD": ("earlyAD", "lateAD"),
}
COLORS = {"nonAD": "#4477AA", "earlyAD": "#EE9944", "lateAD": "#CC6677"}


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype_backend="numpy_nullable")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, na_rep="NA")


def hc3_matrix(X: np.ndarray, Y: np.ndarray, contrast: np.ndarray) -> dict:
    """Vectorized OLS contrast and HC3 SE for many response columns."""
    n, p = X.shape
    df = n - p
    if df < 1 or np.linalg.matrix_rank(X) != p:
        raise np.linalg.LinAlgError("rank deficient design")
    inv = np.linalg.inv(X.T @ X)
    beta = inv @ X.T @ Y
    resid = Y - X @ beta
    hat = np.einsum("ij,jk,ik->i", X, inv, X)
    one_minus_h = np.maximum(1.0 - hat, np.finfo(float).eps)
    influence = X @ inv @ contrast
    variance = (influence[:, None] ** 2 * (resid / one_minus_h[:, None]) ** 2).sum(axis=0)
    se = np.sqrt(np.maximum(variance, 0.0))
    estimate = contrast @ beta
    tstat = np.divide(estimate, se, out=np.full_like(estimate, np.nan), where=se > 0)
    pval = 2.0 * student_t.sf(np.abs(tstat), df)
    critical = student_t.ppf(0.975, df)
    return {"beta": beta, "estimate": estimate, "se": se, "p": pval,
            "low": estimate - critical * se, "high": estimate + critical * se,
            "df": df, "rank": p}


def add_covariates(meta: pd.DataFrame, base: np.ndarray, names: list[str],
                   notes: list[str]) -> tuple[np.ndarray, list[str], dict]:
    """Add complete, nonredundant covariates without sacrificing estimability."""
    X = base.copy()
    included: list[str] = []
    transforms: dict = {}
    for col in ("age_death", "Sex", "pmi"):
        if col not in meta:
            notes.append(f"{col}: unavailable")
            continue
        values = meta[col]
        if values.isna().any():
            notes.append(f"{col}: omitted because donor metadata are missing")
            continue
        if col in ("age_death", "pmi"):
            v = pd.to_numeric(values, errors="coerce").to_numpy(float)
            if not np.isfinite(v).all() or np.std(v) == 0:
                notes.append(f"{col}: omitted because non-finite or constant")
                continue
            mean, sd = float(v.mean()), float(v.std(ddof=0))
            block = ((v - mean) / sd)[:, None]
            block_names = [col]
            transforms[col] = {"mean": mean, "sd": sd}
        else:
            cat = values.astype(str)
            levels = sorted(cat.unique())
            if len(levels) < 2:
                notes.append("Sex: omitted because constant")
                continue
            block = np.column_stack([(cat == lev).to_numpy(float) for lev in levels[1:]])
            block_names = [f"Sex[{lev}]" for lev in levels[1:]]
            transforms["Sex"] = {"levels": levels, "reference": levels[0]}
        trial = np.column_stack([X, block])
        rank = np.linalg.matrix_rank(trial)
        if rank != trial.shape[1]:
            notes.append(f"{col}: omitted because redundant with pathology/previous covariates")
            continue
        if len(meta) - trial.shape[1] < MIN_RESID_DF:
            notes.append(f"{col}: omitted because fewer than {MIN_RESID_DF} residual df")
            continue
        X = trial
        included.extend(block_names)
        notes.append(f"{col}: included")
    names.extend(included)
    return X, included, transforms


def build_design(meta: pd.DataFrame, categorical: bool = False) -> tuple[np.ndarray, list[str], list[str], dict]:
    notes: list[str] = []
    if categorical:
        path = meta["pathology"].astype(str)
        base = np.column_stack([
            np.ones(len(meta)),
            (path == "earlyAD").to_numpy(float),
            (path == "lateAD").to_numpy(float),
        ])
        names = ["intercept", "pathology_earlyAD", "pathology_lateAD"]
    else:
        base = np.column_stack([np.ones(len(meta)), meta["pathology"].map(STAGE).to_numpy(float)])
        names = ["intercept", "pathology_stage"]
    X, included, transforms = add_covariates(meta, base, names, notes)
    return X, names, notes, transforms


def feature_counts(meta: pd.DataFrame) -> dict[str, int]:
    return {p: int((meta.pathology == p).sum()) for p in PATHS}


def main() -> None:
    overwrite = "--overwrite" in sys.argv[1:]
    if (OUT / "config" / "analysis_config.json").exists() and not overwrite:
        raise SystemExit(f"Refusing to overwrite nonempty output directory: {OUT}")
    for sub in ("tables", "figures", "methods", "report", "config", "logs"):
        (OUT / sub).mkdir(parents=True, exist_ok=True)

    selected_path = EXP / "tables" / "main_figure_selected_features.tsv"
    selected = read_tsv(selected_path)
    selected = selected[["region", "source", "target", "interaction_id"]].copy()
    selected["region"] = selected.region.astype(str)
    selected["source"] = selected.source.astype(str)
    selected["target"] = selected.target.astype(str)
    selected["interaction_id"] = selected.interaction_id.astype(str)
    if selected.duplicated(["region", "source", "target", "interaction_id"]).any():
        raise ValueError("Fixed selected-feature list has duplicate rows")

    # Read only donor-level metadata columns in backed mode; expression matrix is not loaded.
    h5 = ad.read_h5ad(SOURCE, backed="r")
    obs = h5.obs[["ROSMAP_IndividualID", "Pathology", "age_death", "Sex", "pmi"]].copy()
    h5.file.close()
    obs = obs.rename(columns={"ROSMAP_IndividualID": "donor_id", "Pathology": "pathology"})
    obs["donor_id"] = obs.donor_id.astype("string").str.strip()
    obs["pathology"] = obs.pathology.astype("string").str.strip()
    donor_meta = obs.drop_duplicates().copy()
    duplicate_conflicts = donor_meta.groupby("donor_id", dropna=False).size()
    bad_ids = duplicate_conflicts[duplicate_conflicts > 1]
    if len(bad_ids):
        # Verify whether conflicts are only representation artifacts; refuse ambiguity.
        raise ValueError(f"Conflicting source covariates/pathology for {len(bad_ids)} donor IDs")
    donor_meta = donor_meta.drop_duplicates("donor_id").reset_index(drop=True)
    donor_meta["age_death"] = pd.to_numeric(donor_meta.age_death, errors="coerce")
    donor_meta["pmi"] = pd.to_numeric(donor_meta.pmi, errors="coerce")
    donor_meta["Sex"] = donor_meta.Sex.astype("string").str.strip()
    donor_meta["pathology"] = donor_meta.pathology.astype(str)
    donor_meta["donor_id"] = donor_meta.donor_id.astype(str)
    write_tsv(donor_meta, OUT / "tables" / "donor_metadata_used.tsv")

    fixed_features = set(map(tuple, selected[["region", "source", "target", "interaction_id"]].itertuples(index=False, name=None)))
    progression_parts, contrast_parts, family_counts, coverage_rows = [], [], [], []
    selected_scores: list[dict] = []
    candidate_meta_all: dict[str, pd.DataFrame] = {}
    selected_models: dict[tuple, dict] = {}
    pair_coverage: list[dict] = []
    network_rows, regional_counts = [], []

    input_paths: list[Path] = [SOURCE, selected_path]
    for region in REGIONS:
        score_path = EXP / "intermediate" / f"{region}_donor_nonzero_probabilities.tsv.gz"
        cov_path = EXP / "tables" / f"{region}_donor_celltype_coverage.tsv"
        cand_path = EXP / "tables" / f"{region}_candidate_interactions.tsv"
        status_path = EXP / "tables" / "donor_network_status.tsv"
        for path in (score_path, cov_path, cand_path, status_path):
            if not path.exists():
                raise FileNotFoundError(path)
        input_paths.extend((score_path, cov_path, cand_path))
        score = read_tsv(score_path)
        cov = read_tsv(cov_path)
        cand = read_tsv(cand_path)
        if cand.interaction_id.astype(str).duplicated().any():
            raise ValueError(f"Duplicate candidate interaction IDs in {region}")
        cand["interaction_id"] = cand.interaction_id.astype(str)
        candidate_meta_all[region] = cand
        for col in ("region", "pathology", "donor_id", "celltype", "donor_network_status"):
            cov[col] = cov[col].astype(str)
        cov["n_cells"] = pd.to_numeric(cov.n_cells, errors="coerce").fillna(0).astype(int)
        status = cov[["region", "pathology", "donor_id", "donor_network_status", "failure_message"]].drop_duplicates()
        network_rows.append(status)
        repeated = status.groupby("donor_id").pathology.nunique()
        if (repeated > 1).any():
            raise ValueError(f"Donor pathology disagreement across cell types in {region}")
        # Every donor × type appears once. Missing celltype rows are not treated as zero cells.
        if cov.duplicated(["donor_id", "celltype"]).any():
            raise ValueError(f"Duplicate donor-celltype coverage rows in {region}")
        donors = status[["region", "pathology", "donor_id", "donor_network_status"]].copy()
        donors = donors.merge(donor_meta, on="donor_id", how="left", suffixes=("", "_source"), validate="one_to_one")
        mismatch = donors.pathology.ne(donors.pathology_source)
        if mismatch.any():
            raise ValueError(f"Pathology differs between CellChat coverage and source metadata in {region}: {int(mismatch.sum())} donors")
        if donors[["age_death", "Sex", "pmi"]].isna().all(axis=None):
            raise ValueError("Source-H5AD donor covariates failed to map")
        source_cov = cov[["donor_id", "celltype", "n_cells", "status"]].copy()
        nmat = source_cov.pivot(index="donor_id", columns="celltype", values="n_cells").reindex(index=donors.donor_id, columns=CELLTYPES).fillna(0).astype(int)
        status_map = donors.set_index("donor_id").donor_network_status
        eligible_network = donors.donor_network_status.eq("completed").to_numpy()
        donor_idx = {d: i for i, d in enumerate(donors.donor_id)}
        cand_ids = cand.interaction_id.tolist()
        feature_idx = {f: i for i, f in enumerate(cand_ids)}
        score["donor_id"] = score.donor_id.astype(str)
        score["source"] = score.source.astype(str)
        score["target"] = score.target.astype(str)
        score["interaction_id"] = score.interaction_id.astype(str)
        score["communication_probability"] = pd.to_numeric(score.communication_probability, errors="coerce")
        if score.duplicated(["donor_id", "source", "target", "interaction_id"]).any():
            raise ValueError(f"Duplicate sparse probability entries in {region}")
        # Scores on failed or cell-ineligible donor pairs cannot be interpreted as measured zeros.
        selected_region = selected[selected.region == region]
        selected_keys = set(map(tuple, selected_region[["source", "target", "interaction_id"]].itertuples(index=False, name=None)))
        for row in selected_region.itertuples(index=False):
            pair_scores = score[(score.source == row.source) & (score.target == row.target) & (score.interaction_id == row.interaction_id)]
            for drow in donors.itertuples(index=False):
                ix = donor_idx[drow.donor_id]
                ok = bool(eligible_network[ix] and nmat.iloc[ix][row.source] > MIN_CELLS and nmat.iloc[ix][row.target] > MIN_CELLS)
                if not ok:
                    continue
                hit = pair_scores[pair_scores.donor_id == drow.donor_id]
                val = float(hit.communication_probability.iloc[0]) if len(hit) else 0.0
                selected_scores.append({"region": region, "donor_id": drow.donor_id,
                                        "pathology": drow.pathology, "source": row.source,
                                        "target": row.target, "interaction_id": row.interaction_id,
                                        "score": val, "score_row_present": bool(len(hit)),
                                        "score_status": "nonzero_recorded" if len(hit) else "evaluated_zero"})
        # Donor counts and failures for QC report.
        for path in PATHS:
            regional_counts.append({"region": region, "pathology": path,
                                    "n_donors_all": int((donors.pathology == path).sum()),
                                    "n_network_completed": int(((donors.pathology == path) & donors.donor_network_status.eq("completed")).sum()),
                                    "n_network_failed_or_skipped": int(((donors.pathology == path) & ~donors.donor_network_status.eq("completed")).sum())})
        # Candidate feature models are fit for each directed source-target pair.
        for source in CELLTYPES:
            for target in CELLTYPES:
                elig = eligible_network & (nmat[source].to_numpy() > MIN_CELLS) & (nmat[target].to_numpy() > MIN_CELLS)
                dm = donors.loc[elig].copy().reset_index(drop=True)
                group_n = feature_counts(dm)
                available = {p: int((dm.pathology == p).sum()) for p in PATHS}
                coverage_record = {"region": region, "source": source, "target": target,
                                   **{f"n_{p}": available[p] for p in PATHS},
                                   "n_eligible_donors": len(dm),
                                   "n_candidate_interactions": len(cand_ids),
                                   "n_eligible_donor_interaction_slots": len(dm) * len(cand_ids),
                                   "n_nonzero_score_rows": 0, "n_evaluated_zero_score_slots": 0,
                                   "n_excluded_missing_donor_interaction_slots": 0,
                                   "n_excluded_failed_or_skipped_donors": 0,
                                   "n_excluded_cell_count_ineligible_donors": 0}
                pair_coverage.append(coverage_record)
                if len(dm) == 0:
                    continue
                # Matrix columns are the complete pathology-blind candidate universe. Sparse absence is zero only here,
                # after restricting to completed networks and donor pairs that passed the >20-cell rule.
                Y = np.zeros((len(dm), len(cand_ids)), dtype=np.float64)
                ss = score[(score.source == source) & (score.target == target)]
                if len(ss):
                    good = ss[ss.donor_id.isin(set(dm.donor_id)) & ss.interaction_id.isin(feature_idx)]
                    rr = good.donor_id.map({d: i for i, d in enumerate(dm.donor_id)}).to_numpy(int)
                    cc = good.interaction_id.map(feature_idx).to_numpy(int)
                    vv = good.communication_probability.to_numpy(float)
                    if not np.isfinite(vv).all():
                        raise ValueError(f"Non-finite probability in {region}/{source}/{target}")
                    if (vv < 0).any():
                        raise ValueError(f"Negative communication probability in {region}/{source}/{target}")
                    Y[rr, cc] = vv
                    n_nonzero = int(np.count_nonzero(vv))
                else:
                    n_nonzero = 0
                excluded_network = int((~eligible_network).sum())
                excluded_cell_count = int((eligible_network & ~((nmat[source].to_numpy() > MIN_CELLS) & (nmat[target].to_numpy() > MIN_CELLS))).sum())
                coverage_record.update({
                    "n_nonzero_score_rows": n_nonzero,
                    "n_evaluated_zero_score_slots": len(dm) * len(cand_ids) - n_nonzero,
                    "n_excluded_failed_or_skipped_donors": excluded_network,
                    "n_excluded_cell_count_ineligible_donors": excluded_cell_count,
                    "n_excluded_missing_donor_interaction_slots": (excluded_network + excluded_cell_count) * len(cand_ids),
                })
                enough = all(available[p] >= MIN_GROUP_N for p in PATHS)
                Xtrend, trend_names, trend_notes, transforms = build_design(dm)
                trend_ok = enough and np.linalg.matrix_rank(Xtrend) == Xtrend.shape[1] and len(dm) - Xtrend.shape[1] >= MIN_RESID_DF
                if trend_ok:
                    trend_fit = hc3_matrix(Xtrend, Y, np.eye(Xtrend.shape[1])[trend_names.index("pathology_stage")])
                    varying = np.nanmax(Y, axis=0) > np.nanmin(Y, axis=0)
                    valid_p = np.isfinite(trend_fit["p"]) & varying
                    for j, iid in enumerate(cand_ids):
                        status_value = "tested" if valid_p[j] else ("no_score_variation" if not varying[j] else "not_estimable")
                        rec = {"region": region, "source": source, "target": target, "interaction_id": iid,
                               **{f"n_{p}": available[p] for p in PATHS},
                               "mean_nonAD": float(Y[dm.pathology.eq("nonAD").to_numpy(), j].mean()) if available["nonAD"] else np.nan,
                               "mean_earlyAD": float(Y[dm.pathology.eq("earlyAD").to_numpy(), j].mean()) if available["earlyAD"] else np.nan,
                               "mean_lateAD": float(Y[dm.pathology.eq("lateAD").to_numpy(), j].mean()) if available["lateAD"] else np.nan,
                               "beta_pathology_stage": float(trend_fit["estimate"][j]) if varying[j] else np.nan,
                               "se_HC3": float(trend_fit["se"][j]) if varying[j] else np.nan,
                               "ci_low": float(trend_fit["low"][j]) if varying[j] else np.nan,
                               "ci_high": float(trend_fit["high"][j]) if varying[j] else np.nan,
                               "p_value": float(trend_fit["p"][j]) if valid_p[j] else np.nan,
                               "status": status_value, "design_rank": trend_fit["rank"],
                               "residual_df": trend_fit["df"], "covariates": ";".join(trend_names[2:]),
                               "covariate_decisions": "; ".join(trend_notes),
                               "model_formula": "score ~ pathology_stage" + (" + " + " + ".join(trend_names[2:]) if trend_names[2:] else "")}
                        # adjusted predictions at the mean/categorical reference covariate profile
                        b = trend_fit["beta"][:, j]
                        mean_row = Xtrend.mean(axis=0); mean_row[0] = 1.0
                        pred = []
                        for stage_value in (0.0, 1.0, 2.0):
                            xrow = mean_row.copy(); xrow[1] = stage_value
                            pred.append(float(xrow @ b))
                        rec["prediction_nonAD_mean_covars"] = pred[0]
                        rec["prediction_earlyAD_mean_covars"] = pred[1]
                        rec["prediction_lateAD_mean_covars"] = pred[2]
                        progression_parts.append(rec)
                        key = (region, source, target, iid)
                        if key in fixed_features:
                            selected_models[key] = {"trend_beta": b.copy(), "trend_names": trend_names,
                                                    "X": Xtrend.copy(), "dm": dm.copy(), "y": Y[:, j].copy(),
                                                    "transforms": transforms}
                else:
                    reason = "insufficient_group_coverage" if not enough else "not_estimable_design"
                    for j, iid in enumerate(cand_ids):
                        progression_parts.append({"region": region, "source": source, "target": target,
                                                  "interaction_id": iid,
                                                  **{f"n_{p}": available[p] for p in PATHS},
                                                  "mean_nonAD": np.nan, "mean_earlyAD": np.nan, "mean_lateAD": np.nan,
                                                  "beta_pathology_stage": np.nan, "se_HC3": np.nan, "ci_low": np.nan,
                                                  "ci_high": np.nan, "p_value": np.nan, "status": reason,
                                                  "design_rank": int(np.linalg.matrix_rank(Xtrend)),
                                                  "residual_df": len(dm) - int(np.linalg.matrix_rank(Xtrend)),
                                                  "covariates": ";".join(trend_names[2:]),
                                                  "covariate_decisions": "; ".join(trend_notes),
                                                  "model_formula": "score ~ pathology_stage"})
                # Pairwise categorical models, adjusted in separate contrast families.
                for cname, (g1, g2) in CONTRASTS.items():
                    use = dm.pathology.isin([g1, g2]).to_numpy()
                    dmp = dm.loc[use].reset_index(drop=True)
                    yp = Y[use, :]
                    n1, n2 = int((dmp.pathology == g1).sum()), int((dmp.pathology == g2).sum())
                    if n1 >= MIN_GROUP_N and n2 >= MIN_GROUP_N:
                        base = np.column_stack([np.ones(len(dmp)), (dmp.pathology == g2).to_numpy(float)])
                        pair_names = ["intercept", f"pathology_{g2}"]
                        pair_notes: list[str] = []
                        Xpair, _, _ = add_covariates(dmp, base, pair_names, pair_notes)
                        pair_rank = int(np.linalg.matrix_rank(Xpair))
                        pair_ok = pair_rank == Xpair.shape[1] and len(dmp) - pair_rank >= MIN_RESID_DF
                    else:
                        Xpair, pair_names, pair_notes, pair_ok = np.empty((len(dmp), 0)), [], [], False
                        pair_rank = 0
                    if pair_ok:
                        cvec = np.zeros(Xpair.shape[1]); cvec[1] = 1.0
                        fit = hc3_matrix(Xpair, yp, cvec)
                        varying = np.ptp(yp, axis=0) > 0
                        for j, iid in enumerate(cand_ids):
                            tested = varying[j] and np.isfinite(fit["p"][j])
                            contrast_parts.append({"region": region, "comparison": cname,
                                                   "condition_1": g1, "condition_2": g2,
                                                   "source": source, "target": target, "interaction_id": iid,
                                                   "n_condition_1": n1, "n_condition_2": n2,
                                                   "effect_condition2_minus_condition1": float(fit["estimate"][j]) if varying[j] else np.nan,
                                                   "se_HC3": float(fit["se"][j]) if varying[j] else np.nan,
                                                   "ci_low": float(fit["low"][j]) if varying[j] else np.nan,
                                                   "ci_high": float(fit["high"][j]) if varying[j] else np.nan,
                                                   "p_value": float(fit["p"][j]) if tested else np.nan,
                                                   "status": "tested" if tested else "no_score_variation",
                                                   "design_rank": fit["rank"], "residual_df": fit["df"],
                                                   "covariates": ";".join(pair_names[2:]),
                                                   "covariate_decisions": "; ".join(pair_notes)})
                    else:
                        reason = "insufficient_donors" if n1 < MIN_GROUP_N or n2 < MIN_GROUP_N else "not_estimable_design"
                        for iid in cand_ids:
                            contrast_parts.append({"region": region, "comparison": cname,
                                                   "condition_1": g1, "condition_2": g2,
                                                   "source": source, "target": target, "interaction_id": iid,
                                                   "n_condition_1": n1, "n_condition_2": n2,
                                                   "effect_condition2_minus_condition1": np.nan, "se_HC3": np.nan,
                                                   "ci_low": np.nan, "ci_high": np.nan, "p_value": np.nan,
                                                   "status": reason, "design_rank": pair_rank,
                                                   "residual_df": len(dmp) - pair_rank,
                                                   "covariates": ";".join(pair_names[2:]),
                                                   "covariate_decisions": "; ".join(pair_notes)})
                del Y
        print(f"Audited and modeled {region}", flush=True)

    progression = pd.DataFrame(progression_parts)
    contrasts = pd.DataFrame(contrast_parts)
    progression["q_value_within_family"] = np.nan
    for _, idx in progression.groupby(["region", "source", "target"], sort=False).groups.items():
        p = progression.loc[idx, "p_value"]
        ok = p.notna() & np.isfinite(pd.to_numeric(p, errors="coerce"))
        if ok.any():
            progression.loc[np.asarray(idx)[ok.to_numpy()], "q_value_within_family"] = pd.Series(p[ok].to_numpy()).rank(method="min").to_numpy()  # overwritten below by exact BH
            vals = p[ok].astype(float).to_numpy()
            order = np.argsort(vals)
            ranks = np.arange(1, len(vals) + 1)
            q = vals[order] * len(vals) / ranks
            q = np.minimum.accumulate(q[::-1])[::-1].clip(0, 1)
            back = np.empty_like(q); back[order] = q
            progression.loc[np.asarray(idx)[ok.to_numpy()], "q_value_within_family"] = back
        family_counts.append({"test": "pathology_trend", "region": progression.loc[idx[0], "region"],
                              "source": progression.loc[idx[0], "source"], "target": progression.loc[idx[0], "target"],
                              "n_eligible_candidates": len(idx), "n_finite_p": int(ok.sum())})
    contrasts["q_value_within_family"] = np.nan
    for key, idx in contrasts.groupby(["region", "comparison", "source", "target"], sort=False).groups.items():
        p = contrasts.loc[idx, "p_value"]
        ok = p.notna() & np.isfinite(pd.to_numeric(p, errors="coerce"))
        if ok.any():
            vals = p[ok].astype(float).to_numpy(); order = np.argsort(vals)
            q = vals[order] * len(vals) / np.arange(1, len(vals) + 1)
            q = np.minimum.accumulate(q[::-1])[::-1].clip(0, 1)
            back = np.empty_like(q); back[order] = q
            contrasts.loc[np.asarray(idx)[ok.to_numpy()], "q_value_within_family"] = back
        family_counts.append({"test": key[1], "region": key[0], "source": key[2], "target": key[3],
                              "n_eligible_candidates": len(idx), "n_finite_p": int(ok.sum())})

    write_tsv(progression, OUT / "tables" / "pathology_progression_tests.tsv")
    write_tsv(contrasts, OUT / "tables" / "pathology_pairwise_contrasts.tsv")
    families = pd.DataFrame(family_counts)
    families["n_q_lt_0_05"] = 0
    for family in families.itertuples(index=False):
        if family.test == "pathology_trend":
            subset = progression[(progression.region == family.region) & (progression.source == family.source) & (progression.target == family.target)]
        else:
            subset = contrasts[(contrasts.region == family.region) & (contrasts.source == family.source) & (contrasts.target == family.target) & (contrasts.comparison == family.test)]
        families.loc[(families.test == family.test) & (families.region == family.region) & (families.source == family.source) & (families.target == family.target), "n_q_lt_0_05"] = int((subset.q_value_within_family < 0.05).sum())
    write_tsv(families, OUT / "tables" / "fdr_family_counts.tsv")

    selected_scores_df = pd.DataFrame(selected_scores)
    if len(selected_scores_df) == 0:
        raise ValueError("No selected display feature has any eligible donor scores")
    selected_scores_df = selected_scores_df.merge(selected, on=["region", "source", "target", "interaction_id"], how="left", validate="many_to_one")
    selected = selected.merge(progression[["region", "source", "target", "interaction_id", "beta_pathology_stage", "se_HC3", "ci_low", "ci_high", "p_value", "q_value_within_family", "status", "covariates", "prediction_nonAD_mean_covars", "prediction_earlyAD_mean_covars", "prediction_lateAD_mean_covars"]], on=["region", "source", "target", "interaction_id"], how="left", validate="one_to_one")
    write_tsv(selected, OUT / "tables" / "fixed_selected_features_refreshed.tsv")
    write_tsv(selected_scores_df, OUT / "tables" / "fixed_selected_feature_donor_scores.tsv")

    # Bootstrap every selected feature and every family-specific trend hit. Donors, not cells, are resampled within stage.
    boot_targets = set(map(tuple, selected[["region", "source", "target", "interaction_id"]].itertuples(index=False, name=None)))
    hits = progression.loc[progression.q_value_within_family < 0.05, ["region", "source", "target", "interaction_id"]]
    boot_targets.update(map(tuple, hits.itertuples(index=False, name=None)))
    rng = np.random.default_rng(SEED)
    boot_rows = []
    for key in sorted(boot_targets):
        model = selected_models.get(key)
        if model is None:
            # Features with inadequate group coverage are intentionally not bootstrapped.
            boot_rows.append({"region": key[0], "source": key[1], "target": key[2], "interaction_id": key[3], "status": "not_bootstrapped_not_estimable"})
            continue
        dm, X, y = model["dm"], model["X"], model["y"]
        group_indices = [np.flatnonzero(dm.pathology.to_numpy() == p) for p in PATHS]
        if any(len(z) < MIN_GROUP_N for z in group_indices):
            boot_rows.append({"region": key[0], "source": key[1], "target": key[2], "interaction_id": key[3], "status": "not_bootstrapped_insufficient_groups"})
            continue
        slopes, categorical = [], {name: [] for name in CONTRASTS}
        Xcat, cat_names, _, _ = build_design(dm, categorical=True)
        cat_contrasts = {}
        for cname, (g1, g2) in CONTRASTS.items():
            c = np.zeros(Xcat.shape[1]); ix = {"nonAD": None, "earlyAD": cat_names.index("pathology_earlyAD"), "lateAD": cat_names.index("pathology_lateAD")}
            if ix[g2] is not None: c[ix[g2]] += 1
            if ix[g1] is not None: c[ix[g1]] -= 1
            cat_contrasts[cname] = c
        base_X = X
        for _ in range(BOOTSTRAPS):
            take = np.concatenate([rng.choice(ix, size=len(ix), replace=True) for ix in group_indices])
            xb, yb = base_X[take], y[take]
            coef, _, rank, _ = np.linalg.lstsq(xb, yb, rcond=None)
            if rank == xb.shape[1]:
                slopes.append(float(coef[1]))
            xcb = Xcat[take]
            bb, _, rank_cat, _ = np.linalg.lstsq(xcb, y[take], rcond=None)
            if rank_cat == xcb.shape[1]:
                for cname, c in cat_contrasts.items(): categorical[cname].append(float(c @ bb))
        rec = {"region": key[0], "source": key[1], "target": key[2], "interaction_id": key[3],
               "bootstrap_replicates": BOOTSTRAPS, "successful_trend_fits": len(slopes),
               "failed_or_rank_deficient_trend_fits": BOOTSTRAPS - len(slopes),
               "trend_ci_low": float(np.quantile(slopes, .025)) if slopes else np.nan,
               "trend_ci_high": float(np.quantile(slopes, .975)) if slopes else np.nan,
               "fraction_positive_trend": float(np.mean(np.asarray(slopes) > 0)) if slopes else np.nan,
               "status": "tested" if slopes else "no_successful_bootstrap_fit"}
        for cname, vals in categorical.items():
            rec[f"{cname}_successful_fits"] = len(vals)
            rec[f"{cname}_failed_or_rank_deficient_fits"] = BOOTSTRAPS - len(vals)
            rec[f"{cname}_ci_low"] = float(np.quantile(vals, .025)) if vals else np.nan
            rec[f"{cname}_ci_high"] = float(np.quantile(vals, .975)) if vals else np.nan
        boot_rows.append(rec)
    bootstrap = pd.DataFrame(boot_rows)
    write_tsv(bootstrap, OUT / "tables" / "stratified_donor_bootstrap.tsv")

    # Coverage tables: preserve all donors including failures; distinguish missing from evaluated zero.
    networks = pd.concat(network_rows, ignore_index=True).drop_duplicates()
    write_tsv(networks, OUT / "tables" / "network_status_audit.tsv")
    regional_counts_df = pd.DataFrame(regional_counts)
    write_tsv(regional_counts_df, OUT / "tables" / "regional_donor_network_counts.tsv")
    pair_coverage_df = pd.DataFrame(pair_coverage)
    # Compute counts of tested interactions after model fitting.
    tested_counts = progression.groupby(["region", "source", "target"]).agg(
        n_tested_interactions=("p_value", lambda x: int(x.notna().sum())),
        n_q_lt_0_05=("q_value_within_family", lambda x: int((x < .05).sum()))).reset_index()
    pair_coverage_df = pair_coverage_df.merge(tested_counts, on=["region", "source", "target"], how="left")
    write_tsv(pair_coverage_df, OUT / "tables" / "directed_pair_coverage.tsv")
    write_tsv(pair_coverage_df[["region", "source", "target", "n_candidate_interactions",
                                "n_eligible_donors", "n_eligible_donor_interaction_slots",
                                "n_nonzero_score_rows", "n_evaluated_zero_score_slots",
                                "n_excluded_missing_donor_interaction_slots",
                                "n_excluded_failed_or_skipped_donors",
                                "n_excluded_cell_count_ineligible_donors"]],
              OUT / "tables" / "score_zero_missing_audit.tsv")
    celltype_coverage = []
    for region in REGIONS:
        cov = read_tsv(EXP / "tables" / f"{region}_donor_celltype_coverage.tsv")
        cov["n_cells"] = pd.to_numeric(cov.n_cells, errors="coerce").fillna(0)
        celltype_coverage.append(cov)
    cellcov = pd.concat(celltype_coverage, ignore_index=True)
    cell_summary = cellcov.groupby(["region", "pathology", "celltype"], dropna=False).agg(
        n_donors_any=("donor_id", "nunique"),
        n_donors_gt20=("n_cells", lambda x: int((x > MIN_CELLS).sum())),
        n_donors_network_completed=("donor_network_status", lambda x: int((x == "completed").sum())),
        median_cells=("n_cells", "median")).reset_index()
    write_tsv(cell_summary, OUT / "tables" / "celltype_coverage_summary.tsv")

    # Figure 1: all fixed selected features with individual donors and adjusted model trend.
    fig, axes = plt.subplots(max(1, int(np.ceil(len(selected) / 3))), 3,
                             figsize=(15, max(4, 3.25 * int(np.ceil(len(selected) / 3)))), squeeze=False)
    score_lookup = selected_scores_df
    for ax, row in zip(axes.flat, selected.itertuples(index=False)):
        key = (row.region, row.source, row.target, row.interaction_id)
        d = score_lookup[(score_lookup.region == row.region) & (score_lookup.source == row.source) &
                         (score_lookup.target == row.target) & (score_lookup.interaction_id == row.interaction_id)]
        for x, path in enumerate(PATHS):
            vals = d[d.pathology == path].score.to_numpy(float)
            jitter = rng.normal(0, .055, len(vals))
            ax.scatter(x + jitter, vals, color=COLORS[path], alpha=.78, s=19, edgecolor="white", linewidth=.25, zorder=3)
            ax.text(x, .98, f"n={len(vals)}", transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=8)
        pred = [row.prediction_nonAD_mean_covars, row.prediction_earlyAD_mean_covars, row.prediction_lateAD_mean_covars]
        if np.isfinite(pred).all():
            ax.plot(range(3), pred, color="#333333", linewidth=1.6, marker="_", markersize=7, zorder=2)
        title = f"{row.region}: {row.source} → {row.target}\n{row.interaction_id}"
        qtxt = "q=NA" if pd.isna(row.q_value_within_family) else f"BH q={row.q_value_within_family:.3g}"
        btxt = "β=NA" if pd.isna(row.beta_pathology_stage) else f"β={row.beta_pathology_stage:.3g} [{row.ci_low:.3g}, {row.ci_high:.3g}]"
        ax.set_title(f"{title}\n{btxt}; {qtxt}", fontsize=8.3)
        ax.set_xticks(range(3), PATHS, rotation=20, ha="right")
        ax.set_ylabel("CellChat probability")
        ax.grid(axis="y", color="#dddddd", linewidth=.5)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes.flat[len(selected):]: ax.set_visible(False)
    fig.suptitle("Fixed selected CellChat features: donor probabilities and adjusted pathology trends", fontsize=14, y=1.002)
    fig.tight_layout()
    fig.savefig(OUT / "figures" / "selected_features.pdf", bbox_inches="tight")
    fig.savefig(OUT / "figures" / "selected_features.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Figure 2: directional ligand-receptor diagrams for every fixed selected row.
    diagrams = selected.merge(pd.concat([candidate_meta_all[r].assign(region=r) for r in REGIONS], ignore_index=True),
                              on=["region", "interaction_id"], how="left")
    write_tsv(diagrams, OUT / "tables" / "selected_interaction_diagram_metadata.tsv")
    fig, axes = plt.subplots(max(1, int(np.ceil(len(diagrams) / 2))), 2,
                             figsize=(14, max(4, 2.35 * int(np.ceil(len(diagrams) / 2)))), squeeze=False)
    for ax, row in zip(axes.flat, diagrams.itertuples(index=False)):
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
        ax.text(.5, .95, f"{row.region} | {row.source} → {row.target}", ha="center", va="top", weight="bold", fontsize=10)
        ax.text(.12, .55, row.source, ha="center", va="center", fontsize=9,
                bbox=dict(boxstyle="round,pad=.45", fc="#dceaf7", ec="#4477AA"))
        ax.text(.43, .55, str(row.ligand) if pd.notna(row.ligand) else "ligand", ha="center", va="center", fontsize=8,
                bbox=dict(boxstyle="round,pad=.4", fc="#f5f0d9", ec="#aa9944"))
        ax.text(.68, .55, str(row.receptor) if pd.notna(row.receptor) else "receptor", ha="center", va="center", fontsize=8,
                bbox=dict(boxstyle="round,pad=.4", fc="#f5f0d9", ec="#aa9944"))
        ax.text(.9, .55, row.target, ha="center", va="center", fontsize=9,
                bbox=dict(boxstyle="round,pad=.45", fc="#e7e1f2", ec="#775da6"))
        for a, b in ((.22, .33), (.54, .58), (.78, .81)):
            ax.annotate("", xy=(b, .55), xytext=(a, .55), arrowprops=dict(arrowstyle="-|>", color="#444444", lw=1.2))
        ax.text(.5, .16, f"{row.interaction_id}\nInferred CellChat communication hypothesis",
                ha="center", va="center", fontsize=7.5, color="#555555")
    for ax in axes.flat[len(diagrams):]: ax.set_visible(False)
    fig.suptitle("Directed sender–ligand–receptor–receiver paths for the fixed feature list", fontsize=14)
    fig.tight_layout()
    fig.savefig(OUT / "figures" / "selected_interaction_diagrams.pdf", bbox_inches="tight")
    fig.savefig(OUT / "figures" / "selected_interaction_diagrams.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Figure 3: overview includes all directed pairs; tiles summarize all LR trends in the family.
    overview = progression.groupby(["region", "source", "target"], dropna=False).agg(
        median_beta=("beta_pathology_stage", "median"),
        n_estimable=("p_value", lambda x: int(x.notna().sum())),
        n_q_lt_0_05=("q_value_within_family", lambda x: int((x < .05).sum()))).reset_index()
    overview["pair"] = overview.source + " → " + overview.target
    write_tsv(overview, OUT / "tables" / "all_pair_overview_summary.tsv")
    pairs_order = [f"{s} → {t}" for s in CELLTYPES for t in CELLTYPES]
    matrix = overview.pivot(index="pair", columns="region", values="median_beta").reindex(index=pairs_order, columns=REGIONS)
    fig, ax = plt.subplots(figsize=(11.5, 17))
    finite = matrix.to_numpy(float); lim = np.nanquantile(np.abs(finite), .98) if np.isfinite(finite).any() else 1
    lim = max(float(lim), 1e-8)
    image = ax.imshow(matrix.to_numpy(float), cmap="coolwarm", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_xticks(np.arange(len(REGIONS)), REGIONS); ax.set_yticks(np.arange(len(pairs_order)), pairs_order, fontsize=7)
    for r, pair in enumerate(pairs_order):
        for c, region in enumerate(REGIONS):
            rr = overview[(overview.region == region) & (overview.pair == pair)]
            if len(rr) == 0 or rr.iloc[0].n_estimable == 0:
                label = "n<3"
            else:
                v = int(rr.iloc[0].n_q_lt_0_05); n = int(rr.iloc[0].n_estimable)
                label = f"{v}/{n}*" if v else f"0/{n}"
            ax.text(c, r, label, ha="center", va="center", fontsize=5.4,
                    color="black" if np.isfinite(matrix.loc[pair, region]) else "#666666")
    ax.set_title("All directed cell-type pairs: median LR-interaction progression slope\nCell labels show within-family q<0.05 hits / estimable LR tests; n<3 = insufficient coverage")
    ax.set_xlabel("Brain region"); ax.set_ylabel("Directed sender → receiver pair")
    fig.colorbar(image, ax=ax, label="Median HC3-adjusted slope per pathology stage")
    fig.tight_layout()
    fig.savefig(OUT / "figures" / "all_pair_overview.pdf", bbox_inches="tight")
    fig.savefig(OUT / "figures" / "all_pair_overview.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Figure 4: pathology donor counts and eligibility for every pair.
    fig = plt.figure(figsize=(17, 19))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 5], hspace=.22)
    ax0 = fig.add_subplot(gs[0])
    bar_df = regional_counts_df
    x = np.arange(len(REGIONS)); width = .24
    for i, path in enumerate(PATHS):
        vals = [int(bar_df[(bar_df.region == r) & (bar_df.pathology == path)].n_donors_all.iloc[0]) for r in REGIONS]
        ax0.bar(x + (i-1)*width, vals, width, color=COLORS[path], label=path)
        for xx, v in zip(x + (i-1)*width, vals): ax0.text(xx, v+.7, str(v), ha="center", fontsize=7)
    ax0.set_xticks(x, REGIONS); ax0.set_ylabel("Donors"); ax0.set_title("Donor numbers by region and pathology (labels give all donor counts)")
    ax0.legend(ncol=3, frameon=False); ax0.spines[["top", "right"]].set_visible(False)
    ax1 = fig.add_subplot(gs[1])
    pair_order = [f"{s} → {t}" for s in CELLTYPES for t in CELLTYPES]
    cols = [(r, p) for r in REGIONS for p in PATHS]
    covmat = np.full((len(pair_order), len(cols)), np.nan)
    for i, pair in enumerate(pair_order):
        s, t = pair.split(" → ")
        for j, (r, p) in enumerate(cols):
            covrow = pair_coverage_df[(pair_coverage_df.region == r) & (pair_coverage_df.source == s) & (pair_coverage_df.target == t)]
            if len(covrow): covmat[i, j] = float(covrow.iloc[0][f"n_{p}"])
    im = ax1.imshow(covmat, cmap="viridis", aspect="auto", vmin=0, vmax=max(1, np.nanmax(covmat)))
    ax1.set_yticks(np.arange(len(pair_order)), pair_order, fontsize=6.2)
    ax1.set_xticks(np.arange(len(cols)), [f"{r}\n{p}" for r,p in cols], fontsize=7)
    ax1.set_title("Eligible donors per directed sender → receiver pair (must pass >20 cells and completed network)")
    for i in range(len(pair_order)):
        for j in range(len(cols)):
            val = covmat[i,j]
            ax1.text(j, i, "NA" if not np.isfinite(val) else str(int(val)), ha="center", va="center", fontsize=4.2,
                     color="white" if np.isfinite(val) and val < np.nanmax(covmat)*.7 else "black")
    fig.colorbar(im, ax=ax1, label="Eligible donor count")
    fig.savefig(OUT / "figures" / "coverage_qc.pdf", bbox_inches="tight")
    fig.savefig(OUT / "figures" / "coverage_qc.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Provenance and exact analysis configuration.
    checksums = []
    for path in dict.fromkeys(input_paths):
        checksums.append({"path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    config = {
        "run_utc": datetime.now(timezone.utc).isoformat(), "script": str(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()), "source_h5ad": str(SOURCE),
        "regions": REGIONS, "pathology_code": STAGE, "celltypes": CELLTYPES,
        "cell_group_eligibility": "more than 20 cells per donor x region x cell type and donor network completed",
        "donor_metadata_covariates": ["age_death", "Sex", "pmi"],
        "model": "donor-level OLS with HC3 standard errors; progression nonAD=0, earlyAD=1, lateAD=2",
        "minimum_donors_per_pathology": MIN_GROUP_N, "minimum_residual_df": MIN_RESID_DF,
        "covariate_selection_order": ["age_death", "Sex", "pmi"],
        "FDR": "BH separately within region x directed source-target x test; finite p-values only; q<0.05",
        "bootstrap": {"replicates": BOOTSTRAPS, "sampling": "donors with replacement within pathology groups, observed group sizes preserved", "seed": SEED},
        "CellChat": {"version": "2.2.0.9001", "database": "Human CellChatDB v1, as labeled in candidate metadata; region-specific LRsig from pooled object", "average_type": "triMean", "population_size": False, "within_donor_permutations": 100},
        "software": {"python": sys.version, "platform": platform.platform(), "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__, "anndata": ad.__version__},
        "checksums": checksums,
    }
    (OUT / "config" / "analysis_config.json").write_text(json.dumps(config, indent=2, default=str) + "\n")
    (OUT / "config" / "input_checksums.tsv").write_text(pd.DataFrame(checksums).to_csv(sep="\t", index=False))
    (OUT / "config" / "python_session.txt").write_text("\n".join(f"{k}: {v}" for k,v in config["software"].items()) + "\n")
    r_session = EXP / "intermediate" / "AG_sessionInfo.txt"
    if r_session.exists():
        (OUT / "config" / "CellChat_R_sessionInfo.txt").write_text(r_session.read_text())
    for source_name, output_name in [
        ("existing_pooled_celltype_differences_descriptive.tsv", "pooled_celltype_differences_descriptive.tsv"),
        ("existing_pooled_pathway_differences_descriptive.tsv", "pooled_pathway_differences_descriptive.tsv"),
    ]:
        source_path = EXP / "tables" / source_name
        if source_path.exists():
            (OUT / "tables" / output_name).write_bytes(source_path.read_bytes())

    # Methods, consistent with executed settings.
    method = """## CellChat analysis

CellChat communication was evaluated from the six region-specific normalized expression inputs derived from the ROSMAP single-cell RNA sequencing dataset. The donor was the biological replicate; communication features were stratified by brain region, sender cell type, receiver cell type, and ligand–receptor interaction. CellChat 2.2.0.9001 used the region-specific, pathology-blind ligand–receptor candidate set from the corresponding pooled regional object (Human CellChatDB v1, as recorded in the candidate metadata), the triMean expression summary, `population.size = FALSE`, and 100 internal permutations per donor. These internal CellChat permutations support its communication inference and are distinct from the donor-level pathology tests.

A donor-by-cell-type group was eligible when it contained more than 20 cells and the donor's CellChat network completed. A score was treated as zero only when the completed donor network and both cell groups were eligible and the interaction belonged to the fixed candidate universe; failed networks and ineligible groups were retained as missing and excluded from that feature's model. Pathology was coded nonAD=0, earlyAD=1, and lateAD=2. The primary model was donor-level ordinary least squares of communication probability on pathology stage, with HC3 heteroskedasticity-robust standard errors. Age at death, sex, and post-mortem interval were considered in that order and included only when complete, nonconstant, nonredundant with the existing design, and compatible with at least five residual degrees of freedom. Three categorical pathology contrasts were secondary and used the same donor covariates. Tests required at least three eligible donors per pathology group, full-rank design, and at least five residual degrees of freedom.

Benjamini–Hochberg adjustment was applied separately within each region × directed sender-to-receiver cell-type pair × prespecified test, across all eligible ligand–receptor candidates in that family. Only finite p-values from estimable, variable-score tests entered BH; constant scores and non-estimable tests retained NA p/q values with explicit statuses. These q-values control error within each family under the usual BH dependence conditions and do not claim global FDR control across families. For selected and FDR-supported features, 1,000 bootstrap replicates resampled donors within pathology groups while preserving each group's observed donor count; percentile intervals and trend direction frequencies describe donor-level uncertainty and stability rather than replacing regression inference.
"""
    (OUT / "methods" / "cellchat_methods.md").write_text(method)
    figure_legends = """# Figure legends

**Figure 1. Donor-level scores for every fixed selected CellChat feature.** Points show eligible donors; the line shows the fitted donor-level progression model at the mean covariate profile. CellChat groups required >20 cells and a completed donor network. The annotated slope is the OLS pathology-stage coefficient (nonAD=0, earlyAD=1, lateAD=2), with HC3 95% confidence interval; age at death, sex, and PMI are included only when the pair-specific design is estimable. Bootstrap uncertainty uses 1,000 donor resamples within pathology group. BH q-values are adjusted within region × directed sender/receiver pair × progression test. Every row of the fixed selected-feature list is shown regardless of new significance.

**Figure 2. Inferred ligand–receptor paths for all fixed selected features.** Arrows show the direction encoded by the CellChat sender → ligand–receptor → receiver result. Ligand and receptor labels come from candidate metadata. These are inferred communication hypotheses, not experimentally established or causal signaling.

**Figure 3. Overview of all directed cell-type pairs across regions.** Tile color is the median HC3-adjusted pathology-stage slope across estimable ligand–receptor interactions in each region × directed pair family. Cell labels give within-family BH q<0.05 interactions / estimable interactions; `n<3` marks insufficient donor coverage. The donor is the analysis unit, and groups require >20 cells plus a completed network.

**Figure 4. Donor and directed-pair coverage.** Top bars show all donors by region and pathology. Lower tiles show eligible donors for each directed pair after the >20-cell and completed-network rules. Failed or ineligible donor pairs are missing, not zero scores.
"""
    (OUT / "methods" / "figure_legends.md").write_text(figure_legends)

    # Self-contained report: figures embedded as base64, full TSVs included as links.
    n_sig = int((progression.q_value_within_family < .05).sum())
    n_tested = int(progression.p_value.notna().sum())
    n_contrast_sig = int((contrasts.q_value_within_family < .05).sum())
    n_selected_sig = int((selected.q_value_within_family < .05).sum())
    sig_table = progression.loc[progression.q_value_within_family < .05].sort_values("q_value_within_family").head(50)
    table_html = sig_table[["region", "source", "target", "interaction_id", "beta_pathology_stage", "ci_low", "ci_high", "p_value", "q_value_within_family", "n_nonAD", "n_earlyAD", "n_lateAD"]].to_html(index=False, float_format=lambda x: f"{x:.3g}", classes="data") if len(sig_table) else "<p>No within-family BH q&lt;0.05 progression hits.</p>"
    imgs = []
    for name, title in [("selected_features.png", "Fixed selected communication features"), ("selected_interaction_diagrams.png", "Directional interaction diagrams"), ("all_pair_overview.png", "All directed-pair overview"), ("coverage_qc.png", "Donor coverage and QC")]:
        raw = base64.b64encode((OUT / "figures" / name).read_bytes()).decode("ascii")
        imgs.append(f"<h2>{html.escape(title)}</h2><img alt='{html.escape(title)}' src='data:image/png;base64,{raw}'>")
    links = [
        ("Full progression tests", "../tables/pathology_progression_tests.tsv"),
        ("Full categorical contrasts", "../tables/pathology_pairwise_contrasts.tsv"),
        ("FDR family counts", "../tables/fdr_family_counts.tsv"),
        ("Bootstrap uncertainty", "../tables/stratified_donor_bootstrap.tsv"),
        ("Donor network audit", "../tables/network_status_audit.tsv"),
        ("Donor metadata used", "../tables/donor_metadata_used.tsv"),
        ("Score zero and missingness audit", "../tables/score_zero_missing_audit.tsv"),
        ("Cell-type coverage summary", "../tables/celltype_coverage_summary.tsv"),
        ("Directed-pair coverage", "../tables/directed_pair_coverage.tsv"),
        ("Pooled cell-type results (descriptive)", "../tables/pooled_celltype_differences_descriptive.tsv"),
        ("Pooled pathway results (descriptive)", "../tables/pooled_pathway_differences_descriptive.tsv"),
        ("Methods subsection", "../methods/cellchat_methods.md"),
        ("Figure legends", "../methods/figure_legends.md"),
        ("Analysis configuration", "../config/analysis_config.json"),
    ]
    report = f"""<!doctype html><html><head><meta charset='utf-8'><title>EXP002 CellChat finalization</title>
<style>body{{font:16px/1.5 system-ui,sans-serif;max-width:1200px;margin:2rem auto;padding:0 1rem;color:#222}}img{{width:100%;height:auto;border:1px solid #ddd}}table.data{{border-collapse:collapse;font-size:13px}}.data td,.data th{{border:1px solid #ccc;padding:4px 7px}}.note{{background:#f4f6f8;padding:1rem;border-left:4px solid #4477aa}}li{{margin:.2rem 0}}</style></head><body>
<h1>EXP002 CellChat pathology analysis: finalized donor-level results</h1>
<p>Generated {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}. Donor is the biological replicate.</p>
<div class='note'><b>Scope:</b> Reanalysis of existing donor-level CellChat scores; no donor CellChat networks were rerun. CellChat results are inferred ligand–receptor communication hypotheses, not demonstrated causal signaling.</div>
<h2>Summary</h2><ul><li>Progression tests with finite p-values: {n_tested:,}</li><li>Progression results with within-family BH q&lt;0.05: {n_sig:,}</li><li>Secondary categorical contrast results with within-family BH q&lt;0.05: {n_contrast_sig:,}</li><li>Fixed display features: {len(selected)}; within-family q&lt;0.05: {n_selected_sig}</li><li>FDR families: {len(families):,}; each is one region × directed sender/receiver pair × test.</li></ul>
<p>Progression uses donor-level OLS with pathology coded nonAD=0, earlyAD=1, lateAD=2 and HC3 standard errors. Age at death, sex, and PMI are used only when the pair-specific design remains estimable. A measured zero is assigned only for an eligible donor pair with a completed network; failures and ineligible groups remain missing.</p>
<h2>Sample coverage and score availability</h2><p>Regional donor counts, network completion, cell-type coverage, and aggregate counts of recorded nonzero scores, evaluated zero slots, and excluded missing donor-interaction slots are available in the linked audit tables. Eligibility requires >20 cells for both sender and receiver and a completed network.</p>
<h2>Top within-family progression results</h2>{table_html}
<h2>Figures</h2>{''.join(imgs)}
<h2>Methods</h2><pre>{html.escape(method)}</pre>
<h2>Downloadable full results and provenance</h2><ul>{''.join(f'<li><a href="{u}">{html.escape(n)}</a></li>' for n,u in links)}</ul>
<h2>Limitations</h2><ul><li>Within-family BH q-values do not control FDR across the full collection of regions, directions, and contrasts.</li><li>BH control uses its usual independence or positive-dependence conditions; correlated ligand–receptor tests may not satisfy these exactly.</li><li>LateAD donor groups are smaller; estimates and q-values are limited by donor coverage.</li><li>CellChat outputs communication hypotheses and do not establish causal signaling.</li></ul>
</body></html>"""
    (OUT / "report" / "cellchat_final_report.html").write_text(report)
    summary = [
        f"Progression features represented: {len(progression):,}",
        f"Finite progression tests: {n_tested:,}",
        f"Within-family progression q<0.05: {n_sig:,}",
        f"Categorical contrast rows: {len(contrasts):,}",
        f"Within-family contrast q<0.05: {n_contrast_sig:,}",
        f"Fixed selected features: {len(selected)}; q<0.05: {n_selected_sig}",
        f"Bootstrap targets: {len(bootstrap):,}; replicates per estimable feature: {BOOTSTRAPS}",
    ]
    (OUT / "report" / "analysis_summary.txt").write_text("\n".join(summary) + "\n")
    print("FINALIZATION_COMPLETE", *summary, sep="\n", flush=True)


if __name__ == "__main__":
    main()
