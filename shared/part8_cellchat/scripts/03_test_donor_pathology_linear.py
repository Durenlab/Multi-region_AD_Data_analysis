#!/usr/bin/env python3
"""Test seven-class donor-level CellChat scores against ordered pathology.

The only inferential test is score ~ pathology_stage with HC3 robust standard
errors. A BH family is one region and directed sender/receiver pair, across
its estimable candidate ligand-receptor interactions.
"""
from pathlib import Path
import os

import numpy as np
import pandas as pd
from scipy.stats import t as student_t


EXP = Path(os.environ["CELLCHAT_WORK_DIR"])
INPUT_TABLES = EXP / "tables"
INTERMEDIATE = EXP / "intermediate"
OUT_TABLES = EXP / "finalized" / "tables"
OUT = EXP / "finalized"
REGIONS = ["AG", "EC", "HC", "MTC", "PFC", "TH"]
CELLTYPES = ["Exc", "Inh", "Ast", "Oli", "OPC", "Mic_Immune", "Vasc"]
PATHOLOGY = ["nonAD", "earlyAD", "lateAD"]
STAGE = {"nonAD": 0.0, "earlyAD": 1.0, "lateAD": 2.0}
MIN_CELLS = 20
MIN_PER_STAGE = 3


def bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values for finite tests only."""
    p = np.asarray(pvalues, dtype=float)
    q = np.full(p.shape, np.nan)
    good = np.isfinite(p)
    if not good.any():
        return q
    pv = p[good]
    order = np.argsort(pv)
    ranked = pv[order] * len(pv) / np.arange(1, len(pv) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1].clip(0, 1)
    qgood = np.empty_like(ranked)
    qgood[order] = ranked
    q[good] = qgood
    return q


def fit_simple_ols_hc3(y: np.ndarray, x: np.ndarray) -> tuple[float, float, float, float, float]:
    """Fit y ~ 1 + pathology stage; return slope, HC3 SE, CI, and p-value."""
    n = len(y)
    df = n - 2
    xbar = float(x.mean())
    x_centered = x - xbar
    sxx = float(np.dot(x_centered, x_centered))
    slope = float(np.dot(x_centered, y - y.mean()) / sxx)
    residual = y - (float(y.mean()) + slope * x_centered)
    leverage = 1.0 / n + (x_centered * x_centered) / sxx
    slope_influence = x_centered / sxx
    variance = float(np.sum((slope_influence * residual / np.maximum(1.0 - leverage, np.finfo(float).eps)) ** 2))
    se = float(np.sqrt(max(variance, 0.0)))
    if not np.isfinite(se) or se <= 0:
        return slope, se, np.nan, np.nan, np.nan
    statistic = slope / se
    pvalue = float(2 * student_t.sf(abs(statistic), df))
    critical = float(student_t.ppf(0.975, df))
    return slope, se, slope - critical * se, slope + critical * se, pvalue


def main() -> None:
    if (OUT_TABLES / "pathology_progression_tests.tsv").exists():
        raise FileExistsError("Progression results already exist; use a new CELLCHAT_WORK_DIR")
    OUT_TABLES.mkdir(parents=True, exist_ok=True)
    all_results = []
    summary = []

    for region in REGIONS:
        print(f"Processing {region}", flush=True)
        candidates = pd.read_csv(INPUT_TABLES / f"{region}_candidate_interactions.tsv", sep="\t")
        interaction_ids = candidates.interaction_id.astype(str).drop_duplicates().tolist()
        coverage = pd.read_csv(INPUT_TABLES / f"{region}_donor_celltype_coverage.tsv", sep="\t")
        if not set(coverage.celltype.dropna()).issubset(CELLTYPES):
            raise ValueError(f"{region} coverage includes an unknown cell class")
        coverage = coverage[coverage.donor_network_status.eq("completed")].copy()
        coverage["donor_id"] = coverage.donor_id.astype(str)
        coverage["pathology"] = coverage.pathology.astype(str)
        if coverage.groupby("donor_id").pathology.nunique().max() > 1:
            raise ValueError(f"Conflicting pathology labels in {region}")
        donors = coverage[["donor_id", "pathology"]].drop_duplicates("donor_id").set_index("donor_id")
        cell_counts = coverage.pivot_table(index="donor_id", columns="celltype", values="n_cells", aggfunc="max")
        cell_counts = cell_counts.reindex(index=donors.index, columns=CELLTYPES, fill_value=0).fillna(0)
        scores = pd.read_csv(INTERMEDIATE / f"{region}_donor_nonzero_probabilities.tsv.gz", sep="\t")
        if not set(scores.source.dropna()).issubset(CELLTYPES) or not set(scores.target.dropna()).issubset(CELLTYPES):
            raise ValueError(f"{region} scores include an unknown cell class")
        scores["donor_id"] = scores.donor_id.astype(str)
        scores["interaction_id"] = scores.interaction_id.astype(str)
        scores["source"] = scores.source.astype(str)
        scores["target"] = scores.target.astype(str)
        scores = scores[scores.interaction_id.isin(interaction_ids)]
        results = []

        for source in CELLTYPES:
            for target in CELLTYPES:
                eligible = (cell_counts[source] > MIN_CELLS) & (cell_counts[target] > MIN_CELLS)
                meta = donors.loc[eligible].copy()
                stage_counts = meta.pathology.value_counts().reindex(PATHOLOGY, fill_value=0)
                subset = scores[(scores.source == source) & (scores.target == target)]
                if subset.duplicated(["donor_id", "interaction_id"]).any():
                    raise ValueError(f"Duplicate donor-feature scores for {region} {source}->{target}")
                matrix = subset.pivot(index="donor_id", columns="interaction_id", values="communication_probability")
                matrix = matrix.reindex(index=meta.index, columns=interaction_ids).fillna(0.0)
                x = meta.pathology.map(STAGE).to_numpy(dtype=float)

                for interaction_id in interaction_ids:
                    y = matrix[interaction_id].to_numpy(dtype=float)
                    row = {
                        "region": region, "source": source, "target": target,
                        "interaction_id": interaction_id,
                        "n_nonAD": int(stage_counts["nonAD"]),
                        "n_earlyAD": int(stage_counts["earlyAD"]),
                        "n_lateAD": int(stage_counts["lateAD"]),
                        "mean_nonAD": float(y[meta.pathology.to_numpy() == "nonAD"].mean()) if stage_counts["nonAD"] else np.nan,
                        "mean_earlyAD": float(y[meta.pathology.to_numpy() == "earlyAD"].mean()) if stage_counts["earlyAD"] else np.nan,
                        "mean_lateAD": float(y[meta.pathology.to_numpy() == "lateAD"].mean()) if stage_counts["lateAD"] else np.nan,
                        "beta_pathology_stage": np.nan, "se_HC3": np.nan,
                        "ci_low": np.nan, "ci_high": np.nan, "p_value": np.nan,
                        "status": "insufficient_group_coverage",
                        "design_rank": 2, "residual_df": max(len(meta) - 2, 0),
                        "covariates": "none",
                        "covariate_decisions": "No covariates included",
                        "model_formula": "score ~ pathology_stage",
                        "prediction_nonAD_mean_covars": np.nan,
                        "prediction_earlyAD_mean_covars": np.nan,
                        "prediction_lateAD_mean_covars": np.nan,
                    }
                    if (stage_counts >= MIN_PER_STAGE).all():
                        if np.unique(y).size < 2:
                            row["status"] = "no_score_variation"
                        else:
                            beta, se, lo, hi, pvalue = fit_simple_ols_hc3(y, x)
                            if np.isfinite(pvalue):
                                row.update({"beta_pathology_stage": beta, "se_HC3": se,
                                            "ci_low": lo, "ci_high": hi, "p_value": pvalue,
                                            "status": "tested"})
                                # Model predictions use E[y] at stage 0/1/2.
                                intercept = float(y.mean() - beta * x.mean())
                                row["prediction_nonAD_mean_covars"] = intercept
                                row["prediction_earlyAD_mean_covars"] = intercept + beta
                                row["prediction_lateAD_mean_covars"] = intercept + 2 * beta
                            else:
                                row["status"] = "not_estimable"
                    results.append(row)

        result = pd.DataFrame(results)
        testable = result.p_value.notna() & np.isfinite(result.p_value)
        result["q_value_within_family"] = np.nan
        for _, family_idx in result.loc[testable].groupby(["source", "target"]).groups.items():
            result.loc[family_idx, "q_value_within_family"] = bh_adjust(result.loc[family_idx, "p_value"].to_numpy())
        # Compatibility field: this mirrors the restored pair-specific q values.
        result["q_value_region_bh"] = result.q_value_within_family
        result["FDR_family"] = result.apply(lambda r: f"{region}: {r.source}->{r.target}", axis=1)
        all_results.append(result)
        summary.append({
            "region": region,
            "n_unique_ligand_receptor_pairs": len(interaction_ids),
            "n_ligand_receptor_tests": int(testable.sum()),
            "n_test_families": int(result.loc[testable, ["source", "target"]].drop_duplicates().shape[0]),
            "FDR_calculation": "Within each region x directed sender->receiver family, BH across finite p-values for unique ligand-receptor interaction IDs",
            "model": "score ~ pathology_stage (nonAD=0, earlyAD=1, lateAD=2); HC3 robust SE",
            "n_no_score_variation": int((result.status == "no_score_variation").sum()),
            "n_insufficient_group_coverage": int((result.status == "insufficient_group_coverage").sum()),
        })
        print(f"{region}: {int(testable.sum())} tests; {len(interaction_ids)} unique ligand-receptor interactions", flush=True)

    combined = pd.concat(all_results, ignore_index=True)
    combined.to_csv(OUT_TABLES / "pathology_progression_tests.tsv", sep="\t", index=False, na_rep="NA")
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(OUT_TABLES / "linear_test_fdr_summary.tsv", sep="\t", index=False)
    families = (combined.assign(significant=combined.q_value_within_family.lt(0.05))
                .groupby(["region", "source", "target"], as_index=False)
                .agg(n_candidates=("interaction_id", "size"),
                     n_tested=("p_value", "count"),
                     n_q_lt_0_05=("significant", "sum")))
    families.to_csv(OUT_TABLES / "fdr_family_counts.tsv", sep="\t", index=False)
    hits = combined.loc[
        combined.status.eq("tested") & combined.q_value_within_family.lt(0.05)
    ].copy()
    hits.to_csv(OUT_TABLES / "significant_progression.tsv", sep="\t", index=False)
    print(f"Saved {len(combined):,} candidate tests and {len(hits):,} FDR-significant rows")


if __name__ == "__main__":
    main()
