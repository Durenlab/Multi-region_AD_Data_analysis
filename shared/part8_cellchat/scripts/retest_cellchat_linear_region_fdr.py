#!/usr/bin/env python3
"""Refit CellChat donor-level progression models without covariates.

One BH family is defined per region and directed sender/receiver pair,
containing the estimable ligand-receptor tests for that pair.
"""
from pathlib import Path
import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.stats import t as student_t


EXP = Path(__file__).resolve().parents[1]
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
    OUT_TABLES.mkdir(parents=True, exist_ok=True)
    all_results = []
    summary = []

    for region in REGIONS:
        print(f"Processing {region}", flush=True)
        candidates = pd.read_csv(INPUT_TABLES / f"{region}_candidate_interactions.tsv", sep="\t")
        interaction_ids = candidates.interaction_id.astype(str).drop_duplicates().tolist()
        coverage = pd.read_csv(INPUT_TABLES / f"{region}_donor_celltype_coverage.tsv", sep="\t")
        coverage = coverage[coverage.donor_network_status.eq("completed")].copy()
        coverage["donor_id"] = coverage.donor_id.astype(str)
        coverage["pathology"] = coverage.pathology.astype(str)
        if coverage.groupby("donor_id").pathology.nunique().max() > 1:
            raise ValueError(f"Conflicting pathology labels in {region}")
        donors = coverage[["donor_id", "pathology"]].drop_duplicates("donor_id").set_index("donor_id")
        cell_counts = coverage.pivot_table(index="donor_id", columns="celltype", values="n_cells", aggfunc="max")
        cell_counts = cell_counts.reindex(index=donors.index, columns=CELLTYPES, fill_value=0).fillna(0)
        scores = pd.read_csv(INTERMEDIATE / f"{region}_donor_nonzero_probabilities.tsv.gz", sep="\t")
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
    # Refresh the legacy per-pair family audit for the simple-model results;
    # leave its categorical contrast rows untouched.
    family_path = OUT_TABLES / "fdr_family_counts.tsv"
    if family_path.exists():
        family_counts = pd.read_csv(family_path, sep="\t")
        trend_counts = (combined.assign(_hit=combined.q_value_within_family < .05)
                        .groupby(["region", "source", "target"], as_index=False)
                        .agg(n_eligible_candidates=("interaction_id", "size"),
                             n_finite_p=("p_value", "count"),
                             n_q_lt_0_05=("_hit", "sum")))
        trend_counts.insert(0, "test", "pathology_trend")
        keep = family_counts.test.ne("pathology_trend")
        family_counts = pd.concat([family_counts.loc[keep], trend_counts], ignore_index=True)
        family_counts.to_csv(family_path, sep="\t", index=False)

    # Keep the fixed-feature estimates and all-pair overview synchronized with
    # the updated model so downstream figures do not retain obsolete q-values.
    keys = ["region", "source", "target", "interaction_id"]
    selected_keys = pd.read_csv(OUT_TABLES / "fixed_selected_features_refreshed.tsv", sep="\t")[keys]
    selected = selected_keys.merge(combined, on=keys, how="left", validate="one_to_one")
    selected_cols = keys + ["beta_pathology_stage", "se_HC3", "ci_low", "ci_high", "p_value",
                            "q_value_within_family", "status", "covariates",
                            "prediction_nonAD_mean_covars", "prediction_earlyAD_mean_covars",
                            "prediction_lateAD_mean_covars"]
    selected[selected_cols].to_csv(OUT_TABLES / "fixed_selected_features_refreshed.tsv", sep="\t", index=False, na_rep="NA")
    pair_summary = (combined.groupby(["region", "source", "target"], as_index=False)
                    .agg(median_beta=("beta_pathology_stage", "median"),
                         n_estimable=("p_value", "count"),
                         n_q_lt_0_05=("q_value_within_family", lambda x: int((x < .05).sum()))))
    pair_summary["pair"] = pair_summary.source + " → " + pair_summary.target
    pair_summary.to_csv(OUT_TABLES / "all_pair_overview_summary.tsv", sep="\t", index=False, na_rep="NA")

    config_path = OUT / "config" / "analysis_config.json"
    config = json.loads(config_path.read_text())
    config["run_utc_linear_model_update"] = datetime.now(timezone.utc).isoformat()
    config["donor_metadata_covariates"] = []
    config["model"] = "donor-level OLS: score ~ pathology_stage (nonAD=0, earlyAD=1, lateAD=2); HC3 robust standard errors; no donor covariates"
    config["FDR"] = "Progression: BH separately within each region x directed sender-receiver pair, over finite p-values for unique ligand-receptor interaction IDs. Secondary pairwise contrasts retain their previously reported analysis."
    config["linear_test_summary"] = str(OUT_TABLES / "linear_test_fdr_summary.tsv")
    config_path.write_text(json.dumps(config, indent=2) + "\n")

    methods = OUT / "methods" / "cellchat_methods.md"
    methods.write_text("""## CellChat analysis\n\nCellChat communication was evaluated from region-specific normalized expression inputs derived from the ROSMAP single-cell RNA sequencing dataset. The donor was the biological replicate; communication features were stratified by brain region, sender cell type, receiver cell type, and ligand–receptor interaction. CellChat 2.2.0.9001 used the region-specific, pathology-blind ligand–receptor candidate set from the corresponding pooled regional object (Human CellChatDB v1, as recorded in the candidate metadata), the triMean expression summary, `population.size = FALSE`, and 100 internal permutations per donor. These internal permutations support CellChat communication inference and are distinct from donor-level pathology tests.\n\nA donor-by-cell-type group was eligible when it contained more than 20 cells and the donor's CellChat network completed. A score was treated as zero only when the completed donor network and both cell groups were eligible and the interaction belonged to the fixed candidate universe; failed networks and ineligible groups were missing and excluded from that feature's model. Pathology was coded nonAD=0, earlyAD=1, and lateAD=2. The updated progression model is donor-level ordinary least squares, `score ~ pathology_stage`, with HC3 heteroskedasticity-robust standard errors and no age, sex, or PMI covariates. Tests required at least three eligible donors in each pathology group.\n\nFor progression analysis, Benjamini–Hochberg adjustment is applied separately within each region × directed sender→receiver pair, across finite p-values for the unique ligand–receptor interaction IDs tested in that pair. Pairs in the same region use separate FDR families. Constant-score and non-estimable tests have missing p/q values and are excluded from BH. The result table reports nominal `p_value` and pair-specific `q_value_within_family`; the compatibility column `q_value_region_bh` mirrors these pair-specific values. The summary table reports distinct ligand–receptor interaction IDs, estimable interaction-by-cell-pair tests, and the number of pair-specific families in each region.\n\nPreviously reported categorical pathology contrasts were not refit in this update and retain their prior covariate-adjusted analyses. Existing bootstrap results also were not recalculated. CellChat outputs remain inferred communication hypotheses and do not establish causal signaling.\n""")
    figure_legend = OUT / "methods" / "figure_legends.md"
    figure_legend.write_text("""# Figure legends\n\n**Figure 1. Donor-level CellChat pathology trends across six regions.** Panel A displays all 49 directed sender→receiver cell-class pairs in each region. Tile color represents the median OLS pathology-stage slope across estimable ligand–receptor interactions; the number in each tile is the count with pair-specific BH q<0.05. Panel B shows fixed selected features with donor scores and simple linear pathology-stage trends. Pathology is coded nonAD=0, earlyAD=1, lateAD=2. Eligible cell groups required >20 cells and a completed donor network. Slopes use HC3 robust standard errors. The progression-test BH family contains all estimable ligand–receptor interaction × directed sender/receiver tests in the region.\n\n**Figure 2. Inferred ligand–receptor paths for fixed selected features.** Arrows show the direction encoded by the CellChat sender→ligand–receptor→receiver result. Ligand and receptor labels come from candidate metadata. These are inferred communication hypotheses, not experimentally established or causal signaling.\n\n**Figure 3. Donor and directed-pair coverage.** Donor totals are shown by region and pathology. Lower tiles show eligible donors for each directed pair after the >20-cell and completed-network rules. Failed or ineligible donor pairs are missing, not zero scores.\n""")

    q_hits = int((combined.q_value_within_family < .05).sum())
    selected_hits = int((selected.q_value_within_family < .05).sum())
    summary_lines = [
        f"Progression feature rows: {len(combined):,}",
        f"Finite progression tests: {int(combined.p_value.notna().sum()):,}",
        f"Pair-specific BH progression q<0.05: {q_hits:,}",
        f"Pair-specific FDR families: {int(summary_df.n_test_families.sum()):,} across six regions",
        f"Fixed selected features: {len(selected)}; pair-specific q<0.05: {selected_hits}",
        "Secondary categorical contrasts and bootstrap summaries were retained from the preceding covariate-adjusted analysis and were not recalculated.",
    ]
    (OUT / "report" / "analysis_summary.txt").write_text("\n".join(summary_lines) + "\n")
    update_note = f"""# CellChat simple linear progression reanalysis\n\n## Updated model and FDR family\n\nDonor scores were refit using `score ~ pathology_stage`, with nonAD=0, earlyAD=1, lateAD=2. Age at death, sex, and PMI were not included. HC3 robust standard errors were retained. Tests required at least three eligible donors per pathology group; eligible donor-celltype groups required >20 cells and a completed donor network.\n\nBH correction was run separately by region, across all finite p-values for interaction ID × directed sender→receiver tests. This is **one FDR family per region**. The test table reports both nominal `p_value` and regional adjusted `q_value_region_bh` (also mirrored in `q_value_within_family` for compatibility).\n\n## Test counts\n\n`linear_test_fdr_summary.tsv` provides the requested region-level counts: distinct ligand–receptor interaction IDs, estimable interaction-by-cell-pair tests, number of FDR families, and the adjustment rule. The count of tests is larger than the number of distinct ligand–receptor IDs because each interaction may be tested in multiple directed cell-type pairs.\n\n## Results\n\nThere are {int(combined.p_value.notna().sum()):,} estimable tests across six regions, with {q_hits} passing pair-specific BH q<0.05. The significant tests are:\n\n{combined.loc[combined.q_value_region_bh < .05, ["region", "source", "target", "interaction_id", "beta_pathology_stage", "p_value", "q_value_region_bh", "n_nonAD", "n_earlyAD", "n_lateAD"]].sort_values("q_value_region_bh").to_markdown(index=False)}\n\n## Scope\n\nThe updated progression table, selected-feature estimates, region-by-pair overview, methods, and configuration now use this model/FDR definition. Previously reported categorical contrasts and bootstrap summaries were not recalculated and remain from the preceding covariate-adjusted analysis.\n"""
    (OUT / "report" / "linear_model_reanalysis_report.md").write_text(update_note)
    report_path = OUT / "report" / "linear_model_reanalysis_report.md"
    report_text = report_path.read_text()
    report_text = report_text.replace(
        "BH correction was run separately by region, across all finite p-values for interaction ID × directed sender→receiver tests. This is **one FDR family per region**. The test table reports both nominal `p_value` and regional adjusted `q_value_region_bh` (also mirrored in `q_value_within_family` for compatibility).",
        "BH correction is run separately for each **region × directed sender→receiver pair**, across finite p-values for unique ligand–receptor interactions tested in that pair. Pairs in the same region have separate FDR families. The test table reports nominal `p_value` and pair-specific `q_value_within_family`; `q_value_region_bh` is a compatibility column mirroring those same values.")
    report_text = report_text.replace("q_value_region_bh", "q_value_within_family")
    report_text = report_text.replace(
        "The test table reports nominal `p_value` and pair-specific `q_value_within_family`; `q_value_within_family` is a compatibility column mirroring those same values.",
        "The test table reports nominal `p_value` and pair-specific `q_value_within_family`. The compatibility column `q_value_region_bh` mirrors the same pair-specific values.")
    report_text = report_text.replace("pair-specific BH q<0.05. The significant tests", "pair-specific BH q<0.05. These significant tests")
    report_path.write_text(report_text)
    figure_text = figure_legend.read_text().replace(
        "The progression-test BH family contains all estimable ligand–receptor interaction × directed sender/receiver tests in the region.",
        "BH correction is performed separately within each region × directed sender/receiver pair, across the estimable ligand–receptor interactions in that pair.")
    figure_legend.write_text(figure_text)

    html_path = OUT / "report" / "cellchat_final_report.html"
    if html_path.exists():
        html = html_path.read_text()
        notice = ("<div style='background:#fff4dc;padding:1rem;border-left:4px solid #c98200'>"
                  "<b>Analysis update:</b> The progression model and pair-specific FDR family were updated after this HTML report was generated. "
                  "Use <a href='linear_model_reanalysis_report.md'>the updated analysis note</a> and current linked TSV tables for progression results. "
                  "The report's embedded summary, hit table, and figures may reflect the previous model.</div>")
        if "Analysis update:</b> The progression model" not in html:
            html = html.replace("<body>", "<body>\n" + notice, 1)
            html_path.write_text(html)
    print("Wrote pathology_progression_tests.tsv and linear_test_fdr_summary.tsv", flush=True)


# pandas.to_markdown uses an optional dependency absent from the cluster env.
# Provide a plain-text fallback for the report table when tabulate is missing.
try:
    import tabulate as _tabulate  # noqa: F401
except ImportError:
    pd.DataFrame.to_markdown = lambda frame, index=False: frame.to_string(index=index)


if __name__ == "__main__":
    main()
