# Paths requiring attention before a fresh rerun

The current model outputs live under `REGA_separate` and `REGA_shared/current`. The publication's five RNA preparation scripts accept explicit paths. Most other scripts retain their source paths and need site-specific path configuration elsewhere.

The following included code still contains an `AD_REGA` or `AD_REGA_region_wise` path. Some references designate unchanged raw/reference inputs; others read historical outputs. Inspect the use in each script before applying it to the current results. In particular, the `exp004` integration code predates the current six-module shared selection and should not be treated as a current-result rebuild without updating its inputs and candidate list.

- `manuscript_analysis/experiments/EXP005_results/scripts/spp1_top_driver_tf_overlap.R`
- `manuscript_analysis/experiments/EXP008_results/scripts/plot_spp1_celltype_expression.R`
- `manuscript_analysis/experiments/EXP011_results/pathology_boxplots/plot_spp1_pathology_boxplots.py`
- `manuscript_analysis/experiments/EXP011_results/scripts/analyze_spp1_tf_activity_correlations.py`
- `manuscript_analysis/experiments/exp004/characterize_regions.py`
- `manuscript_analysis/experiments/exp004/compare_ec_modules.py`
- `manuscript_analysis/experiments/exp004/finalize_audit.py`
- `manuscript_analysis/experiments/exp004/integrate_evidence.py`
- `manuscript_analysis/experiments/exp004/plot_ec_cellchat.py`
- `manuscript_analysis/experiments/exp004/write_report.py`
- `manuscript_analysis/experiments/exp008/scripts/plot_tf_tg_expression.R`
- `manuscript_analysis/experiments/exp008/scripts/prepare_cross_celltype_manifest.py`
- `shared/part8_cellchat/scripts/finalize_cellchat.py`
- `shared/part8_cellchat/scripts/submit_donor_cellchat.sbatch`

The model/analysis code and its provenance hashes are in `SOURCE_MANIFEST.tsv`. For the five consolidated RNA preparation scripts, the source and published hashes differ by design. Configure remaining site-specific paths before a fresh run.
