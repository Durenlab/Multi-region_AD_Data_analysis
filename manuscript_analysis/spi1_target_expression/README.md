# SPI1 target-gene expression

This folder reproduces the module-restricted Oli_EC SPI1 target-gene expression analysis reported in Figure 4i (shared model) and Supplementary Figure 6d (independent model).

1. `01_test_oli_ec_tf_targets.R` reads each model's TF-to-target GRN, gene and TF module assignments, cleaned donor expression, and pathology labels. For every Oli_EC TF with linear-model pathology FDR < 0.05, it keeps nonzero GRN edges only to genes assigned to that TF's module **in the same fitted model**. It then selects the top 10% by signed edge weight from those module-matched genes, standardizes each selected gene across donors, averages those z-scores per donor, and tests the positive ordered-pathology trend. It writes selected targets, donor scores, regression results, and an audit. The within-pair BH correction includes the tested Oli_EC TFs.
2. `02_plot_spi1_target_expression.R` reads the donor scores, checks the trend results for the same TF test family, and plots SPI1. It writes the panel used for Figure 4i or Supplementary Figure 6d.

Only module-restricted targets are analyzed.

For SPI1, the separate-model assignment is Oli_EC_M14; the shared-model assignment is the Oli_EC component of EC_M14. These are distinct fitted modules despite both carrying the label M14. In the audited saved target lists, all 21 separate-model SPI1 genes and all 30 shared-model SPI1 genes match their respective M14 gene assignments.

From the `Alzheimer_project_explore` directory, using the existing R environment:

```bash
ml anaconda3
source activate bench_r
export AD_MANUSCRIPT_WORK_ROOT="$PWD"
export REGA_DATA_ROOT="$(dirname "$PWD")"
export SINGLE_GRN_ROOT="$REGA_DATA_ROOT/AD_Driver_TF/cell_type/GRN"
analysis_dir="REGA_Alzheimer_analysis/manuscript_analysis/spi1_target_expression"
Rscript --vanilla "$analysis_dir/01_test_oli_ec_tf_targets.R" Oli_EC outputs/spi1_target_expression/separate celltype_specific
Rscript --vanilla "$analysis_dir/01_test_oli_ec_tf_targets.R" Oli_EC outputs/spi1_target_expression/shared cross_celltype
Rscript --vanilla "$analysis_dir/02_plot_spi1_target_expression.R" \
  outputs/spi1_target_expression/separate outputs/spi1_target_expression/shared
```

The scripts require the external data roots documented in [`EXTERNAL_INPUTS.md`](../../docs/EXTERNAL_INPUTS.md). Outputs stay outside the publication code folder. No job submission is required for this focused analysis.
