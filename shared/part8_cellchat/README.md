# CellChat analysis: donor-level linear pathology test

This is the intended publication workflow. The input is the study-harmonized
snRNA AnnData with raw UMI counts in `X`, human gene symbols in `var_names`, and
`BrainRegion`, `RNA.Class`, `ROSMAP_IndividualID`, `Sample`, and `Pathology` in
`obs`. The source file and large derived results are supplied separately.
Vasc remains in the CellChat inference and saved results.

The seven inferred classes are `Exc`, `Inh`, `Ast`, `Oli`, `OPC`,
`Mic_Immune`, and `Vasc`. Manuscript figures and narrative focus on the six
nonvascular classes; the all-results and supplementary tables retain Vasc rows.
The six regions are `AG`, `EC`, `HC`, `MTC`, `PFC`, and `TH`.

## Execution order

Run data-heavy stages in a compute allocation. The existing `graph_clean`
environment supplies the Python preprocessing packages, and the existing
`cellchat` environment supplies CellChat and its R dependencies. No new
environment is required. The optional Slurm donor wrapper is a draft and must
be reviewed before submission.

1. **Regional inputs:** `scripts/01_prepare_cellchat_region_inputs.py` reads
   `CELLCHAT_SOURCE_H5AD` and writes six normalized H5AD files in
   `CELLCHAT_REGION_INPUT_DIR`. It samples at most 100
   nuclei per donor × cell class and 5,000 per regional class, removes genes
   with zero counts in that region, then normalizes each nucleus to 10,000
   counts and applies `log1p`. It also writes donor/cell coverage and a
   provenance manifest. The source H5AD is unchanged.
2. **Candidate LR universe:** `scripts/02_define_candidate_lr.R` screens each
   seven-class regional input without pathology labels. It uses human CellChatDB
   interactions in Secreted Signaling, ECM–Receptor, and Cell–Cell Contact and
   the original overexpression thresholds (`thresh.pc=0.10`,
   `thresh.fc=0.05`, `thresh.p=0.05`, `min.cells=10`). It saves only the selected
   interaction set and database under `candidate_universe/`. This is an input
   filter, not a pooled communication or pathology test.
3. **Donor communication scores:** `scripts/run_donor_cellchat_by_region.R`
   fits one CellChat network per donor using the same regional candidate LR
   set. A donor's cell class is eligible only with **more than 20 nuclei**;
   donors with fewer than two eligible classes are skipped. CellChat uses
   `triMean`, `population.size=FALSE`, and 100 internal permutations. Its
   internal p values are not pathology-test p values. The script writes sparse
   nonzero communication probabilities, candidate annotations, and coverage
   and completion status. Missing/failed donor-cell pairs stay missing;
   absent sparse rows count as zero only for completed eligible networks.
4. **Only pathology test:** `scripts/03_test_donor_pathology_linear.py` fits
   `score ~ pathology_stage` for each region × directed sender→receiver × LR
   feature, with nonAD=0, earlyAD=1, lateAD=2 and HC3 robust standard errors.
   It requires at least three eligible donors in each pathology group. Only
   estimable finite p values enter BH correction, separately within each
   region × directed sender→receiver pair. No Spearman, pairwise pathology
   tests, or donor covariates are used.
5. **Tables and figures:** `scripts/04_build_cellchat_supplement.py` joins the
   tested features to LR metadata and writes all linear tests plus the
   significant Supplementary Table 6 rows. Use
   `scripts/make_figure2_style_cellchat_by_region.py` and
   `finalized/figures/region_significant_LR/panel_c_dotplot.py` to generate
   regional figures from `pathology_progression_tests.tsv`; figure code omits
   Vasc for display while retaining its rows in the analysis table.

Set `CELLCHAT_WORK_DIR` to a separate results directory outside the code
release; no large output belongs in the GitHub folder. The preparation script
also takes `CELLCHAT_SOURCE_H5AD` and `CELLCHAT_REGION_INPUT_DIR`. The R
scripts take paths as arguments. For each region, call the candidate script
before the donor script:

```bash
cd shared/part8_cellchat
export CELLCHAT_SOURCE_H5AD=/path/to/study_harmonized_raw_counts.h5ad
export CELLCHAT_WORK_DIR=/path/to/cellchat_results
export CELLCHAT_REGION_INPUT_DIR="$CELLCHAT_WORK_DIR/01_region_inputs"
python scripts/01_prepare_cellchat_region_inputs.py
for REGION in AG EC HC MTC PFC TH; do
  Rscript scripts/02_define_candidate_lr.R "$REGION" "$CELLCHAT_WORK_DIR/01_region_inputs/${REGION}_cellchat_input.h5ad" "$CELLCHAT_WORK_DIR"
  Rscript scripts/run_donor_cellchat_by_region.R "$REGION" "$CELLCHAT_WORK_DIR/01_region_inputs/${REGION}_cellchat_input.h5ad" "$CELLCHAT_WORK_DIR/candidate_universe/${REGION}_candidate_lr.rds" "$CELLCHAT_WORK_DIR" 20 100
done
python scripts/03_test_donor_pathology_linear.py
python scripts/04_build_cellchat_supplement.py
python scripts/make_figure2_style_cellchat_by_region.py
python finalized/figures/region_significant_LR/panel_c_dotplot.py
```

The R steps may be scheduled by region, but the linear test must wait until all
six donor outputs are complete. The commands above describe order; they are
not a request to run the full dataset on a login node.
The source-study filtering that created the historical `AD_scRNA_used.h5ad`
is documented in the original `check_data.ipynb`; supply an equivalently
harmonized source H5AD when starting this workflow.

## Version boundary

The saved CellChat results and manuscript Supplementary Table 6 were produced
from a seven-class input. Omitting Vasc from the written interpretation does
not require refitting the donor-level linear tests. This curated code has not
yet been run end to end against those saved inputs, so its outputs should be
compared with the existing tables before replacing them. Earlier pooled
pathology comparisons, Spearman analyses, categorical contrasts, and bootstrap
summaries are outside the final inferential workflow.
