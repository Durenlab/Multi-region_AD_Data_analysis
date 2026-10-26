# Workflow map

## Independent models (`separate/`)

1. `01_prepare/scripts`: create pseudobulk means from source snRNA data after study selection, excluding Vasc and groups with fewer than 10 nuclei; normalize expression for all 36 cell type–region pairs; align donor metadata; remove covariate effects; and build REGA inputs. Regulatory-prior helpers prepare cCRE, RE–TG, and motif resources. The full `Mic_Immune` label is preserved.
2. `02_run_rega/scripts`: fit 20-module models for each of the 36 nonvascular cell type–region combinations, test module activity against pathology, and select candidates. The current set has 23 pathology-associated modules.
3. `03_process_rega/scripts`: run `08export_stage2_matrices.py` once to export H, V, BV, WBV, and binary sparse W/B from each stage-two file; then assign genes, regulatory elements, and motifs to modules.
4. `04_ldsc/scripts`: prepare selected-module peak annotations, LD scores, and AD heritability enrichment.
5. `05_gsea/scripts` and `05_gsea/figures`: GO Biological Process and KEGG enrichment and summaries.
6. `06_magma/scripts`: GWAS gene-level analysis used for TF genetic evidence.
7. `07_tf/scripts`: TF activity, TF pathology associations, module assignment in candidate contexts, expression/activity comparisons, classification, and inferred GRNs. These scripts also produce the full context-wise classifier plots. The selected OPC_AG ROC panel, TF-count presentation, and overall-classifier panel are in `manuscript_analysis/tf_comparison`. Activity and expression pathology comparisons cover all 36 contexts; lack of a candidate module does not remove a context from those comparisons.

The preparation commands below use `PREP_DIR` as the output directory and a source snRNA AnnData supplied by the user:

```bash
python separate/01_prepare/scripts/01make_pseudobulk.py --source SOURCE_H5AD --output-dir PREP_DIR/pseudobulk_rna
python separate/01_prepare/scripts/02normalize_pseudobulk.py --pseudobulk PREP_DIR/pseudobulk_rna/rna_pseudobulk_celltype_individual_region_mean.h5ad --protein-coding-genes PREP_DIR/protein_coding_genes_list.tsv --output-dir PREP_DIR/original_individuals
python separate/01_prepare/scripts/03make_metadata.py --source SOURCE_H5AD --matrix-dir PREP_DIR/original_individuals --output-dir PREP_DIR/metadata_by_group
Rscript separate/01_prepare/scripts/04remove_batch.R PREP_DIR
bash separate/01_prepare/scripts/merge_bedfiles.sh CCRE_DIR PREP_DIR
bash separate/01_prepare/scripts/05TF_RE_prepare.sh TSS_WINDOWS_BED PREP_DIR
for CELLTYPE in Ast Exc Inh Mic_Immune Oli OPC; do bash separate/01_prepare/scripts/05scan_motifs.sh PREP_DIR MOTIF_FILE "$CELLTYPE"; done
bash separate/01_prepare/scripts/05process_motifs.sh PREP_DIR
Rscript separate/01_prepare/scripts/06prepare_REGA_inputs.r PREP_DIR
```

`CCRE_DIR` contains the subtype BED files, `TSS_WINDOWS_BED` contains the coding-gene TSS windows, and `MOTIF_FILE` is the HOMER motif collection. The `Mic` subtype cCRE files are merged into a `Mic_Immune` input with the full cell-type label. The commands produce `PREP_DIR/rega_inputs/REGION/CELLTYPE` for the 36 independent model fits.

## Region-wise shared models (`shared/`)

1. `part1_prepare_REGA`: align donors across the six nonvascular cell classes within each region, adjust expression, and build REGA inputs and blocks.
2. `part2_run_REGA`: fit a common regional H with cell type-specific regulatory loadings, test pathology, and select modules. The current selection is `AG_M13`, `EC_M14`, `EC_M17`, `MTC_M6`, `MTC_M13`, and `PFC_M19`.
3. `part3_process_REGA_output`: run `09extract_matrices.py` once per region to export H, V, BV, WBV, and binary sparse W/B from the same checkpoint load; then assign genes and REs and summarize the cell type composition of modules.
4. `part4_ldsc`: build and test 36 selected cell type–module annotations, with BH correction across those 36 tests.
5. `part5_gsea`: run GO/KEGG enrichment and aggregate positive hits for the current modules.
6. `part6_Age_related`: compare module gene rankings with healthy-aging signatures.
7. `part7_TF`: infer TF activity, assign TFs, test pathology, and reconstruct regulatory networks.
8. `part8_cellchat`: from a study-harmonized raw-count H5AD, prepare six regional inputs containing seven cell classes (including Vasc), define a pathology-blind ligand–receptor candidate set, infer donor-level scores, fit only the ordered-pathology linear model with HC3 standard errors, and generate tables and figures. The manuscript display focuses on six nonvascular classes while the analysis tables retain Vasc. See [`part8_cellchat/README.md`](../shared/part8_cellchat/README.md). No pooled pathology comparison, Spearman test, or pairwise pathology test belongs to the final inferential workflow.

The shared preparation scripts use the independent-model RNA matrices and align donors within each region before constructing regional model inputs.

## Manuscript integration (`manuscript_analysis/`)

[manuscript_analysis/README.md](../manuscript_analysis/README.md) lists the paper-facing analyses and their inputs in manuscript order. The two-step `spi1_target_expression` workflow tests module-restricted Oli_EC targets and then plots SPI1 for Figure 4i and Supplementary Figure 6d. `spp1_grn_overlap` supports Supplementary Table 7. `tf_comparison` contains the TF activity-versus-expression panels, including a focused OPC_AG ROC notebook; the context-wise fitting script already generates the full PR and ROC summaries. `module_tables` and `driver_tfs` assemble the reported supplementary tables and TF-module figures. Output tables and figures are intentionally absent.
