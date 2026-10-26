# External inputs and path configuration

The code does not bundle restricted or large input files. To rerun it, supply the corresponding datasets and references through environment variables or explicit script arguments. Published code does not embed the original Palmetto project path.

| Resource | Used for |
| --- | --- |
| Processed snRNA-seq AnnData and individual/region/pathology metadata | RNA pseudobulk, model inputs, TF expression, CellChat |
| Cell type cCRE BED files, hg38 coordinates, protein-coding gene annotation | Regulatory priors and REGA input construction |
| HOMER motif collection and motif-to-TF correspondence | RE–motif graph and TF activity/assignment |
| Bellenguez AD GWAS summary statistics and reference genome/liftover files | MAGMA and LDSC |
| 1000 Genomes Phase 3 European reference files, HapMap3 SNPs, baseline-LD resources | LDSC annotations, LD scores, heritability |
| [Enrichr GO Biological Process 2026](https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2026) and [Enrichr KEGG 2026](https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2026) human gene sets | GSEA; retrieved by `separate/05_gsea/scripts/get_KEGG_2026_gmt.sh` |
| Healthy-aging reference gene lists | Aging enrichment |
| CellChat human ligand–receptor database | Cell–cell communication |

For CellChat, supply a study-harmonized raw-count H5AD through
`CELLCHAT_SOURCE_H5AD` and a separate output directory through
`CELLCHAT_WORK_DIR` and `CELLCHAT_REGION_INPUT_DIR`. The CellChat inference
retains all seven source cell classes, including Vasc; manuscript display can
focus on the six nonvascular classes. See
[`shared/part8_cellchat/README.md`](../shared/part8_cellchat/README.md).

## Configure local paths

Set the roots needed by the stage you run. For example:

```bash
export REGA_DATA_ROOT=/path/to/analysis-data
export REGA_SEPARATE_ROOT="$REGA_DATA_ROOT/REGA_separate"
export REGA_SHARED_ROOT="$REGA_DATA_ROOT/REGA_shared"
export AD_MANUSCRIPT_WORK_ROOT=/path/to/manuscript-work
export AD_DRIVER_TF_ROOT=/path/to/driver-tf-inputs
export SINGLE_GRN_ROOT=/path/to/independent-model-grns
```

`REGA_SEPARATE_ROOT` is required by the independent and shared model scripts and figure notebooks. `REGA_SHARED_ROOT` points to the current shared-model results. `AD_MANUSCRIPT_WORK_ROOT` contains the external `manuscript/` tables and TF-list input used by the manuscript builders; `REGA_DATA_ROOT` is the parent of the two REGA result roots when they share a parent. `AD_DRIVER_TF_ROOT` contains the donor pathology and driver-list inputs. Override the individual files with `AD_PATHOLOGY_FILE` and `AD_SNP_TF_LIST` if needed. Set `SINGLE_GRN_ROOT` explicitly to the independent-model GRNs used for the reported results; the shared GRN defaults to `REGA_SHARED_ROOT/part7_TF/05GRN/results` and can be overridden with `SHARED_GRN_ROOT`.

For the reported Oli_EC SPI1 target-gene expression panels, follow the two-step commands in [`manuscript_analysis/spi1_target_expression/README.md`](../manuscript_analysis/spi1_target_expression/README.md). Supply `SINGLE_GRN_ROOT` for the separate-model GRN and use separate output directories for the two fitted models. The broader shared-model target-expression scripts under `shared/part7_TF/06_target_expression` have their own preflight and Slurm array; `TF_TARGET_RESULTS_ROOT` selects that broader workflow's output directory. `SPP1_DRIVER_OUTPUT_ROOT` and `TF_VOLCANO_OUTPUT_ROOT` control the corresponding figure/table outputs.

For LDSC, set `LDSC_CHAIN_FILE`, `LDSC_PHASE3_DIR`, `LDSC_BIM_DIR`, `LDSC_CODE_DIR`, and `BEDTOOLS` for the stages that use them. For MAGMA, set `AD_GWAS_GRCH38_FILE`, `MAGMA_GENE_LOC`, `LDSC_CHAIN_FILE`, `LDSC_BIM_DIR`, and a fresh `MAGMA_OUTPUT_ROOT`; `MAGMA_BIN` can override the executable under `REGA_SEPARATE_ROOT/external_inputs/magma`. CrossMap must be available on `PATH`. These resources are external and are not included in the repository.

The saved GMT inputs used in the separate-model GSEA were `GO_Biological_Process_2026_Human.gmt` (5,952 gene sets; SHA-256 `c9a726f9d4b0a8c8122ce154b633e6d35402785ade640b262be3ec18e7cb276c`) and `KEGG_2026.Human.gmt` (351 gene sets; SHA-256 `ff014e36320a8e4e09962f81ee7a720fb0001d3c349d021c118766b69fd78139`). Enrichr libraries can change at the same URL, so these checksums identify the analyzed files.

The code was run with the Palmetto `graph_clean` environment for Python and `bench_r` for most R stages. CellChat used the existing `cellchat` environment. External tools include bedtools, HOMER, MAGMA, LDSC, and CrossMap. Exact software package versions and dataset accessions should be added to the manuscript and release metadata before a public reproducibility claim is made.
