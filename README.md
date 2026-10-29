# Multicellular regulatory circuits connect genetic risk and aging to regional vulnerability in Alzheimer’s disease

This repository contains the analysis code used for the study. It covers cell-type-specific and cross-cell-type REGA models, their downstream analyses, and code for the tables and figures reported in the manuscript. Large input datasets, fitted models, and result files are supplied separately.

The REGA analyses use six nonvascular cell classes (Exc, Inh, Ast, Oli, OPC, and Mic_Immune) across six brain regions (AG, EC, HC, MTC, PFC, and TH). RNA pseudobulk profiles are retained when a cell type–region–individual group contains at least 10 nuclei. The CellChat inference includes the source data's vascular class; the manuscript reports the six nonvascular classes.

## External data and references

| Resource | Source or location | Used for |
| --- | --- | --- |
| Processed snRNA-seq AnnData (`snRNA_Matrix.2263395_Cells_July7_2025.h5ad`) | [Liu et al. multiregion snRNA data](https://compbio.mit.edu/AD_Multiomic_MultiRegion/02_snRNA_h5ad/) | Individual-level RNA pseudobulk profiles and metadata |
| Cell subtype cCRE BED files | [Liu et al. cCRE catalog](https://compbio.mit.edu/AD_Multiomic_MultiRegion/03_cCREs/) | Regulatory-element catalogs for REGA priors |
| Alzheimer’s disease GWAS summary statistics (Bellenguez et al., 2022; `GCST90027158`) | [NHGRI–EBI GWAS Catalog study](https://www.ebi.ac.uk/gwas/studies/GCST90027158) | LDSC and MAGMA analyses |
| 1000 Genomes Phase 3 European PLINK reference, European allele frequencies, and HapMap3 regression weights excluding the MHC | [S-LDSC reference files](https://zenodo.org/records/7768714) (`1000G_Phase3_plinkfiles.tgz`, `1000G_Phase3_frq.tgz`, `1000G_Phase3_weights_hm3_no_MHC.tgz`) | Annotation construction, LD-score estimation, and heritability regression |
| HapMap3 SNP list | [LDSC reference-data instructions](https://github.com/bulik/ldsc/wiki/Cell%20type-specific-analyses); supplied to this workflow as `LDSC_PHASE3_DIR/list.txt` | SNP selection for LD-score estimation and GWAS munging |
| Baseline-LD v2.2 reference (97 annotations) | [S-LDSC reference files](https://zenodo.org/records/7768714) (`1000G_Phase3_baselineLD_v2.2_ldscores.tgz`) | Baseline annotations in partitioned heritability analysis |
| Human GO Biological Process 2026 and KEGG 2026 gene sets | [Ma’ayan Lab Enrichr](https://maayanlab.cloud/Enrichr/) | Module gene-set enrichment |
| Healthy-aging gene signatures | [Mesecar et al., 2026, Supplementary Data 6](https://doi.org/10.1038/s41514-026-00391-9) | Aging enrichment of module gene rankings |

The workflow also requires protein-coding gene annotation, an hg38-to-hg19 chain file, HOMER motifs and motif-to-TF mappings, and the CellChat human ligand–receptor database. See [external inputs and path configuration](docs/EXTERNAL_INPUTS.md) for the required files and settings. These data and references are not included in the code repository.

## Code layout

| Directory | Analysis |
| --- | --- |
| [`separate/`](separate/) | Cell-type-specific models: pseudobulk and regulatory-prior preparation, 36 independent REGA fits, pathology association, module assignment, LDSC, GSEA, MAGMA, TF activity, and TF–gene network analyses. |
| [`shared/`](shared/) | Cross-cell-type models: one jointly fitted model per region, shared-module pathology and assignment analyses, LDSC, GSEA, aging, TF and network analyses, and CellChat communication analysis. |
| [`manuscript_analysis/`](manuscript_analysis/) | Additional analyses used for reported figures and supplementary tables, including TF activity versus expression, candidate driver TFs, SPI1 target expression, and SPP1 network overlap. |

The stage-by-stage order for the cell-type-specific and cross-cell-type analyses is in [docs/WORKFLOW.md](docs/WORKFLOW.md). The [manuscript analysis guide](manuscript_analysis/README.md) maps downstream code to reported results. [`SOURCE_MANIFEST.tsv`](SOURCE_MANIFEST.tsv) records the origins of the curated code files.
