# External inputs and path configuration

The code does not bundle restricted or large input files. To rerun it, provide the corresponding datasets and references and update the cluster-specific paths retained in the scripts.

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

The principal paths used on the source system were `REGA_separate`, `REGA_shared/current`, and `Alzheimer_project_explore/experiments` under `/project/zduren/durenlab/naqing_palmetto/`. Many scripts have absolute paths inside their source code. Some accept environment variables such as `REGA_INPUTS_ROOT`, `REGA_MATRICES_ROOT`, `REGA_SHARED_ROOT`, or output-root arguments; consult each script before use. Update hard-coded paths in a separate working copy when running elsewhere, while retaining this code snapshot and `SOURCE_MANIFEST.tsv` as the provenance record.

The saved GMT inputs used in the separate-model GSEA were `GO_Biological_Process_2026_Human.gmt` (5,952 gene sets; SHA-256 `c9a726f9d4b0a8c8122ce154b633e6d35402785ade640b262be3ec18e7cb276c`) and `KEGG_2026.Human.gmt` (351 gene sets; SHA-256 `ff014e36320a8e4e09962f81ee7a720fb0001d3c349d021c118766b69fd78139`). Enrichr libraries can change at the same URL, so these checksums identify the analyzed files.

The code was run with the Palmetto `graph_clean` environment for Python and `bench_r` for R. External tools include bedtools, HOMER, MAGMA, LDSC, CrossMap, and CellChat. Exact software package versions and dataset accessions should be added to the manuscript and release metadata before a public reproducibility claim is made.
