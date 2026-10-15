# REGA Alzheimer analysis code

This folder contains the final analysis workflow for the Alzheimer manuscript. It combines **independent cell type–region REGA models** and **region-wise shared REGA models**, followed by pathology, module assignment, LDSC, GSEA, MAGMA, TF, aging, CellChat, and manuscript integration code. The analyzed set contains six nonvascular cell classes across six regions. Each RNA pseudobulk group has at least 10 nuclei.

No raw RNA/ATAC data, model checkpoints, matrices, results, GWAS statistics, reference panels, logs, or figure outputs are included. Notebook outputs and attachments were removed. Original working files were copied and remain in their analysis folders.

## Layout

| Directory | Contents |
| --- | --- |
| `separate/01_prepare` through `separate/07_tf` | Independent-model preparation, fitting, pathology, assignments, LDSC, GSEA, MAGMA, TF analyses, and figure code |
| `shared/part1_prepare_REGA` through `shared/part8_cellchat` | Shared-model workflow across all six regions and downstream CellChat analysis |
| `manuscript_analysis` | Cross-branch integration, manuscript table builders, TF/GRN follow-up, and CellChat figure code |
| `SOURCE_MANIFEST.tsv` | One row per published code file, with the source location and SHA-256 hashes before and after notebook output removal |

## How the current source was chosen

The independent-model code comes from `REGA_separate`. Its RNA preparation scripts express one workflow for all 36 nonvascular cell type–region pairs. The shared-model code combines the selected scripts from `REGA_shared` into one six-region workflow. Manuscript scripts come from `Alzheimer_project_explore/experiments`. Original source files and published hashes are recorded in the manifest.

## Analysis order and external inputs

Read [WORKFLOW.md](docs/WORKFLOW.md) for the stage order and [EXTERNAL_INPUTS.md](docs/EXTERNAL_INPUTS.md) for data and reference resources that must be supplied separately. The original analyses ran on Clemson Palmetto using the existing `graph_clean` Python and `bench_r` R environments, with additional external programs for motif scanning, genomic intervals, MAGMA, LDSC, and CellChat. The `.sbatch` files record cluster requests as run; review and adapt them for another site.

The preparation scripts provide a direct path from source snRNA data to inputs for all 36 independent models; see [WORKFLOW.md](docs/WORKFLOW.md). Other scripts largely preserve the source code and may retain absolute paths to the original cluster directories, so a reader must configure input and output paths for their own data layout before rerunning. [PATH_AUDIT.md](docs/PATH_AUDIT.md) identifies code that still mentions historical `AD_REGA` locations. This is a code release, not a one-command executable workflow. No jobs are submitted by this repository.

## Citation
NA
