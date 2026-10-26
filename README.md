# REGA Alzheimer analysis code

This folder contains the final analysis workflow for the Alzheimer manuscript. It combines **independent cell type–region REGA models** and **region-wise shared REGA models**, followed by pathology, module assignment, LDSC, GSEA, MAGMA, TF, aging, CellChat, and manuscript integration code. REGA models analyze six nonvascular cell classes across six regions, with at least 10 nuclei per RNA pseudobulk group. CellChat inference retains the original seven classes, including Vasc.

No raw RNA/ATAC data, model checkpoints, matrices, results, GWAS statistics, reference panels, logs, or figure outputs are included. Notebook outputs and attachments were removed. Most working files were copied and remain in their analysis folders. The CellChat stage is a curated donor-level linear/HC3 workflow and is documented separately under `shared/part8_cellchat/README.md`. CellChat inference retains seven cell classes, including Vasc; manuscript figures focus on the six nonvascular classes. The curated code has not yet been checked end to end against the saved result tables.

## Layout

| Directory | Contents |
| --- | --- |
| `separate/01_prepare` through `separate/07_tf` | Independent-model preparation, fitting, pathology, assignments, LDSC, GSEA, MAGMA, TF analyses, and figure code |
| `shared/part1_prepare_REGA` through `shared/part8_cellchat` | Shared-model workflow across all six regions and downstream CellChat analysis |
| `manuscript_analysis` | Paper-facing tables and figures, grouped by biological result; see its README for workflow order. |
| `SOURCE_MANIFEST.tsv` | One row per published code file, with the source location and SHA-256 hashes before and after notebook output removal |

## Manuscript analysis

[manuscript_analysis/README.md](manuscript_analysis/README.md) maps each folder to the reported figures and supplementary tables. The original numbered source paths remain in `SOURCE_MANIFEST.tsv` for provenance. Standalone SPP1 expression and SPP1–TF-activity correlation analyses were exploratory and are not included here.
CellChat code formerly filed under `EXP002` is organized in `shared/part8_cellchat`.

## How the current source was chosen

The independent-model code comes from `REGA_separate`. Its RNA preparation scripts express one workflow for all 36 nonvascular cell type–region pairs. The shared-model code combines the selected scripts from `REGA_shared` into one six-region workflow. Manuscript scripts come from `Alzheimer_project_explore/experiments`. Original source files and published hashes are recorded in the manifest where applicable; newly curated CellChat scripts have no unchanged source counterpart.

## Analysis order and external inputs

Read [WORKFLOW.md](docs/WORKFLOW.md) for the stage order and [EXTERNAL_INPUTS.md](docs/EXTERNAL_INPUTS.md) for data and reference resources that must be supplied separately. The original analyses ran on Clemson Palmetto using the existing `graph_clean` Python and `bench_r` R environments; CellChat used the dedicated `cellchat` environment. Other tools include HOMER, bedtools, MAGMA, LDSC, and CrossMap. Some `.sbatch` files are drafts; review requests and paths before submission.

The preparation scripts provide a direct path from source snRNA data to inputs for all 36 independent models; see [WORKFLOW.md](docs/WORKFLOW.md). The CellChat stage has its own [workflow](shared/part8_cellchat/README.md), with seven classes in inference and six in manuscript figures. Configure the external roots and references in [EXTERNAL_INPUTS.md](docs/EXTERNAL_INPUTS.md) before rerunning a stage. [PATH_AUDIT.md](docs/PATH_AUDIT.md) records the remaining configuration boundary. This is a code release, not a one-command executable workflow. No jobs are submitted by this repository.
