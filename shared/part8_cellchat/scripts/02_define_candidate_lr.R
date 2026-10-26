#!/usr/bin/env Rscript

# Define one pathology-blind ligand–receptor universe per region. This step
# screens signaling genes/interactions, but does not test pathology or estimate
# pooled communication probabilities. Only donor-level networks are inferred.

suppressPackageStartupMessages({
  library(CellChat)
  library(zellkonverter)
  library(SummarizedExperiment)
  library(Matrix)
  library(future)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: 02_define_candidate_lr.R REGION REGION_H5AD OUTPUT_ROOT")
}
region <- args[[1L]]
input_h5ad <- normalizePath(args[[2L]], mustWork = TRUE)
output_root <- args[[3L]]
regions <- c("AG", "EC", "HC", "MTC", "PFC", "TH")
celltypes <- c("Exc", "Inh", "Ast", "Oli", "OPC", "Mic_Immune", "Vasc")
if (!(region %in% regions)) stop("Unexpected region: ", region)

candidate_dir <- file.path(output_root, "candidate_universe")
dir.create(candidate_dir, recursive = TRUE, showWarnings = FALSE)
candidate_rds <- file.path(candidate_dir, paste0(region, "_candidate_lr.rds"))
candidate_tsv <- file.path(candidate_dir, paste0(region, "_candidate_lr.tsv"))
if (file.exists(candidate_rds) || file.exists(candidate_tsv)) {
  stop("Candidate output already exists for ", region)
}

sce <- zellkonverter::readH5AD(input_h5ad)
if (!("X" %in% SummarizedExperiment::assayNames(sce))) stop("Missing X assay")
expression <- SummarizedExperiment::assay(sce, "X")
if (!inherits(expression, "dgCMatrix")) expression <- as(expression, "dgCMatrix")
if (is.null(rownames(expression)) || is.null(colnames(expression))) {
  stop("Expression matrix needs gene and nucleus names")
}
meta <- as.data.frame(SummarizedExperiment::colData(sce))
if (!("celltype" %in% names(meta))) stop("Missing celltype metadata")
meta$celltype <- trimws(as.character(meta$celltype))
if (any(!meta$celltype %in% celltypes)) stop("Input contains excluded or unknown cell types")
if (anyDuplicated(rownames(expression))) stop("Duplicated gene names")
if (ncol(expression) != nrow(meta)) stop("Expression/metadata size mismatch")
rownames(meta) <- colnames(expression)
meta$celltype <- factor(meta$celltype, levels = celltypes)

cellchat <- CellChat::createCellChat(
  object = expression, meta = meta, group.by = "celltype",
  datatype = "RNA", do.sparse = TRUE
)
cellchat@DB <- CellChat::subsetDB(
  CellChatDB.human,
  search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact")
)
cellchat <- CellChat::subsetData(cellchat)
workers <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
if (is.na(workers) || workers < 1L) workers <- 1L
future::plan(future::multisession, workers = workers)
set.seed(20260728L)
cellchat <- CellChat::identifyOverExpressedGenes(
  cellchat, thresh.pc = 0.10, thresh.fc = 0.05, thresh.p = 0.05,
  do.fast = requireNamespace("presto", quietly = TRUE), min.cells = 10
)
cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
future::plan(future::sequential)
lr <- cellchat@LR$LRsig
if (is.null(lr) || nrow(lr) == 0L || is.null(rownames(lr)) ||
    anyDuplicated(rownames(lr))) {
  stop("No uniquely identified candidate ligand–receptor interactions")
}

# The donor script needs only this selected interaction table and the database;
# storing them avoids an unnecessary all-donor CellChat object in the release.
saveRDS(
  list(
    LRsig = lr, DB = cellchat@DB,
    region = region,
    celltypes = celltypes,
    source_h5ad = input_h5ad,
    selection = list(thresh.pc = 0.10, thresh.fc = 0.05,
                     thresh.p = 0.05, min.cells = 10,
                     categories = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"),
                     seed = 20260728L)
  ),
  candidate_rds
)
write.table(
  data.frame(interaction_id = rownames(lr), lr, check.names = FALSE),
  candidate_tsv, sep = "\t", quote = FALSE, row.names = FALSE
)
capture.output(sessionInfo(), file = file.path(candidate_dir, paste0(region, "_sessionInfo.txt")))
cat(region, ": ", nrow(lr), " pathology-blind candidate interactions\n", sep = "")
