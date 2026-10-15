#!/usr/bin/env Rscript

# Compute donor-level CellChat probabilities for a fixed, pathology-blind
# candidate interaction set from the existing all-donor regional CellChat run.
# The donor is the unit of replication in downstream pathology comparisons.

suppressPackageStartupMessages({
  library(CellChat)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(future)
})

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 8 * 1024^3)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L || length(args) > 6L) {
  stop(paste(
    "Usage: run_donor_cellchat_by_region.R",
    "<REGION> <REGION_H5AD> <ALL_DONOR_CELLCHAT_RDS> <OUTPUT_ROOT>",
    "[MIN_CELLS=20] [NBOOT=100]"
  ))
}

region <- args[[1]]
input_h5ad <- normalizePath(args[[2]], mustWork = TRUE)
pooled_rds <- normalizePath(args[[3]], mustWork = TRUE)
output_root <- args[[4]]
min_cells <- if (length(args) >= 5L) as.integer(args[[5]]) else 20L
nboot <- if (length(args) >= 6L) as.integer(args[[6]]) else 100L

if (is.na(min_cells) || min_cells < 1L) stop("MIN_CELLS must be positive")
if (is.na(nboot) || nboot < 20L) stop("NBOOT must be at least 20")

region_order <- c("AG", "EC", "HC", "MTC", "PFC", "TH")
pathology_order <- c("nonAD", "earlyAD", "lateAD")
celltype_order <- c("Exc", "Inh", "Ast", "Oli", "OPC", "Mic_Immune", "Vasc")
if (!(region %in% region_order)) stop("Unexpected region: ", region)

intermediate_dir <- file.path(output_root, "intermediate")
table_dir <- file.path(output_root, "tables")
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

score_path <- file.path(
  intermediate_dir,
  paste0(region, "_donor_nonzero_probabilities.tsv.gz")
)
coverage_path <- file.path(
  table_dir,
  paste0(region, "_donor_celltype_coverage.tsv")
)
candidate_path <- file.path(
  table_dir,
  paste0(region, "_candidate_interactions.tsv")
)
if (file.exists(score_path) || file.exists(coverage_path) || file.exists(candidate_path)) {
  stop("Output exists; refusing to overwrite: ", region)
}

cat("Region:", region, "\n")
cat("Input:", input_h5ad, "\n")
cat("CellChat version:", as.character(packageVersion("CellChat")), "\n")
cat("Minimum retained cells per donor-celltype:", min_cells + 1L, "\n")
cat("Permutation count:", nboot, "\n")

cat("Reading regional AnnData...\n")
sce <- zellkonverter::readH5AD(input_h5ad)
assay_names <- SummarizedExperiment::assayNames(sce)
if (!("X" %in% assay_names)) stop("Assay X was not found")
data_input <- SummarizedExperiment::assay(sce, "X")
if (!inherits(data_input, "dgCMatrix")) data_input <- as(data_input, "dgCMatrix")
meta <- as.data.frame(SummarizedExperiment::colData(sce))
rownames(meta) <- colnames(data_input)

required <- c("Pathology", "celltype", "donor_id", "samples")
missing <- setdiff(required, colnames(meta))
if (length(missing)) stop("Missing metadata: ", paste(missing, collapse = ", "))
for (nm in required) meta[[nm]] <- trimws(as.character(meta[[nm]]))
if (anyNA(meta[, required]) || any(vapply(meta[, required], function(x) any(x == ""), logical(1)))) {
  stop("Required metadata contains missing or blank values")
}
if (!all(meta$samples == meta$donor_id)) stop("samples must match donor_id")
if (any(!meta$Pathology %in% pathology_order)) stop("Unexpected pathology label")
if (any(!meta$celltype %in% celltype_order)) stop("Unexpected celltype label")

donor_pathology <- unique(meta[c("donor_id", "Pathology")])
if (anyDuplicated(donor_pathology$donor_id)) {
  stop("A donor has more than one pathology label in this region")
}

cat("Loading pathology-blind interaction universe...\n")
pooled <- readRDS(pooled_rds)
if (!inherits(pooled, "CellChat")) stop("Pooled input is not a CellChat object")
lr_use <- pooled@LR$LRsig
if (is.null(lr_use) || nrow(lr_use) == 0L) stop("No regional candidate interactions found")
if (is.null(rownames(lr_use)) || anyDuplicated(rownames(lr_use))) {
  stop("Candidate interaction row names must be unique")
}
candidate <- data.frame(
  interaction_id = rownames(lr_use),
  lr_use,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
db_use <- pooled@DB
rm(pooled)
gc()

workers <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
if (is.na(workers) || workers < 1L) workers <- 1L
future::plan(future::multisession, workers = workers)
on.exit(future::plan(future::sequential), add = TRUE)

donor_ids <- unique(meta$donor_id)
coverage_rows <- list()
score_chunks <- list()
score_chunk_n <- 0L

for (donor_index in seq_along(donor_ids)) {
  donor <- donor_ids[[donor_index]]
  donor_idx <- which(meta$donor_id == donor)
  donor_meta <- meta[donor_idx, , drop = FALSE]
  pathology <- unique(donor_meta$Pathology)
  if (length(pathology) != 1L) stop("Donor has inconsistent pathology: ", donor)

  count_tab <- table(factor(donor_meta$celltype, levels = celltype_order))
  keep_types <- celltype_order[as.integer(count_tab) > min_cells]
  coverage_start <- length(coverage_rows) + 1L
  for (ct in celltype_order) {
    n_ct <- as.integer(count_tab[[ct]])
    coverage_rows[[length(coverage_rows) + 1L]] <- data.frame(
      region = region,
      pathology = pathology,
      donor_id = donor,
      celltype = ct,
      n_cells = n_ct,
      min_cells_required = min_cells + 1L,
      retained_for_network = n_ct > min_cells,
      donor_network_status = "pending",
      failure_message = "",
      status = if (length(keep_types) < 2L) "skipped_fewer_than_two_celltypes" else
        if (n_ct > min_cells) "retained" else "below_minimum",
      stringsAsFactors = FALSE
    )
  }

  cat(sprintf("[%d/%d] %s (%s): %d cells, %d cell types retained\n",
              donor_index, length(donor_ids), donor, pathology,
              length(donor_idx), length(keep_types)))
  if (length(keep_types) < 2L) {
    for (ii in coverage_start:length(coverage_rows)) {
      coverage_rows[[ii]]$donor_network_status <- "skipped_fewer_than_two_celltypes"
    }
    next
  }

  use_local <- donor_meta$celltype %in% keep_types
  donor_cols <- donor_idx[use_local]
  donor_data <- data_input[, donor_cols, drop = FALSE]
  donor_meta <- donor_meta[use_local, , drop = FALSE]
  donor_meta$celltype <- factor(donor_meta$celltype, levels = keep_types)
  rownames(donor_meta) <- colnames(donor_data)

  seed <- 20260923L + donor_index
  fit <- tryCatch({
    cellchat <- createCellChat(
      object = donor_data,
      meta = donor_meta,
      group.by = "celltype",
      datatype = "RNA",
      do.sparse = TRUE
    )
    cellchat@DB <- db_use
    cellchat <- subsetData(cellchat)

    # Use one common interaction set for every donor in this region. This set
    # comes from the all-donor regional analysis, not pathology-specific objects.
    set.seed(seed)
    cellchat <- computeCommunProb(
      cellchat,
      type = "triMean",
      LR.use = lr_use,
      raw.use = TRUE,
      population.size = FALSE,
      nboot = nboot,
      seed.use = seed
    )
    cellchat <- filterCommunication(
      cellchat,
      min.cells = min_cells,
      min.samples = 1L
    )
    list(prob = cellchat@net$prob, pval = cellchat@net$pval)
  }, error = function(e) e)

  if (inherits(fit, "error")) {
    message("Donor analysis failed for ", donor, ": ", conditionMessage(fit))
    for (ii in coverage_start:length(coverage_rows)) {
      coverage_rows[[ii]]$donor_network_status <- "failed"
      coverage_rows[[ii]]$failure_message <- conditionMessage(fit)
    }
    rm(donor_data, donor_meta, fit)
    gc(verbose = FALSE)
    next
  }
  for (ii in coverage_start:length(coverage_rows)) {
    coverage_rows[[ii]]$donor_network_status <- "completed"
  }
  prob <- fit$prob
  pval <- fit$pval
  idx <- which(prob > 0, arr.ind = TRUE)
  if (nrow(idx) > 0L) {
    score_chunk_n <- score_chunk_n + 1L
    score_chunks[[score_chunk_n]] <- data.frame(
      region = region,
      pathology = pathology,
      donor_id = donor,
      source = dimnames(prob)[[1]][idx[, 1]],
      target = dimnames(prob)[[2]][idx[, 2]],
      interaction_id = dimnames(prob)[[3]][idx[, 3]],
      communication_probability = as.numeric(prob[idx]),
      within_donor_p_value = as.numeric(pval[idx]),
      random_seed = seed,
      stringsAsFactors = FALSE
    )
  }
  rm(fit, donor_data, donor_meta, prob, pval, idx)
  gc(verbose = FALSE)
}

coverage <- do.call(rbind, coverage_rows)
write.table(coverage, coverage_path, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(candidate, candidate_path, sep = "\t", row.names = FALSE, quote = FALSE)

if (length(score_chunks)) {
  scores <- do.call(rbind, score_chunks)
} else {
  scores <- data.frame(
    region = character(), pathology = character(), donor_id = character(),
    source = character(), target = character(), interaction_id = character(),
    communication_probability = numeric(), within_donor_p_value = numeric(),
    random_seed = integer()
  )
}
con <- gzfile(score_path, open = "wt")
write.table(scores, con, sep = "\t", row.names = FALSE, quote = FALSE)
close(con)

sink(file.path(intermediate_dir, paste0(region, "_sessionInfo.txt")))
print(sessionInfo())
sink()

cat("Completed:", region, "donors:", length(donor_ids),
    "nonzero donor scores:", nrow(scores), "\n")
