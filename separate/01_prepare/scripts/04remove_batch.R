#!/usr/bin/env Rscript
# Adjust expression for sex, age at death, and postmortem interval in every pair.

suppressPackageStartupMessages(library(limma))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript 04remove_batch.R PREP_DIR")

prep <- args[[1]]
input_dir <- file.path(prep, "original_individuals")
meta_dir <- file.path(prep, "metadata_by_group")
out_dir <- file.path(prep, "cleaned_expression")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rna_files <- sort(list.files(input_dir, pattern = "\\.tsv$", full.names = TRUE))
rna_files <- rna_files[basename(rna_files) != "summary_table.tsv"]
rna_files <- rna_files[!grepl("^Vasc_", basename(rna_files))]
if (length(rna_files) != 36L) stop("Expected 36 normalized nonvascular matrices")

for (rna_path in rna_files) {
  label <- sub("\\.tsv$", "", basename(rna_path))
  meta_path <- file.path(meta_dir, paste0(label, "_metadata.tsv"))
  if (!file.exists(meta_path)) stop("Missing metadata: ", meta_path)
  rna_mat <- read.table(rna_path, header = TRUE, sep = "\t", row.names = 1,
                        check.names = FALSE)
  meta_df <- read.table(meta_path, header = TRUE, sep = "\t", row.names = 1,
                        check.names = FALSE)
  if (!setequal(colnames(rna_mat), rownames(meta_df))) stop("Metadata mismatch: ", label)
  meta_df <- meta_df[colnames(rna_mat), , drop = FALSE]
  rna_mat <- as.matrix(rna_mat)

  design <- model.matrix(~Pathology, data = meta_df[, c(1, 2, 4, 5, 6, 7)])
  cleaned <- removeBatchEffect(
    rna_mat,
    batch = as.factor(meta_df$Sex),
    covariates = cbind(as.numeric(meta_df$age_death), as.numeric(meta_df$pmi)),
    design = design
  )
  write.table(cleaned, file.path(out_dir, paste0("cleaned_", label, ".tsv")),
              sep = "\t", quote = FALSE, col.names = NA)
  cat("Finished:", label, ncol(cleaned), "donors,", nrow(cleaned), "genes\n")
}
