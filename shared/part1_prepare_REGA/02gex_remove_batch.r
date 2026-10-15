# ============================================================
# Usage:
#   Rscript 02gex_remove_batch.r 1
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
file_idx <- suppressWarnings(as.integer(args[1]))

if (is.na(file_idx) || file_idx < 1L) {
    stop("Please provide a positive file index, e.g. Rscript 02gex_remove_batch.r 1")
}

library(limma)

# ============================================================
# Directories
# ============================================================

stage_dir <- file.path(Sys.getenv("REGA_SHARED_ROOT", normalizePath("..")), "part1_prepare_REGA")
file_list <- sort(list.files(
    file.path(stage_dir, "00region_consistent_individuals"),
    pattern = "_region_consistent\\.tsv$", full.names = TRUE
))
file_list <- file_list[!grepl("^Vasc_", basename(file_list))]
if (length(file_list) != 36L) stop("Expected 36 RNA matrices across six cell types and regions")
if (file_idx > length(file_list)) {
    stop("File index ", file_idx, " exceeds ", length(file_list), " staged RNA matrices")
}
rna_path <- file_list[file_idx]
file_name <- basename(rna_path)

meta_dir <- file.path(stage_dir, "01region_consistent_metadata")
out_dir <- file.path(stage_dir, "02cleaned_region_consistent_gex")

if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

cat("Processing RNA file:", rna_path, "\n")

# ============================================================
# 1. Load RNA matrix
# genes x individuals
# ============================================================

rna_mat <- read.table(
    rna_path,
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE
)

# ============================================================
# 2. Find matching metadata
# Example:
#   Ast_AG_region_consistent.tsv
#   Ast_AG_region_consistent_metadata.tsv
# ============================================================

meta_name <- sub(
    "\\.tsv$",
    "_metadata.tsv",
    file_name
)

meta_path <- file.path(meta_dir, meta_name)

if (!file.exists(meta_path)) {
    stop(paste("Metadata file not found:", meta_path))
}

meta_df <- read.table(
    meta_path,
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE
)

# ============================================================
# 3. Align RNA and metadata
# ============================================================

common_ids <- intersect(colnames(rna_mat), rownames(meta_df))

if (length(common_ids) == 0) {
    stop("No overlapping individuals between RNA matrix and metadata.")
}

rna_mat <- as.matrix(rna_mat[, common_ids, drop = FALSE])
meta_df <- meta_df[common_ids, , drop = FALSE]

# Safety check
stopifnot(all(colnames(rna_mat) == rownames(meta_df)))

cat("Number of individuals:", ncol(rna_mat), "\n")
cat("Number of genes:", nrow(rna_mat), "\n")

# ============================================================
# 4. Clean metadata variables
# ============================================================

meta_df$Pathology <- as.factor(meta_df$Pathology)
meta_df$Sex <- as.factor(meta_df$Sex)

meta_df$age_death <- as.numeric(meta_df$age_death)
meta_df$pmi <- as.numeric(meta_df$pmi)

# Remove samples with missing covariates
keep <- complete.cases(
    meta_df[, c("Pathology", "Sex", "age_death", "pmi")]
)

rna_mat <- rna_mat[, keep, drop = FALSE]
meta_df <- meta_df[keep, , drop = FALSE]

cat("Individuals after removing missing covariates:", ncol(rna_mat), "\n")

# ============================================================
# 5. Design: protect Pathology effect
# ============================================================
# This means we remove Sex, age_death, and pmi effects,
# while preserving variation associated with Pathology.

design <- model.matrix(
    ~ Pathology,
    data = meta_df
)

# ============================================================
# 6. Batch/covariate correction
# ============================================================

cleaned <- removeBatchEffect(
    rna_mat,
    batch = meta_df$Sex,
    covariates = cbind(
        age_death = meta_df$age_death,
        pmi = meta_df$pmi
    ),
    design = design
)

# ============================================================
# 7. Save cleaned matrix
# ============================================================

out_path <- file.path(
    out_dir,
    paste0("cleaned_", file_name)
)

write.table(
    cleaned,
    out_path,
    sep = "\t",
    quote = FALSE,
    col.names = NA
)

cat("Finished:", file_name, "\n")
cat("Saved to:", out_path, "\n")
