#!/usr/bin/env Rscript
shared_root <- Sys.getenv("REGA_SHARED_ROOT")
if (!nzchar(shared_root)) {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(file_arg)) stop("Set REGA_SHARED_ROOT when sourcing this script")
  shared_root <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
  while (!(basename(shared_root) %in% c("REGA_shared", "shared")) && dirname(shared_root) != shared_root) shared_root <- dirname(shared_root)
  if (!(basename(shared_root) %in% c("REGA_shared", "shared"))) stop("Could not find shared analysis root")
}
separate_root <- Sys.getenv("REGA_SEPARATE_ROOT")
if (!nzchar(separate_root)) stop("Set REGA_SEPARATE_ROOT")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

# ============================================================
# Region-wise TF assignment using observed / expected TF-module matrix
#
# Logic:
#   For each brain region:
#     1. Read all celltype-region TF-by-module matrices in that region
#     2. Rename TFs as {celltype_region}_{TF}
#     3. Concatenate rows across cell types
#     4. Compute expected matrix:
#          expected[i,j] = rowSum[i] * colSum[j] / totalSum
#     5. Compute fold change:
#          FC[i,j] = observed[i,j] / expected[i,j]
#     6. Assign each TF to the module with max FC
#     7. Keep only TF with max FC >= FC_THRESHOLD
#
# Input:
#   /tf_by_module/*_TF_by_module.tsv
#
# Output:
#   /tf_assignment/<region>/
# ============================================================

# ============================================================
# User config
# ============================================================

BASE_DIR <- shared_root

INPUT_ROOT <- file.path(
  BASE_DIR,
  "part7_TF/01tf_by_module"
)

OUT_ROOT <- file.path(
  BASE_DIR,
  "part7_TF/02tf_assignment_FC"
)

FC_THRESHOLD <- 2

dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Helper functions
# ============================================================

parse_file_info <- function(fn_base) {
  # Example:
  #   Ast_PFC_TF_by_module.tsv
  #   Mic_Immune_MTC_TF_by_module.tsv
  #
  # Returns:
  #   celltype_region = Ast_PFC
  #   region = PFC
  #
  # Assumes region is the last underscore-separated field
  # before the "_TF_by_module.tsv" suffix.
  x <- sub("_TF_by_module\\.tsv$", "", fn_base)

  parts <- strsplit(x, "_", fixed = TRUE)[[1]]
  if (length(parts) < 2) {
    stop("Cannot parse celltype_region/region from filename: ", fn_base)
  }

  region <- parts[length(parts)]
  celltype_region <- paste(parts, collapse = "_")
  celltype <- paste(parts[-length(parts)], collapse = "_")

  list(
    celltype_region = celltype_region,
    celltype = celltype,
    region = region
  )
}

make_expected_matrix <- function(M) {
  row_sum <- rowSums(M, na.rm = TRUE)
  col_sum <- colSums(M, na.rm = TRUE)
  total_sum <- sum(M, na.rm = TRUE)

  if (!is.finite(total_sum) || total_sum <= 0) {
    stop("Total matrix sum is zero or non-finite; cannot compute expected matrix.")
  }

  expected <- outer(row_sum, col_sum, FUN = "*") / total_sum
  dimnames(expected) <- dimnames(M)

  expected
}

assign_by_fc <- function(M, fc_threshold = 2) {
  M <- as.matrix(M)
  M[!is.finite(M)] <- 0
  M[M < 0] <- 0

  expected <- make_expected_matrix(M)

  fc <- M / expected
  fc[!is.finite(fc)] <- NA_real_

  max_idx <- max.col(fc, ties.method = "first")
  max_fc <- fc[cbind(seq_len(nrow(fc)), max_idx)]
  best_module <- colnames(fc)[max_idx]

  assigned <- ifelse(
    is.finite(max_fc) & max_fc >= fc_threshold,
    best_module,
    NA_character_
  )

  result <- data.frame(
    cell_gene = rownames(M),
    assigned_module = assigned,
    max_fc = max_fc,
    observed_at_assigned = M[cbind(seq_len(nrow(M)), max_idx)],
    expected_at_assigned = expected[cbind(seq_len(nrow(M)), max_idx)],
    stringsAsFactors = FALSE
  )

  list(
    assignment = result,
    fc = fc,
    expected = expected
  )
}

# ============================================================
# Load all TF-by-module files once
# ============================================================

files <- list.files(
  INPUT_ROOT,
  pattern = "_TF_by_module\\.tsv$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No TF-by-module TSV files found in: ", INPUT_ROOT)
}

cat("Found TF-by-module files:", length(files), "\n")

file_info <- lapply(basename(files), parse_file_info)
file_meta <- bind_rows(lapply(file_info, as.data.frame))

file_meta$path <- files

# Group files by region
region_groups <- split(file_meta, file_meta$region)
selected_region <- Sys.getenv("REGA_REGION", unset = "")
if (nzchar(selected_region)) {
  if (!selected_region %in% names(region_groups)) stop("Unknown or unavailable REGA_REGION: ", selected_region)
  region_groups <- region_groups[selected_region]
}

# ============================================================
# Main processing
# ============================================================

region_summaries <- list()

for (region_name in names(region_groups)) {

  cat("\n============================================================\n")
  cat("Region:", region_name, "\n")
  cat("============================================================\n")

  out_dir <- file.path(OUT_ROOT, region_name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  region_files <- region_groups[[region_name]]

  matrix_list <- list()
  row_meta_list <- list()
  module_order <- NULL

  # ------------------------------------------------------------
  # Read and concatenate matrices for this region
  # ------------------------------------------------------------

  for (i in seq_len(nrow(region_files))) {

    f <- region_files$path[i]
    celltype_region <- region_files$celltype_region[i]

    cat("Reading:", basename(f), "->", celltype_region, "\n")

    df <- read_tsv(f, show_col_types = FALSE)

    if (ncol(df) < 2) {
      warning("Skipping file with <2 columns: ", f)
      next
    }

    tf_names <- as.character(df[[1]])

    M <- as.matrix(df[, -1, drop = FALSE])
    storage.mode(M) <- "numeric"

    M[!is.finite(M)] <- 0
    M[M < 0] <- 0

    rownames(M) <- paste0(celltype_region, "_", tf_names)

    if (is.null(module_order)) {
      module_order <- colnames(M)
    }

    # Add missing modules as zeros and reorder
    missing_modules <- setdiff(module_order, colnames(M))
    if (length(missing_modules) > 0) {
      M <- cbind(
        M,
        matrix(
          0,
          nrow = nrow(M),
          ncol = length(missing_modules),
          dimnames = list(rownames(M), missing_modules)
        )
      )
    }

    M <- M[, module_order, drop = FALSE]

    matrix_list[[celltype_region]] <- M

    row_meta_list[[celltype_region]] <- data.frame(
      cell_gene = rownames(M),
      region = region_name,
      celltype_region = celltype_region,
      tf = tf_names,
      stringsAsFactors = FALSE
    )
  }

  if (length(matrix_list) == 0) {
    warning("No valid matrices loaded for region: ", region_name)
    next
  }

  observed <- do.call(rbind, matrix_list)
  row_meta <- bind_rows(row_meta_list)

  cat("Concatenated observed matrix:\n")
  cat("  rows:", nrow(observed), "\n")
  cat("  cols:", ncol(observed), "\n")

  # ------------------------------------------------------------
  # Save observed matrix
  # ------------------------------------------------------------

  observed_out <- file.path(
    out_dir,
    paste0(region_name, "_observed_cell_TF_by_module.tsv")
  )

  write_tsv(
    as.data.frame(observed) |>
      tibble::rownames_to_column("cell_gene"),
    observed_out
  )

  # ------------------------------------------------------------
  # Compute expected and FC, then assign
  # ------------------------------------------------------------

  res <- assign_by_fc(
    observed,
    fc_threshold = FC_THRESHOLD
  )

  assignment_raw <- res$assignment
  fc <- res$fc
  expected <- res$expected

  assignment_df <- assignment_raw %>%
    left_join(row_meta, by = "cell_gene") %>%
    select(
      region,
      celltype_region,
      tf,
      cell_gene,
      assigned_module,
      max_fc,
      observed_at_assigned,
      expected_at_assigned
    ) %>%
    mutate(
      assigned = !is.na(assigned_module)
    )

  assigned_long <- assignment_df %>%
    filter(assigned) %>%
    transmute(
      region = region,
      celltype_region = celltype_region,
      module = assigned_module,
      tf = tf,
      cell_gene = cell_gene,
      max_fc = max_fc,
      observed = observed_at_assigned,
      expected = expected_at_assigned
    ) %>%
    arrange(region, celltype_region, module, desc(max_fc))

  # ------------------------------------------------------------
  # Save assignment outputs
  # ------------------------------------------------------------

  all_assignment_out <- file.path(
    out_dir,
    paste0(region_name, "_ALL_cell_TF_FC_assignment_with_unassigned.tsv")
  )

  assigned_out <- file.path(
    out_dir,
    paste0(region_name, "_ALL_celltype_region_TF_module_assignments_FC.tsv")
  )

  write_tsv(assignment_df, all_assignment_out)
  write_tsv(assigned_long, assigned_out)

  # ------------------------------------------------------------
  # Save expected and FC matrices
  # ------------------------------------------------------------

  expected_out <- file.path(
    out_dir,
    paste0(region_name, "_expected_cell_TF_by_module.tsv")
  )

  fc_out <- file.path(
    out_dir,
    paste0(region_name, "_FC_observed_over_expected_cell_TF_by_module.tsv")
  )

  write_tsv(
    as.data.frame(expected) |>
      tibble::rownames_to_column("cell_gene"),
    expected_out
  )

  write_tsv(
    as.data.frame(fc) |>
      tibble::rownames_to_column("cell_gene"),
    fc_out
  )

  # ------------------------------------------------------------
  # Summary tables
  # ------------------------------------------------------------

  module_counts <- assigned_long %>%
    count(region, celltype_region, module, name = "n_tfs") %>%
    arrange(region, celltype_region, module)

  module_counts_out <- file.path(
    out_dir,
    paste0(region_name, "_FC_assignment_module_counts.tsv")
  )

  write_tsv(module_counts, module_counts_out)

  celltype_summary <- assignment_df %>%
    group_by(region, celltype_region) %>%
    summarise(
      n_tfs_total = n(),
      n_tfs_assigned = sum(assigned),
      assigned_fraction = n_tfs_assigned / n_tfs_total,
      median_max_fc = median(max_fc, na.rm = TRUE),
      mean_max_fc = mean(max_fc, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_tfs_assigned))

  celltype_summary_out <- file.path(
    out_dir,
    paste0(region_name, "_FC_assignment_celltype_summary.tsv")
  )

  write_tsv(celltype_summary, celltype_summary_out)

  module_summary <- assigned_long %>%
    group_by(region, module) %>%
    summarise(
      n_tfs_assigned = n(),
      n_celltype_region = n_distinct(celltype_region),
      median_max_fc = median(max_fc, na.rm = TRUE),
      mean_max_fc = mean(max_fc, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(module)

  module_summary_out <- file.path(
    out_dir,
    paste0(region_name, "_FC_assignment_module_summary.tsv")
  )

  write_tsv(module_summary, module_summary_out)

  cat("\nSaved:\n")
  cat("  Observed matrix:", observed_out, "\n")
  cat("  Expected matrix:", expected_out, "\n")
  cat("  FC matrix:", fc_out, "\n")
  cat("  All assignment:", all_assignment_out, "\n")
  cat("  Assigned TFs:", assigned_out, "\n")
  cat("  Module counts:", module_counts_out, "\n")
  cat("  Celltype summary:", celltype_summary_out, "\n")
  cat("  Module summary:", module_summary_out, "\n")

  cat("\nAssignment summary:\n")
  cat("  Total cell-TFs:", nrow(assignment_df), "\n")
  cat("  Assigned cell-TFs:", nrow(assigned_long), "\n")
  cat("  Assigned fraction:", round(nrow(assigned_long) / nrow(assignment_df), 4), "\n")

  region_summaries[[region_name]] <- data.frame(
    region = region_name,
    n_celltype_region = length(matrix_list),
    n_cell_tf_total = nrow(assignment_df),
    n_assigned = nrow(assigned_long),
    assigned_fraction = nrow(assigned_long) / nrow(assignment_df),
    out_dir = out_dir,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# Overall summary
# ============================================================

if (length(region_summaries) > 0) {
  summary_df <- bind_rows(region_summaries)

  summary_out <- file.path(
    OUT_ROOT,
    if (nzchar(selected_region)) paste0(selected_region, "_FC_TF_assignment_summary.tsv") else "region_FC_TF_assignment_summary.tsv"
  )

  write_tsv(summary_df, summary_out)

  cat("\nSaved overall summary:\n")
  cat("  ", summary_out, "\n")
}

cat("\nAll done.\n")
