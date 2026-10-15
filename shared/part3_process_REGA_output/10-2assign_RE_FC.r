#!/usr/bin/env Rscript
shared_root <- Sys.getenv("REGA_SHARED_ROOT")
if (!nzchar(shared_root)) {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(file_arg)) stop("Set REGA_SHARED_ROOT when sourcing this script")
  shared_root <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
  while (!(basename(shared_root) %in% c("REGA_shared", "shared")) && dirname(shared_root) != shared_root) shared_root <- dirname(shared_root)
  if (!(basename(shared_root) %in% c("REGA_shared", "shared"))) stop("Could not find shared analysis root")
}
separate_root <- Sys.getenv("REGA_SEPARATE_ROOT", "/project/zduren/durenlab/naqing_palmetto/REGA_separate")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

# ============================================================
# Region-wise RE assignment using observed / expected WBV
#
# Logic:
#   For each brain region:
#     1. Read all celltype-region RE matrices
#     2. Rename rows as {celltype_region}__{RE}
#     3. Concatenate rows across cell types
#     4. Compute expected matrix:
#          expected[i,j] = rowSum[i] * colSum[j] / totalSum
#     5. Compute fold change:
#          FC[i,j] = observed[i,j] / expected[i,j]
#     6. Assign each RE to the module with max FC
#     7. Keep only REs with max FC >= FC_THRESHOLD
#
# Input:
#   08extracted_matrices_region_shareH/{REGION}/*_BV_RE_by_module.tsv
#
# Output:
#   09RE_module_assignment_region_shareH_FC/{REGION}/...
# ============================================================

# ============================================================
# User config
# ============================================================

BASE_DIR <- shared_root

INPUT_ROOT <- file.path(
  BASE_DIR,
  "part3_process_REGA_output/09extracted_matrices_region_shareH"
)

OUT_ROOT <- file.path(
  BASE_DIR,
  "part3_process_REGA_output/10-2RE_module_assignment_region_shareH_FC"
)

regions <- c("PFC", "TH", "AG", "EC", "HC", "MTC")
selected_region <- Sys.getenv("REGA_REGION", unset = "")
if (nzchar(selected_region)) {
  if (!selected_region %in% regions) stop("Unknown REGA_REGION: ", selected_region)
  regions <- selected_region
}

FC_THRESHOLD <- 2

dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Helper functions
# ============================================================

parse_cell_from_re_file <- function(fn_base) {
  # Example:
  #   Ast_AG_BV_RE_by_module.tsv -> Ast_AG
  sub("_BV_RE_by_module\\.tsv$", "", fn_base)
}

split_cell_re <- function(cell_re_name) {
  # We name rows as {celltype_region}__{RE}
  # RE itself may contain colon/dash, but not "__"
  parts <- strsplit(cell_re_name, "__", fixed = TRUE)[[1]]

  if (length(parts) < 2) {
    return(list(celltype_region = NA_character_, RE = cell_re_name))
  }

  celltype_region <- paste(parts[-length(parts)], collapse = "__")
  RE <- parts[length(parts)]

  list(celltype_region = celltype_region, RE = RE)
}

make_expected_matrix <- function(M) {
  row_sum <- rowSums(M, na.rm = TRUE)
  col_sum <- colSums(M, na.rm = TRUE)
  total_sum <- sum(M, na.rm = TRUE)

  if (!is.finite(total_sum) || total_sum <= 0) {
    stop("Total WBV sum is zero or non-finite; cannot compute expected matrix.")
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
  fc[!is.finite(fc)] <- NA

  max_idx <- max.col(fc, ties.method = "first")
  max_fc <- fc[cbind(seq_len(nrow(fc)), max_idx)]
  best_module <- colnames(fc)[max_idx]

  assigned <- ifelse(
    is.finite(max_fc) & max_fc >= fc_threshold,
    best_module,
    NA_character_
  )

  result <- data.frame(
    cell_RE = rownames(M),
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
# Main function
# ============================================================

run_one_region <- function(region) {

  cat("\n============================================================\n")
  cat("Region:", region, "\n")
  cat("============================================================\n")

  input_dir <- file.path(INPUT_ROOT, region)
  out_dir <- file.path(OUT_ROOT, region)

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  if (!dir.exists(input_dir)) {
    warning("Missing input dir: ", input_dir)
    return(NULL)
  }

  files <- list.files(
    input_dir,
    pattern = "_BV_RE_by_module\\.tsv$",
    full.names = TRUE
  )

  if (length(files) == 0) {
    warning("No RE-by-module files found in: ", input_dir)
    return(NULL)
  }

  cat("Found RE files:", length(files), "\n")

  matrix_list <- list()
  row_meta_list <- list()
  module_order <- NULL

  # ============================================================
  # Read and concatenate RE matrices
  # ============================================================

  for (f in files) {

    fname <- basename(f)
    celltype_region <- parse_cell_from_re_file(fname)

    cat("Reading:", fname, "->", celltype_region, "\n")

    df <- read_tsv(f, show_col_types = FALSE)

    if (ncol(df) < 2) {
      warning("Skipping file with <2 columns: ", f)
      next
    }

    re_col <- colnames(df)[1]
    res <- as.character(df[[re_col]])

    M <- as.matrix(df[, -1, drop = FALSE])
    mode(M) <- "numeric"

    M[!is.finite(M)] <- 0
    M[M < 0] <- 0

    rownames(M) <- paste0(celltype_region, "__", res)

    if (is.null(module_order)) {
      module_order <- colnames(M)
    }

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
      cell_RE = rownames(M),
      region = region,
      celltype_region = celltype_region,
      RE = res,
      stringsAsFactors = FALSE
    )
  }

  if (length(matrix_list) == 0) {
    warning("No valid matrices loaded for region: ", region)
    return(NULL)
  }

  observed <- do.call(rbind, matrix_list)
  row_meta <- bind_rows(row_meta_list)

  cat("Concatenated observed matrix:\n")
  cat("  rows:", nrow(observed), "\n")
  cat("  cols:", ncol(observed), "\n")

  # ============================================================
  # Save observed matrix
  # ============================================================

  observed_out <- file.path(
    out_dir,
    paste0(region, "_observed_WBV_cellRE_by_module.tsv")
  )

  write_tsv(
    as.data.frame(observed) %>%
      tibble::rownames_to_column("cell_RE"),
    observed_out
  )

  # ============================================================
  # Compute expected and FC, then assign
  # ============================================================

  res <- assign_by_fc(
    observed,
    fc_threshold = FC_THRESHOLD
  )

  assignment_raw <- res$assignment
  fc <- res$fc
  expected <- res$expected

  assignment_df <- assignment_raw %>%
    left_join(row_meta, by = "cell_RE") %>%
    select(
      region,
      celltype_region,
      RE,
      cell_RE,
      assigned_module,
      max_fc,
      observed_at_assigned,
      expected_at_assigned
    ) %>%
    mutate(
      assigned = !is.na(assigned_module)
    )

  # Long table of assigned REs only
  assigned_long <- assignment_df %>%
    filter(assigned) %>%
    transmute(
      region = region,
      celltype_region = celltype_region,
      module = assigned_module,
      RE = RE,
      cell_RE = cell_RE,
      max_fc = max_fc,
      observed = observed_at_assigned,
      expected = expected_at_assigned
    ) %>%
    arrange(region, celltype_region, module, desc(max_fc))

  # ============================================================
  # Save assignment outputs
  # ============================================================

  all_assignment_out <- file.path(
    out_dir,
    paste0(region, "_ALL_cellRE_FC_assignment_with_unassigned.tsv")
  )

  assigned_out <- file.path(
    out_dir,
    paste0(region, "_ALL_celltype_region_RE_module_assignments_FC.tsv")
  )

  write_tsv(assignment_df, all_assignment_out)
  write_tsv(assigned_long, assigned_out)

  # ============================================================
  # Save expected and FC matrices
  # ============================================================

  expected_out <- file.path(
    out_dir,
    paste0(region, "_expected_WBV_cellRE_by_module.tsv")
  )

  fc_out <- file.path(
    out_dir,
    paste0(region, "_FC_observed_over_expected_cellRE_by_module.tsv")
  )

  write_tsv(
    as.data.frame(expected) %>%
      tibble::rownames_to_column("cell_RE"),
    expected_out
  )

  write_tsv(
    as.data.frame(fc) %>%
      tibble::rownames_to_column("cell_RE"),
    fc_out
  )

  # ============================================================
  # Summary tables
  # ============================================================

  module_counts <- assigned_long %>%
    count(region, celltype_region, module, name = "n_REs") %>%
    arrange(region, celltype_region, module)

  module_counts_out <- file.path(
    out_dir,
    paste0(region, "_FC_assignment_module_counts.tsv")
  )

  write_tsv(module_counts, module_counts_out)

  celltype_summary <- assignment_df %>%
    group_by(region, celltype_region) %>%
    summarise(
      n_RE_total = n(),
      n_RE_assigned = sum(assigned),
      assigned_fraction = n_RE_assigned / n_RE_total,
      median_max_fc = median(max_fc, na.rm = TRUE),
      mean_max_fc = mean(max_fc, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_RE_assigned))

  celltype_summary_out <- file.path(
    out_dir,
    paste0(region, "_FC_assignment_celltype_summary.tsv")
  )

  write_tsv(celltype_summary, celltype_summary_out)

  module_summary <- assigned_long %>%
    group_by(region, module) %>%
    summarise(
      n_RE_assigned = n(),
      n_celltype_region = n_distinct(celltype_region),
      median_max_fc = median(max_fc, na.rm = TRUE),
      mean_max_fc = mean(max_fc, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(module)

  module_summary_out <- file.path(
    out_dir,
    paste0(region, "_FC_assignment_module_summary.tsv")
  )

  write_tsv(module_summary, module_summary_out)

  cat("\nSaved:\n")
  cat("  Observed matrix:", observed_out, "\n")
  cat("  Expected matrix:", expected_out, "\n")
  cat("  FC matrix:", fc_out, "\n")
  cat("  All assignment:", all_assignment_out, "\n")
  cat("  Assigned REs:", assigned_out, "\n")
  cat("  Module counts:", module_counts_out, "\n")
  cat("  Celltype summary:", celltype_summary_out, "\n")
  cat("  Module summary:", module_summary_out, "\n")

  cat("\nAssignment summary:\n")
  cat("  Total cell-REs:", nrow(assignment_df), "\n")
  cat("  Assigned cell-REs:", nrow(assigned_long), "\n")
  cat("  Assigned fraction:", round(nrow(assigned_long) / nrow(assignment_df), 4), "\n")

  return(data.frame(
    region = region,
    n_celltype_region = length(matrix_list),
    n_cell_RE_total = nrow(assignment_df),
    n_assigned = nrow(assigned_long),
    assigned_fraction = nrow(assigned_long) / nrow(assignment_df),
    out_file = assigned_out,
    stringsAsFactors = FALSE
  ))
}

# ============================================================
# Run all regions
# ============================================================

summary_list <- list()

for (region in regions) {
  res <- run_one_region(region)

  if (!is.null(res)) {
    summary_list[[region]] <- res
  }
}

if (length(summary_list) > 0) {

  summary_df <- bind_rows(summary_list)

  summary_out <- file.path(
    OUT_ROOT,
    if (nzchar(selected_region)) paste0(selected_region, "_FC_RE_assignment_summary.tsv") else "region_FC_RE_assignment_summary.tsv"
  )

  write_tsv(summary_df, summary_out)

  cat("\nSaved overall summary:\n")
  cat("  ", summary_out, "\n")
}

cat("\nAll done.\n")
