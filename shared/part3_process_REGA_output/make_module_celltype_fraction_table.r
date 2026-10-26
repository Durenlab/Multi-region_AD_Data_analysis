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
})

# ============================================================
# Tabulate cell-type composition of assigned genes per module
#
# For each brain region:
#   one horizontal stacked barplot
#   y-axis: M1-M20
#   fill: cell type
#   value: fraction of assigned genes in that module
# ============================================================

# ============================================================
# Paths and settings
# ============================================================

BASE_DIR <- shared_root

INPUT_ROOT <- file.path(
  BASE_DIR,
  "part3_process_REGA_output/10-1gene_module_assignment_region_shareH_FC"
)

OUT_DIR <- file.path(
  BASE_DIR,
  "part3_process_REGA_output/11module_celltype_fraction_tables"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

regions <- c("PFC", "TH", "AG", "EC", "HC", "MTC")
modules <- paste0("M", 1:20)

celltype_order <- c(
  "Ast",
  "Exc",
  "Inh",
  "Mic_Immune",
  "Oli",
  "OPC"
)

# Fixed colors for 7 cell types
celltype_colors <- c(
    "Exc"= "#3182bd",
    "Inh"="#ff7f0e",
    "Ast"= "#31a354",
    "Oli"="#e41a1c",
    "OPC"= "#9467bd",
    "Mic_Immune"= "#9c755f"
)

# ============================================================
# Helper functions
# ============================================================

parse_celltype <- function(celltype_region, region) {
  # Examples:
  #   Ast_PFC -> Ast
  #   Mic_Immune_PFC -> Mic_Immune
  suffix <- paste0("_", region)
  sub(paste0(suffix, "$"), "", celltype_region)
}

tabulate_one_region <- function(region) {

  cat("\n============================================================\n")
  cat("Region:", region, "\n")
  cat("============================================================\n")

  in_file <- file.path(
    INPUT_ROOT,
    region,
    paste0(region, "_ALL_celltype_region_gene_module_assignments_FC.tsv")
  )

  if (!file.exists(in_file)) {
    warning("Missing assignment file: ", in_file)
    return(NULL)
  }

  df <- read_tsv(in_file, show_col_types = FALSE)

  required_cols <- c("region", "celltype_region", "module", "gene")
  missing_cols <- setdiff(required_cols, colnames(df))

  if (length(missing_cols) > 0) {
    stop("Missing columns in ", in_file, ": ", paste(missing_cols, collapse = ", "))
  }

  df <- df %>%
    filter(region == !!region) %>%
    mutate(
      module = factor(module, levels = modules),
      celltype = parse_celltype(celltype_region, .env$region),
      celltype = factor(celltype, levels = celltype_order)
    ) %>%
    filter(!is.na(module), !is.na(celltype))

  if (nrow(df) == 0) {
    warning("No assigned genes found for region: ", region)
    return(NULL)
  }

  # ------------------------------------------------------------
  # Count assigned genes by module and cell type
  # ------------------------------------------------------------

  count_df <- df %>%
    count(module, celltype, name = "n_genes") %>%
    complete(
      module = factor(modules, levels = modules),
      celltype = factor(celltype_order, levels = celltype_order),
      fill = list(n_genes = 0)
    ) %>%
    group_by(module) %>%
    mutate(
      module_total = sum(n_genes),
      fraction = ifelse(module_total > 0, n_genes / module_total, 0),
      percent = 100 * fraction
    ) %>%
    ungroup()

  # Save table
  table_out <- file.path(
    OUT_DIR,
    paste0(region, "_module_celltype_gene_fraction.tsv")
  )

  write_tsv(count_df, table_out)

  cat("Saved table:", table_out, "\n")
  return(data.frame(
    region = region,
    n_assigned_genes = nrow(df),
    table = table_out,
    stringsAsFactors = FALSE
  ))
}

# ============================================================
# Run all regions
# ============================================================

summary_list <- list()

for (region in regions) {
  res <- tabulate_one_region(region)

  if (!is.null(res)) {
    summary_list[[region]] <- res
  }
}

if (length(summary_list) > 0) {
  summary_df <- bind_rows(summary_list)

  summary_out <- file.path(
    OUT_DIR,
    "module_celltype_fraction_table_summary.tsv"
  )

  write_tsv(summary_df, summary_out)

  cat("\nSaved summary:", summary_out, "\n")
}

cat("\nAll done.\n")