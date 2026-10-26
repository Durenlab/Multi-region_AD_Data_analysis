#!/usr/bin/env Rscript
rega_separate_root <- Sys.getenv("REGA_SEPARATE_ROOT")
if (!nzchar(rega_separate_root)) stop("Set REGA_SEPARATE_ROOT")
suppressPackageStartupMessages({ library(readr); library(dplyr) })

base_dir <- Sys.getenv(
  "REGA_PATHOLOGY_ROOT",
  unset = file.path(rega_separate_root, "02_run_rega/pathology_association")
)
input_file <- file.path(base_dir, "results", "All_four_tests_H_pathology_results.tsv")
if (!file.exists(input_file)) {
  input_file <- file.path(base_dir, "All_four_tests_H_pathology_results.tsv")
}
figure_dir <- file.path(base_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

results <- read_tsv(input_file, col_types = cols(slope = col_double()), show_col_types = FALSE)
required <- c("test", "celltype_region", "celltype", "region", "Module", "Comparison", "pval_adj")
missing <- setdiff(required, names(results))
if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))

candidate_rows <- results %>%
  filter(
    (test == "Linear_Regression" & is.finite(pval_adj) & pval_adj < 0.05) |
    (test == "Pairwise_MannWhitneyU" & Comparison == "nonAD vs lateAD" &
       is.finite(pval_adj) & pval_adj < 0.05)
  )

selected <- candidate_rows %>%
  group_by(celltype_region, celltype, region, Module) %>%
  summarise(best_FDR = min(pval_adj, na.rm = TRUE), .groups = "drop") %>%
  arrange(best_FDR, celltype, region, Module) %>%
  mutate(module_label = ifelse(grepl("^M[0-9]+$", Module), Module, paste0("M", Module)),
         minus_log10_FDR = -log10(best_FDR))

if (!nrow(selected)) stop("No modules pass the requested FDR criteria.")

heatmap <- selected %>% select(celltype, region, module_label, best_FDR, minus_log10_FDR)
write_tsv(heatmap, file.path(figure_dir, "selected_modules_heatmap_plot_data.tsv"))
writeLines(paste(selected$celltype_region, selected$module_label, sep = "_"),
           file.path(figure_dir, "selected_modules.txt"))
writeLines(paste(selected$celltype_region, selected$module_label, sep = "_"),
           file.path(base_dir, "selected_module.txt"))
write_tsv(candidate_rows %>% arrange(celltype_region, Module, test, Comparison),
          file.path(base_dir, "selected_module_significant_test_details.tsv"))
writeLines(c("x", sort(unique(selected$celltype_region))),
           file.path(base_dir, "selected_context.txt"))
