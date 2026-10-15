#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(ggplot2)
  library(ggrepel)
})

project_root <- "/project/zduren/durenlab/naqing_palmetto/Alzheimer_project_explore"
driver_root <- "/project/zduren/durenlab/naqing_palmetto/AD_Driver_TF"
snp_file <- file.path(project_root, "manuscript/TFs_with_SNP.txt")
driver_overlap_file <- file.path(project_root, "manuscript/supp/TFs_with_SNP_driver_TF_overlap.tsv")
driver_tables <- list(
  cell_type = file.path(project_root, "manuscript/supp/Supplementary_Table_celltype_specific_driver_TFs.tsv"),
  cross_cell_type = file.path(project_root, "manuscript/supp/Supplementary_Table_cross_celltype_driver_TFs.tsv")
)
output_root <- Sys.getenv("EXP009_OUTPUT_ROOT", unset = file.path(project_root, "experiments/EXP009_results"))
figures_only <- tolower(Sys.getenv("EXP009_FIGURES_ONLY", "false")) %in% c("1", "true", "yes")
pair_only <- tolower(Sys.getenv("EXP009_PAIR_ONLY", "false")) %in% c("1", "true", "yes")
target_pair <- Sys.getenv("EXP009_CELLTYPE_REGION", unset = "")
target_scope <- Sys.getenv("EXP009_SCOPE", unset = "")
figure_width <- as.numeric(Sys.getenv("EXP009_FIGURE_WIDTH", "6.5"))
figure_height <- as.numeric(Sys.getenv("EXP009_FIGURE_HEIGHT", "7"))
png_dpi <- as.numeric(Sys.getenv("EXP009_PNG_DPI", "300"))
axis_text_size <- as.numeric(Sys.getenv("EXP009_AXIS_TEXT_SIZE", "26.4"))
axis_title_size <- as.numeric(Sys.getenv("EXP009_AXIS_TITLE_SIZE", "22"))
plot_title_size <- as.numeric(Sys.getenv("EXP009_PLOT_TITLE_SIZE", "22"))
plot_subtitle_size <- as.numeric(Sys.getenv("EXP009_PLOT_SUBTITLE_SIZE", "18"))
tf_label_size <- as.numeric(Sys.getenv("EXP009_TF_LABEL_SIZE", "5.4"))
max_labels_per_plot <- as.integer(Sys.getenv("EXP009_MAX_LABELS_PER_PLOT", "2"))
highlight_tfs <- toupper(trimws(strsplit(Sys.getenv("EXP009_HIGHLIGHT_TFS", ""), ",", fixed = TRUE)[[1]]))
highlight_tfs <- highlight_tfs[nzchar(highlight_tfs)]
legend_text_size <- as.numeric(Sys.getenv("EXP009_LEGEND_TEXT_SIZE", as.character(axis_text_size)))
caption_size <- as.numeric(Sys.getenv("EXP009_CAPTION_SIZE", "12"))
empty_label_size <- as.numeric(Sys.getenv("EXP009_EMPTY_LABEL_SIZE", "6"))
inputs <- list(
  cell_type = list(
    activity = file.path(driver_root, "cell_type/TF_driver.tsv"),
    expression = file.path(driver_root, "cell_type/GEX_driver.tsv")
  ),
  cross_cell_type = list(
    activity = file.path(driver_root, "cross_celltype/TF_driver.tsv"),
    expression = file.path(driver_root, "cross_celltype/GEX_driver.tsv")
  )
)
required <- c("celltype_region", "TF", "lm_pval_pathology", "lm_fdr_pathology",
              "mean_nonAD", "mean_lateAD", "mannwhitney_lateAD_vs_nonAD_pval",
              "mannwhitney_lateAD_vs_nonAD_fdr")
all_inputs <- c(vapply(inputs, function(x) x$activity, character(1)),
                vapply(inputs, function(x) x$expression, character(1)),
                unlist(driver_tables, use.names = FALSE), snp_file)
for (path in all_inputs) if (!file.exists(path)) stop("Missing input: ", path)

ad_snp_tfs <- unique(toupper(trimws(read_lines(snp_file))))
ad_snp_tfs <- ad_snp_tfs[!is.na(ad_snp_tfs) & nzchar(ad_snp_tfs)]

driver_overlap <- data.frame()
for (scope_name in names(driver_tables)) {
  driver_data <- read_tsv(driver_tables[[scope_name]], show_col_types = FALSE,
                          col_types = cols(.default = col_character()))
  needed <- c("celltype_region", "TF", "TF_assignment", "significant_test")
  if (!all(needed %in% names(driver_data)))
    stop("Driver table lacks required columns: ", driver_tables[[scope_name]])
  matches <- driver_data[toupper(trimws(driver_data$TF)) %in% ad_snp_tfs, needed, drop = FALSE]
  matches$analysis_scope <- scope_name
  driver_overlap <- rbind(driver_overlap, matches)
}
driver_overlap <- unique(driver_overlap[c("analysis_scope", "celltype_region", "TF",
                                         "TF_assignment", "significant_test")])
if (!figures_only) write_tsv(driver_overlap, driver_overlap_file)
cat("New SNP TF list:", length(ad_snp_tfs), "TFs; overlaps:", nrow(driver_overlap),
    "driver rows across", length(unique(driver_overlap$celltype_region)), "celltype-region pairs\n")

significant_color <- "#D62728"
nonsignificant_color <- "#BDBDBD"

save_plot <- function(plot, stem) {
  ggsave(paste0(stem, ".png"), plot, width = figure_width, height = figure_height,
         dpi = png_dpi, bg = "white")
  ggsave(paste0(stem, ".pdf"), plot, width = figure_width, height = figure_height,
         device = grDevices::cairo_pdf, bg = "white")
}

empty_volcano <- function(title, method) {
  ggplot(data.frame(x = 0, y = 0), aes(x, y)) +
    geom_blank() +
    annotate("text", x = 0, y = 0, label = "No complete rows", size = empty_label_size) +
    labs(title = title, subtitle = method, x = "Log2FC", y = "-log10(minimum BH FDR)") +
    theme_classic(base_size = 11) +
    theme(aspect.ratio = 1,
          axis.text = element_text(size = axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.title = element_text(size = plot_title_size),
          plot.subtitle = element_text(size = plot_subtitle_size, face = "bold"))
}

make_volcano <- function(data, title, method, y_limits = NULL) {
  if (method == "Expression") {
    data$mean_nonAD_plot <- as.numeric(data$expression_mean_nonAD) + 1
    data$mean_lateAD_plot <- as.numeric(data$expression_mean_lateAD) + 1
    data$lm_fdr_plot <- as.numeric(data$expression_lm_fdr)
    data$mw_fdr_plot <- as.numeric(data$expression_mw_fdr)
  } else {
    data$mean_nonAD_plot <- as.numeric(data$mean_nonAD)
    data$mean_lateAD_plot <- as.numeric(data$mean_lateAD)
    data$lm_fdr_plot <- as.numeric(data$lm_fdr_pathology)
    data$mw_fdr_plot <- as.numeric(data$mannwhitney_lateAD_vs_nonAD_fdr)
  }
  data$significant_lm <- is.finite(data$lm_fdr_plot) & data$lm_fdr_plot < 0.05
  data$significant_mw <- is.finite(data$mw_fdr_plot) & data$mw_fdr_plot < 0.05
  data$significant <- data$significant_lm | data$significant_mw
  data$fdr_plot <- pmin(data$lm_fdr_plot, data$mw_fdr_plot, na.rm = TRUE)
  data$fdr_plot[!is.finite(data$fdr_plot)] <- NA_real_
  data <- data[
    is.finite(data$mean_nonAD_plot) & is.finite(data$mean_lateAD_plot) &
      data$mean_nonAD_plot > 0 & data$mean_lateAD_plot > 0 &
      is.finite(data$fdr_plot) & data$fdr_plot >= 0,
    , drop = FALSE
  ]
  if (!nrow(data)) return(empty_volcano(title, method))

  if (method == "Expression") {
    # Retain the plotted expression difference/scale used by the notebook.
    data$x_value <- (data$mean_lateAD_plot - data$mean_nonAD_plot) * 1.442695
  } else {
    data$x_value <- log2(data$mean_lateAD_plot / data$mean_nonAD_plot)
  }
  data$y_value <- -log10(pmax(data$fdr_plot, .Machine$double.xmin))
  data$has_AD_SNP <- trimws(as.character(data$TF)) %in% ad_snp_tfs
  data$plot_color <- ifelse(data$significant, "Significant", "Not significant")
  data$plot_label <- ifelse(data$has_AD_SNP, paste0("\u2605 ", data$TF), data$TF)
  data$x_value <- pmax(pmin(data$x_value, 6), -6)

  # Consolidate points that overlap at the figure's visual resolution.
  # Statistical values remain unchanged; only displayed coordinates are
  # rounded to 0.1 units. Keep all rows in the supplementary table and show
  # no more than two TF labels across the entire plot.
  data$x_plot <- round(data$x_value, 1)
  data$y_plot <- round(data$y_value, 1)
  point_key <- paste(data$x_plot, data$y_plot, sep = "\r")
  point_data <- data[!duplicated(point_key), , drop = FALSE]
  point_data$significant <- vapply(point_key[!duplicated(point_key)], function(key)
    any(data$significant[point_key == key]), logical(1))
  point_data$plot_color <- ifelse(point_data$significant, "Significant", "Not significant")

  p <- ggplot(data, aes(x = x_plot, y = y_plot)) +
    # Keep points at their exact estimates. AD-SNP stars are repelled as part
    # of their labels below and linked back to these points with leader lines.
    geom_point(data = point_data, aes(color = plot_color), size = 2.1, alpha = 0.9) +
    scale_color_manual(values = c("Significant" = significant_color,
                                  "Not significant" = nonsignificant_color),
                       breaks = c("Significant", "Not significant"), drop = TRUE) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    coord_cartesian(xlim = c(-6, 6), ylim = y_limits) +
    labs(title = title, subtitle = method, x = "Log2FC", y = "-log10(minimum BH FDR)", color = NULL,
         caption = "Red = LM or lateAD vs nonAD Mann–Whitney BH FDR < 0.05; ★ = AD-SNP TF; overlapping coordinates grouped at 0.1 units; at most 2 TF labels shown") +
    theme_classic(base_size = 11) +
    theme(aspect.ratio = 1,
          axis.text = element_text(size = axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.title = element_text(size = plot_title_size),
          plot.subtitle = element_text(size = plot_subtitle_size, face = "bold"),
          legend.position = "bottom",
          legend.direction = "horizontal", legend.text = element_text(size = legend_text_size),
          plot.caption = element_text(size = caption_size, hjust = 0))

  # Label at most two TFs per plot, prioritizing significant TFs and then
  # AD-SNP TFs; among comparable rows, show the smallest available BH FDRs.
  # Color follows significance: red for significant and grey for nonsignificant.
  data$highlight_requested <- toupper(trimws(as.character(data$TF))) %in% highlight_tfs
  labels <- data[data$significant | data$has_AD_SNP | data$highlight_requested, , drop = FALSE]
  if (nrow(labels) && max_labels_per_plot > 0) {
    labels <- labels[order(!labels$highlight_requested, !labels$significant,
                           !labels$has_AD_SNP,
                           labels$fdr_plot, labels$TF), , drop = FALSE]
    labels <- head(labels, max_labels_per_plot)
  } else {
    labels <- labels[FALSE, , drop = FALSE]
  }
  repel <- list(
    inherit.aes = FALSE, size = tf_label_size, max.overlaps = Inf,
    box.padding = 0.6, point.padding = 0.55, force = 3,
    min.segment.length = 0, segment.color = "grey50", segment.size = 0.3
  )
  if (nrow(labels)) {
    p <- p + do.call(ggrepel::geom_text_repel, c(list(
      data = labels, aes(x = x_plot, y = y_plot, label = plot_label, color = plot_color)
    ), repel))
  }
  p
}

shared_pair_y_limits <- function(data) {
  fdr_values <- c(
    suppressWarnings(as.numeric(data$lm_fdr_pathology)),
    suppressWarnings(as.numeric(data$mannwhitney_lateAD_vs_nonAD_fdr)),
    suppressWarnings(as.numeric(data$expression_lm_fdr)),
    suppressWarnings(as.numeric(data$expression_mw_fdr))
  )
  fdr_values <- fdr_values[is.finite(fdr_values) & fdr_values > 0]
  if (!length(fdr_values)) return(c(0, -log10(0.05) * 1.1))
  upper <- max(-log10(min(fdr_values)), -log10(0.05)) * 1.08
  c(0, upper)
}

scopes_to_run <- names(inputs)
if (nzchar(target_scope)) scopes_to_run <- intersect(scopes_to_run, target_scope)
for (scope in scopes_to_run) {
  activity_file <- inputs[[scope]]$activity
  expression_file <- inputs[[scope]]$expression
  output_dir <- file.path(output_root, scope)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  activity <- read_tsv(activity_file, show_col_types = FALSE,
                       col_types = cols(.default = col_character()))
  expression_lm <- read_tsv(expression_file, show_col_types = FALSE,
                            col_types = cols(.default = col_character()))
  for (nm in required) {
    if (!nm %in% names(activity)) stop("Missing activity column: ", nm)
    if (!nm %in% names(expression_lm)) stop("Missing expression column: ", nm)
  }

  # Preserve the first expression row per celltype_region/TF while allowing
  # either test to be unavailable when the other has an estimable result.
  expression_lm <- expression_lm[
    complete.cases(expression_lm[c("celltype_region", "TF", "mean_nonAD", "mean_lateAD")]),
    , drop = FALSE
  ]
  expression_join <- expression_lm[
    !duplicated(expression_lm[c("celltype_region", "TF")]),
    c("celltype_region", "TF", "lm_pval_pathology", "lm_fdr_pathology",
      "mean_nonAD", "mean_lateAD", "mannwhitney_lateAD_vs_nonAD_pval",
      "mannwhitney_lateAD_vs_nonAD_fdr"), drop = FALSE
  ]
  names(expression_join)[3:6] <- c(
    "expression_lm_pval", "expression_lm_fdr",
    "expression_mean_nonAD", "expression_mean_lateAD"
  )
  names(expression_join)[7:8] <- c("expression_mw_pval", "expression_mw_fdr")
  joined <- merge(activity, expression_join,
                  by = c("celltype_region", "TF"), all.x = TRUE, sort = FALSE)
  joined <- joined[complete.cases(joined[c("celltype_region", "TF", "mean_nonAD", "mean_lateAD",
                                           "expression_mean_nonAD", "expression_mean_lateAD")]), , drop = FALSE]
  if (!"celltype" %in% names(joined))
    joined$celltype <- sub("_[^_]+$", "", joined$celltype_region)
  if (!"region" %in% names(joined))
    joined$region <- sub("^.*_", "", joined$celltype_region)
  joined$significant_expression_lm <-
    suppressWarnings(as.numeric(joined$expression_lm_fdr)) < 0.05
  joined$significant_expression_mw <-
    suppressWarnings(as.numeric(joined$expression_mw_fdr)) < 0.05
  joined$significant_in_expression <- joined$significant_expression_lm | joined$significant_expression_mw
  joined$significant_activity_lm <-
    suppressWarnings(as.numeric(joined$lm_fdr_pathology)) < 0.05
  joined$significant_activity_mw <-
    suppressWarnings(as.numeric(joined$mannwhitney_lateAD_vs_nonAD_fdr)) < 0.05
  joined$significant_in_activity <- joined$significant_activity_lm | joined$significant_activity_mw
  joined$has_AD_SNP <- trimws(as.character(joined$TF)) %in% ad_snp_tfs
  if (!figures_only) write_tsv(joined, file.path(output_dir, "supplementary_table1.tsv"))

  pairs <- unique(joined[c("celltype_region", "celltype", "region")])
  pairs <- pairs[order(pairs$celltype_region), , drop = FALSE]
  if (nzchar(target_pair))
    pairs <- pairs[pairs$celltype_region == target_pair, , drop = FALSE]
  if (pair_only) {
    target_pairs <- unique(driver_overlap$celltype_region[driver_overlap$analysis_scope == scope])
    pairs <- pairs[pairs$celltype_region %in% target_pairs, , drop = FALSE]
    if (!nrow(pairs)) {
      cat(scope, "has no matching driver/SNP celltype-region pairs; skipping plots\n")
      next
    }
  }
  celltypes <- sort(unique(joined$celltype))
  regions <- sort(unique(joined$region))
  for (method in c("Expression", "Activity")) {
    method_dir <- file.path(output_dir, method)
    dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)
    for (i in seq_len(nrow(pairs))) {
      key <- pairs$celltype_region[[i]]
      pair_data <- joined[joined$celltype_region == key, , drop = FALSE]
      p <- make_volcano(pair_data, key, method,
                        y_limits = shared_pair_y_limits(pair_data))
      save_plot(p, file.path(method_dir, paste0(key, "_volcano")))
    }
    if (!pair_only && !nzchar(target_pair)) {
      for (ct in celltypes) {
        p <- make_volcano(joined[joined$celltype == ct, , drop = FALSE],
                          paste(ct, "across regions"), method)
        save_plot(p, file.path(method_dir, paste0(ct, "_across_regions_volcano")))
      }
      for (rg in regions) {
        p <- make_volcano(joined[joined$region == rg, , drop = FALSE],
                          paste(rg, "across cell types"), method)
        save_plot(p, file.path(method_dir, paste0(rg, "_across_celltypes_volcano")))
      }
    }
  }
  cat(scope, "joined rows:", nrow(joined),
      "AD-SNP TF rows:", sum(joined$has_AD_SNP),
      "celltype-region pairs:", nrow(pairs), "\n")
}
cat("Saved remade plots under:", output_root, "\n")
