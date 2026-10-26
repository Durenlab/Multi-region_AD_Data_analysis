#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(ggplot2)
  library(ggrepel)
})

project_root <- Sys.getenv("AD_MANUSCRIPT_WORK_ROOT", unset = getwd())
data_root <- Sys.getenv("REGA_DATA_ROOT", unset = dirname(project_root))
separate_root <- Sys.getenv("REGA_SEPARATE_ROOT", unset = file.path(data_root, "REGA_separate"))
shared_root <- Sys.getenv("REGA_SHARED_ROOT", unset = file.path(data_root, "REGA_shared"))
snp_file <- Sys.getenv("AD_SNP_TF_LIST", unset = file.path(project_root, "experiments/TF_prioritize/TF_with_SNP/TFs_with_SNP.txt"))
driver_overlap_file <- file.path(project_root, "manuscript/supp/TFs_with_SNP_driver_TF_overlap.tsv")
driver_tables <- list(
  cell_type = file.path(project_root, "manuscript/supp/Supplementary_Table_celltype_specific_driver_TFs.tsv"),
  cross_cell_type = file.path(project_root, "manuscript/supp/Supplementary_Table_cross_celltype_driver_TFs.tsv")
)
output_root <- Sys.getenv("TF_VOLCANO_OUTPUT_ROOT", unset = file.path(project_root, "outputs", "tf_activity_expression_comparison"))
figures_only <- tolower(Sys.getenv("TF_VOLCANO_FIGURES_ONLY", "false")) %in% c("1", "true", "yes")
pair_only <- tolower(Sys.getenv("TF_VOLCANO_PAIR_ONLY", "false")) %in% c("1", "true", "yes")
target_pair <- Sys.getenv("TF_VOLCANO_CELLTYPE_REGION", unset = "")
target_scope <- Sys.getenv("TF_VOLCANO_SCOPE", unset = "")
figure_width <- as.numeric(Sys.getenv("TF_VOLCANO_FIGURE_WIDTH", "6.5"))
figure_height <- as.numeric(Sys.getenv("TF_VOLCANO_FIGURE_HEIGHT", "7"))
png_dpi <- as.numeric(Sys.getenv("TF_VOLCANO_PNG_DPI", "300"))
axis_text_size <- as.numeric(Sys.getenv("TF_VOLCANO_AXIS_TEXT_SIZE", "12"))
axis_title_size <- as.numeric(Sys.getenv("TF_VOLCANO_AXIS_TITLE_SIZE", "14"))
plot_title_size <- as.numeric(Sys.getenv("TF_VOLCANO_PLOT_TITLE_SIZE", "13"))
plot_subtitle_size <- as.numeric(Sys.getenv("TF_VOLCANO_PLOT_SUBTITLE_SIZE", "16"))
tf_label_size <- as.numeric(Sys.getenv("TF_VOLCANO_TF_LABEL_SIZE", "2.8"))
highlight_tfs <- toupper(trimws(strsplit(Sys.getenv("TF_VOLCANO_FOCUS_TFS", ""), ",", fixed = TRUE)[[1]]))
highlight_tfs <- highlight_tfs[nzchar(highlight_tfs)]
legend_text_size <- as.numeric(Sys.getenv("TF_VOLCANO_LEGEND_TEXT_SIZE", "10"))
caption_size <- as.numeric(Sys.getenv("TF_VOLCANO_CAPTION_SIZE", "8.5"))
empty_label_size <- as.numeric(Sys.getenv("TF_VOLCANO_EMPTY_LABEL_SIZE", "6"))
inputs <- list(
  cell_type = list(
    activity = file.path(separate_root, "07_tf/results/activity_pathology/all_celltype_region_TF_pathology_association.tsv"),
    expression = file.path(separate_root, "07_tf/results/expression_pathology/all_celltype_region_TF_pathology_association.tsv")
  ),
  cross_cell_type = list(
    activity = file.path(shared_root, "part7_TF/04TF_pathology/all_celltype_region_TF_pathology_association.tsv"),
    expression = file.path(separate_root, "07_tf/results/expression_pathology/all_celltype_region_TF_pathology_association.tsv")
  )
)
required <- c("celltype_region", "TF", "lm_pval_pathology", "lm_fdr_pathology",
              "mean_nonAD", "mean_lateAD", "mannwhitney_lateAD_vs_nonAD_pval",
              "mannwhitney_lateAD_vs_nonAD_fdr")
needs_driver_tables <- !figures_only || pair_only
all_inputs <- c(vapply(inputs, function(x) x$activity, character(1)),
                vapply(inputs, function(x) x$expression, character(1)),
                if (needs_driver_tables) unlist(driver_tables, use.names = FALSE),
                snp_file)
for (path in all_inputs) if (!file.exists(path)) stop("Missing input: ", path)

ad_snp_tfs <- unique(toupper(trimws(read_lines(snp_file))))
ad_snp_tfs <- ad_snp_tfs[!is.na(ad_snp_tfs) & nzchar(ad_snp_tfs)]

driver_overlap <- data.frame()
if (needs_driver_tables) {
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
}
cat("SNP-listed TFs:", length(ad_snp_tfs), "\n")

significant_color <- "#E34A4A"
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
    labs(title = title, subtitle = method, x = "Log2 fold change (lateAD − nonAD)", y = "-log10(minimum BH FDR)") +
    theme_classic(base_size = 12) +
    theme(aspect.ratio = 1,
          axis.text = element_text(size = axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.title = element_text(size = plot_title_size),
          plot.subtitle = element_text(size = plot_subtitle_size, face = "bold"))
}

make_volcano <- function(data, title, method, y_limits = NULL) {
  if (method == "Expression") {
    data$mean_nonAD_plot <- as.numeric(data$expression_mean_nonAD)
    data$mean_lateAD_plot <- as.numeric(data$expression_mean_lateAD)
    data$lm_fdr_plot <- as.numeric(data$expression_lm_fdr)
    data$mw_fdr_plot <- as.numeric(data$expression_mw_fdr)
    data$lm_p_plot <- as.numeric(data$expression_lm_pval)
    data$mw_p_plot <- as.numeric(data$expression_mw_pval)
  } else {
    data$mean_nonAD_plot <- as.numeric(data$mean_nonAD)
    data$mean_lateAD_plot <- as.numeric(data$mean_lateAD)
    data$lm_fdr_plot <- as.numeric(data$lm_fdr_pathology)
    data$mw_fdr_plot <- as.numeric(data$mannwhitney_lateAD_vs_nonAD_fdr)
    data$lm_p_plot <- as.numeric(data$lm_pval_pathology)
    data$mw_p_plot <- as.numeric(data$mannwhitney_lateAD_vs_nonAD_pval)
  }
  data$significant_lm <- is.finite(data$lm_fdr_plot) & data$lm_fdr_plot < 0.05
  data$significant_mw <- is.finite(data$mw_fdr_plot) & data$mw_fdr_plot < 0.05
  data$significant <- data$significant_lm | data$significant_mw
  data$fdr_plot <- pmin(data$lm_fdr_plot, data$mw_fdr_plot, na.rm = TRUE)
  data$fdr_plot[!is.finite(data$fdr_plot)] <- NA_real_
  data <- data[
    is.finite(data$mean_nonAD_plot) & is.finite(data$mean_lateAD_plot) &
      is.finite(data$fdr_plot) & data$fdr_plot >= 0 &
      data$fdr_plot <= 1,
    , drop = FALSE
  ]
  if (!nrow(data)) return(empty_volcano(title, method))

  data$x_value <- (data$mean_lateAD_plot - data$mean_nonAD_plot) * 1.442695
  data$y_value <- -log10(pmax(data$fdr_plot, .Machine$double.xmin))
  data$has_AD_SNP <- trimws(as.character(data$TF)) %in% ad_snp_tfs
  data$plot_color <- ifelse(data$significant, "Significant", "Not significant")
  data$plot_label <- ifelse(data$has_AD_SNP, paste0("\u2605 ", data$TF), data$TF)
  data$highlight_requested <- toupper(trimws(as.character(data$TF))) %in% highlight_tfs
  data$x_value <- pmax(pmin(data$x_value, 6), -6)

  # Consolidate only points with identical plotted coordinates so the FDR
  # threshold remains an exact visual boundary.
  data$x_plot <- data$x_value
  data$y_plot <- data$y_value
  point_key <- paste(data$x_plot, data$y_plot, sep = "\r")
  point_data <- data[!duplicated(point_key), , drop = FALSE]
  point_data$significant <- vapply(point_key[!duplicated(point_key)], function(key)
    any(data$significant[point_key == key]), logical(1))
  point_data$has_AD_SNP <- vapply(point_key[!duplicated(point_key)], function(key)
    any(data$has_AD_SNP[point_key == key]), logical(1))
  point_data$plot_color <- ifelse(point_data$significant, "Significant", "Not significant")

  p <- ggplot(data, aes(x = x_plot, y = y_plot)) +
    # SNP status is marked with a star before the label, matching the reference.
    geom_point(data = point_data, aes(color = plot_color), size = 1.4, alpha = 0.9) +
    scale_color_manual(values = c("Significant" = significant_color,
                                  "Not significant" = nonsignificant_color),
                       breaks = c("Significant", "Not significant"), drop = TRUE) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    coord_cartesian(xlim = c(-6, 6), ylim = y_limits) +
    labs(title = title, subtitle = method, x = "Log2FC", y = "-log10(minimum BH FDR)", color = NULL,
         caption = "Red = LM or lateAD vs nonAD Mann–Whitney BH FDR < 0.05; ★ = AD-SNP TF") +
    theme_classic(base_size = 12) +
    theme(aspect.ratio = 1,
          axis.text = element_text(size = axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.title = element_text(size = plot_title_size),
          plot.subtitle = element_text(size = plot_subtitle_size, face = "bold"),
          legend.position = "bottom",
          legend.direction = "horizontal", legend.text = element_text(size = legend_text_size),
          plot.caption = element_text(size = caption_size, hjust = 0))

  # Label every significant or SNP-listed TF; label color follows FDR status.
  labels <- data[data$significant | data$has_AD_SNP | data$highlight_requested, , drop = FALSE]
  repel <- list(
    inherit.aes = FALSE, size = tf_label_size, max.overlaps = Inf,
    box.padding = 0.65, point.padding = 0.3, force = 8, force_pull = 0.4,
    direction = "both", max.iter = 12000, max.time = 1.5, seed = 2026,
    min.segment.length = 0, segment.color = "grey50", segment.size = 0.25
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
  fdr_values <- fdr_values[is.finite(fdr_values) & fdr_values >= 0]
  if (!length(fdr_values)) return(c(0, -log10(0.05) * 1.1))
  upper <- max(-log10(max(min(fdr_values), .Machine$double.xmin)), -log10(0.05)) * 1.08
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
