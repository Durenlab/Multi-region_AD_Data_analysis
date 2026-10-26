#!/usr/bin/env Rscript
rega_separate_root <- Sys.getenv("REGA_SEPARATE_ROOT")
if (!nzchar(rega_separate_root)) stop("Set REGA_SEPARATE_ROOT")

# ============================================================
# Compare TF expression vs TF activity for pathology association
# (nonAD vs lateAD) across all celltype-region datasets.
#
# Per TF, two tests are run:
#   1) Linear model: value ~ Pathology
#   2) Mann-Whitney U / Wilcoxon rank-sum test
#
# IMPORTANT FOR THIS VERSION:
#   - Summaries use NUMBER of significant TFs, not fraction.
#   - Expression and activity use the same TF identities after removing
#     all-zero rows in either representation.
#   - For each test separately, a celltype-region is omitted from plots
#     when BOTH Expression and Activity have 0 significant TFs.
#   - Overall boxplots compare significant-TF counts across matched
#     celltype-region datasets.
#   - The paired plot additionally reports a paired Wilcoxon test across
#     matched celltype-region significant-TF counts.
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

# ------------------------------------------------------------
# 1. User settings
# ------------------------------------------------------------

tf_root <- file.path(rega_separate_root, "07_tf/results/activity/")

metadata_dir <- file.path(rega_separate_root, "01_prepare/data/metadata_by_group")

output_dir <- Sys.getenv(
  "REGA_TF_COMPARISON_OUTPUT_DIR",
  unset = file.path(rega_separate_root, "07_tf/results/expression_activity_counts_same_tf")
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

alpha <- 0.05

# Which p-value defines a significant TF:
#   "p_value" = nominal p < 0.05
#   "p_adj"   = BH-FDR < 0.05
significance_column <- "p_adj"

# Minimum usable individuals required in EACH pathology group for a TF.
min_n_per_group <- 3

group_ref <- "nonAD"
group_alt <- "lateAD"

# ------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------

read_tf_matrix <- function(file) {
  x <- read.csv(
    file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (ncol(x) < 2) {
    stop("Matrix has fewer than 2 columns: ", file)
  }

  first_col <- x[[1]]
  first_name <- colnames(x)[1]

  likely_rowname_col <-
    first_name %in% c("", "X", "row.names", "Row.names", "TF", "gene") ||
    (!is.numeric(first_col) && length(unique(first_col)) == nrow(x))

  if (likely_rowname_col) {
    tf_names <- as.character(first_col)
    x <- x[, -1, drop = FALSE]
  } else {
    tf_names <- rownames(x)
  }

  x[] <- lapply(x, function(z) suppressWarnings(as.numeric(z)))
  x <- as.data.frame(x, check.names = FALSE)
  rownames(x) <- make.unique(tf_names)
  colnames(x) <- trimws(colnames(x))

  x
}

read_metadata <- function(file) {
  meta <- read.delim(
    file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  required <- c("ROSMAP_IndividualID", "RNA.Class", "BrainRegion", "Pathology")
  missing_cols <- setdiff(required, colnames(meta))

  if (length(missing_cols) > 0) {
    stop(
      "Metadata missing required column(s) in ", file, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  meta %>%
    transmute(
      individual = trimws(as.character(ROSMAP_IndividualID)),
      celltype = trimws(as.character(RNA.Class)),
      region = trimws(as.character(BrainRegion)),
      Pathology = trimws(as.character(Pathology))
    ) %>%
    filter(Pathology %in% c(group_ref, group_alt)) %>%
    distinct(individual, .keep_all = TRUE)
}

run_tf_tests <- function(mat, meta, data_type, celltype, region, pair_name) {
  common_ids <- intersect(colnames(mat), meta$individual)

  if (length(common_ids) == 0) {
    warning("No overlapping individuals for ", pair_name, " / ", data_type)
    return(tibble())
  }

  meta_sub <- meta %>%
    filter(individual %in% common_ids) %>%
    mutate(Pathology = factor(Pathology, levels = c(group_ref, group_alt)))

  mat_sub <- mat[, meta_sub$individual, drop = FALSE]

  one_tf <- function(tf) {
    y <- as.numeric(mat_sub[tf, ])

    dat <- tibble(
      value = y,
      Pathology = meta_sub$Pathology
    ) %>%
      filter(!is.na(value), !is.na(Pathology))

    n_ref <- sum(dat$Pathology == group_ref)
    n_alt <- sum(dat$Pathology == group_alt)

    mean_ref <- if (n_ref > 0) mean(dat$value[dat$Pathology == group_ref]) else NA_real_
    mean_alt <- if (n_alt > 0) mean(dat$value[dat$Pathology == group_alt]) else NA_real_

    if (n_ref < min_n_per_group || n_alt < min_n_per_group) {
      return(tibble(
        TF = tf,
        celltype = celltype,
        region = region,
        celltype_region = pair_name,
        data_type = data_type,
        n_nonAD = n_ref,
        n_lateAD = n_alt,
        mean_nonAD = mean_ref,
        mean_lateAD = mean_alt,
        delta_lateAD_minus_nonAD = NA_real_,
        lm_beta = NA_real_,
        lm_p = NA_real_,
        mw_p = NA_real_
      ))
    }

    fit <- tryCatch(
      lm(value ~ Pathology, data = dat),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      lm_beta <- NA_real_
      lm_p <- NA_real_
    } else {
      coef_tab <- summary(fit)$coefficients
      coef_row <- grep("^Pathology", rownames(coef_tab), value = TRUE)

      if (length(coef_row) == 1) {
        lm_beta <- unname(coef_tab[coef_row, "Estimate"])
        lm_p <- unname(coef_tab[coef_row, "Pr(>|t|)"])
      } else {
        lm_beta <- NA_real_
        lm_p <- NA_real_
      }
    }

    mw_p <- tryCatch(
      wilcox.test(value ~ Pathology, data = dat, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )

    tibble(
      TF = tf,
      celltype = celltype,
      region = region,
      celltype_region = pair_name,
      data_type = data_type,
      n_nonAD = n_ref,
      n_lateAD = n_alt,
      mean_nonAD = mean_ref,
      mean_lateAD = mean_alt,
      delta_lateAD_minus_nonAD = mean_alt - mean_ref,
      lm_beta = lm_beta,
      lm_p = lm_p,
      mw_p = mw_p
    )
  }

  bind_rows(lapply(rownames(mat_sub), one_tf))
}

adjust_finite_p <- function(p) {
  out <- rep(NA_real_, length(p))
  keep <- is.finite(p)
  if (any(keep)) {
    out[keep] <- p.adjust(
      p[keep],
      method = "BH"
    )
  }
  out
}
                
make_long_results <- function(x) {
  if (nrow(x) == 0) return(tibble())

  bind_rows(
    x %>%
      transmute(
        TF, celltype, region, celltype_region, data_type,
        n_nonAD, n_lateAD,
        mean_nonAD, mean_lateAD, delta_lateAD_minus_nonAD,
        method = "Linear model",
        effect = lm_beta,
        p_value = lm_p
      ),

    x %>%
      transmute(
        TF, celltype, region, celltype_region, data_type,
        n_nonAD, n_lateAD,
        mean_nonAD, mean_lateAD, delta_lateAD_minus_nonAD,
        method = "Mann-Whitney",
        effect = delta_lateAD_minus_nonAD,
        p_value = mw_p
      )
  ) %>%
    group_by(celltype_region, data_type, method) %>%
    mutate(
      p_adj = adjust_finite_p(p_value)
    ) %>%
    ungroup()
}

format_p <- function(p) {
  if (length(p) == 0 || is.na(p)) return("NA")
  if (p < 2.2e-16) return("< 2.2e-16")
  if (p < 0.001) return(formatC(p, format = "e", digits = 2))
  formatC(p, format = "f", digits = 3)
}

# Keep only celltype-region pairs for which at least one of Expression or
# Activity has >=1 significant TF, separately for each test method.
filter_zero_zero_pairs <- function(summary_data, method_name) {
  summary_data %>%
    filter(method == method_name) %>%
    group_by(celltype_region) %>%
    filter(sum(n_significant, na.rm = TRUE) > 0) %>%
    ungroup()
}


# ------------------------------------------------------------
# 3. Discover metadata files
# ------------------------------------------------------------

metadata_files <- list.files(
  metadata_dir,
  pattern = "_metadata\\.tsv$",
  full.names = TRUE
)
metadata_files <- metadata_files[!grepl("^Vasc_", basename(metadata_files))]

if (length(metadata_files) == 0) {
  stop("No metadata files found in: ", metadata_dir)
}

message("Found ", length(metadata_files), " metadata files.")

all_results <- list()
processed_pairs <- character()
skipped_pairs <- character()

# ------------------------------------------------------------
# 4. Loop over celltype-region datasets
# ------------------------------------------------------------

for (meta_file in metadata_files) {
  pair_name <- sub("_metadata\\.tsv$", "", basename(meta_file))
  message("\nProcessing: ", pair_name)

  meta <- tryCatch(
    read_metadata(meta_file),
    error = function(e) {
      warning(conditionMessage(e))
      NULL
    }
  )

  if (is.null(meta) || nrow(meta) == 0) {
    warning("No usable nonAD/lateAD metadata for ", pair_name)
    skipped_pairs <- c(skipped_pairs, pair_name)
    next
  }

  celltypes_here <- unique(meta$celltype)
  regions_here <- unique(meta$region)

  if (length(celltypes_here) != 1 || length(regions_here) != 1) {
    warning(
      "Expected one celltype and one region in metadata for ", pair_name,
      "; found celltypes={", paste(celltypes_here, collapse = ","),
      "}, regions={", paste(regions_here, collapse = ","), "}. Skipping."
    )
    skipped_pairs <- c(skipped_pairs, pair_name)
    next
  }

  celltype <- celltypes_here[1]
  region <- regions_here[1]

  pair_dir <- file.path(tf_root)

  expression_file <- file.path(
    pair_dir,
    paste0(pair_name, "_TF_expression_pcc_cutoff.csv")
  )

  activity_file <- file.path(
    pair_dir,
    paste0(pair_name, "_TF_activity_pcc_cutoff.csv")
  )

  if (!file.exists(expression_file) || !file.exists(activity_file)) {
    warning(
      "Missing expression and/or activity file for ", pair_name,
      "\n  expression: ", expression_file,
      "\n  activity:   ", activity_file
    )
    skipped_pairs <- c(skipped_pairs, pair_name)
    next
  }

  expr00 <- read_tf_matrix(expression_file)
  act00 <- read_tf_matrix(activity_file)
  expr0 <- expr00[rowSums(abs(expr00), na.rm = TRUE) != 0, , drop = FALSE]
  act0 <- act00[rowSums(abs(act00), na.rm = TRUE) != 0, , drop = FALSE]
  common_TF <- intersect(rownames(expr0), rownames(act0))
  if (length(common_TF) == 0) {
    stop("No shared nonzero TFs for ", pair_name)
  }
  expr <- expr0[common_TF, , drop = FALSE]
  act <- act0[common_TF, , drop = FALSE]
  message("  Matched nonzero TF rows: ", length(common_TF))

  expr_res <- run_tf_tests(
    expr, meta,
    data_type = "Expression",
    celltype = celltype,
    region = region,
    pair_name = pair_name
  )

  act_res <- run_tf_tests(
    act, meta,
    data_type = "Activity",
    celltype = celltype,
    region = region,
    pair_name = pair_name
  )

  all_results[[pair_name]] <- bind_rows(expr_res, act_res)
  processed_pairs <- c(processed_pairs, pair_name)
}

if (length(unique(processed_pairs)) != 36 || length(skipped_pairs) != 0) {
  stop("Expected all 36 non-Vasc celltype-region pairs; processed ",
       length(unique(processed_pairs)), "; skipped ", length(unique(skipped_pairs)))
}

raw_results <- bind_rows(all_results)

if (nrow(raw_results) == 0) {
  stop("No celltype-region dataset was successfully processed.")
}

results_long <- make_long_results(raw_results)

# ------------------------------------------------------------
# 5. Define significance and summarize COUNTS
# ------------------------------------------------------------

if (!significance_column %in% c("p_value", "p_adj")) {
  stop("significance_column must be either 'p_value' or 'p_adj'.")
}

results_long <- results_long %>%
  mutate(significant = !is.na(.data[[significance_column]]) & .data[[significance_column]] < alpha)

summary_df <- results_long %>%
  group_by(celltype, region, celltype_region, data_type, method) %>%
  summarise(
    n_TF_total_in_matrix = n(),
    n_TF_tested = sum(!is.na(.data[[significance_column]])),
    n_significant = sum(significant, na.rm = TRUE),
    n_not_significant = n_TF_tested - n_significant,
    .groups = "drop"
  ) %>%
  mutate(
    data_type = factor(data_type, levels = c("Expression", "Activity")),
    method = factor(method, levels = c("Linear model", "Mann-Whitney"))
  )

# A plotting-filtered summary is saved separately. Filtering is performed
# within each statistical method, because a pair can be zero-zero for one
# method but not the other.
summary_plot_df <- summary_df %>%
  group_by(method, celltype_region) %>%
  filter(sum(n_significant, na.rm = TRUE) > 0) %>%
  ungroup()

# ------------------------------------------------------------
# 6. Save result tables
# ------------------------------------------------------------

write_tsv(
  raw_results,
  file.path(output_dir, "TF_expression_activity_tests_wide.tsv")
)

write_tsv(
  results_long,
  file.path(output_dir, "TF_expression_activity_tests_long.tsv")
)

write_tsv(
  summary_df,
  file.path(output_dir, "significant_TF_count_summary_all.tsv")
)

write_tsv(
  summary_plot_df,
  file.path(output_dir, "significant_TF_count_summary_nonzero_pairs.tsv")
)

# ------------------------------------------------------------
# 7. Faceted barplot: NUMBER of significant TFs
# ------------------------------------------------------------

plot_count_facets <- function(summary_data, method_name, output_prefix) {
  dat <- filter_zero_zero_pairs(summary_data, method_name) %>%
    arrange(region, celltype)

  if (nrow(dat) == 0) {
    warning("No nonzero celltype-region pairs for method: ", method_name)
    return(NULL)
  }

  pair_order <- unique(dat$celltype_region)
  dat$celltype_region <- factor(dat$celltype_region, levels = pair_order)

  n_pairs <- length(pair_order)
  ncol_facets <- min(6, max(1, ceiling(sqrt(n_pairs))))
  nrow_facets <- ceiling(n_pairs / ncol_facets)

  y_max <- max(dat$n_significant, na.rm = TRUE)
  y_upper <- max(1, y_max * 1.18 + 0.5)

  p <- ggplot(dat, aes(x = data_type, y = n_significant, fill = data_type)) +
    geom_col(width = 0.68, color = "black", linewidth = 0.25) +
    geom_text(
      aes(label = n_significant),
      vjust = -0.35,
      size = 3.0
    ) +
    facet_wrap(~ celltype_region, ncol = ncol_facets) +
    scale_y_continuous(
      limits = c(0, y_upper),
      breaks = scales::pretty_breaks(n = 5),
      expand = expansion(mult = c(0, 0.01))
    ) +
    labs(
      x = NULL,
      y = "Number of significant TFs",
      title = paste0(method_name, ": TF expression vs TF activity"),
      subtitle = paste0(
        "nonAD vs lateAD; ", significance_column, " < ", alpha,
        "; zero-zero celltype-region pairs removed"
      )
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 9),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 25, hjust = 1),
      plot.title = element_text(face = "bold")
    )

  fig_width <- max(8, ncol_facets * 2.15)
  fig_height <- max(5, nrow_facets * 2.25 + 1.2)

  ggsave(
    file.path(output_dir, paste0(output_prefix, ".pdf")),
    p,
    width = fig_width,
    height = fig_height,
    units = "in"
  )

  ggsave(
    file.path(output_dir, paste0(output_prefix, ".png")),
    p,
    width = fig_width,
    height = fig_height,
    units = "in",
    dpi = 300
  )

  p
}

p_linear <- plot_count_facets(
  summary_df,
  method_name = "Linear model",
  output_prefix = "significant_TF_count_linear_model_all_celltype_region"
)

p_mw <- plot_count_facets(
  summary_df,
  method_name = "Mann-Whitney",
  output_prefix = "significant_TF_count_mann_whitney_all_celltype_region"
)

# ------------------------------------------------------------
# 8. Overall boxplot, UNPAIRED display
#    One distribution = Expression significant-TF counts
#    One distribution = Activity significant-TF counts
# ------------------------------------------------------------

plot_overall_boxplot <- function(summary_data, method_name, output_prefix) {
  dat <- filter_zero_zero_pairs(summary_data, method_name)

  if (nrow(dat) == 0) {
    warning("No nonzero celltype-region pairs for method: ", method_name)
    return(NULL)
  }

  # Use the same matched celltype-region set for Expression and Activity.
  complete_pairs <- dat %>%
    count(celltype_region, name = "n_rows") %>%
    filter(n_rows == 2) %>%
    pull(celltype_region)

  dat <- dat %>%
    filter(celltype_region %in% complete_pairs)

  if (nrow(dat) == 0) {
    warning("No complete Expression/Activity pairs for method: ", method_name)
    return(NULL)
  }

  y_max <- max(dat$n_significant, na.rm = TRUE)
  y_annot <- max(1, y_max * 1.16 + 0.5)
  y_upper <- max(1, y_max * 1.30 + 1)

  p <- ggplot(dat, aes(x = data_type, y = n_significant, fill = data_type)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.65) +
    geom_jitter(width = 0.08, height = 0, size = 2.0, alpha = 0.8) +
    annotate(
      "segment",
      x = 1, xend = 2,
      y = y_annot, yend = y_annot,
      linewidth = 0.45
    ) +
    scale_y_continuous(
      limits = c(0, y_upper),
      breaks = scales::pretty_breaks(n = 6),
      expand = expansion(mult = c(0, 0.01))
    ) +
    labs(
      x = NULL,
      y = "Number of significant TFs",
      title = paste0(method_name, ": significant TF counts across celltype-region pairs"),
      subtitle = paste0(
        "Each point = one celltype-region; zero-zero pairs removed; n = ",
        n_distinct(dat$celltype_region), " pairs"
      )
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )

  ggsave(
    file.path(output_dir, paste0(output_prefix, ".pdf")),
    p,
    width = 5.8,
    height = 5.4,
    units = "in"
  )

  ggsave(
    file.path(output_dir, paste0(output_prefix, ".png")),
    p,
    width = 5.8,
    height = 5.4,
    units = "in",
    dpi = 300
  )

  list(plot = p, data = dat)
}

# ------------------------------------------------------------
# 9. Overall PAIRED boxplot
#    Same boxplots, with matched celltype-region values connected.
# ------------------------------------------------------------

plot_overall_paired_boxplot <- function(summary_data, method_name, output_prefix) {
  dat <- filter_zero_zero_pairs(summary_data, method_name)

  if (nrow(dat) == 0) {
    warning("No nonzero celltype-region pairs for method: ", method_name)
    return(NULL)
  }

  # Keep only truly matched pairs that have both Expression and Activity.
  complete_pairs <- dat %>%
    count(celltype_region, name = "n_rows") %>%
    filter(n_rows == 2) %>%
    pull(celltype_region)

  dat <- dat %>%
    filter(celltype_region %in% complete_pairs)

  wide <- dat %>%
    select(celltype_region, data_type, n_significant) %>%
    pivot_wider(names_from = data_type, values_from = n_significant)

  paired_p <- if (
    nrow(wide) >= 2 &&
    all(c("Expression", "Activity") %in% colnames(wide))
  ) {
    tryCatch(
      wilcox.test(
        wide$Expression,
        wide$Activity,
        paired = TRUE,
        exact = FALSE
      )$p.value,
      error = function(e) NA_real_
    )
  } else {
    NA_real_
  }

  paired_label <- paste0("Paired Wilcoxon p = ", format_p(paired_p))

  y_max <- max(dat$n_significant, na.rm = TRUE)
  y_annot <- max(1, y_max * 1.17 + 0.5)
  y_upper <- max(1, y_max * 1.38 + 1)

  p <- ggplot(dat, aes(x = data_type, y = n_significant)) +
    geom_line(
      aes(group = celltype_region),
      color = "grey65",
      linewidth = 0.45,
      alpha = 0.65
    ) +
    geom_boxplot(
      aes(fill = data_type),
      width = 0.55,
      outlier.shape = NA,
      alpha = 0.65
    ) +
    geom_point(
      aes(fill = data_type),
      shape = 21,
      size = 2.2,
      alpha = 0.85
    ) +
    annotate(
      "segment",
      x = 1, xend = 2,
      y = y_annot, yend = y_annot,
      linewidth = 0.45
    ) +
    annotate(
      "text",
      x = 1.5,
      y = y_annot + 0.015 * y_upper,
      label = paired_label,
      size = 3.6
    ) +
    scale_y_continuous(
      limits = c(0, y_upper),
      breaks = scales::pretty_breaks(n = 6),
      expand = expansion(mult = c(0, 0.01))
    ) +
    labs(
      x = NULL,
      y = "Number of significant TFs",
      title = paste0(method_name, ": paired celltype-region comparison"),
      subtitle = paste0(
        "Lines connect the same celltype-region; zero-zero pairs removed; n = ",
        nrow(wide), " matched pairs"
      )
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )

  ggsave(
    file.path(output_dir, paste0(output_prefix, ".pdf")),
    p,
    width = 5.8,
    height = 5.4,
    units = "in"
  )

  ggsave(
    file.path(output_dir, paste0(output_prefix, ".png")),
    p,
    width = 5.8,
    height = 5.4,
    units = "in",
    dpi = 300
  )

  list(
    plot = p,
    paired_wilcox_p = paired_p,
    paired_wide = wide,
    data = dat
  )
}

# Make both overall figures for each test version.
linear_overall <- plot_overall_boxplot(
  summary_df,
  method_name = "Linear model",
  output_prefix = "overall_boxplot_linear_model_significant_TF_count"
)

mw_overall <- plot_overall_boxplot(
  summary_df,
  method_name = "Mann-Whitney",
  output_prefix = "overall_boxplot_mann_whitney_significant_TF_count"
)

linear_paired <- plot_overall_paired_boxplot(
  summary_df,
  method_name = "Linear model",
  output_prefix = "overall_paired_boxplot_linear_model_significant_TF_count"
)

mw_paired <- plot_overall_paired_boxplot(
  summary_df,
  method_name = "Mann-Whitney",
  output_prefix = "overall_paired_boxplot_mann_whitney_significant_TF_count"
)

# ------------------------------------------------------------
# 10. Save overall statistical-test summaries
# ------------------------------------------------------------

extract_test_row <- function(method_name, overall_obj, paired_obj) {
  if (is.null(overall_obj) || is.null(paired_obj)) return(tibble())

  tibble(
    method = method_name,
    n_celltype_region_nonzero = n_distinct(overall_obj$data$celltype_region),
    n_matched_pairs = nrow(paired_obj$paired_wide),
    paired_wilcoxon_p = paired_obj$paired_wilcox_p
  )
}

overall_test_summary <- bind_rows(
  extract_test_row("Linear model", linear_overall, linear_paired),
  extract_test_row("Mann-Whitney", mw_overall, mw_paired)
)

write_tsv(
  overall_test_summary,
  file.path(output_dir, "overall_expression_vs_activity_statistical_tests.tsv")
)

# Save paired count values explicitly for easy inspection.
if (!is.null(linear_paired)) {
  write_tsv(
    linear_paired$paired_wide,
    file.path(output_dir, "paired_significant_TF_counts_linear_model.tsv")
  )
}

if (!is.null(mw_paired)) {
  write_tsv(
    mw_paired$paired_wide,
    file.path(output_dir, "paired_significant_TF_counts_mann_whitney.tsv")
  )
}

# ------------------------------------------------------------
# 11. Run summary
# ------------------------------------------------------------

zero_zero_removed <- summary_df %>%
  group_by(method, celltype_region) %>%
  summarise(
    total_sig = sum(n_significant, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(total_sig == 0)

writeLines(
  c(
    paste0("Significance column: ", significance_column),
    paste0("Alpha: ", alpha),
    paste0("Minimum n per pathology group: ", min_n_per_group),
    paste0("Processed pairs (", length(unique(processed_pairs)), "): ", paste(unique(processed_pairs), collapse = ", ")),
    paste0("Skipped pairs (", length(unique(skipped_pairs)), "): ", paste(unique(skipped_pairs), collapse = ", ")),
    paste0(
      "Linear-model zero-zero pairs removed from plots: ",
      sum(zero_zero_removed$method == "Linear model")
    ),
    paste0(
      "Mann-Whitney zero-zero pairs removed from plots: ",
      sum(zero_zero_removed$method == "Mann-Whitney")
    )
  ),
  file.path(output_dir, "run_summary.txt")
)

cat("\n============================================================\n")
cat("Finished TF expression vs activity COUNT comparison\n")
cat("============================================================\n")
cat("Processed celltype-region pairs:", length(unique(processed_pairs)), "\n")
cat("Significance:", significance_column, "<", alpha, "\n")
cat("Output directory:\n", output_dir, "\n\n")

cat("Significant-TF count summary:\n")
print(
  summary_df %>%
    select(
      celltype_region, data_type, method,
      n_TF_total_in_matrix, n_TF_tested, n_significant
    ) %>%
    arrange(method, celltype_region, data_type),
  n = Inf
)

cat("\nOverall tests:\n")
print(overall_test_summary, n = Inf)
