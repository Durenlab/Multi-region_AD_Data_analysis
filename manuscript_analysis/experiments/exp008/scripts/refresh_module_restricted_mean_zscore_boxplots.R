#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  stop("Usage: refresh_module_restricted_mean_zscore_boxplots.R <EXP008-results-root> [<results-root> ...]")
}

pathology_levels <- c("nonAD", "earlyAD", "lateAD")
pathology_colors <- c(nonAD = "#66AABB", earlyAD = "#EEAA55", lateAD = "#CC6677")
result_columns <- c(
  "pair", "celltype", "region", "TF", "mode", "TF_module",
  "n_selected_target_genes", "n_expression_matched_target_genes",
  "n_nonAD", "n_earlyAD", "n_lateAD", "status", "test_alternative",
  "beta_pathology", "ci_low", "ci_high", "p_value", "summary_scale",
  "fdr_within_pair_mode"
)

fit_greater <- function(dat) {
  y <- suppressWarnings(as.numeric(dat$mean_target_expression))
  pathology <- factor(dat$Pathology, levels = pathology_levels)
  keep <- is.finite(y) & !is.na(pathology)
  n_group <- table(factor(pathology[keep], levels = pathology_levels))
  out <- list(
    n_nonAD = as.integer(n_group[[1]]),
    n_earlyAD = as.integer(n_group[[2]]),
    n_lateAD = as.integer(n_group[[3]]),
    status = if (grepl("^Vasc_", dat$pair[[1]])) "descriptive_only_vascular" else "not_tested_low_group_n",
    beta_pathology = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_
  )
  if (!grepl("^Vasc_", dat$pair[[1]]) && all(n_group >= 5L)) {
    model_data <- data.frame(
      y = y[keep],
      pathology_code = as.integer(pathology[keep]) - 1L
    )
    fit <- tryCatch(lm(y ~ pathology_code, data = model_data), error = function(e) NULL)
    if (!is.null(fit)) {
      coef <- summary(fit)$coefficients
      if ("pathology_code" %in% rownames(coef)) {
        beta <- unname(coef["pathology_code", "Estimate"])
        t_stat <- unname(coef["pathology_code", "t value"])
        if (is.finite(beta) && is.finite(t_stat)) {
          interval <- tryCatch(confint(fit)["pathology_code", ], error = function(e) c(NA_real_, NA_real_))
          out$beta_pathology <- beta
          out$ci_low <- unname(interval[[1]])
          out$ci_high <- unname(interval[[2]])
          out$p_value <- pt(t_stat, df = df.residual(fit), lower.tail = FALSE)
          out$status <- "tested"
        } else {
          out$status <- "not_estimable"
        }
      } else {
        out$status <- "not_estimable"
      }
    } else {
      out$status <- "not_estimable"
    }
  }
  out
}

format_p <- function(x) {
  if (!is.finite(x)) return("NA")
  if (x < 0.001) sprintf("%.2e", x) else sprintf("%.3f", x)
}

for (root_arg in args) {
  root <- normalizePath(root_arg, mustWork = TRUE)
  mode_root <- file.path(root, "module_restricted")
  score_root <- file.path(mode_root, "mean_target_expression_by_individual")
  score_files <- sort(list.files(
    score_root,
    pattern = "_mean_target_expression_by_individual\\.tsv\\.gz$",
    recursive = TRUE,
    full.names = TRUE
  ))
  if (!length(score_files)) stop("No donor-level score files under ", score_root)

  refreshed_pairs <- character()
  refreshed_rows <- list()
  refreshed_plots <- 0L
  for (score_file in score_files) {
    scores <- read_tsv(score_file, show_col_types = FALSE,
                       col_types = cols(.default = col_character()))
    required <- c("Individual", "Pathology", "pair", "celltype", "region", "TF", "mode",
                  "mean_target_expression", "n_selected_target_genes",
                  "n_expression_matched_target_genes")
    missing <- setdiff(required, names(scores))
    if (length(missing)) stop("Missing score columns in ", score_file, ": ", paste(missing, collapse = ", "))
    pair <- unique(scores$pair)
    if (length(pair) != 1L) stop("Expected one pair per score file: ", score_file)
    scores <- scores |> filter(mode == "module_restricted")
    if (!nrow(scores)) next
    refreshed_pairs <- c(refreshed_pairs, pair)
    result_file <- file.path(mode_root, "regression_results", pair,
                             paste0(pair, "_TF_level_regression.tsv"))
    old <- if (file.exists(result_file)) {
      read_tsv(result_file, show_col_types = FALSE,
               col_types = cols(.default = col_character()))
    } else tibble()
    for (column in result_columns) if (!column %in% names(old)) old[[column]] <- NA_character_
    old_module <- old |> filter(mode == "module_restricted")
    tf_groups <- split(scores, scores$TF)
    pair_rows <- lapply(names(tf_groups), function(tf) {
      dat <- tf_groups[[tf]]
      fit <- fit_greater(dat)
      module_row <- old_module[match(tf, old_module$TF), , drop = FALSE]
      module_value <- if (nrow(module_row) && !is.na(module_row$TF_module[[1]])) {
        module_row$TF_module[[1]]
      } else NA_character_
      data.frame(
        pair = pair,
        celltype = dat$celltype[[1]],
        region = dat$region[[1]],
        TF = tf,
        mode = "module_restricted",
        TF_module = module_value,
        n_selected_target_genes = dat$n_selected_target_genes[[1]],
        n_expression_matched_target_genes = dat$n_expression_matched_target_genes[[1]],
        n_nonAD = fit$n_nonAD,
        n_earlyAD = fit$n_earlyAD,
        n_lateAD = fit$n_lateAD,
        status = fit$status,
        test_alternative = "greater",
        beta_pathology = fit$beta_pathology,
        ci_low = fit$ci_low,
        ci_high = fit$ci_high,
        p_value = fit$p_value,
        summary_scale = "mean_per_gene_zscore_across_individuals",
        fdr_within_pair_mode = NA_real_,
        stringsAsFactors = FALSE
      )
    })
    pair_results <- bind_rows(pair_rows)
    valid <- which(is.finite(pair_results$p_value))
    if (length(valid)) pair_results$fdr_within_pair_mode[valid] <- p.adjust(pair_results$p_value[valid], method = "BH")

    for (tf in names(tf_groups)) {
      dat <- tf_groups[[tf]] |>
        mutate(Pathology = factor(Pathology, levels = pathology_levels),
               mean_target_expression = suppressWarnings(as.numeric(mean_target_expression)))
      test <- pair_results[pair_results$TF == tf, , drop = FALSE]
      subtitle <- if (is.finite(test$p_value[[1]])) {
        sprintf("LM greater-than p = %s; BH FDR = %s | n = %d / %d / %d",
                format_p(test$p_value[[1]]), format_p(test$fdr_within_pair_mode[[1]]),
                test$n_nonAD[[1]], test$n_earlyAD[[1]], test$n_lateAD[[1]])
      } else {
        sprintf("%s | n = %d / %d / %d", test$status[[1]],
                test$n_nonAD[[1]], test$n_earlyAD[[1]], test$n_lateAD[[1]])
      }
      plot <- ggplot(dat, aes(x = Pathology, y = mean_target_expression, fill = Pathology)) +
        geom_boxplot(outlier.shape = NA, width = 0.7) +
        geom_jitter(width = 0.12, size = 0.7, alpha = 0.55) +
        scale_x_discrete(limits = pathology_levels, drop = FALSE) +
        scale_fill_manual(values = pathology_colors, drop = FALSE) +
        labs(title = paste(tf, "target genes in", pair), subtitle = subtitle,
             x = "Pathology group", y = "Mean per-gene z-score across selected targets") +
        theme_classic(base_size = 9) +
        theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1))
      output <- file.path(mode_root, "plots", pair,
                          paste0(tf, "_module_restricted_target_expression_boxplot.pdf"))
      dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
      ggsave(output, plot, width = 3.5, height = 3.5, units = "in", limitsize = FALSE)
      refreshed_plots <- refreshed_plots + 1L
    }

    old_other <- old |> filter(mode != "module_restricted")
    old_module <- old_module |> filter(!TF %in% pair_results$TF)
    if (nrow(old_module)) {
      old_module$test_alternative <- "greater"
      old_module$beta_pathology <- NA_character_
      old_module$ci_low <- NA_character_
      old_module$ci_high <- NA_character_
      old_module$p_value <- NA_character_
      old_module$fdr_within_pair_mode <- NA_character_
    }
    replacement_module <- bind_rows(
      old_module |> select(all_of(result_columns)),
      pair_results |> mutate(across(everything(), as.character)) |> select(all_of(result_columns))
    )
    replacement <- bind_rows(old_other |> select(all_of(result_columns)), replacement_module)
    write_tsv(replacement, result_file, na = "NA")
    refreshed_rows[[length(refreshed_rows) + 1L]] <- replacement_module
    cat(pair, ": refreshed ", length(tf_groups), " module-restricted TF plots\n", sep = "")
  }

  summary_file <- file.path(root, "summary", "TF_level_regression_results.tsv")
  if (file.exists(summary_file) && length(refreshed_pairs)) {
    summary <- read_tsv(summary_file, show_col_types = FALSE,
                         col_types = cols(.default = col_character()))
    for (column in result_columns) if (!column %in% names(summary)) summary[[column]] <- NA_character_
    summary <- summary |>
      filter(mode != "module_restricted" | !pair %in% refreshed_pairs)
    replacements <- bind_rows(refreshed_rows) |> mutate(across(everything(), as.character))
    write_tsv(bind_rows(summary, replacements), summary_file, na = "NA")
  }
  cat("Updated ", refreshed_plots, " boxplots under ", mode_root,
      " using one-sided greater-than regressions; pairs: ", length(unique(refreshed_pairs)), "\n", sep = "")
}
