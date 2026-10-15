#!/usr/bin/env Rscript

# Donor-level pathology tests, sensitivity summaries, main figure, and report.
# Requires the six regional outputs from run_donor_cellchat_by_region.R.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(rmarkdown)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: analyze_donor_cellchat.R <RESULTS_ROOT>")
}
results_root <- normalizePath(args[[1]], mustWork = TRUE)
table_dir <- file.path(results_root, "tables")
intermediate_dir <- file.path(results_root, "intermediate")
figure_dir <- file.path(results_root, "figures")
report_dir <- file.path(results_root, "report")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

regions <- c("AG", "EC", "HC", "MTC", "PFC", "TH")
pathology_order <- c("nonAD", "earlyAD", "lateAD")
pathology_code <- c(nonAD = 0, earlyAD = 1, lateAD = 2)
celltype_order <- c("Exc", "Inh", "Ast", "Oli", "OPC", "Mic_Immune", "Vasc")
min_cells <- 20L
bootstrap_reps <- 1000L
set.seed(20260923)

read_tsv <- function(path) {
  if (is.character(path) && grepl("\\.gz$", path)) {
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
    return(read.delim(con, sep = "\t", check.names = FALSE,
                       stringsAsFactors = FALSE, quote = ""))
  }
  read.delim(path, sep = "\t", check.names = FALSE,
             stringsAsFactors = FALSE, quote = "")
}
adjust_finite_bh <- function(p) {
  out <- rep(NA_real_, length(p))
  keep <- is.finite(p)
  if (any(keep)) out[keep] <- p.adjust(p[keep], method = "BH")
  out
}

trend_results <- list()
pairwise_results <- list()
coverage_all <- list()
candidate_display <- list()
scores_display <- list()
bootstrap_display <- list()

comparisons <- list(
  earlyAD_minus_nonAD = c("nonAD", "earlyAD"),
  lateAD_minus_nonAD = c("nonAD", "lateAD"),
  lateAD_minus_earlyAD = c("earlyAD", "lateAD")
)

for (region in regions) {
  score_path <- file.path(intermediate_dir,
                          paste0(region, "_donor_nonzero_probabilities.tsv.gz"))
  coverage_path <- file.path(table_dir,
                             paste0(region, "_donor_celltype_coverage.tsv"))
  candidate_path <- file.path(table_dir,
                              paste0(region, "_candidate_interactions.tsv"))
  missing <- c(score_path, coverage_path, candidate_path)[
    !file.exists(c(score_path, coverage_path, candidate_path))
  ]
  if (length(missing)) stop("Missing required file(s): ", paste(missing, collapse = ", "))

  cov <- read_tsv(coverage_path)
  coverage_all[[region]] <- cov
  cand <- read_tsv(candidate_path)
  score <- read_tsv(score_path)
  cov$n_cells <- as.integer(cov$n_cells)
  cov$pathology <- as.character(cov$pathology)
  cov$donor_id <- as.character(cov$donor_id)

  cov <- cov[cov$donor_network_status == "completed", , drop = FALSE]
  cov <- cov[!duplicated(cov[c("region", "pathology", "donor_id", "celltype")]), , drop = FALSE]
  donor_meta <- unique(cov[c("region", "pathology", "donor_id")])
  donor_meta$stage <- unname(pathology_code[donor_meta$pathology])
  if (anyNA(donor_meta$stage)) stop("Unexpected pathology in donor coverage for ", region)

  # A source-target pair is testable only when both donor cell groups exceeded
  # CellChat's min.cells threshold. Other pairs remain missing, not zero.
  count_matrix <- matrix(0L, nrow = nrow(donor_meta), ncol = length(celltype_order),
                         dimnames = list(donor_meta$donor_id, celltype_order))
  ri <- match(cov$donor_id, donor_meta$donor_id)
  ci <- match(cov$celltype, celltype_order)
  count_matrix[cbind(ri, ci)] <- cov$n_cells

  feature_grid <- expand.grid(
    interaction_id = cand$interaction_id,
    target = celltype_order,
    source = celltype_order,
    stringsAsFactors = FALSE
  )
  feature_grid$key <- paste(feature_grid$source, feature_grid$target,
                            feature_grid$interaction_id, sep = "\034")
  score$key <- paste(score$source, score$target, score$interaction_id, sep = "\034")
  score_by_feature <- split(score[c("donor_id", "communication_probability")], score$key)

  tr <- vector("list", nrow(feature_grid))
  pw <- vector("list", nrow(feature_grid) * length(comparisons))
  pw_i <- 0L
  for (i in seq_len(nrow(feature_grid))) {
    f <- feature_grid[i, , drop = FALSE]
    i_source <- match(f$source, celltype_order)
    i_target <- match(f$target, celltype_order)
    valid <- count_matrix[, i_source] > min_cells & count_matrix[, i_target] > min_cells
    dm <- donor_meta[valid, , drop = FALSE]
    if (nrow(dm) == 0L) next

    y <- rep(0, nrow(dm))
    sf <- score_by_feature[[f$key]]
    if (!is.null(sf) && nrow(sf)) {
      hit <- match(dm$donor_id, sf$donor_id)
      present <- !is.na(hit)
      y[present] <- as.numeric(sf$communication_probability[hit[present]])
    }
    has_variance <- length(unique(y)) > 1L
    stage_counts <- table(factor(dm$pathology, levels = pathology_order))
    trend_p <- NA_real_
    trend_rho <- 0
    trend_status <- "insufficient_group_coverage"
    if (all(stage_counts >= 3L)) {
      if (has_variance) {
        trend_test <- tryCatch(
          suppressWarnings(cor.test(dm$stage, y, method = "spearman", exact = FALSE)),
          error = function(e) NULL
        )
        if (!is.null(trend_test) && is.finite(trend_test$p.value)) {
          trend_p <- trend_test$p.value
          trend_rho <- unname(trend_test$estimate)
          trend_status <- "tested"
        } else {
          trend_status <- "not_estimable"
        }
      } else {
        trend_p <- 1
        trend_status <- "no_score_variation"
      }
    }
    tr[[i]] <- data.frame(
      region = region, source = f$source, target = f$target,
      interaction_id = f$interaction_id,
      n_nonAD = as.integer(stage_counts[["nonAD"]]),
      n_earlyAD = as.integer(stage_counts[["earlyAD"]]),
      n_lateAD = as.integer(stage_counts[["lateAD"]]),
      spearman_rho = trend_rho, p_value = trend_p,
      status = trend_status, stringsAsFactors = FALSE
    )

    for (comparison in names(comparisons)) {
      groups <- comparisons[[comparison]]
      y1 <- y[dm$pathology == groups[[1]]]
      y2 <- y[dm$pathology == groups[[2]]]
      p <- NA_real_
      status <- "insufficient_donors"
      delta <- NA_real_
      median_difference <- NA_real_
      if (length(y1) >= 3L && length(y2) >= 3L) {
        wt <- tryCatch(
          suppressWarnings(wilcox.test(y2, y1, exact = FALSE, correct = TRUE)),
          error = function(e) NULL
        )
        if (!is.null(wt) && is.finite(wt$p.value)) {
          p <- wt$p.value
          status <- if (length(unique(c(y1, y2))) > 1L) "tested" else "no_score_variation"
        }
        cmp <- outer(y2, y1, "-")
        delta <- (sum(cmp > 0) - sum(cmp < 0)) / length(cmp)
        median_difference <- median(y2) - median(y1)
      }
      pw_i <- pw_i + 1L
      pw[[pw_i]] <- data.frame(
        region = region, comparison = comparison,
        condition_1 = groups[[1]], condition_2 = groups[[2]],
        source = f$source, target = f$target,
        interaction_id = f$interaction_id,
        n_condition_1 = length(y1), n_condition_2 = length(y2),
        median_difference = median_difference,
        rank_biserial_delta = delta,
        p_value = p, status = status, stringsAsFactors = FALSE
      )
    }
  }

  tr <- do.call(rbind, tr[!vapply(tr, is.null, logical(1))])
  pw <- do.call(rbind, pw[seq_len(pw_i)])
  trend_results[[region]] <- tr
  pairwise_results[[region]] <- pw

  # Keep candidate metadata for labels and interpretation in the report.
  candidate_display[[region]] <- cand
  rm(score_by_feature, feature_grid, score, cov, donor_meta, count_matrix)
  gc(verbose = FALSE)
  cat("Tested feature families for", region, "\n")
}

trend <- do.call(rbind, trend_results)
pairwise <- do.call(rbind, pairwise_results)
trend$q_value <- ave(trend$p_value, trend$region, FUN = adjust_finite_bh)
pairwise$q_value <- ave(
  pairwise$p_value,
  interaction(pairwise$region, pairwise$comparison, drop = TRUE),
  FUN = adjust_finite_bh
)
write.table(trend, file.path(table_dir, "donor_pathology_progression_tests.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(pairwise, file.path(table_dir, "donor_pathology_pairwise_tests.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

network_status <- do.call(rbind, lapply(coverage_all, function(x) {
  unique(x[c("region", "pathology", "donor_id", "donor_network_status", "failure_message")])
}))
network_status$n_donors <- 1L
network_status_summary <- aggregate(
  n_donors ~ region + pathology + donor_network_status,
  data = network_status, FUN = sum
)
write.table(network_status, file.path(table_dir, "donor_network_status.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(network_status_summary,
            file.path(table_dir, "donor_network_status_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Select at most two best-ranked progression features per region. If none pass
# BH, the figure labels these as exploratory rankings and shows their q-values.
selected <- do.call(rbind, lapply(regions, function(r) {
  x <- trend[trend$region == r & trend$status == "tested", , drop = FALSE]
  if (!nrow(x)) return(NULL)
  x <- x[order(x$q_value, x$p_value, -abs(x$spearman_rho)), , drop = FALSE]
  head(x, 2L)
}))
if (is.null(selected) || nrow(selected) == 0L) stop("No donor-level trend features could be displayed")

selected$feature_key <- paste(selected$source, selected$target,
                              selected$interaction_id, sep = "\034")
selected$feature_label <- paste0(
  selected$region, ": ", selected$source, " → ", selected$target,
  " | ", selected$interaction_id, " (q=", formatC(selected$q_value, digits = 2, format = "g"), ")"
)
write.table(selected, file.path(table_dir, "main_figure_selected_features.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Donor-level data and bootstrap uncertainty for selected features only.
for (j in seq_len(nrow(selected))) {
  r <- selected$region[[j]]
  score_path <- file.path(intermediate_dir,
                          paste0(r, "_donor_nonzero_probabilities.tsv.gz"))
  coverage_path <- file.path(table_dir,
                             paste0(r, "_donor_celltype_coverage.tsv"))
  score <- read_tsv(score_path)
  cov <- read_tsv(coverage_path)
  f <- selected[j, , drop = FALSE]
  valid <- cov$n_cells[cov$celltype == f$source] > min_cells
  donors_src <- cov$donor_id[cov$celltype == f$source][valid]
  n_src <- cov$n_cells[cov$celltype == f$source][valid]
  valid_t <- cov$n_cells[cov$celltype == f$target] > min_cells
  donors_tgt <- cov$donor_id[cov$celltype == f$target][valid_t]
  donors <- intersect(donors_src, donors_tgt)
  meta <- unique(cov[cov$donor_id %in% donors,
                     c("donor_id", "pathology"), drop = FALSE])
  feature <- score[score$source == f$source & score$target == f$target &
                     score$interaction_id == f$interaction_id, , drop = FALSE]
  value <- rep(0, nrow(meta))
  ii <- match(meta$donor_id, feature$donor_id)
  hit <- !is.na(ii)
  value[hit] <- feature$communication_probability[ii[hit]]
  one <- data.frame(region = r, donor_id = meta$donor_id,
                    pathology = factor(meta$pathology, levels = pathology_order),
                    feature_label = f$feature_label, score = value)
  scores_display[[length(scores_display) + 1L]] <- one

  for (pathology in pathology_order) {
    v <- one$score[one$pathology == pathology]
    if (length(v) < 3L) next
    b <- replicate(bootstrap_reps, median(sample(v, length(v), replace = TRUE)))
    bootstrap_display[[length(bootstrap_display) + 1L]] <- data.frame(
      region = r, feature_label = f$feature_label, pathology = pathology,
      n_donors = length(v), median_score = median(v),
      median_ci_low = unname(quantile(b, 0.025)),
      median_ci_high = unname(quantile(b, 0.975)),
      stringsAsFactors = FALSE
    )
  }
  rm(score, cov)
}
score_plot_data <- do.call(rbind, scores_display)
bootstrap_summary <- do.call(rbind, bootstrap_display)
bootstrap_summary$pathology <- factor(bootstrap_summary$pathology,
                                      levels = pathology_order)
write.table(score_plot_data, file.path(table_dir, "main_figure_donor_scores.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(bootstrap_summary, file.path(table_dir, "main_figure_median_bootstrap_ci.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Donor-resampling sensitivity for selected features. Resample each group at its
# observed size; this does not replace the primary donor-level rank tests.
stability <- list()
for (j in seq_len(nrow(selected))) {
  f <- selected[j, , drop = FALSE]
  x <- score_plot_data[score_plot_data$feature_label == f$feature_label, , drop = FALSE]
  vals <- lapply(pathology_order, function(g) x$score[x$pathology == g])
  names(vals) <- pathology_order
  observed <- f$spearman_rho
  boot_rho <- replicate(bootstrap_reps, {
    yy <- unlist(lapply(pathology_order, function(g) {
      v <- vals[[g]]
      if (!length(v)) return(numeric())
      sample(v, length(v), replace = TRUE)
    }), use.names = FALSE)
    xx <- unlist(lapply(pathology_order, function(g) {
      rep(unname(pathology_code[g]), length(vals[[g]]))
    }), use.names = FALSE)
    if (length(unique(yy)) < 2L || length(unique(xx)) < 2L) return(NA_real_)
    suppressWarnings(cor(xx, yy, method = "spearman"))
  })
  stability[[j]] <- data.frame(
    region = f$region, feature_label = f$feature_label,
    observed_spearman_rho = observed,
    bootstrap_rho_ci_low = unname(quantile(boot_rho, 0.025, na.rm = TRUE)),
    bootstrap_rho_ci_high = unname(quantile(boot_rho, 0.975, na.rm = TRUE)),
    fraction_same_direction = mean(sign(boot_rho[is.finite(boot_rho)]) == sign(observed)),
    bootstrap_replicates = bootstrap_reps,
    stringsAsFactors = FALSE
  )
}
stability <- do.call(rbind, stability)
write.table(stability, file.path(table_dir, "main_figure_donor_resampling_stability.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Coverage panel uses the original donor/celltype audit, including donors with
# no cells in a cell type (which are not included as zero-score observations).
audit <- read_tsv(file.path(table_dir, "ALL_regions_pathology_donor_celltype_counts.tsv"))
audit$n_cells_donor_celltype <- as.integer(audit$n_cells_donor_celltype)
coverage_plot <- aggregate(
  donor_id ~ region + Pathology + celltype,
  data = audit[audit$n_cells_donor_celltype > min_cells, ],
  FUN = function(z) length(unique(z))
)
names(coverage_plot)[4] <- "n_donors_gt20"
coverage_plot$region <- factor(coverage_plot$region, levels = regions)
coverage_plot$Pathology <- factor(coverage_plot$Pathology, levels = pathology_order)
coverage_plot$celltype <- factor(coverage_plot$celltype, levels = celltype_order)
p_coverage <- ggplot(coverage_plot, aes(x = pathology, y = celltype, fill = n_donors_gt20)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = n_donors_gt20), size = 2.7) +
  facet_wrap(~region, nrow = 1) +
  scale_fill_viridis_c(name = "Donors\nwith >20 cells") +
  labs(title = "A. Donor and cell-type coverage", x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1),
        strip.text = element_text(face = "bold"))

p_scores <- ggplot(score_plot_data,
                   aes(x = pathology, y = score, color = pathology)) +
  geom_point(position = position_jitter(width = 0.11, height = 0),
             alpha = 0.72, size = 1.6) +
  geom_pointrange(data = bootstrap_summary,
                  aes(x = pathology, y = median_score,
                      ymin = median_ci_low, ymax = median_ci_high),
                  inherit.aes = FALSE, color = "black", linewidth = 0.35) +
  facet_wrap(~feature_label, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c(nonAD = "#66AABB", earlyAD = "#EEAA55", lateAD = "#CC6677")) +
  labs(title = "B. Donor-level CellChat probabilities",
       subtitle = "Dots are donors; bars show within-group bootstrap 95% CIs for the median",
       x = NULL, y = "Communication probability") +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 35, hjust = 1),
        strip.text = element_text(size = 7.5))

p_trend <- ggplot(selected,
                  aes(x = spearman_rho, y = reorder(feature_label, spearman_rho),
                      color = q_value < 0.05)) +
  geom_vline(xintercept = 0, color = "grey65", linewidth = 0.35) +
  geom_point(size = 2.4) +
  scale_color_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "#999999"),
                     labels = c(`TRUE` = "BH q < 0.05", `FALSE` = "Not FDR-significant"),
                     name = NULL) +
  labs(title = "C. Ordered pathology trend",
       subtitle = "At most two top-ranked features per region",
       x = "Spearman rho (nonAD=0, earlyAD=1, lateAD=2)", y = NULL) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 7))

main_figure <- p_coverage / p_scores / p_trend +
  plot_layout(heights = c(1.2, 2.1, 1.2)) +
  plot_annotation(
    title = "Donor-aware CellChat pathology analysis",
    subtitle = paste0("All eligible donors retained; within-group donor counts vary. ",
                      "Pooled condition networks are descriptive.")
  )
ggsave(file.path(figure_dir, "cellchat_pathology_main_figure.pdf"),
       main_figure, width = 16, height = 14, limitsize = FALSE)
ggsave(file.path(figure_dir, "cellchat_pathology_main_figure.png"),
       main_figure, width = 16, height = 14, dpi = 300, limitsize = FALSE)

summary_lines <- c(
  paste("Trend features tested:", nrow(trend)),
  paste("Trend features with BH q < 0.05:", sum(trend$q_value < 0.05, na.rm = TRUE)),
  paste("Pairwise tests:", nrow(pairwise)),
  paste("Pairwise tests with BH q < 0.05:", sum(pairwise$q_value < 0.05, na.rm = TRUE)),
  paste("Selected figure features passing BH q < 0.05:", sum(selected$q_value < 0.05, na.rm = TRUE))
)
writeLines(summary_lines, file.path(report_dir, "analysis_summary.txt"))

script_arg <- commandArgs()[grep("^--file=", commandArgs())][[1]]
script_path <- sub("^--file=", "", script_arg)
template <- file.path(dirname(normalizePath(script_path, mustWork = TRUE)),
                      "cellchat_pathology_report.Rmd")
if (!file.exists(template)) {
  template <- file.path(results_root, "scripts", "cellchat_pathology_report.Rmd")
}
rmarkdown::render(
  input = template,
  output_file = "cellchat_pathology_report.html",
  output_dir = report_dir,
  params = list(results_root = results_root),
  envir = new.env(parent = globalenv()),
  quiet = TRUE,
  output_options = list(self_contained = TRUE)
)

cat("Analysis complete. Main figure:",
    file.path(figure_dir, "cellchat_pathology_main_figure.pdf"), "\n")
cat("HTML report:", file.path(report_dir, "cellchat_pathology_report.html"), "\n")
