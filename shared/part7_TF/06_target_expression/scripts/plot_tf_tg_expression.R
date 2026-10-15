#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Matrix)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 4L) {
  stop("Usage: plot_tf_tg_expression.R <celltype_region> <EXP008_results-root> [celltype_specific|cross_celltype] [TF|ALL]")
}
pair <- args[[1]]
out_root <- args[[2]]
universe <- if (length(args) >= 3L) args[[3]] else "celltype_specific"
tf_selector <- if (length(args) >= 4L) args[[4]] else "ALL"
if (!universe %in% c("celltype_specific", "cross_celltype")) stop("Unknown analysis universe: ", universe)
adrega <- Sys.getenv("REGA_SEPARATE_ROOT", unset = "/project/zduren/durenlab/naqing_palmetto/REGA_separate")
shared_root <- Sys.getenv("REGA_SHARED_ROOT", unset = normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "../../..")))
adrega_cross <- shared_root
driver_file <- if (universe == "celltype_specific") {
  file.path(adrega, "07_tf/results/activity_pathology/all_celltype_region_TF_pathology_association.tsv")
} else {
  file.path(adrega_cross, "part7_TF/04DriverTF_lm/all_celltype_region_TF_pathology_association.tsv")
}
grn_dir <- Sys.getenv("REGA_SEPARATE_GRN_ROOT", unset = "")
if (universe == "celltype_specific" && !nzchar(grn_dir)) {
  stop("Set REGA_SEPARATE_GRN_ROOT to GRNs derived from the current independent-model matrices")
}
if (universe == "celltype_specific") {
  gex_file <- file.path(adrega, "01_prepare/data/cleaned_expression", paste0("cleaned_", pair, ".tsv"))
  gene_module_file <- file.path(adrega, "03_process_rega/module_assignments", pair, "genes_module_assignments.tsv")
  tf_module_file <- file.path(adrega, "07_tf/results/tf_assignment", paste0(pair, "_TF_module_assignments.tsv"))
  grn_file <- file.path(grn_dir, paste0(pair, "_TG_TF.rds"))
  assignment_entity <- "entity"
  assignment_pair_column <- NULL
  tf_assignment_entity <- "entity"
  module_source_label <- "pair-specific TF assignment table for the driver TF; gene assignment table for targets"
} else {
  gex_file <- file.path(adrega_cross, "part1_prepare_REGA/02cleaned_region_consistent_gex", paste0("cleaned_", pair, "_region_consistent.tsv"))
  region <- sub("^.*_", "", pair)
  gene_module_file <- file.path(adrega_cross, "part3_process_REGA_output/10-1gene_module_assignment_region_shareH_FC", region, paste0(region, "_ALL_cellgene_FC_assignment_with_unassigned.tsv"))
  tf_module_file <- file.path(adrega_cross, "part7_TF/02tf_assignment_FC", region, paste0(region, "_ALL_cell_TF_FC_assignment_with_unassigned.tsv"))
  grn_file <- file.path(shared_root, "part7_TF/05GRN/results", paste0(pair, "_TG_TF.rds"))
  assignment_entity <- "gene"
  assignment_pair_column <- "celltype_region"
  tf_assignment_entity <- "tf"
  module_source_label <- "region-wise TF assignment for the TF; region-wise gene assignment for targets"
}
pathology_file <- Sys.getenv("INDIVIDUAL_PATHOLOGY_FILE", unset = "/project/zduren/durenlab/naqing_palmetto/AD_Driver_TF/individual_pathology.tsv")
pathology_levels <- c("nonAD", "earlyAD", "lateAD")
pathology_colors <- c(nonAD = "#66AABB", earlyAD = "#EEAA55", lateAD = "#CC6677")
test_alternative <- "greater"
trend_p_greater <- function(fit, coefficient = "pathology_code") {
  t_stat <- summary(fit)$coefficients[coefficient, "t value"]
  stats::pt(t_stat, df = stats::df.residual(fit), lower.tail = FALSE)
}
dir.create(file.path(out_root, "audit"), recursive = TRUE, showWarnings = FALSE)

read_feature_matrix <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  lines <- trimws(readLines(path, warn = FALSE))
  lines <- lines[nzchar(lines)]
  parts <- strsplit(lines, "[[:space:]]+")
  header <- parts[[1]]
  rows <- parts[-1]
  if (!length(rows)) stop("No data rows in: ", path)
  row_lengths <- lengths(rows)
  if (all(row_lengths == length(header) + 1L)) {
    ids <- vapply(rows, `[[`, character(1), 1L)
    value_rows <- lapply(rows, `[`, -1L)
    value_names <- header
  } else if (all(row_lengths == length(header)) && header[[1]] %in% c("gene", "peak", "RE", "feature", "entity")) {
    ids <- vapply(rows, `[[`, character(1), 1L)
    value_rows <- lapply(rows, `[`, -1L)
    value_names <- header[-1]
  } else {
    stop("Unexpected field count in ", path, "; expected feature IDs plus sample columns")
  }
  values <- do.call(rbind, lapply(value_rows, function(x) {
    y <- suppressWarnings(as.numeric(x))
    bad <- is.na(y) & !(x %in% c("NA", "NaN", "nan", ""))
    if (any(bad)) stop("Non-numeric expression value in ", path, ": ", x[which(bad)[1]])
    y
  }))
  rownames(values) <- ids
  colnames(values) <- value_names
  if (anyDuplicated(rownames(values)) || anyDuplicated(colnames(values))) stop("Duplicate feature/sample IDs in ", path)
  values
}

read_assignment <- function(path, entity_column = "entity", pair_column = NULL, pair_value = NULL) {
  if (!file.exists(path)) return(list(map = NULL, status = "missing_file"))
  x <- read_tsv(path, show_col_types = FALSE,
                col_types = cols(.default = col_character()), na = c("", "NA", "NaN"))
  if (!all(c(entity_column, "assigned_module") %in% names(x)) || (!is.null(pair_column) && !pair_column %in% names(x))) return(list(map = NULL, status = "invalid_schema"))
  if (!is.null(pair_column)) x <- x |> filter(.data[[pair_column]] == .env$pair_value)
  x <- x |> transmute(entity = .data[[entity_column]], assigned_module) |> filter(!is.na(entity), entity != "")
  conflicts <- x |> group_by(entity) |>
    summarise(n_module = n_distinct(assigned_module[!is.na(assigned_module) & assigned_module != ""]), .groups = "drop") |>
    filter(n_module > 1L)
  if (nrow(conflicts)) return(list(map = NULL, status = paste0("conflicting_assignments:", paste(head(conflicts$entity, 5), collapse = ";"))))
  x <- x |> distinct(entity, .keep_all = TRUE)
  list(map = setNames(x$assigned_module, x$entity), status = "ok")
}

drivers <- read_tsv(driver_file, show_col_types = FALSE,
                    col_types = cols(.default = col_character())) |>
  mutate(lm_fdr_pathology = suppressWarnings(as.numeric(lm_fdr_pathology)),
         lm_beta_pathology = suppressWarnings(as.numeric(lm_beta_pathology)),
         lm_pval_pathology = suppressWarnings(as.numeric(lm_pval_pathology))) |>
  filter(celltype_region == pair, is.finite(lm_fdr_pathology), lm_fdr_pathology < 0.05) |>
  distinct(TF, .keep_all = TRUE)
if (tf_selector != "ALL") drivers <- drivers |> filter(TF == tf_selector)

if (!nrow(drivers)) {
  write_tsv(tibble(pair = pair, status = "skipped_no_significant_driver", n_driver_TF = 0L),
            file.path(out_root, "audit", paste0(pair, "_pair_status.tsv")))
  cat(pair, ": no significant driver TFs\n")
  quit(save = "no", status = 0L)
}
for (path in c(gex_file, pathology_file)) if (!file.exists(path)) stop("Missing required input: ", path)
# Missing module tables make the restricted mode unavailable, but must not
# prevent the unrestricted analysis from running. read_assignment() records
# the missing-file status for the restricted-mode audit below.
relation_file <- grn_file
if (!file.exists(relation_file)) stop("Missing pair-specific TG_TF matrix: ", relation_file)
tg_tf <- readRDS(relation_file)
if (!(is.matrix(tg_tf) || inherits(tg_tf, "Matrix")) || is.null(rownames(tg_tf)) || is.null(colnames(tg_tf))) {
  stop("TG_TF must be a gene-by-TF matrix with row/column names: ", relation_file)
}

gex <- read_feature_matrix(gex_file)
pathology <- read_tsv(pathology_file, show_col_types = FALSE,
                      col_types = cols(.default = col_character()), na = c("", "NA", "NaN"))
if (!all(c("ROSMAP_IndividualID", "Pathology") %in% names(pathology))) stop("Invalid individual pathology schema")
pathology <- pathology |> select(ROSMAP_IndividualID, Pathology) |> distinct()
if (anyDuplicated(pathology$ROSMAP_IndividualID)) stop("Conflicting pathology labels per individual")
donor_meta <- tibble(Individual = colnames(gex)) |>
  left_join(pathology, by = c("Individual" = "ROSMAP_IndividualID"))
if (anyNA(donor_meta$Pathology)) {
  write_tsv(filter(donor_meta, is.na(Pathology)), file.path(out_root, "audit", paste0(pair, "_unmatched_expression_donors.tsv")))
  stop("Expression-matrix donor IDs unmatched to pathology for ", pair)
}
donor_meta$Pathology <- factor(donor_meta$Pathology, levels = pathology_levels)
if (anyNA(donor_meta$Pathology)) stop("Unexpected pathology labels in ", pair)

gene_assignment <- read_assignment(gene_module_file, assignment_entity, assignment_pair_column, pair)
gene_module_map <- gene_assignment$map
tf_assignment <- read_assignment(tf_module_file, tf_assignment_entity, assignment_pair_column, pair)
tf_module_map <- tf_assignment$map
driver_module <- function(tf) {
  if (is.null(tf_module_map) || !(tf %in% names(tf_module_map))) return(NA_character_)
  z <- unname(tf_module_map[[tf]])
  if (!length(z) || is.na(z) || !nzchar(z)) NA_character_ else z
}
candidate_rows <- list(); score_rows <- list(); gene_expr_rows <- list(); tf_test_rows <- list(); tg_test_rows <- list(); audit_rows <- list()

for (i in seq_len(nrow(drivers))) {
  tf <- drivers$TF[[i]]
  drow <- drivers[i, , drop = FALSE]
  if (!(tf %in% colnames(tg_tf))) {
    for (mode in c("unrestricted", "module_restricted")) {
      audit_rows[[length(audit_rows) + 1L]] <- tibble(
        pair = pair, TF = tf, mode = mode, status = "TF_missing_from_TG_TF", module_source = module_source_label,
        TF_module = driver_module(tf), n_candidates = 0L, n_selected = 0L, n_expression_matched = 0L)
      tf_test_rows[[length(tf_test_rows) + 1L]] <- tibble(
        pair = pair, celltype = drow$celltype, region = drow$region, TF = tf, mode = mode,
        TF_module = driver_module(tf), n_selected_target_genes = 0L,
        n_expression_matched_target_genes = 0L, n_nonAD = 0L, n_earlyAD = 0L, n_lateAD = 0L,
        status = "TF_missing_from_TG_TF", test_alternative = test_alternative,
        beta_pathology = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_)
    }
    next
  }
  weights <- as.numeric(tg_tf[, tf]); names(weights) <- rownames(tg_tf)
  weights <- weights[is.finite(weights) & weights != 0]
  tf_module <- driver_module(tf)
  for (mode in c("unrestricted", "module_restricted")) {
    mode_weights <- weights
    status <- "selected"
    if (mode == "module_restricted") {
      if (gene_assignment$status != "ok" || tf_assignment$status != "ok" || is.na(tf_module)) {
        bad_assignment <- if (gene_assignment$status != "ok") gene_assignment$status else tf_assignment$status
        status <- if (gene_assignment$status != "ok" || tf_assignment$status != "ok") paste0("module_unavailable_", bad_assignment) else "TF_gene_module_missing"
        mode_weights <- numeric()
      } else {
        gene_modules <- unname(gene_module_map[names(mode_weights)])
        mode_weights <- mode_weights[!is.na(gene_modules) & gene_modules == tf_module]
        if (!length(mode_weights)) status <- "no_candidates_in_TF_gene_module"
      }
    }
    n_candidate <- length(mode_weights)
    if (n_candidate) {
      order_idx <- order(-mode_weights, names(mode_weights), method = "radix")
      n_select <- max(1L, ceiling(0.10 * n_candidate))
      selected <- names(mode_weights)[order_idx[seq_len(min(n_select, n_candidate))]]
    } else selected <- character()
    expr_genes <- intersect(selected, rownames(gex))
    if (n_candidate && !length(expr_genes)) status <- "selected_genes_absent_from_expression"
    if (length(expr_genes) && status == "selected") status <- "selected_with_expression"
    modules <- if (!is.null(gene_module_map)) unname(gene_module_map[selected]) else rep(NA_character_, length(selected))
    target_tbl <- tibble(pair = pair, celltype = drow$celltype, region = drow$region,
                         TF = tf, mode = mode, TF_module = tf_module,
                         target_gene = selected, signed_TF_TG_weight = unname(mode_weights[selected]),
                         rank = seq_along(selected), target_gene_module = modules,
                         expression_available = selected %in% rownames(gex),
                         driver_beta_pathology = drow$lm_beta_pathology,
                         driver_p_pathology = drow$lm_pval_pathology,
                         driver_fdr_pathology = drow$lm_fdr_pathology)
    if (length(selected)) candidate_rows[[length(candidate_rows) + 1L]] <- target_tbl
    audit_rows[[length(audit_rows) + 1L]] <- tibble(
      pair = pair, TF = tf, mode = mode, status = status,
      module_source = module_source_label,
      TF_module = tf_module, n_candidates = n_candidate,
      n_selected = length(selected), n_expression_matched = length(expr_genes),
      expression_missing_target_genes = length(setdiff(selected, rownames(gex)))
    )
    if (!length(expr_genes)) {
      tf_test_rows[[length(tf_test_rows) + 1L]] <- tibble(
        pair = pair, celltype = drow$celltype, region = drow$region, TF = tf, mode = mode,
        TF_module = tf_module, n_selected_target_genes = length(selected),
        n_expression_matched_target_genes = 0L, n_nonAD = 0L, n_earlyAD = 0L, n_lateAD = 0L,
        status = status, test_alternative = test_alternative,
        beta_pathology = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_)
      next
    }

    expression_matrix <- gex[expr_genes, , drop = FALSE]
    long <- expand.grid(target_gene = rownames(expression_matrix), Individual = colnames(expression_matrix),
                        stringsAsFactors = FALSE) |>
      mutate(expression = as.vector(expression_matrix), pair = pair, celltype = drow$celltype,
             region = drow$region, TF = tf, mode = mode) |>
      left_join(donor_meta, by = "Individual") |>
      left_join(target_tbl |> select(target_gene, rank, signed_TF_TG_weight, target_gene_module), by = "target_gene")
    gene_expr_rows[[length(gene_expr_rows) + 1L]] <- long

    gene_tests <- list()
    for (gene in rownames(expression_matrix)) {
      vals <- expression_matrix[gene, ]
      keep <- is.finite(vals)
      grp <- table(factor(donor_meta$Pathology[keep], levels = pathology_levels))
      pval <- beta <- ci_low <- ci_high <- NA_real_
      test_status <- if (grepl("^Vasc_", pair)) "descriptive_only_vascular" else "not_tested_low_group_n"
      if (!grepl("^Vasc_", pair) && all(grp >= 5L)) {
        dat <- tibble(y = as.numeric(vals[keep]), pathology_code = as.integer(donor_meta$Pathology[keep]) - 1L)
        fit <- lm(y ~ pathology_code, data = dat)
        coef <- summary(fit)$coefficients
        beta <- coef["pathology_code", "Estimate"]
        pval <- trend_p_greater(fit)
        ci <- confint(fit)["pathology_code", ]; ci_low <- ci[1]; ci_high <- ci[2]
        test_status <- "tested"
      }
      gene_tests[[length(gene_tests) + 1L]] <- tibble(
        pair = pair, celltype = drow$celltype, region = drow$region, TF = tf, mode = mode,
        target_gene = gene, target_gene_module = if (!is.null(gene_module_map)) unname(gene_module_map[gene]) else NA_character_,
        n_nonAD = as.integer(grp[1]), n_earlyAD = as.integer(grp[2]), n_lateAD = as.integer(grp[3]),
        status = test_status, test_alternative = test_alternative,
        beta_pathology = beta, ci_low = ci_low, ci_high = ci_high, p_value = pval
      )
    }
    gene_test_tbl <- bind_rows(gene_tests)
    gene_test_tbl$fdr_within_TF_pair_mode <- NA_real_
    valid <- which(is.finite(gene_test_tbl$p_value))
    if (length(valid)) gene_test_tbl$fdr_within_TF_pair_mode[valid] <- p.adjust(gene_test_tbl$p_value[valid], method = "BH")
    tg_test_rows[[length(tg_test_rows) + 1L]] <- gene_test_tbl

    # Standardize each target gene across individuals before taking the donor
    # mean so genes with larger stored expression values cannot dominate it.
    # Genes with fewer than two finite observations or zero variance have no
    # estimable z-score and are excluded from the donor summary.
    expression_z <- expression_matrix
    for (gene in rownames(expression_matrix)) {
      x <- expression_matrix[gene, ]
      ok <- is.finite(x)
      gene_sd <- if (sum(ok) >= 2L) sd(x[ok]) else NA_real_
      if (!is.finite(gene_sd) || gene_sd == 0) {
        expression_z[gene, ] <- NA_real_
      } else {
        z <- rep(NA_real_, length(x))
        z[ok] <- (x[ok] - mean(x[ok])) / gene_sd
        expression_z[gene, ] <- z
      }
    }
    donor_means <- vapply(seq_len(ncol(expression_z)), function(j) {
      x <- expression_z[, j]
      ok <- is.finite(x)
      if (!any(ok)) return(c(mean_expression = NA_real_, n_target_genes = 0))
      c(mean_expression = mean(x[ok]), n_target_genes = sum(ok))
    }, numeric(2))
    score_tbl <- donor_meta |>
      mutate(pair = pair, celltype = drow$celltype, region = drow$region, TF = tf, mode = mode,
             mean_target_expression = as.numeric(donor_means["mean_expression", ]),
             summary_scale = "mean_per_gene_zscore_across_individuals",
             n_target_genes = as.integer(donor_means["n_target_genes", ]),
             n_unstandardizable_target_genes = sum(rowSums(is.finite(expression_z)) == 0L),
             n_selected_target_genes = length(selected), n_expression_matched_target_genes = length(expr_genes))
    score_rows[[length(score_rows) + 1L]] <- score_tbl
    scored <- is.finite(score_tbl$mean_target_expression)
    n_by_group <- table(factor(score_tbl$Pathology[scored], levels = pathology_levels))
    beta <- pval <- ci_low <- ci_high <- NA_real_
    test_status <- if (grepl("^Vasc_", pair)) "descriptive_only_vascular" else "not_tested_low_group_n"
    if (!grepl("^Vasc_", pair) && all(n_by_group >= 5L)) {
      dat <- score_tbl[scored, ] |>
        transmute(y = mean_target_expression, pathology_code = as.integer(Pathology) - 1L)
      fit <- lm(y ~ pathology_code, data = dat)
      coef <- summary(fit)$coefficients
      beta <- coef["pathology_code", "Estimate"]
      pval <- trend_p_greater(fit)
      ci <- confint(fit)["pathology_code", ]; ci_low <- ci[1]; ci_high <- ci[2]
      test_status <- "tested"
    }
    tf_test_rows[[length(tf_test_rows) + 1L]] <- tibble(
      pair = pair, celltype = drow$celltype, region = drow$region, TF = tf, mode = mode,
      TF_module = tf_module, n_selected_target_genes = length(selected),
      n_expression_matched_target_genes = length(expr_genes),
      n_nonAD = as.integer(n_by_group[1]), n_earlyAD = as.integer(n_by_group[2]), n_lateAD = as.integer(n_by_group[3]),
      status = test_status, test_alternative = test_alternative,
      beta_pathology = beta, ci_low = ci_low, ci_high = ci_high, p_value = pval
    )
  }
}

audit_tbl <- bind_rows(audit_rows)
write_tsv(audit_tbl, file.path(out_root, "audit", paste0(pair, "_selection_audit.tsv")))
for (mode in c("unrestricted", "module_restricted")) {
  mode_dir <- file.path(out_root, mode)
  pair_dir <- file.path(mode_dir, pair)
  for (sub in c("target_genes", "expression_by_gene_individual", "mean_target_expression_by_individual", "regression_results", "plots")) {
    dir.create(file.path(mode_dir, sub, pair), recursive = TRUE, showWarnings = FALSE)
  }
  filter_mode <- function(xs, mode_value) {
    x <- bind_rows(xs)
    if (!"mode" %in% names(x)) return(x)
    dplyr::filter(x, .data$mode == .env$mode_value)
  }
  targets <- filter_mode(candidate_rows, mode)
  expr_long <- filter_mode(gene_expr_rows, mode)
  scores <- filter_mode(score_rows, mode)
  tf_results <- filter_mode(tf_test_rows, mode)
  tg_results <- filter_mode(tg_test_rows, mode)
  if (nrow(tf_results)) {
    tf_results$summary_scale <- "mean_per_gene_zscore_across_individuals"
    valid <- which(is.finite(tf_results$p_value))
    tf_results$fdr_within_pair_mode <- NA_real_
    if (length(valid)) tf_results$fdr_within_pair_mode[valid] <- p.adjust(tf_results$p_value[valid], method = "BH")
    for (j in seq_len(nrow(tf_results))) {
      tf <- tf_results$TF[j]
      z <- scores |> filter(TF == tf)
      if (!nrow(z)) next
      pval <- tf_results$p_value[j]; fdr <- tf_results$fdr_within_pair_mode[j]
      subtitle <- if (is.finite(pval)) {
        sprintf("LM greater-than p = %.3g; BH FDR = %.3g | n = %d / %d / %d",
                pval, fdr, tf_results$n_nonAD[j], tf_results$n_earlyAD[j], tf_results$n_lateAD[j])
      } else {
        sprintf("%s | n = %d / %d / %d", tf_results$status[j],
                tf_results$n_nonAD[j], tf_results$n_earlyAD[j], tf_results$n_lateAD[j])
      }
      p <- ggplot(z, aes(x = Pathology, y = mean_target_expression, fill = Pathology)) +
        geom_boxplot(outlier.shape = NA, width = 0.7) +
        geom_jitter(width = 0.12, size = 0.7, alpha = 0.55) +
        scale_x_discrete(limits = pathology_levels, drop = FALSE) +
        scale_fill_manual(values = pathology_colors, drop = FALSE) +
        labs(title = paste(tf, "target genes in", pair), subtitle = subtitle,
             x = "Pathology group", y = "Mean per-gene z-score across selected targets") +
        theme_classic(base_size = 9) + theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1))
      ggsave(file.path(mode_dir, "plots", pair, paste0(tf, "_", mode, "_target_expression_boxplot.pdf")),
             p, width = 3.5, height = 3.5, units = "in", limitsize = FALSE)
    }
  }
  if (nrow(targets)) write_tsv(targets, file.path(mode_dir, "target_genes", pair, paste0(pair, "_target_genes.tsv.gz")))
  if (nrow(expr_long)) write_tsv(expr_long, file.path(mode_dir, "expression_by_gene_individual", pair, paste0(pair, "_target_gene_expression_by_individual.tsv.gz")))
  if (nrow(scores)) write_tsv(scores, file.path(mode_dir, "mean_target_expression_by_individual", pair, paste0(pair, "_mean_target_expression_by_individual.tsv.gz")))
  if (nrow(tf_results)) write_tsv(tf_results, file.path(mode_dir, "regression_results", pair, paste0(pair, "_TF_level_regression.tsv")))
  if (nrow(tg_results)) write_tsv(tg_results, file.path(mode_dir, "regression_results", pair, paste0(pair, "_TG_level_regression.tsv.gz")))
}
write_tsv(tibble(pair = pair, status = "complete", n_driver_TF = nrow(drivers),
                 n_module_restricted_available = sum(audit_tbl$mode == "module_restricted" & audit_tbl$status %in% c("selected_with_expression", "selected"))),
          file.path(out_root, "audit", paste0(pair, "_pair_status.tsv")))
cat(pair, ": processed ", nrow(drivers), " driver TFs; output under ", out_root, "\n", sep = "")
