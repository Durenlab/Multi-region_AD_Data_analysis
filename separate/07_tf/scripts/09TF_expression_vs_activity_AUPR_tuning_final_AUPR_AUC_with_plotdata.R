rega_separate_root <- Sys.getenv("REGA_SEPARATE_ROOT")
if (!nzchar(rega_separate_root)) stop("Set REGA_SEPARATE_ROOT")
# ============================================================
# TF Expression vs TF Activity classification: ALL celltype-region pairs
#
# For each celltype-region:
#   1) Read Expression / Activity TF x individual matrices
#   2) Remove all-zero TFs separately, then intersect TFs
#   3) Align individuals using metadata
#   4) Tune alpha + lambda by pooled out-of-fold (OOF) AUPR from one stratified K-fold CV
#   5) Fix selected alpha + lambda
#   6) Re-train with LOOCV and predict each held-out individual
#   7) Compute BOTH precision-recall (PR) / AUPR and ROC / AUC for Expression and Activity
#
# Final outputs:
#   - one summary TSV for all celltype-region pairs (contains both AUPR and AUC)
#   - one big 36-panel PR figure and one big 36-panel ROC figure
#   - paired AUPR and AUC boxplots across all celltype-region pairs
#   - paired AUPR and AUC connected-point plots
#   - overall AUPR and AUC statistics
#   - plot-ready TSV data frames for regenerating every figure in a notebook
#
# NOTE:
# Hyperparameters are tuned once using the full dataset before LOOCV.
# Therefore this is an exploratory, non-nested estimate and can be
# slightly optimistic. The workflow is intentionally kept consistent
# with the user's current single-pair script.
# ============================================================

suppressPackageStartupMessages({
  library(glmnet)
  library(pROC)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# ============================================================
# 0. USER SETTINGS
# ============================================================

celltypes <- c("Ast", "Exc", "Inh", "Mic_Immune", "OPC", "Oli")
regions   <- c("AG", "EC", "HC", "MTC", "PFC", "TH")
#alpha_grid <- c(0, 0.25, 0.5, 0.75, 1)
#alpha_grid <- c(0.25, 0.5, 0.75)
alpha_grid <- c(0.5)

nfolds_tune <- 36
seed <- 20260902

# Tuning rule: lambda.min philosophy = choose the lambda with MAX pooled OOF AUPR
# (We do NOT use cv_fit$lambda.min, because the tuning metric is our custom pooled OOF AUPR.)
lambda_rule <- "lambda.min"
#lambda_rule <- "lambda.1se"


TF_root <- file.path(rega_separate_root, "07_tf/results/activity")
metadata_dir <- file.path(rega_separate_root, "01_prepare/data/metadata_by_group")

output_dir <- file.path(rega_separate_root, "07_tf/figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Output directory:\n", output_dir, "\n\n")
cat("tuning metric: pooled OOF AUPR\n")
cat("lambda rule:", lambda_rule, "(minimize 1 - pooled OOF AUPR)\n")
cat("alpha grid:", paste(alpha_grid, collapse = ", "), "\n\n")

# ============================================================
# 1. HELPER FUNCTIONS
# ============================================================

read_tf_matrix <- function(file) {
  x <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)

  tf_names <- as.character(x[[1]])
  x <- x[, -1, drop = FALSE]
  x[] <- lapply(x, function(z) suppressWarnings(as.numeric(z)))
  x <- as.matrix(x)

  rownames(x) <- tf_names
  colnames(x) <- trimws(colnames(x))
  x
}

make_stratified_foldid <- function(y01, k = 5, seed = 1) {
  foldid <- integer(length(y01))

  for (cls in c(0, 1)) {
    idx <- which(y01 == cls)

    set.seed(seed + cls + 100)
    idx <- sample(idx)

    foldid[idx] <- rep(seq_len(k), length.out = length(idx))
  }

  foldid
}

# Compute precision-recall curve and AUPR (average precision).
# Positive class is y01 == 1 (lateAD).
# Tied prediction scores are handled together at the same threshold.
compute_pr_metrics <- function(y01, prob) {
  y01 <- as.integer(y01)
  prob <- as.numeric(prob)

  valid <- is.finite(prob) & !is.na(y01)
  y01 <- y01[valid]
  prob <- prob[valid]

  n_pos <- sum(y01 == 1)
  n_neg <- sum(y01 == 0)

  if (length(y01) == 0 || n_pos == 0 || n_neg == 0) {
    return(list(
      aupr = NA_real_,
      curve = data.frame(Recall = numeric(0), Precision = numeric(0)),
      prevalence = if (length(y01) == 0) NA_real_ else mean(y01 == 1)
    ))
  }

  ord <- order(prob, decreasing = TRUE)
  y_sorted <- y01[ord]
  prob_sorted <- prob[ord]

  tp <- cumsum(y_sorted == 1)
  fp <- cumsum(y_sorted == 0)

  # Evaluate only at the end of each tied-score block.
  threshold_end <- c(which(diff(prob_sorted) != 0), length(prob_sorted))

  tp_t <- tp[threshold_end]
  fp_t <- fp[threshold_end]

  recall <- tp_t / n_pos
  precision <- tp_t / (tp_t + fp_t)

  # Include the conventional PR starting point (recall=0, precision=1).
  recall_curve <- c(0, recall)
  precision_curve <- c(1, precision)

  # Stepwise area under the PR curve = average precision (AP).
  aupr <- sum(diff(recall_curve) * precision_curve[-1])

  list(
    aupr = aupr,
    curve = data.frame(
      Recall = recall_curve,
      Precision = precision_curve
    ),
    prevalence = n_pos / length(y01)
  )
}

# Tune alpha and lambda using pooled out-of-fold AUPR.
#
# For each alpha:
#   1) Run ONE stratified K-fold cv.glmnet with keep = TRUE.
#   2) For every lambda, collect each sample's held-out prediction
#      from cv_fit$fit.preval.
#   3) Convert the binomial linear predictor (logit) to disease probability.
#   4) Pool ALL samples across folds and calculate ONE AUPR.
#
# Then choose the (alpha, lambda) pair with the largest pooled OOF AUPR.
# This is the "lambda.min" philosophy here: choose the CV-performance optimum,
# not the more strongly regularized 1-SE solution.
tune_alpha_lambda <- function(x, y01, foldid, alpha_grid, lambda_rule = "lambda.min") {

  if (lambda_rule != "lambda.min") {
    stop("This pooled-OOF-AUPR tuning function is configured for lambda.min only.")
  }

  result_list <- vector("list", length(alpha_grid))
  cvfit_list  <- vector("list", length(alpha_grid))

  for (i in seq_along(alpha_grid)) {

    a <- alpha_grid[i]

    cv_fit <- suppressWarnings(
      cv.glmnet(
        x = x,
        y = y01,
        family = "binomial",
        alpha = a,
        lambda.min.ratio = 1e-3,
        foldid = foldid,
        type.measure = "deviance",  # not used for selection; needed only by cv.glmnet
        standardize = TRUE,
        keep = TRUE
      )
    )

    # fit.preval: n_samples x n_lambda matrix of held-out predictions.
    # For binomial models these are on the logit/link scale.
    oof_eta  <- cv_fit$fit.preval
    oof_prob <- plogis(oof_eta)

    pooled_aupr <- vapply(
      seq_along(cv_fit$lambda),
      function(j) {
        prob_j <- oof_prob[, j]

        if (any(!is.finite(prob_j))) {
          return(NA_real_)
        }

        compute_pr_metrics(y01, prob_j)$aupr
      },
      numeric(1)
    )

    result_list[[i]] <- data.frame(
      alpha = a,
      lambda = cv_fit$lambda,
      pooled_oof_aupr = pooled_aupr,
      pooled_oof_error = 1 - pooled_aupr,
      nzero = cv_fit$nzero
    )

    cvfit_list[[i]] <- cv_fit
  }

  tuning_table <- dplyr::bind_rows(result_list)

  if (all(is.na(tuning_table$pooled_oof_aupr))) {
    stop("All pooled OOF AUPR values are NA during alpha/lambda tuning.")
  }

  # Custom lambda.min rule for pooled AUPR:
  # minimize (1 - pooled OOF AUPR), which is exactly equivalent to
  # maximizing pooled OOF AUPR. This is NOT cv_fit$lambda.min.
  best_index <- which.min(tuning_table$pooled_oof_error)

  tuning_table$selected <- FALSE
  tuning_table$selected[best_index] <- TRUE

  best_alpha  <- tuning_table$alpha[best_index]
  best_lambda <- tuning_table$lambda[best_index]
  best_aupr   <- tuning_table$pooled_oof_aupr[best_index]

  # Recover the tuning OOF disease probabilities for the selected model.
  alpha_index <- which(alpha_grid == best_alpha)[1]
  best_cv_fit <- cvfit_list[[alpha_index]]
  lambda_index <- which.min(abs(log(best_cv_fit$lambda) - log(best_lambda)))
  best_oof_prob <- plogis(best_cv_fit$fit.preval[, lambda_index])

  list(
    table = tuning_table,
    best_alpha = best_alpha,
    best_lambda = best_lambda,
    best_aupr = best_aupr,
    best_oof_prob = best_oof_prob,
    best_cv_fit = best_cv_fit
  )
}

format_p <- function(p) {
  if (length(p) == 0 || is.na(p) || !is.finite(p)) return("NA")
  if (p < 0.001) return(formatC(p, format = "e", digits = 2))
  formatC(p, format = "f", digits = 3)
}

# ============================================================
# 2. LOOP OVER ALL CELLTYPE-REGION PAIRS
# ============================================================

summary_list <- list()
pr_list <- list()
roc_list <- list()
prediction_list <- list()
tuning_list <- list()

pair_counter <- 0

for (region in regions) {
  for (celltype in celltypes) {

    pair_counter <- pair_counter + 1
    pair_name <- paste(celltype, region, sep = "_")

    cat("\n============================================================\n")
    cat(sprintf("[%d/%d] %s\n", pair_counter, length(celltypes) * length(regions), pair_name))
    cat("============================================================\n")

    pair_dir <- TF_root
    expression_file <- file.path(pair_dir, paste0(pair_name, "_TF_expression_pcc_cutoff.csv"))
    activity_file   <- file.path(pair_dir, paste0(pair_name, "_TF_activity_pcc_cutoff.csv"))
    metadata_file   <- file.path(metadata_dir, paste0(pair_name, "_metadata.tsv"))

    needed <- c(expression_file, activity_file, metadata_file)
    if (any(!file.exists(needed))) {
      cat("SKIP: missing input file(s)\n")
      print(needed[!file.exists(needed)])
      next
    }

    # --------------------------------------------------------
    # Read matrices
    # --------------------------------------------------------
    expr_raw <- read_tf_matrix(expression_file)
    act_raw  <- read_tf_matrix(activity_file)

    # Remove all-zero TFs separately.
    expr_keep <- rowSums(abs(expr_raw), na.rm = TRUE) != 0
    act_keep  <- rowSums(abs(act_raw),  na.rm = TRUE) != 0

    expr_nonzero <- expr_raw[expr_keep, , drop = FALSE]
    act_nonzero  <- act_raw[act_keep,  , drop = FALSE]

    # Use exactly the same TFs in both models.
    shared_tfs <- intersect(rownames(expr_nonzero), rownames(act_nonzero))

    if (length(shared_tfs) == 0) {
      cat("SKIP: no shared non-zero TFs\n")
      next
    }

    expr_shared <- expr_nonzero[shared_tfs, , drop = FALSE]
    act_shared  <- act_nonzero[shared_tfs,  , drop = FALSE]

    # --------------------------------------------------------
    # Metadata + individual alignment
    # --------------------------------------------------------
    meta <- read.delim(
      metadata_file,
      check.names = FALSE,
      stringsAsFactors = FALSE
    ) %>%
      transmute(
        individual = trimws(as.character(ROSMAP_IndividualID)),
        Pathology = trimws(as.character(Pathology)),
        celltype_meta = trimws(as.character(RNA.Class)),
        region_meta = trimws(as.character(BrainRegion))
      ) %>%
      filter(
        Pathology %in% c("nonAD", "lateAD"),
        celltype_meta == celltype,
        region_meta == region
      ) %>%
      distinct(individual, .keep_all = TRUE)

    common_ids <- Reduce(
      intersect,
      list(meta$individual, colnames(expr_shared), colnames(act_shared))
    )

    meta_sub <- meta %>%
      filter(individual %in% common_ids)

    ids <- meta_sub$individual

    if (length(ids) < 4) {
      cat("SKIP: too few matched individuals\n")
      next
    }

    expr <- expr_shared[, ids, drop = FALSE]
    act  <- act_shared[, ids, drop = FALSE]

    x_expr <- t(expr)
    x_act  <- t(act)

    y01 <- ifelse(meta_sub$Pathology == "lateAD", 1, 0)
    class_tab <- table(y01)

    if (length(class_tab) < 2 || any(class_tab < 2)) {
      cat("SKIP: too few individuals in one pathology group\n")
      next
    }

    if (any(!is.finite(x_expr)) || any(!is.finite(x_act))) {
      cat("SKIP: NA/Inf detected\n")
      next
    }

    cat("n individuals:", length(ids), "\n")
    cat("nonAD:", sum(y01 == 0), " lateAD:", sum(y01 == 1), "\n")
    cat("shared TFs:", ncol(x_expr), "\n")

    # --------------------------------------------------------
    # One stratified CV split for parameter tuning.
    # Same fold assignment for Expression and Activity.
    # --------------------------------------------------------
    #k_tune <- min(nfolds_tune, min(class_tab))
    k_tune =length(ids)  
    print(paste("k_tune is :",k_tune))
    if (k_tune < 2) {
      cat("SKIP: insufficient samples for tuning CV\n")
      next
    }

    pair_seed <- seed
    foldid_tune <- make_stratified_foldid(y01, k = k_tune, seed = pair_seed)

    # --------------------------------------------------------
    # Tune Expression and Activity separately.
    # --------------------------------------------------------
    expr_tune <- tune_alpha_lambda(
      x = x_expr,
      y01 = y01,
      foldid = foldid_tune,
      alpha_grid = alpha_grid,
      lambda_rule = lambda_rule
    )

    act_tune <- tune_alpha_lambda(
      x = x_act,
      y01 = y01,
      foldid = foldid_tune,
      alpha_grid = alpha_grid,
      lambda_rule = lambda_rule
    )

    cat("Expression best alpha/lambda/pooled OOF AUPR:",
        expr_tune$best_alpha, expr_tune$best_lambda, round(expr_tune$best_aupr, 4), "\n")
    cat("Activity best alpha/lambda/pooled OOF AUPR:  ",
        act_tune$best_alpha, act_tune$best_lambda, round(act_tune$best_aupr, 4), "\n")

    tuning_list[[paste0(pair_name, "_Expression")]] <- expr_tune$table %>%
      mutate(pair = pair_name, Data = "Expression", .before = 1)

    tuning_list[[paste0(pair_name, "_Activity")]] <- act_tune$table %>%
      mutate(pair = pair_name, Data = "Activity", .before = 1)

    # --------------------------------------------------------
    # Fixed-parameter LOOCV.
    # --------------------------------------------------------
    n <- length(ids)
    expr_prob <- rep(NA_real_, n)
    act_prob  <- rep(NA_real_, n)

    for (i in seq_len(n)) {

      train_idx <- setdiff(seq_len(n), i)
      test_idx <- i

      fit_expr <- suppressWarnings(
        glmnet(
          x = x_expr[train_idx, , drop = FALSE],
          y = y01[train_idx],
          family = "binomial",
          alpha = expr_tune$best_alpha,
          lambda = expr_tune$best_lambda,
          standardize = TRUE
        )
      )

      expr_prob[i] <- as.numeric(
        predict(
          fit_expr,
          newx = x_expr[test_idx, , drop = FALSE],
          s = expr_tune$best_lambda,
          type = "response"
        )
      )

      fit_act <- suppressWarnings(
        glmnet(
          x = x_act[train_idx, , drop = FALSE],
          y = y01[train_idx],
          family = "binomial",
          alpha = act_tune$best_alpha,
          lambda = act_tune$best_lambda,
          standardize = TRUE
        )
      )

      act_prob[i] <- as.numeric(
        predict(
          fit_act,
          newx = x_act[test_idx, , drop = FALSE],
          s = act_tune$best_lambda,
          type = "response"
        )
      )
    }

    pred_df <- data.frame(
      celltype = celltype,
      region = region,
      pair = pair_name,
      individual = ids,
      Pathology = meta_sub$Pathology,
      Expression_prob = expr_prob,
      Activity_prob = act_prob
    )

    prediction_list[[pair_name]] <- pred_df

    # --------------------------------------------------------
    # Final evaluation: BOTH PR/AUPR and ROC/AUC
    # --------------------------------------------------------
    # lateAD is the positive class (y = 1).

    # ---- PR / AUPR ----
    pr_expr <- compute_pr_metrics(y01, pred_df$Expression_prob)
    pr_act  <- compute_pr_metrics(y01, pred_df$Activity_prob)

    aupr_expr <- pr_expr$aupr
    aupr_act  <- pr_act$aupr

    # ---- ROC / AUC ----
    roc_expr <- pROC::roc(
      response = y01,
      predictor = pred_df$Expression_prob,
      levels = c(0, 1),
      direction = "<",
      quiet = TRUE
    )

    roc_act <- pROC::roc(
      response = y01,
      predictor = pred_df$Activity_prob,
      levels = c(0, 1),
      direction = "<",
      quiet = TRUE
    )

    auc_expr <- as.numeric(pROC::auc(roc_expr))
    auc_act  <- as.numeric(pROC::auc(roc_act))

    cat("Expression AUPR:", round(aupr_expr, 3), "\n")
    cat("Activity AUPR:  ", round(aupr_act, 3), "\n")
    cat("Expression AUC: ", round(auc_expr, 3), "\n")
    cat("Activity AUC:   ", round(auc_act, 3), "\n")

    summary_list[[pair_name]] <- data.frame(
      celltype = celltype,
      region = region,
      pair = pair_name,
      n_individuals = n,
      n_nonAD = sum(y01 == 0),
      n_lateAD = sum(y01 == 1),
      n_shared_TFs = ncol(x_expr),
      Expression_best_alpha = expr_tune$best_alpha,
      Expression_best_lambda = expr_tune$best_lambda,
      Expression_tuning_pooled_OOF_AUPR = expr_tune$best_aupr,
      Activity_best_alpha = act_tune$best_alpha,
      Activity_best_lambda = act_tune$best_lambda,
      Activity_tuning_pooled_OOF_AUPR = act_tune$best_aupr,
      Expression_AUPR = aupr_expr,
      Activity_AUPR = aupr_act,
      Activity_minus_Expression_AUPR = aupr_act - aupr_expr,
      Expression_AUC = auc_expr,
      Activity_AUC = auc_act,
      Activity_minus_Expression_AUC = auc_act - auc_expr,
      lateAD_prevalence = mean(y01 == 1)
    )

    pr_list[[paste0(pair_name, "_Expression")]] <- pr_expr$curve %>%
      mutate(
        celltype = celltype,
        region = region,
        pair = pair_name,
        Data = "Expression",
        .before = 1
      )

    pr_list[[paste0(pair_name, "_Activity")]] <- pr_act$curve %>%
      mutate(
        celltype = celltype,
        region = region,
        pair = pair_name,
        Data = "Activity",
        .before = 1
      )

    roc_list[[paste0(pair_name, "_Expression")]] <- data.frame(
      celltype = celltype,
      region = region,
      pair = pair_name,
      FPR = 1 - roc_expr$specificities,
      TPR = roc_expr$sensitivities,
      Data = "Expression"
    )

    roc_list[[paste0(pair_name, "_Activity")]] <- data.frame(
      celltype = celltype,
      region = region,
      pair = pair_name,
      FPR = 1 - roc_act$specificities,
      TPR = roc_act$sensitivities,
      Data = "Activity"
    )
  }
}

# ============================================================
# 3. COMBINE + SAVE TABLES
# ============================================================

summary_df <- bind_rows(summary_list)
pr_df <- bind_rows(pr_list)
roc_df <- bind_rows(roc_list)
prediction_df <- bind_rows(prediction_list)
tuning_df <- bind_rows(tuning_list)

if (nrow(summary_df) == 0) {
  stop("No celltype-region pair was successfully analyzed.")
}

write.table(
  summary_df,
  file.path(output_dir, "ALL_celltype_region_PR_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Same summary with a filename that explicitly indicates both final metrics.
write.table(
  summary_df,
  file.path(output_dir, "ALL_celltype_region_AUPR_AUC_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  prediction_df,
  file.path(output_dir, "ALL_celltype_region_LOOCV_predictions.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  tuning_df,
  file.path(output_dir, "ALL_celltype_region_parameter_tuning.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Save plot-ready data frames so figures can be regenerated without rerunning models.
write.table(
  pr_df,
  file.path(output_dir, "ALL_celltype_region_PR_curve_plotdata.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  roc_df,
  file.path(output_dir, "ALL_celltype_region_ROC_curve_plotdata.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ============================================================
# 4. BIG PR FIGURE: ALL CELLTYPE-REGION PANELS
# ============================================================

# Arrange region blocks, then cell types within each region.
pair_levels <- unlist(lapply(regions, function(r) paste(celltypes, r, sep = "_")))
pair_levels <- pair_levels[pair_levels %in% summary_df$pair]

pr_df$pair <- factor(pr_df$pair, levels = pair_levels)
roc_df$pair <- factor(roc_df$pair, levels = pair_levels)
summary_df$pair <- factor(summary_df$pair, levels = pair_levels)

pr_df <- pr_df %>% arrange(pair, Data, Recall, Precision)
roc_df <- roc_df %>% arrange(pair, Data, FPR, TPR)

annotation_df <- summary_df %>%
  mutate(
    label = paste0(
      "Expr=", sprintf("%.2f", Expression_AUPR), "\n",
      "Act=", sprintf("%.2f", Activity_AUPR)
    ),
    x = 0.97,
    y = 0.05
  )

baseline_df <- summary_df %>%
  transmute(
    pair = pair,
    baseline = lateAD_prevalence
  )

p_pr_all <- ggplot(
  pr_df,
  aes(x = Recall, y = Precision, color = Data, linetype = Data)
) +
  geom_hline(
    data = baseline_df,
    aes(yintercept = baseline),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.3
  ) +
  geom_step(linewidth = 0.75) +
  geom_text(
    data = annotation_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 0,
    size = 2.2
  ) +
  facet_wrap(~ pair, ncol = 6) +
  scale_color_manual(
    values = c(
      "Expression" = "#E76F51",
      "Activity" = "#159D9C"
    )
  ) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(
    title = "TF expression vs TF activity: nonAD vs lateAD classification",
    subtitle = paste0(
      "Tune alpha/lambda by pooled OOF AUPR (", lambda_rule,
      "), then fixed-parameter LOOCV; dashed line = lateAD prevalence"
    ),
    x = "Recall",
    y = "Precision",
    color = NULL,
    linetype = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 8),
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p_pr_all)

ggsave(
  file.path(output_dir, "ALL_celltype_region_PR_bigplot.pdf"),
  p_pr_all,
  width = 13,
  height = 13,
  units = "in",
  limitsize = FALSE
)

ggsave(
  file.path(output_dir, "ALL_celltype_region_PR_bigplot.png"),
  p_pr_all,
  width = 13,
  height = 13,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)


# ============================================================
# 4B. BIG ROC FIGURE: ALL CELLTYPE-REGION PANELS
# ============================================================

annotation_auc_df <- summary_df %>%
  mutate(
    label = paste0(
      "Expr=", sprintf("%.2f", Expression_AUC), "\n",
      "Act=", sprintf("%.2f", Activity_AUC)
    ),
    x = 0.97,
    y = 0.05
  )

p_roc_all <- ggplot(
  roc_df,
  aes(x = FPR, y = TPR, color = Data, linetype = Data)
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.3
  ) +
  geom_step(linewidth = 0.75) +
  geom_text(
    data = annotation_auc_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 0,
    size = 2.2
  ) +
  facet_wrap(~ pair, ncol = 6) +
  scale_color_manual(
    values = c(
      "Expression" = "#E76F51",
      "Activity" = "#159D9C"
    )
  ) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(
    title = "TF expression vs TF activity: nonAD vs lateAD classification",
    subtitle = paste0(
      "Parameters tuned by pooled OOF AUPR (", lambda_rule,
      "); final fixed-parameter LOOCV evaluated by ROC/AUC"
    ),
    x = "False positive rate (1 - specificity)",
    y = "True positive rate (sensitivity)",
    color = NULL,
    linetype = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 8),
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p_roc_all)

ggsave(
  file.path(output_dir, "ALL_celltype_region_ROC_bigplot.pdf"),
  p_roc_all,
  width = 13,
  height = 13,
  units = "in",
  limitsize = FALSE
)

ggsave(
  file.path(output_dir, "ALL_celltype_region_ROC_bigplot.png"),
  p_roc_all,
  width = 13,
  height = 13,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)

# ============================================================
# 5. AUPR BOXPLOT ACROSS ALL CELLTYPE-REGION PAIRS
# ============================================================

aupr_long <- summary_df %>%
  select(pair, celltype, region, Expression_AUPR, Activity_AUPR) %>%
  pivot_longer(
    cols = c(Expression_AUPR, Activity_AUPR),
    names_to = "Data",
    values_to = "AUPR"
  ) %>%
  mutate(
    Data = recode(
      Data,
      Expression_AUPR = "Expression",
      Activity_AUPR = "Activity"
    ),
    Data = factor(Data, levels = c("Expression", "Activity"))
  )

write.table(
  aupr_long,
  file.path(output_dir, "ALL_celltype_region_AUPR_summary_plotdata.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# One-sided paired test: Activity AUPR > Expression AUPR.
wilcox_aupr_p <- wilcox.test(
  summary_df$Activity_AUPR,
  summary_df$Expression_AUPR,
  paired = TRUE,
  alternative = "greater",
  exact = FALSE
)$p.value

cat("\n============================================================\n")
cat("Overall paired AUPR comparison\n")
cat("============================================================\n")
cat("Number of matched celltype-region pairs:", nrow(summary_df), "\n")
cat("Mean Expression AUPR:", mean(summary_df$Expression_AUPR), "\n")
cat("Mean Activity AUPR:  ", mean(summary_df$Activity_AUPR), "\n")
cat("Median Expression AUPR:", median(summary_df$Expression_AUPR), "\n")
cat("Median Activity AUPR:  ", median(summary_df$Activity_AUPR), "\n")
cat("Paired Wilcoxon p (Activity > Expression):", format_p(wilcox_aupr_p), "\n")

# Two boxplots + all individual celltype-region points.
p_aupr_box <- ggplot(
  aupr_long,
  aes(x = Data, y = AUPR, color = Data)
) +
  geom_boxplot(
    width = 0.45,
    outlier.shape = NA,
    linewidth = 0.7,
    alpha = 0
  ) +
  geom_jitter(
    width = 0.08,
    height = 0,
    size = 2,
    alpha = 0.7
  ) +
  annotate(
    "text",
    x = 1.5,
    y = 1.02,
    label = paste0("Paired Wilcoxon p = ", format_p(wilcox_aupr_p)),
    size = 4
  ) +
  scale_color_manual(
    values = c(
      "Expression" = "#E76F51",
      "Activity" = "#159D9C"
    )
  ) +
  coord_cartesian(ylim = c(0, 1.08), clip = "off") +
  labs(
    title = "AUPR across all celltype-region datasets",
    subtitle = paste0("n = ", nrow(summary_df), " matched celltype-region pairs"),
    x = NULL,
    y = "AUPR"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

print(p_aupr_box)

ggsave(
  file.path(output_dir, "ALL_celltype_region_AUPR_boxplot_Wilcoxon.pdf"),
  p_aupr_box,
  width = 5.3,
  height = 5.3,
  units = "in"
)

ggsave(
  file.path(output_dir, "ALL_celltype_region_AUPR_boxplot_Wilcoxon.png"),
  p_aupr_box,
  width = 5.3,
  height = 5.3,
  units = "in",
  dpi = 300
)

# Optional paired version: same celltype-region connected by a line.
p_aupr_paired <- ggplot(
  aupr_long,
  aes(x = Data, y = AUPR, group = pair)
) +
  geom_line(alpha = 0.30, linewidth = 0.5) +
  geom_point(aes(color = Data), size = 2, alpha = 0.8) +
  geom_boxplot(
    aes(group = Data, color = Data),
    width = 0.45,
    outlier.shape = NA,
    fill = NA,
    linewidth = 0.7
  ) +
  annotate(
    "text",
    x = 1.5,
    y = 1.02,
    label = paste0("Paired Wilcoxon p = ", format_p(wilcox_aupr_p)),
    size = 4
  ) +
  scale_color_manual(
    values = c(
      "Expression" = "#E76F51",
      "Activity" = "#159D9C"
    )
  ) +
  coord_cartesian(ylim = c(0, 1.08), clip = "off") +
  labs(
    title = "Paired AUPR comparison across celltype-region datasets",
    subtitle = "Each line connects Expression and Activity from the same celltype-region",
    x = NULL,
    y = "AUPR"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

print(p_aupr_paired)

ggsave(
  file.path(output_dir, "ALL_celltype_region_AUPR_paired_Wilcoxon.pdf"),
  p_aupr_paired,
  width = 5.5,
  height = 5.3,
  units = "in"
)

ggsave(
  file.path(output_dir, "ALL_celltype_region_AUPR_paired_Wilcoxon.png"),
  p_aupr_paired,
  width = 5.5,
  height = 5.3,
  units = "in",
  dpi = 300
)


# ============================================================
# 5B. AUC BOXPLOT ACROSS ALL CELLTYPE-REGION PAIRS
# ============================================================

auc_long <- summary_df %>%
  select(pair, celltype, region, Expression_AUC, Activity_AUC) %>%
  pivot_longer(
    cols = c(Expression_AUC, Activity_AUC),
    names_to = "Data",
    values_to = "AUC"
  ) %>%
  mutate(
    Data = recode(
      Data,
      Expression_AUC = "Expression",
      Activity_AUC = "Activity"
    ),
    Data = factor(Data, levels = c("Expression", "Activity"))
  )

write.table(
  auc_long,
  file.path(output_dir, "ALL_celltype_region_AUC_summary_plotdata.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# One-sided paired test: Activity AUC > Expression AUC.
wilcox_auc_p <- wilcox.test(
  summary_df$Activity_AUC,
  summary_df$Expression_AUC,
  paired = TRUE,
  alternative = "greater",
  exact = FALSE
)$p.value

cat("\n============================================================\n")
cat("Overall paired AUC comparison\n")
cat("============================================================\n")
cat("Number of matched celltype-region pairs:", nrow(summary_df), "\n")
cat("Mean Expression AUC:", mean(summary_df$Expression_AUC), "\n")
cat("Mean Activity AUC:  ", mean(summary_df$Activity_AUC), "\n")
cat("Median Expression AUC:", median(summary_df$Expression_AUC), "\n")
cat("Median Activity AUC:  ", median(summary_df$Activity_AUC), "\n")
cat("Paired Wilcoxon p (Activity > Expression):", format_p(wilcox_auc_p), "\n")

p_auc_box <- ggplot(
  auc_long,
  aes(x = Data, y = AUC, color = Data)
) +
  geom_boxplot(
    width = 0.45,
    outlier.shape = NA,
    linewidth = 0.7,
    alpha = 0
  ) +
  geom_jitter(
    width = 0.08,
    height = 0,
    size = 2,
    alpha = 0.7
  ) +
  annotate(
    "text",
    x = 1.5,
    y = 1.02,
    label = paste0("Paired Wilcoxon p = ", format_p(wilcox_auc_p)),
    size = 4
  ) +
  scale_color_manual(
    values = c(
      "Expression" = "#E76F51",
      "Activity" = "#159D9C"
    )
  ) +
  coord_cartesian(ylim = c(0, 1.08), clip = "off") +
  labs(
    title = "AUC across all celltype-region datasets",
    subtitle = paste0("n = ", nrow(summary_df), " matched celltype-region pairs"),
    x = NULL,
    y = "AUC"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

print(p_auc_box)

ggsave(
  file.path(output_dir, "ALL_celltype_region_AUC_boxplot_Wilcoxon.pdf"),
  p_auc_box,
  width = 5.3,
  height = 5.3,
  units = "in"
)

ggsave(
  file.path(output_dir, "ALL_celltype_region_AUC_boxplot_Wilcoxon.png"),
  p_auc_box,
  width = 5.3,
  height = 5.3,
  units = "in",
  dpi = 300
)

p_auc_paired <- ggplot(
  auc_long,
  aes(x = Data, y = AUC, group = pair)
) +
  geom_line(alpha = 0.30, linewidth = 0.5) +
  geom_point(aes(color = Data), size = 2, alpha = 0.8) +
  geom_boxplot(
    aes(group = Data, color = Data),
    width = 0.45,
    outlier.shape = NA,
    fill = NA,
    linewidth = 0.7
  ) +
  annotate(
    "text",
    x = 1.5,
    y = 1.02,
    label = paste0("Paired Wilcoxon p = ", format_p(wilcox_auc_p)),
    size = 4
  ) +
  scale_color_manual(
    values = c(
      "Expression" = "#E76F51",
      "Activity" = "#159D9C"
    )
  ) +
  coord_cartesian(ylim = c(0, 1.08), clip = "off") +
  labs(
    title = "Paired AUC comparison across celltype-region datasets",
    subtitle = "Each line connects Expression and Activity from the same celltype-region",
    x = NULL,
    y = "AUC"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

print(p_auc_paired)

ggsave(
  file.path(output_dir, "ALL_celltype_region_AUC_paired_Wilcoxon.pdf"),
  p_auc_paired,
  width = 5.5,
  height = 5.3,
  units = "in"
)

ggsave(
  file.path(output_dir, "ALL_celltype_region_AUC_paired_Wilcoxon.png"),
  p_auc_paired,
  width = 5.5,
  height = 5.3,
  units = "in",
  dpi = 300
)

# ============================================================
# 6. SAVE OVERALL STATISTICS
# ============================================================

overall_aupr_stats <- data.frame(
  n_pairs = nrow(summary_df),
  mean_Expression_AUPR = mean(summary_df$Expression_AUPR),
  mean_Activity_AUPR = mean(summary_df$Activity_AUPR),
  median_Expression_AUPR = median(summary_df$Expression_AUPR),
  median_Activity_AUPR = median(summary_df$Activity_AUPR),
  mean_Activity_minus_Expression_AUPR = mean(summary_df$Activity_minus_Expression_AUPR),
  median_Activity_minus_Expression_AUPR = median(summary_df$Activity_minus_Expression_AUPR),
  n_Activity_better_AUPR = sum(summary_df$Activity_AUPR > summary_df$Expression_AUPR),
  n_Expression_better_AUPR = sum(summary_df$Expression_AUPR > summary_df$Activity_AUPR),
  paired_Wilcoxon_p_Activity_greater = wilcox_aupr_p,
  tuning_metric = "pooled OOF AUPR",
  lambda_rule = lambda_rule
)

write.table(
  overall_aupr_stats,
  file.path(output_dir, "ALL_celltype_region_overall_AUPR_statistics.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

overall_auc_stats <- data.frame(
  n_pairs = nrow(summary_df),
  mean_Expression_AUC = mean(summary_df$Expression_AUC),
  mean_Activity_AUC = mean(summary_df$Activity_AUC),
  median_Expression_AUC = median(summary_df$Expression_AUC),
  median_Activity_AUC = median(summary_df$Activity_AUC),
  mean_Activity_minus_Expression_AUC = mean(summary_df$Activity_minus_Expression_AUC),
  median_Activity_minus_Expression_AUC = median(summary_df$Activity_minus_Expression_AUC),
  n_Activity_better_AUC = sum(summary_df$Activity_AUC > summary_df$Expression_AUC),
  n_Expression_better_AUC = sum(summary_df$Expression_AUC > summary_df$Activity_AUC),
  paired_Wilcoxon_p_Activity_greater = wilcox_auc_p,
  tuning_metric = "pooled OOF AUPR",
  lambda_rule = lambda_rule
)

write.table(
  overall_auc_stats,
  file.path(output_dir, "ALL_celltype_region_overall_AUC_statistics.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# One compact combined table for convenience.
overall_both_stats <- data.frame(
  n_pairs = nrow(summary_df),
  mean_Expression_AUPR = mean(summary_df$Expression_AUPR),
  mean_Activity_AUPR = mean(summary_df$Activity_AUPR),
  paired_Wilcoxon_AUPR_p_Activity_greater = wilcox_aupr_p,
  mean_Expression_AUC = mean(summary_df$Expression_AUC),
  mean_Activity_AUC = mean(summary_df$Activity_AUC),
  paired_Wilcoxon_AUC_p_Activity_greater = wilcox_auc_p,
  tuning_metric = "pooled OOF AUPR",
  lambda_rule = lambda_rule
)

write.table(
  overall_both_stats,
  file.path(output_dir, "ALL_celltype_region_overall_AUPR_AUC_statistics.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("\n============================================================\n")
cat("DONE\n")
cat("Output directory:\n", output_dir, "\n")
cat("Main outputs:\n")
cat("  ALL_celltype_region_PR_summary.tsv\n")
cat("  ALL_celltype_region_AUPR_AUC_summary.tsv\n")
cat("  ALL_celltype_region_LOOCV_predictions.tsv\n")
cat("  ALL_celltype_region_parameter_tuning.tsv\n")
cat("  ALL_celltype_region_PR_curve_plotdata.tsv\n")
cat("  ALL_celltype_region_ROC_curve_plotdata.tsv\n")
cat("  ALL_celltype_region_AUPR_summary_plotdata.tsv\n")
cat("  ALL_celltype_region_AUC_summary_plotdata.tsv\n")
cat("  ALL_celltype_region_PR_bigplot.pdf\n")
cat("  ALL_celltype_region_ROC_bigplot.pdf\n")
cat("  ALL_celltype_region_AUPR_boxplot_Wilcoxon.pdf\n")
cat("  ALL_celltype_region_AUPR_paired_Wilcoxon.pdf\n")
cat("  ALL_celltype_region_AUC_boxplot_Wilcoxon.pdf\n")
cat("  ALL_celltype_region_AUC_paired_Wilcoxon.pdf\n")
cat("  ALL_celltype_region_overall_AUPR_statistics.tsv\n")
cat("  ALL_celltype_region_overall_AUC_statistics.tsv\n")
cat("  ALL_celltype_region_overall_AUPR_AUC_statistics.tsv\n")
cat("============================================================\n")
