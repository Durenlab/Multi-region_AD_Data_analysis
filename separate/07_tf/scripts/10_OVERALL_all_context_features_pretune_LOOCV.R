#!/usr/bin/env Rscript
rega_separate_root <- Sys.getenv("REGA_SEPARATE_ROOT")
if (!nzchar(rega_separate_root)) stop("Set REGA_SEPARATE_ROOT")

# ============================================================
# OVERALL TF Expression vs TF Activity classifier
#
# Goal:
#   Build ONE overall classifier across all celltype-region contexts.
#
# Feature construction:
#   - Keep only individuals present in ALL 36 celltype-region datasets,
#     in BOTH Expression and Activity matrices, and in metadata.
#   - Within each celltype-region, after restricting to this common cohort:
#       1) remove all-zero TFs separately from Expression and Activity
#       2) intersect TFs so both representations use the same TFs
#       3) remove TFs with zero variance in either representation
#   - Rename each TF as: celltype_region_TFname
#     Thus, the same TF in different cellular/regional contexts is treated
#     as a distinct feature.
#   - Stack all celltype-region feature blocks by rows, then transpose to
#     obtain individual x feature matrices.
#
# Classification logic (same finalized framework):
#   - Penalized logistic regression, alpha fixed at 1 (lasso)
#   - Tune lambda using pooled out-of-fold AUPR
#   - Fix alpha/lambda
#   - Perform LOOCV for final held-out predictions
#   - lateAD is the positive class
#
# Primary evaluation:
#   - PR curve / AUPR
#   - Paired permutation test for Activity > Expression AUPR
#   - Stratified paired bootstrap CI for delta AUPR
#
# Complementary evaluation:
#   - ROC / AUC
#   - Paired DeLong test for Expression vs Activity AUC
#
# NOTE:
# Hyperparameters are selected using all common individuals before the final
# LOOCV. Therefore the resulting LOOCV performance is non-nested and may be
# optimistic in absolute terms. The main purpose is the matched comparison
# of TF Expression versus TF Activity under the same cohort and feature IDs.
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

# Finalized elastic-net setting.
alpha_grid <- c(1)
lambda_rule <- "lambda.min"
seed <- 20260902

# Resampling used only for inference around the FINAL OOF predictions.
n_boot <- 2000
n_perm <- 10000

TF_root <- Sys.getenv(
  "REGA_TF_ACTIVITY_ROOT",
  unset = file.path(rega_separate_root, "07_tf/results/activity")
)
metadata_dir <- Sys.getenv(
  "REGA_TF_METADATA_DIR",
  unset = file.path(rega_separate_root, "01_prepare/data/metadata_by_group")
)

output_dir <- Sys.getenv(
  "REGA_TF_COMPARISON_OUTPUT_DIR",
  unset = file.path(rega_separate_root, "07_tf/results/classification_overall")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

expr_color <- "#E76F51"
act_color  <- "#159D9C"

cat("Output directory:\n", output_dir, "\n\n")
cat("Feature design: celltype_region_TFname\n")
cat("alpha:", paste(alpha_grid, collapse = ", "), "\n")
cat("tuning metric: pooled OOF AUPR\n\n")

# ============================================================
# 1. HELPER FUNCTIONS
# ============================================================

read_tf_matrix <- function(file) {
  x <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)

  if (ncol(x) < 2) {
    stop("TF matrix has fewer than 2 columns: ", file)
  }

  tf_names <- as.character(x[[1]])
  x <- x[, -1, drop = FALSE]
  x[] <- lapply(x, function(z) suppressWarnings(as.numeric(z)))
  x <- as.matrix(x)

  rownames(x) <- make.unique(tf_names)
  colnames(x) <- trimws(colnames(x))
  x
}

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

  threshold_end <- c(which(diff(prob_sorted) != 0), length(prob_sorted))

  tp_t <- tp[threshold_end]
  fp_t <- fp[threshold_end]

  recall <- tp_t / n_pos
  precision <- tp_t / (tp_t + fp_t)

  recall_curve <- c(0, recall)
  precision_curve <- c(1, precision)

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

# True LOOCV fold IDs for pooled OOF tuning.
# Each individual is held out exactly once during the tuning CV.
make_loocv_foldid <- function(n) {
  seq_len(n)
}

# Tune lambda using pooled OOF AUPR.
# alpha is fixed at 1 by alpha_grid, but the function remains general.
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
        type.measure = "deviance",  # not used for selection
        standardize = TRUE,
        keep = TRUE
      )
    )

    oof_eta <- cv_fit$fit.preval
    oof_prob <- plogis(oof_eta)

    pooled_aupr <- vapply(
      seq_along(cv_fit$lambda),
      function(j) {
        prob_j <- oof_prob[, j]
        if (any(!is.finite(prob_j))) return(NA_real_)
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

  tuning_table <- bind_rows(result_list)

  if (all(is.na(tuning_table$pooled_oof_aupr))) {
    stop("All pooled OOF AUPR values are NA during tuning.")
  }

  best_index <- which.min(tuning_table$pooled_oof_error)
  tuning_table$selected <- FALSE
  tuning_table$selected[best_index] <- TRUE

  best_alpha  <- tuning_table$alpha[best_index]
  best_lambda <- tuning_table$lambda[best_index]
  best_aupr   <- tuning_table$pooled_oof_aupr[best_index]

  list(
    table = tuning_table,
    best_alpha = best_alpha,
    best_lambda = best_lambda,
    best_aupr = best_aupr
  )
}

format_p <- function(p) {
  if (length(p) == 0 || is.na(p) || !is.finite(p)) return("NA")
  if (p < 0.001) return(formatC(p, format = "e", digits = 2))
  formatC(p, format = "f", digits = 3)
}

# Paired stratified bootstrap: resample lateAD and nonAD donors separately,
# using the SAME sampled donor indices for Expression and Activity.
bootstrap_metric_difference <- function(
    y01,
    expr_prob,
    act_prob,
    metric = c("AUPR", "AUC"),
    B = 2000,
    seed = 1
) {
  metric <- match.arg(metric)

  idx0 <- which(y01 == 0)
  idx1 <- which(y01 == 1)

  set.seed(seed)
  delta <- rep(NA_real_, B)

  for (b in seq_len(B)) {
    boot_idx <- c(
      sample(idx0, length(idx0), replace = TRUE),
      sample(idx1, length(idx1), replace = TRUE)
    )

    yb <- y01[boot_idx]
    eb <- expr_prob[boot_idx]
    ab <- act_prob[boot_idx]

    if (metric == "AUPR") {
      e <- compute_pr_metrics(yb, eb)$aupr
      a <- compute_pr_metrics(yb, ab)$aupr
    } else {
      e <- as.numeric(pROC::auc(
        pROC::roc(yb, eb, levels = c(0, 1), direction = "<", quiet = TRUE)
      ))
      a <- as.numeric(pROC::auc(
        pROC::roc(yb, ab, levels = c(0, 1), direction = "<", quiet = TRUE)
      ))
    }

    delta[b] <- a - e
  }

  delta <- delta[is.finite(delta)]

  list(
    delta = delta,
    lower = unname(quantile(delta, 0.025, na.rm = TRUE)),
    upper = unname(quantile(delta, 0.975, na.rm = TRUE))
  )
}

# Paired permutation test for AUPR difference.
# Under the null, Expression/Activity prediction labels are exchangeable
# within each donor. Alternative: Activity AUPR > Expression AUPR.
paired_permutation_aupr <- function(
    y01,
    expr_prob,
    act_prob,
    B = 10000,
    seed = 1
) {
  observed <- compute_pr_metrics(y01, act_prob)$aupr -
    compute_pr_metrics(y01, expr_prob)$aupr

  n <- length(y01)
  set.seed(seed)
  null_delta <- rep(NA_real_, B)

  for (b in seq_len(B)) {
    swap <- runif(n) < 0.5

    e_perm <- expr_prob
    a_perm <- act_prob

    e_perm[swap] <- act_prob[swap]
    a_perm[swap] <- expr_prob[swap]

    null_delta[b] <- compute_pr_metrics(y01, a_perm)$aupr -
      compute_pr_metrics(y01, e_perm)$aupr
  }

  p_greater <- (1 + sum(null_delta >= observed, na.rm = TRUE)) /
    (1 + sum(is.finite(null_delta)))

  list(
    observed_delta = observed,
    p_greater = p_greater,
    null_delta = null_delta
  )
}

# ============================================================
# 2. READ ALL 36 DATASETS AND FIND THE COMMON DONOR COHORT
# ============================================================

pair_names <- unlist(lapply(regions, function(r) paste(celltypes, r, sep = "_")))

expr_raw_list <- list()
act_raw_list <- list()
meta_list <- list()
donor_sets <- list()

for (pair_name in pair_names) {

  parts <- strsplit(pair_name, "_")[[1]]

  # Because Mic_Immune contains an underscore, infer celltype/region by suffix.
  region <- sub("^.*_", "", pair_name)
  celltype <- sub(paste0("_", region, "$"), "", pair_name)

  expression_file <- file.path(
    TF_root,
    paste0(pair_name, "_TF_expression_pcc_cutoff.csv")
  )
  activity_file <- file.path(
    TF_root,
    paste0(pair_name, "_TF_activity_pcc_cutoff.csv")
  )
  metadata_file <- file.path(
    metadata_dir,
    paste0(pair_name, "_metadata.tsv")
  )

  needed <- c(expression_file, activity_file, metadata_file)
  if (any(!file.exists(needed))) {
    stop(
      "OVERALL model requires all celltype-region datasets. Missing file(s) for ",
      pair_name, ":\n",
      paste(needed[!file.exists(needed)], collapse = "\n")
    )
  }

  expr_raw <- read_tf_matrix(expression_file)
  act_raw  <- read_tf_matrix(activity_file)

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

  usable_ids <- Reduce(
    intersect,
    list(colnames(expr_raw), colnames(act_raw), meta$individual)
  )

  if (length(usable_ids) == 0) {
    stop("No usable donors in ", pair_name)
  }

  expr_raw_list[[pair_name]] <- expr_raw
  act_raw_list[[pair_name]] <- act_raw
  meta_list[[pair_name]] <- meta
  donor_sets[[pair_name]] <- usable_ids

  cat(
    pair_name,
    ": usable donors = ", length(usable_ids),
    ", Expression TFs = ", nrow(expr_raw),
    ", Activity TFs = ", nrow(act_raw),
    "\n",
    sep = ""
  )
}

common_ids <- Reduce(intersect, donor_sets)

if (length(common_ids) < 4) {
  stop("Fewer than 4 donors are present in all 36 celltype-region datasets.")
}

cat("\n============================================================\n")
cat("COMMON DONOR COHORT\n")
cat("============================================================\n")
cat("Donors present in ALL 36 datasets:", length(common_ids), "\n")

# Check pathology labels are consistent across all metadata files.
label_long <- bind_rows(lapply(names(meta_list), function(pair_name) {
  meta_list[[pair_name]] %>%
    filter(individual %in% common_ids) %>%
    select(individual, Pathology) %>%
    mutate(pair = pair_name)
}))

label_check <- label_long %>%
  group_by(individual) %>%
  summarise(
    n_labels = n_distinct(Pathology),
    Pathology = first(Pathology),
    n_pairs = n_distinct(pair),
    .groups = "drop"
  )

if (any(label_check$n_labels != 1)) {
  bad <- label_check$individual[label_check$n_labels != 1]
  stop(
    "Inconsistent pathology labels across celltype-region metadata for donor(s): ",
    paste(bad, collapse = ", ")
  )
}

if (any(label_check$n_pairs != length(pair_names))) {
  stop("Some common donors do not have metadata records in all 36 datasets.")
}

# Preserve a deterministic donor order.
common_ids <- sort(common_ids)
label_map <- setNames(label_check$Pathology, label_check$individual)
Pathology <- unname(label_map[common_ids])
y01 <- ifelse(Pathology == "lateAD", 1, 0)

if (length(unique(y01)) < 2) {
  stop("The common donor cohort does not contain both nonAD and lateAD.")
}

cat("nonAD:", sum(y01 == 0), "\n")
cat("lateAD:", sum(y01 == 1), "\n\n")

# ============================================================
# 3. BUILD THE STACKED OVERALL FEATURE MATRICES
# ============================================================

expr_blocks <- list()
act_blocks <- list()
feature_map_list <- list()
feature_summary_list <- list()

for (pair_name in pair_names) {

  expr <- expr_raw_list[[pair_name]][, common_ids, drop = FALSE]
  act  <- act_raw_list[[pair_name]][, common_ids, drop = FALSE]

  # Remove all-zero TFs separately within the FINAL common donor cohort.
  expr_nonzero <- rowSums(abs(expr), na.rm = TRUE) != 0
  act_nonzero  <- rowSums(abs(act),  na.rm = TRUE) != 0

  expr <- expr[expr_nonzero, , drop = FALSE]
  act  <- act[act_nonzero,  , drop = FALSE]

  # Fair comparison: exactly the same TF identities in Expression and Activity.
  shared_tfs <- intersect(rownames(expr), rownames(act))

  if (length(shared_tfs) == 0) {
    stop("No shared non-zero TFs after common-donor filtering for ", pair_name)
  }

  expr <- expr[shared_tfs, , drop = FALSE]
  act  <- act[shared_tfs,  , drop = FALSE]

  # Remove zero-variance features in EITHER representation, from BOTH.
  expr_sd <- apply(expr, 1, sd, na.rm = TRUE)
  act_sd  <- apply(act,  1, sd, na.rm = TRUE)

  variable_both <-
    is.finite(expr_sd) & expr_sd > 0 &
    is.finite(act_sd) & act_sd > 0

  shared_tfs_final <- shared_tfs[variable_both]

  if (length(shared_tfs_final) == 0) {
    stop("No shared variable TFs for ", pair_name)
  }

  expr <- expr[shared_tfs_final, , drop = FALSE]
  act  <- act[shared_tfs_final,  , drop = FALSE]

  # Same TF in different contexts becomes a distinct feature.
  feature_names <- paste(pair_name, shared_tfs_final, sep = "_")

  rownames(expr) <- feature_names
  rownames(act)  <- feature_names

  expr_blocks[[pair_name]] <- expr
  act_blocks[[pair_name]]  <- act

  feature_map_list[[pair_name]] <- data.frame(
    overall_feature = feature_names,
    pair = pair_name,
    TF = shared_tfs_final,
    Expression_SD = expr_sd[shared_tfs_final],
    Activity_SD = act_sd[shared_tfs_final],
    stringsAsFactors = FALSE
  )

  feature_summary_list[[pair_name]] <- data.frame(
    pair = pair_name,
    n_expression_nonzero = sum(expr_nonzero),
    n_activity_nonzero = sum(act_nonzero),
    n_shared_nonzero_TFs = length(shared_tfs),
    n_shared_variable_TFs_used = length(shared_tfs_final)
  )

  cat(
    pair_name,
    ": final shared features = ", length(shared_tfs_final),
    "\n",
    sep = ""
  )
}

expr_overall <- do.call(rbind, expr_blocks)
act_overall  <- do.call(rbind, act_blocks)
feature_map  <- bind_rows(feature_map_list)
feature_summary <- bind_rows(feature_summary_list)

if (!identical(rownames(expr_overall), rownames(act_overall))) {
  stop("Internal error: overall Expression and Activity feature IDs differ.")
}
if (!identical(colnames(expr_overall), colnames(act_overall))) {
  stop("Internal error: overall Expression and Activity donor orders differ.")
}
if (any(!is.finite(expr_overall)) || any(!is.finite(act_overall))) {
  stop("NA/Inf detected in the final overall feature matrices.")
}

x_expr <- t(expr_overall)
x_act  <- t(act_overall)

cat("\n============================================================\n")
cat("OVERALL FEATURE MATRICES\n")
cat("============================================================\n")
cat("Individuals:", nrow(x_expr), "\n")
cat("Features per representation:", ncol(x_expr), "\n")
cat("Expression dimensions:", paste(dim(x_expr), collapse = " x "), "\n")
cat("Activity dimensions:  ", paste(dim(x_act), collapse = " x "), "\n\n")

# Save construction diagnostics before modeling.
write.table(
  data.frame(
    individual = common_ids,
    Pathology = Pathology
  ),
  file.path(output_dir, "OVERALL_common_individuals.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  feature_map,
  file.path(output_dir, "OVERALL_feature_map.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  feature_summary,
  file.path(output_dir, "OVERALL_feature_count_by_celltype_region.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ============================================================
# 4. PRE-TUNE LAMBDA USING POOLED OOF AUPR
# ============================================================

# Use identical tuning folds for Expression and Activity.
# This is explicit LOOCV during tuning: each donor is held out once.
foldid_tune <- make_loocv_foldid(length(common_ids))

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

cat("Expression best alpha/lambda/tuning AUPR:",
    expr_tune$best_alpha, expr_tune$best_lambda,
    round(expr_tune$best_aupr, 4), "\n")
cat("Activity best alpha/lambda/tuning AUPR:  ",
    act_tune$best_alpha, act_tune$best_lambda,
    round(act_tune$best_aupr, 4), "\n\n")

expr_tuning_table <- expr_tune$table %>%
  mutate(Data = "Expression", .before = 1)
act_tuning_table <- act_tune$table %>%
  mutate(Data = "Activity", .before = 1)

tuning_df <- bind_rows(expr_tuning_table, act_tuning_table)

write.table(
  tuning_df,
  file.path(output_dir, "OVERALL_parameter_tuning.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ============================================================
# 5. FIX PARAMETERS, THEN FINAL LOOCV
# ============================================================

n <- length(common_ids)
expr_prob <- rep(NA_real_, n)
act_prob  <- rep(NA_real_, n)

loocv_diag <- vector("list", n * 2)
diag_i <- 0

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

  coef_expr <- as.matrix(coef(fit_expr, s = expr_tune$best_lambda))
  beta_expr <- coef_expr[-1, 1]

  diag_i <- diag_i + 1
  loocv_diag[[diag_i]] <- data.frame(
    individual = common_ids[i],
    Pathology = Pathology[i],
    Data = "Expression",
    alpha = expr_tune$best_alpha,
    lambda = expr_tune$best_lambda,
    n_nonzero_features = sum(beta_expr != 0),
    nonzero_fraction = mean(beta_expr != 0),
    coefficient_L1 = sum(abs(beta_expr)),
    coefficient_L2 = sqrt(sum(beta_expr^2))
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

  coef_act <- as.matrix(coef(fit_act, s = act_tune$best_lambda))
  beta_act <- coef_act[-1, 1]

  diag_i <- diag_i + 1
  loocv_diag[[diag_i]] <- data.frame(
    individual = common_ids[i],
    Pathology = Pathology[i],
    Data = "Activity",
    alpha = act_tune$best_alpha,
    lambda = act_tune$best_lambda,
    n_nonzero_features = sum(beta_act != 0),
    nonzero_fraction = mean(beta_act != 0),
    coefficient_L1 = sum(abs(beta_act)),
    coefficient_L2 = sqrt(sum(beta_act^2))
  )
}

loocv_diag_df <- bind_rows(loocv_diag)

prediction_df <- data.frame(
  individual = common_ids,
  Pathology = Pathology,
  y01 = y01,
  Expression_prob = expr_prob,
  Activity_prob = act_prob
)

write.table(
  prediction_df,
  file.path(output_dir, "OVERALL_LOOCV_predictions.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  loocv_diag_df,
  file.path(output_dir, "OVERALL_LOOCV_model_diagnostics.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ============================================================
# 6. FINAL PR / AUPR AND ROC / AUC COMPARISON
# ============================================================

pr_expr <- compute_pr_metrics(y01, expr_prob)
pr_act  <- compute_pr_metrics(y01, act_prob)

aupr_expr <- pr_expr$aupr
aupr_act  <- pr_act$aupr
delta_aupr <- aupr_act - aupr_expr

roc_expr <- pROC::roc(
  response = y01,
  predictor = expr_prob,
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

roc_act <- pROC::roc(
  response = y01,
  predictor = act_prob,
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

auc_expr <- as.numeric(pROC::auc(roc_expr))
auc_act  <- as.numeric(pROC::auc(roc_act))
delta_auc <- auc_act - auc_expr

ci_auc_expr <- as.numeric(pROC::ci.auc(roc_expr))
ci_auc_act  <- as.numeric(pROC::ci.auc(roc_act))

p_delong <- tryCatch(
  pROC::roc.test(
    roc_expr,
    roc_act,
    paired = TRUE,
    method = "delong",
    alternative = "two.sided"
  )$p.value,
  error = function(e) NA_real_
)

# AUPR inference on the paired final OOF predictions.
boot_aupr <- bootstrap_metric_difference(
  y01,
  expr_prob,
  act_prob,
  metric = "AUPR",
  B = n_boot,
  seed = seed + 501
)

perm_aupr <- paired_permutation_aupr(
  y01,
  expr_prob,
  act_prob,
  B = n_perm,
  seed = seed + 502
)

# Also provide a bootstrap CI for delta AUC for descriptive symmetry.
boot_auc <- bootstrap_metric_difference(
  y01,
  expr_prob,
  act_prob,
  metric = "AUC",
  B = n_boot,
  seed = seed + 503
)

cat("============================================================\n")
cat("FINAL OVERALL CLASSIFICATION\n")
cat("============================================================\n")
cat("Expression AUPR:", round(aupr_expr, 4), "\n")
cat("Activity AUPR:  ", round(aupr_act, 4), "\n")
cat("Delta AUPR (Act-Expr):", round(delta_aupr, 4), "\n")
cat("Delta AUPR bootstrap 95% CI:",
    round(boot_aupr$lower, 4), "to", round(boot_aupr$upper, 4), "\n")
cat("Paired permutation p (Activity AUPR > Expression):",
    format_p(perm_aupr$p_greater), "\n\n")
cat("Expression AUC:", round(auc_expr, 4), "\n")
cat("Activity AUC:  ", round(auc_act, 4), "\n")
cat("Delta AUC (Act-Expr):", round(delta_auc, 4), "\n")
cat("Paired DeLong p:", format_p(p_delong), "\n")

# ============================================================
# 7. FIT FULL-DATA REFERENCE MODELS FOR FEATURE INTERPRETATION ONLY
# ============================================================

# These models are NOT used to calculate final performance.
# They are saved only to identify which celltype-region-TF features receive
# non-zero coefficients in a single reference fit using the selected settings.
full_expr_fit <- glmnet(
  x = x_expr,
  y = y01,
  family = "binomial",
  alpha = expr_tune$best_alpha,
  lambda = expr_tune$best_lambda,
  standardize = TRUE
)

full_act_fit <- glmnet(
  x = x_act,
  y = y01,
  family = "binomial",
  alpha = act_tune$best_alpha,
  lambda = act_tune$best_lambda,
  standardize = TRUE
)

expr_coef_mat <- as.matrix(coef(full_expr_fit, s = expr_tune$best_lambda))
act_coef_mat  <- as.matrix(coef(full_act_fit, s = act_tune$best_lambda))

expr_coef_df <- data.frame(
  overall_feature = rownames(expr_coef_mat)[-1],
  Expression_beta = expr_coef_mat[-1, 1],
  stringsAsFactors = FALSE
)

act_coef_df <- data.frame(
  overall_feature = rownames(act_coef_mat)[-1],
  Activity_beta = act_coef_mat[-1, 1],
  stringsAsFactors = FALSE
)

coef_df <- feature_map %>%
  left_join(expr_coef_df, by = "overall_feature") %>%
  left_join(act_coef_df, by = "overall_feature") %>%
  mutate(
    Expression_abs_beta = abs(Expression_beta),
    Activity_abs_beta = abs(Activity_beta),
    Expression_nonzero = Expression_beta != 0,
    Activity_nonzero = Activity_beta != 0
  )

write.table(
  coef_df,
  file.path(output_dir, "OVERALL_full_model_feature_coefficients.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

context_coef_summary <- coef_df %>%
  group_by(pair) %>%
  summarise(
    Expression_nonzero_features = sum(Expression_nonzero, na.rm = TRUE),
    Activity_nonzero_features = sum(Activity_nonzero, na.rm = TRUE),
    Expression_sum_abs_beta = sum(Expression_abs_beta, na.rm = TRUE),
    Activity_sum_abs_beta = sum(Activity_abs_beta, na.rm = TRUE),
    .groups = "drop"
  )

write.table(
  context_coef_summary,
  file.path(output_dir, "OVERALL_full_model_context_contribution.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ============================================================
# 8. SAVE COMPACT SUMMARY
# ============================================================

overall_summary <- data.frame(
  n_individuals = length(common_ids),
  n_nonAD = sum(y01 == 0),
  n_lateAD = sum(y01 == 1),
  n_features_per_representation = ncol(x_expr),
  alpha = alpha_grid[1],
  Expression_best_lambda = expr_tune$best_lambda,
  Activity_best_lambda = act_tune$best_lambda,
  Expression_tuning_pooled_OOF_AUPR = expr_tune$best_aupr,
  Activity_tuning_pooled_OOF_AUPR = act_tune$best_aupr,
  Expression_AUPR = aupr_expr,
  Activity_AUPR = aupr_act,
  Activity_minus_Expression_AUPR = delta_aupr,
  Delta_AUPR_bootstrap_CI_lower = boot_aupr$lower,
  Delta_AUPR_bootstrap_CI_upper = boot_aupr$upper,
  AUPR_paired_permutation_p_Activity_greater = perm_aupr$p_greater,
  Expression_AUC = auc_expr,
  Expression_AUC_CI_lower = ci_auc_expr[1],
  Expression_AUC_CI_upper = ci_auc_expr[3],
  Activity_AUC = auc_act,
  Activity_AUC_CI_lower = ci_auc_act[1],
  Activity_AUC_CI_upper = ci_auc_act[3],
  Activity_minus_Expression_AUC = delta_auc,
  Delta_AUC_bootstrap_CI_lower = boot_auc$lower,
  Delta_AUC_bootstrap_CI_upper = boot_auc$upper,
  paired_DeLong_p_two_sided = p_delong,
  lateAD_prevalence = mean(y01 == 1)
)

write.table(
  overall_summary,
  file.path(output_dir, "OVERALL_classification_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ============================================================
# 9. FIGURE 1: OVERALL PRECISION-RECALL CURVE
# ============================================================

pr_plot_df <- bind_rows(
  pr_expr$curve %>% mutate(Data = "Expression"),
  pr_act$curve  %>% mutate(Data = "Activity")
) %>%
  mutate(Data = factor(Data, levels = c("Expression", "Activity")))

pr_annotation <- paste0(
  "Expression AUPR = ", sprintf("%.3f", aupr_expr), "\n",
  "Activity AUPR = ", sprintf("%.3f", aupr_act), "\n",
  "Delta = ", sprintf("%+.3f", delta_aupr), "\n",
  "95% CI Delta = [", sprintf("%.3f", boot_aupr$lower), ", ",
  sprintf("%.3f", boot_aupr$upper), "]\n",
  "Permutation p = ", format_p(perm_aupr$p_greater)
)

p_pr <- ggplot(
  pr_plot_df,
  aes(x = Recall, y = Precision, color = Data, linetype = Data)
) +
  geom_hline(
    yintercept = mean(y01 == 1),
    linetype = "dashed",
    linewidth = 0.45
  ) +
  geom_step(linewidth = 1.1) +
  annotate(
    "text",
    x = 0.97,
    y = 0.05,
    hjust = 1,
    vjust = 0,
    label = pr_annotation,
    size = 4
  ) +
  scale_color_manual(values = c("Expression" = expr_color, "Activity" = act_color)) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(
    title = "Overall TF expression vs TF activity classification",
    subtitle = paste0(
      "Common donors across all celltype-region datasets; ",
      ncol(x_expr), " matched context-specific TF features"
    ),
    x = "Recall",
    y = "Precision",
    color = NULL,
    linetype = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p_pr)

ggsave(
  file.path(output_dir, "FIG1_OVERALL_PR_AUPR.pdf"),
  p_pr, width = 6.3, height = 6.0, units = "in"
)
ggsave(
  file.path(output_dir, "FIG1_OVERALL_PR_AUPR.png"),
  p_pr, width = 6.3, height = 6.0, units = "in", dpi = 300
)

# ============================================================
# 10. FIGURE 2: OVERALL ROC CURVE
# ============================================================

roc_plot_df <- bind_rows(
  data.frame(
    FPR = 1 - roc_expr$specificities,
    TPR = roc_expr$sensitivities,
    Data = "Expression"
  ),
  data.frame(
    FPR = 1 - roc_act$specificities,
    TPR = roc_act$sensitivities,
    Data = "Activity"
  )
) %>%
  mutate(Data = factor(Data, levels = c("Expression", "Activity"))) %>%
  arrange(Data, FPR, TPR)

roc_annotation <- paste0(
  "Expression AUC = ", sprintf("%.3f", auc_expr), "\n",
  "Activity AUC = ", sprintf("%.3f", auc_act), "\n",
  "Delta = ", sprintf("%+.3f", delta_auc), "\n",
  "DeLong p = ", format_p(p_delong)
)

p_roc <- ggplot(
  roc_plot_df,
  aes(x = FPR, y = TPR, color = Data, linetype = Data)
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.45
  ) +
  geom_step(linewidth = 1.1) +
  annotate(
    "text",
    x = 0.97,
    y = 0.05,
    hjust = 1,
    vjust = 0,
    label = roc_annotation,
    size = 4
  ) +
  scale_color_manual(values = c("Expression" = expr_color, "Activity" = act_color)) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(
    title = "Overall TF expression vs TF activity classification",
    subtitle = "Fixed-parameter LOOCV after pooled-OOF-AUPR tuning",
    x = "False positive rate (1 - specificity)",
    y = "True positive rate (sensitivity)",
    color = NULL,
    linetype = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p_roc)

ggsave(
  file.path(output_dir, "FIG2_OVERALL_ROC_AUC.pdf"),
  p_roc, width = 6.3, height = 6.0, units = "in"
)
ggsave(
  file.path(output_dir, "FIG2_OVERALL_ROC_AUC.png"),
  p_roc, width = 6.3, height = 6.0, units = "in", dpi = 300
)

# ============================================================
# 11. FIGURE 3: INDIVIDUAL HELD-OUT DISEASE PROBABILITIES
# ============================================================

prob_long <- prediction_df %>%
  select(individual, Pathology, Expression_prob, Activity_prob) %>%
  pivot_longer(
    cols = c(Expression_prob, Activity_prob),
    names_to = "Data",
    values_to = "Predicted_probability"
  ) %>%
  mutate(
    Data = recode(
      Data,
      Expression_prob = "Expression",
      Activity_prob = "Activity"
    ),
    Data = factor(Data, levels = c("Expression", "Activity")),
    Pathology = factor(Pathology, levels = c("nonAD", "lateAD"))
  )

p_prob <- ggplot(
  prob_long,
  aes(x = Pathology, y = Predicted_probability, color = Data)
) +
  geom_boxplot(
    aes(group = interaction(Pathology, Data)),
    position = position_dodge(width = 0.72),
    width = 0.58,
    outlier.shape = NA,
    fill = NA,
    linewidth = 0.7
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.72
    ),
    size = 2.2,
    alpha = 0.8
  ) +
  scale_color_manual(values = c("Expression" = expr_color, "Activity" = act_color)) +
  labs(
    title = "Held-out disease probabilities from the overall classifier",
    subtitle = "Each point is one individual's LOOCV prediction",
    x = NULL,
    y = "Predicted probability of lateAD",
    color = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

print(p_prob)

ggsave(
  file.path(output_dir, "FIG3_OVERALL_LOOCV_probabilities.pdf"),
  p_prob, width = 6.3, height = 5.5, units = "in"
)
ggsave(
  file.path(output_dir, "FIG3_OVERALL_LOOCV_probabilities.png"),
  p_prob, width = 6.3, height = 5.5, units = "in", dpi = 300
)

# ============================================================
# 12. FIGURE 4: WHICH CELLTYPE-REGION CONTEXTS CONTRIBUTE FEATURES?
# ============================================================

context_long <- context_coef_summary %>%
  select(pair, Expression_nonzero_features, Activity_nonzero_features) %>%
  pivot_longer(
    cols = c(Expression_nonzero_features, Activity_nonzero_features),
    names_to = "Data",
    values_to = "n_nonzero_features"
  ) %>%
  mutate(
    Data = recode(
      Data,
      Expression_nonzero_features = "Expression",
      Activity_nonzero_features = "Activity"
    ),
    Data = factor(Data, levels = c("Expression", "Activity"))
  )

pair_order <- context_coef_summary %>%
  mutate(total = Expression_nonzero_features + Activity_nonzero_features) %>%
  arrange(total) %>%
  pull(pair)

context_long$pair <- factor(context_long$pair, levels = pair_order)

p_context <- ggplot(
  context_long,
  aes(x = pair, y = n_nonzero_features, fill = Data)
) +
  geom_col(position = "dodge", width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = c("Expression" = expr_color, "Activity" = act_color)) +
  labs(
    title = "Context-specific features retained in the full-data reference model",
    subtitle = "For interpretation only; the reference model is not used for performance estimation",
    x = NULL,
    y = "Number of non-zero TF features",
    fill = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

print(p_context)

ggsave(
  file.path(output_dir, "FIG4_OVERALL_context_feature_contribution.pdf"),
  p_context, width = 7.5, height = 9.0, units = "in"
)
ggsave(
  file.path(output_dir, "FIG4_OVERALL_context_feature_contribution.png"),
  p_context, width = 7.5, height = 9.0, units = "in", dpi = 300
)

# ============================================================
# 13. FINISH
# ============================================================

cat("\n============================================================\n")
cat("DONE: OVERALL ALL-CONTEXT CLASSIFIER\n")
cat("============================================================\n")
cat("Common individuals:", length(common_ids), "\n")
cat("Features per representation:", ncol(x_expr), "\n")
cat("Output directory:\n", output_dir, "\n\n")
cat("Main outputs:\n")
cat("  OVERALL_common_individuals.tsv\n")
cat("  OVERALL_feature_map.tsv\n")
cat("  OVERALL_feature_count_by_celltype_region.tsv\n")
cat("  OVERALL_parameter_tuning.tsv\n")
cat("  OVERALL_LOOCV_predictions.tsv\n")
cat("  OVERALL_LOOCV_model_diagnostics.tsv\n")
cat("  OVERALL_full_model_feature_coefficients.tsv\n")
cat("  OVERALL_full_model_context_contribution.tsv\n")
cat("  OVERALL_classification_summary.tsv\n")
cat("  FIG1_OVERALL_PR_AUPR.pdf\n")
cat("  FIG2_OVERALL_ROC_AUC.pdf\n")
cat("  FIG3_OVERALL_LOOCV_probabilities.pdf\n")
cat("  FIG4_OVERALL_context_feature_contribution.pdf\n")
cat("============================================================\n")
