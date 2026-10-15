#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(tibble)
})

# ============================================================
# Config
# ============================================================

BASE_DIR <- "/project/zduren/durenlab/naqing_palmetto/REGA_separate/"



INPUT_ROOT <- Sys.getenv("REGA_MATRICES_ROOT", unset = file.path(BASE_DIR, "03_process_rega/matrices"))

OUT_ROOT <- Sys.getenv("REGA_ASSIGNMENTS_ROOT", unset = file.path(BASE_DIR, "03_process_rega/module_assignments"))

TOPK_FRAC <- 0.1
P_ABS <- 0.10
TIES_METHOD <- "first"

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# Assignment functions
# ============================================================

normalize_L1L1 <- function(M) {
  col_L1 <- colSums(M)
  col_L1[!is.finite(col_L1) | col_L1 == 0] <- 1

  M1 <- sweep(M, 2, col_L1, "/")

  row_L1 <- rowSums(M1)
  row_L1[!is.finite(row_L1) | row_L1 == 0] <- 1

  sweep(M1, 1, row_L1, "/")
}


assign_topk_sets <- function(M_in, K = 50) {
  mods <- colnames(M_in)

  if (K <= 0) {
    return(setNames(vector("list", length(mods)), mods))
  }

  m <- as.matrix(M_in)
  m[!is.finite(m)] <- -Inf

  out <- setNames(vector("list", ncol(m)), mods)
  rn <- rownames(m)

  for (j in seq_len(ncol(m))) {
    v <- m[, j]
    idx <- which(is.finite(v))

    if (length(idx) == 0) {
      out[[j]] <- character(0)
    } else {
      ord <- order(v[idx], decreasing = TRUE)
      out[[j]] <- rn[head(idx[ord], K)]
    }
  }

  out
}


harden_sets_by_score <- function(M_score, set_list) {
  mods <- colnames(M_score)

  for (mn in mods) {
    if (is.null(set_list[[mn]])) {
      set_list[[mn]] <- character(0)
    }
  }

  set_list <- set_list[mods]

  all_entities <- unique(unlist(set_list, use.names = FALSE))
  entities <- intersect(all_entities[!is.na(all_entities)], rownames(M_score))

  out <- setNames(vector("list", length(mods)), mods)

  if (length(entities) == 0) {
    return(out)
  }

  for (g in entities) {
    cand_mods <- names(Filter(function(v) g %in% v, set_list))

    if (length(cand_mods) == 0) {
      next
    }

    sc <- M_score[g, cand_mods, drop = TRUE]
    sc[!is.finite(sc)] <- -Inf

    if (all(sc == -Inf)) {
      next
    }

    best_mod <- cand_mods[which.max(sc)]
    out[[best_mod]] <- c(out[[best_mod]], g)
  }

  lapply(out, unique)
}


assign_label_L1L1 <- function(M, p = 0.25, ties.method = "first",
                              mode = c("quantile", "absolute")) {
  mode <- match.arg(mode)

  M2 <- normalize_L1L1(M)

  S <- M2
  S[!is.finite(S)] <- -Inf

  idx <- max.col(S, ties.method = ties.method)
  d <- S[cbind(seq_len(nrow(S)), idx)]

  thr <- if (mode == "quantile") {
    stats::quantile(d[is.finite(d)], probs = p, names = FALSE, type = 7)
  } else {
    p
  }

  lab <- integer(nrow(S))
  keep <- which(is.finite(d) & d >= thr)
  lab[keep] <- idx[keep]

  list(labels = lab, M_norm = M2, thr = thr)
}


build_union_module_sets <- function(M, topk = 50, p = 0.25,
                                    ties.method = "first") {
  stopifnot(is.matrix(M), p >= 0, p <= 1)

  if (is.null(colnames(M))) {
    stop("Matrix columns must have module names")
  }

  if (is.null(rownames(M))) {
    stop("Matrix rows must have entity names")
  }

  M_norm <- normalize_L1L1(M)

  module_names <- colnames(M_norm)
  entity_names <- rownames(M_norm)

  if (topk > 0) {
    topk_sets <- assign_topk_sets(M_norm, K = topk)
    topk_hardened <- harden_sets_by_score(M_norm, topk_sets)
  } else {
    topk_hardened <- setNames(vector("list", length(module_names)), module_names)
  }

  lab_res <- assign_label_L1L1(
    M,
    p = p,
    ties.method = ties.method,
    mode = "absolute"
  )

  labels <- lab_res$labels

  label_sets <- setNames(vector("list", ncol(M_norm)), module_names)

  for (k in seq_along(module_names)) {
    entities_k <- entity_names[which(labels == k)]
    label_sets[[k]] <- entities_k
  }

  union_sets <- setNames(vector("list", length(module_names)), module_names)

  for (mn in module_names) {
    g1 <- topk_hardened[[mn]]
    if (is.null(g1)) g1 <- character(0)

    g2 <- label_sets[[mn]]
    if (is.null(g2)) g2 <- character(0)

    union_sets[[mn]] <- unique(c(g1, g2))
  }

  list(
    union_sets = union_sets,
    topk_hardened = topk_hardened,
    label_sets = label_sets,
    M_norm = M_norm,
    params = list(topk = topk, p = p, mode = "absolute")
  )
}


assign_entities_parallel_union_absolute <- function(M, topk = 50, p = 0.25,
                                                    ties.method = "first") {
  res_union <- build_union_module_sets(
    M,
    topk = topk,
    p = p,
    ties.method = ties.method
  )

  hardened_sets <- harden_sets_by_score(
    res_union$M_norm,
    res_union$union_sets
  )

  list(
    hardened_sets = hardened_sets,
    union_sets = res_union$union_sets,
    topk_hardened = res_union$topk_hardened,
    label_sets = res_union$label_sets,
    M_norm = res_union$M_norm,
    params = res_union$params
  )
}


# ============================================================
# I/O helpers
# ============================================================

read_matrix_with_names_tsv <- function(path) {
  df <- read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (ncol(df) < 2) {
    stop("Input file must have at least 2 columns: ", path)
  }

  rn <- as.character(df[[1]])

  mat <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"

  rownames(mat) <- rn
  colnames(mat) <- colnames(df)[-1]

  mat
}


write_sets_and_assignments <- function(sets, M, prefix, out_dir) {
  out_sets_file <- file.path(out_dir, sprintf("%s_sets.txt", prefix))

  con <- file(out_sets_file, "wt")
  combined_lines_local <- list()

  counts <- integer(length(sets))
  names(counts) <- names(sets)

  for (mn in names(sets)) {
    ents <- sets[[mn]]

    if (length(ents) == 0) {
      next
    }

    set_name <- sprintf("%s_%s", prefix, mn)
    line <- paste(c(set_name, ents), collapse = "\t")

    writeLines(line, con)

    combined_lines_local[[length(combined_lines_local) + 1]] <- line
    counts[mn] <- length(ents)
  }

  close(con)

  all_entities <- rownames(M)
  ent2module <- setNames(rep(NA_character_, length(all_entities)), all_entities)

  for (mn in names(sets)) {
    ents <- sets[[mn]]

    if (length(ents) == 0) {
      next
    }

    ent2module[ents] <- mn
  }

  df_assign <- tibble(
    entity = names(ent2module),
    assigned_module = as.character(ent2module)
  )

  out_assign_file <- file.path(
    out_dir,
    sprintf("%s_module_assignments.tsv", prefix)
  )

  write.table(
    df_assign,
    file = out_assign_file,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )

  counts_df <- tibble(
    module = names(counts),
    n = as.integer(counts)
  )

  out_counts_file <- file.path(
    out_dir,
    sprintf("module_%s_counts.csv", prefix)
  )

  write.csv(counts_df, out_counts_file, row.names = FALSE)

  list(
    sets_file = out_sets_file,
    assign_file = out_assign_file,
    counts_file = out_counts_file,
    combined_lines = combined_lines_local,
    counts_df = counts_df
  )
}


plot_counts_simple <- function(counts_df, title, out_png, out_pdf) {
  if (nrow(counts_df) == 0) {
    message("No counts to plot for ", title)
    return(NULL)
  }

  dfc <- counts_df %>%
    arrange(desc(n)) %>%
    mutate(module = factor(module, levels = module))

  p <- ggplot(dfc, aes(x = module, y = n)) +
    geom_col() +
    labs(title = title, x = NULL, y = "count") +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      plot.title = element_text(size = 11, face = "bold")
    )

  ggsave(out_png, p, width = 10, height = 4, dpi = 200)
  ggsave(out_pdf, p, width = 10, height = 4)

  list(png = out_png, pdf = out_pdf)
}


process_one_matrix <- function(matrix_file, entity_type, out_dir, plot_dir,
                               combined_lines_all, counts_all) {
  if (!file.exists(matrix_file)) {
    message("[SKIP] ", entity_type, " file not found: ", matrix_file)
    return(list(
      combined_lines_all = combined_lines_all,
      counts_all = counts_all
    ))
  }

  message("[READ] ", matrix_file)

  M <- read_matrix_with_names_tsv(matrix_file)

  if (is.null(colnames(M))) {
    colnames(M) <- paste0("M", seq_len(ncol(M)))
  }

  topk <- max(1L, ceiling(nrow(M) * TOPK_FRAC))

  message(
    "[RUN] ",
    entity_type,
    ": n_entities = ",
    nrow(M),
    ", n_modules = ",
    ncol(M),
    ", topk = ",
    topk
  )

  res <- assign_entities_parallel_union_absolute(
    M,
    topk = topk,
    p = P_ABS,
    ties.method = TIES_METHOD
  )

  out <- write_sets_and_assignments(
    res$hardened_sets,
    M,
    entity_type,
    out_dir
  )

  combined_lines_all <- c(combined_lines_all, out$combined_lines)
  counts_all[[entity_type]] <- out$counts_df

  plot_counts_simple(
    out$counts_df,
    paste0(entity_type, ": module counts"),
    file.path(plot_dir, paste0(entity_type, "_module_counts.png")),
    file.path(plot_dir, paste0(entity_type, "_module_counts.pdf"))
  )

  message("[DONE] ", entity_type, " processed.")

  list(
    combined_lines_all = combined_lines_all,
    counts_all = counts_all
  )
}


# ============================================================
# Discover every exported celltype-region matrix folder.
celltype_regions <- basename(list.dirs(INPUT_ROOT, recursive = FALSE, full.names = TRUE))
celltype_regions <- celltype_regions[nzchar(celltype_regions)]
celltype_regions <- celltype_regions[!startsWith(celltype_regions, ".")]
celltype_regions <- celltype_regions[!startsWith(celltype_regions, "Vasc_")]

message("Found ", length(celltype_regions), " celltype_region entries.")


# Main loop over celltype_region
# ============================================================

all_group_counts <- list()

for (celltype_region in celltype_regions) {
  message("\n============================================================")
  message("Processing ", celltype_region)
  message("============================================================")

  input_dir <- file.path(INPUT_ROOT, celltype_region)

  if (!dir.exists(input_dir)) {
    message("[SKIP] input directory not found: ", input_dir)
    next
  }

  out_dir <- file.path(OUT_ROOT, celltype_region)
  plot_dir <- file.path(out_dir, "plots")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  wbv_file <- file.path(
    input_dir,
    paste0(celltype_region, "_WBV_genes_by_module.tsv")
  )

  bv_file <- file.path(
    input_dir,
    paste0(celltype_region, "_BV_RE_by_module.tsv")
  )

  v_file <- file.path(
    input_dir,
    paste0(celltype_region, "_V_TF_by_module.tsv")
  )

  h_file <- file.path(
    input_dir,
    paste0(celltype_region, "_H_module_by_sample.tsv")
  )

  combined_lines_all <- list()
  counts_all <- list()

  # ------------------------------------------------------------
  # Genes from WBV: genes x module
  # ------------------------------------------------------------

  tmp <- process_one_matrix(
    matrix_file = wbv_file,
    entity_type = "genes",
    out_dir = out_dir,
    plot_dir = plot_dir,
    combined_lines_all = combined_lines_all,
    counts_all = counts_all
  )

  combined_lines_all <- tmp$combined_lines_all
  counts_all <- tmp$counts_all

  # ------------------------------------------------------------
  # RE from BV: RE x module
  # ------------------------------------------------------------

  tmp <- process_one_matrix(
    matrix_file = bv_file,
    entity_type = "RE",
    out_dir = out_dir,
    plot_dir = plot_dir,
    combined_lines_all = combined_lines_all,
    counts_all = counts_all
  )

  combined_lines_all <- tmp$combined_lines_all
  counts_all <- tmp$counts_all

  # ------------------------------------------------------------
  # TF from V: TF x module
  # ------------------------------------------------------------

  tmp <- process_one_matrix(
    matrix_file = v_file,
    entity_type = "TF",
    out_dir = out_dir,
    plot_dir = plot_dir,
    combined_lines_all = combined_lines_all,
    counts_all = counts_all
  )

  combined_lines_all <- tmp$combined_lines_all
  counts_all <- tmp$counts_all

  # ------------------------------------------------------------
  # Samples from H
  # H file is module x sample, so transpose to sample x module
  # ------------------------------------------------------------

  if (file.exists(h_file)) {
    message("[READ] ", h_file)

    H <- read_matrix_with_names_tsv(h_file)

    H_t <- t(H)

    topk <- max(1L, ceiling(nrow(H_t) * TOPK_FRAC))

    message(
      "[RUN] sample: n_samples = ",
      nrow(H_t),
      ", n_modules = ",
      ncol(H_t),
      ", topk = ",
      topk
    )

    res_s <- assign_entities_parallel_union_absolute(
      H_t,
      topk = topk,
      p = P_ABS,
      ties.method = TIES_METHOD
    )

    out_s <- write_sets_and_assignments(
      res_s$hardened_sets,
      H_t,
      "sample",
      out_dir
    )

    combined_lines_all <- c(combined_lines_all, out_s$combined_lines)
    counts_all$sample <- out_s$counts_df

    plot_counts_simple(
      out_s$counts_df,
      "sample: module counts",
      file.path(plot_dir, "sample_module_counts.png"),
      file.path(plot_dir, "sample_module_counts.pdf")
    )

    message("[DONE] sample processed.")
  } else {
    message("[SKIP] H file not found: ", h_file)
  }

  # ------------------------------------------------------------
  # Write combined sets for this celltype_region
  # ------------------------------------------------------------

  combined_sets_path <- file.path(out_dir, "ALL_types_sets.txt")

  if (length(combined_lines_all) > 0) {
    writeLines(unlist(combined_lines_all), combined_sets_path)
    message("[COMBINED] written: ", combined_sets_path)
  } else {
    message("[COMBINED] No sets were created for ", celltype_region)
  }

  # ------------------------------------------------------------
  # Write combined counts for this celltype_region
  # ------------------------------------------------------------

  if (length(counts_all) > 0) {
    combined_counts <- bind_rows(lapply(names(counts_all), function(nm) {
      df <- counts_all[[nm]]
      df$type <- nm
      df
    })) %>%
      select(type, module, n)

    combined_counts_path <- file.path(out_dir, "module_counts_all_types.csv")

    write.csv(
      combined_counts,
      combined_counts_path,
      row.names = FALSE
    )

    message("[SUMMARY] combined counts saved: ", combined_counts_path)

    combined_counts$celltype_region <- celltype_region
    all_group_counts[[celltype_region]] <- combined_counts
  }

  message("[DONE] Output directory: ", out_dir)
}


# ============================================================
# Global summary across all celltype_regions
# ============================================================

if (length(all_group_counts) > 0) {
  global_counts <- bind_rows(all_group_counts) %>%
    select(celltype_region, type, module, n)

  global_counts_path <- file.path(
    OUT_ROOT,
    "module_counts_all_celltype_regions.csv"
  )

  write.csv(global_counts, global_counts_path, row.names = FALSE)

  message("\n[GLOBAL SUMMARY] written: ", global_counts_path)
}

message("\n[DONE] All output directory: ", OUT_ROOT)
