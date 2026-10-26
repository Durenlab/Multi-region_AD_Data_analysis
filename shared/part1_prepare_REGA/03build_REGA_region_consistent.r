shared_root <- Sys.getenv("REGA_SHARED_ROOT")
if (!nzchar(shared_root)) {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(file_arg)) stop("Set REGA_SHARED_ROOT when sourcing this script")
  shared_root <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
  while (!(basename(shared_root) %in% c("REGA_shared", "shared")) && dirname(shared_root) != shared_root) shared_root <- dirname(shared_root)
  if (!(basename(shared_root) %in% c("REGA_shared", "shared"))) stop("Could not find shared analysis root")
}
separate_root <- Sys.getenv("REGA_SEPARATE_ROOT")
if (!nzchar(separate_root)) stop("Set REGA_SEPARATE_ROOT")
## ============================================================
## Build REGA input matrices for region-wise consistent individuals
## Outputs:
##   {ct}_{rg}.re_tg_matrix.h5
##   {ct}_{rg}.GEX_sorted.csv
##   {ct}_{rg}.peak_motif_dense.csv
## ============================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  library(rhdf5)
})

options(stringsAsFactors = FALSE)

## ============================================================
## Settings
## ============================================================

celltypes <- c("Ast", "Exc", "OPC", "Oli", "Inh", "Mic_Immune")
regions   <- c("PFC", "TH", "AG", "EC", "HC", "MTC")

tau <- 1e5
window_size <- 300000

## Input paths
cleaned_gex_dir <- file.path(shared_root, "part1_prepare_REGA/02cleaned_region_consistent_gex")

re_tg_dir <- file.path(separate_root, "01_prepare/data/re_tg_matrices/")
motif_dir <- file.path(separate_root, "01_prepare/data/motif_scan_results/")

## Output path
out_root <- file.path(shared_root, "part1_prepare_REGA/03REGA_input_region_consistent")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

## Optional summary
summary_list <- list()


## ============================================================
## Step 0: Build global motif list across all celltypes
## This makes peak_motif_dense.csv have identical motif columns
## across all celltype-region inputs.
## ============================================================

cat("Building global motif list...\n")

global_motifs <- c()

for (ct in celltypes) {
  motif_path <- sprintf("%s/%s_peak_motif.txt", motif_dir, ct)

  if (!file.exists(motif_path)) {
    warning("Missing motif file: ", motif_path)
    next
  }

  motifs_tmp <- read.table(
    motif_path,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
  )

  global_motifs <- union(global_motifs, unique(as.character(motifs_tmp$V2)))
}

global_motifs <- sort(global_motifs)

if (length(global_motifs) == 0) {
  stop("Global motif list is empty. Check motif input files.")
}

cat("Global motif count:", length(global_motifs), "\n")


## ============================================================
## Main loop
## ============================================================

for (ct in celltypes) {

  cat("\n==============================\n")
  cat("Celltype:", ct, "\n")
  cat("==============================\n")

  ## ------------------------------------------------------------
  ## Celltype-level RE-TG and motif inputs
  ## ------------------------------------------------------------

  rt_path <- sprintf("%s/%s_RE_TG_matrix.bed", re_tg_dir, ct)
  motif_path <- sprintf("%s/%s_peak_motif.txt", motif_dir, ct)

  if (!file.exists(rt_path)) {
    warning("Missing RE-TG file: ", rt_path)
    next
  }

  if (!file.exists(motif_path)) {
    warning("Missing motif file: ", motif_path)
    next
  }

  re_tg_overlap <- read.table(
    rt_path,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
  )

  motifs <- read.table(
    motif_path,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
  )

  colnames(re_tg_overlap) <- c(
    "g_chr", "g_start", "g_end", "gene",
    "p_chr", "p_start", "p_end", "peak"
  )

  re_tg_overlap$gene <- as.character(re_tg_overlap$gene)
  re_tg_overlap$peak <- as.character(re_tg_overlap$peak)
  motifs$V1 <- as.character(motifs$V1)
  motifs$V2 <- as.character(motifs$V2)

  ## ------------------------------------------------------------
  ## Peak intersection between RE-TG and motif result
  ## ------------------------------------------------------------

  peaks_rt  <- unique(re_tg_overlap$peak)
  peaks_mtf <- unique(motifs$V1)

  filtered_elements_global <- intersect(peaks_rt, peaks_mtf)

  if (length(filtered_elements_global) == 0) {
    warning("No shared peaks between RE-TG and motif for celltype: ", ct)
    next
  }

  motifs_ct <- motifs[motifs$V1 %in% filtered_elements_global, , drop = FALSE]

  re_tg_overlap_ct <- re_tg_overlap %>%
    filter(peak %in% filtered_elements_global)

  ## ============================================================
  ## Region loop
  ## ============================================================

  for (rg in regions) {

    key <- paste0(ct, "_", rg)

    cat("\n---- Processing ", key, " ----\n", sep = "")

    ## ------------------------------------------------------------
    ## Expression input
    ## ------------------------------------------------------------
    gex_path <- sprintf(
      "%s/cleaned_%s_%s_region_consistent.tsv",
      cleaned_gex_dir,
      ct,
      rg
    )
    if (!file.exists(gex_path)) {
      cat("Missing GEX:", gex_path, "-> skipping\n")
      next
    }

    ## Output files expected by the chromosome-block builder
    re_tg_h5      <- file.path(out_root, paste0(key, ".re_tg_matrix.h5"))
    gex_outfile   <- file.path(out_root, paste0(key, ".GEX_sorted.csv"))
    motif_outfile <- file.path(out_root, paste0(key, ".peak_motif_dense.csv"))

    ## ------------------------------------------------------------
    ## Read GEX
    ## Genes x individuals
    ## ------------------------------------------------------------

    gex_csv <- read.table(
      gex_path,
      header = TRUE,
      sep = "\t",
      row.names = 1,
      check.names = FALSE
    )

    gex_mat <- as.matrix(gex_csv)
    storage.mode(gex_mat) <- "numeric"

    ## REGA uses non-negative matrix factorization logic,
    ## so negative values after limma correction are set to 0.
    gex_mat[gex_mat < 0] <- 0

    adjustedExp <- gex_mat
    gene_list <- rownames(adjustedExp)

    ## ------------------------------------------------------------
    ## Gene filtering by available RE-TG links
    ## ------------------------------------------------------------

    rt <- re_tg_overlap_ct

    genes_valid <- intersect(gene_list, unique(rt$gene))

    if (length(genes_valid) == 0) {
      warning("No valid genes after intersecting GEX and RE-TG for: ", key)
      next
    }

    rt <- rt %>%
      filter(gene %in% genes_valid)

    ## Keep GEX in the same order as genes_valid
    gex_filtered <- adjustedExp[genes_valid, , drop = FALSE]

    ## ------------------------------------------------------------
    ## RE-TG weight matrix W
    ## genes x peaks
    ## ------------------------------------------------------------

    TSS <- rt$g_end - window_size
    peak_mid <- (rt$p_start + rt$p_end) / 2
    dist <- abs(peak_mid - TSS)
    weight <- exp(-dist / tau)

    matrix_rows1 <- genes_valid
    matrix_cols1 <- unique(rt$peak)

    i <- match(rt$gene, matrix_rows1)
    j <- match(rt$peak, matrix_cols1)

    keep <- !is.na(i) & !is.na(j)

    i <- i[keep]
    j <- j[keep]
    w <- weight[keep]

    ## Aggregate duplicated gene-peak links by max weight
    if (any(duplicated(paste(i, j, sep = "_")))) {
      dfw <- data.frame(i = i, j = j, w = w) %>%
        group_by(i, j) %>%
        summarise(w = max(w), .groups = "drop")

      i <- dfw$i
      j <- dfw$j
      w <- dfw$w
    }

    re_tg_matrix <- sparseMatrix(
      i = i,
      j = j,
      x = as.numeric(w),
      dims = c(length(matrix_rows1), length(matrix_cols1)),
      dimnames = list(matrix_rows1, matrix_cols1)
    )

    ## ------------------------------------------------------------
    ## Peak-motif matrix B
    ## peaks x motifs
    ## Use global_motifs to force identical motif columns
    ## ------------------------------------------------------------

    matrix_rows_peak <- colnames(re_tg_matrix)

    motifs2 <- motifs_ct %>%
      select(V1, V2) %>%
      rename(peak = V1, motif = V2) %>%
      filter(peak %in% matrix_rows_peak)

    ii <- match(motifs2$peak, matrix_rows_peak)
    jj <- match(motifs2$motif, global_motifs)

    keep2 <- !is.na(ii) & !is.na(jj)

    ii <- ii[keep2]
    jj <- jj[keep2]

    score_matrix <- sparseMatrix(
      i = ii,
      j = jj,
      x = rep.int(1L, length(ii)),
      dims = c(length(matrix_rows_peak), length(global_motifs)),
      dimnames = list(matrix_rows_peak, global_motifs)
    )

    ## ------------------------------------------------------------
    ## Final peak alignment check
    ## ------------------------------------------------------------

    if (!identical(colnames(re_tg_matrix), rownames(score_matrix))) {
      common_peaks <- intersect(colnames(re_tg_matrix), rownames(score_matrix))

      re_tg_matrix <- re_tg_matrix[, common_peaks, drop = FALSE]
      score_matrix <- score_matrix[common_peaks, , drop = FALSE]
    }

    if (nrow(re_tg_matrix) == 0 || ncol(re_tg_matrix) == 0) {
      warning("Empty RE-TG matrix after alignment for: ", key)
      next
    }

    if (nrow(score_matrix) == 0 || ncol(score_matrix) == 0) {
      warning("Empty motif matrix after alignment for: ", key)
      next
    }

    ## Match GEX to final gene order in re_tg_matrix
    gex_filtered <- adjustedExp[rownames(re_tg_matrix), , drop = FALSE]

    ## ------------------------------------------------------------
    ## Diagnostics
    ## ------------------------------------------------------------

    cat(sprintf(
      "[%s] GEX: %d x %d | RE-TG: %d x %d | Motif: %d x %d\n",
      key,
      nrow(gex_filtered), ncol(gex_filtered),
      nrow(re_tg_matrix), ncol(re_tg_matrix),
      nrow(score_matrix), ncol(score_matrix)
    ))

    ## ------------------------------------------------------------
    ## Save RE-TG HDF5
    ## ------------------------------------------------------------

    if (file.exists(re_tg_h5)) {
      file.remove(re_tg_h5)
    }

    trip <- summary(re_tg_matrix)

    h5write(trip$i, re_tg_h5, "i")
    h5write(trip$j, re_tg_h5, "j")
    h5write(as.numeric(trip$x), re_tg_h5, "x")
    h5write(dim(re_tg_matrix), re_tg_h5, "dims")
    h5write(rownames(re_tg_matrix), re_tg_h5, "rownames")
    h5write(colnames(re_tg_matrix), re_tg_h5, "colnames")

    ## ------------------------------------------------------------
    ## Save GEX and motif matrix
    ## ------------------------------------------------------------

    write.csv(gex_filtered, file = gex_outfile, quote = FALSE)

    ## Dense output because your colleague's block builder uses pd.read_csv
    ## and expects a normal RE x motif dataframe.
    score_matrix_dense <- as.matrix(score_matrix)

    write.csv(score_matrix_dense, file = motif_outfile, quote = FALSE)

    cat("Saved outputs for:", key, "\n")

    summary_list[[key]] <- data.frame(
      key = key,
      celltype = ct,
      region = rg,
      n_genes = nrow(gex_filtered),
      n_individuals = ncol(gex_filtered),
      n_peaks = ncol(re_tg_matrix),
      n_motifs = ncol(score_matrix),
      stringsAsFactors = FALSE
    )
  }
}

## ============================================================
## Save summary
## ============================================================

if (length(summary_list) > 0) {
  summary_df <- bind_rows(summary_list)
  summary_out <- file.path(out_root, "input_summary.tsv")
  write.table(
    summary_df,
    file = summary_out,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  cat("\nSaved summary to:", summary_out, "\n")
}

cat("\nDone building shared-H REGA input data.\n")
