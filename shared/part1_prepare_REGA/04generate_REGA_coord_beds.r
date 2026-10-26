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
## Generate gene/peak BED files for region-wise REGA block builder
##
## Input:
##   03REGA_input_region_consistent/{ct}_{rg}.re_tg_matrix.h5
##   05celltype_peaks/05RE_TG_matrixs/{ct}_RE_TG_matrix.bed
##
## Output:
##   03REGA_input_region_consistent/coords/{ct}_{rg}.expr_genes_300k.bed
##   03REGA_input_region_consistent/coords/{ct}_{rg}.PeakCalls_named.bed
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(rhdf5)
})

options(stringsAsFactors = FALSE)

## ============================================================
## Settings
## ============================================================

celltypes <- c("Ast", "Exc", "OPC", "Oli", "Inh", "Mic_Immune")
regions   <- c("PFC", "TH", "AG", "EC", "HC", "MTC")

input_root <- file.path(shared_root, "part1_prepare_REGA/03REGA_input_region_consistent")

re_tg_dir <- file.path(separate_root, "01_prepare/data/re_tg_matrices/")

coord_dir <- file.path(input_root, "coords")
dir.create(coord_dir, recursive = TRUE, showWarnings = FALSE)

summary_list <- list()

## ============================================================
## Helper function
## ============================================================

read_h5_names <- function(h5_path) {
  genes <- as.character(h5read(h5_path, "rownames"))
  peaks <- as.character(h5read(h5_path, "colnames"))
  list(genes = genes, peaks = peaks)
}

## ============================================================
## Main loop
## ============================================================

for (ct in celltypes) {

  cat("\n==============================\n")
  cat("Celltype:", ct, "\n")
  cat("==============================\n")

  rt_path <- file.path(re_tg_dir, paste0(ct, "_RE_TG_matrix.bed"))

  if (!file.exists(rt_path)) {
    warning("Missing RE-TG bed file: ", rt_path)
    next
  }

  ## ------------------------------------------------------------
  ## Read original RE-TG bed
  ## ------------------------------------------------------------

  rt <- read.table(
    rt_path,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
  )

  colnames(rt) <- c(
    "g_chr", "g_start", "g_end", "gene",
    "p_chr", "p_start", "p_end", "peak"
  )

  rt$gene <- as.character(rt$gene)
  rt$peak <- as.character(rt$peak)

  ## ------------------------------------------------------------
  ## Process each region
  ## ------------------------------------------------------------

  for (rg in regions) {

    key <- paste0(ct, "_", rg)

    h5_path <- file.path(input_root, paste0(key, ".re_tg_matrix.h5"))

    if (!file.exists(h5_path)) {
      cat("Missing H5:", h5_path, "-> skipping\n")
      next
    }

    cat("Processing:", key, "\n")

    ## Read final genes and peaks from generated RE-TG matrix
    h5_names <- read_h5_names(h5_path)

    final_genes <- h5_names$genes
    final_peaks <- h5_names$peaks

    ## ------------------------------------------------------------
    ## Gene BED
    ## If duplicated gene coordinates exist, keep the longest interval.
    ## ------------------------------------------------------------

    gene_bed <- rt %>%
      filter(gene %in% final_genes) %>%
      select(g_chr, g_start, g_end, gene) %>%
      distinct() %>%
      mutate(width = g_end - g_start) %>%
      arrange(gene, desc(width)) %>%
      group_by(gene) %>%
      slice(1) %>%
      ungroup() %>%
      select(g_chr, g_start, g_end, gene)

    ## Reorder gene_bed to match H5 rownames
    gene_bed <- gene_bed[match(final_genes, gene_bed$gene), ]

    ## ------------------------------------------------------------
    ## Peak BED
    ## If duplicated peak coordinates exist, keep the longest interval.
    ## ------------------------------------------------------------

    peak_bed <- rt %>%
      filter(peak %in% final_peaks) %>%
      select(p_chr, p_start, p_end, peak) %>%
      distinct() %>%
      mutate(width = p_end - p_start) %>%
      arrange(peak, desc(width)) %>%
      group_by(peak) %>%
      slice(1) %>%
      ungroup() %>%
      select(p_chr, p_start, p_end, peak)

    ## Reorder peak_bed to match H5 colnames
    peak_bed <- peak_bed[match(final_peaks, peak_bed$peak), ]

    ## ------------------------------------------------------------
    ## Strict checks
    ## ------------------------------------------------------------

    if (any(is.na(gene_bed$gene))) {
      missing_genes <- final_genes[is.na(gene_bed$gene)]
      stop(
        key, ": some H5 genes are missing from RE-TG bed. Examples: ",
        paste(head(missing_genes, 10), collapse = ", ")
      )
    }

    if (any(is.na(peak_bed$peak))) {
      missing_peaks <- final_peaks[is.na(peak_bed$peak)]
      stop(
        key, ": some H5 peaks are missing from RE-TG bed. Examples: ",
        paste(head(missing_peaks, 10), collapse = ", ")
      )
    }

    if (!identical(gene_bed$gene, final_genes)) {
      stop(key, ": gene BED order does not match H5 rownames.")
    }

    if (!identical(peak_bed$peak, final_peaks)) {
      stop(key, ": peak BED order does not match H5 colnames.")
    }

    ## ------------------------------------------------------------
    ## Output paths
    ## ------------------------------------------------------------

    gene_bed_out <- file.path(coord_dir, paste0(key, ".expr_genes_300k.bed"))
    peak_bed_out <- file.path(coord_dir, paste0(key, ".PeakCalls_named.bed"))

    write.table(
      gene_bed,
      file = gene_bed_out,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )

    write.table(
      peak_bed,
      file = peak_bed_out,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )

    cat("  Gene BED:", gene_bed_out, "n =", nrow(gene_bed), "\n")
    cat("  Peak BED:", peak_bed_out, "n =", nrow(peak_bed), "\n")

    summary_list[[key]] <- data.frame(
      key = key,
      celltype = ct,
      region = rg,
      n_genes = nrow(gene_bed),
      n_peaks = nrow(peak_bed),
      gene_bed = gene_bed_out,
      peak_bed = peak_bed_out,
      stringsAsFactors = FALSE
    )
  }
}

## ============================================================
## Save summary
## ============================================================

if (length(summary_list) > 0) {
  summary_df <- bind_rows(summary_list)

  summary_out <- file.path(coord_dir, "coord_bed_summary.tsv")

  write.table(
    summary_df,
    file = summary_out,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat("\nSaved summary to:", summary_out, "\n")
}

cat("\nDone generating coordinate BED files.\n")
