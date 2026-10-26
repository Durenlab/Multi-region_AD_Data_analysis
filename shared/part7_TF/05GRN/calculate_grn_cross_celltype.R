#!/usr/bin/env Rscript

# Workflow
# Input:
#   1. One cell-type/region's WBV, BV and V module matrices
#   2. Sparse W (TG-RE) and B (RE-motif) binary edge tables
#   3. Shared motif_to_TF prior and pair-specific M1/PCC matrices.
#   4. H and VH are checked as inputs/provenance; VH generated the PCC upstream.
#
# Calculation:
#   TG_RE = (WBV %*% t(BV)) restricted to W edges
#   RE_TF = (BV %*% t(V)) restricted to B edges, then multiplied by
#           weighted motif-to-TF relationships
#   TG_TF = WBV %*% t(V) %*% weighted motif-to-TF relationships
#
# Output:
#   One output directory containing three sparse RDS matrices:
#     <CELLTYPE_REGION>_TG_RE.rds   # genes x REs
#     <CELLTYPE_REGION>_RE_TF.rds   # REs x TFs
#     <CELLTYPE_REGION>_TG_TF.rds   # genes x TFs
#   and one metadata file:
#     <CELLTYPE_REGION>_GRN_metadata.rds
#
# Example:
#   Rscript calculate_grn_cross_celltype.R \
#     --celltype_region Exc_AG \
#     --output_dir part7_TF/05GRN/results

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
shared_root <- Sys.getenv("REGA_SHARED_ROOT", unset = normalizePath(file.path(script_dir, "../..")))
separate_root <- Sys.getenv("REGA_SEPARATE_ROOT")
if (!nzchar(separate_root)) stop("Set REGA_SEPARATE_ROOT")

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

## -------------------------------------------------------------------------
## Argument parsing
## -------------------------------------------------------------------------

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args)) {
    stop("Required argument: --celltype_region CELLTYPE_REGION")
  }

  result <- list(
    celltype_region = NULL,
    matrix_dir = file.path(shared_root, "part3_process_REGA_output/09extracted_matrices_region_shareH"),
    tf_dir = file.path(shared_root, "part7_TF/03TF_activity"),
    prior_file = file.path(separate_root, "external_inputs/motif_to_TF.txt"),
    output_dir = NULL,
    edge_chunk_size = 100000L,
    overwrite = FALSE
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--overwrite")) {
      result$overwrite <- TRUE
      i <- i + 1L
      next
    }
    if (i == length(args)) {
      stop("Missing value after ", key)
    }
    value <- args[[i + 1L]]
    if (key == "--celltype_region") {
      result$celltype_region <- value
    } else if (key == "--matrix_dir") {
      result$matrix_dir <- value
    } else if (key == "--tf_dir") {
      result$tf_dir <- value
    } else if (key == "--prior_file") {
      result$prior_file <- value
    } else if (key == "--output_dir") {
      result$output_dir <- value
    } else if (key == "--edge_chunk_size") {
      result$edge_chunk_size <- as.integer(value)
    } else {
      stop("Unknown argument: ", key)
    }
    i <- i + 2L
  }

  if (is.null(result$celltype_region) || !nzchar(result$celltype_region)) {
    stop("--celltype_region cannot be empty")
  }
  if (!grepl("^(Exc|Inh|Ast|Oli|OPC|Mic_Immune)_(AG|EC|HC|MTC|PFC|TH)$", result$celltype_region)) {
    stop("Expected one of the 36 non-Vasc cell type/region pairs")
  }
  if (is.na(result$edge_chunk_size) || result$edge_chunk_size < 1L) {
    stop("--edge_chunk_size must be a positive integer")
  }
  if (is.null(result$output_dir)) {
    result$output_dir <- file.path(shared_root, "part7_TF/05GRN/results")
  }
  result
}

## -------------------------------------------------------------------------
## Input readers
## -------------------------------------------------------------------------

read_feature_matrix <- function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }

  lines <- trimws(readLines(path, warn = FALSE))
  lines <- lines[nzchar(lines)]
  if (length(lines) < 2L) {
    stop("File contains no data rows: ", path)
  }

  parts <- strsplit(lines, "[[:space:]]+")
  header <- parts[[1L]]
  rows <- parts[-1L]
  row_lengths <- lengths(rows)

  if (length(unique(row_lengths)) != 1L) {
    stop("Data rows have inconsistent column numbers in ", path)
  }

  n_data_cols <- row_lengths[[1L]]
  n_header_cols <- length(header)

  if (n_data_cols == n_header_cols + 1L) {
    ids <- vapply(rows, `[[`, character(1), 1L)
    value_rows <- lapply(rows, function(x) x[-1L])
    value_names <- header
  } else if (
    n_data_cols == n_header_cols &&
    header[[1L]] %in% c(
      "gene", "peak", "RE", "feature", "entity", "RE_ID",
      "peak_id", "feature_id", "region", "motif", "TF"
    )
  ) {
    ids <- vapply(rows, `[[`, character(1), 1L)
    value_rows <- lapply(rows, function(x) x[-1L])
    value_names <- header[-1L]
  } else if (n_data_cols == n_header_cols) {
    ids <- as.character(seq_along(rows))
    value_rows <- rows
    value_names <- header
  } else {
    stop(
      "Unexpected column count in ", path,
      ". Header columns: ", n_header_cols,
      "; data columns: ", n_data_cols
    )
  }

  all_tokens <- unlist(value_rows, use.names = FALSE)
  all_tokens[all_tokens %in% c("", "NA", "NaN", "NULL", ".")] <- NA_character_
  numeric_values <- suppressWarnings(as.numeric(all_tokens))
  invalid_tokens <- unique(all_tokens[!is.na(all_tokens) & is.na(numeric_values)])
  if (length(invalid_tokens)) {
    stop(
      "Non-numeric values detected in ", path, ": ",
      paste(head(invalid_tokens, 20L), collapse = ", ")
    )
  }

  values <- matrix(
    numeric_values,
    nrow = length(value_rows),
    ncol = length(value_names),
    byrow = TRUE,
    dimnames = list(ids, value_names)
  )

  if (anyNA(values)) {
    stop("NA values detected in feature matrix: ", path)
  }
  values
}

read_edge_matrix <- function(path, row_ids, col_ids, relation_name) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }
  dat <- fread(path, sep = "\t", header = TRUE, check.names = FALSE)
  expected <- switch(relation_name,
                     W = c("gene", "RE", "value"),
                     B = c("RE", "motif", "value"),
                     stop("Unknown edge relation: ", relation_name))
  if (!identical(names(dat), expected)) {
    stop(
      relation_name, " edge file must contain columns ",
      paste(expected, collapse = ", "), ": ", path
    )
  }
  setnames(dat, expected[1:2], c("row", "col"))
  dat[, row := as.character(row)]
  dat[, col := as.character(col)]
  dat[, value := as.numeric(value)]
  if (anyNA(dat$value)) stop("NA/non-numeric edge values in ", path)

  i <- match(dat$row, row_ids)
  j <- match(dat$col, col_ids)
  if (anyNA(i) || anyNA(j)) {
    bad_rows <- unique(dat$row[is.na(i)])
    bad_cols <- unique(dat$col[is.na(j)])
    stop(
      relation_name, " edge IDs do not match matrix IDs. Missing rows: ",
      paste(head(bad_rows, 5L), collapse = ", "),
      "; missing columns: ",
      paste(head(bad_cols, 5L), collapse = ", ")
    )
  }

  result <- sparseMatrix(
    i = i,
    j = j,
    x = dat$value,
    dims = c(length(row_ids), length(col_ids)),
    dimnames = list(row_ids, col_ids)
  )
  drop0(result)
}

read_motif_tf <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)

  # Keep this identical to the original working script. read.table()
  # automatically interprets the first field as row names when the header
  # has one fewer field than the data rows.
  result <- as.matrix(read.table(path, header = TRUE, sep = ""))
  storage.mode(result) <- "numeric"
  if (anyNA(result)) stop("NA/non-numeric values in motif_to_TF: ", path)
  result
}

read_pcc <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  result <- as.matrix(read.csv(path, row.names = 1L, check.names = FALSE))
  storage.mode(result) <- "numeric"
  if (anyNA(result)) stop("NA/non-numeric values in PCC matrix: ", path)
  result
}

## -------------------------------------------------------------------------
## Efficient masked matrix products
## -------------------------------------------------------------------------

edge_table <- function(x) {
  s <- summary(x)
  data.table(
    left = rownames(x)[s$i],
    right = colnames(x)[s$j],
    left_index = s$i,
    right_index = s$j,
    edge_value = as.numeric(s$x)
  )
}

masked_dot_values <- function(left, right, edges, chunk_size) {
  if (!nrow(edges)) return(numeric())
  result <- numeric(nrow(edges))
  starts <- seq.int(1L, nrow(edges), by = chunk_size)

  for (start in starts) {
    end <- min(start + chunk_size - 1L, nrow(edges))
    left_chunk <- left[edges$left_index[start:end], , drop = FALSE]
    right_chunk <- right[edges$right_index[start:end], , drop = FALSE]
    result[start:end] <- rowSums(left_chunk * right_chunk)
    rm(left_chunk, right_chunk)
  }
  result
}

save_matrix <- function(x, path, overwrite) {
  if (file.exists(path) && !overwrite) {
    stop("Output exists; use --overwrite to replace it: ", path)
  }
  saveRDS(x, path, compress = "gzip")
}

## -------------------------------------------------------------------------
## Main workflow
## -------------------------------------------------------------------------

opt <- parse_args()
celltype_region <- opt$celltype_region
region <- sub("^.*_", "", celltype_region)
extracted_dir <- file.path(opt$matrix_dir, region)

wbv_file <- file.path(extracted_dir, paste0(celltype_region, "_WBV_genes_by_module.tsv"))
bv_file <- file.path(extracted_dir, paste0(celltype_region, "_BV_RE_by_module.tsv"))
v_file <- file.path(extracted_dir, paste0(celltype_region, "_V_TF_by_module.tsv"))
w_file <- file.path(extracted_dir, paste0(celltype_region, "_W_binary_sparse.tsv"))
b_file <- file.path(extracted_dir, paste0(celltype_region, "_B_binary_sparse.tsv"))
motif_tf_file <- opt$prior_file
tf_pair_dir <- file.path(opt$tf_dir, celltype_region, "TF_activity_pcc_pos_only")
pcc_file <- file.path(tf_pair_dir, paste0(celltype_region, "_TF_motif_pcc_M2.csv"))
m1_file <- file.path(tf_pair_dir, paste0(celltype_region, "_TF_motif_binary_M1.csv"))
vh_file <- file.path(tf_pair_dir, paste0(celltype_region, "_VH_from_V_times_H.csv"))
h_file <- file.path(extracted_dir, paste0(region, "_H_module_by_sample.tsv"))

required_files <- c(wbv_file, bv_file, v_file, w_file, b_file,
                    motif_tf_file, pcc_file, m1_file, vh_file, h_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop("Missing input files: ", paste(missing_files, collapse = ", "))

prefix <- file.path(opt$output_dir, celltype_region)
output_files <- paste0(prefix, c("_TG_RE.rds", "_RE_TF.rds", "_TG_TF.rds", "_GRN_metadata.rds"))
if (!opt$overwrite && any(file.exists(output_files))) {
  stop("Output exists; use --overwrite to replace this pair: ", paste(output_files[file.exists(output_files)], collapse = ", "))
}

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

message("[", celltype_region, "] Reading module matrices")
WBV <- read_feature_matrix(wbv_file)
BV <- read_feature_matrix(bv_file)
V <- read_feature_matrix(v_file)

if (anyDuplicated(rownames(WBV)) || anyDuplicated(rownames(BV)) || anyDuplicated(rownames(V))) {
  stop("Duplicated row IDs detected in WBV, BV, or V")
}

motif_to_TF <- read_motif_tf(motif_tf_file)
PCC <- read_pcc(pcc_file)
M1 <- read_pcc(m1_file)
if (!all(M1 %in% c(0, 1))) stop("M1 is not binary: ", m1_file)
if (!setequal(rownames(M1), rownames(PCC)) || !setequal(colnames(M1), colnames(PCC))) {
  stop("M1 and PCC TF/motif IDs differ for ", celltype_region)
}
shared_motifs <- intersect(rownames(motif_to_TF), colnames(M1))
shared_TFs <- intersect(colnames(motif_to_TF), rownames(M1))
if (!length(shared_motifs) || !length(shared_TFs)) stop("No shared prior/M1 IDs")
prior_support <- motif_to_TF[shared_motifs, shared_TFs, drop = FALSE] != 0
m1_support <- t(M1[shared_TFs, shared_motifs, drop = FALSE]) != 0
if (!identical(unname(prior_support), unname(m1_support))) {
  stop("Shared prior and pair-specific M1 motif-TF links differ for ", celltype_region)
}

## TG_RE: calculate only the gene-RE edges present in W.
modules_tg_re <- intersect(colnames(WBV), colnames(BV))
if (!length(modules_tg_re)) stop("No shared modules between WBV and BV")
WBV_tg_re <- WBV[, modules_tg_re, drop = FALSE]
BV_tg_re <- BV[, modules_tg_re, drop = FALSE]
W <- read_edge_matrix(w_file, rownames(WBV_tg_re), rownames(BV_tg_re), "W")
W_edges <- edge_table(W)

message("[", celltype_region, "] Computing TG_RE on ", nrow(W_edges), " W edges")
TG_RE_values <- masked_dot_values(
  WBV_tg_re, BV_tg_re, W_edges, opt$edge_chunk_size
) * W_edges$edge_value
TG_RE <- sparseMatrix(
  i = match(W_edges$left, rownames(WBV_tg_re)),
  j = match(W_edges$right, rownames(BV_tg_re)),
  x = TG_RE_values,
  dims = c(nrow(WBV_tg_re), nrow(BV_tg_re)),
  dimnames = list(rownames(WBV_tg_re), rownames(BV_tg_re))
)
TG_RE <- drop0(TG_RE)

## RE_TF: calculate only the RE-motif edges present in B.
modules_re_motif <- intersect(colnames(BV), colnames(V))
if (!length(modules_re_motif)) stop("No shared modules between BV and V")
BV_re_motif <- BV[, modules_re_motif, drop = FALSE]
V_re_motif <- V[, modules_re_motif, drop = FALSE]
B <- read_edge_matrix(b_file, rownames(BV_re_motif), rownames(V_re_motif), "B")
B_edges <- edge_table(B)

message("[", celltype_region, "] Computing RE_motif on ", nrow(B_edges), " B edges")
RE_motif_values <- masked_dot_values(
  BV_re_motif, V_re_motif, B_edges, opt$edge_chunk_size
) * B_edges$edge_value
RE_motif_B <- sparseMatrix(
  i = match(B_edges$left, rownames(BV_re_motif)),
  j = match(B_edges$right, rownames(V_re_motif)),
  x = RE_motif_values,
  dims = c(nrow(BV_re_motif), nrow(V_re_motif)),
  dimnames = list(rownames(BV_re_motif), rownames(V_re_motif))
)
RE_motif_B <- drop0(RE_motif_B)

motifs <- Reduce(intersect, list(
  colnames(RE_motif_B), rownames(motif_to_TF), colnames(PCC)
))
TFs <- intersect(colnames(motif_to_TF), rownames(PCC))
message("[", celltype_region, "] motif_to_TF dimensions: ",
        paste(dim(motif_to_TF), collapse = " x "))
message("[", celltype_region, "] PCC dimensions: ",
        paste(dim(PCC), collapse = " x "))
message("[", celltype_region, "] shared TFs: ", length(TFs))
if (!length(TFs)) {
  message("  motif_to_TF columns: ",
          paste(head(colnames(motif_to_TF), 10L), collapse = ", "))
  message("  PCC rows: ",
          paste(head(rownames(PCC), 10L), collapse = ", "))
}
if (!length(motifs)) stop("No shared motifs across RE_motif, motif_to_TF, and PCC")
if (!length(TFs)) stop("No shared TFs between motif_to_TF and PCC")

weighted_motif_to_TF <- Matrix(
  motif_to_TF[motifs, TFs, drop = FALSE] * t(PCC[TFs, motifs, drop = FALSE]),
  sparse = TRUE
)
dimnames(weighted_motif_to_TF) <- list(motifs, TFs)
RE_motif_B <- RE_motif_B[, motifs, drop = FALSE]

message("[", celltype_region, "] Computing RE_TF")
RE_TF <- drop0(RE_motif_B %*% weighted_motif_to_TF)

## TG_TF: use associativity to avoid constructing the large gene-motif matrix.
modules_tg_tf <- intersect(colnames(WBV), colnames(V))
if (!length(modules_tg_tf)) stop("No shared modules between WBV and V")
WBV_tg_tf <- WBV[, modules_tg_tf, drop = FALSE]
V_tg_tf <- V[, modules_tg_tf, drop = FALSE]
motifs_tg_tf <- Reduce(intersect, list(
  rownames(V_tg_tf), rownames(motif_to_TF), colnames(PCC)
))
if (!length(motifs_tg_tf)) stop("No shared motifs for TG_TF")
V_tg_tf <- V_tg_tf[motifs_tg_tf, , drop = FALSE]
M_tg_tf <- weighted_motif_to_TF[motifs_tg_tf, TFs, drop = FALSE]

message("[", celltype_region, "] Computing TG_TF")
module_to_TF <- t(V_tg_tf) %*% M_tg_tf
TG_TF <- drop0(Matrix(WBV_tg_tf, sparse = TRUE) %*% module_to_TF)

## Save sparse matrices. Dimnames are retained in each RDS file.
save_matrix(TG_RE, paste0(prefix, "_TG_RE.rds"), opt$overwrite)
save_matrix(RE_TF, paste0(prefix, "_RE_TF.rds"), opt$overwrite)
save_matrix(TG_TF, paste0(prefix, "_TG_TF.rds"), opt$overwrite)

metadata <- list(
  celltype_region = celltype_region,
  dimensions = list(
    TG_RE = dim(TG_RE),
    RE_TF = dim(RE_TF),
    TG_TF = dim(TG_TF)
  ),
  nonzero = list(
    TG_RE = nnzero(TG_RE),
    RE_TF = nnzero(RE_TF),
    TG_TF = nnzero(TG_TF)
  ),
  modules = list(
    TG_RE = modules_tg_re,
    RE_motif = modules_re_motif,
    TG_TF = modules_tg_tf
  ),
  motifs = motifs,
  TFs = TFs,
  inputs = c(
    WBV = wbv_file,
    BV = bv_file,
    V = v_file,
    W = w_file,
    B = b_file,
    motif_to_TF = motif_tf_file,
    PCC = pcc_file,
    M1 = m1_file,
    VH = vh_file,
    H = h_file
  )
)
save_matrix(metadata, paste0(prefix, "_GRN_metadata.rds"), opt$overwrite)

message("[", celltype_region, "] Done")
message("  TG_RE: ", paste(dim(TG_RE), collapse = " x "),
        "; nonzero = ", nnzero(TG_RE))
message("  RE_TF: ", paste(dim(RE_TF), collapse = " x "),
        "; nonzero = ", nnzero(RE_TF))
message("  TG_TF: ", paste(dim(TG_TF), collapse = " x "),
        "; nonzero = ", nnzero(TG_TF))
