suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else "."
out_root <- if (length(args) >= 2) args[[2]] else file.path(project_root, "experiments", "SPP1_driver_TF_overlap")

cell_grn_dir <- "/project/zduren/durenlab/naqing_palmetto/AD_Driver_TF/cell_type/GRN"
cell_driver_file <- "/project/zduren/durenlab/naqing_palmetto/AD_Driver_TF/cell_type/TF_driver.tsv"
cross_driver_file <- "/project/zduren/durenlab/naqing_palmetto/AD_Driver_TF/cross_celltype/TF_driver.tsv"
cross_grn_dir <- file.path(project_root, "experiments", "EXP005_results")

cell_drivers <- read.delim(cell_driver_file, check.names = FALSE, stringsAsFactors = FALSE)
cross_drivers <- read.delim(cross_driver_file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c("celltype_region", "celltype", "region", "TF", "lm_fdr_pathology") %in% names(cell_drivers)))
stopifnot(all(c("celltype_region", "celltype", "region", "TF", "lm_fdr_pathology") %in% names(cross_drivers)))

pairs <- data.frame(
  celltype = rep(c("Mic_Immune", "OPC", "Oli"), 2),
  region = c(rep("AG", 3), rep("EC", 3)),
  stringsAsFactors = FALSE
)

read_edges <- function(path, pair, source) {
  if (!file.exists(path)) stop("Missing GRN matrix: ", path)
  x <- readRDS(path)
  if (!(inherits(x, "Matrix") || is.matrix(x) || is.data.frame(x))) stop("Unsupported GRN object: ", path)
  if (is.null(rownames(x)) || is.null(colnames(x))) stop("GRN matrix lacks dimnames: ", path)
  if (!("SPP1" %in% rownames(x))) {
    return(data.frame(celltype_region = pair, celltype = sub("_(AG|EC)$", "", pair),
      region = sub(".*_", "", pair), GRN_source = source, TF = character(),
      target_gene = character(), edge_weight = numeric(), stringsAsFactors = FALSE))
  }
  vals <- as.numeric(x["SPP1", , drop = TRUE])
  names(vals) <- colnames(x)
  keep <- is.finite(vals) & vals != 0
  data.frame(celltype_region = pair, celltype = sub("_(AG|EC)$", "", pair),
    region = sub(".*_", "", pair), GRN_source = source, TF = names(vals)[keep],
    target_gene = "SPP1", edge_weight = vals[keep], stringsAsFactors = FALSE)
}

all_edges <- list()
overlaps <- list()
audits <- list()
for (i in seq_len(nrow(pairs))) {
  celltype <- pairs$celltype[[i]]
  region <- pairs$region[[i]]
  pair <- paste(celltype, region, sep = "_")
  region_dir <- file.path(out_root, region)
  dir.create(region_dir, recursive = TRUE, showWarnings = FALSE)

  grn_specs <- list(
    list(source = "cell_type_specific", path = file.path(cell_grn_dir, paste0(pair, "_TG_TF.rds")), drivers = cell_drivers),
    list(source = "cross_celltype", path = file.path(cross_grn_dir, paste0(pair, "_TG_TF.rds")), drivers = cross_drivers)
  )
  for (spec in grn_specs) {
    e <- read_edges(spec$path, pair, spec$source)
    d <- spec$drivers[spec$drivers$celltype_region == pair &
      is.finite(spec$drivers$lm_fdr_pathology) & spec$drivers$lm_fdr_pathology < 0.05, , drop = FALSE]
    d <- d[!duplicated(d$TF), , drop = FALSE]
    m <- match(e$TF, d$TF)
    e$is_pathology_driver <- !is.na(m)
    e$driver_lm_beta_pathology <- ifelse(e$is_pathology_driver, d$lm_beta_pathology[m], NA_real_)
    e$driver_lm_pval_pathology <- ifelse(e$is_pathology_driver, d$lm_pval_pathology[m], NA_real_)
    e$driver_lm_fdr_pathology <- ifelse(e$is_pathology_driver, d$lm_fdr_pathology[m], NA_real_)
    e$driver_list <- spec$source
    all_edges[[length(all_edges) + 1]] <- e
    overlaps[[length(overlaps) + 1]] <- e[e$is_pathology_driver, , drop = FALSE]
    audits[[length(audits) + 1]] <- data.frame(
      celltype_region = pair, region = region, celltype = celltype,
      GRN_source = spec$source, grn_file = spec$path,
      driver_file = if (spec$source == "cell_type_specific") cell_driver_file else cross_driver_file,
      n_nonzero_SPP1_TF_edges = nrow(e), n_significant_driver_TFs = nrow(d),
      n_overlapping_TFs = sum(e$is_pathology_driver),
      significance_rule = "lm_fdr_pathology < 0.05",
      stringsAsFactors = FALSE
    )
  }
}

write_region <- function(region) {
  dir <- file.path(out_root, region)
  edge_df <- do.call(rbind, all_edges)
  overlap_df <- do.call(rbind, overlaps)
  audit_df <- do.call(rbind, audits)
  edge_df <- edge_df[edge_df$region == region, , drop = FALSE]
  overlap_df <- overlap_df[overlap_df$region == region, , drop = FALSE]
  audit_df <- audit_df[audit_df$region == region, , drop = FALSE]
  edge_df <- edge_df[order(edge_df$celltype, edge_df$GRN_source, -abs(edge_df$edge_weight), edge_df$TF), , drop = FALSE]
  overlap_df <- overlap_df[order(overlap_df$celltype, overlap_df$GRN_source, overlap_df$driver_lm_fdr_pathology, overlap_df$TF), , drop = FALSE]
  write.table(overlap_df, file.path(dir, "SPP1_driver_TF_overlaps.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  write.table(edge_df, file.path(dir, "SPP1_GRN_TF_edges_all.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  write.table(audit_df, file.path(dir, "source_pair_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
write_region("AG")
write_region("EC")

writeLines(c(
  "SPP1 TF driver overlap: method and interpretation",
  "",
  "For each of Mic_Immune, OPC, and Oli in AG and EC, nonzero TF-to-SPP1 edges were extracted separately from the cell type-specific GRN and the EXP005 cross-celltype GRN.",
  "An overlap is a TF with a matching celltype_region and lm_fdr_pathology < 0.05 in the corresponding cell_type/TF_driver.tsv or cross_celltype/TF_driver.tsv list.",
  "The overlap table reports the GRN source, region, driver cell type, edge weight, and pathology association statistics. SPP1_GRN_TF_edges_all.tsv includes every nonzero SPP1 regulator and flags whether it overlaps the driver list.",
  "These are inferred TF-to-target regulatory links and statistical driver-list overlaps; they do not establish causal regulation of SPP1."
), file.path(out_root, "README.txt"))

cat("Wrote AG and EC SPP1 TF overlap tables under", out_root, "\n")
print(do.call(rbind, audits)[, c("celltype_region", "GRN_source", "n_nonzero_SPP1_TF_edges", "n_significant_driver_TFs", "n_overlapping_TFs")], row.names = FALSE)
