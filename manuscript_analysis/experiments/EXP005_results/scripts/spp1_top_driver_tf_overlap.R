suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else "."
out_root <- if (length(args) >= 2) args[[2]] else file.path(project_root, "experiments", "SPP1_top_driver_TF_overlap")

cell_grn_dir <- "/project/zduren/durenlab/naqing_palmetto/AD_Driver_TF/cell_type/GRN"
cell_driver_file <- "/project/zduren/durenlab/naqing_palmetto/AD_Driver_TF/cell_type/TF_driver.tsv"
cross_driver_file <- "/project/zduren/durenlab/naqing_palmetto/AD_Driver_TF/cross_celltype/TF_driver.tsv"
cross_grn_dir <- file.path(project_root, "experiments", "EXP005_results")
specific_module_dir <- "/project/zduren/durenlab/naqing_palmetto/AD_REGA/part7_TF/04tf_assignment"
cross_module_dir <- "/project/zduren/durenlab/naqing_palmetto/AD_REGA_region_wise/part4_TF/02tf_assignment_FC"
candidate_module_list_source <- "docs/project_context.md; experiments/EXP004_integrate_results_request.md"
celltype_specific_candidates <- c(
  "Exc_AG_M8", "Inh_AG_M8", "Ast_AG_M2", "Oli_AG_M6", "Vasc_AG_M11",
  "Exc_EC_M13", "Ast_EC_M17", "Oli_EC_M14", "Oli_EC_M16",
  "Mic_immune_EC_M16", "Vasc_EC_M11", "Ast_HC_M13", "OPC_HC_M13",
  "Exc_MTC_M8", "Exc_MTC_M19", "Ast_MTC_M8", "Oli_MTC_M8",
  "Oli_MTC_M13", "OPC_MTC_M8", "Mic_Immune_MTC_M8", "Exc_PFC_M1",
  "OPC_PFC_M14", "Vasc_PFC_M20"
)
cross_celltype_candidates <- c(
  "AG_M13", "EC_M14", "EC_M17", "MTC_M6", "MTC_M13", "PFC_M14", "PFC_M18"
)

cell_drivers <- read.delim(cell_driver_file, check.names = FALSE, stringsAsFactors = FALSE)
cross_drivers <- read.delim(cross_driver_file, check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- c("celltype_region", "celltype", "region", "TF",
  "lm_beta_pathology", "lm_pval_pathology", "lm_fdr_pathology",
  "mannwhitney_lateAD_vs_nonAD_pval", "mannwhitney_lateAD_vs_nonAD_fdr")
stopifnot(all(required_columns %in% names(cell_drivers)))
stopifnot(all(required_columns %in% names(cross_drivers)))
driver_numeric_columns <- c("lm_beta_pathology", "lm_pval_pathology", "lm_fdr_pathology",
  "mannwhitney_lateAD_vs_nonAD_pval", "mannwhitney_lateAD_vs_nonAD_fdr")
for (column in driver_numeric_columns) {
  cell_drivers[[column]] <- suppressWarnings(as.numeric(cell_drivers[[column]]))
  cross_drivers[[column]] <- suppressWarnings(as.numeric(cross_drivers[[column]]))
}

pairs <- data.frame(
  celltype = rep(c("Mic_Immune", "OPC", "Oli"), 2),
  region = c(rep("AG", 3), rep("EC", 3)),
  stringsAsFactors = FALSE
)

read_edges <- function(path, pair, source) {
  if (!file.exists(path)) stop("Missing GRN matrix: ", path)
  x <- readRDS(path)
  if (!(inherits(x, "Matrix") || is.matrix(x) || is.data.frame(x))) {
    stop("Unsupported GRN object: ", path)
  }
  if (is.null(rownames(x)) || is.null(colnames(x))) {
    stop("GRN matrix lacks dimnames: ", path)
  }
  if (!("SPP1" %in% rownames(x))) {
    return(data.frame(celltype_region = pair,
      celltype = sub("_(AG|EC)$", "", pair), region = sub(".*_", "", pair),
      GRN_source = source, TF = character(), target_gene = character(),
      edge_weight = numeric(), stringsAsFactors = FALSE))
  }
  vals <- as.numeric(x["SPP1", , drop = TRUE])
  names(vals) <- colnames(x)
  keep <- is.finite(vals) & vals != 0
  data.frame(celltype_region = pair, celltype = sub("_(AG|EC)$", "", pair),
    region = sub(".*_", "", pair), GRN_source = source, TF = names(vals)[keep],
    target_gene = "SPP1", edge_weight = vals[keep], stringsAsFactors = FALSE)
}

get_module_assignment <- function(pair, region, source) {
  if (source == "cell_type_specific") {
    path <- file.path(specific_module_dir, paste0(pair, "_TF_module_assignments.tsv"))
    if (!file.exists(path)) return(list(path = path, data = NULL, source = source))
    d <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    stopifnot(all(c("entity", "assigned_module") %in% names(d)))
    d <- data.frame(celltype_region = pair, TF = d$entity,
      TF_assigned_module = as.character(d$assigned_module), stringsAsFactors = FALSE)
  } else {
    path <- file.path(cross_module_dir, region,
      paste0(region, "_ALL_celltype_region_TF_module_assignments_FC.tsv"))
    if (!file.exists(path)) return(list(path = path, data = NULL, source = source))
    d <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    stopifnot(all(c("celltype_region", "tf", "module") %in% names(d)))
    d <- data.frame(celltype_region = as.character(d$celltype_region),
      TF = as.character(d$tf), TF_assigned_module = as.character(d$module),
      stringsAsFactors = FALSE)
  }
  list(path = path, data = d, source = source)
}

annotate_modules <- function(edges, assignment, pair) {
  edges$TF_assigned_module <- NA_character_
  edges$TF_module_assignment_status <- if (is.null(assignment$data)) {
    "assignment_table_missing"
  } else {
    "no_assignment_row"
  }
  if (!is.null(assignment$data) && nrow(edges)) {
    d <- assignment$data[assignment$data$celltype_region == pair, , drop = FALSE]
    for (i in seq_len(nrow(edges))) {
      modules <- unique(d$TF_assigned_module[d$TF == edges$TF[[i]]])
      modules <- sort(modules[!is.na(modules) & nzchar(modules)])
      if (length(modules)) {
        edges$TF_assigned_module[[i]] <- paste(modules, collapse = ";")
        edges$TF_module_assignment_status[[i]] <- if (length(modules) == 1) {
          "assigned"
        } else {
          "multiple_module_assignments"
        }
      }
    }
  }
  edges$TF_module_assignment_source <- assignment$source
  edges$TF_module_assignment_file <- assignment$path
  edges
}

annotate_candidate_modules <- function(edges, pair, source) {
  candidate_set <- if (source == "cell_type_specific") {
    celltype_specific_candidates
  } else {
    cross_celltype_candidates
  }
  module_prefix <- if (source == "cell_type_specific") {
    sub("^Mic_Immune", "Mic_immune", pair)
  } else {
    sub("^.*_", "", pair)
  }
  edges$assigned_module_identifier <- NA_character_
  edges$is_disease_related_candidate_module <- NA
  edges$disease_related_candidate_module_set <- if (source == "cell_type_specific") {
    "cell_type_specific"
  } else {
    "cross_celltype"
  }
  edges$candidate_module_list_source <- candidate_module_list_source
  for (i in seq_len(nrow(edges))) {
    if (edges$TF_module_assignment_status[[i]] == "assigned") {
      full_id <- paste0(module_prefix, "_", edges$TF_assigned_module[[i]])
      edges$assigned_module_identifier[[i]] <- full_id
      edges$is_disease_related_candidate_module[[i]] <- full_id %in% candidate_set
    }
  }
  edges
}

all_edges <- list()
overlaps <- list()
audits <- list()
for (i in seq_len(nrow(pairs))) {
  celltype <- pairs$celltype[[i]]
  region <- pairs$region[[i]]
  pair <- paste(celltype, region, sep = "_")
  dir.create(file.path(out_root, region), recursive = TRUE, showWarnings = FALSE)

  grn_specs <- list(
    list(source = "cell_type_specific",
      path = file.path(cell_grn_dir, paste0(pair, "_TG_TF.rds")),
      drivers = cell_drivers, driver_file = cell_driver_file),
    list(source = "cross_celltype",
      path = file.path(cross_grn_dir, paste0(pair, "_TG_TF.rds")),
      drivers = cross_drivers, driver_file = cross_driver_file)
  )
  for (spec in grn_specs) {
    e <- read_edges(spec$path, pair, spec$source)
    module_assignment <- get_module_assignment(pair, region, spec$source)
    e <- annotate_modules(e, module_assignment, pair)
    e <- annotate_candidate_modules(e, pair, spec$source)
    # Keep every finite, nonzero TF->SPP1 edge. Rank by signed weight for
    # ordering only; the rank does not filter the edge set.
    e <- e[order(-e$edge_weight, e$TF), , drop = FALSE]
    e$edge_weight_rank <- seq_len(nrow(e))
    e$n_nonzero_TF_edges <- nrow(e)
    e$is_nonzero_weight <- is.finite(e$edge_weight) & e$edge_weight != 0

    d <- spec$drivers[spec$drivers$celltype_region == pair, , drop = FALSE]
    d <- d[!duplicated(d$TF), , drop = FALSE]
    m <- match(e$TF, d$TF)
    has_driver_row <- !is.na(m)
    driver_lm_fdr <- d$lm_fdr_pathology[m]
    driver_mw_fdr <- d$mannwhitney_lateAD_vs_nonAD_fdr[m]
    e$is_lm_fdr_driver <- has_driver_row & is.finite(driver_lm_fdr) & driver_lm_fdr < 0.05
    e$is_pairwise_mw_fdr_driver <- has_driver_row & is.finite(driver_mw_fdr) & driver_mw_fdr < 0.05
    e$is_pathology_driver <- e$is_lm_fdr_driver | e$is_pairwise_mw_fdr_driver
    e$driver_lm_beta_pathology <- ifelse(has_driver_row,
      d$lm_beta_pathology[m], NA_real_)
    e$driver_lm_pval_pathology <- ifelse(has_driver_row,
      d$lm_pval_pathology[m], NA_real_)
    e$driver_lm_fdr_pathology <- ifelse(has_driver_row,
      d$lm_fdr_pathology[m], NA_real_)
    e$driver_mannwhitney_lateAD_vs_nonAD_pval <- ifelse(has_driver_row,
      d$mannwhitney_lateAD_vs_nonAD_pval[m], NA_real_)
    e$driver_mannwhitney_lateAD_vs_nonAD_fdr <- ifelse(has_driver_row,
      d$mannwhitney_lateAD_vs_nonAD_fdr[m], NA_real_)
    e$driver_is_lm_fdr_significant <- e$is_lm_fdr_driver
    e$driver_is_pairwise_mw_fdr_significant <- e$is_pairwise_mw_fdr_driver
    e$is_nonzero_driver_overlap <- e$is_nonzero_weight & e$is_pathology_driver
    e$driver_list <- spec$source
    all_edges[[length(all_edges) + 1]] <- e
    overlaps[[length(overlaps) + 1]] <- e[e$is_nonzero_driver_overlap, , drop = FALSE]

    audits[[length(audits) + 1]] <- data.frame(
      celltype_region = pair, region = region, celltype = celltype,
      GRN_source = spec$source, grn_file = spec$path,
      driver_file = spec$driver_file,
      TF_module_assignment_file = module_assignment$path,
      n_nonzero_SPP1_TF_edges = nrow(e),
      n_lm_fdr_driver_TFs = sum(is.finite(d$lm_fdr_pathology) & d$lm_fdr_pathology < 0.05),
      n_pairwise_mw_fdr_driver_TFs = sum(is.finite(d$mannwhitney_lateAD_vs_nonAD_fdr) &
        d$mannwhitney_lateAD_vs_nonAD_fdr < 0.05),
      n_pathology_driver_TFs = sum((is.finite(d$lm_fdr_pathology) & d$lm_fdr_pathology < 0.05) |
        (is.finite(d$mannwhitney_lateAD_vs_nonAD_fdr) & d$mannwhitney_lateAD_vs_nonAD_fdr < 0.05)),
      n_driver_TFs_with_nonzero_SPP1_edge = sum(e$is_pathology_driver),
      n_lm_fdr_driver_TFs_with_nonzero_SPP1_edge = sum(e$is_lm_fdr_driver & e$is_nonzero_weight),
      n_pairwise_mw_fdr_driver_TFs_with_nonzero_SPP1_edge = sum(e$is_pairwise_mw_fdr_driver & e$is_nonzero_weight),
      n_nonzero_driver_overlaps = sum(e$is_nonzero_driver_overlap),
      n_nonzero_overlaps_in_candidate_modules = sum(e$is_nonzero_driver_overlap &
        e$is_disease_related_candidate_module %in% TRUE),
      n_nonzero_overlaps_outside_candidate_modules = sum(e$is_nonzero_driver_overlap &
        e$is_disease_related_candidate_module %in% FALSE),
      edge_selection_rule = "all finite, nonzero TF-to-SPP1 weights; no weight-rank cutoff",
      driver_rule = "matching celltype_region and (lm_fdr_pathology < 0.05 OR mannwhitney_lateAD_vs_nonAD_fdr < 0.05)",
      stringsAsFactors = FALSE
    )
  }
}

write_region <- function(region) {
  region_dir <- file.path(out_root, region)
  edge_df <- do.call(rbind, all_edges)
  overlap_df <- do.call(rbind, overlaps)
  audit_df <- do.call(rbind, audits)
  edge_df <- edge_df[edge_df$region == region, , drop = FALSE]
  overlap_df <- overlap_df[overlap_df$region == region, , drop = FALSE]
  audit_df <- audit_df[audit_df$region == region, , drop = FALSE]
  edge_df <- edge_df[order(edge_df$celltype, edge_df$GRN_source, edge_df$edge_weight_rank), , drop = FALSE]
  overlap_df <- overlap_df[order(overlap_df$celltype, overlap_df$GRN_source,
    overlap_df$edge_weight_rank, overlap_df$TF), , drop = FALSE]
  write.table(overlap_df,
    file.path(region_dir, "SPP1_nonzero_driver_TF_overlaps.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  write.table(edge_df,
    file.path(region_dir, "SPP1_GRN_TF_edges_ranked.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  write.table(audit_df,
    file.path(region_dir, "source_pair_audit.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
write_region("AG")
write_region("EC")

writeLines(c(
  "SPP1 nonzero-weight pathology-driver TF overlap",
  "",
  "For Mic_Immune, OPC, and Oli in AG and EC, nonzero finite TF-to-SPP1 GRN edges were extracted separately from the cell type-specific and EXP005 cross-celltype GRNs.",
  "All finite, nonzero TF-to-SPP1 GRN edges were retained. Edges are sorted by signed weight descending for presentation only; there is no weight-rank cutoff.",
  "An overlap requires a nonzero TF-to-SPP1 edge and a matching celltype_region with either lm_fdr_pathology < 0.05 or mannwhitney_lateAD_vs_nonAD_fdr < 0.05.",
  "TF module assignments come from the celltype-specific TF assignment table for cell_type_specific GRNs and the region-wise FC TF assignment table for cross_celltype GRNs. Missing tables and unmatched TF rows are labeled in the assignment-status column.",
  paste0("Disease-related candidate membership is matched to the source-specific candidate lists documented in ", candidate_module_list_source, ". Cell-type-specific IDs use <celltype_region>_<module>; cross-celltype IDs use <region>_<module>. Missing or ambiguous assignments remain unevaluated (NA)."),
  "SPP1_nonzero_driver_TF_overlaps.tsv contains pathology-driver TFs with nonzero SPP1 edges, both driver-test statistics and criterion flags, assigned modules, and disease-related candidate flags. SPP1_GRN_TF_edges_ranked.tsv reports all nonzero edges, their rank for ordering, driver-list overlap status, both driver-test statistics and criterion flags, assigned module, and candidate flag. source_pair_audit.tsv reports edge, driver, and overlap counts by criterion.",
  "Legacy SPP1_top10_driver_TF_overlaps.tsv files from the previous top-10% analysis may remain in the region folders; they are superseded by SPP1_nonzero_driver_TF_overlaps.tsv.",
  "GRN edges and driver overlaps are computational inferences and associations; they do not establish causal regulation of SPP1."
), file.path(out_root, "README.txt"))

cat("Wrote AG and EC SPP1 nonzero-edge driver overlap tables under", out_root, "\n")
print(do.call(rbind, audits)[, c("celltype_region", "GRN_source",
  "n_nonzero_SPP1_TF_edges", "n_lm_fdr_driver_TFs",
  "n_pairwise_mw_fdr_driver_TFs", "n_pathology_driver_TFs",
  "n_nonzero_driver_overlaps")],
  row.names = FALSE)
