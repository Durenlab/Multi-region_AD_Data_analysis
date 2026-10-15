shared_root <- normalizePath(Sys.getenv("REGA_SHARED_ROOT", unset = normalizePath("..")))
stage_dir <- file.path(shared_root, "part6_Age_related")
gsea_root <- file.path(stage_dir, "GSEA_all_celltype_region")
celltypes <- c("Exc", "Inh", "Ast", "Oli", "OPC", "Mic_Immune")
regions <- c("AG", "EC", "HC", "MTC", "PFC", "TH")

tables <- lapply(regions, function(region) {
  do.call(rbind, lapply(celltypes, function(celltype) {
    context <- paste(celltype, region, sep = "_")
    path <- file.path(gsea_root, context, paste0(context, "_aging_GSEA.tsv"))
    dat <- read.delim(path, check.names = FALSE)
    if (nrow(dat) != 20L || anyDuplicated(dat$module) ||
        !setequal(dat$module, paste0("M", 1:20)) ||
        any(dat$context != context)) {
      stop("Expected 20 distinct modules for ", context)
    }
    expected <- p.adjust(dat$pval, method = "BH")
    if (any(!is.finite(dat$pval)) || any(!is.finite(dat$FDR)) ||
        any(abs(dat$FDR - expected) > 1e-8)) {
      stop("Saved BH q-values do not match 20-module correction for ", context)
    }
    dat$region <- region
    dat$recomputed_q20 <- expected
    dat
  }))
})
all_results <- do.call(rbind, tables)
if (nrow(all_results) != 720L) stop("Expected 36 contexts × 20 modules")

pfc <- all_results$region == "PFC"
selected <- readLines(file.path(shared_root, "part2_run_REGA/08pathology_tests/selected_module.txt"))
selected_rows <- paste(all_results$region, all_results$module, sep = "_") %in% selected
if (sum(pfc) != 120L || sum(selected_rows) != 36L) stop("Unexpected BH family size")

all_results$q120_pfc <- NA_real_
all_results$q120_pfc[pfc] <- p.adjust(all_results$pval[pfc], method = "BH")
all_results$q720_all_contexts <- p.adjust(all_results$pval, method = "BH")
all_results$q36_selected_only <- NA_real_
all_results$q36_selected_only[selected_rows] <-
  p.adjust(all_results$pval[selected_rows], method = "BH")

# Sensitivity analysis: reported zeros elsewhere are bounded at the resolution
# of 1,000 permutations. The PFC_M19 nominal p-values themselves are nonzero.
p_floor <- pmax(all_results$pval, 1 / 1001)
all_results$q720_with_zero_floor <- p.adjust(p_floor, method = "BH")

audit <- all_results[pfc & all_results$module == "M19", ]
audit <- audit[match(celltypes, audit$celltype), ]
audit <- data.frame(
  celltype = audit$celltype,
  region_module = "PFC_M19",
  NES = audit$NES,
  nominal_p = audit$pval,
  saved_q20 = audit$FDR,
  recomputed_q20 = audit$recomputed_q20,
  q36_selected_only_not_used = audit$q36_selected_only,
  q120_all_PFC = audit$q120_pfc,
  q720_all_contexts = audit$q720_all_contexts,
  q720_with_zero_floor = audit$q720_with_zero_floor,
  significant_q20 = audit$NES > 0 & audit$recomputed_q20 < 0.05,
  significant_q120 = audit$NES > 0 & audit$q120_pfc < 0.05,
  significant_q720 = audit$NES > 0 & audit$q720_all_contexts < 0.05
)
if (sum(audit$significant_q20) != 4L ||
    !identical(audit$significant_q20, audit$significant_q120) ||
    !identical(audit$significant_q20, audit$significant_q720)) {
  stop("PFC_M19 significance changes across BH families")
}

out <- file.path(stage_dir, "PFC_M19_aging_multiple_testing_audit.tsv")
write.table(audit, out, sep = "\t", quote = FALSE, row.names = FALSE)
print(audit, row.names = FALSE)
cat("Saved:", out, "\n")
cat("Verified: 36 contexts × 20 modules; max saved/recomputed q20 difference =",
    max(abs(all_results$FDR - all_results$recomputed_q20)), "\n")
