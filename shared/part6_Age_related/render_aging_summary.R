suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
})

shared_root <- normalizePath(Sys.getenv("REGA_SHARED_ROOT", unset = normalizePath("../..")))
stage_dir <- file.path(shared_root, "part6_Age_related")
gsea_root <- file.path(stage_dir, "GSEA_all_celltype_region")
selected_file <- file.path(shared_root, "part2_run_REGA/08pathology_tests/selected_module.txt")
candidate_modules <- trimws(readLines(selected_file))
candidate_modules <- candidate_modules[nzchar(candidate_modules)]
celltypes <- c("Exc", "Inh", "Ast", "Oli", "OPC", "Mic_Immune")
regions <- c("AG", "EC", "HC", "MTC", "PFC", "TH")

if (anyDuplicated(candidate_modules) || !all(grepl("^[A-Z]+_M[0-9]+$", candidate_modules))) {
  stop("Invalid candidate-module selection: ", selected_file)
}

read_context <- function(region, celltype) {
  context <- paste(celltype, region, sep = "_")
  input <- file.path(gsea_root, context, paste0(context, "_aging_GSEA.tsv"))
  if (!file.exists(input)) stop("Missing aging GSEA table: ", input)
  dat <- read_tsv(input, show_col_types = FALSE)
  required <- c("context", "celltype", "geneset", "module", "ES", "NES", "pval", "FDR", "significant_positive")
  if (!all(required %in% names(dat))) stop("Missing aging GSEA columns: ", input)
  if (nrow(dat) != 20L || anyDuplicated(dat$module) ||
      !identical(sort(dat$module), sort(paste0("M", 1:20)))) {
    stop("Expected exactly 20 distinct modules in ", input)
  }
  if (any(dat$context != context) || any(dat$celltype != celltype) ||
      any(dat$geneset != paste0(celltype, "_AGING_UP"))) {
    stop("Context or gene set mismatch in ", input)
  }
  if (any(!is.finite(dat$pval)) || any(!is.finite(dat$FDR)) || any(!is.finite(dat$NES))) {
    stop("Non-finite statistics in ", input)
  }
  if (any(abs(dat$FDR - p.adjust(dat$pval, method = "BH")) > 1e-8)) {
    stop("FDR does not match BH across 20 modules in ", input)
  }
  reported <- tolower(as.character(dat$significant_positive)) %in% c("true", "t", "1")
  expected <- dat$NES > 0 & dat$FDR < 0.05
  if (!identical(reported, expected)) stop("Significance flag mismatch in ", input)
  dat %>% mutate(region = region, source_file = input, significant_positive = expected)
}

all_results <- expand_grid(region = regions, celltype = celltypes) %>%
  rowwise() %>%
  mutate(result = list(read_context(region, celltype))) %>%
  ungroup() %>%
  select(result) %>%
  unnest(result)
if (nrow(all_results) != 720L) stop("Expected 720 aging GSEA rows")

candidate_df <- tibble(candidate_module = candidate_modules) %>%
  separate(candidate_module, into = c("region", "module"), sep = "_", remove = FALSE)
candidate_results <- all_results %>%
  inner_join(candidate_df, by = c("region", "module")) %>%
  select(candidate_module, region, module, celltype, context, geneset,
         ES, NES, pval, FDR, leading_genes, significant_positive, source_file)
if (nrow(candidate_results) != length(candidate_modules) * length(celltypes) ||
    anyDuplicated(candidate_results[c("candidate_module", "celltype")])) {
  stop("Selected candidate-module grid is incomplete or duplicated")
}

significant_results <- candidate_results %>%
  filter(significant_positive) %>%
  arrange(match(candidate_module, candidate_modules), match(celltype, celltypes))
summary_table <- candidate_df %>%
  left_join(candidate_results %>%
              group_by(candidate_module) %>%
              summarise(n_significant_celltypes = sum(significant_positive),
                        significant_celltypes = paste(celltype[significant_positive], collapse = ";"),
                        .groups = "drop"),
            by = "candidate_module") %>%
  select(candidate_module, region, module, n_significant_celltypes, significant_celltypes)

write_tsv(candidate_results, file.path(gsea_root, "cross_celltype_candidate_modules_aging_GSEA_all.tsv"))
write_tsv(significant_results, file.path(gsea_root, "cross_celltype_candidate_modules_aging_GSEA_significant.tsv"))
write_tsv(summary_table, file.path(gsea_root, "cross_celltype_candidate_modules_aging_GSEA_summary.tsv"))
message("Selected modules: ", paste(candidate_modules, collapse = ", "))
message("Significant positive cell-type/module contexts: ", nrow(significant_results))
print(summary_table)

region_colors <- c(AG = "#8DA0CB", EC = "#FC8D62", HC = "#66C2A5",
                   MTC = "#E78AC3", PFC = "#B3B3B3", TH = "#A6D854")
plot_df <- candidate_results %>%
  mutate(x = match(candidate_module, candidate_modules),
         y = match(celltype, rev(celltypes)),
         neglog10_FDR = -log10(pmax(FDR, 1e-6)))
if (!setequal(unique(plot_df$candidate_module), candidate_modules) ||
    nrow(plot_df) != length(candidate_modules) * length(celltypes)) {
  stop("Plot data do not match the current candidate-module selection")
}
write_tsv(tibble(candidate_module = candidate_modules),
          file.path(gsea_root, "plotted_candidate_modules.tsv"))
strip_df <- candidate_df %>%
  mutate(x = match(candidate_module, candidate_modules), y = length(celltypes) + 0.75)

p <- ggplot() +
  geom_hline(yintercept = seq_along(celltypes), color = "grey92", linewidth = 0.4) +
  geom_vline(xintercept = seq_along(candidate_modules), color = "grey94", linewidth = 0.4) +
  geom_tile(data = strip_df, aes(x = x, y = y, fill = region),
            width = 0.82, height = 0.24) +
  geom_point(data = filter(plot_df, !significant_positive),
             aes(x = x, y = y), color = "grey82", size = 2.5) +
  geom_point(data = filter(plot_df, significant_positive),
             aes(x = x, y = y, size = NES, color = neglog10_FDR), alpha = 0.95) +
  scale_size_area(name = "NES", max_size = 10,
                  limits = c(0, max(plot_df$NES[plot_df$significant_positive], 1)),
                  breaks = c(1, 1.5, 2)) +
  scale_color_gradient(name = expression(-log[10]("BH q")),
                       low = "#D9E6F2", high = "#08519C") +
  scale_fill_manual(values = region_colors, name = "Region") +
  scale_x_continuous(breaks = seq_along(candidate_modules), labels = candidate_modules,
                     limits = c(0.4, length(candidate_modules) + 0.6), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_along(celltypes), labels = rev(celltypes),
                     limits = c(0.5, length(celltypes) + 1.05), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  guides(color = guide_colorbar(order = 1, title.position = "top",
                                barheight = grid::unit(2.0, "cm"),
                                barwidth = grid::unit(0.45, "cm")),
         fill = guide_legend(order = 2), size = guide_legend(order = 3)) +
  theme_classic(base_size = 14) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"), legend.position = "right",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))

ggsave(file.path(gsea_root, "cross_celltype_candidate_modules_aging_GSEA_bubbleplot.pdf"),
       p, width = 8, height = 5.5)
ggsave(file.path(gsea_root, "cross_celltype_candidate_modules_aging_GSEA_bubbleplot.png"),
       p, width = 8, height = 5.5, dpi = 300)

# ---- Labeled NES plot (historical notebook's final figure) ----
# Show only positive aging-enrichment hits that pass within-context BH q < 0.05.
# The value printed inside each dot is the normalized enrichment score (NES).
nes_df <- filter(plot_df, significant_positive)
p_nes <- ggplot(nes_df, aes(x = x, y = y)) +
  geom_point(aes(size = NES), color = "#2878A5", alpha = 0.9) +
  geom_text(aes(label = sprintf("%.1f", NES)), color = "white", size = 4.5) +
  scale_size_area(name = "NES", max_size = 12,
                  limits = c(0, max(nes_df$NES, 1)), breaks = c(1, 1.5, 2)) +
  scale_x_continuous(breaks = seq_along(candidate_modules), labels = candidate_modules,
                     limits = c(0.4, length(candidate_modules) + 0.6), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_along(celltypes), labels = rev(celltypes),
                     limits = c(0.5, length(celltypes) + 0.5), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 14) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
        panel.grid.minor = element_blank(), legend.position = "right")

ggsave(file.path(gsea_root, "cross_celltype_candidate_modules_aging_GSEA_bubbleplotNES_Dot.pdf"),
       p_nes, width = 8, height = 5.5)
ggsave(file.path(gsea_root, "cross_celltype_candidate_modules_aging_GSEA_bubbleplotNES_Dot.png"),
       p_nes, width = 8, height = 5.5, dpi = 300)
message("Saved aging summary and both figures in: ", gsea_root)
