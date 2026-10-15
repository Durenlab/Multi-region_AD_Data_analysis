suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else "."
output_dir <- if (length(args) >= 2) args[[2]] else file.path(project_root, "experiments", "EXP008_results", "SPP1_expression")
expression_dir <- "/project/zduren/durenlab/naqing_palmetto/AD_REGA/part1_prepare_REGA/04cleaned_gex"
pathology_file <- "/project/zduren/durenlab/naqing_palmetto/AD_Driver_TF/individual_pathology.tsv"

dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

pathology <- read.delim(pathology_file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c("ROSMAP_IndividualID", "Pathology") %in% names(pathology)))
pathology <- pathology[, c("ROSMAP_IndividualID", "Pathology")]
pathology <- pathology[!duplicated(pathology$ROSMAP_IndividualID), , drop = FALSE]
pathology$Pathology <- factor(pathology$Pathology,
  levels = c("nonAD", "earlyAD", "lateAD"), ordered = TRUE)

celltypes <- c("Oli", "Mic_Immune")
long <- list()
for (celltype in celltypes) {
  pair <- paste(celltype, "EC", sep = "_")
  path <- file.path(expression_dir, paste0("cleaned_", pair, ".tsv"))
  if (!file.exists(path)) stop("Missing expression matrix: ", path)
  expr <- read.delim(path, check.names = FALSE, row.names = 1)
  if (!("SPP1" %in% rownames(expr))) stop("SPP1 is absent from ", path)
  values <- data.frame(
    ROSMAP_IndividualID = colnames(expr),
    expression = as.numeric(expr["SPP1", ]),
    stringsAsFactors = FALSE
  )
  if (anyDuplicated(values$ROSMAP_IndividualID)) stop("Duplicate donor IDs in ", path)
  merged <- merge(values, pathology, by = "ROSMAP_IndividualID", all.x = TRUE, sort = FALSE)
  if (anyNA(merged$Pathology)) stop("Unmatched donor pathology labels for ", pair)
  merged$celltype_region <- pair
  long[[pair]] <- merged
}

plot_data <- do.call(rbind, long)
plot_data$celltype_region <- factor(plot_data$celltype_region,
  levels = c("Oli_EC", "Mic_Immune_EC"))
plot_data$pathology_code <- as.integer(plot_data$Pathology) - 1L
group_n <- aggregate(ROSMAP_IndividualID ~ celltype_region + Pathology,
  data = plot_data, FUN = length)
names(group_n)[[3]] <- "n_donors"
expected <- c(nonAD = 25L, earlyAD = 14L, lateAD = 7L)
for (pair in levels(plot_data$celltype_region)) {
  observed <- setNames(group_n$n_donors[group_n$celltype_region == pair],
    as.character(group_n$Pathology[group_n$celltype_region == pair]))
  if (!identical(as.integer(observed[names(expected)]), as.integer(expected))) {
    stop("Unexpected pathology group counts in ", pair)
  }
}

trend_results <- do.call(rbind, lapply(levels(plot_data$celltype_region), function(pair) {
  pair_data <- plot_data[plot_data$celltype_region == pair, , drop = FALSE]
  fit <- lm(expression ~ pathology_code, data = pair_data)
  coef_table <- summary(fit)$coefficients
  beta <- unname(coef_table["pathology_code", "Estimate"])
  t_stat <- unname(coef_table["pathology_code", "t value"])
  p_value_two_sided <- unname(coef_table["pathology_code", "Pr(>|t|)"])
  p_value_right_tailed <- pt(t_stat, df = df.residual(fit), lower.tail = FALSE)
  ci <- confint(fit, "pathology_code", level = 0.95)
  data.frame(
    celltype_region = pair,
    n_donors = nrow(pair_data),
    pathology_coding = "nonAD=0, earlyAD=1, lateAD=2",
    model = "SPP1 expression ~ pathology_code",
    slope_per_pathology_step = beta,
    slope_ci95_low = ci[1],
    slope_ci95_high = ci[2],
    p_value_right_tailed_beta_gt_0 = p_value_right_tailed,
    p_value_two_sided = p_value_two_sided,
    trend_label = paste0("Linear trend: β=", formatC(beta, format = "f", digits = 2),
      ", p(β>0)=", format.pval(p_value_right_tailed, digits = 2, eps = 0.001)),
    stringsAsFactors = FALSE
  )
}))
trend_results$celltype_region <- factor(trend_results$celltype_region,
  levels = c("Oli_EC", "Mic_Immune_EC"))
facet_labels <- setNames(
  paste0(as.character(trend_results$celltype_region), "\n",
    sub("^Linear trend: ", "", trend_results$trend_label)),
  as.character(trend_results$celltype_region)
)

write.table(plot_data, file.path(output_dir, "tables", "SPP1_expression_by_donor.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(group_n, file.path(output_dir, "tables", "SPP1_expression_group_counts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(trend_results, file.path(output_dir, "tables", "SPP1_expression_pathology_trends.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

pathology_colors <- c(nonAD = "#66AABB", earlyAD = "#EEAA55", lateAD = "#CC6677")
x_labels <- c(nonAD = "nonAD (n=25)", earlyAD = "earlyAD (n=14)", lateAD = "lateAD (n=7)")
p <- ggplot(plot_data, aes(x = Pathology, y = expression, fill = Pathology, color = Pathology)) +
  geom_boxplot(width = 0.64, alpha = 0.42, outlier.shape = NA, color = "#333333", linewidth = 0.55) +
  geom_jitter(width = 0.12, height = 0, size = 1.8, alpha = 0.72) +
  facet_wrap(~celltype_region, nrow = 1, scales = "free_y",
    labeller = labeller(celltype_region = as_labeller(facet_labels))) +
  scale_fill_manual(values = pathology_colors, drop = FALSE) +
  scale_color_manual(values = pathology_colors, drop = FALSE) +
  scale_x_discrete(labels = x_labels, drop = FALSE) +
  labs(
    title = "SPP1 expression across pathology groups",
    subtitle = "Right-tailed trend test: H1 slope > 0; pathology coded 0–2",
    x = NULL,
    y = "SPP1 expression (log1p CPM, batch-adjusted)",
    caption = "Expression values are log1p CPM with sex, age at death, and PMI adjusted; pathology was protected during batch correction."
  ) +
  guides(fill = "none", color = "none") +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(size = 8, hjust = 0),
    axis.text.x = element_text(color = "#333333")
  )

ggsave(file.path(output_dir, "figures", "SPP1_Oli_EC_Mic_Immune_EC_expression_by_pathology.pdf"),
  p, width = 8.2, height = 4.8, device = cairo_pdf)
ggsave(file.path(output_dir, "figures", "SPP1_Oli_EC_Mic_Immune_EC_expression_by_pathology.png"),
  p, width = 8.2, height = 4.8, dpi = 300)
cat("Wrote SPP1 donor-level plot and tables to", output_dir, "\n")
print(group_n, row.names = FALSE)
