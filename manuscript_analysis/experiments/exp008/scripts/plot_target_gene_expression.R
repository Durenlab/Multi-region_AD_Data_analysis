#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L || length(args) > 5L) {
  stop("Usage: plot_target_gene_expression.R <EXP008-results-root> <pair> <TF> <target_gene> [unrestricted|module_restricted]")
}
root <- args[[1]]
pair <- args[[2]]
tf <- args[[3]]
gene <- args[[4]]
mode <- if (length(args) == 5L) args[[5]] else "unrestricted"
if (!mode %in% c("unrestricted", "module_restricted")) stop("Unknown target mode: ", mode)

mode_root <- file.path(root, mode)
expression_file <- file.path(mode_root, "expression_by_gene_individual", pair,
                             paste0(pair, "_target_gene_expression_by_individual.tsv.gz"))
test_file <- file.path(mode_root, "regression_results", pair,
                       paste0(pair, "_TG_level_regression.tsv.gz"))
out_dir <- file.path(mode_root, "plots", pair)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read_tsv(expression_file, show_col_types = FALSE,
                col_types = cols(.default = col_character())) |>
  filter(.data$pair == .env$pair, .data$TF == .env$tf,
         .data$mode == .env$mode, .data$target_gene == .env$gene) |>
  mutate(expression = suppressWarnings(as.numeric(.data$expression)))
if (!nrow(dat)) stop("No selected expression records for ", pair, "/", tf, "/", gene, "/", mode)
dat$Pathology <- factor(dat$Pathology, levels = c("nonAD", "earlyAD", "lateAD"))

test <- read_tsv(test_file, show_col_types = FALSE,
                 col_types = cols(.default = col_character())) |>
  filter(.data$pair == .env$pair, .data$TF == .env$tf,
         .data$mode == .env$mode, .data$target_gene == .env$gene) |>
  mutate(across(c(p_value, fdr_within_TF_pair_mode), as.numeric))
if (nrow(test) != 1L) stop("Expected one target-gene regression row; found ", nrow(test))
n <- table(factor(dat$Pathology[is.finite(dat$expression)],
                  levels = c("nonAD", "earlyAD", "lateAD")))
format_p <- function(x) {
  if (!is.finite(x)) return("NA")
  if (x < 0.001) return(sprintf("%.2e", x))
  sprintf("%.3f", x)
}
subtitle <- sprintf("EXP008 right-tailed trend p = %s\nEXP008 within-TF BH FDR = %s",
                    format_p(test$p_value[[1]]),
                    format_p(test$fdr_within_TF_pair_mode[[1]]))

pal <- c(nonAD = "#66AABB", earlyAD = "#EEAA55", lateAD = "#CC6677")
p <- ggplot(dat, aes(x = Pathology, y = expression, fill = Pathology, color = Pathology)) +
  geom_boxplot(outlier.shape = NA, width = 0.45, alpha = 0.65,
               color = "#333333", linewidth = 0.5) +
  geom_jitter(width = 0.045, height = 0, size = 1.6, alpha = 0.75,
              show.legend = FALSE) +
  scale_x_discrete(limits = c("nonAD", "earlyAD", "lateAD"), drop = FALSE) +
  scale_fill_manual(values = pal, drop = FALSE) +
  scale_color_manual(values = pal, drop = FALSE) +
  labs(title = paste(pair, "-", gene, "target of", tf), subtitle = subtitle,
       x = "Pathology", y = "Target-gene expression") +
  theme_classic(base_size = 9) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, margin = margin(b = 22)),
        plot.subtitle = element_text(hjust = 0.5, size = 7),
        axis.text.x = element_text(angle = 0, hjust = 0.5))

out <- file.path(out_dir, paste0(tf, "_", gene, "_", mode, "_expression_boxplot.pdf"))
ggsave(out, p, width = 2.8, height = 3.5, units = "in", limitsize = FALSE)
cat("Saved ", out, "\n", sep = "")
