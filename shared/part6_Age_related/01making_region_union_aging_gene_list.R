library(readxl)
library(dplyr)

stage_dir <- file.path(Sys.getenv("REGA_SHARED_ROOT", normalizePath("..")), "part6_Age_related")
file <- file.path(stage_dir, "41514_2026_391_MOESM7_ESM.xlsx")

# external cell type -> our cell type
celltype_map <- list(
  Exc = "ExN",
  Inh = "InN",
  Ast = "Astrocyte",
  Oli = "Oligodendrocyte",
  OPC = "OPC",
  Mic_Immune = "Microglia",
  Vasc = c("Endothelial", "Mural")
)

aging_sets <- list()
summary_all <- list()

for(ct in names(celltype_map)) {

  # read corresponding sheet(s)
  dat <- bind_rows(
    lapply(celltype_map[[ct]], function(x)
      read_excel(file, sheet = x))
  ) %>%
    select(gene_name, region, direction) %>%
    distinct()

  # genes in each region
  region_count <- dat %>%
    count(region, name = "n_genes")

  # genes with opposite direction
  opposite_genes <- dat %>%
    group_by(gene_name) %>%
    summarise(n_direction = n_distinct(direction)) %>%
    filter(n_direction > 1) %>%
    pull(gene_name)

  # final unified aging gene set
  aging_genes <- dat %>%
    filter(!gene_name %in% opposite_genes) %>%
    distinct(gene_name, direction)

  aging_sets[[ct]] <- aging_genes

  cat("\n###", ct, "###\n")
  print(region_count)

  cat("Unique genes before merging:",
      n_distinct(dat$gene_name), "\n")

  cat("Removed opposite genes:",
      length(opposite_genes), "\n")

  cat("Final genes:",
      nrow(aging_genes), "\n")

  cat("Up:", sum(aging_genes$direction == 1),
      " Down:", sum(aging_genes$direction == -1), "\n")
}

saveRDS(aging_sets, file.path(stage_dir, "All_celltype_all_direction_aging_sets.rds"))
