library(dplyr)
stage_dir <- file.path(Sys.getenv("REGA_SHARED_ROOT", normalizePath("..")), "part6_Age_related")
aging_sets <- readRDS(file.path(stage_dir, "All_celltype_all_direction_aging_sets.rds"))
# keep only aging-up genes
aging_up_sets <- lapply(aging_sets, function(x) {
  x %>%
    filter(direction == 1) %>%
    pull(gene_name) %>%
    unique()
})

# check gene number
sapply(aging_up_sets, length)

# write GMT
gmt_lines <- sapply(names(aging_up_sets), function(ct) {
  paste(
    c(
      paste0(ct, "_AGING_UP"),
      paste0(ct, "_aging_upregulated_genes"),
      aging_up_sets[[ct]]
    ),
    collapse = "\t"
  )
})

writeLines(
  gmt_lines,
  file.path(stage_dir, "aging_up_gene_sets.gmt")
)
