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
library(data.table)

# ============================================================
# 1. Paths
# ============================================================

# 修改为包含 PFC、AG、EC、MTC、HC、TH 文件夹的目录
BASE_DIR <- file.path(shared_root, "part2_run_REGA/08pathology_tests")

OUT_FILE <- file.path(
  BASE_DIR,
  "ALL_regions_sharedH_All_four_tests_H_pathology_results.tsv"
)

# ============================================================
# 2. Region order
# ============================================================

region_order <- c("PFC", "AG", "EC", "MTC", "HC", "TH")

# ============================================================
# 3. Read one region
# ============================================================

read_region_result <- function(region) {

  input_file <- file.path(
    BASE_DIR,
    region,
    paste0(
      region,
      "_sharedH_All_four_tests_H_pathology_results.tsv"
    )
  )

  if (!file.exists(input_file)) {
    stop("Input file does not exist: ", input_file)
  }

  message("Reading: ", input_file)

  result <- fread(
    input_file,
    sep = "\t",
    header = TRUE,
    na.strings = c("", "NA", "NaN"),
    fill = TRUE,
    data.table = TRUE
  )

  if (nrow(result) == 0) {
    warning("File contains zero rows: ", input_file)
  }

  # Add Region if it is absent
  if (!"Region" %in% colnames(result)) {
    result[, Region := region]
  }

  # Check Region values
  observed_regions <- unique(na.omit(result$Region))

  if (
    length(observed_regions) > 0 &&
    any(observed_regions != region)
  ) {
    warning(
      "Unexpected Region value in ", input_file,
      ". Observed: ",
      paste(observed_regions, collapse = ", "),
      ". The Region column will be set to ", region, "."
    )
  }

  # Ensure the correct region label
  result[, Region := region]

  # Record original file
  result[, source_file := basename(input_file)]

  return(result)
}

# ============================================================
# 4. Read and combine all regions
# ============================================================

result_list <- lapply(region_order, read_region_result)

combined_results <- rbindlist(
  result_list,
  use.names = TRUE,
  fill = TRUE
)

# ============================================================
# 5. Format and sort
# ============================================================

combined_results[
  ,
  Region := factor(
    Region,
    levels = region_order,
    ordered = TRUE
  )
]

# Extract module number for natural sorting: M1, M2, ..., M10
combined_results[
  ,
  module_number := as.integer(
    sub("^M", "", Module)
  )
]

# Optional test order
test_order <- c(
  "Kruskal_Wallis",
  "Linear_Regression",
  "Pairwise_MannWhitneyU",
  "Spearman"
)

combined_results[
  ,
  test := factor(
    test,
    levels = test_order,
    ordered = TRUE
  )
]

setorder(
  combined_results,
  Region,
  module_number,
  test,
  Comparison
)

# Return factor columns to normal character columns before saving
combined_results[, Region := as.character(Region)]
combined_results[, test := as.character(test)]

# Remove temporary sorting column
combined_results[, module_number := NULL]

# Put important columns first
first_columns <- c(
  "Region",
  "test",
  "Module",
  "Comparison"
)

remaining_columns <- setdiff(
  colnames(combined_results),
  first_columns
)

setcolorder(
  combined_results,
  c(first_columns, remaining_columns)
)

# ============================================================
# 6. Save combined results
# ============================================================

fwrite(
  combined_results,
  file = OUT_FILE,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ============================================================
# 7. Summary checks
# ============================================================

cat("\nCombined file written to:\n")
cat(OUT_FILE, "\n\n")

cat("Combined dimensions:\n")
cat(
  nrow(combined_results),
  "rows x",
  ncol(combined_results),
  "columns\n\n"
)

cat("Rows by Region:\n")
print(combined_results[, .N, by = Region])

cat("\nRows by Region and test:\n")
print(
  combined_results[
    ,
    .N,
    by = .(Region, test)
  ][
    order(
      match(Region, region_order),
      match(test, test_order)
    )
  ]
)

cat("\nUnique modules by Region:\n")
print(
  combined_results[
    ,
    .(
      n_modules = uniqueN(Module),
      modules = paste(
        sort(
          unique(Module),
          method = "radix"
        ),
        collapse = ", "
      )
    ),
    by = Region
  ]
)