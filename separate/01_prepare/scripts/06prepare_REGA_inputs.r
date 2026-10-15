## =========================

## 构建 RE–TG 与 Peak–Motif 矩阵（Region × Celltype 批处理版）

## =========================

suppressPackageStartupMessages({
library(dplyr)
library(Matrix)
library(rhdf5)
})

options(stringsAsFactors = FALSE)

## ============ 参数 ============

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript 06prepare_REGA_inputs.r PREP_DIR")
prep <- args[[1]]
celltypes <- c("Ast", "Exc", "Inh", "Mic_Immune", "Oli", "OPC")
regions   <- c("AG", "EC", "HC", "MTC", "PFC", "TH")

tau <- 1e5  # 距离衰减尺度（100kb）

## ============ 主循环 ============

for (ct in celltypes) {
## ---- celltype-level inputs ----

rt_path    <- file.path(prep, "re_tg_matrices", paste0(ct, "_RE_TG_matrix.bed"))
motif_path <- file.path(prep, "motif_scan_results", paste0(ct, "_peak_motif.txt"))
## ---- read invariant inputs ----

re_tg_overlap <- read.table(rt_path, header = FALSE, sep = "\t")
motifs <- read.table(motif_path, header = FALSE, sep = "\t")

## ---- 标准化列名 ----

colnames(re_tg_overlap) <- c(
"g_chr","g_start","g_end","gene",
"p_chr","p_start","p_end","peak"
)

## ---- 峰集合交集（只做一次） ----

peaks_rt  <- unique(re_tg_overlap$peak)
peaks_mtf <- unique(motifs$V1)
filtered_elements_global <- intersect(peaks_rt, peaks_mtf)

motifs <- motifs[motifs$V1 %in% filtered_elements_global, , drop = FALSE]
re_tg_overlap <- re_tg_overlap %>%
filter(peak %in% filtered_elements_global)

## ============ region 循环 ============

for (rg in regions) {
## ---- 路径 ----
gex_path <- file.path(prep, "cleaned_expression", paste0("cleaned_", ct, "_", rg, ".tsv"))

if (!file.exists(gex_path)) { 
    stop("Missing cleaned expression: ", gex_path)
}
    
out_dir <- file.path(prep, "rega_inputs", rg, ct)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

re_tg_h5       <- file.path(out_dir, "re_tg_matrix.h5")
gex_outfile    <- file.path(out_dir, "input_GEX.csv")
motif_outfile  <- file.path(out_dir, "input_binding.csv")

## ---- 读入 GEX ----
gex_csv <- read.table(gex_path)
gex_mat <- as.matrix(gex_csv)
gex_mat[gex_mat < 0] <- 0
adjustedExp <- gex_mat

gene_list <- rownames(adjustedExp)

## ---- 基因过滤 ----
rt <- re_tg_overlap
genes_valid <- intersect(gene_list, unique(rt$gene))
rt <- rt %>% filter(gene %in% genes_valid)

## ---- RE–TG 权重 ----
#TSS      <- rt$g_end
window_size <- 300000
TSS      <- rt$g_end - window_size
peak_mid <- (rt$p_start + rt$p_end) / 2
dist     <- abs(peak_mid - TSS)
weight   <- exp(-dist / tau)

matrix_rows1 <- genes_valid
matrix_cols1 <- unique(rt$peak)

i <- match(rt$gene, matrix_rows1)
j <- match(rt$peak, matrix_cols1)
keep <- !is.na(i) & !is.na(j)

i <- i[keep]; j <- j[keep]; w <- weight[keep]

## ---- 聚合重复 ----
if (any(duplicated(data.frame(i, j)))) {
  dfw <- data.frame(i = i, j = j, w = w) %>%
    group_by(i, j) %>%
    summarise(w = max(w), .groups = "drop")

  i <- dfw$i; j <- dfw$j; w <- dfw$w
}

re_tg_matrix <- sparseMatrix(
  i = i, j = j, x = as.numeric(w),
  dims = c(length(matrix_rows1), length(matrix_cols1)),
  dimnames = list(matrix_rows1, matrix_cols1)
)

## ---- Peak–Motif ----
matrix_rows_peak  <- matrix_cols1
motifs_list       <- sort(unique(motifs$V2))

motifs2 <- motifs %>%
  select(V1, V2) %>%
  rename(peak = V1, motif = V2) %>%
  filter(peak %in% matrix_rows_peak)

ii <- match(motifs2$peak, matrix_rows_peak)
jj <- match(motifs2$motif, motifs_list)

keep2 <- !is.na(ii) & !is.na(jj)
ii <- ii[keep2]; jj <- jj[keep2]

score_matrix_dense <- matrix(
  0L,
  nrow = length(matrix_rows_peak),
  ncol = length(motifs_list),
  dimnames = list(matrix_rows_peak, motifs_list)
)

score_matrix_dense[cbind(ii, jj)] <- 1L

score_matrix <- sparseMatrix(
  i = ii, j = jj, x = rep.int(1L, length(ii)),
  dims = c(length(matrix_rows_peak), length(motifs_list)),
  dimnames = list(matrix_rows_peak, motifs_list)
)

## ---- 对齐检查 ----
if (!identical(colnames(re_tg_matrix), rownames(score_matrix))) {
  common_peaks <- intersect(colnames(re_tg_matrix), rownames(score_matrix))
  re_tg_matrix <- re_tg_matrix[, common_peaks, drop = FALSE]
  score_matrix <- score_matrix[common_peaks, , drop = FALSE]
}
## ---- 诊断 ----
cat(sprintf("[RE–TG] %dx%d | [Motif] %dx%d\n",
            nrow(re_tg_matrix), ncol(re_tg_matrix),
            nrow(score_matrix), ncol(score_matrix)))

score_matrix_dense <- as.matrix(score_matrix)
## ---- 输出 RE–TG (HDF5) ----
if (file.exists(re_tg_h5)) file.remove(re_tg_h5)

trip <- summary(re_tg_matrix)

h5write(trip$i,             re_tg_h5, "i")
h5write(trip$j,             re_tg_h5, "j")
h5write(as.numeric(trip$x), re_tg_h5, "x")
h5write(dim(re_tg_matrix),  re_tg_h5, "dims")
h5write(rownames(re_tg_matrix), re_tg_h5, "rownames")
h5write(colnames(re_tg_matrix), re_tg_h5, "colnames")

## ---- 输出 GEX ----
gex_filtered <- adjustedExp[rownames(adjustedExp) %in% genes_valid, , drop = FALSE]
gex_filtered <- gex_filtered[match(genes_valid, rownames(gex_filtered)), , drop = FALSE]

write.csv(gex_filtered, file = gex_outfile)

## ---- 输出 motif ----
write.csv(score_matrix_dense, file = motif_outfile)

cat("Saved outputs to:", out_dir, "\n")
}
}
