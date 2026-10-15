#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(Matrix))
root <- normalizePath('.')
tab <- file.path(root,'experiments/EXP004_results/tables')
x <- read.delim(file.path(tab,'grn_queries.tsv'),check.names=FALSE,stringsAsFactors=FALSE)
x$edge_value <- NA_real_; x$edge_status <- 'not_read'; x$matrix_source <- NA_character_; x$metadata_source <- NA_character_
x$TG_row_index <- NA_integer_; x$TF_column_index <- NA_integer_; x$within_TF_target_percentile <- NA_real_
meta_rows <- list()
for(pair in unique(x$celltype_region)) {
 p <- file.path(root,'experiments/EXP005_results',paste0(pair,'_TG_TF.rds'))
 mp <- file.path(root,'experiments/EXP005_results',paste0(pair,'_GRN_metadata.rds'))
 ii <- which(x$celltype_region==pair); x$matrix_source[ii] <- p; x$metadata_source[ii] <- mp
 if(!file.exists(p)) {x$edge_status[ii] <- 'matrix_unavailable'; next}
 m <- readRDS(p); md <- readRDS(mp)
 ri <- match(x$TG[ii],rownames(m));ci <- match(x$TF[ii],colnames(m));ok <- !is.na(ri)&!is.na(ci)
 x$TG_row_index[ii] <- ri; x$TF_column_index[ii] <- ci
 x$edge_status[ii] <- 'identifier_absent';x$edge_value[ii[ok]] <- as.numeric(m[cbind(ri[ok],ci[ok])])
 x$edge_status[ii[ok]] <- ifelse(x$edge_value[ii[ok]]>0,'positive_weight','zero_or_negative')
 for(tf in unique(x$TF[ii[ok]])) {
  jj <- ii[ok][x$TF[ii[ok]]==tf]; v <- as.numeric(m[,tf]);percentile <- rank(v,ties.method='average')/length(v)
  x$within_TF_target_percentile[jj] <- percentile[match(x$TG[jj],rownames(m))]
 }
 meta_rows[[pair]] <- data.frame(celltype_region=pair,n_genes=nrow(m),n_TFs=ncol(m),n_nonzero=nnzero(m),density=nnzero(m)/(nrow(m)*ncol(m)),modules=paste(md$modules$TG_TF,collapse=';'),matrix_source=p,metadata_source=mp,WBV_input=unname(md$inputs['WBV']),V_input=unname(md$inputs['V']),PCC_input=unname(md$inputs['PCC']))
 cat(pair,': queried',length(ii),'edges; density',nnzero(m)/(nrow(m)*ncol(m)),'\n')
 rm(m,md);gc(verbose=FALSE)
}
write.table(x,file.path(tab,'grn_edges.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
write.table(do.call(rbind,meta_rows),file.path(tab,'grn_metadata.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
# Reproduce the notebook's exact TF-name handling and denominator.
sig <- read.table(file.path(dirname(root),'AD_REGA/part7_TF/sig_Genes_pval.tsv'))
prior_path <- file.path(dirname(root),'AD_REGA/part7_TF/motif_by_TF_matrix.tsv')
prior_names <- colnames(read.delim(prior_path,nrows=1,header=TRUE))[-1]
notebook_genes <- intersect(prior_names,sig$V1)
write.table(data.frame(TF=notebook_genes,notebook_denominator=length(notebook_genes)),file.path(tab,'magma_exact_R_family.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
