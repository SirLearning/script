#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(qqman)
})

option_list <- list(
  make_option(c("--input"), type="character"),
  make_option(c("--outprefix"), type="character", default="plot")
)
opt <- parse_args(OptionParser(option_list=option_list))

d <- fread(opt$input)
# Try to guess columns
cn <- tolower(names(d))
get_col <- function(keys){
  for(k in keys){
    i <- which(cn==tolower(k))
    if(length(i)>0) return(names(d)[i[1]])
  }
  return(NA)
}
col_snp <- get_col(c("snp","marker","rsid","rs"))
col_chr <- get_col(c("chr","chrom","chromosome"))
col_pos <- get_col(c("pos","position","bp"))
col_p   <- get_col(c("p","pvalue","p.value","p_wald","p_lrt","p_score"))

if(any(is.na(c(col_snp,col_chr,col_pos,col_p)))){
  stop("Cannot determine required columns from input. Need SNP, CHR, POS, P")
}

man <- data.frame(SNP = d[[col_snp]], CHR=as.numeric(d[[col_chr]]), BP=as.numeric(d[[col_pos]]), P=as.numeric(d[[col_p]]))

png(paste0(opt$outprefix, ".manhattan.png"), width=1600, height=900)
try(manhattan(man, genomewideline=-log10(5e-8), suggestiveline=-log10(1e-5), cex=0.6), silent=TRUE)
dev.off()

png(paste0(opt$outprefix, ".qq.png"), width=1000, height=1000)
try(qq(man$P), silent=TRUE)
dev.off()
