#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option(c("--bed"), type="character"),
  make_option(c("--bim"), type="character"),
  make_option(c("--fam"), type="character"),
  make_option(c("--pheno"), type="character"),
  make_option(c("--idcol"), type="character", default="IID"),
  make_option(c("--trait"), type="character"),
  make_option(c("--covar"), type="character", default=NA),
  make_option(c("--models"), type="character", default="GLM,MLM,FarmCPU,Blink"),
  make_option(c("--outdir"), type="character", default="gwas_results")
)
opt <- parse_args(OptionParser(option_list=option_list))

pkg_install <- function(pkgs){
  for(p in pkgs){
    if(!requireNamespace(p, quietly=TRUE)){
      install.packages(p, repos="https://cloud.r-project.org")
    }
    suppressPackageStartupMessages(library(p, character.only=TRUE))
  }
}

pkg_install(c("rMVP","data.table","qqman","bigmemory","bigsnpr","plink2R"))

# Genotype input via PLINK
library(plink2R)
ret <- plink2R::read_plink(bfile = sub("\\.bed$", "", opt$bed))
G <- ret$bed
SNP <- ret$bim
IND <- ret$fam

pheno <- fread(opt$pheno)
setnames(pheno, names(pheno), make.names(names(pheno)))
idcol <- make.names(opt$idcol)
trait <- make.names(opt$trait)
stopifnot(idcol %in% names(pheno), trait %in% names(pheno))

ph <- merge(data.table(IID=IND$fid, IID2=IND$iid)[,.(IID=paste(IID))], setnames(copy(pheno), idcol, "IID"), by="IID", all.x=TRUE)

covar <- NULL
if(!is.na(opt$covar) && file.exists(opt$covar)){
  cov <- fread(opt$covar)
  setnames(cov, names(cov), make.names(names(cov)))
  cov <- merge(data.table(IID=IND$fid, IID2=IND$iid)[,.(IID=paste(IID))], setnames(copy(cov), make.names(opt$idcol), "IID"), by="IID", all.x=TRUE)
  covar <- as.matrix(cov[, setdiff(names(cov), c("IID")), with=FALSE])
}

# Prepare rMVP inputs
library(rMVP)

# MVP requires numeric genotype: convert NA to mean
Gm <- G
Gm[is.na(Gm)] <- mean(Gm, na.rm=TRUE)

M <- data.frame(SNP=SNP$snp, Chromosome=SNP$chr, Position=SNP$pos)

models <- toupper(trimws(strsplit(opt$models, ",")[[1]]))
if(!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive=TRUE)
setwd(opt$outdir)

for(m in models){
  message("Running ", m)
  try({
    MVP(
      phe=as.matrix(ph[[trait]]),
      geno=as.matrix(Gm),
      map=M,
      covariates=covar,
      nPC.MAF=5,
      method=m,
      threshold=0.05/nrow(M),
      file.output=TRUE,
      MEMORY=TRUE,
      drawQQ=FALSE,
      drawManhattan=FALSE,
      out=paste0("rMVP_", m)
    )
  }, silent=TRUE)
}

# Collect results
res <- list.files(pattern="MVP\.Geno\.association\..*\.csv$", full.names=TRUE, recursive=TRUE)
if(length(res)>0){
  all <- rbindlist(lapply(res, fread), fill=TRUE)
  fwrite(all, file="RESULTS.tsv", sep='\t')
}
