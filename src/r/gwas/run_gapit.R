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

# Ensure packages
pkg_install <- function(pkgs){
  for(p in pkgs){
    if(!requireNamespace(p, quietly=TRUE)){
      install.packages(p, repos="https://cloud.r-project.org")
    }
    suppressPackageStartupMessages(library(p, character.only=TRUE))
  }
}

pkg_install(c("GAPIT3","data.table","qqman","plink2R"))

# Read PLINK in R
read_plink <- function(bed,bim,fam){
  library(plink2R)
  ret <- plink2R::read_plink(bfile = sub("\\.bed$", "", bed))
  ret
}

# Load genotype
geno <- read_plink(opt$bed, opt$bim, opt$fam)
G <- geno$bed
SNP <- geno$bim
IND <- geno$fam

# Phenotype
pheno <- fread(opt$pheno)
setnames(pheno, names(pheno), make.names(names(pheno)))
idcol <- make.names(opt$idcol)
trait <- make.names(opt$trait)
stopifnot(idcol %in% names(pheno), trait %in% names(pheno))

# Merge phenotype order with fam
ph <- merge(data.table(IID=IND$fid, IID2=IND$iid)[,.(IID=paste(IID))], setnames(copy(pheno), idcol, "IID"), by="IID", all.x=TRUE)

# Covariates
covar <- NULL
if(!is.na(opt$covar) && file.exists(opt$covar)){
  cov <- fread(opt$covar)
  setnames(cov, names(cov), make.names(names(cov)))
  cov <- merge(data.table(IID=IND$fid, IID2=IND$iid)[,.(IID=paste(IID))], setnames(copy(cov), make.names(opt$idcol), "IID"), by="IID", all.x=TRUE)
  covar <- as.matrix(cov[, setdiff(names(cov), c("IID")), with=FALSE])
}

models <- strsplit(opt$models, ",")[[1]]
models <- toupper(trimws(models))

outdir <- opt$outdir
if(!dir.exists(outdir)) dir.create(outdir, recursive=TRUE)
setwd(outdir)

# Prepare GAPIT input
Y <- data.frame(IID=ph$IID, y=ph[[trait]])
M <- data.frame(SNP=SNP$snp, Chromosome=SNP$chr, Position=SNP$pos)
GD <- G
GD[is.na(GD)] <- mean(GD, na.rm=TRUE)

# Run GAPIT
suppressPackageStartupMessages(library(GAPIT3))

run_model <- function(model){
  message("Running ", model)
  fileprefix <- paste0("GAPIT_",model)
  if(model %in% c("GLM","MLM","FarmCPU","Blink")){
    res <- GAPIT(
      Y=Y,
      GD=GD,
      GM=M,
      CV=covar,
      model=model,
      file.output=TRUE,
      file.output.Manhattan=FALSE,
      file.output.QQ=FALSE,
      name.of.trait=opt$trait,
      output.fn=fileprefix
    )
  }
}

for(m in models){
  try(run_model(m), silent=TRUE)
}

# Collect results into RESULTS.tsv if exists
res_files <- list.files(pattern="GAPIT\..*\.GWAS\.Results\.csv$", recursive=TRUE, full.names=TRUE)
if(length(res_files)>0){
  all <- rbindlist(lapply(res_files, fread), fill=TRUE)
  fwrite(all, file="RESULTS.tsv", sep='\t')
}
