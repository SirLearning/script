#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option(c("--pheno"), type="character"),
  make_option(c("--idcol"), type="character", default="IID"),
  make_option(c("--trait"), type="character"),
  make_option(c("--envcol"), type="character", default=""),
  make_option(c("--repcol"), type="character", default=""),
  make_option(c("--fixed"), type="character", default=""),
  make_option(c("--random"), type="character", default=""),
  make_option(c("--blue"), action="store_true", default=FALSE),
  make_option(c("--blup"), action="store_true", default=FALSE),
  make_option(c("--outdir"), type="character", default=".")
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

pkg_install(c("lme4","emmeans","ggplot2","data.table"))

ph <- fread(opt$pheno)
setnames(ph, names(ph), make.names(names(ph)))
id <- make.names(opt$idcol)
trait <- make.names(opt$trait)
if(!(id %in% names(ph) && trait %in% names(ph))) stop("Missing id or trait column in pheno")

env <- if(nchar(opt$envcol)>0) make.names(opt$envcol) else NA
rep <- if(nchar(opt$repcol)>0) make.names(opt$repcol) else NA
if(!is.na(env) && !(env %in% names(ph))) env <- NA
if(!is.na(rep) && !(rep %in% names(ph))) rep <- NA

# Convert to factors
for(col in c(id, env, rep)) if(!is.na(col)) ph[[col]] <- as.factor(ph[[col]])

# Summary stats
summ <- ph[, .(
  N=.N,
  mean=mean(get(trait), na.rm=TRUE),
  sd=sd(get(trait), na.rm=TRUE),
  min=min(get(trait), na.rm=TRUE),
  q25=quantile(get(trait),0.25,na.rm=TRUE),
  median=median(get(trait), na.rm=TRUE),
  q75=quantile(get(trait),0.75,na.rm=TRUE),
  max=max(get(trait), na.rm=TRUE)
)]
if(!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive=TRUE)

fwrite(summ, file=file.path(opt$outdir, "pheno_summary.tsv"), sep='\t')

# QC plots
suppressPackageStartupMessages(library(ggplot2))

p1 <- ggplot(ph, aes(x=.data[[trait]])) + geom_histogram(bins=30, fill="#69b3a2") + theme_minimal()
ggsave(file.path(opt$outdir, "qc_hist.png"), p1, width=6, height=4, dpi=150)

if(!is.na(env)){
  p2 <- ggplot(ph, aes(x=.data[[env]], y=.data[[trait]])) + geom_boxplot() + theme_minimal() + theme(axis.text.x=element_text(angle=45,hjust=1))
  ggsave(file.path(opt$outdir, "qc_box.png"), p2, width=7, height=4, dpi=150)
}

# BLUE/BLUP using mixed models
get_formula <- function(type=c("fixed","random")){
  type <- match.arg(type)
  fx <- trimws(opt$fixed)
  rx <- trimws(opt$random)
  rhs <- NULL
  if(type=="fixed"){
    rhs <- fx
  } else {
    rhs <- rx
  }
  rhs <- if(nchar(rhs)>0) rhs else NA
  rhs
}

blue_out <- NULL
blup_out <- NULL

if(opt$blue || opt$blup){
  suppressPackageStartupMessages(library(lme4))
  # Construct formula: trait ~ fixed + (1|random)
  fx <- get_formula("fixed")
  rx <- get_formula("random")
  rhs <- "1"
  if(!is.na(fx)) rhs <- paste(rhs, "+", fx)
  if(!is.na(rx)) rhs <- paste(rhs, "+ (1|", rx, ")")
  frm <- as.formula(paste(trait, "~", rhs))
  fit <- lme4::lmer(frm, data=ph, REML=TRUE)
  # BLUE: estimated marginal means per genotype (id)
  if(opt$blue){
    suppressPackageStartupMessages(library(emmeans))
    emm <- emmeans::emmeans(fit, specs = stats::as.formula(paste("~", id)))
    blue <- as.data.table(as.data.frame(emm))
    setnames(blue, old=names(blue)[1], new="IID")
    setnames(blue, old="emmean", new=trait)
    fwrite(blue[,.(IID,get(trait))], file=file.path(opt$outdir, "blue.tsv"), sep='\t')
    blue_out <- blue
  }
  # BLUP: ranef by genotype if genotype is modeled as random; otherwise fallback to conditional modes for id
  if(opt$blup){
    re <- lme4::ranef(fit, condVar=FALSE)
    # try to find id in random effects
    comp <- names(re)
    candidate <- comp[grep(paste0("^", id, "$"), comp)]
    if(length(candidate)==0 && length(comp)>0) candidate <- comp[1]
    bl <- as.data.table(re[[candidate]])
    bl[, IID := rownames(re[[candidate]])]
    setnames(bl, old="(Intercept)", new=trait)
    fwrite(bl[,.(IID,get(trait))], file=file.path(opt$outdir, "blup.tsv"), sep='\t')
    blup_out <- bl
  }
}

# 输出一个 pheno_processed.tsv 的占位，主流程将根据选择复制
fwrite(ph, file=file.path(opt$outdir, "pheno_processed.tsv"), sep='\t')
