#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(argparser)
    library(dplyr)
    library(reshape2)
    library(RColorBrewer)
})

p <- arg_parser("Plot Allele Frequency by Latitude/Group")
p <- add_argument(p, "--geno", help="Genotype file (.012 from vcftools)", required=TRUE)
p <- add_argument(p, "--indv", help="Individual file (.012.indv from vcftools)", required=TRUE)
p <- add_argument(p, "--pos", help="Position file (.012.pos from vcftools)", required=TRUE)
p <- add_argument(p, "--meta", help="Metadata file (Header: SampleID, Latitude, ...)", required=TRUE)
p <- add_argument(p, "--out", help="Output prefix", default="allele_trend")
p <- add_argument(p, "--snp_idx", help="Index of SNP to analyze (1-based, default all if not too many)", default=NULL)

argv <- parse_args(p)

# Read Data
cat("Reading individuals...\n")
indv <- read.table(argv$indv, stringsAsFactors=FALSE)[,1]

cat("Reading positions...\n")
pos <- read.table(argv$pos, stringsAsFactors=FALSE)
# pos columns: CHROM, POS

cat("Reading genotypes...\n")
# .012 file has no header, rows are samples, cols are SNPs. First col is row index (skip it)
geno <- read.table(argv$geno, stringsAsFactors=FALSE)
geno <- geno[,-1] # Remove index column
rownames(geno) <- indv
colnames(geno) <- paste(pos[,1], pos[,2], sep="_")

cat("Reading metadata...\n")
meta <- read.table(argv$meta, header=TRUE, stringsAsFactors=FALSE)
# Assume first column is SampleID
colnames(meta)[1] <- "SampleID"
rownames(meta) <- meta$SampleID

# Filter samples
common_samples <- intersect(indv, meta$SampleID)
cat("Number of common samples:", length(common_samples), "\n")

if(length(common_samples) == 0) stop("No common samples between genotype and metadata!")

geno <- geno[common_samples, , drop=FALSE]
meta <- meta[common_samples, , drop=FALSE]

# Function to process one SNP
process_snp <- function(snp_name, snp_values, meta_df) {
    df <- meta_df
    df$Genotype <- snp_values
    
    # Filter missing
    df <- df[df$Genotype != -1 & !is.na(df$Latitude), ]
    
    # Binning Logic
    df$RegionGroup <- NA
    
    # Try numeric
    lat_numeric <- suppressWarnings(as.numeric(df$Latitude))
    
    # Numeric bins
    df$RegionGroup[which(lat_numeric > 20 & lat_numeric <= 30)] <- "20-30"
    df$RegionGroup[which(lat_numeric > 30 & lat_numeric <= 40)] <- "30-40"
    df$RegionGroup[which(lat_numeric > 40 & lat_numeric <= 50)] <- "40-50"
    df$RegionGroup[which(lat_numeric > 50)] <- ">50"
    
    # Categorical/Special bins (e.g. DD, AABB)
    # If Latitude was not numeric (NA in lat_numeric), use original value
    is_non_numeric <- is.na(lat_numeric) & !is.na(df$Latitude)
    if(any(is_non_numeric)) {
        df$RegionGroup[is_non_numeric] <- df$Latitude[is_non_numeric]
    }
    
    df <- df[!is.na(df$RegionGroup), ]
    
    # Order levels
    # Custom order: Ancestors first, then Latitudes
    known_levels <- c("DD", "AABB", "20-30", "30-40", "40-50", ">50")
    other_levels <- setdiff(unique(df$RegionGroup), known_levels)
    final_levels <- c(known_levels[known_levels %in% unique(df$RegionGroup)], other_levels)
    
    df$RegionGroup <- factor(df$RegionGroup, levels = final_levels)
    
    # Aggregate
    # Genotype is 0, 1, 2. Frequency of Alt allele = Sum(GT) / (2 * N)
    # Or just Mean(GT) / 2
    agg <- df %>% group_by(RegionGroup) %>% 
        summarise(
            Freq = mean(Genotype) / 2,
            Count = n()
        )
    
    return(agg)
}

# Loop over SNPs
pdf(paste0(argv$out, ".pdf"), width=8, height=6)

snps_to_plot <- if(!is.na(argv$snp_idx)) colnames(geno)[as.numeric(argv$snp_idx)] else colnames(geno)

for(snp in snps_to_plot) {
    cat("Processing SNP:", snp, "\n")
    agg_res <- process_snp(snp, geno[,snp], meta)
    
    p <- ggplot(agg_res, aes(x=RegionGroup, y=Freq, group=1)) +
        geom_line(color="steelblue", size=1) +
        geom_point(aes(size=Count), color="red") +
        theme_minimal() +
        labs(title=paste("Allele Frequency Trend:", snp),
             x="Latitude / Group", y="Alt Allele Frequency") +
        ylim(0, 1)
    
    print(p)
}

dev.off()
cat("Done.\n")
