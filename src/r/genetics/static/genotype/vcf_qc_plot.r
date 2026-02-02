#!/usr/bin/env Rscript

library(optparse)
library(ggplot2)
library(reshape2)

option_list <- list(
  make_option(c("--imiss"), type = "character", default = NULL, help = "Input .imiss file (Individual missingness)"),
  make_option(c("--lmiss"), type = "character", default = NULL, help = "Input .lmiss file (Site missingness)"),
  make_option(c("--het"), type = "character", default = NULL, help = "Input .het file (Heterozygosity)"),
  make_option(c("--frq"), type = "character", default = NULL, help = "Input .frq file (Allele frequency)"),
  make_option(c("--depth"), type = "character", default = NULL, help = "Input .idepth file (Individual depth)"),
  make_option(c("--site_depth"), type = "character", default = NULL, help = "Input .ldepth.mean file (Site depth)"),
  make_option(c("-o", "--output"), type = "character", default = "vcf_qc_plots.pdf", help = "Output PDF file")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$imiss) && is.null(opt$lmiss) && is.null(opt$het) && is.null(opt$frq) && is.null(opt$depth) && is.null(opt$site_depth)) {
  print_help(opt_parser)
  stop("At least one input file must be provided", call. = FALSE)
}

pdf(opt$output, width = 10, height = 8)

# 1. Individual Missingness
if (!is.null(opt$imiss)) {
  imiss <- read.table(opt$imiss, header = TRUE, stringsAsFactors = FALSE)
  p1 <- ggplot(imiss, aes(x = F_MISS)) +
    geom_histogram(binwidth = 0.01, fill = "steelblue", color = "black") +
    labs(title = "Individual Missingness", x = "Fraction of Missing Data", y = "Count") +
    theme_minimal()
  print(p1)
}

# 2. Site Missingness
if (!is.null(opt$lmiss)) {
  lmiss <- read.table(opt$lmiss, header = TRUE, stringsAsFactors = FALSE)
  p2 <- ggplot(lmiss, aes(x = F_MISS)) +
    geom_histogram(binwidth = 0.01, fill = "forestgreen", color = "black") +
    labs(title = "Site Missingness", x = "Fraction of Missing Data", y = "Count") +
    theme_minimal()
  print(p2)
}

# 3. Heterozygosity
if (!is.null(opt$het)) {
  het <- read.table(opt$het, header = TRUE, stringsAsFactors = FALSE)
  # Calculate Heterozygous Proportion: (N_SITES - O(HOM)) / N_SITES
  het$HET_PROP <- (het$N_SITES - het$O.HOM.) / het$N_SITES
  
  p3 <- ggplot(het, aes(x = HET_PROP)) +
    geom_histogram(binwidth = 0.01, fill = "purple", color = "black") +
    labs(title = "Individual Heterozygosity", x = "Heterozygous Proportion", y = "Count") +
    theme_minimal()
  print(p3)
  
  p3b <- ggplot(het, aes(y = HET_PROP)) +
    geom_boxplot(fill = "purple", alpha = 0.5) +
    labs(title = "Individual Heterozygosity Boxplot", y = "Heterozygous Proportion") +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  print(p3b)
}

# 4. Allele Frequency (MAF)
if (!is.null(opt$frq)) {
  # VCFtools .frq output format: CHROM POS N_ALLELES N_CHR {ALLELE:FREQ}
  # We need to parse the allele frequencies. Usually the second allele is the minor one if biallelic, 
  # but we should check. For simplicity, we assume biallelic and take the minimum frequency.
  
  # Reading .frq file can be tricky because of variable columns. 
  # We assume standard biallelic output or use readLines if complex.
  # Try reading as table first, skipping header.
  frq <- tryCatch({
      read.table(opt$frq, header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
  }, error = function(e) {
      NULL
  })
  
  if (!is.null(frq)) {
      # Extract frequencies. Columns 5 and 6 usually contain Allele:Freq.
      # Example: A:0.9  G:0.1
      # We parse these.
      
      parse_maf <- function(x) {
          parts <- unlist(strsplit(x, ":"))
          if (length(parts) == 2) return(as.numeric(parts[2]))
          return(NA)
      }
      
      # Assuming biallelic for MAF plot
      if (ncol(frq) >= 6) {
          frq$Freq1 <- sapply(frq[,5], parse_maf)
          frq$Freq2 <- sapply(frq[,6], parse_maf)
          frq$MAF <- pmin(frq$Freq1, frq$Freq2, na.rm = TRUE)
          
          p4 <- ggplot(frq, aes(x = MAF)) +
            geom_histogram(binwidth = 0.01, fill = "orange", color = "black") +
            labs(title = "Minor Allele Frequency (MAF) Distribution", x = "MAF", y = "Count") +
            theme_minimal()
          print(p4)
      }
  }
}

# 5. Individual Depth
if (!is.null(opt$depth)) {
  idepth <- read.table(opt$depth, header = TRUE, stringsAsFactors = FALSE)
  p5 <- ggplot(idepth, aes(x = MEAN_DEPTH)) +
    geom_histogram(binwidth = 1, fill = "firebrick", color = "black") +
    labs(title = "Individual Mean Depth", x = "Mean Depth", y = "Count") +
    theme_minimal()
  print(p5)
}

# 6. Site Mean Depth
if (!is.null(opt$site_depth)) {
  ldepth <- read.table(opt$site_depth, header = TRUE, stringsAsFactors = FALSE)
  p6 <- ggplot(ldepth, aes(x = MEAN_DEPTH)) +
    geom_histogram(binwidth = 1, fill = "darkcyan", color = "black") +
    labs(title = "Site Mean Depth", x = "Mean Depth", y = "Count") +
    theme_minimal()
  print(p6)
  
  if ("VAR_DEPTH" %in% colnames(ldepth)) {
      ldepth$SD_DEPTH <- sqrt(ldepth$VAR_DEPTH)
      p6b <- ggplot(ldepth, aes(x = MEAN_DEPTH, y = SD_DEPTH)) +
        geom_point(alpha = 0.1, color = "darkblue") +
        labs(title = "Site Depth: Mean vs SD", x = "Mean Depth", y = "Standard Deviation") +
        theme_minimal()
      print(p6b)
      
      p6c <- ggplot(ldepth, aes(x = SD_DEPTH)) +
        geom_histogram(binwidth = 1, fill = "cyan", color = "black") +
        labs(title = "Site Depth Standard Deviation", x = "SD Depth", y = "Count") +
        theme_minimal()
      print(p6c)
  }
}

dev.off()
