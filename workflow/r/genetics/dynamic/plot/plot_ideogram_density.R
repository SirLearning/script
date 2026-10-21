#!/usr/bin/env Rscript

library(RIdeogram)
library(optparse)

option_list <- list(
  make_option(c("-k", "--karyotype"), type="character", default=NULL, 
              help="Karyotype file (Chr Start End)", metavar="character"),
  make_option(c("-d", "--density"), type="character", default=NULL, 
              help="Density file (Chr Start End Value)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="ideogram.svg", 
              help="Output SVG file name", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$karyotype) || is.null(opt$density)){
  print_help(opt_parser)
  stop("Both karyotype and density files must be supplied", call.=FALSE)
}

# Read data
karyotype <- read.table(opt$karyotype, header=TRUE, stringsAsFactors=FALSE)
density <- read.table(opt$density, header=TRUE, stringsAsFactors=FALSE)

# Check columns
# RIdeogram expects: Chr, Start, End for karyotype
# Chr, Start, End, Value for density (Value column name can vary but usually 4th)

# Plot
ideogram(karyotype = karyotype, overlaid = density)
convertSVG("chromosome.svg", device = "png", file = opt$output)
# Note: RIdeogram outputs "chromosome.svg" by default in current dir. 
# convertSVG converts it. 
# If we want to rename the svg itself:
if (file.exists("chromosome.svg")) {
    file.rename("chromosome.svg", opt$output)
}

print(paste("Plot saved to", opt$output))
