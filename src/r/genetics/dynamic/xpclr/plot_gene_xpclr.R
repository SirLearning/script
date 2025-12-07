#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(argparser)
    library(ggplot2)
})

p <- arg_parser("Plot XP-CLR scores for specific genes/regions")
p <- add_argument(p, "--xpclr_file", help="XP-CLR output file (with header: Chr, WindowStart, WindowStop, SNPcount, MeanY, Wstat)", required=TRUE)
p <- add_argument(p, "--gene_pos", help="Gene position (bp)", type="numeric", required=TRUE)
p <- add_argument(p, "--gene_name", help="Gene name for title", default="Gene")
p <- add_argument(p, "--chr_name", help="Chromosome name", default="Chr")
p <- add_argument(p, "--threshold", help="Significance threshold line", type="numeric", default=NULL)
p <- add_argument(p, "--centromere", help="Centromere position (bp)", type="numeric", default=NULL)
p <- add_argument(p, "--out", help="Output PDF file", default="xpclr_plot.pdf")

argv <- parse_args(p)

# Read Data
data <- read.table(argv$xpclr_file, header=TRUE, stringsAsFactors=FALSE)

# Plot
pdf(argv$out, width = 10, height = 5)

# Basic plot
# X axis in Mb
x_mb <- data$WindowStart / 1000000
y_val <- data$MeanY

# Determine Y limits to make it look nice
y_max <- max(c(y_val, if(!is.null(argv$threshold)) argv$threshold else 0), na.rm=TRUE) * 1.1
y_min <- min(y_val, na.rm=TRUE)
if (y_min > 0) y_min <- -2 # Give some space at bottom like original script

plot(x_mb, y_val, 
     ylim = c(y_min, y_max), 
     axes = TRUE, 
     col = "skyblue", 
     type = "l", 
     cex = 2, 
     lwd = 3, 
     ann = FALSE)

# Add Gene Line
abline(v = argv$gene_pos / 1000000, col = "red", cex=2, lwd = 1)

# Add Threshold Line
if (!is.null(argv$threshold)) {
    abline(h = argv$threshold, col = "red", lty = 3, cex=2, lwd = 2)
}

# Add Centromere Point
if (!is.null(argv$centromere)) {
    points(argv$centromere / 1000000, 0, pch=20, cex=2, col="grey")
}

# Titles
title(main = paste(argv$gene_name, "on", argv$chr_name), 
      xlab = paste(argv$chr_name, "(Mb)"), 
      ylab = 'XP-CLR Score', 
      font.lab = 1, 
      cex.lab = 1.5)

dev.off()

cat("Plot saved to", argv$out, "\n")
