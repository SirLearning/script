#!/usr/bin/env Rscript
library(optparse)
library(pheatmap)
library(RColorBrewer)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, help = "Input matrix file (heatmap_format2.txt)"),
  make_option(c("-o", "--output"), type = "character", default = "heatmap.pdf", help = "Output plot file")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file is required", call. = FALSE)
}

# Read data
# The user script reads heatmap_format2.txt with header=T, row.names=1
data <- read.table(opt$input, header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
data_mat <- as.matrix(data)

# Clean column names if they contain "type_" (based on user script logic)
# User script logic: x <- c(x,strsplit(colnames(data), "type_")[[i]][2])
if (any(grepl("type_", colnames(data_mat)))) {
    new_cols <- sapply(colnames(data_mat), function(x) {
        parts <- unlist(strsplit(x, "type_"))
        if (length(parts) > 1) return(parts[2]) else return(x)
    })
    colnames(data_mat) <- new_cols
}

# Define colors (from user script)
# colsA = c("#F7FCFD","#8C96C6","#F7FCF5","#74C476","#FFFFCC","#FD8D3C")
colsA <- c("#F7FCFD", "#8C96C6", "#F7FCF5", "#74C476", "#FFFFCC", "#FD8D3C")

# Plot
pdf(opt$output, width = 10, height = 10)
pheatmap(data_mat, 
         cluster_rows = FALSE, 
         cluster_cols = TRUE, 
         border_color = NA, 
         color = colsA, 
         fontsize = 8)
dev.off()
