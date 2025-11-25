#!/usr/bin/env Rscript
library(optparse)
library(reshape2)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, help = "Input table with columns: File, Lineage, Gene"),
  make_option(c("-o", "--output"), type = "character", default = "heatmap_format2.txt", help = "Output matrix file")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file is required", call. = FALSE)
}

# Read input
# Expected format: File <tab> Lineage <tab> Gene
# No header expected in the concatenated file
data <- read.table(opt$input, header = FALSE, stringsAsFactors = FALSE)
colnames(data) <- c("File", "Lineage", "Gene")

# Assign values based on lineage for visualization
# A -> 1, B -> 0.6, D -> 0.2
data$Value <- 0
data$Value[data$Lineage == "A"] <- 1
data$Value[data$Lineage == "B"] <- 0.6
data$Value[data$Lineage == "D"] <- 0.2

# Cast to wide format: Rows = File, Cols = Gene
# If a gene is present, it gets the value. If multiple lineages have it (unlikely for subgenome specific genes but possible), max value is taken.
matrix_data <- dcast(data, File ~ Gene, value.var = "Value", fill = 0, fun.aggregate = max)

write.table(matrix_data, opt$output, sep = "\t", quote = FALSE, row.names = FALSE)
