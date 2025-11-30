#!/usr/bin/env Rscript

library(ggplot2)
library(stringr)
library(optparse)

option_list <- list(
  make_option(c("-i", "--input_dir"), type="character", default=NULL, 
              help="Directory containing FST result files", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="fst_boxplot.pdf", 
              help="Output PDF file name", metavar="character"),
  make_option(c("-p", "--pattern"), type="character", default="*.windowed.weir.fst", 
              help="File pattern to match", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_dir)){
  print_help(opt_parser)
  stop("Input directory must be supplied", call.=FALSE)
}

file_list <- list.files(path = opt$input_dir, pattern = opt$pattern, full.names = TRUE)

if (length(file_list) == 0) {
  stop("No files found matching the pattern in the directory.")
}

data_list <- lapply(file_list, function(x) {
  df <- read.table(x, header=TRUE, stringsAsFactors=FALSE)
  # Assume filename is like group1_group2.suffix
  fname <- basename(x)
  # Remove extension
  fname <- sub(opt$pattern, "", fname, fixed=TRUE) # Simple removal if pattern matches end
  # If pattern is glob, this might not work perfectly for stripping. 
  # Let's just strip known extensions
  fname <- gsub("\\.windowed\\.weir\\.fst$", "", fname)
  fname <- gsub("\\.fst$", "", fname)
  
  df$Comparison <- fname
  return(df)
})

all_data <- do.call(rbind, data_list)

# Check for column names. VCFtools output usually has WEIGHTED_FST and MEAN_FST
if (!"WEIGHTED_FST" %in% colnames(all_data)) {
    # Try to find a column that looks like FST
    fst_col <- grep("fst", colnames(all_data), ignore.case=TRUE, value=TRUE)
    if(length(fst_col) > 0) {
        all_data$WEIGHTED_FST <- all_data[[fst_col[1]]]
    } else {
        stop("Could not find WEIGHTED_FST or similar column in input files.")
    }
}

pdf(file = opt$output, width=10, height=6)
p <- ggplot(all_data, aes(x=Comparison, y=WEIGHTED_FST, fill=Comparison)) + 
  geom_boxplot() + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title="FST Distribution", y="Weighted FST")
print(p)
dev.off()

print(paste("Plot saved to", opt$output))
