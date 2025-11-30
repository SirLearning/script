#!/usr/bin/env Rscript

library(ggplot2)
library(optparse)

option_list <- list(
  make_option(c("-i", "--input_dir"), type="character", default=NULL, 
              help="Directory containing Tajima's D result files", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="tajima_d_distribution.pdf", 
              help="Output PDF file name", metavar="character"),
  make_option(c("-p", "--pattern"), type="character", default="*.Tajima.D", 
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
  # VCFtools .Tajima.D: CHROM BIN_START N_SNPS TajimaD
  df <- read.table(x, header=TRUE, stringsAsFactors=FALSE)
  
  fname <- basename(x)
  fname <- gsub("\\.Tajima\\.D$", "", fname)
  
  df$Group <- fname
  return(df)
})

all_data <- do.call(rbind, data_list)

# Ensure TajimaD column exists
if (!"TajimaD" %in% colnames(all_data)) {
    stop("Column 'TajimaD' not found in input files.")
}

pdf(file = opt$output, width=10, height=6)

# Boxplot
p1 <- ggplot(all_data, aes(x=Group, y=TajimaD, fill=Group)) + 
  geom_boxplot() + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title="Tajima's D Distribution", y="Tajima's D")
print(p1)

# Density plot
p2 <- ggplot(all_data, aes(x=TajimaD, color=Group, fill=Group)) + 
  geom_density(alpha=0.3) + 
  theme_bw() +
  labs(title="Tajima's D Density", x="Tajima's D")
print(p2)

dev.off()

print(paste("Plot saved to", opt$output))
