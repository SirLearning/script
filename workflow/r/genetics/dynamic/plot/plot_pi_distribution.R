#!/usr/bin/env Rscript

library(ggplot2)
library(optparse)

option_list <- list(
  make_option(c("-i", "--input_dir"), type="character", default=NULL, 
              help="Directory containing Pi result files", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="pi_distribution.pdf", 
              help="Output PDF file name", metavar="character"),
  make_option(c("-p", "--pattern"), type="character", default="*.sites.pi", 
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
  # VCFtools .sites.pi: CHROM POS PI
  # VCFtools .windowed.pi: CHROM BIN_START BIN_END N_VARIANTS PI
  df <- read.table(x, header=TRUE, stringsAsFactors=FALSE)
  
  fname <- basename(x)
  fname <- gsub("\\.sites\\.pi$", "", fname)
  fname <- gsub("\\.windowed\\.pi$", "", fname)
  fname <- gsub("\\.pi$", "", fname)
  
  df$Group <- fname
  return(df)
})

all_data <- do.call(rbind, data_list)

# Ensure PI column exists
if (!"PI" %in% colnames(all_data)) {
    stop("Column 'PI' not found in input files.")
}

pdf(file = opt$output, width=10, height=6)

# Boxplot
p1 <- ggplot(all_data, aes(x=Group, y=PI, fill=Group)) + 
  geom_boxplot() + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title="Nucleotide Diversity (Pi) Distribution", y="Pi")
print(p1)

# Density plot
p2 <- ggplot(all_data, aes(x=PI, color=Group, fill=Group)) + 
  geom_density(alpha=0.3) + 
  theme_bw() +
  labs(title="Nucleotide Diversity (Pi) Density", x="Pi")
print(p2)

dev.off()

print(paste("Plot saved to", opt$output))
