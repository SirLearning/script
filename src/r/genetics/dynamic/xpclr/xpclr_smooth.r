#!/usr/bin/env Rscript
library(optparse)
library(GenWin)
library(dplyr)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, help = "Input XP-CLR output file"),
  make_option(c("-o", "--output"), type = "character", default = NULL, help = "Output smoothed file"),
  make_option(c("-p", "--plot"), type = "character", default = NULL, help = "Output PDF plot file (optional)"),
  make_option(c("-s", "--smoothness"), type = "integer", default = 2000, help = "Smoothness parameter for splineAnalyze"),
  make_option(c("-m", "--method"), type = "integer", default = 2, help = "Method for splineAnalyze (default: 2)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("Input and output files are required", call. = FALSE)
}

# Read data
# Assuming space separated based on user script
# Columns: likely V1..Vn. User script uses V6 (score) and V4 (map position)
data <- read.table(opt$input, sep = " ", header = FALSE, stringsAsFactors = FALSE)

# Filter NA and Inf
# User script: Data2<-na.omit(Data1); Data <- filter(Data2, Data2[,6] != Inf)
data <- na.omit(data)
data <- data[data[,6] != Inf, ]

if (nrow(data) == 0) {
    warning("No valid data points found in input file.")
    # Create empty output file with header
    write.table(data.frame(WindowStart=numeric(), WindowStop=numeric(), SNPcount=numeric(), MeanY=numeric(), Wstat=numeric()), 
                opt$output, sep="\t", col.names=TRUE, row.names=FALSE)
    quit(save="no")
}

# Normalize (Z-score)
# User script: Z[j,1]=(Data[j,6]-mean(Data$V6))/sd(Data$V6)
z_scores <- (data[,6] - mean(data[,6])) / sd(data[,6])
z_matrix <- matrix(z_scores, ncol = 1)

# Spline Analyze
# User script: NORM=splineAnalyze(Y=Z,map=Data$V4,smoothness = 2000,plotRaw=F,plotWindows = T,method = 2)
if (!is.null(opt$plot)) {
    pdf(opt$plot, width = 12, height = 8)
    plot_windows <- TRUE
} else {
    plot_windows <- FALSE
}

tryCatch({
    norm <- splineAnalyze(Y = z_matrix, map = data[,4], smoothness = opt$smoothness, plotRaw = FALSE, plotWindows = plot_windows, method = opt$method)
    
    if (!is.null(opt$plot)) {
        dev.off()
    }

    # Write output
    # User script: normScore=NORM$windowData
    write.table(norm$windowData, opt$output, sep = "\t", col.names = TRUE, row.names = FALSE)

}, error = function(e) {
    print(paste("Error in splineAnalyze:", e$message))
    if (!is.null(opt$plot)) { dev.off() }
    stop(e)
})
