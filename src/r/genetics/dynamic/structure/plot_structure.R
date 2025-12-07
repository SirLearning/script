#!/usr/bin/env Rscript

# Script to analyze and plot population structure results (e.g., from ADMIXTURE)
# Based on pophelper package

suppressPackageStartupMessages({
    library(pophelper)
    library(ggplot2)
    library(gridExtra)
    library(reshape2)
    library(pheatmap)
    library(argparser)
    library(RColorBrewer)
})

# Parse arguments
p <- arg_parser("Analyze and plot population structure results (Q matrices)")
p <- add_argument(p, "--input_dir", help="Directory containing Q matrix files", default=".")
p <- add_argument(p, "--pattern", help="Pattern to match Q files (e.g., '*.Q')", default="*.Q")
p <- add_argument(p, "--output_dir", help="Output directory", default=".")
p <- add_argument(p, "--prefix", help="Output prefix", default="structure_analysis")

argv <- parse_args(p)

input_dir <- argv$input_dir
pattern <- argv$pattern
output_dir <- argv$output_dir
prefix <- argv$prefix

if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

cat("Searching for Q files in:", input_dir, "with pattern:", pattern, "\n")
sfiles <- list.files(path = input_dir, pattern = pattern, full.names = TRUE)

if (length(sfiles) == 0) {
    stop("No Q files found!")
}

cat("Found", length(sfiles), "files.\n")

# Read Q files
slist <- readQ(files = sfiles)

# Tabulate Q
cat("Tabulating Q matrices...\n")
tr1 <- tabulateQ(qlist = slist)
write.table(tr1, file = file.path(output_dir, paste0(prefix, ".tabulateQ.txt")), sep = "\t", quote = FALSE, row.names = FALSE)

# Summarise Q
cat("Summarising Q matrices...\n")
sr1 <- summariseQ(tr1)
write.table(sr1, file = file.path(output_dir, paste0(prefix, ".summariseQ.txt")), sep = "\t", quote = FALSE, row.names = FALSE)

# Plotting
cat("Generating plots...\n")

# Basic barplots
# Plotting all runs
pdf(file.path(output_dir, paste0(prefix, ".plots.pdf")), width = 12, height = 8)
tryCatch({
    # Enhanced plotQ with sorting and colors
    # Use Set2 or Set3 for distinct colors
    # sortind="all" sorts individuals by cluster membership
    p1 <- plotQ(slist, exportplot = FALSE, returnplot = TRUE, imgoutput = "join", 
                basesize = 11, showyaxis = TRUE, showticks = TRUE,
                sortind = "all", showindlab = FALSE, sharedindlab = TRUE,
                clustercol = brewer.pal(max(8, length(slist)), "Set2"))
    
    # If multiple plots, arrange them
    if (length(p1$plot) > 1) {
        do.call(grid.arrange, c(p1$plot, ncol=1))
    } else {
        grid.arrange(p1$plot[[1]])
    }
}, error = function(e) {
    cat("Error plotting joined Q-plots: ", e$message, "\n")
})
dev.off()

# Evanno method if possible (needs multiple runs per K usually, but pophelper can try)
# If we have multiple K, we can try to plot Delta K if applicable, 
# but standard ADMIXTURE output usually needs CV error for that.
# Here we just plot the Q matrices.

cat("Analysis complete. Results saved to:", output_dir, "\n")
