#!/usr/bin/env Rscript

# Simple PCA plotter from PLINK eigenvec/eigenval
# Usage:
#   Rscript plot_pca_from_eigenvec.R \
#     --eigenvec path/to/file.eigenvec \
#     --eigenval path/to/file.eigenval \
#     --metadata path/to/metadata.(csv|tsv) \
#     --out-prefix PREFIX
# Notes:
# - metadata is optional; if provided, it will merge by FID/IID when possible
# - output: PREFIX.pca_pc12.pdf

suppressWarnings(suppressMessages({
  library(ggplot2)
}))

# Parse args
args <- commandArgs(trailingOnly = TRUE)
parse_args <- function(a) {
  out <- list()
  i <- 1
  while (i <= length(a)) {
    key <- a[i]
    if (startsWith(key, "--")) {
      val <- if (i + 1 <= length(a)) a[i + 1] else NA
      out[[sub("^--", "", key)]] <- val
      i <- i + 2
    } else i <- i + 1
  }
  out
}
opt <- parse_args(args)
stopifnot(!is.null(opt$eigenvec))
if (is.null(opt$out-prefix)) opt$`out-prefix` <- "pca"

# Helpers
read_eigenvec <- function(path) {
  hdrline <- tryCatch(readLines(path, n = 1L), error = function(e) "")
  has_hdr <- grepl("^FID\\s+IID", hdrline)
  df <- tryCatch(read.table(path, header = has_hdr, stringsAsFactors = FALSE), 
                 error = function(e) read.table(path, header = FALSE, stringsAsFactors = FALSE))
  if (!has_hdr) {
    if (ncol(df) >= 2) {
      colnames(df)[1:2] <- c("FID", "IID")
      if (ncol(df) > 2) colnames(df)[3:ncol(df)] <- paste0("PC", seq_len(ncol(df) - 2))
    }
  } else {
    # Normalize column names to PC1, PC2, ... when possible
    if (ncol(df) > 2) {
      pcn <- paste0("PC", seq_len(ncol(df) - 2))
      colnames(df)[3:ncol(df)] <- pcn
    }
  }
  df
}

read_eigenval <- function(path) {
  if (is.null(path) || is.na(path) || path == "") return(NULL)
  tryCatch(scan(path, quiet = TRUE), error = function(e) NULL)
}

read_metadata <- function(path) {
  if (is.null(path) || is.na(path) || path == "") return(NULL)
  if (!file.exists(path)) return(NULL)
  ext <- tolower(sub("^.*\\.(\\w+)$", "\\1", path))
  if (ext %in% c("csv")) {
    tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  } else {
    tryCatch(read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE), error = function(e) NULL)
  }
}

# Load data
E <- read_eigenvec(opt$eigenvec)
V <- read_eigenval(opt$eigenval)
M <- read_metadata(opt$metadata)

# Merge metadata if present
D <- E
if (!is.null(M)) {
  if (all(c("FID","IID") %in% colnames(M))) {
    D <- merge(E, M, by = c("FID","IID"), all.x = TRUE)
  } else if ("FID" %in% colnames(M)) {
    D <- merge(E, M, by = c("FID"), all.x = TRUE)
  } else if ("IID" %in% colnames(M)) {
    D <- merge(E, M, by = c("IID"), all.x = TRUE)
  }
}

# Axis labels from eigenvalues if available
xlab <- "PC1"
ylab <- "PC2"
if (!is.null(V) && length(V) >= 2) {
  prop <- V / sum(V)
  xlab <- sprintf("PC1 (%.2f%%)", prop[1] * 100)
  ylab <- sprintf("PC2 (%.2f%%)", prop[2] * 100)
}

# Choose color column if present
color_col <- NULL
cand <- c("group", "Group", "population", "Population", "ploidy", "lineage")
for (c in cand) if (c %in% colnames(D)) { color_col <- c; break }

p <- ggplot(D, aes(x = .data[["PC1"]], y = .data[["PC2"]])) +
  geom_point(aes(color = if (!is.null(color_col)) .data[[color_col]] else NULL), alpha = 0.8, size = 1.5, na.rm = TRUE) +
  labs(x = xlab, y = ylab, color = ifelse(is.null(color_col), NULL, color_col)) +
  theme_classic()

outfile <- paste0(opt$`out-prefix`, ".pca_pc12.pdf")
pdf(outfile, width = 6, height = 5)
print(p)
dev.off()

cat("Wrote ", outfile, "\n", sep = "")
