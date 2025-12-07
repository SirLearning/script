#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(argparser)
    library(RColorBrewer)
    library(maps) # For world map borders
})

p <- arg_parser("Plot Sample Distribution on World Map")
p <- add_argument(p, "--meta", help="Metadata file (Header: SampleID, Latitude, Longitude, Region/Group)", required=TRUE)
p <- add_argument(p, "--out", help="Output prefix", default="sample_map")
p <- add_argument(p, "--group_col", help="Column name for grouping/coloring", default="Region")

argv <- parse_args(p)

cat("Reading metadata...\n")
meta <- read.table(argv$meta, header=TRUE, stringsAsFactors=FALSE)

# Check columns
req_cols <- c("Latitude", "Longitude", argv$group_col)
missing_cols <- setdiff(req_cols, colnames(meta))
if(length(missing_cols) > 0) {
    stop(paste("Missing columns in metadata:", paste(missing_cols, collapse=", ")))
}

# Filter invalid coordinates
df <- meta[!is.na(meta$Latitude) & !is.na(meta$Longitude), ]
cat("Plotting", nrow(df), "samples with valid coordinates.\n")

# World map
mapworld <- borders("world", colour = "gray70", fill="white")

# Colors
n_groups <- length(unique(df[[argv$group_col]]))
if(n_groups <= 8) {
    colors <- brewer.pal(max(3, n_groups), "Dark2")[1:n_groups]
} else {
    colors <- rainbow(n_groups)
}

p <- ggplot() + 
    mapworld + 
    geom_point(data=df, aes_string(x="Longitude", y="Latitude", color=argv$group_col), alpha=0.6, size=2) +
    scale_color_manual(values = colors) +
    theme_classic() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.border = element_blank()) +
    labs(title = "Sample Distribution")

pdf(paste0(argv$out, ".pdf"), width=12, height=8)
print(p)
dev.off()

cat("Map saved to", paste0(argv$out, ".pdf"), "\n")
