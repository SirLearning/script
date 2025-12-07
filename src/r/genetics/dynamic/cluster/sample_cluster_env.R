#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(argparser)
    library(cluster)
    library(factoextra)
    library(ggplot2)
    library(maps)
})

p <- arg_parser("Sample Clustering based on Environmental Data")
p <- add_argument(p, "--input", help="Input file (Tab-separated). First column SampleID. Must contain Latitude and Longitude columns.", required=TRUE)
p <- add_argument(p, "--k", help="Number of clusters", type="numeric", default=5)
p <- add_argument(p, "--out", help="Output prefix", default="sample_cluster")
p <- add_argument(p, "--cluster_cols", help="Comma-separated list of columns to use for clustering. If not provided, uses all numeric columns.", default=NULL)

argv <- parse_args(p)

# Read data
# Check if file exists
if (!file.exists(argv$input)) {
    stop(paste("Input file not found:", argv$input))
}

data <- read.table(argv$input, header=TRUE, stringsAsFactors=FALSE, sep="\t")

# Handle ID column (Assume first column is ID)
rownames(data) <- data[,1]
data_clean <- data[,-1]

# Normalize column names for Lat/Lon
colnames(data_clean)[grep("^lat", colnames(data_clean), ignore.case=TRUE)] <- "Latitude"
colnames(data_clean)[grep("^lon|^logi", colnames(data_clean), ignore.case=TRUE)] <- "Longitude"

if(!("Latitude" %in% colnames(data_clean)) | !("Longitude" %in% colnames(data_clean))) {
    stop("Latitude and Longitude columns are required in the input file.")
}

# Select columns for clustering
if(!is.null(argv$cluster_cols)) {
    cols <- strsplit(argv$cluster_cols, ",")[[1]]
    # Check if cols exist
    missing_cols <- setdiff(cols, colnames(data_clean))
    if(length(missing_cols) > 0) {
        stop(paste("Columns not found:", paste(missing_cols, collapse=", ")))
    }
    cluster_data <- data_clean[, cols]
} else {
    # Use all numeric columns
    nums <- unlist(lapply(data_clean, is.numeric))
    cluster_data <- data_clean[, nums]
}

# Remove rows with NA
cluster_data <- na.omit(cluster_data)
common_rows <- rownames(cluster_data)
data_clean <- data_clean[common_rows, ]

cat("Clustering", nrow(cluster_data), "samples using", ncol(cluster_data), "variables.\n")

# Scale
df <- scale(cluster_data, center = TRUE, scale = TRUE)

# Distance
dist_mat <- dist(df, method = "euclidean")

# Hclust
hc <- hclust(d = dist_mat, method = "ward.D2")

# Cut tree
groups <- cutree(hc, k = argv$k)
data_clean$Cluster <- as.factor(groups)

# Calculate mean Lat/Lon for clusters (as in original script)
lat_mean <- tapply(data_clean$Latitude, data_clean$Cluster, mean, na.rm = TRUE)
lon_mean <- tapply(data_clean$Longitude, data_clean$Cluster, mean, na.rm = TRUE)

# Add mean coordinates to data
data_clean$ClusterLat <- lat_mean[data_clean$Cluster]
data_clean$ClusterLon <- lon_mean[data_clean$Cluster]

# Save results
out_df <- cbind(SampleID=rownames(data_clean), data_clean)
write.table(out_df, paste0(argv$out, ".txt"), sep="\t", row.names=FALSE, quote=FALSE)

# Plot Map
mapworld <- borders("world", colour = "gray70", fill="gray90")

p <- ggplot() + 
    mapworld + 
    geom_point(data=data_clean, aes(x=Longitude, y=Latitude, color=Cluster), size=2, alpha=0.7) +
    theme_minimal() +
    labs(title=paste("Sample Clustering (K=", argv$k, ")", sep="")) +
    theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())

pdf(paste0(argv$out, "_map.pdf"), width=10, height=6)
print(p)
dev.off()

# Plot Dendrogram
pdf(paste0(argv$out, "_dendrogram.pdf"), width=12, height=6)
tryCatch({
    print(fviz_dend(hc, k = argv$k, cex = 0.5, 
              k_colors = "jco", rect = TRUE, rect_border = "jco", rect_fill = TRUE,
              main = paste("Dendrogram - Ward.D2 (K=", argv$k, ")", sep="")))
}, error = function(e) {
    cat("Error plotting dendrogram:", e$message, "\n")
})
dev.off()

cat("Analysis complete.\n")
