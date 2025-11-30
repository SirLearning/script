#!/usr/bin/env Rscript
library(optparse)
library(ggplot2)
library(ggExtra)
library(dplyr)
library(aplot)

# 解析命令行参数
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, help = "Input depth file (columns: Depth, SD)"),
  make_option(c("-o", "--output"), type = "character", default = "depth_sd_plot.pdf", help = "Output plot file"),
  make_option(c("-m", "--method"), type = "character", default = "marginal", help = "Plot method: 'marginal' (ggExtra) or 'aplot' (custom layout)"),
  make_option(c("-n", "--sample_size"), type = "integer", default = 10000, help = "Number of points to sample for plotting (default: 10000)"),
  make_option(c("--xlim"), type = "double", default = 15, help = "X-axis limit (Depth)"),
  make_option(c("--ylim"), type = "double", default = 10, help = "Y-axis limit (SD)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file is required", call. = FALSE)
}

# 读取数据
# 假设输入文件没有表头，第一列是 Depth，第二列是 SD
data <- read.table(opt$input, header = FALSE, stringsAsFactors = FALSE)
colnames(data) <- c("Depth", "SD")

# 数据过滤与抽样
data <- data %>%
  filter(Depth > 0 & Depth <= opt$xlim & SD <= opt$ylim)

if (nrow(data) > opt$sample_size) {
  data <- data %>% sample_n(opt$sample_size)
}

# 绘图
if (opt$method == "marginal") {
  # Method 1: ggExtra
  p <- ggplot(data, aes(x = Depth, y = SD)) +
    geom_point(alpha = 0.3, size = 0.8, color = "steelblue") +
    xlab("Mean Depth") + ylab("Standard Deviation (SD)") +
    xlim(0, opt$xlim) + ylim(0, opt$ylim) +
    theme_classic() +
    theme(legend.position = 'none')

  p_final <- ggMarginal(p, type = "density", fill = "grey", alpha = 0.5)
  
  ggsave(opt$output, plot = p_final, width = 6, height = 6)
  
} else if (opt$method == "aplot") {
  # Method 2: aplot
  p1 <- ggplot(data, aes(x = Depth, y = SD)) +
    geom_point(alpha = 0.1, color = "steelblue") +
    theme_bw() +
    xlab("Mean Depth") + ylab("SD") +
    xlim(0, opt$xlim) + ylim(0, opt$ylim)

  p2 <- ggplot(data, aes(Depth)) +
    geom_density(fill = "grey", alpha = 0.5) +
    scale_y_continuous(expand = c(0, 0)) +
    xlim(0, opt$xlim) +
    theme_minimal() +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

  p3 <- ggplot(data, aes(SD)) +
    geom_density(fill = "grey", alpha = 0.5) +
    scale_y_continuous(expand = c(0, 0)) +
    xlim(0, opt$ylim) +
    theme_minimal() +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
    coord_flip()

  p_final <- p1 %>%
    insert_top(p2, height = 0.3) %>%
    insert_right(p3, width = 0.3)
    
  ggsave(opt$output, plot = p_final, width = 6, height = 6)
} else {
  stop("Invalid method. Choose 'marginal' or 'aplot'.")
}

cat(paste("Plot saved to", opt$output, "\n"))
