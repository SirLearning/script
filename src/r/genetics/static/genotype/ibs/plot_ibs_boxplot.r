#!/usr/bin/env Rscript
library(optparse)
library(ggplot2)
library(RColorBrewer)

# 解析命令行参数
option_list <- list(
  make_option(c("-i", "--ibs_files"), type="character", help="Comma-separated list of IBS distance matrix files (e.g. A.mibs,B.mibs)"),
  make_option(c("-l", "--labels"), type="character", default="Genome", help="Comma-separated list of labels for each IBS file (e.g. A,B). If not provided, uses filenames."),
  make_option(c("-d", "--id_files"), type="character", default=NULL, help="Comma-separated list of Sample ID files (e.g. A.mibs.id,B.mibs.id). Required if IBS files have no header. Order must match --ibs_files."),
  make_option(c("-g", "--group_file"), type="character", help="Group information file (Two columns: SampleID Group)"),
  make_option(c("-r", "--ref_sample"), type="character", default="CS", help="Reference sample ID to calculate IBS distance against [default= %default]"),
  make_option(c("-o", "--output"), type="character", default="ibs_boxplot.pdf", help="Output plot filename [default= %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$ibs_files) || is.null(opt$group_file)) {
  stop("Missing required arguments: --ibs_files or --group_file")
}

# 解析文件列表
ibs_file_list <- strsplit(opt$ibs_files, ",")[[1]]
label_list <- if (!is.null(opt$labels) && opt$labels != "Genome") strsplit(opt$labels, ",")[[1]] else basename(ibs_file_list)
id_file_list <- if (!is.null(opt$id_files)) strsplit(opt$id_files, ",")[[1]] else NULL

if (length(label_list) != length(ibs_file_list)) {
    warning("Number of labels does not match number of IBS files. Using filenames as labels.")
    label_list <- basename(ibs_file_list)
}

if (!is.null(id_file_list) && length(id_file_list) != length(ibs_file_list)) {
    stop("Number of ID files must match number of IBS files.")
}

# 2. 读取分组信息
cat("Reading group information...\n")
group_data <- read.table(opt$group_file, header=FALSE, stringsAsFactors=FALSE)
colnames(group_data) <- c("SampleID", "Group")

# 3. 循环读取并合并数据
all_plot_data <- data.frame()

for (i in 1:length(ibs_file_list)) {
    ibs_file <- ibs_file_list[i]
    label <- label_list[i]
    id_file <- if (!is.null(id_file_list)) id_file_list[i] else NULL
    
    cat(paste0("Processing file ", i, ": ", ibs_file, " (Label: ", label, ")...\n"))
    
    # 读取 IBS 矩阵
    if (!is.null(id_file)) {
        ibs_mat <- read.table(ibs_file, header=FALSE)
        ids <- read.table(id_file, header=FALSE, stringsAsFactors=FALSE)
        sample_names <- ids[,2]
        colnames(ibs_mat) <- sample_names
        rownames(ibs_mat) <- sample_names
    } else {
        ibs_mat <- read.table(ibs_file, header=TRUE, row.names=1, check.names=FALSE)
    }
    
    if (!opt$ref_sample %in% colnames(ibs_mat)) {
        warning(paste("Reference sample", opt$ref_sample, "not found in", ibs_file, "- skipping this file."))
        next
    }
    
    # 提取参照样本数据
    ref_vals <- ibs_mat[, opt$ref_sample, drop=FALSE]
    curr_data <- data.frame(SampleID = rownames(ref_vals), Value = ref_vals[,1])
    
    # 合并分组信息
    curr_data <- merge(curr_data, group_data, by="SampleID")
    curr_data$Lineage <- label # 添加 Lineage/Label 列
    
    all_plot_data <- rbind(all_plot_data, curr_data)
}

if (nrow(all_plot_data) == 0) {
    stop("No data available for plotting.")
}

# 4. 绘图
cat("Plotting...\n")

# 设置颜色
n_lineages <- length(unique(all_plot_data$Lineage))
# 使用 Lineage 作为填充色，Group 作为 X 轴
# 如果只有一个文件/Lineage，颜色策略可以回退到按 Group 着色
fill_var <- if (n_lineages > 1) "Lineage" else "Group"

p <- ggplot(all_plot_data, aes(x=Group, y=Value, fill=get(fill_var))) + 
  geom_boxplot(outlier.shape = NA, alpha=0.8) +
  # 如果有多个 Lineage，使用 position_dodge 让它们并排显示
  geom_jitter(position=position_dodge(width=0.75), size=0.5, alpha=0.5) + 
  theme_classic() +
  labs(x = "Group", y = paste("IBS Distance/Similarity with", opt$ref_sample), fill = fill_var) +
  theme(
    axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
    axis.text.y = element_text(size=12)
  )

if (n_lineages > 1) {
    p <- p + scale_fill_brewer(palette = "Set2") # 多 Lineage 使用 Set2
} else {
    n_groups <- length(unique(all_plot_data$Group))
    my_palette <- colorRampPalette(brewer.pal(8, "Accent"))(n_groups)
    p <- p + scale_fill_manual(values = my_palette) + theme(legend.position = "none")
}

ggsave(opt$output, plot=p, width=12, height=6)
cat(paste("Plot saved to", opt$output, "\n"))
