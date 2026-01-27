library(ggplot2)

# 读取数据
# PLINK2 的 .eigenvec 文件可能会有表头（以 # 开头），comment.char="" 防止 # 被当做注释忽略
# 第一列通常是 IID (新版 PLINK2) 或者 FID IID (旧版 PLINK 1.9)
# 根据报错，数据只有 11 列，说明只有 IID 而没有 FID，或者 FID 和 IID 合并了
# pca_data <- read.table("/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.pca.eigenvec", header = FALSE, comment.char = "#")
pca_data <- read.table("/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.pca.raw_sample.eigenvec", header = FALSE, comment.char = "#")

# 根据实际列数来命名
if (ncol(pca_data) == 11) {
    colnames(pca_data) <- c("IID", paste0("PC", 1:10))
} else if (ncol(pca_data) == 12) {
    colnames(pca_data) <- c("FID", "IID", paste0("PC", 1:10))
} else {
    print(paste("Warning: Unexpected number of columns:", ncol(pca_data)))
    # 尝试自动命名
    colnames(pca_data) <- c("ID", paste0("PC", 1:(ncol(pca_data)-1)))
}

# 计算 PVE (提前计算以用于标签)
# eigenvals <- scan("/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.pca.eigenval")
eigenvals <- scan("/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.pca.raw_sample.eigenval")
pve <- eigenvals / sum(eigenvals) * 100
# 查看前两个 PC 的解释率
print(pve[1:4])

# 读取分组信息
group_file <- "/data1/dazheng_tusr1/vmap4.VCF.v1/sample_group.txt"
if (file.exists(group_file)) {
    group_data <- read.table(group_file, header=FALSE, col.names=c("IID", "Group"))
    # 去重：如果样本出现在多个组，保留第一个
    group_data <- group_data[!duplicated(group_data$IID), ]
    
    # 合并 PCA 数据和分组信息
    pca_data <- merge(pca_data, group_data, by="IID", all.x = TRUE)
    # 将未匹配到的样本标记为 "Unknown"
    pca_data$Group[is.na(pca_data$Group)] <- "Unknown"
} else {
    print(paste("Warning: Group file not found:", group_file))
    pca_data$Group <- "Unknown"
}

# 绘制 PC1 vs PC2
p <- ggplot(pca_data, aes(x = PC1, y = -PC2, color = Group)) + # PC2 * -1
  geom_point(alpha = 0.7) + 
  theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5), 
    legend.position = c(0.95, 0.7), # Legend inside top-right
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = NA)
  ) + 
  labs(title = "PCA of Genomic Data",
       x = sprintf("PC1 (%.2f%%)", pve[1]),
       y = sprintf("PC2 * -1 (%.2f%%)", pve[2]),
       color = "Population") # Legend 标题

# 保存为 PNG
ggsave("pca_plot_reform.png", plot = p, width = 8, height = 6, dpi = 300)
print("Plot saved to pca_plot_reform.png")
p1 <- p

# 绘制 PC3 vs PC4
p <- ggplot(pca_data, aes(x = -PC3, y = -PC4, color = Group)) +
  geom_point(alpha = 0.7) +
  theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5), 
    legend.position = c(0.15, 0.95), # Legend inside top-right
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = NA)
  ) + 
  labs(title = "PCA of Genomic Data",
       x = sprintf("PC3 * -1 (%.2f%%)", pve[3]),
       y = sprintf("PC4 * -1 (%.2f%%)", pve[4]),
       color = "Population")

# 保存为 PNG
ggsave("pca_plot_3_4_reform.png", plot = p, width = 8, height = 6, dpi = 300)
print("Plot saved to pca_plot_3_4_reform.png")
p2 <- p

# 组合绘图
# Check if gridExtra is installed, if not, try to use just standard par (not easy with ggplot) or suggest install
# However, in this environment we assume it might be present or we can install it.
if (!require("gridExtra")) install.packages("gridExtra", repos="http://cran.us.r-project.org")
library(gridExtra)

combined_plot <- arrangeGrob(p1, p2, ncol = 2)
ggsave("pca_combined.png", plot = combined_plot, width = 16, height = 6, dpi = 300)
print("Plot saved to pca_combined.png")

# 绘制 PC1 vs PC2 (Zoomed PC1: -0.015 to 0.055)
p_zoom <- ggplot(pca_data, aes(x = PC1, y = -PC2, color = Group)) +
  geom_point(alpha = 0.7) + 
  theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5), 
    legend.position = c(0.95, 0.5), # Legend inside top-right (adjust if needed)
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = NA)
  ) + 
  labs(title = "PCA of Genomic Data (Zoomed)",
       x = sprintf("PC1 (%.2f%%)", pve[1]),
       y = sprintf("PC2 * -1 (%.2f%%)", pve[2]),
       color = "Population") +
  coord_cartesian(xlim = c(-0.015, 0.055))

ggsave("pca_plot_PC1_zoom.png", plot = p_zoom, width = 8, height = 6, dpi = 300)
print("Plot saved to pca_plot_PC1_zoom.png")