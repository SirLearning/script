library(data.table)
library(gridExtra)
library(ggplot2)

# 1. 读取相似度矩阵 (IBS)
ibs_sim <- fread("/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.ibs.mibs")
# ibs_sim <- fread("/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.ibs.raw_sample.mibs")
# 2. 将相似度转为距离 (Distance = 1 - Similarity)
# 检查并替换 NA 值
# 如果 IBS 计算中存在缺失（比如某两个样本没有共同的非缺失位点），IBS 可能为 NA
# 一种处理方法是将 NA 视为完全不相似 (即距离 = 1，或者 IBS = 0)
ibs_sim[is.na(ibs_sim)] <- 0 
dist_matrix <- as.dist(1 - ibs_sim)

# 3. 运行经典 MDS (cmdscale)
# 这将完全基于你 IBS 矩阵里的数值，不涉及任何位点填补
mds_res <- cmdscale(dist_matrix, k = 4)

# 4. 转换为 Dataframe 绘图
mds_df <- data.frame(IID = read.table("/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.ibs.mibs.id")$V2, 
                     C1 = mds_res[,1], 
                     C2 = mds_res[,2],
                     C3 = mds_res[,3],
                     C4 = mds_res[,4])

# 读取分组信息
group_file <- "/data1/dazheng_tusr1/vmap4.VCF.v1/sample_group.txt"
if (file.exists(group_file)) {
    group_data <- read.table(group_file, header=FALSE, col.names=c("IID", "Group"))
    # 去重
    group_data <- group_data[!duplicated(group_data$IID), ]
    
    # 合并 MDS 数据和分组信息
    mds_df <- merge(mds_df, group_data, by="IID", all.x = TRUE)
    # 将未匹配到的样本标记为 "Unknown"
    mds_df$Group[is.na(mds_df$Group)] <- "Unknown"
} else {
    print(paste("Warning: Group file not found:", group_file))
    mds_df$Group <- "Unknown"
}

# 5. 绘图与保存
library(ggplot2)

# Plot C1 vs C2
p1 <- ggplot(mds_df, aes(x = C1, y = C2, color = Group)) +
  geom_point(alpha = 0.7) +
  theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5), 
    legend.position = c(0.95, 0.4), # Legend inside top-right
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = NA)
  ) + 
  labs(title = "MDS based on IBS (C1 vs C2)",
       x = "Coordinate 1",
       y = "Coordinate 2",
       color = "Population")

ggsave("ibs_mds_plot_C1C2_reform.png", plot = p1, width = 8, height = 6, dpi = 300)

# Plot C3 vs C4
p2 <- ggplot(mds_df, aes(x = C3, y = C4, color = Group)) +
  geom_point(alpha = 0.7) +
  theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5), 
    legend.position = c(0.05, 0.4), # Legend inside top-left
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = NA)
  ) + 
  labs(title = "MDS based on IBS (C3 vs C4)",
       x = "Coordinate 3",
       y = "Coordinate 4",
       color = "Population")

ggsave("ibs_mds_plot_C3C4_reform.png", plot = p2, width = 8, height = 6, dpi = 300)

# 组合绘图
combined_plot <- arrangeGrob(p1, p2, ncol = 2)
ggsave("ibs_mds_combined.png", plot = combined_plot, width = 16, height = 6, dpi = 300)
print("Plots saved to ibs_mds_plot_C1C2_reform.png, ibs_mds_plot_C3C4_reform.png and ibs_mds_combined.png")