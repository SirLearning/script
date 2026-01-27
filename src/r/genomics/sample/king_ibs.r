library(ggplot2)
library(scales)

# 读取数据
df <- read.table("/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.king.kin0", header = TRUE, comment.char = "")

# 统计 Mean 和 Median
mean_val <- mean(df$KINSHIP)
median_val <- median(df$KINSHIP)
print(paste("Mean:", mean_val, "Median:", median_val))

# 1. 基础散点图美化
p_base <- ggplot(df, aes(x = IBS0, y = KINSHIP)) +
  geom_point(alpha = 0.4, color = "darkblue", size = 0.5) +
  geom_hline(yintercept = c(0.0442, 0.0884, 0.177, 0.354), 
             linetype = "dashed", color = "red", alpha=0.5) +
  annotate("text", x = 0.05, y = 0.36, label = "Duplicate/MZ Twin", color = "red", hjust=0, size=3) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Kinship vs. IBS0 Scatter Plot",
       x = "Fraction of Zero IBS Sites (IBS0)",
       y = "KING Kinship Coefficient")

ggsave("king_ibs_plot_reform.png", plot = p_base, width = 8, height = 6, dpi = 300)
print("Plot saved to king_ibs_plot_reform.png")


# 2. 截断Y轴图 (-1 到 0.5) + 小于 -1 的点 + 阈值区间标注
# 分离数据
df_main <- df[df$KINSHIP >= -1, ]
df_low <- df[df$KINSHIP < -1, ]

p_cutoff <- ggplot() +
  # 主体点
  geom_point(data = df_main, aes(x = IBS0, y = KINSHIP), 
             alpha = 0.4, color = "darkblue", size = 0.5) +
  # 小于 -1 的点
  geom_point(data = df_low, aes(x = IBS0, y = -1.05),
             alpha = 0.1, color = "darkred", size = 0.2, shape = 124) +
  
  # 辅助线
  geom_hline(yintercept = c(0.0442, 0.0884, 0.177, 0.354), 
             linetype = "dashed", color = "red", alpha=0.5) +
  geom_hline(yintercept = -1, linetype = "solid", color = "gray") +
  
  # 阈值区间标注
  annotate("text", x = 0.002, y = 0.4, label = ">0.354: Duplicate", color = "red", hjust=0, size=2.5) +
  annotate("text", x = 0.002, y = 0.26, label = "0.177-0.354: 1st Degree", color = "red", hjust=0, size=2.5) +
  annotate("text", x = 0.002, y = 0.13, label = "0.0884-0.177: 2nd Degree", color = "red", hjust=0, size=2.5) +
  annotate("text", x = 0.002, y = 0.065, label = "0.0442-0.0884: 3rd Degree", color = "red", hjust=0, size=2.5) +
  
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(-1.1, 0.5), breaks = c(-1, 0, 0.0442, 0.0884, 0.177, 0.354, 0.5)) +
  labs(title = "Kinship vs. IBS0 (Zoomed -1 to 0.5)",
       subtitle = "Intervals annotated; points < -1 shown at y=-1.05",
       x = "Fraction of Zero IBS Sites (IBS0)",
       y = "KING Kinship Coefficient")

# 保存为 PNG (推荐用于大量点数据的压缩)
# 如果数据量巨大，CairoPNG 或 raster graphics 会更快，且文件更小
# 这里使用 type="cairo" 如果可用，可以提升速度和质量
ggsave("king_ibs_plot_cutoff.png", plot = p_cutoff, width = 8, height = 6, dpi = 300)
print("Plot saved to king_ibs_plot_cutoff.png")


# --- Helper to add Mean/Median to Legend ---
add_lines_legend <- function(p) {
  p + geom_vline(aes(xintercept = mean_val, color="Mean"), linetype="dashed", linewidth=0.8) +
      geom_vline(aes(xintercept = median_val, color="Median"), linetype="dotted", linewidth=1.0) +
      scale_color_manual(name = "Statistics", values = c("Mean" = "blue", "Median" = "orange"),
                         labels = c(sprintf("Mean: %.4f", mean_val), sprintf("Median: %.4f", median_val))) +
      theme(legend.position = c(0.1, 0.9), legend.background = element_rect(fill = NA))
}

# 3. KING值分布图 (完整范围)
p_dist_full <- ggplot(df, aes(x = KINSHIP)) +
  geom_histogram(bins = 100, fill = "forestgreen", color = "white", linewidth=0.1) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Distribution of KING Kinship (Full Range)",
       x = "KING Kinship Coefficient",
       y = "Count")
p_dist_full <- add_lines_legend(p_dist_full)

ggsave("king_dist_full.png", plot = p_dist_full, width = 8, height = 6, dpi = 300)
print("Plot saved to king_dist_full.png")

# 4. KING值分布图 (完整范围, Log Scale Y) - 使用 scale_y_log10(labels = scales::comma)
# 解决 log(0) 问题：加一个小常数或者只显示 >0 的 bin
# ggplot 处理 log scalse时会自动丢弃 <=0 的值，这是Warning的原因，通常可以忽略
# labels = scales::comma 需要 library(scales)

p_dist_full_log <- ggplot(df, aes(x = KINSHIP)) +
  geom_histogram(bins = 100, fill = "forestgreen", color = "white", linewidth=0.1) +
  scale_y_log10(labels = comma) + # 只要普通的数字显示
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Distribution of KING Kinship (Full Range, Log Scale)",
       x = "KING Kinship Coefficient",
       y = "Count (Log Scale)")
p_dist_full_log <- add_lines_legend(p_dist_full_log)

ggsave("king_dist_full_log.png", plot = p_dist_full_log, width = 8, height = 6, dpi = 300)
print("Plot saved to king_dist_full_log.png")

# 5. KING值分布图 (-1 到 0.5)
# 必须先过滤数据再绘图，才能重新计算 bin
df_zoom <- df[df$KINSHIP >= -1 & df$KINSHIP <= 0.5, ]
# 如果数据很少，bins=100会导致空隙，但这里点应该很多
p_dist_zoom <- ggplot(df_zoom, aes(x = KINSHIP)) +
  geom_histogram(bins = 100, fill = "forestgreen", color = "white", linewidth=0.1) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  # coord_cartesian(xlim = c(-1, 0.5)) +  # 不需要 coord_cartesian 因为我们已经subset了数据
  labs(title = "Distribution of KING Kinship (-1 to 0.5)",
       x = "KING Kinship Coefficient",
       y = "Count")
p_dist_zoom <- add_lines_legend(p_dist_zoom)

ggsave("king_dist_zoom.png", plot = p_dist_zoom, width = 8, height = 6, dpi = 300)
print("Plot saved to king_dist_zoom.png")

# 6. KING值分布图 (-12 到 0.5)
df_neg12 <- df[df$KINSHIP >= -12 & df$KINSHIP <= 0.5, ]
# 重新计算 bins
p_dist_neg12 <- ggplot(df_neg12, aes(x = KINSHIP)) +
  geom_histogram(bins = 200, fill = "forestgreen", color = "white", linewidth=0.1) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Distribution of KING Kinship (-12 to 0.5)",
       x = "KING Kinship Coefficient",
       y = "Count")
p_dist_neg12 <- add_lines_legend(p_dist_neg12)

ggsave("king_dist_neg12.png", plot = p_dist_neg12, width = 8, height = 6, dpi = 300)
print("Plot saved to king_dist_neg12.png")