# 导入必要的包
library(RIdeogram)
library(ggplot2)

# 创建示例数据，可以根据实际情况替换为你的数据
set.seed(123)
chromosome <- "chr1"
intervals <- 1:100
transposon_counts <- rpois(100, lambda = 5)  # 示例数据

# 将数据放入数据框
data <- data.frame(chromosome, intervals, transposon_counts)
# 创建染色体模式图
data(human_karyotype, package="RIdeogram")
data(gene_density, package="RIdeogram")
data(Random_RNAs_500, package="RIdeogram")
ideogram(karyotype = human_karyotype, overlaid = gene_density, label = Random_RNAs_500, label_type = "marker", colorset1 = c("#ffffff", "#fce6c9", "#ff8000"))

# 创建每个区间内转座子数量的分布图
transposon_plot <- ggplot(data, aes(x = intervals, y = transposon_counts)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Transposon Count Distribution",
       x = "Interval",
       y = "Transposon Count") +
  theme_minimal()
