library(igraph)

# 1. 准备数据
# 注意：KING 的输出文件通常以 #IID1 开头，需要设置 comment.char="" 避免被当做注释忽略
# 且列名实际上是 IID1 和 IID2
king_data <- read.table("/data1/dazheng_tusr1/vmap4.VCF.v1/chr002_king_filter_0.2.kin0", header = TRUE, comment.char = "")

# 修正第一列的列名（去除可能的 # 前缀或 R 自动添加的 X.）
colnames(king_data)[1] <- "IID1"

# 2. 设定阈值（例如 0.0442 以上视为有关系）
# 使用 subset 函数来过滤，它可以自动排除 KINSHIP 为 NA 的行
edges <- subset(king_data, KINSHIP > 0.0442, select = c("IID1", "IID2"))

# 再次确保边数据中没有 NA（因为 graph_from_data_frame 不允许）
edges <- na.omit(edges)

# 3. 创建图形对象
g <- graph_from_data_frame(edges, directed = FALSE)

# 打印基本统计信息，帮助理解图像复杂度
cat("----------------------------------\n")
cat("图的统计信息:\n")
cat("节点数 (Samples):", vcount(g), "\n")
cat("边数 (Relationships):", ecount(g), "\n")
cat("----------------------------------\n")

# 4. 绘图优化
png("kinship_network.png", width = 2400, height = 2400, res = 300)

# 使用 Fruchterman-Reingold 布局算法，这种算法会把有联系的节点拉近，没有联系的推远
# 从而由于亲缘关系形成自然的“簇” (Clusters)
l <- layout_with_fr(g) 

plot(g, 
     layout = l,                  # 使用力导向布局
     vertex.size = 1.2,           # 节点缩小
     vertex.color = rgb(0.8, 0.2, 0.2, 0.5), # 节点红色半透明
     vertex.frame.color = NA,     # 去掉节点外面的圆圈 (border)
     vertex.label = NA,           # 不显示标签
     edge.width = 0.3,            # 边设为细线 (因为全黑如果不细会变成一团黑)
     edge.arrow.size = 0,         # 去除箭头（无向图通常没有，但保险起见）
     edge.color = "black",        # 边使用纯黑色，不透明，确保最突出
     main = paste0("Kinship Network (Kinship > 0.0442)\nNodes: ", vcount(g), "  Edges: ", ecount(g)))

dev.off()
print("Plot saved to kinship_network.png")