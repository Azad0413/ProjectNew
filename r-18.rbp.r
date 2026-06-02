rm(list = ls())
setwd('/home/data/project/54.LZPJ-3001-6/')
file.dir = '18_RBP'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)

gene2rbp <- readxl::read_xlsx('rbp.xlsx')

# === 1. 数据准备 ===
edges <- gene2rbp[, c("RBP", "geneName")]
colnames(edges) <- c("from", "to")

# 构建网络
graph <- graph_from_data_frame(edges, directed = TRUE)

# 提取所有节点名称
nodes <- data.frame(name = V(graph)$name, stringsAsFactors = FALSE)

# === 2. 标注节点类型 ===
nodes$type <- ifelse(nodes$name %in% edges$from, "RBP", "mRNA")

# === 3. 添加度数信息 ===
nodes$deg <- degree(graph, mode = "all")

# 将类型与度数写入图的属性中
V(graph)$type <- nodes$type
V(graph)$deg <- nodes$deg

# === 4. 绘图 ===
set.seed(123)
p <- ggraph(graph, layout = 'kk') +
  geom_edge_link(arrow = arrow(length = unit(4, 'mm')),
                 end_cap = circle(3, 'mm'),
                 edge_alpha = 0.5) +
  geom_node_point(aes(color = type, size = deg)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3,max.overlaps = Inf) +
  scale_color_manual(
    values = c("RBP" = "#E64B35", "mRNA" = "#4DBBD5"),
    name = "Node Type"
  ) +
  scale_size_continuous(name = "Degree", range = c(3, 6)) +
  theme_void() +
  labs(title = "RBP-mRNA Regulation Network")

p
ggsave(filename = '01.mRNA-RBP.pdf', p, width = 6, height = 5)
ggsave(filename = '01.mRNA-RBP.png', p, width = 6, height = 5, units = 'in', dpi = 300, bg = 'white')
