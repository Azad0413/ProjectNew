rm(list = ls())
setwd('/home/data/project/54.LZPJ-3001-6/')
file.dir = '15_network'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)

##mirnet-----
gene2mir <- read.csv('gene2mir.csv')
length(unique(gene2mir$ID))
length(unique(gene2mir$Target))

mirwalk <- read.csv('miRWalk_miRNA_Targets.csv')

gene2mir$link <- paste0(gene2mir$Target,gene2mir$ID)
mirwalk$link <- paste0(mirwalk$genesymbol,mirwalk$mirnaid)

##mirdb
mirdb <- readxl::read_xlsx('mirdb.xlsx')
mirdb$link <- paste0(mirdb$`Gene Symbol`,mirdb$`miRNA Name`)

##取交集
mircor <- gene2mir[which(gene2mir$link %in% mirwalk$link),]
mircor <- mircor[which(mircor$link %in% mirdb$link),]
length(unique(mircor$Target))
##2
length(unique(mircor$ID))
#30

mir <- data.frame(mir = unique(mircor$ID))

write.csv(mir,file = 'intersect.mir.csv',row.names = F,quote = F)

mir2lnc <- read.csv('mir2lnc.csv')
length(unique(mir2lnc$Target))
##415
length(unique(mir2lnc$ID))
##23

encori.all <- read.delim2('../../../files/ENCORI/lncRNA/miRNAcat_hsa.txt',skip = 3)
encori.lnc <- encori.all[which(encori.all$miRNAname %in% unique(mir$mir)),]
table(encori.lnc$geneType)
encori.lnc <- encori.lnc[which(encori.lnc$geneType == 'lncRNA'),]
length(unique(encori.lnc$miRNAid))
##16
length(unique(encori.lnc$geneName))
##18
write.csv(encori.lnc,file = 'encori.lnc.csv',row.names = F)
encori.lnc <- encori.lnc[,c(2,4)]
colnames(mir2lnc)
colnames(encori.lnc) <- c('ID','lncRNA')
mir2lnc <- mir2lnc[,c(1,3)]
encori.lnc$link <- paste0(encori.lnc$ID,encori.lnc$lncRNA)
mir2lnc$link <- paste0(mir2lnc$ID,mir2lnc$Target)
library(tidyverse)
##取交集
mir2lnc <- mir2lnc[which(mir2lnc$link %in% encori.lnc$link),] %>% dplyr::select(-'link')
length(unique(mir2lnc$ID))
##13
length(unique(mir2lnc$Target))
##9
mircor <- mircor[,c(1,3)]
mircor <- mircor[which(mircor$ID %in% mir2lnc$ID),]
length(unique(mircor$ID))
length(unique(mircor$Target))
##2
colnames(mircor) <- c('miRNA','gene')
colnames(mir2lnc) <- c('miRNA','lncRNA')
cornet <- merge(mircor,mir2lnc,by = 'miRNA')
write.csv(cornet,file = 'mRNA.miRNA.lncRNA.csv',row.names = F,quote = F)
length(unique(cornet$gene))
length(unique(cornet$miRNA))
##13
length(unique(cornet$lncRNA))
##9
unique(cornet$gene)

library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)

str(cornet)
# 构建边集：miRNA -> gene，lncRNA -> gene
edges1 <- cornet[, c("miRNA", "gene")]
colnames(edges1) <- c("from", "to")

edges2 <- cornet[, c("lncRNA", "miRNA")]
colnames(edges2) <- c("from", "to")

edges <- distinct(rbind(edges1, edges2))  # 合并后去重

# 构建图
graph <- graph_from_data_frame(edges, directed = TRUE)

# 提取节点名称
nodes <- data.frame(name = V(graph)$name, stringsAsFactors = FALSE)

# 标注节点类型
nodes <- nodes %>%
  mutate(type = case_when(
    name %in% cornet$miRNA ~ "miRNA",
    name %in% cornet$lncRNA ~ "lncRNA",
    name %in% cornet$gene ~ "gene",
    TRUE ~ "other"
  ))

# 加入度数信息
nodes$deg <- degree(graph, mode = "all")

# 添加节点属性
V(graph)$type <- nodes$type
V(graph)$deg <- nodes$deg

set.seed(123)

p <- ggraph(graph, layout = 'fr') +
  geom_edge_link(arrow = arrow(length = unit(2, 'mm')),
                 end_cap = circle(3, 'mm'),
                 edge_alpha = 0.5) +
  geom_node_point(aes(color = type, size = deg)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3, max.overlaps = Inf) +
  scale_color_manual(
    values = c(
      "miRNA" = "#E07B91",  # 红
      "gene" = "#FC8D62",   # 蓝
      "lncRNA" = "#A6D854"  # 绿
    ),
    name = "Node Type"
  ) +
  scale_size_continuous(range = c(3, 6), name = "Degree") +
  theme_void() +
  labs(title = "mRNA-miRNA-lncRNA Regulatory Network")

p

# 保存图像
ggsave(filename = '01.mRNA-miRNA–lncRNA.pdf', p, width = 7, height = 6)
ggsave(filename = '01.mRNA-miRNA–lncRNA.png', p, width = 7, height = 6, dpi = 300, bg = 'white')



##tf-----
gene2tf <- read.csv('gene2tf_jaspar.csv')
length(unique(gene2tf$ID))
# === 1. 数据准备 ===
edges <- gene2tf[, c("ID", "Target")]
colnames(edges) <- c("from", "to")

# 构建网络
graph <- graph_from_data_frame(edges, directed = TRUE)

# 提取所有节点名称
nodes <- data.frame(name = V(graph)$name, stringsAsFactors = FALSE)

# === 2. 标注节点类型 ===
nodes$type <- ifelse(nodes$name %in% edges$from, "TF", "mRNA")

# === 3. 添加度数信息 ===
nodes$deg <- degree(graph, mode = "all")

# 将类型与度数写入图的属性中
V(graph)$type <- nodes$type
V(graph)$deg <- nodes$deg

# === 4. 绘图 ===
set.seed(888)
p <- ggraph(graph, layout = 'fr') +
  geom_edge_link(arrow = arrow(length = unit(4, 'mm')),
                 end_cap = circle(3, 'mm'),
                 edge_alpha = 0.5) +
  geom_node_point(aes(color = type, size = deg)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3,max.overlaps = Inf) +
  scale_color_manual(
    values = c("TF" = "#A6D854", "mRNA" = "#FC8D62"),
    name = "Node Type"
  ) +
  scale_size_continuous(name = "Degree", range = c(3, 6)) +
  theme_void() +
  labs(title = "TF-mRNA Regulation Network")

p
ggsave(filename = '02.mRNA-TF.pdf', p, width = 6, height = 5)
ggsave(filename = '02.mRNA-TF.png', p, width = 6, height = 5, units = 'in', dpi = 300, bg = 'white')

