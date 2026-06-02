rm(list = ls())
setwd('/home/data/project/54.LZPJ-3001-6/')
file.dir = '03_PPI'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
library(igraph)
library(ggraph)
library(tidyverse)
library(tidygraph)

ppi <- read_tsv('string_interactions_short.tsv')
g <- graph_from_data_frame(d=ppi[, c("#node1", "node2")], directed=FALSE)

g <- graph_from_data_frame(d=ppi[, c("#node1", "node2", "combined_score")], directed=FALSE)
tg <- as_tbl_graph(g) %>%
  mutate(deg = centrality_degree())

set.seed(123)
p <- ggraph(tg, layout = 'fr') + 
  geom_edge_link(aes(width = combined_score),
                 edge_alpha = 0.5,
                 lineend = "round",
                 linejoin = "round") +
  geom_node_point(aes(color = deg), size = 4) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_edge_width(range = c(0.2, 2), name = "Combined Score") +
  scale_color_viridis_c(name = "Degree", option = "B", direction = -1) +
  theme_void() +
  labs(title = "PPI Network (Weighted by Combined Score)")

print(p)
ggsave(filename = '01.PPI.pdf',p,width = 6,height = 5)
ggsave(filename = '01.PPI.png',p,width = 6,height = 5,units = 'in',dpi = 300,bg = 'white')
