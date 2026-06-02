rm(list = ls())
setwd('/home/data/project/54.LZPJ-3001-6/')
file.dir = '02_DEPCDs'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
deg <- read.csv('../01_DEGs/DEG_sig.csv',row.names = 1)
geneset <- read.csv('../00_rawdata/19pcd.geneset.csv')
candi <- data.frame(symbol = intersect(rownames(deg),geneset$Gene))
write.csv(candi,file = 'depcds.csv',row.names = F,quote = F)

library(ggvenn)
venn_list <- list('DEGs' = rownames(deg), 'PCDGs' = geneset$Gene)
text <- paste0("**DE-PCDGs:**<br>", paste(candi$symbol, collapse = "<br>"))

library(ggvenn)
library(ggtext)
library(patchwork)
# 绘制 Venn 图
venn <- ggvenn(venn_list, 
               fill_color = c("#33a02c","#EFC000FF"), 
               stroke_size = 0.5, set_name_size = 5) +
  theme(plot.margin = margin(10, 50, 10, 10))

venn
# 文本框
box <- ggplot() +
  geom_textbox(
    aes(x = 0.5, y = 0.5, label = text),  # 中心位置
    box.color = "#EFC000FF",
    fill = "white",
    width = unit(0.7, "npc"),   # 宽度占满整个画布
    height = unit(0.5, "npc"),  # 高度占满整个画布（新加）
    hjust = 0.5, vjust = 0.5, # 居中定位文本框
    halign = 0,               # 文本框内靠左对齐文字
    size = 4,
    lineheight = 1.3
  ) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))
box
# patchwork组合
final_plot <- (venn | box) +
  plot_layout(widths = c(2, 1), heights = c(1, 1)) &
  theme(plot.margin = margin(5, 5, 5, 5))
final_plot
# 保存图像，防止裁切
ggsave("01.venn_plot.pdf", plot = final_plot,
       width = 8, height = 6)
ggsave("01.venn_plot.png", plot = final_plot,
       width = 8, height = 6, units = "in", dpi = 300)




