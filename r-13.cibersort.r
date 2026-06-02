rm(list = ls())
setwd('/home/data/project/54.LZPJ-3001-6/')
file.dir = '13_cibersort'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)

dat <- read.csv('../00_rawdata/dat.GSE65204.csv',row.names = 1,check.names = F)
group <- read.csv('../00_rawdata/group.GSE65204.csv')
str(group)
group$group <- gsub('Asthma','Childhood Asthma',group$group)
dat <- dat[,group$sample]
library(tidyverse)
library(immunedeconv)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
set_cibersort_binary("../../../files/TIIC/CIBERSORT.R")
set_cibersort_mat("../../../files/TIIC/LM22.txt")
res.cibersort <- deconvolute(as.matrix(dat), method = 'cibersort')
write.table(res.cibersort, 'cibersort.txt', sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
save(file = 'cibersort.Rdata', res.cibersort)

display.brewer.all()
mypalette <- colorRampPalette(brewer.pal(9, "Set3"))
plot_cibersort <- function(data, group) {
  data %>%
    gather(sample, fraction, -cell_type) %>%
    merge(group, by = 'sample') %>%
    ggplot(aes(x = sample, y = fraction, fill = cell_type)) +
    geom_bar(position = 'stack', stat = 'identity') +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw() +
    labs(x = '', y = 'Relative Percent', fill = '') +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 'top') +
    scale_fill_manual(values = mypalette(22)) +
    facet_grid(~group, scales = "free", space = "free")
}

pdf('01.cibersort.box.pdf', width = 10, height = 8)
print(plot_cibersort(res.cibersort, group))
dev.off()

png('01.cibersort.box.png', width = 10, height = 8, units = 'in', res = 300)
print(plot_cibersort(res.cibersort, group))
dev.off()

res.cibersort2 <- res.cibersort %>% column_to_rownames(var = 'cell_type')
keep <- rowSums(res.cibersort2 > 0) >= floor(0.25 * ncol(res.cibersort2))
res.cibersort2 <- res.cibersort2[keep, ]
dat.cibersort <- res.cibersort2 %>%
  t %>% as.data.frame() %>%
  rownames_to_column(var = "sample")
dat.cibersort <- merge(group, dat.cibersort, by = "sample")
dat.cibersort2 <- gather(dat.cibersort, ImmuneCell, Score, -c("sample", "group"))

stat_cibersort <- dat.cibersort2 %>%
  group_by(ImmuneCell) %>%
  wilcox_test(Score ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p")

write.csv(stat_cibersort, 'stat.cibersort.csv', row.names = TRUE)
DE.cibersort <- stat_cibersort %>% filter(p < 0.05)
write.csv(DE.cibersort, 'DE.cibersort.csv', row.names = TRUE)

violin.cibersort <- dat.cibersort2 %>% filter(ImmuneCell %in% DE.cibersort$ImmuneCell)
violin.cibersort$group <- factor(violin.cibersort$group, levels = c('Control', 'Childhood Asthma'))
num_cells <- length(unique(violin.cibersort$ImmuneCell))
nrow_val <- ceiling(num_cells / 6)
cibersort_plot <- ggplot(violin.cibersort, aes(x = group, y = Score, fill = group)) +
  geom_violin(trim = FALSE, color = "black") +
  stat_boxplot(geom = "errorbar", width = 0.1, position = position_dodge(0.9)) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.shape = NA, fill = 'white') +
  scale_fill_manual(values = c("#A6D854", "#FC8D62"), name = "Group") +
  labs(title = "Immune Cell", x = "", y = "Score") +
  stat_compare_means(label = "p.signif", method = 'wilcox.test', paired = FALSE, label.x = 1.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        axis.text.x = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_wrap(~ImmuneCell, scales = "free", nrow = nrow_val)
cibersort_plot
ggsave('02.cibersort.plot.pdf', plot = cibersort_plot, width = 12, height = nrow_val*2+3)
ggsave('02.cibersort.plot.png', plot = cibersort_plot, width = 12, height = nrow_val*2+3, dpi = 300,bg = 'white')

#基因相关性-----
library(vegan)
score1<- t(res.cibersort2[which(rownames(res.cibersort2) %in% DE.cibersort$ImmuneCell),])%>%as.data.frame()
hubgene <- read.csv('../05_ML/hubgene.csv')
dat1 <- dat
hub.dat1 <- dat1[hubgene$symbol,]
hub.dat1 <- t(hub.dat1)%>%as.data.frame()
library(Hmisc)
library(psych)
library(tidyverse)
library(ggcor)

# 提取所有样本交集
common_samples <- intersect(rownames(hub.dat1), rownames(score1))

# 保证样本顺序一致
hub.sub <- hub.dat1[common_samples, ]
score.sub <- score1[common_samples, ]

# 计算 Spearman 相关性：gene × cell
corr_result <- corr.test(hub.sub, score.sub, method = "spearman")

# 整理成长格式用于 anno_link
corr_long <- as.data.frame(corr_result$r) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cell", values_to = "r") %>%
  left_join(
    as.data.frame(corr_result$p) %>%
      rownames_to_column("gene") %>%
      pivot_longer(-gene, names_to = "cell", values_to = "p.value"),
    by = c("gene", "cell")
  )

# 提取想强调的基因
genes_of_interest <- hubgene$symbol

mantel_like <- corr_long %>%
  filter(gene %in% genes_of_interest) %>%
  mutate(
    rd = case_when(
      abs(r) < 0.3 ~ "< 0.3 (weak)",
      abs(r) < 0.6 ~ "< 0.6 (moderate)",
      TRUE ~ "> 0.6 (strong)"
    ),
    
    pd = case_when(
      p.value < 0.001 ~ "< 0.001",
      p.value < 0.01 ~ "< 0.01",
      p.value < 0.05 ~ "< 0.05",
      TRUE ~ "ns"
    )
  )
# 绘图：热图为 cell × cell，连线为 gene × cell 的显著 Spearman
pdf('03.heatmap_gene.pdf', w=10, h=7)
quickcor(score1, type = "lower") + 
  geom_circle2() +
  scale_fill_gradient2(low = "#9ECAE1", mid = "white", high = "#984EA3", midpoint = 0) +
  anno_link(aes(colour = pd, size = rd), data = mantel_like) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c('#56B4E9', "#E69F00", "#CC79A7", "#A2A2A288")) +
  guides(size = guide_legend(title = "Spearman |r|", order = 2), 
         colour = guide_legend(title = "Spearman p", order = 3), 
         fill = guide_colorbar(title = "Cell-Cell Spearman", order = 4))
dev.off()

png('03.heatmap_gene.png',w=10,h=7,units = 'in',res = 300,bg = 'white')
quickcor(score1, type = "lower") + 
  geom_circle2() +
  scale_fill_gradient2(low = "#9ECAE1", mid = "white", high = "#984EA3", midpoint = 0) +
  anno_link(aes(colour = pd, size = rd), data = mantel_like) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c('#56B4E9', "#E69F00", "#CC79A7", "#A2A2A288")) +
  guides(size = guide_legend(title = "Spearman |r|", order = 2), 
         colour = guide_legend(title = "Spearman p", order = 3), 
         fill = guide_colorbar(title = "Cell-Cell Spearman", order = 4))
dev.off()




##topcor-------
library(ggcorrplot)
library(corrplot)
library(Hmisc)
library(ggExtra)

# NTRK2 B cell plasma
marker <- as.numeric(hub.dat1[,'NTRK2'])
risk <- as.numeric(score.sub[,'B cell plasma'])
df1 <- cbind(risk,marker)%>%as.data.frame()
p1 <- ggplot(df1,aes(marker,risk))+
  xlab('NTRK2')+ylab('B cell plasma')+
  geom_point()+geom_smooth(method = 'lm',formula = y~x)+
  stat_cor(data=df1, method = "spearman")+
  theme_bw()
p1  
p2 <- ggMarginal(p1,type = 'density',xparams = list(fill='#FC8D62'),
                 yparams = list(fill='#FC8D62'))

p2
ggsave(filename = '04.topgenecor_pos.pdf', p2, w = 5, h = 5)
ggsave(filename = '04.topgenecor_pos.png', p2, w = 5, h = 5,units = 'in',dpi = 600)

# CEACAM5 
marker <- as.numeric(hub.dat1[,'CEACAM5'])
risk <- as.numeric(score.sub[,'T cell follicular helper'])
df1 <- cbind(risk,marker)%>%as.data.frame()
p1 <- ggplot(df1,aes(marker,risk))+
  xlab('CEACAM5')+ylab('T cell follicular helper')+
  geom_point()+geom_smooth(method = 'lm',formula = y~x)+
  stat_cor(data=df1, method = "spearman")+
  theme_bw()
p1  
p2 <- ggMarginal(p1,type = 'density',xparams = list(fill='#A6D854'),
                 yparams = list(fill='#A6D854'))

p2
ggsave(filename = '05.topgenecor_neg.pdf', p2, w = 5, h = 5)
ggsave(filename = '05.topgenecor_neg.png', p2, w = 5, h = 5,units = 'in',dpi = 600)



