rm(list = ls())
setwd('/home/data/project/54.LZPJ-3001-6/')
file.dir = '08_expression'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
library(tidyverse)
dat <- read.csv('../00_rawdata/dat.GSE65204.csv',row.names = 1,check.names = F)
group <- read.csv('../00_rawdata/group.GSE65204.csv')
table(group$group)
hubgene <- read.csv('../05_ML/hubgene.csv')
control.sample <- group$sample[which(group$group=='Control')]
hub_exp<-dat[hubgene$symbol,]
hub_exp2 <- as.data.frame(hub_exp) %>% na.omit()
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'Control','Childhood Asthma')

hub_exp2$Group <- factor(hub_exp2$Group,levels = c('Control','Childhood Asthma'))
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
stat.test<-hub_exp2%>%
  group_by(Symbol)%>%
  wilcox_test(expr ~ Group)%>%
  adjust_pvalue(method = 'fdr')
stat.test$p<-ifelse(stat.test$p<0.0001,"****",
                    ifelse(stat.test$p<0.001,"***",
                           ifelse(stat.test$p<0.01,"**",
                                  ifelse(stat.test$p<0.05,"*",'ns'))))
colnames(hub_exp2)
exp_plot <- ggplot(hub_exp2, aes(x = Group, 
                                 y = expr, 
                                 fill = Group)) +
  # 精致误差线
  stat_boxplot(geom = "errorbar",
               width = 0.15,
               position = position_dodge(0.75)) +
  
  # 箱线图：无离群点，适度宽度
  geom_boxplot(width = 0.5,
               position = position_dodge(0.75),
               outlier.shape = NA,
               color = "black", 
               lwd = 0.3) +
  
  # 添加散点（可选 jitter 更自然）
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
              size = 0.8,
              shape = 21,
              stroke = 0.2,
              color = "black") +
  
  scale_fill_manual(values = c("#33a02c", "#EFC000FF")) +
  
  # 添加显著性标记
  stat_compare_means(aes(group = Group), 
                     label = "p.signif", method = 'wilcox.test',
                     hide.ns = F, 
                     size = 4, 
                     label.y.npc = "top",label.x = 1.4) +
  
  # 标签设置
  labs(title = "Expression Level in Trian Set", x = "", y = "Expression", size = 20) + 
  
  # 美化主题
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5,size = 16),
    axis.text.x = element_text(angle = 75, hjust = 1,  size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey60", size = 0.4)
  )+facet_wrap(~Symbol,scales = 'free',nrow = 1)
exp_plot
ggsave('01.expression_GSE65204.pdf',exp_plot,width = 7,height = 5)
ggsave('01.expression_GSE65204.png',exp_plot,width = 7,height = 5,units = 'in',dpi = 300)


#验证集----------
library(tidyverse)
dat <- read.csv('..//00_rawdata/dat.GSE19187.csv',row.names = 1,check.names = F)
allgene <- data.frame(row.names(dat))
group <- read.csv('../00_rawdata/group.GSE19187.csv')
table(group$group)
control.sample <- group$sample[which(group$group=='Control')]
hub_exp<-dat[hubgene$symbol,]
hub_exp2 <- as.data.frame(hub_exp) %>% na.omit()
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

## 样本分组
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'Control','Childhood Asthma')

hub_exp2$Group <- factor(hub_exp2$Group,levels = c('Control','Childhood Asthma'))
##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
stat.test<-hub_exp2%>%
  group_by(Symbol)%>%
  wilcox_test(expr ~ Group)%>%
  adjust_pvalue(method = 'fdr')
stat.test$p<-ifelse(stat.test$p<0.0001,"****",
                    ifelse(stat.test$p<0.001,"***",
                           ifelse(stat.test$p<0.01,"**",
                                  ifelse(stat.test$p<0.05,"*",'ns'))))
colnames(hub_exp2)
exp_plot <- ggplot(hub_exp2, aes(x = Group, 
                                 y = expr, 
                                 fill = Group)) +
  # 精致误差线
  stat_boxplot(geom = "errorbar",
               width = 0.15,
               position = position_dodge(0.75)) +
  
  # 箱线图：无离群点，适度宽度
  geom_boxplot(width = 0.5,
               position = position_dodge(0.75),
               outlier.shape = NA,
               color = "black", 
               lwd = 0.3) +
  
  # 添加散点（可选 jitter 更自然）
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
              size = 0.8,
              shape = 21,
              stroke = 0.2,
              color = "black") +
  
  scale_fill_manual(values = c("#33a02c", "#EFC000FF")) + 
  
  # 添加显著性标记
  stat_compare_means(aes(group = Group), 
                     label = "p.signif", method = 'anova',
                     hide.ns = F, 
                     size = 4, 
                     label.y.npc = "top",label.x = 1.4) +
  
  # 标签设置
  labs(title = "Expression Level in Test Set", x = "", y = "Expression", size = 20) + 
  
  # 美化主题
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5,size = 16),
    axis.text.x = element_text(angle = 75, hjust = 1,  size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey60", size = 0.4)
  )+facet_wrap(~Symbol,scales = 'free',nrow = 1)
exp_plot
ggsave('02.expression_GSE19187.pdf',exp_plot,width = 7,height = 5)
ggsave('02.expression_GSE19187.png',exp_plot,width = 7,height = 5,units = 'in',dpi = 300)
