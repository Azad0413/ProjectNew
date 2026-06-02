rm(list = ls())
setwd('/home/data/project/54.LZPJ-3001-6/')
file.dir = '14_Immunity'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)

#01 4 种免疫细胞标记基因-------
# （注：四种免疫细胞的标记基因参考文献 PMID：30761122，Th1：T-bet (TBX21)，
# STAT4，STAT1，IFN-γ (IFNG)和 TNF-α (TNF)；Th2：GATA3，STAT6，STAT5A
# 和 IL13；Th17：STAT3 和 IL17A；Treg：FOXP3，CCR8，STAT5B，TGFβ (TGFB1)）。
dat <- read.csv('../00_rawdata/dat.GSE65204.csv',row.names = 1,check.names = F)
group <- read.csv('../00_rawdata/group.GSE65204.csv') 
group$group <- gsub('Asthma','Childhood Asthma',group$group)
dat <- dat[,group$sample]
table(group$group)
control.sample <- group$sample[which(group$group== 'Control')]
allgene <- as.data.frame(rownames(dat))
geneset <- data.frame(symbol = c('TBX21','STAT4','STAT1','IFG','TNF',
                                 'GATA3','STAT6','STAT5A','IL13','STAT3','IL17A',
                                 'FOXP3','CCR8','STAT5B','TGFB1'),
                      Type = c(rep('Th1',5),rep('Th2',4),rep('Th17',2),rep('Treg',4)))
immu.dat <- dat[geneset$symbol,] %>% na.omit()
immu.dat <- immu.dat[-8,]
##加小提琴图
violin.dat <- immu.dat %>% as.data.frame()
violin.dat$Symbol<-rownames(violin.dat)
violin.dat<-gather(violin.dat,
                   key = sample,
                   value = score,
                   -c("Symbol"))

## 样本分组
violin.dat$Group <- group$group[match(violin.dat$sample,group$sample)]
violin.dat$Group <- factor(violin.dat$Group,levels = c('Control','Childhood Asthma'))
table(violin.dat$Group)
# 1. 合并 Type 信息
violin.dat <- left_join(violin.dat, geneset, by = c("Symbol" = "symbol"))

# 2. 分别进行检验
library(rstatix)
table(violin.dat$Type)
stat <- violin.dat %>%
  # filter(Type == "MHC-I") %>%
  group_by(Symbol) %>%
  wilcox_test(score ~ Group) %>%
  adjust_pvalue(method = "fdr") #%>%
# mutate(Type = "MHC-I")

write.csv(stat,file = 'res.wilcox.immunity.csv',row.names = F,quote = F)
# 4. 提取显著差异的基因
de.res <- stat %>% filter(p < 0.05)
length(unique(de.res$Symbol))  # 显著的 Symbol 数目
# 5. 只保留显著的基因用于可视化
violin.plotdata <- violin.dat %>% filter(Symbol %in% de.res$Symbol)

# 6. 绘图：使用 facet_wrap 分面展示
p2 <- ggplot(violin.plotdata, aes(x = Symbol, y = score, fill = Group)) +
  stat_boxplot(geom = "errorbar", width = 0.15, position = position_dodge(0.75)) +
  geom_boxplot(width = 0.55, position = position_dodge(0.75), outlier.shape = NA,
               color = "black", linewidth = 0.4) +
  geom_jitter(shape = 21, size = 1.2, stroke = 0.2, color = "black",
              position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.75),
              alpha = 0.8) +
  stat_compare_means(aes(group = Group), label = "p.signif", hide.ns = TRUE, method = "wilcox.test") +
  scale_fill_manual(values = c("#A6D854", "#FC8D62")) +
  labs(x = "", y = "Expression") +
  facet_wrap(~Type, scales = "free",nrow = 2) +  # 按类型分面展示
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", size = 0.5)
  )

p2

# 7. 保存图像
ggsave("01.Immunity.pdf", p2, width = 6, height = 5)
ggsave("01.Immunity.png", p2, width = 6, height = 5, units = 'in', dpi = 300, bg = 'white')

##相关性---------
library(ggcorrplot)
library(corrplot)
library(Hmisc)
hubgene <- read.csv('../05_ML/hubgene.csv')
hub.exp <- dat[hubgene$symbol,]
detiics <- immu.dat#[de.res$Symbol,]
nc<-t(rbind(hub.exp,detiics))
m=rcorr(nc,type = 'spearman')$r[1:nrow(hub.exp),(ncol(nc)-length(rownames(detiics))+1):ncol(nc)]
m<-t(m)

p=rcorr(nc,type = 'spearman')$P[1:nrow(hub.exp),(ncol(nc)-length(rownames(detiics))+1):ncol(nc)]
p<-t(p)
library(reshape2)
correlation <- reshape2::melt(m)
colnames(correlation) <- c('ImmueGenes','Gene','correlation')
pvalue <- reshape2::melt(p)
correlation$pvalue <- pvalue$value
write.csv(correlation,file = 'correlation_Immunity.csv',row.names = F)

library(reshape2)
library(ggplot2)

# 转换为 long format
corr_long <- reshape2::melt(t(m))
pval_long <- reshape2::melt(t(p))
colnames(corr_long) <- c("Var1", "Var2", "cor")
colnames(pval_long) <- c("Var1", "Var2", "pval")

# 合并数据
df <- merge(corr_long, pval_long, by = c("Var1", "Var2"))

# 添加显著性标记
df$sig <- cut(df$pval,
              breaks = c(-Inf,0.0001, 0.001, 0.01, 0.05, Inf),
              labels = c("****", "***", "**", "*", ""))

# 用 -log10(pval) 控制填色（显著性越高越深）
# df$logp <- -log10(df$pval)
df$logp <- df$cor

pdf(file = '02.correlation_Immunity.pdf',w=6,h=5)
ggplot(df, aes(x = Var1, y = Var2, fill = logp)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sig), size = 5, color = "black") +
  scale_fill_gradientn(colors = c("#A6D854",'white',"#FC8D62"),
                       name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        panel.grid = element_blank()) +
  labs(title = "Spearman Correlation - HubGenes and Marker Genes",
       x = NULL, y = NULL)#+coord_flip()

dev.off()


png(file = '02.correlation_Immunity.png',w=6,h=5,units = 'in',res = 300,bg = 'white')
ggplot(df, aes(x = Var1, y = Var2, fill = logp)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sig), size = 5, color = "black") +
  scale_fill_gradientn(colors = c("#A6D854",'white',"#FC8D62"),
                       name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        panel.grid = element_blank()) +
  labs(title = "Spearman Correlation - HubGenes and Marker Genes",
       x = NULL, y = NULL)#+coord_flip()

dev.off()
