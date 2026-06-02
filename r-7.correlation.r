rm(list = ls())
setwd('/home/data/project/54.LZPJ-3001-6/')
file.dir = '07_correlation'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
library(ggcorrplot)
library(corrplot)
library(Hmisc)
hubgene <- read.csv('../05_ML/hubgene.csv')
dat <- read.csv('../00_rawdata/dat.GSE65204.csv',row.names = 1)

hub.exp <- dat[hubgene$symbol,]
nc<-t(hub.exp)
m=rcorr(nc,type = 'spearman')$r
p=rcorr(nc,type = 'spearman')$P
p<-t(p)
library(reshape2)
correlation <- reshape2::melt(m)
colnames(correlation) <- c('GeneA','GeneB','correlation')
pvalue <- reshape2::melt(p)
correlation$pvalue <- pvalue$value
write.csv(correlation,file = 'train_correlation.csv',row.names = F)

# 准备上下三角矩阵
corr_mat <- m
pval_mat <- p
colnames(corr_mat) <- rownames(corr_mat) <- colnames(nc)
colnames(pval_mat) <- rownames(pval_mat) <- colnames(nc)

# 创建上下三角数据框
lower_mat <- corr_mat
upper_mat <- matrix(NA, nrow = nrow(m), ncol = ncol(m))
upper_mat[upper.tri(upper_mat)] <- 1  # 用于标记上三角

# 转换为 long 格式
library(reshape2)
lower_df <- melt(lower_mat, na.rm = TRUE)
colnames(lower_df) <- c("Var1", "Var2", "cor")

upper_df <- melt(pval_mat, na.rm = TRUE)
colnames(upper_df) <- c("Var1", "Var2", "pval")
upper_df$sig <- cut(upper_df$pval,
                    breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                    labels = c("****", "***", "**", "*", ""))

# 过滤上三角位置（避免下三角也画上去）
upper_df <- upper_df[as.numeric(upper_df$Var1) < as.numeric(upper_df$Var2), ]
lower_df <- lower_df[as.numeric(lower_df$Var1) > as.numeric(lower_df$Var2), ]
levels(lower_df$Var1)
levels(upper_df$Var1)
lower_df$Var1 <- factor(lower_df$Var1,levels = c('CEACAM5','NTRK2','CTSC','CCK','CTSG'))
library(ggplot2)
library(ggnewscale)  # ✅ 支持多个 fill
pdf("01.Train_correlation.pdf", width = 6, height = 5)
ggplot() +
  # 下三角：相关系数填色
  geom_tile(data = lower_df, aes(x = Var1, y = Var2, fill = cor), color = "white") +
  scale_fill_gradientn(colors = c("#377eb8", "white", "#e41a1c"),
                       limits = c(-1, 1),
                       name = "Correlation") +
  ggnewscale::new_scale_fill() +  # 🟢 添加新的 fill 变量映射
  # 上三角：用 p 值着色（-log10(p)）
  geom_tile(data = upper_df, aes(x = Var1, y = Var2, fill = -log10(pval + 1e-20)), color = "white") +
  scale_fill_gradientn(colors = c("white", "#c7e9c0", "#74c476", "#238b45"),
                       name = "-log10(p)") +
  
  # 上三角：显著性标记
  geom_text(data = upper_df, aes(x = Var1, y = Var2, label = sig), size = 4) +
  
  # 下三角：相关性数值
  geom_text(data = lower_df, aes(x = Var1, y = Var2, label = round(cor, 2)), size = 4) +
  
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        panel.grid = element_blank()) +
  labs(title = "Train Cohort: Correlation Between HubGenes",
       x = NULL, y = NULL)

dev.off()

png("01.Train_correlation.png", width = 6, height = 5,units = 'in',res = 300,bg = 'white')
ggplot() +
  # 下三角：相关系数填色
  geom_tile(data = lower_df, aes(x = Var1, y = Var2, fill = cor), color = "white") +
  scale_fill_gradientn(colors = c("#377eb8", "white", "#e41a1c"),
                       limits = c(-1, 1),
                       name = "Correlation") +
  ggnewscale::new_scale_fill() +  # 🟢 添加新的 fill 变量映射
  # 上三角：用 p 值着色（-log10(p)）
  geom_tile(data = upper_df, aes(x = Var1, y = Var2, fill = -log10(pval + 1e-20)), color = "white") +
  scale_fill_gradientn(colors = c("white", "#c7e9c0", "#74c476", "#238b45"),
                       name = "-log10(p)") +
  
  # 上三角：显著性标记
  geom_text(data = upper_df, aes(x = Var1, y = Var2, label = sig), size = 4) +
  
  # 下三角：相关性数值
  geom_text(data = lower_df, aes(x = Var1, y = Var2, label = round(cor, 2)), size = 4) +
  
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        panel.grid = element_blank()) +
  labs(title = "Train Cohort: Correlation Between HubGenes",
       x = NULL, y = NULL)

dev.off()


##验证集相关性-----
library(ggcorrplot)
library(corrplot)
library(Hmisc)
hubgene <- read.csv('../05_ML/hubgene.csv')
dat <- read.csv('../00_rawdata/dat.GSE19187.csv',row.names = 1)

hub.exp <- dat[hubgene$symbol,]
nc<-t(hub.exp)
m=rcorr(nc,type = 'spearman')$r
p=rcorr(nc,type = 'spearman')$P
p<-t(p)
library(reshape2)
correlation <- reshape2::melt(m)
colnames(correlation) <- c('GeneA','GeneB','correlation')
pvalue <- reshape2::melt(p)
correlation$pvalue <- pvalue$value
write.csv(correlation,file = 'test_correlation.csv',row.names = F)

# 准备上下三角矩阵
corr_mat <- m
pval_mat <- p
colnames(corr_mat) <- rownames(corr_mat) <- colnames(nc)
colnames(pval_mat) <- rownames(pval_mat) <- colnames(nc)

# 创建上下三角数据框
lower_mat <- corr_mat
upper_mat <- matrix(NA, nrow = nrow(m), ncol = ncol(m))
upper_mat[upper.tri(upper_mat)] <- 1  # 用于标记上三角

# 转换为 long 格式
library(reshape2)
lower_df <- melt(lower_mat, na.rm = TRUE)
colnames(lower_df) <- c("Var1", "Var2", "cor")

upper_df <- melt(pval_mat, na.rm = TRUE)
colnames(upper_df) <- c("Var1", "Var2", "pval")
upper_df$sig <- cut(upper_df$pval,
                    breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                    labels = c("****", "***", "**", "*", ""))

# 过滤上三角位置（避免下三角也画上去）
upper_df <- upper_df[as.numeric(upper_df$Var1) < as.numeric(upper_df$Var2), ]
lower_df <- lower_df[as.numeric(lower_df$Var1) > as.numeric(lower_df$Var2), ]
levels(lower_df$Var1)
levels(upper_df$Var1)
lower_df$Var1 <- factor(lower_df$Var1,levels = c('CEACAM5','NTRK2','CTSC','CCK','CTSG'))
library(ggplot2)
library(ggnewscale)  # ✅ 支持多个 fill
pdf("02.Test_correlation.pdf", width = 6, height = 5)
ggplot() +
  # 下三角：相关系数填色
  geom_tile(data = lower_df, aes(x = Var1, y = Var2, fill = cor), color = "white") +
  scale_fill_gradientn(colors = c("#377eb8", "white", "#e41a1c"),
                       limits = c(-1, 1),
                       name = "Correlation") +
  ggnewscale::new_scale_fill() +  # 🟢 添加新的 fill 变量映射
  # 上三角：用 p 值着色（-log10(p)）
  geom_tile(data = upper_df, aes(x = Var1, y = Var2, fill = -log10(pval + 1e-20)), color = "white") +
  scale_fill_gradientn(colors = c("white", "#c7e9c0", "#74c476", "#238b45"),
                       name = "-log10(p)") +
  
  # 上三角：显著性标记
  geom_text(data = upper_df, aes(x = Var1, y = Var2, label = sig), size = 4) +
  
  # 下三角：相关性数值
  geom_text(data = lower_df, aes(x = Var1, y = Var2, label = round(cor, 2)), size = 4) +
  
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        panel.grid = element_blank()) +
  labs(title = "Test Cohort: Correlation Between HubGenes",
       x = NULL, y = NULL)

dev.off()

png("02.Test_correlation.png", width = 6, height = 5,units = 'in',res = 300,bg = 'white')
ggplot() +
  # 下三角：相关系数填色
  geom_tile(data = lower_df, aes(x = Var1, y = Var2, fill = cor), color = "white") +
  scale_fill_gradientn(colors = c("#377eb8", "white", "#e41a1c"),
                       limits = c(-1, 1),
                       name = "Correlation") +
  ggnewscale::new_scale_fill() +  # 🟢 添加新的 fill 变量映射
  # 上三角：用 p 值着色（-log10(p)）
  geom_tile(data = upper_df, aes(x = Var1, y = Var2, fill = -log10(pval + 1e-20)), color = "white") +
  scale_fill_gradientn(colors = c("white", "#c7e9c0", "#74c476", "#238b45"),
                       name = "-log10(p)") +
  
  # 上三角：显著性标记
  geom_text(data = upper_df, aes(x = Var1, y = Var2, label = sig), size = 4) +
  
  # 下三角：相关性数值
  geom_text(data = lower_df, aes(x = Var1, y = Var2, label = round(cor, 2)), size = 4) +
  
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        panel.grid = element_blank()) +
  labs(title = "Test Cohort: Correlation Between HubGenes",
       x = NULL, y = NULL)

dev.off()

