rm(list = ls())
setwd('/home/data/project/54.LZPJ-3001-6/')
file.dir = '06_roc'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
#训练集------
library(magrittr)
library(ggplot2)
dat <- read.csv('../00_rawdata/dat.GSE65204.csv',row.names = 1)
group <- read.csv('../00_rawdata/group.GSE65204.csv')
table(group$group)
dat <- dat[,group$sample]
hubgene <- read.csv('../05_ML/hubgene.csv')

group$group = factor(group$group, levels = c("Asthma", "Control"),labels = c('Childhood_Asthma','Control'))
hub_exp<-dat[hubgene$symbol,]
library(ggplot2)
library(pROC)
library(ROCR)
roc_list <- list()
for (gene in rownames(hub_exp)) {
  exp_vec <- as.numeric(hub_exp[gene, ])
  response <- as.numeric(group$group == "Childhood_Asthma")  
  
  roc_obj <- pROC::roc(response, exp_vec, quiet = TRUE)
  roc_list[[gene]] <- roc_obj
}
library(IOBR)
colors <- c("#E69F00", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
pdf(file = '01.Train_roc.pdf',width = 5,height = 5)
plot(roc_list[[1]], col = colors[1], lwd = 2, main = "Train Cohort: ROC curves", legacy.axes = TRUE)

for (i in 2:length(roc_list)) {
  plot(roc_list[[i]], col = colors[i], lwd = 2, add = TRUE)
}

legend_text <- sapply(names(roc_list), function(gene) {
  auc_val <- round(auc(roc_list[[gene]]), 2)
  paste0(gene, " (AUC=", auc_val, ")")
})
legend("bottomright", legend = legend_text, col = colors, lwd = 2, cex = 0.8)
dev.off()

png(file = '01.Train_roc.png',width = 5,height = 5,units = 'in',res = 300,bg = 'white')
plot(roc_list[[1]], col = colors[1], lwd = 2, main = "Train Cohort: ROC curves", legacy.axes = TRUE)

# 添加其余曲线
for (i in 2:length(roc_list)) {
  plot(roc_list[[i]], col = colors[i], lwd = 2, add = TRUE)
}

# 添加图例（显示 AUC）
legend_text <- sapply(names(roc_list), function(gene) {
  auc_val <- round(auc(roc_list[[gene]]), 2)
  paste0(gene, " (AUC=", auc_val, ")")
})
legend("bottomright", legend = legend_text, col = colors, lwd = 2, cex = 0.8)
dev.off()

#验证集------
library(magrittr)
library(ggplot2)
dat <- read.csv('../00_rawdata/dat.GSE19187.csv',row.names = 1)
group <- read.csv('../00_rawdata/group.GSE19187.csv')
table(group$group)
dat <- dat[,group$sample]
hubgene <- read.csv('../05_ML/hubgene.csv')
group$group = factor(group$group, levels = c("Asthma", "Control"))
hub_exp<-dat[hubgene$symbol,]
library(ggplot2)
library(pROC)
library(ROCR)
roc_list <- list()
for (gene in rownames(hub_exp)) {
  exp_vec <- as.numeric(hub_exp[gene, ])
  response <- as.numeric(group$group == "Asthma")  
  
  roc_obj <- pROC::roc(response, exp_vec, quiet = TRUE)
  roc_obj[["auc"]] <- ifelse(roc_obj[["auc"]] > 0.5,roc_obj[["auc"]],1-roc_obj[["auc"]])
  roc_list[[gene]] <- roc_obj
}

library(IOBR)
colors <- c("#E69F00", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
pdf(file = '02.Test_roc.pdf',width = 5,height = 5)
plot(roc_list[[1]], col = colors[1], lwd = 2, main = "Test Cohort: ROC curves", legacy.axes = TRUE)

# 添加其余曲线
for (i in 2:length(roc_list)) {
  plot(roc_list[[i]], col = colors[i], lwd = 2, add = TRUE)
}

legend_text <- sapply(names(roc_list), function(gene) {
  auc_val <- round(auc(roc_list[[gene]]), 2)
  auc_val <- ifelse(auc_val > 0.5,auc_val,1-auc_val)
  paste0(gene, " (AUC=", auc_val, ")")
})
legend("bottomright", legend = legend_text, col = colors, lwd = 2, cex = 0.8)
dev.off()

png(file = '02.Test_roc.png',width = 5,height = 5,units = 'in',res = 300,bg = 'white')
plot(roc_list[[1]], col = colors[1], lwd = 2, main = "Test Cohort: ROC curves", legacy.axes = TRUE)

for (i in 2:length(roc_list)) {
  plot(roc_list[[i]], col = colors[i], lwd = 2, add = TRUE)
}

legend_text <- sapply(names(roc_list), function(gene) {
  auc_val <- round(auc(roc_list[[gene]]), 2)
  auc_val <- ifelse(auc_val > 0.5,auc_val,1-auc_val)
  paste0(gene, " (AUC=", auc_val, ")")
})
legend("bottomright", legend = legend_text, col = colors, lwd = 2, cex = 0.8)
dev.off()
