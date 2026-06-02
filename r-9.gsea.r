rm(list = ls())
setwd('/home/data/project/54.LZPJ-3001-6/')
file.dir = '09_GSEA'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(limma)
library(tidyverse)

dat <- read.csv('../00_rawdata/dat.GSE65204.csv',row.names = 1)
group <- read.csv('../00_rawdata/group.GSE65204.csv')
hubgene <- read.csv('../05_ML/hubgene.csv')
hub.exp<-dat[hubgene$symbol,]
hub.exp<-hub.exp[order(rownames(hub.exp),decreasing = F),] 
kegggmt<-read.gmt("../../../files/GSVA/c2.cp.kegg.v7.5.1.symbols.gmt") 
results <- list()
for (i in c(1:5)){
  train<-t(hub.exp)
  group=as.vector(ifelse(train[,i]>median(train[,i]),'high','low'))  
  group <- factor(group,levels = c("high","low")) 
  condition<-data.frame(group)
  rownames(condition)<-rownames(train)
  design<-model.matrix(~0+group)   
  contrast.matrix<-makeContrasts(contrasts = "grouphigh-grouplow", levels = design) 
  fit <- lmFit(dat,design)     
  fit1 <- contrasts.fit(fit, contrast.matrix)    
  fit2 <- eBayes(fit1)  
  tempOutput <- topTable(fit2, coef=1, n=nrow(fit2),adjust="fdr")  
  genelist<-data.frame(SYMBOL=rownames(tempOutput ),logFC=tempOutput $logFC)
  gene=bitr(genelist$SYMBOL,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") 
  gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
  gene_df <- data.frame(ENTREZID=gene$ENTREZID ,
                        SYMBOL = gene$SYMBOL) 
  gene_df <- merge(gene_df,genelist,by="SYMBOL")
  geneList<-gene_df$logFC 
  names(geneList)=gene_df $SYMBOL 
  geneList=sort(geneList,decreasing = T)
  set.seed(123)
  KEGG<-GSEA(geneList,TERM2GENE = kegggmt,pvalueCutoff = 1) 
  results[[i]] <- KEGG
  sortKEGG<-data.frame(KEGG)
  sortKEGG <- sortKEGG[which(sortKEGG$p.adjust < 0.05),]
  sortKEGG<-sortKEGG[order(sortKEGG$p.adjust, decreasing = F),]
  write.table(sortKEGG,file=paste('0',i,'.',rownames(hub.exp)[i],'.xls',sep=''),sep="\t",quote=F,row.names=F)

  ifelse(length(rownames(sortKEGG) )> 5,
         paths <- rownames(sortKEGG[c(1:5),]),
         paths <- rownames(sortKEGG)
  )
  # 定义调色板
  my_colors <- c("#E69F00", "#009E73", "#0072B2", "#D55E00", "#CC79A7") # 可替换为你喜欢的颜色
  my_colors <- my_colors[1:length(paths)]
  #特定通路作图
  # trace(gseaplot2,edit=T)
  p = gseaplot2(KEGG,paths,base_size =10,color = my_colors,rel_heights = c(1.5, 0.3, 0.5),title=rownames(hub.exp)[i])
  fn1 = paste0("0", i, ".", rownames(hub.exp)[i], ".png")
  fn2 = paste0("0", i, ".", rownames(hub.exp)[i],".pdf")
  png(fn1,w=8,h=6,units = "in",res = 300,bg='white')
  print(p)
  dev.off()
  pdf(fn2,w=8,h=6)
  print(p)
  dev.off()
}

##共同通路------
library(dplyr)
library(ggplot2)
library(readr)

files <- list.files(pattern = "^0[0-9]\\..*\\.xls$")
all_data <- lapply(files, read.delim, stringsAsFactors = FALSE)

# 给每个结果加上基因名标签
gene_names <- gsub("^0[0-9]\\.(.*?)\\.xls$", "\\1", files)
for (i in seq_along(all_data)) {
  all_data[[i]]$Gene <- gene_names[i]
}

df_all <- bind_rows(all_data)

common_pathways <- df_all %>%
  filter(p.adjust < 0.05) %>%
  group_by(ID) %>%
  summarise(n_gene = n_distinct(Gene), .groups = "drop") %>%
  filter(n_gene == length(unique(df_all$Gene))) %>%
  pull(ID)

df_common <- df_all %>%
  filter(ID %in% common_pathways)

# 计算 -log10(p.adjust) 作为显著性
df_common$logP <- -log10(df_common$p.adjust)
write.csv(df_common,file = 'common_pathway.csv',row.names = F)

# 气泡图
p <- ggplot(df_common, aes(x = Gene, y = Description, color = enrichmentScore, size = logP)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = "#1177b0", high = "#c04851") + # 色盲友好
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Hub Gene",
       y = "KEGG Pathway",
       color = "enrichmentScore",
       size = "-log10(p.adjust)",
       title = "Common Enriched Pathways of Hub Genes")
p
ggsave(filename = '06.Common.pathway.pdf',p,width = 8,height = 6)
ggsave(filename = '06.Common.pathway.png',p,width = 8,height = 6,units = 'in',dpi = 300,bg = 'white')
