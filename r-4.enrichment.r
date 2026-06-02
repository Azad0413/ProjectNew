rm(list = ls())
setwd('/home/data/project/54.LZPJ-3001-6/')
file.dir = '04_enrichment'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
candi <- read.csv('../02_DEPCDs/depcds.csv')
gene_transform <- bitr(candi$symbol,
                       fromType = "SYMBOL",
                       toType = c("ENTREZID",'ENSEMBL'),
                       OrgDb = "org.Hs.eg.db")
ego <- enrichGO(gene = gene_transform$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 1,
                readable = TRUE)
write.csv(ego,file = "GO.csv",row.names = F)
go_bar<-barplot(ego, showCategory=5, split="ONTOLOGY",label_format = 80) +
  facet_grid(ONTOLOGY ~ ., scales = "free",space = 'free')+scale_fill_gradient(low = '#EFC000FF',high = '#33a02c')
go_bar
ggsave('01.GO_bar.png',go_bar,width =8,height = 6,units = 'in',dpi = 300,bg = 'white')
ggsave('01.GO_bar.pdf',go_bar,width =8,height = 6)

kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 1)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.csv(kk@result,file = "KEGG.csv",row.names = F)
kk_dot <- dotplot(kk, showCategory=15,label_format = 50)
kk_dot

library(ggrepel)
kk_result <- kk@result[which(kk@result$pvalue < 0.05),]#[c(1:15),]
colors <- c(
  "#F6D6D6", "#F3C1B6", "#F2AC9D", "#EFA58C", "#EC9E7A", "#E89969", "#E4A671",
  "#E0B479", "#DCC281", "#D6CE8B", "#C9D694", "#BBD99F", "#ACDAA9", "#9BDAB4",
  "#8AD5BF", "#78CED1", "#6CBFD9", "#62AEDB", "#5A9CDB", "#5688D2", "#5975C3",
  "#6463AE", "#72579B", "#825088", "#904F77", "#9C5169", "#A5555F"
)

colors <- c(
  "#F6D6D6",  "#EC9E7A", "#E89969", "#E4A671",
  "#E0B479",  "#BBD99F", "#ACDAA9", "#9BDAB4",
  "#8AD5BF", "#78CED1", "#6CBFD9","#5688D2", "#5975C3",
  "#6463AE", "#72579B", "#825088", "#904F77", "#9C5169", "#A5555F"
)
p <- ggplot(kk_result,aes(x=-log10(pvalue),y=Description))+
  geom_point(aes(size=Count,color=Description),alpha=2)+
  scale_size(range = c(10,25))+
  theme_bw()+
  theme(legend.position = c('none'))+
  geom_text_repel(aes(label=Description),size=3,segment.color='black',show.legend = T)+
  theme(axis.text.y = element_blank())  +
  scale_color_manual(values = colors) 
p   
ggsave(filename = '02.KEGG_dot.pdf',p,w=6,h=5)
ggsave(filename = '02.KEGG_dot.png',p,w=6,h=5,dpi=300,units = 'in')
