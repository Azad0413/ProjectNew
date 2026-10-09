rm(list = ls())
library(tidyverse)
library(DESeq2)
library(ggVolcano)
library(ggvenn)

df.count = read.delim2("cleandata/count.tsv", row.names = 1) |> 
  apply(c(1,2), as.integer) |> as.data.frame()
grps = ifelse(str_ends(colnames(df.count), ".1.."),"normal","tumor")
df.meta = data.frame(sample = colnames(df.count), group = grps)
df.meta$group = factor(df.meta$group, levels = c("normal","tumor"))

dds = DESeqDataSetFromMatrix(df.count, df.meta, design = ~group)
dds = DESeq(dds)
res.dds = results(dds) %>% as.data.frame
res.dds = na.omit(res.dds)
res.dds = res.dds[c(2,5,6)]
colnames(res.dds) = c("logFC","P.Value","adj.P.Val")
res.dds = cbind(Gene = rownames(res.dds), res.dds)

res.dds$Direction = ifelse(res.dds$logFC > 1, "Up", "No")
res.dds$Direction = ifelse(res.dds$logFC < -1, "Down", res.dds$Direction)
res.dds$Direction = ifelse(res.dds$adj.P.Val < 0.05, res.dds$Direction, "No")
res.dds.out = res.dds[res.dds$Direction != "No",]
write.table(res.dds.out, "02.DEG/01.TCGA.DEG.xls", quote = F, sep = "\t", 
            row.names = F)

df.plot = add_regulate(res.dds, log2FC_name = "logFC",
                       fdr_name = "adj.P.Val", log2FC = 1, fdr = 0.05)
df.plot$regulate[which(df.plot$regulate == "Normal")] = "No"
df.plot = df.plot[order(df.plot$padj, decreasing = F),]
top.gene = c(
  subset(df.plot, Direction == "Up")[1:10,"Gene"],
  subset(df.plot, Direction == "Down")[1:10,"Gene"]
)
ggvolcano(df.plot, x = "log2FoldChange", y = "padj", log2FC_cut = 1,
          label = "Gene", label_number = 10, legend_position = 'UL',
          output = FALSE, custom_label = top.gene) + 
  scale_color_manual(values = c("deepskyblue","grey30","salmon"))
ggsave("02.DEG/02.TCGA.DEG.png", width = 6, height = 6, units = "in", 
       dpi = 300, bg = "white")
ggsave("02.DEG/02.TCGA.DEG.pdf", width = 6, height = 6, units = "in", 
       dpi = 300, bg = "white")



rm(list = ls())
df.count = read.delim2("cleandata/count.tsv", row.names = 1) |> 
  apply(c(1,2), as.integer) |> as.data.frame()
df.meta = read.delim2("interm/TCGA.cluster.xls")
df.meta = df.meta[df.meta$cluster != "cluster_1",]
df.meta$group = factor(df.meta$cluster, levels = c("cluster_2","cluster_3"))
df.count = df.count[df.meta$sample]

dds = DESeqDataSetFromMatrix(df.count, df.meta, design = ~group)
dds = DESeq(dds)
res.dds = results(dds) %>% as.data.frame
res.dds = na.omit(res.dds)
res.dds = res.dds[c(2,5,6)]
colnames(res.dds) = c("logFC","P.Value","adj.P.Val")
res.dds = cbind(Gene = rownames(res.dds), res.dds)

res.dds$Direction = ifelse(res.dds$logFC > 1, "Up", "No")
res.dds$Direction = ifelse(res.dds$logFC < -1, "Down", res.dds$Direction)
res.dds$Direction = ifelse(res.dds$adj.P.Val < 0.05, res.dds$Direction, "No")
res.dds.out = res.dds[res.dds$Direction != "No",]
write.table(res.dds.out, "02.DEG/03.TCGA.Cluster.DEG.xls", quote = F, sep = "\t", 
            row.names = F)

df.plot = add_regulate(res.dds, log2FC_name = "logFC",
                       fdr_name = "adj.P.Val", log2FC = 1, fdr = 0.05)
df.plot$regulate[which(df.plot$regulate == "Normal")] = "No"
df.plot = df.plot[order(df.plot$padj, decreasing = F),]
top.gene = c(
  subset(df.plot, Direction == "Up")[1:10,"Gene"],
  subset(df.plot, Direction == "Down")[1:10,"Gene"]
)
ggvolcano(df.plot, x = "log2FoldChange", y = "padj", log2FC_cut = 1,
          label = "Gene", label_number = 10, legend_position = 'UR',
          output = FALSE, custom_label = top.gene) + 
  scale_color_manual(values = c("deepskyblue","grey30","salmon"))
ggsave("02.DEG/04.TCGA.Cluster.DEG.png", width = 6, height = 6, units = "in", 
       dpi = 300, bg = "white")
ggsave("02.DEG/04.TCGA.Cluster.DEG.pdf", width = 6, height = 6, units = "in", 
       dpi = 300, bg = "white")




rm(list = ls())
df.1 = read.delim2("02.DEG/01.TCGA.DEG.xls")
df.2 = read.delim2("02.DEG/03.TCGA.Cluster.DEG.xls")

ggvenn(list(`DEG\nbetween groups` = df.1$Gene, `DEG\nbetween clusters` = df.2$Gene)) + 
  ylim(c(-1,1.5))
ggsave("02.DEG/05.Venn.png", width = 6, height = 6, units = "in", dpi = 300, bg = "white")
ggsave("02.DEG/05.Venn.pdf", width = 6, height = 6, units = "in", dpi = 300, bg = "white")

inter = intersect(df.1$Gene, df.2$Gene)
write.table(data.frame(gene = inter), "02.DEG/06.Intersection.xls", quote = F, 
            sep = "\t", row.names = F)




