rm(list = ls())
library(tidyverse)
library(survival)
library(survminer)
library(org.Hs.eg.db)
library(ggsci)
library(ConsensusClusterPlus)
library(msigdbr)
library(GSVA)

df.exp = read.delim2('cleandata/fpkm.tsv', row.names = 1) |>
  type.convert(as.is = T)
gs = c("MFN1","MFN2","ITPR1","HSPA9","VDAC1","FIS1","BCAP31","VAPB","RMDN3","PDZD8","FKBP8")
df.exp = df.exp[gs,]

grp = ifelse(str_ends(colnames(df.exp), "1.."), "Control", "Tumor")
p = apply(df.exp, 1, \(x) wilcox.test(x~grp)$p.value)
df.diff = data.frame(gene = rownames(df.exp),
                     p.val = p,
                     p.adj = p.adjust(p, method = 'BH'))
df.diff = df.diff[order(df.diff$p.adj, decreasing = F),]
write.table(df.diff, "01.Gene/01.Selected.Diff.xls", quote = F, sep = "\t", row.names = F)

df.exp = cbind(gene = rownames(df.exp), df.exp)
df.exp = reshape2::melt(df.exp, id.var = "gene", variable.name = "sample")
df.exp$value = log2(df.exp$value + 1)
df.exp$group = ifelse(str_ends(df.exp$sample, "1.."), "Control", "Tumor")
df.exp$p = df.diff$p.adj[match(df.exp$gene, df.diff$gene)]
df.exp$lab = paste0(df.exp$gene, "\np.adj = ", signif(df.exp$p,3))

ggplot(df.exp, aes(x = group, y = value)) + 
  geom_boxplot(aes(color = group), show.legend = F) + 
  scale_color_aaas() + 
  facet_wrap(~lab, nrow = 2, scale = "free") + 
  theme_bw() + 
  labs(x = "Group", y = "Expression")
ggsave('01.Gene/02.Selected.Diff.png', width = 10, height = 6, units = 'in', dpi = 300, bg = 'white')
ggsave('01.Gene/02.Selected.Diff.pdf', width = 10, height = 6, units = 'in', dpi = 300, bg = 'white')



rm(list = ls())
df.exp = read.delim2('cleandata/fpkm.tsv', row.names = 1) |>
  type.convert(as.is = T) |> 
  apply(c(1,2), \(x) log2(x+1)) |> 
  as.data.frame()
df.surv = read.delim2('cleandata/survival.tsv') |>
  type.convert(as.is = T)
gs = c("MFN1","MFN2","ITPR1","HSPA9","VDAC1","FIS1","BCAP31","VAPB","RMDN3","PDZD8","FKBP8")
df.exp = df.exp[gs,df.surv$sample] |> t() |> as.data.frame()

for (g in gs) {
  df.surv$gene = df.exp[[g]]
  cut = survminer::surv_cutpoint(df.surv, time = "OS.time", event = "OS", variables = "gene")
  cut = cut$cutpoint$cutpoint
  df.surv$group = ifelse(df.surv$gene > cut, "High", "Low")
  surv.fit = survfit(
    Surv(time = OS.time, event = OS) ~ group,
    data = df.surv
  )
  cut.lab = paste0("cutpoint = ", round(cut, 2))
  p.surv = ggsurvplot(
    surv.fit,
    data = df.surv,
    conf.int = F,
    risk.table = TRUE,
    break.x.by = 800,
    pval = T,
    legend.labs = c("Low", "High"),
    palette = pal_aaas()(2)
  ) +
    xlab("Time (Day)")
  p.surv$plot = p.surv$plot +
    guides(color = guide_legend(title = cut.lab)) +
    ggtitle("TCGA",g) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5))
  p.surv$table = p.surv$table + ylab(g)
  
  fn1 = paste0("01.Gene/03.KM/",g,".KM.png")
  fn2 = paste0("01.Gene/03.KM/",g,".KM.pdf")
  
  ggsave(fn1, print(p.surv, newpage = F), width = 6, height = 6, units = 'in', dpi = 300, bg = 'white')
  ggsave(fn2, print(p.surv, newpage = F), width = 6, height = 6, units = 'in', dpi = 300, bg = 'white')
}



rm(list = ls())
df.exp = read.delim2('cleandata/fpkm.tsv', row.names = 1) |>
  type.convert(as.is = T)
df.surv = read.delim2('cleandata/survival.tsv') |>
  type.convert(as.is = T)
gs = c("HSPA9","VDAC1","FIS1","BCAP31","PDZD8")
df.exp = df.exp[gs,df.surv$sample]
df.exp = apply(df.exp, c(1,2), \(x) log2(x+1))
df.exp = t(scale(t(df.exp)))

res = ConsensusClusterPlus(as.matrix(df.exp), 
                           title = "interm/TCGA.CC", 
                           distance = "euclidean",
                           plot = "pdf",
                           innerLinkage = "average",
                           finalLinkage = "average",
                           clusterAlg = "pam",
                           maxK = 6, 
                           seed = 2025)
df.clst = data.frame(sample = colnames(df.exp), cluster = res[[3]]$consensusClass)
df.clst$cluster = paste0("cluster_", df.clst$cluster)
write.table(df.clst, "interm/TCGA.cluster.xls", quote = F, sep = "\t", row.names = F)


rm(list = ls())
df.surv = read.delim2('cleandata/survival.tsv') |>
  type.convert(as.is = T)
df.clst = read.delim2('interm/TCGA.cluster.xls')
df.full = merge(df.surv, df.clst)
df.full = df.full[order(df.full$cluster, decreasing = F),]

surv.fit = survfit(
    Surv(time = OS.time, event = OS) ~ cluster,
    data = df.full
)
p.surv = ggsurvplot(
  surv.fit,
  data = df.full,
  conf.int = F,
  risk.table = TRUE,
  break.x.by = 800,
  pval = T,
  palette = pal_bmj()(5), 
  legend.labs = unique(df.full$cluster)
) +
  xlab("Time (Day)")
p.surv$plot = p.surv$plot +
  ggtitle("TCGA Clusters") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5))
p.surv$table = p.surv$table + ylab("Cluster")
p.surv
ggsave('01.Gene/06.KM.Cluster.png', print(p.surv, newpage = F), width = 6, height = 6, units = 'in', dpi = 300, bg = 'white')
ggsave('01.Gene/06.KM.Cluster.pdf', print(p.surv, newpage = F), width = 6, height = 6, units = 'in', dpi = 300, bg = 'white')



rm(list = ls())
df.exp = read.delim2('cleandata/fpkm.tsv', row.names = 1) |>
  type.convert(as.is = T)
df.grp = read.delim2("interm/TCGA.cluster.xls")
gs = c("HSPA9","VDAC1","FIS1","BCAP31","PDZD8")
df.exp = df.exp[gs,df.grp$sample]

p = apply(df.exp, 1, \(x) kruskal.test(x~df.grp$cluster)$p.value)
df.diff = data.frame(gene = rownames(df.exp),
                     p.val = p,
                     p.adj = p.adjust(p, method = 'BH'))
df.diff = df.diff[order(df.diff$p.adj, decreasing = F),]
write.table(df.diff, '01.Gene/07.Cluster.Gene.xls', sep = '\t', quote = F, row.names = F)

df.exp$gene = rownames(df.exp)
df.exp = reshape2::melt(df.exp, id.var = "gene", variable.name = "sample")
df.exp$cluster = df.grp$cluster[match(df.exp$sample, df.grp$sample)]
df.exp$p = df.diff$p.adj[match(df.exp$gene, df.diff$gene)]
df.exp$lab = paste0(df.exp$gene, "\np.adj = ", signif(df.exp$p, 3))
df.exp$value = log2(df.exp$value+1)

ggplot(df.exp, aes(x = cluster, y = value)) + 
  geom_boxplot(aes(color = cluster)) + 
  scale_color_bmj() + 
  facet_wrap(~lab, nrow = 1, scale = "free") + 
  theme_bw() + 
  labs(x = NULL, y = "Expression") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
ggsave('01.Gene/08.Cluter.Gene.png', width = 10, height = 5, units = 'in', dpi = 300, bg = 'white')
ggsave('01.Gene/08.Cluter.Gene.pdf', width = 10, height = 5, units = 'in', dpi = 300, bg = 'white')


rm(list = ls())
df.exp = read.delim2("cleandata/fpkm.tsv", row.names = 1) |> 
  type.convert(as.is = T) |> 
  apply(c(1,2), \(x) log2(x+1)) |> 
  as.data.frame()
gs = msigdbr()
gs = gs[gs$gs_subcollection %in% c("GO:BP","GO:CC","GO:MF"),]
gs = split(gs$gene_symbol, gs$gs_name)
df = gsva(ssgseaParam(as.matrix(df.exp), gs))
df = as.data.frame(df)
df = cbind(pathway = rownames(df), df)
write.table(df, "interm/tcga.go.tsv", sep = "\t", quote = F, row.names = F)



rm(list = ls())
df.exp = read.delim2('interm/tcga.go.tsv', row.names = 1) |>
  type.convert(as.is = T)
df.grp = read.delim2("interm/TCGA.cluster.xls")
df.exp = df.exp[df.grp$sample]
df.exp = df.exp[str_starts(rownames(df.exp), "GOBP"),]

p = apply(df.exp, 1, \(x) kruskal.test(x~df.grp$cluster)$p.value)
df.diff = data.frame(GOBP = rownames(df.exp),
                     p.val = p,
                     p.adj = p.adjust(p, method = 'BH'))
df.diff = df.diff[order(df.diff$p.adj, decreasing = F),]
write.table(df.diff, '01.Gene/09.Cluster.GOBP.xls', sep = '\t', quote = F, row.names = F)


df.diff = na.omit(df.diff)
df.diff = df.diff[df.diff$p.adj < 0.05, ]
df.diff = df.diff[1:30,]
df.diff$lab = ifelse(df.diff$p.adj <= 0.001, "***", 
                     ifelse(df.diff$p.adj <= 0.01, "**", "*"))
df.diff$lab= paste0(df.diff$GOBP, "(", df.diff$lab, ")")
df.exp = df.exp[df.diff$GOBP,]
rownames(df.exp) = df.diff$lab

df.grp = df.grp[order(df.grp$cluster, decreasing = F),]
rownames(df.grp) = df.grp$sample
df.exp = df.exp[df.grp$sample]

pal_group = pal_bmj()(3)
names(pal_group) = unique(df.grp$cluster)
gaps_col = cumsum(table(df.grp$cluster))

p = pheatmap::pheatmap(
  as.matrix(df.exp),
  color = colorRampPalette(c("blue","white","red"))(100),
  border_color = "white",
  scale = "row",
  show_colnames = F, 
  cluster_cols = F,
  gaps_col = gaps_col,
  annotation_col = df.grp['cluster'],
  annotation_colors = list(cluster = pal_group)
)
p = ggplotify::as.ggplot(p)
ggsave('01.Gene/10.Cluster.GOBP.png', p, width = 12, height = 6, units = 'in', dpi = 300, bg = 'white')
ggsave('01.Gene/10.Cluster.GOBP.pdf', p, width = 12, height = 6, units = 'in', dpi = 300, bg = 'white')



rm(list = ls())
df.exp = read.delim2('interm/tcga.go.tsv', row.names = 1) |>
  type.convert(as.is = T)
df.grp = read.delim2("interm/TCGA.cluster.xls")
df.exp = df.exp[df.grp$sample]
df.exp = df.exp[str_starts(rownames(df.exp), "GOCC"),]

p = apply(df.exp, 1, \(x) kruskal.test(x~df.grp$cluster)$p.value)
df.diff = data.frame(GOCC = rownames(df.exp),
                     p.val = p,
                     p.adj = p.adjust(p, method = 'BH'))
df.diff = df.diff[order(df.diff$p.adj, decreasing = F),]
write.table(df.diff, '01.Gene/11.Cluster.GOCC.xls', sep = '\t', quote = F, row.names = F)


df.diff = na.omit(df.diff)
df.diff = df.diff[df.diff$p.adj < 0.05, ]
df.diff = df.diff[1:30,]
df.diff$lab = ifelse(df.diff$p.adj <= 0.001, "***", 
                     ifelse(df.diff$p.adj <= 0.01, "**", "*"))
df.diff$lab= paste0(df.diff$GOCC, "(", df.diff$lab, ")")
df.exp = df.exp[df.diff$GOCC,]
rownames(df.exp) = df.diff$lab

df.grp = df.grp[order(df.grp$cluster, decreasing = F),]
rownames(df.grp) = df.grp$sample
df.exp = df.exp[df.grp$sample]

pal_group = pal_bmj()(3)
names(pal_group) = unique(df.grp$cluster)
gaps_col = cumsum(table(df.grp$cluster))

p = pheatmap::pheatmap(
  as.matrix(df.exp),
  color = colorRampPalette(c("blue","white","red"))(100),
  border_color = "white",
  scale = "row",
  show_colnames = F, 
  cluster_cols = F,
  gaps_col = gaps_col,
  annotation_col = df.grp['cluster'],
  annotation_colors = list(cluster = pal_group)
)
p = ggplotify::as.ggplot(p)
ggsave('01.Gene/12.Cluster.GOCC.png', p, width = 12, height = 6, units = 'in', dpi = 300, bg = 'white')
ggsave('01.Gene/12.Cluster.GOCC.pdf', p, width = 12, height = 6, units = 'in', dpi = 300, bg = 'white')




rm(list = ls())
df.exp = read.delim2('interm/tcga.go.tsv', row.names = 1) |>
  type.convert(as.is = T)
df.grp = read.delim2("interm/TCGA.cluster.xls")
df.exp = df.exp[df.grp$sample]
df.exp = df.exp[str_starts(rownames(df.exp), "GOMF"),]

p = apply(df.exp, 1, \(x) kruskal.test(x~df.grp$cluster)$p.value)
df.diff = data.frame(GOMF = rownames(df.exp),
                     p.val = p,
                     p.adj = p.adjust(p, method = 'BH'))
df.diff = df.diff[order(df.diff$p.adj, decreasing = F),]
write.table(df.diff, '01.Gene/13.Cluster.GOMF.xls', sep = '\t', quote = F, row.names = F)


df.diff = na.omit(df.diff)
df.diff = df.diff[df.diff$p.adj < 0.05, ]
df.diff = df.diff[1:30,]
df.diff$lab = ifelse(df.diff$p.adj <= 0.001, "***", 
                     ifelse(df.diff$p.adj <= 0.01, "**", "*"))
df.diff$lab= paste0(df.diff$GOMF, "(", df.diff$lab, ")")
df.exp = df.exp[df.diff$GOMF,]
rownames(df.exp) = df.diff$lab

df.grp = df.grp[order(df.grp$cluster, decreasing = F),]
rownames(df.grp) = df.grp$sample
df.exp = df.exp[df.grp$sample]

pal_group = pal_bmj()(3)
names(pal_group) = unique(df.grp$cluster)
gaps_col = cumsum(table(df.grp$cluster))

p = pheatmap::pheatmap(
  as.matrix(df.exp),
  color = colorRampPalette(c("blue","white","red"))(100),
  border_color = "white",
  scale = "row",
  show_colnames = F, 
  cluster_cols = F,
  gaps_col = gaps_col,
  annotation_col = df.grp['cluster'],
  annotation_colors = list(cluster = pal_group)
)
p = ggplotify::as.ggplot(p)
ggsave('01.Gene/14.Cluster.GOMF.png', p, width = 16, height = 6, units = 'in', dpi = 300, bg = 'white')
ggsave('01.Gene/14.Cluster.GOMF.pdf', p, width = 16, height = 6, units = 'in', dpi = 300, bg = 'white')




rm(list = ls())
df.exp = read.delim2("cleandata/fpkm.tsv", row.names = 1) |> 
  type.convert(as.is = T) |> 
  apply(c(1,2), \(x) log2(x+1)) |> 
  as.data.frame()
gs = msigdbr()
gs = gs[gs$gs_subcollection %in% c("CP:KEGG_LEGACY"),]
gs = split(gs$gene_symbol, gs$gs_name)
df = gsva(ssgseaParam(as.matrix(df.exp), gs))
df = as.data.frame(df)
df = cbind(pathway = rownames(df), df)
write.table(df, "interm/tcga.kegg.tsv", sep = "\t", quote = F, row.names = F)




rm(list = ls())
df.exp = read.delim2('interm/tcga.kegg.tsv', row.names = 1) |>
  type.convert(as.is = T)
df.grp = read.delim2("interm/TCGA.cluster.xls")
df.exp = df.exp[df.grp$sample]

p = apply(df.exp, 1, \(x) kruskal.test(x~df.grp$cluster)$p.value)
df.diff = data.frame(KEGG = rownames(df.exp),
                     p.val = p,
                     p.adj = p.adjust(p, method = 'BH'))
df.diff = df.diff[order(df.diff$p.adj, decreasing = F),]
write.table(df.diff, '01.Gene/15.Cluster.KEGG.xls', sep = '\t', quote = F, row.names = F)


df.diff = na.omit(df.diff)
df.diff = df.diff[df.diff$p.adj < 0.05, ]
df.diff = df.diff[1:30,]
df.diff$lab = ifelse(df.diff$p.adj <= 0.001, "***", 
                     ifelse(df.diff$p.adj <= 0.01, "**", "*"))
df.diff$lab= paste0(df.diff$KEGG, "(", df.diff$lab, ")")
df.exp = df.exp[df.diff$KEGG,]
rownames(df.exp) = df.diff$lab

df.grp = df.grp[order(df.grp$cluster, decreasing = F),]
rownames(df.grp) = df.grp$sample
df.exp = df.exp[df.grp$sample]

pal_group = pal_bmj()(3)
names(pal_group) = unique(df.grp$cluster)
gaps_col = cumsum(table(df.grp$cluster))

p = pheatmap::pheatmap(
  as.matrix(df.exp),
  color = colorRampPalette(c("blue","white","red"))(100),
  border_color = "white",
  scale = "row",
  show_colnames = F, 
  cluster_cols = F,
  gaps_col = gaps_col,
  annotation_col = df.grp['cluster'],
  annotation_colors = list(cluster = pal_group)
)
p = ggplotify::as.ggplot(p)
ggsave('01.Gene/16.Cluster.KEGG.png', p, width = 12, height = 6, units = 'in', dpi = 300, bg = 'white')
ggsave('01.Gene/16.Cluster.KEGG.pdf', p, width = 12, height = 6, units = 'in', dpi = 300, bg = 'white')


