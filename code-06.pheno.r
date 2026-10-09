rm(list = ls())
library(tidyverse)
library(ggsci)
library(ggpubr)
library(GSVA)
library(msigdbr)
library(clusterProfiler)
library(GseaVis)

df.exp = read.delim2("cleandata/fpkm.tsv", row.names = 1) |> 
  apply(c(1,2), \(x) log2(as.numeric(x)+1)) |> as.data.frame()
gs = msigdbr(collection = "H")
gs = split(gs$gene_symbol, gs$gs_name)
df = gsva(ssgseaParam(as.matrix(df.exp), gs))
df = as.data.frame(df)
df = cbind(pathway = rownames(df), df)
write.table(df, "interm/tcga.hallmark.tsv", sep = "\t", quote = F, row.names = F)


rm(list = ls())
res.gsva = read.delim2("interm/tcga.hallmark.tsv", row.name = 1) |> 
  type.convert(as.is = T)
rownames(res.gsva) = str_remove(rownames(res.gsva), "HALLMARK_")

df.grp = read.delim2("interm/risk.group.xls")
res.gsva = res.gsva[,df.grp$sample]
df.diff = apply(res.gsva, 1, function(x){
  p = t.test(x~df.grp$risk.group)$p.value
  t = t.test(x~df.grp$risk.group)$statistic[1]
  return(c(p,t))
}) |> t() |> as.data.frame()
df.diff = cbind(rownames(df.diff), df.diff)
colnames(df.diff) = c("HALLMARK","p","t")
df.diff$p.adj = p.adjust(df.diff$p, method = "BH")
df.out = df.diff[c(1,3,2,4)]
df.out = df.out[order(df.out$p.adj, decreasing = F),]
write.table(df.out, "06.Phenotype/01.GSVA.Diff.xls", quote = F, sep = "\t", row.names = F)

df.diff = df.diff[order(df.diff$t, decreasing = F),]
df.diff$HALLMARK = factor(df.diff$HALLMARK, levels = df.diff$HALLMARK)
df.diff$col = ifelse(df.diff$t > 0, "up", "down")
df.diff$col = ifelse(df.diff$p.adj <0.05, df.diff$col, "no")

ggplot(df.diff, aes(x = t, y = HALLMARK, fill = col)) + 
  geom_bar(stat = "identity") + 
  geom_text(data = subset(df.diff, t>0 & p.adj<0.05), mapping = aes(x = -0.3, y = HALLMARK, label = HALLMARK), hjust = 1, size = 3) +
  geom_text(data = subset(df.diff, t<0 & p.adj<0.05), mapping = aes(x = 0.3, y = HALLMARK, label = HALLMARK), hjust = 0, size = 3) +
  geom_text(data = subset(df.diff, t>0 & p.adj>0.05), mapping = aes(x = -0.3, y = HALLMARK, label = HALLMARK), color = "grey", hjust = 1, size = 3) +
  geom_text(data = subset(df.diff, t<0 & p.adj>0.05), mapping = aes(x = 0.3, y = HALLMARK, label = HALLMARK), color = "grey", hjust = 0, size = 3) +
  scale_x_continuous(limits = c(-7,7), breaks = seq(-6,6,2)) +
  ylab(NULL) + xlab("t value of GSVA score") +
  scale_fill_manual(values = c("up"="darkred","no"="grey","down"="darkblue")) +
  theme_bw(base_size = 12) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), legend.position = "none", panel.grid.major = element_line(color = "white"))
ggsave("06.Phenotype/02.GSVA.Diff.png", width = 8, height = 10, units = "in", dpi = 300, bg = "white")
ggsave("06.Phenotype/02.GSVA.Diff.pdf", width = 8, height = 10, units = "in", dpi = 300, bg = "white")


rm(list = ls())
df.risk = read.delim2("interm/risk.group.xls") |> 
  type.convert(as.is = T)
df.exp = read.delim2("cleandata/fpkm.tsv", row.names = 1) |> 
  apply(c(1,2), \(x) log2(as.numeric(x)+1)) |> as.data.frame()
df.exp = df.exp[df.risk$sample] |> t() |> as.data.frame()
df.exp = df.exp[apply(df.exp, 2, sd)!=0]
res = linkET::correlate(df.exp, df.risk['risk'], method = "spearman")
df.r = as.data.frame(res$r)
r = df.r$risk
names(r) = rownames(df.r)
r = r[order(r, decreasing = T)]


df.m = msigdbr()
df.kegg = subset(df.m, gs_subcollection %in% c("CP:KEGG_LEGACY"))[c(9,1)]

res.gsea = GSEA(r, TERM2GENE = df.kegg, pvalueCutoff = 0.05, 
                  seed = 1, pAdjustMethod = "BH", eps = 0)
top.kegg = res.gsea@result
write.table(top.kegg, paste0("06.Phenotype/03.GSEA.KEGG.xls"), quote = F, sep = "\t", row.names = F)

top.kegg = top.kegg[order(top.kegg$p.adjust, decreasing = F),]
up.kegg = subset(top.kegg, NES > 0)$Description[1:10]
p = gseaNb(
  object = res.gsea,
  curveCol = pal_igv()(10),
  geneSetID = up.kegg,
  subPlot = 2,
)
p[[1]] = p[[1]] + ggtitle("Positive correlated KEGG pathways") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave('06.Phenotype/04.GSEA.KEGG.Pos.png', p, width = 8, height = 6, units = 'in', dpi = 300, bg = 'white')
ggsave('06.Phenotype/04.GSEA.KEGG.Pos.pdf', p, width = 8, height = 6, units = 'in', dpi = 300, bg = 'white')

down.kegg = subset(top.kegg, NES < 0)$Description[1:10]
p = gseaNb(
  object = res.gsea,
  curveCol = pal_igv()(10),
  geneSetID = down.kegg,
  subPlot = 2,
)
p[[1]] = p[[1]] + ggtitle("Negative correlated KEGG pathways") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave('06.Phenotype/05.GSEA.KEGG.Neg.png', p, width = 8, height = 6, units = 'in', dpi = 300, bg = 'white')
ggsave('06.Phenotype/05.GSEA.KEGG.Neg.pdf', p, width = 8, height = 6, units = 'in', dpi = 300, bg = 'white')



df.m = msigdbr()
df.kegg = subset(df.m, gs_subcollection %in% c("GO:BP","GO:CC","GO:MF"))[c(9,1)]

res.gsea = GSEA(r, TERM2GENE = df.kegg, pvalueCutoff = 0.05, 
                  seed = 1, pAdjustMethod = "BH", eps = 0)
top.kegg = res.gsea@result
write.table(top.kegg, paste0("06.Phenotype/06.GSEA.GO.xls"), quote = F, sep = "\t", row.names = F)

top.kegg = top.kegg[order(top.kegg$p.adjust, decreasing = F),]
up.kegg = subset(top.kegg, NES > 0)$Description[1:10]
p = gseaNb(
  object = res.gsea,
  curveCol = pal_igv()(10),
  geneSetID = up.kegg,
  subPlot = 2,
)
p[[1]] = p[[1]] + ggtitle("Positive correlated GO pathways") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave('06.Phenotype/07.GSEA.GO.Pos.png', p, width = 8, height = 6, units = 'in', dpi = 300, bg = 'white')
ggsave('06.Phenotype/07.GSEA.GO.Pos.pdf', p, width = 8, height = 6, units = 'in', dpi = 300, bg = 'white')

down.kegg = subset(top.kegg, NES < 0)$Description[1:10]
p = gseaNb(
  object = res.gsea,
  curveCol = pal_igv()(10),
  geneSetID = down.kegg,
  subPlot = 2,
)
p[[1]] = p[[1]] + ggtitle("Negative correlated GO pathways") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave('06.Phenotype/08.GSEA.GO.Neg.png', p, width = 8, height = 6, units = 'in', dpi = 300, bg = 'white')
ggsave('06.Phenotype/08.GSEA.GO.Neg.pdf', p, width = 8, height = 6, units = 'in', dpi = 300, bg = 'white')




rm(list = ls())
library(IOBR)
df.exp = read.delim2("cleandata/fpkm.tsv", row.names = 1) |>
  type.convert(as.is = T)
df.cb = deconvo_cibersort(df.exp, arrays = F, perm = 1)
df.cb = df.cb[1:23]
write.table(df.cb, "interm/cibersort.xls", sep = "\t", quote = F, row.names = F)
df.xc = deconvo_xcell(df.exp, arrays = F)
write.table(df.xc, "interm/xcell.xls", sep = "\t", quote = F, row.names = F)
df.qs = deconvo_quantiseq(df.exp, tumor = T, arrays = F, scale_mrna = F)
write.table(df.qs, "interm/quantiseq.xls", sep = "\t", quote = F, row.names = F)
df.ep = deconvo_epic(df.exp, tumor = T)
write.table(df.ep, "interm/epic.xls", sep = "\t", quote = F, row.names = F)
df.mc = deconvo_mcpcounter(df.exp)
write.table(df.mc, "interm/mcpcounter.xls", sep = "\t", quote = F, row.names = F)


rm(list = ls())
library(ComplexHeatmap)
df.grp = read.delim2("interm/risk.group.xls")

df.cb = read.delim2("interm/cibersort.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> as.data.frame()
df.cb = df.cb[df.grp$sample,]
p.cb = apply(df.cb, 2, \(x){
  wilcox.test(x~df.grp$risk.group)$p.value
})
out.cb = data.frame(cell = colnames(df.cb), 
                    p.val = p.cb, 
                    p.adj = p.adjust(p.cb, method = "BH"))

df.xc = read.delim2("interm/xcell.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> 
  as.data.frame()

imcells = c("aDC_xCell","B.cells_xCell","Basophils_xCell","CD4._memory_T.cells_xCell",
            "CD4._naive_T.cells_xCell","CD4._T.cells_xCell","CD4._Tcm_xCell","CD4._Tem_xCell",
            "CD8._naive_T.cells_xCell","CD8._T.cells_xCell","CD8._Tcm_xCell","CD8._Tem_xCell",                    
            "cDC_xCell","Class.switched_memory_B.cells_xCell","DC_xCell",
            "Eosinophils_xCell","iDC_xCell","Macrophages_xCell","Macrophages_M1_xCell",
            "Macrophages_M2_xCell","Mast_cells_xCell","Megakaryocytes_xCell","Memory_B.cells_xCell",
            "Monocytes_xCell","naive_B.cells_xCell","Neutrophils_xCell","NK_cells_xCell","NKT_xCell",
            "pDC_xCell","Plasma_cells_xCell","pro_B.cells_xCell","Tgd_cells_xCell","Th1_cells_xCell",                   
            "Th2_cells_xCell","Tregs_xCell")
df.xc = df.xc[df.grp$sample,imcells]
p.xc = apply(df.xc, 2, \(x){
  wilcox.test(x~df.grp$risk.group)$p.value
})
out.xc = data.frame(cell = colnames(df.xc), 
                    p.val = p.xc, 
                    p.adj = p.adjust(p.xc, method = "BH"))

df.ep = read.delim2("interm/epic.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> 
  as.data.frame()
df.ep = df.ep[df.grp$sample,1:6]
p.ep = apply(df.ep, 2, \(x){
  wilcox.test(x~df.grp$risk.group)$p.value
})
out.ep = data.frame(cell = colnames(df.ep), 
                    p.val = p.ep, 
                    p.adj = p.adjust(p.ep, method = "BH"))

df.qs = read.delim2("interm/quantiseq.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> 
  as.data.frame()
df.qs = df.qs[df.grp$sample,1:10]
p.qs = apply(df.qs, 2, \(x){
  wilcox.test(x~df.grp$risk.group)$p.value
})
out.qs = data.frame(cell = colnames(df.qs), 
                    p.val = p.qs, 
                    p.adj = p.adjust(p.qs, method = "BH"))

df.mc = read.delim2("interm/mcpcounter.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> 
  as.data.frame()
df.mc = df.mc[df.grp$sample,1:7]
p.mc = apply(df.mc, 2, \(x){
  wilcox.test(x~df.grp$risk.group)$p.value
})
out.mc = data.frame(cell = colnames(df.mc), 
                    p.val = p.mc, 
                    p.adj = p.adjust(p.mc, method = "BH"))
library(writexl)
write_xlsx(list(cibersort=out.cb, xcell=out.xc, epic=out.ep, 
                quantiseq=out.qs, mcpcounter=out.mc),
           "06.Phenotype/09.Immune.Diff.xls")

df.cb = df.cb[out.cb$cell[which(out.cb$p.adj<0.05)]]
df.cb$sample = rownames(df.cb)
df.cb = reshape2::melt(df.cb, id.var = "sample", variable.name = "cell")
df.cb$group = df.grp$risk.group[match(df.cb$sample, df.grp$sample)]
df.cb$p = out.cb$p.adj[match(df.cb$cell, out.cb$cell)]
df.cb$cell = str_remove(df.cb$cell, "_CIBERSORT")
df.cb$lab = paste0(df.cb$cell, "\np.adj = ", signif(df.cb$p,3))
df.cb$group = factor(df.cb$group, levels = c("Low","High"))

ggplot(df.cb, aes(x = group, y = value)) + 
  geom_violin(aes(color = group, fill = group), alpha = 0.1, show.legend = F) +
  geom_boxplot(aes(color = group), width = 0.1, alpha = 0, show.legend = F) +
  facet_wrap(~lab) + 
  scale_fill_futurama() +
  scale_color_futurama() +
  labs(title = "CIBERSORT", x = "Risk", y = "Score") +
  SCP::theme_scp() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave('06.Phenotype/10.CIBERSORT.Diff.png', width = 6, height = 5, units = 'in', dpi = 300, bg = 'white')
ggsave('06.Phenotype/10.CIBERSORT.Diff.pdf', width = 6, height = 5, units = 'in', dpi = 300, bg = 'white')

df.ep = df.ep[out.ep$cell[which(out.ep$p.adj<0.05)]]
df.ep$sample = rownames(df.ep)
df.ep = reshape2::melt(df.ep, id.var = "sample", variable.name = "cell")
df.ep$group = df.grp$risk.group[match(df.ep$sample, df.grp$sample)]
df.ep$p = out.ep$p.adj[match(df.ep$cell, out.ep$cell)]
df.ep$cell = str_remove(df.ep$cell, "_EPIC")
df.ep$lab = paste0(df.ep$cell, "\np.adj = ", signif(df.ep$p,3))
df.ep$group = factor(df.ep$group, levels = c("Low","High"))

ggplot(df.ep, aes(x = group, y = value)) + 
  geom_violin(aes(color = group, fill = group), alpha = 0.1, show.legend = F) +
  geom_boxplot(aes(color = group), width = 0.1, alpha = 0, show.legend = F) +
  facet_wrap(~lab, scale = "free", nrow = 1) + 
  scale_fill_futurama() +
  scale_color_futurama() +
  labs(title = "EPIC", x = "Risk", y = "Score") +
  SCP::theme_scp() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave('06.Phenotype/11.EPIC.Diff.png', width = 10, height = 5, units = 'in', dpi = 300, bg = 'white')
ggsave('06.Phenotype/11.EPIC.Diff.pdf', width = 10, height = 5, units = 'in', dpi = 300, bg = 'white')


rm(list = ls())
df.tide = read.delim2("../database/TIDE/TCGA.COADREAD.RNASeq.norm_subtract.OS_base") |> 
  type.convert(as.is = T)
df.grp = read.delim2("interm/risk.group.xls")
df.grp$sample = str_sub(df.grp$sample, 1, 12)
df.grp = df.grp[!duplicated(df.grp$sample),]
df.tide$sample = rownames(df.tide) |> make.names()
df.full = merge(df.grp, df.tide)
df.full = type.convert(df.full, as.is = T)
df.full$group = factor(df.full$risk.group, levels = c("Low","High"))

p1 = ggplot(df.full, aes(x = group, y = Dysfunction)) + 
  geom_violin(aes(color = group, fill = group), alpha = 0.1, show.legend = F) +
  geom_boxplot(aes(color = group), width = 0.1, alpha = 0, show.legend = F) +
  scale_fill_futurama() +
  scale_color_futurama() +
  stat_compare_means(label.x = 1.5, hjust = 0.5) +
  SCP::theme_scp() + 
  labs(x = "Risk") +
  theme(plot.title = element_text(hjust = 0.5))
p2 = ggplot(df.full, aes(x = group, y = Exclusion)) + 
  geom_violin(aes(color = group, fill = group), alpha = 0.1, show.legend = F) +
  geom_boxplot(aes(color = group), width = 0.1, alpha = 0, show.legend = F) +
  scale_fill_futurama() +
  scale_color_futurama() +
  stat_compare_means(label.x = 1.5, hjust = 0.5) +
  SCP::theme_scp() + 
  labs(x = "Risk") +
  theme(plot.title = element_text(hjust = 0.5))
p1+p2
ggsave('06.Phenotype/12.TIDE.png', width = 6, height = 4, units = 'in', dpi = 300, bg = 'white')
ggsave('06.Phenotype/12.TIDE.pdf', width = 6, height = 4, units = 'in', dpi = 300, bg = 'white')


rm(list = ls())
df.exp = read.delim2("cleandata/fpkm.tsv", row.names = 1) |> 
  type.convert(as.is = T)
gs = c("ADORA2A","BTLA","BTNL2","CD160",
  "CD200","CD200R1","CD244","CD27",
  "CD274","CD276","CD28","CD40",
  "CD40LG","CD44","CD48","CD70",
  "CD80","CD86","CTLA4","HAVCR2",
  "HHLA2","ICOS","ICOSLG","IDO1",
  "IDO2","KIR3DL1","LAG3","LAIR1",
  "LGALS9","NRP1","PDCD1","PDCD1LG2",
  "TIGIT","TMIGD2","TNFRSF9","TNFRSF14",
  "TNFRSF18","TNFRSF25","TNFRSF4",
  "TNFRSF8","TNFSF4","TNFSF9","TNFSF14",
  "TNFSF15","TNFSF18","VTCN1")
df.exp = df.exp[gs,]
df.grp = read.delim2("interm/risk.group.xls")
df.exp = df.exp[df.grp$sample]
df.exp = apply(df.exp, c(1,2), \(x) log2(x+1)) |> as.data.frame()

p = apply(df.exp, 1, \(x){wilcox.test(x~df.grp$risk.group)$p.value})
df.out = data.frame(gene = rownames(df.exp),
                    p.val = p,
                    p.adj = p.adjust(p, method = "BH"))
df.out = df.out[order(df.out$p.adj, decreasing = F),]
write.table(df.out, "06.Phenotype/13.Checkpoint.xls", quote = F, sep = "\t", row.names = F)

df.out = df.out[df.out$p.adj < 0.05, ]
df.exp = df.exp[df.out$gene,]
df.exp$gene = rownames(df.exp)
df.exp = reshape2::melt(df.exp, id.var = "gene", variable.name = "sample")
df.exp$group = df.grp$risk.group[match(df.exp$sample, df.grp$sample)]
df.exp$group = factor(df.exp$group, levels = c("Low","High"))
df.exp$p = df.out$p.adj[match(df.exp$gene, df.out$gene)]
df.exp$lab = paste0(df.exp$gene, "\np.adj = ", signif(df.exp$p, 3))

ggplot(df.exp, aes(x = group, y = value)) +
  geom_boxplot(aes(fill = group), show.legend = F) + 
  facet_wrap(~lab, scale = "free", nrow = 2) +
  scale_fill_futurama(name = "Risk") + 
  labs(x = "Risk", y = "Expression") +
  theme_bw() + 
  theme(legend.position = "top", 
        axis.text.x.bottom = element_text(angle = 45, hjust = 1))
ggsave("06.Phenotype/14.Checkpoint.png", width = 10, height = 6, units = "in", dpi = 300, bg = "white")
ggsave("06.Phenotype/14.Checkpoint.pdf", width = 10, height = 6, units = "in", dpi = 300, bg = "white")


df.exp = read.delim2("cleandata/fpkm.tsv", row.names = 1) |> 
  type.convert(as.is = T)
gs = read.delim2("03.Risk/05.Coefficients.xls")
df.im = df.exp[df.out$gene,] |> t() |> as.data.frame()
df.hub = df.exp[gs$gene,] |> t() |> as.data.frame()
res.cor = linkET::correlate(df.im, df.hub, method = "spearman")
df.r = res.cor$r |> as.data.frame()
df.p = res.cor$p |> as.data.frame()
df.r$gene = rownames(df.r)
df.p$gene = rownames(df.p)
df.r = reshape2::melt(df.r, id.var = "gene", variable.name = "hub", value.name = "r")
df.p = reshape2::melt(df.p, id.var = "gene", variable.name = "hub", value.name = "p")
df.cor = merge(df.r, df.p)
write.table(df.cor, "06.Phenotype/15.Checkpoint.Cor.xls", quote = F, sep = "\t", row.names = F)


df.cor = df.cor[df.cor$p < 0.05,]
ggplot(df.cor, aes(x = gene, y = hub)) + 
  geom_point(aes(size = p, color = r)) + 
  scale_size_continuous(range = c(7,4), transform = "log10") + 
  scale_color_gradient2(low = "darkblue", high = "darkred") + 
  labs(y = "Hub genes", x = "Immune checkpoints") +
  SCP::theme_scp() + 
  theme(axis.text.x.bottom = element_text(hjust = 1, angle = 30))
ggsave('06.Phenotype/16.Checkpoint.Cor.png', width = 7, height = 8, units = 'in', dpi = 300, bg = 'white')
ggsave('06.Phenotype/16.Checkpoint.Cor.pdf', width = 7, height = 8, units = 'in', dpi = 300, bg = 'white')


rm(list = ls())
library(oncoPredict)

df.exp = read.delim2("cleandata/fpkm.tsv", row.names = 1) |> 
  apply(c(1,2), \(x) log2(as.numeric(x)+1)) |> as.data.frame()

train.Expr = readRDS(file="../database/GDSC2_Expr (RMA Normalized and Log Transformed).rds")
train.Res = readRDS(file = "../database/GDSC2_Res.rds")

calcPhenotype(trainingExprData = train.Expr,
              trainingPtype = train.Res,
              testExprData = as.matrix(df.exp), 
              batchCorrect = 'eb',
              powerTransformPhenotype = F,
              minNumSamples = 20,
              printOutput = T,
              removeLowVaryingGenes = 0.4,
              removeLowVaringGenesFrom = "homogenizeData")

df = read.csv("calcPhenotype_Output/DrugPredictions.csv")
colnames(df)[1] = "sample"
write.table(df, "interm/drug.pred.xls", quote = F, sep = "\t", row.names = F)



rm(list = ls())
library(cowplot)
df.im = read.delim2("interm/drug.pred.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> 
  as.data.frame()
df.grp = read.delim2("interm/risk.group.xls")

df.im = df.im[df.grp$sample,]
df.diff = data.frame(
  drug = colnames(df.im),
  p.val = apply(df.im, 2, function(x){
    r = wilcox.test(x~df.grp$risk.group)
    r$p.value
  }))

df.diff$p.adj = p.adjust(df.diff$p.val, method = "BH")
df.diff = df.diff[order(df.diff$p.adj, decreasing = F),]
df.diff$drug = str_replace(df.diff$drug, "X5.Fluorouracil_1073", "5.Fluorouracil_1073")
write.table(df.diff, "06.Phenotype/17.Drug.Diff.xls", quote = F, sep = "\t", row.names = F)

df.diff = df.diff[1:8,]
df.diff$lab = paste0("p.adj = ", signif(df.diff$p.adj, 3))

df.im = df.im[df.diff$drug]
df.im$sample = rownames(df.im)
df.im = reshape2::melt(df.im, id.var = "sample", variable.name = "drug")
df.im$risk = df.grp$risk[match(df.im$sample, df.grp$sample)]
df.im$group = df.grp$risk.group[match(df.im$sample, df.grp$sample)]
df.im = type.convert(df.im, as.is = T)
df.im$group = factor(df.im$group, levels = c("Low","High"))

plot.list = lapply(unique(df.im$drug), \(x){
  d = df.im[df.im$drug == x,]
  
  lab = df.diff[df.diff$drug == x, "lab"]
  ymax = ifelse(max(d$value) > 0, max(d$value) * 1.1, max(d$value) * 0.9)
  p1 = ggplot(d, aes(x = group, y = value)) +
    geom_violin(aes(fill = group), show.legend = F, alpha = 0.5, linewidth = 1) +
    geom_boxplot(aes(fill = group), show.legend = F, width = 0.3) +
    annotate("text", x = 1.5, y = ymax, label = lab, size = 5, hjust = 0.5) + 
    scale_fill_futurama() + 
    theme_bw() + 
    labs(x = "Group", y = "Predicted IC50") + 
    theme(aspect.ratio = 1)
  
  res.cor = cor.test(d$value, d$risk, method = "spearman")
  if(res.cor$p.value == 0) {
    res.cor$p.value = "p < 1e-16"
  } else {
    res.cor$p.value = paste0("p = ", signif(res.cor$p.value, 3))
  }
  lab = paste0("r = ", signif(res.cor$estimate,2), ", ", res.cor$p.value)
  
  p2 = ggplot(d, aes(x = risk, y = value)) + 
    geom_point(aes(color = group), show.legend = F) + 
    scale_colour_futurama() + 
    annotate("text", x = 0, y = max(d$value) * 1.05, label = lab, size = 5, hjust = 0) +
    theme_bw() + 
    labs(x = "Risk", y = "Predicted IC50") + 
    theme(aspect.ratio = 1)
  p = p1 + p2
  main_title = ggdraw() + 
    draw_label(str_replace_all(x, fixed("."), "-"), fontface = 'bold', size = 16, x = 0.5, y = 0.5)
  plot_grid(main_title, p, ncol = 1, rel_heights = c(0.1,1))
})


p = plot_grid(plotlist = plot.list, ncol = 2)
ggsave("06.Phenotype/18.Drug.Diff.png", p, width = 13, height = 14, units = "in", dpi = 300, bg = "white")
ggsave("06.Phenotype/18.Drug.Diff.pdf", p, width = 13, height = 14, units = "in", dpi = 300, bg = "white")






