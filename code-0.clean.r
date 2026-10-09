rm(list = ls())
library(org.Hs.eg.db)
library(tidyverse)
library(GEOquery)

df.count1 = read.delim2("../database/TCGA/TCGA-COAD.star_counts.tsv.gz")
df.fpkm1 = read.delim2("../database/TCGA/TCGA-COAD.star_fpkm.tsv.gz")
df.pheno1 = read.delim2("../database/TCGA/TCGA-COAD.clinical.tsv.gz")
df.survival1 = read.delim2("../database/TCGA/TCGA-COAD.survival.tsv.gz")

df.count2 = read.delim2("../database/TCGA/TCGA-READ.star_counts.tsv.gz")
df.fpkm2 = read.delim2("../database/TCGA/TCGA-READ.star_fpkm.tsv.gz")
df.pheno2 = read.delim2("../database/TCGA/TCGA-READ.clinical.tsv.gz")
df.survival2 = read.delim2("../database/TCGA/TCGA-READ.survival.tsv.gz")

df.count = merge(df.count1, df.count2, by = "Ensembl_ID")
df.fpkm = merge(df.fpkm1, df.fpkm2, by = "Ensembl_ID")

rownames(df.count) = df.count[[1]]; df.count = df.count[-1]
rownames(df.fpkm) = df.fpkm[[1]]; df.fpkm = df.fpkm[-1]

identical(colnames(df.pheno1), colnames(df.pheno2))
df.pheno1 = df.pheno1[colnames(df.pheno2)]
df.pheno = rbind(df.pheno1, df.pheno2)
df.survival = rbind(df.survival1, df.survival2)

df.pheno$sample = make.names(df.pheno$sample)
df.survival$sample = make.names(df.survival$sample)
smps = intersect(df.pheno$sample, df.survival$sample)
smps = smps[!str_ends(smps, "1..")]
smps = smps[smps %in% colnames(df.count)]
df.pheno = subset(df.pheno, sample %in% smps)
df.survival = subset(df.survival, sample %in% smps)

df.count$geneid = str_remove(rownames(df.count), "\\..*")
df.fpkm$geneid = str_remove(rownames(df.fpkm), "\\..*")
df.count = subset(df.count, !duplicated(geneid))
df.fpkm = subset(df.fpkm, !duplicated(geneid))

df.count$gene = mapIds(org.Hs.eg.db, df.count$geneid, "SYMBOL", "ENSEMBL", multiVals = "first")
df.fpkm$gene = mapIds(org.Hs.eg.db, df.fpkm$geneid, "SYMBOL", "ENSEMBL", multiVals = "first")

df.count = na.omit(df.count)
df.fpkm = na.omit(df.fpkm)
df.count = subset(df.count, !duplicated(gene))
df.fpkm = subset(df.fpkm, !duplicated(gene))
rownames(df.count) = df.count$gene
rownames(df.fpkm) = df.fpkm$gene
df.count = df.count[-which(colnames(df.count) %in% c("gene","geneid"))]
df.fpkm = df.fpkm[-which(colnames(df.fpkm) %in% c("gene","geneid"))]

df.count = apply(df.count, c(1,2), function(x){
  x = 2^(as.numeric(x))-1 |> as.integer()
}) |> as.data.frame()

df.fpkm = apply(df.fpkm, c(1,2), function(x){
  x = 2^(as.numeric(x))-1
}) |> as.data.frame()

df.count = cbind(gene = rownames(df.count), df.count)
df.fpkm = cbind(gene = rownames(df.fpkm), df.fpkm)

str_ends(colnames(df.fpkm)[-1], "1..") |> table()
# 640 CAC 51 NC

write.table(df.count, "cleandata/count.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(df.fpkm, "cleandata/fpkm.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(df.survival, "cleandata/survival.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(df.pheno, "cleandata/pheno.tsv", sep = "\t", row.names = F, col.names = T, quote = F)


rm(list = ls())
gs = getGEO(filename = "rawdata/GSE41258_series_matrix.txt.gz", getGPL = F)
df.exp = gs |>  exprs() |>  as.data.frame()
df.pheno = gs |>  pData() |>  as.data.frame()
df.annot = getGEO("GPL96", destdir = "rawdata", AnnotGPL = T) |> Table()
df.annot = subset(df.annot, !str_detect(`Gene symbol`, fixed("/")))
df.exp$gene = df.annot$`Gene symbol`[match(rownames(df.exp), df.annot$ID)]
df.exp = subset(df.exp, gene != "")
df.exp = na.omit(df.exp)
df.exp = aggregate(.~gene, FUN = max, data = df.exp)
df.pheno = df.pheno[df.pheno$`tissue:ch1` == "Primary Tumor",]
df.pheno = df.pheno[!is.na(df.pheno$title),]
df.exp = df.exp[c("gene", df.pheno$geo_accession)]
write.table(df.exp, "cleandata/GSE41258.exp.xls", sep = "\t", quote = F, row.names = F)
write.table(df.pheno, "cleandata/GSE41258.pheno.xls", sep = "\t", quote = F, row.names = F)


rm(list = ls())
gs = getGEO(filename = "rawdata/GSE39582_series_matrix.txt.gz", getGPL = F)
df.exp = gs |>  exprs() |>  as.data.frame()
df.pheno = gs |>  pData() |>  as.data.frame()
df.annot = getGEO("GPL570", destdir = "rawdata", AnnotGPL = T) |> Table()
df.annot = subset(df.annot, !str_detect(`Gene symbol`, fixed("/")))
df.exp$gene = df.annot$`Gene symbol`[match(rownames(df.exp), df.annot$ID)]
df.exp = subset(df.exp, gene != "")
df.exp = na.omit(df.exp)
df.exp = aggregate(.~gene, FUN = max, data = df.exp)
write.table(df.exp, "cleandata/GSE39582.exp.xls", sep = "\t", quote = F, row.names = F)
write.table(df.pheno, "cleandata/GSE39582.pheno.xls", sep = "\t", quote = F, row.names = F)



rm(list = ls())
gs = getGEO(filename = "rawdata/GSE17537_series_matrix.txt.gz", getGPL = F)
df.exp = gs |> exprs() |> as.data.frame()
df.pheno = gs |> pData() |> as.data.frame()
df.annot = getGEO(filename = "rawdata/GPL570.annot.gz", getGPL = F) |> Table()
df.annot = subset(df.annot, !str_detect(`Gene symbol`, "//"))
df.exp$gene = df.annot$`Gene symbol`[match(rownames(df.exp), df.annot$ID)]
df.exp = subset(df.exp, gene != "")
df.exp = na.omit(df.exp)
df.exp = aggregate(.~gene, FUN = max, data = df.exp)
write.table(df.exp, "cleandata/GSE17537.exp.xls", sep = "\t", quote = F, row.names = F)
write.table(df.pheno, "cleandata/GSE17537.pheno.xls", sep = "\t", quote = F, row.names = F)


