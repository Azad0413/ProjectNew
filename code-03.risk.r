rm(list = ls())
library(glmnet)
library(survival)
library(survminer)
library(tidyverse)
library(ggsci)
library(patchwork)
library(survivalROC)

df = read.delim2("cleandata/fpkm.tsv", row.names = 1) |> 
  apply(c(1,2), \(x) log2(as.numeric(x)+1)) |> as.data.frame()
df.surv = read.delim2("cleandata/survival.tsv")
df.surv$OS.time = as.numeric(df.surv$OS.time)
df.surv$OS = as.numeric(df.surv$OS)
df = df[df.surv$sample]

gs = read.delim2("02.DEG/06.Intersection.xls")
gs = gs$gene
df.tmp = read.delim2("cleandata/GSE39582.exp.xls")
gs = gs[gs %in% df.tmp$gene]
df.exp = df[gs,]
df.exp = na.omit(df.exp)

res.univar = apply(df.exp, 1, function(x){
  res.sub = coxph(Surv(time = df.surv$OS.time, event = df.surv$OS) ~ x) |> 
    summary()
  p.value = signif(res.sub$wald["pvalue"], digits=4)
  hr = signif(res.sub$coef[2], digits=3)
  hr.low = signif(res.sub$conf.int[,"lower .95"],digits = 3)
  hr.high = signif(res.sub$conf.int[,"upper .95"],digits = 3)
  return(c(hr, p.value, hr.low, hr.high))
}) |> t() |> as.data.frame()

colnames(res.univar) = c("HR","p.value","CI.low","CI.high")
res.univar$Gene = rownames(res.univar)
res.univar = res.univar[c(5,1,3:4,2)]
res.univar.plot = subset(res.univar, p.value < 0.05)
res.univar.plot = res.univar.plot[order(res.univar.plot$p.value, decreasing = F),]
write.table(res.univar.plot, "03.Risk/01.UnivarCox.xls", sep = "\t", col.names = T, row.names = F, quote = F)

res.univar.plot = res.univar.plot[order(res.univar.plot$HR, decreasing = T),]
res.univar.plot$Gene = factor(res.univar.plot$Gene, levels = res.univar.plot$Gene)
res.univar.plot$plabel = res.univar.plot$p.value %>% signif(3)
res.univar.plot$hlabel = paste0(res.univar.plot$HR,"(", res.univar.plot$CI.low, "-", res.univar.plot$CI.high, ")")

library(aplot)
p = ggplot(data = res.univar.plot) + 
  geom_vline(xintercept = 1, color = "grey60", linetype = 1) + 
  geom_segment(aes(x = CI.low, xend = CI.high, y = Gene, yend = Gene), color = "darkblue", lwd = 1) +
  geom_point(aes(y = Gene, x = HR), size = 3, shape = 15, color = "red") +
  scale_x_continuous(trans = "log10") +
  theme_minimal() + 
  ylab(NULL) + xlab("Harzard Ratio") + 
  scale_color_manual(values = c("green","red"), guide = guide_none()) +
  theme(axis.line.x = element_line(linewidth = 1), 
        axis.title = element_text(size = 14), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_blank(),
        panel.grid = element_blank(), 
        plot.background = element_rect(fill = "white", colour = "transparent"))
dt = ggplot(data = res.univar.plot) + 
  geom_text(aes(y = Gene, x = 1, label = plabel), hjust = 0.5) +
  theme_minimal() + ggtitle("p.value") + xlab(NULL) + ylab(NULL) +
  theme(axis.text = element_blank(), 
        plot.background = element_rect(fill = "white", colour = "transparent"),
        plot.title = element_text(hjust = 0.5), 
        line = element_blank())
dt2 = ggplot(data = res.univar.plot) + 
  geom_text(aes(y = Gene, x = 1, label = hlabel), hjust = 0.5) +
  theme_minimal() + ggtitle("Harzad Ratio") + xlab(NULL) + ylab(NULL) +
  theme(axis.text = element_blank(), 
        plot.background = element_rect(fill = "white", colour = "transparent"),
        plot.title = element_text(hjust = 0.5), 
        line = element_blank())
dt3 = ggplot(data = res.univar.plot) + 
  geom_text(aes(y = Gene, x = 1, label = Gene), hjust = 0) +
  theme_minimal() + xlab(NULL) + ylab(NULL) +
  theme(axis.text= element_blank(), 
        plot.background = element_rect(fill = "white", colour = "transparent"),
        plot.title = element_text(hjust = 0.5), 
        line = element_blank())

p.final = insert_left(p, dt2, width = 0.35) %>% insert_left(dt, width = 0.2) %>% insert_left(dt3, width = 0.4)
p.final
ggsave("03.Risk/02.UnivarCox.png", p.final, width = 8, height = 12, units = "in", dpi = 300, bg = "white")
ggsave("03.Risk/02.UnivarCox.pdf", p.final, width = 8, height = 12, units = "in", dpi = 300, bg = "white")

df.exp = df.exp[as.character(gs),]

set.seed(2025)
res.lasso = cv.glmnet(t(df.exp), Surv(time = df.surv$OS.time, event = df.surv$OS), family = "cox",
                      nfolds = 10, alpha = 1, type.measure = "C")
plot(res.lasso)

ggsave("03.Risk/03.Lasso.CV.png", plot(res.lasso), width = 8, height = 7, dpi = 300, units = "in", bg = "white")
ggsave("03.Risk/03.Lasso.CV.pdf", plot(res.lasso), width = 8, height = 7, dpi = 300, units = "in", bg = "white")
ggsave("03.Risk/04.Lasso.Coef.png", plot(res.lasso$glmnet.fit, xvar = 'lambda'), width = 8, height = 7, dpi = 300, units = "in", bg = "white")
ggsave("03.Risk/04.Lasso.Coef.pdf", plot(res.lasso$glmnet.fit, xvar = 'lambda'), width = 8, height = 7, dpi = 300, units = "in", bg = "white")

df.coef = coef(res.lasso, s = res.lasso$lambda.min)
df.coef = cbind(gene = rownames(df.coef), coefficient = df.coef[,1]) |> as.data.frame()
df.coef = subset(df.coef, coefficient != 0)
write.table(df.coef, "03.Risk/05.Coefficients.xls", sep = "\t", quote = F, col.names = T, row.names = F)

df.surv$risk = predict(res.lasso, s = res.lasso$lambda.min, newx = t(df.exp), type = "link")[,1]
df.surv$risk.group = ifelse(df.surv$risk > median(df.surv$risk), "High", "Low")
df.surv$risk.group = factor(df.surv$risk.group, levels = c("Low","High"))
df.surv = df.surv[order(df.surv$risk, decreasing = F),]
df.surv$rn = 1:nrow(df.surv)
df.surv$status = ifelse(df.surv$OS == 0, "Alive", "Dead")
df.exp = df.exp[df.surv$sample]

df.exp = df.exp[df.coef$gene,]
df.hm = apply(df.exp, 1, scale) |> t() |> as.data.frame()
colnames(df.hm) = colnames(df.exp)
df.hm = cbind(gene = rownames(df.hm), df.hm)
df.hm = reshape2::melt(df.hm, id.var = "gene", variable.name = "sample")

p1 = ggplot(df.surv, aes(x = rn, y = risk, color = risk.group)) + 
  geom_hline(yintercept = median(df.surv$risk), linetype = 2) +
  geom_vline(xintercept = sum(df.surv$risk.group=="Low"), linetype = 2) +
  geom_point() +
  scale_color_aaas(name = "Risk") +
  ggtitle("TCGA") +
  ylab("Risk Score") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
p2 = ggplot(df.surv, aes(x = rn, y = OS.time, color = status)) + 
  geom_vline(xintercept = sum(df.surv$risk.group=="Low"), linetype = 2) +
  geom_point() +
  scale_color_manual(values = c("#82B29A","#F2CC8E"), name = "Status") +
  ylab("Survival Time (Days)") +
  theme_minimal() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())
p3 = ggplot(df.hm, aes(x = sample, y = gene)) + 
  geom_tile(aes(fill = value)) + 
  ylab("Gene") + xlab("Sample with increasing risk score") +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_gradient2(low = "darkblue", high = "darkred", name = "Expression") +
  theme_minimal() + 
  theme(axis.text.x = element_blank())
p = p1/p2/p3 + plot_layout(heights = c(1,1,2))
ggsave("03.Risk/06.TCGA.Dist.png", p, width = 9, height = 9, units = "in", dpi = 300, bg = "white")
ggsave("03.Risk/06.TCGA.Dist.pdf", p, width = 9, height = 9, units = "in", dpi = 300, bg = "white")

surv.fit = survfit(Surv(time = OS.time,event = OS) ~ risk.group, data = df.surv)
p.surv = ggsurvplot(surv.fit, data = df.surv, conf.int = F, 
                    risk.table = TRUE, break.x.by=800,
                    pval = T, legend.labs = c("Low", "High"), 
                    palette = pal_aaas()(2)) + 
  xlab("Time (Day)")
p.surv$plot = p.surv$plot + 
  guides(color = guide_legend(title = "Risk")) + 
  ggtitle("TCGA") +
  theme(plot.title = element_text(hjust = 0.5))
p.surv$table = p.surv$table + ylab("Risk")
p.surv
ggsave("03.Risk/07.TCGA.KM.png", print(p.surv), width = 7, height = 7, units = "in", dpi = 300, bg = "white")
ggsave("03.Risk/07.TCGA.KM.pdf", print(p.surv, newpage = F), width = 7, height = 7, units = "in", dpi = 300, bg = "white")
write.table(df.surv, "interm/risk.group.xls", quote = F, row.names = F, sep = "\t")



rm(list = ls())
df.exp = read.delim2("cleandata/GSE39582.exp.xls", row.names = 1) |> 
  type.convert(as.is = T)
df.surv = read.delim2("cleandata/GSE39582.pheno.xls")

df.coef = read.delim2("03.Risk/05.Coefficients.xls")
coef = as.numeric(df.coef$coefficient)
df.exp = df.exp[df.coef$gene,df.surv$geo_accession]
df.surv$risk = apply(df.exp, 2, function(x){x %*% coef})

df.surv$OS = as.numeric(df.surv$os.event.ch1)
df.surv$OS.time = as.numeric(df.surv$os.delay..months..ch1)
df.surv = df.surv[!is.na(df.surv$OS.time),]
df.surv$risk.group = ifelse(df.surv$risk >= median(df.surv$risk), "High", "Low")
df.surv$risk.group = factor(df.surv$risk.group, levels = c("Low","High"))

df.surv = df.surv[order(df.surv$risk, decreasing = F),]
df.surv$rn = 1:nrow(df.surv)
df.surv$status = ifelse(df.surv$OS == 0, "Alive", "Dead")
df.exp = df.exp[df.surv$geo_accession]

df.hm = apply(df.exp, 1, scale) |> t() |> as.data.frame()
colnames(df.hm) = colnames(df.exp)
df.hm = cbind(gene = rownames(df.hm), df.hm)
df.hm = reshape2::melt(df.hm, id.var = "gene", variable.name = "sample")

p1 = ggplot(df.surv, aes(x = rn, y = risk, color = risk.group)) + 
  geom_hline(yintercept = median(df.surv$risk), linetype = 2) +
  geom_vline(xintercept = sum(df.surv$risk.group=="Low"), linetype = 2) +
  geom_point() +
  scale_color_aaas(name = "Risk") +
  ggtitle("GSE39582") +
  ylab("Risk Score") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
p2 = ggplot(df.surv, aes(x = rn, y = OS.time, color = status)) + 
  geom_vline(xintercept = sum(df.surv$risk.group=="Low"), linetype = 2) +
  geom_point() +
  scale_color_manual(values = c("#82B29A","#F2CC8E"), name = "Status") +
  ylab("Survival Time (Months)") +
  theme_minimal() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())
p3 = ggplot(df.hm, aes(x = sample, y = gene)) + 
  geom_tile(aes(fill = value)) + 
  ylab("Gene") + xlab("Sample with increasing risk score") +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_gradient2(low = "darkblue", high = "darkred", name = "Expression") +
  theme_minimal() + 
  theme(axis.text.x = element_blank())
p = p1/p2/p3 + plot_layout(heights = c(1,1,2))
ggsave("03.Risk/08.GSE39582.Dist.png", p, width = 9, height = 9, units = "in", dpi = 300, bg = "white")
ggsave("03.Risk/08.GSE39582.Dist.pdf", p, width = 9, height = 9, units = "in", dpi = 300, bg = "white")

surv.fit = survfit(Surv(time = OS.time,event = OS) ~ risk.group, data = df.surv)
p.surv = ggsurvplot(surv.fit, data = df.surv, conf.int = F, risk.table = TRUE, 
                    break.x.by=24, palette = c("#3B4992FF","#EE0000FF"),
                    pval = T, legend.labs = c("Low", "High")) + 
  xlab("Time (Months)")
p.surv$plot = p.surv$plot + 
  guides(color = guide_legend(title = "Risk")) + 
  ggtitle("GSE39582") +
  theme(plot.title = element_text(hjust = 0.5))
p.surv$table = p.surv$table + ylab("Risk")
p.surv
ggsave("03.Risk/09.GSE39582.KM.png", print(p.surv), width = 7, height = 7, units = "in", dpi = 300, bg = "white")
ggsave("03.Risk/09.GSE39582.KM.pdf", print(p.surv, newpage = F), width = 7, height = 7, units = "in", dpi = 300, bg = "white")
write.table(df.surv, "interm/GSE39582.risk.xls", quote = F, row.names = F, sep = "\t")



rm(list = ls())
df.exp = read.delim2('cleandata/fpkm.tsv', row.names = 1) |>
  type.convert(as.is = T)
gs = read.delim2("03.Risk/05.Coefficients.xls")
df.exp = df.exp[gs$gene,]

df.exp = cbind(gene = rownames(df.exp), df.exp)
df.exp = reshape2::melt(df.exp, id.var = "gene", variable.name = "sample")
df.exp$value = log2(df.exp$value + 1)
df.exp$group = ifelse(str_ends(df.exp$sample, "1.."), "Control", "Tumor")
df.diff = read.delim2("02.DEG/01.TCGA.DEG.xls") |> type.convert(as.is = T)
df.exp$p = df.diff$adj.P.Val[match(df.exp$gene, df.diff$Gene)]
df.exp$lab = paste0(df.exp$gene, "\np.adj = ", signif(df.exp$p,3))

ggplot(df.exp, aes(x = group, y = value)) + 
  geom_boxplot(aes(color = group), show.legend = F) + 
  scale_color_aaas() + 
  facet_wrap(~lab, nrow = 4, scale = "free") + 
  theme_bw() + 
  labs(x = "Group", y = "Expression")
ggsave('03.Risk/10.Model.Gene.png', width = 10, height = 10, units = 'in', dpi = 300, bg = 'white')
ggsave('03.Risk/10.Model.Gene.pdf', width = 10, height = 10, units = 'in', dpi = 300, bg = 'white')




