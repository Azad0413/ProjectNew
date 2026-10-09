rm(list = ls())
library(tidyverse)
library(org.Hs.eg.db)
library(survival)
library(rms)
library(ggpubr)
library(ggsci)
library(forestploter)
library(grid)
library(timeROC)

df.surv = read.delim2("interm/risk.group.xls")
df.pheno = read.delim2("cleandata/pheno.tsv")
df.pheno = df.pheno[c("sample","age_at_index.demographic",
                      "ajcc_pathologic_t.diagnoses",
                      "ajcc_pathologic_n.diagnoses",
                      "ajcc_pathologic_m.diagnoses",
                      "gender.demographic")]

colnames(df.pheno) = c("sample","Age","T_Stage","N_Stage","M_Stage","Gender")
df.full = merge(df.surv, df.pheno, by = "sample")

df.full$risk = as.numeric(df.full$risk)
res.risk = coxph(Surv(time = OS.time, event = OS) ~ risk, data = df.full) %>% summary
res.risk = c(res.risk$conf.int[-2], res.risk$coefficients[5])

df.full$Age = as.numeric(df.full$Age)
res.age = coxph(Surv(time = OS.time, event = OS) ~ Age, data = df.full) %>% summary
res.age = c(res.age$conf.int[-2], res.age$coefficients[5])

res.gender = coxph(Surv(time = OS.time, event = OS) ~ Gender, data = df.full) %>% summary
res.gender = c(res.gender$conf.int[-2], res.gender$coefficients[5])

df.full$T_Stage = str_match(df.full$T_Stage, "T[1-4]")[,1]
res.t = coxph(Surv(time = OS.time, event = OS) ~ T_Stage, data = df.full) %>% summary
res.t = cbind(res.t$conf.int[,-2], res.t$coefficients[,5])

df.full$N_Stage = str_match(df.full$N_Stage, "N[0-2]")[,1]
res.n = coxph(Surv(time = OS.time, event = OS) ~ N_Stage, data = df.full) %>% summary
res.n = cbind(res.n$conf.int[,-2], res.n$coefficients[,5])

df.full$M_Stage = str_match(df.full$M_Stage, "M[0-1]")[,1]
res.m = coxph(Surv(time = OS.time, event = OS) ~ M_Stage, data = df.full) %>% summary
res.m = c(res.m$conf.int[-2], res.m$coefficients[5])


res.na = rep(NA,4)
res = rbind(res.risk,res.age,
            res.na,res.na,res.gender,
            res.na,res.na,res.t,
            res.na,res.na,res.n,
            res.na,res.na,res.m) %>% as.data.frame()
colnames(res) = c("hr","low","up","pv")

res$Characteristic = c("Risk", "Age", 
                       "Gender", "  Female", "  Male",
                       "T Stage", "  T1", "  T2", "  T3", "  T4", 
                       "N Stage", "  N0", "  N1", "  N2", 
                       "M Stage", "  M0", "  M1"
                       )
res$HR = paste0(round(res$hr,2), " (", round(res$low,2), "-", round(res$up,2), ")")
res$HR[res$Characteristic %in% c("Gender","T Stage","N Stage","M Stage")] = " "
res$HR[res$Characteristic %in% c("  Female","  T1","  N0","  M0")] = "1 (ref)"
res$p = signif(res$pv,3)
res$p[is.na(res$p)] = " "
res$blank = paste(rep(" ",20), collapse = " ")
colnames(res)[ncol(res)] = " "
colnames(res)[which(colnames(res)=="HR")]="HR (95% CI)"

tm = forest_theme(
  base_size = 10,
  ci_pch = 18,
  ci_col = "black",
  ci_lty = 1,
  ci_lwd = 1.8,
  ci_Theight = 0.2,
  refline_lwd = 1,
  refline_lty = 2,
  refline_col = "grey20"
)


p = forest(res[c(5,6,8,7)], 
       est = res$hr,
       lower = res$low,
       upper = res$up,
       ci_column = 3,
       ref_line = 1,
       size = 1,
       xlim = c(-1,9),
       theme = tm,
       ticks_at = c(1,3,6,9))
p = edit_plot(p, part = "header", gp = gpar(fontface = "bold"))
p = add_border(p, part = "header", where = "bottom")
p = add_border(p, part = "header", where = "top")
p = edit_plot(p, which = "background", gp = gpar(fill = "white"))
sigrows = which(res$pv < 0.05)
p = edit_plot(p, row = sigrows, col = 2, gp = gpar(fontface = "bold"))
p = edit_plot(p, row = sigrows, col = 3, which = "ci", gp = gpar(col = "green3", lwd = 2))
p

png("05.Prognostic/01.Univar.Cox.png", width = 7, height = 7, units = "in", bg = "white", res = 300)
plot(p)
dev.off()
pdf("05.Prognostic/01.Univar.Cox.pdf", width = 7, height = 7)
plot(p)
dev.off()



res.mul = coxph(Surv(time = OS.time, event = OS) ~ risk + Age + T_Stage + N_Stage + M_Stage, data = df.full) %>% summary
res.mul = cbind(res.mul$conf.int[,-2], res.mul$coefficients[,5])

res.na = rep(NA,4)
res.mul = rbind(res.mul[1:2,],
                res.na,res.na,res.mul[3:5,],
                res.na,res.na,res.mul[6:7,],
                res.na,res.na,res.mul[8,]) %>% as.data.frame()
colnames(res.mul) = c("hr","low","up","pv")
res.mul$Characteristic = c("Risk", "Age",
                       "T Stage", "  T1", "  T2", "  T3", "  T4", 
                       "N Stage", "  N0", "  N1", "  N2", 
                       "M Stage", "  M0", "  M1"
)
res.mul$HR = paste0(round(res.mul$hr,2), " (", round(res.mul$low,2), "-", round(res.mul$up,2), ")")
res.mul$HR[res.mul$Characteristic %in% c("Gender","Age","T Stage","N Stage","M Stage")] = " "
res.mul$HR[res.mul$Characteristic %in% c("  Female","  T1","  N0","  M0")] = "1 (ref)"
res.mul$p = signif(res.mul$pv,3)
res.mul$p[is.na(res.mul$p)] = " "
res.mul$blank = paste(rep(" ",20), collapse = " ")
colnames(res.mul)[ncol(res.mul)] = " "
colnames(res.mul)[which(colnames(res.mul)=="HR")]="HR (95% CI)"

p = forest(res.mul[c(5,6,8,7)], 
           est = res.mul$hr,
           lower = res.mul$low,
           upper = res.mul$up,
           ci_column = 3,
           ref_line = 1,
           size = 1,
           xlim = c(-1,9),
           theme = tm,
           ticks_at = c(1,3,6,9))
p = edit_plot(p, part = "header", gp = gpar(fontface = "bold"))
p = add_border(p, part = "header", where = "bottom")
p = add_border(p, part = "header", where = "top")
p = edit_plot(p, which = "background", gp = gpar(fill = "white"))
sigrows = which(res.mul$pv < 0.05)
p = edit_plot(p, row = sigrows, col = 2, gp = gpar(fontface = "bold"))
p = edit_plot(p, row = sigrows, col = 3, which = "ci", gp = gpar(col = "green3", lwd = 2))
p

png("05.Prognostic/02.Multivar.Cox.png", width = 7, height = 7, units = "in", bg = "white", res = 300)
plot(p)
dev.off()
pdf("05.Prognostic/02.Multivar.Cox.pdf", width = 7, height = 7)
plot(p)
dev.off()


colnames(df.full)[which(colnames(df.full)=="risk")] = "Risk"
data.dist = datadist(df.full)
options(datadist = "data.dist")
model4nomo = rms::cph(Surv(time = OS.time, event = OS)~Risk+Age+T_Stage+N_Stage+M_Stage,
                      data = df.full, surv = T, x = T, y = T)

surv <- Survival(model4nomo)
surv1 <- function(x)surv(365*3,lp=x)
surv2 <- function(x)surv(365*2,lp=x)
surv3 <- function(x)surv(365*1,lp=x)

nomoplot = nomogram(model4nomo, fun = list(surv3,surv2,surv1), lp = F, funlabel = c("1-year Survival", "2-year Survival", "3-year Survival"),
                    maxscale = 100, fun.at=c('0.95', '0.9','0.8','0.7','0.5','0.3','0.1'))
png("05.Prognostic/03.Nomogram.png", width = 10, height = 8, res = 1200, bg = "white", units = "in")
plot(nomoplot, xfrac=0.2)
dev.off()
pdf("05.Prognostic/03.Nomogram.pdf", width = 10, height = 8)
plot(nomoplot, xfrac=0.2)
dev.off()

form = Surv(time = OS.time, event = OS) ~ Risk+Age+T_Stage+N_Stage+M_Stage
model = cph(form, data = df.full, time.inc = 365, surv = T, x = T, y = T)
set.seed(123456)
cal1 = calibrate(model, cmethod="KM", method="boot", u=365, m=50, B=500)
model = cph(form, data = df.full, time.inc = 365*2, surv = T, x = T, y = T)
set.seed(123456)
cal2 = calibrate(model, cmethod="KM", method="boot", u=365*2, m=50, B=500)
model = cph(form, data = df.full, time.inc = 365*3, surv = T, x = T, y = T)
set.seed(123456)
cal3 = calibrate(model, cmethod="KM", method="boot", u=365*3, m=50, B=500)

df.1 = data.frame(Predicted = cal1[,'mean.predicted'], Observed = cal1[,'KM'], 
                  upper = cal1[,'KM'] + cal1[,'std.err'],
                  lower = cal1[,'KM'] - cal1[,'std.err'],
                  year = "1-year")
df.2 = data.frame(Predicted = cal2[,'mean.predicted'], Observed = cal2[,'KM'], 
                  upper = cal2[,'KM'] + cal2[,'std.err'],
                  lower = cal2[,'KM'] - cal2[,'std.err'], 
                  year = "2-year")
df.3 = data.frame(Predicted = cal3[,'mean.predicted'], Observed = cal3[,'KM'],  
                  upper = cal3[,'KM'] + cal3[,'std.err'],
                  lower = cal3[,'KM'] - cal3[,'std.err'],
                  year = "3-year")

ggplot(mapping = aes(x = Predicted, y = Observed, color = year, ymin = lower, ymax = upper)) + 
  geom_line(data = df.1) + 
  geom_line(data = df.2) + 
  geom_line(data = df.3) +
  geom_abline(slope = 1, color = "grey40") + 
  ylim(0,1) + xlim(0,1) +
  scale_color_d3(name = NULL) +
  SCP::theme_scp() +
  theme(legend.position = "inside", legend.position.inside = c(0.75, 0.25))
ggsave("05.Prognostic/04.Calibration.png", width = 5, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("05.Prognostic/04.Calibration.pdf", width = 5, height = 5, units = "in", dpi = 300, bg = "white")


library(ggDCA)
res.mul.risk = coxph(Surv(OS.time,event=OS)~Risk, data = df.full)
res.mul.age = coxph(Surv(OS.time,event=OS)~Age, data = df.full)
res.mul.t = coxph(Surv(OS.time,event=OS)~T_Stage, data = df.full)
res.mul.n = coxph(Surv(OS.time,event=OS)~N_Stage, data = df.full)
res.mul.m = coxph(Surv(OS.time,event=OS)~M_Stage, data = df.full)
res.mul.full = coxph(Surv(OS.time,event=OS)~Risk+Age+T_Stage+N_Stage+M_Stage, data = df.full)

res.dca = dca(res.mul.full, res.mul.age, res.mul.risk, res.mul.t, res.mul.n, res.mul.m, 
              model.names = c("Full Model","Risk","Age","T Stage","N Stage","M Stage"))
ggplot(res.dca, lwd = 0.7, linetype = F) + 
  scale_color_nejm()
ggsave("05.Prognostic/05.DCA.png", width = 7, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("05.Prognostic/05.DCA.pdf", width = 7, height = 5, units = "in", dpi = 300, bg = "white")

res.mul = coxph(
  Surv(time = OS.time, event = OS) ~ Risk + Age + T_Stage + N_Stage + M_Stage,
  data = df.full
)

df.test = read.delim2("interm/GSE39582.risk.xls") |> type.convert(as.is = T)
df.test$Risk = df.test$risk
df.test$Age = as.numeric(df.test$age.at.diagnosis..year..ch1)
df.test$Gender = df.test$Sex.ch1 |> str_to_lower()
df.test$T_Stage = str_match(df.test$tnm.t.ch1, "T[1-4]")[,1]
df.test$N_Stage = str_match(df.test$tnm.n.ch1, "N[0-2]")[,1]
df.test$M_Stage = str_match(df.test$tnm.m.ch1, "M[0-1]")[,1]
df.test$pred = predict(res.mul, df.test, type = "lp")
df.full$pred = predict(res.mul, df.full, type = "lp")


roc.train = timeROC(
  T = df.full$OS.time,
  delta = df.full$OS,
  marker = df.full$pred,
  cause = 1,
  weighting = "marginal",
  times = 365 * c(1:3),
  iid = TRUE
)

labels = roc.train$AUC |> signif(3)
labels = paste0(c("1-year","2-year","3-year"), " (AUC = ", labels, ")")

png("05.Prognostic/06.ROC.Train.png", width = 6, height = 6, units = "in", res = 300)
plot(roc.train, time = 365, col = "red", title = F)
plot(roc.train, time = 365*2, col = "green", title = F, add = T)
plot(roc.train, time = 365*3, col = "blue", title = F, add = T)
legend("bottomright",
       labels,
       col=c("red","green","blue"),
       lty=1,
       lwd=1)
title("TCGA")
dev.off()

pdf("05.Prognostic/06.ROC.Train.pdf", width = 6, height = 6)
plot(roc.train, time = 365, col = "red", title = F)
plot(roc.train, time = 365*2, col = "green", title = F, add = T)
plot(roc.train, time = 365*3, col = "blue", title = F, add = T)
legend("bottomright",
       labels,
       col=c("red","green","blue"),
       lty=1,
       lwd=1)
title("TCGA")
dev.off()


df.test = df.test[!is.na(df.test$pred),]
roc.test = timeROC(
  T = df.test$OS.time,
  delta = df.test$OS,
  marker = df.test$pred,
  cause = 1,
  weighting = "marginal",
  times = 12 * c(1:3),
  iid = TRUE
)
labels = roc.test$AUC |> signif(3)
labels = paste0(c("1-year","2-year","3-year"), " (AUC = ", labels, ")")


png("05.Prognostic/07.ROC.Test.png", width = 6, height = 6, units = "in", res = 300)
plot(roc.test, time = 12, col = "red", title = F)
plot(roc.test, time = 12*2, col = "green", title = F, add = T)
plot(roc.test, time = 12*3, col = "blue", title = F, add = T)
legend("bottomright",
       labels,
       col=c("red","green","blue"),
       lty=1,
       lwd=1)
title("GSE39582")
dev.off()

pdf("05.Prognostic/07.ROC.Test.pdf", width = 6, height = 6)
plot(roc.test, time = 12, col = "red", title = F)
plot(roc.test, time = 12*2, col = "green", title = F, add = T)
plot(roc.test, time = 12*3, col = "blue", title = F, add = T)
legend("bottomright",
       labels,
       col=c("red","green","blue"),
       lty=1,
       lwd=1)
title("GSE39582")
dev.off()


rm(list = ls())
library(ggalluvial)
df.full = read.delim2("interm/risk.group.xls")
df.clst = read.delim2("interm/TCGA.cluster.xls")
df.pheno = read.delim2("cleandata/pheno.tsv")

df.full$stage = df.pheno$ajcc_pathologic_stage.diagnoses[match(df.full$sample, df.pheno$sample)]
df.full$stage = str_match(df.full$stage, "Stage I+V?")[,1]
df.full$cluster = df.clst$cluster[match(df.full$sample, df.clst$sample)]
df.plot = df.full[c("cluster","risk.group","stage")]
rownames(df.plot) = df.full$sample
df.plot = na.omit(df.plot)
df.plot = to_lodes_form(df.plot)
df.plot$stratum = factor(
  df.plot$stratum,
  levels = c(unique(df.clst$cluster), "Low", "High", unique(df.full$stage))
)
cols = c(pal_bmj()(3), pal_futurama()(2), viridis::viridis(4))
names(cols) = levels(df.plot$stratum)

ggplot(
  df.plot,
  aes(
    x = x,
    stratum = stratum,
    alluvium = alluvium,
    fill = stratum,
    label = stratum
  )
) +
  scale_x_discrete(expand = c(.1, 0)) +
  geom_flow(width = 0.55) +
  geom_stratum(alpha = .5, width = 0.55) +
  ggfittext::geom_fit_text(stat = "stratum", width = 0.55, min.size = 3) +
  scale_fill_manual(values = cols) +
  scale_x_discrete(labels = c("Cluster", "Risk", "Stage")) +
  scale_y_continuous(expand = c(0.01, 0)) +
  theme_void() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "in")
  )
ggsave('05.Prognostic/08.Sankey.png', width = 6, height = 6, units = 'in', dpi = 300, bg = 'white')
ggsave('05.Prognostic/08.Sankey.pdf', width = 6, height = 6, units = 'in', dpi = 300, bg = 'white')


