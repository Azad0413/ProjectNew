rm(list = ls())
library(tidyverse)
library(ggsci)
library(ggpubr)

df.risk = read.delim2("interm/risk.group.xls")
df.pheno = read.delim2("cleandata/pheno.tsv")
df.full = merge(df.risk, df.pheno, by = "sample")
df.full$risk = as.numeric(df.full$risk)
df.full$risk.group = factor(df.full$risk.group, levels = c("Low","High"))

df.full$t = str_match(df.full$ajcc_pathologic_t.diagnoses, "T[1-4]") %>% as.vector
ggplot(subset(df.full, !is.na(t)), aes(x = t, y = risk)) +
  geom_violin(aes(fill = t)) +
  geom_boxplot(fill = "white", width = 0.1) +
  stat_compare_means(label.x = 2.5, hjust = 0.5) +
  scale_fill_bmj() +
  guides(fill = guide_none()) +
  xlab(NULL) + ylab("Risk") + theme_classic2()
ggsave("04.Clinic/01.Tstage.png", width = 4, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("04.Clinic/01.Tstage.pdf", width = 4, height = 5, units = "in", dpi = 300, bg = "white")

df.full$n = str_match(df.full$ajcc_pathologic_n.diagnoses, "N[0-3]")
ggplot(subset(df.full, !is.na(n)), aes(x = n, y = risk)) +
  geom_violin(aes(fill = n)) +
  geom_boxplot(fill = "white", width = 0.3) +
  stat_compare_means(label.x = 2, hjust = 0.5) +
  scale_fill_bmj() +
  guides(fill = guide_none()) +
  xlab(NULL) + ylab("Risk") + theme_classic2()
ggsave("04.Clinic/02.Nstage.png", width = 4, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("04.Clinic/02.Nstage.pdf", width = 4, height = 5, units = "in", dpi = 300, bg = "white")

df.full$m = str_match(df.full$ajcc_pathologic_m.diagnoses, "M[0-1]")
ggplot(subset(df.full, !is.na(m)), aes(x = m, y = risk)) +
  geom_violin(aes(fill = m)) +
  geom_boxplot(fill = "white", width = 0.3) +
  stat_compare_means(label.x = 1.5, hjust = 0.5) +
  scale_fill_bmj() +
  guides(fill = guide_none()) +
  xlab(NULL) + ylab("Risk") + theme_classic2()
ggsave("04.Clinic/03.Mstage.png", width = 4, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("04.Clinic/03.Mstage.pdf", width = 4, height = 5, units = "in", dpi = 300, bg = "white")

ggplot(subset(df.full, !is.na(gender.demographic)), 
       aes(x = str_to_sentence(gender.demographic), y = risk)) +
  geom_violin(aes(fill = str_to_sentence(gender.demographic))) +
  geom_boxplot(fill = "white", width = 0.3) +
  stat_compare_means(label.x = 1.5, hjust = 0.5) +
  scale_fill_bmj() +
  guides(fill = guide_none()) +
  xlab(NULL) + ylab("Risk") + theme_classic2()
ggsave("04.Clinic/04.Gender.png", width = 4, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("04.Clinic/04.Gender.pdf", width = 4, height = 5, units = "in", dpi = 300, bg = "white")

df.full$Age = as.numeric(df.full$age_at_index.demographic)
df.full = subset(df.full, !is.na(Age))
lab.cor = cor(df.full$risk, df.full$Age, method = "spearman") %>% round(2)
lab.p = cor.test(df.full$risk, df.full$Age, method = "spearman")$p.value %>% signif(3)
lab = paste0("R = ",lab.cor,", ","p = ",lab.p)
ggscatter(df.full, x = "risk", y = "Age",
          color = "risk.group", shape = 20, size = 3, 
          add = "reg.line",
          add.params = list(color = "blue", fill = "lightgray"), 
          conf.int = TRUE) +
  scale_color_aaas(guide = guide_none()) + 
  ggtitle(lab) + xlab("Risk") + ylab("Age") +
  theme_bw(base_size = 14) + 
  theme(plot.title = element_text(hjust = 0.5, face = "italic", size = 12), aspect.ratio = 1)
ggsave("04.Clinic/05.Age.png", width = 5, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("04.Clinic/05.Age.pdf", width = 5, height = 5, units = "in", dpi = 300, bg = "white")


