rm(list = ls())
setwd('/home/data/project/54.LZPJ-3001-6/')
file.dir = '21_Disease'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
library(readr)
file <- list.files()

gene <- read.csv('../05_ML/hubgene.csv')$symbol
stopifnot(length(file) == length(gene))
all.dat <- mapply(function(f, g) {
  dat <- read_tsv(f)
  dat$Gene <- g
  return(dat)
}, f = file, g = gene, SIMPLIFY = FALSE)


names(all.dat) <- gene

ctd <- do.call(rbind,all.dat)
# Diabetic Nephropathies。Diabetic Retinopathy  Diabetic Foot  Diabetic Ketoacidosis
ctd <- ctd[which(ctd$`Disease Name` %in% c('Diabetic Nephropathies','Diabetic Retinopathy','Diabetic Foot','Diabetic Ketoacidosis')),]
colnames(ctd)
library(dplyr)
library(ggplot2)

bar.dat <- ctd[,c(7,1,5)] %>% as.data.frame()
bar.dat$Group <- bar.dat$`Disease Name`
bar.dat$`Inference Score` <- as.numeric(bar.dat$`Inference Score`)
colnames(bar.dat) <- c('Symbol','Disease','Score','Group')

# 每个基因保留 Score 最高的前 5 个疾病
bar.dat.top5 <- bar.dat %>%
  group_by(Symbol) %>%
  slice_max(order_by = Score, n = 5) %>%
  ungroup()
unique(bar.dat.top5$Disease)
# 画图
p <- ggplot(bar.dat.top5, aes(x = Score, y = Symbol, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c(
    "Lung Neoplasms" = "#3A85A8",
    "Lung Diseasesy" = "#BC8FA7",
    "Lung Injury" = '#629363',
    "Asthma" = '#A24057',
    'Pulmonary Fibrosis' = '#A6D854',
    'Pulmonary Edema' = '#FC8D62',
    'Pneumonia' = '#E78AC3'
  )) +
  labs(x = "Inference Score", y = "Hub Gene", fill = "") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 12),
    legend.position = "right",
    panel.grid.major.y = element_blank()
  )

p
ggsave(filename = '01.ctd.pdf',p,w=7,h=4)
ggsave(filename = '01.ctd.png',p,w=7,h=4,units = 'in',dpi = 300,bg = 'white')
