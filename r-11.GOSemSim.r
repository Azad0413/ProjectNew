rm(list = ls())
setwd('/home/data/project/54.LZPJ-3001-6/')
file.dir = '11_GOSemSim'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)

hubgene <- read.csv('../05_ML/hubgene.csv')
library(GOSemSim)
library(org.Hs.eg.db)
library(clusterProfiler)
bp<-godata('org.Hs.eg.db',ont = "BP",computeIC = F)
cc<-godata('org.Hs.eg.db',ont = "CC",computeIC = F)
mf<-godata('org.Hs.eg.db',ont = "MF",computeIC = F)
geneid2<-hubgene$symbol
gene_transform2<-bitr(geneid2,fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)
genelist<-gene_transform2
simbp<-mgeneSim(genelist$ENTREZID,
                semData = bp,
                measure = "Wang",
                drop = NULL,
                combine = "BMA")
simcc<-mgeneSim(genelist$ENTREZID,
                semData = cc,
                measure = "Wang",
                drop = NULL,
                combine = "BMA")

simmf<-mgeneSim(genelist$ENTREZID,
                semData = mf,
                measure = "Wang",
                drop = NULL,
                combine = "BMA")

fsim<-(simcc*simbp)^(1/3)
View(fsim)
colnames(fsim)<-genelist$SYMBOL
rownames(fsim)<-genelist$SYMBOL

for (i in 1:ncol(fsim)) {
  fsim[i,i] <- NA
}
dat_sim<-reshape2::melt(fsim)
dat_sim<-dat_sim[!is.na(dat_sim$value),]
dat_sim<-dat_sim[,c(1,3)]
head(dat_sim)
dat.median<-aggregate(value~Var1,dat_sim,median)
dat.mean<-aggregate(value~Var1,dat_sim,mean)
m<-dat.mean$value
m2<-dat.median$value
names(m)<-dat.mean$Var1
names(m2)<-dat.median$Var1
dat_sim$Var1<-factor(dat_sim$Var1,
                     levels = names(sort(m)))
str(dat_sim)
write.csv(fsim,file = 'GOSemSim.csv',row.names = T,quote = F)
pdf(file = '01.GOSemSim.pdf',w=5,h=3)
ggplot(dat_sim,
       aes(x=Var1,y=value,fill=factor(Var1)))+
  scale_fill_brewer(palette = "Set2")+
  geom_boxplot()+
  coord_flip()+       # 坐标轴互换
  xlab("")+ylab("")+
  theme_bw()+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 18),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13))+
  guides(fill='none')
dev.off()

png(file = '01.GOSemSim.png',w=5,h=3,units = 'in',res = 300,bg = 'white')
ggplot(dat_sim,
       aes(x=Var1,y=value,fill=factor(Var1)))+
  scale_fill_brewer(palette = "Set2")+
  geom_boxplot()+
  coord_flip()+       # 坐标轴互换
  xlab("")+ylab("")+
  theme_bw()+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 18),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13))+
  guides(fill='none')
dev.off()

