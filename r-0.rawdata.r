rm(list = ls())
setwd('/home/data/project/54.LZPJ-3001-6/')
file.dir = '00_rawdata'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)

library(GEOquery)
library(Biobase)
library(GEOquery)
library(Biobase)
library(tidyverse)

##鼻上皮-------
# GSE65204 
gset<-getGEO("GSE65204",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL14550",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)
colnames(gpl)
probe2symobl<-gpl %>%
  dplyr::select('ID','GENE_SYMBOL')%>%
  filter('GENE_SYMBOL'!='')%>%
  separate('GENE_SYMBOL',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symobl=probe2symobl[probe2symobl$symbol!='',]
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
dat<-dat %>%
  inner_join(probe2symobl,by='ID')%>%
  dplyr::select(-ID)%>%     
  dplyr::select(symbol,everything())%>%   
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    
  arrange(desc(rowMean))%>%       
  distinct(symbol,.keep_all = T)%>%      
  dplyr::select(-rowMean)%>%    
  tibble::column_to_rownames(colnames(.)[1])   
pd<-pData(a)
group <- data.frame(sample = pd$geo_accession,group = pd$`asthma:ch1`)
table(group$group)
group$group <- ifelse(group$group == 'TRUE','Asthma','Control')
group <- group[order(group$group),]
dat <- dat[,group$sample]
range(dat)
write.csv(dat,file = 'dat.GSE65204.csv',row.names = T)
write.csv(group,file = 'group.GSE65204.csv',row.names = F)


# GSE19187-----
gset<-getGEO("GSE19187",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL6244",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)
colnames(gpl)
probe2symobl<-gpl %>%
  dplyr::select('ID','gene_assignment')%>%
  filter('gene_assignment'!='')%>%
  separate('gene_assignment',c('drop','symbol'),sep = ' // ')%>%
  dplyr::select(-drop)
probe2symobl=probe2symobl[probe2symobl$symbol!='',]
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
dat<-dat %>%
  inner_join(probe2symobl,by='ID')%>%
  dplyr::select(-ID)%>%     
  dplyr::select(symbol,everything())%>%   
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    
  arrange(desc(rowMean))%>%       
  distinct(symbol,.keep_all = T)%>%      
  dplyr::select(-rowMean)%>%    
  tibble::column_to_rownames(colnames(.)[1])   
pd<-pData(gset[[1]])
group <- data.frame(sample = pd$geo_accession,group = pd$`disease state:ch1`)
table(group$group)
group <- group[which(group$group %in% c('Healthy','Uncontrolled asthma','Rhinitis with controlled asthma')),]
group$group <- ifelse(group$group == 'Healthy','Control','Asthma')
group <- group[order(group$group),]
dat <- dat[,group$sample]
range(dat)
write.csv(dat,file = 'dat.GSE19187.csv',row.names = T)
write.csv(group,file = 'group.GSE19187.csv',row.names = F)

