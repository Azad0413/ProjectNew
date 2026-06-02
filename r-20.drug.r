rm(list = ls())
setwd('/home/data/project/54.LZPJ-3001-6/')
file.dir = '20_drug'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)

library(enrichR)
dbs <- listEnrichrDbs()
##DsigDB
dbs$libraryName
dbs <- c("DSigDB")
hubgene <- read.csv('../05_ML/hubgene.csv')
symbol <- hubgene$symbol
enrichr <- enrichr(symbol,dbs)
result <- data.frame(enrichr$DSigDB)
result <- result[which(result$P.value<0.05),]

write.csv(result,file = 'drug.res.csv',row.names = F)

pdf(file = '01.drug.bar.pdf',width = 8,height = 6)
enrichR::plotEnrich(result,y = 'Ratio',orderBy = "Score",title = '')
dev.off()

png(file = '01.drug.bar.png',width = 8,height = 6,units = 'in',res = 300,bg = 'white')
enrichR::plotEnrich(result,y = 'Ratio',orderBy = "Score",title = '')
dev.off()
