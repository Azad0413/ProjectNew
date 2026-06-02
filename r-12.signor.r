rm(list = ls())
setwd('/home/data/project/54.LZPJ-3001-6/')
file.dir = '12_SIGNOR2'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
# 加载库
library(rsvg)
# 转换 SVG 到 PDF
rsvg_pdf("drawing_pure.svg", "01.signor.pdf")
rsvg_png("drawing_pure.svg", "01.signor.png")
