rm(list = ls())
setwd('/home/data/project/54.LZPJ-3001-6/')
file.dir = '17_BioGRID'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
png_files <- list.files(pattern = "\\.png$", full.names = TRUE)

# 遍历所有 PNG 文件并转换为 PDF
for (png_file in png_files) {
  # 获取文件名并构造新的 PDF 文件名
  pdf_file <- sub("\\.png$", ".pdf", png_file)
  
  # 打开一个新的 PDF 文件
  pdf(pdf_file, width = 6, height = 6)  # 可以根据需要调整图像大小
  
  # 读取 PNG 图像并在 PDF 中显示
  img <- png::readPNG(png_file)  # 读取 PNG 图像
  grid::grid.raster(img)  # 在 PDF 中显示图像
  
  # 关闭 PDF 设备
  dev.off()
  
  # 保留原图：无需额外操作，文件保持在原目录下
}

cat("All PNG files have been converted to PDF!\n")
