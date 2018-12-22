# 主题：Affy原始芯片数据读入、质量控制、标准化
# 此脚本为主线脚本（重在整体分析流程）
# 同目录下有一支线脚本（其他必须的处理及相关自定义函数所在）_supply.r$
# 自定义函数也可能位于同目录一单独文件,^function_
# 目录详解见Readme.txt
# -------------------------------
getwd() # 酌情设置工作目录
rm(list = ls())
options(stringsAsFactors = F)
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
chooseBioCmirror()
4
if (!requireNamespace("GEOquery")) BiocManager::install("GEOquery")
n
if (!requireNamespace("simpleaffy")) BiocManager::install("simpleaffy")
n
library(GEOquery)

accession <- 'GSE75560'
# 1.通过R下载GEO原始芯片数据，并解打包、解压
getGEOSuppFiles(accession)
fpath <- paste0('./origin/', accession, '_RAW.tar')

exdir_path = paste0('origin/', accession)
dir.create(exdir_path) # 创建文件夹

untar(fpath, exdir = exdir_path)
cels <- list.files(exdir_path, pattern = "[gz]")
sapply(paste0(exdir_path, '/', cels), gunzip)

# 1.处理自己下好的数据
dir.create(exdir_path)

untar(fpath, exdir = exdir_path)
cels <- list.files(exdir_path, pattern = "[gz]")
sapply(paste0(exdir_path, '/', cels), gunzip) ## 循环处理数据

# 2.处理样本分类信息(见supply文件)

# 3. 使用simleaffy:read.affy()函数读入cel文件(affy芯片原始数据文件)
library(simpleaffy)
celfiles <- read.affy(covdesc="phen_gse75560.txt", path = exdir_path)
celfiles
phenoData(celfiles)

# 保存重要的中间数据（读入后的affy原始芯片数据）#
save(celfiles,file = '69925affybatch.Rdata')

#---------------------------------
# 清空数据，载入数据
library(simpleaffy)
load('69925affybatch.Rdata')

# -----------
# 标准化前可视化
# -----------
# 载入色彩包
library(RColorBrewer)
library(affyPLM)
# 回归计算
Pset <- fitPLM(celfiles)
# 设置调色板
colors <- brewer.pal(12,"Set3") 
cols <- brewer.pal(8, "Set1")

# 密度图
pdf(file = "hist_cel.pdf",width=12,height=8)
hist(celfiles, main="Original_Hist", col=cols)
dev.off()

# R&N图&deg图
pdf(file = "R&N&_cel.pdf",width=40,height=7)
# RLE箱线图
Mbox(Pset,ylim = c(-1,1),col = colors,main="Original_RLE",las = 3)
# NUSE 箱线图
boxplot(Pset,ylim = c(0.95,1.22),col = colors,main="Original_NUSE",las = 3)
dev.off()

# 标准化

cel_mas5 <- mas5(celfiles)
# cel_gcrma <- gcrma(celfiles)
cel_rma <- rma(celfiles)

# -----------------
# 标准化后可视化
# -----------------
# 标准化前已经将所需包及颜色准备

# 密度图
pdf(file = "hist_rma.pdf",width=12,height=8)
hist(cel_rma, main="Original_Hist", col=cols)
# hist(cel_mas5, main="mas5_Hist", col=cols)
dev.off()

# R&N图&deg图
pdf(file = "R&N&_rma.pdf",width=40,height=7)
# RLE箱线图
Mbox(cel_rma,ylim = c(-1,1),col = colors,main="rma_RLE",las = 3)
# NUSE 箱线图
boxplot(cel_rma,ylim = c(0.5,1.5),col = colors,main="rma_NUSE",las = 3)
dev.off()

# 提取表达矩阵
eset_rma <- exprs(cel_rma)
# eset_mas5 <- exprs(cel_mas5)
eset_rma[1:4,1:4]
dimnames(eset_rma)[[2]] <- gsub('\\.CEL','',dimnames(eset_rma)[[2]])

# 输出结果 #
write.table(eset_rma, file = 'GSE69925_rma.txt', sep = '\t', row.names = T, col.names = T)


# -----------------------------
# 30784 ------------------自动化，只需要更改gse_id即可-------已写成函数----------
# -----------------------------

# 1.处理自己下好的数据
gse_id <- 'GSE30784'  # 输入

# 进入亚工作目录（分析完后切记要跳回父目录！！！）
setwd(gse_id)

file <- paste(gse_id, 'RAW.tar', sep = '_') #'GSE30784_RAW.tar'
cel_data <- paste('data', gse_id, sep = '_') # data_GSE30784

untar(file, exdir= cel_data)
cels <- list.files(cel_data, pattern = "[gz]")
## 循环处理数据
sapply(paste(cel_data, cels, sep="/"), gunzip)

# 2.处理样本分类信息(见supply文件)
## 若能在shell中或excel中处理好，则可不用R处理，
## 处理好后统一命名phenodata.txt置于cel_data目录下

# 3. 使用simleaffy:read.affy()函数读入cel文件(affy芯片原始数据文件)
library(simpleaffy)
celfiles <- read.affy(covdesc="phenodata.txt", path = cel_data)
celfiles
phenoData(celfiles)

# 保存重要的中间数据（读入后的affy原始芯片数据）#
out_bacth <- paste(gse_id, 'affybatch.Rdata', sep = '_')
save(celfiles,file = out_bacth)

# --------------------
# 清空数据，载入数据
# --------------------
library(simpleaffy)
# load(out_bacth)

# -----------
# 标准化前可视化
# -----------
# 载入色彩包
library(RColorBrewer)
library(affyPLM)
# 回归计算
Pset <- fitPLM(celfiles)
# 设置调色板
colors <- brewer.pal(12,"Set3") 
cols <- brewer.pal(8, "Set1")

# 密度图
pdf(file = "hist_cel.pdf",width=12,height=8)   #g
hist(celfiles, main="Original_Hist", col=cols)
dev.off()

# R&N图&deg图
pdf(file = "R&N&_cel.pdf",width=40,height=7)
# RLE箱线图
Mbox(Pset,ylim = c(-1,1),col = colors,main="Original_RLE",las = 3)
# NUSE 箱线图
boxplot(Pset,ylim = c(0.95,1.22),col = colors,main="Original_NUSE",las = 3)
dev.off()

# 标准化

# cel_mas5 <- mas5(celfiles)
# cel_gcrma <- gcrma(celfiles)
cel_rma <- rma(celfiles)

# -----------------
# 标准化后可视化
# -----------------
# 标准化前已经将所需包及颜色准备

# 密度图
pdf(file = "hist_rma.pdf",width=12,height=8)
hist(cel_rma, main="Original_Hist", col=cols)
# hist(cel_mas5, main="mas5_Hist", col=cols)
dev.off()

# R&N图&deg图
pdf(file = "R&N&_rma.pdf",width=40,height=7)
# RLE箱线图
Mbox(cel_rma,ylim = c(-1,1),col = colors,main="rma_RLE",las = 3)
# NUSE 箱线图
boxplot(cel_rma,ylim = c(0.5,1.5),col = colors,main="rma_NUSE",las = 3)
dev.off()

# 提取表达矩阵
eset_rma <- exprs(cel_rma)
# eset_mas5 <- exprs(cel_mas5)

eset_rma[1:4,1:4]
dimnames(eset_rma)[[2]] <- gsub('\\.CEL','',dimnames(eset_rma)[[2]])

# 输出结果 #
out_txt <- paste(gse_id, 'rma.txt', sep = '_')
write.table(eset_rma, file = out_txt, sep = '\t', row.names = T, col.names = T)

# 还原工作目录（重要！！！）
setwd('../')

# ----------------------------------------
# 


