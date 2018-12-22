read_affy_norma <- function(gse) {
	# 此函数的需要两层工作目录，总工作目录，次级工作目录（保存有原始芯片RAW.tar数据）
	# 启动R时应将工作目录设置为次级工作目录的父目录（总目录）
	# 若如此设置目录关系，可以多次运行函数，批量读入，写出数据，互不干扰
	# 读入原始数据需要小号大量计算机资源，内存危机，若内存不够可先设置虚拟内存硬上
	
	# phenodata.txt 文件需要预先做好，第一列(无列名):cel文件名,后续列(有列名)：其他信息【建议3列】
	# gse格式为：'GSE30784'

	# 此函数仅返回cel_rma(RMA标准化后的ExpressionSet对象)
	# 余重要数据或结果写入次级工作目录gse_id
	# 1.affybatch对象（原始cel读入后的返回对象）
	# 2.RMA标准化后的表达矩阵文件及Rdata（包括cel_rma）
	# 3.标准化前后的质量控制相关可视化图片
	# ----------------------------

	# 1.处理自己下好的数据
	gse_id <- gse  # 输入
	
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
	out_rma <- paste(gse_id, 'rma.Rdata', sep = '_')
	save(cel_rma,eset_rma, file = out_rma)
	write.table(eset_rma, file = out_txt, sep = '\t', row.names = T, col.names = T)
	
	# 返回结果
	return(cel_rma)
	# 还原工作目录（重要！！！）
	setwd('../')
}