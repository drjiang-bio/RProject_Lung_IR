rm(list = ls());gc()
library(GEOquery)
library(limma)
library(Vennerable)
options(stringsAsFactors = F)

BiocManager::install(c("RBGL","graph"))
install.packages("devtools");library(devtools)
install_github("js229/Vennerable")


eset <- read.table('./origin/GSE75560_series_matrix.txt', header = T, 
                   sep = '\t', comment.char = '!', row.names = 1)
range(eset[,1])
targets <- read.table('./origin/GSE75560/phen_gse75560.txt', 
                      sep = '\t', header = T)
type <- targets$class
names(eset) <- targets$class

gene_diff <- function(data_methy) {
  # 本函数需对beta矩阵 事先排序, 将正常样本至于矩阵前列
  # 本函数需 事先计算 对照（正常）与实验组（癌症）的样本数目
  # dara_methy:处理好的基因总体平均甲基化beta值矩阵，列为样本，行为基因;
  # normalNum: 正常样品的数目; tumorNum: 癌症样品的数目
  # grade: 样本分类信息，使用数字向量(1,1,2,2,2)，对应矩阵样本，须正常样本在前列
  data = data_methy
  res_l <- apply(data, 1, function(x) {
    rt <- rbind(expression=x, grade=c(1,2,2))
    rt <- as.matrix(t(rt))
    wilcoxTest <- wilcox.test(expression ~ grade, data=rt)
    logFC = log2(x[1] + 1)-log2(x[2] + 1)
    return(c(logFC, wilcoxTest$p.value))
    })
  res_df <- data.frame(t(do.call(cbind, res_l)))
  names(res_df) <- c('logFC', 'pvalue')
  #对p值进行矫正
  fdr <- p.adjust(as.numeric(as.vector(res_df[,'pvalue'])), method="fdr")
  res_df <- cbind(res_df, FDR=fdr)
  return(res_df)
}

test <- gene_diff(eset)

table(sapply(1:10000, function(i) {
  wilcox.test(ex~grade, t(rbind(ex=eset[i,], grade=c(1,1,2))))$p.value
}) )

eset[, c(3,2)][1,],cbind()
rbind(ex=eset[, c(3,2)][1,], grade=c(1,2))




targets$class <- c("RNA1","RNA2","RNA3")

f <- factor(targets$class, levels=c("RNA1","RNA2","RNA3"))
f
design <- model.matrix(~0+f)
colnames(design) <- targets$class

fit <- lmFit(esetl, design)
contrast.matrix <- makeContrasts(RNA2-RNA1, RNA3-RNA2, RNA3-RNA1, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, coef=1, adjust="BH")
results <- decideTests(fit2)
vennDiagram(results)
topTableF(fit2, number=30)




l_limma_diff <- function(eset_l,type_l,FC = 2,p=0.05) {
  # 参数均数据框，注意检查数据格式
  # 构建分类信息：1.核对分表表,2.提取实验分类信息
  type <- type_l[match(names(eset_l),type_l[, 1]),] # g
  type <- factor(type[, 2])
  
  # 设计实验矩阵design及对比模型contrast.matrix
  design <- model.matrix(~-1 + type)
  
  duibi <- paste(colnames(design), collapse = ' - ')
  
  contrast.matrix <- makeContrasts(contrasts = duibi,levels = design)
  # 线性模型拟合，据对比模型行差值计算，贝叶斯检验
  fit <- lmFit(eset_l, design)                          # g
  fit1 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit1)
  # 检验结果报表及筛选
  dif <- topTable(fit2, coef = duibi, n = nrow(fit2), lfc = log2(FC))
  if (FC != 1) {
    dif <- dif[dif[, 'P.Value'] < p, ]
  }
  return(dif)
}












# eso data GSE20347, GSE23400, GSE38129
# ----------------------------------------------------
# 1.批量下载数据并保存数据列表
# ----------------------------------------------------
get_gse <- function(gse_id_number) {
  # gse_id_number 格式为 20347
  # 下载gse20347，并返回ExpressionSet对象，并储存在当前目录
  gse_id <- paste('GSE', gse_id_number, sep = '')
  gse <- getGEO(gse_id, destdir = '.', getGPL = F)
  gse <- gse[[1]]
  save(gse, file = paste(gse_id, 'Rdata', sep = '.'))
  return(gse)
}

id_list <- c(20347, 38129, 23400, 23036, 9844)

gse_list <- list()

for (i in 1:length(id_list)) {
  gse_list[[i]] <- get_gse(id_list[i])
}

names(gse_list) <- c(20347, 38129, 23400, 23036, 9844)

# 保存批量下载的gse matrix 列表
save(gse_list, file = 'gse_list.Rdata')
# --------------------------------------------------------------------
# 2.提取、整理exprs,pData数据，并做差异分析
# --------------------------------------------------------------------
load('gse_list.Rdata')

# 提取exprs、处理
eset_list <- lapply(gse_list, exprs)
eset_list <- lapply(eset_list, data.frame)
eset_list[[1]][1:5,1:5]

# 提取pData，处理 
pData_list <- lapply(gse_list, pData)
type_list <- lapply(pData_list, function(x) x[,c(2,1)])

View(pData_list[[1]])
View(type_list[[1]])

type_l_f <- list()
for (i in 1:5) {
  type = type_list[[i]]
  names(type) <- c('GSM_ID', 'class')
  type <- data.frame(apply(type, 2, as.character))
  type[grep('normal|N$', type$class), 2] <- 'normal'  ## g
  type[grep('[^normal]', type$class), 2] <- 'cancer'
  type_l_f[[i]] <- type
}


l_limma_diff <- function(eset_l,type_l,FC = 2,p=0.05) {
  # 参数均数据框，注意检查数据格式
  # 构建分类信息：1.核对分表表,2.提取实验分类信息
  type <- type_l[match(names(eset_l),type_l[, 1]),] # g
  type <- factor(type[, 2])
  
  # 设计实验矩阵design及对比模型contrast.matrix
  design <- model.matrix(~-1 + type)
  
  duibi <- paste(colnames(design), collapse = ' - ')
  
  contrast.matrix <- makeContrasts(contrasts = duibi,levels = design)
  # 线性模型拟合，据对比模型行差值计算，贝叶斯检验
  fit <- lmFit(eset_l, design)                          # g
  fit1 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit1)
  # 检验结果报表及筛选
  dif <- topTable(fit2, coef = duibi, n = nrow(fit2), lfc = log2(FC))
  if (FC != 1) {
    dif <- dif[dif[, 'P.Value'] < p, ]
  }
  return(dif)
}


diff_log21 <- lapply(1:5, function(i) {
  l_limma_diff(eset_list[[i]], type_l_f[[i]])
  })

diff_log20 <- lapply(1:5, function(i) {
  l_limma_diff(eset_list[[i]], type_l_f[[i]], FC = 1)
})

save(diff_log21, file = 'diff_log21.Rdata')
save(diff_log20, file = 'diff_log20.Rdata')

# ------------------------------------------------------------
# 3.注释表达矩阵、差异分析结果
# ------------------------------------------------------------

load('diff_log21.Rdata')
load('gse_list.Rdata')

annot_file <- list.files('.', pattern = '^GPL', full.names = T)
annot_l <- lapply(annot_file, function(x) read.table(x, sep = '\t', header = T))
names(annot_l) <- c(570, 571, 96)


eset_list <- lapply(gse_list, exprs)
eset_list <- lapply(eset_list, data.frame)
eset_list[[1]][1:5,1:5]
names(gse_list) <- c('20347_571', '38129_571', '23400_96', '23036_571', '9844_570')
names(eset_list) <- c('20347_571', '38129_571', '23400_96', '23036_571', '9844_570')

names(diff_log21) <- c('20347_571', '38129_571', '23400_96', '23036_571', '9844_570')

# 注释整个癌症芯片(调用自定义函数geo_pi_annot)

geo_pi_annot <- function(annot_list,exprs_list) {
  # 两参数均要求为list，annot:注释，exprs表达矩阵
  # 此函数判断annot与exprs为行数，需注意此局限性，仔细检查
  
  exprs_ann <- list()
  for (i in 1:length(exprs_list)) {
    # 从annot_list匹配exprs的对应注释数据
    n <- grep(nrow(exprs_list[[i]]),sapply(annot_list,nrow))
    annot <- annot_list[[n]]
    exprs <- exprs_list[[i]]
    # 
    keep_probe <- row.names(exprs) %in% annot[,1]
    exprs <- exprs[keep_probe,]
    exprs <- cbind(annot[match(row.names(exprs),annot[,1]),3], exprs)
    names(exprs)[1] <- 'symbol'
    #
    exprs_ann[[i]] <- exprs
  }
  return(exprs_ann)
}
eset_list_an <- geo_pi_annot(annot_l,eset_list)
names(eset_list_an) <- c('20347_571', '38129_571', '23400_96', '23036_571', '9844_570')

eset_list_an[[2]][1:5,1:5]

geo_pi_annot2 <- function(annot_list,exprs_list) {
  # 两参数均要求为list，annot:注释，exprs表达矩阵
  # 此函数判断annot与exprs为行数，需注意此局限性，仔细检查
  
  exprs_ann <- list()
  for (i in 1:length(exprs_list)) {
    # 从annot_list匹配exprs的对应注释数据
    n <- grep(strsplit(names(exprs_list[i]),'_')[[1]][2], names(annot_list))
    # n <- grep(nrow(exprs_list[[i]]),sapply(annot_list,nrow))
    annot <- annot_list[[n]]
    exprs <- exprs_list[[i]]
    # 
    keep_probe <- row.names(exprs) %in% annot[,1]
    exprs <- exprs[keep_probe,]
    exprs <- cbind(annot[match(row.names(exprs),annot[,1]),3], exprs)
    names(exprs)[1] <- 'symbol'
    #
    exprs_ann[[i]] <- exprs
  }
  return(exprs_ann)
}

diff_log21_an <- geo_pi_annot2(annot_l,diff_log21)
names(diff_log21_an) <- c('20347_571', '38129_571', '23400_96', '23036_571', '9844_570')

save(diff_log21_an, eset_list_an, file = 'diff_eset_annot.Rdata')

# ----------------------------------------------------------
# 4.集合运算，得到食管癌（3张）、头颈部癌（2张）各自的公共差异基因，以及它们两个的交集与差集
# # 可以叫此步看做建模的第一步，即得到待选特征基因
# ----------------------------------------------------------
load('diff_eset_annot.Rdata')
# 集合运算

eso_diff <- diff_log21_an[1:3]
head_diff <- diff_log21_an[4:5]

eso_jiao <- Venn(lapply(eso_diff, function(x) x[,1]))
eso_jiao <- eso_jiao@IntersectionSets$`111`

head_jiao <- Venn(lapply(head_diff, function(x) x[,1]))
head_jiao <- head_jiao@IntersectionSets$`11`

eso_head <- Venn(list(eso = eso_jiao, head = head_jiao))
eso_head_cha <- list(eso_dan = eso_head@IntersectionSets$`10`,
                      head_dan = eso_head@IntersectionSets$`01`)
eso_head_jiao <- eso_head@IntersectionSets$`11`

save(eso_head_cha, eso_head_jiao, head_jiao, eso_jiao, file = 'gene_Set.Rdata')

# ---------------------------------------------
# 5. 训练集构建
# ---------------------------------------------
load('gene_Set.Rdata')
load('diff_eset_annot.Rdata')

# 5.1.1 根据差集提取集合（基因子集）
cha <- as.character(unlist(eso_head_cha))
eset_zi <- lapply(eset_list_an, function(x) x[x[,1]%in% cha,])
# cha_probe_id <- lapply(diff_log21_an, function(x) x[x[,1]%in% cha,])
unique(eset_zi[[5]]$symbol)

# 5.1.2 根据交集提取集合（基因子集）
eset_zi <- lapply(eset_list_an, function(x) x[x[,1]%in% eso_head_jiao,])

rmDupID <- function(a) {
  # 求每行均值排序后去重，设置行名
  exprSet <- a[,-1]
  rowMeans <- apply(exprSet, 1, function(x) mean(as.numeric(x), na.rm=T))
  a <- a[order(rowMeans, decreasing=T),]
  exprSet <- a[!duplicated(a[, 1]),]
  #
  exprSet <- exprSet[!is.na(exprSet[, 1]),]
  rownames(exprSet) = exprSet[, 1]
  exprSet <- exprSet[, -1]
  return(exprSet)
}

# 对集合去重
rm_eset_zi <- lapply(eset_zi, rmDupID)

# 5.2 根据pData提取子集（样本子集）
# ---- from 2 part
load('gse_list.Rdata')

# 提取pData，处理 
pData_list <- lapply(gse_list, pData)
type_list <- lapply(pData_list, function(x) x[,c(2,1)])

View(pData_list[[1]])
View(type_list[[1]])

type_l_f <- list()
for (i in 1:5) {
  type = type_list[[i]]
  names(type) <- c('GSM_ID', 'class')
  type <- data.frame(apply(type, 2, as.character))
  type[grep('normal|N$', type$class), 2] <- 'normal'  ## g
  type[grep('[^normal]', type$class), 2] <- 'cancer'
  type_l_f[[i]] <- type
}

# 取子集
eset_cancer <- lapply(1:5, function(i) {
  nc <- grep('cancer', type_l_f[[i]]$class)
  rm_eset_zi[[i]][, nc]
})

# 5.3 分别合并数据

lapply(eset_cancer, function(x) grep('^$',row.names(x)))

# 食管数据转置、排序
eso_trian <- lapply(eset_cancer[1:3], function(x) {
  a <- data.frame(t(x))
  nc <- order(names(a))
  a[,nc]
  })
# 合并
all(names(eso_trian[[1]]) == names(eso_trian[[3]]))
eso_trian_he <- do.call(rbind, eso_trian) 
eso_trian_he$class <- 'eso'

# 头颈部转置、排序
head_trian <- lapply(eset_cancer[4:5], function(x) {
  a <- data.frame(t(x))
  nc <- order(names(a))
  a[,nc]
})
# 合并
all(names(head_trian[[1]]) == names(head_trian[[2]]))
head_trian_he <- do.call(rbind, head_trian) 
head_trian_he$class <- 'head'

# 食管、头颈部合并
all(names(head_trian_he) == names(head_trian_he))
train <- rbind(eso_trian_he, head_trian_he)
train$class <- as.factor(train$class)

# 保存数据时一定小心到底是差集训练集还是交集训练集，避免错误覆盖本地文件

# 保存结果（差集训练集） # 5.1.1结果
save(train, file = 'train.Rdata')

# 保存数据（交集训练集） # 5.1.2结果
save(train, file = 'train_cha.Rdata')

# --------------------------
load('train.Rdata')







