rm(list = ls());gc()
library(GEOquery)
library(limma)
library(Vennerable)
options(stringsAsFactors = F)

eset <- read.table('./origin/GSE75560_series_matrix.txt', header = T, 
                   sep = '\t', comment.char = '!', row.names = 1)
range(eset[,1])
targets <- read.table('./origin/GSE75560/phen_gse75560.txt', 
                      sep = '\t', header = T)

annot <- read.csv('origin/Porcine.na36.annot.csv', comment.char = '#')
annot2 <- annot[,c(1,15)]
annot2 <- annot2[!grepl('^---$', annot2$Gene.Symbol),]

# 方差分析
res <- apply(eset, 1, function(x) {
  rt <- data.frame(t(rbind(expression=x, grade=c(1,2,3))))
  fit <- aov(expression ~ grade, data=rt)
  summ <- summary(fit)
  pvalue <- summ[[1]]$`Pr(>F)`[1]
  return(pvalue)
})

genedf <- data.frame(res, check.names = F)
genedf$symbol <- annot2[match(rownames(genedf), annot2$Probe.Set.ID), 2]
genedf <- na.omit(genedf)
names(genedf)[1] <- 'pval'
genes <- genedf[genedf$pval<0.05,]
genes <- dplyr::arrange(genes, pval)
genes <- genes[!duplicated(genes$symbol),]
write.csv(genes, file = './result/genes_aov_p05.csv')

# 注释
# 1
library(dplyr)
eset$symbol <- annot2[match(rownames(eset), annot2$Probe.Set.ID), 2]
eset <- na.omit(eset)
library(data.table)
eset <- data.table(eset)
eset <- eset[, lapply(.SD, median), by=symbol]
eset <- data.frame(eset)
eset$bi32 <- eset[,3] / eset[,2]
eset$bi43 <- eset[,4] / eset[,3]

index1 <- eset$bi32 < 1.2 & eset$bi32 > 0.8
sum(index1, na.rm = T)
index2 <- eset$bi43 > 2 | eset$bi43 < 0.5
index <- index1 & index2
sum(index, na.rm = T)

eset_s <- eset[index, ]
eset_s <- eset_s[eset_s$bi32 != 0,]
eset_s <- eset_s[eset_s$bi32 != Inf,]
eset_s <- na.omit(eset_s)
eset_s <- arrange(eset_s, desc(bi43))
write.csv(eset_s, file='./result/eset_double_filter.csv', row.names = F)
eset_sig <- read.csv('./result/eset_double_filter.csv')
library(clusterProfiler)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
library(org.Ss.eg.db)
keytypes(org.Ss.eg.db)
ego <- enrichGO(gene = eset_sig$symbol, OrgDb = org.Ss.eg.db,
                keyType = 'SYMBOL', ont = "ALL",
                pAdjustMethod = "BH", pvalueCutoff = 0.05,
                qvalueCutoff = 0.5, minGSSize = 1)
ego@result

ids <- bitr(eset_sig$symbol, fromType="SYMBOL", toType="UNIPROT", OrgDb="org.Ss.eg.db")
kk <- enrichKEGG(gene = ids,
                 organism = 'ssc',
                 pvalueCutoff = 0.05)
file <- list.files('./result/', pattern = '^david', full.names = T)
david <- lapply(file[1:4], function(x) read.table(x, header = T, sep = '\t'))
file[1:4]
names(david) <- c('BP', 'CC', 'KEGG', 'MF')
david <- do.call(rbind, david)
david <- dplyr::arrange(david, Category, PValue)
dav <- dplyr::filter(david, PValue < 0.05)

write.csv(dav, file = './result/david_all.csv', row.names = F)

library(ggplot2)
ggplot(data = dav, aes(x=reorder(Term, PValue), y=-log10(PValue), fill=Category)) + 
  geom_bar(stat = 'identity', position = 'dodge', width = 0.7) +
  scale_y_continuous(expand = c(0,0)) + 
  coord_flip() +
  theme_bw() + # ggtitle('Biologiacl Progress') +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(face = 'italic', colour = 'black', 
                                 size = 12),
        axis.title = element_text(face = 'italic', colour = 'black', 
                                  size = 12),
        plot.title = element_text(vjust = -7, hjust = -0.42),
        panel.border = element_blank(), 
        axis.line = element_line(colour = 'black', size = 1, 
                                 lineend = 'square'))# +
  facet_grid(. ~ Category, scales = 'free') +
  theme(strip.text = element_text(face = 'bold', size = rel(1.1)),
        strip.background = element_rect(fill = 'lightblue', 
                                        colour = 'black', size = 1))
ggsave('./result/david_all.pdf', width = 12, height = 8, units = 'in')














eset$bi32 <- eset[,3] / eset[,2]
eset$ca32 <- eset[,3] - eset[,2]


eset_or <- eset[order(eset$bi32), ]
eset_q <- dplyr::filter(eset_or, bi32 != 0 & bi32 != Inf)
write.csv(genes, file = './result/genes_aov_p05.csv')

aov_eset <- merge(genes, eset_q, by='symbol')
library(data.table)
aov_eset <- data.table(aov_eset)
tt <- aov_eset[, lapply(.SD, median), by=symbol]
write.csv(tt, file = 'result/genes_anov_p05_eset_unique.csv')
sum(tt$bi32  > 1.5 | tt$bi32 < 0.5)
sum(tt$bi32 < 0.8)
range(tt$bi32)

eset_sig <- rbind(head(eset_q, 200), tail(eset_q, 200))
eset_sig <- dplyr::arrange(eset_sig, bi32)
eset_sig2 <- with(eset_q, eset_q[bi32 > 1 | bi32 < 0.5,])

write.csv(eset_sig, file = './result/eset_sig400.csv')
eset_sig <- read.csv('./result/eset_sig400.csv', row.names = 1)
library(clusterProfiler)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
library(org.Ss.eg.db)
ego <- enrichGO(gene = eset_sig$symbol, OrgDb = org.Ss.eg.db,
                 keyType = 'SYMBOL', ont = "ALL",
                 pAdjustMethod = "BH", pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
?enrichGO()
ego@result
BiocManager::install('org.Ss.eg.db')


mb <- c('ADAR1', 'Eif2ak2', 'PKR', 'Mavs', 'MAVS', 'Ifih1', 
        'MDA5', 'Ddx58', 'RIG-I', 'NTRK1')
mb <- read.csv('origin/g.csv', header = F)
mbl <- unlist(as.list(mb))
mbl <- mbl[!grepl('^$', mbl)]
MB <- toupper(mbl)
MB <- as.character(MB)

have <- MB[MB %in% annot2$Gene.Symbol]
'ADAR' %in% annot2$Gene.Symbol
library(dplyr)
unique(annot2$Gene.Symbol) %>% length()

annot <- read.csv('origin/Porcine.na36.annot.csv', comment.char = '#')
annot2 <- annot[,c(1,15)]
annot2 <- annot2[!grepl('^---$', annot2$Gene.Symbol),]


eset$symbol <- annot2[match(rownames(eset), annot2$Probe.Set.ID), 2]
eseta <- na.omit(eset)
eset_s <- eseta[eseta$symbol %in% have, ]
eset_s$bi23 <- eset_s[,2] / eset_s[,3]
eset_s$cha23 <- eset_s[,2] - eset_s[,3]

eset_s1 <- eseta[eseta$symbol == 'ADAR', ]
eset_s1$bi23 <- eset_s1[,3] / eset_s1[,2]

gene_n_s <- na.omit(annot2[match(gene_n, annot2$Probe.Set.ID), 2])

intersect(mb, annot2$Gene.Symbol)
intersect(MB, annot2$Gene.Symbol) %in% gene_n_s


data(litter, package = 'multcomp')
litter[1:4,1:4]


res_l <- apply(eset[,2:3], 1, function(x) {
  rt <- rbind(expression=x, grade=c(1,2))
  rt <- as.matrix(t(rt))
  wilcoxTest <- wilcox.test(expression ~ grade, data=rt)
  logFC = log2(x[1] + 1)-log2(x[2] + 1)
  return(c(logFC, wilcoxTest$p.value))
})
res_df <- data.frame(t(res_l))
names(res_df) <- c('logFC', 'pvalue')
res_df <- with(res_df, res_df[order(logFC),])
sig <- with(res_df, res_df[abs(logFC) > 1, ])

annot <- read.csv('origin/Porcine.na36.annot.csv', comment.char = '#')
annot2 <- annot[,c(1,15)]
annot2 <- annot2[!grepl('^---$', annot2$Gene.Symbol),]
intersect(annot2$Probe.Set.ID, rownames(eset))
names(annot)
mb <- c('ADAR1', 'Eif2ak2', 'PKR', 'Mavs', 'MAVS', 'Ifih1', 'MDA5', 'Ddx58', 'RIG-I', 'NTRK1')
MB <- toupper(mb)
sig$symbol <- annot2[match(rownames(sig), annot2$Probe.Set.ID), 2]
MB %in% annot2$Gene.Symbol
MB %in% sig$symbol

eset2 <- eset
eset2$symbol <- annot2[match(rownames(eset2), annot2$Probe.Set.ID), 2]
eset2 <- na.omit(eset2)
eset2$d23 <- eset2[,2] / eset2[,3]
eset2 <- eset2[eset2$d23 > 2 | eset2$d23 < 0.5,]

intersect(eset2$symbol, sig$symbol)

2^-1
mb %in% eset$symbol
e <- eset
table(res_df$pvalue)

test <- gene_diff(eset)
i=1
table(sapply(1:10000, function(ii) {
  dt <- t(rbind(ex=eset[,2:3][ii,], grade=c(1,2)))
  wilcox.test(ex~grade, dt)$p.value
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







