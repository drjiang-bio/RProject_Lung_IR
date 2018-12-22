# 支线Affy_cel_read.r
# -------------------------------

# 2.处理样本分类信息
# 目标：第一列为文件名，第二列为文件名简写（常用GSM_ID），第三列为表型
# shell中ls > xx.txt, 打开文件检查，在GEO网页复制等方法
pd <- read.table('./origin/GSE75560/phen_gse75560.txt', sep = '\t', header = T)
head(pd)[1:3, ]
names(pd) <- c('gsm','id', 'class')
# sum(grepl('_$', substr(pd[,1],1,11)))
# pd[,1] <- substr(pd[,1],1,10)
#pd$Target <- 'ESCC'
#pd[,4] <- NULL
write.table(pd, './origin/GSE75560/phen_gse75560.txt', 
            sep = '\t', row.names = F)
