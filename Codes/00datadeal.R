rm(list = ls())

foldpath <- paste0("D:/mydata/GEO/colorectal_cancer/")

if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)
library(GEOquery)
library(idmap3)
options(timeout = 99999999)
library(readxl)

# ----------------------------------GSE17536------------------
rm(list = ls())
GSEID <- "GSE17536"
gset <- getGEO(filename = paste0(GSEID,"_series_matrix.txt.gz"),AnnotGPL = T, destdir = "./")
# gset <- getGEO("GSE106090",AnnotGPL = T, destdir = "./",getGPL = T)

gset[1]

# 提取表达量
exp <- exprs(gset)

# 临床信息
cli <- pData(gset)  #获取临床信息
save(cli,file = paste0(GSEID,"_cliall.RData"))

head(cli)

gpl<- fData(gset)

# 执行下面代码后，便通过只取第一个Gene symbol来达到去重效果，这时再打开gpl发现已经没有对应多个Gene symbol的情况了
gpl$`Gene symbol`<-data.frame(sapply(gpl$`Gene symbol`,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]

# 接下来对表达矩阵exp进行探针ID转化，首先将其转换为数据框的形式，否则一会会报错
exp<-as.data.frame(exp)

# 此时虽然我们能看到exp的探针ID号（最左侧），但这个ID号并不是真正在列中，所以我们还需要将探针ID放到数据框里面，执行下面命令
exp$ID<-rownames(exp)	# 增加新的一列（最后一列），存放基因ID信息

gpl<-gpl[,c(1,3)]
# 接下来通过将gpl合并到exp中，从而在exp中将芯片ID号转化为Gene symbol
exp_symbol<-merge(exp,gpl,by="ID")

# 首先查看有多少重复的Gene symbol
table(duplicated(exp_symbol$`ID`))
# 输入以下命令进行去重操作同时将矩阵转换成标准的表达矩阵形式
# exp_unique<-avereps(exp_symbol[,-c(1,ncol(exp_symbol))],ID=exp_symbol$`ID`)

exp_symbol <- exp_symbol[,-1]

# 删除重复的
exp_symbol <- exp_symbol[!duplicated(exp_symbol$`Gene symbol`),]
exp_symbol <- na.omit(exp_symbol)
rownames(exp_symbol) <- exp_symbol$`Gene symbol`

exp_symbol <- exp_symbol[,-ncol(exp_symbol)]

table(cli$characteristics_ch1)

GSE17536 <- as.data.frame(t(exp_symbol))
GSE17536 <- GSE17536[rownames(cli),]
GSE17536$ID <- rownames(GSE17536)

surv <- cli[,c("geo_accession",
               "overall survival follow-up time:ch1",
               "overall_event (death from any cause):ch1")]
colnames(surv) <- c("ID","OS.time","OS")
surv$OS.time <- as.numeric(surv$OS.time)
surv$OS.time <- round(surv$OS.time*30)
surv$OS <- ifelse(surv$OS == "no death",0,1)

GSE17536 <- merge(surv,GSE17536,by = "ID")

save(GSE17536,file = "GSE17536.RData")
write.csv(exp_symbol,paste0("dat_",GSEID,".csv"))

# ----------------------------------GSE17537------------------------------
rm(list = ls())
GSEID <- "GSE17537"
gset <- getGEO(filename = paste0(GSEID,"_series_matrix.txt.gz"),AnnotGPL = T, destdir = "./")
# gset <- getGEO("GSE106090",AnnotGPL = T, destdir = "./",getGPL = T)

gset[1]

# 提取表达量
exp <- exprs(gset)

# 临床信息
cli <- pData(gset)  #获取临床信息
save(cli,file = paste0(GSEID,"_cliall.RData"))

head(cli)

gpl<- fData(gset)

# 执行下面代码后，便通过只取第一个Gene symbol来达到去重效果，这时再打开gpl发现已经没有对应多个Gene symbol的情况了
gpl$`Gene symbol`<-data.frame(sapply(gpl$`Gene symbol`,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]

# 接下来对表达矩阵exp进行探针ID转化，首先将其转换为数据框的形式，否则一会会报错
exp<-as.data.frame(exp)

# 此时虽然我们能看到exp的探针ID号（最左侧），但这个ID号并不是真正在列中，所以我们还需要将探针ID放到数据框里面，执行下面命令
exp$ID<-rownames(exp)	# 增加新的一列（最后一列），存放基因ID信息

gpl<-gpl[,c(1,3)]
# 接下来通过将gpl合并到exp中，从而在exp中将芯片ID号转化为Gene symbol
exp_symbol<-merge(exp,gpl,by="ID")

# 首先查看有多少重复的Gene symbol
table(duplicated(exp_symbol$`ID`))
# 输入以下命令进行去重操作同时将矩阵转换成标准的表达矩阵形式
# exp_unique<-avereps(exp_symbol[,-c(1,ncol(exp_symbol))],ID=exp_symbol$`ID`)

exp_symbol <- exp_symbol[,-1]

# 删除重复的
exp_symbol <- exp_symbol[!duplicated(exp_symbol$`Gene symbol`),]
exp_symbol <- na.omit(exp_symbol)
rownames(exp_symbol) <- exp_symbol$`Gene symbol`

exp_symbol <- exp_symbol[,-ncol(exp_symbol)]

GSE17537 <- as.data.frame(t(exp_symbol))
GSE17537 <- GSE17537[rownames(cli),]
GSE17537$ID <- rownames(GSE17537)

surv <- cli[,c("geo_accession",
               "overall survival follow-up time:ch1",
               "overall_event (death from any cause):ch1")]
colnames(surv) <- c("ID","OS.time","OS")
surv$OS.time <- as.numeric(surv$OS.time)
range(surv$OS.time)
surv$OS.time <- round(surv$OS.time*30)
surv$OS <- ifelse(surv$OS == "no death",0,1)

GSE17537 <- merge(surv,GSE17537,by = "ID")

save(GSE17537,file = paste0(GSEID,".RData"))
write.csv(exp_symbol,paste0("dat_",GSEID,".csv"))

# ----------------------------------GSE29621------------------------------
rm(list = ls())
GSEID <- "GSE29621"
gset <- getGEO(filename = paste0(GSEID,"_series_matrix.txt.gz"),AnnotGPL = T, destdir = "./")
# gset <- getGEO("GSE106090",AnnotGPL = T, destdir = "./",getGPL = T)

gset[1]

# 提取表达量
exp <- exprs(gset)

# 临床信息
cli <- pData(gset)  #获取临床信息
save(cli,file = paste0(GSEID,"_cliall.RData"))

head(cli)

gpl<- fData(gset)

# 执行下面代码后，便通过只取第一个Gene symbol来达到去重效果，这时再打开gpl发现已经没有对应多个Gene symbol的情况了
gpl$`Gene symbol`<-data.frame(sapply(gpl$`Gene symbol`,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]

# 接下来对表达矩阵exp进行探针ID转化，首先将其转换为数据框的形式，否则一会会报错
exp<-as.data.frame(exp)

# 此时虽然我们能看到exp的探针ID号（最左侧），但这个ID号并不是真正在列中，所以我们还需要将探针ID放到数据框里面，执行下面命令
exp$ID<-rownames(exp)	# 增加新的一列（最后一列），存放基因ID信息

gpl<-gpl[,c(1,3)]
# 接下来通过将gpl合并到exp中，从而在exp中将芯片ID号转化为Gene symbol
exp_symbol<-merge(exp,gpl,by="ID")

# 首先查看有多少重复的Gene symbol
table(duplicated(exp_symbol$`ID`))
# 输入以下命令进行去重操作同时将矩阵转换成标准的表达矩阵形式
# exp_unique<-avereps(exp_symbol[,-c(1,ncol(exp_symbol))],ID=exp_symbol$`ID`)

exp_symbol <- exp_symbol[,-1]

# 删除重复的
exp_symbol <- exp_symbol[!duplicated(exp_symbol$`Gene symbol`),]
exp_symbol <- na.omit(exp_symbol)
rownames(exp_symbol) <- exp_symbol$`Gene symbol`

exp_symbol <- exp_symbol[,-ncol(exp_symbol)]

GSE29621 <- as.data.frame(t(exp_symbol))
GSE29621$ID <- rownames(GSE29621)

surv <- cli[,c("geo_accession",
               "characteristics_ch1.2",
               "characteristics_ch1.3",
               "characteristics_ch1.4",
               "characteristics_ch1.6",
               "os event:ch1",
               "overall survival (os):ch1")]
colnames(surv) <- c("ID","T","N","M","AJCC","OS","OS.time")
table(surv$T)
surv$T <- ifelse(surv$T=="t stage: 2",2,
                 ifelse(surv$T=="t stage: 3",3,4))

save(surv, file = paste0(GSEID,"_stage.RData"))

tdiag <- surv[,c("ID","T")]
surv <- surv[,c("ID","OS.time","OS")]
surv$OS.time <- as.numeric(surv$OS.time)
surv$OS.time <- surv$OS.time*30
table(surv$OS)
surv$OS <- ifelse(surv$OS == "alive",0,1)

t_GSE29621 <- merge(tdiag,GSE29621,by = "ID")
colnames(t_GSE29621)[2] <- "T_stage"
save(t_GSE29621, file = paste0("t_",GSEID,".RData"))

GSE29621 <- merge(surv,GSE29621, by = "ID")

save(GSE29621,file = "GSE29621.RData")
write.csv(exp_symbol,paste0("dat_",GSEID,".csv"))

# ----------------------------------GSE72970------------------------------
rm(list = ls())
GSEID <- "GSE72970"
gset <- getGEO(filename = paste0(GSEID,"_series_matrix.txt.gz"),AnnotGPL = T, destdir = "./")
# gset <- getGEO("GSE106090",AnnotGPL = T, destdir = "./",getGPL = T)

gset[1]

# 提取表达量
exp <- exprs(gset)

# 临床信息
cli <- pData(gset)  #获取临床信息
save(cli,file = paste0(GSEID,"_cliall.RData"))

head(cli)

gpl<- fData(gset)

# 执行下面代码后，便通过只取第一个Gene symbol来达到去重效果，这时再打开gpl发现已经没有对应多个Gene symbol的情况了
gpl$`Gene symbol`<-data.frame(sapply(gpl$`Gene symbol`,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]

# 接下来对表达矩阵exp进行探针ID转化，首先将其转换为数据框的形式，否则一会会报错
exp<-as.data.frame(exp)

# 此时虽然我们能看到exp的探针ID号（最左侧），但这个ID号并不是真正在列中，所以我们还需要将探针ID放到数据框里面，执行下面命令
exp$ID<-rownames(exp)	# 增加新的一列（最后一列），存放基因ID信息

gpl<-gpl[,c(1,3)]
# 接下来通过将gpl合并到exp中，从而在exp中将芯片ID号转化为Gene symbol
exp_symbol<-merge(exp,gpl,by="ID")
exp_symbol <- na.omit(exp_symbol)
# 首先查看有多少重复的Gene symbol
table(duplicated(exp_symbol$`ID`))
# 输入以下命令进行去重操作同时将矩阵转换成标准的表达矩阵形式
# exp_unique<-avereps(exp_symbol[,-c(1,ncol(exp_symbol))],ID=exp_symbol$`ID`)

exp_symbol <- exp_symbol[,-1]

# 删除重复的
exp_symbol <- exp_symbol[!duplicated(exp_symbol$`Gene symbol`),]
exp_symbol <- na.omit(exp_symbol)
rownames(exp_symbol) <- exp_symbol$`Gene symbol`

exp_symbol <- exp_symbol[,-ncol(exp_symbol)]

GSE72970 <- exp_symbol
GSE72970 <- as.data.frame(t(GSE72970))
GSE72970$ID <- rownames(GSE72970)

surv <- cli[,c("geo_accession",
               "os:ch1",
               "os censored:ch1",
               "pn:ch1",
               "pt:ch1")]

colnames(surv) <- c("ID","OS.time","OS","N_stage","T_stage")
table(surv$T_stage)
surv <- surv[!surv$T_stage=="pTX",]
surv$T_stage <- factor(surv$T_stage,levels = c("pT1",
                                               "pT2",
                                               "pT3",
                                               "pT4"),
                       label = c(1,
                                 2,
                                 3,
                                 4))
surv$T_stage <- as.numeric(surv$T_stage)
surv$OS <- as.numeric(surv$OS)
surv$OS.time <- as.numeric(surv$OS.time)
range(surv$OS.time)
surv$OS.time <- surv$OS.time*30

save(surv, file = paste0(GSEID,"_stage.RData"))

tdiag <- surv[,c("ID","T_stage")]
surv <- surv[,c("ID","OS.time","OS")]

t_GSE72970 <- merge(tdiag,GSE72970,by = "ID")
save(t_GSE72970, file = paste0("t_",GSEID,".RData"))

GSE72970 <- merge(surv,GSE72970, by = "ID")

save(GSE72970,file = "GSE72970.RData")
write.csv(exp_symbol,paste0("dat_",GSEID,".csv"))

# ----------------------------------GSE69657------------------------------
rm(list = ls())
GSEID <- "GSE69657"
gset <- getGEO(filename = paste0(GSEID,"_series_matrix.txt.gz"),AnnotGPL = T, destdir = "./")
# gset <- getGEO("GSE106090",AnnotGPL = T, destdir = "./",getGPL = T)

gset[1]

# 提取表达量
exp <- exprs(gset)

# 临床信息
cli <- pData(gset)  #获取临床信息
save(cli,file = paste0(GSEID,"_cliall.RData"))

head(cli)

gpl<- fData(gset)

# 执行下面代码后，便通过只取第一个Gene symbol来达到去重效果，这时再打开gpl发现已经没有对应多个Gene symbol的情况了
gpl$`Gene symbol`<-data.frame(sapply(gpl$`Gene symbol`,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]

# 接下来对表达矩阵exp进行探针ID转化，首先将其转换为数据框的形式，否则一会会报错
exp<-as.data.frame(exp)

# 此时虽然我们能看到exp的探针ID号（最左侧），但这个ID号并不是真正在列中，所以我们还需要将探针ID放到数据框里面，执行下面命令
exp$ID<-rownames(exp)	# 增加新的一列（最后一列），存放基因ID信息

gpl<-gpl[,c(1,3)]
# 接下来通过将gpl合并到exp中，从而在exp中将芯片ID号转化为Gene symbol
exp_symbol<-merge(exp,gpl,by="ID")

# 首先查看有多少重复的Gene symbol
table(duplicated(exp_symbol$`ID`))
# 输入以下命令进行去重操作同时将矩阵转换成标准的表达矩阵形式
# exp_unique<-avereps(exp_symbol[,-c(1,ncol(exp_symbol))],ID=exp_symbol$`ID`)

exp_symbol <- exp_symbol[,-1]

# 删除重复的
exp_symbol <- exp_symbol[!duplicated(exp_symbol$`Gene symbol`),]
exp_symbol <- na.omit(exp_symbol)
rownames(exp_symbol) <- exp_symbol$`Gene symbol`

exp_symbol <- exp_symbol[,-ncol(exp_symbol)]

surv <- cli[,c("geo_accession",
               "tumor stage:ch1")]

surv$T_stage <- substr(surv$`tumor stage:ch1`, 1, 2)
surv$N_stage <- substr(surv$`tumor stage:ch1`, 3, 4)
surv$M_stage <- substr(surv$`tumor stage:ch1`, 5, 6)
surv$AJCC_stage <- substr(surv$`tumor stage:ch1`, 8, 9)
surv <- surv[,-2]
colnames(surv)[1] <- "ID"

table(surv$T_stage)
surv <- surv[!surv$T_stage=="pTX",]
surv$T_stage <- factor(surv$T_stage,levels = c("T1",
                                               "T2",
                                               "T3",
                                               "T4"),
                       label = c(1,
                                 2,
                                 3,
                                 4))
surv$T_stage <- as.numeric(surv$T_stage)
# surv$OS <- as.numeric(surv$OS)
# surv$OS.time <- as.numeric(surv$OS.time)
# range(surv$OS.time)
# surv$OS.time <- surv$OS.time*30

save(surv, file = paste0(GSEID,"_stage.RData"))

tdiag <- surv[,c("ID","T_stage")]
# surv <- surv[,c("ID","OS.time","OS")]

GSE69657 <- exp_symbol
GSE69657 <- as.data.frame(t(GSE69657))
GSE69657$ID <- rownames(GSE69657)
t_GSE69657 <- merge(tdiag,GSE69657,by = "ID")
save(t_GSE69657, file = "t_GSE69657.RData")

# =========================================================================
# TCGA
# =========================================================================
rm(list = ls())
foldpath <- paste0("D:/mydata/处理完成/CRC/")

if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)

load("tcga_crc.RData")
load("tcga_crc_cli.RData")

table(cliall$tstage)
cliall <- cliall[!cliall$tstage=="",]
cliall$tstage <- factor(cliall$tstage ,levels = c("T1","Tis",
                                               "T2",
                                               "T3",
                                               "T4","T4a","T4b"),
                       label = c(0,0,
                                 0,
                                 1,
                                 1,1,1))
cliall$tstage <- as.character(cliall$tstage)
cliall$tstage <- as.numeric(cliall$tstage)
table(cliall$tstage)

tdiag <- cliall[,c("ID","tstage")]
colnames(tdiag)[2] <- "T_stage"

TCGA_count <- merge(tdiag,crc_count,by = "ID")
TCGA_count <- TCGA_count[TCGA_count$group=="Tumor",]
TCGA_count <- TCGA_count[,-3]

TCGA_fpkm <- merge(tdiag,crc_fpkm,by = "ID")
TCGA_fpkm <- TCGA_fpkm[TCGA_fpkm$group=="Tumor",]
TCGA_fpkm <- TCGA_fpkm[,-3]

save(TCGA_count,TCGA_fpkm,file = "TCGA_tstage.RData")

# =========================================================================
# 整合
# =========================================================================
rm(list = ls())
dir ="D:/mydata/GEO/colorectal_cancer/"
# samples=list.files(dir,pattern = "^GSE.*RData$")

pattern <- "^t_GSE[0-9]{4,6}\\.RData"

# 列出所有匹配的文件
samples <- list.files(path = dir, pattern = pattern, full.names = FALSE)
# 
# index <- c()
# for (sample in samples){
#   if (nchar(sample)<17){
#     index <- c(index,sample)
#   }
# }

for (i in samples){
  load(paste0(dir,"/",i))
}

load("TCGA_tstage.RData")

table(t_GSE69657$T_stage)
t_GSE29621$T_stage <- ifelse(t_GSE29621$T_stage <3,0,1)
t_GSE69657$T_stage <- ifelse(t_GSE69657$T_stage <3,0,1)
t_GSE72970$T_stage <- ifelse(t_GSE72970$T_stage <3,0,1)

data <- list("TCGA"=TCGA_fpkm,
             "GSE29621" = t_GSE29621,
             "GSE69657" = t_GSE69657,
             "GSE72970" = t_GSE72970)

# # 获取所有数据框的列名
all_colnames <- lapply(data, colnames)

# 找到所有列名的交集
common_colnames <- Reduce(intersect, all_colnames)

# 使用交集列名统一所有数据框的列名
data <- lapply(data, function(df) {
  df <- df[common_colnames]
  setNames(df, common_colnames)
})

save(data, file = "t_data.RData")

for (i in 1:length(data)){
  dataname = names(data[i])
  print(dataname)
  print(table(data[[i]]$T_stage))
}
