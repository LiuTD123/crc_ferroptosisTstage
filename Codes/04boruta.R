rm(list = ls())
foldpath <- paste0("D:/workdir/23/04boruta")

if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)
set.seed(125)
library(glmnet)
library(dplyr)
library(ggplot2)
library(tidyverse)
# ---------------------------------------机器学习部分
# 
library(xgboost)
library(randomForest)
library(caret)
library(DALEX)
library(stats)
library(e1071)
#  # # ------------------Boruta -------

library(Boruta)
# load("D:/mydata/处理完成/CRC/t_data.RData")
load("../08batch/data_all.RData")
load("../03venn/inter.RData")

for (i in 1:4){
  colnames(data[[i]])[2] <- "group"
}

dat <- data[[1]]
dat <- dat[order(dat$group,decreasing = F),]

rownames(dat) <- dat$ID
dat1 <- dat[-1]

candigenes <- inter
dat1 <- dat1[,c("group",candigenes)]

dat1$group <- as.factor(dat1$group)

# num <- length(colnames(dat1))-1
# dat1[,1:num] <- as.data.frame(lapply(dat1[,1:num],as.numeric))
dat1 <- na.omit(dat1)

res.Boruta<-Boruta(x=dat1[,2: ncol(dat1)], y=as.factor(dat1[,1]), pValue=0.01, mcAdj=T,
                   maxRuns=1000)
Boruta<-attStats(res.Boruta) #给出Boruta算法的结果
write.csv(Boruta,"09.Boruta.csv") #,head = T
table(res.Boruta$finalDecision)
# Tentative Confirmed  Rejected 
# 7        59        28 

Boruta <- Boruta[order(Boruta$medianImp, decreasing = T),]
boruta_geneids<-Boruta[Boruta$decision=='Confirmed',]%>%rownames(.)
boruta_geneids

##定义一个函数提取每个变量对应的重要性值。

library(dplyr)
boruta.imp <- function(x){
  imp <- reshape2::melt(x$ImpHistory, na.rm=T)[,-1]
  colnames(imp) <- c("Variable","Importance")
  imp <- imp[is.finite(imp$Importance),]
  
  variableGrp <- data.frame(Variable=names(x$finalDecision),
                            finalDecision=x$finalDecision)
  
  showGrp <- data.frame(Variable=c("shadowMax", "shadowMean", "shadowMin"),
                        finalDecision=c("shadowMax", "shadowMean", "shadowMin"))
  
  variableGrp <- rbind(variableGrp, showGrp)
  
  boruta.variable.imp <- merge(imp, variableGrp, all.x=T)
  
  # sortedVariable <- boruta.variable.imp %>% group_by(Variable) %>%
  #   summarise(median=median(Importance)) %>% arrange(median)
  # sortedVariable <- as.vector(sortedVariable$Variable)
  # 
  # 
  # boruta.variable.imp$Variable <- factor(boruta.variable.imp$Variable, levels=sortedVariable)
  
  invisible(boruta.variable.imp)
}

boruta.variable.imp <- boruta.imp(res.Boruta)

head(boruta.variable.imp)
boruta.variable.imp <- boruta.variable.imp[]

# 绘制Boruta算法运行过程中各个变量的重要性得分的变化 （绿色是重要的变量，红色是不重要的变量，蓝色是影子变量，黄色是Tentative变量）
# 用来查看是否有必要增加迭代的次数以便再次确认Tentative变量中是否有一部分为有意义的特征变量。
# ，黄色变量部分随着迭代还是有部分可能高于最高值，可以继续尝试增加迭代次数。

pdf("boruta_imphistory.pdf",w=5,h=4)
Boruta::plotImpHistory(res.Boruta)
dev.off()

# 提取重要变量
boruta.finalVars <- data.frame(Item=getSelectedAttributes(res.Boruta, withTentative = F), Type="Boruta")

# devtools::install_github("Tong-Chen/YSX")
# library(YSX)
# devtools::install_github("Tong-Chen/ImageGP")
library(ImageGP)

p <- sp_boxplot(boruta.variable.imp, melted=T, xvariable = "Variable", yvariable = "Importance",
                legend_variable = "finalDecision", legend_variable_order = c("Tentative", "Confirmed", "Rejected", "shadowMax", "shadowMean", "shadowMin"),
                xtics_angle = 90,
                xvariable_order = c(rownames(Boruta)))

# ggsave(p,file = "09.Boruta.pdf",w = 18, h = 8,family='Times')
# ggsave(p,file = "09.Boruta.png",w = 18, h = 8,family='Times')
pdf("09.Boruta.pdf",w = 4, h = 3,family='Times')
p
dev.off()

write.table(boruta_geneids,"10.boruta_geneids.csv",sep = "\t",row.names = F,col.names = T,quote = F)
save(boruta_geneids, file = "borutagenes.RData")
