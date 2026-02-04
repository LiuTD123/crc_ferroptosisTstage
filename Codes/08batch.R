rm(list = ls())

foldpath <- paste0("D:/workdir/23/08batch")

if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(factoextra)
library(FactoMineR)
library(tidyverse)
library(scatterplot3d)
library(sva)
library(limma)

# 训练集合并数据-----
# --------------------------------substantia nigra------------------------------
load("D:/mydata/处理完成/CRC/t_data.RData")

for (i in 1:4){
  colnames(data[[i]])[2] <- "group"
}

groupall <- data.frame("ID" = character(0),
                       "group" = character(0),
                       "dataset" = character(0)
)
comdata <- data.frame()

normdata <- list()
for (i in 1:3){
  d <- data[[i]]
  
  rownames(d) <- d$ID
  group <- d[,c(1,2)]
  group$dataset <- names(data)[i]
  groupall <- rbind(groupall,
                    group)
  
  d <- d[,-c(1,2)]
  
  dimnames = list(rownames(d),
                  colnames(d))
  
  d <- matrix(as.numeric(as.matrix(d)),
              nrow = nrow(d),
              dimnames = dimnames)
  normdata <- normalizeBetweenArrays(d)
  
  # normdata <- as.data.frame(normdata)
  # normdata$ID <- rownames(normdata)
  # normdata <- merge(group,normdata)
  # normdata[[i]] <- normdata
  
  comdata <- rbind(comdata,
                   normdata)
}

# save(normdata, file = "normadata.RData")
# ----------保存------------
save(comdata,groupall, file = "comdata_pre.RData")

load("comdata_pre.RData")

# batchType <- c(rep(1,29),rep(2,35),rep(3,29))
batchType <- as.numeric(as.factor(groupall$dataset))

talgroup_list=factor(groupall$group,levels = c("Control","SLE"))

modType=groupall$group
mod = model.matrix(~modType)
# PH 去批次效应-----
# 去批次前
comexp <- as.data.frame(t(na.omit(comdata)))

table(groupall$dataset)
colors <- c(
  rep("#B0C4DE",65),
  rep("#FF69B4", 30),
  # rep("#FFB6C1",95), 
  rep("#FFA07A",633)
  # rep("#D8BFD8",306)
)


# 去批次后------------
# comexp_ComBat=ComBat(dat=comexp, batch=batchType, #使用ComBat法去批次
#                      mod=mod, par.prior=TRUE, ref.batch = 1)
comexp_ComBat=ComBat(dat=comexp, batch=groupall$dataset, #使用ComBat法去批次
                     par.prior=TRUE, ref.batch = "TCGA")

# ------2D PCA------------
# before
dat.pca <- PCA(as.data.frame(t(comexp)), graph = FALSE)
pca_plot <- fviz_pca_ind(dat.pca,
                         # axes = c(1, 2),
                         geom.ind = "points",#仅显示"点(points)"（但不是“文本(text)”）（show "points" only (but not "text")）
                         # col.ind = talgroup_list,
                         col.ind = groupall$dataset,
                         palette = c("#B0C4DE",
                                     "#FF69B4", 
                                     "#FFB6C1", 
                                     "#FFA07A"
                                     # "#ADD8E6"
                                     ),
                         addEllipses = TRUE, 
                         legend.title = "Groups")
pca_plot
ggsave(plot = pca_plot,filename ="PCA_befor.pdf", height = 5, width = 7)

# after
# PCA
batchType <- as.factor(batchType)

dat.pca2 <- PCA(as.data.frame(t(comexp_ComBat)), graph = FALSE)
pca_plot2 <- fviz_pca_ind(dat.pca2,
                          geom.ind = "points",#仅显示"点(points)"（但不是“文本(text)”）（show "points" only (but not "text")）
                          # col.ind = talgroup_list,
                          col.ind = groupall$dataset,
                          palette = c("#B0C4DE",
                                      "#FF69B4", 
                                      "#FFB6C1",
                                      "#FFA07A"
                                      # "#FFFFE0", 
                                      # "#ADD8E6"
                                      ),
                          addEllipses = TRUE, 
                          legend.title = "Groups")
pca_plot2
ggsave(plot = pca_plot2,filename ="PCA_after.pdf", height = 5, width = 7)

# -----------------------
load("comdata.RData")

data <- as.data.frame(t(comexp_ComBat))
data$ID <- rownames(data)
data_all <- merge(groupall,data,by = "ID")

table(data_all$dataset)

table(data_all$group)
# data_all$group <- ifelse(data_all$group == "DKD",1,0)

data = list("TCGA" = data_all[data_all$dataset == "TCGA",],
            "GSE29621" = data_all[data_all$dataset == "GSE29621",],
            "GSE69657" = data_all[data_all$dataset == "GSE69657",],
            "GSE72970" = data_all[data_all$dataset == "GSE72970",])

data <- lapply(data, 
               function(x) {
                 subset(x, select = -3)
               })

group <- groupall
save(data,group,
     file = "data_all.RData")

