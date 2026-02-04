options(timeout = Inf)

rm(list = ls())
foldpath <- "D:/workdir/23/14wgcna_image"
if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)

enableWGCNAThreads()

library(GSVA)
library(magrittr)
library(tidyr)
library(ggpubr)
library(tibble)
library(dplyr)
library(ggplot2)
library(limma)
library(WGCNA)
library(stringr)
library(rstatix)

load("../08batch/data_all.RData")
genedata <- data[[1]]
genedata$ID <- substr(genedata$ID,1,12)

imagedata <- read_excel("D:/radiomicworkdir/CRC/TCGA.xlsx")
colnames(imagedata)[1] <- "ID"
copatients <- intersect(genedata$ID,imagedata$ID)
genedata <- genedata[genedata$ID %in% copatients,]
imagedata <- as.data.frame(imagedata[imagedata$ID %in% copatients,])

rownames(genedata) <- genedata$ID
rownames(imagedata) <- imagedata$ID
genedata <- genedata[,-1]
imagedata <- imagedata[,-1]
# 
# genes <- read.csv("../04boruta/10.boruta_geneids.csv")$x
# data <- as.data.frame(t(genedata[,genes]))

data <- as.data.frame(t(imagedata))
# data格式：行名为基因，列名为样本

# step 1 ：输入数据的准备 ---------
RNAseq_voom   <- data
datExpr <- t(RNAseq_voom[order(apply(RNAseq_voom,1,mad), decreasing = T),])
dim(datExpr)
method <- "ward.D"
datExpr_tree <- hclust(dist(datExpr),method = method)

pdf("sample_clusting.pdf")
plot(datExpr_tree,main = "Sample clusting",sub = "",
     xlab = "", cex.lab = 2,
     labels=FALSE,
     cex.axis = 1,cex.main = 1,cex.lab =1)

dev.off()

pdf(file = paste0("00.clust_ALL.pdf"),width = 16,height = 5)
a <- dev.cur()
png(file = paste0("00.clust_ALL.png"),width= 16, height= 5, units="in", res=300)
dev.control("enable")
par(mar = c(0,5,2,0))
plot(datExpr_tree,main = "Sample clusting",sub = "",
     labels=FALSE,
     xlab = "", cex.lab = 2,
     cex.axis = 1,cex.main = 1,cex.lab =1)

dev.copy(which = a)
dev.off()
dev.off()

#####-------------------- 没有离群值就不要cut
clust = cutreeStatic(datExpr_tree, cutHeight = 1.5e12, minSize = 10)
table(clust)
# clust
# 1  2  3  4  5
# 99 33 18 13 10
keepSamples = (clust==1)
datExpr = datExpr[keepSamples,]

# step 2 ：确定最佳ß -------
if(T){
  cutline <- 0.5
  powers = c(1:20)
  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,RsquaredCut = cutline)
  sft$powerEstimate# 最佳软阈值
  #同时输出PDF与PNG
  pdf(file = paste0("01.Threshold.pdf"),width = 12,height = 7)
  a <- dev.cur()   #记录pdf设备
  png(file = paste0("01.Threshold.png"),width=12, height=7, units="in", res=300)
  dev.control("enable")
  # Plot the results:
  par(mfrow = c(1,2))
  cex1 = 1
  
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  abline(h = cutline,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(
    sft$fitIndices[,1], sft$fitIndices[,5],
    xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
    main = paste("Mean connectivity")
  )
  
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.copy(which = a)  #复制来自png设备的图片到pdf
  dev.off()
  dev.off()
  sft$powerEstimate
}

power <- 10
# 软阈值powerEstimate 为0.85：7
save.image("sft_0.5.Rdata")

# load("sft.Rdata")
# datExpr0 <- lapply(datExpr, as.numeric)
datExpr0 <- datExpr %>% as.data.frame()
datExpr0 <- as.data.frame(lapply(datExpr0, as.numeric))

# step 3 ：构建加权共表达网络 --------

if(T){
  cor <- WGCNA::cor
  corType <-  "pearson"
  net = blockwiseModules(
    datExpr0,
    power = power,
    maxBlockSize = ncol(datExpr0),
    corType =  corType,
    deepSplit = 2,
    minModuleSize = 30,
    mergeCutHeight = 0.3,
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs = F,
    verbose = 5 #整数级别的详细程度。
  )
  cor<-stats::cor
  table(net$colors)
}

save.image("sft_0.5.Rdata")

# 0   1   2 
# 854 131  66

# load("net.Rdata")
# step 4 ：模块可视化  --------
library(RColorBrewer)
mergedColors = labels2colors(net$colors)
table(mergedColors)
moduleColors <- mergedColors
pdf(file = paste0("02.Cluster_Dendrogram.pdf"),width = 12,height = 6)
a <- dev.cur()
png(file = paste0("02.Cluster_Dendrogram.png"),width=12, height= 6, units="in", res=300)
dev.control("enable")
plotDendroAndColors(net$dendrograms[[1]],
                    mergedColors[net$blockGenes[[1]]],
                    "Module",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    cex.rowText = 2,
                    cex.colorLabels = 1.5, cex.dendroLabels = 1.5,
                    addGuide = TRUE,
                    guideHang = 0.2)
dev.copy(which = a)  #复制来自png设备的图片到pdf
dev.off()
dev.off()
# step 5 ：模块和性状的关系 -------

ssgsea <- read.csv("../13corimage/GSVA_Score.csv")
rownames(ssgsea) <- ssgsea$X
ssgsea <- ssgsea[,-1]
# train_group <- train_group[!duplicated(train_group),]
# rownames(train_group)<- train_group$sample

if(T){
  moduleColors <- labels2colors(net$colors)
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  # MEs0 <- MEs0[,-which(colnames(MEs0) == "MEgrey")]
  MEs = orderMEs(MEs0)
  # save(MEs,file = "MEs.Rdata")
  # design <- train_group[match(rownames(MEs),rownames(train_group)),,drop = F]
  design <- ssgsea[match(rownames(MEs),rownames(ssgsea)),,drop = F]
  
  # design <- ifelse(design$group == "Keloid", 1,0) %>% as.data.frame()
  design$Group <- ifelse(design$Group == "Tumor", 1,0)
  
  if (corType=="pearson") {
    moduleTraitCor = cor(MEs, design$Score, use = "p")
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(design)) %>% signif(3)
    modTraitPadj = p.adjust(moduleTraitPvalue, method = "BH")
  } else {
    modTraitCorP = bicorAndPvalue(MEs, design, robustY=F)
    moduleTraitCor = modTraitCorP$bicor
    moduleTraitPvalue   = modTraitCorP$p
    modTraitPadj = p.adjust(moduleTraitPvalue, method = "BH")
  }
  sizeGrWindow(10,6)
  textMatrix = paste(signif(moduleTraitCor, 4), "\n(",
                     signif(modTraitPadj, 3), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  
  pdf(file = paste0("03.Module-relation.pdf"),width = 8,height = 14)
  a <- dev.cur()   #记录pdf设备
  png(file = paste0("03.Module-relation.png"),width = 8, height=14, units="in", res=300)
  dev.control("enable")
  par(mar = c(5,15,3,5));#项的作用是调整绘图区域距离外围框线的距离。下、左、上、右
  labeledHeatmap(
    Matrix = moduleTraitCor,
    #test
    xLabels = c("Gene"),
    # xLabels = c("H_HCC","L_HCC"),
    xLabelsPosition = "bottom",
    naColor = "grey",
    yLabels = names(MEs),
    ySymbols = names(MEs),
    colorLabels = FALSE,
    colors = blueWhiteRed(50),
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 1.3,
    cex.lab = 2,
    zlim = c(-1,1),
    xLabelsAngle = 0,
    xLabelsAdj = 0.5,
    main = paste("Module-trait relationships"))
  dev.copy(which = a)
  dev.off()
  dev.off()
}

# save.image("modulTrait_0.85.Rdata")
# step 6 ：提取指定模块的基因并绘制热图 ------------s

library(tibble)
moduleTraitCor <- tibble::rownames_to_column(as.data.frame(moduleTraitCor))
moduleTraitPvalue <- tibble::rownames_to_column(as.data.frame(moduleTraitPvalue))
moduleTrait <- merge(moduleTraitCor,moduleTraitPvalue,by.x = 1,by.y = 1,keep_all = T)
colnames(moduleTrait) <- c('module','scores.cor','scores.p')

# 提取p<0.05，同时|cor|>0.5
moduleTrait_sig <- moduleTrait[which(moduleTrait$scores.p < 0.05 & abs(moduleTrait$scores.cor) > 0.5),]
moduleTrait_sig
module_name1 <- grep(max(abs(moduleTrait_sig$scores.cor)),moduleTrait_sig$scores.cor) %>%
  moduleTrait_sig[.,1] %>% str_sub(.,3,15)
module_name1
# 模块颜色：0.8:black; 0.85:red; 0.9: yellow
# 直接提取模块内的基因 ------
if(T){
  # Select module
  module1 = "blue"
  # Select module probes
  probes = colnames(datExpr)
  
  inModule1 = (moduleColors %in% c(module1))
  modProbes1 = probes[inModule1]
  modProbes1 <- c(modProbes1) %>% as.data.frame() %>% distinct()
  colnames(modProbes1) <- 'symbol'
  write.csv(modProbes1,"04.WGCNA_genes.csv",row.names = F)
}
dim(modProbes1)
# [1] 424   1

save.image("modProbes1.rds")
#

