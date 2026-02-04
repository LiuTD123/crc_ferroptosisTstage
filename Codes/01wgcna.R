rm(list = ls())

foldpath <- paste0("D:/workdir/23/01wgcna")

if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)
# options(stringsAsFactors = F)
library(dplyr)
library(ggplot2)
library(limma)
library(WGCNA)
library(stringr)
library(rstatix)
library(GSVA)
library(magrittr)
library(tidyr)
library(ggpubr)
library(tibble)

enableWGCNAThreads()

# ---------------------------------------------WGCNA--------------------
load("D:/mydata/处理完成/CRC/t_data.RData")
load("D:/mydata/处理完成/CRC/TCGA_tstage.RData")
# load("D:/mydata/处理完成/CRC/tcga_crc.RData")
load("D:/mydata/处理完成/CRC/tcga_crc_cli.RData")
cliall <- cliall[!cliall$tstage=="",]
cliall$tstage <- factor(cliall$tstage ,levels = c("T1","Tis",
                                                  "T2",
                                                  "T3",
                                                  "T4","T4a","T4b"),
                        label = c("T1/0","T1/0",
                                  "T2",
                                  "T3",
                                  "T4","T4","T4"))
table(cliall$tstage)
group <- cliall[,c("ID","tstage")]
colnames(group)[2]<-"group"
group$group <- as.character(group$group)
group <- as.data.frame(na.omit(group))

dataExpr <- data[[1]]
dataExpr <- TCGA_count
group <- group[group$ID %in% dataExpr$ID,]

group <- TCGA_count[,c(1,2)]
group <- dataExpr[,c(1,2)]
colnames(group)[2]<- "group"
group$group <- ifelse(group$group == 1, "High","Low")

rownames(dataExpr) <- dataExpr$ID
dataExpr <- as.data.frame(dataExpr[,-c(1,2)])

dataExpr <- as.data.frame(t(dataExpr[,-c(1,2)]))
dataExpr <- na.omit(dataExpr)
# data格式：行名为基因，列名为样本
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
# 对二元变量，如样本性状信息计算相关性时，
# 或基因表达严重依赖于疾病状态时，需设置下面参数
maxPOutliers = ifelse(corType=="pearson",1,0.05)
# 关联样品性状的二元变量时，设置
robustY = ifelse(corType=="pearson",T,F)
type = "unsigned"

# 数据筛选-----------
## 筛选中位绝对偏差前75%的基因，至少MAD大于0.01
## 筛选后会降低运算量，也会失去部分信息
## 也可不做筛选，使MAD大于0即可
m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad > 
                                max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]

## 转换为样品在行，基因在列的矩阵
dataExpr <- as.data.frame(t(dataExprVar))

## 检测缺失值
gsg = goodSamplesGenes(dataExpr, verbose = 3)

##  Flagging genes and samples with too many missing values...
##   ..step 1

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)

dim(dataExpr)

# "#ffb2b2","#b2e7cb","#7570B3","orange","#e67e22"
## 查看是否有离群样品------------------
datExpr_tree = hclust(dist(dataExpr), method = "ward.D2")
# plot(datExpr_tree, main = "Sample clustering to detect outliers", sub="", xlab="")

tcga_group <- group

colors <- ifelse(tcga_group$group=="T1/0","#029149",ifelse(
  tcga_group$group=="T2","#6E568C",ifelse(
    tcga_group$group=="T3","#E0367A","#D8D155"
  )
))

pdf(file="2.Sample_dendrogram_and_trait_heatmap.pdf",width=20,height=12) 
par(mai=c(1,1,1,1))
plotDendroAndColors(datExpr_tree, colors = ifelse(tcga_group$group=="High","#E0367A","#D8D155"), 
                    groupLabels = names(tcga_group),
                    dendroLabels =F, #不显示样本标签
                    main = "Sample dendrogram and trait heatmap",cex.colorLabels = 1.5, cex.dendroLabels = 1, cex.rowText = 2) 
legend(x= -0.15,y = 0.6,xpd = TRUE,
       c("High\n","Low\n"),
       col=c("#E0367A","#D8D155"),
       lty=1, lwd=2, cex = 0.7)
dev.off() 
# "#EB4B17", "#2775AB", '#4C8045',"#D8D155"
# --------软阈值-----------
# 软阈值的筛选原则是使构建的网络更符合无标度网络特征。
# powers = c(c(1:10), seq(from = 12, to=30, by=2))
powers = c(c(1:20))
cutvalue <- 0.75
sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                        RsquaredCut = cutvalue,
                        networkType=type, verbose=5)
sft$powerEstimate

# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
# 网络越符合无标度特征 (non-scale)
pdf(file = paste0("01.Threshold.pdf"),width = 12,height = 7)
par(mfrow = c(1,2))
cex1 = 1
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# 筛选标准。R-square=0.85
abline(h=cutvalue,col="red")

# Soft threshold与平均连通性
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")
power = sft$powerEstimate
power <- 9
# 16
dev.off()

if(T){
  cor <- WGCNA::cor
  corType <-  "pearson"
  exprMat = "TOMfilebase"
  net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                         TOMType = type, minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs=TRUE, corType = corType, 
                         maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                         saveTOMFileBase = paste0(exprMat, ".tom"),
                         verbose = 3)
  cor<-stats::cor
  table(net$colors)
}

# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22   23   24   25   26 
# 5559 4740 1584 1308 1268  892  699  652  612  395  337  321  149  147  141  136  111  110  103   71   63   59   59   48   39   38   38 
# 27 
# 33
save.image("net.Rdata")

load("net.Rdata")
# ---层级聚类树展示模块---------
## 灰色的为**未分类**到模块的基因。
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
# Plot the dendrogram and the module colors underneath
# 如果对结果不满意，还可以recutBlockwiseTrees，节省计算时间
pdf(file = "02.Cluster_Dendrogram.pdf",width = 12,height = 6)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# 绘制模块之间相关性图----------
# module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
MEs = net$MEs

### 不需要重新计算，改下列名字就好
### 官方教程是重新计算的，起始可以不用这么麻烦
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)

# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距
pdf(file = "eigene adjacency heatmap.pdf",w = 5, h = 8)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap",
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T,
                      xLabelsAngle = 90)
dev.off()

## 如果有表型数据，也可以跟ME数据放一起，一起出图
MEs_colpheno = orderMEs(cbind(MEs_col, factor(group$group)))
plotEigengeneNetworks(MEs_colpheno, "Eigengene adjacency heatmap",
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T,
                      xLabelsAngle = 90)

# 可视化基因网络 (TOM plot)----------------
# 如果采用分步计算，或设置的blocksize>=总基因数，直接load计算好的TOM结果
# 否则需要再计算一遍，比较耗费时间
# TOM = TOMsimilarityFromExpr(dataExpr, power=power, corType=corType, networkType=type)
load(net$TOMFiles[1], verbose=T)

## Loading objects:
##   TOM

TOM <- as.matrix(TOM)

dissTOM = 1-TOM
# Transform dissTOM with a power to make moderately strong 
# connections more visible in the heatmap
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function

# 这一部分特别耗时，行列同时做层级聚类
pdf(file = "cluster_dendrograms",w= 8, h= 8)
TOMplot(plotTOM, net$dendrograms, moduleColors, 
        main = "Network heatmap plot, all genes")
dev.off()

# 导出网络用于Cytoscape-----------
probes = colnames(dataExpr)
dimnames(TOM) <- list(probes, probes)

# Export the network into edge and node list files Cytoscape can read
# threshold 默认为0.5, 可以根据自己的需要调整，也可以都导出后在
# cytoscape中再调整
cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = paste(exprMat, ".edges.txt", sep=""),
                               nodeFile = paste(exprMat, ".nodes.txt", sep=""),
                               weighted = TRUE, threshold = 0,
                               nodeNames = probes, nodeAttr = moduleColors)

# step 5 ：模块和性状的关系 -------
# ssgsea <- read.csv("GSVA_Score.csv")
# rownames(ssgsea) <- ssgsea$X
# ssgsea <- ssgsea[,-1]
train_group <- tcga_group
train_group <- group

rownames(train_group)<- train_group$ID

dataExpr <- dataExpr[train_group$ID,]
if(T){
  moduleColors <- labels2colors(net$colors)
  MEs0 = moduleEigengenes(dataExpr, moduleColors)$eigengenes
  # MEs0 <- MEs0[,-which(colnames(MEs0) == "MEgrey")]
  MEs = orderMEs(MEs0)
  # save(MEs,file = "MEs.Rdata")
  # design <- train_group[match(rownames(MEs),rownames(train_group)),,drop = F]
  # design <- ssgsea[match(rownames(MEs),rownames(ssgsea)),,drop = F]
  design <- model.matrix(~0+factor(train_group$group))
  # colnames(design) <- c("T1/0","T2","T3","T4")
  colnames(design) <- c("High","Low")
  design <- design[,c(2,1)]
  # design$group <- ifelse(design$group == "PD", 1,0)
  nGenes = ncol(dataExpr)
  nSamples = nrow(dataExpr)
  
  if (corType=="pearson") {
    moduleTraitCor = cor(MEs, design, use = "p")
    # moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(design)) %>% signif(3)
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples) %>% signif(3)
    modTraitPadj = p.adjust(moduleTraitPvalue, method = "BH")
  } else {
    modTraitCorP = bicorAndPvalue(MEs, design, robustY=F)
    moduleTraitCor = modTraitCorP$bicor
    moduleTraitPvalue   = modTraitCorP$p
    modTraitPadj = p.adjust(moduleTraitPvalue, method = "BH")
  }
  sizeGrWindow(10,6)
  textMatrix = paste(signif(moduleTraitCor, 4), "\n(",
                     signif(moduleTraitPvalue, 3), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  
  pdf(file = paste0("03.Module-relation.pdf"),width = 8,height = 17)
  par(mar = c(5,15,3,5));#项的作用是调整绘图区域距离外围框线的距离。下、左、上、右
  labeledHeatmap(
    Matrix = moduleTraitCor,
    #test
    # xLabels = c("Gene"),
    xLabels = colnames(design),
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
  dev.off()
}

save.image("modulTrait.Rdata")
# step 6 ：提取指定模块的基因并绘制热图 ------------
load("modulTrait.Rdata")

library(tibble)
moduleTraitCor <- tibble::rownames_to_column(as.data.frame(moduleTraitCor))
moduleTraitPvalue <- tibble::rownames_to_column(as.data.frame(moduleTraitPvalue))
moduleTrait <- merge(moduleTraitCor,moduleTraitPvalue,by.x = 1,by.y = 1,keep_all = T)
colnames(moduleTrait) <- c('module','scores.cor','scores.p')

# 提取p<0.05，同时|cor|>0.5
moduleTrait_sig <- moduleTrait[which(moduleTrait$scores.p < 0.05 & abs(moduleTrait$scores.cor) > 0.3),]
moduleTrait_sig
module_name1 <- grep(max(abs(moduleTrait_sig$scores.cor)),moduleTrait_sig$scores.cor) %>%
  moduleTrait_sig[.,1] %>% str_sub(.,3,15)
module_name1
# 模块颜色：0.8:black; 0.85:red; 0.9: yellow
# 直接提取模块内的基因 ------

# 选择一个最正相关和一个最负相关
module_name_p4 <- "purple"
module_name_p3 <- "midnightblue"
module_name_p2 <- "blue"
module_name_p1 <- "brown"

if(T){
  probes = colnames(dataExpr)
  
  inModule1 = (moduleColors %in% c(module_name_p1,
                                   module_name_p2,
                                   module_name_p3,
                                   module_name_p4))
  modProbes1 = probes[inModule1]
  modProbes1 <- c(modProbes1) %>% as.data.frame() %>% distinct()
  colnames(modProbes1) <- 'symbol'
  # write.csv(modProbes1,"04.WGCNA_genes.csv",row.names = F)
  wgcna_genes <- modProbes1$symbol
  save(wgcna_genes, file ="wgcna_genes.RData")
}
length(wgcna_genes)
# [1] 2795

save.image("modProbes1.RData")
#

