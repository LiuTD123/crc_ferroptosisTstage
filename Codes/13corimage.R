options(timeout = Inf)

rm(list = ls())
foldpath <- "D:/workdir/23/13corimage"
if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)

library(readxl)
library(readxl)
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
# 有两个矩阵：genedata 和 imagedata
# genedata: 行名为样本，列名为基因，值为表达量
# imagedata: 行名为样本，列名为影像特征，值为特征值
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

# 初始化结果存储
cor_results <- list()
sig_features <- list()

genes <- read.csv("../04boruta/10.boruta_geneids.csv")$x
genedata <- genedata[,genes]
# 
corcut = 0.7
pcut = 0.05
# # 计算每个基因与每个影像特征的相关性
# for(gene in colnames(genedata)) {
#   # 提取当前基因的表达向量
#   gene_expr <- genedata[, gene]
#   
#   # 计算与所有影像特征的相关性
#   cor_values <- apply(imagedata, 2, function(img_feature) {
#     cor_test <- cor.test(gene_expr, img_feature, method = "pearson")
#     return(c(cor = cor_test$estimate, p.value = cor_test$p.value))
#   })
#   
#   # 转置结果
#   cor_df <- as.data.frame(t(cor_values))
#   colnames(cor_df) <- c("correlation", "p_value")
#   
#   # 存储结果
#   cor_results[[gene]] <- cor_df
#   
#   # 筛选满足条件的影像特征 (|r| > 0.6 且 p < 0.05)
#   sig_indices <- which(abs(cor_df$correlation) > corcut & cor_df$p_value < pcut)
#   if(length(sig_indices) > 0) {
#     sig_features[[gene]] <- rownames(cor_df)[sig_indices]
#   }
# }
# 
# # 创建相关性矩阵
# cor_matrix <- matrix(NA, nrow = ncol(genedata), ncol = ncol(imagedata),
#                      dimnames = list(colnames(genedata), colnames(imagedata)))
# pvalue_matrix <- matrix(NA, nrow = ncol(genedata), ncol = ncol(imagedata),
#                         dimnames = list(colnames(genedata), colnames(imagedata)))
# 
# # 填充相关性矩阵和p值矩阵
# for(gene in names(cor_results)) {
#   cor_matrix[gene, ] <- cor_results[[gene]]$correlation
#   pvalue_matrix[gene, ] <- cor_results[[gene]]$p_value
# }
# 
# # 输出结果
# print("基因与影像特征的相关性矩阵:")
# print(cor_matrix[1:5, 1:5])  # 只显示前5行和前5列
# 
# print("显著性p值矩阵:")
# print(pvalue_matrix[1:5, 1:5])  # 只显示前5行和前5列
# 
# print("与基因显著相关的影像特征 (|r| > 0.6, p < 0.05):")
# for(gene in names(sig_features)) {
#   cat(paste0("基因 ", gene, " 相关的影像特征: ", 
#              paste(sig_features[[gene]], collapse = ", "), "\n"))
# }

# ======================================
method <- "spearman"
cor_r <- cor(genedata,imagedata,method = method) 
#cor_p <- WGCNA::corPvalueStudent(cor_r,length(rownames(dat_exp_diff)))
library(psych)
d <- corr.test(genedata,imagedata,use="complete",method = method)
cor_p <- d$p
cor_r2 <- cor_r %>% as.data.frame %>% tibble::rownames_to_column(var = "Gene") %>%
  tidyr::gather(., ImageFeature,Correlation,-Gene)#转换数据长短 cor
cor_p2 <- cor_p %>% as.data.frame %>% tibble::rownames_to_column(var = "Gene") %>% 
  tidyr::gather(., ImageFeature, Pvalue, -Gene)#转换数据长短 p

cor_dat <- cbind(cor_r2, cor_p2)[,c("Gene","ImageFeature","Correlation","Pvalue")]
#cor_dat <- cor_dat[cor_dat$Cell %in% stat_res$Cell,]
write.csv(cor_dat,"04.correlation_cor.csv")

# 相关性热图带显著性----
data <- read.csv('04.correlation_cor.csv',row.names = 1,check.names = F)

data <- cor_dat
data <- data %>%
  mutate(text = case_when( #设置label，并加入判断，当P值符合特定条件就显示"\n"外加特定数量的*号
    Pvalue <= 0.001 ~ "***", #P<0.001就显示回车加三个星号
    between(Pvalue, 0.001, 0.01) ~ "**", #P为0.001-0.01 显示回车加两个*号
    between(Pvalue, 0.01, 0.05) ~ "*",  #P为0.01-0.05 显示回车加一个星号
    T ~ ""))
# data <- data %>%
#   mutate(text = case_when(  # 一定要 get 到 case_when() 函数奥秘
#     Pvalue > 0 ~ paste(round(Pvalue, 2), "+"), # round() 只保留两位小数
#     Pvalue < 0 ~ paste(round(Pvalue, 2), "-")))
library(ComplexHeatmap)
library(circlize)
library(dendextend)

col_fun <- colorRamp2(c(-1, 0, 1), c("#4DBBD5FF", "white", "#E64B35FF"))

ht <- Heatmap(
  cor_r_sig2,
  name = "Correlation",
  col = col_fun,
  rect_gp = gpar(col = "grey", lwd = 1),
  
  # 聚类设置
  cluster_rows = T,
  cluster_columns = F,
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  show_column_names = T,
  
  # 行列名称设置
  row_names_gp = gpar(fontsize = 12, fontface = "bold"),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  # row_names_max_width = max_text_width(rownames(sig_cor_matrix)),
  # column_names_max_height = max_text_width(colnames(sig_cor_matrix)),
  
  # 热图参数
  heatmap_legend_param = list(
    title = "Correlation",
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 10),
    legend_height = unit(4, "cm")
    # at = c(-1, -0.5, 0, 0.5, 1),
    # labels = c("-1.0", "-0.5", "0", "0.5", "1.0")
  )
)

pdf(paste0('correlation_Gene_ImageFeatures_',corcut,'2.pdf'), height = 3, width = 5)
ht
dev.off()

corcut = 0.7
pcut = 0.05

data_sig <- data[abs(data$Correlation) > corcut & data$Pvalue<pcut,]

sig_ImageFeatures <- unique(data_sig$ImageFeature)
save(sig_ImageFeatures, file = paste0("sig_ImageFeatures_",corcut,".RData"))

cor_r_sig <- cor_r[,unique(data_sig$ImageFeature)]

write.csv(cor_r_sig,file = "cor_r_sig.csv")

featurenames <- as.data.frame(colnames(cor_r_sig))
colnames(featurenames) <- "Name"
featurenames$Simp <- ""
for (i in 1:nrow(featurenames)){
  featurenames$Simp[i] <- paste0("F",{i})
}
rownames(featurenames) <- featurenames$Name
write.csv(featurenames,file = "featurenames.csv")

featurenames <- as.data.frame(t(featurenames))

cor_r_sig2 <- as.data.frame(cor_r_sig)
cor_r_sig2 <- rbind(cor_r_sig2,featurenames)

cor_r_sig2 <- cor_r_sig2[c(15,1:14),]
write.csv(cor_r_sig2,file = "cor_r_sig2.csv")

colnames(cor_r_sig2) <- cor_r_sig2[15,]
cor_r_sig2 <- cor_r_sig2[-c(14,15),]
for (i in 1:ncol(cor_r_sig2)){
  cor_r_sig2[[i]] <- as.numeric(cor_r_sig2[[i]])
}
cor_r_sig2 <- as.matrix(cor_r_sig2)
