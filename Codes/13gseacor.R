options(timeout = Inf)

rm(list = ls())
foldpath <- "D:/workdir/23/13corimage"
if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)

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
enableWGCNAThreads()

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
imagedata <- imagedata[imagedata$ID %in% copatients,]

rownames(genedata) <- genedata$ID
rownames(imagedata) <- imagedata$ID
genedata <- genedata[,-1]
imagedata <- imagedata[,-1]

# ========================================
# GSEA
# ========================================

# 计算ssGSEA-score -------

genes <- read.csv("../04boruta/10.boruta_geneids.csv")$x
gene_expr <- genedata[,genes]

gene_expr <- na.omit(gene_expr)
gene_expr <- as.data.frame(t(gene_expr))
gene_select <- rownames(gene_expr)
gene_list <- list(`gene_select` = unique(gene_select))
gene_list <- as.list(gene_list)

genedata <- as.data.frame(t(genedata[,-1]))
sp <- ssgseaParam(as.matrix(genedata),  gene_list)
ssgsea_score <-  gsva(sp) # method = "ssgsea"，"gsva"，"zscore"，"plage"无显著差异
# ssgsea_score_save <- rownames_to_column(as.data.frame(ssgsea_score),var = "DE_NMRGS")
# ssgsea_score_save <- t(as.data.frame(ssgsea_score))
# ssgsea_score_save <- ssgsea_score_save %>% as.data.frame()
# 差异箱线图----
score <- as.data.frame(t(ssgsea_score))
HL <- as.data.frame(score)

rownames(HL) <- gsub('\\.','\\-',rownames(HL))

# HL <- HL[train_group$sample,] %>% as.data.frame()
HL$type <- group[copatients,]$group
# rownames(HL) <- colnames(data)

HL <- HL[!duplicated(HL),]
HL$Group <- ifelse(HL$Group==1,"High","Low")
colnames(HL) <- c("Score","Group")
write.table(HL, file = "GSVA_Score.csv",sep = ",",row.names = T,col.names = NA,quote = F)

color_mapping <- c("High" = "#D20A13","Low" = "#58CDD9")

library(gghalves)
pdf("GSVA_score barplot.pdf", width = 7, height = 5)
p <- ggplot(HL, aes(x = Group, y = Score)) +
  stat_compare_means(aes(group = Group),
                     method = "anova",
                     label = "p.signif",
                     label.x = 1.5,
                     # label.y = 10,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  geom_point(aes(fill =Group, color = Group),
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = .1,
                                             dodge.width = 0),
             size = 3, shape = 20, alpha = 1)+
  geom_boxplot(position = position_nudge(x = 0.15),
               width = .1, alpha = 1, outlier.shape = NA)+
  geom_half_violin(aes(fill =Group, color = Group),position = position_nudge(x=0.2,y=0),
                   adjust = 0.8, trim = T, alpha = 1, colour = NA, side = "r")+
  scale_fill_manual(values = color_mapping) + # 指定填充颜色
  scale_color_manual(values = color_mapping) +   # 指定轮廓颜色
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_flip() + theme(legend.position = c(0.85, 0.3))+
  theme_bw()+
  theme(panel.grid=element_blank())+
  ggtitle("GSVA_score in PD vs Normal") +
  theme(plot.title    = element_text(color = 'black', size   = 16, hjust = 0.5),
        plot.subtitle = element_text(color = 'black', size   = 16,hjust = 0.5),
        plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
        axis.text.x   = element_text(color = 'black', size = 16, angle = 0),
        axis.text.y   = element_text(color = 'black', size = 16, angle = 0),
        axis.title.x  = element_text(color = 'black', size = 16, angle = 0),
        axis.title.y  = element_text(color = 'black', size = 16, angle = 90),
        legend.title  = element_text(color = 'black', size  = 16),
        legend.text   = element_text(color = 'black', size   = 16),
        axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
        axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
        panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA)) # 图四周框起来
p

dev.off()



