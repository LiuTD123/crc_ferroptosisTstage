rm(list = ls())

foldpath <- paste0("D:/workdir/23/03venn")

if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)
#韦恩图（VennDiagram 包，适用样本数 2-5）
library(VennDiagram)
library(ggvenn)
library(readr)
library(readxl)
library(tidyverse)

driverg <- read.csv("../FerroptosisGenes/ferroptosis_driver.csv")$symbol
markerg <- read.csv("../FerroptosisGenes/ferroptosis_marker.csv")$symbol
suppressorg <- read.csv("../FerroptosisGenes/ferroptosis_suppressor.csv")$symbol
unclassifiedg <- read.csv("../FerroptosisGenes/ferroptosis_unclassified.csv")$symbol

ferroptosis_genes <- c(driverg,
                       markerg,
                       suppressorg,
                       unclassifiedg)

# ================================================================
degs <- readRDS("../02degs/tcga_filter_diff2.rds")

inter <- Reduce(intersect,list(rownames(degs),ferroptosis_genes))
# inter <- Reduce(intersect,list(rownames(DEGs1),DEGs2$X,cluster$genesymbol))

write.csv(inter,"ferroptosis_genes.csv",row.names = T)
save(inter, file = "inter.RData")
vennlist <- list("DEGs"=rownames(degs),"Ferroptosis"=ferroptosis_genes)

pdf('04.venn.pdf',w=6,h=8)
ggvenn(vennlist,
       c("DEGs","Ferroptosis"),
       fill_color = c("#58CDD9","#FFD121"),
       show_percentage = T,
       fill_alpha=0.5,
       stroke_alpha=1,
       stroke_size=0.4,
       text_size=5,
       stroke_color=NA,
       stroke_linetype="solid",
       set_name_color=c("#58CDD9","#FFD121"),
       set_name_size=6,
       text_color="black")
dev.off()
