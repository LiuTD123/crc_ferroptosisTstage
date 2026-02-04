options(timeout = Inf)

rm(list = ls())
foldpath <- "D:/workdir/23//10gsea"
if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)
#  09_GSEA----
library(enrichplot)
library(psych)
# BiocManager::install("GSEABase",force = TRUE)
library(GSEABase)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(stringr)
library(ggplot2)
require(DOSE)
library(Hmisc)
library(future)
library(future.apply)
library(limma)
library(psych)
library(RColorBrewer)
# BiocManager::install("GSVA")
library(GSVA)
# devtools::install_github("junjunlab/GseaVis",force = TRUE)
library(GseaVis)
# devtools::install_github("satijalab/seurat-data")

options(timeout = 9999999)

# KEGG_df = msigdbr(species = "Homo sapiens",collection = "C2")%>% 
#   dplyr::select(gs_description,gene_symbol)

# subcollection	
# Sub-collection abbreviation, such as "CGP" or "BP". Use msigdbr_collections() for the available options.
print(n=25,msigdbr_collections())

KEGG_df = msigdbr(species = "Homo sapiens",category = "C2", subcategory ="CP")%>% 
  dplyr::select(gs_description,gene_symbol)
head(KEGG_df)
# 
KEGG_df = msigdbr(species = "Homo sapiens",subcollection ="IMMUNESIGDB")%>%
  dplyr::select(gs_description,gene_symbol)
# head(KEGG_df)

# Human Phenotype Ontology
KEGG_df = msigdbr(species = "Homo sapiens",subcollection ="HPO")%>%
  dplyr::select(gs_description,gene_symbol)

genes <- c("SLC2A3",
           "PTGS2",
           "ANGPTL7",
           "FABP4",
           "RGS4")
hub_gene <- genes

load("../08batch/data_all.RData")
exp <- data[[1]]
rownames(exp) <- exp$ID
exp <- exp[,-c(1:2)]
exp <- as.data.frame(t(exp))
genes2 <- hub_gene

# ii <- "A2M"
# # ii <- "CFH"
# # tar.exp <- t(exp[ii,])
for (ii in genes2){
  tar.exp <- t(exp[ii,])
  cor_df <- cor(t(exp), tar.exp, method="pearson")
  # print(head(cor_df))
  data1 = data.frame(gene=rownames(cor_df),cor=cor_df[, 1])
  # print(head(data1))
  ge = data1$cor
  names(ge) = data1$gene
  ge = sort(ge,decreasing = T)
  print(head(ge))
  em_kegg <- GSEA(ge, TERM2GENE = KEGG_df,pAdjustMethod = "none")
  print(head(em_kegg))

  p1 <- GseaVis::gseaNb(object = em_kegg,
                        # rank.gene = ii,
                        # rank.gene.nudgey= 8,#the gene label nudge y on rank plot, defalut is 2.
                        # addGene = T, #whether add gene name on the curve, defalut is FALSE.
                        # geneSize = 6,# gene label text size, defalut is 4.
                        geneSetID = em_kegg@result$ID[1:8],
                        # legend.position = c(x = 1,y = 1),
                        addPval = F,
                        pvalX = 1,pvalY = 1,
                        ght.facet = T,
                        termWidth = 45, #the width or the term name, defalut is 40.
                        curveCol = brewer.pal(8,'Paired'))+
    labs(title = ii)+
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14,vjust=65),
          plot.subtitle = element_text(hjust = 0.5, color = "grey50"),
          plot.title.position = "plot")

  p1
  write.csv(em_kegg ,file=paste0("GSEA",ii, "_KEGG_gsea.csv"))
  ggsave(paste0("GSEA",ii,'.KEGG.pdf'),width = 10,height = 5, p1) #GSEA分析
  ggsave(paste0("GSEA",ii,'.KEGG.png'),width = 10,height = 5, p1)
}

# for (ii in genes2){
#   tar.exp <- t(exp[ii,])
#   cor_df <- cor(t(exp), tar.exp, method="pearson")
#   # print(head(cor_df))
#   data1 = data.frame(gene=rownames(cor_df),cor=cor_df[, 1])
#   # print(head(data1))
#   ge = data1$cor
#   names(ge) = data1$gene
#   ge = sort(ge,decreasing = T)
#   print(head(ge))
#   em_kegg <- GSEA(ge, TERM2GENE = KEGG_df,pAdjustMethod = "none")
#   print(head(em_kegg))
#   
#   p1 <- GseaVis::gseaNb(object = em_kegg,
#                         # rank.gene = ii,
#                         # rank.gene.nudgey= 8,#the gene label nudge y on rank plot, defalut is 2.
#                         # addGene = T, #whether add gene name on the curve, defalut is FALSE.
#                         # geneSize = 6,# gene label text size, defalut is 4.
#                         geneSetID = em_kegg@result$ID[1:5], 
#                         addPval = T,
#                         pvalX = 1,pvalY = 1,
#                         ght.facet = T,
#                         termWidth = 80, #the width or the term name, defalut is 40.
#                         curveCol = brewer.pal(5,'Paired'))+
#     labs(title = ii)+
#     theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18,vjust=55),
#           plot.subtitle = element_text(hjust = 0.5, color = "grey50"),
#           plot.title.position = "plot")
#   
#   p1
#   write.csv(em_kegg ,file=paste0("GSEA",ii, "_KEGG_gsea.csv"))
#   ggsave(paste0("GSEA",ii,'.KEGG.pdf'),width = 8,height = 5, p1) #GSEA分析
#   ggsave(paste0("GSEA",ii,'.KEGG.png'),width = 8,height = 5, p1)
# }
