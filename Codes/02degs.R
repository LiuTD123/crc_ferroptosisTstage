rm(list = ls())

foldpath <- paste0("D:/workdir/23/02degs")

if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)

mycol=rep(c("#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)

library(DESeq2)
library(tibble)
library(limma)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(stringr)
library(pheatmap)
library(ComplexHeatmap)

deg_cutoff = 0.5

load("D:/mydata/处理完成/CRC/t_data.RData")
load("D:/mydata/处理完成/CRC/TCGA_tstage.RData")

if(T){
  tcga_count <- TCGA_count
  # colnames(tcga_count) <- substr(colnames(tcga_count) ,1,15)
  
  tcga_group <- tcga_count[,c(1,2)]
  # tcga_group$sample <- gsub("\\.1","",tcga_group$sample)
  
  tcga_group$T_stage <- as.factor(tcga_group$T_stage)
  tcga_group$T_stage <- ifelse(tcga_group$T_stage == 1,"High","Low")
  tcga_group$T_stage <- as.factor(tcga_group$T_stage)
  
  rownames(tcga_count) <- tcga_count$ID
  tcga_count <- tcga_count[,-c(1,2)]
  tcga_count <- as.data.frame(t(tcga_count))
  dds <- DESeqDataSetFromMatrix(countData = tcga_count, colData = tcga_group, design = ~ T_stage)
  dds <- DESeq(dds)
  
  resultsNames(dds)
  res <- results(dds,contrast = c("T_stage","High","Low")) 
  
  DE_RNA <- res %>% as.data.frame()
  logfc <- tem1 <- deg_cutoff
  # logfc <- tem1 <- 1
  # pSig <- "p.adj"
  pSig <- "pvalue"
  str(DE_RNA)
  
  if(pSig == "pvalue"){
    DE_RNA[which(DE_RNA$log2FoldChange > logfc & DE_RNA$pvalue < 0.05),'sig'] <- 'up'
    DE_RNA[which(DE_RNA$log2FoldChange < - logfc & DE_RNA$pvalue < 0.05),'sig'] <- 'down'
    if(logfc > 0){
      DE_RNA[which(DE_RNA$pvalue >= 0.05 & abs(DE_RNA$log2FoldChange) < logfc) ,'sig'] <- 'none'
    }else{
      DE_RNA[which(DE_RNA$pvalue >= 0.05 ) ,'sig'] <- 'none'
    }
    result_up <- subset(DE_RNA, sig == 'up')
    result_down <- subset(DE_RNA, sig == 'down')
    select_DERNA <- subset(DE_RNA,sig %in% c('up','down'))
  }else{
    DE_RNA[which(DE_RNA$log2FoldChange > logfc & DE_RNA$padj < 0.05),'sig'] <- 'up'
    DE_RNA[which(DE_RNA$log2FoldChange < - logfc & DE_RNA$padj < 0.05),'sig'] <- 'down'
    if(logfc > 0){
      DE_RNA[which(DE_RNA$padj >= 0.05 & abs(DE_RNA$log2FoldChange) < logfc) ,'sig'] <- 'none'
    }else{
      DE_RNA[which(DE_RNA$padj >= 0.05 ) ,'sig'] <- 'none'
    }
    result_up <- subset(DE_RNA, sig == 'up')
    result_down <- subset(DE_RNA, sig == 'down')
    select_DERNA <- subset(DE_RNA,sig %in% c('up','down'))
  }
  
  dim(select_DERNA) 
  
  paste0("共",dim(select_DERNA)[[1]],"个，其中上调",dim(result_up)[[1]],"个，下调",dim(result_down)[[1]],"个")
  # [1] "共1619个，其中上调772个，下调847个"
  write.table(select_DERNA, file = "02.DEG.txt",sep="\t",quote=F,col.names=T)
  write.table(result_down,file="03.DEG_down.txt",sep="\t",quote=F,col.names=T)
  write.table(result_up, file="04.DEG_up.txt",sep="\t",quote=F,col.names=T)
  ##所有基因信息
  write.table(res, file="01.DEG_ALL.txt",sep="\t",quote=F,col.names=T)
  #差异基因名
  saveRDS(res,"tcga_DEG2_ALL.rds")
  saveRDS(select_DERNA,"tcga_filter_diff2.rds")
  saveRDS(result_down,"TCGA_deg2_dow.rds")
  saveRDS(result_up,"TCGA_deg2_up.rds")
}
# volcanoplot ----
if(T){
  tem1 = deg_cutoff
  
  result_up <- readRDS("TCGA_deg2_up.rds") %>% as.data.frame()
  result_down <- readRDS("TCGA_deg2_dow.rds") %>% as.data.frame()
  deg_all <- readRDS("tcga_DEG2_ALL.rds") %>% as.data.frame() %>% na.omit(.)
  
  tem_up <- result_up[order(result_up$log2FoldChange,decreasing = T),]
  tem_down <- result_down[order(result_down$log2FoldChange,decreasing = T),]
  
  # 标注top10变异基因
  data_repel <- rbind(head(tem_up,10),tail(tem_down,10))
  
  data_repel$sig <- ifelse(str_detect(data_repel$log2FoldChange,"-") == TRUE,"down","up")
  data_repel <- data_repel[,-1]
  colnames(data_repel)[1] <- "logFC"
  result <- deg_all[,-1]
  colnames(result)[1] = "logFC"
  # logFC=log2(result$logFC)   #数据预处理
  cut_off_pvalue = 0.05  
  cut_off_logFC = tem1    
  result$change = ifelse(result$pvalue < cut_off_pvalue & abs(result$logFC) >= cut_off_logFC, 
                         ifelse(result$logFC> cut_off_logFC ,'up','down'),
                         'stable')
  p1 <- ggplot(
    result, aes(x = logFC, y = -log10(pvalue), colour=change)) +
    geom_point(alpha=0.6, size=2) +
    # "#223D6C","#D20A13"
    scale_color_manual(values=c("#223D6C", "gray","#D20A13"))+
    geom_vline(xintercept=c(-tem1,tem1),lty=4,col="#828282",lwd=0.8) +
    geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="#828282",lwd=0.8) +
    labs(x="log2FC",
         y="-log10 (pvalue)")+
    theme_bw()+
    # scale_x_continuous(limits = c(-15,18))+
    # scale_y_continuous(limits = c(0,250))+
    theme(legend.justification = c(1,1),
          legend.position.inside = c(1 ,1),
          legend.background = element_rect(fill = "white", color = "black", size = 0.2),
          panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(face="bold",color="black",family = "Times",size=12),
          plot.title = element_text(hjust = 0.5, face = "bold",color = "black",family = "Times",size = 18),
          axis.text.x = element_text(face = "bold",color = "black",size = 12),
          axis.text.y = element_text(face = "bold",color = "black",size = 12),
          axis.title.x = element_text(face = "bold",color = "black",family = "Times",size = 18),
          axis.title.y = element_text(face = "bold",color = "black",family = "Times",size = 18),
          plot.subtitle = element_text(hjust = 0.5,family = "Times", size = 12, face = "italic", colour = "black")
    )+ 
    geom_label_repel(
      data = data_repel[which(data_repel$sig=="up"),],
      aes(label = rownames((data_repel[which(data_repel$sig=="up"),]))),
      fontface = "italic",
      force = 1,  # 设置标签重叠时的排斥力
      size = 5, #size
      color = "black",
      family = "Times",
      segment.color = "grey",
      segment.size =1,
      segment.alpha = 0.5,
      show.legend = F,
      direction = "y",    # 按y轴调整标签位置方向，若想水平对齐则为x
      # hjust = 1,   # 对齐标签：0右对齐，1左对齐，0.5居中
      # xlim = c(5, 13),
      # ylim = c(5,90)
      max.overlaps = getOption("ggrepel.max.overlaps", default = 10)               # 设置排斥重叠过多标签的阈值为无穷大，保持始终显示所有标签
    ) +
    geom_label_repel(
      data = data_repel[which(data_repel$sig=="down"),],
      aes(label = rownames((data_repel[which(data_repel$sig=="down"),]))),
      fontface = "italic",
      size = 5,
      color = "black",
      family = "Times",
      segment.color = "grey",
      segment.size = 1,
      segment.alpha = 0.5,
      show.legend = F,
      direction = "y",
      # hjust = 1,
      force = 1,   # 设置标签重叠时的排斥力
      # xlim = c(-5,-10),
      # ylim = c(5,90)
      max.overlaps = getOption("ggrepel.max.overlaps", default = 10)                # 设置排斥重叠过多标签的阈值为无穷大，保持始终显示所有标签
    ) 
  ggsave(filename = "03.DEGs_volcano.pdf", width = 7, height = 8,p1)
}

# heatmap-----

if(T){
  rm(list = ls())
  
  count_data <- readRDS("..\\00data\\tcga_count.rds") %>% na.omit()
  tcga_group <-readRDS("..\\00data\\tcga_group.rds")
  tcga_group <- tcga_group[order(tcga_group$group,decreasing = T),]
  count_data <- count_data[,tcga_group$sample]
  
  select_DERNA <- readRDS("tcga_filter_diff2.rds")%>% as.data.frame()%>% na.omit()
  DE_RNA <- select_DERNA[order(select_DERNA$log2FoldChange,decreasing = T),] 
  
  data_repel2 <-  rbind(DE_RNA,tail(DE_RNA,10))
  data_repel2$Symbols <- rownames(data_repel2)
  args_type <- "count"
  if (args_type != "GEO2R"){
    if ( args_type == "count"){
      exprSet <- as.data.frame(log(edgeR::cpm(count_data)+0.1))
      mat <- exprSet[rownames(data_repel2),]
    } else {
      mat  <- expr[(rownames(data_repel2)),]
    }}
  mat <- t(scale(t(mat)))
  mat[mat > 2] <- 2
  mat[mat < -2] <- -2
  group <- tcga_group
  group$group <- factor(group$group, levels = c("Normal", "Tumor"))
  group <- group[!duplicated(group$sample),]
  annotation_col <- data.frame(
    Group = group$group,
    row.names = group$sample
  )
  # mat <- mat[data_repel2[order(data_repel2$log2FoldChange, decreasing = T),]$Symbols,] %>% as.data.frame()
  mat <- mat[,group$sample]
  gene_col <- data.frame(Expression=c(rep("Up",10), rep("Down", 10)))
  # rownames(gene_col) <- rownames(mat)
  
  newnames <- lapply(
    rownames(mat),
    function(x) bquote(italic(.(x))))
  
  ##  https://www.jianshu.com/p/671ea93108e4
  p <- densityHeatmap(mat ,title = "Distribution as heatmap", ylab = " ",height = unit(8, "cm")) %v%
    HeatmapAnnotation(Group = annotation_col$Group, col = list(Group = c("Tumor" = "#D20A13", "Normal" = "#223D6C"))) %v%
    Heatmap(mat, row_names_gp = gpar(fontsize = 9),
            cluster_rows = F,
            show_column_names = F,
            show_row_names = F,
            name = "mat", height = unit(6, "cm"),
            col = colorRampPalette(c("#58CDD9", "#D20A13"))(100))
  
  pdf(file = paste0("04.DEGs_pheatmap.pdf"),width = 8,height = 8)
  par(mar = c(2,2,2,2),cex=1.5,family="Times")
  print(p)
  dev.off()
}

# -----heatmap------
rt <- count_data
rt <- rt[rownames(DE_RNA),]
#group <- colData

heat<-rt[rownames(rt)%in%
           rownames(rbind(DE_RNA[order(DE_RNA$log2FoldChange,decreasing = T),],
                          DE_RNA[order(DE_RNA$log2FoldChange,decreasing = F),])),]

x <- heat #log2(heat+1)

a<-tcga_group[order(tcga_group$group),]
x<-x[,a$sample]

mat <- t(scale(t(x)))#归一化
df1 <- as.data.frame(mat)
# mat[mat < (-2)] <- (-2)
# mat[mat > 2] <- 2

pdf('03.heatmap_all.pdf',  w=6,h=6,family='Times')
densityHeatmap(mat ,title = "Distribution as heatmap", ylab = " ",height = unit(6, "cm")) %v%
  HeatmapAnnotation(Group = a$group, col = list(Group = c("Tumor" = "#B72230", "Normal" = "#104680"))) %v%
  Heatmap(mat, 
          row_names_gp = gpar(fontsize = 9),
          show_column_names = F,
          show_row_names = F,
          ###show_colnames = FALSE,
          name = "expression", 
          ##cluster_cols = F,
          cluster_rows = T,
          height = unit(6, "cm"),
          #cluster_columns = FALSE,
          ###cluster_rows = FALSE,
          col = colorRampPalette(c("#20B2AA", "white","red"))(100))
dev.off()
