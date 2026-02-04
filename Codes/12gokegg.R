options(timeout = Inf)
rm(list = ls())

foldpath <- "D:/workdir/20DM009/07gokegg"
if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)
library(tibble)
library(GOplot)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(GOplot)

load("D:/mydata/处理完成/胃癌/00data/diagnostic_data.RData")
lasso_geneids <- read.csv("../04cox/01.lasso_genes.csv")

load("../02degs/inter_genes.RData")
inputgenes <- inter
# inputgenes <- lasso_geneids$symbol

data_train <- data[[1]]

GO <- enrichGO(gene = inputgenes,
               OrgDb  ="org.Hs.eg.db",
               keyType = "SYMBOL",
               pAdjustMethod = "BH",
               pvalueCutoff =0.05,
               minGSSize = 10,
               ont="all",
               readable =T)
save(GO,file = 'GO.Rdata')
write.table(GO,file = "GO.xls",sep = "\t",quote = F,row.names = F)
# 弦图----
go_result<-read.table('GO.xls',header = T,sep = '\t',check.names = F)
go=data.frame(Category = go_result$ONTOLOGY,
              ID = go_result$ID,
              Term = go_result$Description,
              Genes = gsub("/", ", ", go_result$geneID),
              adj_pval = go_result$p.adjust)

load("../02degs/TraindataDEGs.RData")
diffSig <- diffSig %>% as.data.frame()

id.fc <-diffSig
genelist <- data.frame(ID = rownames(id.fc), logFC = id.fc$logFC)
rownames(genelist)=genelist[,1]
circ_GO<-GOplot::circle_dat(go,genelist)
pdf('01.ChordalChart.pdf',w=9,h=6)
GOCircle(data=circ_GO,
         nsub=10, ###指定显示通路个数
         rad1 = 2, rad2 = 3, ##rad1和rad2分别代表内圈和外圈的大小
         zsc.col=c("#EE0000B2",'#5F97D2'),  ##内圈Z-score颜色设置
         lfc.col = c('#EE0000B2','#5F97D2'),  ##logFC颜色设置
         label.size = 4,  ##标签字号设置
         label.fontface='bold',   ##标签字体设置
         table.legend = T   ##是否显示右侧表格，若为F则不显示
)
dev.off()

# KEGG
gene_transform <- bitr(inputgenes,
                       fromType = "SYMBOL",
                       toType = c("ENTREZID",'ENSEMBL'),
                       OrgDb = "org.Hs.eg.db")
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 organism = "hsa",
                 keyType = "kegg",
                 pvalueCutoff =0.05,
                 qvalueCutoff =1)

kk <- setReadable(kk, #前面分析的结果
                  OrgDb = "org.Hs.eg.db", #人类数据库
                  keyType = "ENTREZID") #要转换的基因类型

kegg_result <- kk@result
hh <- as.data.frame(kegg_result)
rownames(hh) <- 1:nrow(hh)
hh <- hh[hh$pvalue <= 0.05,]
hh <- hh[order(hh$Count,decreasing = T),]
dim(hh) # 44  14
write.table(hh,file = "00.KEGG.txt",sep = "\t",quote = F,row.names = F)

rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
kk <- hh[1:10,]
paste0(kk$Description,collapse = "（）；")
kegg <-  data.frame(Category = hh[,1],ID = hh[,'ID'],Term = hh[,'Description'], 
                    Genes = gsub("/", ", ", hh[,'geneID']), adj_pval = hh[,'pvalue'])

# diffSig <- readRDS("../01DEG/tcga_DEG2_ALL.rds")
# diffSig <- diffSig %>% as.data.frame()
# id.fc <- diffSig
#基因变化倍数
head(id.fc)
genelist <- data.frame(ID = rownames(id.fc), logFC = id.fc$logFC)             

#把富集分析和倍数整合在一起
circ <- circle_dat(kegg, genelist)
head(circ)
circ.gsym <- circ
#参数设置
n = 10 #
chord <- chord_dat(circ, genelist, kegg$Term[1:n])
head(chord)
chord.gsym <- chord
head(chord.gsym)
# 更改函数
# trace(GOChord,edit = T)
# 绘图
p1 <- GOChord(chord.gsym, 
              space = 0.01, #基因方块间隙
              gene.order = 'logFC', 
              lfc.col = c('#EE0000B2', 'white', '#5F97D2'),
              gene.space = 0.3,
              gene.size = 6, 
              border.size = 0.1, 
              ribbon.col=brewer.pal(10, "Set3"),
              process.label = 12) 
p1#
# 保存结果

pdf(file = paste0("02.KEGG_chord.pdf"),width = 12,height = 13)

# par(mar = c(2,2,2,5),cex=1,family="Times")
p1

dev.off()
untrace(GOChord)
