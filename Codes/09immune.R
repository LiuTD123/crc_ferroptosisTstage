options(timeout = Inf)

rm(list = ls())
foldpath <- "D:/workdir/23/09immune"
if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)
#####免疫浸润分析
library(psych)
library(ggcorrplot)
library(corrplot)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(tibble)
library(reshape2)
library(magrittr)
library(GSVA)
library(ggplot2)
library(pheatmap)
library(immunedeconv)
library(tidyr)
library(IOBR)
library(reshape2)
library(rstatix)
library(patchwork)
library(CIBERSORT)
# rforge <- "http://r-forge.r-project.org"
# install.packages("estimate", repos = rforge, dependencies = TRUE)
# install.packages("D:/workdir/Rpackage/estimate_1.0.13.tar.gz")

# 输入参考
# OvarianCancerExpr <- system.file("extdata", "sample_input.txt", package = "estimate")
# exampleinput <- read.table(OvarianCancerExpr)[1:4,1:4]
# s516      s518      s519      s520
# C9orf152  4.881540  4.575656  3.739469  3.695996
# ELMO2     7.298054  7.555440  7.533202  7.382355
# CREB3L1   5.569164  5.700406  5.959730  5.770007
# RPS11    13.389937 13.848820 13.642862 13.654622

library(estimate)

load("../08batch/data_all.RData")
data_train <- data[[1]]
rownames(data_train) <- data_train$ID
data_train <- data_train[,-c(1,2)]
data_train <- as.data.frame(t(data_train))

head(data_train)[1:4,1:4]

write.table(data_train,file = "comexp.txt",sep = "\t", quote =F)
input_file_dir <- './comexp.txt'
output_file_dir <- './genes.gct'
output_estimate <- './estimate_score.gct'

which(duplicated(data_train$symbol))

# 根据表达谱生成gct
filterCommonGenes(input.f = input_file_dir, 
                  output.f = "genes.gct", 
                  id = "GeneSymbol")

estimateScore(input.ds = output_file_dir,
              output.ds = "estimate_score.gct", 
              platform = "affymetrix")

scores <- read.table("estimate_score.gct",skip = 2,header = T)
rownames(scores) <- scores[,1]
scores <- t(scores[,3:ncol(scores)])

# plotPurity函数可以根据保存好的文件来挑选对应的样本进行可视化：
# plotPurity(scores = "estimate_score.gct", samples = "GSM1192715", platform = "affymetrix")


# -------------------------
# =========================
# lasso_geneids <- read.csv("../04cox/01.lasso_genes.csv")
# genes <- lasso_geneids$symbol

load("../04boruta/borutagenes.RData")

genes <- c("SLC2A3",
           "PTGS2",
           "ANGPTL7",
           "FABP4",
           "RGS4")

genes <- boruta_geneids
# --------------------
data_train <- data[[1]]
# exprSet <- data.frame(scores,
#                       CSF2RB = data_train[,"CSF2RB"],
#                       COL15A1 = data_train[,"COL15A1"],
#                       MME = data_train[,"MME"],
#                       NEFM = data_train[,"NEFM"],
#                       CYP24A1 = data_train[,"CYP24A1"]
# )

exprSet <- data.frame(scores,
                      value1 = data_train[,genes[1]],
                      value2 = data_train[,genes[2]],
                      value3 = data_train[,genes[3]],
                      value4 = data_train[,genes[4]],
                      value5 = data_train[,genes[5]]
                      # value5 = data_train[,genes[6]]
)

names(exprSet)[5:(length(genes)+4)] <- genes

library(ggstatsplot)
library(rlang)

methods <- c("StromalScore",
             "ImmuneScore",
             "ESTIMATEScore",
             "TumorPurity"
)

for (method in methods){
  plot_list <- list()
  
  for (i in 1:length(genes)){
    cor_test <- cor.test(
      x = exprSet[[method]],
      y = exprSet[[genes[i]]],
      method = "spearman",
      exact = FALSE  # 近似计算提高大样本效率
    )
    
    # 格式化输出：rho保留3位小数，p值科学计数法
    rho <- formatC(cor_test$estimate, digits = 3, format = "f")
    p_value <- ifelse(
      cor_test$p.value < 0.001,
      "< 0.001",
      formatC(cor_test$p.value, digits = 3, format = "f")
    )
    
    # 创建自定义副标题
    custom_subtitle <- bquote(
      rho == .(rho) * ", " * italic(p) * .(paste0(" = ", p_value))
    )
    
    # 生成图形（移除公式和统计标签）
    plot_list[[i]] <- ggscatterstats(
      data = exprSet,
      x = !!sym(method),
      y = !!sym(genes[i]),
      marginal = TRUE,
      type = "nonparametric",
      title = "log2TPM",
      subtitle = custom_subtitle,  # 使用自定义副标题
      label.args = list(label = NULL),  # 移除图中统计标签
      results.subtitle = FALSE,  # 禁用默认副标题
      bf.message = FALSE  # 移除贝叶斯因子信息
    )
  }
  
  p <- Reduce(`+`, plot_list) + patchwork::plot_layout(ncol = length(plot_list))
  
  ggsave(p,file = paste0(method,".pdf"),w = 4*length(plot_list),h=4)
}

# =======================================================
# # ----------------------------cibersort评估细胞丰度-免疫细胞比例--------------------------
a = data_train

k = !duplicated(rownames(a));table(k)
exp = a[k,]

rownames(exp) = exp$ID
colnames(exp) = str_remove(colnames(exp),"TPM")

exp[1:4,1:4]
exp2 = as.data.frame(exp)
# exp2 = rownames_to_column(exp2)
# exp2 <- as.data.frame(t(exp))
write.table(exp2,file = "exp.txt",row.names = F,quote = F,sep = "\t")

# 内置数据LM22.txt，记录了22种免疫细胞的基因表达特征数据。
lm22f = system.file("extdata", "LM22.txt", package = "CIBERSORT")

exp2 <- na.omit(exp2)
TME.results = cibersort(lm22f,
                        "comexp.txt",
                        QN = F
)

# TME.results <- TME.results[-1,]
save(TME.results,file = "TME.result.RData")
TME.results[1:4,1:4]
# B cells naive B cells memory Plasma cells T cells CD8
# GSM404005   0.000000000    0.006760713    0.6669458  0.01723630
# GSM404006   0.000000000    0.031748910    0.7168789  0.00000000
# GSM404007   0.006649664    0.008262864    0.3980714  0.00000000
# GSM404008   0.000000000    0.071516561    0.5579762  0.01617254

re <- TME.results[,-(23:25)]

# 堆积柱状图
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
Group = str_sub(colnames(exp),1,str_length(colnames(exp)))

# 肿瘤正常--------
# tcga_group <-readRDS("..\\00data\\tcga_group.rds")
# group <- tcga_group
# 高低风险--------
group <- data[[1]][,c(1,2)]
table(group$group)
group$group <- ifelse(group$group==1,"High","Low")

re2 <- re[group$ID,]
exp2 <- as.data.frame(t(exp2[group$ID,-c(1,2)]))

# ---------
dat <- re2 %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  mutate(group = group$group) %>%
  gather(key = Cell_type,value = Proportion,-ID,-group) %>%
  arrange(group)

dat$ID = factor(dat$ID,ordered = T,levels = unique(dat$ID)) #定横坐标顺序

write.csv(dat,file = "cell_prop_risk.csv")

# 先把group排序，然后将sample设为了因子，确定排序后的顺序为水平，所以两图的顺序是对应的。
dat2 = data.frame(a = 1:ncol(exp2),
                  b = 1,
                  group = sort(group$group))

p1 = ggplot(dat2,aes(x = a, y = b)) + 
  geom_tile(aes(fill = group)) + 
  scale_fill_manual(values = c("Low"="#4DBBD5FF","High"="#E64B35FF")) +
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        axis.title = element_blank()) + 
  scale_x_continuous(expand = c(0, 0)) +
  labs(fill = "Group")
p2 = ggplot(dat,aes(ID, Proportion,fill = Cell_type)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
  ) + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(22))

library(patchwork)

p <- p1 / p2 + plot_layout(heights = c(1,10),guides = "collect" ) &
  theme(legend.position = "bottom")
p
ggsave(filename = '04.histon_plot_risk.pdf',p,w=15,h=8)
# 
# # #------------------------- 箱线图--------------------------
# # # 全是0的行去掉
k = colSums(re)>0;table(k)
# #
# # ## k
# # ## FALSE  TRUE
# # ##     1    21
# #
re2 = re2[,k]
library(tinyarray)

# group <- t(group)
# factor(group[2,])
# #

# p <- draw_boxplot(t(re2)%>%as.data.frame(),factor(group$group,levels = c("High","Low"),
#                                                   labels = c("High","Low")),
#                   drop = F,
#                   method = "wilcox.test",
#                   color = mypalette(length(unique(Group))))+
#   scale_fill_manual(values = c("Low"="#4DBBD5FF","High"="#E64B35FF"))+
#   labs(x = "Cell Type", y = "Estimated Proportion")

p <- draw_boxplot(t(re) %>% as.data.frame(), factor(group$group),
                  drop = F,
                  method = "wilcox.test",
                  color = mypalette(length(unique(Group)))) +
  scale_fill_manual(values = c("Low"="#4DBBD5FF","High"="#E64B35FF")) +
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_text(
    size = 16,
    angle = 60,     # 倾斜45度防标签重叠
    hjust = 1,      # 右对齐
    vjust = 1       # 垂直对齐
  ))

p
ggsave(filename = '05.box_plot2_risk.pdf',p,w=8,h=6)

# 细胞相关性----
save.image("statistic_over.Rdata")

load("statistic_over.Rdata")

tiics_result <- re2

tes <- p$plot_env$p$data
pvalue <- compare_means(exp ~ group, data = tes, group.by = "rows",
                        symnum.args = list(cutpoints = c(0, 0.01, 0.05, 1), 
                                           symbols = c("***", "*", "ns")))
print(pvalue[pvalue$p<0.05,]$rows)

tiics_result <- t(tiics_result) %>% as.data.frame()

diffcell <- as.character(pvalue[pvalue$p<0.05,]$rows)

res3 <- tiics_result[diffcell,]
res3 <- t(res3)%>%as.data.frame()

# #过滤掉表达为0的
# res3 <- res3[,which(colSums(res3) > 0)]
# res3 <- res3[,order(colnames(res3),decreasing = F)]
library(ggcorrplot)
cor_data <- cor(res3,method="spearman")
corp <- cor_pmat(res3)

write.csv(cor_data,'03.cell_cor_r.csv',quote=F)
write.csv(corp,'03.cell_cor_p.csv',quote=F)
env.cor <- round(cor((res3),method="spearman"), 3)
# env.p <-round(cor_pmat((gene_dat),method = "spearman"),3) 
cor_p <- WGCNA::corPvalueStudent(env.cor,nrow(res3))

pdf("cor.pdf", width = 7, height = 7)
cor.plot<-corrplot(corr =env.cor,p.mat = cor_p,type="upper",
                   col = colorRampPalette(c("#4DBBD5FF", "white", "#E64B35FF"))(50),
                   tl.pos="lt",tl.col="black", 
                   insig = "label_sig", sig.level = c(.001,.01, .05),
                   pch.cex=1.5,pch.col = "black",order = "AOE")
cor.plot<-corrplot(corr = env.cor,type="lower",add=TRUE,method="number",
                   col = colorRampPalette(c("#4DBBD5FF", "white", "#E64B35FF"))(50),
                   tl.pos="n",tl.col="black",tl.cex=1.2,
                   diag=FALSE, cl.pos="n",pch.col = "black",
                   number.cex = 1,order = "AOE")
dev.off()

#  生物标志物与关键免疫细胞相关性----
exprlog <- exp2
# group <- comgroup
# multi <- modelgenes

gene <- genes

expr <- t(exprlog) %>% as.data.frame()
expr <- expr[,gene]
colnames(group) <- c('sample', 'group')

tiics_result <- re2
tiics_result <- t(tiics_result) %>% as.data.frame()
tiics_result <- tiics_result[, group$sample] %>% as.matrix()
tiics_result <- t(tiics_result) %>% as.data.frame()
tiics_result <- tiics_result[rownames(expr),]

diff <- diffcell
tiics_result <- tiics_result[,diff]


tem <- intersect(rownames(expr), rownames(tiics_result))
expr <- expr[tem,]
tiics_result <- tiics_result[tem,]
identical(rownames(expr), rownames(tiics_result))
cor_r <- cor(expr,tiics_result,method = "spearman") 
#cor_p <- WGCNA::corPvalueStudent(cor_r,length(rownames(dat_exp_diff)))
library(psych)
d <- corr.test(expr,tiics_result,use="complete",method = 'spearman')
cor_p <- d$p
cor_r2 <- cor_r %>% as.data.frame %>% tibble::rownames_to_column(var = "gene") %>%
  tidyr::gather(., cell,Correlation,-gene)#转换数据长短 cor
cor_p2 <- cor_p %>% as.data.frame %>% tibble::rownames_to_column(var = "gene") %>% 
  tidyr::gather(., cell, Pvalue, -gene)#转换数据长短 p

cor_dat <- cbind(cor_r2, cor_p2)[,c("gene","cell","Correlation","Pvalue")]
#cor_dat <- cor_dat[cor_dat$Cell %in% stat_res$Cell,]
write.csv(cor_dat,"04.correlation_cor.csv")

# 相关性热图带显著性----
data <- read.csv('04.correlation_cor.csv',row.names = 1,check.names = F)
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
p <- 
  ggplot(data, aes(gene, cell)) + 
  geom_tile(aes(fill = Correlation), colour = "grey", size = 1)+
  scale_fill_gradient2(low = "#4DBBD5FF",mid = "white",high = "#E64B35FF") + # 这里可以用 windowns 小工具 takecolor 取色，看中哪个文章就吸哪个文章
  # 比如这篇 https://www.nature.com/articles/nmeth.1902 
  geom_text(aes(label = text),col ="black",size = 5) +
  theme_minimal() + # 不要背景
  theme(axis.title.x=element_blank(), # 去掉 title
        axis.ticks.x=element_blank(), # 去掉x 轴
        axis.title.y=element_blank(), # 去掉 y 轴
        axis.text.x = element_text(hjust = 0.5, size = 14, face = "bold",angle = 90,vjust = 0), # 调整x轴文字，字体加粗
        axis.text.y = element_text(size = 14, face = "bold")) + #调整y轴文字
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation")) +   # 修改 legend 内容
  scale_x_discrete(position = "top") #
p

ggsave(file=paste0('correlation_biomarker_','.pdf'), height = 4, width = 7, p)

