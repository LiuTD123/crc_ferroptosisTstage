rm(list = ls())
foldpath <- "D:/workdir/23/11expr"
if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)
library(pROC)
library(glmnet)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tibble)
library(tidyr)
# -------------------------------------------训练集表达量验证
load("../08batch/data_all.RData")
lasso_geneids <- read.csv("../04boruta/10.boruta_geneids.csv")
hub_gene <- lasso_geneids$x

data_train <- data[[1]]

data_train <- as.data.frame(data_train[,c("group",hub_gene)])

table(data_train$group)

aa <- gather(data = data_train,key = gene,value = exp,-c(,"group")) 

# -------------------做一个差异基因列表
# ------------------------------绘图

violin_dat1 <- read.table("../02degs/01.DEG_ALL.txt",sep = "\t",header = T)
violin_dat1$gene <- rownames(violin_dat1)
violin_dat1 <- violin_dat1[gene,]

data <- violin_dat1 %>%
  mutate(text = case_when( #设置label，并加入判断，当P值符合特定条件就显示"\n"外加特定数量的*号
    pvalue <= 0.001 ~ "***", #P<0.001就显示回车加三个星号
    between(pvalue, 0.001, 0.01) ~ "**", #P为0.001-0.01 显示回车加两个*号
    between(pvalue, 0.01, 0.05) ~ "*",  #P为0.01-0.05 显示回车加一个星号
    T ~ "ns"))
data <- data[data$gene %in% hub_gene,]

aa$group <- factor(aa$group, levels = c(1,0),
                   labels = c("High","Low"))

pdf("01.train_gene.pdf",w=6,h=4)
ggplot(aa, aes(x=gene, y=exp, fill = group)) +
  scale_fill_manual(values = c("Low"="#4DBBD5FF","High"="#E64B35FF")) + #设置颜色
  # ylim(3,15)+
  geom_boxplot()+
  # stat_compare_means(aes(group = group),method = 'wilcox.test',label = "p.signif",cex=4)+
  annotate(geom = "text", x = data$gene, y = 13, size = 5,#(修改显著问题)
           label =as.character(data$text)) +
  xlab("") + 
  ylab("Expression") + 
  labs(title = "")+ 
  theme(panel.background=element_blank(),panel.border=element_rect(linetype="solid",fill=NA))+#去除背景加上黑色外边框
  theme(axis.text.x=element_text(face = "bold", color="black",angle = 70,vjust = 1, hjust = 1 ))#angle = 70表示横轴标题倾斜70度
dev.off()

# ------------------循环-----------
for (i in 2:length(data)){
  data_train <- data[[i]]
  
  data_train <- as.data.frame(data_train[,c("group",hub_gene)])
  
  aa <- gather(data = data_train,key = gene,value = exp,-c(,"group")) 

  aa$group <- factor(aa$group, levels = c(1,0),
                     labels = c("High","Low"))
  
  p <- ggplot(aa, aes(x=gene, y=exp, fill = group)) +
    scale_fill_manual(values = c("Low"="#4DBBD5FF","High"="#E64B35FF")) + 
    geom_boxplot() +
    stat_compare_means(aes(group = group), 
                       method = 'wilcox.test',
                       label = "p.signif",
                       cex = 4) +
    # annotate(geom = "text", 
    #          x = data$gene, 
    #          y = 13, 
    #          size = 5,
    #          label = as.character(data$text)) +
    xlab("") + 
    ylab("Expression") + 
    labs(title = paste0(names(data)[i]," Expression")) +  # 添加标题
    theme(
      panel.background = element_blank(),
      panel.border = element_rect(linetype = "solid", fill = NA),
      axis.text.x = element_text(face = "bold", 
                                 color = "black",
                                 angle = 70, 
                                 vjust = 1, 
                                 hjust = 1),
      plot.title = element_text(hjust = 0.5,  # 标题居中
                                face = "bold",  # 粗体
                                size = 16)     # 增大字号
    )
  
  ggsave(p, file = paste0(names(data)[i],"_expr.pdf"), w=6,h=4)
}
