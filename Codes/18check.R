rm(list = ls())

foldpath <- "D:/workdir/23/18check"

if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)

options(timeout = Inf)

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggpubr)

# check----
check <- read.csv("../../basedata/key_checkpoint_gene.csv", header=T)

load("../08batch/data_all.RData")
# data <- comdata
for (i in 1:4){
  colnames(data[[i]])[2] <- "group"
}

TrainExp <- data[[1]]

checkpoint <- check$Symbol[which(check$Symbol%in%colnames(TrainExp))]
TrainExp <- TrainExp[,c("ID","group",checkpoint)]
TrainExp <- na.omit(TrainExp)
TrainExp$group <- ifelse(TrainExp$group==1,"High","Low")
aa <- gather(data = TrainExp,key = gene,value = exp,-c("ID","group")) 

# -------------------做一个差异基因列表
library(rstatix)
stat_res <- aa %>% 
  group_by(gene) %>% 
  wilcox_test(exp ~ group) %>% 
  adjust_pvalue(method = "BH") %>% 
  add_significance("p")
stat_res
# write.csv(stat_res, file = "wilcoxon_res_diffgene.csv", row.names = F, quote = F)

deffgenes <- c()

for (i in 1:nrow(stat_res)){
  p <- stat_res[i, 8]
  if (p<0.05){
    deffgenes<-append(deffgenes, stat_res[i,1])
  }
}

stat_res_diff <- subset(stat_res,p<0.05)
write.csv(stat_res, file = "wilcoxon_res.csv", row.names = F, quote = F)
write.csv(stat_res_diff, file = "wilcoxon_res_diffgene.csv", row.names = F, quote = F)
# ------------------------------绘图

# pdf("01.check.pdf",w=12,h=6)
# ggplot(aa, aes(x=gene, y=exp, fill = group)) +
#   scale_fill_manual(values = c("High" = "#D20A13","Low" = "#58CDD9")) + #设置颜色
#   ylim(0,12)+
#   geom_boxplot()+
#   stat_compare_means(aes(group = group),method = 'wilcox.test',label = "p.signif",cex=4)+
#   # annotate(geom = "text", x = violin_dat1$gene, y = 13, size = 5,#(修改显著问题)
#   #          label =as.character( violin_dat1$P.Value.signif)) +
#   xlab("") + 
#   ylab("Expression") + 
#   labs(title = "")+ 
#   theme(panel.background=element_blank(),panel.border=element_rect(linetype="solid",fill=NA))+#去除背景加上黑色外边框
#   theme(axis.text.x=element_text(face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 ))#angle = 70表示横轴标题倾斜70度
# dev.off()

rownames(TrainExp) <- TrainExp$ID
TrainExp <- TrainExp[,-c(1,2)]
group <- group[rownames(TrainExp),]
group$group <- ifelse(group$group == 1, "High","Low")
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
Group = str_sub(colnames(TrainExp),1,str_length(colnames(TrainExp)))

pdf("01.check.pdf",w=12,h=6)
draw_boxplot(t(TrainExp) %>% as.data.frame(), factor(group$group),
                  drop = F,
                  method = "wilcox.test",
                  color = mypalette(length(unique(Group)))
                  ) +
  scale_fill_manual(values = c("High" = "#D20A13","Low" = "#58CDD9")) +
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_text(
    size = 16,
    angle = 60,     # 倾斜45度防标签重叠
    hjust = 1,      # 右对齐
    vjust = 1       # 垂直对齐
  ))
dev.off()

# png("01.check_old.png",w=12,h=6)
# ggplot(aa, aes(x=gene, y=exp, fill = risk)) +
#   scale_fill_manual(values = c("High" = "#D20A13","Low" = "#58CDD9")) + #设置颜色
#   ylim(0,12)+
#   geom_boxplot()+
#   stat_compare_means(aes(group = risk),method = 'wilcox.test',label = "p.signif",cex=4)+
#   # annotate(geom = "text", x = violin_dat1$gene, y = 13, size = 5,#(修改显著问题)
#   #          label =as.character( violin_dat1$P.Value.signif)) +
#   xlab("") + 
#   ylab("Expression") + 
#   labs(title = "")+ 
#   theme(panel.background=element_blank(),panel.border=element_rect(linetype="solid",fill=NA))+#去除背景加上黑色外边框
#   theme(axis.text.x=element_text(face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 ))#angle = 70表示横轴标题倾斜70度
# dev.off()

# ----------差异免疫检查点和风险评分基因相关性分析------------
library(dplyr)
library(psych)
library(Hmisc)

TrainExp <- readRDS("../00data/tcga_fpkm.rds")
check <- read.csv("wilcoxon_res_diffgene.csv", header=T)

riskScore <- read.table("../06COX/10.all_train_riskScore.txt",head = T)

TrainExp <- subset(TrainExp, rownames(TrainExp) %in% check$gene)
# TrainExp <- TrainExp[check$genes,]
TrainExp <- na.omit(TrainExp)
# riskScore <- na.omit(riskScore[,2])
riskScore$sample <-rownames(riskScore)
riskScore <- riskScore[,c(9,7)]

intersect_result <- intersect(rownames(riskScore),colnames(TrainExp))

TrainExp <- TrainExp[,intersect_result]


# 相关性分析
TrainExp <- t(TrainExp) %>% data.frame()
identical(rownames(TrainExp), riskScore$sample)
# riskScore <- riskScore[,1,drop = F] %>% rownames_to_column(.,var = "sample")
riskScore <- data.frame(riskScore)

riskScore <- na.omit(riskScore)
cor_r <- cor(TrainExp,riskScore[,2],method = "spearman") 

d <- corr.test(TrainExp,riskScore[,2],use="complete",method = 'spearman')
cor_p <- d$p
cor_r2 <- cor_r %>% as.data.frame %>% tibble::rownames_to_column(var = "gene") %>%
  tidyr::gather(., cell,Correlation,-gene)#转换数据长短 cor
cor_p2 <- cor_p %>% as.data.frame %>% tibble::rownames_to_column(var = "gene") %>% 
  tidyr::gather(., cell, Pvalue, -gene)#转换数据长短 p

cor_dat <- cbind(cor_r2, cor_p2)[,c("gene","cell","Correlation","Pvalue")]
cor_dat$cell <- "riskscore"

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

data$cell <- "riskscore"
p <- 
  ggplot(data, aes(cell,gene)) + 
  geom_tile(aes(fill = Correlation), colour = "black", size = 0.5)+
  scale_fill_gradient2(low = "#5C5DAF",mid = "white",high = "#EA2E2D") + # 这里可以用 windowns 小工具 takecolor 取色，看中哪个文章就吸哪个文章
  # 比如这篇 https://www.nature.com/articles/nmeth.1902 
  geom_text(aes(label = text),col ="black") +
  theme_minimal() + # 不要背景
  theme(axis.title.x=element_blank(), # 去掉 title
        axis.ticks.x=element_blank(), # 去掉x 轴
        axis.title.y=element_blank(), # 去掉 y 轴
        axis.text.x = element_text(hjust = 0.5,  face = "bold"), # 调整x轴文字，字体加粗
        axis.text.y = element_text(hjust = 0.5,face = "bold")) + #调整y轴文字
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation")) +   # 修改 legend 内容
  scale_x_discrete(position = "top") + #
  # scale_y_discrete(position = "top") +
  coord_fixed(ratio = 0.2)
p
ggsave(file=paste0('correlation_biomarker_','.png'), height = 6, width = 4 ,p)
ggsave(file=paste0('correlation_biomarker_','.pdf'), height = 6, width = 4, p)



# ---------------------------预后和免疫检查点相关性-------------

library(ggplot2)
library(magrittr)
library(tibble)
library(reshape2)
library(dplyr)
library(psych)
library(Hmisc)

TrainExp <- readRDS("../00data/tcga_fpkm.rds")
multi <- read.table("../06COX/Lasso_Coefficients.xls", header = T) %>% as.data.frame()#cox_result_step2.rds;mul_cox_result.rds
hubgene <- read.csv("../06COX/hub_gene.csv")
multi <- multi[multi$gene %in% hubgene$x,]

intersect_gene <- multi$gene

# TrainExp <- t(TrainExp) %>% data.frame()

checkExp <- TrainExp[check$gene,]
checkExp <- na.omit(checkExp)
preExp <- TrainExp[multi$gene,]
preExp <- na.omit(preExp)

# 相关性分析
checkExp <- t(checkExp) %>% data.frame()
preExp <- t(preExp) %>% data.frame()

cor_r <- cor(preExp,checkExp,method = "spearman") 

d <- corr.test(preExp,checkExp,use="complete",method = 'spearman')
cor_p <- d$p
cor_r2 <- cor_r %>% as.data.frame %>% tibble::rownames_to_column(var = "gene") %>%
  tidyr::gather(., cell,Correlation,-gene)#转换数据长短 cor
cor_p2 <- cor_p %>% as.data.frame %>% tibble::rownames_to_column(var = "gene") %>% 
  tidyr::gather(., cell, Pvalue, -gene)#转换数据长短 p

cor_dat <- cbind(cor_r2, cor_p2)[,c("gene","cell","Correlation","Pvalue")]
#cor_dat <- cor_dat[cor_dat$Cell %in% stat_res$Cell,]
write.csv(cor_dat,"04.pre_check_correlation_cor.csv")


# 相关性热图带显著性----
data <- read.csv('04.pre_check_correlation_cor.csv',row.names = 1,check.names = F)
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

p <- ggplot(data, aes(gene, cell)) + 
  geom_tile(aes(fill = Correlation), colour = "grey", size = 0.5) +
  scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") +
  geom_text(aes(label = text), col = "black", size = 5) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),# 去掉 title
    axis.ticks.x = element_blank(),# 去掉x 轴
    axis.title.y = element_blank(),# 去掉 y 轴
    axis.text.x = element_text(hjust = 0.5, size = 10, face = "bold"),
    axis.text.y = element_text(hjust = 0.5, size = 10, face = "bold"),
    legend.title = element_text(hjust = 0.5), # 调整图例标题对齐
    legend.text = element_text(hjust = 0.5) # 调整图例文本对齐
  ) +
  labs(
    fill = paste0(
      " * p < 0.05\n",
      " ** p < 0.01\n",
      " *** p < 0.001\n",
      "Correlation"
    )
  ) +
  scale_x_discrete(position = "top")

p
ggsave(file=paste0('pre_correlation','.png'), height = 6, width = 5,p)
ggsave(file=paste0('pre_correlation','.pdf'), height = 6, width = 5, p)
