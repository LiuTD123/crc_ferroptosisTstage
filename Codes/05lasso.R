rm(list = ls())
foldpath <- paste0("D:/workdir/23/05lasso")

if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)
library(magrittr)
library(ggplot2)
library(glmnet)
library(reshape2)
library(ggsci)
library(tibble)
library(caret)

# PH 机器学习----
### lasso -----
library(dplyr)
library(tidyverse)
library(xgboost)
library(randomForest)
library(DALEX)
library(stats)
library(e1071)
library(pROC)
library(scatterplot3d)

# candigene <- read.csv("../04boruta/10.boruta_geneids.csv")
candigene <- read.csv("../03venn/ferroptosis_genes.csv")

load("../08batch/data_all.RData")

for (i in 1:4){
  colnames(data[[i]])[2] <- "group"
}

dat <- data[[1]]

group <- dat[,c(1,2)]
rownames(dat) <- dat$ID
# dat <- dat[,-c(1,2)]

dat<-dat[group$ID,c("group",candigene$x)]
# dat$group<-factor(dat$group,levels = c('PD','Healthy'))
dat$group<- as.factor(dat$group)

set.seed(888) #设置种子

res.lasso <- cv.glmnet(as.matrix(dat[-1]), dat$group, family = "binomial", 
                       type.measure = "auc",nfolds = 10)

# 提取系数并准备数据进行美化
coef_df <- as.data.frame(as.matrix(coef(res.lasso$glmnet.fit)))
coef_df$coef <- row.names(coef_df)
coef_df <- melt(coef_df, id = "coef")
coef_df$variable <- as.numeric(gsub("s", "", coef_df$variable))
coef_df$lambda <- res.lasso$glmnet.fit$lambda[coef_df$variable + 1]
coef_df$norm <- apply(abs(coef(res.lasso$glmnet.fit)[-1,]), 2, sum)[coef_df$variable + 1]
coef_df <- coef_df[!coef_df$coef == "(Intercept)",]

# 美化Lasso逻辑系数惩罚图
pdf("lambda_beautified_PH_1.pdf", w= 8,  h =6)
ggplot(coef_df, aes(log(lambda), value, color = coef)) + 
  geom_vline(xintercept = log(res.lasso$lambda.min), size = 0.8, color = 'grey60', alpha = 0.8, linetype = 2) +
  geom_vline(xintercept = log(res.lasso$lambda.1se), size = 0.8, color = 'grey60', alpha = 0.8, linetype = 2) +
  geom_line(size = 1) + 
  xlab("log Lambda") + 
  ylab('Coefficients') + 
  theme_bw(base_rect_size = 2) + 
  scale_color_manual(values = c(pal_npg()(10), pal_d3()(10), pal_jco()(3), pal_lancet()(10), pal_aaas()(10), pal_simpsons()(10), pal_gsea()(10), pal_jama()(10))) + 
  scale_x_continuous(expand = c(0.1, 0.1)) +  # 调整坐标轴扩展
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 15, color = 'black'), 
        axis.text = element_text(size = 12, color = 'black'), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12, color = 'black'), 
        legend.position = 'right') + 
  annotate('text', x = log(res.lasso$lambda.min)-0.5, y = max(coef_df$value) - 0.1, label = paste('lambda.min =', round(res.lasso$lambda.min, 4)), color = 'black',size = 5) +  # 调整文本位置
  annotate('text', x = log(res.lasso$lambda.1se)-0, y = max(coef_df$value) - 0.3, label = paste('lambda.1se =', round(res.lasso$lambda.1se, 4)), color = 'black',size = 5) + 
  guides(col = guide_legend(ncol = 2))
dev.off()

# 交叉验证Lasso图
cv_df <- data.frame(lambda = res.lasso$lambda,
                    cvm = res.lasso$cvm,
                    cvsd = res.lasso$cvsd, 
                    cvup = res.lasso$cvup,
                    cvlo = res.lasso$cvlo,
                    nzero = res.lasso$nzero) 
cv_df$ll <- log(cv_df$lambda) 
cv_df$NZERO <- paste0(cv_df$nzero, ' vars')

pdf("cvlasso_beautified_PH_1.pdf", w= 6, h = 6)
ggplot(cv_df, aes(ll, cvm, color = NZERO)) + 
  geom_errorbar(aes(x = ll, ymin = cvlo, ymax = cvup), width = 0.05, size = 1) + 
  geom_vline(xintercept = log(res.lasso$lambda.min), size = 0.8, color = 'grey60', alpha = 0.8, linetype = 2) + 
  geom_vline(xintercept = log(res.lasso$lambda.1se), size = 0.8, color = 'grey60', alpha = 0.8, linetype = 2) + 
  geom_point(size = 2) + 
  xlab("Log Lambda") +
  ylab('Partial Likelihood Deviance') + 
  theme_bw(base_rect_size = 1.5) + 
  scale_color_manual(values = c(pal_npg()(10), pal_d3()(10), pal_lancet()(10), pal_aaas()(10))) + 
  scale_x_continuous(expand = c(0.1, 0.1)) +  # 调整坐标轴扩展
  scale_y_continuous(expand = c(0.02, 0.02)) + 
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 15, color = 'black'), 
        axis.text = element_text(size = 12, color = 'black'), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12, color = 'black'), 
        legend.position = 'bottom') + 
  annotate('text', x = log(res.lasso$lambda.min), y = max(cv_df$cvm)-0.07, label = paste('lambda.min =', round(res.lasso$lambda.min, 4)), color = 'black') +  # 调整文本位置
  annotate('text', x = log(res.lasso$lambda.1se)-0.5, y = max(cv_df$cvm)-0.12, label = paste('lambda.1se =', round(res.lasso$lambda.1se, 4)), color = 'black') + 
  guides(col = guide_legend(nrow = 3))
dev.off()

# # plot(res.lasso)
# # plot(res.lasso$glmnet.fit, xvar = 'lambda')
# ggsave("01.lasso_PH.CV.pdf", plot(res.lasso), width = 6, height = 5, dpi = 300, units = "in", bg = "white")
# ggsave("02.lasso_PH.Coef.pdf", plot(res.lasso$glmnet.fit, xvar = 'lambda'), width = 6, height = 5, dpi = 300, units = "in", bg = "white")
# # dev.off()
l.coef<-coef(res.lasso$glmnet.fit,s=res.lasso$lambda.min,exact= F)
l.coef
coef.min = coef(res.lasso, s = "lambda.min")  ## lambda.min & lambda.1se 取一个
res.lasso$lambda.1se #  0.07986446
res.lasso$lambda.min #  0.007802837

# 找出那些回归系数没有被惩罚为0的
active.min = which(coef.min@i != 0)
# 提取基因名称
lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1]
lasso_geneids <- lasso_geneids[-1]
lasso_geneids
# [1] 12个
lasso_geneids <- data.frame(lasso_geneids)
colnames(lasso_geneids) <- 'symbol'
write.csv(lasso_geneids,file = '03.lasso.gene_PH_1.csv',row.names = T)

# LASSO特征系数
df.coef = cbind(gene = rownames(coef.min), coefficient = coef.min[,1]) %>% as.data.frame()
df.coef = subset(df.coef, coefficient != 0) %>% as.data.frame
write.table(df.coef, "Lasso_Coefficients.xls", sep = "\t", quote = F, col.names = T, row.names = F)  #6个
