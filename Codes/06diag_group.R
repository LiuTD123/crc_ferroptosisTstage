rm(list = ls())

foldpath <- "D:/workdir/23/06diag"

if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)
library(glmnet)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(xgboost)
library(randomForest)
library(caret)
library(DALEX)
library(stats)
library(e1071)
library(readxl)
library(gplots)
library(tableone)

# set.seed(9999999)
# ----数据集整理---------
# load("D:/mydata/处理完成/CRC/t_data.RData")
load("../08batch/data_all.RData")
# data <- comdata
for (i in 1:4){
  colnames(data[[i]])[2] <- "group"
}

# modelgenes <- read.csv("../05lasso/03.lasso.gene_PH_1.csv")
modelgenes <- read.csv("../04boruta/10.boruta_geneids.csv")

for (i in 1:4){
  data[[i]]$group <- as.factor(data[[i]]$group)
  data[[i]] <- data[[i]][,-1]
  data[[i]] <- data[[i]][,c("group",modelgenes$x)]
}

library(dplyr)
input <- data[[1]]

testset <- rbind(data[[2]],
                 data[[3]]
                 # data[[4]]
)
data[[5]] <- testset
names(data)[5] <- "Test"

test <- data[[5]]
colnames(input)[14] <- "BLOC1S5_TXNDC5"
colnames(test)[14] <- "BLOC1S5_TXNDC5"
# 构建随机森林模型
# 注意：outcome ~ . 表示使用除outcome外的所有其他变量作为特征
rf_model <- randomForest(group ~ .,
                         data = input,
                         ntree = 500,       # 树的数量，通常500-1000足够
                         mtry = sqrt(ncol(input) - 1), # 默认值，特征数开方
                         importance = TRUE, # 计算变量重要性
                         proximity = TRUE)  # 计算样本 proximity，可用于聚类

# 查看模型基本信息
print(rf_model)

table(input$group)

# 对测试集进行预测 (type = "prob" 会得到概率)
test_predictions_prob <- predict(rf_model, newdata = test, type = "prob")
test_predictions_class <- predict(rf_model, newdata = test, type = "response")

# 查看预测的类别（基于默认0.5截断值）
table(Predicted = test_predictions_class, Actual = test$group)

# 计算准确率等指标
confusionMatrix(test_predictions_class, test$group)
# Confusion Matrix and Statistics
# 
# Reference
# Prediction  0  1
# 0  1  2
# 1 10 82
# 
# Accuracy : 0.8737         
# 95% CI : (0.7897, 0.933)
# No Information Rate : 0.8842         
# P-Value [Acc > NIR] : 0.69578        
# 
# Kappa : 0.0981         
# 
# Mcnemar's Test P-Value : 0.04331        
#                                          
#             Sensitivity : 0.09091        
#             Specificity : 0.97619        
#          Pos Pred Value : 0.33333        
#          Neg Pred Value : 0.89130        
#              Prevalence : 0.11579        
#          Detection Rate : 0.01053        
#    Detection Prevalence : 0.03158        
#       Balanced Accuracy : 0.53355        
#                                          
#        'Positive' Class : 0 


# 绘制ROC曲线并计算AUC
roc_obj <- roc(test$group, test_predictions_prob[, 2]) # 取高风险(1)的概率
plot(roc_obj, main = paste0("ROC Curve (AUC = ", auc(roc_obj), ")"))

# --------------
# 区分分组
# -----------------
# 方法一：使用中位数划分（最常用、最简单）
# 对全体数据（或测试集）进行预测，得到每个样本属于高风险的概率
all_predictions_prob <- predict(rf_model, newdata = test, type = "prob")
test$risk_score <- all_predictions_prob[, 2] # 获取属于1（高风险）的概率

# 根据风险评分的中位数进行分组
risk_median <- median(test$risk_score)
test$risk_group <- ifelse(test$risk_score > risk_median, "High", "Low")
test$risk_group <- factor(test$risk_group, levels = c("Low", "High"))
# ===============

# 方法二：使用ROC曲线 Youden指数确定最佳截断点（更精确）
# 使用训练集或全体数据来找到最佳截断点
# 这里以训练集为例
train_predictions_prob <- predict(rf_model, newdata = input, type = "prob")
input$risk_score <- train_predictions_prob[, 2] # 获取属于1（高风险）的概率
roc_train <- roc(input$group, train_predictions_prob[, 2])

# 找到Youden指数最大的点（敏感度+特异度-1）
best_threshold <- coords(roc_train, "best", transpose = FALSE, best.method = "youden")$threshold

# 使用这个最佳截断点对全体数据进行分组
input$risk_group_youden <- ifelse(input$risk_score > best_threshold, "High", "Low")
input$risk_group_youden <- factor(input$risk_group_youden, levels = c("Low", "High"))

# 查看分组情况
table(input$risk_group_youden)
