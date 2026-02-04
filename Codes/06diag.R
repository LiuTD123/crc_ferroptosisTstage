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
  # data[[i]]$group <- ifelse(data[[i]]$group == "PD",1,0)
  data[[i]] <- data[[i]][,-1]
  data[[i]] <- data[[i]][,c("group",modelgenes$x)]
}

library(dplyr)
input <- data[[1]]

ctrl <- trainControl(method="cv", number=100)
classif_knn <- train(group~., data = input,
                     method = "knn", trControl=ctrl, tuneLength=10)

classif_pls <- train(group~., data = input,
                     method = "pls",validation = 'CV')

# 支持向量机模型
classif_svm <- train(group~., data = input, 
                     method = "svmRadial")

# 随机森林模型
classif_rf <- train(group~., data = input, 
                    method = "rf",
                    ntree = 300)
# 广义线性模型
classif_glm <- train(group~., data = input,
                     method = 'glm')


# # 极限梯度提升模型
x = model.matrix(group~.,input)
model_martix_train<-model.matrix(group~., input)

data_train <- xgb.DMatrix(x, label =as.numeric(as.factor(input$group)))

params <- list(
  objective = "reg:squarederror"
)
# 
classif_xgboost <- xgb.train(params, data_train, nrounds = 100)
# 
# save(classif_xgboost,classif_glm,classif_rf,classif_svm,file="fourModel_classif.RData")
explainer_knn<-explain(classif_knn,label = "KNN",
                       data = input,
                       y = input$group)

explainer_pls<-explain(classif_pls,label = "PLS",
                       data = input,
                       y = input$group)

explainer_svm<-explain(classif_svm,label = "SVM",
                       data = input,
                       y = input$group)

explainer_rf<-explain(classif_rf,label = "RF",
                      data = input,
                      y = input$group)

explainer_glm<-explain(classif_glm,label = "GLM",
                       data = input,
                       y = input$group)

# -------------------------------xgboost
predict_logit <- function(model,x){
  raw_x <-predict(model,x)
  exp(raw_x)/(1+exp(raw_x))
}

logit <- function(x){
  exp(x)/(1+exp(x))
}

explainer_xgboost<-explain(classif_xgboost,
                           label = "xgboost",
                           data = x,
                           y = as.numeric(input$group),
                           predict_function = predict_logit,
                           link = logit
)

# save(explainer_xgboost,explainer_glm,explainer_rf,explainer_svm,file="fourModel_explainer.RData")

# model performance
# per_knn<-model_performance(explainer_knn)
library(pROC)
per_pls<-model_performance(explainer_pls, measure = "auc")
per_knn<-model_performance(explainer_knn, measure = "auc")
per_svm<-model_performance(explainer_svm, measure = "auc")
per_rf<-model_performance(explainer_rf, measure = "auc")
per_glm<-model_performance(explainer_glm, measure = "auc")
per_xgboost<-model_performance(explainer_xgboost, measure = "auc")

auc_pls <- auc(per_pls$residuals$observed,per_pls$residuals$predicted)
auc_glm <- auc(per_glm$residuals$observed,per_glm$residuals$predicted)
auc_rf <- auc(per_rf$residuals$observed,per_rf$residuals$predicted)
auc_svm<- auc(per_svm$residuals$observed,per_svm$residuals$predicted)
auc_knn <- auc(per_knn$residuals$observed,per_knn$residuals$predicted)
auc_xgboost <- auc(per_xgboost$residuals$observed,per_xgboost$residuals$predicted)

# ROC曲线
psvm <- plot(per_svm, geom = "roc")+
  annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_svm, 3)))

pknn <- plot(per_knn, geom = "roc")+
  annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_knn, 3)))

prf <- plot(per_rf, geom = "roc")+
  annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_rf, 3)))

pglm <- plot(per_glm, geom = "roc")+
  annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_glm, 3)))

ppls <- plot(per_pls, geom = "roc")+
  annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_pls, 3)))

pxgboost <- plot(per_xgboost, geom = "roc")+
  annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_xgboost, 3)))

# library(cowplot)
# pdf("Train_ML_all_ROC2.pdf",h = 8, w=14)
# plot_grid(pknn, pglm, ppls, prf, psvm, pxgboost, ncol=3)
# dev.off()

# ------------ROC曲线放大一个图------
plotlist <- list(pknn, pglm, ppls, prf, psvm, pxgboost)

method <- c("KNN","GLM","PLS","RF","SVM","XGBoost")
col=c("#EB4B17", "#2775AB", "#91612D",'#4C8045',"#D8D155","#E0867B","#35112D","#E0367A")

auclist <- c(auc_knn,auc_glm,auc_pls,auc_rf,auc_svm,auc_xgboost)

pdf("Train_diagnostic_auc.pdf",w = 6, h = 6)
for (i in 1:length(plotlist)){
  if (i == 1){
    plot(x = plotlist[[i]]$data$fpr, y = plotlist[[i]]$data$tpr, lwd=3, type = "l", col = col[i],
         xlab = "False positive rate", ylab = "True positive rate",
         main = paste0("Gastric Cancer Diagnostic"),
         bty = "l", xaxt = "n")
  } else {
    lines(plotlist[[i]]$data$fpr, y = plotlist[[i]]$data$tpr,
          lwd=3, col = col[i])
  }
}
lines(x=c(0,1),
      y=c(0,1),
      lwd=3, col = "gray")
legend("bottomright",
       legend=c(paste0(method[1]," AUC=",round(auclist[1],3),"\n"),
                paste0(method[2]," AUC=",round(auclist[2],3),"\n"),
                paste0(method[3]," AUC=",round(auclist[3],3),"\n"),
                paste0(method[4]," AUC=",round(auclist[4],3),"\n"),
                paste0(method[5]," AUC=",round(auclist[5],3),"\n"),
                paste0(method[6]," AUC=",round(auclist[6],3)
                )),
       col=col,
       lwd=2,
       title="Curve")
dev.off()

# # # ===============================
# # # ------------------------------------------------------
# # 
# # plot model performance
# pdf("04.machine_learning_residuals_line2.pdf",height = 4,width = 6)
# plot(per_knn,per_glm,per_pls,per_rf,per_svm,per_xgboost)
# dev.off()
# 
# pdf("04.machine_learning_residuals_box2.pdf",height = 4,width = 6)
# plot(per_knn,per_glm,per_pls,per_rf,per_svm,per_xgboost,geom = "boxplot")
# dev.off()
# 
# # importance
# importance_knn<-variable_importance(
#   explainer_knn,
#   loss_function = loss_root_mean_square
# )
# importance_glm<-variable_importance(
#   explainer_glm,
#   loss_function = loss_root_mean_square
# )
# importance_pls<-variable_importance(
#   explainer_pls,
#   loss_function = loss_root_mean_square
# )
# importance_rf<-variable_importance(
#   explainer_rf,
#   loss_function = loss_root_mean_square
# )
# importance_svm<-variable_importance(
#   explainer_svm,
#   loss_function = loss_root_mean_square
# )
# importance_xgboost<-variable_importance(
#   explainer_xgboost,
#   loss_function = loss_root_mean_square
# )
# write.csv(importance_knn,file = "importance_knn.csv")
# write.csv(importance_glm,file = "importance_glm.csv")
# write.csv(importance_pls,file = "importance_pls.csv")
# write.csv(importance_rf,file = "importance_rf.csv")
# write.csv(importance_xgboost,file = "importance_xgboost.csv")
# write.csv(importance_svm,file = "importance_svm.csv")
# 
# pdf("04.machine_learning_importance2.pdf",height = 11,width = 22)
# # par(mfrow = c(2,3), xpd=TRUE)
# # plot(importance_knn,importance_glm,importance_pls,importance_rf,importance_svm,importance_xgboost)
# knn <- plot(importance_knn)
# glm <- plot(importance_glm)
# pls <- plot(importance_pls)
# rf <- plot(importance_rf)
# svm <- plot(importance_svm)
# xgboost <- plot(importance_xgboost)
# # par(2,2,2,2)
# plot_grid(knn, glm, pls, rf, svm, xgboost, ncol=3)
# dev.off()
# 
# # save.image("modeltrained.RData")
# 
# # # 筛选6个模型中RSME值前50%的基因的交集，确定为特征基因----
# # quantile(importance_glm$dropout_loss, 0.5)
# glm <- importance_glm[importance_glm$dropout_loss < quantile(importance_glm$dropout_loss, 0.5),]
# pls <- importance_pls[importance_pls$dropout_loss < quantile(importance_pls$dropout_loss, 0.5),]
# svm <- importance_svm[importance_svm$dropout_loss < quantile(importance_svm$dropout_loss, 0.5),]
# rf <- importance_rf[importance_rf$dropout_loss < quantile(importance_rf$dropout_loss, 0.5),]
# knn <- importance_rf[importance_knn$dropout_loss < quantile(importance_knn$dropout_loss, 0.5),]
# xgboost <- importance_rf[importance_xgboost$dropout_loss < quantile(importance_xgboost$dropout_loss, 0.5),]
# 
# # xgboost <- importance_xgboost[importance_xgboost$dropout_loss < 0.345,]
# 
# glm_gene <- unique(glm$variable) %>% as.data.frame()
# glm_gene <- subset(glm_gene,glm_gene$. != "_full_model_" & glm_gene$. != "group")
# #
# pls_gene <- unique(pls$variable) %>% as.data.frame()
# pls_gene <- subset(pls_gene,pls_gene$. != "_full_model_" & pls_gene$. != "group")
# #
# svm_gene <- unique(svm$variable) %>% as.data.frame()
# svm_gene <-subset(svm_gene,svm_gene$. != "_full_model_" & svm_gene$. != "group")
# 
# rf_gene <- unique(rf$variable) %>% as.data.frame()
# rf_gene <- subset(rf_gene,rf_gene$. != "_full_model_" & rf_gene$. != "group")
# 
# knn_gene <- unique(knn$variable) %>% as.data.frame()
# knn_gene <- subset(knn_gene,knn_gene$. != "_full_model_" & knn_gene$. != "group")
# 
# xgboost_gene <- unique(xgboost$variable) %>% as.data.frame()
# xgboost_gene <- subset(xgboost_gene,xgboost_gene$. != "_full_model_" & xgboost_gene$. != "group")
# 
# # 交集基因
# hub_gene <- Reduce(intersect, list(glm_gene$.,rf_gene$.,pls_gene$.,svm_gene$.,knn_gene$.,xgboost_gene$.)) # 4个
# write.csv(hub_gene,'hub_gene.csv',row.names = T)
# 
# # 绘图upset图
# library(limma)
# library(VennDiagram)
# library(ggVennDiagram)
# library(VennDiagram)
# datalist <- list("XGboost"=xgboost_gene$.,
#                  "RF"=rf_gene$.,
#                  "KNN"=knn_gene$.,
#                  "GLM"=glm_gene$.,
#                  "PLS"=pls_gene$.,
#                  "SVM" = svm_gene$.
# )
# # symbol交集（少于4个时用）
# library(grid)
# library(ggvenn)
# 
# veengene<-hub_gene
# 
# library(readxl)
# Venn <- cbind(unique(xgboost_gene$.),
#               unique(rf_gene$.),
#               unique(knn_gene$.),
#               unique(glm_gene$.),
#               unique(pls_gene$.),
#               unique(svm_gene$.)) %>%
#   as.data.frame()
# colnames(Venn) <- c("XGboost",
#                     "RF",
#                     "KNN",
#                     "GLM",
#                     "PLS",
#                     "SVM")
# 
# Upsetdata <- function(data){
#   require(tidyverse)
#   Genes <- data %>%
#     do.call(c,.) %>%
#     unique()
#   upset0 = matrix(NA, nrow = length(Genes), ncol = length(data))%>%
#     as.data.frame() %>%
#     `rownames<-`(Genes) %>%
#     rownames_to_column("Gene")
#   colnames(upset0)<-c("Gene",names(data))
#   for (i in 2:c(length(data)+1)) {
#     upset0[match(data[[i-1]],upset0[["Gene"]]),i] = 1
#     upset0[is.na(upset0[[i]]),i] <- 0
#   }
#   return(upset0)
# }
# upset_data <- Upsetdata(Venn)
# library(UpSetR)
# p1 <- upset(upset_data,
#             nsets =6,           # 绘制全部13个基因集 默认绘制前5个基因集
#             keep.order=TRUE,     # 保持输入文件列表中基因顺序
#             point.size =2,
#             matrix.color=c("#BC80BD"),  # 点线图颜色
#             line.size =0.3,
#             matrix.dot.alpha = 0.5,     # 非交集点透明度
#             shade.color = "blue",       #背景颜色
#             shade.alpha = 0.1)          #背景透明度
# library(ggplotify)
# require(ggplotify)
# # g1 <- as.ggplot(p1)
# pdf('upsetgenes.pdf',width = 6,height = 6)
# p1
# dev.off()
# 
# # -0---------------------
# # ---------------验证-----------------
# # ====================================

testset <- rbind(data[[2]],
                 data[[3]]
                 # data[[4]]
                 )
data[[5]] <- testset
names(data)[5] <- "Test"

for (i in 2:5){
  datasets <- names(data[i])
  
  testset <- data[[i]]
  
  explainer_knn<-explain(classif_knn,label = "KNN",
                         data = testset,
                         y = testset$group)
  
  explainer_pls<-explain(classif_pls,label = "PLS",
                         data = testset,
                         y = testset$group)
  
  explainer_svm<-explain(classif_svm,label = "SVM",
                         data = testset,
                         y = testset$group)
  
  explainer_rf<-explain(classif_rf,label = "RF",
                        data = testset,
                        y = testset$group)
  
  explainer_glm<-explain(classif_glm,label = "GLM",
                         data = testset,
                         y = testset$group)
  
  # -------------------------------xgboost
  x = model.matrix(group~.,testset)
  
  params <- list(
    objective = "reg:squarederror"
  )
  # 
  # classif_xgboost <- xgb.train(params, data_train, nrounds = 100)
  predict_logit <- function(model,x){
    raw_x <-predict(model,x)
    exp(raw_x)/(1+exp(raw_x))
  }
  
  logit <- function(x){
    exp(x)/(1+exp(x))
  }
  
  explainer_xgboost<-explain(classif_xgboost,
                             label = "xgboost",
                             data = x,
                             y = as.numeric(testset$group),
                             predict_function = predict_logit,
                             link = logit
  )
  
  
  library(pROC)
  per_pls<-model_performance(explainer_pls, measure = "auc")
  per_knn<-model_performance(explainer_knn, measure = "auc")
  per_svm<-model_performance(explainer_svm, measure = "auc")
  per_rf<-model_performance(explainer_rf, measure = "auc")
  per_glm<-model_performance(explainer_glm, measure = "auc")
  per_xgboost<-model_performance(explainer_xgboost, measure = "auc")
  
  auc_pls <- auc(per_pls$residuals$observed,per_pls$residuals$predicted)
  auc_glm <- auc(per_glm$residuals$observed,per_glm$residuals$predicted)
  auc_rf <- auc(per_rf$residuals$observed,per_rf$residuals$predicted)
  auc_svm<- auc(per_svm$residuals$observed,per_svm$residuals$predicted)
  auc_knn <- auc(per_knn$residuals$observed,per_knn$residuals$predicted)
  auc_xgboost <- auc(per_xgboost$residuals$observed,per_xgboost$residuals$predicted)
  
  # ROC曲线
  psvm <- plot(per_svm, geom = "roc")+
    annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_svm, 3)))
  
  pknn <- plot(per_knn, geom = "roc")+
    annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_knn, 3)))
  
  prf <- plot(per_rf, geom = "roc")+
    annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_rf, 3)))
  
  pglm <- plot(per_glm, geom = "roc")+
    annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_glm, 3)))
  
  ppls <- plot(per_pls, geom = "roc")+
    annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_pls, 3)))
  
  pxgboost <- plot(per_xgboost, geom = "roc")+
    annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_xgboost, 3)))
  
  # ---------画图-----------
  plotlist <- list(pknn, pglm, ppls, prf, psvm, pxgboost)
  method <- c("KNN","GLM","PLS","RF","SVM","XGBoost")
  col=c("#EB4B17", "#2775AB", "#91612D",'#4C8045',"#D8D155","#E0867B","#35112D","#E0367A")
  
  auclist <- c(auc_knn,auc_glm,auc_pls,auc_rf,auc_svm,auc_xgboost)
  
  pdf(paste0(datasets,"_diagnostic_auc.pdf"),w = 6, h = 6)
  for (i in 1:length(plotlist)){
    if (i == 1){
      plot(x = plotlist[[i]]$data$fpr, y = plotlist[[i]]$data$tpr, lwd=3, type = "l", col = col[i],
           xlab = "False positive rate", ylab = "True positive rate",
           main = paste0(datasets," Diagnostic"),
           bty = "l", xaxt = "n")
    } else {
      lines(plotlist[[i]]$data$fpr, y = plotlist[[i]]$data$tpr,
            lwd=3, col = col[i])
    }
  }
  lines(x=c(0,1),
        y=c(0,1),
        lwd=3, col = "gray")
  legend("bottomright",
         legend=c(paste0(method[1]," AUC=",round(auclist[1],3),"\n"),
                  paste0(method[2]," AUC=",round(auclist[2],3),"\n"),
                  paste0(method[3]," AUC=",round(auclist[3],3),"\n"),
                  paste0(method[4]," AUC=",round(auclist[4],3),"\n"),
                  paste0(method[5]," AUC=",round(auclist[5],3),"\n"),
                  paste0(method[6]," AUC=",round(auclist[6],3)
                  )),
         col=col,
         lwd=2,
         title="Curve")
  dev.off()
}
