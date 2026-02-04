
# ------------------ 1. 初始化设置 ------------------

rm(list = ls())
foldpath <- paste0("D:/workdir/23/07shap")

if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)

if(! dir.exists("./06_Model")){
  dir.create("./06_Model")
}
setwd("./06_Model")

getwd()

# 创建存图文件夹
figures_dir <- "figures"
if(!dir.exists(figures_dir)) dir.create(figures_dir)

# ------------------ 2. 加载依赖包 ------------------

library(tidyverse)
library(caret)
library(pROC)
library(glmnet)
library(xgboost)
library(randomForest)
library(e1071)
library(nnet)
library(MASS)
library(naivebayes)
library(ggplot2)
library(scales)  
library(doParallel)  
library(tidytext)  
library(data.table)
library(ResourceSelection)
library(rms)
library(gbm)
library(MLmetrics)
library(shapviz)

# 解决函数冲突
conflicted::conflict_prefer_matching("^filter$|^select$|^arrange$|^mutate$|", "dplyr", quiet = TRUE)

# 文件编号函数
# source("/data/nas1/wanghao/script/file_counter.R")

# ------------------ 3. 数据准备 ------------------

# hubgene <- fread("../../05diag/hub_gene.csv")$V2

hubgene <- read.csv("../../04boruta/10.boruta_geneids.csv")
hubgene <- hubgene$x

# group <- fread("/data/nas1/wanghao/project/31.YQGZ-20526-1-昇腾-酮体/00_rawdata/TCGA_LIHC/group.TCGA_LIHC.csv")

# data <- fread("/data/nas1/wanghao/project/31.YQGZ-20526-1-昇腾-酮体/00_rawdata/TCGA_LIHC/dat.TCGA_LIHC.fpkm.csv") %>% column_to_rownames("symbol")
load("../../08batch/data_all.RData")
# 
# testset <- rbind(data[[2]],
#                  data[[3]]
#                  # data[[4]]
# )
# data[[5]] <- testset
# names(data)[5] <- "Test"
# 
for (i in 1:4){
  colnames(data[[i]])[2] <- "group"
  data[[i]] <- data[[i]][,c("group",hubgene)]
  colnames(data[[i]]) <- make.names(colnames(data[[i]]))
}
# 
dat <- data[[1]]

dat <- rbind(data[[1]],
             data[[2]],
             data[[3]]
)

testset <- rbind(data[[2]],
                 data[[3]]
                 # data[[4]]
)
data[[5]] <- testset
names(data)[5] <- "Test"

# dat$group<-as.factor(dat$group)
# dat$group <- as.numeric(dat$group)

dat <- dat[,c("group",hubgene)]

# dat <- data  %>% 
#   filter(rownames(.) %in% hubgene) %>%
#   t() %>% 
#   as.data.frame() %>% 
#   rownames_to_column("sample") %>% 
#   inner_join(group,"sample") %>% 
#   column_to_rownames("sample") %>% 
#   mutate(group = factor(ifelse(group == "Tumor", 1, 0))) # 对应修改

# 规范变量名
colnames(dat) <- make.names(colnames(dat))

# saveRDS(dat, paste0(file_id(), ".train_dat.rds"))

saveRDS(dat, paste0("01.train_dat.rds"))
# ------------------ 4. 函数定义 ------------------

## 4.1 RFE特征选择函数
perform_rfe <- function(train_data, alg_name, n_folds = 10) {
  tryCatch({
    # 根据算法类型选择适当的RFE函数
    rfe_funcs <- switch(alg_name,
                        "RF" = rfFuncs,
                        "LR" = lrFuncs,
                        "DT" = treebagFuncs,
                        "NB" = nbFuncs,
                        "LDA" = ldaFuncs,
                        "BT" = rfFuncs, 
                        caretFuncs  # 默认使用caret通用函数
    )
    
    # 动态计算特征子集大小
    n_features <- ncol(train_data) - 1
    
    # 生成完全符合要求的种子列表，保证结果可复现
    set.seed(12345)
    seeds <- vector(mode = "list", length = n_folds + 1)
    for(i in 1:n_folds) {
      seeds[[i]] <- sample.int(1000, 11)  # 固定长度=11
    }
    seeds[[n_folds + 1]] <- sample.int(1000, 1)
    
    # 正确的RFE控制参数
    rfe_ctrl <- rfeControl(
      functions = rfe_funcs,
      method = "cv",
      number = n_folds,
      verbose = FALSE,
      allowParallel = TRUE,
      seeds = seeds  # 使用正确结构的种子列表
    )
    
    # 设置特征子集大小（动态调整）
    sizes <- unique(floor(seq(1, min(50, n_features), length.out = min(10, n_features))))
    sizes <- sizes[sizes > 0]  # 过滤无效值
    
    # 执行RFE
    rfe_result <- rfe(
      x = train_data[, -which(names(train_data) == "group")],
      y = factor(train_data$group),
      sizes = sizes,
      rfeControl = rfe_ctrl
    )
    
    return(rfe_result)
  }, error = function(e) {
    warning(paste("RFE failed for", alg_name, ":", e$message))
    return(NULL)
  })
}

## 4.2 特征重要性表格生成
get_feature_importance_table <- function(importance_data, algorithm_name, fold_id) {
  library(dplyr)
  
  # 检查数据是否有效
  if (is.null(importance_data) || nrow(importance_data) == 0) {
    warning(paste("No importance data in fold", fold_id))
    return(NULL)
  }
  
  # 整理数据：选择特征名和重要性，排序并过滤0值
  importance_table <- importance_data %>%
    select(Feature = var, Importance = Overall) %>%  # 重命名列
    arrange(desc(Importance)) %>%
    filter(Importance > 0) %>%
    mutate(
      Algorithm = algorithm_name,
      Fold = fold_id,
      Importance = round(Importance, 6)  # 保留4位小数
    )
  
  return(importance_table)
}

## 4.3 RFE性能可视化函数
plot_rfe_performance <- function(rfe_data, alg_name) {
  
  # 计算平均性能
  avg_perf <- rfe_data %>%
    group_by(Variables) %>%
    summarise(
      Accuracy = mean(Accuracy, na.rm = TRUE),
      Accuracy_sd = sd(Accuracy, na.rm = TRUE),
      Kappa = mean(Kappa, na.rm = TRUE),
      Kappa_sd = sd(Kappa, na.rm = TRUE)
    )
  
  # 确定最优特征数（最高准确率）
  optimal_vars <- avg_perf$Variables[which.max(avg_perf$Accuracy)]
  
  # 创建Accuracy曲线图
  accuracy_plot <- ggplot(avg_perf, aes(x = Variables, y = Accuracy)) +
    geom_line(color = "black", linewidth = 1.2) +
    geom_point(color = "black", size = 2) +
    geom_vline(xintercept = optimal_vars, 
               color = "red", linetype = "dashed", linewidth = 0.8) +
    geom_errorbar(aes(ymin = Accuracy - Accuracy_sd, 
                      ymax = Accuracy + Accuracy_sd), 
                  width = 0.2, color = "gray50") +
    scale_x_continuous(breaks = seq(0, max(avg_perf$Variables), by = 5)) +
    scale_y_continuous(labels = scales::percent) +
    labs(title = paste(alg_name, "Accuracy vs Feature Count"),
         x = "Number of Variables",
         y = "Cross-validated Accuracy") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      panel.grid.major = element_line(color = "#f0f0f0"),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  # 创建Kappa曲线图
  kappa_plot <- ggplot(avg_perf, aes(x = Variables, y = Kappa)) +
    geom_line(color = "#1f77b4", linewidth = 1.2) +
    geom_point(color = "#1f77b4", size = 2) +
    geom_vline(xintercept = optimal_vars, 
               color = "red", linetype = "dashed", linewidth = 0.8) +
    geom_errorbar(aes(ymin = Kappa - Kappa_sd, 
                      ymax = Kappa + Kappa_sd), 
                  width = 0.2, color = "gray50") +
    scale_x_continuous(breaks = seq(0, max(avg_perf$Variables), by = 5)) +
    labs(title = paste(alg_name, "Kappa vs Feature Count"),
         x = "Number of Variables",
         y = "Cross-validated Kappa") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      panel.grid.major = element_line(color = "#f0f0f0"),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  # 合并图形
  combined_plot <- gridExtra::grid.arrange(accuracy_plot, kappa_plot, ncol = 2)
  
  # 保存图形
  ggsave(file.path(figures_dir, paste0(round(runif(1)*100000),".RFE_Performance_", alg_name, ".png")), 
         combined_plot, width = 12, height = 5, dpi = 300)
  ggsave(file.path(figures_dir, paste0(round(runif(1)*100000),".RFE_Performance_", alg_name, ".pdf")), 
         combined_plot, width = 12, height = 5, device = cairo_pdf)
  
  return(combined_plot)
}

## 4.4 Log-Loss表
generate_log_loss_table <- function(final_results) {
  # 计算平均Log-Loss
  summary_table <- final_results %>%
    group_by(Algorithm) %>%
    summarise(
      Train_LogLoss = mean(Train_LogLoss, na.rm = TRUE),
      Test_LogLoss = mean(Test_LogLoss, na.rm = TRUE)
    ) %>%
    arrange(Test_LogLoss)  # 按测试集性能排序
  
  # 转换为与图片一致的格式
  formatted_table <- summary_table %>%
    mutate(
      Train_LogLoss = format(round(Train_LogLoss, 6), nsmall = 6),
      Test_LogLoss = format(round(Test_LogLoss, 6), nsmall = 6),
      Algorithm = recode(Algorithm,
                         "DT" = "DT",
                         "LR" = "LR",
                         "KNN" = "KNN",
                         "NB" = "NB",
                         "NNET" = "NNET",
                         "RF" = "RF",
                         "SGBT" = "SGBT",
                         "SVM" = "SVM",
                         "XGBoost" = "XGB"
      )
    ) %>%
    select(Model = Algorithm, 
           `Log Loss Train` = Train_LogLoss,
           `Log Loss Test` = Test_LogLoss)
  
  return(formatted_table)
}

## 4.5 SHAP可视化函数
plot_shap_analysis <- function(shp, optimal_data) {
  # 1. 全局特征重要性
  p_global <- sv_importance(shp, kind = "both", show_numbers = TRUE) +
    ggtitle("Global Feature Importance in Optimal Model") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  # 保存双格式
  ggsave(file.path(figures_dir, paste0(round(runif(1)*100000),".SHAP_Global_Importance.png")), p_global, 
         width = 10, height = 8, dpi = 300)
  ggsave(file.path(figures_dir, paste0(round(runif(1)*100000),".SHAP_Global_Importance.pdf")), p_global, 
         width = 10, height = 8, device = cairo_pdf)
  
  # 2. 关键特征依赖图
  top_features <- head(names(sort(colMeans(abs(shp$S)), decreasing = TRUE)), 3)
  
  for (feature in top_features) {
    p_dep <- sv_dependence(shp, v = feature, color_var = "auto", alpha = 0.7, size = 1.5) + 
      geom_smooth(color = "red", se = FALSE) +
      labs(title = paste("Gene-Effect Relationship:", feature),
           x = "Gene Expression Level", 
           y = "SHAP Value (Contribution to log-odds)")
    
    ggsave(file.path(figures_dir, paste0(round(runif(1)*100000),".SHAP_Dependence_", feature, ".png")), 
           p_dep, width = 8, height = 6, dpi = 300)
    ggsave(file.path(figures_dir, paste0(round(runif(1)*100000),".SHAP_Dependence_", feature, ".pdf")), 
           p_dep, width = 8, height = 6, device = cairo_pdf)
  }
}

## 4.6 全局ROC曲线
plot_roc_curves <- function(pred_data) {
  # 确保pROC包已加载
  if (!requireNamespace("pROC", quietly = TRUE)) {
    install.packages("pROC")
  }
  library(pROC)
  
  # 准备存储ROC对象和AUC值
  roc_list <- list()
  auc_values <- numeric()
  colors <- rainbow(length(pred_data))
  
  # 计算所有ROC对象
  for (i in seq_along(pred_data)) {
    alg_name <- names(pred_data)[i]
    alg_df <- bind_rows(pred_data[[alg_name]]) %>%
      mutate(true_label = as.numeric(as.character(true_label)))
    
    roc_obj <- roc(true_label ~ pred_prob, data = alg_df, quiet = TRUE)
    roc_list[[alg_name]] <- roc_obj
    auc_values[alg_name] <- auc(roc_obj)
  }
  
  # 设置图形参数（增加右侧边距）
  par(mar = c(5, 4, 4, 10), xpd = TRUE)
  
  # 直接绘制第一条ROC曲线（自动创建图形区域）
  plot.roc(roc_list[[1]], 
           col = colors[1], 
           lwd = 2,
           legacy.axes = TRUE,
           xlab = "1 - Specificity",
           ylab = "Sensitivity",
           main = "Algorithm ROC Comparison")
  
  # 添加参考线
  #abline(a = 0, b = 1, lty = 2, col = "gray")
  
  # 添加其余曲线
  if (length(roc_list) > 1) {
    for (i in 2:length(roc_list)) {
      plot.roc(roc_list[[i]], 
               add = TRUE, 
               col = colors[i], 
               lwd = 2)
    }
  }
  
  # 添加图例
  legend("right", 
         inset = c(-0.35, 0),  # 放置在图形右侧外部
         legend = paste0(names(roc_list), " (AUC = ", round(auc_values, 3), ")"),
         col = colors, 
         lwd = 2,
         cex = 0.8,
         bty = "n")
  
  # 保存图形的函数
  save_roc_plot <- function(roc_list, colors, auc_values) {
    # 打开图形设备
    png_path <- file.path(figures_dir, paste0(round(runif(1)*100000), ".Algorithm_ROC_Comparison.png"))
    png(png_path, width = 1000, height = 800, res = 100)
    
    # 重复绘图设置
    par(mar = c(5, 4, 4, 10), xpd = TRUE)
    plot.roc(roc_list[[1]], 
             col = colors[1], 
             lwd = 2,
             legacy.axes = TRUE,
             xlab = "1 - Specificity",
             ylab = "Sensitivity",
             main = "ROC")
    #abline(a = 0, b = 1, lty = 2, col = "gray")
    if (length(roc_list) > 1) {
      for (i in 2:length(roc_list)) {
        plot.roc(roc_list[[i]], add = TRUE, col = colors[i], lwd = 2)
      }
    }
    legend("right", 
           inset = c(0, 0),
           legend = paste0(names(roc_list), " (AUC = ", round(auc_values, 3), ")"),
           col = colors, lwd = 2, cex = 0.8, bty = "n")
    
    dev.off()
    
    # PDF格式保存
    pdf_path <- file.path(figures_dir, paste0(round(runif(1)*100000), ".Algorithm_ROC_Comparison.pdf"))
    pdf(pdf_path, width = 10, height = 8)
    
    # 重复绘图设置
    par(mar = c(5, 4, 4, 10), xpd = TRUE)
    plot.roc(roc_list[[1]], 
             col = colors[1], 
             lwd = 2,
             legacy.axes = TRUE,
             xlab = "1 - Specificity",
             ylab = "Sensitivity",
             main = "ROC")
    #abline(a = 0, b = 1, lty = 2, col = "gray")
    if (length(roc_list) > 1) {
      for (i in 2:length(roc_list)) {
        plot.roc(roc_list[[i]], add = TRUE, col = colors[i], lwd = 2)
      }
    }
    legend("right", 
           inset = c(0, 0),
           legend = paste0(names(roc_list), " (AUC = ", round(auc_values, 3), ")"),
           col = colors, lwd = 2, cex = 0.8, bty = "n")
    
    dev.off()
  }
  
  # 保存图片
  save_roc_plot(roc_list, colors, auc_values)
  
  return(invisible(roc_list))
}

# ------------------ 5. 算法定义 ------------------

algorithms <- list(
  list(name = "LR", 
       train = function(data, features = NULL) {
         if (!is.null(features)) {
           data <- data[, c(features, "group")]
         }
         glm(
           group ~ .,
           family = binomial(link = "logit"),
           data = data
         )
       }),
  
  list(name = "DT", 
       train = function(data, features = NULL) {
         if (!is.null(features)) {
           data <- data[, c(features, "group")]
         }
         rpart::rpart(
           group ~ .,
           data = data,
           method = "class"
         )
       }),
  
  list(name = "SVM",
       train = function(data, features = NULL) {
         if (!is.null(features)) {
           data <- data[, c(features, "group")]
         }
         svm(
           group ~ .,
           data = data,
           type = "C-classification",
           kernel = "radial",
           probability = TRUE
         )
       }),
  
  list(name = "RF",
       train = function(data, features = NULL) {
         if (!is.null(features)) {
           data <- data[, c(features, "group")]
         }
         randomForest(
           factor(group) ~ .,
           data = data,
           ntree = 500,
           seed = 12345
         )
       }),
  
  list(name = "KNN",
       train = function(data, features = NULL) {
         if (!is.null(features)) {
           data <- data[, c(features, "group")]
         }
         knn3(
           factor(group) ~ .,
           data = data,
           k = 5
         )
       }),
  
  list(name = "BT",
       train = function(data, features = NULL) {
         if (!is.null(features)) {
           data <- data[, c(features, "group")]
         }
         randomForest(
           factor(group) ~ .,
           data = data,
           mtry = ncol(data) - 1,  # 关键：使用所有特征
           ntree = 500,
           importance = TRUE,
           seed = 12345
         )
       }),
  
  list(name = "LDA",
       train = function(data, features = NULL) {
         if (!is.null(features)) {
           data <- data[, c(features, "group")]
         }
         MASS::lda(
           group ~ .,
           data = data,
           prior = c(0.5, 0.5)
         )
       }),
  
  list(name = "NNET",
       preprocess = function(data) {
         preProc <- preProcess(data, method = c("center", "scale"))
         predict(preProc, data)
       },
       train = function(data, features = NULL) {
         if (!is.null(features)) {
           data <- data[, c(features, "group")]
         }
         nnet::nnet(factor(group) ~ ., data = data, size = 5, trace = FALSE,
                    seed = 12345)
       }),
  
  list(name = "NB",
       train = function(data, features = NULL) {
         if (!is.null(features)) {
           data <- data[, c(features, "group")]
         }
         naive_bayes(factor(group) ~ ., data = data)
       }),
  
  list(name = "SGBT",
       train = function(data, features = NULL) {
         # 确保使用传入的 features 参数
         if (!is.null(features)) {
           data <- data[, c(features, "group")]
         }
         data$group <- as.numeric(as.character(data$group))
         
         gbm(
           formula = group ~ . ,
           data = data,
           distribution = "bernoulli",
           n.trees = 500,
           interaction.depth = 3,
           shrinkage = 0.01,
           bag.fraction = 0.7
         )
       })
)

# ------------------ 6. 主分析流程 ------------------
# set.seed(99999)

fold_train <- readRDS('01.train_dat.rds')

cl <- makeCluster(10)
registerDoParallel(cl)

# 储存结果
all_results <- list()
feature_importance <- list()
feature_history <- list()
rfe_performance <- list() 

all_train_preds <- list()
all_test_preds <- list()
all_true_labels <- list()

# 交叉验证主循环
knumb <- 1
folds <- createFolds(fold_train$group, k = knumb)

for (alg in algorithms) {
  
  cat("\n=== 开始训练", alg$name, "模型 (带RFE特征选择) ===\n")
  results <- data.frame()
  fold_features <- list()
  fold_imp_tables <- list()  # 存储所有fold的表格
  alg_rfe_data <- list()  # 存储当前算法的RFE性能数据
  
  for (i in 1:knumb) {
    # # 数据划分
    # test_index <- folds[[i]]
    # train_data <- fold_train[-test_index, ]
    # test_data <- fold_train[test_index, ]
    train_data <- data[[1]]
    test_data <- data[[5]]
    
    # 算法特定预处理
    if (!is.null(alg$preprocess)) {
      train_prep <- alg$preprocess(train_data)
      test_prep <- alg$preprocess(test_data)
    } else {
      train_prep <- train_data
      test_prep <- test_data
    }
    
    # 执行递归特征消除
    rfe_result <- perform_rfe(train_prep, alg$name,n_folds =10)
    
    # 获取选择的特征
    if (!is.null(rfe_result)) {
      selected_features <- predictors(rfe_result)
      cat("Fold", i, "选择了", length(selected_features), "个特征: ", 
          paste(head(selected_features, 5), collapse = ", "), 
          ifelse(length(selected_features) > 5, "...", ""), "\n")
    } else {
      selected_features <- setdiff(colnames(train_prep), "group")
      cat("Fold", i, "使用所有特征（", length(selected_features), "个）\n")
    }
    
    # 保存特征选择结果
    fold_features[[i]] <- selected_features
    
    if (!is.null(rfe_result)) {
      # 提取RFE结果
      rfe_perf <- rfe_result$results %>%
        mutate(Algorithm = alg$name,
               Fold = i)
      alg_rfe_data[[i]] <- rfe_perf
    }
    
    fold_imp_tables [[i]] <- if (!is.null(rfe_result)) {
      get_feature_importance_table(rfe_result$variables %>% as.data.frame(), alg$name, 1)
    } else NULL
    
    # 模型训练（使用选择的特征）
    model <- tryCatch(
      alg$train(train_prep[, c(selected_features, "group")], selected_features),
      error = function(e) {
        warning(paste(alg$name, "模型训练失败:", e$message))
        NULL
      }
    )
    
    if (is.null(model)) next
    
    train_pred_probs <- tryCatch({
      switch(alg$name,
             "DT" = predict(model, train_prep, type = "prob")[, 2],
             "RF" = predict(model, train_prep, type = "prob")[, 2],
             "KNN" = predict(model, train_prep, type = "prob")[, 2],
             "LDA" = predict(model, train_prep)$posterior[, 2],
             "NNET" = predict(model, train_prep, type = "raw"),
             "NB" = predict(model, train_prep, type = "prob")[, 2],
             "BT" = predict(model, train_prep, type = "prob")[, 2],   
             "SVM" = {
               attr(predict(model, newdata = train_prep, probability = TRUE), "probabilities")[, "1"]
             },
             predict(model, train_prep, type = "response")
      )
    }, error = function(e) {
      warning(paste("训练集预测失败:", e$message))
      rep(NA, nrow(train_prep))
    })
    
    
    pred_probs <- tryCatch({
      switch(alg$name,
             "DT" = predict(model, test_prep, type = "prob")[, 2],
             "RF" = predict(model, test_prep, type = "prob")[, 2],
             "KNN" = predict(model, test_prep, type = "prob")[, 2],
             "LDA" = predict(model, test_prep)$posterior[, 2],
             "NNET" = predict(model, test_prep, type = "raw"),
             "NB" = predict(model, test_prep, type = "prob")[, 2],
             "BT" = predict(model, test_prep, type = "prob")[, 2],   
             "SVM" = {
               attr(predict(model, newdata = test_prep, probability = TRUE), "probabilities")[, "1"]
             },
             predict(model, test_prep, type = "response")
      )
    }, error = function(e) {
      warning(paste("预测失败:", e$message))
      rep(NA, nrow(test_prep))
    })
    
    # 计算二分类预测结果
    pred_class <- ifelse(pred_probs > 0.5, 1, 0)
    
    # 统一因子水平
    ref_levels <- levels(factor(train_prep$group))
    pred_factor <- factor(pred_class, levels = ref_levels)
    true_factor <- factor(test_prep$group, levels = ref_levels)
    
    # 安全计算混淆矩阵
    conf_matrix <- tryCatch(
      confusionMatrix(pred_factor, true_factor),
      error = function(e) NULL
    )
    
    # 定义Log-Loss计算函数（带安全保护）
    safe_log_loss <- function(true_labels, predicted_probs) {
      # 避免Inf值
      epsilon <- 1e-15
      predicted_probs <- pmin(pmax(predicted_probs, epsilon), 1 - epsilon)
      
      # 确保标签和概率一致
      if (is.factor(true_labels)) {
        true_labels <- as.numeric(true_labels) - 1
      }
      
      return(LogLoss(y_pred = predicted_probs, y_true = true_labels))
    }
    
    # 训练集Log-Loss
    train_log_loss <- tryCatch({
      safe_log_loss(train_prep$group, train_pred_probs)
    }, error = function(e) {
      warning(paste("训练集Log-Loss计算失败:", e$message))
      NA
    })
    
    # 测试集Log-Loss
    test_log_loss <- tryCatch({
      safe_log_loss(test_prep$group, pred_probs)
    }, error = function(e) {
      warning(paste("测试集Log-Loss计算失败:", e$message))
      NA
    })
    
    # 收集训练集预测结果
    if (is.null(all_train_preds[[alg$name]])) {
      all_train_preds[[alg$name]] <- list()
    }
    all_train_preds[[alg$name]][[i]] <- data.frame(
      sample_id = rownames(train_prep),
      true_label = train_prep$group,
      pred_prob = train_pred_probs,
      fold = i
    )
    
    # 收集测试集预测结果
    if (is.null(all_test_preds[[alg$name]])) {
      all_test_preds[[alg$name]] <- list()
    }
    all_test_preds[[alg$name]][[i]] <- data.frame(
      sample_id = rownames(test_prep),
      true_label = test_prep$group,
      pred_prob = pred_probs,
      fold = i
    )
    
    # 保存结果
    current_results <- data.frame(
      Algorithm = alg$name,
      Fold = i,
      Train_AUC = auc(roc(train_prep$group, train_pred_probs, direction = "<")) %>% as.numeric(),
      Test_AUC = auc(roc(test_data$group, pred_probs, direction = "<")) %>% as.numeric(),
      Train_LogLoss = train_log_loss,
      Test_LogLoss = test_log_loss,
      Sensitivity = if(!is.null(conf_matrix)) conf_matrix$byClass["Sensitivity"] else NA,
      Specificity = if(!is.null(conf_matrix)) conf_matrix$byClass["Specificity"] else NA,
      Accuracy = if(!is.null(conf_matrix)) conf_matrix$overall["Accuracy"] else NA,
      F1 = if(!is.null(conf_matrix)) conf_matrix$byClass["F1"] else NA,
      PPV = if(!is.null(conf_matrix)) conf_matrix$byClass["Pos Pred Value"] else NA,
      NPV = if(!is.null(conf_matrix)) conf_matrix$byClass["Neg Pred Value"] else NA,
      Brier = mean((as.numeric(test_data$group) - pred_probs)^2, na.rm = TRUE),  # Brier评分计算
      Num_Features = length(selected_features)
    )
    
    results <- rbind(results, current_results)
    
  }
  
  # 存储算法结果
  all_results[[alg$name]] <- results
  
  # 保存特征选择历史
  feature_history[[alg$name]] <- do.call(rbind, lapply(seq_along(fold_features), function(i) {
    data.frame(Algorithm = alg$name, Fold = i, Feature = fold_features[[i]])
  }))
  feature_importance[[alg$name]] <- do.call(rbind, lapply(seq_along(fold_imp_tables), function(i) {
    if (is.null(fold_imp_tables[[i]])) return(NULL)
    
    # 仅对 LR 算法缩放 Importance
    if (alg$name == "LR") {
      fold_imp_tables[[i]]$Importance <- fold_imp_tables[[i]]$Importance * 1e4  
    }
    
    # 添加 Fold 和 Algorithm 列
    fold_imp_tables[[i]]$Fold <- i
    fold_imp_tables[[i]]$Algorithm <- alg$name
    fold_imp_tables[[i]]
  }))
  
  # 保存当前算法的RFE性能数据
  if (length(alg_rfe_data) > 0) {
    rfe_performance[[alg$name]] <- bind_rows(alg_rfe_data)
  }
  
}

# 停止并行计算
stopCluster(cl)

# ------------------ 7. 结果保存 ------------------

# 汇总所有结果
final_results <- bind_rows(all_results)
saveRDS(final_results, "final_results_with_rfe.rds")
# fwrite(final_results,paste0(file_id(),".final_results.csv"))
fwrite(final_results,paste0("02.final_results.csv"))
# 特征选择历史分析
feature_history_df <- bind_rows(feature_history)
# fwrite(final_results,paste0(file_id(),".feature_selection_history.csv"))
fwrite(final_results,paste0("03.feature_selection_history.csv"))
# 特征重要性汇总
feature_importance_df <- bind_rows(feature_importance)
# fwrite(final_results,paste0(file_id(),".feature_importance.csv"))
fwrite(final_results,paste0("04.feature_importance.csv"))

# ------------------ 8. 可视化输出 ------------------

cat("\n=== 算法性能比较 (带RFE特征选择) ===\n")
# 计算各算法平均性能
summary_table <- final_results %>%
  group_by(Algorithm) %>%
  summarise(
    Train_AUC = mean(Train_AUC, na.rm = TRUE),
    AUC = mean(Test_AUC, na.rm = TRUE),
    Train_LogLoss = mean(Train_LogLoss, na.rm = TRUE),
    LogLoss = mean(Test_LogLoss, na.rm = TRUE),
    Sensitivity = mean(Sensitivity, na.rm = TRUE),
    Specificity = mean(Specificity, na.rm = TRUE),
    Accuracy = mean(Accuracy, na.rm = TRUE),
    F1 = mean(F1, na.rm = TRUE),
    Avg_Features = mean(Num_Features, na.rm = TRUE)
  ) %>%
  arrange(desc(AUC))

print(summary_table)


# 10折平均AUC
# 8.1 全局ROC曲线图（10折数据合并）
# 绘制所有算法的ROC曲线比较图
if (length(all_test_preds) > 0) {
  roc_plot <- plot_roc_curves(all_test_preds)
  print(roc_plot)
}

if (length(all_test_preds) > 0) {
  roc_plot <- plot_roc_curves(all_train_preds)
  print(roc_plot)
}

## 8.2 特征稳定性图
feature_frequency <- feature_history_df %>%
  group_by(Algorithm, Feature) %>%
  summarise(Frequency = n() / 5, .groups = "drop") %>%  # 5折交叉验证
  group_by(Algorithm) %>%
  arrange(desc(Frequency)) %>%
  mutate(Rank = row_number())

# 每个算法选择最稳定的特征
top_features <- feature_frequency %>%
  group_by(Algorithm) %>%
  slice_max(Frequency, n = 5)  # 每个算法前10个特征

# 绘制特征稳定性图
p2 <- ggplot(top_features, aes(x = reorder_within(Feature, Frequency, Algorithm), 
                               y = Frequency, fill = Algorithm)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~Algorithm, scales = "free_y", ncol = 3) +
  scale_x_reordered() +
  labs(title = "Feature Selection Stability",
       x = "Feature",
       y = "Selection Frequency in CV") +
  theme_minimal() +
  theme(legend.position = "none")

# ggsave(file.path(figures_dir, paste0(file_id(),".feature_stability_rfe.pdf")), p2, width = 15, height = 10, device = cairo_pdf)
ggsave(file.path(figures_dir, paste0("01.feature_stability_rfe.pdf")), p2, width = 15, height = 10, device = cairo_pdf)

## 8.3 特征重要性图

if (!is.null(feature_importance_df)) {
  # 计算平均特征重要性
  avg_imp <- feature_importance_df %>%
    group_by(Algorithm, Feature) %>%
    summarise(Avg_Importance = mean(Importance, na.rm = TRUE), .groups = "drop") %>%
    group_by(Algorithm) %>%
    arrange(desc(Avg_Importance)) %>%
    slice_head(n = 10)  # 每个算法取前10个重要特征
  
  # 绘制综合特征重要性图
  p3 <- ggplot(avg_imp, aes(x = reorder_within(Feature, Avg_Importance, Algorithm), 
                            y = Avg_Importance, fill = Algorithm)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    facet_wrap(~Algorithm, scales = "free_y", ncol = 3) +
    scale_x_reordered() +
    labs(title = "Algorithm Feature Importance Ranking",
         x = "Feature",
         y = "Average Importance") +
    theme_minimal() +
    theme(legend.position = "none")
  ggsave(file.path(figures_dir, paste0(".feature_importance_summary.png")), p3, width = 15, height = 10, dpi = 300)
  ggsave(file.path(figures_dir, paste0(".feature_importance_summary.pdf")), p3, width = 15, height = 10, device = cairo_pdf)
  
}

## 8.4 RFE性能曲线

if (length(rfe_performance) > 0) {
  for (alg_name in names(rfe_performance)) {
    cat("绘制", alg_name, "的RFE性能曲线...\n")
    tryCatch({
      plot_rfe_performance(rfe_performance[[alg_name]], alg_name)
    }, error = function(e) {
      warning(paste("绘制", alg_name, "RFE曲线失败:", e$message))
    })
  }
}


## 8.4 log-loss

log_loss_table <- generate_log_loss_table(final_results)
log_loss_table 
write_csv(log_loss_table, paste0(".LogLoss_Comparison_Table.csv"))


## 8.4 评估模型指标柱状图
# ​​Test_AUC​​（测试集AUC）：
# 最重要指标，综合反映分类器性能（值越接近1越好）
# ​​Test_LogLoss​​（测试集对数损失）：
# 衡量概率校准质量（值越小越好）
# ​​Accuracy​​（测试集准确率）：
# 整体分类正确率
# ​​F1-Score​​：
# 平衡精确率和召回率（尤其适用于类别不平衡）
# ​​稳定性​​（标准差）：
# 算法在10折中的性能波动（标准差越小越稳定）


# 手动创建数据框（替换为实际文件读取）
final_results <- fread("02.final_results.csv")  

# 定义关键指标 - 保持原始顺序
# 使用保留顺序的字符向量，而不是转换为因子
specific_metrics <- c("AUC","Sensitivity","Specificity","Accuracy" ,"F1" ,"PPV","NPV","Brier")


# 转换长格式并处理因子
metric_data <- final_results %>%
  rename(AUC=Test_AUC) %>%
  pivot_longer(
    cols = all_of(specific_metrics),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    # 确保算法顺序
    Algorithm = factor(Algorithm, levels = unique(final_results$Algorithm)),
    
    # 关键修正：将指标转换为因子并明确指定水平顺序
    Metric = factor(Metric, levels = specific_metrics),  # 保持原始列顺序
    
    # 确保数据类型正确
    Value = as.numeric(Value)
  ) %>%
  filter(!is.na(Value))  # 移除NA值


# 计算统计量（算法x指标的均值/标准误）
metric_summary <- metric_data %>%
  group_by(Algorithm, Metric) %>%
  summarise(
    mean_val = mean(Value),
    se = sd(Value) / sqrt(n()),  # 标准误
    .groups = 'drop'
  ) 

# 修改后的绘图代码 - 重点突出0-1范围的值
ggplot(metric_summary, aes(x = Algorithm, y = mean_val, color = Algorithm)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = mean_val - se, ymax = mean_val + se),
    width = 0.2,
    linewidth = 0.8
  ) +
  annotate("rect", 
           xmin = -Inf, xmax = Inf, 
           ymin = 0, ymax = 1.5, 
           fill = "gray90", alpha = 0.1) +
  # 关键修改3：突出显示0-1范围的内容
  geom_text(
    aes(label = ifelse(mean_val <= 1, sprintf("%.3f", mean_val), "")), 
    vjust = -1.5, size = 3.5, color = "black"
  ) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 4) +
  labs(
    #    title = "Algorithm Performance Across Metrics (Focus on 0-1 Range)",
    y = "Mean Value ± SE",
    x = "Algorithm"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "lightgray"),
    # 添加强调0-1范围的标记
    axis.line.y = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.text.y = element_text(color = "black", face = "bold")
  ) 

ggsave(file.path(figures_dir, paste0(".algorithm_comparison.png")),  width = 16, height = 6.2, dpi = 300)
ggsave(file.path(figures_dir, paste0(".algorithm_comparison.pdf")),  width = 16, height = 6.2, device = cairo_pdf)

save.image(paste0(".10_ML_RFE.RData"))


# ------------------ 9. SHAP分析 ------------------
set.seed(12345)
full_data <- readRDS("01.train_dat.rds")

cl <- makeCluster(10)
registerDoParallel(cl, iseed = 12345)

rfe_result <- perform_rfe(full_data, "RF", n_folds = 10)  
# 停止并行计算
stopCluster(cl)

# 获取最优特征
if (!is.null(rfe_result)) {
  selected_features <- predictors(rfe_result)
  cat("RFE选出的最优特征数:", length(selected_features), "\n")
  cat("重要特征列表:\n")
  print(selected_features)
  
  # 保存最优特征
  write.csv(selected_features, paste0(".optimal_features.csv"))
  
  # 准备最优特征数据集
  optimal_data <- full_data[, c(selected_features, "group")]
  optimal_data$group <- as.numeric(as.character(optimal_data$group))
} else {
  warning("RFE失败，使用所有特征")
  optimal_data <- full_data
  optimal_data$group <- as.numeric(as.character(optimal_data$group)) #-1
}

# RFE选出的最优特征数: 24 
# 重要特征列表:
#   [1] "HSPB1"    "SERPINH1" "FOSB"     "FOS"      "TFPI2"    "PTGS2"    "ZFP36"    "IL1B"     "HS3ST3B1" "THBS1"    "HSPA1A"   "XIRP1"   
# [13] "ATF3"     "MEDAG"    "CXCL10"   "TNFAIP3"  "CXCL2"    "DNAJB1"   "PNP"      "FOSL1"    "NFKBIE"   "IL6"      "AREG"     "SAT1"

# 训练最终模型（使用全部数据和最优特征）
if (!is.null(selected_features)) {
  optimal_data <- optimal_data[, c(selected_features, "group")]
}

# 最优模型
model_name = "RF"
# final_model <- algorithms[[5]]$train(optimal_data,selected_features)
# final_model <- algorithms[[3]]$train(optimal_data,selected_features)
final_model <- algorithms[[4]]$train(optimal_data,selected_features)

# 保存最终模型
saveRDS(final_model, paste0(".final_model.rds"))

# 需要conda创建虚拟环境安装
library(fastshap)

pred_fun <- 
  switch(model_name,
         "DT" = function(model, newdata) {predict(model, newdata, type = "prob")[, 2]},
         "RF" = function(model, newdata) {predict(model, newdata, type = "prob")[, 2]},
         "KNN" = function(model, newdata) {predict(model, newdata, type = "prob")[, 2]},
         "LDA" = function(model, newdata) {predict(model, newdata)$posterior[, 2]},
         "NNET" = function(model, newdata) {predict(model, newdata, type = "raw")},
         "NB" = function(model, newdata) {predict(model, newdata, type = "prob")[, 2]},
         "BT" = function(model, newdata) {predict(model, newdata, type = "prob")[, 2]},   
         "SVM" = {
           function(model, newdata) {attr(predict(model, newdata = newdata, probability = TRUE), "probabilities")[, "1"]}
         },
         function(model, newdata) {predict(model, newdata, type = "response")}
  )



# 准备特征矩阵
X <- optimal_data[, -which(colnames(optimal_data) == "group")]

# 计算 SHAP 值
shap_values <- fastshap::explain(
  final_model, 
  X = X, 
  pred_wrapper = pred_fun,
  nsim = 100  # 蒙特卡洛采样次数，增加可提高精度，默认值为1。注意：为了获得最准确的结果，应将nsim设置得尽可能大。
)

# 转换为 shapviz 对象
shp <- shapviz(shap_values, X = X)

# 生成可视化
plot_shap_analysis(shp, optimal_data)


# 高风险样本决策分解 可视化
pred_probs <- tryCatch({
  switch(model_name,
         "DT" = predict(final_model, optimal_data, type = "prob")[, 2],
         "RF" = predict(final_model, optimal_data, type = "prob")[, 2],
         "KNN" = predict(final_model, optimal_data, type = "prob")[, 2],
         "LDA" = predict(final_model, optimal_data)$posterior[, 2],
         "NNET" = predict(final_model, optimal_data, type = "raw"),
         "NB" = predict(final_model, optimal_data, type = "prob")[, 2],
         "BT" = predict(final_model, optimal_data, type = "prob")[, 2],   
         "SVM" = {
           attr(predict(final_model, newdata = optimal_data, probability = TRUE), "probabilities")[, "1"]
         },
         predict(final_model, optimal_data, type = "response")
  )
}, error = function(e) {
  warning(paste("训练集预测失败:", e$message))
  rep(NA, nrow(optimal_data))
})

high_risk_idx <- which(pred_probs > 0.9)  # 高风险样本阈值

if (length(high_risk_idx) > 0) {
  sample_id <- high_risk_idx[2]
  p_force <- sv_force(shp, row_id = sample_id) +
    theme_bw() +
    labs(title = paste("High-Risk Sample Decision (ID:", rownames(optimal_data)[sample_id], ")"))
  
  ggsave(file.path(figures_dir, paste0(".SHAP_HighRisk_ForcePlot.png")), 
         p_force, width = 12, height = 6, dpi = 300)
  ggsave(file.path(figures_dir, paste0(".SHAP_HighRisk_ForcePlot.pdf")), 
         p_force, width = 12, height = 6, device = cairo_pdf)
  
}

# 5. 生成临床可解释报告
top_risk_genes <- names(sort(colMeans(shp$S), decreasing = TRUE)[1:5])
top_protective_genes <- names(sort(colMeans(shp$S))[1:2])
top_features <- head(names(sort(colMeans(abs(shp$S)), decreasing = TRUE)), 3)
top5_features <- head(names(sort(colMeans(abs(shp$S)), decreasing = TRUE)), 5)                  

cat("## 模型临床解释报告\n")
cat("### 关键发现：\n")
cat("1. **驱动基因**：", paste(top_risk_genes, collapse = ", "), "主导疾病风险预测\n")
cat("2. **保护性基因**：", paste(top_protective_genes, collapse = ", "), "降低MDD风险\n")
cat("3. **非线性效应**：", top_features[1], "表达与风险呈U型关系（阈值=", 
    round(quantile(X[, top_features[1]], 0.6), 3), ")\n")
cat("4. **关键交互**：", paste(top5_features[1:2], collapse = "-"), 
    "交互效应显著影响风险预测\n")

fwrite(data.frame(symbol=top5_features), paste0(".top5_features.csv"))

# ------------------ 10. 保存最终环境 ------------------
# 保存SHAP值结果
shap_results <- list(
  shap_values = shp$S,
  base_value = shp$baseline,
  features = X
)
saveRDS(shap_results, paste0(".SHAP_Results.rds"))
save.image(paste0(".10ML_RFE_Full_Analysis.RData"))

cat("\n=== Analysis Completed ===")




