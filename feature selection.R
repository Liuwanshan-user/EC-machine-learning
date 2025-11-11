# ===================================================================
# 蛋白组学数据特征选择和多算法模型评估 - 完整版
# ===================================================================

# 检查并安装所需的包
packages_needed <- c(
  "ggplot2", "cowplot", "dplyr", "pROC", "glmnet", 
  "randomForest", "e1071", "nnet", "caret", "PRROC",
  "scales", "readr", "gridExtra"
)

new_packages <- packages_needed[!(packages_needed %in% installed.packages()[, "Package"])]
if (length(new_packages)) {
  install.packages(new_packages)
}

# 加载所需包
library(ggplot2)
library(cowplot)
library(dplyr)
library(pROC)
library(glmnet)
library(randomForest)
library(e1071)
library(nnet)
library(caret)
library(PRROC)
library(scales)
library(readr)
library(gridExtra)

# ------------------------
# 统一设置字体和大小
# ------------------------
if (.Platform$OS.type == "windows") {
  windowsFonts(
    Arial = windowsFont("Arial"),
    Times = windowsFont("Times New Roman"),
    Helvetica = windowsFont("Helvetica")
  )
}

font_family <- "Arial"
axis_text_size <- 10
font_size <- 10
line_width <- 1.0
plot_width <- 7.5
plot_height <- 7

# 设置随机种子
set.seed(123)

# 设置工作路径
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ===================================================================
# 1. 数据读取和预处理
# ===================================================================

cat("==================== 数据读取 ====================\n")
data <- read.csv("significant_metabolite_data.csv", check.names = FALSE)

# 提取训练集和测试集
train_data <- data %>% filter(cohort == "train")
test_data <- data %>% filter(cohort == "test")

# 提取特征和标签
train_x <- train_data %>% select(-c(id, label, cohort))
train_y <- as.numeric(as.character(train_data$label))
test_x <- test_data %>% select(-c(id, label, cohort))
test_y <- as.numeric(as.character(test_data$label))

cat(sprintf("训练集样本数: %d\n", nrow(train_data)))
cat(sprintf("测试集样本数: %d\n", nrow(test_data)))
cat(sprintf("特征数: %d\n", ncol(train_x)))
cat(sprintf("训练集标签分布: 0=%d, 1=%d\n", sum(train_y==0), sum(train_y==1)))
cat(sprintf("测试集标签分布: 0=%d, 1=%d\n", sum(test_y==0), sum(test_y==1)))

# Log2转换和标准化
cat("\n==================== 数据预处理 ====================\n")
train_x_log <- log2(train_x + 1)
test_x_log <- log2(test_x + 1)

# 标准化 (使用训练集的均值和标准差)
train_mean <- apply(train_x_log, 2, mean)
train_sd <- apply(train_x_log, 2, sd)

train_x_scaled <- scale(train_x_log, center = train_mean, scale = train_sd)
test_x_scaled <- scale(test_x_log, center = train_mean, scale = train_sd)

cat("预处理完成: Log2转换 + Z-score标准化\n")

# ===================================================================
# 2. 特征选择前的性能评估
# ===================================================================

cat("\n==================== 特征选择前性能评估 ====================\n")

# 定义配色方案
model_colors <- c(
  "#3AB5B3",  # 青色 - MLP
  "#7B6C9F",  # 紫色 - LASSO
  "#A188BD",  # 浅紫 - Logistic
  "#BBC5DE",  # 浅蓝 - Random Forest
  "#E7777F"   # 粉色 - SVM
)

model_names <- c("MLP", "LASSO", "Logistic", "Random Forest", "SVM")

# 简化的预测函数用于快速评估
quick_predict_cv <- function(X, y, model_type, n_folds = 5) {
  folds <- createFolds(as.factor(y), k = n_folds)
  cv_probs <- numeric(length(y))
  
  for (i in 1:n_folds) {
    test_idx <- folds[[i]]
    train_idx <- -test_idx
    
    X_train_cv <- X[train_idx, , drop = FALSE]
    y_train_cv <- y[train_idx]
    X_test_cv <- X[test_idx, , drop = FALSE]
    
    if (model_type == "LASSO") {
      model <- cv.glmnet(X_train_cv, y_train_cv, family = "binomial", alpha = 1)
      cv_probs[test_idx] <- predict(model, X_test_cv, s = "lambda.min", type = "response")[,1]
    } else if (model_type == "Random Forest") {
      model <- randomForest(X_train_cv, as.factor(y_train_cv), ntree = 500)
      cv_probs[test_idx] <- predict(model, X_test_cv, type = "prob")[, 2]
    }
  }
  return(cv_probs)
}

quick_predict <- function(X_train, y_train, X_test, model_type) {
  if (model_type == "LASSO") {
    model <- cv.glmnet(X_train, y_train, family = "binomial", alpha = 1)
    probs <- predict(model, X_test, s = "lambda.min", type = "response")[,1]
  } else if (model_type == "Random Forest") {
    model <- randomForest(X_train, as.factor(y_train), ntree = 500)
    probs <- predict(model, X_test, type = "prob")[, 2]
  }
  return(probs)
}

# 计算特征选择前的性能
cat("计算特征选择前的性能...\n")
before_train_lasso <- quick_predict_cv(train_x_scaled, train_y, "LASSO")
before_test_lasso <- quick_predict(train_x_scaled, train_y, test_x_scaled, "LASSO")
before_train_rf <- quick_predict_cv(train_x_scaled, train_y, "Random Forest")
before_test_rf <- quick_predict(train_x_scaled, train_y, test_x_scaled, "Random Forest")

# 计算ROC
roc_before_train_lasso <- roc(train_y, before_train_lasso, quiet = TRUE)
roc_before_test_lasso <- roc(test_y, before_test_lasso, quiet = TRUE)
roc_before_train_rf <- roc(train_y, before_train_rf, quiet = TRUE)
roc_before_test_rf <- roc(test_y, before_test_rf, quiet = TRUE)

# 计算PR
pr_before_train_lasso <- pr.curve(
  scores.class0 = before_train_lasso[train_y == 1],
  scores.class1 = before_train_lasso[train_y == 0],
  curve = TRUE
)
pr_before_test_lasso <- pr.curve(
  scores.class0 = before_test_lasso[test_y == 1],
  scores.class1 = before_test_lasso[test_y == 0],
  curve = TRUE
)
pr_before_train_rf <- pr.curve(
  scores.class0 = before_train_rf[train_y == 1],
  scores.class1 = before_train_rf[train_y == 0],
  curve = TRUE
)
pr_before_test_rf <- pr.curve(
  scores.class0 = before_test_rf[test_y == 1],
  scores.class1 = before_test_rf[test_y == 0],
  curve = TRUE
)

# ===================================================================
# 3. 特征选择: LASSO + 随机森林交集
# ===================================================================

cat("\n==================== 特征选择 ====================\n")

# LASSO特征选择
cv_lasso <- cv.glmnet(
  as.matrix(train_x_scaled), 
  train_y, 
  family = "binomial", 
  alpha = 1,
  nfolds = 5
)

lasso_coef <- coef(cv_lasso, s = "lambda.min")
lasso_features <- rownames(lasso_coef)[which(lasso_coef != 0)][-1]

cat(sprintf("LASSO选择特征数: %d\n", length(lasso_features)))

# 随机森林特征重要性
rf_model <- randomForest(
  x = train_x_scaled, 
  y = as.factor(train_y), 
  importance = TRUE,
  ntree = 500
)

rf_importance <- importance(rf_model)
rf_importance_sorted <- sort(rf_importance[, "MeanDecreaseGini"], decreasing = TRUE)
top_n <- min(30, length(rf_importance_sorted))
rf_features <- names(rf_importance_sorted)[1:top_n]

cat(sprintf("随机森林TOP%d特征\n", top_n))

# 取交集
selected_features <- intersect(lasso_features, rf_features)
cat(sprintf("\n特征交集数量: %d\n", length(selected_features)))

if (length(selected_features) == 0) {
  cat("警告: 交集为空,使用LASSO特征\n")
  selected_features <- lasso_features
}

# 保存选择的特征
write.csv(
  data.frame(Feature = selected_features),
  "selected_features.csv",
  row.names = FALSE
)

# 更新数据集为选择的特征
train_x_final <- train_x_scaled[, selected_features, drop = FALSE]
test_x_final <- test_x_scaled[, selected_features, drop = FALSE]

# ===================================================================
# 4. 定义模型训练函数
# ===================================================================

# 5折交叉验证函数 - 返回概率和fold编号
cv_predict <- function(X, y, model_type, n_folds = 5) {
  folds <- createFolds(as.factor(y), k = n_folds)
  cv_probs <- numeric(length(y))
  fold_ids <- numeric(length(y))
  
  for (i in 1:n_folds) {
    test_idx <- folds[[i]]
    train_idx <- -test_idx
    
    X_train_cv <- X[train_idx, , drop = FALSE]
    y_train_cv <- y[train_idx]
    X_test_cv <- X[test_idx, , drop = FALSE]
    
    # 记录fold编号
    fold_ids[test_idx] <- i
    
    # 训练模型
    if (model_type == "MLP") {
      df_train <- data.frame(y = as.factor(y_train_cv), X_train_cv)
      model <- nnet(y ~ ., data = df_train, size = 10, maxit = 500, trace = FALSE)
      df_test <- data.frame(X_test_cv)
      pred_matrix <- predict(model, df_test, type = "raw")
      
      if(is.matrix(pred_matrix)) {
        cv_probs[test_idx] <- pred_matrix[, ncol(pred_matrix)]
      } else {
        cv_probs[test_idx] <- as.numeric(pred_matrix)
      }
      
    } else if (model_type == "LASSO") {
      model <- cv.glmnet(X_train_cv, y_train_cv, family = "binomial", alpha = 1)
      cv_probs[test_idx] <- predict(model, X_test_cv, s = "lambda.min", type = "response")[,1]
      
    } else if (model_type == "Logistic") {
      df_train <- data.frame(y = y_train_cv, X_train_cv)
      model <- glm(y ~ ., data = df_train, family = binomial)
      df_test <- data.frame(X_test_cv)
      cv_probs[test_idx] <- predict(model, df_test, type = "response")
      
    } else if (model_type == "Random Forest") {
      model <- randomForest(X_train_cv, as.factor(y_train_cv), ntree = 500)
      cv_probs[test_idx] <- predict(model, X_test_cv, type = "prob")[, 2]
      
    } else if (model_type == "SVM") {
      model <- svm(X_train_cv, as.factor(y_train_cv), kernel = "radial", probability = TRUE)
      pred_probs <- attr(predict(model, X_test_cv, probability = TRUE), "probabilities")
      if("1" %in% colnames(pred_probs)) {
        cv_probs[test_idx] <- pred_probs[, "1"]
      } else {
        cv_probs[test_idx] <- pred_probs[, 2]
      }
    }
  }
  return(list(probs = cv_probs, folds = fold_ids))
}

# 测试集预测函数
test_predict <- function(X_train, y_train, X_test, model_type) {
  if (model_type == "MLP") {
    df_train <- data.frame(y = as.factor(y_train), X_train)
    model <- nnet(y ~ ., data = df_train, size = 10, maxit = 500, trace = FALSE)
    df_test <- data.frame(X_test)
    pred_matrix <- predict(model, df_test, type = "raw")
    
    if(is.matrix(pred_matrix)) {
      probs <- pred_matrix[, ncol(pred_matrix)]
    } else {
      probs <- as.numeric(pred_matrix)
    }
    
  } else if (model_type == "LASSO") {
    model <- cv.glmnet(X_train, y_train, family = "binomial", alpha = 1)
    probs <- predict(model, X_test, s = "lambda.min", type = "response")[,1]
    
  } else if (model_type == "Logistic") {
    df_train <- data.frame(y = y_train, X_train)
    model <- glm(y ~ ., data = df_train, family = binomial)
    df_test <- data.frame(X_test)
    probs <- predict(model, df_test, type = "response")
    
  } else if (model_type == "Random Forest") {
    model <- randomForest(X_train, as.factor(y_train), ntree = 500)
    probs <- predict(model, X_test, type = "prob")[, 2]
    
  } else if (model_type == "SVM") {
    model <- svm(X_train, as.factor(y_train), kernel = "radial", probability = TRUE)
    pred_probs <- attr(predict(model, X_test, probability = TRUE), "probabilities")
    if("1" %in% colnames(pred_probs)) {
      probs <- pred_probs[, "1"]
    } else {
      probs <- pred_probs[, 2]
    }
  }
  return(probs)
}

# ===================================================================
# 5. 训练所有模型并获取预测概率
# ===================================================================

cat("\n==================== 模型训练 ====================\n")

# 存储结果
cv_results <- list()
cv_folds <- list()
test_results <- list()

for (model in model_names) {
  cat(sprintf("训练 %s...\n", model))
  
  # 5折交叉验证
  cv_output <- cv_predict(train_x_final, train_y, model)
  cv_results[[model]] <- cv_output$probs
  cv_folds[[model]] <- cv_output$folds
  
  # 测试集预测
  test_probs <- test_predict(train_x_final, train_y, test_x_final, model)
  test_results[[model]] <- test_probs
}

# ===================================================================
# 6. 导出概率值
# ===================================================================

cat("\n==================== 导出概率值 ====================\n")

# 训练集5折交叉验证概率 (包含fold列)
cv_prob_df <- data.frame(
  id = train_data$id,
  label = train_data$label,
  fold = cv_folds[["MLP"]],  # 所有模型的fold一致
  MLP = cv_results[["MLP"]],
  LASSO = cv_results[["LASSO"]],
  Logistic = cv_results[["Logistic"]],
  RandomForest = cv_results[["Random Forest"]],
  SVM = cv_results[["SVM"]]
)

write.csv(cv_prob_df, "train_cv_probabilities.csv", row.names = FALSE)

# 测试集概率
test_prob_df <- data.frame(
  id = test_data$id,
  label = test_data$label,
  MLP = test_results[["MLP"]],
  LASSO = test_results[["LASSO"]],
  Logistic = test_results[["Logistic"]],
  RandomForest = test_results[["Random Forest"]],
  SVM = test_results[["SVM"]]
)

write.csv(test_prob_df, "test_probabilities.csv", row.names = FALSE)

# ===================================================================
# 7. 计算性能指标
# ===================================================================

cat("\n==================== 计算性能指标 ====================\n")

# 计算性能指标函数 (90%置信区间)
calculate_metrics <- function(true_labels, pred_probs) {
  # ROC曲线和AUC
  roc_obj <- roc(true_labels, pred_probs, quiet = TRUE, levels = c(0, 1), direction = "<")
  auc_val <- as.numeric(auc(roc_obj))
  auc_ci <- ci.auc(roc_obj, conf.level = 0.90)  # 90%置信区间
  
  # PR曲线和AUC
  pr_obj <- pr.curve(
    scores.class0 = pred_probs[true_labels == 1],
    scores.class1 = pred_probs[true_labels == 0],
    curve = TRUE
  )
  pr_auc <- pr_obj$auc.integral
  
  # 最优阈值
  coords_obj <- coords(roc_obj, "best", best.method = "youden")
  threshold <- coords_obj$threshold
  
  # 预测标签
  pred_labels <- ifelse(pred_probs > threshold, 1, 0)
  
  # 混淆矩阵
  cm <- confusionMatrix(
    as.factor(pred_labels), 
    as.factor(true_labels),
    positive = "1"
  )
  
  sensitivity <- cm$byClass["Sensitivity"]
  specificity <- cm$byClass["Specificity"]
  accuracy <- cm$overall["Accuracy"]
  
  return(list(
    AUC_ROC = auc_val,
    AUC_ROC_CI_lower = auc_ci[1],
    AUC_ROC_CI_upper = auc_ci[3],
    AUC_PR = pr_auc,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Accuracy = accuracy,
    roc_obj = roc_obj,
    pr_obj = pr_obj
  ))
}

# 计算所有模型的性能
metrics_list_cv <- list()
metrics_list_test <- list()

for (model in model_names) {
  cat(sprintf("计算 %s 性能指标...\n", model))
  metrics_list_cv[[model]] <- calculate_metrics(train_y, cv_results[[model]])
  metrics_list_test[[model]] <- calculate_metrics(test_y, test_results[[model]])
}

# 导出性能指标表格
metrics_df <- data.frame(
  Model = rep(model_names, 2),
  Dataset = c(rep("Train_CV", 5), rep("Test", 5)),
  AUC_ROC = c(
    sapply(metrics_list_cv, function(x) x$AUC_ROC),
    sapply(metrics_list_test, function(x) x$AUC_ROC)
  ),
  AUC_ROC_CI_lower = c(
    sapply(metrics_list_cv, function(x) x$AUC_ROC_CI_lower),
    sapply(metrics_list_test, function(x) x$AUC_ROC_CI_lower)
  ),
  AUC_ROC_CI_upper = c(
    sapply(metrics_list_cv, function(x) x$AUC_ROC_CI_upper),
    sapply(metrics_list_test, function(x) x$AUC_ROC_CI_upper)
  ),
  AUC_PR = c(
    sapply(metrics_list_cv, function(x) x$AUC_PR),
    sapply(metrics_list_test, function(x) x$AUC_PR)
  ),
  Sensitivity = c(
    sapply(metrics_list_cv, function(x) x$Sensitivity),
    sapply(metrics_list_test, function(x) x$Sensitivity)
  ),
  Specificity = c(
    sapply(metrics_list_cv, function(x) x$Specificity),
    sapply(metrics_list_test, function(x) x$Specificity)
  ),
  Accuracy = c(
    sapply(metrics_list_cv, function(x) x$Accuracy),
    sapply(metrics_list_test, function(x) x$Accuracy)
  )
)

write.csv(metrics_df, "model_performance_metrics.csv", row.names = FALSE)

# ===================================================================
# 8. 绘制ROC和PR曲线 - 4张图，每张图包含5个算法
# ===================================================================

cat("\n==================== 绘制ROC和PR曲线 ====================\n")

# 训练集ROC曲线 (5折交叉验证) - 特征选择前后对比
plot_train_roc_before <- ggplot() +
  theme_cowplot() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = line_width) +
  labs(x = "1 - Specificity", y = "Sensitivity", title = "Train CV ROC (Before Feature Selection)") +
  theme(
    text = element_text(family = font_family, size = font_size, face = "bold"),
    axis.text = element_text(size = axis_text_size, color = "black", face = "bold"),
    panel.border = element_blank(),
    axis.line.x.bottom = element_line(color = "black", linewidth = line_width),
    axis.line.y.left = element_line(color = "black", linewidth = line_width),
    axis.ticks = element_line(color = "black", linewidth = line_width),
    legend.position = c(0.65, 0.25),
    legend.text = element_text(size = 8),
    legend.title = element_blank(),
    legend.key.size = unit(0.4, "cm"),
    plot.title = element_text(hjust = 0.5, size = font_size)
  )

# 添加LASSO和RF的特征选择前曲线
auc_before_train_lasso <- as.numeric(auc(roc_before_train_lasso))
ci_before_train_lasso <- ci.auc(roc_before_train_lasso, conf.level = 0.90)
label_before_lasso <- sprintf("LASSO %.2f (%.2f-%.2f)", 
                              auc_before_train_lasso, 
                              ci_before_train_lasso[1], 
                              ci_before_train_lasso[3])

auc_before_train_rf <- as.numeric(auc(roc_before_train_rf))
ci_before_train_rf <- ci.auc(roc_before_train_rf, conf.level = 0.90)
label_before_rf <- sprintf("Random Forest %.2f (%.2f-%.2f)", 
                           auc_before_train_rf, 
                           ci_before_train_rf[1], 
                           ci_before_train_rf[3])

plot_train_roc_before <- plot_train_roc_before +
  geom_line(
    data = data.frame(x = 1 - roc_before_train_lasso$specificities, 
                      y = roc_before_train_lasso$sensitivities),
    aes(x = x, y = y, color = label_before_lasso),
    linewidth = line_width
  ) +
  geom_line(
    data = data.frame(x = 1 - roc_before_train_rf$specificities, 
                      y = roc_before_train_rf$sensitivities),
    aes(x = x, y = y, color = label_before_rf),
    linewidth = line_width
  ) +
  scale_color_manual(values = c(model_colors[2], model_colors[4]))

# 训练集ROC曲线 (5折交叉验证) - 特征选择后
plot_train_roc <- ggplot() +
  theme_cowplot() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = line_width) +
  labs(x = "1 - Specificity", y = "Sensitivity", title = "Train CV ROC (After Feature Selection)") +
  theme(
    text = element_text(family = font_family, size = font_size, face = "bold"),
    axis.text = element_text(size = axis_text_size, color = "black", face = "bold"),
    panel.border = element_blank(),
    axis.line.x.bottom = element_line(color = "black", linewidth = line_width),
    axis.line.y.left = element_line(color = "black", linewidth = line_width),
    axis.ticks = element_line(color = "black", linewidth = line_width),
    legend.position = c(0.65, 0.25),
    legend.text = element_text(size = 8),
    legend.title = element_blank(),
    legend.key.size = unit(0.4, "cm"),
    plot.title = element_text(hjust = 0.5, size = font_size)
  )

for (i in seq_along(model_names)) {
  model <- model_names[i]
  roc_data <- metrics_list_cv[[model]]$roc_obj
  auc_val <- metrics_list_cv[[model]]$AUC_ROC
  ci_lower <- metrics_list_cv[[model]]$AUC_ROC_CI_lower
  ci_upper <- metrics_list_cv[[model]]$AUC_ROC_CI_upper
  
  label_text <- sprintf("%s %.2f (%.2f-%.2f)", model, auc_val, ci_lower, ci_upper)
  
  plot_train_roc <- plot_train_roc +
    geom_line(
      data = data.frame(x = 1 - roc_data$specificities, y = roc_data$sensitivities),
      aes(x = x, y = y, color = label_text),
      linewidth = line_width
    )
}

plot_train_roc <- plot_train_roc +
  scale_color_manual(values = model_colors[1:5])

# 测试集ROC曲线 - 特征选择前后对比
plot_test_roc_before <- ggplot() +
  theme_cowplot() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = line_width) +
  labs(x = "1 - Specificity", y = "Sensitivity", title = "Test ROC (Before Feature Selection)") +
  theme(
    text = element_text(family = font_family, size = font_size, face = "bold"),
    axis.text = element_text(size = axis_text_size, color = "black", face = "bold"),
    panel.border = element_blank(),
    axis.line.x.bottom = element_line(color = "black", linewidth = line_width),
    axis.line.y.left = element_line(color = "black", linewidth = line_width),
    axis.ticks = element_line(color = "black", linewidth = line_width),
    legend.position = c(0.65, 0.25),
    legend.text = element_text(size = 8),
    legend.title = element_blank(),
    legend.key.size = unit(0.4, "cm"),
    plot.title = element_text(hjust = 0.5, size = font_size)
  )

auc_before_test_lasso <- as.numeric(auc(roc_before_test_lasso))
ci_before_test_lasso <- ci.auc(roc_before_test_lasso, conf.level = 0.90)
label_before_test_lasso <- sprintf("LASSO %.2f (%.2f-%.2f)", 
                                   auc_before_test_lasso, 
                                   ci_before_test_lasso[1], 
                                   ci_before_test_lasso[3])

auc_before_test_rf <- as.numeric(auc(roc_before_test_rf))
ci_before_test_rf <- ci.auc(roc_before_test_rf, conf.level = 0.90)
label_before_test_rf <- sprintf("Random Forest %.2f (%.2f-%.2f)", 
                                auc_before_test_rf, 
                                ci_before_test_rf[1], 
                                ci_before_test_rf[3])

plot_test_roc_before <- plot_test_roc_before +
  geom_line(
    data = data.frame(x = 1 - roc_before_test_lasso$specificities, 
                      y = roc_before_test_lasso$sensitivities),
    aes(x = x, y = y, color = label_before_test_lasso),
    linewidth = line_width
  ) +
  geom_line(
    data = data.frame(x = 1 - roc_before_test_rf$specificities, 
                      y = roc_before_test_rf$sensitivities),
    aes(x = x, y = y, color = label_before_test_rf),
    linewidth = line_width
  ) +
  scale_color_manual(values = c(model_colors[2], model_colors[4]))

# 测试集ROC曲线 - 特征选择后
plot_test_roc <- ggplot() +
  theme_cowplot() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = line_width) +
  labs(x = "1 - Specificity", y = "Sensitivity", title = "Test ROC (After Feature Selection)") +
  theme(
    text = element_text(family = font_family, size = font_size, face = "bold"),
    axis.text = element_text(size = axis_text_size, color = "black", face = "bold"),
    panel.border = element_blank(),
    axis.line.x.bottom = element_line(color = "black", linewidth = line_width),
    axis.line.y.left = element_line(color = "black", linewidth = line_width),
    axis.ticks = element_line(color = "black", linewidth = line_width),
    legend.position = c(0.65, 0.25),
    legend.text = element_text(size = 8),
    legend.title = element_blank(),
    legend.key.size = unit(0.4, "cm"),
    plot.title = element_text(hjust = 0.5, size = font_size)
  )

for (i in seq_along(model_names)) {
  model <- model_names[i]
  roc_data <- metrics_list_test[[model]]$roc_obj
  auc_val <- metrics_list_test[[model]]$AUC_ROC
  ci_lower <- metrics_list_test[[model]]$AUC_ROC_CI_lower
  ci_upper <- metrics_list_test[[model]]$AUC_ROC_CI_upper
  
  label_text <- sprintf("%s %.2f (%.2f-%.2f)", model, auc_val, ci_lower, ci_upper)
  
  plot_test_roc <- plot_test_roc +
    geom_line(
      data = data.frame(x = 1 - roc_data$specificities, y = roc_data$sensitivities),
      aes(x = x, y = y, color = label_text),
      linewidth = line_width
    )
}

plot_test_roc <- plot_test_roc +
  scale_color_manual(values = model_colors[1:5])

# 训练集PR曲线 - 特征选择前
plot_train_pr_before <- ggplot() +
  theme_cowplot() +
  labs(x = "Recall", y = "Precision", title = "Train CV PR (Before Feature Selection)") +
  theme(
    text = element_text(family = font_family, size = font_size, face = "bold"),
    axis.text = element_text(size = axis_text_size, color = "black", face = "bold"),
    panel.border = element_blank(),
    axis.line.x.bottom = element_line(color = "black", linewidth = line_width),
    axis.line.y.left = element_line(color = "black", linewidth = line_width),
    axis.ticks = element_line(color = "black", linewidth = line_width),
    legend.position = c(0.35, 0.25),
    legend.text = element_text(size = 8),
    legend.title = element_blank(),
    legend.key.size = unit(0.4, "cm"),
    plot.title = element_text(hjust = 0.5, size = font_size)
  )

auc_pr_before_train_lasso <- pr_before_train_lasso$auc.integral
label_pr_before_lasso <- sprintf("LASSO %.2f", auc_pr_before_train_lasso)

auc_pr_before_train_rf <- pr_before_train_rf$auc.integral
label_pr_before_rf <- sprintf("Random Forest %.2f", auc_pr_before_train_rf)

plot_train_pr_before <- plot_train_pr_before +
  geom_line(
    data = data.frame(x = pr_before_train_lasso$curve[, 1], 
                      y = pr_before_train_lasso$curve[, 2]),
    aes(x = x, y = y, color = label_pr_before_lasso),
    linewidth = line_width
  ) +
  geom_line(
    data = data.frame(x = pr_before_train_rf$curve[, 1], 
                      y = pr_before_train_rf$curve[, 2]),
    aes(x = x, y = y, color = label_pr_before_rf),
    linewidth = line_width
  ) +
  scale_color_manual(values = c(model_colors[2], model_colors[4]))

# 训练集PR曲线 - 特征选择后
plot_train_pr <- ggplot() +
  theme_cowplot() +
  labs(x = "Recall", y = "Precision", title = "Train CV PR (After Feature Selection)") +
  theme(
    text = element_text(family = font_family, size = font_size, face = "bold"),
    axis.text = element_text(size = axis_text_size, color = "black", face = "bold"),
    panel.border = element_blank(),
    axis.line.x.bottom = element_line(color = "black", linewidth = line_width),
    axis.line.y.left = element_line(color = "black", linewidth = line_width),
    axis.ticks = element_line(color = "black", linewidth = line_width),
    legend.position = c(0.35, 0.25),
    legend.text = element_text(size = 8),
    legend.title = element_blank(),
    legend.key.size = unit(0.4, "cm"),
    plot.title = element_text(hjust = 0.5, size = font_size)
  )

for (i in seq_along(model_names)) {
  model <- model_names[i]
  pr_data <- metrics_list_cv[[model]]$pr_obj
  pr_auc <- metrics_list_cv[[model]]$AUC_PR
  
  label_text <- sprintf("%s %.2f", model, pr_auc)
  
  plot_train_pr <- plot_train_pr +
    geom_line(
      data = data.frame(x = pr_data$curve[, 1], y = pr_data$curve[, 2]),
      aes(x = x, y = y, color = label_text),
      linewidth = line_width
    )
}

plot_train_pr <- plot_train_pr +
  scale_color_manual(values = model_colors[1:5])

# 测试集PR曲线 - 特征选择前
plot_test_pr_before <- ggplot() +
  theme_cowplot() +
  labs(x = "Recall", y = "Precision", title = "Test PR (Before Feature Selection)") +
  theme(
    text = element_text(family = font_family, size = font_size, face = "bold"),
    axis.text = element_text(size = axis_text_size, color = "black", face = "bold"),
    panel.border = element_blank(),
    axis.line.x.bottom = element_line(color = "black", linewidth = line_width),
    axis.line.y.left = element_line(color = "black", linewidth = line_width),
    axis.ticks = element_line(color = "black", linewidth = line_width),
    legend.position = c(0.35, 0.25),
    legend.text = element_text(size = 8),
    legend.title = element_blank(),
    legend.key.size = unit(0.4, "cm"),
    plot.title = element_text(hjust = 0.5, size = font_size)
  )

auc_pr_before_test_lasso <- pr_before_test_lasso$auc.integral
label_pr_before_test_lasso <- sprintf("LASSO %.2f", auc_pr_before_test_lasso)

auc_pr_before_test_rf <- pr_before_test_rf$auc.integral
label_pr_before_test_rf <- sprintf("Random Forest %.2f", auc_pr_before_test_rf)

plot_test_pr_before <- plot_test_pr_before +
  geom_line(
    data = data.frame(x = pr_before_test_lasso$curve[, 1], 
                      y = pr_before_test_lasso$curve[, 2]),
    aes(x = x, y = y, color = label_pr_before_test_lasso),
    linewidth = line_width
  ) +
  geom_line(
    data = data.frame(x = pr_before_test_rf$curve[, 1], 
                      y = pr_before_test_rf$curve[, 2]),
    aes(x = x, y = y, color = label_pr_before_test_rf),
    linewidth = line_width
  ) +
  scale_color_manual(values = c(model_colors[2], model_colors[4]))

# 测试集PR曲线 - 特征选择后
plot_test_pr <- ggplot() +
  theme_cowplot() +
  labs(x = "Recall", y = "Precision", title = "Test PR (After Feature Selection)") +
  theme(
    text = element_text(family = font_family, size = font_size, face = "bold"),
    axis.text = element_text(size = axis_text_size, color = "black", face = "bold"),
    panel.border = element_blank(),
    axis.line.x.bottom = element_line(color = "black", linewidth = line_width),
    axis.line.y.left = element_line(color = "black", linewidth = line_width),
    axis.ticks = element_line(color = "black", linewidth = line_width),
    legend.position = c(0.35, 0.25),
    legend.text = element_text(size = 8),
    legend.title = element_blank(),
    legend.key.size = unit(0.4, "cm"),
    plot.title = element_text(hjust = 0.5, size = font_size)
  )

for (i in seq_along(model_names)) {
  model <- model_names[i]
  pr_data <- metrics_list_test[[model]]$pr_obj
  pr_auc <- metrics_list_test[[model]]$AUC_PR
  
  label_text <- sprintf("%s %.2f", model, pr_auc)
  
  plot_test_pr <- plot_test_pr +
    geom_line(
      data = data.frame(x = pr_data$curve[, 1], y = pr_data$curve[, 2]),
      aes(x = x, y = y, color = label_text),
      linewidth = line_width
    )
}

plot_test_pr <- plot_test_pr +
  scale_color_manual(values = model_colors[1:5])

# 保存8张独立图片
ggsave("train_cv_roc_before.jpg", plot_train_roc_before, width = plot_width, height = plot_height, units = "cm", dpi = 300)
ggsave("train_cv_roc_after.jpg", plot_train_roc, width = plot_width, height = plot_height, units = "cm", dpi = 300)
ggsave("test_roc_before.jpg", plot_test_roc_before, width = plot_width, height = plot_height, units = "cm", dpi = 300)
ggsave("test_roc_after.jpg", plot_test_roc, width = plot_width, height = plot_height, units = "cm", dpi = 300)
ggsave("train_cv_pr_before.jpg", plot_train_pr_before, width = plot_width, height = plot_height, units = "cm", dpi = 300)
ggsave("train_cv_pr_after.jpg", plot_train_pr, width = plot_width, height = plot_height, units = "cm", dpi = 300)
ggsave("test_pr_before.jpg", plot_test_pr_before, width = plot_width, height = plot_height, units = "cm", dpi = 300)
ggsave("test_pr_after.jpg", plot_test_pr, width = plot_width, height = plot_height, units = "cm", dpi = 300)

# 保存组合图 (特征选择后的4张)
combined_plot_after <- grid.arrange(
  plot_train_roc, plot_test_roc,
  plot_train_pr, plot_test_pr,
  ncol = 2
)

ggsave("combined_roc_pr_after.jpg", combined_plot_after, width = 15, height = 14, units = "cm", dpi = 300)

# 保存组合图 (特征选择前的4张 - 仅LASSO和RF)
combined_plot_before <- grid.arrange(
  plot_train_roc_before, plot_test_roc_before,
  plot_train_pr_before, plot_test_pr_before,
  ncol = 2
)

ggsave("combined_roc_pr_before.jpg", combined_plot_before, width = 15, height = 14, units = "cm", dpi = 300)

# ===================================================================
# 完成
# ===================================================================

cat("\n==================== 分析完成 ====================\n")
cat("生成文件:\n")
cat("1. selected_features.csv - 选择的特征\n")
cat("2. train_cv_probabilities.csv - 训练集5折交叉验证概率(含fold列)\n")
cat("3. test_probabilities.csv - 测试集概率\n")
cat("4. model_performance_metrics.csv - 模型性能指标(90%置信区间)\n")
cat("5-8. 特征选择前的ROC和PR曲线 (仅LASSO和RF)\n")
cat("9-12. 特征选择后的ROC和PR曲线 (5个算法)\n")
cat("13. combined_roc_pr_before.jpg - 特征选择前组合图\n")
cat("14. combined_roc_pr_after.jpg - 特征选择后组合图\n")
cat("==================================================\n")