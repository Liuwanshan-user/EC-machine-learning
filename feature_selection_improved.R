# ===================================================================
# 蛋白组学数据特征选择和多算法模型评估 - 改进版
# 特征选择前后完整对比，95%置信区间，每折概率导出
# ===================================================================

# 检查并安装所需的包
packages_needed <- c(
  "ggplot2", "cowplot", "dplyr", "pROC", "glmnet",
  "randomForest", "e1071", "neuralnet", "caret", "PRROC",
  "scales", "readr", "gridExtra", "boot"
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
library(neuralnet)
library(caret)
library(PRROC)
library(scales)
library(readr)
library(gridExtra)
library(boot)

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
# 2. 定义模型训练函数
# ===================================================================

# 定义配色方案
model_colors <- c(
  "#3AB5B3",  # 青色 - NN
  "#7B6C9F",  # 紫色 - LASSO
  "#A188BD",  # 浅紫 - Logistic
  "#BBC5DE",  # 浅蓝 - Random Forest
  "#E7777F"   # 粉色 - SVM
)

model_names <- c("NN", "LASSO", "Logistic", "Random Forest", "SVM")

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
    if (model_type == "NN") {
      # 两层神经网络，每层50个神经元
      # 使用简单的数字列名避免解析问题
      X_train_clean <- as.data.frame(X_train_cv)
      X_test_clean <- as.data.frame(X_test_cv)

      # 使用 X1, X2, X3... 命名
      n_features <- ncol(X_train_clean)
      simple_names <- paste0("X", 1:n_features)
      colnames(X_train_clean) <- simple_names
      colnames(X_test_clean) <- simple_names

      df_train <- data.frame(y = y_train_cv, X_train_clean)
      df_test <- X_test_clean

      # 创建公式
      formula_str <- paste("y ~", paste(simple_names, collapse = " + "))
      formula_obj <- as.formula(formula_str)

      # 训练两层神经网络
      model <- neuralnet(formula_obj, data = df_train,
                        hidden = c(50, 50),
                        linear.output = FALSE,
                        threshold = 0.01,
                        stepmax = 1e6,
                        rep = 1,
                        err.fct = "ce",
                        act.fct = "logistic",
                        likelihood = TRUE)

      # 预测
      pred_result <- compute(model, df_test)
      cv_probs[test_idx] <- pred_result$net.result[, 1]

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
  if (model_type == "NN") {
    # 两层神经网络，每层50个神经元
    # 使用简单的数字列名避免解析问题
    X_train_clean <- as.data.frame(X_train)
    X_test_clean <- as.data.frame(X_test)

    # 使用 X1, X2, X3... 命名
    n_features <- ncol(X_train_clean)
    simple_names <- paste0("X", 1:n_features)
    colnames(X_train_clean) <- simple_names
    colnames(X_test_clean) <- simple_names

    df_train <- data.frame(y = y_train, X_train_clean)
    df_test <- X_test_clean

    # 创建公式
    formula_str <- paste("y ~", paste(simple_names, collapse = " + "))
    formula_obj <- as.formula(formula_str)

    # 训练两层神经网络
    model <- neuralnet(formula_obj, data = df_train,
                      hidden = c(50, 50),
                      linear.output = FALSE,
                      threshold = 0.01,
                      stepmax = 1e6,
                      rep = 1,
                      err.fct = "ce",
                      act.fct = "logistic",
                      likelihood = TRUE)

    # 预测
    pred_result <- compute(model, df_test)
    probs <- pred_result$net.result[, 1]

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

# Bootstrap函数计算PR AUC的置信区间
bootstrap_pr_auc <- function(y_true, y_pred, n_boot = 1000, conf_level = 0.95) {
  auc_boot <- numeric(n_boot)
  n <- length(y_true)

  for (i in 1:n_boot) {
    idx <- sample(1:n, n, replace = TRUE)
    y_boot <- y_true[idx]
    pred_boot <- y_pred[idx]

    # 确保至少有两个类别
    if (length(unique(y_boot)) < 2) {
      next
    }

    pr_obj <- tryCatch({
      pr.curve(
        scores.class0 = pred_boot[y_boot == 1],
        scores.class1 = pred_boot[y_boot == 0],
        curve = FALSE
      )
    }, error = function(e) NULL)

    if (!is.null(pr_obj)) {
      auc_boot[i] <- pr_obj$auc.integral
    }
  }

  # 移除NA值
  auc_boot <- auc_boot[auc_boot > 0]

  alpha <- 1 - conf_level
  ci_lower <- quantile(auc_boot, alpha/2, na.rm = TRUE)
  ci_upper <- quantile(auc_boot, 1 - alpha/2, na.rm = TRUE)

  return(c(lower = ci_lower, upper = ci_upper))
}

# 计算性能指标函数 (95%置信区间)
calculate_metrics <- function(true_labels, pred_probs) {
  # ROC曲线和AUC
  roc_obj <- roc(true_labels, pred_probs, quiet = TRUE, levels = c(0, 1), direction = "<")
  auc_val <- as.numeric(auc(roc_obj))
  auc_ci <- ci.auc(roc_obj, conf.level = 0.95)  # 95%置信区间

  # PR曲线和AUC
  pr_obj <- pr.curve(
    scores.class0 = pred_probs[true_labels == 1],
    scores.class1 = pred_probs[true_labels == 0],
    curve = TRUE
  )
  pr_auc <- pr_obj$auc.integral

  # PR AUC的95%置信区间 (使用bootstrap)
  pr_ci <- bootstrap_pr_auc(true_labels, pred_probs, n_boot = 1000, conf_level = 0.95)

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
    AUC_PR_CI_lower = pr_ci[1],
    AUC_PR_CI_upper = pr_ci[2],
    Sensitivity = sensitivity,
    Specificity = specificity,
    Accuracy = accuracy,
    roc_obj = roc_obj,
    pr_obj = pr_obj
  ))
}

# ===================================================================
# 3. 特征选择前的完整分析
# ===================================================================

cat("\n==================== 特征选择前 - 所有算法 ====================\n")

# 存储结果
before_cv_results <- list()
before_cv_folds <- list()
before_test_results <- list()

for (model in model_names) {
  cat(sprintf("特征选择前 - 训练 %s...\n", model))

  # 5折交叉验证
  cv_output <- cv_predict(train_x_scaled, train_y, model)
  before_cv_results[[model]] <- cv_output$probs
  before_cv_folds[[model]] <- cv_output$folds

  # 测试集预测
  test_probs <- test_predict(train_x_scaled, train_y, test_x_scaled, model)
  before_test_results[[model]] <- test_probs
}

# 导出特征选择前的概率
cat("\n导出特征选择前的概率...\n")

# 训练集5折交叉验证概率
before_cv_prob_df <- data.frame(
  id = train_data$id,
  label = train_data$label,
  fold = before_cv_folds[["NN"]],
  NN = before_cv_results[["NN"]],
  LASSO = before_cv_results[["LASSO"]],
  Logistic = before_cv_results[["Logistic"]],
  RandomForest = before_cv_results[["Random Forest"]],
  SVM = before_cv_results[["SVM"]]
)

write.csv(before_cv_prob_df, "before_fs_train_cv_probabilities.csv", row.names = FALSE)

# 测试集概率
before_test_prob_df <- data.frame(
  id = test_data$id,
  label = test_data$label,
  NN = before_test_results[["NN"]],
  LASSO = before_test_results[["LASSO"]],
  Logistic = before_test_results[["Logistic"]],
  RandomForest = before_test_results[["Random Forest"]],
  SVM = before_test_results[["SVM"]]
)

write.csv(before_test_prob_df, "before_fs_test_probabilities.csv", row.names = FALSE)

# 计算特征选择前的性能指标
cat("\n计算特征选择前的性能指标...\n")

before_metrics_list_cv <- list()
before_metrics_list_test <- list()

for (model in model_names) {
  cat(sprintf("特征选择前 - 计算 %s 性能指标...\n", model))
  before_metrics_list_cv[[model]] <- calculate_metrics(train_y, before_cv_results[[model]])
  before_metrics_list_test[[model]] <- calculate_metrics(test_y, before_test_results[[model]])
}

# 导出特征选择前的性能指标
before_metrics_df <- data.frame(
  Model = rep(model_names, 2),
  Dataset = c(rep("Train_CV", 5), rep("Test", 5)),
  AUC_ROC = c(
    sapply(before_metrics_list_cv, function(x) x$AUC_ROC),
    sapply(before_metrics_list_test, function(x) x$AUC_ROC)
  ),
  AUC_ROC_CI_lower = c(
    sapply(before_metrics_list_cv, function(x) x$AUC_ROC_CI_lower),
    sapply(before_metrics_list_test, function(x) x$AUC_ROC_CI_lower)
  ),
  AUC_ROC_CI_upper = c(
    sapply(before_metrics_list_cv, function(x) x$AUC_ROC_CI_upper),
    sapply(before_metrics_list_test, function(x) x$AUC_ROC_CI_upper)
  ),
  AUC_PR = c(
    sapply(before_metrics_list_cv, function(x) x$AUC_PR),
    sapply(before_metrics_list_test, function(x) x$AUC_PR)
  ),
  AUC_PR_CI_lower = c(
    sapply(before_metrics_list_cv, function(x) x$AUC_PR_CI_lower),
    sapply(before_metrics_list_test, function(x) x$AUC_PR_CI_lower)
  ),
  AUC_PR_CI_upper = c(
    sapply(before_metrics_list_cv, function(x) x$AUC_PR_CI_upper),
    sapply(before_metrics_list_test, function(x) x$AUC_PR_CI_upper)
  ),
  Sensitivity = c(
    sapply(before_metrics_list_cv, function(x) x$Sensitivity),
    sapply(before_metrics_list_test, function(x) x$Sensitivity)
  ),
  Specificity = c(
    sapply(before_metrics_list_cv, function(x) x$Specificity),
    sapply(before_metrics_list_test, function(x) x$Specificity)
  ),
  Accuracy = c(
    sapply(before_metrics_list_cv, function(x) x$Accuracy),
    sapply(before_metrics_list_test, function(x) x$Accuracy)
  )
)

write.csv(before_metrics_df, "before_fs_performance_metrics.csv", row.names = FALSE)

# ===================================================================
# 4. 特征选择: LASSO + 随机森林交集
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
# 5. 特征选择后的完整分析
# ===================================================================

cat("\n==================== 特征选择后 - 所有算法 ====================\n")

# 存储结果
after_cv_results <- list()
after_cv_folds <- list()
after_test_results <- list()

for (model in model_names) {
  cat(sprintf("特征选择后 - 训练 %s...\n", model))

  # 5折交叉验证
  cv_output <- cv_predict(train_x_final, train_y, model)
  after_cv_results[[model]] <- cv_output$probs
  after_cv_folds[[model]] <- cv_output$folds

  # 测试集预测
  test_probs <- test_predict(train_x_final, train_y, test_x_final, model)
  after_test_results[[model]] <- test_probs
}

# 导出特征选择后的概率
cat("\n导出特征选择后的概率...\n")

# 训练集5折交叉验证概率
after_cv_prob_df <- data.frame(
  id = train_data$id,
  label = train_data$label,
  fold = after_cv_folds[["NN"]],
  NN = after_cv_results[["NN"]],
  LASSO = after_cv_results[["LASSO"]],
  Logistic = after_cv_results[["Logistic"]],
  RandomForest = after_cv_results[["Random Forest"]],
  SVM = after_cv_results[["SVM"]]
)

write.csv(after_cv_prob_df, "after_fs_train_cv_probabilities.csv", row.names = FALSE)

# 测试集概率
after_test_prob_df <- data.frame(
  id = test_data$id,
  label = test_data$label,
  NN = after_test_results[["NN"]],
  LASSO = after_test_results[["LASSO"]],
  Logistic = after_test_results[["Logistic"]],
  RandomForest = after_test_results[["Random Forest"]],
  SVM = after_test_results[["SVM"]]
)

write.csv(after_test_prob_df, "after_fs_test_probabilities.csv", row.names = FALSE)

# 计算特征选择后的性能指标
cat("\n计算特征选择后的性能指标...\n")

after_metrics_list_cv <- list()
after_metrics_list_test <- list()

for (model in model_names) {
  cat(sprintf("特征选择后 - 计算 %s 性能指标...\n", model))
  after_metrics_list_cv[[model]] <- calculate_metrics(train_y, after_cv_results[[model]])
  after_metrics_list_test[[model]] <- calculate_metrics(test_y, after_test_results[[model]])
}

# 导出特征选择后的性能指标
after_metrics_df <- data.frame(
  Model = rep(model_names, 2),
  Dataset = c(rep("Train_CV", 5), rep("Test", 5)),
  AUC_ROC = c(
    sapply(after_metrics_list_cv, function(x) x$AUC_ROC),
    sapply(after_metrics_list_test, function(x) x$AUC_ROC)
  ),
  AUC_ROC_CI_lower = c(
    sapply(after_metrics_list_cv, function(x) x$AUC_ROC_CI_lower),
    sapply(after_metrics_list_test, function(x) x$AUC_ROC_CI_lower)
  ),
  AUC_ROC_CI_upper = c(
    sapply(after_metrics_list_cv, function(x) x$AUC_ROC_CI_upper),
    sapply(after_metrics_list_test, function(x) x$AUC_ROC_CI_upper)
  ),
  AUC_PR = c(
    sapply(after_metrics_list_cv, function(x) x$AUC_PR),
    sapply(after_metrics_list_test, function(x) x$AUC_PR)
  ),
  AUC_PR_CI_lower = c(
    sapply(after_metrics_list_cv, function(x) x$AUC_PR_CI_lower),
    sapply(after_metrics_list_test, function(x) x$AUC_PR_CI_lower)
  ),
  AUC_PR_CI_upper = c(
    sapply(after_metrics_list_cv, function(x) x$AUC_PR_CI_upper),
    sapply(after_metrics_list_test, function(x) x$AUC_PR_CI_upper)
  ),
  Sensitivity = c(
    sapply(after_metrics_list_cv, function(x) x$Sensitivity),
    sapply(after_metrics_list_test, function(x) x$Sensitivity)
  ),
  Specificity = c(
    sapply(after_metrics_list_cv, function(x) x$Specificity),
    sapply(after_metrics_list_test, function(x) x$Specificity)
  ),
  Accuracy = c(
    sapply(after_metrics_list_cv, function(x) x$Accuracy),
    sapply(after_metrics_list_test, function(x) x$Accuracy)
  )
)

write.csv(after_metrics_df, "after_fs_performance_metrics.csv", row.names = FALSE)

# ===================================================================
# 6. 绘制ROC和PR曲线
# ===================================================================

cat("\n==================== 绘制ROC和PR曲线 ====================\n")

# 函数：创建ROC图
create_roc_plot <- function(metrics_list) {
  plot_obj <- ggplot() +
    theme_cowplot() +
    labs(x = "1 - Specificity", y = "Sensitivity") +
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
      legend.key.size = unit(0.4, "cm")
    )

  for (i in seq_along(model_names)) {
    model <- model_names[i]
    roc_data <- metrics_list[[model]]$roc_obj
    auc_val <- metrics_list[[model]]$AUC_ROC

    label_text <- sprintf("%s %.2f", model, auc_val)

    plot_obj <- plot_obj +
      geom_line(
        data = data.frame(x = 1 - roc_data$specificities, y = roc_data$sensitivities),
        aes(x = x, y = y, color = label_text),
        linewidth = line_width
      )
  }

  plot_obj <- plot_obj + scale_color_manual(values = model_colors[1:5])
  return(plot_obj)
}

# 函数：创建PR图
create_pr_plot <- function(metrics_list) {
  plot_obj <- ggplot() +
    theme_cowplot() +
    labs(x = "Recall", y = "Precision") +
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
      legend.key.size = unit(0.4, "cm")
    )

  for (i in seq_along(model_names)) {
    model <- model_names[i]
    pr_data <- metrics_list[[model]]$pr_obj
    pr_auc <- metrics_list[[model]]$AUC_PR

    label_text <- sprintf("%s %.2f", model, pr_auc)

    plot_obj <- plot_obj +
      geom_line(
        data = data.frame(x = pr_data$curve[, 1], y = pr_data$curve[, 2]),
        aes(x = x, y = y, color = label_text),
        linewidth = line_width
      )
  }

  plot_obj <- plot_obj + scale_color_manual(values = model_colors[1:5])
  return(plot_obj)
}

# 特征选择前的图
plot_before_train_roc <- create_roc_plot(before_metrics_list_cv)
plot_before_test_roc <- create_roc_plot(before_metrics_list_test)
plot_before_train_pr <- create_pr_plot(before_metrics_list_cv)
plot_before_test_pr <- create_pr_plot(before_metrics_list_test)

# 特征选择后的图
plot_after_train_roc <- create_roc_plot(after_metrics_list_cv)
plot_after_test_roc <- create_roc_plot(after_metrics_list_test)
plot_after_train_pr <- create_pr_plot(after_metrics_list_cv)
plot_after_test_pr <- create_pr_plot(after_metrics_list_test)

# 保存特征选择前的图
ggsave("before_fs_train_cv_roc.jpg", plot_before_train_roc, width = plot_width, height = plot_height, units = "cm", dpi = 300)
ggsave("before_fs_test_roc.jpg", plot_before_test_roc, width = plot_width, height = plot_height, units = "cm", dpi = 300)
ggsave("before_fs_train_cv_pr.jpg", plot_before_train_pr, width = plot_width, height = plot_height, units = "cm", dpi = 300)
ggsave("before_fs_test_pr.jpg", plot_before_test_pr, width = plot_width, height = plot_height, units = "cm", dpi = 300)

# 保存特征选择后的图
ggsave("after_fs_train_cv_roc.jpg", plot_after_train_roc, width = plot_width, height = plot_height, units = "cm", dpi = 300)
ggsave("after_fs_test_roc.jpg", plot_after_test_roc, width = plot_width, height = plot_height, units = "cm", dpi = 300)
ggsave("after_fs_train_cv_pr.jpg", plot_after_train_pr, width = plot_width, height = plot_height, units = "cm", dpi = 300)
ggsave("after_fs_test_pr.jpg", plot_after_test_pr, width = plot_width, height = plot_height, units = "cm", dpi = 300)

# ===================================================================
# 7. 完成
# ===================================================================

cat("\n==================== 分析完成 ====================\n")
cat("\n特征选择前生成的文件:\n")
cat("1. before_fs_train_cv_probabilities.csv - 训练集5折CV概率(含fold列)\n")
cat("2. before_fs_test_probabilities.csv - 测试集概率\n")
cat("3. before_fs_performance_metrics.csv - 性能指标(95%CI)\n")
cat("4. before_fs_train_cv_roc.jpg - 训练集ROC曲线\n")
cat("5. before_fs_test_roc.jpg - 测试集ROC曲线\n")
cat("6. before_fs_train_cv_pr.jpg - 训练集PR曲线\n")
cat("7. before_fs_test_pr.jpg - 测试集PR曲线\n")
cat("\n特征选择:\n")
cat("8. selected_features.csv - 选择的特征列表\n")
cat("\n特征选择后生成的文件:\n")
cat("9. after_fs_train_cv_probabilities.csv - 训练集5折CV概率(含fold列)\n")
cat("10. after_fs_test_probabilities.csv - 测试集概率\n")
cat("11. after_fs_performance_metrics.csv - 性能指标(95%CI)\n")
cat("12. after_fs_train_cv_roc.jpg - 训练集ROC曲线\n")
cat("13. after_fs_test_roc.jpg - 测试集ROC曲线\n")
cat("14. after_fs_train_cv_pr.jpg - 训练集PR曲线\n")
cat("15. after_fs_test_pr.jpg - 测试集PR曲线\n")
cat("==================================================\n")
