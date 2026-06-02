rm(list = ls())
setwd('/home/data/project/54.LZPJ-3001-6/')
file.dir = '05_ML'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)

library(rsample)
library(tidyverse)
library(tidymodels)
library(MLmetrics)
library(PRROC)
library(pROC)
library(caret)
library(bonsai)

dat <- read.csv('../00_rawdata/dat.GSE65204.csv',row.names = 1)
group <- read.csv('../00_rawdata/group.GSE65204.csv')
candigene <- read.csv('../02_DEPCDs/depcds.csv')
hub.dat <- dat[candigene$symbol,]
hub.dat <- t(hub.dat) %>% as.data.frame()
hub.dat$group <- group$group[match(group$sample,rownames(hub.dat))]
hub.dat$group <- factor(hub.dat$group , levels = c("Control", "Asthma"))
data <- hub.dat
data_test <- read.csv('../00_rawdata/dat.GSE19187.csv',row.names = 1)
group_test <- read.csv('../00_rawdata/group.GSE19187.csv')
allgene <- data.frame(rownames(data_test))
data_test <- t(data_test[candigene$symbol,]) %>% as.data.frame()
data_test$group <- group_test$group[match(group_test$sample,rownames(data_test))]
data_test$group <- factor(data_test$group,levels = c('Control','Asthma'))
data_train <- data
# 交叉验证
set.seed(666)
data_folds <- vfold_cv(data_train,v = 5, strata = group)
library(corrplot)
M <- cor(hub.dat[,-ncol(hub.dat)],method = 'spearman')
corrplot(M, type = "lower",tl.col = "black",tl.cex = 0.8,
         order = "hclust")

library(themis)

data_rec <- recipe(group ~ ., data = data_train) %>%
  step_corr(all_predictors(), threshold = 0.8,method = 'spearman')
rec_prepped <- prep(data_rec, training = data_train)
tidy(rec_prepped, number = 1)  # 1 是指 step_corr 是第1步；改为对应步骤编号
rec_prepped <- prep(data_rec, training = data_train)

original_vars <- colnames(dplyr::select(data_train, -group))

processed_data <- juice(rec_prepped)  # 获取处理后的数据
retained_vars <- colnames(dplyr::select(processed_data, -group))

removed_vars <- setdiff(original_vars, retained_vars)

# 打印
cat("保留的变量数量:", length(retained_vars), "\n")
cat("被去掉的变量数量:", length(removed_vars), "\n")
removed_vars
lr_spec <- logistic_reg() %>% set_engine("glm") %>% set_mode("classification")

rf_spec <- rand_forest(mtry = tune(), trees = tune()) %>%
  set_engine("ranger",importance = "permutation") %>%
  set_mode("classification")


enet_spec <- logistic_reg(penalty = tune(), mixture = tune()) %>%
  set_engine("glmnet") %>%
  set_mode("classification")

library(discrim)      
nb_spec <- naive_Bayes(smoothness = tune()) %>%
  set_engine("klaR") %>%
  set_mode("classification")
xgb_spec <- boost_tree(trees = tune(), tree_depth = tune(), learn_rate = tune()) %>%
  set_engine("xgboost") %>%
  set_mode("classification")

adaboost_spec <- boost_tree(
  trees = tune()
) %>%
  set_engine("C5.0") %>%
  set_mode("classification")


# 创建 workflow
models <- list(
  LR = workflow() %>% add_model(lr_spec) %>% add_recipe(data_rec),
  # SVM = workflow() %>% add_model(svm_spec) %>% add_recipe(data_rec),
  RF = workflow() %>% add_model(rf_spec) %>% add_recipe(data_rec),
  ENet = workflow() %>% add_model(enet_spec) %>% add_recipe(data_rec),
  XGB = workflow() %>% add_model(xgb_spec) %>% add_recipe(data_rec),
  # MLP = workflow() %>% add_model(mlp_spec) %>% add_recipe(data_rec),
  naive_Bayes = workflow() %>% add_model(nb_spec) %>% add_recipe(data_rec),
  # lightGBM = workflow() %>% add_model(lightgbm_spec) %>% add_recipe(data_rec),
  Adaboost =  workflow() %>% add_model(adaboost_spec) %>% add_recipe(data_rec)
)

rf_param <- parameters(
  mtry(range = c(1,2)),  # 视具体变量个数定
  trees(range = c(1700,2000))
)
nb_param <- parameters(nb_spec)
nb_param <- parameters(smoothness(range = c(0.1, 2)))

print(nb_param)
# 
ctrl <- control_bayes(verbose = TRUE,seed = 123456)
results <- list()
# name <- 'RF'
for (name in names(models)) {
  if (name == "LR") {
    set.seed(987654321)
    results[[name]] <- models[[name]] %>%
      fit_resamples(resamples = data_folds)
  } else if (name == "RF") {
    set.seed(987654321)
    results[[name]] <- models[[name]] %>%
      tune_bayes(resamples = data_folds, iter = 25, param_info = rf_param, control = ctrl)
  } else if (name == "naive_Bayes") {
    set.seed(987654321)
    results[[name]] <- models[[name]] %>%
      tune_bayes(resamples = data_folds, iter = 25, param_info = nb_param, control = ctrl)
  } else if (name == "lightGBM") {
    set.seed(987654321)
    results[[name]] <- models[[name]] %>%
      tune_bayes(resamples = data_folds, iter = 25, param_info = lightgbm_params, control = ctrl)
  } else {
    set.seed(123456)
    results[[name]] <- models[[name]] %>%
      tune_bayes(resamples = data_folds, iter = 25, control = ctrl)
  }
}
saveRDS(results,file = 'bayes.rds')
# 选择最优参数------
best_params <- lapply(results, function(res) {
  if ("tune_results" %in% class(res)) {
    select_best(res, metric = "roc_auc")
  } else {
    NULL
  }
})
show_best(results[["RF"]], metric = "roc_auc", n = 10)

library(dplyr)

best_params_df <- bind_rows(
  LR = best_params$LR,
  RF = best_params$RF,
  ENet = best_params$ENet,
  XGB = best_params$XGB,
  naive_Bayes = best_params$naive_Bayes,
  Adaboost = best_params$Adaboost,
  .id = "Model"
)

print(best_params_df)
write.csv(best_params_df,file = 'best_params.csv',row.names = F,quote = F)

final_models <- list()
for (name in names(models)) {
  if (name == "LR") {
    # 获取 LR 模型对象
    set.seed(123)
    lr_model <- lr_spec 
    # 使用数据训练 LR 模型
    final_models[[name]] <- lr_model %>% fit(group ~ ., data = data_train)
  } else {
    final_models[[name]] <- models[[name]] %>%
      finalize_workflow(best_params[[name]]) %>%
      # last_fit(data_split)
      fit(data_train)
  }
}

saveRDS(final_models,file = 'final.models.rds')

##训练集--------
model_results <- data.frame(
  Model = character(),
  AUC = numeric(),
  Accuracy = numeric(),
  Sensitivity = numeric(),
  Precision = numeric(),
  F1_Score = numeric(),
  PRAUC = numeric(),
  stringsAsFactors = FALSE
)
confusion <- list()

# 循环遍历每个模型并计算预测结果和评估指标
for (model_name in names(final_models)) {
  # 提取模型的预测结果
  pred_prob <- predict(final_models[[model_name]], data_train, type = "prob")[,2]
  predictions <- predict(final_models[[model_name]], data_train)
  predictions_class <- predictions$.pred_class
  aucvalue <- auc(roc(data_train$group, pred_prob$.pred_Asthma))
  # 计算混淆矩阵
  conf_matrix <- ConfusionMatrix(as.factor(predictions_class), as.factor(data_train$group))
  confusion[[model_name]] <- conf_matrix
  # 计算其他指标
  accuracy <- Accuracy(predictions_class, data_train$group)
  sensitivity <- Sensitivity(predictions_class, data_train$group)
  precision <- Precision(predictions_class, data_train$group)
  f1 <- F1_Score(predictions_class, data_train$group)
  # 计算 PR 曲线并捕获可能的错误
  prauc_value <- tryCatch({
    pr <- pr.curve(scores.class0 = predictions_class[data_train$group == 'Asthma'],
                   scores.class1 = predictions_class[data_train$group == 'Control'],
                   curve = TRUE)
    pr$auc.integral
  }, error = function(e) {
    # 如果计算失败，返回 0
    0
  })
  # 将结果存储到数据框中
  model_results <- rbind(model_results, data.frame(
    Model = model_name,
    AUC = aucvalue,
    Accuracy = accuracy,
    Sensitivity = sensitivity,
    Precision = precision,
    F1_Score = f1,
    PRAUC = prauc_value
  ))
  data_train[[paste0(model_name, "_prob")]] <- pred_prob$.pred_Asthma
}
write.csv(model_results, "model_evaluation_results.csv", row.names = FALSE,quote = F)

library(ggplot2)
library(reshape2)

model_results[, c(2:7)] <- lapply(model_results[, c(2:7)], function(x) as.numeric(unlist(x)))
results_long <-  melt(model_results, id.vars = "Model")
colnames(results_long)
str(model_results)
p <- ggplot(results_long, aes(x = variable, y = value, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  labs(x = "", y = "", title = "Train Cohort: Comparison of Model Performance") +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  scale_fill_manual(values = c("LR" = "#F781BF",
                               "RF" = "#b5aa0f",
                               "ENet" = "#984EA3",
                               'XGB' = '#CD534CFF',
                               'naive_Bayes' = '#3B4992FF',
                               'Adaboost' = '#8DD3C7'))
p
ggsave(filename = '01.MLPerformance_train.pdf',p,width = 7,height = 4)
ggsave(filename = '01.MLPerformance_train.png',p,width = 7,height = 4,units = 'in',dpi = 300)

library(ggplot2)
library(pROC)

colors <- c("LR" = "#F781BF",
            "RF" = "#b5aa0f",
            "ENet" = "#984EA3",
            'XGB' = '#CD534CFF',
            'naive_Bayes' = '#3B4992FF',
            'Adaboost' = '#8DD3C7')

pdf(file = '02.ROC_train.pdf', width = 5, height = 5.5, onefile = FALSE)
plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1), xlab = "1 - Specificity", ylab = "Sensitivity", 
     main = "Train Cohort: ROC Curves for Models")
legend_text <- c()
for (model_name in names(final_models)) {
  pred_prob <- data_train[[paste0(model_name, "_prob")]]
  roc_obj <- roc(data_train$group, pred_prob)
  auc_value <- auc(roc_obj)
  
  lines(1 - roc_obj$specificities, roc_obj$sensitivities, col = colors[model_name], lwd = 2)
  
  legend_text <- c(legend_text, paste0(model_name, " (AUC = ", round(auc_value, 3), ")"))
}

legend("bottomright", legend = legend_text, col = colors[names(final_models)], lwd = 2, bty = "n")
dev.off()

png(file = '02.ROC_train.png', width = 5, height = 5.5, units = 'in', bg = 'white', res = 300)
plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1), xlab = "1 - Specificity", ylab = "Sensitivity",
     main = "Train Cohort: ROC Curves for Models")
legend_text <- c()
for (model_name in names(final_models)) {
  pred_prob <- data_train[[paste0(model_name, "_prob")]]
  roc_obj <- roc(data_train$group, pred_prob)
  auc_value <- auc(roc_obj)
  
  lines(1 - roc_obj$specificities, roc_obj$sensitivities, col = colors[model_name], lwd = 2)
  
  legend_text <- c(legend_text, paste0(model_name, " (AUC = ", round(auc_value, 3), ")"))
}

legend("bottomright", legend = legend_text, col = colors[names(final_models)], lwd = 2, bty = "n")
dev.off()


##验证集--------
model_results <- data.frame(
  Model = character(),
  AUC = numeric(),
  Accuracy = numeric(),
  Sensitivity = numeric(),
  Precision = numeric(),
  F1_Score = numeric(),
  PRAUC = numeric(),
  stringsAsFactors = FALSE
)
confusion <- list()
# 循环遍历每个模型并计算预测结果和评估指标
for (model_name in names(final_models)) {
  # 提取模型的预测结果
  pred_prob <- predict(final_models[[model_name]], data_test, type = "prob")[,2]
  predictions <- predict(final_models[[model_name]], data_test)
  predictions_class <- predictions$.pred_class
  aucvalue <- auc(roc(data_test$group, pred_prob$.pred_Asthma))
  # 计算混淆矩阵
  conf_matrix <- ConfusionMatrix(as.factor(predictions_class), as.factor(data_test$group))
  confusion[[model_name]] <- conf_matrix
  # 计算其他指标
  accuracy <- Accuracy(predictions_class, data_test$group)
  sensitivity <- Sensitivity(predictions_class, data_test$group)
  precision <- Precision(predictions_class, data_test$group)
  f1 <- F1_Score(predictions_class, data_test$group)
  # 计算 PR 曲线并捕获可能的错误
  prauc_value <- tryCatch({
    pr <- pr.curve(scores.class0 = predictions_class[data_test$group == 'Asthma'],
                   scores.class1 = predictions_class[data_test$group == 'Control'],
                   curve = TRUE)
    pr$auc.integral
  }, error = function(e) {
    0
  })
  # 将结果存储到数据框中
  model_results <- rbind(model_results, data.frame(
    Model = model_name,
    AUC = aucvalue,
    Accuracy = accuracy,
    Sensitivity = sensitivity,
    Precision = precision,
    F1_Score = f1,
    PRAUC = prauc_value
  ))
  data_test[[paste0(model_name, "_prob")]] <- pred_prob$.pred_Asthma
}



write.csv(model_results, "test_model_evaluation_results.csv", row.names = FALSE,quote = F)

library(ggplot2)
library(reshape2)

model_results[, c(2:7)] <- lapply(model_results[, c(2:7)], function(x) as.numeric(unlist(x)))
results_long <-  melt(model_results, id.vars = "Model")
colnames(results_long)
str(model_results)
p <- ggplot(results_long, aes(x = variable, y = value, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  labs(x = "", y = "", title = "Test Cohort: Comparison of Model Performance") +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  scale_fill_manual(values = c("LR" = "#F781BF",
                               "RF" = "#b5aa0f",
                               "ENet" = "#984EA3",
                               'XGB' = '#CD534CFF',
                               'naive_Bayes' = '#3B4992FF',
                               'Adaboost' = '#8DD3C7'))
p
ggsave(filename = '03.MLPerformance_Test.pdf',p,width = 7,height = 4)
ggsave(filename = '03.MLPerformance_Test.png',p,width = 7,height = 4,units = 'in',dpi = 600)

library(ggplot2)
library(pROC)

colors <- c("LR" = "#F781BF",
            "RF" = "#b5aa0f",
            "ENet" = "#984EA3",
            'XGB' = '#CD534CFF',
            'naive_Bayes' = '#3B4992FF',
            'Adaboost' = '#8DD3C7')

pdf(file = '04.ROC_test.pdf', width = 5, height = 5.5, onefile = FALSE)
plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1), xlab = "1 - Specificity", ylab = "Sensitivity",
     main = "Test Cohort: ROC Curves for Models")
legend_text <- c()
for (model_name in names(final_models)) {
  pred_prob <- data_test[[paste0(model_name, "_prob")]]
  roc_obj <- roc(data_test$group, pred_prob)
  auc_value <- auc(roc_obj)
  
  lines(1 - roc_obj$specificities, roc_obj$sensitivities, col = colors[model_name], lwd = 2)
  
  legend_text <- c(legend_text, paste0(model_name, " (AUC = ", round(auc_value, 3), ")"))
}

legend("bottomright", legend = legend_text, col = colors[names(final_models)], lwd = 2, bty = "n")
dev.off()

png(file = '04.ROC_test.png', width = 5, height = 5.5, units = 'in', bg = 'white', res = 600)
plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1), xlab = "1 - Specificity", ylab = "Sensitivity",
     main = "Test Cohort: ROC Curves for Models")
legend_text <- c()
for (model_name in names(final_models)) {
  pred_prob <- data_test[[paste0(model_name, "_prob")]]
  roc_obj <- roc(data_test$group, pred_prob)
  auc_value <- auc(roc_obj)
  
  lines(1 - roc_obj$specificities, roc_obj$sensitivities, col = colors[model_name], lwd = 2)
  
  legend_text <- c(legend_text, paste0(model_name, " (AUC = ", round(auc_value, 3), ")"))
}

legend("bottomright", legend = legend_text, col = colors[names(final_models)], lwd = 2, bty = "n")
dev.off()

##SHAP---------
fit_rf = extract_fit_parsnip(final_models$RF) 
rf_feats <- names(fit_rf$fit$variable.importance)

library(fastshap)
X_train_df <- as.data.frame(data_train[, rf_feats])  # 特征


pfun <- function(object, newdata) {
  preds <- predict(object, data = newdata, type = "response")
  return(preds$predictions[, 2])
}


shap_values <- fastshap::explain(
  object = fit_rf$fit,
  X = X_train_df,
  pred_wrapper = pfun,
  nsim = 100,  # 可调整模拟次数
  feature_names = colnames(X_train_df)
)

boxplot(shap_values, las = 2, main = "SHAP Values by Feature")
library(reshape2)
library(ggplot2)

shap_long <- melt(shap_values)
colnames(shap_long) <- c("Sample", "Feature", "SHAP")

shap_mean <- aggregate(SHAP ~ Feature, data = shap_long, FUN = function(x) mean(abs(x)))
shap_long$Feature <- factor(shap_long$Feature, levels = shap_mean$Feature[order(shap_mean$SHAP)])


library(fmsb)

shap_means <- apply(shap_values, 2, function(x) mean(abs(x)))

top_n <- 11
top_genes <- sort(shap_means, decreasing = TRUE)[1:top_n]

radar_data <- rbind(rep(max(top_genes), top_n),   # max range
                    rep(0, top_n),                # min range
                    as.numeric(top_genes))        # 实际数据
colnames(radar_data) <- names(top_genes)
rownames(radar_data) <- c("max", "min", "SHAP")
class(radar_data)
radar_data <- as.data.frame(radar_data)
pdf(file = '05.imp.radar.pdf',width = 6,height = 6)
radarchart(radar_data,
           axistype = 1,
           pcol = "#F0027F", pfcol = rgb(1, 0, 0, 0.3), plwd = 2,
           cglcol = "grey", cglty = 1, axislabcol = "black",
           vlcex = 0.8, title = "Features SHAP Importance")
dev.off()

png(file = '05.imp.radar.png',width = 6,height = 6,units = 'in',res = 300,bg = 'white')
radarchart(radar_data,
           axistype = 1,
           pcol = "#F0027F", pfcol = rgb(1, 0, 0, 0.3), plwd = 2,
           cglcol = "grey", cglty = 1, axislabcol = "black",
           vlcex = 0.8, title = "Features SHAP Importance")
dev.off()

p <- ggplot(shap_long, aes(x = Feature, y = SHAP)) +
  geom_boxplot(fill = "#8DD3C7") +
  coord_flip() +
  labs(title = "SHAP Value Distribution by Gene",
       y = "SHAP value", x = "Gene") +
  theme_minimal()
p
ggsave(filename = '06.imp_shap_rf.pdf',p,w=6,h=5)
ggsave(filename = '06.imp_shap_rf.png',p,w=6,h=5,units = 'in',dpi = 300,bg = 'white')


top_genes <- sort(shap_means, decreasing = TRUE)[1:5]
hubgene <- data.frame(symbol = names(top_genes))
write.csv(hubgene,file = 'hubgene.csv',row.names = F,quote = F)

shapvalue <- as.data.frame(shap_mean)
write.csv(shapvalue,file = 'shapvalue.rf.csv',row.names = F,quote = F)
