# ============================================================
# Figure S10：XGBoost Kernel SHAP 蜂群图
# 用保存的分组模型解释测试集，从训练集抽取背景数据
# ============================================================

if (.Platform$OS.type == "windows") {
  invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
}

run_dir <- paste0(
  "D:/环境/meta数据csv/EDC_Overall_AnalysisA_SVG/三水平出图/",
  "补充分组验证_run_20261107_200813"
)
output_dir <- run_dir

project_dir <- "E:/codex 1edc/manuscript_review/ml_minimal_revision"
local_libs <- c(
  file.path(project_dir, "Rlib"),
  file.path(dirname(project_dir), "revision_existing_data", "Rlib")
)
.libPaths(c(local_libs[dir.exists(local_libs)], .libPaths()))

for (pkg in c(
  "caret", "xgboost", "ggplot2", "kernelshap", "ggbeeswarm", "svglite"
)) {
  loadNamespace(pkg)
}

suppressPackageStartupMessages({
  library(caret)
  library(ggplot2)
})

# 1. 读取这次保存的模型、效应量和测试集预测
stopifnot(identical(
  trimws(readLines(file.path(run_dir, "run_status.txt"), warn = FALSE)),
  "COMPLETED"
))

bundle <- readRDS(file.path(run_dir, "models_and_recipes.rds"))
model_path <- file.path(run_dir, bundle$xgboost_file)
stopifnot(file.exists(model_path))

input <- read.csv(
  file.path(run_dir, "analysis_input_used.csv"),
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM",
  stringsAsFactors = FALSE
)

test_pred <- read.csv(
  file.path(run_dir, "grouped_predictions.csv"),
  fileEncoding = "UTF-8-BOM",
  stringsAsFactors = FALSE
)
test_pred <- test_pred[test_pred$Model == "XGBoost", , drop = FALSE]

train_idx <- match(bundle$train_effect_ids, input$Effect_ID)
test_idx <- match(test_pred$Effect_ID, input$Effect_ID)

stopifnot(
  !anyDuplicated(input$Effect_ID),
  !anyDuplicated(test_pred$Effect_ID),
  length(train_idx) == 174L,
  nrow(test_pred) == 37L,
  !anyNA(train_idx),
  !anyNA(test_idx),
  !length(intersect(train_idx, test_idx)),
  setequal(test_pred$Effect_ID, bundle$test_effect_ids),
  max(abs(test_pred$Observed - input$ln_aOR[test_idx])) < 1e-10
)

# 2. 应用模型训练时保存的预处理规则
apply_recipe <- function(data, recipe) {
  data <- as.data.frame(data[, recipe$features, drop = FALSE])
  
  for (feature in intersect(names(recipe$medians), names(data))) {
    values <- suppressWarnings(as.numeric(data[[feature]]))
    values[!is.finite(values)] <- recipe$medians[[feature]]
    data[[feature]] <- values
  }
  
  for (feature in recipe$categorical) {
    values <- as.character(data[[feature]])
    values[is.na(values) | values == ""] <- "Missing"
    values[!values %in% recipe$levels[[feature]]] <- "Missing"
    data[[feature]] <- factor(values, levels = recipe$levels[[feature]])
  }
  
  x <- as.data.frame(
    predict(recipe$dummy, newdata = data),
    check.names = FALSE
  )
  names(x) <- make.names(names(x), unique = TRUE)
  
  for (feature in setdiff(recipe$columns, names(x))) {
    x[[feature]] <- 0
  }
  
  x <- as.data.frame(
    predict(recipe$scale, x[, recipe$columns, drop = FALSE])
  )
  stopifnot(all(is.finite(as.matrix(x))))
  x
}

xgb <- xgboost::xgb.load(model_path)

predict_raw <- function(object, X, ...) {
  x <- apply_recipe(X, bundle$recipe)
  dm <- xgboost::xgb.DMatrix(as.matrix(x), nthread = 2L)
  as.numeric(predict(object, dm))
}

features <- bundle$features
train_raw <- input[train_idx, features, drop = FALSE]
test_raw <- input[test_idx, features, drop = FALSE]

# 核对载入的模型是否复现已保存的测试预测
check_test <- predict_raw(xgb, test_raw)
stopifnot(max(abs(check_test - test_pred$Predicted)) < 1e-6)

# 3. 从训练集抽取 100 条背景记录，计算 37 条测试记录的 SHAP
set.seed(2026)
bg_index <- sample(seq_len(nrow(train_raw)), min(100L, nrow(train_raw)))
bg <- train_raw[bg_index, , drop = FALSE]

set.seed(2026)
shap <- kernelshap::kernelshap(
  object = xgb,
  X = test_raw,
  bg_X = bg,
  pred_fun = predict_raw
)

S <- as.matrix(shap$S)

# 核对背景预测均值与 SHAP 贡献之和能还原模型预测
stopifnot(
  nrow(S) == 37L,
  setequal(colnames(S), features),
  all(is.finite(S)),
  max(abs(as.numeric(shap$baseline) + rowSums(S) - check_test)) < 1e-6
)

# 4. 准备重要性和蜂群图数据
importance <- data.frame(
  Feature = colnames(S),
  Importance = colMeans(abs(S))
)
importance <- importance[
  order(importance$Importance, decreasing = TRUE), , drop = FALSE
]

feature_order <- rev(importance$Feature)
numeric_features <- names(bundle$recipe$medians)

# 数值变量按各自的测试集范围映射到 0–1。
# 分类变量和缺失值用灰色。
feature_colour <- function(feature) {
  if (!feature %in% numeric_features) {
    return(rep(NA_real_, nrow(test_raw)))
  }
  
  values <- suppressWarnings(as.numeric(test_raw[[feature]]))
  usable <- is.finite(values)
  
  if (!any(usable)) {
    return(rep(NA_real_, length(values)))
  }
  
  extent <- range(values[usable])
  values[!usable] <- NA_real_
  
  if (diff(extent) == 0) {
    values[usable] <- 0.5
  } else {
    values[usable] <- (values[usable] - extent[1]) / diff(extent)
  }
  
  values
}

plot_data <- do.call(rbind, lapply(colnames(S), function(feature) {
  data.frame(
    Effect_ID = test_pred$Effect_ID,
    Study_ID = test_pred$Study_ID,
    Feature = feature,
    SHAP = S[, feature],
    Feature_value = feature_colour(feature),
    row.names = NULL
  )
}))

plot_data$Feature <- factor(plot_data$Feature, levels = feature_order)
importance$y_num <- match(importance$Feature, feature_order)

extent <- range(c(0, plot_data$SHAP))
span <- diff(extent)
stopifnot(span > 0, max(importance$Importance) > 0)

bar_start <- extent[1] - 0.08 * span
bar_scale <- 0.95 * span / max(importance$Importance)
importance$bar_end <- bar_start + importance$Importance * bar_scale

# 5. 按原 S10 布局绘制蜂群图
p_shap <- ggplot() +
  geom_rect(
    data = importance,
    aes(
      xmin = bar_start, xmax = bar_end,
      ymin = y_num - 0.275, ymax = y_num + 0.275
    ),
    fill = "#87B1D4",
    alpha = 0.70
  ) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.6) +
  ggbeeswarm::geom_quasirandom(
    data = plot_data,
    aes(x = SHAP, y = Feature, colour = Feature_value),
    orientation = "y",
    method = "pseudorandom",
    width = 0.35,
    alpha = 0.85,
    size = 1.6,
    stroke = 0
  ) +
  scale_colour_gradient(
    low = "#4575B4",
    high = "#D73027",
    na.value = "#666666",
    limits = c(0, 1),
    breaks = c(0, 1),
    labels = c("Low", "High"),
    name = "Numeric feature value"
  ) +
  scale_x_continuous(
    name = "SHAP value (impact on model output)",
    limits = c(bar_start, max(extent[2], importance$bar_end)),
    expand = expansion(mult = c(0, 0.05)),
    sec.axis = sec_axis(
      ~ (. - bar_start) / bar_scale,
      name = "Mean (|SHAP value|)"
    )
  ) +
  labs(y = "Feature", title = "XGBoost SHAP Summary") +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(
      colour = "grey92", linetype = "dashed", linewidth = 0.5
    ),
    axis.line.x.bottom = element_line(colour = "black", linewidth = 0.8),
    axis.line.x.top = element_line(colour = "black", linewidth = 0.8),
    axis.line.y.left = element_line(colour = "black", linewidth = 0.8),
    axis.text.y = element_text(colour = "black", size = 13),
    axis.text.x = element_text(colour = "black", size = 12, face = "bold"),
    axis.ticks = element_line(colour = "black"),
    legend.position = "right",
    legend.margin = margin(l = 10),
    axis.title.y = element_text(
      margin = margin(r = 15), face = "bold", size = 14
    ),
    axis.title.x.bottom = element_text(
      margin = margin(t = 15), face = "bold", size = 14
    ),
    axis.title.x.top = element_text(
      margin = margin(b = 15), face = "bold", size = 14
    ),
    plot.title = element_text(
      hjust = 0.5, face = "bold", size = 22, colour = "black"
    ),
    plot.margin = margin(t = 20, r = 50, b = 20, l = 120)
  )

# 6. 保存 SVG 和计算结果
output_file <- file.path(output_dir, "FigureS10_XGBoost_kernelshap_grouped.svg")

set.seed(2026)
ggsave(
  output_file,
  plot = p_shap,
  width = 12.5,
  height = 7.5,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

write.csv(
  importance[, c("Feature", "Importance")],
  file.path(output_dir, "FigureS10_XGBoost_mean_abs_SHAP.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  plot_data,
  file.path(output_dir, "FigureS10_XGBoost_SHAP_values.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  data.frame(Effect_ID = bundle$train_effect_ids[bg_index]),
  file.path(output_dir, "FigureS10_XGBoost_background_ids.csv"),
  row.names = FALSE
)

saveRDS(shap, file.path(output_dir, "FigureS10_XGBoost_kernelshap.rds"))

if (interactive()) {
  set.seed(2026)
  print(p_shap)
}

print(importance[, c("Feature", "Importance")], row.names = FALSE)
cat("Figure S10 已保存到：", output_file, "\n")

