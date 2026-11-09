# ============================================================
# Figure S9：XGBoost 预测与残差图
# 读取单次分组验证保存的模型，沿用原 S9 的灰蓝配色和布局
# ============================================================

if (.Platform$OS.type == "windows") {
  invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
}

run_dir <- paste0(
  "D:/环境/meta数据csv/EDC_Overall_AnalysisA_SVG/三水平出图/",
  "补充分组验证_run_20261107_200813"
)

project_dir <- "E:/codex 1edc/manuscript_review/ml_minimal_revision"
local_libs <- c(
  file.path(project_dir, "Rlib"),
  file.path(dirname(project_dir), "revision_existing_data", "Rlib")
)
.libPaths(c(local_libs[dir.exists(local_libs)], .libPaths()))

for (pkg in c("caret", "xgboost", "ggplot2", "ggExtra", "cowplot", "svglite")) {
  loadNamespace(pkg)
}

suppressPackageStartupMessages({
  library(caret)
  library(ggplot2)
})

# 1. 读取已保存的模型、预处理规则和测试预测
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

# 2. 应用训练时保存的预处理规则
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

x_train <- apply_recipe(input[train_idx, , drop = FALSE], bundle$recipe)
x_test <- apply_recipe(input[test_idx, , drop = FALSE], bundle$recipe)

xgb <- xgboost::xgb.load(model_path)

train_pred <- as.numeric(predict(
  xgb, xgboost::xgb.DMatrix(as.matrix(x_train))
))

check_test <- as.numeric(predict(
  xgb, xgboost::xgb.DMatrix(as.matrix(x_test))
))

# 确认载入的 XGBoost 模型与这次保存的留出预测一致
stopifnot(max(abs(check_test - test_pred$Predicted)) < 1e-6)

# 3. 准备绘图数据和指标
plot_data <- rbind(
  data.frame(
    Observed = input$ln_aOR[train_idx],
    Predicted = train_pred,
    Set = "Train data"
  ),
  data.frame(
    Observed = test_pred$Observed,
    Predicted = test_pred$Predicted,
    Set = "Test data"
  )
)

plot_data$Set <- factor(plot_data$Set, levels = c("Train data", "Test data"))
plot_data$Residual <- plot_data$Observed - plot_data$Predicted

train_mae <- mean(abs(plot_data$Residual[plot_data$Set == "Train data"]))
test_error <- test_pred$Observed - test_pred$Predicted
test_mae <- mean(abs(test_error))
test_rmse <- sqrt(mean(test_error^2))
test_r2 <- 1 - sum(test_error^2) /
  sum((test_pred$Observed - mean(test_pred$Observed))^2)

colours <- c("Train data" = "#B8BCC2", "Test data" = "#6F9FCF")

theme_ml <- theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_line(
      colour = "grey85", linetype = "dashed", linewidth = 0.45
    ),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.1),
    axis.text = element_text(colour = "black", size = 14),
    axis.title = element_text(colour = "black", face = "bold", size = 15),
    legend.title = element_blank(),
    legend.text = element_text(colour = "black", size = 13),
    legend.background = element_rect(
      colour = "black", fill = "white", linewidth = 0.6
    ),
    plot.margin = margin(t = 10, r = 28, b = 35, l = 38)
  )

extent <- range(c(plot_data$Observed, plot_data$Predicted))
span <- diff(extent)
limits <- extent + c(-0.08, 0.08) * span

# 4. 上图：观测值、预测值和边缘分布
p_prediction <- ggplot(plot_data, aes(Observed, Predicted)) +
  geom_point(
    aes(colour = Set, fill = Set),
    shape = 21, size = 2.8, alpha = 0.78, stroke = 0.35
  ) +
  geom_abline(
    intercept = 0, slope = 1,
    linetype = "dashed", colour = "black", linewidth = 0.85
  ) +
  scale_colour_manual(values = colours) +
  scale_fill_manual(values = colours) +
  coord_fixed(ratio = 1, xlim = limits, ylim = limits) +
  labs(title = "XGBoost", x = NULL, y = NULL) +
  theme_ml +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
    legend.position = "inside",
    legend.position.inside = c(0.78, 0.18)
  ) +
  annotate(
    "text",
    x = limits[1] + span * 0.05,
    y = limits[2] - span * 0.06,
    label = sprintf("Test R² = %.3f", test_r2),
    hjust = 0, vjust = 1, size = 5.5,
    fontface = "bold", colour = "#A64040"
  ) +
  annotate(
    "text",
    x = limits[1] + span * 0.05,
    y = limits[2] - span * 0.15,
    label = sprintf("Test RMSE = %.3f", test_rmse),
    hjust = 0, vjust = 1, size = 5.5,
    fontface = "bold", colour = "#A64040"
  )

p_marginal <- ggExtra::ggMarginal(
  p_prediction,
  groupColour = TRUE,
  groupFill = TRUE,
  type = "densigram",
  alpha = 0.35,
  size = 4
)

# 5. 下图：残差
residual_extent <- range(plot_data$Residual)
residual_span <- diff(residual_extent)
right <- max(plot_data$Predicted)

p_residual <- ggplot(plot_data, aes(Predicted, Residual)) +
  geom_point(
    aes(colour = Set, fill = Set),
    shape = 21, size = 2.6, alpha = 0.78, stroke = 0.35
  ) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.9) +
  scale_colour_manual(values = colours) +
  scale_fill_manual(values = colours) +
  labs(x = "Predicted ln(aOR)", y = "Residual") +
  theme_ml +
  theme(legend.position = "none") +
  annotate(
    "text",
    x = right,
    y = residual_extent[2] - residual_span * 0.05,
    label = sprintf("Train MAE = %.3f", train_mae),
    hjust = 1, vjust = 1, size = 4.8, fontface = "bold"
  ) +
  annotate(
    "text",
    x = right,
    y = residual_extent[2] - residual_span * 0.25,
    label = sprintf("Test MAE = %.3f", test_mae),
    hjust = 1, vjust = 1, size = 4.8, fontface = "bold"
  )

# 6. 拼图并按原尺寸保存为 SVG
figure_s9 <- cowplot::plot_grid(
  p_marginal,
  p_residual,
  ncol = 1,
  rel_heights = c(4, 1.25),
  align = "v",
  axis = "lr"
)

output_file <- file.path(
  run_dir,
  "FigureS9_XGBoost_prediction_residual_grouped.svg"
)

ggsave(
  output_file,
  plot = figure_s9,
  width = 8.5,
  height = 11,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

if (interactive()) print(figure_s9)

cat(
  sprintf(
    "\nXGBoost Test R² = %.3f，RMSE = %.3f，MAE = %.3f\n",
    test_r2, test_rmse, test_mae
  )
)
cat("Figure S9 已保存到：", output_file, "\n")

