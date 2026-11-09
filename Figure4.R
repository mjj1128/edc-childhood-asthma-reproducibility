# Figure 4：原样式的四模型观测值与预测值散点图
# 读取 20261107_200813 的分组验证结果，无需重新训练模型。

if (.Platform$OS.type == "windows") {
  invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
}

run_dir <- paste0(
  "D:/环境/meta数据csv/EDC_Overall_AnalysisA_SVG/三水平出图/",
  "补充分组验证_run_20261107_200813"
)
output_dir <- file.path(run_dir, "Figure4_分组验证")

project_dir <- "E:/codex 1edc/manuscript_review/ml_minimal_revision"
local_libs <- c(
  file.path(project_dir, "Rlib"),
  file.path(dirname(project_dir), "revision_existing_data", "Rlib")
)
.libPaths(c(local_libs[dir.exists(local_libs)], .libPaths()))

for (pkg in c("caret", "ranger", "glmnet", "xgboost", "ggplot2", "svglite")) {
  loadNamespace(pkg)
}

suppressPackageStartupMessages({
  library(caret)
  library(ggplot2)
})

# 1. 读取已完成运行的模型、预处理规则和预测值
stopifnot(identical(
  trimws(readLines(file.path(run_dir, "run_status.txt"), warn = FALSE)),
  "COMPLETED"
))

bundle <- readRDS(file.path(run_dir, "models_and_recipes.rds"))

input <- read.csv(
  file.path(run_dir, "analysis_input_used.csv"),
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM",
  stringsAsFactors = FALSE
)

test_predictions <- read.csv(
  file.path(run_dir, "grouped_predictions.csv"),
  fileEncoding = "UTF-8-BOM",
  stringsAsFactors = FALSE
)

stopifnot(
  !anyDuplicated(input$Effect_ID),
  !anyDuplicated(bundle$train_effect_ids),
  !anyDuplicated(bundle$test_effect_ids),
  !length(intersect(bundle$train_effect_ids, bundle$test_effect_ids))
)

train_idx <- match(bundle$train_effect_ids, input$Effect_ID)
test_idx <- match(bundle$test_effect_ids, input$Effect_ID)

stopifnot(
  length(train_idx) == 174L,
  length(test_idx) == 37L,
  !anyNA(train_idx),
  !anyNA(test_idx)
)

train <- input[train_idx, , drop = FALSE]
test <- input[test_idx, , drop = FALSE]

# 2. 应用已保存的预处理规则
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

x_train <- apply_recipe(train, bundle$recipe)
x_test <- apply_recipe(test, bundle$recipe)

models <- bundle$models
models$XGBoost <- xgboost::xgb.load(file.path(run_dir, bundle$xgboost_file))

predict_saved <- function(model, x) {
  fit <- models[[model]]
  
  if (model == "RandomForest") {
    return(as.numeric(predict(fit, data = x)$predictions))
  }
  
  if (model == "XGBoost") {
    return(as.numeric(predict(fit, xgboost::xgb.DMatrix(as.matrix(x)))))
  }
  
  as.numeric(predict(
    fit,
    newx = as.matrix(x),
    s = bundle$parameters[[model]]$lambda
  ))
}

metrics_label <- function(label, observed, predicted) {
  r2 <- 1 - sum((observed - predicted)^2) /
    sum((observed - mean(observed))^2)
  rmse <- sqrt(mean((observed - predicted)^2))
  
  sprintf("%s\nR² = %.2f\nRMSE = %.2f", label, r2, rmse)
}

# 3. 按原 Figure 4 的样式画单个模型面板
make_panel <- function(model, display_name) {
  saved_test <- test_predictions[
    test_predictions$Model == model, , drop = FALSE
  ]
  
  stopifnot(
    nrow(saved_test) == nrow(test),
    !anyDuplicated(saved_test$Effect_ID),
    setequal(saved_test$Effect_ID, test$Effect_ID)
  )
  
  saved_test <- saved_test[
    match(test$Effect_ID, saved_test$Effect_ID), , drop = FALSE
  ]
  
  # 核对当前读取的模型确实对应这份留出预测。
  stopifnot(
    max(abs(saved_test$Observed - test$ln_aOR)) < 1e-10,
    max(abs(saved_test$Predicted - predict_saved(model, x_test))) < 1e-6
  )
  
  train_pred <- predict_saved(model, x_train)
  
  plot_data <- rbind(
    data.frame(
      Observed = train$ln_aOR, Predicted = train_pred, Set = "Train"
    ),
    data.frame(
      Observed = saved_test$Observed,
      Predicted = saved_test$Predicted,
      Set = "Test"
    )
  )
  plot_data$Set <- factor(plot_data$Set, levels = c("Train", "Test"))
  
  extent <- range(c(plot_data$Observed, plot_data$Predicted))
  span <- diff(extent)
  limits <- extent + c(-0.10, 0.10) * span
  colours <- c(Train = "black", Test = "#E64B35")
  
  ggplot(plot_data, aes(x = Observed, y = Predicted)) +
    geom_abline(
      intercept = 0, slope = 1,
      linetype = "dashed", colour = "#8A2BE2", linewidth = 0.85
    ) +
    geom_point(
      aes(colour = Set, fill = Set),
      shape = 21, size = 2.5, alpha = 0.82, stroke = 0.35
    ) +
    scale_colour_manual(values = colours) +
    scale_fill_manual(values = colours) +
    coord_fixed(ratio = 1, xlim = limits, ylim = limits, clip = "off") +
    annotate(
      "label",
      x = limits[1] + span * 0.05,
      y = limits[2] - span * 0.08,
      label = metrics_label("Train", train$ln_aOR, train_pred),
      hjust = 0, vjust = 1, size = 3.2, linewidth = 0.25,
      fill = "white", colour = "black"
    ) +
    annotate(
      "label",
      x = limits[2] - span * 0.05,
      y = limits[1] + span * 0.08,
      label = metrics_label("Test", saved_test$Observed, saved_test$Predicted),
      hjust = 1, vjust = 0, size = 3.2, linewidth = 0.25,
      fill = "white", colour = "black"
    ) +
    annotate(
      "label",
      x = limits[1] + span * 0.07,
      y = limits[1] + span * 0.72,
      label = display_name,
      hjust = 0, vjust = 0.5, size = 4.2, fontface = "bold",
      linewidth = 0, fill = "#1F5EA8", colour = "white"
    ) +
    labs(x = "Observed ln(aOR)", y = "Predicted ln(aOR)") +
    theme_classic(base_size = 15, base_family = "sans") +
    theme(
      axis.title = element_text(face = "bold", colour = "black", size = 16),
      axis.text = element_text(colour = "black", size = 13),
      axis.line = element_line(colour = "black", linewidth = 0.7),
      panel.border = element_rect(fill = NA, colour = "#555555", linewidth = 0.85),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(colour = "black", size = 12),
      legend.background = element_blank(),
      legend.key = element_blank(),
      plot.margin = margin(10, 12, 10, 12),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    ) +
    guides(
      colour = guide_legend(
        nrow = 1, override.aes = list(size = 4, alpha = 0.95)
      ),
      fill = "none"
    )
}

# 4. 按原稿 A–D 的顺序输出四张 SVG
panel_names <- c(
  RandomForest = "Random Forest",
  XGBoost = "XGBoost",
  Ridge = "Ridge",
  ElasticNet = "Elastic Net"
)

panels <- setNames(lapply(names(panel_names), function(model) {
  make_panel(model, panel_names[[model]])
}), names(panel_names))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(dir.exists(output_dir))

for (i in seq_along(panels)) {
  file <- file.path(
    output_dir,
    paste0("Figure4", LETTERS[i], "_", names(panels)[i], "_grouped.svg")
  )
  
  ggsave(
    filename = file,
    plot = panels[[i]],
    width = 5.6,
    height = 5.5,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )
}

if (interactive()) print(panels$RandomForest)
cat("Figure 4 的四张 SVG 已保存到：", output_dir, "\n")

