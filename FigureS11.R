# ============================================================
# Figure S11：XGBoost 部分依赖图
# A. Mean age
# B. Male proportion
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

for (pkg in c("caret", "xgboost", "ggplot2", "svglite")) {
  loadNamespace(pkg)
}

suppressPackageStartupMessages({
  library(caret)
  library(ggplot2)
})

# 1. 读取这次保存的模型、原始输入和留出预测
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

xgb <- xgboost::xgb.load(model_path)

predict_raw <- function(data) {
  x <- apply_recipe(data, bundle$recipe)
  dm <- xgboost::xgb.DMatrix(as.matrix(x), nthread = 2L)
  as.numeric(predict(xgb, dm))
}

# 先核对载入的模型能复现这次保存的测试预测
check_test <- predict_raw(input[test_idx, , drop = FALSE])
stopifnot(max(abs(check_test - test_pred$Predicted)) < 1e-6)

test_raw <- input[test_idx, bundle$features, drop = FALSE]
train_raw <- input[train_idx, bundle$features, drop = FALSE]

# 3. 计算单个变量的 PDP 和逐条预测分位带
calc_pdp <- function(feature) {
  # NR 等非数值记录不参与确定横轴范围。
  observed_test <- suppressWarnings(as.numeric(test_raw[[feature]]))
  observed_test <- observed_test[is.finite(observed_test)]
  
  observed_train <- suppressWarnings(as.numeric(train_raw[[feature]]))
  observed_train <- observed_train[is.finite(observed_train)]
  
  stopifnot(length(observed_test) >= 5L, length(observed_train) >= 5L)
  
  # 测试集第 2.5–97.5 百分位，限制在训练集实际观察范围内。
  limits <- as.numeric(quantile(observed_test, c(0.025, 0.975)))
  limits <- c(
    max(limits[1], min(observed_train)),
    min(limits[2], max(observed_train))
  )
  stopifnot(limits[1] < limits[2])
  
  grid <- seq(limits[1], limits[2], length.out = 60L)
  
  # 每个横轴值依次替换目标变量，保留每条记录的其余变量。
  prediction_matrix <- vapply(
    grid,
    function(value) {
      modified <- test_raw
      modified[[feature]] <- value
      predict_raw(modified)
    },
    numeric(nrow(test_raw))
  )
  stopifnot(all(is.finite(prediction_matrix)))
  
  interval <- apply(
    prediction_matrix,
    2L,
    quantile,
    probs = c(0.10, 0.90)
  )
  
  data.frame(
    Feature = feature,
    Value = grid,
    mean_yhat = colMeans(prediction_matrix),
    lower = interval[1, ],
    upper = interval[2, ],
    Reference_effects = nrow(test_raw),
    Observed_test_values = length(observed_test)
  )
}

# 4. 沿用原图的蓝色曲线、圆点和浅蓝色带
plot_pdp <- function(curve, title, axis_label) {
  ggplot(curve, aes(x = Value)) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      fill = "#DCE9F4",
      alpha = 0.78
    ) +
    geom_line(
      aes(y = mean_yhat),
      colour = "#3F7FB3",
      linewidth = 1.35
    ) +
    geom_point(
      aes(y = mean_yhat),
      colour = "white",
      fill = "#4F8FC2",
      shape = 21,
      size = 2.7,
      stroke = 0.35,
      alpha = 0.95
    ) +
    labs(
      title = title,
      x = axis_label,
      y = "Model-predicted ln(aOR)"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(
        face = "bold", size = 21, hjust = 0.5, colour = "black",
        margin = margin(b = 12)
      ),
      axis.title = element_text(face = "bold", colour = "black", size = 15),
      axis.text = element_text(colour = "black", size = 12),
      axis.line = element_line(colour = "black", linewidth = 0.65),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.45),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        colour = "#555555", fill = NA, linewidth = 0.85
      ),
      plot.margin = margin(t = 18, r = 28, b = 22, l = 34),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

# 5. 生成 S11A 和 S11B
age <- calc_pdp("Mean age")
male <- calc_pdp("Male proportion")

p_age <- plot_pdp(age, "Mean age", "Mean age (years)")
p_male <- plot_pdp(male, "Male proportion", "Male proportion (0–1)")

save_panel <- function(plot, name) {
  ggsave(
    filename = file.path(output_dir, name),
    plot = plot,
    width = 7.4,
    height = 5.2,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )
}

save_panel(p_age, "FigureS11A_XGBoost_PDP_Mean_age.svg")
save_panel(p_male, "FigureS11B_XGBoost_PDP_Male_proportion.svg")

write.csv(
  rbind(age, male),
  file.path(output_dir, "FigureS11_XGBoost_PDP_values.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

if (interactive()) {
  print(p_age)
  print(p_male)
}

cat("Figure S11 已保存到：", output_dir, "\n")

