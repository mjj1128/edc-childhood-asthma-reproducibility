# ============================================================
# 单次研究分组验证
# 28 个研究组训练，7 个研究组测试
# 训练集内部 10 折分组验证
# ============================================================

if (.Platform$OS.type == "windows") {
  invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
  invisible(Sys.setlocale("LC_COLLATE", "Chinese (Simplified)_China.utf8"))
}

run_grouped_holdout <- function() {
  
  # ----------------------------------------------------------
  # 1. 路径、依赖和固定设置
  # ----------------------------------------------------------
  
  project_dir <- "E:/codex 1edc/manuscript_review/ml_minimal_revision"
  output_base <- "D:/环境/meta数据csv/EDC_Overall_AnalysisA_SVG/三水平出图"
  
  local_libraries <- c(
    file.path(project_dir, "Rlib"),
    file.path(dirname(project_dir), "revision_existing_data", "Rlib")
  )
  .libPaths(c(local_libraries[dir.exists(local_libraries)], .libPaths()))
  
  packages <- c("caret", "ranger", "glmnet", "xgboost", "metafor")
  missing_packages <- setdiff(
    packages,
    rownames(utils::installed.packages())
  )
  
  if (length(missing_packages)) {
    stop(
      "缺少 R 包：", paste(missing_packages, collapse = ", "),
      "。请先在重启后的 R 会话中安装所缺的包。",
      call. = FALSE
    )
  }
  
  for (package in packages) {
    tryCatch(
      loadNamespace(package),
      error = function(e) {
        stop(
          "无法加载 R 包 ", package, "：", conditionMessage(e),
          call. = FALSE
        )
      }
    )
  }
  
  # caret 的 dummyVars 在处理对比编码时需要加载到搜索路径。
  suppressPackageStartupMessages(library(caret))
  
  ml_path <- file.path(project_dir, "inputs", "input_ml.csv")
  meta_path <- file.path(project_dir, "inputs", "input_meta.csv")
  population_path <- file.path(project_dir, "inputs", "population_groups.csv")
  
  if (!file.exists(ml_path)) {
    ml_path <- paste0(
      "E:/作业缪娟/2026-09-14/",
      "referenced-chatgpt-conversation-this-is-an/outputs/",
      "机器学习_修改_with_Exposure_Contrast.csv"
    )
  }
  
  if (!file.exists(meta_path)) {
    meta_path <- "D:/环境/meta数据csv/三水平_meta_AnalysisA_最终版.csv"
  }
  
  if (!file.exists(ml_path) || !file.exists(meta_path)) {
    stop("找不到输入文件。请核对 ml_path 和 meta_path。", call. = FALSE)
  }
  
  split_seed <- 42L
  cv_seed <- 43L
  model_seed <- 42L
  inner_k <- 10L
  threads <- 2L
  conservative_CEAS <- TRUE
  
  output_dir <- file.path(
    output_base,
    paste0("补充分组验证_run_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  )
  
  if (dir.exists(output_dir)) {
    stop("输出目录已存在，请稍后再运行。", call. = FALSE)
  }
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (!dir.exists(output_dir)) {
    stop("无法创建输出目录：", output_dir, call. = FALSE)
  }
  
  status_path <- file.path(output_dir, "run_status.txt")
  completed <- FALSE
  on.exit(
    if (!completed) writeLines("FAILED: partial outputs", status_path),
    add = TRUE
  )
  writeLines("RUNNING", status_path)
  
  save_csv <- function(data, filename) {
    utils::write.csv(
      data,
      file.path(output_dir, filename),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )
  }
  
  save_csv(
    data.frame(
      Input = c("ML", "Meta"),
      Path = c(ml_path, meta_path),
      MD5 = unname(tools::md5sum(c(ml_path, meta_path)))
    ),
    "input_checksums.csv"
  )
  
  # ----------------------------------------------------------
  # 2. 读取机器学习数据，并与 Meta 数据核对效应量
  # ----------------------------------------------------------
  
  numeric_features <- c(
    "Effective sample size",
    "Mean age",
    "Male proportion"
  )
  
  categorical_features <- c(
    "EDC class",
    "Chemical form",
    "Exposure window",
    "Study design",
    "Region",
    "Outcome group",
    "Exposure contrast"
  )
  
  features <- c(numeric_features, categorical_features)
  
  data <- read.csv(
    ml_path,
    check.names = FALSE,
    fileEncoding = "UTF-8-BOM",
    na.strings = c("", "NA")
  )
  names(data) <- trimws(names(data))
  
  ready_column <- if ("Analysis_Ready" %in% names(data)) {
    "Analysis_Ready"
  } else {
    "Analysis Ready"
  }
  
  required <- c(features, "ln_aOR", "SE", "Study-ID", "Effect_ID", ready_column)
  absent <- setdiff(required, names(data))
  if (length(absent)) {
    stop("输入数据缺少列：", paste(absent, collapse = ", "), call. = FALSE)
  }
  
  data$ln_aOR <- as.numeric(data$ln_aOR)
  data <- data[
    which(data[[ready_column]] == 1 & is.finite(data$ln_aOR)),
    , drop = FALSE
  ]
  data$Study_ID <- trimws(as.character(data[["Study-ID"]]))
  
  stopifnot(
    nrow(data) == 211L,
    length(unique(data$Study_ID)) == 36L,
    !anyDuplicated(data$Effect_ID),
    !anyNA(data$Effect_ID),
    !anyNA(data$Study_ID),
    all(nzchar(data$Study_ID)),
    all(is.finite(data$SE)),
    all(data$SE > 0)
  )
  
  meta <- read.csv(
    meta_path,
    check.names = FALSE,
    fileEncoding = "UTF-8-BOM",
    na.strings = c("", "NA", "NR", "N/A")
  )
  
  stopifnot(
    all(c("Effect_ID", "ln_aOR", "SE", "Sample_ID") %in% names(meta)),
    !anyDuplicated(meta$Effect_ID)
  )
  
  matched_meta <- match(data$Effect_ID, meta$Effect_ID)
  stopifnot(
    !anyNA(matched_meta),
    max(abs(data$ln_aOR - meta$ln_aOR[matched_meta])) < 1e-12,
    max(abs(data$SE - meta$SE[matched_meta])) < 1e-12
  )
  data$Sample_ID <- as.character(meta$Sample_ID[matched_meta])
  
  # ----------------------------------------------------------
  # 3. 外层分组留出与内层分组折
  # ----------------------------------------------------------
  
  group <- paste0("REPORT::", data$Study_ID)
  
  # 可选文件用于合并已核实的跨报告共同人群。
  # 文件列名为 Study_ID、Population_ID。
  if (file.exists(population_path)) {
    population <- read.csv(
      population_path,
      fileEncoding = "UTF-8-BOM",
      colClasses = "character"
    )
    
    stopifnot(
      all(c("Study_ID", "Population_ID") %in% names(population)),
      !anyDuplicated(population$Study_ID),
      !anyNA(population$Study_ID),
      !anyNA(population$Population_ID),
      all(nzchar(population$Study_ID)),
      all(nzchar(population$Population_ID)),
      all(population$Study_ID %in% data$Study_ID)
    )
    
    matched_population <- match(data$Study_ID, population$Study_ID)
    has_mapping <- !is.na(matched_population)
    group[has_mapping] <- paste0(
      "POPULATION::",
      population$Population_ID[matched_population[has_mapping]]
    )
  }
  
  # CEAS 报告按可能存在人群重叠保守合并。
  if (conservative_CEAS) {
    ceas_groups <- unique(group[grepl("CEAS", data$Sample_ID, ignore.case = TRUE)])
    if (length(ceas_groups) > 1L) {
      group[group %in% ceas_groups] <- "CEAS_potential_overlap_block"
    }
  }
  
  groups <- sort(unique(group))
  
  set.seed(split_seed)
  training_groups <- sample(groups, size = floor(0.8 * length(groups)))
  
  train_index <- which(group %in% training_groups)
  test_index <- which(!group %in% training_groups)
  train_group <- group[train_index]
  test_group <- group[test_index]
  
  stopifnot(
    length(train_index) > 0L,
    length(test_index) > 0L,
    length(training_groups) >= inner_k,
    !length(intersect(train_group, test_group))
  )
  
  train <- data[train_index, , drop = FALSE]
  test <- data[test_index, , drop = FALSE]
  
  set.seed(cv_seed)
  shuffled_groups <- sample(sort(unique(train_group)))
  fold_map <- setNames(
    rep(seq_len(inner_k), length.out = length(shuffled_groups)),
    shuffled_groups
  )
  inner_fold <- as.integer(fold_map[train_group])
  
  save_csv(
    data.frame(
      Effect_ID = data$Effect_ID,
      Study_ID = data$Study_ID,
      Sample_ID = data$Sample_ID,
      Validation_group = group,
      Set = ifelse(group %in% training_groups, "Train", "Test")
    ),
    "outer_groups.csv"
  )
  
  save_csv(
    data.frame(
      Effect_ID = train$Effect_ID,
      Study_ID = train$Study_ID,
      Validation_group = train_group,
      Inner_fold = inner_fold
    ),
    "inner_groups.csv"
  )
  
  cat(
    "外层分组完成。训练组：", length(unique(train_group)),
    "，测试组：", length(unique(test_group)),
    "，训练效应量：", nrow(train),
    "，测试效应量：", nrow(test), "\n",
    sep = ""
  )
  
  # ----------------------------------------------------------
  # 4. 预处理函数
  #    每个内层折的中位数、编码和缩放都只用该折训练数据拟合
  # ----------------------------------------------------------
  
  clean_data <- function(data, recipe) {
    result <- as.data.frame(data[, recipe$features, drop = FALSE])
    
    for (feature in numeric_features) {
      if (feature %in% names(result)) {
        value <- suppressWarnings(as.numeric(result[[feature]]))
        value[!is.finite(value)] <- recipe$medians[[feature]]
        result[[feature]] <- value
      }
    }
    
    for (feature in recipe$categorical) {
      value <- as.character(result[[feature]])
      value[is.na(value) | value == ""] <- "Missing"
      value[!value %in% recipe$levels[[feature]]] <- "Missing"
      result[[feature]] <- factor(value, levels = recipe$levels[[feature]])
    }
    
    result
  }
  
  fit_recipe <- function(data) {
    recipe <- list(
      features = features,
      categorical = categorical_features,
      medians = list(),
      levels = list()
    )
    
    for (feature in numeric_features) {
      value <- suppressWarnings(as.numeric(data[[feature]]))
      middle <- median(value[is.finite(value)], na.rm = TRUE)
      recipe$medians[[feature]] <- if (is.finite(middle)) middle else 0
    }
    
    constant_categorical <- character()
    for (feature in categorical_features) {
      value <- as.character(data[[feature]])
      value[is.na(value) | value == ""] <- "Missing"
      recipe$levels[[feature]] <- unique(c(value, "Missing"))
      
      if (length(unique(value)) < 2L) {
        constant_categorical <- c(constant_categorical, feature)
      }
    }
    
    recipe$features <- setdiff(recipe$features, constant_categorical)
    recipe$categorical <- setdiff(recipe$categorical, constant_categorical)
    
    formula <- as.formula(
      paste("~", paste(sprintf("`%s`", recipe$features), collapse = " + "))
    )
    environment(formula) <- baseenv()
    
    cleaned <- clean_data(data, recipe)
    recipe$dummy <- caret::dummyVars(formula, data = cleaned)
    
    matrix_data <- as.data.frame(
      predict(recipe$dummy, newdata = cleaned),
      check.names = FALSE
    )
    names(matrix_data) <- make.names(names(matrix_data), unique = TRUE)
    
    near_zero <- caret::nearZeroVar(matrix_data)
    if (length(near_zero)) {
      matrix_data <- matrix_data[, -near_zero, drop = FALSE]
    }
    
    if (ncol(matrix_data) < 2L) {
      stop("某个训练折保留下来的预测变量少于 2 个。", call. = FALSE)
    }
    
    recipe$columns <- names(matrix_data)
    recipe$scale <- caret::preProcess(
      matrix_data,
      method = c("center", "scale")
    )
    recipe
  }
  
  apply_recipe <- function(data, recipe) {
    matrix_data <- as.data.frame(
      predict(recipe$dummy, newdata = clean_data(data, recipe)),
      check.names = FALSE
    )
    names(matrix_data) <- make.names(names(matrix_data), unique = TRUE)
    
    for (feature in setdiff(recipe$columns, names(matrix_data))) {
      matrix_data[[feature]] <- 0
    }
    
    matrix_data <- as.data.frame(
      predict(recipe$scale, matrix_data[, recipe$columns, drop = FALSE])
    )
    stopifnot(all(is.finite(as.matrix(matrix_data))))
    matrix_data
  }
  
  outer_recipe <- fit_recipe(train)
  x_train <- apply_recipe(train, outer_recipe)
  x_test <- apply_recipe(test, outer_recipe)
  y_train <- train$ln_aOR
  
  fold_data <- lapply(seq_len(inner_k), function(fold) {
    in_train <- inner_fold != fold
    in_validation <- inner_fold == fold
    
    stopifnot(
      !length(intersect(train_group[in_train], train_group[in_validation]))
    )
    
    fold_train <- train[in_train, , drop = FALSE]
    fold_validation <- train[in_validation, , drop = FALSE]
    recipe <- fit_recipe(fold_train)
    
    list(
      x = apply_recipe(fold_train, recipe),
      y = fold_train$ln_aOR,
      validation_x = apply_recipe(fold_validation, recipe),
      validation_y = fold_validation$ln_aOR,
      recipe = recipe
    )
  })
  
  # ----------------------------------------------------------
  # 5. 评分函数与原模型参数范围
  # ----------------------------------------------------------
  
  score <- function(observed, predicted, weights = rep(1, length(observed))) {
    stopifnot(
      length(observed) == length(predicted),
      all(is.finite(observed)),
      all(is.finite(predicted))
    )
    
    residual <- observed - predicted
    total_ss <- sum(weights * (observed - weighted.mean(observed, weights))^2)
    
    c(
      RMSE = sqrt(weighted.mean(residual^2, weights)),
      MAE = weighted.mean(abs(residual), weights),
      R2 = if (total_ss > 0) {
        1 - sum(weights * residual^2) / total_ss
      } else {
        NA_real_
      }
    )
  }
  
  predictors <- ncol(x_train)
  mtry_values <- unique(pmax(
    1, c(floor(sqrt(predictors)), floor(predictors / 3), floor(predictors / 2))
  ))
  
  grids <- list(
    Ridge = expand.grid(
      lambda = 10^seq(-4, 2, length.out = 40)
    ),
    ElasticNet = expand.grid(
      alpha = c(0.25, 0.50, 0.75),
      lambda = 10^seq(-4, 1, length.out = 30)
    ),
    RandomForest = expand.grid(
      mtry = mtry_values,
      min.node.size = c(1, 2, 5, 10)
    ),
    XGBoost = expand.grid(
      max_depth = c(2, 3, 4),
      eta = c(0.03, 0.05, 0.10),
      subsample = c(0.7, 0.9),
      colsample_bytree = c(0.7, 0.9),
      min_child_weight = c(1, 3)
    )
  )
  
  fit_regular <- function(model, x, y, parameters, importance = "none") {
    if (model == "RandomForest") {
      return(ranger::ranger(
        x = x,
        y = y,
        num.trees = 1000,
        num.threads = threads,
        mtry = min(parameters$mtry, ncol(x)),
        min.node.size = parameters$min.node.size,
        seed = model_seed,
        importance = importance
      ))
    }
    
    glmnet::glmnet(
      as.matrix(x),
      y,
      alpha = if (model == "Ridge") 0 else parameters$alpha,
      lambda = parameters$lambda,
      standardize = FALSE
    )
  }
  
  predict_regular <- function(model, fit, x, parameters) {
    if (model == "RandomForest") {
      return(as.numeric(predict(fit, data = x)$predictions))
    }
    as.numeric(predict(fit, newx = as.matrix(x), s = parameters$lambda))
  }
  
  xgb_parameters <- function(parameters) {
    list(
      objective = "reg:squarederror",
      eval_metric = "rmse",
      max_depth = parameters$max_depth,
      eta = parameters$eta,
      subsample = parameters$subsample,
      colsample_bytree = parameters$colsample_bytree,
      min_child_weight = parameters$min_child_weight,
      seed = model_seed,
      nthread = threads
    )
  }
  
  xgb_fold_rmse <- function(fold, parameters) {
    train_matrix <- xgboost::xgb.DMatrix(as.matrix(fold$x), label = fold$y)
    validation_matrix <- xgboost::xgb.DMatrix(
      as.matrix(fold$validation_x),
      label = fold$validation_y
    )
    
    args <- list(
      params = xgb_parameters(parameters),
      data = train_matrix,
      nrounds = 300L,
      verbose = 0
    )
    
    evaluation <- list(validation = validation_matrix)
    if ("evals" %in% names(formals(xgboost::xgb.train))) {
      args$evals <- evaluation
    } else {
      args$watchlist <- evaluation
    }
    
    fit <- do.call(xgboost::xgb.train, args)
    evaluation_log <- fit$evaluation_log
    if (is.null(evaluation_log)) {
      evaluation_log <- attr(fit, "evaluation_log")
    }
    
    evaluation_log <- as.data.frame(evaluation_log)
    rmse_column <- grep("validation.*rmse", names(evaluation_log), value = TRUE)
    if (length(rmse_column) != 1L) {
      stop("未找到 XGBoost 的内层验证 RMSE。", call. = FALSE)
    }
    
    as.numeric(evaluation_log[[rmse_column]])
  }
  
  # ----------------------------------------------------------
  # 6. 只用内层分组验证调参
  # ----------------------------------------------------------
  
  best_parameters <- list()
  inner_cv_rmse <- numeric()
  
  for (model in names(grids)) {
    grid <- grids[[model]]
    records <- vector("list", nrow(grid))
    
    for (candidate in seq_len(nrow(grid))) {
      parameters <- grid[candidate, , drop = FALSE]
      
      if (model == "XGBoost") {
        curves <- lapply(fold_data, xgb_fold_rmse, parameters = parameters)
        average_rmse <- rowMeans(do.call(cbind, curves))
        
        best_so_far <- Inf
        rounds_without_improvement <- 0L
        last_round <- length(average_rmse)
        
        for (round in seq_along(average_rmse)) {
          if (average_rmse[round] < best_so_far) {
            best_so_far <- average_rmse[round]
            rounds_without_improvement <- 0L
          } else {
            rounds_without_improvement <- rounds_without_improvement + 1L
          }
          
          if (rounds_without_improvement >= 20L) {
            last_round <- round
            break
          }
        }
        
        best_round <- which.min(average_rmse[seq_len(last_round)])
        records[[candidate]] <- cbind(
          parameters,
          best_nrounds = best_round,
          CV_RMSE = average_rmse[best_round],
          CV_R2 = NA_real_
        )
      } else {
        fold_metrics <- sapply(fold_data, function(fold) {
          fit <- fit_regular(model, fold$x, fold$y, parameters)
          predicted <- predict_regular(model, fit, fold$validation_x, parameters)
          score(fold$validation_y, predicted)
        })
        
        records[[candidate]] <- cbind(
          parameters,
          CV_RMSE = mean(fold_metrics["RMSE", ], na.rm = TRUE),
          CV_R2 = mean(fold_metrics["R2", ], na.rm = TRUE)
        )
      }
      
      if (candidate == 1L || candidate %% 12L == 0L || candidate == nrow(grid)) {
        cat(model, "候选参数", candidate, "/", nrow(grid), "\n")
        flush.console()
      }
    }
    
    tuning_table <- do.call(rbind, records)
    
    # 沿用原标准。Ridge 和 Elastic Net 按平均折 R² 调参，
    # 随机森林和 XGBoost 按平均折 RMSE 调参。
    if (model %in% c("Ridge", "ElasticNet")) {
      if (!any(is.finite(tuning_table$CV_R2))) {
        stop("模型 ", model, " 的内层 R² 均不可计算。", call. = FALSE)
      }
      best_row <- which.max(tuning_table$CV_R2)
    } else {
      best_row <- which.min(tuning_table$CV_RMSE)
    }
    
    best_parameters[[model]] <- tuning_table[best_row, , drop = FALSE]
    inner_cv_rmse[model] <- tuning_table$CV_RMSE[best_row]
    
    tuning_table$Selected <- seq_len(nrow(tuning_table)) == best_row
    save_csv(tuning_table, paste0(model, "_inner_tuning.csv"))
  }
  
  selected_model <- names(which.min(inner_cv_rmse))
  
  save_csv(
    data.frame(
      Model = names(inner_cv_rmse),
      Internal_CV_RMSE = unname(inner_cv_rmse),
      Selected = names(inner_cv_rmse) == selected_model
    ),
    "algorithm_selection.csv"
  )
  
  # ----------------------------------------------------------
  # 7. 在完整训练集拟合模型，再预测外层测试集
  # ----------------------------------------------------------
  
  models <- list()
  predictions <- list()
  
  for (model in names(grids)) {
    parameters <- best_parameters[[model]]
    
    if (model == "XGBoost") {
      models[[model]] <- xgboost::xgb.train(
        params = xgb_parameters(parameters),
        data = xgboost::xgb.DMatrix(as.matrix(x_train), label = y_train),
        nrounds = as.integer(parameters$best_nrounds),
        verbose = 0
      )
      
      predictions[[model]] <- as.numeric(predict(
        models[[model]],
        xgboost::xgb.DMatrix(as.matrix(x_test))
      ))
    } else {
      models[[model]] <- fit_regular(
        model, x_train, y_train, parameters,
        importance = if (model == "RandomForest") "permutation" else "none"
      )
      
      predictions[[model]] <- predict_regular(
        model, models[[model]], x_test, parameters
      )
    }
  }
  
  # 简化 Meta 回归。两个线性协变量为平均年龄和男性比例。
  meta_features <- c("Mean age", "Male proportion")
  
  meta_train_x <- as.matrix(
    clean_data(train, outer_recipe)[, meta_features, drop = FALSE]
  )
  meta_test_x <- as.matrix(
    clean_data(test, outer_recipe)[, meta_features, drop = FALSE]
  )
  
  meta_train_data <- data.frame(
    y = y_train,
    vi = train$SE^2,
    Study_ID = train$Study_ID
  )
  
  models$MetaRegression <- metafor::rma.mv(
    yi = y,
    V = vi,
    mods = meta_train_x,
    random = ~ 1 | Study_ID,
    data = meta_train_data,
    method = "REML"
  )
  
  predictions$MetaRegression <- as.numeric(
    predict(models$MetaRegression, newmods = meta_test_x)$pred
  )
  predictions$Mean <- rep(mean(y_train), nrow(test))
  predictions$SelectedML <- predictions[[selected_model]]
  
  # ----------------------------------------------------------
  # 8. 导出性能、逐条预测及可复核材料
  # ----------------------------------------------------------
  
  report_counts <- table(test$Study_ID)
  weights_by_report <- 1 / as.numeric(report_counts[test$Study_ID])
  
  performance <- do.call(rbind, lapply(names(predictions), function(model) {
    do.call(rbind, lapply(c("effect_equal", "report_equal"), function(weighting) {
      weights <- if (weighting == "effect_equal") {
        rep(1, nrow(test))
      } else {
        weights_by_report
      }
      
      metrics <- score(test$ln_aOR, predictions[[model]], weights)
      data.frame(
        Model = model,
        Weighting = weighting,
        RMSE = unname(metrics["RMSE"]),
        MAE = unname(metrics["MAE"]),
        R2 = unname(metrics["R2"]),
        Test_reports = length(unique(test$Study_ID)),
        Test_groups = length(unique(test_group))
      )
    }))
  }))
  
  prediction_table <- do.call(rbind, lapply(names(predictions), function(model) {
    data.frame(
      Effect_ID = test$Effect_ID,
      Study_ID = test$Study_ID,
      Validation_group = test_group,
      Model = model,
      Observed = test$ln_aOR,
      Predicted = predictions[[model]]
    )
  }))
  
  save_csv(performance, "grouped_performance.csv")
  save_csv(prediction_table, "grouped_predictions.csv")
  save_csv(data, "analysis_input_used.csv")
  
  save_csv(
    data.frame(
      Setting = c(
        "split_seed", "cv_seed", "model_seed", "inner_folds",
        "conservative_CEAS", "training_groups", "test_groups",
        "train_effects", "test_effects"
      ),
      Value = c(
        split_seed, cv_seed, model_seed, inner_k,
        conservative_CEAS, length(unique(train_group)), length(unique(test_group)),
        nrow(train), nrow(test)
      )
    ),
    "run_settings.csv"
  )
  
  xgboost::xgb.save(
    models$XGBoost,
    file.path(output_dir, "xgboost_model.json")
  )
  models$XGBoost <- NULL
  
  saveRDS(
    list(
      models = models,
      xgboost_file = "xgboost_model.json",
      recipe = outer_recipe,
      inner_recipes = lapply(fold_data, `[[`, "recipe"),
      parameters = best_parameters,
      selected_by_internal_CV = selected_model,
      features = features,
      train_effect_ids = train$Effect_ID,
      test_effect_ids = test$Effect_ID
    ),
    file.path(output_dir, "models_and_recipes.rds")
  )
  
  capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
  
  writeLines(
    c(
      "Single grouped holdout with 10-fold grouped inner cross-validation.",
      "Preprocessing was fitted separately within each training fold.",
      "Models were selected by training-set cross-validation.",
      "CEAS reports were conservatively grouped for potential overlap."
    ),
    file.path(output_dir, "analysis_notes.txt")
  )
  
  writeLines("COMPLETED", status_path)
  completed <- TRUE
  
  cat("\n按训练集内部 CV_RMSE 选择的模型：", selected_model, "\n")
  print(performance, row.names = FALSE)
  cat("\n结果保存位置：", output_dir, "\n")
  
  invisible(list(
    output_dir = output_dir,
    selected_model = selected_model,
    performance = performance
  ))
}

# 执行完整训练与验证
supplement_result <- run_grouped_holdout()
