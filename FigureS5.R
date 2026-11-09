############################################################
# Metals and asthma associations
# Metals类别总效应：纳入所有符合Analysis A条件的Metals效应量
# 具体化合物：k < 2不单独显示，但仍保留在Metals类别总效应中
# n >= 3：研究聚类随机效应 + CR2
# 输出SVG；无底部备注
############################################################

suppressPackageStartupMessages({
  library(metafor)
  library(clubSandwich)
  library(svglite)
})

# 1. 路径 ---------------------------------------------------

file_path <- "D:/环境/meta数据csv/三水平_meta_AnalysisA_最终版.csv"
out_path <- "D:/环境/meta数据csv/EDC_Overall_AnalysisA_SVG/三水平出图"

if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

output_svg <- file.path(
  out_path,
  "FigureS4_Metals_compound_estimates.svg"
)

# 2. 读取及筛选数据 -----------------------------------------

raw_dat <- read.csv(
  file_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)

is_yes <- function(x) {
  trimws(as.character(x)) %in%
    c("1", "TRUE", "True", "true", "YES", "Yes", "yes", "Y", "y")
}

find_col <- function(dat, candidates) {
  found <- candidates[candidates %in% colnames(dat)]
  
  if (length(found) == 0) {
    stop(
      paste0(
        "找不到变量列：",
        paste(candidates, collapse = " / ")
      )
    )
  }
  
  found[1]
}

class_col <- find_col(
  raw_dat,
  c("EDC_class_Clean", "EDC_class", "EDC_Class")
)

compound_col <- find_col(
  raw_dat,
  c("Compound_Clean", "Compound", "Chemical", "Exposure")
)

# 先筛选全部符合Analysis A条件的Metals记录
class_dat <- raw_dat[
  is_yes(raw_dat$Analysis_A) &
    is_yes(raw_dat$Analysis_Ready) &
    is_yes(raw_dat$Outcome_Eligible) &
    is_yes(raw_dat$Main_Effect) &
    trimws(as.character(raw_dat[[class_col]])) == "Metals",
  ,
  drop = FALSE
]

class_dat$Study_ID <- as.character(class_dat[["Study-ID"]])
class_dat$Sample_ID <- as.character(class_dat[["Sample_ID"]])
class_dat$Compound <- trimws(as.character(class_dat[[compound_col]]))
class_dat$OR <- suppressWarnings(as.numeric(class_dat$OR))
class_dat$yi <- suppressWarnings(as.numeric(class_dat$ln_aOR))
class_dat$SE <- suppressWarnings(as.numeric(class_dat$SE))
class_dat$vi <- class_dat$SE^2

# 去除无效记录
class_dat <- class_dat[
  !is.na(class_dat$Study_ID) &
    class_dat$Study_ID != "" &
    !is.na(class_dat$Sample_ID) &
    class_dat$Sample_ID != "" &
    !is.na(class_dat$Compound) &
    class_dat$Compound != "" &
    is.finite(class_dat$OR) &
    class_dat$OR > 0 &
    is.finite(class_dat$yi) &
    is.finite(class_dat$vi) &
    class_dat$vi > 0,
  ,
  drop = FALSE
]

if (nrow(class_dat) == 0) {
  stop("没有符合条件的Metals数据。")
}

# 只用于具体化合物结果的展示筛选。
# 不能用此数据集计算Metals类别总效应。
compound_k <- table(class_dat$Compound)
eligible_compounds <- names(compound_k)[compound_k >= 2]

compound_dat <- class_dat[
  class_dat$Compound %in% eligible_compounds,
  ,
  drop = FALSE
]

# 3. 模型函数 -----------------------------------------------

format_p <- function(p) {
  if (!is.finite(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

# 三水平随机效应模型 + 以Study_ID为聚类单位的CR2稳健推断
fit_cr2 <- function(d) {
  
  n_studies <- length(unique(d$Study_ID))
  
  if (n_studies < 3 || nrow(d) < 2) {
    return(NULL)
  }
  
  model <- tryCatch(
    rma.mv(
      yi = yi,
      V = vi,
      random = ~ 1 | Study_ID / Sample_ID,
      method = "REML",
      data = d
    ),
    error = function(e) NULL
  )
  
  if (is.null(model)) return(NULL)
  
  cr2 <- tryCatch(
    coef_test(
      model,
      vcov = "CR2",
      cluster = d$Study_ID,
      test = "Satterthwaite"
    ),
    error = function(e) NULL
  )
  
  if (is.null(cr2)) return(NULL)
  
  beta <- as.numeric(cr2[1, "beta"])
  se <- as.numeric(cr2[1, "SE"])
  df <- as.numeric(cr2[1, "df_Satt"])
  
  if (!is.finite(beta) || !is.finite(se) || !is.finite(df) || df <= 0) {
    return(NULL)
  }
  
  list(
    OR = exp(beta),
    CI_lower = exp(beta - qt(0.975, df) * se),
    CI_upper = exp(beta + qt(0.975, df) * se),
    P_value = as.numeric(cr2[1, "p_Satt"])
  )
}

# 具体化合物研究数不足3时，仅提供探索性95%CI，不报告P值
fit_exploratory <- function(d) {
  
  model <- tryCatch(
    rma.uni(
      yi = yi,
      vi = vi,
      method = "REML",
      data = d
    ),
    error = function(e) NULL
  )
  
  if (is.null(model)) return(NULL)
  
  beta <- as.numeric(model$b[1])
  se <- sqrt(as.numeric(vcov(model)[1, 1]))
  
  if (!is.finite(beta) || !is.finite(se)) {
    return(NULL)
  }
  
  list(
    OR = exp(beta),
    CI_lower = exp(beta - 1.96 * se),
    CI_upper = exp(beta + 1.96 * se),
    P_value = NA_real_
  )
}

# 4. Metals类别及具体化合物结果 -----------------------------

# 类别总效应使用全部合格Metals记录，包括仅有一个效应量的化合物
class_result <- fit_cr2(class_dat)

if (is.null(class_result)) {
  stop("Metals总体模型无法拟合。")
}

plot_data <- data.frame(
  Label = "Metals (class)",
  Studies = length(unique(class_dat$Study_ID)),
  Effect_sizes = nrow(class_dat),
  OR = class_result$OR,
  CI_lower = class_result$CI_lower,
  CI_upper = class_result$CI_upper,
  P_value = class_result$P_value,
  Row_type = "Class",
  stringsAsFactors = FALSE
)

# 具体化合物结果只使用k >= 2的化合物
for (compound_name in sort(unique(compound_dat$Compound))) {
  
  d_sub <- compound_dat[
    compound_dat$Compound == compound_name,
    ,
    drop = FALSE
  ]
  
  n_studies <- length(unique(d_sub$Study_ID))
  
  if (n_studies >= 3) {
    result <- fit_cr2(d_sub)
    row_type <- "Pooled"
  } else {
    result <- fit_exploratory(d_sub)
    row_type <- "Exploratory"
  }
  
  if (!is.null(result)) {
    plot_data <- rbind(
      plot_data,
      data.frame(
        Label = compound_name,
        Studies = n_studies,
        Effect_sizes = nrow(d_sub),
        OR = result$OR,
        CI_lower = result$CI_lower,
        CI_upper = result$CI_upper,
        P_value = result$P_value,
        Row_type = row_type,
        stringsAsFactors = FALSE
      )
    )
  }
}

plot_data$Y <- rev(seq_len(nrow(plot_data)))

# 5. 绘图位置与尺度 -----------------------------------------

n_rows <- nrow(plot_data)
figure_height <- max(5.4, 2.9 + n_rows * 0.52)

all_values <- c(plot_data$CI_lower, plot_data$CI_upper, 1)
log_values <- log(all_values[is.finite(all_values) & all_values > 0])

log_min <- min(log_values)
log_max <- max(log_values)
log_range <- log_max - log_min

if (!is.finite(log_range) || log_range <= 0) {
  log_range <- 1
}

log_limits <- c(
  log_min - 0.15 * log_range,
  log_max + 0.15 * log_range
)

label_x <- 0.010
graph_left_x <- 0.285
graph_right_x <- 0.455
result_x <- 0.475
count_x <- 0.665
p_x <- 0.805

map_or_to_x <- function(x) {
  graph_left_x +
    (log(x) - log_limits[1]) /
    (log_limits[2] - log_limits[1]) *
    (graph_right_x - graph_left_x)
}

# 6. 绘制SVG ------------------------------------------------

svglite::svglite(
  file = output_svg,
  width = 12.4,
  height = figure_height,
  bg = "white"
)

par(
  mar = c(1.1, 0.15, 0.2, 0.15),
  family = "sans",
  xaxs = "i",
  yaxs = "i"
)

title_y <- n_rows + 2.05
header_y <- n_rows + 0.90
header_line_y <- n_rows + 0.48
axis_y <- -0.05

plot.new()

plot.window(
  xlim = c(0, 1),
  ylim = c(-1.55, n_rows + 2.35)
)

text(
  0.5, title_y, "Metals",
  font = 2, cex = 1.20, col = "#202020"
)

text(
  label_x, header_y, "Compound",
  pos = 4, font = 2, cex = 0.94
)

text(
  mean(c(graph_left_x, graph_right_x)),
  header_y,
  "Adjusted odds ratio",
  font = 2, cex = 0.94
)

text(
  result_x, header_y, "aOR [95% CI]",
  pos = 4, font = 2, cex = 0.94
)

text(
  count_x, header_y, "Studies / ES",
  font = 2, cex = 0.84
)

text(
  p_x, header_y, "P",
  pos = 2, font = 2, cex = 0.94
)

segments(
  0, header_line_y,
  0.84, header_line_y,
  col = "#303030", lwd = 1.05
)

for (i in seq_len(n_rows)) {
  if (i %% 2 == 0) {
    rect(
      0, plot_data$Y[i] - 0.43,
      0.84, plot_data$Y[i] + 0.43,
      col = "#F3F5F6",
      border = NA
    )
  }
}

segments(
  map_or_to_x(1), axis_y,
  map_or_to_x(1), header_line_y,
  col = "#C51B2D",
  lty = 2,
  lwd = 1.05
)

for (i in seq_len(n_rows)) {
  
  is_class <- plot_data$Row_type[i] == "Class"
  is_exploratory <- plot_data$Row_type[i] == "Exploratory"
  
  line_colour <- if (is_class) {
    "#1B9E77"
  } else if (is_exploratory) {
    "#8A96A3"
  } else {
    "#5B8CCB"
  }
  
  segments(
    map_or_to_x(plot_data$CI_lower[i]),
    plot_data$Y[i],
    map_or_to_x(plot_data$CI_upper[i]),
    plot_data$Y[i],
    col = line_colour,
    lwd = ifelse(is_class, 1.45, 1.10)
  )
  
  if (is_class) {
    
    polygon(
      x = c(
        map_or_to_x(plot_data$CI_lower[i]),
        map_or_to_x(plot_data$OR[i]),
        map_or_to_x(plot_data$CI_upper[i]),
        map_or_to_x(plot_data$OR[i])
      ),
      y = c(
        plot_data$Y[i],
        plot_data$Y[i] + 0.20,
        plot_data$Y[i],
        plot_data$Y[i] - 0.20
      ),
      col = "#1B9E77",
      border = "#1B9E77"
    )
    
  } else {
    
    points(
      map_or_to_x(plot_data$OR[i]),
      plot_data$Y[i],
      pch = 22,
      cex = 1.10,
      bg = ifelse(is_exploratory, "#B5C0CA", "#6A9BD4"),
      col = "#303030",
      lwd = 1
    )
  }
  
  text(
    label_x, plot_data$Y[i],
    plot_data$Label[i],
    pos = 4,
    font = ifelse(is_class, 2, 1),
    cex = 0.96,
    col = "#202020"
  )
  
  text(
    result_x, plot_data$Y[i],
    sprintf(
      "%.2f [%.2f, %.2f]",
      plot_data$OR[i],
      plot_data$CI_lower[i],
      plot_data$CI_upper[i]
    ),
    pos = 4,
    font = ifelse(is_class, 2, 1),
    cex = 0.86,
    col = "#202020"
  )
  
  text(
    count_x, plot_data$Y[i],
    paste0(
      "n=", plot_data$Studies[i],
      "; k=", plot_data$Effect_sizes[i]
    ),
    cex = 0.80,
    font = ifelse(is_class, 2, 1),
    col = "#202020"
  )
  
  significant <- is.finite(plot_data$P_value[i]) &&
    plot_data$P_value[i] < 0.05
  
  text(
    p_x, plot_data$Y[i],
    format_p(plot_data$P_value[i]),
    pos = 2,
    cex = 0.84,
    font = ifelse(significant, 2, 1),
    col = ifelse(significant, "#B2182B", "#303030")
  )
}

candidate_ticks <- c(0.10, 0.20, 0.50, 1, 2, 5, 10, 20)

axis_ticks <- candidate_ticks[
  log(candidate_ticks) >= log_limits[1] &
    log(candidate_ticks) <= log_limits[2]
]

if (!1 %in% axis_ticks) {
  axis_ticks <- sort(unique(c(axis_ticks, 1)))
}

tick_x <- map_or_to_x(axis_ticks)

segments(
  graph_left_x, axis_y,
  graph_right_x, axis_y,
  lwd = 0.90
)

segments(
  tick_x, axis_y,
  tick_x, axis_y - 0.10,
  lwd = 0.80
)

text(
  tick_x,
  axis_y - 0.33,
  format(axis_ticks, trim = TRUE, scientific = FALSE),
  cex = 0.72
)

text(
  mean(c(graph_left_x, graph_right_x)),
  axis_y - 0.75,
  "Adjusted odds ratio",
  cex = 0.80
)

dev.off()

cat(
  "\n完成：Metals SVG已保存至：\n",
  normalizePath(output_svg, winslash = "/"),
  "\n",
  sep = ""
)

