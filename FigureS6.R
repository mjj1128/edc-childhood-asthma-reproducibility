############################################################
# PAHs and asthma associations
#
# 顶部 PAHs (class)：
# 使用全部合格 PAH 效应量进行三水平 + CR2 汇总
#
# 具体物质行：
# k < 2 不显示；但不影响顶部 PAHs 类别总体汇总
# n < 3 时仅显示探索性CI，P列为空
#
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
  "FigureS5_PAHs_compound_estimates.svg"
)

# 2. 读取数据 -----------------------------------------------

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
    stop(paste0("找不到变量列：", paste(candidates, collapse = " / ")))
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

# 3. 筛选全部合格 PAH 效应量 -------------------------------

all_pah_dat <- raw_dat[
  is_yes(raw_dat$Analysis_A) &
    is_yes(raw_dat$Analysis_Ready) &
    is_yes(raw_dat$Outcome_Eligible) &
    is_yes(raw_dat$Main_Effect) &
    trimws(as.character(raw_dat[[class_col]])) == "PAHs",
  ,
  drop = FALSE
]

all_pah_dat$Study_ID <- as.character(all_pah_dat[["Study-ID"]])
all_pah_dat$Sample_ID <- as.character(all_pah_dat$Sample_ID)
all_pah_dat$Compound <- trimws(as.character(all_pah_dat[[compound_col]]))

all_pah_dat$yi <- suppressWarnings(as.numeric(all_pah_dat$ln_aOR))
all_pah_dat$vi <- suppressWarnings(as.numeric(all_pah_dat$SE))^2

all_pah_dat <- all_pah_dat[
  !is.na(all_pah_dat$Study_ID) &
    all_pah_dat$Study_ID != "" &
    !is.na(all_pah_dat$Sample_ID) &
    all_pah_dat$Sample_ID != "" &
    !is.na(all_pah_dat$Compound) &
    all_pah_dat$Compound != "" &
    is.finite(all_pah_dat$yi) &
    is.finite(all_pah_dat$vi) &
    all_pah_dat$vi > 0,
  ,
  drop = FALSE
]

if (nrow(all_pah_dat) == 0) {
  stop("PAHs筛选后没有可用于分析的数据。")
}

# 仅用于展示具体物质行：
# k < 2 的物质不显示
compound_k <- table(all_pah_dat$Compound)

display_compounds <- names(compound_k)[compound_k >= 2]

display_pah_dat <- all_pah_dat[
  all_pah_dat$Compound %in% display_compounds,
  ,
  drop = FALSE
]

# 4. 模型函数 -----------------------------------------------

format_p <- function(p) {
  if (!is.finite(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

# 正式模型：
# 优先三水平 Study_ID / Sample_ID + CR2
# 若该层级不可识别，退回研究聚类随机效应 + CR2
fit_cr2 <- function(d) {
  
  n_studies <- length(unique(d$Study_ID))
  
  if (nrow(d) < 2 || n_studies < 3) {
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
  
  if (is.null(model)) {
    model <- tryCatch(
      rma.mv(
        yi = yi,
        V = vi,
        random = ~ 1 | Study_ID,
        method = "REML",
        data = d
      ),
      error = function(e) NULL
    )
  }
  
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
  
  list(
    OR = exp(beta),
    CI_lower = exp(beta - qt(0.975, df) * se),
    CI_upper = exp(beta + qt(0.975, df) * se),
    P_value = as.numeric(cr2[1, "p_Satt"])
  )
}

# n < 3 的具体物质：探索性CI，不显示P
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
  
  list(
    OR = exp(beta),
    CI_lower = exp(beta - 1.96 * se),
    CI_upper = exp(beta + 1.96 * se),
    P_value = NA_real_
  )
}

# 5. 顶部 PAHs 类别结果：使用全部合格 PAH 数据 ------------

class_result <- fit_cr2(all_pah_dat)

if (is.null(class_result)) {
  stop("PAHs总体三水平 + CR2 模型无法拟合。")
}

plot_data <- data.frame(
  Label = "PAHs (class)",
  Studies = length(unique(all_pah_dat$Study_ID)),
  Effect_sizes = nrow(all_pah_dat),
  OR = class_result$OR,
  CI_lower = class_result$CI_lower,
  CI_upper = class_result$CI_upper,
  P_value = class_result$P_value,
  Row_type = "Class",
  stringsAsFactors = FALSE
)

# 6. 具体 PAH 物质：只显示 k >= 2 --------------------------

if (nrow(display_pah_dat) > 0) {
  
  for (compound_name in sort(unique(display_pah_dat$Compound))) {
    
    d_sub <- display_pah_dat[
      display_pah_dat$Compound == compound_name,
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
}

plot_data$Y <- rev(seq_len(nrow(plot_data)))

# 7. 图形参数 -----------------------------------------------

n_rows <- nrow(plot_data)

figure_height <- max(5.4, 2.9 + n_rows * 0.52)

all_values <- c(
  plot_data$CI_lower,
  plot_data$CI_upper,
  1
)

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

# 紧凑列间距；字体及行距不变
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

# 8. 绘制 SVG -----------------------------------------------

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
  0.5, title_y,
  "PAHs",
  font = 2,
  cex = 1.20,
  col = "#202020"
)

text(
  label_x, header_y,
  "Compound",
  pos = 4,
  font = 2,
  cex = 0.94
)

text(
  mean(c(graph_left_x, graph_right_x)),
  header_y,
  "Adjusted odds ratio",
  font = 2,
  cex = 0.94
)

text(
  result_x, header_y,
  "aOR [95% CI]",
  pos = 4,
  font = 2,
  cex = 0.94
)

text(
  count_x, header_y,
  "Studies / ES",
  font = 2,
  cex = 0.84
)

text(
  p_x, header_y,
  "P",
  pos = 2,
  font = 2,
  cex = 0.94
)

segments(
  0, header_line_y,
  0.84, header_line_y,
  col = "#303030",
  lwd = 1.05
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
    "#E66101"
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
      col = "#E66101",
      border = "#E66101"
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
  "\n完成：PAHs SVG已保存至：\n",
  normalizePath(output_svg, winslash = "/"),
  "\n",
  sep = ""
)

