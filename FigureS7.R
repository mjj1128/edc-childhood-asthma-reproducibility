############################################################
# Phenols and asthma associations
# 小菱形、小方块；标题相对表格内容居中
############################################################

suppressPackageStartupMessages({
  library(metafor)
  library(clubSandwich)
  library(svglite)
})

file_path <- "D:/环境/meta数据csv/三水平_meta_AnalysisA_最终版.csv"
out_path <- "D:/环境/meta数据csv/EDC_Overall_AnalysisA_SVG/三水平出图"

if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)

output_svg <- file.path(out_path, "FigureS6_Phenols_compound_estimates.svg")

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
  found <- candidates[candidates %in% names(dat)]
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

matrix_col <- find_col(
  raw_dat,
  c("Matrix_Clean", "Matrix", "Exposure_matrix")
)

dat <- raw_dat[
  is_yes(raw_dat$Analysis_A) &
    is_yes(raw_dat$Analysis_Ready) &
    is_yes(raw_dat$Outcome_Eligible) &
    is_yes(raw_dat$Main_Effect) &
    trimws(as.character(raw_dat[[class_col]])) == "Phenols",
  ,
  drop = FALSE
]

dat$Study_ID <- as.character(dat[["Study-ID"]])
dat$Sample_ID <- as.character(dat$Sample_ID)
dat$Compound <- trimws(as.character(dat[[compound_col]]))
dat$Matrix <- trimws(as.character(dat[[matrix_col]]))
dat$yi <- suppressWarnings(as.numeric(dat$ln_aOR))
dat$vi <- suppressWarnings(as.numeric(dat$SE))^2

dat <- dat[
  !is.na(dat$Study_ID) & dat$Study_ID != "" &
    !is.na(dat$Sample_ID) & dat$Sample_ID != "" &
    !is.na(dat$Compound) & dat$Compound != "" &
    is.finite(dat$yi) & is.finite(dat$vi) & dat$vi > 0,
  ,
  drop = FALSE
]

if (nrow(dat) == 0) stop("Phenols筛选后没有可用数据。")

dat$Exposure_group <- ifelse(
  grepl("urine|urinary", tolower(dat$Matrix)),
  "Urine",
  "Environment"
)

format_p <- function(p) {
  if (!is.finite(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

fit_cr2 <- function(d) {
  
  if (nrow(d) < 2 || length(unique(d$Study_ID)) < 3) {
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
    Lower = exp(beta - qt(0.975, df) * se),
    Upper = exp(beta + qt(0.975, df) * se),
    P = as.numeric(cr2[1, "p_Satt"])
  )
}

fit_exploratory <- function(d) {
  
  model <- tryCatch(
    rma.uni(yi = yi, vi = vi, method = "REML", data = d),
    error = function(e) NULL
  )
  
  if (is.null(model)) return(NULL)
  
  beta <- as.numeric(model$b[1])
  se <- sqrt(as.numeric(vcov(model)[1, 1]))
  
  list(
    OR = exp(beta),
    Lower = exp(beta - 1.96 * se),
    Upper = exp(beta + 1.96 * se),
    P = NA_real_
  )
}

# 顶部类别汇总
class_fit <- fit_cr2(dat)

if (is.null(class_fit)) {
  stop("Phenols总体模型无法拟合。")
}

plot_data <- data.frame(
  Label = "Phenols (class)",
  Studies = length(unique(dat$Study_ID)),
  Effect_sizes = nrow(dat),
  OR = class_fit$OR,
  Lower = class_fit$Lower,
  Upper = class_fit$Upper,
  P = class_fit$P,
  Type = "Class",
  stringsAsFactors = FALSE
)

# 具体物质：同一物质、同一基质中 k >= 2
dat$Display_key <- paste(dat$Exposure_group, dat$Compound, sep = " | ")
key_k <- table(dat$Display_key)
eligible_key <- names(key_k)[key_k >= 2]

for (key in sort(eligible_key)) {
  
  d_sub <- dat[dat$Display_key == key, , drop = FALSE]
  n_studies <- length(unique(d_sub$Study_ID))
  
  if (n_studies >= 3) {
    result <- fit_cr2(d_sub)
    result_type <- "Pooled"
  } else {
    result <- fit_exploratory(d_sub)
    result_type <- "Exploratory"
  }
  
  if (!is.null(result)) {
    
    group_name <- unique(d_sub$Exposure_group)[1]
    compound_name <- unique(d_sub$Compound)[1]
    
    plot_data <- rbind(
      plot_data,
      data.frame(
        Label = paste0(
          compound_name,
          ifelse(group_name == "Urine", " (urine)", " (environment)")
        ),
        Studies = n_studies,
        Effect_sizes = nrow(d_sub),
        OR = result$OR,
        Lower = result$Lower,
        Upper = result$Upper,
        P = result$P,
        Type = result_type,
        stringsAsFactors = FALSE
      )
    )
  }
}

plot_data$Y <- rev(seq_len(nrow(plot_data)))

n_rows <- nrow(plot_data)
figure_height <- max(5.4, 2.9 + n_rows * 0.52)

all_values <- c(plot_data$Lower, plot_data$Upper, 1)
log_values <- log(all_values[is.finite(all_values) & all_values > 0])

log_min <- min(log_values)
log_max <- max(log_values)
log_range <- log_max - log_min

if (!is.finite(log_range) || log_range <= 0) log_range <- 1

log_limits <- c(
  log_min - 0.15 * log_range,
  log_max + 0.15 * log_range
)

# 紧凑列间距；字体、行距不变
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

svglite::svglite(
  output_svg,
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

# 标题相对内容区 0 到 0.84 居中
text(
  0.42,
  title_y,
  "Phenols",
  font = 2,
  cex = 1.20,
  col = "#202020"
)

text(label_x, header_y, "Compound", pos = 4, font = 2, cex = 0.94)

text(
  mean(c(graph_left_x, graph_right_x)),
  header_y,
  "Adjusted odds ratio",
  font = 2,
  cex = 0.94
)

text(result_x, header_y, "aOR [95% CI]", pos = 4, font = 2, cex = 0.94)
text(count_x, header_y, "Studies / ES", font = 2, cex = 0.84)
text(p_x, header_y, "P", pos = 2, font = 2, cex = 0.94)

segments(0, header_line_y, 0.84, header_line_y, col = "#303030", lwd = 1.05)

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
  
  is_class <- plot_data$Type[i] == "Class"
  is_exploratory <- plot_data$Type[i] == "Exploratory"
  
  line_colour <- if (is_class) {
    "#7570B3"
  } else if (is_exploratory) {
    "#8A96A3"
  } else {
    "#5B8CCB"
  }
  
  # 完整95%CI线段
  segments(
    map_or_to_x(plot_data$Lower[i]),
    plot_data$Y[i],
    map_or_to_x(plot_data$Upper[i]),
    plot_data$Y[i],
    col = line_colour,
    lwd = ifelse(is_class, 1.45, 1.10)
  )
  
  # 小菱形：只标记汇总点，不改变CI线
  if (is_class) {
    
    diamond_x <- map_or_to_x(plot_data$OR[i])
    diamond_half_width <- 0.008
    
    polygon(
      x = c(
        diamond_x - diamond_half_width,
        diamond_x,
        diamond_x + diamond_half_width,
        diamond_x
      ),
      y = c(
        plot_data$Y[i],
        plot_data$Y[i] + 0.10,
        plot_data$Y[i],
        plot_data$Y[i] - 0.10
      ),
      col = "#7570B3",
      border = "#7570B3"
    )
    
  } else {
    
    points(
      map_or_to_x(plot_data$OR[i]),
      plot_data$Y[i],
      pch = 22,
      cex = 0.60,
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
      plot_data$Lower[i],
      plot_data$Upper[i]
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
  
  significant <- is.finite(plot_data$P[i]) &&
    plot_data$P[i] < 0.05
  
  text(
    p_x, plot_data$Y[i],
    format_p(plot_data$P[i]),
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

segments(graph_left_x, axis_y, graph_right_x, axis_y, lwd = 0.90)
segments(tick_x, axis_y, tick_x, axis_y - 0.10, lwd = 0.80)

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
  "\n完成：Phenols SVG已保存至：\n",
  normalizePath(output_svg, winslash = "/"),
  "\n",
  sep = ""
)

