############################################################
# Figure S1 replot
# Compact layout with no label/forest overlap
############################################################

required_packages <- c("readr", "dplyr", "svglite")

missing_packages <- required_packages[
  !(required_packages %in% rownames(installed.packages()))
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(svglite)
})

############################################################
# 1. Paths
############################################################

output_dir <- "D:/环境/meta数据csv/EDC_Overall_AnalysisA_SVG"

plot_data_path <- file.path(
  output_dir,
  "10_FigureS1_plot_data.csv"
)

output_svg <- file.path(
  output_dir,
  "FigureS1_Exposure_timing_asthma_outcomes_compact.svg"
)

if (!file.exists(plot_data_path)) {
  stop("Cannot find plot data: ", plot_data_path)
}

############################################################
# 2. Read plot-data file
############################################################

plot_data <- readr::read_csv(
  plot_data_path,
  show_col_types = FALSE,
  progress = FALSE
) %>%
  mutate(
    Label = as.character(Label),
    Row_type = as.character(Row_type),
    Time_group = as.character(Time_group),
    OR = as.numeric(OR),
    CI_lower = as.numeric(CI_lower),
    CI_upper = as.numeric(CI_upper),
    P_value = as.numeric(P_value),
    Studies = as.numeric(Studies),
    Samples = as.numeric(Samples),
    Effect_sizes = as.numeric(Effect_sizes)
  ) %>%
  filter(
    Status == "OK",
    is.finite(OR),
    is.finite(CI_lower),
    is.finite(CI_upper),
    OR > 0,
    CI_lower > 0,
    CI_upper > 0
  )

if (nrow(plot_data) == 0) {
  stop("No valid rows available for Figure S1.")
}

############################################################
# 3. Formatting functions
############################################################

format_p_value <- function(p) {
  if (!is.finite(p)) return("—")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

get_row_color <- function(row_type, time_group) {
  if (row_type == "Overall") return("#B2182B")
  if (time_group == "Childhood") return("#2166AC")
  if (time_group == "Prenatal") return("#E66101")
  "#4D4D4D"
}

############################################################
# 4. Draw compact SVG
############################################################

draw_figureS1_svg <- function(plot_data, output_file) {
  
  n_rows <- nrow(plot_data)
  y_positions <- rev(seq_len(n_rows))
  
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
  
  # 原有列位置和间距
  label_x <- 0.010
  graph_left_x <- 0.315
  graph_right_x <- 0.500
  result_x <- 0.525
  count_x <- 0.750
  p_x <- 0.910
  
  map_or_to_x <- function(or_value) {
    graph_left_x +
      (log(or_value) - log_limits[1]) /
      (log_limits[2] - log_limits[1]) *
      (graph_right_x - graph_left_x)
  }
  
  svglite::svglite(
    file = output_file,
    width = 12.4,
    height = figure_height,
    bg = "white"
  )
  
  on.exit(grDevices::dev.off(), add = TRUE)
  
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
  
  # 恢复原标题
  text(
    x = 0.5,
    y = title_y,
    labels = "Exposure timing and asthma outcome definitions",
    font = 2,
    cex = 1.20,
    col = "#202020"
  )
  
  text(
    x = label_x,
    y = header_y,
    labels = "Exposure timing / outcome",
    pos = 4,
    font = 2,
    cex = 0.94
  )
  
  text(
    x = mean(c(graph_left_x, graph_right_x)),
    y = header_y,
    labels = "Adjusted odds ratio",
    font = 2,
    cex = 0.94
  )
  
  text(
    x = result_x,
    y = header_y,
    labels = "aOR [95% CI]",
    pos = 4,
    font = 2,
    cex = 0.94
  )
  
  text(
    x = count_x,
    y = header_y,
    labels = "Studies / samples / ES",
    font = 2,
    cex = 0.84
  )
  
  text(
    x = p_x,
    y = header_y,
    labels = "P",
    pos = 2,
    font = 2,
    cex = 0.94
  )
  
  segments(
    x0 = 0,
    y0 = header_line_y,
    x1 = 0.95,
    y1 = header_line_y,
    col = "#303030",
    lwd = 1.05
  )
  
  # 交替灰色背景
  for (i in seq_len(n_rows)) {
    if (i %% 2 == 0) {
      rect(
        xleft = 0,
        ybottom = y_positions[i] - 0.43,
        xright = 0.95,
        ytop = y_positions[i] + 0.43,
        col = "#F3F5F6",
        border = NA
      )
    }
  }
  
  # OR = 1 参考线
  segments(
    x0 = map_or_to_x(1),
    y0 = axis_y,
    x1 = map_or_to_x(1),
    y1 = header_line_y,
    col = "#C51B2D",
    lty = 2,
    lwd = 1.05
  )
  
  # 每一行结果
  for (i in seq_len(n_rows)) {
    
    current_color <- get_row_color(
      plot_data$Row_type[i],
      plot_data$Time_group[i]
    )
    
    segments(
      x0 = map_or_to_x(plot_data$CI_lower[i]),
      y0 = y_positions[i],
      x1 = map_or_to_x(plot_data$CI_upper[i]),
      y1 = y_positions[i],
      col = current_color,
      lwd = 1.45
    )
    
    polygon(
      x = c(
        map_or_to_x(plot_data$CI_lower[i]),
        map_or_to_x(plot_data$OR[i]),
        map_or_to_x(plot_data$CI_upper[i]),
        map_or_to_x(plot_data$OR[i])
      ),
      y = c(
        y_positions[i],
        y_positions[i] + 0.20,
        y_positions[i],
        y_positions[i] - 0.20
      ),
      col = current_color,
      border = current_color
    )
    
    text(
      x = label_x,
      y = y_positions[i],
      labels = plot_data$Label[i],
      pos = 4,
      font = 2,
      cex = 0.96,
      col = "#202020"
    )
    
    text(
      x = result_x,
      y = y_positions[i],
      labels = sprintf(
        "%.2f [%.2f, %.2f]",
        plot_data$OR[i],
        plot_data$CI_lower[i],
        plot_data$CI_upper[i]
      ),
      pos = 4,
      cex = 0.86,
      col = "#202020"
    )
    
    count_label <- if (plot_data$Row_type[i] == "Overall") {
      paste0(
        "n=", plot_data$Studies[i],
        "; s=", plot_data$Samples[i],
        "; k=", plot_data$Effect_sizes[i]
      )
    } else {
      paste0(
        "n=", plot_data$Studies[i],
        "; k=", plot_data$Effect_sizes[i]
      )
    }
    
    text(
      x = count_x,
      y = y_positions[i],
      labels = count_label,
      cex = 0.80,
      col = "#202020"
    )
    
    text(
      x = p_x,
      y = y_positions[i],
      labels = format_p_value(plot_data$P_value[i]),
      pos = 2,
      cex = 0.84,
      font = ifelse(
        is.finite(plot_data$P_value[i]) &&
          plot_data$P_value[i] < 0.05,
        2,
        1
      ),
      col = ifelse(
        is.finite(plot_data$P_value[i]) &&
          plot_data$P_value[i] < 0.05,
        "#B2182B",
        "#303030"
      )
    )
  }
  
  # 横轴
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
    x0 = graph_left_x,
    y0 = axis_y,
    x1 = graph_right_x,
    y1 = axis_y,
    lwd = 0.90
  )
  
  segments(
    x0 = tick_x,
    y0 = axis_y,
    x1 = tick_x,
    y1 = axis_y - 0.10,
    lwd = 0.80
  )
  
  text(
    x = tick_x,
    y = axis_y - 0.33,
    labels = format(axis_ticks, trim = TRUE, scientific = FALSE),
    cex = 0.72
  )
  
  # 删除 (log scale)
  text(
    x = mean(c(graph_left_x, graph_right_x)),
    y = axis_y - 0.75,
    labels = "Adjusted odds ratio",
    cex = 0.80
  )
}

############################################################
# 5. Generate SVG
############################################################

draw_figureS1_svg(
  plot_data = plot_data,
  output_file = output_svg
)

cat(
  "\nCompact SVG saved to:\n",
  normalizePath(output_svg, winslash = "/"),
  "\n",
  sep = ""
)

