make_matrix_plot <- function(matrix_dat, matrix_label, output_file) {
  
  matrix_dat <- matrix_dat[
    is.finite(matrix_dat$yi) &
      is.finite(matrix_dat$vi) &
      matrix_dat$vi > 0,
    ,
    drop = FALSE
  ]
  
  if (nrow(matrix_dat) == 0) {
    stop(paste0(matrix_label, "没有可用于绘图的数据。"))
  }
  
  class_result <- fit_three_level_cr2(matrix_dat)
  
  if (is.null(class_result)) {
    stop(paste0(matrix_label, "总体三水平模型无法拟合。"))
  }
  
  class_row <- data.frame(
    Label = paste0(matrix_label, " phthalates (class)"),
    Studies = length(unique(matrix_dat$Study_ID)),
    Effect_sizes = nrow(matrix_dat),
    OR = class_result$OR,
    CI_lower = class_result$CI_lower,
    CI_upper = class_result$CI_upper,
    P_value = class_result$P_value,
    Row_type = "Class",
    stringsAsFactors = FALSE
  )
  
  compounds <- sort(unique(matrix_dat$Compound))
  compound_rows <- list()
  
  for (compound_name in compounds) {
    
    d_sub <- matrix_dat[
      matrix_dat$Compound == compound_name,
      ,
      drop = FALSE
    ]
    
    n_studies <- length(unique(d_sub$Study_ID))
    
    if (n_studies >= 3) {
      result <- fit_three_level_cr2(d_sub)
      result_type <- "Pooled"
    } else {
      result <- fit_small_study_ci(d_sub)
      result_type <- "Small_study"
    }
    
    if (!is.null(result)) {
      compound_rows[[compound_name]] <- data.frame(
        Label = compound_name,
        Studies = n_studies,
        Effect_sizes = nrow(d_sub),
        OR = result$OR,
        CI_lower = result$CI_lower,
        CI_upper = result$CI_upper,
        P_value = result$P_value,
        Row_type = result_type,
        stringsAsFactors = FALSE
      )
    }
  }
  
  plot_data <- do.call(
    rbind,
    c(list(class_row), compound_rows)
  )
  
  plot_data$Y <- rev(seq_len(nrow(plot_data)))
  
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
  
  # 仅压缩列间距；字体和行距完全不变
  label_x <- 0.010
  graph_left_x <- 0.285
  graph_right_x <- 0.455
  result_x <- 0.475
  count_x <- 0.665
  p_x <- 0.805
  
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
  
  text(
    x = 0.5,
    y = title_y,
    labels = ifelse(
      matrix_label == "Dust",
      "Phthalates (Dust)",
      "Phthalates (Urine)"
    ),
    font = 2,
    cex = 1.20,
    col = "#202020"
  )
  
  text(label_x, header_y, "Compound",
       pos = 4, font = 2, cex = 0.94)
  
  text(mean(c(graph_left_x, graph_right_x)), header_y,
       "Adjusted odds ratio", font = 2, cex = 0.94)
  
  text(result_x, header_y, "aOR [95% CI]",
       pos = 4, font = 2, cex = 0.94)
  
  text(count_x, header_y, "Studies / ES",
       font = 2, cex = 0.84)
  
  text(p_x, header_y, "P",
       pos = 2, font = 2, cex = 0.94)
  
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
    is_small <- plot_data$Row_type[i] == "Small_study"
    
    line_colour <- if (is_class) {
      "#E66101"
    } else if (is_small) {
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
        bg = ifelse(is_small, "#B5C0CA", "#6A9BD4"),
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
}