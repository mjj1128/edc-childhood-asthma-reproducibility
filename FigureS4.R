# Recreate Figure S4 from archived CR2 results and analysis inputs.
# Run from the repository root.

if (.Platform$OS.type == "windows") {
  invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
}

results <- read.csv(
  file.path("reference", "meta", "bisphenols_forest_input.csv"),
  fileEncoding = "UTF-8"
)

effects <- read.csv(
  file.path("data", "meta_effect_estimates.csv"),
  fileEncoding = "UTF-8",
  check.names = FALSE
)

results <- results[order(results$Y, decreasing = TRUE), ]

# The archived results table leaves BPF unpooled. Figure S4 shows
# a descriptive inverse-variance summary using the original BPF inputs.
bisphenols <- subset(effects, EDC_class_Clean == "Bisphenols")
bpf <- subset(bisphenols, Compound_Clean == "BPF")
bpf_row <- results$Compound == "BPF"

stopifnot(
  identical(results$Compound, c("Bisphenols (class)", "BPA", "BPS", "BPF")),
  nrow(bisphenols) == results$K_effect_sizes[results$Is_class],
  length(unique(bisphenols$Study_ID)) == results$N_studies[results$Is_class],
  nrow(bpf) == results$K_effect_sizes[bpf_row],
  length(unique(bpf$Study_ID)) == results$N_studies[bpf_row],
  all(is.finite(bpf$yi)),
  all(is.finite(bpf$vi) & bpf$vi > 0)
)

weights <- 1 / bpf$vi
log_or <- weighted.mean(bpf$yi, weights)
log_se <- sqrt(1 / sum(weights))

results[bpf_row, c("OR", "CI_low", "CI_high")] <-
  as.list(exp(log_or + qnorm(c(0.5, 0.025, 0.975)) * log_se))

results$Estimate_text[bpf_row] <- with(
  results[bpf_row, ],
  sprintf("%.2f [%.2f, %.2f]", OR, CI_low, CI_high)
)
results$P_text[bpf_row] <- ""

draw_bisphenol_forest <- function(results, output_svg) {
  svg(output_svg, width = 12.4, height = 5.4, bg = "white", family = "sans")
  on.exit(dev.off(), add = TRUE)
  
  plot_bounds <- c(0.315, 0.500)
  log_limits <- range(log(c(results$CI_low, results$CI_high, 1)))
  log_limits <- log_limits + c(-0.15, 0.15) * diff(log_limits)
  
  map_or <- function(x) {
    plot_bounds[1] +
      diff(plot_bounds) * (log(x) - log_limits[1]) / diff(log_limits)
  }
  
  label_x <- 0.010
  estimate_x <- 0.525
  count_x <- 0.750
  p_x <- 0.910
  plot_right <- 0.950
  header_y <- nrow(results) + 0.90
  header_line_y <- nrow(results) + 0.48
  axis_y <- -0.05
  
  par(mar = c(1.1, 0.15, 0.2, 0.15), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(
    xlim = c(0, 1),
    ylim = c(-1.55, nrow(results) + 2.35)
  )
  
  text(
    0.5, nrow(results) + 2.05, "Bisphenols",
    font = 2, cex = 1.20, col = "#202020"
  )
  text(label_x, header_y, "Compound", pos = 4, font = 2, cex = 0.94)
  text(mean(plot_bounds), header_y, "Adjusted odds ratio", font = 2, cex = 0.94)
  text(estimate_x, header_y, "aOR [95% CI]", pos = 4, font = 2, cex = 0.94)
  text(count_x, header_y, "Studies / ES", font = 2, cex = 0.84)
  text(p_x, header_y, "P", pos = 2, font = 2, cex = 0.94)
  
  segments(
    0, header_line_y, plot_right, header_line_y,
    col = "#303030", lwd = 1.05
  )
  
  shaded_y <- results$Y[seq(2, nrow(results), by = 2)]
  rect(
    0, shaded_y - 0.43, plot_right, shaded_y + 0.43,
    col = "#F3F5F6", border = NA
  )
  
  segments(
    map_or(1), axis_y, map_or(1), header_line_y,
    col = "#C51B2D", lty = 2, lwd = 1.05
  )
  
  for (i in seq_len(nrow(results))) {
    row <- results[i, ]
    descriptive <- row$Compound == "BPF"
    ci_x <- map_or(c(row$CI_low, row$CI_high))
    or_x <- map_or(row$OR)
    
    line_color <- if (row$Is_class) {
      "#7B3294"
    } else if (descriptive) {
      "#8A96A3"
    } else {
      "#5B8CCB"
    }
    row_font <- if (row$Is_class) 2 else 1
    
    segments(
      ci_x[1], row$Y, ci_x[2], row$Y,
      col = line_color,
      lwd = if (row$Is_class) 1.45 else 1.10
    )
    
    if (row$Is_class) {
      polygon(
        c(ci_x[1], or_x, ci_x[2], or_x),
        row$Y + c(0, 0.20, 0, -0.20),
        col = line_color,
        border = line_color
      )
    } else {
      points(
        or_x, row$Y,
        pch = 22, cex = 0.60,
        bg = if (descriptive) "#B5C0CA" else "#6A9BD4",
        col = "#303030"
      )
    }
    
    text(
      label_x, row$Y, row$Compound,
      pos = 4, font = row_font, cex = 0.96, col = "#202020"
    )
    text(
      estimate_x, row$Y, row$Estimate_text,
      pos = 4, font = row_font, cex = 0.86, col = "#202020"
    )
    text(
      count_x, row$Y, row$Count_text,
      font = row_font, cex = 0.80, col = "#202020"
    )
    
    significant <- !is.na(row$P_value) && row$P_value < 0.05
    text(
      p_x, row$Y, row$P_text,
      pos = 2, cex = 0.84,
      font = if (significant) 2 else 1,
      col = if (significant) "#B2182B" else "#303030"
    )
  }
  
  ticks <- c(1, 2)
  tick_x <- map_or(ticks)
  
  segments(plot_bounds[1], axis_y, plot_bounds[2], axis_y, lwd = 0.90)
  segments(tick_x, axis_y, tick_x, axis_y - 0.10, lwd = 0.80)
  text(tick_x, axis_y - 0.33, ticks, cex = 0.72)
  text(mean(plot_bounds), axis_y - 0.75, "Adjusted odds ratio", cex = 0.80)
}

output_dir <- file.path("generated", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_svg <- file.path(output_dir, "FigureS4_bisphenols.svg")
draw_bisphenol_forest(results, output_svg)
message("Figure S4 saved to ", output_svg)