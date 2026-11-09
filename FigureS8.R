# Figure S8：从已有结果重绘方差分解图

if (.Platform$OS.type == "windows") {
  invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
}

local_lib <- "E:/codex 1edc/manuscript_review/ml_minimal_revision/Rlib"
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))
suppressPackageStartupMessages(library(ggplot2))

base_dir <- "D:/环境/meta数据csv/EDC_Overall_AnalysisA_SVG"
input_csv <- file.path(base_dir, "48_EDC_class_two_level_vs_three_level_CR2.csv")
output_dir <- file.path(base_dir, "三水平出图")

results <- read.csv(
  input_csv,
  fileEncoding = "UTF-8-BOM",
  stringsAsFactors = FALSE
)

# 指定从上到下的类别顺序。
keep <- c("Phthalates", "Phenols", "Bisphenols")
selected <- results[results$EDC_class %in% keep & results$status == "OK", ]

stopifnot(
  nrow(selected) == length(keep),
  setequal(selected$EDC_class, keep),
  all(is.finite(selected$total_model_variance) & selected$total_model_variance > 0)
)

# 沿用原图的计算方式：各方差成分除以两个成分之和。
components <- c("Between-study variance", "Within-study variance")
plot_data <- data.frame(
  EDC_class = factor(rep(selected$EDC_class, 2), levels = rev(keep)),
  Variance_component = rep(components, each = nrow(selected)),
  Proportion = c(selected$tau2_study_level3, selected$tau2_within_study_level2) /
    rep(selected$total_model_variance, 2)
)

p <- ggplot(plot_data, aes(x = EDC_class, y = Proportion, fill = Variance_component)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.45) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Between-study variance" = "#2166AC",
      "Within-study variance" = "#F4A261"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = c("0%", "25%", "50%", "75%", "100%"),
    oob = scales::squish
  ) +
  labs(
    title = "Variance decomposition in three-level EDC meta-analysis",
    x = NULL,
    y = "Proportion of model variance",
    fill = NULL
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(
      hjust = 0.5, face = "bold", size = 19, margin = margin(b = 14)
    ),
    axis.text.y = element_text(size = 14, face = "bold", colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    legend.position = "bottom",
    legend.text = element_text(size = 13),
    plot.margin = margin(15, 20, 18, 20)
  )

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_svg <- file.path(output_dir, "FigureS8_EDC_class_variance_decomposition.svg")

ggsave(output_svg, plot = p, width = 11.5, height = 5.5, device = svglite::svglite)
cat("Figure S8 已保存到 ", output_svg, "\n", sep = "")

