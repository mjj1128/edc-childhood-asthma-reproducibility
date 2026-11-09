if (.Platform$OS.type == "windows") {
  invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
  invisible(Sys.setlocale("LC_COLLATE", "Chinese (Simplified)_China.utf8"))
}

local_lib <- "E:/codex 1edc/manuscript_review/ml_minimal_revision/Rlib"
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))
suppressPackageStartupMessages(library(ggplot2))

plot_dir <- "D:/环境/meta数据csv/EDC_Overall_AnalysisA_SVG/三水平出图"
run_dir <- file.path(plot_dir, "补充分组验证_run_20261107_200813")
output_svg <- file.path(plot_dir, "Figure5_Spearman_training_no_stars.svg")

read_input <- function(name) {
  read.csv(
    file.path(run_dir, name),
    check.names = FALSE,
    fileEncoding = "UTF-8-BOM",
    stringsAsFactors = FALSE
  )
}

data <- read_input("analysis_input_used.csv")
split <- read_input("outer_groups.csv")

stopifnot(
  !anyDuplicated(data$Effect_ID),
  !anyDuplicated(split$Effect_ID),
  setequal(data$Effect_ID, split$Effect_ID)
)

train <- split[which(split$Set == "Train"), ]
test <- split[which(split$Set == "Test"), ]
stopifnot(
  nrow(train) == 174L,
  nrow(test) == 37L,
  length(unique(train$Validation_group)) == 28L,
  !length(intersect(train$Validation_group, test$Validation_group))
)

numeric_features <- c("Effective sample size", "Mean age", "Male proportion")
categorical_features <- c(
  "EDC class", "Chemical form", "Exposure window", "Study design",
  "Region", "Outcome group", "Exposure contrast"
)
vars <- c(numeric_features, categorical_features)
features <- data[match(train$Effect_ID, data$Effect_ID), vars, drop = FALSE]

# 沿用原训练集相关图的缺失值填补与类别编码。
for (name in numeric_features) {
  x <- suppressWarnings(as.numeric(features[[name]]))
  x[!is.finite(x)] <- NA_real_
  x[is.na(x)] <- median(x, na.rm = TRUE)
  features[[name]] <- x
}

for (name in categorical_features) {
  x <- trimws(as.character(features[[name]]))
  x[is.na(x) | !nzchar(x)] <- "Missing"
  features[[name]] <- as.integer(factor(x))
}

# 只计算 Spearman 相关系数。
rho <- cor(features, method = "spearman")
stopifnot(all(is.finite(rho)))

cells <- as.data.frame(as.table(rho), responseName = "r")
upper <- match(cells$Var1, vars) < match(cells$Var2, vars)
cells$Var1 <- factor(cells$Var1, levels = vars)
cells$Var2 <- factor(cells$Var2, levels = rev(vars))
cells$label <- sprintf("%.2f", cells$r)
cells$text_col <- ifelse(abs(cells$r) > 0.55, "white", "black")

p <- ggplot(cells, aes(x = Var1, y = Var2)) +
  geom_tile(
    data = cells[!upper, ], aes(fill = r),
    color = "white", linewidth = 1.15
  ) +
  geom_text(
    data = cells[!upper, ], aes(label = label, color = text_col),
    size = 3.6, fontface = "bold"
  ) +
  geom_point(
    data = cells[upper, ], aes(fill = r, size = abs(r)),
    shape = 21, color = "grey35", stroke = 0.55, alpha = 0.95
  ) +
  scale_fill_gradient2(
    low = "#6F9FCF", mid = "#FFFFFF", high = "#E64B35",
    midpoint = 0, limits = c(-1, 1), name = "Spearman\nCorrelation"
  ) +
  scale_color_identity() +
  scale_size_continuous(range = c(2.6, 9.2), guide = "none") +
  labs(title = "Feature Correlation Matrix", x = NULL, y = NULL) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22, color = "black"),
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1,
      face = "bold", size = 11, color = "black"
    ),
    axis.text.y = element_text(face = "bold", size = 11, color = "black"),
    axis.line = element_blank(),
    panel.border = element_rect(fill = NA, color = "#555555", linewidth = 0.85),
    panel.grid = element_blank(),
    legend.title = element_text(face = "bold", size = 11, color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    legend.key.height = grid::unit(1.15, "cm"),
    plot.margin = ggplot2::margin(16, 26, 16, 26),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  output_svg, plot = p,
  width = 12.2, height = 10.3, units = "in",
  device = svglite::svglite, bg = "white"
)
cat("Figure 5 已保存到 ", output_svg, "\n", sep = "")

