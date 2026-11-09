if (.Platform$OS.type == "windows") {
  invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
}

plot_dir <- "D:/环境/meta数据csv/EDC_Overall_AnalysisA_SVG/三水平出图"
source_svg <- file.path(
  plot_dir, "三水平最终出图",
  "Figure1_Overall_EDC_class_associations_compact.svg"
)
output_svg <- file.path(plot_dir, "Figure3_descriptive.svg")

svg <- readLines(source_svg, encoding = "UTF-8")

# 按表头的位置识别 P 值列。
p_header <- grep("^<text .*?>P</text>$", svg)
stopifnot(length(p_header) == 1L)
p_x <- sub("^<text x='([^']+)'.*$", "\\1", svg[p_header])
p_column <- startsWith(svg, paste0("<text x='", p_x, "' "))
p_rows <- integer()

for (edc in c("Parabens", "OPEs")) {
  label_row <- which(endsWith(svg, paste0(">", edc, "</text>")))
  stopifnot(length(label_row) == 1L, label_row + 3L <= length(svg))
  
  # 原 SVG 每行依次写入类别名、aOR [95% CI]、计数和 P 值。
  p_row <- label_row + 3L
  stopifnot(
    grepl(">n=2; k=[0-9]+</text>$", svg[label_row + 2L]),
    p_column[p_row],
    grepl(">0\\.[0-9]+</text>$", svg[p_row])
  )
  
  # 释放原标签的固定宽度，使新增文字正常排版。
  label <- sub(" textLength='[^']*' lengthAdjust='[^']*'", "", svg[label_row])
  svg[label_row] <- sub("</text>$", " (descriptive)</text>", label)
  p_rows <- c(p_rows, p_row)
}

writeLines(svg[-p_rows], output_svg, useBytes = TRUE)
cat("Figure 3 已保存到 ", output_svg, "\n", sep = "")

