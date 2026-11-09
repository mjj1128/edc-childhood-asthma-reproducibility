# ============================================================
# EDC 研究特征弦图：大字号版
# 仅保留弦图和类别标签
#
# 输出：
# 1. EDC 类别 ↔ 地区
# 2. EDC 类别 ↔ 研究设计
# 3. EDC 类别 ↔ 检测基质
# 4. EDC 类别 ↔ 暴露时间窗
# 5. 四联图
#
# 无标题、无面板编号、无底部图注
# 不使用 readr/vroom
# ============================================================


# 1. 文件路径 -------------------------------------------------

input_csv <- "D:/环境/meta数据csv/三水平_meta_AnalysisA_最终版.csv"

output_dir <- paste0(
  "D:/环境/meta数据csv/",
  "EDC_Overall_AnalysisA_SVG/三水平出图"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# 2. 安装和加载包 ---------------------------------------------

options(timeout = max(1200, getOption("timeout", 60)))

if (!requireNamespace("circlize", quietly = TRUE)) {
  install.packages(
    "circlize",
    repos = "https://cloud.r-project.org",
    type = "binary",
    method = "libcurl"
  )
}

if (!requireNamespace("circlize", quietly = TRUE)) {
  stop(
    "circlize 未成功安装或无法加载，请查看上方报错。",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(circlize)
})


# 3. CSV 读写函数 ---------------------------------------------

read_csv_base <- function(path) {
  
  if (!file.exists(path)) {
    stop("找不到文件：", path, call. = FALSE)
  }
  
  read.csv(
    file = path,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fileEncoding = "UTF-8-BOM",
    na.strings = c("", "NA")
  )
}

write_csv_base <- function(data, path) {
  
  con <- file(
    path,
    open = "wt",
    encoding = "UTF-8"
  )
  
  on.exit(close(con), add = TRUE)
  
  cat("\ufeff", file = con)
  
  write.table(
    data,
    file = con,
    sep = ",",
    row.names = FALSE,
    col.names = TRUE,
    quote = TRUE,
    na = "",
    qmethod = "double"
  )
}


# 4. 读取并筛选数据 -------------------------------------------

raw <- read_csv_base(input_csv)

required_columns <- c(
  "Study-ID",
  "Analysis_A",
  "EDC_class_Clean",
  "Region_Clean",
  "Design_Clean",
  "Matrix_Clean",
  "Exposure_wind"
)

missing_columns <- setdiff(
  required_columns,
  names(raw)
)

if (length(missing_columns) > 0) {
  stop(
    "数据表缺少字段：",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

keep <- !is.na(raw$Analysis_A) & raw$Analysis_A == 1
dat <- raw[keep, , drop = FALSE]

if (nrow(dat) == 0) {
  stop(
    "没有 Analysis_A = 1 的记录。",
    call. = FALSE
  )
}

text_columns <- setdiff(
  required_columns,
  "Analysis_A"
)

for (column in text_columns) {
  dat[[column]] <- trimws(
    as.character(dat[[column]])
  )
}

if (any(
  is.na(dat[["Study-ID"]]) |
  dat[["Study-ID"]] == ""
)) {
  stop(
    "正式纳入的数据存在缺失 Study-ID。",
    call. = FALSE
  )
}

feature_columns <- c(
  "EDC_class_Clean",
  "Region_Clean",
  "Design_Clean",
  "Matrix_Clean",
  "Exposure_wind"
)

# 保留缺失类别，不静默删除报告
for (column in feature_columns) {
  
  missing_value <- (
    is.na(dat[[column]]) |
      dat[[column]] == ""
  )
  
  dat[[column]][missing_value] <- "Not reported"
}

n_reports <- length(
  unique(dat[["Study-ID"]])
)

message(
  "Analysis A：",
  n_reports,
  " 篇报告；",
  nrow(dat),
  " 条效应量记录。"
)

if (
  n_reports != 37 ||
  nrow(dat) != 212
) {
  warning(
    "当前数量与此前核对的 37 篇报告、212 个效应量不同，请确认版本。"
  )
}

# 保留原始分类，不自动合并 Childhood、Childhood (2) 等标签。
# 如需合并，须先确认这些编码的实际含义。


# 5. EDC 类别顺序和配色 ---------------------------------------

edc_order_preferred <- c(
  "Phthalates",
  "Bisphenols",
  "Phenols",
  "Parabens",
  "PAHs",
  "Metals",
  "BFRs",
  "OPEs",
  "PFAS",
  "Pesticides",
  "DINCH",
  "VOCs",
  "Mixed EDCs",
  "Not reported"
)

observed_edc <- unique(dat$EDC_class_Clean)

edc_order <- c(
  edc_order_preferred[
    edc_order_preferred %in% observed_edc
  ],
  sort(
    setdiff(
      observed_edc,
      edc_order_preferred
    )
  )
)

edc_colours <- c(
  "Phthalates"   = "#4C88B5",
  "Bisphenols"   = "#E88973",
  "Phenols"      = "#63B4A5",
  "Parabens"     = "#B49ACB",
  "PAHs"         = "#D8B15A",
  "Metals"       = "#C86B78",
  "BFRs"         = "#7A91B3",
  "OPEs"         = "#A5B66C",
  "PFAS"         = "#6BBCC4",
  "Pesticides"   = "#C99567",
  "DINCH"        = "#9A82A8",
  "VOCs"         = "#91A6A6",
  "Mixed EDCs"   = "#786F65",
  "Not reported" = "#BDBDBD"
)

extra_classes <- setdiff(
  edc_order,
  names(edc_colours)
)

if (length(extra_classes) > 0) {
  
  extra_colours <- setNames(
    grDevices::hcl.colors(
      length(extra_classes),
      palette = "Dark 3"
    ),
    extra_classes
  )
  
  edc_colours <- c(
    edc_colours,
    extra_colours
  )
}


# 6. 生成去重后的连接计数 -------------------------------------
# 同一报告在一个“EDC 类别 × 特征类别”组合内只计一次。
# 一篇报告可贡献多个组合，连接总数可能超过报告总数。

make_links <- function(feature_column) {
  
  records <- unique(
    data.frame(
      Study_ID = dat[["Study-ID"]],
      EDC_class = dat$EDC_class_Clean,
      Feature = dat[[feature_column]],
      stringsAsFactors = FALSE
    )
  )
  
  links <- aggregate(
    x = list(
      n = rep(1L, nrow(records))
    ),
    by = list(
      EDC_class = records$EDC_class,
      Feature = records$Feature
    ),
    FUN = sum
  )
  
  links <- links[
    order(
      match(
        links$EDC_class,
        edc_order
      ),
      links$Feature
    ),
    ,
    drop = FALSE
  ]
  
  rownames(links) <- NULL
  
  links
}

links_region <- make_links("Region_Clean")
links_design <- make_links("Design_Clean")
links_matrix <- make_links("Matrix_Clean")
links_window <- make_links("Exposure_wind")


# 7. 保存计数表 -----------------------------------------------

count_tables <- list(
  "Chord_EDC_Region_Counts.csv" = links_region,
  "Chord_EDC_Design_Counts.csv" = links_design,
  "Chord_EDC_Matrix_Counts.csv" = links_matrix,
  "Chord_EDC_Window_Counts.csv" = links_window
)

for (filename in names(count_tables)) {
  write_csv_base(
    count_tables[[filename]],
    file.path(
      output_dir,
      filename
    )
  )
}


# 8. 字号设置与标签换行 ---------------------------------------
# 字号越大，文字越大。
# 单图原为 0.78，现为 1.15。
# 四联图原为 0.74，现为 1.10。
# 如还需要更大，可尝试分别改为 1.25 和 1.20。

label_cex_single <- 1.15
label_cex_combined <- 1.10

format_label <- function(x) {
  paste(
    strwrap(
      x,
      width = 18
    ),
    collapse = "\n"
  )
}


# 9. 弦图函数 -------------------------------------------------
# 不添加总标题、面板标题、字母编号或底部说明。

draw_chord <- function(
    links,
    feature_order = NULL,
    label_cex = label_cex_single) {
  
  circlize::circos.clear()
  
  on.exit(
    circlize::circos.clear(),
    add = TRUE
  )
  
  left_present <- edc_order[
    edc_order %in% links$EDC_class
  ]
  
  right_present <- unique(links$Feature)
  
  if (is.null(feature_order)) {
    
    right_totals <- aggregate(
      n ~ Feature,
      data = links,
      FUN = sum
    )
    
    right_totals <- right_totals[
      order(
        -right_totals$n,
        right_totals$Feature
      ),
      ,
      drop = FALSE
    ]
    
    right_order <- right_totals$Feature
    
  } else {
    
    right_order <- c(
      feature_order[
        feature_order %in% right_present
      ],
      sort(
        setdiff(
          right_present,
          feature_order
        )
      )
    )
  }
  
  # 两端使用不同内部名称，防止同名类别合并
  left_ids <- paste0(
    "EDC__",
    left_present
  )
  
  right_ids <- paste0(
    "FEATURE__",
    right_order
  )
  
  sector_order <- c(
    left_ids,
    right_ids
  )
  
  label_lookup <- c(
    setNames(
      left_present,
      left_ids
    ),
    setNames(
      right_order,
      right_ids
    )
  )
  
  right_colours <- setNames(
    grDevices::hcl.colors(
      length(right_order),
      palette = "Set 2"
    ),
    right_ids
  )
  
  grid_colours <- c(
    setNames(
      unname(
        edc_colours[left_present]
      ),
      left_ids
    ),
    right_colours
  )
  
  plot_links <- data.frame(
    from = paste0(
      "EDC__",
      links$EDC_class
    ),
    to = paste0(
      "FEATURE__",
      links$Feature
    ),
    value = links$n,
    stringsAsFactors = FALSE
  )
  
  link_colours <- unname(
    edc_colours[links$EDC_class]
  )
  
  # 两组类别之间留较大的空隙
  gap_after <- c(
    rep(
      2,
      length(left_ids) - 1
    ),
    13,
    rep(
      3,
      length(right_ids) - 1
    ),
    13
  )
  
  circlize::circos.par(
    start.degree = 90,
    gap.after = gap_after,
    track.margin = c(0.005, 0.005),
    cell.padding = c(0, 0, 0, 0),
    
    # 扩大画布范围，为大字号标签留出空间
    canvas.xlim = c(-1.25, 1.25),
    canvas.ylim = c(-1.25, 1.25),
    
    points.overflow.warning = FALSE
  )
  
  circlize::chordDiagram(
    x = plot_links,
    order = sector_order,
    grid.col = grid_colours,
    col = link_colours,
    transparency = 0.62,
    annotationTrack = "grid",
    annotationTrackHeight = 0.035,
    
    # 加宽标签轨道，原为 0.17
    preAllocateTracks = list(
      track.height = 0.23
    ),
    
    directional = 0,
    link.border = NA,
    link.sort = TRUE,
    link.largest.ontop = TRUE,
    reduce = 0
  )
  
  # 外圈类别名称
  circlize::circos.trackPlotRegion(
    track.index = 1,
    bg.border = NA,
    panel.fun = function(x, y) {
      
      sector_id <- circlize::get.cell.meta.data(
        "sector.index"
      )
      
      xlim <- circlize::get.cell.meta.data(
        "xlim"
      )
      
      ylim <- circlize::get.cell.meta.data(
        "ylim"
      )
      
      sector_label <- format_label(
        unname(
          label_lookup[sector_id]
        )
      )
      
      circlize::circos.text(
        x = mean(xlim),
        y = mean(ylim),
        labels = sector_label,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0.5, 0.5),
        cex = label_cex,
        col = "#253D4D"
      )
    }
  )
}


# 10. 各图顺序设置 --------------------------------------------

region_order <- c(
  "Asia",
  "Americas",
  "Europe",
  "Africa",
  "Not reported"
)

design_order <- c(
  "Cross-sectional",
  "Cohort",
  "Case-control",
  "Nested case-control",
  "Case-cohort",
  "Not reported"
)

panels <- list(
  list(
    links = links_region,
    file = "Chord_A_EDC_Region.svg",
    order = region_order
  ),
  list(
    links = links_design,
    file = "Chord_B_EDC_Study_Design.svg",
    order = design_order
  ),
  list(
    links = links_matrix,
    file = "Chord_C_EDC_Matrix.svg",
    order = NULL
  ),
  list(
    links = links_window,
    file = "Chord_D_EDC_Exposure_Window.svg",
    order = NULL
  )
)


# 11. 保存四张单独 SVG ----------------------------------------

save_single <- function(panel) {
  
  output_file <- file.path(
    output_dir,
    panel$file
  )
  
  grDevices::svg(
    filename = output_file,
    width = 9,
    height = 9,
    family = "sans",
    bg = "white",
    onefile = TRUE
  )
  
  on.exit(
    grDevices::dev.off(),
    add = TRUE
  )
  
  graphics::par(
    mar = c(0.5, 0.5, 0.5, 0.5),
    cex = 1,
    xpd = NA
  )
  
  draw_chord(
    links = panel$links,
    feature_order = panel$order,
    label_cex = label_cex_single
  )
  
  message(
    "已保存：",
    output_file
  )
}

for (panel in panels) {
  save_single(panel)
}


# 12. 保存四联 SVG --------------------------------------------
# 面板顺序：
# 左上：地区
# 右上：研究设计
# 左下：检测基质
# 右下：暴露时间窗

save_combined <- function() {
  
  output_file <- file.path(
    output_dir,
    "Chord_EDC_Four_Panels.svg"
  )
  
  grDevices::svg(
    filename = output_file,
    width = 16,
    height = 16,
    family = "sans",
    bg = "white",
    onefile = TRUE
  )
  
  on.exit(
    grDevices::dev.off(),
    add = TRUE
  )
  
  graphics::par(
    mfrow = c(2, 2),
    
    # 抵消四联图布局对基础字号的自动缩小
    cex = 1,
    
    mar = c(0.5, 0.5, 0.5, 0.5),
    oma = c(0, 0, 0, 0),
    xpd = NA
  )
  
  for (panel in panels) {
    draw_chord(
      links = panel$links,
      feature_order = panel$order,
      label_cex = label_cex_combined
    )
  }
  
  message(
    "已保存四联图：",
    output_file
  )
}

save_combined()


# 13. 完成 ----------------------------------------------------

message(
  "\n出图完成，保存目录：\n",
  output_dir,
  "\n\n四联图文件：Chord_EDC_Four_Panels.svg",
  "\n已放大类别标签，并加宽标签区域。",
  "\n已去除所有标题、面板编号和底部说明。"
)

