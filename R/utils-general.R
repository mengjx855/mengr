#### Jin-Xin Meng, jinxmeng@zju.edu.cn, 20230610, 20260923 ####

# 20231101: update function: get_freq
# 20250223: undate functions with new grammar.
# 20250223: add function: get_text_color to decide the text 'black' or 'white'.
# 20250309: add function: write_xlsx_with_comment()
# 20260502: delete get_freq(), intersect_multiple(), stat_vec(), calcu_adjusted_r2(),
#           update other functions.
# 20260519: add function 'pairwise_cluster()'.
# 20260527: add function 'set_calcu()'.
# 20260902: add function 'log_message()'.
# 20260916: standardize script metadata, function sections, documentation, and naming style.
# 20260923: remove Chinese Roxygen text, correct export tags, and fix variable references.

## 日志记录函数，支持不同类型的日志信息

#### log_message ####

log_message <- function(
    ..., type = c("info", "success", "warning", "error", "debug")
) {
  type <- match.arg(type)
  time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  text <- paste0(...)
  label <- switch(
    type,
    info    = cli::col_blue("INFO"),
    success = cli::col_green("SUCCESS"),
    warning = cli::col_yellow("WARNING"),
    error   = cli::col_red("ERROR"),
    debug   = cli::col_magenta("DEBUG")
  )

  message(cli::col_grey(time), " | ", label, ": ", text)

  invisible(NULL)
}


#### floor_n ####
#' Floor n utility
#'
#'
#'
#' @param x Primary vector or object supplied to the utility.
#' @param n Requested number of values, features, or results.
#' @return A numeric vector rounded downward at the requested decimal position.
#' @export
floor_n <- function(x, n = 2) {
  base <- 10^n
  floor(x / base) * base
}

#' Ceiling n utility
#'
#'
#'
#' @param x Primary vector or object supplied to the utility.
#' @param n Requested number of values, features, or results.
#' @return A numeric vector rounded upward at the requested decimal position.
#' @export
#### ceiling_n ####

ceiling_n <- function(x, n = 2) {
  base <- 10^n
  ceiling(x / base) * base
}

#### set_calcu ####
# 多个 vector 的集合运算
# ...: 多个 vector，或者一个 lis
# method:
#   intersect: 多个 vector 的交集
#   union: 多个 vector 的并集
#   setdiff: 第一个 vector 减去后面所有 vector 的并集
# unique_out: 是否对结果去重
# sort_out: 是否排序输出

#' Set Calcu utility
#'
#'
#'
#' @param ... Vectors to combine; a single list of vectors is also accepted.
#' @param method Analysis or summary method; supported values are shown in the usage.
#' @param unique_out Whether duplicate values are removed from the set-operation result.
#' @param sort_out Whether the returned set-operation result is sorted.
#' @return A vector containing the requested intersection, union, or sequential set difference.
#' @export
set_calcu <- function(..., method = c("intersect", "union", "setdiff"),
                      unique_out = TRUE, sort_out = FALSE) {
  method <- match.arg(method)

  x <- list(...)

  ## 支持直接输入一个 lis
  if (length(x) == 1 && is.list(x[[1]])) {
    x <- x[[1]]
  }

  ## 去掉 NULL
  x <- x[!vapply(x, is.null, logical(1))]

  if (length(x) == 0) {
    return(character())
  }

  x <- lapply(x, as.vector)

  if (isTRUE(unique_out)) {
    x <- lapply(x, unique)
  }

  if (method == "intersect") {
    out <- Reduce(intersect, x)
  } else if (method == "union") {
    out <- Reduce(union, x)
  } else if (method == "setdiff") {
    if (length(x) == 1) {
      out <- x[[1]]
    } else {
      out <- setdiff(x[[1]], Reduce(union, x[-1]))
    }
  }

  if (isTRUE(sort_out)) {
    out <- sort(out)
  }

  return(out)
}

#### get_text_color ####
# 计算亮度（Luminance）公式：0.299*R + 0.587*G + 0.114*B
#' Get Text Color utility
#'
#'
#'
#' @param color Color specification for `color`.
#' @param threshold Luminance threshold separating dark and light text.
#' @param dark Text color returned for a sufficiently light background.
#' @param light Text color returned for a sufficiently dark background.
#' @return A character vector choosing `dark` or `light` for each input color.
#' @export
get_text_color <- function(color, threshold = 0.5, dark = "black", light = "white") {
  rgb_mat <- grDevices::col2rgb(color) # 将颜色转换为 RGB 值
  luminance <- (0.299 * rgb_mat[1, ] + 0.587 * rgb_mat[2, ] + 0.114 * rgb_mat[3, ]) / 255
  out <- ifelse(luminance >= threshold, dark, light)
  names(out) <- color
  return(out)
}

#### write_xlsx_with_comment ####
#' Write xlsx with Comment utility
#'
#'
#'
#' @param data An input data frame or compatible object.
#' @param filename Path of the workbook or output file.
#' @param comment Character marker written before comment rows in the worksheet.
#' @param sheet Worksheet name used for Excel input or output.
#' @param comment_color Color specification for `comment_color`.
#' @param overwrite Whether an existing workbook may be replaced.
#' @param append Whether to append a worksheet to an existing workbook.
#' @param replace_sheet Whether an existing worksheet with the same name is replaced.
#' @return The output workbook path, invisibly.
#' @export
write_xlsx_with_comment <- function(data, filename, comment = "###", sheet = "Sheet1",
                                    comment_color = "red", overwrite = TRUE,
                                    append = TRUE, replace_sheet = FALSE) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required.")
  }

  if (is.null(comment)) {
    comment <- character(0)
  }

  if (is.list(comment)) {
    comment <- unlist(comment, use.names = FALSE)
  }

  # 1. 如果文件存在并且允许追加，则读取已有 workbook
  if (file.exists(filename) && append) {
    wb <- openxlsx::loadWorkbook(filename)
  } else {
    wb <- openxlsx::createWorkbook()
  }

  # 2. 如果 sheet 已存在，决定是否替换
  existing_sheets <- names(wb)
  if (sheet %in% existing_sheets) {
    if (replace_sheet) {
      openxlsx::removeWorksheet(wb, sheet = sheet)
    } else {
      stop(sprintf(
        "Sheet '%s' already exists. Use replace_sheet = TRUE if you want to overwrite it.",
        sheet
      ))
    }
  }

  # 3. 新增 shee
  openxlsx::addWorksheet(wb, sheet)

  comment_style <- openxlsx::createStyle(fontColour = comment_color)

  # 4. 写入顶部注释
  if (length(comment) > 0) {
    for (i in seq_along(comment)) {
      openxlsx::writeData(
        wb,
        sheet = sheet, x = comment[i], startRow = i,
        startCol = 1, colNames = FALSE
      )
      openxlsx::addStyle(wb,
        sheet = sheet, style = comment_style,
        rows = i, cols = 1, gridExpand = TRUE
      )
    }
    start_row <- length(comment) + 2
  } else {
    start_row <- 1
  }

  # 5. 写入数据
  openxlsx::writeData(
    wb,
    sheet = sheet, x = data, startRow = start_row, startCol = 1
  )

  # 6. 保存
  openxlsx::saveWorkbook(wb, file = filename, overwrite = overwrite)
}

#### read_xlsx_multiple ####
#' Read xlsx Multiple utility
#'
#'
#'
#' @param file Path to an input file.
#' @param sheets Worksheet names to read; `NULL` reads every worksheet.
#' @param ... Additional arguments passed to `openxlsx::read.xlsx()`.
#' @return A named list of data frames, one for each selected worksheet.
#' @export
read_xlsx_multiple <- function(file, sheets = NULL, ...) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required.")
  }

  if (is.null(sheets)) {
    sheets <- openxlsx::getSheetNames(file)
  }

  data <- lapply(
    sheets,
    \(x) openxlsx::read.xlsx(file, sheet = x, ...)
  )

  names(data) <- sheets

  return(data)
}

#### pairwise_cluster ####
# 根据 feature1-feature2-value 三列表构建距离矩阵，并进行层次聚类。
# data:
#   默认前三列分别为 feature1, feature2, value；
#   也可以通过 feature1_col / feature2_col / value_col 指定列名或列号。
#
# 参数:
#   one_minus:
#     FALSE: value 本身就是 distance，例如 1-ANI、1-|cor|、Bray distance
#     TRUE : value 是 similarity，例如 ANI、correlation，需要转换为 distance = 1 - value
#
#   cutoff:
#     用于 cutree(h = cutoff) 的距离阈值。
#     如果 ANI cutoff = 0.95 且 one_minus = TRUE，则 cutoff 应写为 1 - 0.95 = 0.05。
#
#   abs_value:
#     TRUE 时先对 value 取绝对值，适合 correlation，例如 distance = 1 - abs(cor)。
#
# 输出:
#   simplify = TRUE : named vector，名字为 feature，值为 cluster
#   simplify = FALSE: data.frame(name, cluster)
#
# 用法：
#   1. fastANI
#     cl <- pairwise_cluster(
#      data = ani_df, feature1_col = "ref", feature2_col = "query",
#      value_col = "ani", one_minus = TRUE, cutoff = 1 - 0.95,
#      linkage_method = "average", simplify = FALSE
#    )
#   2. 相关性聚类
#     cl <- pairwise_cluster(
#       data = cor_df, feature1_col = "feature1", feature2_col = "feature2",
#       value_col = "cor", one_minus = TRUE, abs_value = TRUE,
#       cutoff = 1 - 0.7, simplify = FALSE
#     )
#   3. value 本身就是距离
#     cl <- pairwise_cluster(
#       data = dist_df, one_minus = FALSE, cutoff = 0.2, simplify = TRUE
#     )

#' Pairwise Cluster utility
#'
#'
#'
#' @param data An input data frame or compatible object.
#' @param feature1_col Name of the `feature1_col` input column.
#' @param feature2_col Name of the `feature2_col` input column.
#' @param value_col Name of the `value_col` input column.
#' @param cutoff Threshold used to form clusters or significance calls, depending on the function.
#' @param one_minus Whether pairwise values are converted to `1 - value` before clustering.
#' @param abs_value Whether to cluster using absolute pairwise values.
#' @param linkage_method Hierarchical-clustering linkage method passed to `stats::hclust()`.
#' @param fill_missing Distance assigned to feature pairs absent from a pairwise table.
#' @param duplicate_fun Function used to combine duplicated feature-pair values.
#' @param simplify Whether to return the simplified tabular result instead of intermediate objects.
#' @return A cluster-membership vector or data frame with `hclust` and distance-matrix attributes.
#' @export
pairwise_cluster <- function(data, feature1_col = NULL, feature2_col = NULL,
                             value_col = NULL, cutoff = 0.05, one_minus = FALSE,
                             abs_value = FALSE, linkage_method = "average",
                             fill_missing = 1, duplicate_fun = mean,
                             simplify = TRUE) {
  data <- data.frame(data, check.names = FALSE)

  ## 默认使用前三列
  if (is.null(feature1_col)) feature1_col <- colnames(data)[1]
  if (is.null(feature2_col)) feature2_col <- colnames(data)[2]
  if (is.null(value_col)) value_col <- colnames(data)[3]

  ## 支持列号指定
  if (is.numeric(feature1_col)) feature1_col <- colnames(data)[feature1_col]
  if (is.numeric(feature2_col)) feature2_col <- colnames(data)[feature2_col]
  if (is.numeric(value_col)) value_col <- colnames(data)[value_col]

  if (!all(c(feature1_col, feature2_col, value_col) %in% colnames(data))) {
    stop("feature1_col / feature2_col / value_col not found in data.")
  }

  linkage_method <- match.arg(
    linkage_method,
    c(
      "average", "complete", "single", "median", "centroid",
      "ward.D", "ward.D2", "mcquitty", "weighted", "ward"
    )
  )

  ## 兼容 scipy 的命名
  if (linkage_method == "weighted") linkage_method <- "mcquitty"
  if (linkage_method == "ward") linkage_method <- "ward.D2"

  ## 整理三列表
  data <- data.frame(
    feature1 = as.character(data[[feature1_col]]),
    feature2 = as.character(data[[feature2_col]]),
    value = suppressWarnings(as.numeric(data[[value_col]])),
    stringsAsFactors = FALSE
  )

  data <- data |>
    dplyr::filter(
      !is.na(feature1),
      !is.na(feature2),
      feature1 != "",
      feature2 != "",
      is.finite(value)
    )

  if (nrow(data) == 0) {
    stop("No valid pairwise records.")
  }

  ## 如果是相关性，可以先取绝对值
  if (isTRUE(abs_value)) {
    data$value <- abs(data$value)
  }

  ## value 转 distance
  data$distance <- if (isTRUE(one_minus)) {
    1 - data$value
  } else {
    data$value
  }

  if (any(data$distance < 0, na.rm = TRUE)) {
    warning("Some distances are < 0. Please check value scale or one_minus setting.")
  }

  ## A-B 和 B-A 统一成一个方向
  ## pmin(feature1, feature2)：一行一行地比对，一行中的值相比，谁小谁就去 node1
  data <- data |>
    dplyr::mutate(
      node1 = pmin(feature1, feature2),
      node2 = pmax(feature1, feature2)
    )

  ## 重复 pair 合并，默认取平均 distance
  data2 <- data |>
    dplyr::group_by(node1, node2) |>
    dplyr::summarise(
      distance = duplicate_fun(distance, na.rm = TRUE),
      .groups = "drop"
    )

  ## 构建距离矩阵
  features <- sort(unique(c(data2$node1, data2$node2)))

  dist_mat <- matrix(
    fill_missing,
    nrow = length(features),
    ncol = length(features),
    dimnames = list(features, features)
  )

  diag(dist_mat) <- 0

  for (i in seq_len(nrow(data2))) {
    a <- data2$node1[i]
    b <- data2$node2[i]
    d <- data2$distance[i]

    dist_mat[a, b] <- d
    dist_mat[b, a] <- d
  }

  ## 层次聚类
  hc <- stats::hclust(
    stats::as.dist(dist_mat),
    method = linkage_method
  )

  cluster <- stats::cutree(hc, h = cutoff)

  if (isTRUE(simplify)) {
    out <- cluster
  } else {
    out <- data.frame(
      name = names(cluster),
      cluster = as.integer(cluster),
      row.names = NULL,
      check.names = FALSE
    )
  }

  attr(out, "hclust") <- hc
  attr(out, "dist_matrix") <- dist_mat

  return(out)
}


#### theme_bw_clean ####
#' Theme Bw Clean utility
#'
#'
#'
#' @param base_size Base font size for the plot theme.
#' @return A ggplot2 theme object.
#' @export
theme_bw_clean <- function(base_size = 12) {
  ## 返回可继续用“+”叠加修改的 ggplot2 theme 对象
  ggplot2::theme(
    axis.ticks = ggplot2::element_line(linewidth = 0.5),
    axis.ticks.length = grid::unit(2, "mm"),
    axis.text = ggplot2::element_text(size = base_size, color = "black"),
    axis.title = ggplot2::element_text(size = base_size, color = "black"),
    strip.background = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(
      size = base_size + 1, color = "black", face = "bold", hjust = 0.5
    ),
    plot.background = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(linewidth = 0.5),
    panel.grid.minor = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(
      fill = NA, linewidth = 0.5, color = "black"
    ),
    panel.background = ggplot2::element_blank(),
    panel.spacing = grid::unit(0, "mm")
  )
}
