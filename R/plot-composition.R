#### Jin-Xin Meng, jinxmeng@zju.edu.cn, 20231113, 20260925 ####

# 20260916: standardize script metadata, function sections, documentation, and naming style.
# 20260923: remove Chinese text from Roxygen documentation.
# 20260925: add plot_circlepack(), document it, and rename this composition-plot file.

#### plot_pie ####

#' Plot a pie chart
#'
#'
#' @param data An input data frame or compatible object.
#' @param name_col Name of the `name_col` input column.
#' @param value_col Name of the `value_col` input column.
#' @param name_level Optional order for `name_level`.
#' @param decreasing Whether categories are sorted in decreasing abundance order.
#' @param top_n Number of highest-ranking features or categories retained.
#' @param other_name Label assigned to features combined into the residual category.
#' @param other_last Whether the residual `Other` category is placed last.
#' @param add_count Whether to include counts in slice labels.
#' @param add_percent Whether to include percentages in slice labels.
#' @param circular_label Whether category labels should follow the circular plotting direction.
#' @param flip_label Whether to reverse labels on the left half of the pie for readability.
#' @param title Optional plot or result title.
#' @param border_color Color specification for `border_color`.
#' @param fill Color specification for fill.
#' @param font_size Base text size used in the plot.
#' @param hemisphere Portion of the circle used for the pie layout; supported values are shown in Usage.
#' @param start Starting angle, in radians, used for the circular coordinate system.
#' @param percent_digits Number of decimal places shown in percentage labels.
#' @return A ggplot-compatible plot object; computed data or fitted objects are retained as attributes when applicable.
#' @export
plot_pie <- function(
  data, name_col = "name", value_col = "n", name_level = NULL,
  decreasing = TRUE, top_n = NULL, other_name = "Other",
  other_last = TRUE, add_count = FALSE, add_percent = TRUE,
  circular_label = FALSE, flip_label = FALSE,
  title = NULL, border_color = "white", fill = "auto",
  font_size = 2, hemisphere = FALSE, start = 0,
  percent_digits = 1
) {
  ## 1. 统一类别和值，并合并重复类别
  plot_df <- .as_df(data)
  .check_columns(plot_df, c(name_col, value_col), object = "data")
  plot_df <- data.frame(
    name = as.character(plot_df[[name_col]]),
    value = as.numeric(plot_df[[value_col]]),
    check.names = FALSE
  ) |>
    dplyr::filter(!is.na(name), is.finite(value), value > 0) |>
    dplyr::group_by(name) |>
    dplyr::summarise(value = sum(value), .groups = "drop")
  if (!nrow(plot_df)) stop("No valid rows remained after filtering.")

  ## 2. 排序，并把 top_n 之外的类别合并为 Other
  if (!is.null(name_level)) {
    order_vec <- match(plot_df$name, name_level)
    plot_df <- plot_df[order(order_vec, na.last = TRUE), , drop = FALSE]
  } else if (isTRUE(decreasing)) {
    plot_df <- dplyr::arrange(plot_df, dplyr::desc(value))
  } else {
    plot_df <- dplyr::arrange(plot_df, value)
  }
  if (!is.null(top_n) && nrow(plot_df) > top_n) {
    if (top_n < 2L) stop("top_n should be at least 2 when Other is required.")
    keep_name <- utils::head(plot_df$name, top_n - 1L)
    plot_df$name <- ifelse(plot_df$name %in% keep_name, plot_df$name, other_name)
    plot_df <- plot_df |>
      dplyr::group_by(name) |>
      dplyr::summarise(value = sum(value), .groups = "drop")
  }
  if (isTRUE(other_last) && other_name %in% plot_df$name) {
    plot_df <- dplyr::arrange(plot_df, name == other_name)
  }

  ## 3. 计算比例、标签位置和 circular label 角度
  plot_df <- plot_df |>
    dplyr::mutate(
      percent = value / sum(value),
      y_position = cumsum(percent) - 0.5 * percent,
      angle_position = (y_position + start / 360) %% 1,
      angle = ifelse(
        angle_position < 0.5,
        360 * angle_position + 180, 360 * angle_position
      )
    )
  if (isTRUE(add_count) && isTRUE(add_percent)) {
    plot_df$label <- paste0(
      plot_df$name, ", ", prettyNum(plot_df$value, big.mark = ","), ", ",
      round(plot_df$percent * 100, percent_digits), "%"
    )
  } else if (isTRUE(add_count)) {
    plot_df$label <- paste0(
      plot_df$name, ", ", prettyNum(plot_df$value, big.mark = ",")
    )
  } else if (isTRUE(add_percent)) {
    plot_df$label <- paste0(
      plot_df$name, ", ",
      round(plot_df$percent * 100, percent_digits), "%"
    )
  } else {
    plot_df$label <- plot_df$name
  }

  ## 4. 生成与类别一一对应的填充色
  if (length(fill) == 1L && fill == "auto") {
    palette <- c(
      "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
      "#a6d854", "#ffd92f", "#8dd3c7", "#fdb462",
      "#80b1d3", "#fccde5", "#d9ef8b", "#fee391"
    )
    fill_color <- rep(palette, length.out = nrow(plot_df))
  } else if (length(fill) == 1L && fill == "hue") {
    fill_color <- scales::hue_pal()(nrow(plot_df))
  } else if (!is.null(names(fill)) && all(plot_df$name %in% names(fill))) {
    fill_color <- unname(fill[plot_df$name])
  } else {
    fill_color <- rep(fill, length.out = nrow(plot_df))
  }

  ## 5. 绘图并按需旋转标签
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = 3, y = percent)) +
    ggplot2::geom_col(
      width = 1, color = border_color, fill = fill_color,
      linewidth = 0.5, show.legend = FALSE
    ) +
    ggplot2::coord_polar("y", start = start * pi / 180) +
    ggplot2::theme_void() +
    ggplot2::theme(
      aspect.ratio = 1,
      plot.title = ggplot2::element_text(
        color = "black", face = "bold", size = 8 + font_size, hjust = 0.5
      )
    )
  if (!is.null(title)) p <- p + ggplot2::labs(title = as.character(title))

  if (isTRUE(circular_label)) {
    plot_df$label_angle <- if (isTRUE(flip_label)) {
      270 - plot_df$angle
    } else {
      -plot_df$angle
    }
    p <- p + ggplot2::geom_text(
      data = plot_df,
      ggplot2::aes(x = 3.5, y = y_position, label = label, angle = label_angle),
      size = font_size, hjust = 0.5
    )
  } else {
    p <- p + ggplot2::geom_text(
      data = plot_df,
      ggplot2::aes(x = 3.5, y = y_position, label = label),
      size = font_size, hjust = 0.5
    )
  }
  if (isTRUE(hemisphere)) p <- p + ggplot2::lims(x = c(0, 3.5))
  attr(p, "plot_df") <- plot_df
  p
}

#### plot_circlepack ####

#' Plot categorical composition with circle packing
#'
#' Draw a circle-packing composition plot in which each category is represented
#' by a circle with area proportional to its aggregated value.
#'
#' @param data A data frame containing category names and numeric values.
#' @param name_col Name of the category column in `data`.
#' @param value_col Name of the numeric-value column in `data`.
#' @param name_level Optional character vector defining category order. Listed
#'   categories are placed first, followed by categories not present in the
#'   vector. When supplied, this argument takes precedence over `decreasing`.
#' @param decreasing Whether categories are ordered by decreasing value when
#'   `name_level` is `NULL`.
#' @param top_n Optional maximum number of displayed categories, including the
#'   residual category. When more categories are present, the first
#'   `top_n - 1` categories are retained and the remainder are combined.
#' @param other_name Label assigned to categories combined by `top_n`.
#' @param other_last Whether the residual category is placed last.
#' @param add_count Whether category values are included in circle labels.
#' @param add_percent Whether percentages are included in circle labels.
#' @param title Optional plot title.
#' @param fill Fill colors. Supply `"auto"` for the built-in palette, `"hue"`
#'   for a hue palette, one color, an unnamed color vector, or a named vector
#'   containing every displayed category.
#' @param border_color Circle-border color.
#' @param border_size Circle-border line width.
#' @param fill_alpha Circle-fill opacity.
#' @param font_size Label text size.
#' @param font_color Label text color.
#' @param font_face Label font face, such as `"plain"`, `"bold"`, or `"italic"`.
#' @param min_label_value Minimum category value required to display
#'   a label. Categories below this threshold remain visible without labels.
#' @param show_legend Whether the category legend is displayed.
#' @param root_name Name assigned to the internal root node. It must not match
#'   any displayed category name.
#' @param percent_digits Number of decimal places shown in percentage labels.
#' @param seed Random seed used to calculate the circle-packing layout. The
#'   caller's random-number state is restored after layout calculation.
#'
#' @details
#' Missing or empty category names, non-finite values, and non-positive values
#' are removed. Duplicate categories are combined by summing their values.
#'
#' The function constructs a two-level graph containing one root node and one
#' leaf node per displayed category. Only leaf nodes are drawn. A fixed
#' coordinate ratio and symmetric plotting limits preserve circular geometry.
#'
#' Category ordering may affect the spatial arrangement of circles but does
#' not change their values or relative areas. Circle positions do not represent
#' statistical significance, similarity or biological relationships.
#'
#' Percentages use the sum of all valid input values as the denominator,
#' including values combined into the residual category.
#'
#' @return A `ggraph` and `ggplot` object. Its `plot_df` attribute contains the
#'   processed categories in columns `name`, `value`, `percent`, and `label`.
#'
#' @examples
#' example_df <- data.frame(
#'   category = c("Bile acids", "Aromatic compounds", "Amines", "Organic acids"),
#'   n = c(7, 5, 2, 1)
#' )
#'
#' if (
#'   requireNamespace("ggraph", quietly = TRUE) &&
#'     requireNamespace("tidygraph", quietly = TRUE) &&
#'     requireNamespace("withr", quietly = TRUE)
#' ) {
#'   plot_circlepack(
#'     example_df,
#'     name_col = "category",
#'     value_col = "n",
#'     title = "Metabolite composition",
#'     add_percent = TRUE
#'   )
#' }
#'
#' @seealso [ggraph::ggraph()] and [ggraph::create_layout()]
#'
#' @export
plot_circlepack <- function(
  data, name_col = "name", value_col = "n", name_level = NULL,
  decreasing = TRUE, top_n = NULL, other_name = "Other",
  other_last = TRUE, add_count = TRUE, add_percent = FALSE,
  title = NULL, fill = "auto", border_color = "white",
  border_size = 0.5, fill_alpha = 0.9, font_size = 3,
  font_color = "black", font_face = "plain", min_label_value = 1,
  show_legend = FALSE, root_name = "Total", percent_digits = 1,
  seed = 2026
) {
  ## 1. 统一类别和值，并合并重复类别
  plot_df <- as.data.frame(data)
  if (!all(c(name_col, value_col) %in% names(plot_df))) {
    stop("name_col or value_col not found in data.")
  }

  plot_df <- data.frame(
    name = as.character(plot_df[[name_col]]),
    value = as.numeric(plot_df[[value_col]]),
    check.names = FALSE
  ) |>
    dplyr::filter(!is.na(name), name != "", is.finite(value), value > 0) |>
    dplyr::group_by(name) |>
    dplyr::summarise(value = sum(value), .groups = "drop")

  if (!nrow(plot_df)) stop("No valid rows remained after filtering.")

  ## 2. 排序，并把 top_n 之外的类别合并为 Other
  ## 只要指定了 name_level，就不再按照数量升序或降序排列
  ## 排序在 Circle Packing 中主要影响节点进入布局算法的顺序，进而影响圆的位置，不会改变圆的面积
  if (!is.null(name_level)) {
    order_vec <- match(plot_df$name, name_level)
    plot_df <- plot_df[order(order_vec, na.last = TRUE), , drop = FALSE]
  } else if (isTRUE(decreasing)) {
    plot_df <- dplyr::arrange(plot_df, dplyr::desc(value))
  } else {
    plot_df <- dplyr::arrange(plot_df, value)
  }

  if (!is.null(top_n) && nrow(plot_df) > top_n) {
    if (top_n < 2L) stop("top_n should be at least 2 when Other is required.")
    keep_name <- utils::head(plot_df$name, top_n - 1L)
    plot_df$name <- ifelse(plot_df$name %in% keep_name, plot_df$name, other_name)
    plot_df <- plot_df |>
      dplyr::group_by(name) |>
      dplyr::summarise(value = sum(value), .groups = "drop")
  }

  if (isTRUE(other_last) && other_name %in% plot_df$name) {
    plot_df <- dplyr::arrange(plot_df, name == other_name)
  }

  ## 3. 计算比例及标签
  plot_df <- plot_df |>
    dplyr::mutate(percent = value / sum(value))

  if (isTRUE(add_count) && isTRUE(add_percent)) {
    plot_df$label <- paste0(
      plot_df$name, "\n(n=", prettyNum(plot_df$value, big.mark = ","), ")",
      "\n", round(plot_df$percent * 100, percent_digits), "%"
    )
  } else if (isTRUE(add_count)) {
    plot_df$label <- paste0(plot_df$name, "\n(n=", prettyNum(plot_df$value, big.mark = ","), ")")
  } else if (isTRUE(add_percent)) {
    plot_df$label <- paste0(
      plot_df$name, "\n", round(plot_df$percent * 100, percent_digits), "%"
    )
  } else {
    plot_df$label <- plot_df$name
  }

  ## 4. 生成类别颜色
  n_category <- nrow(plot_df)

  if (length(fill) == 1L && fill == "auto") {
    palette <- c(
      "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
      "#a6d854", "#ffd92f", "#8dd3c7", "#fdb462",
      "#80b1d3", "#fccde5", "#d9ef8b", "#fee391"
    )
    fill_color <- rep(palette, length.out = n_category)
  } else if (length(fill) == 1L && fill == "hue") {
    fill_color <- scales::hue_pal()(n_category)
  } else if (!is.null(names(fill)) && all(plot_df$name %in% names(fill))) {
    fill_color <- unname(fill[plot_df$name])
  } else {
    fill_color <- rep(fill, length.out = n_category)
  }

  names(fill_color) <- plot_df$name

  ## 5. 构建节点和边
  ## 构建层级图
  ##   Total (root node)
  ## ├── Bile acids (leaf node)
  ## ├── Aromatic compounds (leaf node)
  ## ├── Nucleotide metabolism (leaf node)
  ## ├── Amines (leaf node)
  ## └── Organic acids (leaf node)
  if (root_name %in% plot_df$name) stop("root_name must not match a category name.")

  nodes <- dplyr::bind_rows(
    tibble::tibble(
      name = root_name, value = 0,
      category = NA_character_, label = NA_character_
    ),
    plot_df |>
      dplyr::transmute(
        name = name, value = as.numeric(value),
        category = name, label = label
      )
  )

  edges <- tibble::tibble(
    from = rep(1L, nrow(plot_df)),
    to = seq.int(2L, nrow(nodes))
  )

  ## 将节点和连接关系转换成 tbl_graph 对象，
  ## 这里的 graph 只是数据结构，还没有计算任何圆的位置
  graph <- tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = TRUE)

  ## 6. 计算 Circle Packing 布局
  ## 计算每个圆的位置和半径
  ## x    | 圆心的横坐标
  ## y    | 圆心的纵坐标
  ## r    | 圆的半径
  ## leaf | 是否为叶节点
  layout <- withr::with_seed(
    seed, ggraph::create_layout(graph, layout = "circlepack", weight = value)
  )

  ## 根据圆的实际边界计算统一坐标范围
  leaf_df <- layout[layout$leaf, ]
  ## 计算所有圆实际占据的空间边界
  xlim <- range(c(leaf_df$x - leaf_df$r, leaf_df$x + leaf_df$r))
  ylim <- range(c(leaf_df$y - leaf_df$r, leaf_df$y + leaf_df$r))
  center_x <- mean(xlim)
  center_y <- mean(ylim)
  ## 计算整个布局的中心，并使用较大的坐标跨度构建正方形显示范围。
  ## 1.05 表示额外保留约5%的半径余量
  radius <- max(diff(xlim), diff(ylim)) / 2 * 1.05

  ## 7. 绘图
  p <- ggraph::ggraph(layout) +
    ggraph::geom_node_circle(
      ggplot2::aes(
        filter = .data[["leaf"]],
        fill = .data[["category"]]
      ),
      color = border_color, linewidth = border_size,
      alpha = fill_alpha, show.legend = show_legend
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(
        filter = .data[["leaf"]] & .data[["value"]] >= min_label_value,
        label = .data[["label"]]
      ),
      size = font_size, color = font_color, fontface = font_face,
      lineheight = 0.9, show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = fill_color, drop = FALSE) +
    ggplot2::coord_fixed(
      ratio = 1, expand = FALSE, clip = "on",
      xlim = center_x + c(-radius, radius),
      ylim = center_y + c(-radius, radius)
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      aspect.ratio = 1,
      legend.position = if (show_legend) "right" else "none",
      plot.title = ggplot2::element_text(
        color = "black", face = "bold", size = 8 + font_size, hjust = 0.5
      ),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )

  if (!is.null(title)) p <- p + ggplot2::labs(title = as.character(title), fill = NULL)

  attr(p, "plot_df") <- plot_df
  p
}
