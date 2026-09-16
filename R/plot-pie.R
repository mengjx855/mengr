#### Jin-Xin Meng, jinxmeng@zju.edu.cn, 20231113, 20260916 ####

# 20260916: standardize script metadata, function sections, documentation, and naming style.

#### plot_pie ####

#' Plot a pie chart
#'
#' 从指定名称列和值列绘制饼图或环形图。
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
