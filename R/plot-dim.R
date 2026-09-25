#### Jin-Xin Meng, jinxmeng@zju.edu.cn, 20250417, 20260923 ####

# 20250417: 创建函数 plot_dim, 嵌套所有降维分析的可视化，以统一画图格式；
# 20250617: 升级函数，默认输入的文件第二列和第三列为坐标位置
# 20260526: update some function.
# 20260916: standardize script metadata, function sections, documentation, and naming style.
# 20260923: remove Chinese text from Roxygen documentation.



#### plot_dim ####
# data: 至少包含 sample, X1, X2, group 四列；
#       如果第 2、3 列不是 X1/X2，会自动重命名为 X1/X2。
# display_type: 'line' 表示样本点连到组中心；'point' 表示普通散点图。
# conf_type: 'ellipse', 'encircle', 'none'。
# ellipse_level pass to ggplot::stat_ellipse()
# group_level: 指定分组顺序；默认按数据中出现顺序。
# group_color: 指定分组颜色；默认自动生成。
# expand = 0, s_shape = 1 pass to ggalt::geom_encircle()
# title, subtitle, xlab, ylab, legend_title pass to ggplot::labs()
# show_xxx, show item with feature.
# aspect_ratio and theme pass to ggplot2::theme.

#' Plot Dim utility
#'
#'
#'
#' @param data An input data frame or compatible object.
#' @param group_level Optional order of group levels.
#' @param group_color Optional colors aligned to `group_level`.
#' @param display_type How samples are displayed, such as points alone or points joined to centroids.
#' @param conf_type Confidence-region geometry: ellipse, encircle, or none.
#' @param ellipse_level Optional order for `ellipse_level`.
#' @param expand Numeric expansion applied around groups drawn with `ggforce::geom_mark_ellipse()`.
#' @param s_shape Shape or smoothness parameter used for group-enclosure geometry.
#' @param x_col Name of the `x_col` input column.
#' @param y_col Name of the `y_col` input column.
#' @param title Optional plot or result title.
#' @param subtitle Optional plot subtitle.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param legend_title Legend title; `NULL` uses a context-dependent default.
#' @param add_group_label Whether to label group centroids.
#' @param add_sample_label Whether to label individual samples.
#' @param label_size Text size for sample or group labels.
#' @param point_size Point size used for samples or observations.
#' @param line_width Line width used for plotted paths.
#' @param show_legend Whether to display the plot legend.
#' @param show_grid Whether to draw panel grid lines.
#' @param show_line Whether to draw horizontal and vertical reference lines at zero.
#' @param aspect_ratio Panel aspect ratio passed to `ggplot2::theme()`.
#' @param theme Plot theme preset; supported values are shown in Usage.
#' @param ... Additional arguments passed to `ggplot2::stat_ellipse()` when confidence regions are drawn.
#' @return A ggplot-compatible plot object; computed data or fitted objects are retained as attributes when applicable.
#' @export
plot_dim <- function(
  data, group_level = NULL, group_color = NULL,
  display_type = c("line", "point"),
  conf_type = c("ellipse", "encircle", "none"),
  ellipse_level = .75, expand = 0, s_shape = 1,
  x_col = NULL, y_col = NULL, title = NULL, subtitle = NULL,
  xlab = "PC1", ylab = "PC2", legend_title = NULL,
  add_group_label = FALSE, add_sample_label = FALSE,
  label_size = 2, point_size = 1.5, line_width = .5,
  show_legend = TRUE, show_grid = FALSE, show_line = TRUE,
  aspect_ratio = NULL, theme = c("default", "pubr"), ...
) {
  display_type <- match.arg(display_type)
  conf_type <- match.arg(conf_type)
  theme <- match.arg(theme)

  data <- data.frame(data, check.names = FALSE)

  if (!all(c("sample", "group") %in% colnames(data))) {
    stop("data should contain columns: sample and group.")
  }

  ## 坐标列处理
  ## 默认第 2、3 列为二维坐标；也可以用 x_col/y_col 指定
  if (!is.null(x_col) && !is.null(y_col)) {
    if (!all(c(x_col, y_col) %in% colnames(data))) {
      stop("x_col/y_col not found in data.")
    }

    data <- data |>
      dplyr::rename(
        X1 = dplyr::all_of(x_col),
        X2 = dplyr::all_of(y_col)
      )
  } else {
    if (!all(colnames(data)[2:3] == c("X1", "X2"))) {
      data <- dplyr::rename_with(data, ~ c("X1", "X2"), .cols = 2:3)
    }
  }

  if (conf_type == "encircle" && !requireNamespace("ggalt", quietly = TRUE)) {
    stop("Package 'ggalt' is required when conf_type = 'encircle'.")
  }

  if ((add_group_label || add_sample_label) && !requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required when add_group_label or add_sample_label = TRUE.")
  }

  if (!theme %in% c("default", "pubr")) {
    stop("theme option: default | pubr")
  }

  ## 分组顺序和颜色
  if (is.null(group_level)) {
    group_level <- unique(as.character(data$group))
  }

  data <- dplyr::mutate(data, group = factor(group, group_level))

  if (is.null(group_color)) {
    group_color <- scales::hue_pal()(length(group_level))
    names(group_color) <- group_level
  } else {
    group_color <- rep(group_color, length.out = length(group_level))
    names(group_color) <- group_level
  }

  ## 组中心坐标
  centroid_df <- data |>
    dplyr::group_by(group) |>
    dplyr::summarise(
      X1mean = mean(X1, na.rm = TRUE),
      X2mean = mean(X2, na.rm = TRUE),
      .groups = "drop"
    )

  ## 展示类型为连线：样本点连到组中心
  if (display_type == "line") {
    plot_df <- dplyr::left_join(data, centroid_df, by = "group")

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = X1, y = X2, color = group))

    ## 展示 0-0 基准线
    if (isTRUE(show_line)) {
      p <- p +
        ggplot2::geom_vline(
          xintercept = 0, lty = "longdash",
          linewidth = .4, color = "grey50"
        ) +
        ggplot2::geom_hline(
          yintercept = 0, lty = "longdash",
          linewidth = .4, color = "grey50"
        )
    }

    p <- p +
      ggplot2::geom_segment(
        ggplot2::aes(xend = X1mean, yend = X2mean, color = group),
        linewidth = line_width, show.legend = FALSE
      )

    if (conf_type == "ellipse") {
      p <- p +
        ggplot2::stat_ellipse(
          ggplot2::aes(fill = group, color = group),
          geom = "polygon", alpha = .05, level = ellipse_level,
          linetype = 2, linewidth = line_width, show.legend = FALSE
        )
    }

    if (conf_type == "encircle") {
      p <- p +
        ggalt::geom_encircle(
          ggplot2::aes(fill = group),
          alpha = .05, expand = expand,
          s_shape = s_shape, show.legend = FALSE
        ) +
        ggalt::geom_encircle(
          ggplot2::aes(color = group),
          expand = expand, lty = 2, size = line_width,
          s_shape = s_shape, show.legend = FALSE
        )
    }

    p <- p +
      ggplot2::geom_point(size = point_size, show.legend = FALSE) +
      ggplot2::geom_point(
        data = centroid_df,
        ggplot2::aes(x = X1mean, y = X2mean, color = group),
        size = 2 * point_size, inherit.aes = FALSE
      )
  }

  if (display_type == "point") {
    p <- ggplot2::ggplot(data, ggplot2::aes(x = X1, y = X2, color = group))

    ## 展示 0-0 基准线
    if (isTRUE(show_line)) {
      p <- p +
        ggplot2::geom_vline(
          xintercept = 0, lty = "longdash",
          linewidth = .4, color = "grey50"
        ) +
        ggplot2::geom_hline(
          yintercept = 0, lty = "longdash",
          linewidth = .4, color = "grey50"
        )
    }

    if (conf_type == "ellipse") {
      p <- p +
        ggplot2::stat_ellipse(
          ggplot2::aes(fill = group, color = group),
          geom = "polygon", alpha = .05, level = ellipse_level,
          linetype = 2, linewidth = line_width, show.legend = FALSE
        )
    }

    if (conf_type == "encircle") {
      p <- p +
        ggalt::geom_encircle(
          ggplot2::aes(fill = group),
          alpha = .05, expand = expand,
          s_shape = s_shape, show.legend = FALSE
        ) +
        ggalt::geom_encircle(
          ggplot2::aes(color = group),
          expand = expand, lty = 2, size = line_width,
          s_shape = s_shape, show.legend = FALSE
        )
    }

    p <- p +
      ggplot2::geom_point(size = point_size)
  }

  ## 颜色、标题和坐标轴
  p <- p +
    ggplot2::scale_color_manual(values = group_color, drop = FALSE) +
    ggplot2::scale_fill_manual(values = group_color, drop = FALSE) +
    ggplot2::labs(
      x = xlab, y = ylab,
      title = title, subtitle = subtitle, color = legend_title
    )

  if (isTRUE(add_group_label)) {
    p <- p +
      ggrepel::geom_label_repel(
        data = centroid_df,
        ggplot2::aes(x = X1mean, y = X2mean, label = group, fill = group),
        color = "black", size = label_size, show.legend = FALSE,
        min.segment.length = 10, inherit.aes = FALSE
      )
    show_legend <- FALSE
  }

  ## 添加样本标签
  if (isTRUE(add_sample_label)) {
    p <- p +
      ggrepel::geom_text_repel(
        data = data,
        ggplot2::aes(x = X1, y = X2, label = sample),
        color = "black", size = label_size,
        show.legend = FALSE, inherit.aes = FALSE
      )
  }

  ## 主题
  if (theme == "pubr") {
    p <- p +
      ggpubr::theme_pubr() +
      ggplot2::theme(
        aspect.ratio = aspect_ratio,
        plot.margin = grid::unit(c(2, 2, 2, 2), "mm"),
        plot.title = ggplot2::element_text(hjust = .5, size = 12, face = "bold"),
        legend.position = "right"
      )
  } else {
    p <- p +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.ticks = ggplot2::element_line(linewidth = .5, color = "black"),
        axis.ticks.length = grid::unit(2, "mm"),
        axis.title = ggplot2::element_text(size = 12, color = "black"),
        axis.text = ggplot2::element_text(size = 12, color = "black"),
        axis.line = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = .5, size = 12, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = .5, size = 12, color = "black"),
        plot.margin = grid::unit(c(2, 2, 2, 2), "mm"),
        panel.border = ggplot2::element_rect(linewidth = .5, color = "black", fill = NA),
        panel.background = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        legend.background = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 10, color = "black"),
        legend.title = ggplot2::element_text(size = 10, color = "black"),
        aspect.ratio = aspect_ratio
      )
  }

  ## 是否显示图例
  if (isFALSE(show_legend)) {
    p <- p + ggplot2::guides(color = "none", fill = "none")
  }

  ## 是否显示网格线
  if (isTRUE(show_grid)) {
    p <- p +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(linewidth = .5, color = "grey90"),
        panel.grid.minor = ggplot2::element_blank()
      )
  }

  attr(p, "plot_df") <- data
  attr(p, "centroid_df") <- centroid_df

  return(p)
}
