#### ggradar pair plot ####
# data: 需要包含 proj、cf_col、random_col 和 p_col
# cf_col/random_col: 两组要比较的数值列
# p_col: 用于标记显著性的 p 值列
# star_radius: 星号离中心的距离；如果星号太靠内或太靠外，可以调这个值
#' Plot Ggradar Pair utility
#'
#' `plot_ggradar_pair()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 将成对比较数据整理为 ggradar 所需格式并绘制雷达图。
#'
#' @param data An input data frame or compatible object.
#' @param cf_col Name of the `cf_col` input column.
#' @param random_col Name of the `random_col` input column.
#' @param p_col Name of the `p_col` input column.
#' @param title Optional plot or result title.
#' @param proj_col Name of the `proj_col` input column.
#' @param cf_name Display or identifier name for cf name.
#' @param random_name Display or identifier name for random name.
#' @param colors Color specification for `colors`.
#' @param grid_max Upper limit of the radar grid; `NULL` derives it from the data.
#' @param star_radius Inner radius used to position the radar-chart star polygons.
#' @param axis.label.size Numeric setting for `axis.label.size`.
#' @param grid.label.size Numeric setting for `grid.label.size`.
#' @param legend.position Legend position passed to `ggplot2::theme()`.
#' @return A plot object; analysis data or models may also be stored as attributes.
#' @export
plot_ggradar_pair <- function(data, cf_col, random_col, p_col,
                              title = NULL,
                              proj_col = "proj",
                              cf_name = "CF",
                              random_name = "Random",
                              colors = c(CF = "#E64B35", Random = "#4DBBD5"),
                              grid_max = NULL,
                              star_radius = 1.16,
                              axis.label.size = 3,
                              grid.label.size = 3,
                              legend.position = "bottom") {
  data <- data.frame(data, check.names = FALSE)
  data[[proj_col]] <- as.character(data[[proj_col]])

  ## 1. 整理成 ggradar 需要的格式：
  ##    第一列是 group，后面每一列是一个雷达轴
  radar_data <- data |>
    dplyr::select(
      proj = dplyr::all_of(proj_col),
      CF = dplyr::all_of(cf_col),
      Random = dplyr::all_of(random_col)
    ) |>
    tidyr::pivot_longer(
      cols = c("CF", "Random"),
      names_to = "group",
      values_to = "value"
    ) |>
    tidyr::pivot_wider(
      names_from = proj,
      values_from = value
    ) |>
    dplyr::mutate(group = factor(group, levels = c("CF", "Random"))) |>
    data.frame(check.names = FALSE)

  ## 2. 自动设置雷达图最大刻度
  value_mat <- as.matrix(radar_data[, -1])

  if (is.null(grid_max)) {
    grid_max <- max(value_mat, na.rm = TRUE) * 1.15
    grid_max <- signif(grid_max, 2)
  }

  grid_mid <- grid_max / 2

  ## 3. p 值星号位置
  axis_names <- colnames(radar_data)[-1]

  p_data <- data |>
    dplyr::transmute(
      proj = .data[[proj_col]],
      plab = .add_plab(.data[[p_col]], format = 2)
    ) |>
    dplyr::filter(plab != "") |>
    dplyr::mutate(
      idx = match(proj, axis_names),
      n_axis = length(axis_names),
      angle = pi / 2 - 2 * pi * (idx - 1) / n_axis,
      x = star_radius * cos(angle),
      y = star_radius * sin(angle)
    )

  ## 4. ggradar 作图
  p <- ggradar::ggradar(
    radar_data,
    grid.min = 0,
    grid.mid = grid_mid,
    grid.max = grid_max,
    values.radar = c("0", signif(grid_mid, 2), signif(grid_max, 2)),
    group.colours = colors,
    group.line.width = 0.9,
    group.point.size = 2.2,
    background.circle.colour = "white",
    gridline.min.colour = "grey88",
    gridline.mid.colour = "grey82",
    gridline.max.colour = "grey70",
    axis.label.size = axis.label.size,
    grid.label.size = grid.label.size,
    legend.position = legend.position
  ) +
    ggplot2::labs(title = title) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 13),
      legend.text = ggplot2::element_text(size = 10)
    )

  ## 5. 加显著性星号
  if (nrow(p_data) > 0) {
    p <- p +
      ggplot2::geom_text(
        data = p_data,
        ggplot2::aes(x = x, y = y, label = plab),
        inherit.aes = FALSE,
        size = 5,
        fontface = "bold",
        color = "black"
      )
  }

  return(p)
}

#### radar plot ####
#' Plot Pair Radar utility
#'
#' `plot_pair_radar()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 使用 ggplot2 绘制成对比较的极坐标雷达图。
#'
#' @param data An input data frame or compatible object.
#' @param cf_col Name of the `cf_col` input column.
#' @param random_col Name of the `random_col` input column.
#' @param p_col Name of the `p_col` input column.
#' @param title Optional plot or result title.
#' @param cf_name Display or identifier name for cf name.
#' @param random_name Display or identifier name for random name.
#' @param colors Color specification for `colors`.
#' @param fill_alpha Numeric setting for `fill_alpha`.
#' @param line_width Numeric setting for `line_width`.
#' @param point_size Numeric setting for `point_size`.
#' @param label_size Numeric setting for `label_size`.
#' @param star_size Numeric setting for `star_size`.
#' @param grid_n Number of grid intervals.
#' @return A plot object; analysis data or models may also be stored as attributes.
#' @export
plot_pair_radar <- function(data, cf_col, random_col, p_col,
                            title = NULL,
                            cf_name = "CF",
                            random_name = "Random",
                            colors = c(CF = "#D73027", Random = "#4575B4"),
                            fill_alpha = 0.16,
                            line_width = 0.75,
                            point_size = 2.2,
                            label_size = 3.4,
                            star_size = 5,
                            grid_n = 4) {
  data <- data.frame(data, check.names = FALSE)
  data$proj <- factor(data$proj, levels = data$proj)

  plot_df <- data |>
    dplyr::select(
      proj,
      CF = dplyr::all_of(cf_col),
      Random = dplyr::all_of(random_col)
    ) |>
    tidyr::pivot_longer(
      cols = c("CF", "Random"),
      names_to = "type",
      values_to = "value"
    ) |>
    dplyr::mutate(
      type = factor(type, levels = c("CF", "Random"))
    )

  max_y <- max(plot_df$value, na.rm = TRUE)
  min_y <- 0
  y_breaks <- pretty(c(min_y, max_y), n = grid_n)
  y_max <- max(y_breaks) * 1.18

  p_data <- data |>
    dplyr::transmute(
      proj,
      pval = .data[[p_col]],
      plab = .add_plab(.data[[p_col]], format = 2),
      y = y_max * 0.96
    )

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = proj, y = value, group = type, color = type, fill = type)
  ) +
    ggplot2::geom_polygon(alpha = fill_alpha, linewidth = line_width) +
    ggplot2::geom_line(linewidth = line_width) +
    ggplot2::geom_point(size = point_size, stroke = 0.3) +
    ggplot2::geom_text(
      data = p_data,
      ggplot2::aes(x = proj, y = y, label = plab),
      inherit.aes = FALSE,
      size = star_size,
      fontface = "bold",
      color = "black"
    ) +
    ggplot2::coord_polar(clip = "off") +
    ggplot2::scale_color_manual(
      values = colors,
      labels = c(CF = cf_name, Random = random_name)
    ) +
    ggplot2::scale_fill_manual(
      values = colors,
      labels = c(CF = cf_name, Random = random_name)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, y_max),
      breaks = y_breaks,
      labels = round(y_breaks, 2)
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      color = NULL,
      fill = NULL,
      title = title
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(linewidth = 0.35, color = "grey82"),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = label_size * 3, color = "black"),
      axis.text.y = ggplot2::element_text(size = 8, color = "grey35"),
      axis.title = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 10, color = "black"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 13),
      plot.margin = ggplot2::margin(10, 12, 10, 12)
    )
}
