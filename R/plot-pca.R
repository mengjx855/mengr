#### Jinxin Meng, 20230915, 20250526, v0.4.2 ####

# 20250107: add parameter add_group_label, show_legend, lab_size, show_grid in plot_PCoA function
# 20250115: add sample labels and update group-label parameters.
# 20250317: update function.
# 20250417: plot function pass to plot_dim()
# 20250419: add options sub_sample, sub_group for plot_PCA()
# 20250526: update functions.


#### calcu_PCA ####
# 计算 PCA 坐标和解释度
# profile: 行为 feature，列为 sample 的丰度表
# dim: 输出前几个 PCA 轴
# cumulative_eig: 如果设置，则输出累积解释度达到该阈值所需的轴数
# prefix: 修改 PCA 坐标列名前缀
# add_eig: 是否在列名中加入解释度

#' Calcu PCA utility
#'
#' `calcu_PCA()` provides a reusable mengR workflow with input validation and standardized
#'   output.
#'
#' Chinese summary: 对 feature × sample profile 执行 PCA 并返回坐标和解释方差。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param dim Number of ordination dimensions retained in the returned coordinate table.
#' @param cumulative_eig Cumulative explained-variance threshold used to retain ordination axes.
#' @param prefix Prefix used when naming derived coordinates or labels.
#' @param add_eig Logical control for `add_eig`.
#' @param remove_zero_var Logical control for `remove_zero_var`.
#' @param na_fill Value used to replace missing observations before analysis.
#' @return A result object described in the Details section.
#' @export
calcu_PCA <- function(
    profile, dim = 2, cumulative_eig = NULL, prefix = NULL, add_eig = FALSE,
    remove_zero_var = TRUE, na_fill = 0
  ) {
  
  profile <- data.frame(profile, check.names = FALSE)
  
  if (!is.null(na_fill)) {
    profile[is.na(profile)] <- na_fill
  }
  
  profile <- profile[
    rowSums(profile, na.rm = TRUE) != 0,
    , 
    drop = FALSE
  ]
  
  ## scale. = TRUE 时，零方差 feature 会导致 prcomp 报错
  if (isTRUE(remove_zero_var)) {
    keep <- apply(profile, 1, stats::sd, na.rm = TRUE) > 0
    profile <- profile[keep, , drop = FALSE]
  }
  
  if (nrow(profile) == 0) {
    stop('No valid features remained for PCA.')
  }

  if (ncol(profile) < 2) {
    stop('PCA requires at least two samples.')
  }
  
  PCA <- stats::prcomp(t(profile), scale. = TRUE, center = TRUE)
  PCA_sum <- summary(PCA)
  PCA_eig <- round(PCA_sum$importance[2, ] * 100, 2)
  
  max_dim <- min(ncol(PCA$x), length(PCA_eig))
  
  if (!is.null(cumulative_eig)) {
    dim <- which(cumsum(PCA_eig) >= cumulative_eig)[1]
  }
  
  dim <- min(dim, max_dim)
  
  PCA_points <- data.frame(PCA$x[, seq_len(dim), drop = FALSE]) |>
    tibble::rownames_to_column(var = 'sample')
  
  if (!is.null(prefix)) {
    colnames(PCA_points)[2:(dim + 1)] <- paste0(prefix, '_', seq_len(dim))
  } else {
    colnames(PCA_points)[2:(dim + 1)] <- paste0('PCA', seq_len(dim))
  }
  
  if (isTRUE(add_eig)) {
    colnames(PCA_points)[2:(dim + 1)] <- paste0(
      colnames(PCA_points)[2:(dim + 1)],
      ' (',
      PCA_eig[seq_len(dim)],
      '%)'
    )
  }
  
  out <- list(
    points = PCA_points,
    dim = dim,
    eig = PCA_eig,
    eig_ = paste0(
      colnames(PCA_points)[2:(dim + 1)],
      ' (',
      PCA_eig[seq_len(dim)],
      '%)'
    )
  )
  
  return(out)
}

#### plot_PCA ####
# 绘制 PCA 散点图
# profile: 行为 feature，列为 sample 的丰度表
# group: 样本分组信息表
# sample_col: group 中样本列名，默认 sample
# group_col: group 中分组列名，默认 group
# sub_sample: 指定样本子集
# sub_group: 指定分组子集
# remove_zero_var: 是否删除零方差 feature，避免 prcomp 报错

#' Plot PCA utility
#'
#' `plot_PCA()` provides a reusable mengR workflow with input validation and standardized
#'   output.
#'
#' Chinese summary: 对齐 profile 与分组后执行并绘制 PCA。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param group A sample metadata table containing sample and group columns.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param group_level Optional order of group levels.
#' @param group_color Optional colors aligned to `group_level`.
#' @param sub_sample Optional sample identifiers retained before analysis.
#' @param sub_group Optional group values retained before analysis.
#' @param display_type How samples are displayed, such as points alone or points joined to centroids.
#' @param conf_type Confidence-region geometry: ellipse, encircle, or none.
#' @param ellipse_level Optional order for `ellipse_level`.
#' @param title Optional plot or result title.
#' @param subtitle Optional plot subtitle.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param legend_title Legend title; `NULL` uses a context-dependent default.
#' @param add_group_label Logical control for `add_group_label`.
#' @param add_sample_label Logical control for `add_sample_label`.
#' @param label_size Numeric setting for `label_size`.
#' @param point_size Numeric setting for `point_size`.
#' @param show_legend Logical control for `show_legend`.
#' @param show_grid Logical control for `show_grid`.
#' @param show_line Logical control for `show_line`.
#' @param aspect_ratio Panel aspect ratio passed to `ggplot2::theme()`.
#' @param theme Plot theme preset; supported values are shown in Usage.
#' @param remove_zero_var Logical control for `remove_zero_var`.
#' @param na_fill Value used to replace missing observations before analysis.
#' @param ... Additional arguments passed to the underlying function.
#' @return A plot object; analysis data or models may also be stored as attributes.
#' @export
plot_PCA <- function(
    profile, group, sample_col = 'sample', group_col = 'group', 
    group_level = NULL, group_color = NULL,
    sub_sample = NULL, sub_group = NULL,
    display_type = 'line', conf_type = 'ellipse', 
    ellipse_level = .75, title = NULL, subtitle = NULL, 
    xlab = NULL, ylab = NULL, legend_title = NULL, 
    add_group_label = FALSE, add_sample_label = FALSE, 
    label_size = 3, point_size = 1.5,
    show_legend = TRUE, show_grid = FALSE, show_line = TRUE, 
    aspect_ratio = 3/4, theme = 'default',
    remove_zero_var = TRUE, na_fill = 0, ...
  ) {

  profile <- data.frame(profile, check.names = FALSE)
  group <- data.frame(group, check.names = FALSE)
  
  if (!all(c(sample_col, group_col) %in% colnames(group))) {
    stop('group should contain columns: ', sample_col, ' | ', group_col)
  }
  
  if (!is.null(sub_sample)) {
    group <- dplyr::filter(group, .data[[sample_col]] %in% sub_sample)
  }
  
  if (!is.null(sub_group)) {
    group <- dplyr::filter(group, .data[[group_col]] %in% sub_group)
  }
  
  sample_use <- intersect(group[[sample_col]], colnames(profile))
  
  if (length(sample_use) < 2) {
    stop('PCA requires at least two matched samples.')
  }
  
  group <- group |>
    dplyr::filter(.data[[sample_col]] %in% sample_use)
  
  if (!is.null(na_fill)) {
    profile[is.na(profile)] <- na_fill
  }
  
  profile <- profile[
    rowSums(profile, na.rm = TRUE) != 0,
    ,
    drop = FALSE
  ]
  
  ## scale. = TRUE 时，零方差 feature 会导致 prcomp 报错
  if (isTRUE(remove_zero_var)) {
    keep <- apply(profile, 1, stats::sd, na.rm = TRUE) > 0
    profile <- profile[keep, , drop = FALSE]
  }
  
  if (nrow(profile) == 0) {
    stop('No valid features remained for PCA.')
  }
  
  if (is.null(group_level)) {
    group_level <- as.character(unique(group[[group_col]]))
  }
  
  group_color <- .resolve_group_colors(group_level, group_color)
  
  PCA <- stats::prcomp(t(profile), scale. = TRUE, center = TRUE)
  PCA_sum <- summary(PCA)
  
  PCA_points <- data.frame(PCA$x[, 1:2, drop = FALSE]) |>
    dplyr::rename_with(~ c('X1', 'X2')) |>
    tibble::rownames_to_column('sample')
  
  plot_df <- PCA_points |>
    dplyr::left_join(
      dplyr::select(
        group,
        sample = dplyr::all_of(sample_col),
        group = dplyr::all_of(group_col)
      ),
      by = 'sample'
    ) |>
    dplyr::mutate(group = factor(group, levels = group_level))
  
  if (is.null(xlab)) {
    xlab <- paste0('PC1 (', round(PCA_sum$importance[2, 1] * 100, 2), '%)')
  }
  
  if (is.null(ylab)) {
    ylab <- paste0('PC2 (', round(PCA_sum$importance[2, 2] * 100, 2), '%)')
  }
  
  if (is.null(legend_title)) {
    legend_title <- 'Group'
  }
  
  if (is.null(title)) {
    title <- 'Principal Components Analysis'
  }

  p <- plot_dim(
    data = plot_df, group_level = group_level, group_color = group_color,
    display_type = display_type, conf_type = conf_type, ellipse_level = ellipse_level,
    title = title, subtitle = subtitle, xlab = xlab, ylab = ylab,
    add_group_label = add_group_label, add_sample_label = add_sample_label,
    label_size = label_size, point_size = point_size,
    legend_title = legend_title, show_legend = show_legend,
    show_grid = show_grid, show_line = show_line,
    aspect_ratio = aspect_ratio, theme = theme,
    ...
  )
  
  return(p)
}
