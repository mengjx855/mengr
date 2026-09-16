#### Jin-Xin Meng, jinxmeng@zju.edu.cn, 20250328, 20260916 ####

# 20250416: update function.
# 20250417: plot function pass to plot_dim()
# 20250419: add options sub_sample, sub_group for plot_pls().
# 20260622: update functional style.
# 20260916: rename PLS/OPLS functions to lowercase and standardize documentation and naming.


#### plot_pls ####
# 绘制 PLS-DA 散点图
# profile: 行为 feature，列为 sample 的丰度表
# group: 样本分组信息表
# sample_col: group 中样本列名，默认 sample
# group_col: group 中分组列名，默认 group
# sub_sample: 指定样本子集
# sub_group: 指定分组子集
# remove_zero_var: 是否删除零方差 feature
# predI: PLS-DA predictive component 数量，默认 3

#' Plot PLS utility
#'
#'
#' 使用 ropls 执行 PLS-DA，并绘制前两个 predictive components。
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
#' @param add_group_label Whether to label group centroids.
#' @param add_sample_label Whether to label individual samples.
#' @param label_size Text size for sample or group labels.
#' @param point_size Point size used for samples or observations.
#' @param show_legend Whether to display the plot legend.
#' @param show_grid Whether to draw panel grid lines.
#' @param show_line Whether to draw horizontal and vertical reference lines at zero.
#' @param aspect_ratio Panel aspect ratio passed to `ggplot2::theme()`.
#' @param theme Plot theme preset; supported values are shown in Usage.
#' @param remove_zero_var Whether to remove zero-variance features before analysis.
#' @param na_fill Value used to replace missing observations before analysis.
#' @param predI Number of predictive components fitted by ropls.
#' @param ... Additional arguments passed to `ropls::opls()` and `plot_dim()`.
#' @return A ggplot-compatible plot object; computed data or fitted objects are retained as attributes when applicable.
#' @export
plot_pls <- function(
  profile, group,
  sample_col = "sample", group_col = "group",
  group_level = NULL, group_color = NULL,
  sub_sample = NULL, sub_group = NULL,
  display_type = "line",
  conf_type = "ellipse", ellipse_level = .75,
  title = NULL, subtitle = "default",
  xlab = NULL, ylab = NULL, legend_title = NULL,
  add_group_label = FALSE, add_sample_label = FALSE,
  label_size = 3, point_size = 1.5,
  show_legend = TRUE, show_grid = FALSE, show_line = TRUE,
  aspect_ratio = 3 / 4, theme = "default",
  remove_zero_var = TRUE, na_fill = 0,
  predI = 3, ...
) {
  profile <- data.frame(profile, check.names = FALSE)
  group <- data.frame(group, check.names = FALSE)

  ## 检查分组列
  if (!all(c(sample_col, group_col) %in% colnames(group))) {
    stop(
      "group should contain columns: ",
      sample_col, " | ", group_col
    )
  }

  ## 过滤样本
  if (!is.null(sub_sample)) {
    group <- dplyr::filter(
      group,
      .data[[sample_col]] %in% sub_sample
    )
  }

  ## 过滤分组
  if (!is.null(sub_group)) {
    group <- dplyr::filter(
      group,
      .data[[group_col]] %in% sub_group
    )
  }

  ## 匹配样本，保持 group 中的样本顺序
  group <- group |>
    dplyr::filter(.data[[sample_col]] %in% colnames(profile))

  sample_use <- as.character(group[[sample_col]])

  if (length(sample_use) < 3) {
    stop("PLS-DA requires at least three matched samples.")
  }

  if (length(unique(group[[group_col]])) < 2) {
    stop("PLS-DA requires at least two groups.")
  }

  ## 避免重复样本名
  if (any(duplicated(sample_use))) {
    stop("Duplicated sample names found in group: ", sample_col)
  }

  ## group level
  if (is.null(group_level)) {
    group_level <- as.character(unique(group[[group_col]]))
  }

  ## group color
  if (is.null(group_color)) {
    default_colors <- c(
      "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
      "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3"
    )

    group_color <- rep(
      default_colors,
      times = ceiling(length(group_level) / length(default_colors))
    )[seq_along(group_level)]

    group_color <- structure(
      group_color,
      names = group_level
    )
  } else {
    if (is.null(names(group_color))) {
      group_color <- structure(
        group_color[seq_along(group_level)],
        names = group_level
      )
    } else {
      group_color <- group_color[group_level]
    }
  }

  ## 整理 profile
  profile <- profile[
    ,
    sample_use,
    drop = FALSE
  ]

  if (!is.null(na_fill)) {
    profile[is.na(profile)] <- na_fill
  }

  ## 删除全 0 feature
  profile <- profile[
    rowSums(profile, na.rm = TRUE) != 0, ,
    drop = FALSE
  ]

  ## 删除零方差 feature
  if (isTRUE(remove_zero_var)) {
    keep <- apply(profile, 1, stats::sd, na.rm = TRUE) > 0

    profile <- profile[
      keep, ,
      drop = FALSE
    ]
  }

  if (nrow(profile) == 0) {
    stop("No valid features remained for PLS-DA.")
  }

  ## ropls 输入：sample × feature
  pls <- suppressWarnings(
    ropls::opls(
      x = data.frame(t(profile), check.names = FALSE),
      y = as.vector(group[[group_col]]),
      orthoI = 0,
      predI = predI,
      info.txtC = "none",
      fig.pdfC = "none"
    )
  )

  if (nrow(pls@scoreMN) == 0) {
    stop("No PLS-DA score matrix was generated.")
  }

  if (ncol(pls@scoreMN) < 2) {
    stop("PLS-DA generated fewer than two components.")
  }

  ## 提取坐标
  pls_points <- data.frame(
    pls@scoreMN[, 1:2, drop = FALSE],
    check.names = FALSE
  ) |>
    dplyr::rename_with(~ c("X1", "X2")) |>
    tibble::rownames_to_column("sample")

  plot_df <- pls_points |>
    dplyr::left_join(
      dplyr::select(
        group,
        sample = dplyr::all_of(sample_col),
        group = dplyr::all_of(group_col)
      ),
      by = "sample"
    ) |>
    dplyr::mutate(
      group = factor(group, levels = group_level)
    )

  ## x/y lab
  if (is.null(xlab)) {
    x_var <- suppressWarnings(as.numeric(pls@modelDF[1, 1]) * 100)

    if (is.finite(x_var)) {
      xlab <- paste0("PC1 (", round(x_var, 2), "%)")
    } else {
      xlab <- "PC1"
    }
  }

  if (is.null(ylab)) {
    y_var <- suppressWarnings(as.numeric(pls@modelDF[2, 1]) * 100)

    if (is.finite(y_var)) {
      ylab <- paste0("PC2 (", round(y_var, 2), "%)")
    } else {
      ylab <- "PC2"
    }
  }

  if (is.null(legend_title)) {
    legend_title <- "Group"
  }

  ## subtitle
  if (!is.null(subtitle)) {
    if (identical(subtitle, "default")) {
      subtitle <- substitute(
        R^2 * X == a ~ ~ R^2 * Y == b ~ ~ Q^2 == c ~ ~ RMSEE == d,
        list(
          a = round(pls@summaryDF[1, 1], 3),
          b = round(pls@summaryDF[1, 2], 3),
          c = round(pls@summaryDF[1, 3], 3),
          d = round(pls@summaryDF[1, 4], 3)
        )
      )
    }
  }

  if (is.null(title)) {
    title <- "Partial least squares discriminant analysis"
  }

  p <- plot_dim(
    data = plot_df,
    group_level = group_level,
    group_color = group_color,
    display_type = display_type,
    conf_type = conf_type,
    ellipse_level = ellipse_level,
    title = title,
    subtitle = subtitle,
    xlab = xlab,
    ylab = ylab,
    legend_title = legend_title,
    add_group_label = add_group_label,
    add_sample_label = add_sample_label,
    label_size = label_size,
    point_size = point_size,
    show_legend = show_legend,
    show_grid = show_grid,
    show_line = show_line,
    aspect_ratio = aspect_ratio,
    theme = theme,
    ...
  )

  attr(p, "model") <- pls
  attr(p, "plot_df") <- plot_df

  return(p)
}

#### plot_opls ####
# 绘制 OPLS-DA 散点图
# profile: 行为 feature，列为 sample 的丰度表
# group: 样本分组信息表
# OPLS-DA 只适用于二分类
# top_frac: 选择方差最大的 top 比例 feature 进行分析，默认 0.2

#' Plot OPLS utility
#'
#'
#' 对二分类数据执行 OPLS-DA 并绘制 predictive/orthogonal components。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param group A sample metadata table containing sample and group columns.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param group_level Optional order of group levels.
#' @param group_color Optional colors aligned to `group_level`.
#' @param sub_sample Optional sample identifiers retained before analysis.
#' @param sub_group Optional group values retained before analysis.
#' @param top_frac Fraction of highest-ranking features retained.
#' @param display_type How samples are displayed, such as points alone or points joined to centroids.
#' @param conf_type Confidence-region geometry: ellipse, encircle, or none.
#' @param ellipse_level Optional order for `ellipse_level`.
#' @param title Optional plot or result title.
#' @param subtitle Optional plot subtitle.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param legend_title Legend title; `NULL` uses a context-dependent default.
#' @param add_group_label Whether to label group centroids.
#' @param add_sample_label Whether to label individual samples.
#' @param label_size Text size for sample or group labels.
#' @param point_size Point size used for samples or observations.
#' @param show_legend Whether to display the plot legend.
#' @param show_grid Whether to draw panel grid lines.
#' @param show_line Whether to draw horizontal and vertical reference lines at zero.
#' @param aspect_ratio Panel aspect ratio passed to `ggplot2::theme()`.
#' @param theme Plot theme preset; supported values are shown in Usage.
#' @param remove_zero_var Whether to remove zero-variance features before analysis.
#' @param na_fill Value used to replace missing observations before analysis.
#' @param predI Number of predictive components fitted by ropls.
#' @param orthoI Number of orthogonal components fitted by OPLS-DA.
#' @param ... Additional arguments passed to `ropls::opls()` and `plot_dim()`.
#' @return A ggplot-compatible plot object; computed data or fitted objects are retained as attributes when applicable.
#' @export
plot_opls <- function(
  profile, group,
  sample_col = "sample", group_col = "group",
  group_level = NULL, group_color = NULL,
  sub_sample = NULL, sub_group = NULL,
  top_frac = .2,
  display_type = "line", conf_type = "ellipse",
  ellipse_level = .75,
  title = NULL, subtitle = "default",
  xlab = NULL, ylab = NULL, legend_title = NULL,
  add_group_label = FALSE, add_sample_label = FALSE,
  label_size = 3, point_size = 1.5,
  show_legend = TRUE, show_grid = FALSE, show_line = TRUE,
  aspect_ratio = 3 / 4, theme = "default",
  remove_zero_var = TRUE, na_fill = 0,
  predI = 1, orthoI = NA,
  ...
) {
  profile <- data.frame(profile, check.names = FALSE)
  group <- data.frame(group, check.names = FALSE)

  ## 检查分组列
  if (!all(c(sample_col, group_col) %in% colnames(group))) {
    stop(
      "group should contain columns: ",
      sample_col, " | ", group_col
    )
  }

  ## 过滤样本
  if (!is.null(sub_sample)) {
    group <- dplyr::filter(
      group,
      .data[[sample_col]] %in% sub_sample
    )
  }

  ## 过滤分组
  if (!is.null(sub_group)) {
    group <- dplyr::filter(
      group,
      .data[[group_col]] %in% sub_group
    )
  }

  ## 匹配样本，保持 group 中的样本顺序
  group <- group |>
    dplyr::filter(.data[[sample_col]] %in% colnames(profile))

  sample_use <- as.character(group[[sample_col]])

  if (length(sample_use) < 4) {
    stop("OPLS-DA requires at least four matched samples.")
  }

  if (any(duplicated(sample_use))) {
    stop("Duplicated sample names found in group: ", sample_col)
  }

  ## OPLS-DA 二分类检查
  group_n <- length(unique(as.character(group[[group_col]])))

  if (group_n != 2) {
    stop(
      "OPLS-DA only supports binary classification. ",
      "Current number of groups: ", group_n,
      ". Use plot_pls() for multiple classes."
    )
  }

  ## group level
  if (is.null(group_level)) {
    group_level <- as.character(unique(group[[group_col]]))
  }

  if (length(group_level) != 2) {
    stop("group_level should contain exactly two groups for OPLS-DA.")
  }

  ## group color
  if (is.null(group_color)) {
    default_colors <- c(
      "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
      "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3"
    )

    group_color <- rep(
      default_colors,
      times = ceiling(length(group_level) / length(default_colors))
    )[seq_along(group_level)]

    group_color <- structure(
      group_color,
      names = group_level
    )
  } else {
    if (is.null(names(group_color))) {
      group_color <- structure(
        group_color[seq_along(group_level)],
        names = group_level
      )
    } else {
      group_color <- group_color[group_level]
    }
  }

  ## 整理 profile
  profile <- profile[
    ,
    sample_use,
    drop = FALSE
  ]

  if (!is.null(na_fill)) {
    profile[is.na(profile)] <- na_fill
  }

  ## 删除全 0 feature
  profile <- profile[
    rowSums(profile, na.rm = TRUE) != 0, ,
    drop = FALSE
  ]

  ## 删除零方差 feature
  if (isTRUE(remove_zero_var)) {
    keep <- apply(profile, 1, stats::sd, na.rm = TRUE) > 0

    profile <- profile[
      keep, ,
      drop = FALSE
    ]
  }

  if (nrow(profile) == 0) {
    stop("No valid features remained before top_frac filtering.")
  }

  ## top_frac 选择方差最大的 top 比例 feature
  if (!is.null(top_frac)) {
    if (!is.numeric(top_frac) || length(top_frac) != 1) {
      stop("top_frac should be a numeric value between 0 and 1.")
    }

    if (top_frac <= 0 || top_frac > 1) {
      stop("top_frac should be > 0 and <= 1.")
    }

    feature_sd <- apply(profile, 1, stats::sd, na.rm = TRUE)

    n_top <- max(
      2,
      ceiling(length(feature_sd) * top_frac)
    )

    feature_keep <- names(
      sort(feature_sd, decreasing = TRUE)
    )[seq_len(n_top)]

    profile <- profile[
      feature_keep, ,
      drop = FALSE
    ]
  }

  if (nrow(profile) < 2) {
    stop("OPLS-DA requires at least two valid features.")
  }

  ## ropls 输入：sample × feature
  opls <- suppressWarnings(
    ropls::opls(
      x = data.frame(t(profile), check.names = FALSE),
      y = as.character(group[[group_col]]),
      orthoI = orthoI,
      predI = predI,
      info.txtC = "none",
      fig.pdfC = "none"
    )
  )

  if (nrow(opls@modelDF) == 0) {
    stop(
      "No OPLS-DA model was built. ",
      "The first predictive component may be not significant."
    )
  }

  if (ncol(opls@scoreMN) < 1) {
    stop("No predictive score was generated by OPLS-DA.")
  }

  if (ncol(opls@orthoScoreMN) < 1) {
    stop("No orthogonal score was generated by OPLS-DA.")
  }

  ## 提取 predictive score 和 first orthogonal score
  opls_points <- data.frame(
    X1 = opls@scoreMN[, 1],
    X2 = opls@orthoScoreMN[, 1],
    check.names = FALSE
  ) |>
    tibble::rownames_to_column("sample")

  plot_df <- opls_points |>
    dplyr::left_join(
      dplyr::select(
        group,
        sample = dplyr::all_of(sample_col),
        group = dplyr::all_of(group_col)
      ),
      by = "sample"
    ) |>
    dplyr::mutate(
      group = factor(group, levels = group_level)
    )

  ## x/y lab
  if (is.null(xlab)) {
    x_var <- suppressWarnings(as.numeric(opls@modelDF[1, 1]) * 100)

    if (is.finite(x_var)) {
      xlab <- paste0("Predictive component (", round(x_var, 2), "%)")
    } else {
      xlab <- "Predictive component"
    }
  }

  if (is.null(ylab)) {
    y_var <- suppressWarnings(as.numeric(opls@modelDF[2, 1]) * 100)

    if (is.finite(y_var)) {
      ylab <- paste0("Orthogonal component (", round(y_var, 2), "%)")
    } else {
      ylab <- "Orthogonal component"
    }
  }

  if (is.null(legend_title)) {
    legend_title <- "Group"
  }

  ## subtitle
  if (!is.null(subtitle)) {
    if (identical(subtitle, "default")) {
      subtitle <- substitute(
        R^2 * X == a ~ ~ R^2 * Y == b ~ ~ Q^2 == c ~ ~ RMSEE == d,
        list(
          a = round(opls@summaryDF[1, 1], 3),
          b = round(opls@summaryDF[1, 2], 3),
          c = round(opls@summaryDF[1, 3], 3),
          d = round(opls@summaryDF[1, 4], 3)
        )
      )
    }
  }

  if (is.null(title)) {
    title <- "Orthogonal partial least squares discriminant analysis"
  }

  p <- plot_dim(
    data = plot_df,
    group_level = group_level,
    group_color = group_color,
    display_type = display_type,
    conf_type = conf_type,
    ellipse_level = ellipse_level,
    title = title,
    subtitle = subtitle,
    xlab = xlab,
    ylab = ylab,
    legend_title = legend_title,
    add_group_label = add_group_label,
    add_sample_label = add_sample_label,
    label_size = label_size,
    point_size = point_size,
    show_legend = show_legend,
    show_grid = show_grid,
    show_line = show_line,
    aspect_ratio = aspect_ratio,
    theme = theme,
    ...
  )

  attr(p, "model") <- opls
  attr(p, "plot_df") <- plot_df
  attr(p, "top_frac") <- top_frac

  return(p)
}
