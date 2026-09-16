#### Jin-Xin Meng, jinxmeng@zju.edu.cn, 20221001, 20260916 ####

# 20221001: 添加选择不同距离尺度的参数 dist_method
# 20231204: update function, check_file_name was deprecated.
# 20250107: add parameter add_lab_to_plot, show_legend, lab_size, show_grid in plot_pcoa function
# 20250417: plot function pass to plot_dim()
# 20250418: add plot_pcoa_box().
# 20250617: add new function calcu_adonis_r2().
# 20250701: update function calcu_betadisper().
# 20251111: 计算距离前进行hellinger转换
# 20251116: labels() 获取矩阵行列名
# 20260527: update some functions.
# 20260819: keep pcoa-related functions, others had been moved to diversity.R
# 20260916: standardize documentation and rename internal data-frame variables to the `*_df` style.


#### calcu_pcoa ####
# profile: 行为 feature，列为 sample 的丰度表
# distance: 已经计算好的 dist 对象；如果提供 distance，则忽略 profile
# group: 样本分组表；adonis2 = TRUE 时需要
# sample_col/group_col: group 中样本列和分组列
# dim: 输出前几个 pcoa 轴
# cumulative_eig: 累积解释度阈值；如果设置，则自动决定输出轴数
# dist_method: 距离方法，传给 calcu_distance()
# prefix: pcoa 轴名前缀
# adonis2: 是否计算 PERMANOVA

#' Calcu pcoa utility
#'
#'
#' 从 profile 或距离对象计算 pcoa 坐标和特征值。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param distance A precomputed distance object; when supplied, it takes precedence over `profile`.
#' @param group A sample metadata table containing sample and group columns.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param dim Number of ordination dimensions retained in the returned coordinate table.
#' @param cumulative_eig Cumulative explained-variance threshold used to retain ordination axes.
#' @param dist_method Distance method; available values are validated with `match.arg()`.
#' @param prefix Prefix used when naming derived coordinates or labels.
#' @param adonis2 Whether to run PERMANOVA with `vegan::adonis2()` and attach its result.
#' @param permutations Number of permutations used by the significance test.
#' @param ... Additional arguments passed to distance calculation or `vegan::adonis2()`.
#' @return A list containing the PCoA object, sample coordinates, eigenvalues, and explained variance.
#' @importFrom rlang .data
#' @export
calcu_pcoa <- function(
  profile = NULL, distance = NULL, group = NULL,
  sample_col = "sample", group_col = "group", dim = 2,
  cumulative_eig = NULL,
  dist_method = c("bray", "jaccard", "euclidean", "manhattan", "unifrac"),
  prefix = NULL, adonis2 = FALSE, permutations = 999, ...
) {
  ## 1. 检查距离方法并计算或整理距离矩阵
  dist_method <- .match_distance_method(dist_method)
  ## 如果没有提供 distance，就用 calcu_distance() 计算
  if (is.null(distance)) {
    if (is.null(profile)) {
      stop("Need either profile or distance.")
    }

    distance <- calcu_distance(
      profile = profile,
      dist_method = dist_method,
      ...
    )
  }

  distance <- stats::as.dist(distance)

  ## 2. pcoa
  if (is.null(cumulative_eig)) {
    pcoa <- stats::cmdscale(distance, k = dim, eig = TRUE)
  } else {
    k_all <- attr(distance, "Size") - 1
    pcoa <- suppressWarnings(stats::cmdscale(distance, k = k_all, eig = TRUE))

    eig_tmp <- pcoa$eig
    eig_pct_tmp <- round(eig_tmp / sum(eig_tmp[eig_tmp > 0]) * 100, 2)
    dim <- which(cumsum(eig_pct_tmp) >= cumulative_eig)[1]
  }

  ## 3. eigenvalue percentage
  pcoa_eig <- round(pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]) * 100, 2)

  dim <- min(dim, ncol(pcoa$points))

  pcoa_points <- data.frame(
    pcoa$points[, seq_len(dim), drop = FALSE],
    check.names = FALSE
  ) |>
    tibble::rownames_to_column("sample")

  if (is.null(prefix)) {
    axis_names <- paste0("PC", seq_len(dim))
  } else {
    axis_names <- paste0(prefix, seq_len(dim))
  }

  colnames(pcoa_points)[2:(dim + 1)] <- axis_names

  result <- list(
    point = pcoa_points,
    dist = distance,
    eigenvalue = pcoa_eig,
    eig_lab = paste0(axis_names, "(", pcoa_eig[seq_len(dim)], "%)")
  )

  ## 4. adonis2
  if (isTRUE(adonis2)) {
    if (is.null(group)) {
      stop("Need group information when adonis2 = TRUE.")
    }

    group <- data.frame(group, check.names = FALSE)

    if (!all(c(sample_col, group_col) %in% colnames(group))) {
      stop("group should contain columns: ", sample_col, " | ", group_col)
    }

    sample_order <- labels(distance)

    meta <- group |>
      dplyr::select(
        sample = dplyr::all_of(sample_col),
        group = dplyr::all_of(group_col)
      ) |>
      dplyr::filter(sample %in% sample_order) |>
      dplyr::arrange(match(sample, sample_order))

    if (!identical(meta$sample, sample_order)) {
      stop("Some samples in distance object are missing in group.")
    }

    meta$group <- factor(meta$group)

    adonis <- vegan::adonis2(
      distance ~ group,
      data = meta, permutations = permutations
    )

    r2 <- round(adonis$R2[1], 4)
    r2adj <- calcu_adjusted_r2(adonis)
    pval <- adonis$`Pr(>F)`[1]

    p_lab <- dplyr::case_when(
      pval < 0.001 ~ "***",
      pval < 0.01 ~ "**",
      pval < 0.05 ~ "*",
      TRUE ~ "ns"
    )

    label <- paste0("R2 = ", r2, ", p ", p_lab)

    result[["group"]] <- meta
    result[["adonis2"]] <- adonis
    result[["adonis2_r2"]] <- r2
    result[["adonis2_r2adj"]] <- round(r2adj, 4)
    result[["adonis2_pval"]] <- pval
    result[["adonis2_lab"]] <- label
  }

  result
}

#### plot_pcoa ####
# 绘制 pcoa 图
# profile: 行为 feature，列为 sample 的丰度表
# distance: 已经计算好的距离矩阵或 dist 对象；如果提供 distance，则优先使用 distance
# group: 样本分组信息表
# sample_col/group_col: group 中样本列和分组列
# sub_sample: 指定样本子集
# sub_group: 指定分组子集
# dist_method: 距离方法，默认 bray
# transform: 计算距离前是否转换，默认 hellinger；NULL 表示不转换
# adonis2: 是否添加 PERMANOVA 结果作为 subtitle

#' Plot pcoa utility
#'
#'
#' 计算并绘制 pcoa，可附加 PERMANOVA 统计结果。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param group A sample metadata table containing sample and group columns.
#' @param distance A precomputed distance object; when supplied, it takes precedence over `profile`.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param group_level Optional order of group levels.
#' @param group_color Optional colors aligned to `group_level`.
#' @param sub_sample Optional sample identifiers retained before analysis.
#' @param sub_group Optional group values retained before analysis.
#' @param dist_method Distance method; available values are validated with `match.arg()`.
#' @param transform Optional transformation applied before analysis.
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
#' @param adonis2 Whether to run PERMANOVA with `vegan::adonis2()` and attach its result.
#' @param permutations Number of permutations used by the significance test.
#' @param ... Additional arguments passed to `vegan::vegdist()` and `plot_dim()`.
#' @return A ggplot-compatible plot object; computed data or fitted objects are retained as attributes when applicable.
#' @export
plot_pcoa <- function(
  profile = NULL, group, distance = NULL,
  sample_col = "sample", group_col = "group",
  group_level = NULL, group_color = NULL,
  sub_sample = NULL, sub_group = NULL,
  dist_method = c("bray", "jaccard", "euclidean", "manhattan", "unifrac"),
  transform = c("hellinger", "total", "pa", "clr"),
  display_type = c("line", "point"),
  conf_type = c("ellipse", "encircle", "none"),
  ellipse_level = .75, title = NULL, subtitle = NULL,
  xlab = NULL, ylab = NULL, legend_title = NULL,
  add_group_label = FALSE, add_sample_label = FALSE,
  label_size = 3, point_size = 1.5,
  show_legend = TRUE, show_grid = FALSE, show_line = TRUE,
  aspect_ratio = 3 / 4, theme = c("default", "pubr"),
  adonis2 = TRUE, permutations = 999, ...
) {
  dist_method <- .match_distance_method(dist_method)
  transform <- .match_transform_method(transform)
  display_type <- match.arg(display_type)
  conf_type <- match.arg(conf_type)
  theme <- match.arg(theme)

  group <- data.frame(group, check.names = FALSE)

  if (!all(c(sample_col, group_col) %in% colnames(group))) {
    stop("group should contain columns: ", sample_col, " | ", group_col)
  }

  ## 1. 根据样本或分组筛选 metadata
  if (!is.null(sub_sample)) {
    group <- dplyr::filter(group, .data[[sample_col]] %in% sub_sample)
  }

  if (!is.null(sub_group)) {
    group <- dplyr::filter(group, .data[[group_col]] %in% sub_group)
  }

  ## 2. 计算或整理距离矩阵
  if (is.null(distance)) {
    if (is.null(profile)) {
      stop("Need either profile or distance.")
    }

    profile <- data.frame(profile, check.names = FALSE)
    profile <- as.matrix(profile)
    suppressWarnings(storage.mode(profile) <- "numeric")
    profile[!is.finite(profile)] <- 0

    sample_use <- intersect(group[[sample_col]], colnames(profile))

    if (length(sample_use) < 2) {
      stop("NMDS requires at least two matched samples.")
    }

    group <- group |>
      dplyr::filter(.data[[sample_col]] %in% sample_use)

    profile <- profile[, group[[sample_col]], drop = FALSE]

    ## 删除全 0 feature 和全 0 sample
    profile <- profile[
      rowSums(profile, na.rm = TRUE) != 0,
      colSums(profile, na.rm = TRUE) != 0,
      drop = FALSE
    ]

    if (nrow(profile) == 0) {
      stop("No valid features remained for NMDS.")
    }

    if (ncol(profile) < 2) {
      stop("NMDS requires at least two valid samples.")
    }

    profile_t <- t(profile)

    if (!is.null(transform)) {
      profile_t <- vegan::decostand(profile_t, method = transform)
    }

    profile_t[!is.finite(profile_t)] <- 0

    ## 转换后再次删除全 0 sample，避免 jaccard / bray 出现 empty rows
    profile_t <- profile_t[rowSums(profile_t, na.rm = TRUE) != 0, , drop = FALSE]

    if (nrow(profile_t) < 2) {
      stop("NMDS requires at least two valid samples after transformation.")
    }

    distance <- vegan::vegdist(profile_t, method = dist_method, na.rm = TRUE, ...)
  } else {
    distance <- stats::as.dist(distance)

    distance_sample <- labels(distance)
    sample_use <- distance_sample[distance_sample %in% group[[sample_col]]]

    if (length(sample_use) < 2) {
      stop("NMDS requires at least two matched samples.")
    }

    group <- group |>
      dplyr::filter(.data[[sample_col]] %in% sample_use)

    ## 只有样本确实被筛选时，才把 dist 转成 matrix 后 subset
    ## 如果 distance 和 group 本来完全匹配，就保留原 dist，避免 as.matrix() 的巨大开销
    if (!identical(sample_use, distance_sample)) {
      distance <- stats::as.dist(as.matrix(distance)[sample_use, sample_use])
    }
  }

  ## 3. 对齐 group 和 distance 中的样本顺序
  sample_order <- labels(distance)

  group <- group |>
    dplyr::filter(.data[[sample_col]] %in% sample_order) |>
    dplyr::arrange(match(.data[[sample_col]], sample_order))

  if (!identical(group[[sample_col]], sample_order)) {
    stop("Some samples in distance object are missing in group.")
  }

  ## 4. 设置分组顺序和颜色
  if (is.null(group_level)) {
    group_level <- as.character(unique(group[[group_col]]))
  }

  group_level <- group_level[group_level %in% unique(group[[group_col]])]

  if (length(group_level) == 0) {
    stop("No valid group remained for plotting.")
  }

  group_color <- .resolve_group_colors(group_level, group_color)

  ## 5. pcoa 降维
  pcoa <- stats::cmdscale(distance, k = 2, eig = TRUE)

  pcoa_points <- data.frame(
    pcoa$points[, 1:2, drop = FALSE],
    check.names = FALSE
  ) |>
    dplyr::rename_with(~ c("X1", "X2")) |>
    tibble::rownames_to_column("sample")

  pcoa_eig <- round(
    pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]) * 100,
    digits = 2
  )

  plot_df <- pcoa_points |>
    dplyr::left_join(
      dplyr::select(
        group,
        sample = dplyr::all_of(sample_col),
        group = dplyr::all_of(group_col)
      ),
      by = "sample"
    ) |>
    dplyr::mutate(group = factor(group, levels = group_level))

  ## 6. 坐标轴、标题和图例
  if (is.null(xlab)) {
    xlab <- paste0("pcoa1 (", pcoa_eig[1], "%)")
  }

  if (is.null(ylab)) {
    ylab <- paste0("pcoa2 (", pcoa_eig[2], "%)")
  }

  if (is.null(legend_title)) {
    legend_title <- "Group"
  }

  if (is.null(title)) {
    title <- paste0(stringr::str_to_sentence(dist_method), "-distance pcoa")
  }

  ## 7. PERMANOVA 统计
  if (isTRUE(adonis2) && is.null(subtitle)) {
    adonis_group <- group |>
      dplyr::select(
        sample = dplyr::all_of(sample_col),
        group = dplyr::all_of(group_col)
      ) |>
      dplyr::arrange(match(sample, sample_order))

    adonis_group$group <- factor(adonis_group$group, levels = group_level)

    adonis <- vegan::adonis2(
      distance ~ group,
      data = adonis_group,
      permutations = permutations
    )

    subtitle <- substitute(
      "PERMANOVA: " * R^2 == a ~ ", " ~ italic(p) < b,
      list(
        a = round(adonis$R2[1], 4),
        b = ifelse(adonis$`Pr(>F)`[1] < 0.001, 0.001, adonis$`Pr(>F)`[1])
      )
    )
  }

  ## 8. 传入 plot_dim() 作图
  p <- plot_dim(
    data = plot_df, group_level = group_level, group_color = group_color,
    display_type = display_type, conf_type = conf_type, ellipse_level = ellipse_level,
    title = title, subtitle = subtitle,
    xlab = xlab, ylab = ylab, legend_title = legend_title,
    add_group_label = add_group_label, add_sample_label = add_sample_label,
    label_size = label_size, point_size = point_size,
    show_legend = show_legend, show_grid = show_grid,
    show_line = show_line, aspect_ratio = aspect_ratio,
    theme = theme
  )

  attr(p, "plot_df") <- plot_df
  attr(p, "distance") <- distance
  attr(p, "eigenvalue") <- pcoa_eig

  p
}

#### plot_pcoa_box ####
# 绘制 pcoa 主图，并在上方和右侧添加 PC1/PC2 分布箱线图
# profile/distance 二选一
# group: 样本分组表
# sample_col/group_col: group 中样本列和分组列

#' Plot pcoa Box utility
#'
#'
#' 同时展示 pcoa 散点图及主要坐标轴的边际箱线图。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param group A sample metadata table containing sample and group columns.
#' @param distance A precomputed distance object; when supplied, it takes precedence over `profile`.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param group_level Optional order of group levels.
#' @param group_color Optional colors aligned to `group_level`.
#' @param dist_method Distance method; available values are validated with `match.arg()`.
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
#' @param show_legend Whether to display the plot legend.
#' @param show_grid Whether to draw panel grid lines.
#' @param show_line Whether to draw horizontal and vertical reference lines at zero.
#' @param theme Plot theme preset; supported values are shown in Usage.
#' @param permutations Number of permutations used by the significance test.
#' @param ... Additional arguments passed to `plot_pcoa()`.
#' @return A ggplot-compatible plot object; computed data or fitted objects are retained as attributes when applicable.
#' @export
plot_pcoa_box <- function(
  profile = NULL, group, distance = NULL, sample_col = "sample",
  group_col = "group", group_level = NULL, group_color = NULL,
  dist_method = c("bray", "jaccard", "euclidean", "manhattan", "unifrac"),
  display_type = c("line", "point"),
  conf_type = c("ellipse", "encircle", "none"),
  ellipse_level = .75, title = NULL, subtitle = NULL, xlab = NULL,
  ylab = NULL, legend_title = NULL, add_group_label = FALSE,
  add_sample_label = FALSE, label_size = 1.5, show_legend = TRUE,
  show_grid = FALSE, show_line = TRUE, theme = c("default", "pubr"),
  permutations = 999, ...
) {
  dist_method <- .match_distance_method(dist_method)
  display_type <- match.arg(display_type)
  conf_type <- match.arg(conf_type)
  theme <- match.arg(theme)

  group <- data.frame(group, check.names = FALSE)

  if (!all(c(sample_col, group_col) %in% colnames(group))) {
    stop("group should contain columns: ", sample_col, " | ", group_col)
  }

  res <- calcu_pcoa(
    profile = profile,
    distance = distance,
    group = group,
    sample_col = sample_col,
    group_col = group_col,
    dim = 2,
    dist_method = dist_method,
    adonis2 = TRUE,
    permutations = permutations,
    ...
  )

  meta <- group |>
    dplyr::select(
      sample = dplyr::all_of(sample_col),
      group = dplyr::all_of(group_col)
    )

  plot_df <- res$point
  colnames(plot_df)[2:3] <- c("X1", "X2")

  plot_df <- plot_df |>
    dplyr::left_join(meta, by = "sample") |>
    dplyr::filter(!is.na(group))

  if (is.null(group_level)) {
    group_level <- unique(plot_df$group)
  }

  group_level <- group_level[group_level %in% unique(plot_df$group)]

  if (length(group_level) == 0) {
    stop("No valid group remained for plotting.")
  }

  group_color <- .resolve_group_colors(group_level, group_color)

  plot_df <- plot_df |>
    dplyr::mutate(group = factor(group, levels = group_level))

  if (is.null(xlab)) {
    xlab <- paste0("pcoa1 (", res$eigenvalue[1], "%)")
  }

  if (is.null(ylab)) {
    ylab <- paste0("pcoa2 (", res$eigenvalue[2], "%)")
  }

  if (is.null(legend_title)) {
    legend_title <- "Group"
  }

  if (is.null(title)) {
    title <- paste0(stringr::str_to_sentence(dist_method), "-distance pcoa")
  }

  if (is.null(subtitle)) {
    subtitle <- substitute(
      "PERMANOVA: " * R^2 == a ~ ", " ~ italic(p) < b,
      list(
        a = round(res$adonis2_r2, 4),
        b = ifelse(res$adonis2_pval < 0.001, 0.001, res$adonis2_pval)
      )
    )
  }

  ## 1. pcoa 主图
  p_main <- plot_dim(
    data = plot_df,
    group_level = group_level,
    group_color = group_color,
    display_type = display_type,
    conf_type = conf_type,
    ellipse_level = ellipse_level,
    xlab = xlab,
    ylab = ylab,
    legend_title = legend_title,
    add_group_label = add_group_label,
    add_sample_label = add_sample_label,
    label_size = label_size,
    show_legend = FALSE,
    show_grid = show_grid,
    show_line = show_line,
    theme = theme
  ) +
    ggplot2::scale_x_continuous(expand = c(.01, .01)) +
    ggplot2::scale_y_continuous(expand = c(.01, .01))

  x_limits <- ggplot2::ggplot_build(p_main)$layout$panel_params[[1]]$x.range
  y_limits <- ggplot2::ggplot_build(p_main)$layout$panel_params[[1]]$y.range

  ## 2. PC1 顶部箱线图
  top_df <- plot_df |>
    dplyr::select(group, X1)

  top_diff <- top_df |>
    rstatix::pairwise_wilcox_test(X1 ~ group) |>
    dplyr::mutate(comparison = paste0(group1, "-", group2)) |>
    dplyr::pull(p, name = comparison) |>
    multcompView::multcompLetters()

  top_label <- data.frame(label = top_diff$Letters) |>
    tibble::rownames_to_column("group") |>
    dplyr::left_join(
      top_df |>
        dplyr::group_by(group) |>
        dplyr::slice_min(order_by = X1, n = 1),
      by = "group"
    )

  p_top <- ggplot2::ggplot(top_df, ggplot2::aes(X1, group, color = group)) +
    ggplot2::geom_boxplot(fill = NA, outlier.shape = NA, show.legend = FALSE, width = .61) +
    ggplot2::geom_jitter(size = 1.2, height = .2, show.legend = FALSE) +
    ggplot2::geom_text(
      data = top_label,
      ggplot2::aes(x = X1, y = group, label = label),
      inherit.aes = FALSE,
      vjust = .5,
      size = 4
    ) +
    ggplot2::scale_color_manual(values = group_color, drop = FALSE) +
    ggplot2::scale_x_continuous(expand = c(.01, .01), limits = x_limits) +
    ggplot2::labs(title = title, subtitle = subtitle, y = "", x = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.ticks = ggplot2::element_line(linewidth = .5, color = "black"),
      axis.ticks.length = grid::unit(2, "mm"),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 12, color = "black"),
      axis.text.x = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 12, color = "black", face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 12, color = "black"),
      plot.margin = grid::unit(c(0, 0, 0, 0), "mm"),
      panel.border = ggplot2::element_rect(linewidth = .5, color = "black", fill = NA),
      panel.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "none"
    )

  ## 3. PC2 右侧箱线图
  right_df <- plot_df |>
    dplyr::select(group, X2)

  right_diff <- right_df |>
    rstatix::pairwise_wilcox_test(X2 ~ group) |>
    dplyr::mutate(comparison = paste0(group1, "-", group2)) |>
    dplyr::pull(p, name = comparison) |>
    multcompView::multcompLetters()

  right_label <- data.frame(label = right_diff$Letters) |>
    tibble::rownames_to_column("group") |>
    dplyr::left_join(
      right_df |>
        dplyr::group_by(group) |>
        dplyr::slice_max(order_by = X2, n = 1),
      by = "group"
    )

  p_right <- ggplot2::ggplot(right_df, ggplot2::aes(group, X2, color = group)) +
    ggplot2::geom_boxplot(fill = NA, outlier.shape = NA, show.legend = FALSE, width = .61) +
    ggplot2::geom_jitter(size = 1.2, width = .2, show.legend = show_legend) +
    ggplot2::geom_text(
      data = right_label,
      ggplot2::aes(x = group, y = X2, label = label),
      inherit.aes = FALSE,
      vjust = .5,
      size = 4
    ) +
    ggplot2::scale_color_manual(values = group_color, drop = FALSE) +
    ggplot2::scale_y_continuous(expand = c(.01, .01), limits = y_limits) +
    ggplot2::labs(y = "", x = "", color = legend_title) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.ticks = ggplot2::element_line(linewidth = .5, color = "black"),
      axis.ticks.length = grid::unit(2, "mm"),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 12, color = "black"),
      axis.text.y = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      plot.margin = grid::unit(c(0, 0, 0, 0), "mm"),
      panel.border = ggplot2::element_rect(linewidth = .5, color = "black", fill = NA),
      panel.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 10, color = "black"),
      legend.title = ggplot2::element_text(size = 10, color = "black")
    )

  p <- ggpubr::ggarrange(
    p_top, NULL,
    p_main, p_right,
    ncol = 2,
    nrow = 2,
    heights = c(1.7, 3),
    widths = c(3, 1.5),
    align = "hv"
  )

  attr(p, "pcoa") <- res

  p
}
