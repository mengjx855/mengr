#### Jin-Xin Meng, 20250307, 20260827, v0.2.1 ####

# 20250417: plot function pass to plot_dim()
# 20251129: 未来考虑 NMDS 中加入 vegan::anosim 分析
# 20260419: 修改anosim输入矩阵为距离矩阵
# 20260827: update function.


#### plot_NMDS ####
# 绘制 NMDS 图
# profile: 行为 feature，列为 sample 的丰度表
# distance: 已经计算好的距离矩阵或 dist 对象；如果提供 distance，则优先使用 distance
# group: 样本分组信息表
# sample_col/group_col: group 中样本列和分组列
# sub_sample: 指定样本子集
# sub_group: 指定分组子集
# dist_method: 距离方法，默认 bray
# transform: 计算距离前是否转换，默认 total；NULL 表示不转换
# anosim: 是否添加 ANOSIM 结果作为 subtitle

#' Plot NMDS utility
#'
#' `plot_NMDS()` provides a reusable mengR workflow with input validation and standardized
#'   output.
#'
#' Chinese summary: 从 profile 或距离对象执行 NMDS，可附加 ANOSIM 结果。
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
#' @param add_group_label Logical control for `add_group_label`.
#' @param add_sample_label Logical control for `add_sample_label`.
#' @param label_size Numeric setting for `label_size`.
#' @param point_size Numeric setting for `point_size`.
#' @param show_legend Logical control for `show_legend`.
#' @param show_grid Logical control for `show_grid`.
#' @param show_line Logical control for `show_line`.
#' @param aspect_ratio Panel aspect ratio passed to `ggplot2::theme()`.
#' @param theme Plot theme preset; supported values are shown in Usage.
#' @param anosim Whether to run ANOSIM and annotate its statistic and P value.
#' @param permutations Number of permutations used by the significance test.
#' @param trymax Maximum number of random starts attempted by NMDS.
#' @param ... Additional arguments passed to the underlying function.
#' @return A plot object; analysis data or models may also be stored as attributes.
#' @export
plot_NMDS <- function(
    profile = NULL, group, distance = NULL,
    sample_col = 'sample', group_col = 'group',
    group_level = NULL, group_color = NULL,
    sub_sample = NULL, sub_group = NULL,
    dist_method = c('bray', 'jaccard', 'euclidean', 'manhattan'),
    transform = c('hellinger', 'total', 'pa', 'clr'),
    display_type = c('line', 'point'),
    conf_type = c('ellipse', 'encircle', 'none'),
    ellipse_level = .75, title = NULL, subtitle = NULL,
    xlab = NULL, ylab = NULL, legend_title = NULL,
    add_group_label = FALSE, add_sample_label = FALSE,
    label_size = 2, point_size = 1.5,
    show_legend = TRUE, show_grid = FALSE, show_line = TRUE,
    aspect_ratio = 3/4, theme = c('default', 'pubr'),
    anosim = TRUE, permutations = 999, trymax = 100, ...
) {
  dist_method <- .match_distance_method(dist_method)
  transform <- .match_transform_method(transform)
  display_type <- match.arg(display_type)
  conf_type <- match.arg(conf_type)
  theme <- match.arg(theme)

  group <- data.frame(group, check.names = FALSE)
  
  if (!all(c(sample_col, group_col) %in% colnames(group))) {
    stop('group should contain columns: ', sample_col, ' | ', group_col)
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
      stop('Need either profile or distance.')
    }
    
    profile <- data.frame(profile, check.names = FALSE)
    profile <- as.matrix(profile)
    suppressWarnings(storage.mode(profile) <- 'numeric')
    profile[!is.finite(profile)] <- 0
    
    sample_use <- intersect(group[[sample_col]], colnames(profile))
    if (length(sample_use) < 2) {
      stop('NMDS requires at least two matched samples.')
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
      stop('No valid features remained for NMDS.')
    }
    
    if (ncol(profile) < 2) {
      stop('NMDS requires at least two valid samples.')
    }
    
    profile_t <- t(profile)
    
    if (!is.null(transform)) {
      profile_t <- vegan::decostand(profile_t, method = transform)
    }
    
    profile_t[!is.finite(profile_t)] <- 0
    
    ## 转换后再次删除全 0 sample，避免 jaccard / bray 出现 empty rows
    profile_t <- profile_t[rowSums(profile_t, na.rm = TRUE) != 0, , drop = FALSE]
    
    if (nrow(profile_t) < 2) {
      stop('NMDS requires at least two valid samples after transformation.')
    }
    
    distance <- vegan::vegdist(profile_t, method = dist_method, na.rm = TRUE, ...)
    
  } else {
    
    distance <- stats::as.dist(distance)
    
    distance_sample <- labels(distance)
    sample_use <- distance_sample[distance_sample %in% group[[sample_col]]]
    if (length(sample_use) < 2) {
      stop('NMDS requires at least two matched samples.')
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
    stop('Some samples in distance object are missing in group.')
  }
  
  ## 4. 设置分组顺序和颜色
  if (is.null(group_level)) {
    group_level <- as.character(unique(group[[group_col]]))
  }
  
  group_level <- group_level[group_level %in% unique(group[[group_col]])]
  
  if (length(group_level) == 0) {
    stop('No valid group remained for plotting.')
  }
  
  group_color <- .resolve_group_colors(group_level, group_color)
  
  ## 5. NMDS 降维
  NMDS <- vegan::metaMDS(
    distance, k = 2, trace = 0, trymax = trymax
  )
  
  NMDS_points <- data.frame(
    NMDS$points[, 1:2, drop = FALSE],
    check.names = FALSE
  ) |>
    dplyr::rename_with(~ c('X1', 'X2')) |>
    tibble::rownames_to_column('sample')
  
  NMDS_stress <- round(NMDS$stress, digits = 4)
  
  plot_df <- NMDS_points |>
    dplyr::left_join(
      dplyr::select(
        group,
        sample = dplyr::all_of(sample_col),
        group = dplyr::all_of(group_col)
      ),
      by = 'sample'
    ) |>
    dplyr::mutate(group = factor(group, levels = group_level))
  
  ## 6. 坐标轴、标题和图例
  if (is.null(xlab)) {
    xlab <- 'MDS1'
  }
  
  if (is.null(ylab)) {
    ylab <- 'MDS2'
  }
  
  if (is.null(legend_title)) {
    legend_title <- 'Group'
  }
  
  if (is.null(title)) {
    title <- paste0(stringr::str_to_sentence(dist_method), '-distance NMDS')
  }
  
  ## 7. ANOSIM 统计
  anosim_result <- NULL
  if (isTRUE(anosim) && is.null(subtitle)) {
    
    anosim_group <- group |>
      dplyr::select(
        sample = dplyr::all_of(sample_col),
        group = dplyr::all_of(group_col)
      ) |>
      dplyr::arrange(match(sample, sample_order))
    
    anosim_group$group <- factor(anosim_group$group, levels = group_level)
    
    anosim_result <- vegan::anosim(
      distance,
      as.vector(anosim_group$group),
      permutations = permutations
    )
    
    subtitle <- substitute(
      'Stress =' ~ a ~ ', ANOSIM ' * R == b ~ ', ' ~ italic(p) < c,
      list(
        a = NMDS_stress,
        b = round(anosim_result$statistic, 4),
        c = ifelse(
          anosim_result$signif < 0.001, 0.001, anosim_result$signif
        )
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
    theme = theme, ...
  )
  
  attr(p, 'plot_df') <- plot_df
  attr(p, 'distance') <- distance
  attr(p, 'NMDS') <- NMDS
  
  if (!is.null(anosim_result)) attr(p, 'anosim') <- anosim_result
  
  return(p)
}

#### calcu_pairwise_anosim ####
# 两两 ANOSIM
# profile: 行为 feature，列为 sample 的丰度表
# group: 样本分组表
# sample_col/group_col: group 中样本列和分组列
# group_level: 指定分组顺序
# transform: 计算距离前是否转换，默认 total；NULL 表示不转换

#' Calcu Pairwise Anosim utility
#'
#' `calcu_pairwise_anosim()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 对各分组组合执行 pairwise ANOSIM。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param group A sample metadata table containing sample and group columns.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param group_level Optional order of group levels.
#' @param dist_method Distance method; available values are validated with `match.arg()`.
#' @param transform Optional transformation applied before analysis.
#' @param permutations Number of permutations used by the significance test.
#' @param add_plab Logical control for `add_plab`.
#' @param ... Additional arguments passed to the underlying function.
#' @return A result object described in the Details section.
#' @export
calcu_pairwise_anosim <- function(
    profile, group, sample_col = 'sample', group_col = 'group',
    group_level = NULL,
    dist_method = c('bray', 'jaccard', 'euclidean', 'manhattan'),
    transform = c('total', 'hellinger', 'pa', 'clr'),
    permutations = 999, add_plab = TRUE, ...
) {
  dist_method <- .match_distance_method(dist_method)
  transform <- .match_transform_method(transform)
  
  profile <- data.frame(profile, check.names = FALSE)
  group <- data.frame(group, check.names = FALSE)
  
  if (!all(c(sample_col, group_col) %in% colnames(group))) {
    stop('group should contain columns: ', sample_col, ' | ', group_col)
  }
  
  if (is.null(group_level)) {
    group_level <- as.character(unique(group[[group_col]]))
  }
  
  group_level <- group_level[group_level %in% unique(group[[group_col]])]
  
  if (length(group_level) < 2) {
    stop('At least two groups are required.')
  }
  
  sample_use <- intersect(group[[sample_col]], colnames(profile))
  
  if (length(sample_use) < 2) {
    stop('At least two matched samples are required.')
  }
  
  group <- group |>
    dplyr::filter(.data[[sample_col]] %in% sample_use)
  
  profile <- profile[, group[[sample_col]], drop = FALSE]
  
  if (!is.null(transform)) {
    profile <- vegan::decostand(
      profile, method = transform, MARGIN = 2
    )
  }
  
  data <- purrr::map_dfr(
    utils::combn(group_level, m = 2, simplify = FALSE), \(x) {
      
      meta <- group |>
        dplyr::filter(.data[[group_col]] %in% x) |>
        dplyr::filter(.data[[sample_col]] %in% colnames(profile)) |>
        dplyr::select(
          sample = dplyr::all_of(sample_col),
          group = dplyr::all_of(group_col)
        )
      
      if (nrow(meta) < 3 || length(unique(meta$group)) < 2) {
        return(
          data.frame(
            comparison = paste0(x, collapse = '_vs_'),
            r = NA_real_,
            pval = NA_real_,
            check.names = FALSE
          )
        )
      }
      
      profile_x <- profile[, meta$sample, drop = FALSE]
      
      meta <- meta |>
        dplyr::arrange(match(sample, colnames(profile_x)))
      
      meta$group <- factor(meta$group, levels = x)
      
      anosim_res <- vegan::anosim(
        t(profile_x),
        as.vector(meta$group),
        permutations = permutations,
        distance = dist_method,
        ...
      )
      
      data.frame(
        comparison = paste0(x, collapse = '_vs_'),
        r = round(anosim_res$statistic, 4),
        pval = anosim_res$signif,
        check.names = FALSE
      )
    }
  )
  
  if (isTRUE(add_plab)) {
    data <- data |>
      dplyr::mutate(
        plab = cut(
          pval,
          breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
          labels = c('***', '**', '*', 'ns')
        )
      )
  }
  
  return(data)
}

#### plot_pairwise_anosim ####
# 绘制两两 ANOSIM 结果
# data: calcu_pairwise_anosim() 输出结果
# group_level: 指定分组顺序

#' Plot Pairwise Anosim utility
#'
#' `plot_pairwise_anosim()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 将 pairwise ANOSIM 结果绘制为矩阵式结果图。
#'
#' @param data An input data frame or compatible object.
#' @param group_level Optional order of group levels.
#' @return A plot object; analysis data or models may also be stored as attributes.
#' @export
plot_pairwise_anosim <- function(data, group_level = NULL) {
  
  data <- data.frame(data, check.names = FALSE)
  
  if (!all(c('comparison', 'r', 'pval') %in% colnames(data))) {
    stop('data should contain columns: comparison | r | pval')
  }
  
  if (is.null(group_level)) {
    group_level <- unique(c(
      stringr::str_split_i(data$comparison, '_vs_', 1),
      stringr::str_split_i(data$comparison, '_vs_', 2)
    ))
  }
  
  plot_df <- data.frame(
    x = stringr::str_split_i(data$comparison, '_vs_', 1),
    y = stringr::str_split_i(data$comparison, '_vs_', 2),
    r = data$r,
    pval = data$pval,
    check.names = FALSE
  ) |>
    dplyr::mutate(
      x = factor(x, levels = group_level),
      y = factor(y, levels = rev(group_level)),
      r_size = abs(r),
      plab = dplyr::case_when(
        pval <= 0.001 ~ 'p≤0.001',
        pval < 0.01 ~ 'p<0.01',
        pval < 0.05 ~ 'p<0.05',
        TRUE ~ 'p≥0.05'
      ),
      plab = factor(
        plab,
        levels = c('p≤0.001', 'p<0.01', 'p<0.05', 'p≥0.05')
      )
    )
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x, y)) +
    ggplot2::geom_tile(
      fill = 'transparent',
      color = 'black',
      width = 1,
      height = 1,
      linewidth = .4
    ) +
    ggplot2::geom_point(
      ggplot2::aes(size = r_size, fill = plab),
      shape = 21,
      color = 'black',
      stroke = .4
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        'p≤0.001' = '#f46d43',
        'p<0.01' = '#fee08b',
        'p<0.05' = '#abdda4',
        'p≥0.05' = '#3288bd'
      ),
      breaks = c('p≤0.001', 'p<0.01', 'p<0.05', 'p≥0.05')
    ) +
    ggplot2::scale_size_continuous(range = c(6, 12)) +
    ggplot2::labs(x = '', y = '') +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 10, color = 'black'),
      axis.title = ggplot2::element_text(size = 10, color = 'black'),
      plot.title = ggplot2::element_text(size = 10, color = 'black'),
      panel.border = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      aspect.ratio = 1
    ) +
    ggplot2::guides(
      size = ggplot2::guide_legend(title = 'ANOSIM |R|', order = 1),
      fill = ggplot2::guide_legend(
        title = 'Significance',
        order = 2,
        override.aes = list(size = 4)
      )
    )
  
  return(p)
}
