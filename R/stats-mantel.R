#### Jin-Xin Meng, 20260606, 20260607, 0.1.0 ####

# 20260606: add calcu_mantel()
# 20260607: add make_curve_path(), plot_mantel_lower(), plot_mantel_upper(),

#### calcu_mantel ####
# 根据 feature 分组计算 Mantel test
# profile: 行为 feature，列为 sample 的丰度表
# metadata: feature 与 category/taxa/pathway class 的对应关系
# envs: 环境因子表，包含 sample_col 和环境变量列，
#   多个环境变量合成一个环境距离矩阵：需要 scale
#   每个环境变量单独做 Mantel：不需要 scale
# feature_col/category_col: metadata 中 feature 和分组列
# sample_col: envs 中样本列
# category_keep: 指定保留哪些 category；NULL 表示全部保留
# trans_ra: 是否将 profile 转换为百分比相对丰度
# remove_empty: 是否删除全 0 feature/sample/env variable

#' Calcu Mantel utility
#'
#' `calcu_mantel()` provides a reusable mengR workflow with input validation and standardized
#'   output.
#'
#' Chinese summary: 在多个 profile 或环境变量之间批量执行 Mantel test。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param metadata A metadata or annotation data frame.
#' @param envs Environmental or host-variable table used in Mantel analyses.
#' @param feature_col Name of the feature-identifier column.
#' @param category_col Name of the `category_col` input column.
#' @param sample_col Name of the sample-identifier column.
#' @param category_keep Optional environmental categories retained for Mantel tests.
#' @param trans_ra Whether abundances are converted to relative abundance before analysis.
#' @param base Numeric setting for `base`.
#' @param remove_empty Logical control for `remove_empty`.
#' @param na_fill Value used to replace missing observations before analysis.
#' @param parallel Whether to enable the parallel behavior.
#' @param ... Additional arguments passed to the underlying function.
#' @return A result object described in the Details section.
#' @export
calcu_mantel <- function(
    profile, metadata, envs,
    feature_col = 'name', category_col = 'category',
    sample_col = 'sample', category_keep = NULL,
    trans_ra = FALSE, base = 100, remove_empty = TRUE,
    na_fill = NULL, parallel = 1, ...
) {
  
  profile <- data.frame(profile, check.names = FALSE)
  metadata <- data.frame(metadata, check.names = FALSE)
  envs <- data.frame(envs, check.names = FALSE)
  
  if (!all(c(feature_col, category_col) %in% colnames(metadata))) {
    stop('metadata should contain columns: ', feature_col, ' | ', category_col)
  }
  
  if (!sample_col %in% colnames(envs)) {
    stop('envs should contain sample_col: ', sample_col)
  }
  
  ## 1. 对齐样本
  sample_use <- intersect(envs[[sample_col]], colnames(profile))
  
  if (length(sample_use) < 3) {
    stop('Mantel test requires at least three matched samples.')
  }
  
  profile <- profile[, sample_use, drop = FALSE]
  
  envs <- envs |>
    dplyr::filter(.data[[sample_col]] %in% sample_use) |>
    dplyr::arrange(match(.data[[sample_col]], sample_use))
  
  if (!identical(envs[[sample_col]], sample_use)) {
    stop('Sample order between profile and envs is not matched.')
  }
  
  envs <- envs |>
    dplyr::select(-dplyr::all_of(sample_col))
  
  ## 2. 转换为数值矩阵
  profile <- as.matrix(profile)
  suppressWarnings(storage.mode(profile) <- 'numeric')
  profile[!is.finite(profile)] <- 0
  
  envs <- as.matrix(envs)
  suppressWarnings(storage.mode(envs) <- 'numeric')
  
  if (any(!is.finite(envs))) {
    stop('envs contains NA/NaN/Inf values. Please filter or impute them before Mantel test.')
  }
  
  ## 3. 删除空 feature / sample / env variable
  if (isTRUE(remove_empty)) {
    
    profile <- profile[
      rowSums(profile, na.rm = TRUE) != 0,
      colSums(profile, na.rm = TRUE) != 0,
      drop = FALSE
    ]
    
    envs <- envs[
      ,
      colSums(abs(envs), na.rm = TRUE) != 0,
      drop = FALSE
    ]
  }
  
  if (nrow(profile) == 0) {
    stop('No valid features remained in profile.')
  }
  
  if (ncol(profile) < 3) {
    stop('Mantel test requires at least three valid samples.')
  }
  
  if (ncol(envs) == 0) {
    stop('No valid environmental variables remained in envs.')
  }
  
  ## 4. 是否转换为相对丰度
  if (isTRUE(trans_ra)) {
    
    col_sum <- colSums(profile, na.rm = TRUE)
    col_sum[col_sum == 0] <- NA_real_
    
    profile <- sweep(profile, 2, col_sum, '/') * base
    profile[!is.finite(profile)] <- 0
  }
  
  ## 5. 匹配 feature 与 category
  meta <- data.frame(
    feature = rownames(profile),
    check.names = FALSE
  ) |>
    dplyr::left_join(
      metadata |>
        dplyr::select(
          feature = dplyr::all_of(feature_col),
          category = dplyr::all_of(category_col)
        ),
      by = 'feature'
    )
  
  if (!is.null(na_fill)) {
    meta$category[is.na(meta$category) | meta$category == ''] <- na_fill
  }
  
  meta <- meta |>
    dplyr::filter(!is.na(category), category != '')
  
  if (!is.null(category_keep)) {
    meta <- meta |>
      dplyr::filter(category %in% category_keep)
  }
  
  if (nrow(meta) == 0) {
    stop('No features remained after matching metadata.')
  }
  
  profile <- profile[meta$feature, , drop = FALSE]
  
  ## 6. 转为 sample × feature
  spec <- t(profile) |>
    data.frame(check.names = FALSE)
  
  envs <- data.frame(envs, check.names = FALSE)
  
  ## 7. 生成 spec_select
  spec_info <- data.frame(
    feature = colnames(spec),
    category = meta$category,
    index = seq_len(ncol(spec)),
    check.names = FALSE
  )
  
  spec_select <- split(spec_info$index, spec_info$category)
  
  ## 8. Mantel test
  test <- linkET::mantel_test(
    spec, envs, spec_select = spec_select,
    parallel = parallel, ...
  ) |>
    data.frame(check.names = FALSE)
  
  attr(test, 'spec') <- spec
  attr(test, 'envs') <- envs
  attr(test, 'metadata') <- meta
  attr(test, 'spec_select') <- spec_select
  
  return(test)
}

#### make_curve_path ####
# 生成每条 edge 的二次 Bezier 曲线坐标
# 不依赖 linkET
# data:
#   已经包含起点和终点坐标的数据框。
# x, y:
#   曲线起点坐标列名。
#   在当前 Mantel lower 图中，一般是 spec_x / spec_y。
# xend, yend:
#   曲线终点坐标列名。
#   在当前 Mantel lower 图中，一般是 anchor_x / anchor_y。
# to:
#   用于决定曲率方向的变量，一般是 env。
#   这里不直接根据 x/y 坐标大小判断曲率方向，而是根据 env 在矩阵中的顺序决定。
# to_levels:
#   env 的完整顺序。建议传入 env_names。
# bend:
#   曲线弯曲强度。值越大，曲线越弯。
# n:
#   每条曲线插值点数。越大越平滑。
# reverse:
#   是否反转曲率方向。如果曲线方向和预期相反，改成 TRUE。

#' Make Curve Path utility
#'
#' `make_curve_path()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 为 Mantel 图中的节点连线生成平滑曲线路径坐标。
#'
#' @param data An input data frame or compatible object.
#' @param x Primary vector or object supplied to the utility.
#' @param y Secondary vector or object supplied to the utility.
#' @param xend Ending x coordinate of the generated connection path.
#' @param yend Ending y coordinate of the generated connection path.
#' @param to Target taxonomy rank or identifier type produced by the conversion.
#' @param to_levels Numeric setting for `to_levels`.
#' @param bend Signed curvature of the generated connection path.
#' @param n Requested number of values, features, or results.
#' @param reverse Logical control for `reverse`.
#' @return A result object described in the Details section.
#' @export
make_curve_path <- function(
    data, x = 'spec_x', y = 'spec_y', xend = 'anchor_x', yend = 'anchor_y',
    to = 'env', to_levels = NULL, bend = 0.18, n = 80, reverse = FALSE
) {
  
  data <- data.frame(data, check.names = FALSE)
  
  ## 避免重复生成 edge_id 后出现 edge_id...1 / edge_id...2
  data <- data[, !grepl('^edge_id', colnames(data)), drop = FALSE]
  data$edge_id <- seq_len(nrow(data))
  
  if (is.null(to_levels)) {
    to_levels <- unique(as.character(data[[to]]))
  }
  
  ## 根据 env 在矩阵顺序中的位置决定曲率正负
  ## 这个逻辑类似 linkET::geom_couple() 中 nice_curvature(by = 'to') 的思路
  id <- stats::setNames(seq_along(to_levels), rev(to_levels))
  half <- length(to_levels) / 2
  
  data$.curve_sign <- ifelse(
    id[as.character(data[[to]])] >= half,
    1,
    -1
  )
  
  data$.curve_sign[is.na(data$.curve_sign)] <- 1
  
  if (isTRUE(reverse)) {
    data$.curve_sign <- -data$.curve_sign
  }
  
  res <- lapply(seq_len(nrow(data)), function(i) {
    
    x0 <- data[[x]][i]
    y0 <- data[[y]][i]
    x1 <- data[[xend]][i]
    y1 <- data[[yend]][i]
    
    dx <- x1 - x0
    dy <- y1 - y0
    len <- sqrt(dx^2 + dy^2)
    
    if (!is.finite(len) || len == 0) {
      len <- 1
    }
    
    ## 起点和终点中点
    mx <- (x0 + x1) / 2
    my <- (y0 + y1) / 2
    
    ## 法向量，用于确定 Bezier 控制点偏移方向
    nx <- -dy / len
    ny <-  dx / len
    
    ## 每条边自己的控制点
    cx <- mx + data$.curve_sign[i] * nx * len * abs(bend)
    cy <- my + data$.curve_sign[i] * ny * len * abs(bend)
    
    tt <- seq(0, 1, length.out = n)
    
    curve <- data.frame(
      edge_id = data$edge_id[i],
      t = tt,
      x = (1 - tt)^2 * x0 + 2 * (1 - tt) * tt * cx + tt^2 * x1,
      y = (1 - tt)^2 * y0 + 2 * (1 - tt) * tt * cy + tt^2 * y1,
      check.names = FALSE
    )
    
    meta <- data[rep(i, n), , drop = FALSE]
    meta <- meta[
      ,
      setdiff(colnames(meta), c(x, y, xend, yend, 'edge_id')),
      drop = FALSE
    ]
    
    cbind(curve, meta)
  })
  
  out <- do.call(rbind, res)
  rownames(out) <- NULL
  
  return(out)
}


#### plot_mantel_lower ####
# 绘制 lower 版本 Mantel test + 环境因子相关矩阵图
# test: Mantel test 结果表，至少包含 spec、env、r、p 四列。
# envs: 环境因子表。需要包含 sample_col 指定的样本列。
# spec_range:
#   控制 spec 节点在平行线上的分布范围。
#   例如 c(0.45, 0.95) 表示 y 轴 45% ~ 95%。
# spec_offset: 控制 spec 节点所在平行线相对矩阵对角线向外偏移多少。
# curve_bend: 控制 Mantel 连线弯曲程度。
# curve_reverse: 是否反转曲线方向。

#' Plot Mantel Lower utility
#'
#' `plot_mantel_lower()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 在环境相关矩阵下三角区域叠加 Mantel 关联曲线。
#'
#' @param test Statistical-test result table used for plotting.
#' @param envs Environmental or host-variable table used in Mantel analyses.
#' @param sample_col Name of the sample-identifier column.
#' @param env_cols Environmental-variable columns selected from `envs`.
#' @param spec_col Name of the `spec_col` input column.
#' @param env_col Name of the `env_col` input column.
#' @param r_col Name of the `r_col` input column.
#' @param p_col Name of the `p_col` input column.
#' @param cor_method Correlation method: Pearson, Spearman, or Kendall.
#' @param scale_env Whether environmental variables are standardized before correlation analysis.
#' @param p_filter Maximum P value retained for plotting or downstream analysis.
#' @param r_breaks Numeric breakpoints used to categorize correlation or Mantel r.
#' @param r_labels Labels corresponding to intervals defined by `r_breaks`.
#' @param r_size Numeric setting for `r_size`.
#' @param r_cut_abs Whether absolute r values are used when assigning line-width categories.
#' @param p_breaks Numeric breakpoints used to categorize P values.
#' @param p_labels Labels corresponding to intervals defined by `p_breaks`.
#' @param p_color Color specification for `p_color`.
#' @param fill_color Color specification for `fill_color`.
#' @param grid_col Name of the `grid_col` input column.
#' @param grid_linewidth Numeric setting for `grid_linewidth`.
#' @param square_size_range Numeric setting for `square_size_range`.
#' @param spec_range Relative range reserved for the species-label region.
#' @param spec_offset Offset between species labels and the heatmap boundary.
#' @param spec_hjust Horizontal justification applied to species-side labels.
#' @param curve_bend Signed curvature used for Mantel connection paths.
#' @param curve_n Number of interpolation points used for each connection path.
#' @param curve_reverse Logical control for `curve_reverse`.
#' @param node_fill Color specification for node fill.
#' @param node_size Numeric setting for `node_size`.
#' @param title Optional plot or result title.
#' @param env_label_size Numeric setting for `env_label_size`.
#' @param spec_label_size Numeric setting for `spec_label_size`.
#' @return A plot object; analysis data or models may also be stored as attributes.
#' @export
plot_mantel_lower <- function(
    test, envs, sample_col = 'sample', env_cols = NULL,
    spec_col = 'spec', env_col = 'env', r_col = 'r', p_col = 'p',
    cor_method = c('pearson', 'spearman', 'kendall'),
    scale_env = TRUE, p_filter = 0.05,
    ## Mantel r 分段
    r_breaks = c(0.1, 0.2, 0.3, Inf),
    r_labels = c('0.1 ~ 0.2', '0.2 ~ 0.3', '> 0.3'),
    r_size = c(0.7, 1.4, 2.1),
    r_cut_abs = FALSE, 
    ## Mantel p 分段
    p_breaks = c(-Inf, 0.001, 0.01, 0.05),
    p_labels = c('<= 0.001', '0.001 ~ 0.01', '0.01 ~ 0.05'),
    p_color = c(
      '<= 0.001' = '#d95f02',
      '0.001 ~ 0.01' = '#1f9f79',
      '0.01 ~ 0.05' = 'grey77'
    ),
    ## Pearson 相关矩阵
    fill_color = RColorBrewer::brewer.pal(7, 'PRGn'),
    grid_col = 'black', grid_linewidth = 0.35, square_size_range = c(2, 6),
    ## spec 节点位置
    spec_range = c(0.45, 0.95), spec_offset = 0.5, spec_hjust = 0,
    ## 曲线
    curve_bend = 0.12, curve_n = 80, curve_reverse = TRUE,
    ## node 样式
    node_fill = 'blue', node_size = 2, 
    ## text
    title = 'Mantel test', env_label_size = 3, spec_label_size = 3
) {
  cor_method <- match.arg(cor_method)
  
  test <- data.frame(test, check.names = FALSE)
  envs <- data.frame(envs, check.names = FALSE)
  
  if (!all(c(spec_col, env_col, r_col, p_col) %in% colnames(test))) {
    stop(
      'test should contain columns: ',
      paste(c(spec_col, env_col, r_col, p_col), collapse = ' | ')
    )
  }
  
  if (!sample_col %in% colnames(envs)) {
    stop('envs should contain sample_col: ', sample_col)
  }
  
  if (length(r_labels) != length(r_breaks) - 1) {
    stop('length(r_labels) should be length(r_breaks) - 1.')
  }
  
  if (length(p_labels) != length(p_breaks) - 1) {
    stop('length(p_labels) should be length(p_breaks) - 1.')
  }
  
  if (is.null(names(r_size))) {
    r_size <- stats::setNames(
      rep(r_size, length.out = length(r_labels)),
      r_labels
    )
  }
  
  if (is.null(names(p_color))) {
    p_color <- stats::setNames(
      rep(p_color, length.out = length(p_labels)),
      p_labels
    )
  }
  
  ## 1. 整理环境变量矩阵
  if (is.null(env_cols)) {
    env_cols <- setdiff(colnames(envs), sample_col)
  }
  
  if (!all(env_cols %in% colnames(envs))) {
    stop('Some env_cols are not found in envs.')
  }
  
  env_mat <- envs[, env_cols, drop = FALSE]
  env_mat <- as.matrix(env_mat)
  suppressWarnings(storage.mode(env_mat) <- 'numeric')
  
  keep_env <- apply(env_mat, 2, function(x) {
    x <- x[is.finite(x)]
    length(x) > 2 && stats::sd(x, na.rm = TRUE) > 0
  })
  
  env_mat <- env_mat[, keep_env, drop = FALSE]
  
  if (ncol(env_mat) < 2) {
    stop('Need at least two valid environmental variables.')
  }
  
  if (isTRUE(scale_env)) {
    env_mat <- scale(env_mat)
    env_mat[!is.finite(env_mat)] <- 0
  }
  
  env_names <- colnames(env_mat)
  n_env <- length(env_names)
  
  ## 2. 计算环境变量 Pearson 相关矩阵
  cor_mat <- stats::cor(
    env_mat,
    method = cor_method,
    use = 'pairwise.complete.obs'
  )
  
  cor_long_df <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE)
  
  cor_data <- data.frame(
    env_y = cor_long_df$Var1,
    env_x = cor_long_df$Var2,
    cor = cor_long_df$Freq,
    check.names = FALSE
  )
  
  cor_data$row_id <- match(cor_data$env_y, env_names)
  cor_data$col_id <- match(cor_data$env_x, env_names)
  
  ## lower triangle，不画对角线
  cor_data <- cor_data[
    cor_data$col_id < cor_data$row_id,
    ,
    drop = FALSE
  ]
  
  cor_data$x <- n_env - cor_data$row_id + 1
  cor_data$y <- cor_data$col_id
  
  ## 3. 整理 Mantel 结果
  plot_df <- data.frame(
    spec = test[[spec_col]],
    env = test[[env_col]],
    mantel_r = test[[r_col]],
    mantel_p = test[[p_col]],
    check.names = FALSE
  )
  
  plot_df <- plot_df[
    plot_df$env %in% env_names &
      is.finite(plot_df$mantel_r) &
      is.finite(plot_df$mantel_p),
    ,
    drop = FALSE
  ]
  
  if (!is.null(p_filter)) {
    plot_df <- plot_df[
      plot_df$mantel_p <= p_filter,
      ,
      drop = FALSE
    ]
  }
  
  if (nrow(plot_df) == 0) {
    stop('No Mantel links remained after filtering.')
  }
  
  ## 4. rcut / pcut
  plot_df$.r_for_cut <- if (isTRUE(r_cut_abs)) {
    abs(plot_df$mantel_r)
  } else {
    plot_df$mantel_r
  }
  
  plot_df$rcut <- cut(
    plot_df$.r_for_cut,
    breaks = r_breaks,
    labels = r_labels,
    include.lowest = TRUE,
    right = TRUE
  )
  
  plot_df$pcut <- cut(
    plot_df$mantel_p,
    breaks = p_breaks,
    labels = p_labels,
    include.lowest = TRUE,
    right = TRUE
  )
  
  plot_df$rcut <- factor(plot_df$rcut, levels = r_labels)
  plot_df$pcut <- factor(plot_df$pcut, levels = p_labels)
  
  plot_df <- plot_df[
    !is.na(plot_df$rcut) & !is.na(plot_df$pcut),
    ,
    drop = FALSE
  ]
  
  if (nrow(plot_df) == 0) {
    stop('No Mantel links remained after rcut/pcut.')
  }
  
  ## 5. spec 分类节点
  ## spec 节点放在与对角线平行的一段线上
  spec_names <- unique(as.character(plot_df$spec))
  n_spec <- length(spec_names)
  
  k <- n_env + 1 + ceiling(spec_offset * n_env)
  spec_y_max <- n_env * spec_range[2]
  spec_y_min <- n_env * spec_range[1]
  spec_x_max <- k - spec_y_max
  spec_x_min <- k - spec_y_min
  
  spec_data <- data.frame(
    spec = spec_names,
    spec_x = seq(spec_x_max, spec_x_min, length.out = n_spec),
    spec_y = seq(spec_y_max, spec_y_min, length.out = n_spec),
    spec_hjust = spec_hjust,
    check.names = FALSE
  )
  
  ## 6. 环境变量对角线锚点
  anchor_data <- data.frame(
    env = env_names,
    anchor_x = seq_len(n_env),
    anchor_y = n_env - seq_len(n_env) + 1,
    check.names = FALSE
  ) |>
    dplyr::filter(env %in% plot_df$env)
  
  ## 7. 生成 Mantel 曲线路径
  curve_df <- plot_df |>
    dplyr::left_join(spec_data, by = 'spec') |>
    dplyr::left_join(anchor_data, by = 'env') |>
    make_curve_path(
      x = 'spec_x',
      y = 'spec_y',
      xend = 'anchor_x',
      yend = 'anchor_y',
      to = 'env',
      to_levels = env_names,
      bend = curve_bend,
      n = curve_n,
      reverse = curve_reverse
    )
  
  ## 8. 环境变量标签
  x_label_df <- data.frame(
    env = env_names,
    x = seq_len(n_env),
    y = -0.01,
    check.names = FALSE
  )
  
  y_label_df <- data.frame(
    env = env_names,
    x = -0.01,
    y = n_env - seq_len(n_env) + 1,
    check.names = FALSE
  )
  
  ## 9. 坐标范围
  x_min <- min(c(0, spec_data$spec_x - 1.8), na.rm = TRUE)
  x_max <- max(c(n_env + 0.5, spec_data$spec_x + 3), na.rm = TRUE)
  y_min <- min(c(0, spec_data$spec_y - 1), na.rm = TRUE)
  y_max <- max(c(n_env + 0.5, spec_data$spec_y + 0.5), na.rm = TRUE)
  
  ## 10. 绘图
  p <- ggplot2::ggplot() +
    ## Mantel 曲线
    ggplot2::geom_path(
      data = curve_df,
      ggplot2::aes( 
        x = x, y = y, group = edge_id, colour = pcut, linewidth = rcut
      ), 
      alpha = 0.9, lineend = 'round'
    ) +
    ## Mantel p / r 图例
    ggplot2::scale_colour_manual(
      values = p_color, drop = FALSE, name = "Mantel's p"
    ) +
    ggplot2::scale_linewidth_manual(
      values = r_size, breaks = r_labels, drop = FALSE, name = "Mantel's r"
    ) +
    ## 固定外框格子
    ggplot2::geom_tile(
      data = cor_data, ggplot2::aes(x = x, y = y),
      width = 1, height = 1, fill = 'white', colour = grid_col,
      linewidth = grid_linewidth
    ) +
    ## 对角线锚点
    ggplot2::geom_point(
      data = anchor_data, ggplot2::aes(x = anchor_x, y = anchor_y),
      shape = 21, fill = node_fill, colour = 'black', 
      size = node_size, stroke = 0.35
    ) +
    ## spec 分类节点
    ggplot2::geom_point(
      data = spec_data, ggplot2::aes(x = spec_x, y = spec_y),
      shape = 21, fill = node_fill, colour = 'black',
      size = node_size, stroke = 0.35
    ) +
    ## spec 标签
    ggplot2::geom_text(
      data = spec_data, 
      ggplot2::aes(x = spec_x, y = spec_y, label = spec, hjust = spec_hjust),
      nudge_x = 0.5, size = spec_label_size
    ) +
    ## 内部相关方块
    ggplot2::geom_point(
      data = cor_data,
      ggplot2::aes(
        x = x, y = y, fill = cor, size = abs(cor)
      ),
      shape = 22, colour = 'white', stroke = 0.25
    ) +
    ggplot2::scale_size_continuous(
      range = square_size_range, guide = 'none'
    ) +
    ## Pearson r 色阶
    ggplot2::scale_fill_gradientn(
      colours = fill_color, limits = c(-1, 1), name = "Pearson's r"
    ) +
    ## x 标签
    ggplot2::geom_text(
      data = x_label_df, ggplot2::aes(x = x, y = y, label = env),
      angle = 90, hjust = 1, vjust = 0.5, nudge_y = .3,
      size = env_label_size, color = 'black'
    ) +
    ## y 标签
    ggplot2::geom_text(
      data = y_label_df, ggplot2::aes(x = x, y = y, label = env),
      hjust = 1, size = env_label_size, color = 'black', nudge_x = .3,
    ) +
    ggplot2::labs(title = title) +
    ggplot2::coord_fixed(
      xlim = c(x_min, x_max), ylim = c(y_min, y_max), clip = 'off'
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = 'bold', hjust = 0.5, size = 13, color = 'black'
      )
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        title = "Mantel's p", override.aes = list(size = 3), order = 1
      ),
      linewidth = ggplot2::guide_legend(
        title = "Mantel's r", override.aes = list(colour = 'black'), order = 2
      ),
      fill = ggplot2::guide_colorbar(
        title = "Pearson's r", order = 3
      )
    )
  
  attr(p, 'cor_data') <- cor_data
  attr(p, 'plot_df') <- plot_df
  attr(p, 'spec_data') <- spec_data
  attr(p, 'anchor_data') <- anchor_data
  attr(p, 'curve_df') <- curve_df
  
  return(p)
}

#### plot_mantel_upper ####
# 绘制 upper 版本 Mantel test + 环境因子相关矩阵图
# test: Mantel test 结果表，至少包含 spec、env、r、p 四列。
# envs: 环境因子表。需要包含 sample_col 指定的样本列。
# spec_range:
#   控制 spec 节点在左侧平行线上的分布范围。
#   例如 c(0.45, 0.95) 表示 y 轴 45% ~ 95%。
# spec_offset:
#   控制 spec 节点所在平行线相对矩阵对角线向外偏移多少。
#   upper 图中，spec_offset 越大，spec 节点越靠左。
# curve_bend: 控制 Mantel 连线弯曲程度。
# curve_reverse: 是否反转曲线方向。

#' Plot Mantel Upper utility
#'
#' `plot_mantel_upper()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 在环境相关矩阵上三角区域叠加 Mantel 关联曲线。
#'
#' @param test Statistical-test result table used for plotting.
#' @param envs Environmental or host-variable table used in Mantel analyses.
#' @param sample_col Name of the sample-identifier column.
#' @param env_cols Environmental-variable columns selected from `envs`.
#' @param spec_col Name of the `spec_col` input column.
#' @param env_col Name of the `env_col` input column.
#' @param r_col Name of the `r_col` input column.
#' @param p_col Name of the `p_col` input column.
#' @param cor_method Correlation method: Pearson, Spearman, or Kendall.
#' @param scale_env Whether environmental variables are standardized before correlation analysis.
#' @param p_filter Maximum P value retained for plotting or downstream analysis.
#' @param r_breaks Numeric breakpoints used to categorize correlation or Mantel r.
#' @param r_labels Labels corresponding to intervals defined by `r_breaks`.
#' @param r_size Numeric setting for `r_size`.
#' @param r_cut_abs Whether absolute r values are used when assigning line-width categories.
#' @param p_breaks Numeric breakpoints used to categorize P values.
#' @param p_labels Labels corresponding to intervals defined by `p_breaks`.
#' @param p_color Color specification for `p_color`.
#' @param fill_color Color specification for `fill_color`.
#' @param grid_col Name of the `grid_col` input column.
#' @param grid_linewidth Numeric setting for `grid_linewidth`.
#' @param square_size_range Numeric setting for `square_size_range`.
#' @param spec_range Relative range reserved for the species-label region.
#' @param spec_offset Offset between species labels and the heatmap boundary.
#' @param spec_hjust Horizontal justification applied to species-side labels.
#' @param curve_bend Signed curvature used for Mantel connection paths.
#' @param curve_n Number of interpolation points used for each connection path.
#' @param curve_reverse Logical control for `curve_reverse`.
#' @param node_fill Color specification for node fill.
#' @param node_size Numeric setting for `node_size`.
#' @param title Optional plot or result title.
#' @param env_label_size Numeric setting for `env_label_size`.
#' @param spec_label_size Numeric setting for `spec_label_size`.
#' @return A plot object; analysis data or models may also be stored as attributes.
#' @export
plot_mantel_upper <- function(
    test, envs, sample_col = 'sample', env_cols = NULL,
    spec_col = 'spec', env_col = 'env', r_col = 'r', p_col = 'p',
    cor_method = c('pearson', 'spearman', 'kendall'),
    scale_env = TRUE, p_filter = 0.05,
    ## Mantel r 分段
    r_breaks = c(0.1, 0.2, 0.3, Inf),
    r_labels = c('0.1 ~ 0.2', '0.2 ~ 0.3', '> 0.3'),
    r_size = c(0.7, 1.4, 2.1),
    r_cut_abs = FALSE, 
    ## Mantel p 分段
    p_breaks = c(-Inf, 0.001, 0.01, 0.05),
    p_labels = c('<= 0.001', '0.001 ~ 0.01', '0.01 ~ 0.05'),
    p_color = c(
      '<= 0.001' = '#d95f02',
      '0.001 ~ 0.01' = '#1f9f79',
      '0.01 ~ 0.05' = 'grey77'
    ),
    ## Pearson 相关矩阵
    fill_color = RColorBrewer::brewer.pal(7, 'PRGn'),
    grid_col = 'black', grid_linewidth = 0.35, square_size_range = c(2, 6),
    ## spec 节点位置
    spec_range = c(.2, .6), spec_offset = 0.5, spec_hjust = 1,
    ## 曲线
    curve_bend = 0.12, curve_n = 80, curve_reverse = FALSE,
    ## node 样式
    node_fill = 'blue', node_size = 2, 
    ## text
    title = 'Mantel test', env_label_size = 3, spec_label_size = 3
) {
  cor_method <- match.arg(cor_method)
  
  test <- data.frame(test, check.names = FALSE)
  envs <- data.frame(envs, check.names = FALSE)
  
  if (!all(c(spec_col, env_col, r_col, p_col) %in% colnames(test))) {
    stop(
      'test should contain columns: ',
      paste(c(spec_col, env_col, r_col, p_col), collapse = ' | ')
    )
  }
  
  if (!sample_col %in% colnames(envs)) {
    stop('envs should contain sample_col: ', sample_col)
  }
  
  if (length(r_labels) != length(r_breaks) - 1) {
    stop('length(r_labels) should be length(r_breaks) - 1.')
  }
  
  if (length(p_labels) != length(p_breaks) - 1) {
    stop('length(p_labels) should be length(p_breaks) - 1.')
  }
  
  if (is.null(names(r_size))) {
    r_size <- stats::setNames(
      rep(r_size, length.out = length(r_labels)),
      r_labels
    )
  }
  
  if (is.null(names(p_color))) {
    p_color <- stats::setNames(
      rep(p_color, length.out = length(p_labels)),
      p_labels
    )
  }
  
  ## 1. 整理环境变量矩阵
  if (is.null(env_cols)) {
    env_cols <- setdiff(colnames(envs), sample_col)
  }
  
  if (!all(env_cols %in% colnames(envs))) {
    stop('Some env_cols are not found in envs.')
  }
  
  env_mat <- envs[, env_cols, drop = FALSE]
  env_mat <- as.matrix(env_mat)
  suppressWarnings(storage.mode(env_mat) <- 'numeric')
  
  keep_env <- apply(env_mat, 2, function(x) {
    x <- x[is.finite(x)]
    length(x) > 2 && stats::sd(x, na.rm = TRUE) > 0
  })
  
  env_mat <- env_mat[, keep_env, drop = FALSE]
  
  if (ncol(env_mat) < 2) {
    stop('Need at least two valid environmental variables.')
  }
  
  if (isTRUE(scale_env)) {
    env_mat <- scale(env_mat)
    env_mat[!is.finite(env_mat)] <- 0
  }
  
  env_names <- colnames(env_mat)
  n_env <- length(env_names)
  
  ## 2. 计算环境变量 Pearson 相关矩阵
  cor_mat <- stats::cor(
    env_mat,
    method = cor_method,
    use = 'pairwise.complete.obs'
  )
  
  cor_long_df <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE)
  
  cor_data <- data.frame(
    env_y = cor_long_df$Var1,
    env_x = cor_long_df$Var2,
    cor = cor_long_df$Freq,
    check.names = FALSE
  )
  
  cor_data$row_id <- match(cor_data$env_y, env_names)
  cor_data$col_id <- match(cor_data$env_x, env_names)
  
  ## upper triangle，不画对角线
  cor_data <- cor_data[
    cor_data$col_id > cor_data$row_id,
    ,
    drop = FALSE
  ]
  
  cor_data$x <- cor_data$col_id
  cor_data$y <- n_env - cor_data$row_id + 1
  
  ## 3. 整理 Mantel 结果
  plot_df <- data.frame(
    spec = test[[spec_col]],
    env = test[[env_col]],
    mantel_r = test[[r_col]],
    mantel_p = test[[p_col]],
    check.names = FALSE
  )
  
  plot_df <- plot_df[
    plot_df$env %in% env_names &
      is.finite(plot_df$mantel_r) &
      is.finite(plot_df$mantel_p),
    ,
    drop = FALSE
  ]
  
  if (!is.null(p_filter)) {
    plot_df <- plot_df[
      plot_df$mantel_p <= p_filter,
      ,
      drop = FALSE
    ]
  }
  
  if (nrow(plot_df) == 0) {
    stop('No Mantel links remained after filtering.')
  }
  
  ## 4. rcut / pcut
  plot_df$.r_for_cut <- if (isTRUE(r_cut_abs)) {
    abs(plot_df$mantel_r)
  } else {
    plot_df$mantel_r
  }
  
  plot_df$rcut <- cut(
    plot_df$.r_for_cut,
    breaks = r_breaks,
    labels = r_labels,
    include.lowest = TRUE,
    right = TRUE
  )
  
  plot_df$pcut <- cut(
    plot_df$mantel_p,
    breaks = p_breaks,
    labels = p_labels,
    include.lowest = TRUE,
    right = TRUE
  )
  
  plot_df$rcut <- factor(plot_df$rcut, levels = r_labels)
  plot_df$pcut <- factor(plot_df$pcut, levels = p_labels)
  
  plot_df <- plot_df[
    !is.na(plot_df$rcut) & !is.na(plot_df$pcut),
    ,
    drop = FALSE
  ]
  
  if (nrow(plot_df) == 0) {
    stop('No Mantel links remained after rcut/pcut.')
  }
  
  ## 5. spec 分类节点
  ## upper 图中，spec 节点放在左侧，与对角线平行
  spec_names <- unique(as.character(plot_df$spec))
  n_spec <- length(spec_names)
  
  k <- n_env + 1 - ceiling(spec_offset * n_env)
  spec_y_max <- n_env * spec_range[2]
  spec_y_min <- n_env * spec_range[1]
  spec_x_min <- k - spec_y_max
  spec_x_max <- k - spec_y_min
  
  spec_data <- data.frame(
    spec = spec_names,
    spec_x = seq(spec_x_min, spec_x_max, length.out = n_spec),
    spec_y = seq(spec_y_max, spec_y_min, length.out = n_spec),
    spec_hjust = spec_hjust,
    check.names = FALSE
  )
  
  ## 6. 环境变量对角线锚点
  anchor_data <- data.frame(
    env = env_names,
    anchor_x = seq_len(n_env),
    anchor_y = n_env - seq_len(n_env) + 1,
    check.names = FALSE
  ) |>
    dplyr::filter(env %in% plot_df$env)
  
  ## 7. 生成 Mantel 曲线路径
  curve_df <- plot_df |>
    dplyr::left_join(spec_data, by = 'spec') |>
    dplyr::left_join(anchor_data, by = 'env') |>
    make_curve_path(
      x = 'spec_x',
      y = 'spec_y',
      xend = 'anchor_x',
      yend = 'anchor_y',
      to = 'env',
      to_levels = env_names,
      bend = curve_bend,
      n = curve_n,
      reverse = curve_reverse
    )
  
  ## 8. 环境变量标签
  x_label_df <- data.frame(
    env = env_names,
    x = seq_len(n_env),
    y = n_env + 0.75,
    check.names = FALSE
  )
  
  y_label_df <- data.frame(
    env = env_names,
    x = n_env + 0.65,
    y = n_env - seq_len(n_env) + 1,
    check.names = FALSE
  )
  
  ## 9. 坐标范围
  x_min <- min(c(0.2, spec_data$spec_x - 1.8), na.rm = TRUE)
  x_max <- max(c(n_env + 1.2, spec_data$spec_x + 1), na.rm = TRUE)
  y_min <- min(c(0.4, spec_data$spec_y - 1), na.rm = TRUE)
  y_max <- max(c(n_env + 1.5, spec_data$spec_y + 0.5), na.rm = TRUE)
  
  ## 10. 绘图
  p <- ggplot2::ggplot() +
    
    ## Mantel 曲线
    ggplot2::geom_path(
      data = curve_df,
      ggplot2::aes(
        x = x,
        y = y,
        group = edge_id,
        colour = pcut,
        linewidth = rcut
      ),
      alpha = 0.9,
      lineend = 'round'
    ) +
    
    ## Mantel p / r 图例
    ggplot2::scale_colour_manual(
      values = p_color,
      drop = FALSE,
      name = "Mantel's p"
    ) +
    ggplot2::scale_linewidth_manual(
      values = r_size,
      breaks = r_labels,
      drop = FALSE,
      name = "Mantel's r"
    ) +
    
    ## 固定外框格子
    ggplot2::geom_tile(
      data = cor_data,
      ggplot2::aes(x = x, y = y),
      width = 1,
      height = 1,
      fill = 'white',
      colour = grid_col,
      linewidth = grid_linewidth
    ) +
    
    ## 对角线锚点
    ggplot2::geom_point(
      data = anchor_data,
      ggplot2::aes(x = anchor_x, y = anchor_y),
      shape = 21,
      fill = node_fill,
      colour = 'black',
      size = node_size,
      stroke = 0.35
    ) +
    
    ## spec 分类节点
    ggplot2::geom_point(
      data = spec_data,
      ggplot2::aes(x = spec_x, y = spec_y),
      shape = 21,
      fill = node_fill,
      colour = 'black',
      size = node_size,
      stroke = 0.35
    ) +
    
    ## spec 标签
    ggplot2::geom_text(
      data = spec_data,
      ggplot2::aes(
        x = spec_x,
        y = spec_y,
        label = spec,
        hjust = spec_hjust
      ),
      nudge_x = -0.5,
      size = spec_label_size
    ) +
    
    ## 内部相关方块
    ggplot2::geom_point(
      data = cor_data,
      ggplot2::aes(
        x = x,
        y = y,
        fill = cor,
        size = abs(cor)
      ),
      shape = 22,
      colour = 'white',
      stroke = 0.25
    ) +
    ggplot2::scale_size_continuous(
      range = square_size_range,
      guide = 'none'
    ) +
    
    ## Pearson r 色阶
    ggplot2::scale_fill_gradientn(
      colours = fill_color,
      limits = c(-1, 1),
      name = "Pearson's r"
    ) +
    
    ## x 标签，上方
    ggplot2::geom_text(
      data = x_label_df,
      ggplot2::aes(x = x, y = y, label = env),
      angle = 90,
      hjust = 0,
      vjust = 0.5,
      size = env_label_size
    ) +
    
    ## y 标签，右侧
    ggplot2::geom_text(
      data = y_label_df,
      ggplot2::aes(x = x, y = y, label = env),
      hjust = 0,
      size = env_label_size
    ) +
    
    ggplot2::labs(title = title) +
    ggplot2::coord_fixed(
      xlim = c(x_min, x_max),
      ylim = c(y_min, y_max),
      clip = 'off'
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = 'bold',
        hjust = 0.5,
        size = 13
      )
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        title = "Mantel's p",
        override.aes = list(size = 3),
        order = 1
      ),
      linewidth = ggplot2::guide_legend(
        title = "Mantel's r",
        override.aes = list(colour = 'grey35'),
        order = 2
      ),
      fill = ggplot2::guide_colorbar(
        title = "Pearson's r",
        order = 3
      )
    )
  
  attr(p, 'cor_data') <- cor_data
  attr(p, 'plot_df') <- plot_df
  attr(p, 'spec_data') <- spec_data
  attr(p, 'anchor_data') <- anchor_data
  attr(p, 'curve_df') <- curve_df
  
  return(p)
}
