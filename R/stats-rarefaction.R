#### Jin-Xin Meng, jinxmeng@zju.edu.cn, 20220529, 20260923 ####

# 20260916: standardize script metadata, function sections, documentation, and naming style.
# 20260923: clarify metadata argument names and remove Chinese text from Roxygen documentation.


#### calcu_specaccum ####

#' Calculate a species accumulation curve
#'
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param permutations Number of permutations used by the significance test.
#' @param method Analysis or summary method; supported values are shown in the usage.
#' @return A data frame containing accumulation effort, richness, standard error, and confidence limits.
#' @export
calcu_specaccum <- function(
  profile, permutations = 99,
  method = c("random", "collector", "exact", "rarefaction", "coleman")
) {
  ## 1. 整理为 sample × feature 矩阵
  method <- match.arg(method)
  profile_df <- .as_profile_df(profile, numeric = TRUE)

  ## 2. 计算累积曲线并整理输出
  specaccum_obj <- vegan::specaccum(
    t(profile_df),
    method = method, permutations = permutations
  )
  sd_vec <- specaccum_obj$sd
  if (is.null(sd_vec)) sd_vec <- rep(NA_real_, length(specaccum_obj$sites))
  data.frame(
    sample_n = specaccum_obj$sites,
    richness = specaccum_obj$richness,
    sd = sd_vec,
    check.names = FALSE
  )
}

#### plot_specaccum ####

#' Plot a species accumulation curve
#'
#'
#' @param data An input data frame or compatible object.
#' @param sample_n_col Name of the `sample_n_col` input column.
#' @param richness_col Name of the `richness_col` input column.
#' @param sd_col Name of the `sd_col` input column.
#' @param color Color specification for `color`.
#' @param add_errorbar Whether to draw pointwise error bars.
#' @param add_ribbon Whether to draw an uncertainty ribbon.
#' @param fill Color specification for fill.
#' @param linetype Line type used for the accumulation or rarefaction curve.
#' @param aspect_ratio Panel aspect ratio passed to `ggplot2::theme()`.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param title Optional plot or result title.
#' @return A ggplot-compatible plot object; computed data or fitted objects are retained as attributes when applicable.
#' @export
plot_specaccum <- function(
  data, sample_n_col = "sample_n", richness_col = "richness", sd_col = "sd",
  color = "#de2726", add_errorbar = TRUE, add_ribbon = FALSE,
  fill = "#fcbba1", linetype = "solid", aspect_ratio = 1,
  xlab = "Cumulative samples", ylab = "Cumulative features",
  title = "Rarefaction curve analysis"
) {
  ## 1. 统一作图列
  plot_df <- .as_df(data)
  .check_columns(
    plot_df, c(sample_n_col, richness_col, sd_col),
    object = "data"
  )
  plot_df <- data.frame(
    sample_n = plot_df[[sample_n_col]],
    richness = plot_df[[richness_col]],
    sd = plot_df[[sd_col]], check.names = FALSE
  )
  if (isTRUE(add_errorbar) && isTRUE(add_ribbon)) add_errorbar <- FALSE

  ## 2. 按需要添加误差线或 ribbon
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(sample_n, richness))
  if (isTRUE(add_errorbar)) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = richness - sd, ymax = richness + sd),
      width = 0.6, linewidth = 0.4, color = color, show.legend = FALSE
    )
  }
  if (isTRUE(add_ribbon)) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = richness - sd, ymax = richness + sd),
      fill = fill, show.legend = FALSE
    )
  }

  ## 3. 统一主题并返回
  p +
    ggplot2::geom_line(
      linewidth = 0.4, color = color, linetype = linetype,
      show.legend = FALSE
    ) +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    ggplot2::scale_x_continuous(expand = c(0.02, 0.02)) +
    ggplot2::scale_y_continuous(expand = c(0.02, 0.02)) +
    ggpubr::theme_pubr() +
    ggplot2::theme(aspect.ratio = aspect_ratio)
}

#### calcu_specaccum_by_group ####

#' Calculate species accumulation curves by group
#'
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param sample_meta A sample metadata table containing sample and group columns.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param group_level Optional order of group levels.
#' @param permutations Number of permutations used by the significance test.
#' @param method Analysis or summary method; supported values are shown in the usage.
#' @return A data frame of accumulation curves with a group column.
#' @export
calcu_specaccum_by_group <- function(
  profile, sample_meta, sample_col = "sample", group_col = "group",
  group_level = NULL, permutations = 99,
  method = c("random", "collector", "exact", "rarefaction", "coleman")
) {
  group <- sample_meta
  ## 1. 对齐 profile 与 group
  method <- match.arg(method)
  aligned <- .align_profile_group(
    profile = profile, sample_meta = group,
    sample_col = sample_col, group_col = group_col,
    group_level = group_level
  )
  profile_df <- data.frame(t(aligned$profile_df), check.names = FALSE)
  group_df <- aligned$group_df

  ## 2. 每组独立计算累积曲线
  purrr::map_dfr(aligned$group_level, function(group_name) {
    sample_vec <- group_df[[sample_col]][
      as.character(group_df[[group_col]]) == group_name
    ]
    specaccum_obj <- vegan::specaccum(
      profile_df[sample_vec, , drop = FALSE],
      method = method, permutations = permutations
    )
    sd_vec <- specaccum_obj$sd
    if (is.null(sd_vec)) sd_vec <- rep(NA_real_, length(specaccum_obj$sites))
    data.frame(
      sample_n = specaccum_obj$sites,
      richness = specaccum_obj$richness,
      sd = sd_vec,
      group = group_name, check.names = FALSE
    )
  })
}

#### plot_specaccum_by_group ####

#' Plot species accumulation curves by group
#'
#'
#' @param data An input data frame or compatible object.
#' @param sample_n_col Name of the `sample_n_col` input column.
#' @param richness_col Name of the `richness_col` input column.
#' @param sd_col Name of the `sd_col` input column.
#' @param group_col Name of the grouping column.
#' @param group_level Optional order of group levels.
#' @param group_color Optional colors aligned to `group_level`.
#' @param add_errorbar Whether to draw pointwise error bars.
#' @param add_group_label Whether to label group centroids.
#' @param add_ribbon Whether to draw an uncertainty ribbon.
#' @param fill Color specification for fill.
#' @param aspect_ratio Panel aspect ratio passed to `ggplot2::theme()`.
#' @param linetype Line type used for the accumulation or rarefaction curve.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param title Optional plot or result title.
#' @return A ggplot-compatible plot object; computed data or fitted objects are retained as attributes when applicable.
#' @export
plot_specaccum_by_group <- function(
  data, sample_n_col = "sample_n", richness_col = "richness", sd_col = "sd",
  group_col = "group", group_level = NULL, group_color = NULL,
  add_errorbar = TRUE, add_group_label = FALSE, add_ribbon = FALSE,
  fill = "grey85", aspect_ratio = 1, linetype = "solid",
  xlab = "Number of samples", ylab = "Number of features",
  title = "Rarefaction curve analysis"
) {
  ## 1. 统一作图列和分组顺序
  plot_df <- .as_df(data)
  .check_columns(
    plot_df,
    c(sample_n_col, richness_col, sd_col, group_col),
    object = "data"
  )
  plot_df <- data.frame(
    sample_n = plot_df[[sample_n_col]], richness = plot_df[[richness_col]],
    sd = plot_df[[sd_col]], group = plot_df[[group_col]],
    check.names = FALSE
  )
  if (is.null(group_level)) group_level <- unique(as.character(plot_df$group))
  plot_df$group <- factor(plot_df$group, levels = group_level)
  group_color <- .resolve_group_colors(group_level, group_color)
  if (isTRUE(add_errorbar) && isTRUE(add_ribbon)) add_errorbar <- FALSE

  ## 2. 绘制曲线及误差范围
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(sample_n, richness, color = group, group = group)
  )
  if (isTRUE(add_errorbar)) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = richness - sd, ymax = richness + sd),
      width = 0.3, linewidth = 0.4, show.legend = FALSE
    )
  }
  if (isTRUE(add_ribbon)) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = richness - sd, ymax = richness + sd),
      fill = fill, color = NA, show.legend = FALSE
    )
  }
  p <- p +
    ggplot2::geom_line(linewidth = 0.4, linetype = linetype) +
    ggplot2::scale_color_manual(values = group_color) +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    ggplot2::scale_x_continuous(expand = c(0.02, 0.02)) +
    ggplot2::scale_y_continuous(expand = c(0.02, 0.02)) +
    ggpubr::theme_pubr() +
    ggplot2::theme(aspect.ratio = aspect_ratio)

  ## 3. 可选地在曲线末端添加组名
  if (isTRUE(add_group_label)) {
    label_df <- plot_df |>
      dplyr::group_by(group) |>
      dplyr::slice_max(sample_n, n = 1, with_ties = FALSE) |>
      dplyr::ungroup()
    p <- p +
      ggplot2::geom_label(
        data = label_df,
        ggplot2::aes(sample_n, richness, label = group),
        inherit.aes = FALSE, size = 2
      ) +
      ggplot2::guides(color = "none")
  }
  p
}

#### calcu_specaccum_by_depth ####

#' Calculate sample rarefaction curves by sequencing depth
#'
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param step Increment in sequencing depth between rarefaction points.
#' @param seed Optional random seed for reproducibility.
#' @return A data frame of sample, subsampling depth, and observed feature richness.
#' @export
calcu_specaccum_by_depth <- function(profile, step = 1000, seed = NULL) {
  ## 1. 检查深度并整理 count matrix
  if (!is.numeric(step) || length(step) != 1L || step <= 0) {
    stop("step should be a single positive number.")
  }
  if (!is.null(seed)) set.seed(seed)
  profile_mat <- t(as.matrix(.as_profile_df(profile, numeric = TRUE)))
  if (any(profile_mat < 0) || any(profile_mat %% 1 != 0)) {
    stop("profile should contain non-negative integer counts.")
  }
  depth_vec <- rowSums(profile_mat)
  depth_grid <- seq(step, max(depth_vec), by = step)
  if (!length(depth_grid)) depth_grid <- max(depth_vec)
  if (max(depth_grid) < max(depth_vec)) depth_grid <- c(depth_grid, max(depth_vec))

  ## 2. 每个深度只抽平测序量足够的样本
  result_list <- vector("list", length(depth_grid))
  for (depth_idx in seq_along(depth_grid)) {
    sample_depth <- depth_grid[depth_idx]
    keep_sample <- depth_vec >= sample_depth
    rarefied_mat <- vegan::rrarefy(
      profile_mat[keep_sample, , drop = FALSE], sample_depth
    )
    richness_vec <- vegan::estimateR(rarefied_mat)[1, ]
    result_list[[depth_idx]] <- data.frame(
      sample = names(richness_vec), richness = as.numeric(richness_vec),
      depth = sample_depth, check.names = FALSE
    )
  }
  dplyr::bind_rows(
    data.frame(
      sample = rownames(profile_mat), richness = 0, depth = 0,
      check.names = FALSE
    ),
    dplyr::bind_rows(result_list)
  )
}

#### plot_specaccum_by_depth ####

#' Plot sample rarefaction curves by sequencing depth
#'
#'
#' @param data An input data frame or compatible object.
#' @param sample_col Name of the sample-identifier column.
#' @param depth_col Name of the `depth_col` input column.
#' @param richness_col Name of the `richness_col` input column.
#' @return A ggplot-compatible plot object; computed data or fitted objects are retained as attributes when applicable.
#' @export
plot_specaccum_by_depth <- function(
  data, sample_col = "sample", depth_col = "depth",
  richness_col = "richness"
) {
  ## 统一作图列
  plot_df <- .as_df(data)
  .check_columns(plot_df, c(sample_col, depth_col, richness_col), "data")
  plot_df <- data.frame(
    sample = plot_df[[sample_col]], depth = plot_df[[depth_col]],
    richness = plot_df[[richness_col]], check.names = FALSE
  )
  color_n <- length(unique(plot_df$sample))

  ## 绘制每个样本的稀释曲线
  ggplot2::ggplot(plot_df, ggplot2::aes(depth, richness, color = sample)) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(
      x = "Sequences per sample", y = "Observed features", color = "Sample"
    ) +
    ggplot2::scale_color_manual(values = palc("Rainbow5", n = color_n)) +
    ggplot2::scale_x_continuous(
      labels = scales::label_number(scientific = FALSE)
    ) +
    ggpubr::theme_pubr() +
    ggplot2::theme(
      aspect.ratio = 3 / 4, axis.line = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.ticks.length = grid::unit(2, "mm"),
      legend.position = "right",
      panel.border = ggplot2::element_rect(
        fill = NA, linewidth = 0.5, colour = "black"
      )
    )
}

#### calcu_rankabund ####

#' Calculate rank-abundance data for each sample
#'
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @return A data frame containing sample-wise abundance ranks, abundances, and log abundances.
#' @export
calcu_rankabund <- function(profile) {
  ## 每列作为一个样本传给 BiodiversityR::rankabundance()
  profile_df <- .as_profile_df(profile, numeric = TRUE)
  purrr::map2_dfr(profile_df, colnames(profile_df), function(value_vec, sample_name) {
    sample_mat <- matrix(
      value_vec,
      nrow = 1,
      dimnames = list(sample_name, rownames(profile_df))
    )
    BiodiversityR::rankabundance(sample_mat) |>
      data.frame(row.names = NULL, check.names = FALSE) |>
      dplyr::filter(abundance != 0) |>
      dplyr::transmute(
        sample = sample_name, rank = rank, log_abundance = logabun
      )
  })
}

#### plot_rankabund ####

#' Plot rank-abundance curves
#'
#'
#' @param data An input data frame or compatible object.
#' @param sample_col Name of the sample-identifier column.
#' @param rank_col Name of the `rank_col` input column.
#' @param log_abundance_col Name of the `log_abundance_col` input column.
#' @return A ggplot-compatible plot object; computed data or fitted objects are retained as attributes when applicable.
#' @export
plot_rankabund <- function(
  data, sample_col = "sample", rank_col = "rank",
  log_abundance_col = "log_abundance"
) {
  ## 统一作图列
  plot_df <- .as_df(data)
  .check_columns(
    plot_df, c(sample_col, rank_col, log_abundance_col),
    object = "data"
  )
  plot_df <- data.frame(
    sample = plot_df[[sample_col]], rank = plot_df[[rank_col]],
    log_abundance = plot_df[[log_abundance_col]], check.names = FALSE
  )
  color_n <- length(unique(plot_df$sample))

  ## 绘制 rank-abundance 曲线
  ggplot2::ggplot(
    plot_df, ggplot2::aes(rank, log_abundance, color = sample)
  ) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::labs(
      x = "Feature rank", y = "log10 relative abundance (%)", color = NULL
    ) +
    ggplot2::scale_color_manual(values = palc("Rainbow5", n = color_n)) +
    ggplot2::theme(
      aspect.ratio = 3 / 4, axis.line = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 11, color = "black"),
      axis.ticks.length = grid::unit(2, "mm"),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        fill = NA, linewidth = 0.5, colour = "black"
      ),
      panel.background = ggplot2::element_rect(
        fill = "transparent", color = "black"
      ),
      legend.key = ggplot2::element_rect(fill = NA, color = NA)
    )
}
