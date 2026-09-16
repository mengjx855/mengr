#### Jin-Xin Meng, jinxmeng@zju.edu.cn, 20241204, 20260916 ####

# 20260916: standardize script metadata, function sections, documentation, and naming style.

#### plot_umap ####

#' Plot a UMAP ordination
#'
#' 对齐样本分组后执行 UMAP，并返回统一风格降维图。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param group A sample metadata table containing sample and group columns.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param group_level Optional order of group levels.
#' @param group_color Optional colors aligned to `group_level`.
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
#' @param ... Additional arguments passed to `umap::umap()`.
#' @return A ggplot-compatible plot object; computed data or fitted objects are retained as attributes when applicable.
#' @export
plot_umap <- function(
  profile, group, sample_col = "sample", group_col = "group",
  group_level = NULL, group_color = NULL,
  display_type = c("line", "point"),
  conf_type = c("ellipse", "encircle", "none"), ellipse_level = 0.75,
  title = NULL, subtitle = NULL, xlab = "UMAP_1", ylab = "UMAP_2",
  legend_title = "Group", add_group_label = FALSE,
  add_sample_label = FALSE, label_size = 1.5, point_size = 1.5,
  show_legend = TRUE, show_grid = FALSE, show_line = TRUE,
  aspect_ratio = 3 / 4, theme = c("default", "pubr"), ...
) {
  display_type <- match.arg(display_type)
  conf_type <- match.arg(conf_type)
  theme <- match.arg(theme)

  aligned <- .align_profile_group(
    profile = profile,
    group = group,
    sample_col = sample_col,
    group_col = group_col,
    group_level = group_level
  )
  profile_df <- aligned$profile_df
  group_level <- aligned$group_level
  group_color <- .resolve_group_colors(group_level, group_color)

  umap_obj <- umap::umap(t(as.matrix(profile_df)), ...)
  group_key_df <- data.frame(
    sample = aligned$group_df[[sample_col]],
    group = aligned$group_df[[group_col]],
    check.names = FALSE
  )
  plot_df <- data.frame(
    sample = aligned$sample_vec,
    X1 = umap_obj$layout[, 1],
    X2 = umap_obj$layout[, 2],
    check.names = FALSE
  ) |>
    dplyr::left_join(group_key_df, by = "sample") |>
    dplyr::mutate(group = factor(group, levels = group_level))

  if (is.null(title)) {
    title <- "Uniform manifold approximation and projection analysis"
  }

  plot_dim(
    data = plot_df, group_level = group_level, group_color = group_color,
    display_type = display_type, conf_type = conf_type,
    ellipse_level = ellipse_level, title = title, subtitle = subtitle,
    xlab = xlab, ylab = ylab, legend_title = legend_title,
    add_group_label = add_group_label,
    add_sample_label = add_sample_label, label_size = label_size,
    point_size = point_size, show_legend = show_legend,
    show_grid = show_grid, show_line = show_line,
    aspect_ratio = aspect_ratio, theme = theme
  )
}
