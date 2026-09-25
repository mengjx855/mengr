#### Jin-Xin Meng, jinxmeng@zju.edu.cn, 20220927, 20260923 ####

# 20260916: rename `plot_dbRDA()` to `plot_dbrda()` and standardize documentation and naming.
# 20260923: clarify metadata argument names and remove Chinese text from Roxygen documentation.


#### plot_dbrda ####

#' Plot distance-based redundancy analysis
#'
#' Supply either a feature-by-sample `profile` or a precomputed `distance`.
#' The columns in `constraint_cols` are used as explanatory variables, while
#' `group_col` controls the plot colour.
#'
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param sample_meta A sample metadata table containing sample and group columns.
#' @param distance A precomputed distance object; when supplied, it takes precedence over `profile`.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param constraint_cols Metadata columns included as constrained variables in dbRDA.
#' @param group_level Optional order of group levels.
#' @param group_color Optional colors aligned to `group_level`.
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
#' @param show_variable Whether to draw significant constraining-variable arrows.
#' @param show_line Whether to draw horizontal and vertical reference lines at zero.
#' @param aspect_ratio Panel aspect ratio passed to `ggplot2::theme()`.
#' @param theme Plot theme preset; supported values are shown in Usage.
#' @param permutations Number of permutations used by the significance test.
#' @param ... Additional arguments passed to distance calculation and `plot_dim()`.
#' @return A ggplot-compatible plot object; computed data or fitted objects are retained as attributes when applicable.
#' @export
plot_dbrda <- function(
  profile = NULL, sample_meta, distance = NULL,
  sample_col = "sample", group_col = "group",
  constraint_cols = group_col, group_level = NULL, group_color = NULL,
  dist_method = c(
    "bray", "jaccard", "euclidean", "manhattan", "canberra",
    "kulczynski", "gower", "altGower", "morisita", "horn",
    "mountford", "raup", "binomial", "chao", "cao",
    "mahalanobis", "unifrac"
  ),
  transform = c(
    "hellinger", "total", "max", "frequency", "normalize",
    "range", "rank", "rrank", "standardize", "pa",
    "chi.square", "log", "clr", "rclr", "alr"
  ),
  display_type = c("line", "point"),
  conf_type = c("ellipse", "encircle", "none"), ellipse_level = 0.75,
  title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
  legend_title = "Group", add_group_label = FALSE,
  add_sample_label = FALSE, label_size = 3, point_size = 1.5,
  show_legend = TRUE, show_grid = FALSE, show_variable = TRUE,
  show_line = TRUE, aspect_ratio = 3 / 4,
  theme = c("default", "pubr"), permutations = 999, ...
) {
  group <- sample_meta
  dist_method <- .match_distance_method(dist_method)
  transform <- .match_transform_method(transform)
  display_type <- match.arg(display_type)
  conf_type <- match.arg(conf_type)
  theme <- match.arg(theme)

  group_df <- .as_df(group)
  .check_columns(
    group_df,
    unique(c(sample_col, group_col, constraint_cols)),
    object = "group"
  )
  group_df[[sample_col]] <- as.character(group_df[[sample_col]])
  if (anyDuplicated(group_df[[sample_col]])) {
    stop("Duplicated sample identifiers found in group[[sample_col]].")
  }

  if (is.null(distance)) {
    if (is.null(profile)) stop("Supply either profile or distance.")
    aligned <- .align_profile_group(
      profile = profile, sample_meta = group_df,
      sample_col = sample_col, group_col = group_col,
      group_level = group_level
    )
    group_df <- aligned$group_df
    distance <- calcu_distance(
      profile = aligned$profile_df,
      dist_method = dist_method,
      transform = transform,
      ...
    )
  } else {
    dist_mat <- as.matrix(distance)
    sample_vec <- intersect(rownames(dist_mat), group_df[[sample_col]])
    if (length(sample_vec) < 3L) {
      stop("dbRDA requires at least three matched samples.")
    }
    group_df <- group_df[
      match(sample_vec, group_df[[sample_col]]), ,
      drop = FALSE
    ]
    distance <- stats::as.dist(dist_mat[sample_vec, sample_vec, drop = FALSE])
  }

  if (is.null(group_level)) {
    group_level <- if (is.factor(group_df[[group_col]])) {
      levels(droplevels(group_df[[group_col]]))
    } else {
      unique(as.character(group_df[[group_col]]))
    }
  }
  group_color <- .resolve_group_colors(group_level, group_color)

  model_df <- group_df[, constraint_cols, drop = FALSE]
  rownames(model_df) <- group_df[[sample_col]]
  dbRDA_obj <- vegan::capscale(distance ~ ., data = model_df)
  anova_df <- stats::anova(dbRDA_obj, permutations = permutations)

  site_mat <- vegan::scores(dbRDA_obj, display = "sites", scaling = 1)
  if (ncol(site_mat) < 2L) {
    site_mat <- cbind(site_mat, Axis2 = 0)
  }
  axis_names <- colnames(site_mat)[1:2]
  plot_df <- data.frame(
    sample = rownames(site_mat),
    X1 = site_mat[, 1],
    X2 = site_mat[, 2],
    check.names = FALSE
  ) |>
    dplyr::left_join(
      data.frame(
        sample = group_df[[sample_col]],
        group = group_df[[group_col]],
        check.names = FALSE
      ),
      by = "sample"
    ) |>
    dplyr::mutate(group = factor(group, levels = group_level))

  eig_vec <- dbRDA_obj$CCA$eig
  eig_pct <- if (length(eig_vec)) {
    100 * eig_vec / sum(c(dbRDA_obj$CCA$eig, dbRDA_obj$CA$eig), na.rm = TRUE)
  } else {
    numeric()
  }
  if (is.null(xlab)) {
    xlab <- if (length(eig_pct) >= 1L) {
      sprintf("%s (%.2f%%)", axis_names[1], eig_pct[1])
    } else {
      axis_names[1]
    }
  }
  if (is.null(ylab)) {
    ylab <- if (length(eig_pct) >= 2L) {
      sprintf("%s (%.2f%%)", axis_names[2], eig_pct[2])
    } else {
      axis_names[2]
    }
  }
  if (is.null(title)) {
    title <- paste0(stringr::str_to_sentence(dist_method), "-distance dbRDA")
  }
  if (is.null(subtitle)) {
    adj_r2 <- vegan::RsquareAdj(dbRDA_obj)$adj.r.squared
    p_value <- anova_df[1, "Pr(>F)"]
    subtitle <- sprintf("Adjusted R2 = %.3f, p = %.3g", adj_r2, p_value)
  }

  p <- plot_dim(
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

  if (isTRUE(show_variable)) {
    centroid_mat <- tryCatch(
      vegan::scores(dbRDA_obj, display = "cn", choices = 1:2, scaling = 1),
      error = function(e) NULL
    )
    if (!is.null(centroid_mat) && nrow(centroid_mat)) {
      variable_df <- data.frame(
        variable = rownames(centroid_mat),
        X1 = centroid_mat[, 1],
        X2 = centroid_mat[, 2],
        check.names = FALSE
      )
      p <- p +
        ggplot2::geom_segment(
          data = variable_df,
          ggplot2::aes(x = 0, y = 0, xend = X1, yend = X2),
          inherit.aes = FALSE,
          arrow = grid::arrow(length = grid::unit(0.8, "mm")),
          linewidth = 0.4, color = "black"
        ) +
        ggrepel::geom_text_repel(
          data = variable_df,
          ggplot2::aes(x = X1, y = X2, label = variable),
          inherit.aes = FALSE, size = label_size, color = "black"
        )
    }
  }

  attr(p, "dbRDA") <- dbRDA_obj
  attr(p, "anova") <- anova_df
  p
}
