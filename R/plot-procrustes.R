#### Jinxin Meng, 20231215, 20260522, v0.1.1 ####

# 20260522 v0.1.1: update functions.


#### plot_procrustes ####
# Procrustes 分析并绘图
#
# profile_x/profile_y:
#   profile 表，行为 feature，列为 sample。
#   如果输入 profile，会先根据样本名取交集，再计算距离矩阵。
# dist_x/dist_y:
#   已经计算好的距离矩阵，通常为 dist 对象。
#   如果不输入 profile，则使用 dist_x/dist_y。
# dist_method:
#   profile 转距离矩阵时使用的方法，默认 bray。
# symmetric:
#   是否使用对称 Procrustes 分析。
#   TRUE 更适合比较两个矩阵整体一致性；
#   FALSE 更强调将 profile_y 旋转/缩放到 profile_x。
#   默认是 y 旋转变换到 x。
# permutations:
#   PROTEST 置换检验次数。
# colors:
#   x 为目标矩阵颜色，y 为旋转后矩阵颜色。
# 统计量M2:
#   M2 是 Procrustes 统计量，度量两个形状间差异程度的一个指标。
#   M2 是源数据集的点经过旋转、缩放和/或平移后与目标数据集中对应点的平方距离和。
#   它代表了变换后源数据集的点与目标数据集点之间的不匹配程度。
#   M2 的数值越小，表明两组数据集的形状越相似；反之，则表明它们之间的差异越大
# 返回:
#   ggplot 对象；其中 attributes 保存 proc、proc_test 和 plot_df。

#' Plot Procrustes utility
#'
#' `plot_procrustes()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 比较两个 profile 或距离空间，执行 Procrustes 和 permutation test。
#'
#' @param profile_x Feature-by-sample profile used for model training or the first data space.
#' @param profile_y Feature-by-sample profile used for validation or the second data space.
#' @param dist_x Distance object for the first data space.
#' @param dist_y Distance object for the second data space.
#' @param dist_method Distance method; available values are validated with `match.arg()`.
#' @param symmetric Whether to enable the symmetric behavior.
#' @param permutations Number of permutations used by the significance test.
#' @param seed Optional random seed for reproducibility.
#' @param colors Color specification for `colors`.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param title Optional plot or result title.
#' @param subtitle Optional plot subtitle.
#' @param show_grid Logical control for `show_grid`.
#' @param show_line Logical control for `show_line`.
#' @param show_rotation_axis Logical control for `show_rotation_axis`.
#' @param aspect_ratio Panel aspect ratio passed to `ggplot2::theme()`.
#' @param theme Plot theme preset; supported values are shown in Usage.
#' @param ... Additional arguments passed to the underlying function.
#' @return A plot object; analysis data or models may also be stored as attributes.
#' @export
plot_procrustes <- function(profile_x = NULL, profile_y = NULL,
                            dist_x = NULL, dist_y = NULL,
                            dist_method = c(
                              "bray", "jaccard", "euclidean", "manhattan"
                            ), symmetric = TRUE,
                            permutations = 999, seed = 2026,
                            colors = c(x = "#9BBB59", y = "#957DB1"),
                            xlab = NULL, ylab = NULL, title = NULL,
                            subtitle = NULL, show_grid = FALSE,
                            show_line = TRUE, show_rotation_axis = TRUE,
                            aspect_ratio = 3 / 4,
                            theme = c("default", "pubr"),
                            ...) {
  dist_method <- .match_distance_method(dist_method)
  theme <- match.arg(theme)

  ## 1. 准备距离矩阵
  if (!is.null(profile_x) && !is.null(profile_y)) {
    profile_x <- data.frame(profile_x, check.names = FALSE)
    profile_y <- data.frame(profile_y, check.names = FALSE)

    sample_common <- intersect(colnames(profile_x), colnames(profile_y))

    if (length(sample_common) < 3) {
      stop("At least 3 shared samples are required between profile_x and profile_y.")
    }

    profile_x <- profile_x[, sample_common, drop = FALSE]
    profile_y <- profile_y[, sample_common, drop = FALSE]

    dist_x <- vegan::vegdist(t(profile_x), method = dist_method)
    dist_y <- vegan::vegdist(t(profile_y), method = dist_method)

  } else if (!is.null(dist_x) && !is.null(dist_y)) {
    sample_x <- attr(dist_x, "Labels")
    sample_y <- attr(dist_y, "Labels")
    sample_common <- intersect(sample_x, sample_y)

    if (length(sample_common) < 3) {
      stop("At least 3 shared samples are required between dist_x and dist_y.")
    }

    dist_x <- stats::as.dist(as.matrix(dist_x)[sample_common, sample_common])
    dist_y <- stats::as.dist(as.matrix(dist_y)[sample_common, sample_common])

  } else {
    stop("Please provide profile_x/profile_y or dist_x/dist_y.")
  }

  ## 2. PCoA 降维
  PCoA_x <- stats::cmdscale(dist_x, k = 2)
  PCoA_y <- stats::cmdscale(dist_y, k = 2)

  colnames(PCoA_x) <- c("Dim1", "Dim2")
  colnames(PCoA_y) <- c("Dim1", "Dim2")

  ## 3. Procrustes 分析
  ## profile_x 是目标矩阵，profile_y 会被旋转匹配到 profile_x
  proc <- vegan::procrustes(PCoA_x, PCoA_y, symmetric = symmetric)

  set.seed(seed)
  proc_test <- vegan::protest(
    PCoA_x, PCoA_y,
    permutations = permutations,
    symmetric = symmetric
  )

  ## 4. 提取绘图数据
  proc_point <- data.frame(
    sample = rownames(proc$X),
    X1_target = proc$X[, 1],
    X2_target = proc$X[, 2],
    X1_rotated = proc$Yrot[, 1],
    X2_rotated = proc$Yrot[, 2],
    check.names = FALSE
  )

  proc_coord <- data.frame(proc$rotation, check.names = FALSE)

  # plot(proc, kind = 2)
  # residuals(proc)

  ## 5. 标签
  if (is.null(xlab)) xlab <- "Dim 1"
  if (is.null(ylab)) ylab <- "Dim 2"

  if (is.null(title)) {
    title <- paste0(tools::toTitleCase(dist_method), " distance-based Procrustes analysis")
  }

  if (is.null(subtitle)) {
    subtitle <- substitute(
      M^2 == a ~ ", " ~ italic(p) == b,
      list(
        a = round(proc_test$ss, 4),
        b = signif(proc_test$signif, 3)
      )
    )
  }

  col_x <- unname(colors["x"])
  col_y <- unname(colors["y"])

  if (is.na(col_x)) col_x <- colors[1]
  if (is.na(col_y)) col_y <- colors[2]

  ## 6. 绘图
  p <- ggplot2::ggplot(proc_point) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = X1_rotated, y = X2_rotated,
        xend = (X1_rotated + X1_target) / 2,
        yend = (X2_rotated + X2_target) / 2
      ),
      color = col_y,
      linewidth = .4
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = (X1_rotated + X1_target) / 2,
        y = (X2_rotated + X2_target) / 2,
        xend = X1_target,
        yend = X2_target
      ),
      arrow = grid::arrow(length = grid::unit(0.15, "cm")),
      color = col_x,
      linewidth = .4
    ) +
    ggplot2::geom_point(
      ggplot2::aes(X1_rotated, X2_rotated),
      color = col_y,
      size = 1.8
    ) +
    ggplot2::geom_point(
      ggplot2::aes(X1_target, X2_target),
      color = col_x,
      size = 1.8
    ) +
    ggplot2::labs(
      x = xlab,
      y = ylab,
      title = title,
      subtitle = subtitle
    )

  ## 旋转坐标轴辅助线
  if (isTRUE(show_rotation_axis)) {
    
    slope1 <- proc_coord[1, 2] / proc_coord[1, 1]
    slope2 <- proc_coord[2, 2] / proc_coord[2, 1]

    if (is.finite(slope1)) {
      p <- p + ggplot2::geom_abline(intercept = 0, slope = slope1, linewidth = .4)
    }

    if (is.finite(slope2)) {
      p <- p + ggplot2::geom_abline(intercept = 0, slope = slope2, linewidth = .4)
    }
  }

  ## 中心辅助线
  if (isTRUE(show_line)) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = 0,
        color = "gray70",
        linetype = "longdash",
        linewidth = .4
      ) +
      ggplot2::geom_hline(
        yintercept = 0,
        color = "gray70",
        linetype = "longdash",
        linewidth = .4
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
        plot.subtitle = ggplot2::element_text(size = 12, color = "black"),
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

  if (isTRUE(show_grid)) {
    p <- p +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(linewidth = .4, color = "grey90"),
        panel.grid.minor = ggplot2::element_blank()
      )
  }

  attr(p, "proc") <- proc
  attr(p, "proc_test") <- proc_test
  attr(p, "plot_df") <- proc_point

  return(p)
}
