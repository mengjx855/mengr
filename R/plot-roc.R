#### Jin-Xin Meng, jinxmeng@zju.edu.cn, 20230816, 20260923 ####

# 20250820: update some function.
# 20260519: simplify code, add helper functions, fix CI label and ribbon order.
# 20260916: standardize script metadata, function sections, documentation, and naming style.
# 20260923: remove Chinese text from Roxygen documentation.



#### theme_roc ####
# ROC 图主题
# style:
#   default : 黑框、无网格，适合常规论文图
#   grid    : 黑框、带浅灰虚线网格
#   classic : 只有坐标轴，无外框，比较简洁
#   minimal : 极简风格，适合展示型图
#   nature  : 类 Nature 风格，坐标轴更突出
#   lancet  : 类 Lancet 风格，黑框较明显
#' Theme Roc utility
#'
#'
#'
#' @param base_size Base font size for the plot theme.
#' @param show_grid Whether to draw panel grid lines.
#' @param style Visual style preset; supported values are shown in Usage.
#' @return A ggplot2 theme object.
#' @export
theme_roc <- function(base_size = 12, show_grid = NULL,
                      style = c(
                        "default", "grid", "classic",
                        "minimal", "nature", "lancet"
                      )) {
  style <- match.arg(style)

  if (is.null(show_grid)) {
    show_grid <- style %in% c("grid", "minimal")
  }

  p <- switch(style,
    default = ggplot2::theme_bw(),
    grid = ggplot2::theme_bw(),
    classic = ggplot2::theme_classic(),
    minimal = ggplot2::theme_minimal(),
    nature = ggplot2::theme_classic(),
    lancet = ggplot2::theme_bw()
  )

  p <- p +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = base_size - 1, color = "black"),
      axis.title = ggplot2::element_text(size = base_size, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = .4, color = "black"),
      axis.ticks.length = grid::unit(2, "mm"),
      axis.line = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        face = "bold", hjust = .5, size = base_size, color = "black"
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = .5, size = base_size - 1, color = "black"
      ),
      legend.text = ggplot2::element_text(size = base_size - 1, color = "black"),
      legend.title = ggplot2::element_text(size = base_size - 1, color = "black"),
      aspect.ratio = 1
    )

  ## 不同主题的细节
  if (style %in% c("default", "grid")) {
    p <- p +
      ggplot2::theme(
        panel.border = ggplot2::element_rect(linewidth = .4, color = "black", fill = NA)
      )
  }

  if (style == "classic") {
    p <- p +
      ggplot2::theme(
        axis.line = ggplot2::element_line(linewidth = .4, color = "black"),
        panel.border = ggplot2::element_blank()
      )
  }

  if (style == "minimal") {
    p <- p +
      ggplot2::theme(
        panel.border = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(linewidth = .3, color = "black")
      )
  }

  if (style == "nature") {
    p <- p +
      ggplot2::theme(
        axis.line = ggplot2::element_line(linewidth = .5, color = "black"),
        axis.ticks = ggplot2::element_line(linewidth = .5, color = "black"),
        panel.border = ggplot2::element_blank(),
        legend.key = ggplot2::element_blank()
      )
  }

  if (style == "lancet") {
    p <- p +
      ggplot2::theme(
        panel.border = ggplot2::element_rect(linewidth = .6, color = "black", fill = NA),
        axis.ticks = ggplot2::element_line(linewidth = .5, color = "black"),
        legend.key = ggplot2::element_blank()
      )
  }

  ## 是否添加网格
  if (isTRUE(show_grid)) {
    p <- p +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(
          linewidth = .35, linetype = "longdash", color = "grey88"
        ),
        panel.grid.minor = ggplot2::element_blank()
      )
  }

  return(p)
}


#### roc_auc_label ####
#' Roc Auc Label utility
#'
#'
#'
#' @param roc A single object returned by `pROC::roc()`.
#' @param digits Optional number of decimal digits retained.
#' @param prefix Prefix used when naming derived coordinates or labels.
#' @return A single character label containing the AUC and its confidence interval.
#' @export
roc_auc_label <- function(roc, digits = 3, prefix = "AUC") {
  auc_value <- as.numeric(pROC::auc(roc))

  auc_ci <- tryCatch(
    as.numeric(pROC::ci.auc(roc)),
    error = function(e) c(NA_real_, NA_real_, NA_real_)
  )

  paste0(
    prefix, ": ", round(auc_value, digits),
    "\n(95% CI: ",
    paste(round(auc_ci[c(1, 3)], digits), collapse = " ~ "),
    ")"
  )
}

#### roc_se_data ####
#' Roc Se Data utility
#'
#'
#'
#' @param roc A single object returned by `pROC::roc()`.
#' @param by Spacing between successive specificity values used to calculate the ROC confidence band.
#' @param conf.level Confidence level used for ROC sensitivity intervals.
#' @return A data frame with specificity, sensitivity, and lower and upper confidence limits.
#' @export
roc_se_data <- function(roc, by = 0.01, conf.level = 0.95) {
  roc_se <- tryCatch(
    pROC::ci.se(
      roc,
      specificities = seq(0, 1, by),
      conf.level = conf.level
    ),
    error = function(e) NULL
  )

  if (is.null(roc_se)) {
    return(NULL)
  }

  roc_se <- data.frame(roc_se, check.names = FALSE) |>
    dplyr::rename(
      lower = dplyr::all_of("2.5%"),
      median = dplyr::all_of("50%"),
      upper = dplyr::all_of("97.5%")
    ) |>
    tibble::rownames_to_column("spec") |>
    dplyr::mutate(spec = as.numeric(spec))

  return(roc_se)
}


#### plot_roc ####
# 针对一个 pROC::roc 对象绘制 ROC 曲线。
# roc: pROC::roc() 返回的对象。
# color: ROC 曲线颜色。
# plot_se: 是否添加 sensitivity 的 95% CI ribbon。
# label_pos: AUC 标签位置，格式为 c(x, y)。

#' Plot Roc utility
#'
#'
#'
#' @param roc A single object returned by `pROC::roc()`.
#' @param color Color specification for `color`.
#' @param plot_se Whether to draw the sensitivity confidence band.
#' @param label_pos Length-two numeric vector giving the x and y position of a label.
#' @param title Optional plot or result title.
#' @param subtitle Optional plot subtitle.
#' @param linewidth Width of the ROC curve.
#' @param se_alpha Opacity of the confidence ribbon.
#' @param theme_style ROC theme preset; supported values are shown in Usage.
#' @param base_size Base font size for the plot theme.
#' @param show_grid Whether to draw panel grid lines.
#' @return A ggplot-compatible plot object; computed data or fitted objects are retained as attributes when applicable.
#' @export
plot_roc <- function(roc, color = "#238443", plot_se = FALSE,
                     label_pos = c(0.75, 0.125), title = NULL,
                     subtitle = NULL, linewidth = .6, se_alpha = .12,
                     theme_style = c(
                       "default", "grid", "classic", "minimal",
                       "nature", "lancet"
                     ), base_size = 12,
                     show_grid = FALSE) {
  theme_style <- match.arg(theme_style)
  label <- roc_auc_label(roc)

  ## 基础 ROC 曲线
  p <- pROC::ggroc(
    roc,
    legacy.axes = TRUE,
    linewidth = linewidth,
    color = color
  ) +
    ggplot2::annotate(
      "segment",
      x = 0, xend = 1, y = 0, yend = 1,
      color = "black",
      linetype = "dashed",
      linewidth = .25
    ) +
    ggplot2::annotate(
      "text",
      x = label_pos[1],
      y = label_pos[2],
      label = label,
      size = 3
    ) +
    ggplot2::scale_x_continuous(expand = c(.01, .01)) +
    ggplot2::scale_y_continuous(expand = c(.01, .01)) +
    ggplot2::labs(
      x = "1 - Specificity",
      y = "Sensitivity",
      title = title,
      subtitle = subtitle
    ) +
    theme_roc(
      base_size = base_size,
      show_grid = show_grid,
      style = theme_style
    )

  ## 添加 sensitivity 置信区间
  ## 注意这里用 p$layers 调整顺序，让 ribbon 位于 ROC 曲线下方。
  if (isTRUE(plot_se)) {
    roc_se <- roc_se_data(roc)

    if (!is.null(roc_se)) {
      ribbon_layer <- ggplot2::geom_ribbon(
        data = roc_se,
        ggplot2::aes(x = 1 - spec, ymin = lower, ymax = upper),
        fill = color,
        alpha = se_alpha,
        inherit.aes = FALSE
      )

      p <- p + ribbon_layer
      p$layers <- c(
        p$layers[length(p$layers)],
        p$layers[-length(p$layers)]
      )
    }
  }

  return(p)
}

#### plot_roc_multiple ####
# 针对多个 pROC::roc 对象绘制 ROC 曲线。
# roc_list: roc 对象列表，建议命名。
# colors: 曲线颜色向量；如果为空则自动生成。
# plot_se: 是否添加 sensitivity 95% CI ribbon
#' Plot Roc Multiple utility
#'
#'
#'
#' @param roc_list Named list of objects returned by `pROC::roc()`.
#' @param colors Color specification for `colors`.
#' @param plot_se Whether to draw the sensitivity confidence band.
#' @param title Optional plot or result title.
#' @param subtitle Optional plot subtitle.
#' @param linewidth Width of the ROC curve.
#' @param se_alpha Opacity of the confidence ribbon.
#' @param theme_style ROC theme preset; supported values are shown in Usage.
#' @param base_size Base font size for the plot theme.
#' @param show_grid Whether to draw panel grid lines.
#' @return A ggplot-compatible plot object; computed data or fitted objects are retained as attributes when applicable.
#' @export
plot_roc_multiple <- function(roc_list, colors = NULL, plot_se = FALSE,
                              title = NULL, subtitle = NULL,
                              linewidth = .6, se_alpha = .10,
                              theme_style = c(
                                "default", "grid", "classic", "minimal",
                                "nature", "lancet"
                              ), base_size = 12,
                              show_grid = TRUE) {
  theme_style <- match.arg(theme_style)
  if (is.null(names(roc_list)) || any(names(roc_list) == "")) {
    names(roc_list) <- paste0("ROC_", seq_along(roc_list))
  }

  if (is.null(colors)) {
    base_colors <- c(
      "#80b1d3", "#b3de69", "#fdb462", "#8dd3c7", "#bc80bd",
      "#fb8072", "#ffed6f", "#fccde5", "#bebada", "#ccebc5"
    )
    colors <- rep(base_colors, length.out = length(roc_list))
  }

  names(colors) <- names(roc_list)

  labels <- purrr::map_chr(
    roc_list, \(x) {
      auc_value <- as.numeric(pROC::auc(x))
      auc_ci <- tryCatch(
        as.numeric(pROC::ci.auc(x)),
        error = function(e) c(NA_real_, NA_real_, NA_real_)
      )

      paste0(
        round(auc_value, 3),
        " [95% CI ",
        paste(round(auc_ci[c(1, 3)], 3), collapse = "~"),
        "]"
      )
    }
  )

  labels <- paste0(names(roc_list), " ", labels)

  p <- pROC::ggroc(
    roc_list,
    legacy.axes = TRUE,
    linewidth = linewidth
  ) +
    ggplot2::annotate(
      "segment",
      x = 0, y = 0, xend = 1, yend = 1,
      color = "black",
      linetype = "longdash",
      linewidth = .35
    ) +
    ggplot2::scale_color_manual(
      values = colors,
      breaks = names(roc_list),
      labels = labels
    ) +
    ggplot2::scale_x_continuous(expand = c(.01, .01)) +
    ggplot2::scale_y_continuous(expand = c(.01, .01)) +
    ggplot2::labs(
      x = "1 - Specificity",
      y = "Sensitivity",
      color = "",
      title = title,
      subtitle = subtitle
    ) +
    theme_roc(
      base_size = base_size,
      show_grid = show_grid,
      style = theme_style
    )

  ## 添加多个 ROC 的 sensitivity 置信区间
  if (isTRUE(plot_se)) {
    roc_se <- purrr::imap_dfr(
      roc_list,
      \(x, y) {
        se_df <- roc_se_data(x)
        if (is.null(se_df)) {
          return(NULL)
        }
        dplyr::mutate(se_df, class = y)
      }
    )

    if (nrow(roc_se) > 0) {
      p <- p +
        ggplot2::geom_ribbon(
          data = roc_se,
          ggplot2::aes(x = 1 - spec, ymin = lower, ymax = upper, fill = class),
          alpha = se_alpha,
          inherit.aes = FALSE,
          show.legend = FALSE
        ) +
        ggplot2::scale_fill_manual(values = colors)

      ## 让 ribbon 位于曲线下面
      p$layers <- c(
        p$layers[length(p$layers)],
        p$layers[-length(p$layers)]
      )
    }
  }

  return(p)
}
