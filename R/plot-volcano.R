#### Jin-Xin Meng, 20240328, 20260820, v0.2.0 ####

#' Plot a volcano plot
#'
#' Chinese summary: 根据 effect、P 值和富集方向绘制 volcano plot。
#'
#' @param data An input data frame or compatible object.
#' @param feature_col Name of the feature-identifier column.
#' @param effect_col Name of the `effect_col` input column.
#' @param p_col Name of the `p_col` input column.
#' @param enriched_col Name of the `enriched_col` input column.
#' @param group_level Optional order of group levels.
#' @param group_color Optional colors aligned to `group_level`.
#' @param p_cutoff Numeric setting for `p_cutoff`.
#' @param effect_cutoff Numeric setting for `effect_cutoff`.
#' @param x_limit Optional numeric limits for the x axis.
#' @param title Optional plot or result title.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param aspect_ratio Panel aspect ratio passed to `ggplot2::theme()`.
#' @param subtitle_keywords Keywords whose occurrences are highlighted in the plot subtitle.
#' @return A plot object; analysis data or models may also be stored as attributes.
#' @export
plot_volcano <- function(
    data, feature_col = 'name', effect_col = 'log2FC', p_col = 'pval',
    enriched_col = 'enriched',
    group_level = c('Up', 'None', 'Down'), group_color = NULL,
    p_cutoff = 0.05, effect_cutoff = 1, x_limit = NULL,
    title = 'Volcano plot',
    xlab = expression(log[2] * FoldChange),
    ylab = expression(-log[10] * pvalue),
    aspect_ratio = 3 / 4,
    subtitle_keywords = c('FoldChange', 'pvalue')) {
  ## 1. 统一列名并检查阈值
  plot_df <- .as_df(data)
  .check_columns(
    plot_df, c(feature_col, effect_col, p_col), object = 'data'
  )
  if (!is.numeric(p_cutoff) || length(p_cutoff) != 1L ||
      p_cutoff <= 0 || p_cutoff >= 1) {
    stop('p_cutoff should be a single number between 0 and 1.')
  }
  if (!is.numeric(effect_cutoff) || length(effect_cutoff) != 1L ||
      effect_cutoff < 0) {
    stop('effect_cutoff should be a single non-negative number.')
  }
  plot_df <- dplyr::rename(
    plot_df,
    name = dplyr::all_of(feature_col),
    effect = dplyr::all_of(effect_col),
    pvalue = dplyr::all_of(p_col)
  )
  plot_df$effect <- suppressWarnings(as.numeric(plot_df$effect))
  plot_df$pvalue <- suppressWarnings(as.numeric(plot_df$pvalue))

  ## 2. 根据阈值生成或读取富集方向
  if (enriched_col %in% colnames(plot_df)) {
    plot_df$enriched <- as.character(plot_df[[enriched_col]])
  } else {
    plot_df$enriched <- ifelse(
      plot_df$pvalue < p_cutoff & plot_df$effect >= effect_cutoff,
      'Up',
      ifelse(
        plot_df$pvalue < p_cutoff & plot_df$effect <= -effect_cutoff,
        'Down', 'None'
      )
    )
  }
  plot_df <- plot_df[
    is.finite(plot_df$effect) & is.finite(plot_df$pvalue) &
      plot_df$pvalue >= 0,
    , drop = FALSE
  ]
  if (!nrow(plot_df)) stop('No valid rows remained for plotting.')

  ## p = 0 时用最小正 p 值的十分之一替代，避免 -log10(p) 为 Inf
  positive_p <- plot_df$pvalue[plot_df$pvalue > 0]
  zero_fill <- if (length(positive_p)) min(positive_p) / 10 else .Machine$double.xmin
  plot_df$plot_pvalue <- pmax(plot_df$pvalue, zero_fill)

  ## 3. 可选地截断横轴极端值
  if (!is.null(x_limit)) {
    if (length(x_limit) == 1L) x_limit <- c(-abs(x_limit), abs(x_limit))
    if (length(x_limit) != 2L || x_limit[1] >= x_limit[2]) {
      stop('x_limit should be a positive scalar or an increasing length-2 vector.')
    }
    plot_df$effect <- pmin(pmax(plot_df$effect, x_limit[1]), x_limit[2])
  }

  group_level <- group_level[group_level %in% unique(plot_df$enriched)]
  extra_level <- setdiff(unique(plot_df$enriched), group_level)
  group_level <- c(group_level, extra_level)
  plot_df$enriched <- factor(plot_df$enriched, levels = group_level)
  if (is.null(group_color)) {
    default_color <- c(Up = '#f46d43', None = 'grey75', Down = '#66bd63')
    group_color <- default_color[group_level]
    missing_color <- is.na(group_color)
    group_color[missing_color] <- scales::hue_pal()(sum(missing_color))
  } else {
    group_color <- .resolve_group_colors(group_level, group_color)
  }
  group_label <- stats::setNames(
    paste0(
      group_level, ' (n=',
      vapply(group_level, \(x) sum(plot_df$enriched == x), integer(1)),
      ')'
    ),
    group_level
  )

  ## 4. 绘图并返回 ggplot 对象
  subtitle <- sprintf(
    '%s >= %.3g; %s < %.3g',
    subtitle_keywords[1], effect_cutoff,
    subtitle_keywords[2], p_cutoff
  )
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = effect, y = -log10(plot_pvalue), color = enriched)
  ) +
    ggplot2::geom_point(size = 0.7) +
    ggplot2::scale_color_manual(values = group_color, labels = group_label) +
    ggplot2::scale_y_continuous(expand = c(0.01, 0.01)) +
    ggplot2::geom_vline(
      xintercept = c(-effect_cutoff, effect_cutoff),
      linetype = 2, linewidth = 0.4, color = 'black'
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(p_cutoff),
      linetype = 2, linewidth = 0.4, color = 'black'
    ) +
    ggplot2::labs(
      x = xlab, y = ylab, title = title,
      subtitle = subtitle, color = 'Enriched in'
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.ticks = ggplot2::element_line(linewidth = 0.4, color = 'black'),
      axis.text = ggplot2::element_text(size = 8, color = 'black'),
      axis.title = ggplot2::element_text(size = 8, color = 'black'),
      axis.line = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 10, color = 'black'),
      panel.background = ggplot2::element_rect(
        linewidth = 0.4, color = 'black'
      ),
      panel.grid = ggplot2::element_blank(),
      legend.position.inside = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background = ggplot2::element_blank(),
      aspect.ratio = aspect_ratio
    )
  p
}
