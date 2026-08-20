#### Jin-Xin Meng, 20240307, 20260820, v0.2.0 ####

#' Fit Sloan's neutral community model
#'
#' Chinese summary: 拟合 neutral community model，返回摘要、拟合数据和模型。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param name Display or identifier name for name.
#' @param conf_level Optional order for `conf_level`.
#' @param max_iter Maximum number of iterations allowed during model fitting.
#' @return A result object described in the Details section.
#' @export
calcu_NCM <- function(
    profile, name = NULL, conf_level = 0.95, max_iter = 500) {
  ## 1. 整理 sample × feature count matrix
  profile_mat <- t(as.matrix(.as_profile_df(profile, numeric = TRUE)))
  if (any(!is.finite(profile_mat)) || any(profile_mat < 0)) {
    stop('profile should contain finite non-negative values.')
  }
  sample_n <- nrow(profile_mat)
  mean_depth <- mean(rowSums(profile_mat))
  if (sample_n < 3L || mean_depth <= 0) {
    stop('NCM requires at least three non-empty samples.')
  }

  ## 2. 计算 mean relative abundance 和 occurrence frequency
  mean_count <- colMeans(profile_mat)
  frequency <- colMeans(profile_mat > 0)
  keep_feature <- mean_count > 0 & frequency > 0
  fit_df <- data.frame(
    feature = colnames(profile_mat)[keep_feature],
    mean_abundance = mean_count[keep_feature] / mean_depth,
    frequency = frequency[keep_feature],
    check.names = FALSE
  )
  if (nrow(fit_df) < 3L) stop('NCM requires at least three observed features.')
  fit_df <- fit_df[order(fit_df$mean_abundance), , drop = FALSE]
  d_value <- 1 / mean_depth

  ## 3. 用 nonlinear least squares 估计 migration parameter m
  model_obj <- minpack.lm::nlsLM(
    frequency ~ stats::pbeta(
      d_value,
      mean_depth * m * mean_abundance,
      mean_depth * m * (1 - mean_abundance),
      lower.tail = FALSE
    ),
    data = fit_df, start = list(m = 0.1),
    lower = c(m = 1e-6), upper = c(m = 1e6),
    control = minpack.lm::nls.lm.control(maxiter = max_iter)
  )
  m_value <- as.numeric(stats::coef(model_obj)[['m']])
  fit_df$predicted <- stats::pbeta(
    d_value,
    mean_depth * m_value * fit_df$mean_abundance,
    mean_depth * m_value * (1 - fit_df$mean_abundance),
    lower.tail = FALSE
  )

  ## 4. 计算 Wilson confidence interval 和模型拟合指标
  alpha <- 1 - conf_level
  ci_df <- Hmisc::binconf(
    fit_df$predicted * sample_n, sample_n,
    alpha = alpha, method = 'wilson', return.df = TRUE
  )
  fit_df$lower <- ci_df$Lower
  fit_df$upper <- ci_df$Upper
  fit_df$class <- ifelse(
    fit_df$frequency < fit_df$lower, 'below',
    ifelse(fit_df$frequency > fit_df$upper, 'above', 'neutral')
  )
  r2 <- 1 - sum((fit_df$frequency - fit_df$predicted)^2) /
    sum((fit_df$frequency - mean(fit_df$frequency))^2)
  summary_df <- data.frame(
    r2 = r2, m = m_value, mN = m_value * mean_depth,
    sample_n = sample_n, feature_n = nrow(fit_df), model = 'NCM',
    check.names = FALSE
  )
  if (!is.null(name)) summary_df$name <- name

  list(summary_df = summary_df, fit_df = fit_df, model = model_obj)
}

#' Plot Sloan's neutral community model
#'
#' Chinese summary: 绘制 neutral community model 的拟合曲线和置信边界。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param ncm_result Optional result returned by `calcu_NCM()` to avoid refitting the model.
#' @param title Optional plot or result title.
#' @param class_color Color specification for `class_color`.
#' @return A plot object; analysis data or models may also be stored as attributes.
#' @export
plot_NCM <- function(
    profile = NULL, ncm_result = NULL, title = NULL,
    class_color = c(
      neutral = '#000000', below = '#a52a2a', above = '#29a6a6'
    )) {
  ## 1. 复用已有计算结果；未提供时才重新拟合
  if (is.null(ncm_result)) {
    if (is.null(profile)) stop('Supply either profile or ncm_result.')
    ncm_result <- calcu_NCM(profile)
  }
  if (!all(c('summary_df', 'fit_df', 'model') %in% names(ncm_result))) {
    stop('ncm_result should be returned by calcu_NCM().')
  }
  plot_df <- ncm_result$fit_df
  summary_df <- ncm_result$summary_df
  if (is.null(title)) title <- 'Neutral community model analysis'
  subtitle <- sprintf(
    'R2 = %.4f, m = %.4f, mN = %.4f',
    summary_df$r2[1], summary_df$m[1], summary_df$mN[1]
  )

  ## 2. 绘制观测频率、预测曲线和 confidence interval
  ggplot2::ggplot(plot_df) +
    ggplot2::geom_point(
      ggplot2::aes(log10(mean_abundance), frequency, color = class),
      size = 1.2, show.legend = FALSE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(log10(mean_abundance), predicted),
      linewidth = 0.6, color = 'blue'
    ) +
    ggplot2::geom_line(
      ggplot2::aes(log10(mean_abundance), upper),
      linewidth = 0.6, color = 'grey60', linetype = 'dashed'
    ) +
    ggplot2::geom_line(
      ggplot2::aes(log10(mean_abundance), lower),
      linewidth = 0.6, color = 'grey60', linetype = 'dashed'
    ) +
    ggplot2::scale_color_manual(values = class_color) +
    ggplot2::labs(
      x = 'Mean relative abundance (log10)',
      y = 'Occurrence frequency', title = title, subtitle = subtitle
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(linewidth = 0.4, color = 'black'),
      axis.text = ggplot2::element_text(size = 8, color = 'black'),
      axis.title = ggplot2::element_text(size = 8, color = 'black'),
      plot.title = ggplot2::element_text(size = 10, color = 'black'),
      plot.subtitle = ggplot2::element_text(size = 8, color = 'black'),
      panel.grid = ggplot2::element_blank(), aspect.ratio = 3 / 4
    )
}
