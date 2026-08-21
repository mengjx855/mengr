#### Jinxin Meng, 20250308, 20260502 v0.2.1 ####

# 20250327: update some parameter.
# 20250327: add mean abundance of feature as output in calcu_stamp() with option 'method=wilcox'.
# 20250517: rename file 'calcu_stamp' as 'plot_stamp'.
# 20260301: update some parameter. rename padj as qval.
# 20260301: add new function with 'multiple' suffix for multivariable analysis.
# 20260502: update functions.


#### help function ####
.check_group <- function(
  group, sample_col = "sample", group_col = "group",
  comparison = NULL, two_group = FALSE
) {
  ## 统一为内部使用的 sample/group 两列
  group_df <- .as_df(group)
  .check_columns(group_df, c(sample_col, group_col), object = "group")
  group_df <- data.frame(
    sample = as.character(group_df[[sample_col]]),
    group = group_df[[group_col]], check.names = FALSE
  )
  if (anyDuplicated(group_df$sample)) {
    stop("Duplicated sample identifiers found in group.")
  }
  if (is.null(comparison)) comparison <- unique(as.character(group_df$group))
  group_df <- group_df[group_df$group %in% comparison, , drop = FALSE]
  if (isTRUE(two_group) && length(comparison) != 2L) {
    stop("This function only supports two-group comparison.")
  }
  group_df$group <- factor(group_df$group, levels = comparison)
  list(group_df = group_df, comparison = comparison)
}


.add_p_q_label <- function(data) {
  
  data |>
    dplyr::rename(pval = p) |>
    rstatix::add_significance(
      p.col = "pval",
      output.col = "plab",
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols = c("***", "**", "*", "ns")
    ) |>
    dplyr::mutate(qval = stats::p.adjust(pval, method = "BH")) |>
    rstatix::add_significance(
      p.col = "qval",
      output.col = "qlab",
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols = c("***", "**", "*", "ns")
    ) |>
    dplyr::relocate(plab, qval, qlab, .after = pval) |>
    dplyr::select(-dplyr::any_of(".y.")) |>
    dplyr::rename_with(~ gsub("\\.", "_", .x))
}

#### calcu_stamp ####
#' Calcu Stamp utility
#'
#' `calcu_stamp()` provides a reusable mengR workflow with input validation and standardized
#'   output.
#'
#' Chinese summary: 计算两组 STAMP 风格均值差、置信区间和显著性结果。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param group A sample metadata table containing sample and group columns.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param comparison Two outcome levels ordered as case and control.
#' @param method Analysis or summary method; supported values are shown in the usage.
#' @param var_equal Whether two-sample t tests assume equal group variances.
#' @param exact Whether an exact test is requested when supported by the selected method.
#' @param qvalue Maximum adjusted P value retained in the result.
#' @param pvalue Maximum raw P value retained in the result.
#' @param add_enriched Logical control for `add_enriched`.
#' @return A result object described in the Details section.
#' @export
calcu_stamp <- function(
  profile, group, sample_col = "sample", group_col = "group",
  comparison = NULL, method = c("wilcox", "t"),
  var_equal = FALSE, exact = NULL, qvalue = 0.2,
  pvalue = NULL, add_enriched = TRUE
) {
  ## 1. 检查方法、分组并对齐样本
  method <- match.arg(method)
  group_info <- .check_group(
    group = group, sample_col = sample_col, group_col = group_col,
    comparison = comparison, two_group = TRUE
  )
  comparison <- group_info$comparison
  aligned <- .align_profile_group(
    profile = profile, group = group_info$group_df,
    sample_col = "sample", group_col = "group",
    group_level = comparison
  )
  profile_df <- aligned$profile_df
  group_df <- aligned$group_df

  ## 2. 转为 long format，并计算各组 mean abundance
  test_df <- profile_df |>
    tibble::rownames_to_column("name") |>
    tidyr::pivot_longer(
      cols = -name,
      names_to = "sample",
      values_to = "value"
    ) |>
    dplyr::left_join(group_df, by = "sample") |>
    dplyr::group_by(name)

  mean_df <- test_df |>
    dplyr::ungroup() |>
    dplyr::group_by(name, group) |>
    dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(
      names_from = group,
      values_from = value
    )

  if (method == "t") {
    result_df <- rstatix::t_test(
      test_df, value ~ group,
      var.equal = var_equal, detailed = TRUE
    ) |>
      .add_p_q_label()

    if (all(c("estimate1", "estimate2") %in% colnames(result_df))) {
      colnames(result_df)[
        match(c("estimate1", "estimate2"), colnames(result_df))
      ] <- comparison
    }

    result_df$method <- ifelse(
      isTRUE(var_equal), "Two Sample t-test", "Welch Two Sample t-test"
    )
  }

  if (method == "wilcox") {
    result_df <- rstatix::wilcox_test(
      test_df, value ~ group,
      exact = exact, detailed = TRUE
    ) |>
      .add_p_q_label() |>
      dplyr::left_join(mean_df, by = "name") |>
      dplyr::relocate(dplyr::all_of(comparison), .after = "estimate")
  }

  ## 3. 根据 pvalue 或 qvalue 判断富集方向
  if (isTRUE(add_enriched)) {
    if (is.numeric(pvalue)) {
      result_df <- result_df |>
        dplyr::mutate(
          enriched = dplyr::case_when(
            estimate > 0 & pval < pvalue ~ comparison[1],
            estimate < 0 & pval < pvalue ~ comparison[2],
            TRUE ~ "none"
          )
        )
    } else {
      result_df <- result_df |>
        dplyr::mutate(
          enriched = dplyr::case_when(
            estimate > 0 & qval < qvalue ~ comparison[1],
            estimate < 0 & qval < qvalue ~ comparison[2],
            TRUE ~ "none"
          )
        )
    }
  }
  result_df
}

#### plot_stamp ####
#' Plot Stamp utility
#'
#' `plot_stamp()` provides a reusable mengR workflow with input validation and standardized
#'   output.
#'
#' Chinese summary: 将 `calcu_stamp()` 结果绘制为 STAMP 风格多面板图。
#'
#' @param data An input data frame or compatible object.
#' @param top_n Number of highest-ranking features or categories retained.
#' @param comparison Two outcome levels ordered as case and control.
#' @param palette Color palette name or vector supplied to the plot.
#' @param left_xlab Optional label used for left xlab.
#' @param left_title Optional label used for left title.
#' @param mid_xlab Optional label used for mid xlab.
#' @param mid_title Optional label used for mid title.
#' @param right_xlab Optional label used for right xlab.
#' @param right_title Optional label used for right title.
#' @param show_label Logical control for `show_label`.
#' @return A plot object; analysis data or models may also be stored as attributes.
#' @export
plot_stamp <- function(
  data, top_n = 10, comparison = NULL, palette = c("#E69F00", "#56B4E9"),
  left_xlab = "Mean proportion (%)", left_title = "Relative abundance (%)",
  mid_xlab = "Difference in mean proportions (%)",
  mid_title = "95% confidence intervals",
  right_xlab = "", right_title = "FDR",
  show_label = "qval"
) {
  if (is.null(comparison)) {
    comparison <- as.character(unlist(data[1, c("group1", "group2")]))
  }

  if (!all(comparison %in% colnames(data))) {
    stop("comparison columns are not found in data.")
  }

  if (!show_label %in% colnames(data)) {
    stop("show_label is not found in data.")
  }

  if (!"enriched" %in% colnames(data)) {
    data$enriched <- "none"
  }

  palette <- structure(palette, names = comparison)
  palette <- c(palette, none = "grey80")

  plot_df <- .as_df(data) |>
    dplyr::mutate(abs_estimate = abs(estimate)) |>
    dplyr::arrange(dplyr::desc(abs_estimate)) |>
    utils::head(top_n) |>
    dplyr::arrange(dplyr::desc(estimate)) |>
    dplyr::mutate(
      name = factor(name, name),
      pval = signif(pval, 5),
      qval = signif(qval, 5),
      enriched = ifelse(is.na(enriched), "none", enriched)
    )

  abundance_df <- plot_df |>
    dplyr::select(name, dplyr::all_of(comparison)) |>
    tidyr::pivot_longer(
      cols = -name,
      names_to = "group",
      values_to = "value"
    )

  add_bg_stripe <- function(p, n) {
    if (n > 1) {
      for (i in seq_len(n - 1)) {
        p <- p + ggplot2::annotate(
          "rect",
          ymin = i + 0.5,
          ymax = i + 1.5,
          xmin = -Inf,
          xmax = Inf,
          fill = ifelse(i %% 2 == 0, "white", "gray95")
        )
      }
    }
    p
  }

  p1 <- ggplot2::ggplot(
    abundance_df, ggplot2::aes(value, name, fill = group)
  ) +
    ggplot2::scale_y_discrete(
      limits = levels(plot_df$name),
      labels = function(x) stringr::str_wrap(x, width = 35)
    ) +
    ggplot2::labs(y = "", x = left_xlab, title = left_title) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "transparent"),
      panel.grid = ggplot2::element_blank(),
      axis.ticks.length = grid::unit(0.4, "lines"),
      axis.ticks = ggplot2::element_line(color = "black"),
      axis.line = ggplot2::element_line(color = "black"),
      axis.title.x = ggplot2::element_text(color = "black", size = 12),
      axis.text = ggplot2::element_text(color = "black", size = 12),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(
        size = 10, color = "black", vjust = .5, hjust = .5,
        margin = ggplot2::margin(t = 0, unit = "cm")
      ),
      legend.position = "bottom",
      legend.background = ggplot2::element_rect(fill = "transparent"),
      legend.direction = "horizontal",
      legend.key.width = grid::unit(0.7, "cm"),
      legend.key.height = grid::unit(0.4, "cm"),
      plot.title = ggplot2::element_text(
        size = 12, color = "black", hjust = 0.5,
        margin = ggplot2::margin(0, 0, .3, 0, "cm")
      )
    )

  p1 <- add_bg_stripe(p1, nrow(plot_df))

  p1 <- p1 +
    ggplot2::geom_bar(
      stat = "identity",
      position = "dodge",
      width = 0.7,
      color = "black",
      show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(values = palette)

  p2 <- ggplot2::ggplot(
    plot_df, ggplot2::aes(estimate, name, fill = enriched)
  ) +
    ggplot2::scale_y_discrete(limits = levels(plot_df$name)) +
    ggplot2::labs(y = "", x = mid_xlab, title = mid_title) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "transparent"),
      panel.grid = ggplot2::element_blank(),
      axis.ticks.length = grid::unit(0.4, "lines"),
      axis.ticks = ggplot2::element_line(color = "black"),
      axis.line = ggplot2::element_line(color = "black"),
      axis.title.x = ggplot2::element_text(color = "black", size = 12),
      axis.text = ggplot2::element_text(color = "black", size = 12),
      axis.text.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = "none",
      plot.title = ggplot2::element_text(
        size = 12, color = "black", hjust = 0.5,
        margin = ggplot2::margin(0, 0, .3, 0, "cm")
      )
    )

  p2 <- add_bg_stripe(p2, nrow(plot_df))

  p2 <- p2 +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = conf_low, xmax = conf_high),
      position = ggplot2::position_dodge(0.8),
      width = 0.3,
      linewidth = 0.3
    ) +
    ggplot2::geom_point(shape = 21, size = 4) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "black",
      linewidth = .3
    ) +
    ggplot2::scale_fill_manual(values = palette)

  label_vec <- plot_df[[show_label]]

  if (is.numeric(label_vec)) {
    plot_df$show_label <- sprintf("%.4f", label_vec)
  } else {
    plot_df$show_label <- as.character(label_vec)
  }

  p3 <- ggplot2::ggplot(plot_df, ggplot2::aes(1, name)) +
    ggplot2::geom_text(
      ggplot2::aes(x = 0, y = name, label = show_label),
      color = "black",
      vjust = .5,
      fontface = "bold",
      inherit.aes = FALSE,
      size = 4
    ) +
    ggplot2::labs(x = right_xlab, title = right_title) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 12, color = "black", hjust = 0.5,
        margin = ggplot2::margin(0, 0, .3, 0, "cm")
      )
    )

  p <- p1 + p2 + p3 + patchwork::plot_layout(widths = c(5, 5, 1))

  return(p)
}

#### calcu_stamp_multiple ####
#' Calcu Stamp Multiple utility
#'
#' `calcu_stamp_multiple()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 对多组数据计算 STAMP 风格的整体差异统计。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param group A sample metadata table containing sample and group columns.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param group_level Optional order of group levels.
#' @param method Analysis or summary method; supported values are shown in the usage.
#' @return A result object described in the Details section.
#' @export
calcu_stamp_multiple <- function(
  profile, group, sample_col = "sample", group_col = "group",
  group_level = NULL, method = c("kruskal", "anova")
) {
  ## 1. 检查方法、分组并对齐样本
  method <- match.arg(method)
  group_info <- .check_group(
    group = group, sample_col = sample_col, group_col = group_col,
    comparison = group_level, two_group = FALSE
  )
  aligned <- .align_profile_group(
    profile = profile, group = group_info$group_df,
    sample_col = "sample", group_col = "group",
    group_level = group_info$comparison
  )
  profile_df <- aligned$profile_df
  group_df <- aligned$group_df

  ## 2. 转为 long format，并计算各组 mean abundance
  test_df <- profile_df |>
    tibble::rownames_to_column("name") |>
    tidyr::pivot_longer(
      cols = -name,
      names_to = "sample",
      values_to = "value"
    ) |>
    dplyr::left_join(group_df, by = "sample") |>
    dplyr::group_by(name)

  mean_df <- test_df |>
    dplyr::ungroup() |>
    dplyr::group_by(name, group) |>
    dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(
      names_from = group,
      values_from = value
    )

  ## 3. 进行多组检验并整理输出
  if (method == "kruskal") {
    test_result_df <- rstatix::kruskal_test(test_df, value ~ group) |>
      .add_p_q_label()
    result_df <- dplyr::left_join(mean_df, test_result_df, by = "name")
  }

  if (method == "anova") {
    test_result_df <- rstatix::anova_test(test_df, value ~ group) |>
      data.frame(check.names = FALSE) |>
      dplyr::select(-dplyr::any_of(c("p<.05", "Effect"))) |>
      .add_p_q_label()
    result_df <- dplyr::left_join(mean_df, test_result_df, by = "name") |>
      dplyr::mutate(method = "ANOVA-test")
  }

  result_df
}
