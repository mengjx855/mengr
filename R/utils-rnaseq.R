#### Jinxin Meng, 20251229, 20260527, v0.2 ####

# 20251229: create script, add function 'plot_gsea_barcode()'.
# 20260527: add function 'plot_GO_bar()', 'plot_GO_circular_bar()'.


#### plot_gsea_barcode ####
#' Plot Gsea Barcode utility
#'
#' `plot_gsea_barcode()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 绘制 GSEA barcode/running-score 风格结果图。
#'
#' @param gseaResult A GSEA result object containing ranked genes and enrichment results.
#' @param set_ID Gene-set identifier selected from a GSEA result.
#' @param bar_color Color specification for `bar_color`.
#' @param bar_width Numeric setting for `bar_width`.
#' @param core_enrichment Gene identifiers in the leading-edge or core-enrichment subset.
#' @param add_table Logical control for `add_table`.
#' @return A plot object; analysis data or models may also be stored as attributes.
#' @export
plot_gsea_barcode <- function(
  gseaResult, set_ID, bar_color = c(up = "#de77ae", down = "#35978f"),
  bar_width = 0.8, core_enrichment = FALSE, add_table = TRUE
) {
  if (!inherits(gseaResult, "gseaResult")) {
    stop("Need a gseaResult object")
  }

  if (!all(set_ID %in% gseaResult@result$ID)) {
    stop("IDs not all in gseaResult object")
  }

  if (length(bar_color) != 2) {
    bar_color <- c(up = "#de77ae", down = "#35978f")
  }

  if (is.null(names(bar_color))) {
    bar_color <- structure(bar_color, names = c("up", "down"))
  }

  if (isFALSE(core_enrichment)) {
    plot_bar <- purrr::map_dfr(
      set_ID, ~
        data.frame(
          TERM = .x,
          GENE = gseaResult@geneSets[[.x]]
        )
    )
  } else {
    plot_bar <- dplyr::filter(gseaResult@result, ID %in% set_ID) |>
      dplyr::group_by(ID) |>
      dplyr::group_modify(~ data.frame(GENE = stringr::str_split_1(.x$core_enrichment, "/"))) |>
      dplyr::ungroup() |>
      dplyr::rename(TERM = ID)
  }

  plot_enriched <- dplyr::filter(gseaResult@result, ID %in% set_ID) |>
    dplyr::select(ID, NES) |>
    dplyr::mutate(enriched = ifelse(NES > 0, "up", "down"))

  rank_data <- data.frame(
    GENE = names(gseaResult@geneList),
    RANK = seq_along(gseaResult@geneList)
  ) |>
    dplyr::left_join(plot_bar, by = "GENE") |>
    dplyr::filter(!is.na(TERM)) |>
    tibble::add_column(value = 1) |>
    dplyr::left_join(dplyr::select(plot_enriched, TERM = ID, enriched), by = "TERM")

  p <- ggplot2::ggplot(
    rank_data, ggplot2::aes(x = RANK, y = value)
  ) +
    ggplot2::geom_bar(
      ggplot2::aes(fill = enriched),
      color = NA, stat = "identity",
      width = bar_width, position = ggplot2::position_identity()
    ) +
    ggplot2::facet_grid(rows = ggplot2::vars(TERM)) +
    ggplot2::scale_fill_manual(values = bar_color) +
    ggplot2::scale_x_continuous(
      limits = c(1, length(gseaResult@geneList)), expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(
      x = "Rank in Ordered Dataset", y = "",
      fill = "Enriched in"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.ticks.length = grid::unit(1.8, "mm"),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = "black", size = 10),
      axis.text.y = ggplot2::element_blank(),
      legend.position = "none",
      panel.spacing.y = grid::unit(0, "cm"),
      strip.text.y = ggplot2::element_text(angle = 0, size = 12, hjust = 0, face = "italic"),
      strip.background.y = ggplot2::element_blank(),
      aspect.ratio = 1 / 10
    )

  if (isTRUE(add_table)) {
    label_df <- dplyr::filter(gseaResult@result, ID %in% set_ID) |>
      dplyr::select(ID, Description, NES, FDR = p.adjust) |>
      dplyr::mutate(
        NES = signif(NES, 4),
        FDR = signif(FDR, 4)
      ) |>
      data.frame(row.names = NULL)

    p_add <- ggpubr::ggtexttable(
      label_df,
      theme = ggpubr::ttheme("blank"), rows = NULL
    ) |>
      ggpubr::tab_add_hline(
        at.row = c(1, 2), row.side = "top",
        linewidth = 3, linetype = 1
      ) |>
      ggpubr::tab_add_hline(
        at.row = nrow(label_df) + 1, row.side = "bottom",
        linewidth = 3, linetype = 1
      )

    p <- cowplot::plot_grid(p, p_add, ncol = 1, align = "v")
  }
  return(p)
}

#### plot_GO_bar ####
# GO 富集分析普通横向柱状图
# data: GO enrichment 结果表
# term_col: GO term 名称列，默认 Description
# group_col: GO 分类列，默认 ONTOLOGY
# value_col: 柱子高度列，默认 FoldEnrichment
# sort_col: top_n 筛选前的排序列；NULL 表示按输入顺序
# top_n: 每个 group 保留前 top_n 个 term；NULL 表示不筛选
# sort_decreasing: sort_col 是否降序排序
# group_level: 指定 group 顺序
# group_name: 是否将 BP/CC/MF 改成完整名称
# facet: 是否按 group 分面

#' Plot GO Bar utility
#'
#' `plot_GO_bar()` provides a reusable mengR workflow with input validation and standardized
#'   output.
#'
#' Chinese summary: 绘制 GO enrichment 柱状图，并支持分组、排序和标签。
#'
#' @param data An input data frame or compatible object.
#' @param term_col Name of the `term_col` input column.
#' @param group_col Name of the grouping column.
#' @param value_col Name of the `value_col` input column.
#' @param sort_col Name of the `sort_col` input column.
#' @param top_n Number of highest-ranking features or categories retained.
#' @param sort_decreasing Whether sorting is performed in decreasing order.
#' @param group_level Optional order of group levels.
#' @param group_name Display or identifier name for group name.
#' @param palette Color palette name or vector supplied to the plot.
#' @param title Optional plot or result title.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param bar_width Numeric setting for `bar_width`.
#' @param bar_alpha Numeric setting for `bar_alpha`.
#' @param show_legend Logical control for `show_legend`.
#' @param facet Whether enrichment results are separated into facets.
#' @return A plot object; analysis data or models may also be stored as attributes.
#' @export
plot_GO_bar <- function(
  data,
  term_col = "Description", group_col = "ONTOLOGY",
  value_col = "FoldEnrichment", sort_col = NULL,
  top_n = 10, sort_decreasing = FALSE,
  group_level = NULL,
  group_name = c(
    BP = "Biological Process",
    CC = "Cellular Component",
    MF = "Molecular Function"
  ),
  palette = c("#66c2a5", "#fc8d62", "#8da0cb"),
  title = "GO Enrichment Analysis",
  xlab = "FoldEnrichment", ylab = "",
  bar_width = .75, bar_alpha = .85,
  show_legend = FALSE,
  facet = TRUE
) {
  data <- data.frame(data, check.names = FALSE)

  need_col <- c(term_col, group_col, value_col)

  if (!is.null(sort_col)) {
    need_col <- c(need_col, sort_col)
  }

  if (!all(need_col %in% colnames(data))) {
    stop(
      "data should contain columns: ",
      paste(need_col, collapse = " | ")
    )
  }

  ## 1. 整理数据
  data$.input_order <- seq_len(nrow(data))

  if (is.null(sort_col)) {
    data$.sort <- data$.input_order
  } else {
    data$.sort <- data[[sort_col]]
  }

  data <- data |>
    dplyr::transmute(
      term = .data[[term_col]],
      group = .data[[group_col]],
      value = .data[[value_col]],
      .sort = .data[[".sort"]]
    ) |>
    dplyr::filter(!is.na(term), !is.na(group), !is.na(value)) |>
    dplyr::mutate(
      term = as.character(term),
      group = as.character(group),
      value = suppressWarnings(as.numeric(value)),
      .sort = suppressWarnings(as.numeric(.sort))
    ) |>
    dplyr::filter(is.finite(value))

  if (nrow(data) == 0) {
    stop("No valid enrichment terms remained for plotting.")
  }

  ## 2. GO ontology 名称转换
  if (!is.null(group_name)) {
    data <- data |>
      dplyr::mutate(
        group = dplyr::case_when(
          group %in% names(group_name) ~ unname(group_name[group]),
          TRUE ~ group
        )
      )

    if (!is.null(group_level)) {
      group_level <- ifelse(
        group_level %in% names(group_name),
        unname(group_name[group_level]),
        group_level
      )
    }
  }

  if (is.null(group_level)) {
    group_level <- unique(data$group)
  }

  group_level <- group_level[group_level %in% unique(as.character(data$group))]

  if (length(group_level) == 0) {
    stop("No valid group remained for plotting.")
  }

  data <- data |>
    dplyr::mutate(group = factor(group, levels = group_level)) |>
    dplyr::filter(!is.na(group))

  ## 3. 每个 ontology 取前 top_n 个 term
  if (!is.null(top_n)) {
    if (isTRUE(sort_decreasing)) {
      data <- data |>
        dplyr::arrange(group, dplyr::desc(.sort))
    } else {
      data <- data |>
        dplyr::arrange(group, .sort)
    }

    data <- data |>
      dplyr::group_by(group) |>
      dplyr::slice_head(n = top_n) |>
      dplyr::ungroup()
  }

  ## 4. 按 value 排序，用唯一 ID 避免不同 group 中 term 重名
  data <- data |>
    dplyr::arrange(group, value) |>
    dplyr::mutate(
      .term_id = paste0(dplyr::row_number(), "__", term),
      .term_id = factor(.term_id, levels = .term_id)
    )

  term_label <- stats::setNames(data$term, data$.term_id)

  ## 5. 颜色
  palette <- rep(palette, length.out = length(group_level))
  names(palette) <- group_level

  ## 6. 作图
  p <- ggplot2::ggplot(data, ggplot2::aes(x = value, y = .term_id, fill = group)) +
    ggplot2::geom_col(width = bar_width, alpha = bar_alpha, color = NA) +
    ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
    ggplot2::scale_y_discrete(labels = term_label) +
    ggplot2::labs(x = xlab, y = ylab, title = title, fill = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 10, color = "black"),
      axis.title = ggplot2::element_text(size = 12, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = .5, color = "black"),
      axis.ticks.length = grid::unit(2, "mm"),
      axis.line = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(linewidth = .5, color = "black", fill = NA),
      panel.grid.major = ggplot2::element_line(color = "grey88", linewidth = .5),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = .5, size = 13, face = "bold"),
      strip.text = ggplot2::element_text(face = "bold", size = 11),
      strip.background = ggplot2::element_blank(),
      legend.position = ifelse(show_legend, "right", "none")
    )

  if (isTRUE(facet)) {
    p <- p +
      ggplot2::facet_grid(
        rows = ggplot2::vars(group),
        scales = "free_y",
        space = "free_y"
      ) +
      ggplot2::theme(
        panel.spacing.y = grid::unit(1, "mm")
      )
  }

  return(p)
}

#### plot_GO_circular_bar ####
# GO 富集分析环状柱状图
# data: GO enrichment 结果表
# term_col: GO term 名称列，默认 Description
# group_col: GO 分类列，默认 ONTOLOGY
# value_col: 柱子高度列，默认 FoldEnrichment
# sort_col: top_n 筛选前的排序列；NULL 表示按输入顺序
# top_n: 每个 group 保留前 top_n 个 term；NULL 表示不筛选
# sort_decreasing: sort_col 是否降序排序
# group_level: 指定 group 顺序
# group_name: 是否将 BP/CC/MF 改成完整名称
# empty_bar: 每个 group 后面添加的空白 bar 数量
# grid_breaks: 环形辅助线刻度；NULL 自动生成
# grid_n: 自动生成 grid_breaks 时的参考数量
# show_grid: 是否显示 grid line
# inner_size: 中间空白区域大小，数值越大中间越空
# label_pad: term 标签距离柱子顶端的距离
# label_wrap: term 标签换行宽度；NULL 表示不换行

#' Plot GO Circular Bar utility
#'
#' `plot_GO_circular_bar()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 将 GO enrichment 结果绘制为 circular barplot。
#'
#' @param data An input data frame or compatible object.
#' @param term_col Name of the `term_col` input column.
#' @param group_col Name of the grouping column.
#' @param value_col Name of the `value_col` input column.
#' @param sort_col Name of the `sort_col` input column.
#' @param top_n Number of highest-ranking features or categories retained.
#' @param sort_decreasing Whether sorting is performed in decreasing order.
#' @param group_level Optional order of group levels.
#' @param group_name Display or identifier name for group name.
#' @param empty_bar Angular width reserved as an empty separator in the circular bar plot.
#' @param palette Color palette name or vector supplied to the plot.
#' @param grid_breaks Reference values used to draw circular or radar grid lines.
#' @param grid_n Number of grid intervals.
#' @param show_grid Logical control for `show_grid`.
#' @param inner_size Numeric setting for `inner_size`.
#' @param label_pad Radial padding between circular bars and their labels.
#' @param label_wrap Maximum label width before circular labels are wrapped.
#' @param bar_width Numeric setting for `bar_width`.
#' @param bar_alpha Numeric setting for `bar_alpha`.
#' @param label_size Numeric setting for `label_size`.
#' @param group_label_size Numeric setting for `group_label_size`.
#' @param title Optional plot or result title.
#' @param show_legend Logical control for `show_legend`.
#' @return A plot object; analysis data or models may also be stored as attributes.
#' @export
plot_GO_circular_bar <- function(
  data,
  term_col = "Description", group_col = "ONTOLOGY",
  value_col = "FoldEnrichment", sort_col = NULL,
  top_n = 10, sort_decreasing = FALSE,
  group_level = NULL,
  group_name = c(
    BP = "Biological Process",
    CC = "Cellular Component",
    MF = "Molecular Function"
  ),
  empty_bar = 2, palette = c("#66c2a5", "#fc8d62", "#8da0cb"),
  grid_breaks = NULL, grid_n = 4, show_grid = TRUE,
  inner_size = 1, label_pad = NULL, label_wrap = NULL,
  bar_width = .75, bar_alpha = .65,
  label_size = 3, group_label_size = 3,
  title = "GO Enrichment Analysis",
  show_legend = FALSE
) {
  data <- data.frame(data, check.names = FALSE)

  need_col <- c(term_col, group_col, value_col)

  if (!is.null(sort_col)) {
    need_col <- c(need_col, sort_col)
  }

  if (!all(need_col %in% colnames(data))) {
    stop(
      "data should contain columns: ",
      paste(need_col, collapse = " | ")
    )
  }

  ## 1. 整理数据
  data$.input_order <- seq_len(nrow(data))

  if (is.null(sort_col)) {
    data$.sort <- data$.input_order
  } else {
    data$.sort <- data[[sort_col]]
  }

  data <- data |>
    dplyr::transmute(
      term = .data[[term_col]],
      group = .data[[group_col]],
      value = .data[[value_col]],
      .sort = .data[[".sort"]]
    ) |>
    dplyr::filter(!is.na(term), !is.na(group), !is.na(value)) |>
    dplyr::mutate(
      term = as.character(term),
      group = as.character(group),
      value = suppressWarnings(as.numeric(value)),
      .sort = suppressWarnings(as.numeric(.sort))
    ) |>
    dplyr::filter(is.finite(value))

  if (nrow(data) == 0) {
    stop("No valid enrichment terms remained for plotting.")
  }

  ## 2. GO ontology 名称转换
  if (!is.null(group_name)) {
    data <- data |>
      dplyr::mutate(
        group = dplyr::case_when(
          group %in% names(group_name) ~ unname(group_name[group]),
          TRUE ~ group
        )
      )

    ## 如果 group_level 输入的是 BP/CC/MF，也同步转换
    if (!is.null(group_level)) {
      group_level <- ifelse(
        group_level %in% names(group_name),
        unname(group_name[group_level]),
        group_level
      )
    }
  }

  if (is.null(group_level)) {
    group_level <- unique(data$group)
  }

  group_level <- group_level[group_level %in% unique(as.character(data$group))]

  if (length(group_level) == 0) {
    stop("No valid group remained for plotting.")
  }

  data <- data |>
    dplyr::mutate(group = factor(group, levels = group_level)) |>
    dplyr::filter(!is.na(group))

  ## 3. 每个 ontology 取前 top_n 个 term
  if (!is.null(top_n)) {
    if (isTRUE(sort_decreasing)) {
      data <- data |>
        dplyr::arrange(group, dplyr::desc(.sort))
    } else {
      data <- data |>
        dplyr::arrange(group, .sort)
    }

    data <- data |>
      dplyr::group_by(group) |>
      dplyr::slice_head(n = top_n) |>
      dplyr::ungroup()
  }

  ## 4. 组内按 value 排序，并记录顺序
  data <- data |>
    dplyr::arrange(group, value) |>
    dplyr::group_by(group) |>
    dplyr::mutate(
      .order = dplyr::row_number(),
      .empty = FALSE
    ) |>
    dplyr::ungroup()

  ## 5. 添加 empty bar
  if (empty_bar > 0) {
    to_add <- data.frame(
      term = NA_character_,
      group = rep(group_level, each = empty_bar),
      value = NA_real_,
      .sort = NA_real_,
      .order = rep(seq_len(empty_bar), times = length(group_level)) + 1e6,
      .empty = TRUE,
      check.names = FALSE
    )

    to_add$group <- factor(to_add$group, levels = group_level)

    data <- dplyr::bind_rows(data, to_add)
  }

  data <- data |>
    dplyr::arrange(group, .empty, .order) |>
    dplyr::mutate(id = dplyr::row_number())

  ## 6. label 角度
  label_df <- data |>
    dplyr::filter(!.empty, !is.na(term), is.finite(value))

  number_of_bar <- nrow(data)

  label_df <- label_df |>
    dplyr::mutate(
      angle = 90 - 360 * (id - 0.5) / number_of_bar,
      hjust = ifelse(angle < -90, 1, 0),
      angle = ifelse(angle < -90, angle + 180, angle)
    )

  if (!is.null(label_wrap)) {
    label_df$term <- stringr::str_wrap(label_df$term, width = label_wrap)
  }

  ## 7. group baseline 数据
  base_data <- data |>
    dplyr::filter(!.empty) |>
    dplyr::group_by(group) |>
    dplyr::summarise(
      start = min(id),
      end = max(id),
      .groups = "drop"
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(title = mean(c(start, end))) |>
    dplyr::ungroup()

  ## 8. 坐标范围
  max_value <- max(data$value, na.rm = TRUE)

  if (!is.finite(max_value) || max_value <= 0) {
    stop("value_col should contain positive numeric values.")
  }

  if (is.null(label_pad)) {
    label_pad <- max_value * .05
  }

  if (is.null(grid_breaks)) {
    grid_breaks <- pretty(c(0, max_value), n = grid_n)
    grid_breaks <- grid_breaks[grid_breaks > 0 & grid_breaks < max_value]
  }

  y_min <- -max_value * inner_size
  y_max <- max_value + label_pad * 4

  ## 9. grid line 数据
  ## 只在每个 group 后面的 empty_bar 区域画 grid
  ## 最后一个 group 后面的 empty_bar 也保留

  grid_data <- base_data |>
    dplyr::mutate(
      gap_start = end + 1,
      gap_end = dplyr::lead(start) - 1
    )

  grid_data$gap_end[nrow(grid_data)] <- nrow(data)

  grid_data <- grid_data |>
    dplyr::filter(gap_end >= gap_start)

  ## 10. 颜色
  palette <- rep(palette, length.out = length(group_level))
  names(palette) <- group_level

  ## 11. 作图
  p <- ggplot2::ggplot(data, ggplot2::aes(x = id, y = value, fill = group))

  ## grid lines 放在前面，让柱子盖在上面
  ## 这里只在 empty_bar 空白区域画 grid

  if (isTRUE(show_grid) && nrow(grid_data) > 0 && length(grid_breaks) > 0) {
    for (g in grid_breaks) {
      p <- p +
        ggplot2::geom_segment(
          data = grid_data, ggplot2::aes(x = gap_start, xend = gap_end),
          y = g, yend = g, colour = "grey80", linewidth = .3,
          linetype = "longdash", inherit.aes = FALSE
        )
    }

    p <- p +
      ggplot2::annotate(
        "text",
        x = max(data$id), y = grid_breaks, label = grid_breaks,
        color = "grey50", size = 3, fontface = "bold", hjust = 1
      )
  }

  p <- p +
    ggplot2::geom_col(width = bar_width, alpha = bar_alpha, color = NA, na.rm = TRUE) +
    ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
    ggplot2::scale_y_continuous(limits = c(y_min, y_max), expand = c(0, 0)) +
    ggplot2::coord_polar() +
    ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(
        x = id, y = value + label_pad,
        label = term, hjust = hjust, angle = angle
      ),
      color = "black", size = label_size, alpha = .75,
      fontface = "plain", inherit.aes = FALSE,
      na.rm = TRUE
    ) +
    ggplot2::geom_segment(
      data = base_data,
      ggplot2::aes(x = start, y = y_min * .08, xend = end, yend = y_min * .08),
      colour = "black", alpha = .8,
      linewidth = .6, inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = base_data,
      ggplot2::aes(x = title, y = y_min * .22, label = group),
      colour = "black", alpha = .8, size = group_label_size,
      fontface = "bold", inherit.aes = FALSE
    ) +
    ggplot2::labs(title = title, fill = "") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = ifelse(show_legend, "right", "none"),
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = .5, size = 13, face = "bold"),
      plot.margin = grid::unit(c(5, 5, 5, 5), "mm")
    )

  return(p)
}
