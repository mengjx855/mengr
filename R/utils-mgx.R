#### Jin-Xin Meng, 20250418, 20260619, v0.1.5 ####

# 20250423: update some function.
# 20250502: rename 'profile_KEGG_convert' as 'profile_KEGG_trans'
# 20260502: update some function.
# 20260504: add new function plot_maaslin3_multi()
# 20260619: add new function plot_maaslin3_abundance(), plot_maaslin3_both()


#### 20260504 plot_maaslin3_multi ####
#' Plot multi-level MaAsLin3 taxonomic effects
#'
#' This function visualizes MaAsLin3 association results across multiple
#' taxonomic ranks. It is designed for results in which features are taxonomic
#' strings such as:
#'
#'   k__Bacteria|p__Bacillota|c__Clostridia|o__Eubacteriales|...
#'
#' For a selected metadata variable and value, the function shows the MaAsLin3
#' coefficient at different taxonomic levels. Positive coefficients are plotted
#' in the upper panel and negative coefficients in the lower panel. The middle
#' strip indicates taxonomic ranks.
#'
#' Interpretation:
#'   - coef > 0: the taxon is positively associated with `value`
#'               relative to `ref`.
#'   - coef < 0: the taxon is negatively associated with `value`
#'               and is relatively enriched in the reference group.
#'
#' Example:
#'   metadata = "group", value = "PPM", ref = "CN"
#'   coef > 0 means the taxon is enriched/positively associated with PPM.
#'   coef < 0 means the taxon is enriched/positively associated with CN.
#'
#' Example:
#'   metadata = "gender", value = "M", ref = "F"
#'   coef > 0 means the taxon is enriched/positively associated with males.
#'   coef < 0 means the taxon is enriched/positively associated with females.
#'
#'
#' Chinese summary: 绘制 MaAsLin3 多分类层级的正负效应及 taxonomy 色带。
#'
#' @param maaslin3_result
#'   MaAsLin3 result table. Required columns usually include:
#'   feature, metadata, value, coef, qval_individual.
#'   If columns `model` and `error` exist, the function keeps
#'   model == "abundance" and removes rows with non-NA error.
#'
#' @param metadata
#'   Metadata variable to plot, for example "group" or "gender".
#'
#' @param value
#'   The level of metadata to plot, for example "PPM" or "M".
#'   This corresponds to the `value` column in MaAsLin3 output.
#'
#' @param ref
#'   Reference level used only for labeling and interpretation.
#'   It does not affect the data calculation. For example "CN" or "F".
#'
#' @param feature_col
#'   Column name containing taxonomic features.
#'
#' @param coef_col
#'   Column name containing MaAsLin3 coefficients.
#'
#' @param p_col
#'   Column name used for significance. Common choices:
#'   "qval_individual", "qval_joint", or "pval".
#'
#' @param p_cutoff
#'   Significance threshold used to color points.
#'   Points with p_col < p_cutoff are colored as significant.
#'
#' @param levels_keep
#'   Taxonomic ranks to show. Use rank letters:
#'   "p" = Phylum, "c" = Class, "o" = Order,
#'   "f" = Family, "g" = Genus, "s" = Species.
#'
#' @param point_top_n
#'   Number of taxa to plot per taxonomic rank.
#'   If point_side = "both", top N positive and top N negative taxa are
#'   selected separately within each rank.
#'   If NULL, all selected taxa are plotted.
#'
#' @param point_by
#'   Criterion for selecting plotted taxa:
#'   "p" = select by smallest p/q value;
#'   "abs_coef" = select by largest absolute coefficient.
#'
#' @param point_side
#'   Which coefficient direction to plot:
#'   "both" = select positive and negative taxa separately;
#'   "positive" = only coef > 0;
#'   "negative" = only coef < 0;
#'   "all" = ignore direction and select top taxa per rank.
#'
#' @param label_top_n
#'   Number of taxa to label per taxonomic rank among already plotted taxa.
#'   If label_side = "both", top N positive and top N negative labels are
#'   selected separately within each rank.
#'   If NULL, all plotted taxa are labeled.
#'
#' @param label_by
#'   Criterion for selecting labels:
#'   "p" = label taxa with smallest p/q value;
#'   "abs_coef" = label taxa with largest absolute coefficient.
#'
#' @param label_side
#'   Which coefficient direction to label:
#'   "both", "positive", "negative", or "all".
#'
#' @param manual_label_taxa
#'   Character vector of short taxon names to label manually, for example
#'   c("Lactobacillus", "Muribaculaceae").
#'   These labels are added only if the taxa are already included in
#'   the plotted points.
#'
#' @param manual_label_features
#'   Character vector of full feature strings to label manually.
#'   These labels are added only if the features are already included in
#'   the plotted points.
#'
#' @param point_size
#'   Size of plotted points.
#'
#' @param jitter_width
#'   Horizontal jitter width. The jitter coordinates are generated once and
#'   shared by points and labels, so labels match their corresponding points.
#'
#' @param strip_height
#'   Height of the middle taxonomic-rank strip. If NULL, it is determined
#'   automatically from the coefficient range.
#'
#' @param strip_palette
#'   Color palette for the middle taxonomic-rank strip. Colors are interpolated
#'   automatically according to the number of displayed ranks.
#'
#' @param sig_colors
#'   Named vector for point fill colors. Must contain names "sig" and "ns".
#'   Example: c(sig = "#D53E4F", ns = "#0099C7").
#'
#' @param title
#'   Plot title. If NULL, a default title is generated.
#'
#' @param verbose
#'   Whether to print summary information and recommended ggsave command.
#'
#' @return
#'   A patchwork ggplot object. The selected plotting data and label data are
#'   stored as attributes:
#'     attr(p, "plot_df")
#'     attr(p, "label_df")
#'
#' @examples
#' p <- plot_maaslin3_multi(
#'   maaslin3_result = maaslin_res,
#'   metadata = "group",
#'   value = "PPM",
#'   ref = "CN",
#'   p_col = "qval_individual",
#'   p_cutoff = 0.05,
#'   point_top_n = 100,
#'   label_top_n = 5,
#'   point_side = "both",
#'   label_side = "both",
#'   title = "PPM exposure-associated taxa"
#' )
#'
#' ggsave("PPM_associated_taxa.pdf", p, width = 9, height = 6)
#' @export
plot_maaslin3_multi <- function(
    maaslin3_result,
    metadata = "group", value = "Case", ref = NULL,
    feature_col = "feature", coef_col = "coef",
    p_col = "qval_individual", p_cutoff = 0.05,
    levels_keep = c("p", "c", "o", "f", "g", "s"),
    point_top_n = 20, point_by = c("p", "abs_coef"),
    point_side = c("both", "positive", "negative", "all"),
    label_top_n = 5, label_by = c("p", "abs_coef"),
    label_side = c("both", "positive", "negative", "all"),
    manual_label_taxa = NULL, manual_label_features = NULL,
    point_size = 2.5, jitter_width = 0.3, strip_height = NULL,
    strip_palette = c(
      "#F46D43","#FDAE61","#FEE08B","#FFFFBF",
      "#E6F598","#ABDDA4","#66C2A5","#3288BD"
    ),
    sig_colors = c(sig = "#D53E4F", ns = "#0099C7"),
    title = NULL, verbose = TRUE) {
  
  ## 0. 参数匹配与结果表检查
  ## point_by/label_by 控制按显著性或效应大小排序；
  ## point_side/label_side 控制选择正向、负向或双向结果。
  point_by <- match.arg(point_by)
  label_by <- match.arg(label_by)
  point_side <- match.arg(point_side)
  label_side <- match.arg(label_side)
  
  ## 1. 整理 MaAsLin3 结果
  ## 仅保留 abundance 模型结果，并去除报错行。
  ## 如果没有 model/error 列，则自动跳过。
  maaslin <- data.frame(maaslin3_result, check.names = FALSE)
  
  if ("model" %in% colnames(maaslin)) {
    maaslin <- dplyr::filter(maaslin, model == "abundance")
  }
  if ("error" %in% colnames(maaslin)) {
    maaslin <- dplyr::filter(maaslin, is.na(error))
  }
  
  need_cols <- c(feature_col, "metadata", "value", coef_col, p_col)
  miss_cols <- setdiff(need_cols, colnames(maaslin))
  if (length(miss_cols) > 0) {
    stop("Missing columns in MaAsLin3 result: ", paste(miss_cols, collapse = ", "))
  }
  
  ## 2. 解析 feature 中的分类层级
  ## 从形如 k__Bacteria|p__...|g__... 的字符串中提取最后一级分类。
  ## .rank_letter 用于判断层级，.taxon/.label 用于后续标注。
  rank_map <- c(
    d = "Domain", k = "Kingdom", p = "Phylum", c = "Class",
    o = "Order", f = "Family", g = "Genus", s = "Species", t = "Strain"
  )
  rank_order <- c(
    "Domain", "Kingdom", "Phylum", "Class",
    "Order", "Family", "Genus", "Species", "Strain"
  )
  
  parse_taxa <- function(x) {
    parts <- stringr::str_extract_all(as.character(x), "[dkpcofgst]__[^|;]+")[[1]]
    if (length(parts) == 0) {
      return(data.frame(rank_letter = NA, rank = NA, taxon = x, label = x))
    }
    last <- parts[length(parts)]
    rank_letter <- substr(last, 1, 1)
    taxon <- sub("^[dkpcofgst]__", "", last)
    data.frame(
      rank_letter = rank_letter,
      rank = unname(rank_map[rank_letter]),
      taxon = taxon,
      label = taxon
    )
  }
  
  taxa_info <- dplyr::bind_rows(lapply(maaslin[[feature_col]], parse_taxa))
  maaslin$.rank_letter <- taxa_info$rank_letter
  maaslin$.rank <- taxa_info$rank
  maaslin$.taxon <- taxa_info$taxon
  maaslin$.label <- taxa_info$label
  
  ## 3. 提取目标变量的效应结果
  ## 例如 metadata = "group", value = "PPM" 表示提取 PPM 相对于参考组的效应。
  ## coef > 0 表示与 value 正相关；coef < 0 表示与参考组方向相关。
  data <- maaslin |>
    dplyr::filter(
      .data[["metadata"]] == {{ metadata }},
      .data[["value"]] == {{ value }}
    ) |>
    dplyr::mutate(
      .coef = as.numeric(.data[[coef_col]]),
      .p = as.numeric(.data[[p_col]]),
      .abs_coef = abs(.coef),
      .direction = dplyr::case_when(
        .coef > 0 ~ "positive",
        .coef < 0 ~ "negative",
        TRUE ~ "zero"
      ),
      .sig_group = ifelse(.p < p_cutoff, "sig", "ns")
    ) |>
    dplyr::filter(!is.na(.coef), !is.na(.p), !is.na(.rank))
  
  ## 4. 根据分类层级筛选结果
  ## levels_keep 使用分类层级字母：
  ## p = Phylum, c = Class, o = Order, f = Family, g = Genus, s = Species。
  if (!is.null(levels_keep)) {
    data <- dplyr::filter(data, .rank_letter %in% levels_keep)
  }
  if (nrow(data) == 0) {
    stop("No rows remained after filtering metadata/value/taxonomic levels.")
  }
  
  data$.rank <- factor(data$.rank, levels = rank_order)
  
  ## 5. 定义 top 选择函数
  ## 选择每个分类层级中的 top 结果
  ## by = "p" 时优先选择 p/q 值最小的结果；
  ## by = "abs_coef" 时优先选择绝对效应值最大的结果。
  ## side = "both" 时，正向和负向分别筛选，避免只展示单一方向。
  select_top <- function(x, top_n, by, side) {
    if (side == "positive") x <- dplyr::filter(x, .coef > 0)
    if (side == "negative") x <- dplyr::filter(x, .coef < 0)
    
    if (by == "p") {
      x <- dplyr::arrange(x, .rank, .p, dplyr::desc(.abs_coef))
    } else {
      x <- dplyr::arrange(x, .rank, dplyr::desc(.abs_coef), .p)
    }
    
    if (is.null(top_n)) return(x)
    
    if (side == "both") {
      x |>
        dplyr::filter(.direction %in% c("positive", "negative")) |>
        dplyr::group_by(.rank, .direction) |>
        dplyr::slice_head(n = top_n) |>
        dplyr::ungroup()
    } else {
      x |>
        dplyr::group_by(.rank) |>
        dplyr::slice_head(n = top_n) |>
        dplyr::ungroup()
    }
  }
  
  ## 6. 选择要展示的点
  ## point_top_n 控制每个分类层级展示多少个点。
  ## 如果 point_side = "both"，则每个层级中正向和负向各取 point_top_n 个。
  plot_df <- select_top(data, point_top_n, point_by, point_side)
  plot_df <- plot_df[!duplicated(plot_df[[feature_col]]), , drop = FALSE]
  if (nrow(plot_df) == 0) stop("No taxa selected for plotting.")
  
  ## 7. 选择要标注的标签
  ## 标签只从已经展示的点中选择，不会额外增加新的点。
  ## label_top_n 控制每个分类层级标注多少个标签。
  label_df <- select_top(plot_df, label_top_n, label_by, label_side)
  
  ## 8. 加入手动指定标签
  ## manual_label_taxa 使用短分类名匹配，例如 Lactobacillus；
  ## manual_label_features 使用完整 feature 字符串匹配。
  ## 注意：手动标签也只会在已展示点中寻找。
  manual_label_df <- NULL
  if (!is.null(manual_label_taxa)) {
    manual_label_df <- plot_df |>
      dplyr::filter(.taxon %in% manual_label_taxa | .label %in% manual_label_taxa)
  }
  if (!is.null(manual_label_features)) {
    selected_df <- plot_df |>
      dplyr::filter(.data[[feature_col]] %in% manual_label_features)
    manual_label_df <- dplyr::bind_rows(manual_label_df, selected_df)
  }
  
  label_df <- dplyr::bind_rows(label_df, manual_label_df)
  label_df <- label_df[!duplicated(label_df[[feature_col]]), , drop = FALSE]
  
  if (!is.null(manual_label_taxa)) {
    missing_taxa <- setdiff(manual_label_taxa, unique(plot_df$.taxon))
    if (length(missing_taxa) > 0) {
      message(
        "Manual label taxa not shown because they are not in plotted taxa: ",
        paste(missing_taxa, collapse = ", ")
      )
    }
  }
  
  if (!is.null(manual_label_features)) {
    missing_features <- setdiff(manual_label_features, unique(plot_df[[feature_col]]))
    if (length(missing_features) > 0) {
      message(
        "Manual label features not shown because they are not in plotted taxa: ",
        paste(missing_features, collapse = ", ")
      )
    }
  }
  
  ## 9. 确定实际展示的分类层级
  ## 不能直接使用 levels(plot_df$.rank)，因为 factor levels 可能包含未展示层级。
  ## 这里根据 plot_df 中真实存在的层级重新确定 x 轴顺序。
  rank_levels <- rank_order[rank_order %in% unique(as.character(plot_df$.rank))]
  n_rank <- length(rank_levels)
  
  plot_df$.rank <- factor(plot_df$.rank, levels = rank_levels)
  plot_df$.x_id <- as.numeric(plot_df$.rank)
  
  ## 10. 生成 x 轴位置和 jitter 坐标
  ## jitter 坐标只生成一次，并同时赋给点和标签。
  ## 这样标签可以准确对应到被 jitter 后的点。
  set.seed(123)
  plot_df$.x_jit <- plot_df$.x_id +
    stats::runif(nrow(plot_df), -jitter_width, jitter_width)
  
  label_df$.x_jit <- plot_df$.x_jit[
    match(label_df[[feature_col]], plot_df[[feature_col]])
  ]
  label_df$.x_id <- plot_df$.x_id[
    match(label_df[[feature_col]], plot_df[[feature_col]])
  ]
  
  ## 11. 将正向和负向效应拆成上下两个面板
  ## coef > 0 的结果绘制在上方面板；## coef < 0 的结果绘制在下方面板；
  ## 中间面板仅用于显示分类层级色条，避免点和色条重叠。
  plot_top_df <- dplyr::filter(plot_df, .coef > 0)
  plot_bottom_df <- dplyr::filter(plot_df, .coef < 0)
  label_df_top <- dplyr::filter(label_df, .coef > 0)
  label_df_bottom <- dplyr::filter(label_df, .coef < 0)
  
  ## 12. 生成每个分类层级的灰色背景块
  ## 灰色背景高度根据该分类层级内点的最大/最小 coef 自动确定，
  ## 而不是铺满整个坐标系。
  coef_range <- range(plot_df$.coef, na.rm = TRUE)
  bg_pad <- max(0.3, diff(coef_range) * 0.04)
  
  ## 根据每个分类层级的点范围生成灰色背景块
  ## top 面板从 0 延伸到该层级最大正效应；
  ## bottom 面板从该层级最小负效应延伸到 0。
  make_bg <- function(x, side = c("top", "bottom")) {
    side <- match.arg(side)
    if (nrow(x) == 0) {
      return(data.frame(xmin = numeric(), xmax = numeric(), ymin = numeric(), ymax = numeric()))
    }
    
    out <- x |>
      dplyr::group_by(.rank, .x_id) |>
      dplyr::summarise(
        xmin = dplyr::first(.x_id) - 0.45,
        xmax = dplyr::first(.x_id) + 0.45,
        ymin = if (side == "top") 0 else min(.coef, na.rm = TRUE) - bg_pad,
        ymax = if (side == "top") max(.coef, na.rm = TRUE) + bg_pad else 0,
        .groups = "drop"
      )
    out
  }
  
  bg_df_top <- make_bg(plot_top_df, "top")
  bg_df_bottom <- make_bg(plot_bottom_df, "bottom")
  
  ## 13. 生成中间分类层级色条
  ## 色条颜色使用 Spectral 风格调色板，并根据实际展示的层级数量自动插值。
  if (is.null(strip_height)) {
    strip_half_height <- max(0.20, diff(coef_range) * 0.03)
  } else {
    strip_half_height <- strip_height / 2
  }
  
  strip_df <- data.frame(
    .rank = factor(rank_levels, levels = rank_levels),
    xmin = seq_len(n_rank) - 0.45,
    xmax = seq_len(n_rank) + 0.45,
    ymin = -strip_half_height,
    ymax = strip_half_height,
    fill = grDevices::colorRampPalette(strip_palette)(n_rank)
  )
  
  strip_label_df <- data.frame(
    x = seq_len(n_rank),
    y = 0,
    label = rank_levels
  )
  
  ## 14. 添加标签函数
  ## 优先使用 ggrepel 避免标签重叠；
  ## 如果未安装 ggrepel，则退回到普通 geom_text。
  if (is.null(ref)) ref <- "Reference"
  if (is.null(title)) title <- paste0("MaAsLin3 effect: ", metadata, " = ", value)
  
  p_prefix <- dplyr::if_else(grepl('qval', p_col), 'P.adjust', 'P.value')
  
  sig_labels <- c(
    paste0(p_prefix, " < ", p_cutoff), paste0(p_prefix, " >= ", p_cutoff)
  )
  
  ## 给指定面板添加标签
  ## 标签数据已经继承了点的 jitter 坐标，因此标签会对准对应点。
  add_label <- function(p, label_df) {
    if (nrow(label_df) == 0) return(p)
    
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p + ggrepel::geom_text_repel(
        data = label_df,
        ggplot2::aes(x = .x_jit, y = .coef, label = .label),
        size = 3.5, max.overlaps = Inf, box.padding = 0.35,
        point.padding = 0.25, min.segment.length = 0, 
        seed = 123, show.legend = FALSE
      )
    } else {
      p + ggplot2::geom_text(
        data = label_df,
        ggplot2::aes(x = .x_jit, y = .coef, label = .label),
        size = 3.2, vjust = -0.6, show.legend = FALSE
      )
    }
  }
  
  ## 15. 定义点图层
  ## 统一点图层样式
  ## sig/ns 由 p_col 和 p_cutoff 判断，并通过 color 颜色区分。
  point_layer <- function(data) {
    ggplot2::geom_point(
      data = data,
      ggplot2::aes(x = .x_jit, y = .coef, color = .sig_group),
      size = point_size, alpha = 0.95
    )
  }
  
  base_x <- ggplot2::scale_x_continuous(
    breaks = seq_len(n_rank), labels = rank_levels, expand = c(0.02, 0.02)
  )
  
  color_scale <- ggplot2::scale_color_manual(
    values = sig_colors, breaks = c("sig", "ns"), labels = sig_labels, name = NULL
  )
  
  ## 16. 绘制上方面板
  ## 仅展示 coef > 0 的菌，表示与 value 方向正相关。
  p_top <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = bg_df_top,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = "grey92", color = NA
    ) +
    point_layer(plot_top_df) +
    color_scale + 
    base_x +
    ggplot2::labs(
      title = title,
      subtitle = paste0(
        metadata, " = ", value, " vs ", ref,
        "; point top = ", point_top_n,
        "; label top = ", label_top_n,
        "; p cutoff = ", p_cutoff
      ),
      x = NULL, y = NULL
    ) +
    ggplot2::coord_cartesian(
      ylim = c(0, max(c(bg_df_top$ymax, 0.5), na.rm = TRUE))
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      axis.line.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.length.y = grid::unit(2, "mm"),
      axis.ticks.y = ggplot2::element_line(linewidth = .5, colour = 'black'),
      axis.line.y = ggplot2::element_line(linewidth = .5, colour = 'black'),
      axis.text.y = ggplot2::element_text(size = 12, color = "black"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10),
      plot.margin = ggplot2::margin(t = 5.5, r = 5.5, b = 0, l = 5.5),
      legend.position = "top"
    )
  
  p_top <- add_label(p_top, label_df_top)
  
  ## 17. 绘制中间分类层级色条
  ## 该面板只显示分类层级名称，不展示散点。
  p_mid <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = strip_df,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = .rank),
      inherit.aes = FALSE, color = "black", linewidth = 0.4 
    ) +
    ggplot2::scale_fill_manual(
      values = stats::setNames(strip_df$fill, strip_df$.rank),
      guide = "none"
    ) +
    ggplot2::geom_text(
      data = strip_label_df,
      ggplot2::aes(x = x, y = y, label = label),
      inherit.aes = FALSE, size = 4, fontface = "bold"
    ) +
    base_x +
    ggplot2::coord_cartesian(ylim = c(-strip_half_height, strip_half_height)) +
    ggplot2::labs(y = "MaAsLin3 coefficient") +
    ggplot2::theme_void() +
    ggplot2::theme(
      axis.title.y = ggplot2::element_text(size = 12, color = "black", angle = 90),
      plot.margin = ggplot2::margin(t = 0, r = 5.5, b = 0, l = 5.5)
    )
  
  ## 18. 绘制下方面板
  ## 仅展示 coef < 0 的菌，表示与参考组方向相关。
  p_bottom <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = bg_df_bottom,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = "grey92", color = NA
    ) +
    point_layer(plot_bottom_df) +
    color_scale +
    base_x +
    ggplot2::labs(x = "Taxonomic level", y = NULL) +
    ggplot2::coord_cartesian(
      ylim = c(min(c(bg_df_bottom$ymin, -0.5), na.rm = TRUE), 0)
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      axis.line.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.length.y = grid::unit(2, "mm"),
      axis.ticks.y = ggplot2::element_line(linewidth = .5, colour = 'black'),
      axis.line.y = ggplot2::element_line(linewidth = .5, colour = 'black'),
      axis.text.y = ggplot2::element_text(size = 12, color = "black"),
      plot.margin = ggplot2::margin(t = 0, r = 5.5, b = 5.5, l = 5.5),
      legend.position = "none"
    )
  
  p_bottom <- add_label(p_bottom, label_df_bottom)
  
  ## 19. 拼接三个面板
  ## 上：正向效应；中：分类层级色条；下：负向效应。
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required.")
  }
  
  p <- p_top / p_mid / p_bottom +
    patchwork::plot_layout(heights = c(4.8, 0.6, 4.8))
  
  ## 20. 保存作图数据
  ## plot_df 为最终展示的点；
  ## label_df 为最终标注的标签。
  ## 可用 attr(p, "plot_df") 和 attr(p, "label_df") 取出。
  attr(p, "plot_df") <- plot_df
  attr(p, "label_df") <- label_df
  
  ## 21. 输出当前绘图参数和数据量，便于检查筛选结果
  if (isTRUE(verbose)) {
    message(
      "\nMaAsLin3 multi-level effect plot\n",
      "  metadata: ", metadata, "\n",
      "  value: ", value, "\n",
      "  reference: ", ref, "\n",
      "  p column: ", p_col, "\n",
      "  p cutoff: ", p_cutoff, "\n",
      "  plotted levels: ", paste(rank_levels, collapse = ", "), "\n",
      "  total plotted taxa: ", nrow(plot_df), "\n",
      "  labeled taxa: ", nrow(label_df)
    )
  }
  
  return(p)
}


#### 20260619 plot_maaslin3_abundance ####
#' Plot Maaslin3 Abundance utility
#'
#' `plot_maaslin3_abundance()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 绘制 MaAsLin3 abundance model 的 effect 结果。
#'
#' @param data An input data frame or compatible object.
#' @param p_col Name of the `p_col` input column.
#' @return A plot object; analysis data or models may also be stored as attributes.
#' @export
plot_maaslin3_abundance <- function(
    data, p_col = 'qval_joint'
) {
  
  data <- data.frame(data, check.names = FALSE)
  
  ## 必须手动指定 p_col
  if (missing(p_col) || is.null(p_col) || length(p_col) != 1) {
    stop(
      'Please specify p_col, for example: p_col = "qval_joint".'
    )
  }
  
  if (!p_col %in% colnames(data)) {
    stop("p_col not found in data: ", p_col)
  }
  
  ## 检查必要列
  need_cols <- c(
    "feature", "coef", "null_hypothesis", "model"
  )
  
  miss_cols <- setdiff(need_cols, colnames(data))
  
  if (length(miss_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(miss_cols, collapse = ", ")
    )
  }
  
  ## 整理数据
  plot_df <- data |>
    dplyr::mutate(
      feature = as.character(feature),
      model = tolower(as.character(model)),
      coef = as.numeric(coef),
      null_hypothesis = as.numeric(null_hypothesis),
      p_value = as.numeric(.data[[p_col]])
    ) |>
    dplyr::filter(
      model == "abundance",
      is.finite(coef),
      is.finite(null_hypothesis),
      is.finite(p_value)
    ) |> 
    dplyr::mutate(
      model_label = dplyr::case_when(
        model == "abundance" ~ "Abundance",
        TRUE ~ model
      )
    )
  
  if (nrow(plot_df) == 0) {
    stop("No valid abundance rows remained.")
  }
  
  ## feature 顺序：输入表格第一行显示在图最上方
  feature_levels <- dplyr::arrange(plot_df, coef) |>
    dplyr::pull(feature)
  
  plot_df <- plot_df |>
    dplyr::mutate(
      feature = factor(feature, levels = feature_levels),
      y = match(feature, feature_levels),
      sig = -log10(pmax(p_value, .Machine$double.xmin))
    )
  
  ## 颜色刻度：pretty() 自动生成
  sig_range <- range(plot_df$sig, na.rm = TRUE)
  
  if (diff(sig_range) == 0) {
    fill_breaks <- sig_range
  } else {
    fill_breaks <- pretty(sig_range, n = 4)
    fill_breaks <- fill_breaks[
      fill_breaks >= sig_range[1] & fill_breaks <= sig_range[2]
    ]
    if (length(fill_breaks) == 0) {
      fill_breaks <- sig_range
    }
  }
  
  fill_labels <- formatC(10^(-fill_breaks), format = "g", digits = 2)
  
  fill_title <- if (
    grepl("qval|fdr", p_col, ignore.case = TRUE)
  ) {
    expression(Abundance~P[FDR])
  } else {
    expression(Abundance~P)
  }
  
  ## Null hypothesis 线：每种 model 一个或多个 null 值
  null_df <- plot_df |>
    dplyr::select(model_label, null_hypothesis) |>
    dplyr::distinct()
  
  ## 作图
  p <- ggplot2::ggplot() +
    ## 0 线：只是方向参考，不是 MaAsLin3 的 null
    ggplot2::geom_vline(
      xintercept = 0, linetype = "longdash", linewidth = 0.45, color = "black"
    ) +
    ## null_hypothesis：abundance 实线
    ggplot2::geom_vline(
      data = null_df,
      ggplot2::aes(xintercept = null_hypothesis, linetype = model_label),
      linewidth = 0.55, color = "black"
    ) +
    ## association：abundance 圆形
    ggplot2::geom_point(
      ggplot2::aes(x = coef, y = y, fill = sig), data = plot_df,
      size = 3.4, stroke = 0.8, color = "black", shape = 21
    ) +
    ggplot2::scale_linetype_manual(
      name = "Null hypothesis", values = c(Abundance = "solid")
    ) +
    ggplot2::scale_fill_gradient(
      low = "#f2e6f2", high = "#a200a6", breaks = fill_breaks,
      labels = fill_labels, name = fill_title
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq_along(feature_levels),
      labels = feature_levels,
      expand = ggplot2::expansion(add = c(0.5, 0.5))
    ) +
    ggplot2::labs(
      x = expression(beta~coefficient), y = "Feature"
    ) +
    ggplot2::guides(
      linetype = ggplot2::guide_legend(order = 1),
      fill = ggplot2::guide_colorbar(order = 2,reverse = TRUE)
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "grey88", linewidth = 0.4),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = NA, color = 'black', linewidth = .4),
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(color = "black"),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = .4),
      axis.ticks.length = grid::unit(1.6, 'mm'),
      legend.position = "right",
      legend.key = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(color = "black"),
      legend.text = ggplot2::element_text(color = "black")
    )
  
  attr(p, "plot_df") <- plot_df
  
  return(p)
}

#### 20260619 plot_maaslin3_both ####
#' Plot Maaslin3 both utility
#'
#' `plot_maaslin3_both()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 联合展示 MaAsLin3 abundance 和 prevalence 模型结果。
#'
#' @param data An input data frame or compatible object.
#' @param p_col Name of the `p_col` input column.
#' @param title Optional plot or result title.
#' @param point_size Numeric setting for `point_size`.
#' @param aspect_ratio Panel aspect ratio passed to `ggplot2::theme()`.
#' @return A plot object; analysis data or models may also be stored as attributes.
#' @export
plot_maaslin3_both <- function(
  data, p_col = 'qval_individual', title = NULL, point_size = 3.8,
  aspect_ratio = 2
) {
  
  data <- data.frame(data, check.names = FALSE)
  
  ## 指定 p_col
  if (is.null(p_col) || length(p_col) != 1) {
    stop(
      'Please specify p_col, for example: p_col = "qval_joint".'
    )
  }
  
  if (!p_col %in% colnames(data)) {
    stop("p_col not found in data: ", p_col)
  }
  
  ## 检查必要列
  need_cols <- c(
    "feature", "coef", "null_hypothesis", "model"
  )
  
  miss_cols <- setdiff(need_cols, colnames(data))
  
  if (length(miss_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(miss_cols, collapse = ", ")
    )
  }
  
  ## 整理数据
  plot_df <- data |>
    dplyr::mutate(
      feature = as.character(feature),
      model = tolower(as.character(model)),
      coef = as.numeric(coef),
      null_hypothesis = as.numeric(null_hypothesis),
      p_value = as.numeric(.data[[p_col]])
    ) |>
    dplyr::filter(
      model %in% c("abundance", "prevalence"),
      is.finite(coef),
      is.finite(null_hypothesis),
      is.finite(p_value)
    ) |>
    dplyr::mutate(
      model_label = dplyr::case_when(
        model == "abundance" ~ "Abundance",
        model == "prevalence" ~ "Prevalence",
        TRUE ~ model
      )
    )
  
  if (nrow(plot_df) == 0) {
    stop("No valid abundance/prevalence rows remained.")
  }
  
  ## feature 顺序：输入表格第一行显示在图最上方
  feature_levels <- dplyr::group_by(plot_df, feature) |>
    dplyr::slice_max(order_by = abs(coef),n = 1) |> 
    # dplyr::summarise(coef = mean(coef)) |>
    dplyr::arrange(coef) |>
    dplyr::pull(feature)
  
  plot_df <- plot_df |>
    dplyr::mutate(
      feature = factor(feature, levels = feature_levels),
      y = match(feature, feature_levels),
      sig = -log10(pmax(p_value, .Machine$double.xmin))
    )
  
  df_abun <- plot_df |>
    dplyr::filter(model == "abundance")
  
  df_prev <- plot_df |>
    dplyr::filter(model == "prevalence")
  
  ## Null hypothesis 线：每种 model 一个或多个 null 值
  null_df <- plot_df |>
    dplyr::select(model_label, null_hypothesis) |>
    dplyr::distinct()
  
  ## abundance 颜色刻度
  if (nrow(df_abun) > 0) {
    abun_range <- range(df_abun$sig, na.rm = TRUE)
    if (diff(abun_range) == 0) {
      abun_breaks <- abun_range
    } else {
      abun_breaks <- pretty(abun_range, n = 4)
      abun_breaks <- abun_breaks[
        abun_breaks >= abun_range[1] & abun_breaks <= abun_range[2]
      ]
      if (length(abun_breaks) == 0) {
        abun_breaks <- abun_range
      }
    }
    abun_labels <- formatC(10^(-abun_breaks), format = "g", digits = 2)
  }
  
  ## prevalence 颜色刻度
  if (nrow(df_prev) > 0) {
    prev_range <- range(df_prev$sig, na.rm = TRUE)
    if (diff(prev_range) == 0) {
      prev_breaks <- prev_range
    } else {
      prev_breaks <- pretty(prev_range, n = 4)
      prev_breaks <- prev_breaks[
        prev_breaks >= prev_range[1] & prev_breaks <= prev_range[2]
      ]
      if (length(prev_breaks) == 0) {
        prev_breaks <- prev_range
      }
    }
    prev_labels <- formatC(10^(-prev_breaks), format = "g", digits = 2)
  }
  
  abun_title <- ifelse(
    grepl("qval|fdr", p_col, ignore.case = TRUE), 
    expression(Abundance~P[FDR]),
    expression(Abundance~P)
  )
  
  prev_title <- ifelse(
    grepl("qval|fdr", p_col, ignore.case = TRUE),
    expression(Prevalence~P[FDR]),
    expression(Prevalence~P)
  )
  
  ## 基础图
  p <- ggplot2::ggplot() +
    
    ## MaAsLin3 null_hypothesis 线
    ggplot2::geom_vline(
      data = null_df,
      ggplot2::aes(xintercept = null_hypothesis, linetype = model_label),
      linewidth = 0.55, color = "black"
    )
  
  ## abundance 点：圆形，紫色
  if (nrow(df_abun) > 0) {
    
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(x = coef, y = y, shape = model_label, fill = sig),
        data = df_abun, size = point_size, stroke = 0.8, color = "black"
      ) +
      ggplot2::scale_fill_gradient(
        low = "#f2e6f2", high = "#a200a6", breaks = abun_breaks,
        labels = abun_labels, name = abun_title
      )
  }
  
  ## prevalence 点：三角形，青绿色
  if (nrow(df_prev) > 0) {

    p <- p +
      ggnewscale::new_scale_fill() +
      ggplot2::geom_point(
        ggplot2::aes(x = coef, y = y, shape = model_label, fill = sig),
        data = df_prev, size = point_size, stroke = 0.8, color = "black"
      ) +
      ggplot2::scale_fill_gradient(
        low = "#e6f2f2", high = "#0f9d9a", breaks = prev_breaks,
        labels = prev_labels, name = prev_title
      )
  }
  
  p <- p +
    ggplot2::scale_linetype_manual(
      name = "Null hypothesis", drop = FALSE,
      values = c(Abundance = "solid", Prevalence = "dashed")
    ) +
    ggplot2::scale_shape_manual(
      name = "Association", drop = FALSE,
      values = c(Abundance = 21, Prevalence = 24)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq_along(feature_levels), labels = feature_levels,
      expand = ggplot2::expansion(add = c(0.5, 0.5))
    ) +
    ggplot2::labs(x = expression(beta~coefficient), y = "Feature", title = title) +
    ggplot2::guides(
      linetype = ggplot2::guide_legend(order = 1),
      shape = ggplot2::guide_legend(
        order = 2, override.aes = list(fill = NA, color = "black")
      )
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "grey88", linewidth = 0.4),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = NA, color = 'black', linewidth = .4),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(color = 'black', face = 'bold', hjust = .5, size = 13),
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(color = "black"),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = .4),
      axis.ticks.length = grid::unit(1.6, 'mm'),
      legend.position = "right",
      legend.key = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(color = "black"),
      legend.text = ggplot2::element_text(color = "black"),
      aspect.ratio = aspect_ratio
    )
  
  attr(p, "plot_df") <- plot_df
  
  return(p)
}

#### 20260502 load_KEGG_info ####
#' Load KEGG Info utility
#'
#' `load_KEGG_info()` provides a reusable mengR workflow with input validation and standardized
#'   output.
#'
#' Chinese summary: 读取并整理本地 KEGG 注释表。
#'
#' @param database Numeric setting for `database`.
#' @param level Numeric setting for `level`.
#' @param relation Relationship table or relation type used to connect database identifiers.
#' @param keep_desc Whether KEGG descriptions are retained in the returned table.
#' @param add_prefix Logical control for `add_prefix`.
#' @param distinct Whether duplicated KEGG mapping records are removed.
#' @return A result object described in the Details section.
#' @export
load_KEGG_info <- function(
    database = .mengR_db_file(
      'KEGG', 'KO', 'KO_level_A_B_C_D_Description'
    ),
    level = c('A', 'B', 'C', 'D'), relation = NULL, 
    keep_desc = TRUE, add_prefix = TRUE, distinct = TRUE
  ) {
  if (!file.exists(database)) stop('KEGG annotation file not found: ', database)

  kegg_db <- utils::read.delim(
    database, sep = '\t', quote = '', check.names = FALSE, stringsAsFactors = FALSE
    )
  
  lv_map <- c(A = 'lvA', B = 'lvB', C = 'lvC', D = 'lvD')
  
  level <- toupper(level)
  
  if (!all(level %in% names(lv_map))) {
    stop('`level` should be selected from: A, B, C, D')
  }
  
  ## 如果指定 relation，则输出两个层级之间的对应关系
  ## 例如 relation = c('D', 'C') 表示 KO 到 pathway C 的关系
  if (!is.null(relation)) {
    relation <- toupper(relation)
    
    if (length(relation) != 2) {
      stop("`relation` should be length 2, e.g. c('D', 'C')")
    }
    
    if (!all(relation %in% names(lv_map))) {
      stop('`relation` should be selected from: A, B, C, D')
    }
    
    cols <- lv_map[relation]
    
    if (isTRUE(keep_desc)) {
      desc_cols <- paste0(cols, 'des')
      cols <- as.vector(rbind(cols, desc_cols))
      cols <- cols[cols %in% colnames(kegg_db)]
    }
    
    result <- kegg_db[, cols, drop = FALSE]
    
    if (isTRUE(distinct)) {
      result <- dplyr::distinct(result)
    }
    
    return(result)
  }
  
  ## 否则输出指定 level 的信息
  cols <- lv_map[level]
  
  if (isTRUE(keep_desc)) {
    desc_cols <- paste0(cols, 'des')
    cols <- as.vector(rbind(cols, desc_cols))
    cols <- cols[cols %in% colnames(kegg_db)]
  }
  
  result <- kegg_db[, cols, drop = FALSE]
  
  if (isTRUE(distinct)) {
    result <- dplyr::distinct(result)
  }
  
  result <- dplyr::mutate(
    result, 
    dplyr::across(dplyr::where(is.character), ~ gsub('\"', '', .x) ) 
    )
  
  return(result)
}

#### 20260502 profile_KEGG_trans ####
#' Profile KEGG Trans utility
#'
#' `profile_KEGG_trans()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 将 KO profile 映射并聚合到 KEGG A/B/C 层级。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param to Target taxonomy rank or identifier type produced by the conversion.
#' @param trans_ra Whether abundances are converted to relative abundance before analysis.
#' @param base Numeric setting for `base`.
#' @param rownames_fmt Format used to construct row names after KEGG conversion.
#' @param remove_unknown Logical control for `remove_unknown`.
#' @param filter Optional expression or criterion used to retain converted records.
#' @param database Numeric setting for `database`.
#' @param split_multi Whether records containing multiple identifiers are expanded into separate rows.
#' @param ... Additional arguments passed to the underlying function.
#' @return A result object described in the Details section.
#' @export
profile_KEGG_trans <- function(
    profile, to = c('A', 'B', 'C'), trans_ra = FALSE, base = 100,
    rownames_fmt = c('both', 'id', 'name'), 
    remove_unknown = FALSE, filter = NULL,
    database = .mengR_db_file(
      'KEGG', 'KO', 'KO_level_A_B_C_D_Description'
    ),
    split_multi = TRUE, ...) {
  
  to <- match.arg(to)
  rownames_fmt <- match.arg(rownames_fmt)
  
  lv_col <- c(A = 'lvA', B = 'lvB', C = 'lvC')[[to]]
  desc_col <- paste0(lv_col, 'des')
  
  profile <- data.frame(profile, check.names = FALSE)

  if (!file.exists(database)) stop('KEGG annotation file not found: ', database)
  
  sample_cols <- colnames(profile)
  
  kegg_db <- utils::read.delim(
    database, sep = '\t', quote = '', check.names = FALSE, 
    stringsAsFactors = FALSE, ...
  )
  
  if (!is.null(filter)) {
    kegg_db <- kegg_db |>
      dplyr::filter(lvAdes %in% filter)
  }
  
  ## KO -> KEGG level map
  kegg_map <- kegg_db |>
    dplyr::select(lvD, name = dplyr::all_of(lv_col)) |>
    dplyr::distinct() |>
    dplyr::filter(lvD %in% rownames(profile))
  
  ## KEGG ID -> KEGG description map
  kegg_name <- kegg_db |>
    dplyr::select(name = dplyr::all_of(lv_col),
                  des = dplyr::all_of(desc_col)) |>
    dplyr::distinct()
  
  ## profile long mapping
  data <- profile |>
    tibble::rownames_to_column('lvD') |>
    dplyr::left_join(kegg_map, by = 'lvD') |>
    dplyr::mutate(name = ifelse(is.na(name), 'Unknown', name))
  
  if (isTRUE(remove_unknown)) {
    data <- data |>
      dplyr::filter(name != 'Unknown')
  }
  
  ## 如果一个 KO 对应多个 KEGG 分类，平均分配丰度
  if (isTRUE(split_multi)) {
    data <- data |>
      dplyr::group_by(lvD) |>
      dplyr::mutate(.n = dplyr::n()) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(sample_cols),
          ~ .x / .n
        )
      ) |>
      dplyr::select(-.n)
  }
  
  ## 汇总到 KEGG A/B/C 层级
  data <- data |>
    dplyr::select(-lvD) |>
    dplyr::group_by(name) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(sample_cols),
        ~ sum(.x, na.rm = TRUE)
      ),
      .groups = 'drop'
    )
  
  ## 修改行名格式
  if (rownames_fmt %in% c('name', 'both')) {
    
    data <- data |>
      dplyr::left_join(kegg_name, by = 'name')
    
    data$des[is.na(data$des)] <- data$name[is.na(data$des)]
    
    if (rownames_fmt == 'name') {
      data$rowname <- gsub('"', '', data$des)
    }
    
    if (rownames_fmt == 'both') {
      data$rowname <- ifelse(
        data$name == 'Unknown',
        'Unknown',
        paste0(data$name, '|', gsub('"', '', data$des) )
      )
    }
    
    data <- data |>
      dplyr::select(rowname, dplyr::all_of(sample_cols))
    
  } else {
    
    data <- data |>
      dplyr::rename(rowname = name)
  }
  
  ## 如果不同 ID 对应同一个 description，再合并一次
  data <- data |>
    dplyr::group_by(rowname) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(sample_cols),
        ~ sum(.x, na.rm = TRUE)
      ),
      .groups = 'drop'
    ) |>
    tibble::column_to_rownames('rowname') |>
    data.frame(check.names = FALSE)
  
  ## 转成相对丰度百分比
  if (isTRUE(trans_ra)) {
    cs <- colSums(data, na.rm = TRUE)
    data <- (sweep(data, 2, cs, '/') * base) |> 
      as.data.frame(check.names = FALSE)
  }
  
  rownames(data) <- gsub('\'', '', rownames(data))
  
  return(data)
}

#### 20250423 tidy_LEfSe ####
#' Tidy LEfSe utility
#'
#' `tidy_LEfSe()` provides a reusable mengR workflow with input validation and standardized
#'   output.
#'
#' Chinese summary: 统一 LEfSe 结果列和 taxonomy 字段，便于后续筛选绘图。
#'
#' @param data An input data frame or compatible object.
#' @return A result object described in the Details section.
#' @export
tidy_LEfSe <- function(data) {
  RENAMES <- data.frame(
    full = c('domain','kingdom','phylum','class','order','family','genus','species', 'strain'),
    abbr = c('d__','k__','p__','c__','o__', 'f__','g__','s__','t__')
    )
  
  data <- data.frame(data, row.names = NULL) |> 
    dplyr::mutate(name = purrr::map_vec(Taxa, ~ {
      full_name = unlist(strsplit(.x, split = '\\|'), use.names = FALSE)
      length_name = length(full_name)
      full_name[[length_name]] } ),
      taxa_level = RENAMES$full[match(substr(name, 1, 3), RENAMES$abbr)]) |> 
    dplyr::select(name, taxa_level, full_name = Taxa, comparison = Comparison, 
                  LDA, pval = P.unadj, qval = P.adj, qlab = Significance, 
                  method = Method, enriched = Group)
  return(data)
}

#### 20250423 tidy_CAZyme_profile ####
#' Tidy CAZyme Profile utility
#'
#' `tidy_CAZyme_profile()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 整理 CAZyme profile 名称并按目标分类层级聚合。
#'
#' @param data An input data frame or compatible object.
#' @return A result object described in the Details section.
#' @export
tidy_CAZyme_profile <- function(data) {
  ## 1. 去除 CAZyme 编号后缀，并合并同名条目
  profile_df <- tibble::rownames_to_column(.as_df(data), 'name')
  profile_df$name <- gsub('_\\d+', '', profile_df$name)
  profile_df <- stats::aggregate(. ~ name, data = profile_df, FUN = sum)

  ## 2. 用“|”拆分复合名称，并把丰度平均分配给各名称
  result_list <- lapply(seq_len(nrow(profile_df)), function(row_idx) {
    name_vec <- strsplit(profile_df$name[row_idx], '\\|')[[1]]
    value_vec <- as.numeric(profile_df[row_idx, -1, drop = TRUE]) /
      length(name_vec)
    value_mat <- matrix(
      rep(value_vec, times = length(name_vec)),
      nrow = length(name_vec), byrow = TRUE,
      dimnames = list(NULL, colnames(profile_df)[-1])
    )
    data.frame(new_name = name_vec, value_mat, check.names = FALSE)
  })
  result_df <- dplyr::bind_rows(result_list)
  result_df <- stats::aggregate(. ~ new_name, data = result_df, FUN = sum)
  tibble::column_to_rownames(result_df, 'new_name')
}

