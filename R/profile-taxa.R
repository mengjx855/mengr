#### Jin-Xin Meng, 20220918, 20260523, v0.1.6 ####

# 20250206: fix some bug.
# 20250804: update function.
# 20251217: tidy function.
# 20260502: add profile_mpa() to generate metaphlan-like profile.
# 20260523: update taxa_trans(), plot_compos(), plot_compos_multiple(),
#           plot_compos_manual(), and plot_taxa_boxplot().


#### 20260502 profile_mpa ####
# 生成 MetaPhlAn-like profile 表
# profile: 行为 feature，列为 sample 的丰度表
# taxa: feature 注释表，至少包含 feature_col，以及 d__/p__/c__/o__/f__/g__/s__ 等分类列
# feature_col: taxa 中用于匹配 profile 行名的列，默认 name
# mode:
#   all     输出所有分类层级路径
#   lowest  输出每个 feature 可获得的最低层级路径
#   species 只输出含 species 层级的路径
# normalize:
#   none     不转换
#   relative 转为相对丰度
#   percent  转为百分比

#' Profile Mpa utility
#'
#' `profile_mpa()` provides a reusable mengR workflow with input validation and standardized
#'   output.
#'
#' Chinese summary: 处理 MetaPhlAn 风格 profile，拆分层级并输出指定分类水平。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param taxa Taxon names or taxonomy data used by the operation.
#' @param feature_col Name of the feature-identifier column.
#' @param mode Processing mode; supported values are shown in Usage.
#' @param normalize Whether to enable the normalize behavior.
#' @return A result object described in the Details section.
#' @export
profile_mpa <- function(
  profile, taxa, feature_col = "name", mode = c("all", "lowest", "species"),
  normalize = c("none", "percent", "relative")
) {
  mode <- match.arg(mode)
  normalize <- match.arg(normalize)

  profile <- as.data.frame(profile, check.names = FALSE)
  taxa <- as.data.frame(taxa, check.names = FALSE)

  if (!feature_col %in% colnames(taxa)) {
    stop("tax should contain feature column: ", feature_col)
  }

  sample_cols <- colnames(profile)

  ## 1. 可选：把 profile 转成相对丰度或百分比
  if (normalize != "none") {
    if (normalize == "relative") {
      profile <- profile_trans_ra(profile, base = 1, digits = NULL)
    }

    if (normalize == "percent") {
      profile <- profile_trans_ra(profile, base = 100, digits = NULL)
    }
  }

  profile$.feature <- rownames(profile)

  ## 2. 解析 taxonomy
  ##    不依赖列名，而是根据 d__/p__/c__/o__/f__/g__/s__ 自动识别
  taxa_long <- tidyr::pivot_longer(
    taxa,
    cols = -dplyr::all_of(feature_col),
    names_to = "rank_col", values_to = "taxon"
  )

  taxa_long <- taxa_long |>
    dplyr::mutate(taxon = trimws(as.character(taxon))) |>
    dplyr::filter(!is.na(taxon), taxon != "", grepl("^[dkpcofgs]__", taxon)) |>
    dplyr::mutate(
      rank = substr(taxon, 1, 1),
      rank = dplyr::case_when(
        rank == "d" ~ "k",
        TRUE ~ rank
      ),
      taxon = sub("^[dkpcofgs]__", "", taxon),
      taxon = paste0(rank, "__", taxon)
    ) |>
    dplyr::select(.feature = dplyr::all_of(feature_col), rank, taxon) |>
    dplyr::distinct(.feature, rank, .keep_all = TRUE)

  taxa_wide <- tidyr::pivot_wider(
    taxa_long,
    names_from = rank, values_from = taxon
  )

  ranks <- c("k", "p", "c", "o", "f", "g", "s")

  for (r in ranks) {
    if (!r %in% colnames(taxa_wide)) {
      taxa_wide[[r]] <- NA_character_
    }
  }

  data <- dplyr::left_join(
    profile, taxa_wide,
    by = ".feature"
  )

  ## 3. 构建 clade_name
  make_path <- function(x) {
    x <- x[!is.na(x) & x != ""]
    if (length(x) == 0) {
      return(NA_character_)
    }
    paste(x, collapse = "|")
  }

  ## 4. 生成 MetaPhlAn-like 表
  if (mode == "all") {
    result_list <- list()

    for (i in seq_along(ranks)) {
      rr <- ranks[seq_len(i)]

      clade <- apply(data[, rr, drop = FALSE], 1, make_path)

      path_df <- data.frame(
        clade_name = clade,
        data[, sample_cols, drop = FALSE],
        check.names = FALSE
      )

      path_df <- path_df[!is.na(path_df$clade_name), , drop = FALSE]

      path_df <- path_df |>
        dplyr::group_by(clade_name) |>
        dplyr::summarise(
          dplyr::across(
            dplyr::all_of(sample_cols),
            ~ sum(.x, na.rm = TRUE)
          ),
          .groups = "drop"
        )

      result_list[[i]] <- path_df
    }

    result <- dplyr::bind_rows(result_list) |>
      dplyr::distinct(clade_name, .keep_all = TRUE)
  }

  if (mode == "lowest") {
    clade <- apply(data[, ranks, drop = FALSE], 1, make_path)

    result <- data.frame(
      clade_name = clade,
      data[, sample_cols, drop = FALSE],
      check.names = FALSE
    )

    result <- result[!is.na(result$clade_name), , drop = FALSE]

    result <- result |>
      dplyr::group_by(clade_name) |>
      dplyr::summarise(
        dplyr::across(
          dplyr::all_of(sample_cols),
          ~ sum(.x, na.rm = TRUE)
        ),
        .groups = "drop"
      )
  }

  if (mode == "species") {
    data <- data[!is.na(data$s), , drop = FALSE]

    clade <- apply(data[, ranks, drop = FALSE], 1, make_path)

    result <- data.frame(
      clade_name = clade,
      data[, sample_cols, drop = FALSE],
      check.names = FALSE
    )

    result <- result[!is.na(result$clade_name), , drop = FALSE]

    result <- result |>
      dplyr::group_by(clade_name) |>
      dplyr::summarise(
        dplyr::across(
          dplyr::all_of(sample_cols),
          ~ sum(.x, na.rm = TRUE)
        ),
        .groups = "drop"
      )
  }

  result <- tibble::column_to_rownames(result, "clade_name")

  return(result)
}

#### 20220918 taxa_split ####
# taxonomy: data.frame containing name + taxonomy columns
# taxonomy_rename: optional, e.g. c(name = "OTU_ID", taxonomy = "Taxonomy")
#' Taxa Split utility
#'
#' `taxa_split()` provides a reusable mengR workflow with input validation and standardized
#'   output.
#'
#' Chinese summary: 将完整 taxonomy 字符串拆分为 domain 至 species 等分类列。
#'
#' @param taxonomy Feature taxonomy table used for annotation or aggregation.
#' @param taxonomy_rename Display or identifier name for taxonomy rename.
#' @param sep Field separator used when reading or writing a text file.
#' @param from Source taxonomy rank or identifier type to be converted.
#' @param to Target taxonomy rank or identifier type produced by the conversion.
#' @param na_fill Value used to replace missing observations before analysis.
#' @param rm_suffix Whether suffix text is removed from split taxonomy labels.
#' @param ... Additional arguments passed to the underlying function.
#' @return A result object described in the Details section.
#' @export
taxa_split <- function(
  taxonomy, taxonomy_rename = NULL, sep = ";", from = "d", to = "s",
  na_fill = "Unknown", rm_suffix = FALSE, ...
) {
  taxonomy <- data.frame(taxonomy, check.names = FALSE)

  ## 1. Optional rename
  if (!is.null(taxonomy_rename)) {
    if (is.null(names(taxonomy_rename))) {
      stop("`taxonomy_rename` should be a named vector, e.g. c(name = 'ID', taxonomy = 'Taxon').")
    }
    taxonomy <- dplyr::rename(taxonomy, dplyr::all_of(taxonomy_rename))
  }

  ## 2. Check required columns
  if (!all(c("name", "taxonomy") %in% colnames(taxonomy))) {
    stop("`taxonomy` must contain two columns: `name` and `taxonomy`.")
  }

  ## 3. Taxonomic rank definition
  taxa_info <- c(
    domain = "d", phylum = "p", class = "c", order = "o",
    family = "f", genus = "g", species = "s", strain = "t"
  )

  ## allow full rank names, e.g. from = "phylum"
  normalize_rank <- function(x) {
    x <- tolower(x)
    if (x %in% names(taxa_info)) {
      return(unname(taxa_info[x]))
    }
    if (x %in% taxa_info) {
      return(x)
    }
    if (x == "k") {
      return("d")
    } # compatible with k__Bacteria
    stop("Unknown taxonomic rank: ", x)
  }

  from <- normalize_rank(from)
  to <- normalize_rank(to)

  from_idx <- match(from, taxa_info)
  to_idx <- match(to, taxa_info)

  if (from_idx > to_idx) {
    stop("`from` should be an upstream rank of `to`.")
  }

  rank_idx <- from_idx:to_idx
  rank_names <- names(taxa_info)[rank_idx]
  rank_letters <- unname(taxa_info[rank_idx])
  rank_prefix <- paste0(rank_letters, "__")

  unknown_values <- c(
    "", "NA", "NaN", "NULL",
    "Unassigned", "Unknown", "Uncultured",
    "unassigned", "unknown", "uncultured"
  )

  ## 4. Remove GTDB-style suffix, e.g. Firmicutes_A -> Firmicutes
  remove_suffix <- function(x) {
    x <- sub("_[A-Z]$", "", x)
    x <- sub("_[A-Z](?=\\s)", "", x, perl = TRUE)
    x
  }

  ## 5. Clean one taxon
  clean_taxon <- function(x, prefix) {
    x <- trimws(as.character(x))
    ## remove existing rank prefix, e.g. p__Bacteroidota -> Bacteroidota
    x <- sub("^[a-zA-Z]__", "", x)
    x <- trimws(x)
    if (is.na(x) || x %in% unknown_values) x <- na_fill
    x <- paste0(prefix, x)
    if (isTRUE(rm_suffix)) x <- remove_suffix(x)
    x
  }

  ## 6. Split one taxonomy string
  split_one <- function(x) {
    out <- stats::setNames(paste0(rank_prefix, na_fill), rank_names)

    if (is.na(x) || trimws(x) %in% unknown_values) {
      return(out)
    }

    parts <- trimws(strsplit(as.character(x), split = sep, fixed = TRUE)[[1]])

    ## Case A: taxonomy has prefixes, e.g. d__Bacteria;p__Bacteroidota
    has_prefix <- grepl("^[a-zA-Z]__", parts)

    if (any(has_prefix)) {
      for (p in parts[has_prefix]) {
        r <- tolower(substr(p, 1, 1))
        if (r == "k") r <- "d"

        if (r %in% rank_letters) {
          idx <- match(r, rank_letters)
          out[idx] <- clean_taxon(p, rank_prefix[idx])
        }
      }
    } else {
      ## Case B: taxonomy has no prefixes, e.g. Bacteria;Bacteroidota;Bacteroidia
      full_rank_names <- names(taxa_info)
      full_rank_prefix <- paste0(unname(taxa_info), "__")

      full_out <- stats::setNames(
        paste0(full_rank_prefix, na_fill),
        full_rank_names
      )

      n <- min(length(parts), length(full_out))

      for (i in seq_len(n)) {
        full_out[i] <- clean_taxon(parts[i], full_rank_prefix[i])
      }

      out <- full_out[rank_names]
    }

    out
  }

  ## 7. Apply to all rows
  tax_mat <- t(vapply(
    taxonomy$taxonomy,
    split_one,
    FUN.VALUE = character(length(rank_names))
  ))

  result <- data.frame(
    name = taxonomy$name,
    tax_mat,
    check.names = FALSE,
    row.names = NULL
  )

  colnames(result) <- c("name", rank_names)

  return(result)
}

#### 20220918 taxa_trans ####
# 按 taxonomy 注释将 profile 汇总到指定分类层级
# profile: 行为 feature，列为 sample 的丰度表
# taxonomy: feature 注释表
# feature_col: taxonomy 中用于匹配 profile 行名的列，默认 name
# taxa_col: taxonomy 中用于汇总的分类列，默认 family
# sample_col/group_col: group 中样本列和分组列
# top_n: 保留丰度最高的前 top_n - 1 个分类，其余合并为 Other
# top_list: 指定保留的分类名称，其余合并为 Other
# out_all: 是否输出全部分类；TRUE 时忽略 top_n
# remove_unknown: 是否在计算相对丰度前删除 unknown/unclassified 等分类
# trans_ra: 是否转换为百分比相对丰度
# collapse_group: 是否按 group 合并样本

#' Taxa Trans utility
#'
#' `taxa_trans()` provides a reusable mengR workflow with input validation and standardized
#'   output.
#'
#' Chinese summary: 按 taxonomy 层级聚合 profile，选择 top taxa，并可按组或相对丰度转换。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param taxonomy Feature taxonomy table used for annotation or aggregation.
#' @param group A sample metadata table containing sample and group columns.
#' @param feature_col Name of the feature-identifier column.
#' @param taxa_col Name of the `taxa_col` input column.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param top_n Number of highest-ranking features or categories retained.
#' @param top_list Optional feature or category names retained regardless of abundance rank.
#' @param other_name Label assigned to features combined into the residual category.
#' @param out_all Whether all hierarchy levels or intermediate results are returned.
#' @param na_fill Value used to replace missing observations before analysis.
#' @param trans_ra Whether abundances are converted to relative abundance before analysis.
#' @param base Numeric setting for `base`.
#' @param digits Optional number of decimal digits retained.
#' @param collapse_group Whether samples within each group are aggregated before plotting.
#' @param method Analysis or summary method; supported values are shown in the usage.
#' @param remove_unknown Logical control for `remove_unknown`.
#' @param unknown_pattern Regular expression identifying unknown or unclassified annotations.
#' @return A result object described in the Details section.
#' @export
taxa_trans <- function(
  profile, taxonomy, group = NULL, feature_col = "name", taxa_col = "family",
  sample_col = "sample", group_col = "group", top_n = 12, top_list = NULL,
  other_name = "Other", out_all = FALSE, na_fill = "Unknown",
  trans_ra = FALSE, base = 100, digits = 8, collapse_group = FALSE,
  method = c("mean", "median", "sum"), remove_unknown = FALSE,
  unknown_pattern = "unknown|unclassified|unassigned|uncultured"
) {
  if (is.function(method)) {
    summary_fun <- method
  } else {
    method <- match.arg(method)
    summary_fun <- match.fun(method)
  }
  profile <- data.frame(profile, check.names = FALSE)
  taxonomy <- data.frame(taxonomy, check.names = FALSE)

  if (!all(c(feature_col, taxa_col) %in% colnames(taxonomy))) {
    stop(
      "taxonomy should contain columns: ",
      feature_col, " | ", taxa_col
    )
  }

  if (isTRUE(out_all)) {
    top_n <- 0
  }

  ## 1. 匹配 taxonomy 注释
  taxa <- taxonomy[[taxa_col]][match(rownames(profile), taxonomy[[feature_col]])]
  taxa <- as.character(taxa)
  taxa[is.na(taxa) | taxa == ""] <- na_fill

  data <- profile |>
    dplyr::mutate(taxa = taxa) |>
    dplyr::filter(!is.na(taxa), taxa != "")

  ## 2. 删除 unknown / unclassified 等未明确注释的分类
  ## 默认匹配 unknown, unclassified, unassigned, uncultured
  ## 如果想自定义，可修改 unknown_pattern
  if (isTRUE(remove_unknown)) {
    data <- data |>
      dplyr::filter(
        !grepl(unknown_pattern, taxa, ignore.case = TRUE)
      )
  }

  if (nrow(data) == 0) {
    warning("No features remained after filtering taxonomy.")
    return(data.frame())
  }

  ## 3. 按分类层级汇总
  data <- data |>
    stats::aggregate(. ~ taxa, data = _, FUN = sum, na.rm = TRUE) |>
    tibble::column_to_rownames("taxa") |>
    data.frame(check.names = FALSE)

  ## 4. top_n 合并
  if (top_n > 0 && is.null(top_list)) {
    top_taxa <- data.frame(
      taxa = rownames(data),
      value = rowSums(data, na.rm = TRUE)
    ) |>
      dplyr::arrange(dplyr::desc(value)) |>
      utils::head(n = top_n - 1) |>
      dplyr::pull(taxa)

    data <- data |>
      tibble::rownames_to_column("name") |>
      dplyr::mutate(name = ifelse(name %in% top_taxa, name, other_name)) |>
      stats::aggregate(. ~ name, data = _, FUN = sum, na.rm = TRUE) |>
      tibble::column_to_rownames("name")
  }

  ## 5. top_list 合并
  if (!is.null(top_list)) {
    data <- data |>
      tibble::rownames_to_column("name") |>
      dplyr::mutate(name = ifelse(name %in% top_list, name, other_name)) |>
      stats::aggregate(. ~ name, data = _, FUN = sum, na.rm = TRUE) |>
      tibble::column_to_rownames("name")
  }

  ## 6. 是否合并到 group 水平
  if (isTRUE(collapse_group)) {
    if (is.null(group)) {
      stop("missing group file")
    }

    group <- data.frame(group, check.names = FALSE)

    if (!all(c(sample_col, group_col) %in% colnames(group))) {
      stop(
        "group should contain columns: ",
        sample_col, " | ", group_col
      )
    }

    data <- data.frame(t(data), check.names = FALSE)
    sample_vec <- rownames(data)
    data <- data |>
      dplyr::mutate(
        group = group[[group_col]][match(sample_vec, group[[sample_col]])]
      ) |>
      dplyr::filter(!is.na(group)) |>
      stats::aggregate(. ~ group, data = _, FUN = summary_fun) |>
      tibble::column_to_rownames("group") |>
      t() |>
      data.frame(check.names = FALSE)
  }

  ## 7. 转换为百分比相对丰度
  if (isTRUE(trans_ra)) {
    data <- profile_trans_ra(profile = data, base = base, digits = digits)
  }

  data <- data[
    rowSums(data, na.rm = TRUE) != 0,
    colSums(data, na.rm = TRUE) != 0,
    drop = FALSE
  ]

  return(data)
}

#### 20220918 plot_compos ####
# 绘制分类组成堆叠柱状图
# display: group 或 sample
# to: 显示的分类层级
# top_n/top_list: 控制显示分类数量

#' Plot Compos utility
#'
#' `plot_compos()` provides a reusable mengR workflow with input validation and standardized
#'   output.
#'
#' Chinese summary: 绘制单个分类层级的样本或分组组成堆叠柱状图。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param taxonomy Feature taxonomy table used for annotation or aggregation.
#' @param group A sample metadata table containing sample and group columns.
#' @param display Whether composition is summarized by sample or by group.
#' @param feature_col Name of the feature-identifier column.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param to Target taxonomy rank or identifier type produced by the conversion.
#' @param top_n Number of highest-ranking features or categories retained.
#' @param top_list Optional feature or category names retained regardless of abundance rank.
#' @param group_level Optional order of group levels.
#' @param sample_level Optional order for `sample_level`.
#' @param width Numeric setting for `width`.
#' @param taxa_level Optional order for `taxa_level`.
#' @param taxa_color Color specification for `taxa_color`.
#' @param plot_title Logical control for `plot_title`.
#' @param x_text_angle Rotation angle, in degrees, for x-axis text.
#' @param remove_unknown Logical control for `remove_unknown`.
#' @param unknown_pattern Regular expression identifying unknown or unclassified annotations.
#' @return A plot object; analysis data or models may also be stored as attributes.
#' @export
plot_compos <- function(
  profile, taxonomy, group = NULL, display = c("group", "sample"),
  feature_col = "name", sample_col = "sample", group_col = "group",
  to = "family", top_n = 12, top_list = NULL, group_level = NULL,
  sample_level = NULL, width = .75, taxa_level = NULL, taxa_color = NULL,
  plot_title = NULL, x_text_angle = 90, remove_unknown = FALSE,
  unknown_pattern = "unknown|unclassified|unassigned|uncultured"
) {
  display <- match.arg(display)

  if (missing(profile) || missing(taxonomy)) {
    stop("missing profile or taxonomy.")
  }

  if (display == "group" && is.null(group)) {
    stop("missing group info.")
  }

  profile <- data.frame(profile, check.names = FALSE)
  taxonomy <- data.frame(taxonomy, check.names = FALSE)

  colors <- c(
    "#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", "#59A14F",
    "#8CD17D", "#B6992D", "#F1CE63", "#499894", "#86BCB6",
    "#E15759", "#FF9D9A", "#79706E", "#BAB0AC", "#D37295",
    "#FABFD2", "#B07AA1", "#D4A6C8", "#9D7660", "#D7B5A6"
  )

  other_name <- paste0(stringr::str_to_lower(stringr::str_sub(to, 1, 1)), "__Other")

  data <- taxa_trans(
    profile = profile, taxonomy = taxonomy, group = group,
    feature_col = feature_col, taxa_col = to,
    sample_col = sample_col, group_col = group_col,
    top_n = top_n, top_list = top_list, other_name = other_name,
    collapse_group = display == "group", trans_ra = TRUE,
    remove_unknown = remove_unknown, unknown_pattern = unknown_pattern
  )

  if (nrow(data) == 0) {
    stop("No taxa remained for plotting.")
  }

  if (is.null(taxa_level)) {
    taxa_level <- names(sort(rowSums(data, na.rm = TRUE), decreasing = TRUE)) |> rev()
  }

  if (is.null(taxa_color)) {
    taxa_color <- rep(colors, time = ceiling(nrow(data) / length(colors)))[seq_len(nrow(data))]
  }

  if (is.null(plot_title)) {
    plot_title <- paste0(stringr::str_to_sentence(to), " level composition")
  }

  fill_title <- paste0(stringr::str_to_sentence(to), " taxa")

  if (display == "group") {
    if (is.null(group_level)) group_level <- colnames(data)

    plot_df <- data |>
      tibble::rownames_to_column("name") |>
      tidyr::pivot_longer(-name, names_to = "group", values_to = "value") |>
      dplyr::mutate(
        group = factor(group, levels = group_level),
        name = factor(name, levels = taxa_level)
      )

    p <- ggpubr::ggbarplot(
      plot_df,
      x = "group", y = "value", fill = "name",
      color = "#000000", position = ggplot2::position_stack(),
      linewidth = .4, width = width, palette = taxa_color,
      x.text.angle = x_text_angle, legend = "right"
    )
  } else {
    if (is.null(sample_level)) sample_level <- colnames(data)

    plot_df <- data |>
      tibble::rownames_to_column("name") |>
      tidyr::pivot_longer(-name, names_to = "sample", values_to = "value") |>
      dplyr::mutate(
        sample = factor(sample, levels = sample_level),
        name = factor(name, levels = taxa_level)
      )

    p <- ggpubr::ggbarplot(
      plot_df,
      x = "sample", y = "value", fill = "name",
      color = "#000000", position = ggplot2::position_stack(),
      linewidth = .4, width = width, palette = taxa_color,
      x.text.angle = x_text_angle, legend = "right"
    )
  }

  p <- p +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(x = "", y = "Relative Abundance (%)", title = plot_title, fill = fill_title) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = .4, color = "#000000"),
      axis.ticks = ggplot2::element_line(linewidth = .4, color = "#000000"),
      axis.text = ggplot2::element_text(size = 10, color = "#000000"),
      axis.title = ggplot2::element_text(size = 10, color = "#000000"),
      plot.title = ggplot2::element_text(size = 10, color = "#000000", face = "bold", hjust = .5),
      legend.text = ggplot2::element_text(size = 10, color = "#000000", face = "italic"),
      legend.title = ggplot2::element_text(size = 10, color = "#000000"),
      legend.key.spacing.y = grid::unit(1, "mm"),
      panel.grid = ggplot2::element_blank()
    )

  return(p)
}

#### 20220918 plot_compos_multiple ####
# 一次绘制多个分类层级的组成图
# profile: 行为 feature，列为 sample 的丰度表
# taxonomy: feature 注释表
# group: 样本分组信息
# feature_col: taxonomy 中用于匹配 profile 行名的列
# sample_col/group_col: group 中样本列和分组列
# taxa_levels: 指定要绘制的分类层级；默认自动识别 taxonomy 中除 feature_col 外的列

#' Plot Compos Multiple utility
#'
#' `plot_compos_multiple()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 对多个 taxonomy 层级批量绘制组成图并组合输出。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param taxonomy Feature taxonomy table used for annotation or aggregation.
#' @param group A sample metadata table containing sample and group columns.
#' @param display Whether composition is summarized by sample or by group.
#' @param feature_col Name of the feature-identifier column.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param top_n Number of highest-ranking features or categories retained.
#' @param top_list Optional feature or category names retained regardless of abundance rank.
#' @param group_level Optional order of group levels.
#' @param sample_level Optional order for `sample_level`.
#' @param taxa_color Color specification for `taxa_color`.
#' @param width Numeric setting for `width`.
#' @param taxa_levels Numeric setting for `taxa_levels`.
#' @param nrow Number of rows used when arranging multiple plots.
#' @param ... Additional arguments passed to the underlying function.
#' @return A plot object; analysis data or models may also be stored as attributes.
#' @export
plot_compos_multiple <- function(
  profile, taxonomy, group = NULL, display = "group", feature_col = "name",
  sample_col = "sample", group_col = "group", top_n = 12, top_list = NULL,
  group_level = NULL, sample_level = NULL, taxa_color = NULL, width = .75,
  taxa_levels = NULL, nrow = 2, ...
) {
  if (missing(profile) || missing(taxonomy)) {
    stop("missing profile or taxonomy.")
  }

  taxonomy <- data.frame(taxonomy, check.names = FALSE)

  if (is.null(taxa_levels)) {
    ## 自动识别 taxonomy 中的分类层级列
    ## 支持大小写不一致，例如 Phylum / phylum / PHYLUM
    ## 也支持简写，例如 p / c / o / f / g / s
    rank_patterns <- list(
      # domain  = c("^domain$", "^kingdom$", "^superkingdom$", "^d$", "^k$", "^d__$", "^k__$"),
      phylum  = c("^phylum$", "^p$", "^p__$"),
      class   = c("^class$", "^c$", "^c__$"),
      order   = c("^order$", "^o$", "^o__$"),
      family  = c("^family$", "^f$", "^f__$"),
      genus   = c("^genus$", "^g$", "^g__$"),
      species = c("^species$", "^s$", "^s__$"),
      # strain  = c("^strain$", "^t$", "^t__$")
    )

    col_raw <- colnames(taxonomy)
    col_std <- tolower(col_raw)
    col_std <- gsub("\\s+", "_", col_std)
    col_std <- gsub("[.-]+", "_", col_std)
    col_std <- gsub("_+$", "", col_std)

    taxa_levels <- purrr::map_chr(
      rank_patterns, \(x) {
        hit <- which(grepl(paste(x, collapse = "|"), col_std))
        if (length(hit) == 0) NA_character_ else col_raw[hit[1]]
      }
    )

    taxa_levels <- taxa_levels[!is.na(taxa_levels)]
  }

  if (length(taxa_levels) == 0) {
    stop("No taxonomy rank columns were detected. Please specify taxa_levels manually.")
  }

  p <- purrr::map(
    taxa_levels, \(x) {
      plot_compos(
        profile = profile,
        taxonomy = taxonomy,
        group = group,
        display = display,
        feature_col = feature_col,
        sample_col = sample_col,
        group_col = group_col,
        top_n = top_n,
        top_list = top_list,
        to = x,
        group_level = group_level,
        sample_level = sample_level,
        taxa_color = taxa_color,
        width = width,
        ...
      ) |>
        suppressMessages()
    }
  ) |>
    cowplot::plot_grid(plotlist = _, nrow = nrow, align = "v")

  return(p)
}


#### 20220918 plot_compos_manual ####
# 手动输入已经汇总好的 profile 绘制组成图
# profile: 行为 taxa，列为 sample 或 group
# display:
#   group  先按 group 合并，再画组成图
#   sample 直接按 sample 画组成图
# sample_col/group_col: group 中样本列和分组列
# method: display = group 时，样本合并方法，默认 mean

#' Plot Compos Manual utility
#'
#' `plot_compos_manual()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 对已经整理好的组成数据绘制可精细控制的堆叠柱状图。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param group A sample metadata table containing sample and group columns.
#' @param display Whether composition is summarized by sample or by group.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param top_n Number of highest-ranking features or categories retained.
#' @param out_all Whether all hierarchy levels or intermediate results are returned.
#' @param group_level Optional order of group levels.
#' @param sample_level Optional order for `sample_level`.
#' @param taxa_level Optional order for `taxa_level`.
#' @param taxa_color Color specification for `taxa_color`.
#' @param width Numeric setting for `width`.
#' @param title Optional plot or result title.
#' @param fill_title Color specification for fill title.
#' @param aspect_ratio Panel aspect ratio passed to `ggplot2::theme()`.
#' @param x_text_angle Rotation angle, in degrees, for x-axis text.
#' @param other_name Label assigned to features combined into the residual category.
#' @param method Analysis or summary method; supported values are shown in the usage.
#' @param remove_unknown Logical control for `remove_unknown`.
#' @param unknown_pattern Regular expression identifying unknown or unclassified annotations.
#' @param base Numeric setting for `base`.
#' @param digits Optional number of decimal digits retained.
#' @param ... Additional arguments passed to the underlying function.
#' @return A plot object; analysis data or models may also be stored as attributes.
#' @export
plot_compos_manual <- function(
  profile, group = NULL, display = c("group", "sample"),
  sample_col = "sample", group_col = "group",
  top_n = 12, out_all = FALSE, group_level = NULL, sample_level = NULL,
  taxa_level = NULL, taxa_color = "category20", width = .75,
  title = NULL, fill_title = NULL, aspect_ratio = 1, x_text_angle = 0,
  other_name = "Others", method = "mean", remove_unknown = FALSE,
  unknown_pattern = "unknown|unclassified|unassigned|uncultured",
  base = 100, digits = 8, ...
) {
  display <- match.arg(display)

  if (display == "group" && is.null(group)) {
    stop("missing group info.")
  }

  if (isTRUE(out_all)) {
    top_n <- Inf
  }

  profile <- data.frame(profile, check.names = FALSE)
  profile <- profile[
    rowSums(profile, na.rm = TRUE) != 0,
    colSums(profile, na.rm = TRUE) != 0,
    drop = FALSE
  ]

  pals <- list(
    category20 = c(
      "#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", "#59A14F",
      "#8CD17D", "#B6992D", "#F1CE63", "#499894", "#86BCB6",
      "#E15759", "#FF9D9A", "#79706E", "#BAB0AC", "#D37295",
      "#FABFD2", "#B07AA1", "#D4A6C8", "#9D7660", "#D7B5A6"
    ),
    set3 = c(
      "#80b1d3", "#b3de69", "#fdb462", "#8dd3c7", "#bc80bd",
      "#fb8072", "#ffed6f", "#fccde5", "#bebada", "#ccebc5",
      "#ffffb3", "#d9d9d9"
    ),
    npj = c(
      "#3C5488", "#00A087", "#E64B35", "#4DBBD5", "#F39B7F",
      "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"
    ),
    custom1 = c(
      "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
      "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"
    ),
    custom2 = c(
      "#6B8E8D", "#C27D38", "#B85450", "#8E6C8A", "#7A9E59",
      "#D0A85C", "#7F8C8D", "#A67C52", "#C08A80", "#A4B494"
    ),
    custom3 = c(
      "#3B5BA5", "#57A0D3", "#6CC3A0", "#A1C349", "#F2C14E",
      "#F78154", "#D8576B", "#8D5A97", "#5D576B", "#B8B8B8"
    )
  )

  if (length(taxa_color) == 1) {
    if (!taxa_color %in% names(pals)) {
      stop(
        "taxa_color palette not match reference: ",
        paste(names(pals), collapse = ", ")
      )
    }
    colors <- pals[[taxa_color]]
  } else {
    colors <- taxa_color
  }

  # if (is.null(title)) {
  #   title <- 'Composition Analysis'
  # }

  ## 1. 删除 unknown / unclassified 等未明确注释分类
  ## 注意：如果删除 unknown 后需重新计算相对丰度，应先删除再转换
  if (isTRUE(remove_unknown)) {
    profile <- profile[
      !grepl(unknown_pattern, rownames(profile), ignore.case = TRUE), ,
      drop = FALSE
    ]
  }

  if (nrow(profile) == 0) {
    stop("No taxa remained after removing unknown taxa.")
  }

  ## 2. 准备数据
  if (display == "group") {
    data <- profile_collapse(
      profile = profile,
      group = group,
      sample_col = sample_col,
      group_col = group_col,
      method = method
    ) |>
      profile_trans_ra(base = base, digits = digits, remove_empty = TRUE)
  } else {
    data <- profile_trans_ra(
      profile = profile,
      base = base,
      digits = digits,
      remove_empty = TRUE
    )
  }

  ## 3. top taxa 合并
  ## 3. top taxa 合并
  if (is.finite(top_n) && top_n > 0) {
    taxa_sum <- rowSums(data, na.rm = TRUE)

    select_taxa <- names(sort(taxa_sum, decreasing = TRUE))
    select_taxa <- select_taxa[
      !grepl("unknown|unclass", select_taxa, ignore.case = TRUE)
    ]
    select_taxa <- utils::head(select_taxa, top_n - 1)

    data <- data |>
      tibble::rownames_to_column("name") |>
      dplyr::mutate(
        name = ifelse(name %in% select_taxa, name, other_name)
      ) |>
      dplyr::group_by(name) |>
      dplyr::summarise(
        dplyr::across(
          dplyr::everything(),
          \(x) sum(x, na.rm = TRUE)
        ),
        .groups = "drop"
      ) |>
      tibble::column_to_rownames("name")
  }

  if (is.null(taxa_level)) {
    taxa_level <- names(sort(rowSums(data, na.rm = TRUE), decreasing = TRUE)) |> rev()
  }

  taxa_color <- rep(
    colors,
    time = ceiling(nrow(data) / length(colors))
  )[seq_len(nrow(data))] |>
    rev()

  ## 4. 整理作图数据
  if (display == "group") {
    if (is.null(group_level)) {
      group_level <- colnames(data)
    }

    plot_df <- data |>
      tibble::rownames_to_column("name") |>
      tidyr::pivot_longer(-name, names_to = "group", values_to = "value") |>
      dplyr::mutate(
        group = factor(group, levels = group_level),
        name = factor(name, levels = taxa_level)
      )

    p <- ggpubr::ggbarplot(
      plot_df, "group", "value",
      fill = "name", color = "#000000",
      position = ggplot2::position_stack(), linewidth = .4, width = width,
      palette = taxa_color, x.text.angle = x_text_angle, legend = "right"
    )
  } else {
    if (is.null(sample_level)) {
      sample_level <- colnames(data)
    }

    plot_df <- data |>
      tibble::rownames_to_column("name") |>
      tidyr::pivot_longer(-name, names_to = "sample", values_to = "value") |>
      dplyr::mutate(
        sample = factor(sample, levels = sample_level),
        name = factor(name, levels = taxa_level)
      )

    p <- ggpubr::ggbarplot(
      plot_df, "sample", "value",
      fill = "name", color = "#000000",
      position = ggplot2::position_stack(), linewidth = .4, width = width,
      palette = taxa_color, x.text.angle = x_text_angle, legend = "right"
    )
  }

  p <- p +
    ggplot2::scale_y_continuous(expand = c(.02, 0)) +
    ggplot2::scale_x_discrete(expand = c(0, .5)) +
    ggplot2::labs(x = "", y = "Relative Abundance (%)", title = title, fill = fill_title) +
    ggplot2::theme(
      axis.ticks.length = grid::unit(2, "mm"),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(linewidth = .5, color = "#000000"),
      axis.text = ggplot2::element_text(size = 12, color = "#000000"),
      axis.title = ggplot2::element_text(size = 12, color = "#000000"),
      panel.background = ggplot2::element_rect(fill = NA, linewidth = .5, color = "#000000"),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = .5, size = 12, face = "bold"),
      legend.text = ggplot2::element_text(size = 12, color = "#000000", face = "italic"),
      legend.title = ggplot2::element_text(size = 12, color = "#000000"),
      legend.key.spacing.y = grid::unit(1, "mm"),
      aspect.ratio = aspect_ratio
    )

  return(p)
}

#### 20220918 plot_taxa_boxplot ####
# 对 profile 中每个 taxa/feature 绘制分组箱线图
# profile: 行为 taxa，列为 sample
# group: 样本分组信息，默认包含 sample 和 group
# trans: NULL / LOG10 / LOG2 / SQRT
# method: 差异检验方法，传递给 calcu_diff()

#' Plot Taxa Boxplot utility
#'
#' `plot_taxa_boxplot()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 对每个 taxa 绘制分组箱线图，并添加显著性比较。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param group A sample metadata table containing sample and group columns.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param group_level Optional order of group levels.
#' @param group_color Optional colors aligned to `group_level`.
#' @param trans Transformation applied to abundance values before plotting.
#' @param method Analysis or summary method; supported values are shown in the usage.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param aspect_ratio Panel aspect ratio passed to `ggplot2::theme()`.
#' @param legend_title Legend title; `NULL` uses a context-dependent default.
#' @param show_legend Logical control for `show_legend`.
#' @param x_text_angle Rotation angle, in degrees, for x-axis text.
#' @param ... Additional arguments passed to the underlying function.
#' @return A plot object; analysis data or models may also be stored as attributes.
#' @export
plot_taxa_boxplot <- function(
  profile, group,
  sample_col = "sample", group_col = "group",
  group_level = NULL, group_color = NULL,
  trans = NULL, method = c("wilcox", "anova", "t"),
  xlab = "", ylab = "Relative Abundance (%)",
  aspect_ratio = 1, legend_title = "group",
  show_legend = FALSE, x_text_angle = 0, ...
) {
  method <- match.arg(method)
  if (!all(c(sample_col, group_col) %in% colnames(group))) {
    stop("group should contain columns: ", sample_col, " | ", group_col)
  }

  group <- data.frame(group, check.names = FALSE)
  profile <- data.frame(profile, check.names = FALSE)
  profile <- profile[
    rowSums(profile, na.rm = TRUE) != 0,
    colSums(profile, na.rm = TRUE) != 0,
    drop = FALSE
  ]

  if (!is.null(trans)) {
    trans <- toupper(trans)

    if (trans == "LOG10") {
      profile <- profile_trans_log10(profile)
    } else if (trans == "LOG2") {
      profile <- profile_trans_log2(profile)
    } else if (trans %in% c("SQRT", "SQER")) {
      profile <- profile_trans_sqrt(profile)
    } else {
      stop("trans option: NULL | LOG10 | LOG2 | SQRT")
    }
  }

  feature_names <- rownames(profile)

  data <- data.frame(t(profile), check.names = FALSE) |>
    tibble::rownames_to_column("sample") |>
    dplyr::filter(sample %in% group[[sample_col]]) |>
    dplyr::left_join(
      dplyr::select(group,
        sample = dplyr::all_of(sample_col),
        group = dplyr::all_of(group_col)
      ),
      by = "sample"
    )

  if (is.null(group_level)) {
    group_level <- unique(data$group)
  }

  color <- c(
    "#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#ffff99",
    "#b15928", "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6"
  )

  if (is.null(group_color)) {
    group_color <- rep(
      color,
      times = ceiling(length(group_level) / length(color))
    )[seq_along(group_level)]
  }

  p_list <- purrr::map(feature_names, \(x) {
    plot_df <- data |>
      dplyr::select(sample, group, value = dplyr::all_of(x)) |>
      dplyr::mutate(group = factor(group, levels = group_level))

    comparisons <- calcu_diff(plot_df, value ~ group, method = method, ...) |>
      dplyr::filter(pval < 0.05) |>
      dplyr::pull(comparison) |>
      stringr::str_split("_vs_")

    p <- ggpubr::ggboxplot(
      plot_df, "group", "value",
      fill = "group",
      linewidth = .4, width = .6, outlier.shape = NA,
      palette = group_color, legend = "right",
      x.text.angle = x_text_angle
    ) +
      ggplot2::geom_jitter(
        ggplot2::aes(fill = group),
        size = 2.5, width = .25, shape = 21, height = 0
      ) +
      ggplot2::labs(x = xlab, y = ylab, color = legend_title, title = x) +
      ggplot2::theme(
        axis.ticks.length = grid::unit(1.7, "mm"),
        axis.line = ggplot2::element_line(linewidth = .4, color = "#000000"),
        axis.ticks = ggplot2::element_line(linewidth = .4, color = "#000000"),
        axis.text = ggplot2::element_text(size = 12, color = "#000000"),
        axis.title = ggplot2::element_text(size = 12, color = "#000000"),
        plot.title = ggplot2::element_text(hjust = .5, size = 14, face = "bold.italic"),
        legend.text = ggplot2::element_text(size = 12, color = "#000000"),
        legend.title = ggplot2::element_text(size = 12, color = "#000000"),
        panel.grid = ggplot2::element_blank(),
        aspect.ratio = aspect_ratio
      ) +
      ggplot2::guides(color = "none")

    if (length(comparisons) > 0) {
      p <- p +
        ggsignif::geom_signif(
          comparisons = comparisons,
          step_increase = .05,
          textsize = 5,
          tip_length = .02,
          vjust = .7,
          linewidth = .4,
          map_signif_level = TRUE
        )
    }

    if (isFALSE(show_legend)) {
      p <- p + ggplot2::guides(fill = "none")
    }

    return(p)
  })

  return(p_list)
}
