#### Jin-Xin Meng, jinxmeng@zju.edu.cn, 20211029, 20260916 ####

# 20230101: update function 'calcu_alpha()'.
# 20231204: update function 'check_file_name()' was deprecated.
# 20240304: fix some bug
# 20250404: 修改函数的某些参数名称，plot_alpha() 函数中 添加 add_ref_line 参数
# 20260527: update function.
# 20260819: add functions about beta-diversity from the plot_PCoA.R script.
# 20260916: standardize documentation and rename internal data-frame variables to the `*_df` style without changing returned component names.


#### calcu_alpha ####
# 计算 alpha diversity
# profile: 行为 feature，列为 sample 的丰度表
# method: richness / observed / chao1 / ace / shannon / simpson / pielou / gc / pd
# tree: method = 'pd' 时需要的系统发育树
# base: shannon 和 pielou 的对数底数
# value_col: 输出 alpha 指标列名

#' Calcu Alpha utility
#'
#'
#' 计算 richness、Shannon、Simpson、Chao1、ACE 等 alpha diversity 指标。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param method Analysis or summary method; supported values are shown in the usage.
#' @param tree Phylogenetic tree required for UniFrac distances.
#' @param base Scaling constant used for relative-abundance output.
#' @param value_col Name of the `value_col` input column.
#' @return A sample-level data frame containing the requested alpha-diversity indices.
#' @export
calcu_alpha <- function(
  profile, method = c(
    "richness", "observed", "chao1", "ace",
    "shannon", "simpson", "pielou", "gc", "pd"
  ),
  tree = NULL, base = exp(1), value_col = "value"
) {
  method <- match.arg(method)

  profile <- data.frame(profile, check.names = FALSE)
  profile <- as.matrix(profile)
  suppressWarnings(storage.mode(profile) <- "numeric")

  ## vegan / picante 默认行为 sample，列为 feature
  profile <- t(profile)

  if (method == "richness") {
    result <- rowSums(profile > 0, na.rm = TRUE)
  } else if (method == "observed") {
    result <- vegan::estimateR(ceiling(profile))[1, ]
  } else if (method == "chao1") {
    result <- vegan::estimateR(ceiling(profile))[2, ]
  } else if (method == "ace") {
    result <- vegan::estimateR(ceiling(profile))[4, ]
  } else if (method == "shannon") {
    result <- vegan::diversity(profile, index = "shannon", base = base)
  } else if (method == "simpson") {
    result <- vegan::diversity(profile, index = "simpson")
  } else if (method == "pielou") {
    observed <- vegan::estimateR(ceiling(profile))[1, ]
    shannon <- vegan::diversity(profile, index = "shannon", base = base)
    result <- shannon / log(observed, base = base)
    result[!is.finite(result)] <- NA_real_
  } else if (method == "gc") {
    result <- 1 - rowSums(profile == 1, na.rm = TRUE) / rowSums(profile, na.rm = TRUE)
    result[!is.finite(result)] <- NA_real_
  } else if (method == "pd") {
    if (is.null(tree)) {
      stop("tree should be provided when method = 'pd'.")
    }

    pd <- picante::pd(profile, tree, include.root = FALSE)
    result <- pd[, 1]
    names(result) <- rownames(pd)
  }

  data <- data.frame(
    sample = names(result),
    value = as.numeric(result),
    row.names = NULL,
    check.names = FALSE
  )

  colnames(data)[2] <- value_col

  return(data)
}

#### plot_alpha ####
# 绘制 alpha diversity 箱线图
# data: alpha diversity 结果表，默认包含 sample / value
# group: 样本分组表，默认包含 sample / group
# sample_col: 样本列名
# value_col: alpha 指标列名
# group_col: 分组列名
# group_level: 指定分组顺序
# group_color: 指定分组颜色
# sort_value: NULL / asc / desc，是否按组内 median 排序
# add_ref_line: 指定某一组的 median 作为参考虚线
# show_diff: 是否显示组间差异
# method: wilcox / t.test / anova 等，传给 calcu_diff/stat_compare_means

#' Plot Alpha utility
#'
#'
#' 绘制 alpha diversity 分组图，并支持排序、显著性及参考线。
#'
#' @param data An input data frame or compatible object.
#' @param group A sample metadata table containing sample and group columns.
#' @param sample_col Name of the sample-identifier column.
#' @param value_col Name of the `value_col` input column.
#' @param group_col Name of the grouping column.
#' @param group_level Optional order of group levels.
#' @param group_color Optional colors aligned to `group_level`.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param title Optional plot or result title.
#' @param aspect_ratio Panel aspect ratio passed to `ggplot2::theme()`.
#' @param show_grid Whether to draw panel grid lines.
#' @param show_jitter Whether to overlay jittered sample points.
#' @param rotate_x_text Whether to rotate x-axis text by 45 degrees.
#' @param coord_flip Whether to exchange the x and y axes with `ggplot2::coord_flip()`.
#' @param show_diff Whether to add pairwise significance comparisons.
#' @param method Analysis or summary method; supported values are shown in the usage.
#' @param sort_value Optional ascending or descending ordering of plotted values.
#' @param add_ref_line Whether to draw the reference line specified by `ref_line`.
#' @param ... Additional arguments passed to the boxplot layer or significance test.
#' @return A ggplot-compatible plot object; computed data or fitted objects are retained as attributes when applicable.
#' @export
plot_alpha <- function(
  data, group, sample_col = "sample", value_col = "value",
  group_col = "group", group_level = NULL, group_color = NULL,
  xlab = "", ylab = "", title = "", aspect_ratio = 1, show_grid = TRUE,
  show_jitter = TRUE, rotate_x_text = FALSE, coord_flip = FALSE,
  show_diff = TRUE, method = c("wilcox", "t"), sort_value = NULL,
  add_ref_line = NULL, ...
) {
  method <- match.arg(method)

  data <- data.frame(data, check.names = FALSE)
  group <- data.frame(group, check.names = FALSE)

  if (!all(c(sample_col, value_col) %in% colnames(data))) {
    stop("data should contain columns: ", sample_col, " | ", value_col)
  }

  if (!all(c(sample_col, group_col) %in% colnames(group))) {
    stop("group should contain columns: ", sample_col, " | ", group_col)
  }

  data <- data |>
    dplyr::select(
      sample = dplyr::all_of(sample_col),
      value = dplyr::all_of(value_col)
    ) |>
    dplyr::mutate(value = as.numeric(value))

  group <- group |>
    dplyr::select(
      sample = dplyr::all_of(sample_col),
      group = dplyr::all_of(group_col)
    )

  plot_df <- dplyr::left_join(data, group, by = "sample") |>
    dplyr::filter(!is.na(group), is.finite(value))

  if (nrow(plot_df) == 0) {
    stop("No matched samples between data and group.")
  }

  if (is.null(group_level)) {
    group_level <- unique(plot_df$group)
  }

  group_level <- group_level[group_level %in% unique(plot_df$group)]

  if (length(group_level) == 0) {
    stop("No valid group remained for plotting.")
  }

  if (!is.null(sort_value)) {
    if (!sort_value %in% c("asc", "desc")) {
      message("sort_value should be 'asc', 'desc' or NULL.")
      sort_value <- NULL
    }

    if (sort_value == "asc") {
      group_level <- stats::aggregate(value ~ group, plot_df, median) |>
        dplyr::arrange(value) |>
        dplyr::pull(group) |>
        as.character()
    }

    if (sort_value == "desc") {
      group_level <- stats::aggregate(value ~ group, plot_df, median) |>
        dplyr::arrange(dplyr::desc(value)) |>
        dplyr::pull(group) |>
        as.character()
    }
  }

  if (!is.null(add_ref_line) && !add_ref_line %in% group_level) {
    message("add_ref_line not existing.")
    add_ref_line <- NULL
  }

  group_color <- .resolve_group_colors(group_level, group_color)

  plot_df <- plot_df |>
    dplyr::mutate(group = factor(group, levels = group_level))

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(group, value, fill = group)) +
    ggplot2::geom_boxplot(
      width = .618, linewidth = .4,
      outlier.shape = NA, show.legend = FALSE, ...
    ) +
    ggplot2::scale_fill_manual(values = group_color, drop = FALSE) +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    ggpubr::theme_pubr() +
    ggplot2::theme(
      aspect.ratio = aspect_ratio,
      axis.ticks.length = grid::unit(2, "mm"),
      axis.ticks = ggplot2::element_line(linewidth = .4, color = "black"),
      plot.title = ggplot2::element_text(hjust = .5, size = 12, face = "bold")
    )

  if (!is.null(add_ref_line)) {
    ref_value <- stats::median(
      plot_df$value[plot_df$group == add_ref_line],
      na.rm = TRUE
    )

    p <- p +
      ggplot2::geom_hline(
        yintercept = ref_value,
        linetype = "dashed",
        linewidth = .4,
        color = "#000000"
      )
  }

  if (isTRUE(show_jitter)) {
    p <- p +
      ggplot2::geom_jitter(
        ggplot2::aes(color = group),
        size = .7, width = .2,
        show.legend = FALSE
      ) +
      ggplot2::scale_color_manual(values = group_color, drop = FALSE)
  }

  if (isTRUE(show_diff)) {
    diff <- calcu_diff(
      data = plot_df,
      formula = value ~ group,
      method = method
    )

    comparisons <- diff |>
      dplyr::filter(pval < 0.05) |>
      dplyr::pull(comparison) |>
      strsplit(split = "_vs_")

    if (length(comparisons) > 0) {
      plot_method <- dplyr::case_when(
        method == "wilcox" ~ "wilcox.test",
        method == "t" ~ "t.test"
      )

      p <- p +
        ggsignif::geom_signif(
          comparisons = comparisons,
          step_increase = .09, textsize = 2.5, test = plot_method,
          tip_length = .02, vjust = .1, size = .4, parse = TRUE,
          map_signif_level = \(p) {
            p_lab <- ifelse(
              p < 0.001,
              scales::scientific(p, digits = 2),
              sprintf("%.3f", p)
            )
            sprintf("italic(p)~'='~'%s'", p_lab)
          }
        )
    }
  }

  if (isTRUE(show_grid)) {
    p <- p +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(color = "grey88", linewidth = .4)
      )
  }

  if (isTRUE(rotate_x_text)) {
    p <- p +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = .5))
  }

  if (isTRUE(coord_flip)) {
    p <- p + ggplot2::coord_flip()
  }

  return(p)
}

#### calcu_distance ####
# 计算样本间距离
# profile: 行为 feature，列为 sample 的丰度表
# dist_method: 距离方法，如 bray / jaccard / euclidean / unifrac
# tree: dist_method = 'unifrac' 时需要系统发育树
# weighted: UniFrac 是否加权
# transform: 计算普通距离前是否转换；NULL 表示不转换
# remove_empty: 是否删除全 0 feature 和 sample

#' Calcu Distance utility
#'
#'
#' 计算 Bray、Jaccard、Euclidean、UniFrac 等样本距离，可先转换 profile。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param dist_method Distance method; available values are validated with `match.arg()`.
#' @param tree Phylogenetic tree required for UniFrac distances.
#' @param weighted Whether to calculate weighted rather than unweighted UniFrac.
#' @param transform Optional transformation applied before analysis.
#' @param remove_empty Whether to remove samples whose total abundance is zero.
#' @param ... Additional arguments passed to the selected distance function.
#' @return A `dist` object containing pairwise sample distances.
#' @export
calcu_distance <- function(
  profile,
  dist_method = c(
    "bray", "jaccard", "euclidean", "manhattan", "canberra",
    "kulczynski", "gower", "altGower", "morisita", "horn",
    "mountford", "raup", "binomial", "chao", "cao",
    "mahalanobis", "unifrac"
  ),
  tree = NULL, weighted = TRUE,
  transform = c(
    "hellinger", "total", "max", "frequency", "normalize",
    "range", "rank", "rrank", "standardize", "pa",
    "chi.square", "log", "clr", "rclr", "alr"
  ),
  remove_empty = TRUE, ...
) {
  dist_method <- .match_distance_method(dist_method)
  transform <- .match_transform_method(transform)

  profile_mat <- as.matrix(.as_profile_df(profile, numeric = TRUE))

  profile_mat[!is.finite(profile_mat)] <- 0

  ## 1. 删除全 0 feature 和 sample
  if (isTRUE(remove_empty)) {
    profile_mat <- profile_mat[
      rowSums(profile_mat, na.rm = TRUE) != 0,
      colSums(profile_mat, na.rm = TRUE) != 0,
      drop = FALSE
    ]
  }

  if (nrow(profile_mat) == 0) {
    stop("No valid features remained for distance calculation.")
  }

  if (ncol(profile_mat) < 2) {
    stop("Distance calculation requires at least two valid samples.")
  }

  ## 2. UniFrac 距离
  if (dist_method == "unifrac") {
    if (is.null(tree)) {
      stop("Need phylogenetic tree for UniFrac distance.")
    }

    ps <- phyloseq::phyloseq(
      phyloseq::otu_table(profile_mat, taxa_are_rows = TRUE),
      phyloseq::phy_tree(tree)
    )

    distance <- phyloseq::UniFrac(
      ps,
      weighted = weighted, ...
    )
  } else {
    ## 3. 普通距离：转成 sample × feature
    sample_mat <- t(profile_mat)

    if (!is.null(transform)) {
      sample_mat <- vegan::decostand(sample_mat, method = transform)
    }

    sample_mat[!is.finite(sample_mat)] <- 0

    ## 转换后再删一次全 0 sample，避免 jaccard 报 empty rows
    if (isTRUE(remove_empty)) {
      sample_mat <- sample_mat[
        rowSums(sample_mat, na.rm = TRUE) != 0, ,
        drop = FALSE
      ]
    }

    if (nrow(sample_mat) < 2) {
      stop("Distance calculation requires at least two valid samples.")
    }

    distance <- vegan::vegdist(
      sample_mat,
      method = dist_method, na.rm = TRUE, ...
    )
  }
  return(distance)
}

#### calcu_beta ####
# 从距离矩阵/vegan dist 对象中提取每个 group 内部样本两两距离
# distance: dist 对象或距离矩阵
# metadata: 样本分组信息
# sample_col: metadata 中样本列名
# group_col: metadata 中分组列名
# drop_na_group: 是否删除 group 为 NA 或空字符的样本

#' Calcu Beta utility
#'
#'
#' 基于距离矩阵执行整体 beta-diversity 组间检验。
#'
#' @param distance A precomputed distance object; when supplied, it takes precedence over `profile`.
#' @param metadata A metadata or annotation data frame.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param drop_na_group Whether to remove samples with missing or empty group labels.
#' @return A data frame of within-group sample pairs and their distances.
#' @export
calcu_beta <- function(
  distance, metadata, sample_col = "sample",
  group_col = "group", drop_na_group = TRUE
) {
  dist_mat <- as.matrix(distance)

  metadata <- data.frame(metadata, check.names = FALSE)

  meta <- metadata[, c(sample_col, group_col)]
  colnames(meta) <- c("sample", "group")

  meta$sample <- as.character(meta$sample)
  meta$group <- as.character(meta$group)

  ## 只保留距离矩阵中存在的样本
  meta <- meta[
    meta$sample %in% rownames(dist_mat) &
      !is.na(meta$group) &
      meta$group != "", ,
    drop = FALSE
  ]

  group_level <- unique(meta$group)

  beta_list <- list()

  for (g in group_level) {
    sample_g <- meta$sample[meta$group == g]

    if (length(sample_g) < 2) {
      next
    }

    pair_g <- utils::combn(sample_g, 2)

    idx_x <- match(pair_g[1, ], rownames(dist_mat))
    idx_y <- match(pair_g[2, ], colnames(dist_mat))

    beta_list[[g]] <- data.frame(
      sample_x = pair_g[1, ],
      sample_y = pair_g[2, ],
      value = dist_mat[cbind(idx_x, idx_y)],
      group = g,
      check.names = FALSE
    )
  }

  beta_df <- dplyr::bind_rows(beta_list)

  rownames(beta_df) <- NULL

  return(beta_df)
}

#### calcu_adjusted_r2 ####
# 根据 adonis2 结果计算 adjusted R2

#' Calcu Adjusted R2 utility
#'
#'
#' 从 adonis 类结果中计算或提取 adjusted R-squared。
#'
#' @param adonis_object Object returned by a PERMANOVA/adonis calculation.
#' @return A numeric adjusted R-squared value.
#' @export
calcu_adjusted_r2 <- function(adonis_object) {
  n_observations <- adonis_object$Df[nrow(adonis_object)] + 1
  d_freedom <- adonis_object$Df[1]
  r2 <- adonis_object$R2[1]

  adjusted_r2 <- vegan::RsquareAdj(
    r2, n_observations, d_freedom
  )

  return(adjusted_r2)
}

#### calcu_pairwise_adonis ####
# 两两 PERMANOVA
# profile: 行为 feature，列为 sample 的丰度表
# group: 样本分组表
# sample_col/group_col: group 中样本列和分组列
# group_level: 指定分组顺序

#' Calcu Pairwise Adonis utility
#'
#'
#' 对各组组合执行 pairwise PERMANOVA，并校正 P 值。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param group A sample metadata table containing sample and group columns.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param group_level Optional order of group levels.
#' @param dist_method Distance method; available values are validated with `match.arg()`.
#' @param permutations Number of permutations used by the significance test.
#' @param add_plab Whether to add formatted significance labels.
#' @param ... Additional arguments passed to `vegan::adonis2()`.
#' @return A data frame containing pairwise PERMANOVA statistics, P values, adjusted P values, and optional labels.
#' @export
calcu_pairwise_adonis <- function(
  profile, group, sample_col = "sample", group_col = "group",
  group_level = NULL, dist_method = "bray", permutations = 999,
  add_plab = TRUE, ...
) {
  profile <- data.frame(profile, check.names = FALSE)
  group <- data.frame(group, check.names = FALSE)

  if (!all(c(sample_col, group_col) %in% colnames(group))) {
    stop("group should contain columns: ", sample_col, " | ", group_col)
  }

  if (is.null(group_level)) {
    group_level <- unique(group[[group_col]])
  }

  group_level <- group_level[group_level %in% unique(group[[group_col]])]

  if (length(group_level) < 2) {
    stop("At least two groups are required.")
  }

  profile <- profile |>
    vegan::decostand(MARGIN = 2, method = "total")

  data <- purrr::map_dfr(
    utils::combn(group_level, m = 2, simplify = FALSE), \(x) {
      meta <- group |>
        dplyr::filter(.data[[group_col]] %in% x) |>
        dplyr::filter(.data[[sample_col]] %in% colnames(profile)) |>
        dplyr::select(
          sample = dplyr::all_of(sample_col),
          group = dplyr::all_of(group_col)
        )

      if (nrow(meta) < 3 || length(unique(meta$group)) < 2) {
        return(data.frame(
          comparison = paste0(x, collapse = "_vs_"),
          r2 = NA_real_,
          r2adj = NA_real_,
          pval = NA_real_,
          check.names = FALSE
        ))
      }

      profile_x <- profile[, meta$sample, drop = FALSE]

      meta <- meta |>
        dplyr::arrange(match(sample, colnames(profile_x)))

      meta$group <- factor(meta$group, levels = x)

      adonis <- vegan::adonis2(
        t(profile_x) ~ group,
        data = meta,
        permutations = permutations,
        distance = dist_method,
        ...
      )

      r2adj <- calcu_adjusted_r2(adonis)

      data.frame(
        comparison = paste0(x, collapse = "_vs_"),
        r2 = adonis$R2[1],
        r2adj = r2adj,
        pval = adonis$`Pr(>F)`[1],
        check.names = FALSE
      )
    }
  )

  if (isTRUE(add_plab)) {
    data <- data |>
      dplyr::mutate(
        plab = cut(
          pval,
          breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
          labels = c("***", "**", "*", "ns")
        )
      )
  }

  return(data)
}

#### plot_pairwise_adonis ####
# 绘制两两 PERMANOVA 结果
# data: calcu_pairwise_adonis() 输出结果
# group_level: 指定分组顺序

#' Plot Pairwise Adonis utility
#'
#'
#' 将 pairwise PERMANOVA 结果绘制为矩阵或热图式结果图。
#'
#' @param data An input data frame or compatible object.
#' @param group_level Optional order of group levels.
#' @return A ggplot-compatible plot object; computed data or fitted objects are retained as attributes when applicable.
#' @export
plot_pairwise_adonis <- function(data, group_level = NULL) {
  data <- data.frame(data, check.names = FALSE)

  if (!all(c("comparison", "r2", "pval") %in% colnames(data))) {
    stop("data should contain columns: comparison | r2 | pval")
  }

  if (is.null(group_level)) {
    group_level <- unique(c(
      stringr::str_split_i(data$comparison, "_vs_", 1),
      stringr::str_split_i(data$comparison, "_vs_", 2)
    ))
  }

  plot_df <- data.frame(
    x = stringr::str_split_i(data$comparison, "_vs_", 1),
    y = stringr::str_split_i(data$comparison, "_vs_", 2),
    r2 = data$r2,
    pval = data$pval
  ) |>
    dplyr::mutate(
      x = factor(x, levels = group_level),
      y = factor(y, levels = rev(group_level)),
      plab = dplyr::case_when(
        pval <= 0.001 ~ "p\u22640.001",
        pval < 0.01 ~ "p<0.01",
        pval < 0.05 ~ "p<0.05",
        TRUE ~ "p\u22650.05"
      ),
      plab = factor(
        plab,
        levels = c("p\u22640.001", "p<0.01", "p<0.05", "p\u22650.05")
      )
    )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x, y)) +
    ggplot2::geom_tile(
      fill = "transparent",
      color = "black",
      width = 1,
      height = 1,
      linewidth = .4
    ) +
    ggplot2::geom_point(
      ggplot2::aes(size = r2, fill = plab),
      shape = 21,
      color = "black",
      stroke = .4
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "p\u22640.001" = "#f46d43",
        "p<0.01"  = "#fee08b",
        "p<0.05"  = "#abdda4",
        "p\u22650.05"  = "#3288bd"
      ),
      breaks = c("p\u22640.001", "p<0.01", "p<0.05", "p\u22650.05")
    ) +
    ggplot2::scale_size_continuous(range = c(6, 12)) +
    ggplot2::labs(x = "", y = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 10, color = "black"),
      axis.title = ggplot2::element_text(size = 10, color = "black"),
      panel.border = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      aspect.ratio = 1
    ) +
    ggplot2::guides(
      size = ggplot2::guide_legend(title = stats::as.formula("'Adonis'~R^2"), order = 1),
      fill = ggplot2::guide_legend(title = "Significance", order = 2, override.aes = list(size = 4))
    )

  return(p)
}

#### calcu_betadisper ####
# 计算组内 beta dispersion，并检验不同组之间离散程度是否存在差异
# profile: 行为 feature，列为 sample 的丰度表
# group: 样本分组信息
# distance: 已计算的距离对象；NULL 时根据 profile 计算
# sample_col: group 中样本名称所在列
# group_col: group 中分组信息所在列
# group_level: 分组顺序
# dist_method: 距离方法，如 bray / jaccard / euclidean / unifrac
# permutations: permutation test 的置换次数
# type: 组中心类型，可选 median / centroid
# bias_adjust: 是否对小样本组的 dispersion 进行偏倚校正
# ...: 传入 calcu_distance()

#' Calcu Betadisper utility
#'
#'
#' 检验各组到中心或中位数中心的 multivariate dispersion。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param group A sample metadata table containing sample and group columns.
#' @param distance A precomputed distance object; when supplied, it takes precedence over `profile`.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param group_level Optional order of group levels.
#' @param dist_method Distance method; available values are validated with `match.arg()`.
#' @param permutations Number of permutations used by the significance test.
#' @param type Analysis or value type; supported values are shown in Usage.
#' @param bias_adjust Whether to apply the small-sample bias correction in `vegan::betadisper()`.
#' @param ... Additional arguments forwarded to `calcu_distance()`.
#' @return A list containing the betadisper object, group metadata, permutation/ANOVA/Tukey tests, and sample distances to group centers.
#' @export
calcu_betadisper <- function(
  profile = NULL, group, distance = NULL, sample_col = "sample",
  group_col = "group", group_level = NULL, dist_method = "bray",
  permutations = 999, type = c("median", "centroid"),
  bias_adjust = FALSE, ...
) {
  type <- match.arg(type)

  group <- data.frame(group, check.names = FALSE)
  if (!all(c(sample_col, group_col) %in% colnames(group))) {
    stop("group should contain columns: ", sample_col, " | ", group_col)
  }

  # Calculate distance if distance matrix is not provided
  if (is.null(distance)) {
    if (is.null(profile)) {
      stop("Need either profile or distance.")
    }
    distance <- calcu_distance(profile = profile, dist_method = dist_method, ...)
  }

  distance <- stats::as.dist(distance)
  sample_order <- labels(distance)

  # Match sample information with distance matrix
  meta <- group |>
    dplyr::select(
      sample = dplyr::all_of(sample_col),
      group = dplyr::all_of(group_col)
    ) |>
    dplyr::filter(sample %in% sample_order) |>
    dplyr::arrange(match(sample, sample_order))

  if (!identical(meta$sample, sample_order)) {
    stop(
      "Samples in distance and group do not match, ",
      "or duplicated samples exist in group."
    )
  }

  if (anyNA(meta$group)) {
    stop("Missing values are not allowed in group.")
  }

  # Set group levels
  if (is.null(group_level)) {
    group_level <- unique(meta$group)
  } else {
    if (!all(unique(meta$group) %in% group_level)) {
      stop("group_level should contain all groups in the data.")
    }
    group_level <- group_level[group_level %in% unique(meta$group)]
  }

  meta$group <- factor(meta$group, levels = group_level)

  if (length(group_level) < 2) {
    stop("At least two groups are required for betadisper.")
  }

  # Calculate beta dispersion
  betadisper <- vegan::betadisper(
    distance, meta$group,
    type = type, bias.adjust = bias_adjust
  )

  # Statistical tests
  permutest <- vegan::permutest(
    betadisper,
    permutations = permutations, pairwise = TRUE
  )
  anova_test <- stats::anova(betadisper)
  tukey_test <- stats::TukeyHSD(betadisper)

  # Distance of each sample to its group center
  dist_df <- data.frame(
    sample = names(betadisper$distances),
    value = as.numeric(betadisper$distances),
    check.names = FALSE
  ) |>
    dplyr::left_join(meta, by = "sample") |>
    dplyr::mutate(
      group = factor(group, levels = group_level)
    )

  # Output
  result <- list(
    object = betadisper,
    group = meta,
    permutest = permutest,
    anova = anova_test,
    tukey = tukey_test,
    dist_data = dist_df,
    type = type,
    bias_adjust = bias_adjust
  )

  attr(result$dist_data, "label") <- glue::glue(
    "The distance between each sample and the multivariate center of the group
    to which it belongs, calculated using the '{type}' method."
  )

  return(result)
}

#### plot_betadisper ####
# 绘制 calcu_betadisper 计算得到的组内 beta dispersion
# result: calcu_betadisper 返回的结果
# group_color: 各分组对应的颜色
# title: 图标题
# subtitle: 图副标题；NULL 时自动显示 permutation test 和 ANOVA 的 p 值
# aspect_ratio: 图形纵横比
# show_grid: 是否显示主网格线
# x_text_angle: x 轴文字旋转角度

#' Plot Betadisper utility
#'
#'
#' 绘制组内离散度及组间比较结果。
#'
#' @param result Result object or table to summarize.
#' @param group_color Optional colors aligned to `group_level`.
#' @param title Optional plot or result title.
#' @param subtitle Optional plot subtitle.
#' @param aspect_ratio Panel aspect ratio passed to `ggplot2::theme()`.
#' @param show_grid Whether to draw panel grid lines.
#' @param x_text_angle Rotation angle, in degrees, for x-axis text.
#' @return A ggplot-compatible plot object; computed data or fitted objects are retained as attributes when applicable.
#' @export
plot_betadisper <- function(
  result, group_color = NULL, title = NULL, subtitle = NULL,
  aspect_ratio = NULL, show_grid = FALSE, x_text_angle = 0
) {
  dist_df <- result$dist_data
  group_level <- levels(dist_df$group)

  # Set group colors
  if (is.null(group_color)) {
    palette <- c(
      "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
      "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3"
    )
    group_color <- structure(
      rep(palette, length.out = length(group_level)),
      names = group_level
    )
  } else {
    if (!is.null(names(group_color)) && all(group_level %in% names(group_color))) {
      group_color <- group_color[group_level]
    } else {
      group_color <- rep(group_color, length.out = length(group_level))
      names(group_color) <- group_level
    }
  }

  # Add statistical results to subtitle
  if (is.null(subtitle)) {
    p_perm <- result$permutest$tab[1, "Pr(>F)"]
    p_anova <- result$anova[1, "Pr(>F)"]

    p_perm <- ifelse(
      p_perm < .001, "italic(p) < 0.001",
      glue::glue('italic(p) == {sprintf("%.3f", p_perm)}')
    )

    p_anova <- ifelse(
      p_anova < .001, "italic(p) < 0.001",
      glue::glue('italic(p) == {sprintf("%.3f", p_anova)}')
    )

    subtitle <- glue::glue("'Permutation test:'~{p_perm}*','~~'ANOVA:'~{p_anova}")
    subtitle <- parse(text = subtitle)
  }

  # Y-axis label according to center type
  ylab <- if (identical(result$type, "centroid")) {
    "Distance to centroid"
  } else {
    "Distance to spatial median"
  }

  plt <- ggpubr::ggboxplot(
    dist_df,
    x = "group", y = "value", fill = "group", legend = "none",
    palette = group_color, xlab = "", ylab = ylab, outlier.shape = NA,
    x.text.angle = x_text_angle, title = title, subtitle = subtitle
  ) +
    ggplot2::geom_jitter(
      size = 2.8, width = .25, fill = "white", shape = 21, color = "black",
      show.legend = FALSE
    ) +
    ggplot2::theme(
      aspect.ratio = aspect_ratio,
      axis.ticks.length = grid::unit(2, "mm"),
      plot.title = ggplot2::element_text(hjust = .5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = .5)
    )

  if (isTRUE(show_grid)) {
    plt <- plt +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(linewidth = .5, color = "grey90"),
        panel.grid.minor = ggplot2::element_blank()
      )
  }

  return(plt)
}

#### calcu_adonis_r2 ####
# 根据距离矩阵和分组手动计算 PERMANOVA R2
# R² 的计算公式
# R² = SS_between / SS_total
# SS_between：组间平方和（Between-group Sum of Squares）
# SS_total：总平方和（Total Sum of Squares）
# SS_between = SS_total - SS_within
# SS_within：累加各组中的平方和

#' Calcu Adonis R2 utility
#'
#'
#' 根据距离对象和分组标签计算 PERMANOVA R-squared。
#'
#' @param dist A distance object used by the analysis.
#' @param group_labels Labels corresponding to group labels.
#' @return A numeric PERMANOVA R-squared value.
#' @export
calcu_adonis_r2 <- function(dist, group_labels) {
  dist <- as.matrix(dist)

  if (nrow(dist) != length(group_labels)) {
    stop("length of group_labels should be equal to sample number in dist.")
  }

  group_labels <- as.character(group_labels)

  SS_total <- sum(dist^2, na.rm = TRUE) / (2 * nrow(dist))

  SS_within <- purrr::map_vec(
    unique(group_labels), \(x) {
      sample_index <- which(group_labels == x)
      sub_dist <- dist[sample_index, sample_index, drop = FALSE]

      sum(sub_dist^2, na.rm = TRUE) / (2 * length(sample_index))
    }
  ) |>
    sum(na.rm = TRUE)

  r2 <- (SS_total - SS_within) / SS_total

  return(r2)
}
