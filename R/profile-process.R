#### Jinxin Meng, 20240308, 20260718, v0.4.3 ####

# 20231101: add all_group parameter in profile_filter().
# 20231219: add function profile_replace().
# 20250227: add function profile_top_n/frac().
# 20250605: add function profile_aggregate().
# 20250729: add function CLR() and the CLR profile transformation.
# 20250823: add digits parameter in profile_trans* functions
# 20260419: rename profile_smp2grp() as profile_collapse(), add updated for improving performance
#           remove profile_smp2grp_1(), update profile_replace(), add Hellinger transformation.
# 20260501 v0.4.1: update and re-coding by GPT-5.5.
# 20260526 v0.4.2: add parameter 'remove_unknown' in profile_aggregate().
#                   update function with removing *_rename() parameters.
# 20260718 v0.4.3: add zero-variance feature filtering.


#### LOG ####
# LOG transformation method in MaAsLin2
# use_half_minimum = TRUE: 0 替换为最小正值的一半
# use_half_minimum = FALSE: 0 替换为 pseudocount

#' LOG2 utility
#'
#' `LOG2()` provides a reusable mengR workflow with input validation and standardized output.
#'
#' Chinese summary: 对一个数值向量进行 log2 转换，并处理零值和缺失值。
#'
#' @param x Primary vector or object supplied to the utility.
#' @param pseudocount Positive value used to replace or offset zeros before logarithmic operations.
#' @param use_half_minimum Whether non-positive values are replaced by half the smallest positive value.
#' @return A result object described in the Details section.
#' @export
LOG2 <- function(x, pseudocount = 1e-6, use_half_minimum = TRUE) {
  
  x <- as.numeric(x)
  
  if (isTRUE(use_half_minimum)) {
    min_pos <- min(x[x > 0], na.rm = TRUE)
    if (!is.finite(min_pos)) min_pos <- pseudocount * 2
    x[x <= 0 | is.na(x)] <- min_pos / 2
  } else {
    x[x <= 0 | is.na(x)] <- pseudocount
  }
  
  log2(x)
}

#' LOG10 utility
#'
#' `LOG10()` provides a reusable mengR workflow with input validation and standardized output.
#'
#' Chinese summary: 对一个数值向量进行 log10 转换，并处理零值和缺失值。
#'
#' @param x Primary vector or object supplied to the utility.
#' @param pseudocount Positive value used to replace or offset zeros before logarithmic operations.
#' @param use_half_minimum Whether non-positive values are replaced by half the smallest positive value.
#' @return A result object described in the Details section.
#' @export
LOG10 <- function(x, pseudocount = 1e-6, use_half_minimum = TRUE) {
  
  x <- as.numeric(x)
  
  if (isTRUE(use_half_minimum)) {
    min_pos <- min(x[x > 0], na.rm = TRUE)
    if (!is.finite(min_pos)) min_pos <- pseudocount * 2
    x[x <= 0 | is.na(x)] <- min_pos / 2
  } else {
    x[x <= 0 | is.na(x)] <- pseudocount
  }
  
  log10(x)
}

#' CLR utility
#'
#' `CLR()` provides a reusable mengR workflow with input validation and standardized output.
#'
#' Chinese summary: 对组成型数值向量执行 centered log-ratio 转换。
#'
#' @param x Primary vector or object supplied to the utility.
#' @param pseudocount Positive value used to replace or offset zeros before logarithmic operations.
#' @param use_half_minimum Whether non-positive values are replaced by half the smallest positive value.
#' @return A result object described in the Details section.
#' @export
CLR <- function(x, pseudocount = 1e-6, use_half_minimum = FALSE) {
  
  x <- as.numeric(x)
  
  if (isTRUE(use_half_minimum)) {
    min_pos <- min(x[x > 0], na.rm = TRUE)
    if (!is.finite(min_pos)) min_pos <- pseudocount * 2
    x[x <= 0 | is.na(x)] <- min_pos / 2
  } else {
    x[x <= 0 | is.na(x)] <- pseudocount
  }
  
  gm <- exp(mean(log(x), na.rm = TRUE))  # 计算每行的几何均值
  log(x) - log(gm)                       # 执行CLR转换：log(x) - log(几何均值)
}

#### profile transformations ####
# profile: 行为 feature，列为 sample 的丰度表

#' Profile Trans Log2 utility
#'
#' `profile_trans_log2()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 对 feature × sample profile 按 feature 执行 log2 转换。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param pseudocount Positive value used to replace or offset zeros before logarithmic operations.
#' @param use_half_minimum Whether non-positive values are replaced by half the smallest positive value.
#' @param digits Optional number of decimal digits retained.
#' @return A result object described in the Details section.
#' @export
profile_trans_log2 <- function(profile, pseudocount = 1e-6,
                              use_half_minimum = FALSE, digits = NULL) {
  
  profile <- data.frame(profile, check.names = FALSE)
  profile <- apply(profile, 1, \(x) {
    LOG2(x, pseudocount = pseudocount, use_half_minimum = use_half_minimum)
  }) |>
    t() |>
    data.frame(check.names = FALSE)
  
  if (!is.null(digits)) {
    profile <- round(profile, digits = digits)
  }
  
  return(profile)
}

#' Profile Trans Log10 utility
#'
#' `profile_trans_log10()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 对 feature × sample profile 按 feature 执行 log10 转换。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param pseudocount Positive value used to replace or offset zeros before logarithmic operations.
#' @param use_half_minimum Whether non-positive values are replaced by half the smallest positive value.
#' @param digits Optional number of decimal digits retained.
#' @return A result object described in the Details section.
#' @export
profile_trans_log10 <- function(profile, pseudocount = 1e-6,
                               use_half_minimum = FALSE, digits = NULL) {
  
  profile <- data.frame(profile, check.names = FALSE)
  profile <- apply(profile, 1, \(x) {
    LOG10(x, pseudocount = pseudocount, use_half_minimum = use_half_minimum)
  }) |>
    t() |>
    data.frame(check.names = FALSE)
  
  if (!is.null(digits)) {
    profile <- round(profile, digits = digits)
  }
  
  return(profile)
}

#' Profile Trans Clr utility
#'
#' `profile_trans_clr()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 对每个样本执行 CLR 转换，保持 feature × sample 方向。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param pseudocount Positive value used to replace or offset zeros before logarithmic operations.
#' @param use_half_minimum Whether non-positive values are replaced by half the smallest positive value.
#' @param digits Optional number of decimal digits retained.
#' @return A result object described in the Details section.
#' @export
profile_trans_clr <- function(profile, pseudocount = 1e-6,
                             use_half_minimum = FALSE, digits = NULL) {
  
  profile <- data.frame(profile, check.names = FALSE)
  profile <- apply(profile, 2, \(x) {
    CLR(x, pseudocount = pseudocount, use_half_minimum = use_half_minimum)
  }) |>
    data.frame(check.names = FALSE)
  
  if (!is.null(digits)) {
    profile <- round(profile, digits = digits)
  }
  
  return(profile)
}

#' Profile Trans Sqrt utility
#'
#' `profile_trans_sqrt()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 对 profile 执行平方根转换。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param digits Optional number of decimal digits retained.
#' @return A result object described in the Details section.
#' @export
profile_trans_sqrt <- function(profile, digits = NULL) {
  
  profile <- data.frame(profile, check.names = FALSE)
  profile <- sqrt(profile)
  
  if (!is.null(digits)) {
    profile <- round(profile, digits = digits)
  }
  
  return(profile)
}

#### relative abundance ####
# 将 profile 转换为百分比或相对丰度
# profile: 行为 feature，列为 sample 的丰度表
# base: 转换基数；base = 100 为百分比，base = 1 为相对丰度
# digits: 保留小数位；NULL 表示不四舍五入
# remove_empty: 是否删除全 0 行和全 0 列

#' Profile Trans Ra utility
#'
#' `profile_trans_ra()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 将每个样本转换为相对丰度或百分比，并处理空样本。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param base Numeric setting for `base`.
#' @param digits Optional number of decimal digits retained.
#' @param remove_empty Logical control for `remove_empty`.
#' @return A result object described in the Details section.
#' @export
profile_trans_ra <- function(
    profile, base = 100, digits = 8, remove_empty = FALSE
  ) {
  
  profile <- data.frame(profile, check.names = FALSE)
  profile <- as.matrix(profile)
  suppressWarnings(storage.mode(profile) <- 'numeric')
  
  col_sum <- colSums(profile, na.rm = TRUE)
  col_sum[col_sum == 0] <- NA_real_
  
  profile <- sweep(profile, 2, col_sum, '/') * base
  profile[!is.finite(profile)] <- 0
  
  if (!is.null(digits)) {
    profile <- round(profile, digits = digits)
  }
  
  profile <- data.frame(profile, check.names = FALSE)
  
  if (isTRUE(remove_empty)) {
    profile <- profile[
      rowSums(profile, na.rm = TRUE) != 0,
      colSums(profile, na.rm = TRUE) != 0,
      drop = FALSE
    ]
  }
  
  return(profile)
}

#### Hellinger transformation ####
# Hellinger 转换，返回方向仍然是 feature × sample

#' Profile Trans Hellinger utility
#'
#' `profile_trans_hellinger()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 对 profile 执行 Hellinger 转换，适合部分生态距离分析。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param digits Optional number of decimal digits retained.
#' @return A result object described in the Details section.
#' @export
profile_trans_hellinger <- function(profile, digits = NULL) {
  
  profile <- data.frame(profile, check.names = FALSE)
  profile <- vegan::decostand(t(profile), method = 'hellinger') |>
    t() |>
    data.frame(check.names = FALSE)
  
  if (!is.null(digits)) {
    profile <- round(profile, digits = digits)
  }
  
  return(profile)
}

#### profile_collapse ####
# 按样本分组对 profile 进行合并
# profile: 行为 feature，列为 sample 的丰度表
# group: 样本分组表
# sample_col: group 中样本列名
# group_col: group 中分组列名
# method: mean / median / sum 或自定义函数
# na_fill: 是否填充 profile 中 NA；NULL 表示不填充

#' Profile Collapse utility
#'
#' `profile_collapse()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 按样本分组对 profile 求 mean、median、sum 或自定义统计量。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param group A sample metadata table containing sample and group columns.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param method Analysis or summary method; supported values are shown in the usage.
#' @param na_fill Value used to replace missing observations before analysis.
#' @param group_level Optional order of group levels.
#' @return A result object described in the Details section.
#' @export
profile_collapse <- function(
    profile, group, sample_col = "sample", group_col = "group",
    method = c('mean', 'median', 'sum'), na_fill = 0, group_level = NULL
  ) {
  if (is.function(method)) {
    stat_fun <- method
  } else {
    method <- match.arg(method)
    stat_fun <- match.fun(method)
  }

  profile <- data.frame(profile, check.names = FALSE)
  group <- data.frame(group, check.names = FALSE)
  
  if (!all(c(sample_col, group_col) %in% colnames(group))) {
    stop('group should contain columns: ', sample_col, ' | ', group_col)
  }
  
  if (anyDuplicated(group[[sample_col]])) {
    stop('Duplicated sample names found in group table.')
  }
  
  sample_use <- intersect(colnames(profile), group[[sample_col]])
  
  if (length(sample_use) == 0) {
    stop('No matched samples between profile and group table.')
  }
  
  profile <- profile[, sample_use, drop = FALSE]
  
  if (!is.null(na_fill)) {
    profile[is.na(profile)] <- na_fill
  }
  
  meta <- group |>
    dplyr::select(
      sample = dplyr::all_of(sample_col),
      group = dplyr::all_of(group_col)
    )
  
  profile_long_df <- data.frame(t(profile), check.names = FALSE) |>
    tibble::rownames_to_column('sample') |>
    dplyr::left_join(meta, by = 'sample') |>
    dplyr::filter(!is.na(group))
  
  if (is.null(group_level)) {
    group_level <- unique(profile_long_df$group)
  }

  profile_long_df <- profile_long_df |>
    dplyr::mutate(group = factor(group, levels = group_level))
  
  profile <- profile_long_df |>
    dplyr::select(-sample) |>
    dplyr::group_by(group) |>
    dplyr::summarise(
      dplyr::across(dplyr::where(is.numeric), \(x) stat_fun(x, na.rm = TRUE)),
      .groups = 'drop'
    ) |>
    tibble::column_to_rownames('group') |>
    t() |>
    data.frame(check.names = FALSE)
  
  return(profile)
}

#### profile_filter ####
# 根据 prevalence 或出现样本数过滤 profile
# by_group = FALSE: 在所有样本中筛选
# by_group = TRUE: 在每个分组内筛选
# all_group = TRUE: 所有组都必须通过
# n_group: 至少几个组通过

#' Profile Filter utility
#'
#' `profile_filter()` provides a reusable mengR workflow with input validation and standardized
#'   output.
#'
#' Chinese summary: 按总体或组内 prevalence、出现样本数和最低丰度筛选 feature。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param group A sample metadata table containing sample and group columns.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param by_group Whether calculations are performed separately within each group.
#' @param all_group Whether every group, rather than at least `n_group`, must pass filtering.
#' @param n_group Minimum number of groups that must pass the filtering criterion.
#' @param min_prevalence Minimum fraction of samples in which a feature must be present.
#' @param min_n Minimum number of samples in which a feature must be present.
#' @param min_abundance Minimum abundance above which a feature is considered present.
#' @return A result object described in the Details section.
#' @export
profile_filter <- function(
    profile, group = NULL, sample_col = 'sample', group_col = 'group',
    by_group = FALSE, all_group = FALSE, n_group = 1,
    min_prevalence = 0.1, min_n = NULL, min_abundance = 0
  ) {
  
  profile <- data.frame(profile, check.names = FALSE)
  profile <- as.matrix(profile)
  suppressWarnings(storage.mode(profile) <- 'numeric')
  
  present <- !is.na(profile) & profile > min_abundance
  
  ## 不分组过滤
  if (isFALSE(by_group)) {
    
    present_n <- rowSums(present, na.rm = TRUE)
    
    if (!is.null(min_n)) {
      keep <- present_n >= min_n
    } else {
      keep <- present_n / ncol(profile) >= min_prevalence
    }
    
    profile <- profile[keep, , drop = FALSE]
    return(data.frame(profile, check.names = FALSE))
  }
  
  ## 分组过滤
  if (is.null(group)) {
    stop('if by_group = TRUE, group should be provided.')
  }
  
  group <- data.frame(group, check.names = FALSE)
  
  if (!all(c(sample_col, group_col) %in% colnames(group))) {
    stop('group should contain columns: ', sample_col, ' | ', group_col)
  }
  
  if (anyDuplicated(group[[sample_col]])) {
    stop('Duplicated sample names found in group table.')
  }
  
  group <- group[match(colnames(profile), group[[sample_col]]), ]
  
  if (any(is.na(group[[sample_col]])) || any(is.na(group[[group_col]]))) {
    miss <- colnames(profile)[is.na(group[[sample_col]]) | is.na(group[[group_col]])]
    stop('Some profile samples are not found in group table: ',
         paste(utils::head(miss, 10), collapse = ', '))
  }
  
  group_vec <- as.character(group[[group_col]])
  group_level <- unique(group_vec)
  
  count_mat <- sapply(group_level, \(x) {
    rowSums(present[, group_vec == x, drop = FALSE], na.rm = TRUE)
  })
  
  if (is.null(dim(count_mat))) {
    count_mat <- matrix(
      count_mat,
      nrow = nrow(profile),
      ncol = length(group_level),
      dimnames = list(rownames(profile), group_level)
    )
  }
  
  if (!is.null(min_n)) {
    pass_mat <- count_mat >= min_n
  } else {
    group_size <- as.numeric(table(factor(group_vec, levels = group_level)))
    pass_mat <- sweep(count_mat, 2, group_size, '/') >= min_prevalence
  }
  
  if (isTRUE(all_group)) {
    keep <- rowSums(pass_mat, na.rm = TRUE) == ncol(pass_mat)
  } else {
    keep <- rowSums(pass_mat, na.rm = TRUE) >= n_group
  }
  
  profile <- profile[keep, , drop = FALSE]
  
  return(data.frame(profile, check.names = FALSE))
}

#### profile_top_n ####
# 选择丰度最高的前 n 个 feature
# profile: 行为 feature，列为 sample 的丰度表
# n: 保留 feature 数
# out_other: 是否将未进入 top n 的 feature 合并为 Other
# other_name: 合并后的名称
# sort_method: 排序方法，可选 mean / sum / median，也可以输入自定义函数

#' Profile Top n utility
#'
#' `profile_top_n()` provides a reusable mengR workflow with input validation and standardized
#'   output.
#'
#' Chinese summary: 保留总体丰度最高的 n 个 feature，可把其余 feature 合并为 Other。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param n Requested number of values, features, or results.
#' @param out_other Whether discarded features are combined into an `Other` category.
#' @param other_name Label assigned to features combined into the residual category.
#' @param sort_method Summary statistic or function used to rank features.
#' @return A result object described in the Details section.
#' @export
profile_top_n <- function(profile, n = 12, out_other = FALSE,
                          other_name = 'Other',
                          sort_method = c('mean', 'sum', 'median')) {
  
  profile <- data.frame(profile, check.names = FALSE)
  profile <- as.matrix(profile)
  suppressWarnings(storage.mode(profile) <- 'numeric')
  
  if (n <= 0) {
    stop('n should be > 0.')
  }
  
  if (nrow(profile) == 0) {
    return(data.frame())
  }
  
  n <- min(n, nrow(profile))
  
  if (is.function(sort_method)) {
    stat <- apply(profile, 1, \(x) sort_method(x, na.rm = TRUE))
  } else {
    sort_method <- match.arg(sort_method)
  }

  if (!is.function(sort_method) && sort_method == 'mean') {
    stat <- rowMeans(profile, na.rm = TRUE)
  } else if (!is.function(sort_method) && sort_method == 'sum') {
    stat <- rowSums(profile, na.rm = TRUE)
  } else if (!is.function(sort_method) && sort_method == 'median') {
    stat <- apply(profile, 1, stats::median, na.rm = TRUE)
  }
  
  feats <- names(sort(stat, decreasing = TRUE))[seq_len(n)]
  
  if (isFALSE(out_other)) {
    
    profile <- profile[feats, , drop = FALSE]
    
  } else {
    
    group <- rownames(profile)
    group[!group %in% feats] <- other_name
    
    profile <- rowsum(profile, group = group, reorder = FALSE, na.rm = TRUE)
    
    row_order <- c(feats[feats %in% rownames(profile)], other_name)
    row_order <- row_order[row_order %in% rownames(profile)]
    
    profile <- profile[row_order, , drop = FALSE]
  }
  
  return(data.frame(profile, check.names = FALSE))
}

#### profile_top_frac ####
# 按比例选择丰度最高的 feature
# frac: 保留比例，范围为 (0, 1]
# out_other: 是否将未进入 top fraction 的 feature 合并为 Other

#' Profile Top Frac utility
#'
#' `profile_top_frac()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 按比例保留丰度最高的 feature，可把其余部分合并为 Other。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param frac Fraction of the highest-ranking features retained.
#' @param out_other Whether discarded features are combined into an `Other` category.
#' @param other_name Label assigned to features combined into the residual category.
#' @param sort_method Summary statistic or function used to rank features.
#' @return A result object described in the Details section.
#' @export
profile_top_frac <- function(profile, frac = 0.1, out_other = FALSE,
                             other_name = 'Other',
                             sort_method = c('mean', 'sum', 'median')) {
  
  profile <- data.frame(profile, check.names = FALSE)
  
  if (!is.numeric(frac) || frac <= 0 || frac > 1) {
    stop('frac should be in (0, 1].')
  }
  
  if (nrow(profile) == 0) {
    return(data.frame())
  }
  
  n <- max(1, floor(nrow(profile) * frac))
  
  profile <- profile_top_n(
    profile = profile,
    n = n,
    out_other = out_other,
    other_name = other_name,
    sort_method = sort_method
  )
  
  return(profile)
}

#### profile_replace ####
# replace values <= min_value with fill_value
# trans_ra = TRUE: 先转为相对丰度/百分比

#' Profile Replace utility
#'
#' `profile_replace()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 按阈值替换低丰度值，并可在替换前转换为相对丰度。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param min_value Numeric limit controlling min value.
#' @param fill_value Value assigned to observations that fail the replacement threshold.
#' @param trans_ra Whether abundances are converted to relative abundance before analysis.
#' @param base Numeric setting for `base`.
#' @param remove_empty Logical control for `remove_empty`.
#' @return A result object described in the Details section.
#' @export
profile_replace <- function(profile, min_value = 1, fill_value = 0,
                            trans_ra = FALSE, base = 100,
                            remove_empty = TRUE) {
  
  profile <- data.frame(profile, check.names = FALSE)
  profile <- as.matrix(profile)
  suppressWarnings(storage.mode(profile) <- 'numeric')
  
  if (isTRUE(trans_ra)) {
    col_sum <- colSums(profile, na.rm = TRUE)
    col_sum[col_sum == 0] <- NA_real_
    profile <- sweep(profile, 2, col_sum, '/') * base
  }
  
  profile[!is.finite(profile)] <- fill_value
  profile[profile <= min_value] <- fill_value
  
  if (isTRUE(remove_empty)) {
    profile <- profile[
      rowSums(profile != fill_value, na.rm = TRUE) > 0,
      colSums(profile != fill_value, na.rm = TRUE) > 0,
      drop = FALSE
    ]
  }
  
  return(data.frame(profile, check.names = FALSE))
}

#### profile_adjacency ####
# 将 profile 转换为 feature × sample 的存在/缺失矩阵
# logical = TRUE: 输出 TRUE/FALSE
# logical = FALSE: 输出 1/0
# min_abundance: 大于该阈值视为存在

#' Profile Adjacency utility
#'
#' `profile_adjacency()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 将丰度表转换为 feature × sample 的存在/缺失矩阵。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param logical Whether the adjacency matrix contains logical presence/absence values.
#' @param min_abundance Minimum abundance above which a feature is considered present.
#' @return A result object described in the Details section.
#' @export
profile_adjacency <- function(profile, logical = FALSE, min_abundance = 0) {
  
  profile <- data.frame(profile, check.names = FALSE)
  profile <- as.matrix(profile)
  suppressWarnings(storage.mode(profile) <- 'numeric')
  
  profile <- !is.na(profile) & profile > min_abundance
  
  if (isFALSE(logical)) {
    storage.mode(profile) <- 'integer'
  }
  
  profile <- data.frame(profile, check.names = FALSE)
  
  return(profile)
}

#### profile_prevalence ####
# 计算 feature 在全部样本或各分组中的流行率
# profile: 行为 feature，列为 sample 的丰度表
# group: 样本分组表
# by_group: 是否按分组计算
# min_abundance: 大于该值视为存在
# count: TRUE 输出出现样本数；FALSE 输出流行率
# base: count = FALSE 时的基数；100 表示百分比，1 表示比例

#' Profile Prevalence utility
#'
#' `profile_prevalence()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 计算 feature 在总体或各组中的出现样本数或 prevalence。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param group A sample metadata table containing sample and group columns.
#' @param by_group Whether calculations are performed separately within each group.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param min_abundance Minimum abundance above which a feature is considered present.
#' @param count Whether prevalence is returned as sample counts instead of proportions.
#' @param base Numeric setting for `base`.
#' @return A result object described in the Details section.
#' @export
profile_prevalence <- function(profile, group = NULL, by_group = TRUE,
                               sample_col = 'sample', group_col = 'group',
                               min_abundance = 0, count = FALSE, base = 100) {
  
  profile <- data.frame(profile, check.names = FALSE)
  profile <- as.matrix(profile)
  suppressWarnings(storage.mode(profile) <- 'numeric')
  
  present <- !is.na(profile) & profile > min_abundance
  
  ## 全部样本
  if (isFALSE(by_group)) {
    
    value <- rowSums(present, na.rm = TRUE)
    
    if (isFALSE(count)) {
      value <- value / ncol(profile) * base
    }
    
    out <- data.frame(
      name = rownames(profile),
      prevalence = value,
      check.names = FALSE
    )
    
    return(out)
  }
  
  ## 分组样本
  if (is.null(group)) {
    stop('if by_group = TRUE, group should be provided.')
  }
  
  group <- data.frame(group, check.names = FALSE)
  
  if (!all(c(sample_col, group_col) %in% colnames(group))) {
    stop('group should contain columns: ', sample_col, ' | ', group_col)
  }
  
  sample_use <- intersect(colnames(profile), group[[sample_col]])
  
  if (length(sample_use) == 0) {
    stop('No matched samples between profile and group table.')
  }
  
  present <- present[, sample_use, drop = FALSE]
  
  group_vec <- group[[group_col]][match(sample_use, group[[sample_col]])]
  group_level <- unique(group_vec)
  
  out <- sapply(group_level, \(x) {
    rowSums(present[, group_vec == x, drop = FALSE], na.rm = TRUE)
  })
  
  if (is.null(dim(out))) {
    out <- matrix(
      out,
      nrow = nrow(profile),
      ncol = length(group_level),
      dimnames = list(rownames(profile), group_level)
    )
  }
  
  if (isFALSE(count)) {
    group_size <- as.numeric(table(factor(group_vec, levels = group_level)))
    out <- sweep(out, 2, group_size, '/') * base
  }
  
  out <- data.frame(out, check.names = FALSE)
  
  return(out)
}

#### profile_statistics ####
# 计算 profile 的基础统计量
# by_group = TRUE: 每组分别计算 mean/sd/median/prevalence
# by_group = FALSE: 所有样本整体计算

#' Profile Statistics utility
#'
#' `profile_statistics()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 计算 feature 的 mean、SD、median 和 prevalence，可按组汇总。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param group A sample metadata table containing sample and group columns.
#' @param by_group Whether calculations are performed separately within each group.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param min_abundance Minimum abundance above which a feature is considered present.
#' @param base Numeric setting for `base`.
#' @return A result object described in the Details section.
#' @export
profile_statistics <- function(profile, group = NULL, by_group = TRUE,
                               sample_col = 'sample', group_col = 'group',
                               min_abundance = 0, base = 100) {
  
  profile <- data.frame(profile, check.names = FALSE)
  profile <- as.matrix(profile)
  suppressWarnings(storage.mode(profile) <- 'numeric')
  
  ## 不分组统计
  if (isFALSE(by_group)) {
    
    out <- data.frame(
      name = rownames(profile),
      mean = rowMeans(profile, na.rm = TRUE),
      sd = apply(profile, 1, stats::sd, na.rm = TRUE),
      median = apply(profile, 1, stats::median, na.rm = TRUE),
      prevalence = rowSums(profile > min_abundance, na.rm = TRUE) / ncol(profile) * base,
      check.names = FALSE
    )
    
    return(out)
  }
  
  ## 分组统计
  if (is.null(group)) {
    stop('if by_group = TRUE, group should be provided.')
  }
  
  group <- data.frame(group, check.names = FALSE)
  
  if (!all(c(sample_col, group_col) %in% colnames(group))) {
    stop('group should contain columns: ', sample_col, ' | ', group_col)
  }
  
  sample_use <- intersect(colnames(profile), group[[sample_col]])
  
  if (length(sample_use) == 0) {
    stop('No matched samples between profile and group table.')
  }
  
  profile <- profile[, sample_use, drop = FALSE]
  group_vec <- group[[group_col]][match(sample_use, group[[sample_col]])]
  group_level <- unique(group_vec)
  
  out <- purrr::map_dfc(group_level, \(x) {
    
    group_profile_mat <- profile[, group_vec == x, drop = FALSE]
    
    data.frame(
      mean = rowMeans(group_profile_mat, na.rm = TRUE),
      sd = apply(group_profile_mat, 1, stats::sd, na.rm = TRUE),
      median = apply(group_profile_mat, 1, stats::median, na.rm = TRUE),
      prevalence = rowSums(
        group_profile_mat > min_abundance, na.rm = TRUE
      ) / ncol(group_profile_mat) * base,
      check.names = FALSE
    ) |>
      dplyr::rename_with(~ paste0(x, '_', .x))
  })
  
  out <- out |>
    tibble::add_column(name = rownames(profile), .before = 1) |>
    data.frame(check.names = FALSE)
  
  return(out)
}

#### profile_aggregate ####
# 按行特征的注释信息进行 profile 汇总
# profile: 行为 feature，列为 sample 的丰度表
# metadata: feature 注释表；默认第 1 列为 feature，第 2 列为 group
# feature_col: metadata 中用于匹配 profile 行名的列，默认第 1 列
# group_col: metadata 中用于汇总分类的列，默认第 2 列，可以是多列
# method: 汇总方法，可选 sum / mean / median，也可以输入自定义函数
# unknown: 未注释或空分组的填充值
# remove_unknown: 是否删除 unknown 分组的 feature
#' Profile Aggregate utility
#'
#' `profile_aggregate()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 根据 feature metadata 聚合 profile，并处理未知分类。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param metadata A metadata or annotation data frame.
#' @param feature_col Name of the feature-identifier column.
#' @param group_col Name of the grouping column.
#' @param method Analysis or summary method; supported values are shown in the usage.
#' @param unknown Label assigned to unknown or unclassified taxonomy entries.
#' @param remove_unknown Logical control for `remove_unknown`.
#' @param unknown_pattern Regular expression identifying unknown or unclassified annotations.
#' @param sep Field separator used when reading or writing a text file.
#' @return A result object described in the Details section.
#' @export
profile_aggregate <- function(
    profile, metadata, feature_col = NULL, group_col = NULL, 
    method = c('sum', 'mean', 'median'), unknown = 'unknown',
    remove_unknown = FALSE, unknown_pattern = "unknown|unclassified|unassigned",
    sep = '|'
  ) {
  
  profile <- data.frame(profile, check.names = FALSE) |>
    tibble::rownames_to_column('.feature')
  
  metadata <- data.frame(metadata, check.names = FALSE)
  
  if (ncol(metadata) < 2) {
    stop('metadata should contain at least two columns: feature and group.')
  }
  
  ## 默认 metadata 第 1 列是 feature，第 2 列是 group
  if (is.null(feature_col)) feature_col <- colnames(metadata)[1]
  if (is.null(group_col)) group_col <- colnames(metadata)[2]
  
  ## 支持用列号指定
  if (is.numeric(feature_col)) feature_col <- colnames(metadata)[feature_col]
  if (is.numeric(group_col)) group_col <- colnames(metadata)[group_col]
  
  if (!feature_col %in% colnames(metadata)) {
    stop('feature_col not found in metadata: ', feature_col)
  }
  
  if (!all(group_col %in% colnames(metadata))) {
    stop('group_col not found in metadata: ',
         paste(setdiff(group_col, colnames(metadata)), collapse = ', '))
  }
  
  sample_cols <- setdiff(colnames(profile), '.feature')
  
  if (is.function(method)) {
    fun <- \(x) method(x, na.rm = TRUE)
  } else {
    method <- match.arg(method)
    fun <- switch(
      method,
      sum = \(x) sum(x, na.rm = TRUE),
      mean = \(x) mean(x, na.rm = TRUE),
      median = \(x) stats::median(x, na.rm = TRUE)
    )
  }
  
  ## 匹配注释
  data <- profile |>
    dplyr::left_join(
      dplyr::select(metadata, dplyr::all_of(c(feature_col, group_col))),
      by = c('.feature' = feature_col)
    ) |>
    dplyr::mutate(
      dplyr::across(dplyr::all_of(sample_cols), as.numeric),
      dplyr::across(
        dplyr::all_of(group_col),
        \(x) {
          x <- as.character(x)
          x[is.na(x) | x == ''] <- unknown
          x
        }
      )
    )
  
  ## 删除 unknown / unclassified 等未明确注释的 feature
  ## 多列 group_col 时，只要任意一列匹配 unknown_pattern 就删除
  if (isTRUE(remove_unknown)) {
    data <- data |>
      dplyr::filter(
        dplyr::if_all(
          dplyr::all_of(group_col),
          \(x) !grepl(unknown_pattern, x, ignore.case = TRUE)
        )
      )
  }
  
  if (nrow(data) == 0) {
    warning('No features remained after filtering metadata.')
    return(data.frame())
  }
  
  out <- data |>
    dplyr::select(dplyr::all_of(c(group_col, sample_cols))) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_col))) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(sample_cols), fun),
      .groups = 'drop'
    )
  
  if (length(group_col) == 1) {
    group_name <- out[[group_col]]
  } else {
    group_name <- apply(out[, group_col, drop = FALSE], 1, paste, collapse = sep)
  }
  
  out <- out |>
    dplyr::select(-dplyr::all_of(group_col)) |>
    data.frame(check.names = FALSE)
  
  rownames(out) <- group_name
  
  out <- out[
    rowSums(out, na.rm = TRUE) != 0,
    colSums(out, na.rm = TRUE) != 0,
    drop = FALSE
  ]
  
  return(out)
}

#### remove zero-variance features ####
# 删除丰度表中在所有样本间无变异的 feature
# profile: 行为 feature、列为 sample 的丰度表或数值矩阵
# 判定方法: 按行计算标准差，仅保留标准差有限且大于 0 的 feature
# NA处理: 计算标准差时忽略 NA；全为 NA 或有效值不足的 feature 将被删除
# 返回值: 删除零方差及无法计算方差的 feature 后的丰度表
#' Profile Remove Zero Var utility
#'
#' `profile_remove_zero_var()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 删除跨样本零方差的 feature，避免降维或模型拟合失败。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @return A result object described in the Details section.
#' @export
profile_remove_zero_var <- function(profile) {
  row_sd <- apply(profile, 1, stats::sd, na.rm = TRUE)
  keep <- is.finite(row_sd) & row_sd > 0
  profile[keep, , drop = FALSE]
}
