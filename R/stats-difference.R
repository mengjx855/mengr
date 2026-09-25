#### Jin-Xin Meng, jinxmeng@zju.edu.cn, 20220101, 20260923 ####

# 20220601: 可选择'wilcox rank-sum','one-way anova','student's t test'三种方法做差异分析；
# 20230117: diff_test_profile函数对feature进行差异分析，输入的是标准otu表和group表
# 20231204: 修改diff_test_profile函数中的for循环，使用purrr::map_dfr，速度上稍微快一丢丢。
# 20231204: 修改diff_test函数中的rbind()，使用tibble::add_column()，速度上稍微快一丢丢。
# 20250109: 添加round参数，格式化pval， sig 参数，只输出显著的
# 20250404: 修改部分参数名称，sig <- filter_by_pval, xx_colnames <- xx_rename，使用 map_xx 代替 for循环，加速计算速度
# 20250417: 修改一些BUG
# 20250908: 合并difference_analysis, 改本脚本名称为calcu_difference.R
# 20250915: 更新部分函数，简化函数，使用formula代替分组文件
# 20251013: 更新difference_analysis，计算P与丰度估计使用的数据是否转换要分开；
# 20251105: 修改一些BUG
# 20260523: add new function 'calcu_empirical_p()'
# 20260624: 更新difference_analysis
# 20260916: standardize script metadata, function sections, documentation, and naming style.
# 20260923: clarify metadata argument names and remove Chinese text from Roxygen documentation.



#### .add_plab ####
# 根据 p 值生成显著性标签
# format:
#   1: '***', '**', '*', 'ns'
#   2: '***', '**', '*', ''
#   3: '#', '+', '*', ''
#   4: '+', '**', '*', ''
#   5: 'p<0.001', 'p<0.01', 'p<0.05', ''
#   6: round(pvalue, 3)

.add_plab <- function(pvalue, format = 2) {
  if (missing(pvalue)) {
    cat("  format-1: '***', '**', '*', 'ns'\n")
    cat("  format-2: '***', '**', '*',   ''\n")
    cat("  format-3:   '#',  '+', '*',   ''\n")
    cat("  format-4:   '+', '**', '*',   ''\n")
    cat("  format-5: 'p<0.001', 'p<0.01', 'p<0.05', 'p>0.05'\n")
    cat("  format-6: round(pvalue, digits = 3)\n")
    return(invisible(NULL))
  }

  if (format == 6) {
    return(round(pvalue, 3))
  }

  labs <- switch(as.character(format),
    "1" = c("***", "**", "*", "ns"),
    "2" = c("***", "**", "*", ""),
    "3" = c("#", "+", "*", ""),
    "4" = c("+", "**", "*", ""),
    "5" = c("p<0.001", "p<0.01", "p<0.05", "p>0.05"),
    c("***", "**", "*", "")
  )

  label <- cut(
    pvalue,
    include.lowest = TRUE,
    breaks = c(0, 0.001, 0.01, 0.05, Inf),
    labels = labs
  ) |>
    as.character()

  label[is.na(label)] <- ""

  return(label)
}


## 将非正值和非有限值替换为按 feature 计算的伪计数
#### .replace_nonpositive ####

.replace_nonpositive <- function(x, pseudo_factor = 0.5) {
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  out <- x

  global_min <- suppressWarnings(
    min(out[is.finite(out) & out > 0], na.rm = TRUE)
  )
  if (!is.finite(global_min) || global_min <= 0) global_min <- 1e-12

  for (i in seq_len(nrow(out))) {
    value_vec <- out[i, ]
    positive_vec <- value_vec[is.finite(value_vec) & value_vec > 0]
    pseudo <- if (length(positive_vec)) {
      min(positive_vec, na.rm = TRUE) * pseudo_factor
    } else {
      global_min * pseudo_factor
    }
    if (!is.finite(pseudo) || pseudo <= 0) pseudo <- 1e-12

    value_vec[!is.finite(value_vec) | value_vec <= 0] <- pseudo
    out[i, ] <- value_vec
  }

  return(out)
}


#### calcu_empirical_p ####
# 计算经验 P 值和标准化效应量 SES
# obs: 观测值
# random: 随机 null 分布
# alternative:
#   auto: obs >= random mean 时右尾，否则左尾
#   greater: P(random >= obs)
#   less:  P(random <= obs)
#   two.sided: 双尾经验 P 值

#' Calcu Empirical p utility
#'
#'
#'
#' @param obs Observed statistic compared with the randomization distribution.
#' @param random Randomized statistics forming the empirical null distribution.
#' @param alternative Alternative hypothesis used to calculate the empirical P value.
#' @param simplify Whether to return the simplified tabular result instead of intermediate objects.
#' @return A numeric empirical P value, or a data frame of observed values and empirical P values when `simplify = FALSE`.
#' @export
calcu_empirical_p <- function(
  obs, random, alternative = c("auto", "greater", "less", "two.sided"),
  simplify = FALSE
) {
  alternative <- match.arg(alternative)

  random <- random[is.finite(random)]
  n <- length(random)

  if (n == 0 || !is.finite(obs)) {
    out <- list(
      observed = obs,
      random_mean = NA_real_,
      SES = NA_real_,
      pval = NA_real_,
      side = NA_character_
    )

    if (isTRUE(simplify)) {
      out <- data.frame(out, check.names = FALSE)
    }

    return(out)
  }

  random_mean <- mean(random, na.rm = TRUE)

  if (alternative == "auto") {
    alternative <- ifelse(obs >= random_mean, "greater", "less")
  }

  if (alternative == "greater") {
    p <- (1 + sum(random >= obs, na.rm = TRUE)) / (1 + n)
    side <- "greater"
  }

  if (alternative == "less") {
    p <- (1 + sum(random <= obs, na.rm = TRUE)) / (1 + n)
    side <- "less"
  }

  if (alternative == "two.sided") {
    p1 <- (1 + sum(random >= obs, na.rm = TRUE)) / (1 + n)
    p2 <- (1 + sum(random <= obs, na.rm = TRUE)) / (1 + n)
    p <- min(1, 2 * min(p1, p2))
    side <- "two.sided"
  }

  sd_random <- stats::sd(random, na.rm = TRUE)
  ses <- ifelse(sd_random > 0, (obs - random_mean) / sd_random, NA_real_)

  out <- list(
    observed = obs,
    random_mean = random_mean,
    SES = ses,
    pval = p,
    side = side
  )

  if (isTRUE(simplify)) {
    out <- data.frame(out, check.names = FALSE)
  }

  return(out)
}

#### calcu_diff ####
# 根据 formula 进行两两差异检验
# data: 数据框
# formula: value ~ group
# method: wilcox / anova / t
# var_equal: t 检验是否假设方差齐性
# add_plab: 是否添加显著性标签

#' Calcu Diff utility
#'
#'
#'
#' @param data An input data frame or compatible object.
#' @param formula Model formula defining the response and grouping variables.
#' @param method Analysis or summary method; supported values are shown in the usage.
#' @param var_equal Whether two-sample t tests assume equal group variances.
#' @param add_plab Whether to add formatted significance labels.
#' @param plab_fmt Format string or function used to display significance labels.
#' @param ... Additional arguments passed to the selected rstatix test.
#' @return A data frame of pairwise comparisons, effect summaries, P values, and optional significance labels.
#' @export
calcu_diff <- function(data, formula, method = c("wilcox", "anova", "t"),
                       var_equal = FALSE, add_plab = FALSE,
                       plab_fmt = 2, ...) {
  method <- match.arg(method)

  terms <- stats::terms(formula)
  response <- all.vars(terms)[1]
  group <- all.vars(terms)[2]

  if (!response %in% names(data)) {
    stop("variable ", response, " not in data")
  }

  if (!group %in% names(data)) {
    stop("variable ", group, " not in data")
  }

  data <- data |>
    dplyr::filter(!is.na(.data[[response]]), !is.na(.data[[group]]))

  if (is.factor(data[[group]])) {
    group_level <- levels(droplevels(data[[group]]))
  } else {
    group_level <- unique(as.character(data[[group]]))
  }

  if (length(group_level) < 2) {
    stop("at least two groups are required.")
  }

  comparison <- utils::combn(group_level, m = 2, simplify = FALSE)

  difference <- purrr::map_dfr(
    comparison, \(x) {
      pair_df <- dplyr::filter(data, .data[[group]] %in% x)

      test <- tryCatch(
        {
          if (method == "wilcox") {
            stats::wilcox.test(formula, pair_df, ...) |> suppressWarnings()
          } else if (method == "anova") {
            stats::oneway.test(formula, pair_df, ...)
          } else if (method == "t") {
            stats::t.test(
              formula, pair_df,
              var.equal = var_equal, ...
            ) |>
              suppressWarnings()
          }
        },
        error = function(e) NULL
      )

      method_name <- dplyr::case_when(
        method == "wilcox" ~ "Wilcoxon rank-sum test",
        method == "anova" ~ "One-way ANOVA test",
        method == "t" & isTRUE(var_equal) ~ "Student's t-test",
        method == "t" & isFALSE(var_equal) ~ "Welch t-test"
      )

      data.frame(
        comparison = paste0(x, collapse = "_vs_"),
        pval = ifelse(is.null(test), NA_real_, test$p.value),
        method = method_name
      )
    }
  )

  if (nrow(difference) >= 3) {
    difference <- difference |>
      dplyr::mutate(padj = stats::p.adjust(pval, method = "BH"), .after = "pval")
  }

  if (isTRUE(add_plab)) {
    if ("padj" %in% colnames(difference)) {
      difference <- difference |>
        dplyr::mutate(plab = .add_plab(padj, plab_fmt), .after = "padj")
    } else {
      difference <- difference |>
        dplyr::mutate(plab = .add_plab(pval, plab_fmt), .after = "pval")
    }
  }

  return(difference)
}

#### calcu_diff_profile ####
# 对 profile 中每个 feature 进行差异检验
# profile: 行为 feature，列为 sample
# group: 样本分组信息，至少包含 sample 和 group_by
# comparison: 指定比较组，例如 c('IBD', 'HC')；NULL 时自动两两比较

#' Calcu Diff Profile utility
#'
#'
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param sample_meta A sample metadata table containing sample and group columns.
#' @param group_by Metadata column or grouping definition used for aggregation.
#' @param comparison Two outcome levels ordered as case and control.
#' @param method Analysis or summary method; supported values are shown in the usage.
#' @param add_plab Whether to add formatted significance labels.
#' @param plab_fmt Format string or function used to display significance labels.
#' @param var_equal Whether two-sample t tests assume equal group variances.
#' @param progress Whether progress information is printed during repeated analyses.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param ... Additional arguments forwarded to `calcu_diff()`.
#' @return A feature-by-comparison data frame of test results, adjusted P values, and optional significance labels.
#' @export
calcu_diff_profile <- function(
  profile, sample_meta, group_by = NULL, comparison = NULL,
  method = c("wilcox", "anova", "t"), add_plab = FALSE,
  plab_fmt = 2, var_equal = FALSE, progress = TRUE,
  sample_col = "sample", group_col = "group", ...
) {
  group <- sample_meta
  method <- match.arg(method)

  ## Backward compatibility for the former group_by argument.
  if (!is.null(group_by)) group_col <- group_by

  aligned <- .align_profile_group(
    profile = profile,
    sample_meta = group,
    sample_col = sample_col,
    group_col = group_col
  )

  profile_df <- aligned$profile_df
  group_df <- aligned$group_df
  group_vec <- as.character(group_df[[group_col]])

  if (is.null(comparison)) {
    comparison <- utils::combn(unique(group_vec), m = 2, simplify = FALSE)
  } else if (is.vector(comparison) && !is.list(comparison)) {
    comparison <- list(comparison)
  }

  test_df <- data.frame(t(profile_df), check.names = FALSE) |>
    tibble::rownames_to_column("sample") |>
    dplyr::mutate(group = group_vec)

  feature_vec <- rownames(profile_df)

  result_df <- purrr::map_dfr(
    comparison, \(x) {
      purrr::map_dfr(
        feature_vec, \(feature) {
          test_df |>
            dplyr::select(
              group,
              value = dplyr::all_of(feature)
            ) |>
            dplyr::filter(group %in% x) |>
            calcu_diff(
              value ~ group,
              method = method,
              var_equal = var_equal,
              ...
            ) |>
            tibble::add_column(name = feature, .before = 1)
        },
        .progress = progress
      )
    }
  )


  result_df <- result_df |>
    dplyr::mutate(padj = stats::p.adjust(pval, method = "BH"), .after = "pval")

  if (isTRUE(add_plab)) {
    result_df <- result_df |>
      dplyr::mutate(plab = .add_plab(padj, plab_fmt), .after = "padj")
  }

  return(result_df)
}


#### difference_analysis ####
# 对 profile 做两组差异分析，同时输出均值、FC、log2FC、prevalence 和 p/q 值
#
# profile:
#   行为 feature / metabolite / taxon，列为 sample
#
# group:
#   样本分组信息表
#
# sample_col:
#   group 中样本列名，默认 sample
#
# group_col:
#   group 中分组列名，默认 group
#
# comparison:
#   长度为 2 的比较组，例如 c('Treat', 'Control')
#   FC = comparison[1] / comparison[2]
#   log2FC > 0 表示 comparison[1] 更高
#   log2FC < 0 表示 comparison[2] 更高
#
# input_scale:
#   输入 profile 当前的数据尺度
#   'raw'   : 原始丰度 / 相对丰度 / peak area / intensity
#   'log10' : 已经是 log10 转换后的数据
#   'log2'  : 已经是 log2 转换后的数据
#
# test_trans:
#   差异检验时使用的数据转换方式
#   'none'  : 直接使用输入 profile 做检验
#   'log10' : 使用 raw 尺度数据 log10 转换后做检验
#   'log2'  : 使用 raw 尺度数据 log2 转换后做检验
#   'sqrt'  : 使用 raw 尺度数据 sqrt 转换后做检验
#   'ra'    : 使用 raw 尺度数据转相对丰度后做检验
#   'clr'   : 使用 raw 尺度数据 CLR 转换后做检验
#
# fc_method:
#   'arithmetic':
#     使用 raw 尺度算术均值比值
#     FC = mean_raw_group1 / mean_raw_group2
#     适合微生物相对丰度常规描述
#
#   'geometric':
#     使用几何均值比值
#     FC = geometric_mean_group1 / geometric_mean_group2
#     适合代谢组 peak area / intensity，尤其是 log 转换数据
#
# method:
#   'wilcox': Wilcoxon rank-sum test
#   't'     : t-test；var_equal = FALSE 时为 Welch t-test
#
# exact:
#   传递给 wilcox.test()
#
# var_equal:
#   传递给 t.test()
#   FALSE 为 Welch t-test
#   TRUE 为 Student's t-test
#
# min_abundance:
#   prevalence 判断阈值
#
# fc_pseudo:
#   arithmetic FC 计算时加入的伪值，避免分母为 0
#   微生物相对丰度常用 1e-6
#
# log_pseudo_factor:
#   log 转换或 geometric FC 计算时，0 或非正值替换为该 feature 最小正值的倍数
#   代谢组常用 0.5
#
# 输出：
#   name        : feature 名称
#   comparison  : 比较名称
#   avg_ab1     : group1 的均值；arithmetic 为算术均值，geometric 为几何均值
#   avg_ab2     : group2 的均值；arithmetic 为算术均值，geometric 为几何均值
#   prev1       : group1 中 feature > min_abundance 的样本比例
#   prev2       : group2 中 feature > min_abundance 的样本比例
#   FC          : group1 / group2
#   log2FC      : log2(FC)
#   pval        : p 值
#   padj        : BH 校正 q 值
#   method      : 检验方法
# 输出结果主要列：
#   name        : feature 名称
#   comparison  : 比较名称，例如 Treat_vs_Control
#   avg_ab1     : group1 的丰度估计值
#                 arithmetic 时为算术均值
#                 geometric 时为几何均值
#   avg_ab2     : group2 的丰度估计值
#   prev1       : group1 中 feature prevalence
#   prev2       : group2 中 feature prevalence
#   mean_input1 : group1 在输入尺度下的均值
#   mean_input2 : group2 在输入尺度下的均值
#   mean_raw1   : group1 反变换到 raw 尺度后的算术均值
#   mean_raw2   : group2 反变换到 raw 尺度后的算术均值
#   FC          : group1 / group2 的 fold change
#   log2FC      : log2(FC)
#   pval        : 差异检验 p 值
#   padj        : BH 校正后的 q 值
#   method      : 使用的检验方法
#
# 推荐用法示例：
#
# 1. 微生物相对丰度，Wilcoxon 检验
# res_micro <- difference_analysis(
#   profile = profile_ra,
#   group = group,
#   comparison = c('Treat', 'Control'),
#   input_scale = 'raw',
#   test_trans = 'none',
#   fc_method = 'arithmetic',
#   method = 'wilcox',
#   fc_pseudo = 1e-6
# )
#
# 2. 微生物相对丰度，log10 转换后 Welch t-test
# res_micro_t <- difference_analysis2(
#   profile = profile_ra,
#   group = group,
#   comparison = c('Treat', 'Control'),
#   input_scale = 'raw',
#   test_trans = 'log10',
#   fc_method = 'arithmetic',
#   method = 't',
#   var_equal = FALSE,
#   fc_pseudo = 1e-6
# )
#
# 3. 代谢组原始 peak area / intensity，log10 转换后 Welch t-test，
#    FC 使用几何均值倍数
# res_meta_raw <- difference_analysis2(
#   profile = profile_raw,
#   group = group,
#   comparison = c('DSS', 'Control'),
#   input_scale = 'raw',
#   test_trans = 'log10',
#   fc_method = 'geometric',
#   method = 't',
#   var_equal = FALSE,
#   log_pseudo_factor = 0.5
# )
#
# 4. 代谢组已经 log10 转换后的 intensity
# res_meta_log10 <- difference_analysis2(
#   profile = profile_log10,
#   group = group,
#   comparison = c('DSS', 'Control'),
#   input_scale = 'log10',
#   test_trans = 'none',
#   fc_method = 'geometric',
#   method = 't',
#   var_equal = FALSE
# )
#
# 注意事项：
# 1. comparison 的顺序决定 FC 方向：
#    comparison = c('A', 'B') 时，FC = A / B
#
# 2. 微生物相对丰度是组成型数据，mean RA FC 适合作为描述性指标，
#    不应过度解释为绝对丰度变化。
#    如果需要更严格的组成型差异分析，可结合 ANCOM-BC、ALDEx2、
#    MaAsLin 或 LinDA 等方法。
#
# 3. 代谢组数据如果已经 log10 转换，不能直接用 log10 均值相除计算 FC。
#    推荐使用：
#    FC = 10^(mean_log10_group1 - mean_log10_group2)
#
# 4. 差异检验使用的数据尺度和 FC 使用的数据尺度可以不同。
#    例如代谢组常用 log10 数据做 t-test，但 FC 用几何均值倍数解释。
#
# 5. z-score 数据不适合用于 FC 计算，也不建议作为差异检验的原始输入。

#' Difference Analysis utility
#'
#'
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param sample_meta A sample metadata table containing sample and group columns.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param comparison Two outcome levels ordered as case and control.
#' @param input_scale Scale of the supplied abundances: raw, log10, or log2.
#' @param test_trans Transformation applied to feature values before statistical testing.
#' @param fc_method Method used to calculate fold change; supported values are shown in Usage.
#' @param method Analysis or summary method; supported values are shown in the usage.
#' @param exact Whether an exact test is requested when supported by the selected method.
#' @param var_equal Whether two-sample t tests assume equal group variances.
#' @param min_abundance Minimum abundance above which a feature is considered present.
#' @param fc_pseudo Pseudocount added before calculating fold change.
#' @param log_pseudo_factor Multiplier used to derive a log-scale pseudocount from the minimum positive value.
#' @param progress Whether progress information is printed during repeated analyses.
#' @param digits Optional number of decimal digits retained.
#' @param ... Additional arguments forwarded to `calcu_diff_profile()`.
#' @return A feature-level data frame combining abundance, prevalence, fold-change, and hypothesis-test results.
#' @export
difference_analysis <- function(
  profile, sample_meta, sample_col = "sample", group_col = "group", comparison = NULL,
  input_scale = c("raw", "log10", "log2"),
  test_trans = c("none", "log10", "log2", "sqrt", "ra", "clr"),
  fc_method = c("arithmetic", "geometric"),
  method = c("wilcox", "t"), exact = NULL, var_equal = FALSE,
  min_abundance = 0, fc_pseudo = 1e-6, log_pseudo_factor = 0.5,
  progress = TRUE, digits = 6, ...
) {
  group <- sample_meta
  input_scale <- match.arg(input_scale)
  test_trans <- match.arg(test_trans)
  fc_method <- match.arg(fc_method)
  method <- match.arg(method)

  profile <- data.frame(profile, check.names = FALSE)
  group <- data.frame(group, check.names = FALSE)

  ## 1. group 检查和内部统一
  if (!all(c(sample_col, group_col) %in% colnames(group))) {
    stop("group should contain columns: ", sample_col, " | ", group_col)
  }

  if (is.null(comparison) || length(comparison) != 2) {
    stop("comparison should be c(group1, group2).")
  }

  group_use <- group |>
    dplyr::select(
      sample = dplyr::all_of(sample_col),
      group = dplyr::all_of(group_col)
    ) |>
    dplyr::mutate(
      sample = as.character(.data$sample),
      group = as.character(.data$group)
    ) |>
    dplyr::filter(
      .data$group %in% comparison,
      .data$sample %in% colnames(profile)
    )

  if (nrow(group_use) == 0) {
    stop("No matched samples between group and profile.")
  }

  group_use$group <- factor(group_use$group, levels = comparison)
  group_use <- group_use |> dplyr::arrange(.data$group)

  sample_n <- table(factor(group_use$group, levels = comparison))

  if (any(sample_n == 0)) {
    stop(
      "Each comparison group should contain samples. Current sample size: ",
      paste0(names(sample_n), "=", sample_n, collapse = ", ")
    )
  }

  message(
    "[", format(Sys.time()), "] Sample: ",
    paste0(paste0(names(sample_n), " (n=", sample_n, ")"), collapse = ", ")
  )

  ## 2. 提取并整理矩阵
  profile <- profile[, group_use$sample, drop = FALSE]

  mat_input <- as.matrix(profile)
  storage.mode(mat_input) <- "numeric"

  if (is.null(rownames(mat_input))) {
    rownames(mat_input) <- paste0("feature_", seq_len(nrow(mat_input)))
  }

  ## 根据输入尺度过滤 feature
  if (input_scale == "raw") {
    keep <- rowSums(is.finite(mat_input) & mat_input != 0) > 0
  } else {
    keep <- rowSums(is.finite(mat_input)) > 0
  }

  mat_input <- mat_input[keep, , drop = FALSE]

  if (nrow(mat_input) == 0) {
    stop("No features left after filtering.")
  }

  feature_names <- rownames(mat_input)

  ## 3. 根据 input_scale 反推 raw 尺度矩阵
  mat_raw <- switch(input_scale,
    raw = mat_input,
    log10 = 10^mat_input,
    log2 = 2^mat_input
  )

  ## 4. 为差异检验准备转换矩阵
  mat_test <- mat_input
  if (test_trans != "none") {
    mat_test <- mat_raw

    if (test_trans == "ra") {
      mat_test[!is.finite(mat_test)] <- 0
      library_size <- colSums(mat_test, na.rm = TRUE)
      library_size[library_size == 0] <- 1
      mat_test <- sweep(mat_test, 2, library_size, "/")
    }

    if (test_trans == "sqrt") {
      mat_test[!is.finite(mat_test)] <- 0
      mat_test[mat_test < 0] <- 0
      mat_test <- sqrt(mat_test)
    }

    if (test_trans %in% c("log10", "log2", "clr")) {
      mat_test <- .replace_nonpositive(
        mat_test,
        pseudo_factor = log_pseudo_factor
      )
      mat_test <- switch(test_trans,
        log10 = log10(mat_test),
        log2 = log2(mat_test),
        clr = {
          log_mat <- log(mat_test)
          sweep(log_mat, 2, colMeans(log_mat, na.rm = TRUE), "-")
        }
      )
    }
  }

  ## 6. 逐 feature 差异检验
  message("[", format(Sys.time()), "] Assessing P-value.")

  difference <- purrr::map_dfr(
    feature_names,
    function(f) {
      df <- data.frame(
        value = as.numeric(mat_test[f, ]),
        group = group_use$group
      )

      test <- tryCatch(
        {
          if (method == "wilcox") {
            stats::wilcox.test(
              value ~ group,
              data = df, exact = exact, ...
            ) |>
              suppressWarnings()
          } else {
            stats::t.test(
              value ~ group,
              data = df, var.equal = var_equal, ...
            ) |>
              suppressWarnings()
          }
        },
        error = function(e) NULL
      )

      data.frame(
        name = f,
        pval = if (is.null(test)) NA_real_ else test$p.value,
        method = ifelse(
          method == "wilcox",
          "Wilcoxon rank-sum test",
          ifelse(isTRUE(var_equal), "Student's t-test", "Welch t-test")
        ),
        check.names = FALSE
      )
    },
    .progress = progress
  ) |>
    dplyr::mutate(
      padj = stats::p.adjust(.data$pval, method = "BH"),
      .after = "pval"
    )

  ## 7. 计算均值、FC 和 log2FC
  message("[", format(Sys.time()), "] Mean abundance / intensity and FC.")

  idx1 <- group_use$group == comparison[1]
  idx2 <- group_use$group == comparison[2]

  if (fc_method == "arithmetic") {
    avg_ab1 <- rowMeans(mat_raw[, idx1, drop = FALSE], na.rm = TRUE)
    avg_ab2 <- rowMeans(mat_raw[, idx2, drop = FALSE], na.rm = TRUE)

    FC <- (avg_ab1 + fc_pseudo) / (avg_ab2 + fc_pseudo)
    log2FC <- log2(FC)
  }

  if (fc_method == "geometric") {
    mat_log2_for_fc <- log2(
      .replace_nonpositive(mat_raw, pseudo_factor = log_pseudo_factor)
    )

    mean_log2_1 <- rowMeans(mat_log2_for_fc[, idx1, drop = FALSE], na.rm = TRUE)
    mean_log2_2 <- rowMeans(mat_log2_for_fc[, idx2, drop = FALSE], na.rm = TRUE)

    avg_ab1 <- 2^mean_log2_1
    avg_ab2 <- 2^mean_log2_2

    log2FC <- mean_log2_1 - mean_log2_2
    FC <- 2^log2FC
  }

  ## 8. Prevalence
  message("[", format(Sys.time()), "] Prevalence.")

  prev_mat <- is.finite(mat_raw) & mat_raw > min_abundance

  prevalence <- data.frame(
    name = feature_names,
    prev1 = rowMeans(prev_mat[, idx1, drop = FALSE]),
    prev2 = rowMeans(prev_mat[, idx2, drop = FALSE]),
    check.names = FALSE
  )

  ## 9. 合并结果
  message("[", format(Sys.time()), "] Output result.")

  result <- data.frame(
    name = feature_names,
    comparison = paste0(comparison, collapse = "_vs_"),
    avg_ab1 = avg_ab1,
    avg_ab2 = avg_ab2,
    FC = FC,
    log2FC = log2FC,
    input_scale = input_scale,
    test_trans = test_trans,
    fc_method = fc_method,
    check.names = FALSE
  ) |>
    dplyr::left_join(prevalence, by = "name") |>
    dplyr::left_join(difference, by = "name") |>
    dplyr::relocate(prev1, prev2, .after = avg_ab2) |>
    dplyr::mutate(
      dplyr::across(
        c(avg_ab1, avg_ab2, prev1, prev2, FC, log2FC),
        \(x) round(x, digits)
      )
    )

  message("[", format(Sys.time()), "] end ...")

  return(result)
}
