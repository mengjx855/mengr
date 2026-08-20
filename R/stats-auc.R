#### Jinxin Meng, 20260520, 20260520, v0.0.1 ####

# 20260520 v0.0.1: add functions 'calcu_feature_auc()'

#### calcu_feature_auc ####
# 计算每个 feature 的单特征 AUC，只用一个 feature 的数值作为 predictor，
# 查看它区分两组样本的能力，这表示 pROC 会自动选择方向，使 AUC 通常大于 0.5。
# 因此这里的 AUC 更偏向表示区分能力强弱，而不是方向

# profile:
#   行为 feature，列为 sample 的丰度表、表达矩阵或其他数值矩阵。
# label:
#   可以是命名向量，例如：
#     label = c(S1 = "HC", S2 = "IBD", S3 = "HC")
#   也可以是数据框，默认包含 sample_col 和 group_col 两列。
# levels:
#   二分类顺序，格式为 c(control, case)。
#   例如 levels = c("HC", "IBD")。
# direction:
#   传递给 pROC::roc()。
#   "auto" 表示自动判断方向，使 AUC 通常 >= 0.5。
# signal:
#   signal = 2 * abs(AUC - 0.5)
#   表示单 feature 对分组的区分强度，不关心方向。
# 输出:
#   feature, auc, ci_low, ci_high, signal, mean_control, mean_case, enriched

#' Calcu Feature Auc utility
#'
#' `calcu_feature_auc()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 逐 feature 计算二分类 ROC AUC、置信区间、方向和区分强度。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param label Outcome labels corresponding to the profile samples.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param levels Numeric setting for `levels`.
#' @param direction Direction used for prediction, ROC calculation, or network extraction.
#' @param ci Whether to calculate confidence intervals for AUC estimates.
#' @param quiet Whether to suppress progress messages.
#' @return A result object described in the Details section.
#' @export
calcu_feature_auc <- function(profile, label, sample_col = 'sample', 
                              group_col = 'group', levels = NULL,
                              direction = c('auto', '<', '>'), ci = TRUE,
                              quiet = TRUE) {
  
  direction <- match.arg(direction)
  profile <- data.frame(profile, check.names = FALSE)
  
  ## 1. 整理 label
  if (is.data.frame(label)) {
    
    label <- data.frame(label, check.names = FALSE)
    
    if (!all(c(sample_col, group_col) %in% colnames(label))) {
      stop('label data.frame should contain columns: ', sample_col, ' | ', group_col)
    }
    
    label_vec <- label[[group_col]]
    names(label_vec) <- label[[sample_col]]
    
  } else {
    
    label_vec <- label
    
    if (is.null(names(label_vec))) {
      stop('label should be a named vector or a data.frame with sample/group columns.')
    }
  }
  
  ## 2. 对齐样本
  common <- intersect(colnames(profile), names(label_vec))
  
  if (length(common) < 3) {
    stop('Too few matched samples between profile and label.')
  }
  
  profile <- profile[, common, drop = FALSE]
  label_vec <- label_vec[common]
  
  ## 3. 设置二分类水平
  if (is.null(levels)) {
    levels <- unique(as.character(label_vec))
  }
  
  if (length(levels) != 2) {
    stop('levels should contain exactly two groups, such as c("HC", "IBD").')
  }
  
  label_vec <- factor(label_vec, levels = levels)
  
  ## 4. 逐个 feature 计算 AUC
  res <- lapply(rownames(profile), \(f) {
    
    x <- suppressWarnings(as.numeric(profile[f, ]))
    y <- label_vec
    
    ok <- is.finite(x) & !is.na(y)
    x <- x[ok]
    y <- y[ok]
    
    ## 如果只有一个分组，或 feature 没有变化，则无法计算 ROC
    if (length(unique(y)) < 2 || length(unique(x)) < 2) {
      return(data.frame(
        feature = f,
        auc = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        signal = NA_real_,
        mean_control = mean(x[y == levels[1]], na.rm = TRUE),
        mean_case = mean(x[y == levels[2]], na.rm = TRUE),
        enriched = NA_character_
      ))
    }
    
    roc_obj <- tryCatch(
      pROC::roc(
        response = y,
        predictor = x,
        levels = levels,
        direction = direction,
        quiet = quiet
      ),
      error = function(e) NULL
    )
    
    if (is.null(roc_obj)) {
      return(data.frame(
        feature = f,
        auc = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        signal = NA_real_,
        mean_control = mean(x[y == levels[1]], na.rm = TRUE),
        mean_case = mean(x[y == levels[2]], na.rm = TRUE),
        enriched = NA_character_
      ))
    }
    
    auc_val <- as.numeric(pROC::auc(roc_obj))
    
    if (isTRUE(ci)) {
      auc_ci <- tryCatch(
        as.numeric(pROC::ci.auc(roc_obj)),
        error = function(e) c(NA_real_, NA_real_, NA_real_)
      )
    } else {
      auc_ci <- c(NA_real_, NA_real_, NA_real_)
    }
    
    mean_control <- mean(x[y == levels[1]], na.rm = TRUE)
    mean_case <- mean(x[y == levels[2]], na.rm = TRUE)
    
    enriched <- dplyr::case_when(
      mean_case > mean_control ~ levels[2],
      mean_case < mean_control ~ levels[1],
      TRUE ~ 'none'
    )
    
    data.frame(
      feature = f,
      auc = auc_val,
      ci_low = auc_ci[1],
      ci_high = auc_ci[3],
      signal = 2 * abs(auc_val - 0.5),
      mean_control = mean_control,
      mean_case = mean_case,
      enriched = enriched,
      stringsAsFactors = FALSE
    )
  })
  
  res <- dplyr::bind_rows(res) |>
    dplyr::arrange(dplyr::desc(signal))
  
  return(res)
}

