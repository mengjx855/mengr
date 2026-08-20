#### Jin-Xin Meng, 20220529, 20260820, v0.4.0 ####

.prepare_rf_data <- function(
    profile, group, sample_col = 'sample', group_col = 'group',
    group_level = NULL, na_fill = NULL, remove_zero_var = TRUE) {
  ## 对齐样本并转为 sample × feature 数值表
  aligned <- .align_profile_group(
    profile = profile, group = group,
    sample_col = sample_col, group_col = group_col,
    group_level = group_level
  )
  x_df <- data.frame(
    t(.as_profile_df(aligned$profile_df, numeric = TRUE)),
    check.names = FALSE
  )
  if (anyNA(x_df)) {
    if (is.null(na_fill)) stop('NA values found in profile; set na_fill if needed.')
    x_df[is.na(x_df)] <- na_fill
  }
  if (isTRUE(remove_zero_var)) {
    keep_feature <- vapply(x_df, stats::var, numeric(1), na.rm = TRUE) > 0
    x_df <- x_df[, keep_feature, drop = FALSE]
  }
  if (!ncol(x_df)) stop('No valid features remained for random forest.')

  ## randomForest 对复杂 feature 名不稳定，内部改名并保留映射
  feature_df <- data.frame(
    feature = colnames(x_df),
    model_feature = make.names(colnames(x_df), unique = TRUE),
    check.names = FALSE
  )
  colnames(x_df) <- feature_df$model_feature
  y_factor <- droplevels(aligned$group_df[[group_col]])
  if (nlevels(y_factor) < 2L) {
    stop('Random forest classification requires at least two groups.')
  }
  list(
    x_df = x_df, y_factor = y_factor,
    sample_vec = aligned$sample_vec, feature_df = feature_df
  )
}

.rf_probability_df <- function(model_obj, test_x, level_vec) {
  probability_df <- data.frame(
    stats::predict(model_obj, test_x, type = 'prob'),
    check.names = FALSE
  )
  missing_level <- setdiff(level_vec, colnames(probability_df))
  for (level in missing_level) probability_df[[level]] <- 0
  probability_df[, level_vec, drop = FALSE]
}

#' Generate stratified K-fold random-forest predictions
#'
#' Chinese summary: 执行分层 random-forest k-fold cross-validation。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param group A sample metadata table containing sample and group columns.
#' @param k Number of folds or clusters, according to the analysis performed.
#' @param seed Optional random seed for reproducibility.
#' @param ntree Number of trees fitted by the random-forest model.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param group_level Optional order of group levels.
#' @param remove_zero_var Logical control for `remove_zero_var`.
#' @param na_fill Value used to replace missing observations before analysis.
#' @param progress Whether progress information is printed during repeated analyses.
#' @param ... Additional arguments passed to the underlying function.
#' @return A result object described in the Details section.
#' @export
rf_kfold <- function(
    profile, group, k = 5, seed = 2025, ntree = 1000,
    sample_col = 'sample', group_col = 'group', group_level = NULL,
    remove_zero_var = TRUE, na_fill = NULL, progress = TRUE, ...) {
  ## 1. 统一输入并创建 stratified folds
  prepared <- .prepare_rf_data(
    profile, group, sample_col, group_col, group_level,
    na_fill, remove_zero_var
  )
  if (k < 2L || k > nrow(prepared$x_df)) {
    stop('k should be between 2 and the number of matched samples.')
  }
  set.seed(seed)
  fold_list <- caret::createFolds(
    prepared$y_factor, k = k, returnTrain = FALSE
  )
  result_list <- vector('list', length(fold_list))
  if (isTRUE(progress)) {
    progress_bar <- utils::txtProgressBar(
      min = 0, max = length(fold_list), style = 3, width = 50, char = '#'
    )
    on.exit(close(progress_bar), add = TRUE)
  }

  ## 2. 每折独立训练，并预测留出样本
  for (fold_idx in seq_along(fold_list)) {
    test_idx <- fold_list[[fold_idx]]
    train_idx <- setdiff(seq_len(nrow(prepared$x_df)), test_idx)
    train_y <- droplevels(prepared$y_factor[train_idx])
    if (nlevels(train_y) < 2L) {
      stop('A training fold contains fewer than two groups; reduce k.')
    }
    set.seed(seed + fold_idx - 1L)
    model_obj <- randomForest::randomForest(
      x = prepared$x_df[train_idx, , drop = FALSE],
      y = train_y, ntree = ntree,
      importance = FALSE, proximity = FALSE, ...
    )
    probability_df <- .rf_probability_df(
      model_obj, prepared$x_df[test_idx, , drop = FALSE],
      levels(prepared$y_factor)
    )
    result_list[[fold_idx]] <- data.frame(
      sample = prepared$sample_vec[test_idx],
      actual = as.character(prepared$y_factor[test_idx]),
      predicted = colnames(probability_df)[max.col(probability_df)],
      fold = fold_idx, probability_df, check.names = FALSE
    )
    if (isTRUE(progress)) utils::setTxtProgressBar(progress_bar, fold_idx)
  }
  dplyr::bind_rows(result_list)
}

#' Repeat stratified K-fold random-forest validation
#'
#' Chinese summary: 多次重复 random-forest k-fold，并标记每次重复和 fold。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param group A sample metadata table containing sample and group columns.
#' @param k Number of folds or clusters, according to the analysis performed.
#' @param repeats Number of repeated cross-validation runs.
#' @param seed Optional random seed for reproducibility.
#' @param ntree Number of trees fitted by the random-forest model.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param group_level Optional order of group levels.
#' @param remove_zero_var Logical control for `remove_zero_var`.
#' @param na_fill Value used to replace missing observations before analysis.
#' @param progress Whether progress information is printed during repeated analyses.
#' @param ... Additional arguments passed to the underlying function.
#' @return A result object described in the Details section.
#' @export
rf_repeated_kfold <- function(
    profile, group, k = 5, repeats = 5, seed = 2026, ntree = 1000,
    sample_col = 'sample', group_col = 'group', group_level = NULL,
    remove_zero_var = TRUE, na_fill = NULL, progress = TRUE, ...) {
  ## 1. 重复执行 K-fold，并记录每次预测
  result_list <- vector('list', repeats)
  for (repeat_idx in seq_len(repeats)) {
    repeat_df <- rf_kfold(
      profile = profile, group = group, k = k,
      seed = seed + repeat_idx - 1L, ntree = ntree,
      sample_col = sample_col, group_col = group_col,
      group_level = group_level, remove_zero_var = remove_zero_var,
      na_fill = na_fill, progress = FALSE, ...
    )
    repeat_df$repeat_id <- repeat_idx
    result_list[[repeat_idx]] <- repeat_df
  }
  all_df <- dplyr::bind_rows(result_list)
  info_col <- c('sample', 'actual', 'predicted', 'fold', 'repeat_id')
  probability_col <- setdiff(colnames(all_df), info_col)

  ## 2. 对每个样本的重复预测概率取均值
  result_df <- all_df |>
    dplyr::group_by(sample, actual) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(probability_col), \(x) mean(x, na.rm = TRUE)
      ),
      .groups = 'drop'
    )
  probability_mat <- as.matrix(result_df[, probability_col, drop = FALSE])
  result_df$predicted <- probability_col[max.col(probability_mat)]
  result_df
}

#' Validate random forests across datasets
#'
#' Chinese summary: 按 dataset 留一验证，评估跨队列 random-forest 泛化能力。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param group A sample metadata table containing sample and group columns.
#' @param dataset_col Name of the `dataset_col` input column.
#' @param dataset_level Optional order for `dataset_level`.
#' @param seed Optional random seed for reproducibility.
#' @param ntree Number of trees fitted by the random-forest model.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param positive_class Outcome level treated as the positive class for binary metrics.
#' @param na_fill Value used to replace missing observations before analysis.
#' @param ... Additional arguments passed to the underlying function.
#' @return A result object described in the Details section.
#' @export
rf_cross_dataset_validate <- function(
    profile, group, dataset_col = 'dataset', dataset_level = NULL,
    seed = 2026, ntree = 1000, sample_col = 'sample',
    group_col = 'group', positive_class = NULL, na_fill = 0, ...) {
  ## 1. 检查 metadata 并确定数据集顺序
  profile_df <- .as_profile_df(profile, numeric = TRUE)
  group_df <- .as_df(group)
  .check_columns(
    group_df, c(sample_col, group_col, dataset_col), object = 'group'
  )
  if (is.null(dataset_level)) {
    dataset_level <- unique(as.character(group_df[[dataset_col]]))
  }
  result_list <- list()
  result_idx <- 1L

  ## 2. 依次使用一个数据集训练，其他数据集验证
  for (train_dataset in dataset_level) {
    train_group_df <- group_df[
      group_df[[dataset_col]] == train_dataset, , drop = FALSE
    ]
    train <- .prepare_rf_data(
      profile_df, train_group_df, sample_col, group_col,
      na_fill = na_fill, remove_zero_var = TRUE
    )
    if (nlevels(train$y_factor) != 2L) next
    if (is.null(positive_class)) {
      current_positive <- levels(train$y_factor)[2]
    } else {
      current_positive <- positive_class
    }
    if (!current_positive %in% levels(train$y_factor)) next

    set.seed(seed)
    model_obj <- randomForest::randomForest(
      x = train$x_df, y = train$y_factor, ntree = ntree,
      importance = FALSE, proximity = FALSE, ...
    )
    for (test_dataset in setdiff(dataset_level, train_dataset)) {
      test_group_df <- group_df[
        group_df[[dataset_col]] == test_dataset, , drop = FALSE
      ]
      test <- .prepare_rf_data(
        profile_df, test_group_df, sample_col, group_col,
        group_level = levels(train$y_factor),
        na_fill = na_fill, remove_zero_var = FALSE
      )
      feature_vec <- train$feature_df$feature
      feature_vec <- feature_vec[feature_vec %in% test$feature_df$feature]
      if (!length(feature_vec)) next
      train_col <- train$feature_df$model_feature[
        match(feature_vec, train$feature_df$feature)
      ]
      test_col <- test$feature_df$model_feature[
        match(feature_vec, test$feature_df$feature)
      ]
      test_x <- test$x_df[, test_col, drop = FALSE]
      colnames(test_x) <- train_col
      test_x <- test_x[, colnames(train$x_df), drop = FALSE]

      probability_df <- .rf_probability_df(
        model_obj, test_x, levels(train$y_factor)
      )
      auc_value <- NA_real_
      if (length(unique(test$y_factor)) == 2L) {
        roc_obj <- pROC::roc(
          test$y_factor, probability_df[[current_positive]],
          levels = levels(train$y_factor), quiet = TRUE
        )
        auc_value <- as.numeric(pROC::auc(roc_obj))
      }
      result_list[[result_idx]] <- data.frame(
        train_dataset = train_dataset, test_dataset = test_dataset,
        seed = seed, ntree = ntree, auc = auc_value,
        positive_class = current_positive, check.names = FALSE
      )
      result_idx <- result_idx + 1L
    }
  }
  dplyr::bind_rows(result_list)
}

#' Calculate random-forest feature importance
#'
#' Chinese summary: 提取并整理 random-forest feature importance。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param group A sample metadata table containing sample and group columns.
#' @param seed Optional random seed for reproducibility.
#' @param ntree Number of trees fitted by the random-forest model.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param group_level Optional order of group levels.
#' @param na_fill Value used to replace missing observations before analysis.
#' @param ... Additional arguments passed to the underlying function.
#' @return A result object described in the Details section.
#' @export
rf_importance <- function(
    profile, group, seed = 2026, ntree = 1000,
    sample_col = 'sample', group_col = 'group', group_level = NULL,
    na_fill = 0, ...) {
  ## 1. 预处理并拟合全数据模型
  prepared <- .prepare_rf_data(
    profile, group, sample_col, group_col, group_level,
    na_fill, remove_zero_var = TRUE
  )
  set.seed(seed)
  model_obj <- randomForest::randomForest(
    x = prepared$x_df, y = prepared$y_factor, ntree = ntree,
    importance = TRUE, proximity = FALSE, ...
  )

  ## 2. 把内部安全列名还原成原始 feature 名
  result_df <- randomForest::importance(model_obj) |>
    data.frame(check.names = FALSE) |>
    tibble::rownames_to_column('model_feature') |>
    dplyr::left_join(prepared$feature_df, by = 'model_feature') |>
    dplyr::select(feature, dplyr::everything(), -model_feature)
  if ('MeanDecreaseAccuracy' %in% colnames(result_df)) {
    result_df <- dplyr::arrange(result_df, dplyr::desc(MeanDecreaseAccuracy))
  }
  attr(result_df, 'model') <- model_obj
  result_df
}

#' Leave-one-out random-forest predictions
#'
#' Chinese summary: 执行 leave-one-out random-forest 预测。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param group A sample metadata table containing sample and group columns.
#' @param seed Optional random seed for reproducibility.
#' @param ntree Number of trees fitted by the random-forest model.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param group_level Optional order of group levels.
#' @param na_fill Value used to replace missing observations before analysis.
#' @param progress Whether progress information is printed during repeated analyses.
#' @param ... Additional arguments passed to the underlying function.
#' @return A result object described in the Details section.
#' @export
rf_leave_one_out <- function(
    profile, group, seed = 2026, ntree = 1000,
    sample_col = 'sample', group_col = 'group', group_level = NULL,
    na_fill = 0, progress = TRUE, ...) {
  ## 1. 预处理并逐个留出样本
  prepared <- .prepare_rf_data(
    profile, group, sample_col, group_col, group_level,
    na_fill, remove_zero_var = TRUE
  )
  result_list <- vector('list', nrow(prepared$x_df))
  if (isTRUE(progress)) {
    progress_bar <- utils::txtProgressBar(
      min = 0, max = nrow(prepared$x_df), style = 3, width = 50, char = '#'
    )
    on.exit(close(progress_bar), add = TRUE)
  }
  for (sample_idx in seq_len(nrow(prepared$x_df))) {
    train_idx <- setdiff(seq_len(nrow(prepared$x_df)), sample_idx)
    train_y <- droplevels(prepared$y_factor[train_idx])
    if (nlevels(train_y) < 2L) {
      probability_df <- as.data.frame(matrix(
        NA_real_, nrow = 1, ncol = nlevels(prepared$y_factor),
        dimnames = list(NULL, levels(prepared$y_factor))
      ))
    } else {
      set.seed(seed + sample_idx - 1L)
      model_obj <- randomForest::randomForest(
        x = prepared$x_df[train_idx, , drop = FALSE],
        y = train_y, ntree = ntree,
        importance = FALSE, proximity = FALSE, ...
      )
      probability_df <- .rf_probability_df(
        model_obj, prepared$x_df[sample_idx, , drop = FALSE],
        levels(prepared$y_factor)
      )
    }
    result_list[[sample_idx]] <- data.frame(
      sample = prepared$sample_vec[sample_idx],
      actual = as.character(prepared$y_factor[sample_idx]),
      probability_df, check.names = FALSE
    )
    if (isTRUE(progress)) utils::setTxtProgressBar(progress_bar, sample_idx)
  }
  result_df <- dplyr::bind_rows(result_list)
  probability_col <- levels(prepared$y_factor)
  result_df$predicted <- probability_col[
    max.col(as.matrix(result_df[, probability_col, drop = FALSE]))
  ]
  missing_pred <- rowSums(is.finite(
    as.matrix(result_df[, probability_col, drop = FALSE])
  )) == 0
  result_df$predicted[missing_pred] <- NA_character_
  result_df
}

#' Train a random forest in one dataset and validate it in another
#'
#' Chinese summary: 用训练数据拟合 random forest，并在独立数据中验证。
#'
#' @param profile_x Feature-by-sample profile used for model training or the first data space.
#' @param profile_y Feature-by-sample profile used for validation or the second data space.
#' @param group_x Sample metadata associated with `profile_x`.
#' @param group_y Sample metadata associated with `profile_y`.
#' @param seed Optional random seed for reproducibility.
#' @param ntree Number of trees fitted by the random-forest model.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param group_level Optional order of group levels.
#' @param na_fill Value used to replace missing observations before analysis.
#' @param ... Additional arguments passed to the underlying function.
#' @return A result object described in the Details section.
#' @export
rf_next_validate <- function(
    profile_x, profile_y, group_x, group_y,
    seed = 2026, ntree = 1000,
    sample_col = 'sample', group_col = 'group',
    group_level = NULL, na_fill = 0, ...) {
  ## 1. 仅保留训练集和验证集共有 feature
  feature_vec <- intersect(rownames(profile_x), rownames(profile_y))
  if (!length(feature_vec)) {
    stop('No matched features between profile_x and profile_y.')
  }
  train <- .prepare_rf_data(
    profile_x[feature_vec, , drop = FALSE], group_x,
    sample_col, group_col, group_level, na_fill, TRUE
  )
  test <- .prepare_rf_data(
    profile_y[feature_vec, , drop = FALSE], group_y,
    sample_col, group_col, levels(train$y_factor), na_fill, FALSE
  )

  ## 2. 按训练集保留的 feature 对齐验证矩阵
  selected_feature <- train$feature_df$feature
  test_col <- test$feature_df$model_feature[
    match(selected_feature, test$feature_df$feature)
  ]
  if (anyNA(test_col)) stop('Validation data are missing selected features.')
  test_x <- test$x_df[, test_col, drop = FALSE]
  colnames(test_x) <- train$feature_df$model_feature

  ## 3. 拟合模型并输出各类别概率
  set.seed(seed)
  model_obj <- randomForest::randomForest(
    x = train$x_df, y = train$y_factor, ntree = ntree,
    importance = FALSE, proximity = FALSE, ...
  )
  probability_df <- .rf_probability_df(
    model_obj, test_x, levels(train$y_factor)
  )
  result_df <- data.frame(
    sample = test$sample_vec, actual = as.character(test$y_factor),
    predicted = colnames(probability_df)[max.col(probability_df)],
    probability_df, check.names = FALSE
  )
  attr(result_df, 'model') <- model_obj
  result_df
}
