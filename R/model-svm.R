#### Jin-Xin Meng, 20241024, 20260820, v0.2.0 ####

.prepare_svm_data <- function(
    profile, group, sample_col = 'sample', group_col = 'group',
    group_level = NULL) {
  aligned <- .align_profile_group(
    profile = profile, group = group,
    sample_col = sample_col, group_col = group_col,
    group_level = group_level
  )
  group_factor <- droplevels(aligned$group_df[[group_col]])
  if (nlevels(group_factor) != 2L) {
    stop('SVM functions currently require exactly two groups.')
  }
  list(
    x_mat = t(as.matrix(.as_profile_df(aligned$profile_df, numeric = TRUE))),
    group_factor = group_factor,
    sample_vec = aligned$sample_vec,
    profile_df = aligned$profile_df
  )
}

.fit_tuned_svm <- function(
    x_mat, y, kernel, scale, cost_grid, gamma_grid,
    tune_boot = 10, probability = TRUE) {
  tune_args <- list(
    x = x_mat, y = y, kernel = kernel,
    cost = cost_grid,
    tunecontrol = e1071::tune.control(
      sampling = 'bootstrap', nboot = tune_boot
    )
  )
  if (kernel != 'linear') tune_args$gamma <- gamma_grid
  tune_obj <- do.call(e1071::tune.svm, tune_args)

  fit_args <- list(
    x = x_mat, y = y, kernel = kernel, scale = scale,
    probability = probability,
    cost = tune_obj$best.parameters[['cost']]
  )
  if (kernel != 'linear') {
    fit_args$gamma <- tune_obj$best.parameters[['gamma']]
  }
  list(
    tune_obj = tune_obj,
    model_obj = do.call(e1071::svm, fit_args)
  )
}

.svm_prediction_df <- function(model_obj, x_mat, sample_vec, positive_class) {
  class_vec <- stats::predict(model_obj, x_mat, probability = TRUE)
  prob_mat <- attr(class_vec, 'probabilities')
  if (is.null(prob_mat) || !positive_class %in% colnames(prob_mat)) {
    stop('Probability column not found for positive_class: ', positive_class)
  }
  data.frame(
    sample = sample_vec,
    predicted = as.character(class_vec),
    probability = prob_mat[, positive_class],
    check.names = FALSE
  )
}

.binary_metrics <- function(actual, predicted, positive_class) {
  level_vec <- levels(actual)
  negative_class <- setdiff(level_vec, positive_class)
  confusion_mat <- table(
    predicted = factor(predicted, levels = level_vec),
    actual = factor(actual, levels = level_vec)
  )
  list(
    confusion_matrix = confusion_mat,
    accuracy = sum(diag(confusion_mat)) / sum(confusion_mat),
    sensitivity = confusion_mat[positive_class, positive_class] /
      sum(confusion_mat[, positive_class]),
    specificity = confusion_mat[negative_class, negative_class] /
      sum(confusion_mat[, negative_class])
  )
}

#' Repeated holdout validation for a binary SVM
#'
#' Chinese summary: 拟合带参数搜索的基础 SVM，并返回样本预测和性能。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param group A sample metadata table containing sample and group columns.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param group_level Optional order of group levels.
#' @param positive_class Outcome level treated as the positive class for binary metrics.
#' @param rep Number of repeated train/test runs.
#' @param seed Optional random seed for reproducibility.
#' @param train_prop Proportion of matched samples assigned to the training split.
#' @param kernel SVM kernel; supported values are shown in Usage.
#' @param scale Logical control for `scale`.
#' @param cost_grid Candidate SVM cost values used during tuning.
#' @param gamma_grid Candidate SVM gamma values used during tuning.
#' @param tune_boot Number of bootstrap resamples used during SVM tuning.
#' @return A result object described in the Details section.
#' @export
svm_base <- function(
    profile, group, sample_col = 'sample', group_col = 'group',
    group_level = NULL, positive_class = NULL,
    rep = 5, seed = 2024, train_prop = 0.8,
    kernel = c('radial', 'linear', 'polynomial', 'sigmoid'),
    scale = TRUE, cost_grid = 10^(-1:3), gamma_grid = 10^(-3:1),
    tune_boot = 10
) {
  kernel <- match.arg(kernel)
  prepared <- .prepare_svm_data(
    profile, group, sample_col, group_col, group_level
  )
  if (is.null(positive_class)) {
    positive_class <- levels(prepared$group_factor)[2]
  }
  if (!positive_class %in% levels(prepared$group_factor)) {
    stop('positive_class is not present in group.')
  }
  if (train_prop <= 0 || train_prop >= 1) {
    stop('train_prop should be between 0 and 1.')
  }

  best_auc <- -Inf
  out <- NULL
  for (rep_idx in seq_len(rep)) {
    model_seed <- seed + rep_idx - 1L
    set.seed(model_seed)
    train_idx <- caret::createDataPartition(
      prepared$group_factor, p = train_prop, list = FALSE
    )
    test_idx <- setdiff(seq_len(nrow(prepared$x_mat)), train_idx)
    if (!length(test_idx)) next

    fitted <- .fit_tuned_svm(
      prepared$x_mat[train_idx, , drop = FALSE],
      prepared$group_factor[train_idx], kernel, scale,
      cost_grid, gamma_grid, tune_boot
    )
    pred_df <- .svm_prediction_df(
      fitted$model_obj,
      prepared$x_mat[test_idx, , drop = FALSE],
      prepared$sample_vec[test_idx], positive_class
    )
    pred_df$actual <- as.character(prepared$group_factor[test_idx])
    roc_obj <- pROC::roc(
      response = factor(pred_df$actual, levels = levels(prepared$group_factor)),
      predictor = pred_df$probability,
      levels = levels(prepared$group_factor), quiet = TRUE
    )
    current_auc <- as.numeric(pROC::auc(roc_obj))
    if (current_auc > best_auc) {
      best_auc <- current_auc
      metrics <- .binary_metrics(
        factor(pred_df$actual, levels = levels(prepared$group_factor)),
        pred_df$predicted, positive_class
      )
      out <- c(list(
        seed = model_seed, model = fitted$model_obj,
        tune_model = fitted$tune_obj, pred = pred_df,
        roc = roc_obj, roc_plot = plot_roc(roc_obj),
        kernel = kernel
      ), metrics)
    }
  }
  if (is.null(out)) stop('No valid SVM model was produced.')
  out
}

#' K-fold cross-validated predictions for a binary SVM
#'
#' Chinese summary: 执行分层 SVM k-fold cross-validation。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param group A sample metadata table containing sample and group columns.
#' @param k Number of folds or clusters, according to the analysis performed.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param group_level Optional order of group levels.
#' @param positive_class Outcome level treated as the positive class for binary metrics.
#' @param seed Optional random seed for reproducibility.
#' @param kernel SVM kernel; supported values are shown in Usage.
#' @param scale Logical control for `scale`.
#' @param cost_grid Candidate SVM cost values used during tuning.
#' @param gamma_grid Candidate SVM gamma values used during tuning.
#' @param tune_boot Number of bootstrap resamples used during SVM tuning.
#' @return A result object described in the Details section.
#' @export
svm_kfold <- function(
    profile, group, k = 5, sample_col = 'sample', group_col = 'group',
    group_level = NULL, positive_class = NULL, seed = 2024,
    kernel = c('radial', 'linear', 'polynomial', 'sigmoid'),
    scale = TRUE, cost_grid = 10^(-1:3), gamma_grid = 10^(-3:1),
    tune_boot = 10
) {
  kernel <- match.arg(kernel)
  prepared <- .prepare_svm_data(
    profile, group, sample_col, group_col, group_level
  )
  if (is.null(positive_class)) {
    positive_class <- levels(prepared$group_factor)[2]
  }
  if (k < 2L || k > nrow(prepared$x_mat)) {
    stop('k should be between 2 and the number of matched samples.')
  }

  set.seed(seed)
  fold_list <- caret::createFolds(
    prepared$group_factor, k = k, returnTrain = FALSE
  )
  result_list <- vector('list', length(fold_list))
  for (fold_idx in seq_along(fold_list)) {
    test_idx <- fold_list[[fold_idx]]
    train_idx <- setdiff(seq_len(nrow(prepared$x_mat)), test_idx)
    set.seed(seed + fold_idx - 1L)
    fitted <- .fit_tuned_svm(
      prepared$x_mat[train_idx, , drop = FALSE],
      prepared$group_factor[train_idx], kernel, scale,
      cost_grid, gamma_grid, tune_boot
    )
    pred_df <- .svm_prediction_df(
      fitted$model_obj,
      prepared$x_mat[test_idx, , drop = FALSE],
      prepared$sample_vec[test_idx], positive_class
    )
    pred_df$actual <- as.character(prepared$group_factor[test_idx])
    pred_df$fold <- fold_idx
    result_list[[fold_idx]] <- pred_df
  }
  dplyr::bind_rows(result_list)
}

#' Train an SVM in one dataset and validate it in another
#'
#' Chinese summary: 在训练集调参拟合 SVM，并在独立数据中验证。
#'
#' @param profile_x Feature-by-sample profile used for model training or the first data space.
#' @param profile_y Feature-by-sample profile used for validation or the second data space.
#' @param group_x Sample metadata associated with `profile_x`.
#' @param group_y Sample metadata associated with `profile_y`.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param group_level Optional order of group levels.
#' @param positive_class Outcome level treated as the positive class for binary metrics.
#' @param seed Optional random seed for reproducibility.
#' @param kernel SVM kernel; supported values are shown in Usage.
#' @param scale Logical control for `scale`.
#' @param cost_grid Candidate SVM cost values used during tuning.
#' @param gamma_grid Candidate SVM gamma values used during tuning.
#' @param tune_boot Number of bootstrap resamples used during SVM tuning.
#' @return A result object described in the Details section.
#' @export
svm_next_validate <- function(
    profile_x, profile_y, group_x, group_y,
    sample_col = 'sample', group_col = 'group',
    group_level = NULL, positive_class = NULL, seed = 2024,
    kernel = c('radial', 'linear', 'polynomial', 'sigmoid'),
    scale = TRUE, cost_grid = 10^(-1:3), gamma_grid = 10^(-3:1),
    tune_boot = 20
) {
  kernel <- match.arg(kernel)
  train <- .prepare_svm_data(
    profile_x, group_x, sample_col, group_col, group_level
  )
  if (is.null(positive_class)) positive_class <- levels(train$group_factor)[2]
  test <- .prepare_svm_data(
    profile_y, group_y, sample_col, group_col,
    group_level = levels(train$group_factor)
  )
  if (!identical(levels(train$group_factor), levels(test$group_factor))) {
    stop('Training and validation datasets should contain the same two groups.')
  }

  feature_vec <- intersect(
    rownames(train$profile_df), rownames(test$profile_df)
  )
  if (!length(feature_vec)) {
    stop('No matched features between profile_x and profile_y.')
  }
  train_x <- t(as.matrix(train$profile_df[feature_vec, , drop = FALSE]))
  test_x <- t(as.matrix(test$profile_df[feature_vec, , drop = FALSE]))

  set.seed(seed)
  fitted <- .fit_tuned_svm(
    train_x, train$group_factor, kernel, scale,
    cost_grid, gamma_grid, tune_boot
  )
  pred_df <- .svm_prediction_df(
    fitted$model_obj, test_x, test$sample_vec, positive_class
  )
  pred_df$actual <- as.character(test$group_factor)
  roc_obj <- pROC::roc(
    response = test$group_factor,
    predictor = pred_df$probability,
    levels = levels(train$group_factor), quiet = TRUE
  )
  metrics <- .binary_metrics(
    test$group_factor, pred_df$predicted, positive_class
  )
  c(list(
    tune_model = fitted$tune_obj,
    model = fitted$model_obj, pred = pred_df,
    roc = roc_obj, roc_plot = plot_roc(roc_obj),
    method = 'svm', kernel = kernel, feature = feature_vec
  ), metrics)
}
