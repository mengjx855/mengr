#### Jin-Xin Meng, jinxmeng@zju.edu.cn, 20241023, 20260916 ####

# 20260916: standardize script metadata, function sections, documentation, and naming style.

#### .lasso_outcome ####

.lasso_outcome <- function(group_vec, family, positive_class = NULL) {
  family <- match.arg(family, c("binomial", "gaussian", "poisson"))

  if (family == "binomial") {
    group_vec <- factor(group_vec)
    if (nlevels(group_vec) != 2L) {
      stop("family = 'binomial' requires exactly two groups.")
    }
    if (is.null(positive_class)) positive_class <- levels(group_vec)[2]
    if (!positive_class %in% levels(group_vec)) {
      stop("positive_class is not present in group.")
    }
    negative_class <- setdiff(levels(group_vec), positive_class)
    group_vec <- factor(group_vec, levels = c(negative_class, positive_class))
  } else {
    group_vec <- suppressWarnings(as.numeric(as.character(group_vec)))
    if (any(!is.finite(group_vec))) {
      stop("group should be numeric for gaussian or poisson models.")
    }
  }

  list(y = group_vec, family = family, positive_class = positive_class)
}

#### lasso_kfold ####

#' Cross-validated LASSO predictions
#'
#' 使用 glmnet 和分层 folds 执行 LASSO k-fold cross-validation。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param group A sample metadata table containing sample and group columns.
#' @param k Number of folds or clusters, according to the analysis performed.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param family Model family; supported values are shown in Usage.
#' @param positive_class Outcome level treated as the positive class for binary metrics.
#' @param seed Optional random seed for reproducibility.
#' @param inner_folds Number of internal folds used to tune the regularization parameter.
#' @param ... Additional arguments passed to `glmnet::cv.glmnet()`.
#' @return A data frame with sample IDs, observed outcomes, predictions, and predicted classes for binomial models.
#' @export
lasso_kfold <- function(
  profile, group, k = 5, sample_col = "sample", group_col = "group",
  family = c("binomial", "gaussian", "poisson"),
  positive_class = NULL, seed = 2024, inner_folds = 5, ...
) {
  family <- match.arg(family)
  aligned <- .align_profile_group(
    profile = profile, group = group,
    sample_col = sample_col, group_col = group_col
  )
  x_mat <- t(as.matrix(.as_profile_df(aligned$profile_df, numeric = TRUE)))
  outcome <- .lasso_outcome(
    aligned$group_df[[group_col]], family, positive_class
  )
  y_vec <- outcome$y

  if (k < 2L || k > nrow(x_mat)) {
    stop("k should be between 2 and the number of matched samples.")
  }
  set.seed(seed)
  fold_list <- caret::createFolds(y_vec, k = k, returnTrain = FALSE)
  pred_vec <- rep(NA_real_, nrow(x_mat))

  for (fold_idx in seq_along(fold_list)) {
    test_idx <- fold_list[[fold_idx]]
    train_idx <- setdiff(seq_len(nrow(x_mat)), test_idx)
    nfolds <- min(inner_folds, length(train_idx))
    if (nfolds < 3L) stop("Each training split requires at least three samples.")

    set.seed(seed + fold_idx - 1L)
    cv_fit <- glmnet::cv.glmnet(
      x = x_mat[train_idx, , drop = FALSE],
      y = y_vec[train_idx], family = family,
      alpha = 1, nfolds = nfolds, ...
    )
    pred_vec[test_idx] <- drop(stats::predict(
      cv_fit,
      newx = x_mat[test_idx, , drop = FALSE],
      s = "lambda.min", type = "response"
    ))
  }

  result_df <- data.frame(
    sample = aligned$sample_vec,
    actual = if (is.factor(y_vec)) as.character(y_vec) else y_vec,
    pred = pred_vec,
    check.names = FALSE
  )
  if (family == "binomial") {
    negative_class <- setdiff(levels(y_vec), outcome$positive_class)
    result_df$predicted <- ifelse(
      result_df$pred >= 0.5, outcome$positive_class, negative_class
    )
  }
  result_df
}

#### lasso_next_validate ####

#' Train a LASSO model in one dataset and validate it in another
#'
#' 在独立数据集上应用已训练的 LASSO 模型并评估预测。
#'
#' @param profile_x Feature-by-sample profile used for model training or the first data space.
#' @param profile_y Feature-by-sample profile used for validation or the second data space.
#' @param group_x Sample metadata associated with `profile_x`.
#' @param group_y Sample metadata associated with `profile_y`.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param family Model family; supported values are shown in Usage.
#' @param positive_class Outcome level treated as the positive class for binary metrics.
#' @param seed Optional random seed for reproducibility.
#' @param inner_folds Number of internal folds used to tune the regularization parameter.
#' @param ... Additional arguments passed to `glmnet::cv.glmnet()`.
#' @return A list containing the fitted cross-validated model, shared features, validation predictions, and binary metrics when applicable.
#' @export
lasso_next_validate <- function(
  profile_x, profile_y, group_x, group_y,
  sample_col = "sample", group_col = "group",
  family = c("binomial", "gaussian", "poisson"),
  positive_class = NULL, seed = 2024, inner_folds = 5, ...
) {
  family <- match.arg(family)
  aligned_x <- .align_profile_group(
    profile_x, group_x,
    sample_col = sample_col, group_col = group_col
  )
  aligned_y <- .align_profile_group(
    profile_y, group_y,
    sample_col = sample_col, group_col = group_col
  )
  feature_vec <- intersect(
    rownames(aligned_x$profile_df), rownames(aligned_y$profile_df)
  )
  if (!length(feature_vec)) {
    stop("No matched features between profile_x and profile_y.")
  }

  train_x <- t(as.matrix(.as_profile_df(
    aligned_x$profile_df[feature_vec, , drop = FALSE],
    numeric = TRUE
  )))
  test_x <- t(as.matrix(.as_profile_df(
    aligned_y$profile_df[feature_vec, , drop = FALSE],
    numeric = TRUE
  )))
  train_outcome <- .lasso_outcome(
    aligned_x$group_df[[group_col]], family, positive_class
  )
  test_outcome <- .lasso_outcome(
    aligned_y$group_df[[group_col]], family, train_outcome$positive_class
  )
  if (family == "binomial" && !identical(levels(train_outcome$y), levels(test_outcome$y))) {
    stop("Training and validation datasets should contain the same two groups.")
  }

  nfolds <- min(inner_folds, nrow(train_x))
  if (nfolds < 3L) stop("Training requires at least three matched samples.")
  set.seed(seed)
  cv_fit <- glmnet::cv.glmnet(
    x = train_x, y = train_outcome$y,
    family = family, alpha = 1, nfolds = nfolds, ...
  )
  model_obj <- glmnet::glmnet(
    x = train_x, y = train_outcome$y,
    family = family, alpha = 1, lambda = cv_fit$lambda.min
  )
  pred_vec <- drop(stats::predict(
    model_obj,
    newx = test_x, type = "response"
  ))
  result_df <- data.frame(
    sample = aligned_y$sample_vec,
    actual = if (is.factor(test_outcome$y)) {
      as.character(test_outcome$y)
    } else {
      test_outcome$y
    },
    pred = pred_vec,
    check.names = FALSE
  )

  out <- list(
    model = model_obj, cv_fit = cv_fit,
    pred = result_df, method = "lasso",
    feature = feature_vec
  )
  if (family == "binomial") {
    positive_class <- train_outcome$positive_class
    negative_class <- setdiff(levels(train_outcome$y), positive_class)
    result_df$predicted <- ifelse(
      result_df$pred >= 0.5, positive_class, negative_class
    )
    confusion_mat <- table(
      predicted = factor(result_df$predicted, levels = levels(train_outcome$y)),
      actual = factor(result_df$actual, levels = levels(train_outcome$y))
    )
    total_n <- sum(confusion_mat)
    out$pred <- result_df
    out$confusion_matrix <- confusion_mat
    out$accuracy <- sum(diag(confusion_mat)) / total_n
    out$sensitivity <- confusion_mat[positive_class, positive_class] /
      sum(confusion_mat[, positive_class])
    out$specificity <- confusion_mat[negative_class, negative_class] /
      sum(confusion_mat[, negative_class])
  }
  out
}
