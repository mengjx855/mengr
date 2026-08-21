#### Jin-Xin Meng, 20231028, 20260820, v0.2.0 ####

#' Run feature-wise Fisher exact tests
#'
#' Chinese summary: 对每个 feature 的 2×2 计数表执行 Fisher exact test。
#'
#' @param data An input data frame or compatible object.
#' @param feature_col Name of the feature-identifier column.
#' @param x_pos_col Name of the `x_pos_col` input column.
#' @param y_pos_col Name of the `y_pos_col` input column.
#' @param x_neg_col Name of the `x_neg_col` input column.
#' @param y_neg_col Name of the `y_neg_col` input column.
#' @param add_plab Logical control for `add_plab`.
#' @param add_padj Logical control for `add_padj`.
#' @param p_adjust_method Multiple-testing correction method passed to `stats::p.adjust()`.
#' @param add_enriched Logical control for `add_enriched`.
#' @param enriched_by Column used to identify which comparison group is enriched.
#' @param cutoff Numeric setting for `cutoff`.
#' @param x_name Display or identifier name for x name.
#' @param y_name Display or identifier name for y name.
#' @param ... Additional arguments passed to the underlying function.
#' @return A result object described in the Details section.
#' @export
calcu_fisher <- function(
  data, feature_col = "name", x_pos_col = "x_pos", y_pos_col = "y_pos",
  x_neg_col = "x_neg", y_neg_col = "y_neg",
  add_plab = FALSE, add_padj = TRUE,
  p_adjust_method = c(
    "BH", "holm", "hochberg", "hommel", "bonferroni",
    "BY", "fdr", "none"
  ),
  add_enriched = FALSE, enriched_by = c("padj", "pval"),
  cutoff = 0.05, x_name = "x", y_name = "y", ...
) {
  ## 1. 统一输入列并检查计数
  p_adjust_method <- match.arg(p_adjust_method)
  enriched_by <- match.arg(enriched_by)
  test_df <- .as_df(data)
  required_cols <- c(
    feature_col, x_pos_col, y_pos_col, x_neg_col, y_neg_col
  )
  .check_columns(test_df, required_cols, object = "data")
  test_df <- data.frame(
    name = as.character(test_df[[feature_col]]),
    x_pos = as.numeric(test_df[[x_pos_col]]),
    y_pos = as.numeric(test_df[[y_pos_col]]),
    x_neg = as.numeric(test_df[[x_neg_col]]),
    y_neg = as.numeric(test_df[[y_neg_col]]),
    check.names = FALSE
  )
  count_mat <- as.matrix(test_df[, -1, drop = FALSE])
  if (any(!is.finite(count_mat)) || any(count_mat < 0)) {
    stop("Count columns should contain finite non-negative values.")
  }

  ## 2. 逐个 feature 进行 Fisher exact test
  pval_vec <- vapply(seq_len(nrow(test_df)), function(row_idx) {
    contingency_mat <- matrix(
      c(
        test_df$x_pos[row_idx], test_df$y_pos[row_idx],
        test_df$x_neg[row_idx], test_df$y_neg[row_idx]
      ),
      nrow = 2, byrow = TRUE
    )
    stats::fisher.test(contingency_mat, ...)$p.value
  }, numeric(1))

  ## 3. 计算发生率、校正 P 值和富集方向
  result_df <- data.frame(
    name = test_df$name,
    x_occur = test_df$x_pos / (test_df$x_pos + test_df$x_neg),
    y_occur = test_df$y_pos / (test_df$y_pos + test_df$y_neg),
    pval = pval_vec,
    check.names = FALSE
  )
  result_df$x_occur[!is.finite(result_df$x_occur)] <- NA_real_
  result_df$y_occur[!is.finite(result_df$y_occur)] <- NA_real_
  if (isTRUE(add_padj) || enriched_by == "padj") {
    result_df$padj <- stats::p.adjust(
      result_df$pval,
      method = p_adjust_method
    )
  }
  if (isTRUE(add_plab)) {
    label_source <- if ("padj" %in% colnames(result_df)) {
      result_df$padj
    } else {
      result_df$pval
    }
    result_df$plab <- .add_plab(label_source)
  }
  if (isTRUE(add_enriched)) {
    signif_vec <- result_df[[enriched_by]] < cutoff
    result_df$enriched <- ifelse(
      signif_vec & result_df$x_occur > result_df$y_occur,
      x_name,
      ifelse(
        signif_vec & result_df$x_occur < result_df$y_occur,
        y_name, "none"
      )
    )
  }
  colnames(result_df)[colnames(result_df) == "x_occur"] <- paste0(x_name, "_occur")
  colnames(result_df)[colnames(result_df) == "y_occur"] <- paste0(y_name, "_occur")
  result_df
}
