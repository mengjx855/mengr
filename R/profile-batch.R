#### Jin-Xin Meng, jinxmeng@zju.edu.cn, 20221102, 20260923 ####

# 20260916: standardize script metadata, function sections, documentation, and naming style.
# 20260923: clarify metadata argument names and remove Chinese text from Roxygen documentation.


#### .prepare_batch_data ####

.prepare_batch_data <- function(
  profile, sample_meta, sample_col, batch_col,
  batch2_col = NULL, covariate_cols = NULL
) {
  metadata_df <- .as_df(sample_meta)
  required_cols <- unique(c(
    sample_col, batch_col, batch2_col, covariate_cols
  ))
  required_cols <- required_cols[!is.na(required_cols) & nzchar(required_cols)]
  .check_columns(metadata_df, required_cols, object = "sample_meta")

  aligned <- .align_profile_group(
    profile = profile,
    sample_meta = metadata_df,
    sample_col = sample_col,
    require_group = FALSE
  )
  metadata_df <- aligned$group_df
  metadata_df[[batch_col]] <- factor(metadata_df[[batch_col]])
  if (!is.null(batch2_col)) {
    metadata_df[[batch2_col]] <- factor(metadata_df[[batch2_col]])
  }

  design_mat <- NULL
  if (length(covariate_cols)) {
    formula_obj <- stats::reformulate(covariate_cols)
    design_mat <- stats::model.matrix(formula_obj, data = metadata_df)
  }

  list(
    profile_mat = as.matrix(aligned$profile_df),
    metadata_df = metadata_df,
    design_mat = design_mat
  )
}

#### remove_batch_combat ####

#' Remove batch effects with sva::ComBat
#'
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param sample_meta A sample metadata table containing sample identifiers,
#'   batch assignments, and optional covariates.
#' @param sample_col Name of the sample-identifier column.
#' @param batch_col Name of the `batch_col` input column.
#' @param covariate_cols Metadata columns whose effects should be retained during batch correction.
#' @param par_prior Whether parametric empirical-Bayes priors are used by ComBat.
#' @param prior_plots Whether ComBat diagnostic prior plots are produced.
#' @param mean_only Whether ComBat adjusts batch-specific means without adjusting variances.
#' @param ref_batch Optional reference-batch level passed to ComBat.
#' @param ... Additional arguments passed to `sva::ComBat()`.
#' @return A feature-by-sample numeric data frame after ComBat correction.
#' @export
remove_batch_combat <- function(
  profile, sample_meta, sample_col = "sample", batch_col = "batch",
  covariate_cols = NULL, par_prior = TRUE, prior_plots = FALSE,
  mean_only = FALSE, ref_batch = NULL, ...
) {
  metadata <- sample_meta
  prepared <- .prepare_batch_data(
    profile = profile,
    sample_meta = metadata,
    sample_col = sample_col,
    batch_col = batch_col,
    covariate_cols = covariate_cols
  )

  result_mat <- sva::ComBat(
    dat = prepared$profile_mat,
    batch = prepared$metadata_df[[batch_col]],
    mod = prepared$design_mat,
    par.prior = par_prior,
    prior.plots = prior_plots,
    mean.only = mean_only,
    ref.batch = ref_batch,
    ...
  )
  data.frame(result_mat, check.names = FALSE)
}

#### remove_batch_limma ####

#' Remove batch effects with limma::removeBatchEffect
#'
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param sample_meta A sample metadata table containing sample identifiers,
#'   one or two batch columns, and optional covariates.
#' @param sample_col Name of the sample-identifier column.
#' @param batch_col Name of the `batch_col` input column.
#' @param batch2_col Name of the `batch2_col` input column.
#' @param covariate_cols Metadata columns whose effects should be retained during batch correction.
#' @param ... Additional arguments passed to `limma::removeBatchEffect()`.
#' @return A feature-by-sample numeric data frame after limma batch-effect removal.
#' @export
remove_batch_limma <- function(
  profile, sample_meta, sample_col = "sample", batch_col = "batch",
  batch2_col = NULL, covariate_cols = NULL, ...
) {
  metadata <- sample_meta
  prepared <- .prepare_batch_data(
    profile = profile,
    sample_meta = metadata,
    sample_col = sample_col,
    batch_col = batch_col,
    batch2_col = batch2_col,
    covariate_cols = covariate_cols
  )
  batch2_vec <- if (is.null(batch2_col)) {
    NULL
  } else {
    prepared$metadata_df[[batch2_col]]
  }

  result_mat <- limma::removeBatchEffect(
    x = prepared$profile_mat,
    batch = prepared$metadata_df[[batch_col]],
    batch2 = batch2_vec,
    design = prepared$design_mat,
    ...
  )
  data.frame(result_mat, check.names = FALSE)
}
