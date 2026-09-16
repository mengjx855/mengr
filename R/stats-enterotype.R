#### Jin-Xin Meng, jinxmeng@zju.edu.cn, 20250617, 20260916 ####

# 20260916: standardize script metadata, function sections, documentation, and naming style.

#### calcu_jsd_dist ####

#' Calculate Jensen-Shannon distances between samples
#'
#' 对组成型 profile 计算 Jensen-Shannon divergence 距离。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param pseudocount Positive value used to replace or offset zeros before logarithmic operations.
#' @param normalize Whether to normalize each sample to unit sum before calculation.
#' @return A `dist` object containing pairwise Jensen-Shannon distances.
#' @export
calcu_jsd_dist <- function(profile, pseudocount = 1e-6, normalize = TRUE) {
  ## 1. 整理 feature × sample 矩阵
  profile_mat <- as.matrix(.as_profile_df(profile, numeric = TRUE))
  if (any(!is.finite(profile_mat)) || any(profile_mat < 0)) {
    stop("profile should contain finite non-negative values.")
  }
  if (!is.numeric(pseudocount) || length(pseudocount) != 1L ||
    pseudocount <= 0) {
    stop("pseudocount should be a single positive number.")
  }
  profile_mat[profile_mat == 0] <- pseudocount
  if (isTRUE(normalize)) {
    sample_sum <- colSums(profile_mat)
    if (any(sample_sum <= 0)) stop("Every sample should have a positive sum.")
    profile_mat <- sweep(profile_mat, 2, sample_sum, "/")
  }

  ## 2. 定义 KLD 和 JSD；局部函数不加点前缀
  kld <- function(x, y) sum(x * log(x / y))
  jsd <- function(x, y) {
    midpoint <- (x + y) / 2
    sqrt(0.5 * kld(x, midpoint) + 0.5 * kld(y, midpoint))
  }

  ## 3. 仅计算上三角，随后转为 dist 对象
  sample_n <- ncol(profile_mat)
  distance_mat <- matrix(
    0,
    nrow = sample_n, ncol = sample_n,
    dimnames = list(colnames(profile_mat), colnames(profile_mat))
  )
  if (sample_n > 1L) {
    for (i in seq_len(sample_n - 1L)) {
      for (j in seq.int(i + 1L, sample_n)) {
        distance_mat[i, j] <- jsd(profile_mat[, i], profile_mat[, j])
        distance_mat[j, i] <- distance_mat[i, j]
      }
    }
  }
  distance <- stats::as.dist(distance_mat)
  attr(distance, "method") <- "Jensen-Shannon"
  distance
}

#### pam_clustering ####

#' Partition samples around medoids
#'
#' 对距离对象执行 partitioning around medoids 聚类。
#'
#' @param distance A precomputed distance object; when supplied, it takes precedence over `profile`.
#' @param k Number of folds or clusters, according to the analysis performed.
#' @return A `pam` clustering object from `cluster::pam()`.
#' @export
pam_clustering <- function(distance, k) {
  ## cluster::pam() 在 diss = TRUE 时直接接收距离对象
  if (!is.numeric(k) || length(k) != 1L || k < 2) {
    stop("k should be a single integer greater than or equal to 2.")
  }
  as.vector(cluster::pam(stats::as.dist(distance), k = k, diss = TRUE)$clustering)
}
