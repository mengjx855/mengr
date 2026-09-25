#### Jin-Xin Meng, jinxmeng@zju.edu.cn, 20220529, 20260923 ####

# 20260916: standardize script metadata, function sections, documentation, and naming style.
# 20260923: remove Chinese text from Roxygen documentation.


#### .align_gene_length ####

.align_gene_length <- function(
  profile, gene_length, feature_col = "name", length_col = "length"
) {
  ## 按 profile 行名对齐 gene length，不能依赖输入顺序
  profile_df <- .as_profile_df(profile, numeric = TRUE)
  length_df <- .as_df(gene_length)
  .check_columns(length_df, c(feature_col, length_col), object = "gene_length")
  if (anyDuplicated(length_df[[feature_col]])) {
    stop("Duplicated feature identifiers found in gene_length.")
  }
  length_vec <- as.numeric(length_df[[length_col]][
    match(rownames(profile_df), length_df[[feature_col]])
  ])
  if (any(!is.finite(length_vec)) || any(length_vec <= 0)) {
    stop("Every profile feature should have a finite positive gene length.")
  }
  list(profile_mat = as.matrix(profile_df), length_vec = length_vec)
}

#### rc2tpm ####

#' Convert read counts to TPM
#'
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param gene_length Feature-length table used for TPM or FPKM/RPKM normalization.
#' @param feature_col Name of the feature-identifier column.
#' @param length_col Name of the `length_col` input column.
#' @return A feature-by-sample numeric matrix of transcripts per million.
#' @export
rc2tpm <- function(
  profile, gene_length, feature_col = "name", length_col = "length"
) {
  ## 1. 先计算 reads per kilobase
  aligned <- .align_gene_length(
    profile, gene_length, feature_col, length_col
  )
  rpk_mat <- aligned$profile_mat / (aligned$length_vec / 1000)

  ## 2. 每个样本标准化到一百万
  scale_vec <- colSums(rpk_mat, na.rm = TRUE) / 1e6
  if (any(scale_vec <= 0)) stop("Every sample should have a positive RPK sum.")
  data.frame(sweep(rpk_mat, 2, scale_vec, "/"), check.names = FALSE)
}

#### rc2fpkm ####

#' Convert fragment or read counts to FPKM/RPKM
#'
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param gene_length Feature-length table used for TPM or FPKM/RPKM normalization.
#' @param feature_col Name of the feature-identifier column.
#' @param length_col Name of the `length_col` input column.
#' @return A feature-by-sample numeric matrix of FPKM/RPKM values.
#' @export
rc2fpkm <- function(
  profile, gene_length, feature_col = "name", length_col = "length"
) {
  ## FPKM 与 RPKM 公式相同，区别在输入是 fragment counts 还是 read counts
  aligned <- .align_gene_length(
    profile, gene_length, feature_col, length_col
  )
  library_size <- colSums(aligned$profile_mat, na.rm = TRUE)
  if (any(library_size <= 0)) stop("Every sample should have a positive count sum.")
  result_mat <- aligned$profile_mat * 1e9 / aligned$length_vec
  result_mat <- sweep(result_mat, 2, library_size, "/")
  data.frame(result_mat, check.names = FALSE)
}

#### rc2rpm ####

#' Convert read counts to reads per million
#'
#' Convert a feature-by-sample count matrix to reads per million (RPM). When
#' `library_size` is `NULL`, column sums of `profile` are used. A two-column
#' sample table, a delimited text-file path, or a named numeric vector can
#' instead supply externally calculated library sizes.
#'
#'
#' @param profile A numeric feature-by-sample data frame or matrix.
#' @param library_size `NULL`, a two-column data frame, a path to a CSV/TSV
#'   file, or a named numeric vector. `NULL` uses `colSums(profile)`.
#'   `colSums(profile)`.
#' @param sample_col Sample identifier column in `library_size`.
#' @param library_size_col Library-size column in `library_size`.
#' @param library_size_sep Separator used when `library_size` is a file path.
#'
#' @return A numeric matrix with the same dimensions and dimnames as `profile`.
#' @export
rc2rpm <- function(
  profile, library_size = NULL, sample_col = "sample",
  library_size_col = "library_size", library_size_sep = NULL
) {
  ## 1. 整理 profile；默认用各列之和作为 library size
  profile_mat <- as.matrix(.as_profile_df(profile, numeric = TRUE))
  if (is.null(library_size)) {
    library_size_vec <- colSums(profile_mat, na.rm = TRUE)
  } else if (is.numeric(library_size) && !is.null(names(library_size))) {
    ## 2a. 带名称的数值向量按 sample name 对齐
    library_size_vec <- as.numeric(
      library_size[match(colnames(profile_mat), names(library_size))]
    )
  } else {
    ## 2b. 文件路径先读取为两列表；否则直接整理输入的数据框
    if (is.character(library_size) && length(library_size) == 1) {
      if (!file.exists(library_size)) {
        stop("library_size file was not found: ", library_size)
      }
      if (is.null(library_size_sep)) {
        library_size_sep <- if (grepl("[.]csv$", library_size, ignore.case = TRUE)) {
          ","
        } else {
          "\t"
        }
      }
      library_size_df <- utils::read.delim(
        library_size,
        sep = library_size_sep, check.names = FALSE
      )
    } else {
      library_size_df <- .as_df(library_size)
    }

    ## 2c. 两列表通过显式列名按 sample name 对齐
    .check_columns(
      library_size_df, c(sample_col, library_size_col),
      object = "library_size"
    )
    if (anyDuplicated(library_size_df[[sample_col]])) {
      stop("Duplicated sample identifiers found in library_size.")
    }
    library_size_vec <- as.numeric(library_size_df[[library_size_col]][
      match(colnames(profile_mat), library_size_df[[sample_col]])
    ])
  }

  ## 3. 检查对齐结果并转换为 RPM
  if (any(!is.finite(library_size_vec)) || any(library_size_vec <= 0)) {
    stop("Every profile sample should have a finite positive library size.")
  }
  names(library_size_vec) <- colnames(profile_mat)
  sweep(profile_mat, 2, library_size_vec / 1e6, "/")
}
