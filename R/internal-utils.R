#### Jin-Xin Meng, 20260820, 20260820, v0.0.1 ####

.as_df <- function(x) {
  data.frame(x, check.names = FALSE)
}

.as_profile_df <- function(profile, numeric = FALSE) {
  profile_df <- data.frame(profile, check.names = FALSE)

  if (is.null(rownames(profile_df)) || is.null(colnames(profile_df))) {
    stop('profile should have feature row names and sample column names.')
  }

  if (isTRUE(numeric)) {
    profile_mat <- as.matrix(profile_df)
    suppressWarnings(storage.mode(profile_mat) <- 'numeric')
    profile_df <- data.frame(profile_mat, check.names = FALSE)
  }

  profile_df
}

.check_columns <- function(data, columns, object = 'data') {
  missing_cols <- setdiff(columns, colnames(data))
  if (length(missing_cols)) {
    stop(
      object, ' is missing required columns: ',
      paste(missing_cols, collapse = ', ')
    )
  }
  invisible(TRUE)
}

.match_distance_method <- function(method) {
  match.arg(
    method,
    c(
      'bray', 'jaccard', 'euclidean', 'manhattan', 'canberra',
      'kulczynski', 'gower', 'altGower', 'morisita', 'horn',
      'mountford', 'raup', 'binomial', 'chao', 'cao',
      'mahalanobis', 'unifrac'
    )
  )
}

.match_transform_method <- function(method) {
  if (is.null(method)) return(NULL)

  match.arg(
    method,
    c(
      'hellinger', 'total', 'max', 'frequency', 'normalize',
      'range', 'rank', 'rrank', 'standardize', 'pa',
      'chi.square', 'log', 'clr', 'rclr', 'alr'
    )
  )
}

.default_group_colors <- c(
  '#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3',
  '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3'
)

.resolve_group_colors <- function(group_level, group_color = NULL) {
  if (is.null(group_color)) {
    group_color <- rep(
      .default_group_colors,
      length.out = length(group_level)
    )
  } else if (!is.null(names(group_color))) {
    if (!all(group_level %in% names(group_color))) {
      stop('Named group_color should contain every value in group_level.')
    }
    group_color <- group_color[group_level]
  } else {
    group_color <- rep(group_color, length.out = length(group_level))
  }

  stats::setNames(group_color, group_level)
}

.align_profile_group <- function(
    profile, group, sample_col = 'sample', group_col = 'group',
    group_level = NULL, require_group = TRUE) {

  profile_df <- .as_profile_df(profile)
  group_df <- .as_df(group)

  required_cols <- sample_col
  if (isTRUE(require_group)) required_cols <- c(required_cols, group_col)
  .check_columns(group_df, required_cols, object = 'group')

  group_df[[sample_col]] <- as.character(group_df[[sample_col]])
  if (anyDuplicated(group_df[[sample_col]])) {
    stop('Duplicated sample identifiers found in group[[sample_col]].')
  }

  sample_vec <- intersect(colnames(profile_df), group_df[[sample_col]])
  if (!length(sample_vec)) {
    stop('No matched samples between profile and group.')
  }

  group_df <- group_df[
    match(sample_vec, group_df[[sample_col]]),
    , drop = FALSE
  ]
  profile_df <- profile_df[, sample_vec, drop = FALSE]

  if (isTRUE(require_group)) {
    if (is.null(group_level)) {
      group_level <- if (is.factor(group_df[[group_col]])) {
        levels(droplevels(group_df[[group_col]]))
      } else {
        unique(as.character(group_df[[group_col]]))
      }
    }
    group_df[[group_col]] <- factor(group_df[[group_col]], levels = group_level)
    if (anyNA(group_df[[group_col]])) {
      stop('group_level does not contain every observed group value.')
    }
  }

  list(
    profile_df = profile_df,
    group_df = group_df,
    sample_vec = sample_vec,
    group_level = group_level
  )
}

.profile_long_df <- function(
    profile, group, sample_col = 'sample', group_col = 'group',
    group_level = NULL, feature_col = 'name', value_col = 'value') {

  aligned <- .align_profile_group(
    profile = profile,
    group = group,
    sample_col = sample_col,
    group_col = group_col,
    group_level = group_level
  )

  group_key_df <- data.frame(
    sample = aligned$group_df[[sample_col]],
    group = aligned$group_df[[group_col]],
    check.names = FALSE
  )

  long_df <- aligned$profile_df |>
    tibble::rownames_to_column(feature_col) |>
    tidyr::pivot_longer(
      cols = -dplyr::all_of(feature_col),
      names_to = 'sample',
      values_to = value_col
    ) |>
    dplyr::left_join(group_key_df, by = 'sample')

  list(
    profile_df = aligned$profile_df,
    group_df = group_key_df,
    long_df = long_df,
    group_level = aligned$group_level
  )
}
