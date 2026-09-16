#### Jin-Xin Meng, jinxmeng@zju.edu.cn, 20260820, 20260916 ####

# 20260916: standardize script metadata, function sections, documentation, and naming style.

#### .heatmap_annotation_palette ####

.heatmap_annotation_palette <- function(n) {
  if (n <= 8) {
    return(pald("Set2", n = n))
  }
  if (n <= 12) {
    return(pald("Set3", n = n))
  }
  if (n <= 20) {
    return(pald("Paired2", n = n))
  }
  scales::hue_pal()(n)
}

#### .prepare_heatmap_annotation ####

.prepare_heatmap_annotation <- function(
  annotation, target_names, show_cols, annotation_name
) {
  if (is.null(annotation)) {
    return(NULL)
  }

  annotation_df <- .as_df(annotation)
  .check_columns(annotation_df, "name", annotation_name)
  if (anyDuplicated(annotation_df$name)) {
    stop(annotation_name, "$name should not contain duplicated values.")
  }

  missing_names <- setdiff(target_names, annotation_df$name)
  if (length(missing_names) > 0) {
    stop(
      annotation_name, " is missing names: ",
      paste(utils::head(missing_names, 10), collapse = ", ")
    )
  }

  if (is.null(show_cols)) show_cols <- setdiff(colnames(annotation_df), "name")
  .check_columns(annotation_df, c("name", show_cols), annotation_name)

  annotation_df <- annotation_df[
    match(target_names, annotation_df$name), c("name", show_cols),
    drop = FALSE
  ]
  rownames(annotation_df) <- annotation_df$name
  annotation_df$name <- NULL
  annotation_df
}

#### plot_heatmap ####

#' Draw a heatmap with aligned row and column annotations
#'
#' Draw a `ComplexHeatmap::pheatmap()` heatmap from a numeric matrix. Annotation
#' rows are matched by the explicit `name` column, so their input order does not
#' need to match the profile. Missing categorical palettes are generated
#' automatically, while user-supplied palettes take precedence.
#'
#' 基于数值矩阵绘制 `ComplexHeatmap::pheatmap()` 热图。行、列注释表通过
#' `name` 列与 profile 对齐，因此输入顺序可以不同。分类变量缺失配色时会
#' 自动生成，用户在 `annotation_colors` 中提供的配色优先。
#'
#' @param profile A numeric matrix-like object with features in rows and samples
#'   in columns.
#' @param scale Scaling direction: `"none"`, `"row"`, or `"column"`.
#' @param border_color Cell-border color. Use `NA` to hide borders.
#' @param title Heatmap legend title.
#' @param cellwidth,cellheight Cell dimensions passed to
#'   `ComplexHeatmap::pheatmap()`.
#' @param cluster_rows,cluster_cols Whether to cluster rows or columns; an
#'   `hclust` object is also accepted by the underlying function.
#' @param treeheight_row,treeheight_col Dendrogram sizes.
#' @param show_rownames,show_colnames Whether to show matrix names.
#' @param row_annotation,col_annotation Optional data frames whose `name` column
#'   contains row or column names.
#' @param annotation_legend Whether to show annotation legends.
#' @param show_row_anno,show_col_anno Annotation columns to retain. `NULL` keeps
#'   all columns other than `name`.
#' @param annotation_colors Optional named list of annotation color vectors.
#' @param ... Additional arguments passed to `ComplexHeatmap::pheatmap()`.
#'
#' @return A ComplexHeatmap heatmap object, invisibly drawn according to the
#'   behavior of `ComplexHeatmap::pheatmap()`.
#' @export
plot_heatmap <- function(
  profile, scale = c("row", "none", "column"), border_color = NA,
  title = "Scaled value", cellwidth = 3, cellheight = 3,
  cluster_rows = TRUE, cluster_cols = TRUE, treeheight_row = 15,
  treeheight_col = 15, show_rownames = FALSE, show_colnames = FALSE,
  row_annotation = NULL, col_annotation = NULL, annotation_legend = TRUE,
  show_row_anno = NULL, show_col_anno = NULL, annotation_colors = NULL, ...
) {
  scale <- match.arg(scale)
  profile_mat <- as.matrix(profile)
  suppressWarnings(storage.mode(profile_mat) <- "numeric")
  if (anyNA(profile_mat)) {
    warning("profile contains NA or non-numeric values; they are shown as NA.")
  }
  if (is.null(rownames(profile_mat)) || is.null(colnames(profile_mat))) {
    stop("profile should have both row names and column names.")
  }

  # 1. 根据 name 列对齐注释，并仅保留指定字段
  row_annotation_df <- .prepare_heatmap_annotation(
    row_annotation, rownames(profile_mat), show_row_anno, "row_annotation"
  )
  col_annotation_df <- .prepare_heatmap_annotation(
    col_annotation, colnames(profile_mat), show_col_anno, "col_annotation"
  )

  # 2. 只为 character/factor 注释自动补充缺失的分类配色
  if (is.null(annotation_colors)) annotation_colors <- list()
  annotation_list <- Filter(Negate(is.null), list(row_annotation_df, col_annotation_df))
  for (annotation_df in annotation_list) {
    for (annotation_col in colnames(annotation_df)) {
      value_vec <- annotation_df[[annotation_col]]
      if (!is.numeric(value_vec) && is.null(annotation_colors[[annotation_col]])) {
        level_vec <- unique(as.character(value_vec[!is.na(value_vec)]))
        annotation_colors[[annotation_col]] <- stats::setNames(
          .heatmap_annotation_palette(length(level_vec)), level_vec
        )
      }
    }
  }
  if (length(annotation_colors) == 0) annotation_colors <- NA

  ComplexHeatmap::pheatmap(
    profile_mat,
    scale = scale,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    treeheight_row = treeheight_row,
    treeheight_col = treeheight_col,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    cellwidth = cellwidth,
    cellheight = cellheight,
    border_color = border_color,
    annotation_row = if (is.null(row_annotation_df)) NA else row_annotation_df,
    annotation_col = if (is.null(col_annotation_df)) NA else col_annotation_df,
    annotation_colors = annotation_colors,
    annotation_legend = annotation_legend,
    heatmap_legend_param = list(
      border = "black",
      title = title,
      title_gp = grid::gpar(fontface = "plain", fontsize = 10),
      title_position = "topleft",
      legend_direction = "vertical",
      legend_width = grid::unit(4, "cm"),
      labels_gp = grid::gpar(fontsize = 8)
    ),
    ...
  )
}
