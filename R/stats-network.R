#### Jin-Xin Meng, jinxmeng@zju.edu.cn, 20240305, 20260916 ####

# 20260916: rename `calcu_MEN()` to `calcu_men()` and standardize documentation and naming.

#### calcu_men ####

#' Calculate and plot a microbial ecological network
#'
#' 使用 ggClusterNet 流程构建 microbial ecological network。
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param metadata A metadata or annotation data frame.
#' @param feature_col Name of the feature-identifier column.
#' @param sample_metadata Sample metadata required by normalization methods that model library size.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param p_threshold Maximum adjusted P value retained for a network edge.
#' @param r_threshold Minimum absolute correlation retained for a network edge.
#' @param scale Whether to scale values as described by the selected method.
#' @param scale_method Normalization method applied before network construction.
#' @param method Analysis or summary method; supported values are shown in the usage.
#' @param title Optional plot or result title.
#' @param n_hub Number of hub nodes to highlight; `FALSE` disables hub selection.
#' @param seed Optional random seed for reproducibility.
#' @return A list containing aligned inputs, correlations, graph and network tables, layout data, and the network plot.
#' @export
calcu_men <- function(
  profile, metadata, feature_col = NULL, sample_metadata = NULL,
  sample_col = "sample", group_col = "group",
  p_threshold = 0.05, r_threshold = 0.5,
  scale = FALSE,
  scale_method = c("rela", "sampling", "log", "TMM", "RLE", "upperquartile"),
  method = c("spearman", "pearson", "kendall"),
  title = "", n_hub = FALSE, seed = 2024
) {
  ## 1. 对齐 profile 与 feature taxonomy
  scale_method <- match.arg(scale_method)
  method <- match.arg(method)
  profile_df <- .as_profile_df(profile, numeric = TRUE)
  metadata_df <- .as_df(metadata)
  if (!is.null(feature_col)) {
    .check_columns(metadata_df, feature_col, object = "metadata")
    rownames(metadata_df) <- as.character(metadata_df[[feature_col]])
    metadata_df[[feature_col]] <- NULL
  }
  feature_vec <- intersect(rownames(profile_df), rownames(metadata_df))
  if (!length(feature_vec)) {
    stop("No matched features between profile and metadata.")
  }
  profile_df <- profile_df[feature_vec, , drop = FALSE]
  metadata_df <- metadata_df[feature_vec, , drop = FALSE]
  colnames(metadata_df) <- stringr::str_to_title(colnames(metadata_df))

  ## 2. 构建 phyloseq 对象；edgeR 类标准化需要样本分组
  otu_obj <- phyloseq::otu_table(profile_df, taxa_are_rows = TRUE)
  tax_obj <- phyloseq::tax_table(as.matrix(metadata_df))
  if (isTRUE(scale) && scale_method %in% c("TMM", "RLE", "upperquartile")) {
    if (is.null(sample_metadata)) {
      stop("sample_metadata is required for TMM, RLE or upperquartile scaling.")
    }
    aligned <- .align_profile_group(
      profile = profile_df, group = sample_metadata,
      sample_col = sample_col, group_col = group_col
    )
    profile_df <- aligned$profile_df
    otu_obj <- phyloseq::otu_table(profile_df, taxa_are_rows = TRUE)
    sample_df <- data.frame(
      Group = aligned$group_df[[group_col]],
      row.names = aligned$sample_vec, check.names = FALSE
    )
    ps_obj <- phyloseq::phyloseq(
      otu_obj, tax_obj, phyloseq::sample_data(sample_df)
    )
  } else {
    ps_obj <- phyloseq::phyloseq(otu_obj, tax_obj)
  }

  ## 3. 计算相关矩阵并转换为 igraph
  cor_obj <- ggClusterNet::cor_Big_micro2(
    ps = ps_obj, N = 0,
    p.threshold = p_threshold, r.threshold = r_threshold,
    scale = scale, method = method,
    met.scale = scale_method, p.adj = "BH"
  )
  correlation_mat <- cor_obj[[1]]
  edge_node_list <- ggClusterNet::nodeEdge(corr = correlation_mat)
  graph_obj <- igraph::graph_from_data_frame(
    edge_node_list[[1]],
    directed = FALSE, vertices = edge_node_list[[2]]
  )
  network_stat_df <- ggClusterNet::net_properties.2(
    graph_obj,
    n.hub = n_hub
  ) |>
    data.frame(check.names = FALSE) |>
    tibble::rownames_to_column("name")

  ## 4. 计算模块布局、节点属性和边属性
  layout_list <- ggClusterNet::model_igraph2(
    cor = correlation_mat, method = "cluster_fast_greedy", seed = seed
  )
  abundance_df <- t(ggClusterNet::vegan_otu(ps_obj)) |>
    data.frame(check.names = FALSE)
  taxonomy_df <- ggClusterNet::vegan_tax(ps_obj) |>
    data.frame(check.names = FALSE)
  node_df <- layout_list[[1]] |>
    ggClusterNet::nodeadd(
      otu_table = abundance_df, tax_table = taxonomy_df
    )
  model_df <- layout_list[[2]]
  node_color <- model_df |>
    dplyr::distinct(model, .keep_all = TRUE) |>
    dplyr::pull(color, name = model)

  edge_df_raw <- ggClusterNet::edgeBuild(
    cor = correlation_mat, node = layout_list[[1]]
  )
  edge_df <- dplyr::inner_join(
    model_df |>
      dplyr::select(OTU, model, color) |>
      dplyr::right_join(edge_df_raw, by = c("OTU" = "OTU_1")) |>
      dplyr::rename(OTU_1 = OTU, model1 = model, color1 = color),
    model_df |>
      dplyr::select(OTU, model, color) |>
      dplyr::right_join(edge_df_raw, by = c("OTU" = "OTU_2")) |>
      dplyr::rename(OTU_2 = OTU, model2 = model, color2 = color)
  ) |>
    dplyr::mutate(
      edge_class = ifelse(model1 == model2, model1, "across"),
      edge_color = ifelse(model1 == model2, color1, "#C1C1C1")
    )
  edge_color <- edge_df |>
    dplyr::distinct(edge_class, .keep_all = TRUE) |>
    dplyr::pull(edge_color, name = edge_class)

  ## 5. 绘制网络并返回计算过程中的可复用对象
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edge_df,
      ggplot2::aes(X1, Y1, xend = X2, yend = Y2, color = edge_class),
      linewidth = 1, show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = edge_color) +
    ggnewscale::new_scale_color() +
    ggplot2::geom_point(
      data = model_df, ggplot2::aes(X1, X2, color = model),
      size = 4, show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = node_color) +
    ggplot2::labs(title = title) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 10, hjust = 0.5),
      aspect.ratio = 1,
      plot.margin = grid::unit(c(1, 1, 1, 1), "cm")
    )

  list(
    profile_df = profile_df, metadata_df = metadata_df,
    correlation = cor_obj, graph = graph_obj,
    network_stat_df = network_stat_df, model_df = model_df,
    node_df = node_df, edge_df = edge_df, plot = p
  )
}
