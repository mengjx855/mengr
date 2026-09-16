#### Jin-Xin Meng, jinxmeng@zju.edu.cn, 20220425, 20260916 ####

# 20260720 v0.1.2: add functions: calcu_correlation(), tidy_correlation()
# 20260916: standardize script metadata, function sections, documentation, and naming style.

#### calcu_correlation ####
# 计算相关性矩阵或 long-format 相关性结果
# x: 行为 sample，列为 variable 的数值矩阵或 data.frame
# y: 可选；如果提供，则计算 x variables 与 y variables 之间的相关性
# method: spearman / pearson / kendall
# adjust: p 值校正方法，默认 BH
# output: matrix / long / both
# remove_zero_var: 是否删除零方差 variable
# na_fill: 是否填充 NA；NULL 表示不填充

#' Calcu Correlation utility
#'
#'
#' 计算 feature 间 Pearson、Spearman 或 Kendall 相关及校正 P 值。
#'
#' @param x Primary vector or object supplied to the utility.
#' @param y Secondary vector or object supplied to the utility.
#' @param method Analysis or summary method; supported values are shown in the usage.
#' @param adjust Multiple-testing correction method for correlation P values.
#' @param output Requested output representation or output path.
#' @param remove_zero_var Whether to remove zero-variance features before analysis.
#' @param na_fill Value used to replace missing observations before analysis.
#' @param ... Additional arguments passed to `Hmisc::rcorr()` or `psych::corr.test()`.
#' @return A list of correlation/P-value matrices, a long edge table, or both, according to `output`.
#' @export
calcu_correlation <- function(
  x, y = NULL, method = c("spearman", "pearson", "kendall"),
  adjust = c(
    "BH", "holm", "hochberg", "hommel", "bonferroni",
    "BY", "fdr", "none"
  ),
  output = c("both", "long", "matrix"),
  remove_zero_var = TRUE, na_fill = NULL, ...
) {
  method <- match.arg(method)
  adjust <- match.arg(adjust)
  output <- match.arg(output)

  ## 1. 整理输入矩阵
  x <- data.frame(x, check.names = FALSE)
  x <- as.matrix(x)
  suppressWarnings(storage.mode(x) <- "numeric")

  if (is.null(colnames(x))) {
    colnames(x) <- paste0("X", seq_len(ncol(x)))
  }

  if (!is.null(y)) {
    y <- data.frame(y, check.names = FALSE)
    y <- as.matrix(y)
    suppressWarnings(storage.mode(y) <- "numeric")

    if (is.null(colnames(y))) {
      colnames(y) <- paste0("Y", seq_len(ncol(y)))
    }
  }

  x[is.infinite(x)] <- NA_real_
  if (!is.null(y)) y[is.infinite(y)] <- NA_real_

  if (!is.null(na_fill)) {
    x[is.na(x)] <- na_fill
    if (!is.null(y)) y[is.na(y)] <- na_fill
  }

  ## 2. 对齐样本
  if (!is.null(y)) {
    if (!is.null(rownames(x)) && !is.null(rownames(y))) {
      sample_use <- intersect(rownames(x), rownames(y))

      if (length(sample_use) < 3) {
        stop("Correlation analysis requires at least three matched samples.")
      }

      x <- x[sample_use, , drop = FALSE]
      y <- y[sample_use, , drop = FALSE]
    } else {
      if (nrow(x) != nrow(y)) {
        stop("x and y should have the same sample number.")
      }
    }
  }

  if (nrow(x) < 3) {
    stop("Correlation analysis requires at least three samples.")
  }

  ## 3. 删除零方差 variable
  if (isTRUE(remove_zero_var)) {
    keep_x <- apply(x, 2, \(v) {
      sum(is.finite(v)) >= 3 && stats::sd(v, na.rm = TRUE) > 0
    })

    x <- x[, keep_x, drop = FALSE]

    if (!is.null(y)) {
      keep_y <- apply(y, 2, \(v) {
        sum(is.finite(v)) >= 3 && stats::sd(v, na.rm = TRUE) > 0
      })

      y <- y[, keep_y, drop = FALSE]
    }
  }

  if (ncol(x) == 0) {
    stop("No valid variables remained in x.")
  }

  if (!is.null(y) && ncol(y) == 0) {
    stop("No valid variables remained in y.")
  }

  ## 4. 计算相关性
  ## 注意：这里固定 adjust = 'none'，后面自己统一计算 padj
  if (is.null(y)) {
    test <- psych::corr.test(
      x,
      method = method, adjust = "none",
      ci = FALSE, ...
    )
  } else {
    test <- psych::corr.test(
      x, y,
      method = method, adjust = "none",
      ci = FALSE, ...
    )
  }

  r_matrix <- test$r
  pval_matrix <- test$p

  ## 5. 统一计算 padj
  if (is.null(y)) {
    ## 单矩阵：只对非重复 pair 做 p.adjust
    p_vec <- pval_matrix[lower.tri(pval_matrix)]
    padj_vec <- stats::p.adjust(p_vec, method = adjust)

    padj_matrix <- matrix(
      NA_real_,
      nrow = nrow(pval_matrix),
      ncol = ncol(pval_matrix),
      dimnames = dimnames(pval_matrix)
    )

    padj_matrix[lower.tri(padj_matrix)] <- padj_vec
    padj_matrix[upper.tri(padj_matrix)] <- t(padj_matrix)[upper.tri(padj_matrix)]

    ## pval 也整理成对称矩阵
    pval_matrix[upper.tri(pval_matrix)] <- t(pval_matrix)[upper.tri(pval_matrix)]

    diag(pval_matrix) <- 0
    diag(padj_matrix) <- 0
  } else {
    ## 双矩阵：对所有 x-y pair 做 p.adjust
    padj_matrix <- matrix(
      stats::p.adjust(as.vector(pval_matrix), method = adjust),
      nrow = nrow(pval_matrix),
      ncol = ncol(pval_matrix),
      dimnames = dimnames(pval_matrix)
    )
  }

  ## 6. long-format 输出
  if (is.null(y)) {
    long <- data.frame(
      name_x = rownames(r_matrix)[row(r_matrix)[lower.tri(r_matrix)]],
      name_y = colnames(r_matrix)[col(r_matrix)[lower.tri(r_matrix)]],
      r = r_matrix[lower.tri(r_matrix)],
      pval = pval_matrix[lower.tri(pval_matrix)],
      padj = padj_matrix[lower.tri(padj_matrix)],
      check.names = FALSE
    )
  } else {
    long <- as.data.frame(as.table(r_matrix), stringsAsFactors = FALSE) |>
      dplyr::rename(name_x = Var1, name_y = Var2, r = Freq) |>
      dplyr::left_join(
        as.data.frame(as.table(pval_matrix), stringsAsFactors = FALSE) |>
          dplyr::rename(name_x = Var1, name_y = Var2, pval = Freq),
        by = c("name_x", "name_y")
      ) |>
      dplyr::left_join(
        as.data.frame(as.table(padj_matrix), stringsAsFactors = FALSE) |>
          dplyr::rename(name_x = Var1, name_y = Var2, padj = Freq),
        by = c("name_x", "name_y")
      )
  }

  long <- long |>
    dplyr::mutate(
      abs_r = abs(r),
      cor_type = dplyr::case_when(
        r > 0 ~ "positive",
        r < 0 ~ "negative",
        TRUE ~ "zero"
      )
    ) |>
    dplyr::arrange(padj, dplyr::desc(abs_r))

  ## 7. 输出
  result <- list(
    r = r_matrix,
    pval = pval_matrix,
    padj = padj_matrix,
    long = long,
    method = method,
    adjust = adjust
  )

  if (output == "matrix") {
    return(result[c("r", "pval", "padj")])
  }

  if (output == "long") {
    return(long)
  }

  return(result)
}

#### tidy_correlation ####
# 根据相关系数和 P 值筛选相关性矩阵
# r_mat: 相关系数矩阵，行为一组变量，列为另一组变量
# p_mat: 与 r_mat 对应的 P 值矩阵
# r:     相关系数绝对值阈值，默认 |r| >= 0
# pval:  P 值阈值，默认 p < 0.05
# only_signif:
#   FALSE = 根据显著相关性筛选行和列，并保留子矩阵中的全部原始 r 和 p
#   TRUE  = 仅保留显著相关性，不显著的 r 设为 0，p 设为 1
# 返回值:
# 一个 list，包括：
#   r_mat：筛选后的相关系数矩阵
#   p_mat：与 r_mat 对应的 P 值矩阵
#' Tidy Correlation utility
#'
#'
#' 将相关矩阵整理为 edge long table，并按相关性和显著性筛选。
#'
#' @param r_mat Matrix supplying r mat.
#' @param p_mat Matrix supplying p mat.
#' @param r Correlation-coefficient column or matrix to be converted to tidy form.
#' @param pval P-value column or matrix paired with the correlation coefficients.
#' @param only_signif Whether only statistically significant correlations are retained.
#' @return A filtered long data frame containing feature pairs, correlations, P values, and adjusted P values.
#' @export
tidy_correlation <- function(
  r_mat, p_mat, r = 0, pval = 0.05, only_signif = FALSE
) {
  r_mat <- as.matrix(r_mat)
  p_mat <- as.matrix(p_mat)

  # 检查是否存在行名和列名
  if (
    is.null(rownames(r_mat)) || is.null(colnames(r_mat)) ||
      is.null(rownames(p_mat)) || is.null(colnames(p_mat))
  ) {
    stop("r_mat and p_mat should both have row and column names.")
  }

  common_rows <- rownames(r_mat)[rownames(r_mat) %in% rownames(p_mat)]
  common_cols <- colnames(r_mat)[colnames(r_mat) %in% colnames(p_mat)]
  r_mat <- r_mat[common_rows, common_cols, drop = FALSE]
  p_mat <- p_mat[common_rows, common_cols, drop = FALSE]

  # 判断每个位置是否满足相关系数和 P 值阈值
  signif_mat <- (abs(r_mat) >= r & p_mat < pval & !is.na(r_mat) & !is.na(p_mat))

  # 保留至少存在一个显著相关性的行/列
  keep_rows <- rowSums(signif_mat, na.rm = TRUE) > 0
  keep_cols <- colSums(signif_mat, na.rm = TRUE) > 0

  # 先根据显著相关性筛选行和列
  r_mat <- r_mat[keep_rows, keep_cols, drop = FALSE]
  p_mat <- p_mat[keep_rows, keep_cols, drop = FALSE]

  signif_mat <- signif_mat[keep_rows, keep_cols, drop = FALSE]

  # only_signif = TRUE 时，仅保留显著位置
  if (only_signif) {
    r_mat[!signif_mat] <- 0
    p_mat[!signif_mat] <- 1
  }

  return(
    list(
      r_mat = r_mat,
      p_mat = p_mat
    )
  )
}

#### get_nwk_attr ####
#' Get Nwk Attr utility
#'
#'
#' 从相关 edge table 构建 igraph 网络，整理 edge/node 属性并可选导出。
#'
#' @param adjacency Adjacency matrix used to create or summarize a network.
#' @param suffix Suffix removed from or appended to derived identifiers.
#' @param mode Processing mode; supported values are shown in Usage.
#' @param weighted Whether to calculate weighted rather than unweighted UniFrac.
#' @param interaction Interaction class or direction retained in network results.
#' @param prefix_length Number of leading characters used to derive node prefixes.
#' @param export Optional export format; `NULL` keeps the operation in memory.
#' @param output_dir Existing directory used for optional file export.
#' @return A list containing the igraph object, edge and node attribute data frames, and paths of any exported files.
#' @export
get_nwk_attr <- function(
  adjacency, suffix = NULL,
  mode = c("undirected", "directed", "max", "min", "upper", "lower", "plus"),
  weighted = TRUE,
  interaction = c("all", "between_prefix"), prefix_length = 1,
  export = c("none", "cytoscape", "gephi", "both"),
  output_dir = "."
) {
  ## 1. 检查参数并建立网络
  mode <- match.arg(mode)
  interaction <- match.arg(interaction)
  export <- match.arg(export)
  adjacency_mat <- as.matrix(adjacency)
  if (nrow(adjacency_mat) != ncol(adjacency_mat)) {
    stop("adjacency should be a square matrix.")
  }
  graph_obj <- igraph::graph_from_adjacency_matrix(
    adjacency_mat,
    mode = mode,
    weighted = if (isTRUE(weighted)) "weight" else NULL,
    diag = FALSE
  )

  ## 2. 保留相关性的正负号，同时用绝对值作为网络权重
  if (isTRUE(weighted)) {
    igraph::E(graph_obj)$association <- igraph::E(graph_obj)$weight
    igraph::E(graph_obj)$weight <- abs(igraph::E(graph_obj)$weight)
  }
  edge_mat <- igraph::as_edgelist(graph_obj, names = TRUE)
  edge_df <- data.frame(
    source = edge_mat[, 1], target = edge_mat[, 2],
    weight = if (isTRUE(weighted)) igraph::E(graph_obj)$weight else 1,
    association = if (isTRUE(weighted)) {
      igraph::E(graph_obj)$association
    } else {
      1
    },
    check.names = FALSE
  )

  ## 3. 可选：仅保留不同前缀类型之间的连边
  if (interaction == "between_prefix") {
    source_prefix <- substr(edge_df$source, 1, prefix_length)
    target_prefix <- substr(edge_df$target, 1, prefix_length)
    edge_df <- edge_df[source_prefix != target_prefix, , drop = FALSE]
    graph_obj <- igraph::graph_from_data_frame(
      edge_df,
      directed = mode == "directed"
    )
  }
  node_df <- data.frame(
    node = igraph::V(graph_obj)$name,
    degree = igraph::degree(graph_obj),
    check.names = FALSE
  )

  ## 4. 默认不写文件；仅在明确指定 export 时导出
  file_list <- character()
  if (export != "none") {
    if (is.null(suffix) || !nzchar(suffix)) {
      stop("suffix is required when export is not none.")
    }
    if (!dir.exists(output_dir)) stop("output_dir does not exist: ", output_dir)
    edge_file <- file.path(output_dir, paste0(suffix, "-network-edge.tsv"))
    node_file <- file.path(output_dir, paste0(suffix, "-network-node.tsv"))
    utils::write.table(
      edge_df, edge_file,
      sep = "\t", row.names = FALSE, quote = FALSE
    )
    utils::write.table(
      node_df, node_file,
      sep = "\t", row.names = FALSE, quote = FALSE
    )
    file_list <- c(edge = edge_file, node = node_file)
    if (export %in% c("cytoscape", "both")) {
      gml_file <- file.path(output_dir, paste0(suffix, "-network.gml"))
      igraph::write_graph(graph_obj, gml_file, format = "gml")
      file_list <- c(file_list, cytoscape = gml_file)
    }
    if (export %in% c("gephi", "both")) {
      graphml_file <- file.path(output_dir, paste0(suffix, "-network.graphml"))
      igraph::write_graph(graph_obj, graphml_file, format = "graphml")
      file_list <- c(file_list, gephi = graphml_file)
    }
  }

  list(graph = graph_obj, edge_df = edge_df, node_df = node_df, files = file_list)
}

#### get_nwk_stat ####
#' Get Nwk Stat utility
#'
#'
#' 计算网络节点的 degree、strength、centrality 等拓扑统计量。
#'
#' @param graph An igraph object to summarize.
#' @param prefix_pattern Regular expression used to extract prefixes from node identifiers.
#' @return A node-level data frame of network topology statistics.
#' @export
get_nwk_stat <- function(
  graph, prefix_pattern = c(bacteria = "^b_", fungi = "^f_")
) {
  ## 1. 读取基础网络属性
  association_vec <- igraph::edge_attr(graph, "association")
  if (is.null(association_vec)) {
    association_vec <- igraph::edge_attr(graph, "weight")
  }
  if (is.null(association_vec)) {
    association_vec <- rep(NA_real_, igraph::ecount(graph))
  }
  node_vec <- igraph::V(graph)$name

  ## 2. 按给定前缀统计节点类型
  prefix_count <- vapply(
    prefix_pattern,
    \(pattern) sum(grepl(pattern, node_vec)),
    numeric(1)
  )

  ## 3. 整理为单行结果
  result_df <- data.frame(
    edge_n = igraph::ecount(graph),
    vertex_n = igraph::vcount(graph),
    positive_n = sum(association_vec > 0, na.rm = TRUE),
    negative_n = sum(association_vec < 0, na.rm = TRUE),
    connectance = igraph::edge_density(graph, loops = FALSE),
    average_degree = mean(igraph::degree(graph)),
    check.names = FALSE
  )
  for (name in names(prefix_count)) {
    result_df[[paste0(name, "_n")]] <- prefix_count[[name]]
  }
  result_df
}
