#### Jin-Xin Meng, jinxmeng@zju.edu.cn, 20241023, 20260916 ####

# 20260916: standardize script metadata, function sections, documentation, and naming style.

# 由于reactome没有对应的层级结构，获取某个ID的上级关系需要用此函数进行确定
#### reactome_longest_path ####

#' Reactome Longest Path utility
#'
#'
#' 根据本地 Reactome 父子关系表提取目标通路的最长层级路径。
#'
#' @param path_ID Reactome pathway identifier used as the starting node.
#' @param organism Reactome organism name used to restrict pathway records.
#' @param path_file Path to the Reactome parent-child relation file.
#' @param path_info_file Path to the Reactome pathway-information file.
#' @return A data frame describing the longest root-to-target Reactome path for each requested pathway.
#' @export
reactome_longest_path <- function(
  path_ID, organism = NULL,
  path_file = .mengr_db_file("Reactome", "ReactomePathwaysRelation.txt"),
  path_info_file = .mengr_db_file("Reactome", "ReactomePathways.txt")
) {
  org <- list(
    "BTA" = "Bos taurus", "GGA" = "Gallus gallus", "HSA" = "Homo sapiens",
    "MMU" = "Mus musculus", "SSC" = "Sus scrofa"
  )

  if (!file.exists(path_file)) stop("Reactome relation file not found: ", path_file)
  if (!file.exists(path_info_file)) stop("Reactome pathway file not found: ", path_info_file)

  organism <- ifelse(is.null(organism), unlist(strsplit(path_ID, "-"))[2], organism)
  if (!organism %in% names(org)) {
    stop("Unsupported organism code: ", organism)
  }
  path <- utils::read.delim(path_file, header = FALSE, col.names = c("from", "to")) |>
    dplyr::filter(grepl(organism, to)) |>
    igraph::graph_from_data_frame(directed = TRUE)
  path_info <- utils::read.delim(
    path_info_file,
    header = FALSE, col.names = c("id", "name", "tax")
  ) |>
    dplyr::filter(tax == org[[organism]])

  # 查找从所有节点出发的所有路径, 遍历每个节点作为起点
  all_paths <- purrr::map(igraph::V(path), \(x) igraph::all_simple_paths(path, from = x)) |>
    purrr::list_flatten()

  paths_with_target <- base::Filter(\(x) path_ID %in% names(igraph::V(path)[x]), all_paths)
  if (!length(paths_with_target)) {
    stop("path_ID not found in Reactome hierarchy: ", path_ID)
  }
  longest_path <- paths_with_target[[which.max(purrr::map_int(paths_with_target, length))]]
  out <- data.frame(id = names(igraph::V(path))[longest_path]) |>
    dplyr::mutate(name = path_info$name[match(id, path_info$id)])
  out <- tibble::add_column(out, seq_id = seq_len(nrow(out)), .before = 1)

  return(out)
}
