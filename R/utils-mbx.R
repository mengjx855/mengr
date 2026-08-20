##### Jinxin Meng, 20241121, 20241121, v0.1 ####


#' Summarize chromatograms in an MChromatograms object
#'
#' Extract retention time and intensity vectors from every chromatogram and
#' report retention-time limits plus intensity quantiles and the retention time
#' at maximum intensity.
#'
#' 从 MChromatograms 对象的每条色谱中提取保留时间和强度，汇总保留时间
#' 范围、强度分位数和最高强度所在的保留时间。
#'
#' @param data An MChromatograms-like object supported by `xcms::rtime()` and
#'   `xcms::intensity()`.
#' @param label Retained for compatibility as a descriptive input label.
#'
#' @return A data frame with one row per chromatogram, sorted by the retention
#'   time of maximum intensity.
#' @export
parse_MChromatograms <- function(data, label = "MChromatograms object") {
  if (length(data) == 0) stop(label, ' is empty.')

  # 1. 将每条色谱整理为独立数据框
  chromatogram_list <- lapply(seq_along(data), \(index) {
    data.frame(
      rtime = xcms::rtime(data[index]),
      intensity = xcms::intensity(data[index])
    )
  })
  names(chromatogram_list) <- names(data)
  if (is.null(names(chromatogram_list))) {
    names(chromatogram_list) <- as.character(seq_along(data))
  }

  # 2. 提取保留时间范围、最强峰位置和强度分位数
  result_df <- purrr::imap_dfr(chromatogram_list, \(chromatogram_df, name) {
    if (nrow(chromatogram_df) == 0) return(NULL)
    intensity_quantile <- stats::quantile(
      chromatogram_df$intensity, na.rm = TRUE, names = FALSE
    )
    max_index <- which.max(chromatogram_df$intensity)
    data.frame(
      name = name,
      rt_min = min(chromatogram_df$rtime, na.rm = TRUE),
      rt_max = max(chromatogram_df$rtime, na.rm = TRUE),
      rt_intensity_max = chromatogram_df$rtime[max_index],
      min_intensity = intensity_quantile[1],
      q25_intensity = intensity_quantile[2],
      median_intensity = intensity_quantile[3],
      q75_intensity = intensity_quantile[4],
      max_intensity = intensity_quantile[5],
      check.names = FALSE
    )
  })
  dplyr::arrange(result_df, .data$rt_intensity_max)
}

#### MBT_eKEGG #### 
#' MBT eKEGG utility
#'
#' `MBT_eKEGG()` provides a reusable mengR workflow with input validation and standardized
#'   output.
#'
#' Chinese summary: 将代谢物列表与本地 eKEGG 注释关联并整理 pathway 结果。
#'
#' @param cpd_list Metabolite or compound identifiers to annotate.
#' @param database Numeric setting for `database`.
#' @return A result object described in the Details section.
#' @export
MBT_eKEGG <- function(
    cpd_list,
    database = .mengR_db_file(
      'KEGG', 'enrichment_analysis', 'cpd2path_enrichment.tsv'
    )) {
  if (!file.exists(database)) {
    stop('KEGG enrichment mapping file not found: ', database)
  }

  path_info <- utils::read.delim(database)
  
  eKEGG <- clusterProfiler::enricher(cpd_list, TERM2GENE = path_info, minGSSize = 1, 
                    pvalueCutoff = 1, qvalueCutoff = 1) |> 
    data.frame()
  
  return(eKEGG)
}
