#### Jin-Xin Meng, jinxmeng@zju.edu.cn, 20241121, 20260923 ####

# 20260916: rename exported functions to lowercase and standardize their documentation.
# 20260923: remove Chinese text from Roxygen documentation.



#### parse_mchromatograms ####

#' Summarize chromatograms in an MChromatograms object
#'
#' Extract retention time and intensity vectors from every chromatogram and
#' report retention-time limits plus intensity quantiles and the retention time
#' at maximum intensity.
#'
#'
#' @param data An MChromatograms-like object supported by `xcms::rtime()` and
#'   `xcms::intensity()`.
#' @param label Retained for compatibility as a descriptive input label.
#'
#' @return A data frame with one row per chromatogram, sorted by the retention
#'   time of maximum intensity.
#' @export
parse_mchromatograms <- function(data, label = "MChromatograms object") {
  if (length(data) == 0) stop(label, " is empty.")

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
    if (nrow(chromatogram_df) == 0) {
      return(NULL)
    }
    intensity_quantile <- stats::quantile(
      chromatogram_df$intensity,
      na.rm = TRUE, names = FALSE
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

#### mbt_ekegg ####
#' MBT eKEGG utility
#'
#'
#'
#' @param cpd_list Metabolite or compound identifiers to annotate.
#' @param database Path to the local annotation database file.
#' @return A data frame containing compound-set enrichment results.
#' @export
mbt_ekegg <- function(
  cpd_list,
  database = .mengr_db_file(
    "KEGG", "enrichment_analysis", "cpd2path_enrichment.tsv"
  )
) {
  if (!file.exists(database)) {
    stop("KEGG enrichment mapping file not found: ", database)
  }

  path_info <- utils::read.delim(database)

  eKEGG <- clusterProfiler::enricher(cpd_list,
    TERM2GENE = path_info, minGSSize = 1,
    pvalueCutoff = 1, qvalueCutoff = 1
  ) |>
    data.frame()

  return(eKEGG)
}
