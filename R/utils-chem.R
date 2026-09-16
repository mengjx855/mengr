#### Jin-Xin Meng, jinxmeng@zju.edu.cn, 20260101, 20260916 ####

# 20260916: rename `calcu_DBE()` to `calcu_dbe()` and standardize documentation.

#### calcu_dbe ####

#' Calculate double-bond equivalents from molecular formulas
#'
#' 从分子式计算 double-bond equivalents，并正确处理卤素元素。
#'
#' @param formula Model formula defining the response and grouping variables.
#' @return A numeric vector of double-bond equivalents, preserving missing values for unparseable formulas.
#' @export
calcu_dbe <- function(formula) {
  ## 每个分子式独立解析；Cl、Br 等双字母元素不能按单个大写字母拆分
  vapply(formula, function(formula_one) {
    element_token <- stringr::str_extract_all(
      formula_one, "[A-Z][a-z]?[0-9]*"
    )[[1]]
    if (!length(element_token) || paste0(element_token, collapse = "") != formula_one) {
      stop("Invalid molecular formula: ", formula_one)
    }

    ## 元素后没有数字时按 1 处理；同一元素重复出现时合并计数
    element_vec <- stringr::str_extract(element_token, "^[A-Z][a-z]?")
    count_text <- stringr::str_extract(element_token, "[0-9]+$")
    count_vec <- ifelse(is.na(count_text), 1, as.numeric(count_text))
    atom_count <- rowsum(count_vec, group = element_vec, reorder = FALSE)[, 1]
    get_count <- function(element) {
      if (element %in% names(atom_count)) atom_count[[element]] else 0
    }

    ## 常用有机分子公式：DBE = C + 1 + (N - H - X) / 2
    halogen_n <- sum(vapply(c("F", "Cl", "Br", "I"), get_count, numeric(1)))
    get_count("C") + 1 +
      (get_count("N") - get_count("H") - halogen_n) / 2
  }, numeric(1), USE.NAMES = TRUE)
}
