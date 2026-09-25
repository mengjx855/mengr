#### Jin-Xin Meng, jinxmeng@zju.edu.cn, 20220918, 20260923 ####

# 20250206: fix some bug.
# 20250804: update function.
# 20251217: tidy function.
# 20260502: add profile_mpa() to generate metaphlan-like profile.
# 20260523: update taxa_trans(), plot_compos(), plot_compos_multiple(),
#           plot_compos_manual(), and plot_taxa_boxplot().
# 20260828: move plot-related functions to plot-taxa.R.
# 20260916: standardize script metadata, function sections, documentation, and naming style.
# 20260923: clarify metadata argument names and remove Chinese text from Roxygen documentation.


#### profile_mpa ####
# 生成 MetaPhlAn-like profile 表
# profile: 行为 feature，列为 sample 的丰度表
# taxa: feature 注释表，至少包含 feature_col，以及 d__/p__/c__/o__/f__/g__/s__ 等分类列
# feature_col: taxa 中用于匹配 profile 行名的列，默认 name
# mode:
#   all     输出所有分类层级路径
#   lowest  输出每个 feature 可获得的最低层级路径
#   species 只输出含 species 层级的路径
# normalize:
#   none     不转换
#   relative 转为相对丰度
#   percent  转为百分比

#' @title Build a MetaPhlAn-like profile
#'
#' @description
#' Build a Metaphlan-style profile, split the taxonomic 
#' hierarchy and output the profile of the specified classification level.
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param taxa Taxon names or taxonomy data used by the operation.
#' @param feature_col Name of the feature-identifier column.
#' @param mode Processing mode; supported values are shown in Usage.
#' @param normalize Whether to normalize each sample to unit sum before calculation.
#' @return A list of taxonomic profiles, one for each requested rank.
#' @export
profile_mpa <- function(
  profile, taxa, feature_col = "name", 
  mode = c("all", "lowest", "species"),
  normalize = c("none", "percent", "relative")
) {
  mode <- match.arg(mode)
  normalize <- match.arg(normalize)

  profile <- as.data.frame(profile, check.names = FALSE)
  taxa <- as.data.frame(taxa, check.names = FALSE)

  if (!feature_col %in% colnames(taxa)) {
    stop("tax should contain feature column: ", feature_col)
  }

  sample_cols <- colnames(profile)

  ## 1. 可选：把 profile 转成相对丰度或百分比
  if (normalize != "none") {
    if (normalize == "relative") {
      profile <- profile_trans_ra(profile, base = 1, digits = NULL)
    }

    if (normalize == "percent") {
      profile <- profile_trans_ra(profile, base = 100, digits = NULL)
    }
  }

  profile$.feature <- rownames(profile)

  ## 2. 解析 taxonomy
  ##    不依赖列名，而是根据 d__/p__/c__/o__/f__/g__/s__ 自动识别
  taxa_long <- tidyr::pivot_longer(
    taxa,
    cols = -dplyr::all_of(feature_col),
    names_to = "rank_col", values_to = "taxon"
  )

  taxa_long <- taxa_long |>
    dplyr::mutate(taxon = trimws(as.character(taxon))) |>
    dplyr::filter(!is.na(taxon), taxon != "", grepl("^[dkpcofgs]__", taxon)) |>
    dplyr::mutate(
      rank = substr(taxon, 1, 1),
      rank = dplyr::case_when(
        rank == "d" ~ "k",
        TRUE ~ rank
      ),
      taxon = sub("^[dkpcofgs]__", "", taxon),
      taxon = paste0(rank, "__", taxon)
    ) |>
    dplyr::select(.feature = dplyr::all_of(feature_col), rank, taxon) |>
    dplyr::distinct(.feature, rank, .keep_all = TRUE)

  taxa_wide <- tidyr::pivot_wider(
    taxa_long,
    names_from = rank, values_from = taxon
  )

  ranks <- c("k", "p", "c", "o", "f", "g", "s")

  for (r in ranks) {
    if (!r %in% colnames(taxa_wide)) {
      taxa_wide[[r]] <- NA_character_
    }
  }

  data <- dplyr::left_join(
    profile, taxa_wide,
    by = ".feature"
  )

  ## 3. 构建 clade_name
  make_path <- function(x) {
    x <- x[!is.na(x) & x != ""]
    if (length(x) == 0) {
      return(NA_character_)
    }
    paste(x, collapse = "|")
  }

  ## 4. 生成 MetaPhlAn-like 表
  if (mode == "all") {
    result_list <- list()

    for (i in seq_along(ranks)) {
      rr <- ranks[seq_len(i)]

      clade <- apply(data[, rr, drop = FALSE], 1, make_path)

      path_df <- data.frame(
        clade_name = clade,
        data[, sample_cols, drop = FALSE],
        check.names = FALSE
      )

      path_df <- path_df[!is.na(path_df$clade_name), , drop = FALSE]

      path_df <- path_df |>
        dplyr::group_by(clade_name) |>
        dplyr::summarise(
          dplyr::across(
            dplyr::all_of(sample_cols),
            ~ sum(.x, na.rm = TRUE)
          ),
          .groups = "drop"
        )

      result_list[[i]] <- path_df
    }

    result <- dplyr::bind_rows(result_list) |>
      dplyr::distinct(clade_name, .keep_all = TRUE)
  }

  if (mode == "lowest") {
    clade <- apply(data[, ranks, drop = FALSE], 1, make_path)

    result <- data.frame(
      clade_name = clade,
      data[, sample_cols, drop = FALSE],
      check.names = FALSE
    )

    result <- result[!is.na(result$clade_name), , drop = FALSE]

    result <- result |>
      dplyr::group_by(clade_name) |>
      dplyr::summarise(
        dplyr::across(
          dplyr::all_of(sample_cols),
          ~ sum(.x, na.rm = TRUE)
        ),
        .groups = "drop"
      )
  }

  if (mode == "species") {
    data <- data[!is.na(data$s), , drop = FALSE]

    clade <- apply(data[, ranks, drop = FALSE], 1, make_path)

    result <- data.frame(
      clade_name = clade,
      data[, sample_cols, drop = FALSE],
      check.names = FALSE
    )

    result <- result[!is.na(result$clade_name), , drop = FALSE]

    result <- result |>
      dplyr::group_by(clade_name) |>
      dplyr::summarise(
        dplyr::across(
          dplyr::all_of(sample_cols),
          ~ sum(.x, na.rm = TRUE)
        ),
        .groups = "drop"
      )
  }

  result <- tibble::column_to_rownames(result, "clade_name")

  return(result)
}

#### taxa_split ####
# taxonomy: data.frame containing name + taxonomy columns
# taxonomy_rename: optional, e.g. c(name = "OTU_ID", taxonomy = "Taxonomy")
#' Taxa Split utility
#'
#'
#'
#' @param taxonomy Feature taxonomy table used for annotation or aggregation.
#' @param taxonomy_rename Optional replacement name for the parsed taxonomy column.
#' @param sep Field separator used when reading or writing a text file.
#' @param from Source taxonomy rank or identifier type to be converted.
#' @param to Target taxonomy rank or identifier type produced by the conversion.
#' @param na_fill Value used to replace missing observations before analysis.
#' @param rm_suffix Whether suffix text is removed from split taxonomy labels.
#' @param ... Additional arguments passed to `stringr::str_split_fixed()`.
#' @return A data frame with parsed taxonomy ranks and the original taxonomy identifier.
#' @export
taxa_split <- function(
  taxonomy, taxonomy_rename = NULL, sep = ";", from = "d", to = "s",
  na_fill = "Unknown", rm_suffix = FALSE, ...
) {
  taxonomy <- data.frame(taxonomy, check.names = FALSE)

  ## 1. Optional rename
  if (!is.null(taxonomy_rename)) {
    if (is.null(names(taxonomy_rename))) {
      stop("`taxonomy_rename` should be a named vector, e.g. c(name = 'ID', taxonomy = 'Taxon').")
    }
    taxonomy <- dplyr::rename(taxonomy, dplyr::all_of(taxonomy_rename))
  }

  ## 2. Check required columns
  if (!all(c("name", "taxonomy") %in% colnames(taxonomy))) {
    stop("`taxonomy` must contain two columns: `name` and `taxonomy`.")
  }

  ## 3. Taxonomic rank definition
  taxa_info <- c(
    domain = "d", phylum = "p", class = "c", order = "o",
    family = "f", genus = "g", species = "s", strain = "t"
  )

  ## allow full rank names, e.g. from = "phylum"
  normalize_rank <- function(x) {
    x <- tolower(x)
    if (x %in% names(taxa_info)) {
      return(unname(taxa_info[x]))
    }
    if (x %in% taxa_info) {
      return(x)
    }
    if (x == "k") {
      return("d")
    } # compatible with k__Bacteria
    stop("Unknown taxonomic rank: ", x)
  }

  from <- normalize_rank(from)
  to <- normalize_rank(to)

  from_idx <- match(from, taxa_info)
  to_idx <- match(to, taxa_info)

  if (from_idx > to_idx) {
    stop("`from` should be an upstream rank of `to`.")
  }

  rank_idx <- from_idx:to_idx
  rank_names <- names(taxa_info)[rank_idx]
  rank_letters <- unname(taxa_info[rank_idx])
  rank_prefix <- paste0(rank_letters, "__")

  unknown_values <- c(
    "", "NA", "NaN", "NULL",
    "Unassigned", "Unknown", "Uncultured",
    "unassigned", "unknown", "uncultured"
  )

  ## 4. Remove GTDB-style suffix, e.g. Firmicutes_A -> Firmicutes
  remove_suffix <- function(x) {
    x <- sub("_[A-Z]$", "", x)
    x <- sub("_[A-Z](?=\\s)", "", x, perl = TRUE)
    x
  }

  ## 5. Clean one taxon
  clean_taxon <- function(x, prefix) {
    x <- trimws(as.character(x))
    ## remove existing rank prefix, e.g. p__Bacteroidota -> Bacteroidota
    x <- sub("^[a-zA-Z]__", "", x)
    x <- trimws(x)
    if (is.na(x) || x %in% unknown_values) x <- na_fill
    x <- paste0(prefix, x)
    if (isTRUE(rm_suffix)) x <- remove_suffix(x)
    x
  }

  ## 6. Split one taxonomy string
  split_one <- function(x) {
    out <- stats::setNames(paste0(rank_prefix, na_fill), rank_names)

    if (is.na(x) || trimws(x) %in% unknown_values) {
      return(out)
    }

    parts <- trimws(strsplit(as.character(x), split = sep, fixed = TRUE)[[1]])

    ## Case A: taxonomy has prefixes, e.g. d__Bacteria;p__Bacteroidota
    has_prefix <- grepl("^[a-zA-Z]__", parts)

    if (any(has_prefix)) {
      for (p in parts[has_prefix]) {
        r <- tolower(substr(p, 1, 1))
        if (r == "k") r <- "d"

        if (r %in% rank_letters) {
          idx <- match(r, rank_letters)
          out[idx] <- clean_taxon(p, rank_prefix[idx])
        }
      }
    } else {
      ## Case B: taxonomy has no prefixes, e.g. Bacteria;Bacteroidota;Bacteroidia
      full_rank_names <- names(taxa_info)
      full_rank_prefix <- paste0(unname(taxa_info), "__")

      full_out <- stats::setNames(
        paste0(full_rank_prefix, na_fill),
        full_rank_names
      )

      n <- min(length(parts), length(full_out))

      for (i in seq_len(n)) {
        full_out[i] <- clean_taxon(parts[i], full_rank_prefix[i])
      }

      out <- full_out[rank_names]
    }

    out
  }

  ## 7. Apply to all rows
  tax_mat <- t(vapply(
    taxonomy$taxonomy,
    split_one,
    FUN.VALUE = character(length(rank_names))
  ))

  result <- data.frame(
    name = taxonomy$name,
    tax_mat,
    check.names = FALSE,
    row.names = NULL
  )

  colnames(result) <- c("name", rank_names)

  return(result)
}

#### taxa_trans ####
# 按 taxonomy 注释将 profile 汇总到指定分类层级
# profile: 行为 feature，列为 sample 的丰度表
# taxonomy: feature 注释表
# feature_col: taxonomy 中用于匹配 profile 行名的列，默认 name
# taxa_col: taxonomy 中用于汇总的分类列，默认 family
# sample_col/group_col: group 中样本列和分组列
# top_n: 保留丰度最高的前 top_n - 1 个分类，其余合并为 Other
# top_list: 指定保留的分类名称，其余合并为 Other
# out_all: 是否输出全部分类；TRUE 时忽略 top_n
# remove_unknown: 是否在计算相对丰度前删除 unknown/unclassified 等分类
# trans_ra: 是否转换为百分比相对丰度
# collapse_group: 是否按 group 合并样本

#' Taxa Trans utility
#'
#'
#'
#' @param profile A feature-by-sample numeric matrix-like object.
#' @param taxonomy Feature taxonomy table used for annotation or aggregation.
#' @param sample_meta A sample metadata table containing sample and group columns.
#' @param feature_col Name of the feature-identifier column.
#' @param taxa_col Name of the `taxa_col` input column.
#' @param sample_col Name of the sample-identifier column.
#' @param group_col Name of the grouping column.
#' @param top_n Number of highest-ranking features or categories retained.
#' @param top_list Optional feature or category names retained regardless of abundance rank.
#' @param other_name Label assigned to features combined into the residual category.
#' @param out_all Whether all hierarchy levels or intermediate results are returned.
#' @param na_fill Value used to replace missing observations before analysis.
#' @param trans_ra Whether abundances are converted to relative abundance before analysis.
#' @param base Scaling constant used for relative-abundance output.
#' @param digits Optional number of decimal digits retained.
#' @param collapse_group Whether samples within each group are aggregated before plotting.
#' @param method Analysis or summary method; supported values are shown in the usage.
#' @param remove_unknown Whether to discard unknown or unclassified entries.
#' @param unknown_pattern Regular expression identifying unknown or unclassified annotations.
#' @return An aggregated taxonomic profile, or a list of profiles when multiple ranks are requested.
#' @export
taxa_trans <- function(
  profile, taxonomy, sample_meta = NULL, feature_col = "name", taxa_col = "family",
  sample_col = "sample", group_col = "group", top_n = 12, top_list = NULL,
  other_name = "Other", out_all = FALSE, na_fill = "Unknown",
  trans_ra = FALSE, base = 100, digits = 8, collapse_group = FALSE,
  method = c("mean", "median", "sum"), remove_unknown = FALSE,
  unknown_pattern = "unknown|unclassified|unassigned|uncultured"
) {
  group <- sample_meta
  if (is.function(method)) {
    summary_fun <- method
  } else {
    method <- match.arg(method)
    summary_fun <- match.fun(method)
  }
  profile <- data.frame(profile, check.names = FALSE)
  taxonomy <- data.frame(taxonomy, check.names = FALSE)

  if (!all(c(feature_col, taxa_col) %in% colnames(taxonomy))) {
    stop(
      "taxonomy should contain columns: ",
      feature_col, " | ", taxa_col
    )
  }

  if (isTRUE(out_all)) {
    top_n <- 0
  }

  ## 1. 匹配 taxonomy 注释
  taxa <- taxonomy[[taxa_col]][match(rownames(profile), taxonomy[[feature_col]])]
  taxa <- as.character(taxa)
  taxa[is.na(taxa) | taxa == ""] <- na_fill

  data <- profile |>
    dplyr::mutate(taxa = taxa) |>
    dplyr::filter(!is.na(taxa), taxa != "")

  ## 2. 删除 unknown / unclassified 等未明确注释的分类
  ## 默认匹配 unknown, unclassified, unassigned, uncultured
  ## 如果想自定义，可修改 unknown_pattern
  if (isTRUE(remove_unknown)) {
    data <- data |>
      dplyr::filter(
        !grepl(unknown_pattern, taxa, ignore.case = TRUE)
      )
  }

  if (nrow(data) == 0) {
    warning("No features remained after filtering taxonomy.")
    return(data.frame())
  }

  ## 3. 按分类层级汇总
  data <- data |>
    stats::aggregate(. ~ taxa, data = _, FUN = sum, na.rm = TRUE) |>
    tibble::column_to_rownames("taxa") |>
    data.frame(check.names = FALSE)

  ## 4. top_n 合并
  if (top_n > 0 && is.null(top_list)) {
    top_taxa <- data.frame(
      taxa = rownames(data),
      value = rowSums(data, na.rm = TRUE)
    ) |>
      dplyr::arrange(dplyr::desc(value)) |>
      utils::head(n = top_n - 1) |>
      dplyr::pull(taxa)

    data <- data |>
      tibble::rownames_to_column("name") |>
      dplyr::mutate(name = ifelse(name %in% top_taxa, name, other_name)) |>
      stats::aggregate(. ~ name, data = _, FUN = sum, na.rm = TRUE) |>
      tibble::column_to_rownames("name")
  }

  ## 5. top_list 合并
  if (!is.null(top_list)) {
    data <- data |>
      tibble::rownames_to_column("name") |>
      dplyr::mutate(name = ifelse(name %in% top_list, name, other_name)) |>
      stats::aggregate(. ~ name, data = _, FUN = sum, na.rm = TRUE) |>
      tibble::column_to_rownames("name")
  }

  ## 6. 是否合并到 group 水平
  if (isTRUE(collapse_group)) {
    if (is.null(group)) {
      stop("missing group file")
    }

    group <- data.frame(group, check.names = FALSE)

    if (!all(c(sample_col, group_col) %in% colnames(group))) {
      stop(
        "group should contain columns: ",
        sample_col, " | ", group_col
      )
    }

    data <- data.frame(t(data), check.names = FALSE)
    sample_vec <- rownames(data)
    data <- data |>
      dplyr::mutate(
        group = group[[group_col]][match(sample_vec, group[[sample_col]])]
      ) |>
      dplyr::filter(!is.na(group)) |>
      stats::aggregate(. ~ group, data = _, FUN = summary_fun) |>
      tibble::column_to_rownames("group") |>
      t() |>
      data.frame(check.names = FALSE)
  }

  ## 7. 转换为百分比相对丰度
  if (isTRUE(trans_ra)) {
    data <- profile_trans_ra(profile = data, base = base, digits = digits)
  }

  data <- data[
    rowSums(data, na.rm = TRUE) != 0,
    colSums(data, na.rm = TRUE) != 0,
    drop = FALSE
  ]

  return(data)
}
