#### Jin-Xin Meng, 20260501, 20260501, v0.0.1 ####

# 20260501: create this scripts to parse a variety of database.


#### get_HMDB2KEGG ####
#' Extract HMDB Xrefs utility
#'
#' `extract_HMDB_xrefs()` provides a reusable mengR workflow with input validation and
#'   standardized output.
#'
#' Chinese summary: 流式解析 HMDB XML，提取代谢物名称、标识符和外部数据库交叉引用。
#'
#' @param xml_file Path to the HMDB XML file.
#' @param ids Identifiers to query or retain.
#' @param relation Relationship table or relation type used to connect database identifiers.
#' @param output Requested output representation or output path.
#' @param remove_empty Logical control for `remove_empty`.
#' @param progresults Whether XML parsing progress is reported.
#' @param write_file Whether to enable the write file behavior.
#' @return A result object described in the Details section.
#' @export
extract_HMDB_xrefs <- function(
  xml_file, ids = c("hmdb_id", "name", "kegg_id"),
  relation = NULL, output = c("wide", "long", "relation"),
  remove_empty = TRUE, progresults = TRUE,
  write_file = NULL
) {
  output <- match.arg(output)

  if (!file.exists(xml_file)) {
    stop("xml_file not found: ", xml_file)
  }

  ## 常用字段
  common_xrefs <- c(
    "hmdb_id", "name", "kegg_id", "pubchem_compound_id",
    "chebi_id", "chemspider_id", "drugbank_id", "foodb_id",
    "biocyc_id", "bigg_id", "metlin_id", "knapsack_id",
    "phenol_explorer_compound_id", "vmh_id", "wikipedia_id"
  )

  common_chem <- c(
    "chemical_formula", "average_molecular_weight",
    "monisotopic_molecular_weight", "smiles",
    "inchi", "inchikey", "iupac_name",
    "traditional_iupac", "cas_registry_number"
  )

  if (length(ids) == 1 && ids == "common") {
    ids <- common_xrefs
  }

  if (length(ids) == 1 && ids == "chem") {
    ids <- c("hmdb_id", "name", common_chem)
  }

  if (length(ids) == 1 && ids == "all_common") {
    ids <- c(common_xrefs, common_chem)
  }

  ## 输出关系时，自动把 relation 里的字段加入 ids
  if (!is.null(relation)) {
    relation <- as.character(relation)
    ids <- unique(c("hmdb_id", "name", ids, relation))
    output <- "relation"
  }

  ## 输出列名 -> HMDB XML 字段名
  ## hmdb_id 在 XML 中叫 accession
  field_map <- c(
    hmdb_id = "accession",
    name = "name"
  )

  xml_fields <- ids
  xml_fields[ids %in% names(field_map)] <- field_map[ids[ids %in% names(field_map)]]
  names(xml_fields) <- ids

  message("Reading HMDB XML file...")
  x <- xml2::read_xml(xml_file)

  metabolites <- xml2::xml_find_all(x, ".//*[local-name()='metabolite']")
  n <- length(metabolites)

  message("Found ", n, " metabolite entries.")

  get_text <- function(node, field) {
    value <- xml2::xml_text(
      xml2::xml_find_first(node, paste0("./*[local-name()='", field, "']"))
    )

    if (length(value) == 0 || is.na(value) || trimws(value) == "") {
      return(NA_character_)
    }

    trimws(value)
  }

  result_list <- vector("list", n)

  if (isTRUE(progresults)) {
    pb <- utils::txtProgresultsBar(min = 0, max = n, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  for (i in seq_along(metabolites)) {
    m <- metabolites[[i]]

    result_list[[i]] <- as.data.frame(
      as.list(vapply(xml_fields, function(f) get_text(m, f), character(1))),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    if (isTRUE(progresults) && (i %% 100 == 0 || i == n)) {
      utils::setTxtProgresultsBar(pb, i)
    }
  }

  result <- dplyr::bind_rows(result_list)

  ## 去掉完全没有外部 ID 的行
  if (isTRUE(remove_empty)) {
    id_cols <- setdiff(colnames(result), c("hmdb_id", "name"))

    if (length(id_cols) > 0) {
      result <- result |>
        dplyr::filter(
          rowSums(!is.na(dplyr::across(dplyr::all_of(id_cols))) &
            dplyr::across(dplyr::all_of(id_cols)) != "") > 0
        )
    }
  }

  result <- dplyr::distinct(result)

  ## wide 格式
  if (output == "wide") {
    out <- result
  }

  ## long 格式：hmdb_id / name / database / db_id
  if (output == "long") {
    id_cols <- setdiff(colnames(result), c("hmdb_id", "name"))

    out <- result |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(id_cols),
        names_to = "database",
        values_to = "db_id"
      )

    if (isTRUE(remove_empty)) {
      out <- out |>
        dplyr::filter(!is.na(db_id), db_id != "")
    }

    out <- dplyr::distinct(out)
  }

  ## relation 格式：任意两个或多个字段之间的对应关系
  if (output == "relation") {
    if (is.null(relation)) {
      stop("`relation` is required when output = 'relation'.")
    }

    if (!all(relation %in% colnames(result))) {
      stop(
        "Some relation fields are not found in extracted resultult: ",
        paste(setdiff(relation, colnames(result)), collapse = ", ")
      )
    }

    out <- result |>
      dplyr::select(dplyr::all_of(relation))

    if (isTRUE(remove_empty)) {
      for (cc in relation) {
        out <- out |>
          dplyr::filter(!is.na(.data[[cc]]), .data[[cc]] != "")
      }
    }

    out <- dplyr::distinct(out)
  }

  if (!is.null(write_file)) {
    utils::write.table(
      out,
      file = write_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  return(out)
}
