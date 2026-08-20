## Add a Chinese semantic summary to public roxygen blocks that lack one.
## This script prints an apply_patch document and does not edit source files.

args <- commandArgs(trailingOnly = TRUE)
package_dir <- if (length(args)) args[[1]] else '.'
target_file <- if (length(args) >= 2) args[[2]] else NULL
r_dir <- file.path(package_dir, 'R')

reference_lines <- readLines(
  file.path(package_dir, 'FUNCTION_REFERENCE.md'), encoding = 'UTF-8'
)
table_rows <- grep('^\\| `[^`]+\\(\\)` \\|', reference_lines, value = TRUE)
reference_name <- sub('^\\| `([^`]+)\\(\\)`.*$', '\\1', table_rows)
reference_desc <- sub('^\\| `[^`]+\\(\\)` \\| (.*) \\|$', '\\1', table_rows)
description_map <- stats::setNames(reference_desc, reference_name)

source_files <- list.files(r_dir, pattern = '[.]R$', full.names = TRUE)
if (!is.null(target_file)) {
  target_names <- strsplit(target_file, ',', fixed = TRUE)[[1]]
  source_files <- source_files[basename(source_files) %in% target_names]
  if (!length(source_files)) stop('Unknown R source file: ', target_file)
}

wrap_summary <- function(text) {
  strwrap(
    paste0('Chinese summary: ', text), width = 96,
    prefix = "#'   ", initial = "#' "
  )
}

contains_han <- function(text) {
  any(vapply(text, function(value) {
    code_points <- utf8ToInt(enc2utf8(value))
    any(code_points >= 0x4E00 & code_points <= 0x9FFF)
  }, logical(1)))
}

patch_lines <- '*** Begin Patch'
for (source_file in source_files) {
  source_lines <- readLines(source_file, warn = FALSE, encoding = 'UTF-8')
  definitions <- grep(
    '^[A-Za-z][A-Za-z0-9_.]*[[:space:]]*<-[[:space:]]*function',
    source_lines
  )
  if (!length(definitions)) next

  file_hunks <- character()
  for (line_no in definitions) {
    function_name <- sub(
      '^([A-Za-z][A-Za-z0-9_.]*)[[:space:]]*<-.*$', '\\1',
      source_lines[line_no]
    )
    if (startsWith(function_name, '.') ||
        !function_name %in% names(description_map)) next

    block_end <- line_no - 1L
    if (block_end < 1L || !grepl("^#'", source_lines[block_end])) next
    block_start <- block_end
    while (block_start > 1L && grepl("^#'", source_lines[block_start - 1L])) {
      block_start <- block_start - 1L
    }
    block <- source_lines[block_start:block_end]
    if (contains_han(block)) next

    first_tag <- grep("^#' @", block)[1]
    if (is.na(first_tag)) next
    insertion_at <- block_start + first_tag - 1L
    anchor_start <- max(block_start, insertion_at - 2L)
    old_context <- source_lines[anchor_start:insertion_at]
    new_context <- c(
      source_lines[anchor_start:(insertion_at - 1L)],
      "#'",
      wrap_summary(unname(description_map[[function_name]])),
      "#'",
      source_lines[insertion_at]
    )
    file_hunks <- c(
      file_hunks, '@@', paste0('-', old_context), paste0('+', new_context)
    )
  }
  if (length(file_hunks)) {
    patch_lines <- c(
      patch_lines,
      paste0('*** Update File: ', normalizePath(source_file, winslash = '\\')),
      file_hunks
    )
  }
}
writeLines(c(patch_lines, '*** End Patch'), useBytes = TRUE)
