## Audit generated roxygen2 documentation without building or installing mengR.

package_dir <- commandArgs(trailingOnly = TRUE)
package_dir <- if (length(package_dir)) package_dir[[1]] else '.'

namespace_lines <- readLines(
  file.path(package_dir, 'NAMESPACE'), warn = FALSE, encoding = 'UTF-8'
)
exports <- sub('^export\\((.*)\\)$', '\\1', grep('^export\\(', namespace_lines, value = TRUE))

rd_files <- list.files(
  file.path(package_dir, 'man'), pattern = '[.]Rd$', full.names = TRUE
)
invisible(lapply(rd_files, tools::parse_Rd, encoding = 'UTF-8'))

aliases <- unique(unlist(lapply(rd_files, function(rd_file) {
  lines <- readLines(rd_file, warn = FALSE, encoding = 'UTF-8')
  sub('^\\\\alias\\{(.*)\\}$', '\\1', grep('^\\\\alias\\{', lines, value = TRUE))
})))

missing_docs <- setdiff(exports, aliases)
if (length(missing_docs)) {
  stop('Exported functions without an Rd alias: ', paste(missing_docs, collapse = ', '))
}

contains_han <- function(text) {
  any(vapply(text, function(value) {
    code_points <- utf8ToInt(enc2utf8(value))
    any(code_points >= 0x4E00 & code_points <= 0x9FFF)
  }, logical(1)))
}

rd_by_alias <- lapply(rd_files, function(rd_file) {
  lines <- readLines(rd_file, warn = FALSE, encoding = 'UTF-8')
  file_aliases <- sub(
    '^\\\\alias\\{(.*)\\}$', '\\1', grep('^\\\\alias\\{', lines, value = TRUE)
  )
  list(aliases = file_aliases, has_han = contains_han(lines))
})
missing_chinese <- exports[!vapply(exports, function(name) {
  any(vapply(rd_by_alias, function(item) {
    name %in% item$aliases && item$has_han
  }, logical(1)))
}, logical(1))]
if (length(missing_chinese)) {
  stop(
    'Exported functions without a Chinese summary: ',
    paste(missing_chinese, collapse = ', ')
  )
}

cat('Exports:', length(exports), '\n')
cat('Rd files:', length(rd_files), '\n')
cat('Missing Rd aliases:', length(missing_docs), '\n')
cat('Missing Chinese summaries:', length(missing_chinese), '\n')
cat('All Rd files parsed successfully.\n')
