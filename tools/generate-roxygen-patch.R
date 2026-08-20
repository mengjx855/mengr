## Generate an apply_patch document for missing/incomplete roxygen blocks.
## This script only prints a patch; it never edits package source files.

args <- commandArgs(trailingOnly = TRUE)
package_dir <- if (length(args)) args[[1]] else '.'
target_file <- if (length(args) >= 2) args[[2]] else NULL
r_dir <- file.path(package_dir, 'R')

## 1. Read Chinese summaries from FUNCTION_REFERENCE.md
reference_lines <- readLines(
  file.path(package_dir, 'FUNCTION_REFERENCE.md'), encoding = 'UTF-8'
)
table_rows <- grep('^\\| `[^`]+\\(\\)` \\|', reference_lines, value = TRUE)
reference_name <- sub('^\\| `([^`]+)\\(\\)`.*$', '\\1', table_rows)
reference_desc <- sub('^\\| `[^`]+\\(\\)` \\| (.*) \\|$', '\\1', table_rows)
description_map <- stats::setNames(reference_desc, reference_name)

## 2. Source functions so formal arguments can be inspected
all_source_files <- list.files(r_dir, pattern = '[.]R$', full.names = TRUE)
source_env <- new.env(parent = globalenv())
for (source_file in all_source_files) sys.source(source_file, envir = source_env)
source_files <- all_source_files
if (!is.null(target_file)) {
  target_names <- strsplit(target_file, ',', fixed = TRUE)[[1]]
  source_files <- source_files[basename(source_files) %in% target_names]
  if (!length(source_files)) stop('Unknown R source file: ', target_file)
}

roxygen_line <- function(text, initial = "#' ", subsequent = "#'   ") {
  strwrap(
    text, width = 96, prefix = subsequent, initial = initial
  )
}

param_description <- function(param) {
  exact <- c(
    profile = 'A feature-by-sample numeric matrix-like object.',
    data = 'An input data frame or compatible object.',
    group = 'A sample metadata table containing sample and group columns.',
    metadata = 'A metadata or annotation data frame.',
    sample_col = 'Name of the sample-identifier column.',
    group_col = 'Name of the grouping column.',
    feature_col = 'Name of the feature-identifier column.',
    group_level = 'Optional order of group levels.',
    group_color = 'Optional colors aligned to `group_level`.',
    group_colors = 'Optional named colors for groups.',
    dist_method = 'Distance method; available values are validated with `match.arg()`.',
    transform = 'Optional transformation applied before analysis.',
    method = 'Analysis or summary method; supported values are shown in the usage.',
    title = 'Optional plot or result title.',
    subtitle = 'Optional plot subtitle.',
    xlab = 'Optional x-axis label.',
    ylab = 'Optional y-axis label.',
    seed = 'Optional random seed for reproducibility.',
    digits = 'Optional number of decimal digits retained.',
    quiet = 'Whether to suppress progress messages.',
    verbose = 'Whether to print progress or diagnostic messages.',
    `...` = 'Additional arguments passed to the underlying function.'
  )
  if (param %in% names(exact)) return(unname(exact[[param]]))
  if (grepl('_col$', param)) {
    return(paste0('Name of the `', param, '` input column.'))
  }
  if (grepl('_level$', param)) {
    return(paste0('Optional order for `', param, '`.'))
  }
  if (grepl('color|colour', param)) {
    return(paste0('Color specification for `', param, '`.'))
  }
  if (grepl('show_|add_|remove_|plot_|simplify|overwrite|append|reverse|scale$', param)) {
    return(paste0('Logical control for `', param, '`.'))
  }
  if (grepl('size|width|height|alpha|linewidth|level|cutoff|threshold|base$', param)) {
    return(paste0('Numeric setting for `', param, '`.'))
  }
  paste0('Input controlling `', param, '`.')
}

return_description <- function(name) {
  if (startsWith(name, 'plot_')) {
    return('A plot object; analysis data or models may also be stored as attributes.')
  }
  if (startsWith(name, 'theme_')) {
    return('A ggplot2 theme object.')
  }
  if (grepl('_(show)$', name)) {
    return('A palette preview or a character summary of available palettes.')
  }
  if (name %in% c('cellchat_run')) {
    return('A processed model object.')
  }
  'A result object described in the Details section.'
}

english_title <- function(name) {
  words <- gsub('_', ' ', name, fixed = TRUE)
  paste0(tools::toTitleCase(words), ' utility')
}

## 3. Create patch hunks without changing files directly
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
        !exists(function_name, envir = source_env, inherits = FALSE)) next

    function_obj <- get(function_name, envir = source_env)
    formal_names <- names(formals(function_obj))

    previous_line <- line_no - 1L
    while (previous_line > 0 && !nzchar(trimws(source_lines[previous_line]))) {
      previous_line <- previous_line - 1L
    }
    has_roxygen <- previous_line > 0 &&
      startsWith(trimws(source_lines[previous_line]), "#'")

    additions <- character()
    if (!has_roxygen) {
      chinese_desc <- description_map[[function_name]]
      if (is.null(chinese_desc) || is.na(chinese_desc)) {
        chinese_desc <- paste0('See `', function_name, '()` for this operation.')
      }
      additions <- c(
        roxygen_line(english_title(function_name)),
        "#'",
        roxygen_line(paste0(
          '`', function_name, '()` provides a reusable mengR workflow with ',
          'input validation and standardized output.'
        )),
        "#'",
        roxygen_line(paste0('Chinese summary: ', chinese_desc)),
        "#'"
      )
    } else {
      block_start <- previous_line
      while (block_start > 1 &&
             startsWith(trimws(source_lines[block_start - 1L]), "#'")) {
        block_start <- block_start - 1L
      }
      block_lines <- source_lines[block_start:previous_line]
    }

    if (!has_roxygen) block_lines <- character()
    param_lines <- grep('@param[[:space:]]+', block_lines, value = TRUE)
    documented_params <- character()
    if (length(param_lines)) {
      param_specs <- sub(
        '^.*@param[[:space:]]+([^[:space:]]+).*$','\\1', param_lines
      )
      documented_params <- unlist(
        strsplit(param_specs, ',', fixed = TRUE), use.names = FALSE
      )
    }
    for (param in formal_names) {
      if (!param %in% documented_params) {
        additions <- c(additions, roxygen_line(
          param_description(param),
          initial = paste0("#' @param ", param, " "),
          subsequent = "#'   "
        ))
      }
    }
    if (!any(grepl("@return", block_lines, fixed = TRUE))) {
      additions <- c(additions, roxygen_line(
        return_description(function_name), initial = "#' @return ",
        subsequent = "#'   "
      ))
    }
    if (!any(grepl("@export", block_lines, fixed = TRUE))) {
      additions <- c(additions, "#' @export")
    }
    if (!length(additions)) next

    replacement <- c(additions, source_lines[line_no])
    file_hunks <- c(
      file_hunks,
      '@@',
      paste0('-', source_lines[line_no]),
      paste0('+', replacement)
    )
  }

  if (length(file_hunks)) {
    normalized_file <- normalizePath(source_file, winslash = '\\', mustWork = TRUE)
    patch_lines <- c(
      patch_lines,
      paste0('*** Update File: ', normalized_file),
      file_hunks
    )
  }
}
patch_lines <- c(patch_lines, '*** End Patch')
writeLines(patch_lines, useBytes = TRUE)
