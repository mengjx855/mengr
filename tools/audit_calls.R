args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1]] else "."
files <- list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE)

calls <- character()
definitions <- character()

for (file in files) {
  parsed <- parse(file, keep.source = TRUE)
  data <- getParseData(parsed)
  tokens <- data[
    data$terminal & data$token == "SYMBOL_FUNCTION_CALL",
    , drop = FALSE
  ]
  lines <- readLines(file, warn = FALSE, encoding = "UTF-8")

  if (nrow(tokens)) {
    for (i in seq_len(nrow(tokens))) {
      token <- tokens[i, ]
      before <- if (token$col1 > 1) {
        substr(lines[token$line1], 1, token$col1 - 1)
      } else {
        ""
      }
      if (!grepl("(::|:::)[[:space:]]*$", before)) {
        calls <- c(calls, token$text)
      }
    }
  }

  for (expr in parsed) {
    if (
      is.call(expr) && length(expr) >= 3 &&
      identical(expr[[1]], as.name("<-")) &&
      is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))
    ) {
      definitions <- c(definitions, as.character(expr[[2]]))
    }
  }
}

known <- unique(c(ls(baseenv(), all.names = TRUE), definitions))
unresolved <- sort(setdiff(unique(calls), known))
cat(paste(unresolved, collapse = "\n"), "\n", sep = "")
