#### Jin-Xin Meng, jinxmeng@zju.edu.cn, 20260820, 20260925 ####

# 20260916: rename the public configuration function to `mengr_config()` and standardize documentation.
# 20260923: remove Chinese text from Roxygen documentation.
# 20260925: update package and option names from mengR to mengr.


#### mengr_config ####

#' Configure paths used by mengr
#'
#' @param database Root directory containing shared reference databases.
#' @return A named list containing the active mengr configuration, invisibly.
#' @export
mengr_config <- function(database = NULL) {
  if (!is.null(database)) {
    database <- normalizePath(database, winslash = "/", mustWork = FALSE)
    options(mengr.database = database)
  }

  config <- list(
    database = getOption(
      "mengr.database",
      Sys.getenv("MENGR_DATABASE", unset = "F:/database")
    )
  )

  invisible(config)
}

#### .mengr_db_file ####

.mengr_db_file <- function(...) {
  file.path(mengr_config()$database, ...)
}
