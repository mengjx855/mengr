#### Jin-Xin Meng, jinxmeng@zju.edu.cn, 20260820, 20260916 ####

# 20260916: rename the public configuration function to `mengr_config()` and standardize documentation.

#### mengr_config ####

#' Configure paths used by mengR
#'
#'
#' 设置或读取 mengR 配置；目前主要管理公共数据库根目录。
#'
#' @param database Root directory containing shared reference databases.
#' @return The active mengR configuration, invisibly.
#' @export
mengr_config <- function(database = NULL) {
  if (!is.null(database)) {
    database <- normalizePath(database, winslash = "/", mustWork = FALSE)
    options(mengR.database = database)
  }

  config <- list(
    database = getOption(
      "mengR.database",
      Sys.getenv("MENGR_DATABASE", unset = "F:/database")
    )
  )

  invisible(config)
}

#### .mengr_db_file ####

.mengr_db_file <- function(...) {
  file.path(mengr_config()$database, ...)
}
