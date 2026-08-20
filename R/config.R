#### Jin-Xin Meng, 20260820, 20260820, v0.0.1 ####

#' Configure paths used by mengR
#'
#'
#' Chinese summary: 设置或读取 mengR 配置；目前主要管理公共数据库根目录。
#'
#' @param database Root directory containing shared reference databases.
#' @return The active mengR configuration, invisibly.
#' @export
mengR_config <- function(database = NULL) {
  if (!is.null(database)) {
    database <- normalizePath(database, winslash = '/', mustWork = FALSE)
    options(mengR.database = database)
  }

  config <- list(
    database = getOption(
      'mengR.database',
      Sys.getenv('MENGR_DATABASE', unset = 'F:/database')
    )
  )

  invisible(config)
}

.mengR_db_file <- function(...) {
  file.path(mengR_config()$database, ...)
}
