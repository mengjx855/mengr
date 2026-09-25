#### Jin-Xin Meng, jinxmeng@zju.edu.cn, 20260911, 20260923 ####

# 20260911: create new file for ggplot2-related utility functions.
# 20260911: add a function `add_strip` to add alternating strip backgrounds to ggplot2 plots.
# 20260911: add a function `format_stat_label` to format p-values and q-values for ggplot2 plots.
# 20260916: validate inputs, document both functions, and prevent `add_strip()` from replacing existing fill scales.
# 20260923: correct Roxygen export tags.


#### format_stat_label ####
#' Format statistical values as plotmath labels
#'
#' Format p-values or adjusted p-values as character strings that can be parsed
#' by ggplot2 when `parse = TRUE`. Values below `10^-digits` are displayed as a
#' threshold unless scientific notation is requested. Missing values remain
#' `NA_character_`.
#'
#' @param p Numeric vector of statistical values. Non-missing values must lie
#'   between 0 and 1.
#' @param digits Non-negative integer giving the number of decimal places or
#'   significant digits used for scientific notation.
#' @param scientific Whether values below the display threshold should use
#'   scientific notation instead of a less-than label.
#' @param label Label shown before the value. The value `"p"` is italicized;
#'   other supported labels are displayed literally.
#'
#' @return A character vector with the same length as `p`, suitable for
#'   `ggplot2::geom_text(parse = TRUE)` or `ggplot2::annotate(parse = TRUE)`.
#' @export
format_stat_label <- function(
    p, digits = 3, scientific = FALSE,
    label = c(
      "p", "q", "pvalue", "qvalue", "padjust",
      "p.value", "p.adjust", "FDR"
    )
) {
  label <- match.arg(label)

  if (!is.numeric(p)) {
    stop("p must be a numeric vector.")
  }
  if (any(p < 0 | p > 1, na.rm = TRUE)) {
    stop("Non-missing p values must be between 0 and 1.")
  }
  if (length(digits) != 1L || is.na(digits) || digits < 0 || digits %% 1 != 0) {
    stop("digits must be one non-negative integer.")
  }
  if (length(scientific) != 1L || is.na(scientific)) {
    stop("scientific must be TRUE or FALSE.")
  }

  label_expr <- if (label == "p") {
    "italic(p)"
  } else {
    paste0("'", label, "'")
  }

  threshold <- 10^(-digits)

  decimal_value <- sprintf(paste0("%.", digits, "f"), p)
  scientific_value <- format(p, scientific = TRUE, digits = digits, trim = TRUE)
  cutoff <- sprintf(paste0("%.", digits, "f"), threshold)

  value <- ifelse(
    p < threshold,
    if (scientific) scientific_value else cutoff,
    decimal_value
  )

  operator <- ifelse(p < threshold & !scientific, "<", "=")

  ifelse(
    is.na(p),
    NA_character_,
    sprintf("%s~'%s'~'%s'", label_expr, operator, value)
  )
}


#### add_strip ####
#' Add alternating backgrounds to discrete ggplot2 panels
#'
#' Insert alternating rectangular background layers behind an existing ggplot.
#' Discrete x axes receive vertical bands and discrete y axes receive horizontal
#' bands. Because the colors are set outside `aes()`, an existing fill scale in
#' the plot is left unchanged.
#'
#' @param p A ggplot object.
#' @param fill Character vector of exactly two colors used alternately.
#' @param alpha Opacity of the background bands, between 0 and 1.
#' @param axis Axis to inspect. `"auto"` and `"both"` add bands to every
#'   discrete axis; `"x"` or `"y"` restrict the operation to one axis.
#'
#' @details Axis breaks are obtained from the first panel of the built plot.
#' Plots with free facet scales may therefore require adding backgrounds to
#' individual panels manually.
#'
#' @return The input ggplot with background rectangle layers prepended.
#' @export
add_strip <- function(
    p, fill = c("white", "grey95"), alpha = 1,
    axis = c("auto", "x", "y", "both")
) {
  axis <- match.arg(axis)

  if (!inherits(p, "ggplot")) {
    stop("p must be a ggplot object.")
  }

  if (!is.character(fill) || length(fill) != 2L || anyNA(fill)) {
    stop("fill must contain exactly two colors.")
  }
  invisible(lapply(fill, grDevices::col2rgb))

  if (length(alpha) != 1L || is.na(alpha) || alpha < 0 || alpha > 1) {
    stop("alpha must be in range of [0, 1].")
  }

  build <- ggplot2::ggplot_build(p)
  panel <- build$layout$panel_params[[1]]

  x_breaks <- panel$x$get_breaks()
  y_breaks <- panel$y$get_breaks()

  x_discrete <- inherits(panel$x$scale, "ScaleDiscrete")
  y_discrete <- inherits(panel$y$scale, "ScaleDiscrete")

  if (axis == "x") {
    y_discrete <- FALSE
  } else if (axis == "y") {
    x_discrete <- FALSE
  }

  strip_layers <- list()

  if (x_discrete && length(x_breaks) > 0) {
    x_df <- data.frame(
      xmin = seq_along(x_breaks) - 0.5,
      xmax = seq_along(x_breaks) + 0.5,
      fill_index = rep(seq_along(fill), length.out = length(x_breaks))
    )

    for (fill_index in seq_along(fill)) {
      strip_layers[[length(strip_layers) + 1L]] <- ggplot2::geom_rect(
        data = x_df[x_df$fill_index == fill_index, , drop = FALSE],
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE, fill = fill[[fill_index]], alpha = alpha,
        color = NA, show.legend = FALSE
      )
    }
  }

  if (y_discrete && length(y_breaks) > 0) {
    y_df <- data.frame(
      ymin = seq_along(y_breaks) - 0.5,
      ymax = seq_along(y_breaks) + 0.5,
      fill_index = rep(seq_along(fill), length.out = length(y_breaks))
    )

    for (fill_index in seq_along(fill)) {
      strip_layers[[length(strip_layers) + 1L]] <- ggplot2::geom_rect(
        data = y_df[y_df$fill_index == fill_index, , drop = FALSE],
        ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE, fill = fill[[fill_index]], alpha = alpha,
        color = NA, show.legend = FALSE
      )
    }
  }

  if (length(strip_layers) > 0) {
    p$layers <- c(strip_layers, p$layers)
  }

  p
}
