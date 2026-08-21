#### Jinxin Meng, 20241029, 20250903 v0.2 ####
# 20241029: create functions.
# 20250903: fix bug.


#### plot_msa ####
# input a msa with fa-formated
# muscle/mafft and trimal
#' Plot Msa utility
#'
#' `plot_msa()` provides a reusable mengR workflow with input validation and standardized
#'   output.
#'
#' Chinese summary: 读取 multiple sequence alignment 并按氨基酸类别分块绘图。
#'
#' @param path Path to the required input file.
#' @param width Numeric setting for `width`.
#' @param font_size Numeric setting for `font_size`.
#' @return A plot object; analysis data or models may also be stored as attributes.
#' @export
plot_msa <- function(path, width = 80, font_size = 2) {
  ## 1. 读取 alignment，并把每条序列拆成单字符列
  sequence_list <- Biostrings::readAAMultipleAlignment(path) |>
    as.character() |>
    as.list()
  msa <- purrr::map2_dfc(
    sequence_list, names(sequence_list),
    function(sequence, sequence_name) {
      data.frame(aa = unlist(stringr::str_split(sequence, ""))) |>
        dplyr::rename(!!sequence_name := aa)
    }
  ) |>
    tibble::rownames_to_column("site")

  seq_lens <- nrow(msa)
  num_seqs <- ceiling(seq_lens / width)
  total_lens <- num_seqs * width
  num_ids <- ncol(msa) - 1

  data <- purrr::map2(
    seq(0, width * (num_seqs - 1), width),
    seq(width, width * num_seqs, width), \(x, y)
    dplyr::filter(msa, site %in% seq(x + 1, y, 1)) |>
      tidyr::gather(key = "name", value = "aa", -site) |>
      dplyr::mutate(site = factor(site, seq(x + 1, y, 1)))
  ) |>
    rlang::set_names(paste0("seq_", 1:num_seqs))

  if (seq_lens < total_lens) {
    data[[num_seqs]] <- rbind(
      data[[num_seqs]],
      data.frame(
        site = rep(seq(seq_lens + 1, total_lens, 1), each = num_ids),
        name = rep(unique(data[[num_seqs]]$name), times = (total_lens - seq_lens))
      ) |>
        tibble::add_column(aa = " ")
    ) |>
      dplyr::mutate(site = factor(site, seq(((num_seqs - 1) * width + 1), total_lens, 1)))
  }

  aa_col <- c(
    "E" = "#ff6d6d", "D" = "#ff6d6d", "P" = "#f2be3c", "A" = "#f2be3c", "V" = "#f2be3c",
    "M" = "#f2be3c", "L" = "#f2be3c", "I" = "#f2be3c", "G" = "#f2be3c", "K" = "#769dcc",
    "R" = "#769dcc", "H" = "#769dcc", "N" = "#74ce98", "T" = "#74ce98", "C" = "#74ce98",
    "Q" = "#74ce98", "S" = "#74ce98", "F" = "#ffff66", "Y" = "#ffff66", "W" = "#ffff66",
    "-" = "#ffffff", " " = "#ffffff"
  )

  plots <- purrr::map(data, \(x)
  ggplot2::ggplot(x, ggplot2::aes(site, name, fill = aa)) +
    ggplot2::geom_tile(color = "white", linewidth = .4, show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = aa), size = font_size) +
    ggplot2::labs(x = "", y = "") +
    ggplot2::scale_x_discrete(
      breaks = c(
        as.character(min(as.numeric(as.character(x$site)))),
        as.character(max(as.numeric(as.character(x$site))))
      ),
      labels = c(
        as.character(min(as.numeric(as.character(x$site)))),
        as.character(max(as.numeric(as.character(x$site))))
      )
    ) +
    ggplot2::scale_fill_manual(values = aa_col) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 8),
      axis.text.y = ggplot2::element_text(size = 8, face = "italic", hjust = 1),
      axis.ticks = ggplot2::element_blank()
    ))
  cowplot::plot_grid(plotlist = plots, ncol = 1)
}
