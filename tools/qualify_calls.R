args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1]] else "."
r_dir <- file.path(root, "R")

package_functions <- list(
  dplyr = c(
    "across", "all_of", "any_of", "arrange", "case_when", "desc",
    "distinct", "everything", "filter", "first", "group_by", "group_keys",
    "group_map", "group_modify", "if_all", "if_else", "inner_join",
    "left_join", "mutate", "n", "pull", "relocate", "rename", "rename_all",
    "rename_with", "right_join", "row_number", "rowwise", "select",
    "slice_head", "slice_max", "slice_min", "summarise", "summarise_all",
    "transmute", "ungroup", "where"
  ),
  tibble = c(
    "add_column", "add_row", "as_tibble", "column_to_rownames",
    "rownames_to_column", "tibble"
  ),
  tidyr = c(
    "complete", "drop_na", "fill", "gather", "nest", "pivot_longer",
    "pivot_wider", "replace_na", "separate", "spread", "unite", "unnest"
  ),
  purrr = c(
    "compact", "imap_dfr", "list_flatten", "map", "map_chr", "map_dbl",
    "map_df", "map_dfc", "map_dfr", "map_int", "map_vec", "map2",
    "map2_df", "map2_dfc", "map2_dfr", "map2_vec", "pmap", "reduce"
  ),
  rlang = c("is_null", "set_names"),
  stringr = c(
    "str_extract_all", "str_flatten", "str_replace_all", "str_split",
    "str_split_1", "str_split_i", "str_sub", "str_to_lower",
    "str_to_sentence", "str_to_title", "str_wrap"
  ),
  ggplot2 = c(
    "aes", "annotate", "coord_cartesian", "coord_fixed", "coord_flip",
    "coord_polar", "element_blank", "element_line", "element_rect",
    "element_text", "expansion", "facet_grid", "facet_wrap", "geom_abline",
    "geom_bar", "geom_boxplot", "geom_col", "geom_errorbar", "geom_hline",
    "geom_jitter", "geom_label", "geom_line", "geom_path", "geom_point",
    "geom_polygon", "geom_rect", "geom_ribbon", "geom_segment", "geom_text",
    "geom_tile", "geom_violin", "geom_vline", "ggplot", "ggplot_build", "ggsave",
    "guide_colorbar", "guide_legend", "guides", "labs", "margin",
    "position_dodge", "position_identity", "position_stack",
    "scale_color_manual", "scale_colour_manual", "scale_fill_gradient",
    "scale_fill_gradientn", "scale_fill_manual", "scale_linetype_manual",
    "scale_shape_manual", "scale_size_continuous", "scale_x_continuous",
    "scale_x_discrete", "scale_y_continuous", "scale_y_discrete",
    "stat_ellipse", "theme", "theme_bw", "theme_classic", "theme_minimal",
    "theme_void", "vars"
  ),
  ggpubr = c(
    "ggarrange", "ggbarplot", "ggboxplot", "ggtext", "ggtexttable",
    "stat_compare_means", "tab_add_hline", "theme_pubr", "ttheme"
  ),
  vegan = c(
    "adonis2", "anosim", "betadisper", "capscale", "decostand", "diversity",
    "estimateR", "metaMDS", "permutest", "procrustes", "protest", "rrarefy",
    "RsquareAdj", "specaccum", "specnumber", "vegdist"
  ),
  igraph = c(
    "all_simple_paths", "as_edgelist", "degree", "edge_density", "E",
    "graph_from_adjacency_matrix", "graph_from_data_frame",
    "graph_from_edgelist", "V", "write.graph"
  ),
  ggClusterNet = c(
    "cor_Big_micro2", "edgeBuild", "model1", "model2", "model_igraph2",
    "net_properties.2",
    "nodeadd", "nodeEdge", "vegan_otu", "vegan_tax"
  ),
  phyloseq = c("otu_table", "phy_tree", "phyloseq", "tax_table", "UniFrac"),
  Seurat = "FeaturePlot",
  cluster = "pam",
  clusterProfiler = "enricher",
  Hmisc = "binconf",
  pROC = "roc",
  minpack.lm = c("nlsLM", "nls.lm.control"),
  xcms = c("intensity", "rtime"),
  randomForest = c("importance", "randomForest"),
  e1071 = c("svm", "tune.control", "tune.svm"),
  glmnet = c("cv.glmnet", "glmnet"),
  ComplexHeatmap = "pheatmap",
  grid = c("arrow", "gpar", "unit"),
  graphics = c("legend", "par", "plot", "text"),
  grDevices = c("colorRampPalette", "dev.off", "pdf"),
  stats = c(
    "aggregate", "anova", "as.dist", "as.formula", "cmdscale", "coef",
    "confint", "cor", "cutree", "density", "dist", "fisher.test",
    "formula", "hclust",
    "median", "model.matrix", "na.omit", "oneway.test", "p.adjust",
    "pbeta", "prcomp", "predict", "quantile", "sd", "setNames", "t.test",
    "terms", "TukeyHSD", "var", "wilcox.test"
  ),
  utils = c(
    "combn", "head", "read.delim", "setTxtProgressBar", "tail",
    "txtProgressBar", "write.table"
  )
)

mapping <- unlist(lapply(names(package_functions), function(pkg) {
  stats::setNames(rep(pkg, length(package_functions[[pkg]])), package_functions[[pkg]])
}))

qualify_file <- function(file) {
  parsed <- parse(file, keep.source = TRUE)
  tokens <- getParseData(parsed)
  tokens <- tokens[
    tokens$terminal & tokens$token == "SYMBOL_FUNCTION_CALL" &
      tokens$text %in% names(mapping),
    , drop = FALSE
  ]

  if (!nrow(tokens)) return(invisible(FALSE))

  lines <- readLines(file, warn = FALSE, encoding = "UTF-8")
  tokens <- tokens[order(tokens$line1, tokens$col1, decreasing = TRUE), ]

  for (i in seq_len(nrow(tokens))) {
    token <- tokens[i, ]
    line <- lines[[token$line1]]
    before <- if (token$col1 > 1) substr(line, 1, token$col1 - 1) else ""

    # Calls already using pkg::fun() have an NS_GET token immediately before
    # the function symbol. Checking the source prefix preserves their spelling.
    if (grepl("(::|:::)[[:space:]]*$", before)) next

    replacement <- paste0(mapping[[token$text]], "::", token$text)
    after <- if (token$col2 < nchar(line)) {
      substr(line, token$col2 + 1, nchar(line))
    } else {
      ""
    }
    lines[[token$line1]] <- paste0(before, replacement, after)
  }

  writeLines(lines, file, useBytes = TRUE)
  invisible(TRUE)
}

files <- list.files(r_dir, pattern = "[.]R$", full.names = TRUE)
invisible(lapply(files, qualify_file))
