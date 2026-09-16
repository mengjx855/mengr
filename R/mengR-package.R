#### Jin-Xin Meng, jinxmeng@zju.edu.cn, 20260820, 20260916 ####

# 20260916: standardize package documentation and declare data-masking variables used across functions.

#' mengR: Personal bioinformatics analysis utilities
#'
#' Reusable functions for omics data processing, statistical analysis,
#' machine learning, and publication-oriented visualization.
#'
#' @author Jin-Xin Meng \email{jinxmeng@zju.edu.cn}
#' @importFrom rlang .data
#' @keywords internal
"_PACKAGE"

# Names evaluated non-standardly by dplyr, data.table, ggplot2, and modelling
# helpers. Declaring them keeps R CMD check focused on genuine undefined objects.
utils::globalVariables(c(
  ".", ".abs_coef", ".coef", ".direction", ".empty", ".feature", ".I",
  ".label", ".n", ".N", ".order", ".p", ".rank", ".rank_letter",
  ".sig_group", ".sort", ".taxon", ".term_id", ".x_id", ".x_jit", ":=",
  "aa", "abs_estimate", "abs_r", "abundance", "actual", "anchor_x",
  "anchor_y", "angle", "angle_position", "background_fraction",
  "background_n", "category", "cell_id", "celltype", "clade_name", "coef",
  "color", "color1", "comparison", "Comparison", "conf_high", "conf_low",
  "cor", "db_id", "depth", "Description", "distance", "dx", "dy",
  "edge_class", "edge_id", "effect", "eligible", "end", "enriched", "env",
  "error", "estimate", "expected_target", "FDR", "feature", "feature1",
  "feature2", "Freq", "frequency", "gap_end", "gap_start", "group", "Group",
  "group1", "group2", "hjust", "i.cell_id", "id", "ID", "idx", "index",
  "ix", "ix_latest", "iy", "iy_latest", "label", "label_angle", "LDA",
  "log_abundance", "log2_enrichment", "logabun", "lower", "lvAdes", "lvD",
  "Mean", "mean_abundance", "MeanDecreaseAccuracy", "median", "Method",
  "model", "model_feature", "model_label", "model1", "model2", "N", "n_axis",
  "n_cells", "n_target", "name", "NES", "NMDS", "node1", "node2",
  "null_hypothesis", "OTU", "p", "P.adj", "p.adjust", "P.unadj", "p_value",
  "padj", "pcut", "percent", "plab", "plot_pvalue", "poly_xmax", "poly_xmin",
  "poly_ymax", "poly_ymin", "predicted", "prev1", "prev2", "proj", "pval",
  "qval", "r", "r_size", "r2", "RANK", "rcut", "richness", "roi_id",
  "rowname", "sample_n", "score", "sd", "Sd", "sig", "signal",
  "Significance", "site", "spec", "spec_x", "spec_y", "start",
  "target_fraction", "target_pass", "tax", "Taxa", "taxa_level", "taxon",
  "term", "TERM", "to", "type", "upper", "value", "Var1", "Var2",
  "variable", "vertex_x", "vertex_y", "vi", "x", "X", "x_Mean", "x_N",
  "x_Sd", "X1", "X1_rotated", "X1_target", "X1mean", "X2", "X2_rotated",
  "X2_target", "X2mean", "xmax", "xmin", "y", "Y", "y_Mean", "y_N",
  "y_position", "y_Sd", "Y1", "Y2", "yi", "ymax", "ymin"
))
