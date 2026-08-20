## Replace generic generated parameter text with semantic English descriptions.
## This script prints an apply_patch document and does not edit source files.

args <- commandArgs(trailingOnly = TRUE)
package_dir <- if (length(args)) args[[1]] else '.'
target_file <- if (length(args) >= 2) args[[2]] else NULL
r_dir <- file.path(package_dir, 'R')
source_files <- list.files(r_dir, pattern = '[.]R$', full.names = TRUE)
if (!is.null(target_file)) {
  target_names <- strsplit(target_file, ',', fixed = TRUE)[[1]]
  source_files <- source_files[basename(source_files) %in% target_names]
  if (!length(source_files)) stop('Unknown R source file: ', target_file)
}

describe_param <- function(param) {
  exact <- c(
    profile_x = 'Feature-by-sample profile used for model training or the first data space.',
    profile_y = 'Feature-by-sample profile used for validation or the second data space.',
    group_x = 'Sample metadata associated with `profile_x`.',
    group_y = 'Sample metadata associated with `profile_y`.',
    distance = 'A precomputed distance object; when supplied, it takes precedence over `profile`.',
    dist = 'A distance object used by the analysis.',
    dist_x = 'Distance object for the first data space.',
    dist_y = 'Distance object for the second data space.',
    k = 'Number of folds or clusters, according to the analysis performed.',
    n = 'Requested number of values, features, or results.',
    ntree = 'Number of trees fitted by the random-forest model.',
    repeats = 'Number of repeated cross-validation runs.',
    rep = 'Number of repeated train/test runs.',
    family = 'Model family; supported values are shown in Usage.',
    positive_class = 'Outcome level treated as the positive class for binary metrics.',
    inner_folds = 'Number of internal folds used to tune the regularization parameter.',
    train_prop = 'Proportion of matched samples assigned to the training split.',
    kernel = 'SVM kernel; supported values are shown in Usage.',
    cost_grid = 'Candidate SVM cost values used during tuning.',
    gamma_grid = 'Candidate SVM gamma values used during tuning.',
    tune_boot = 'Number of bootstrap resamples used during SVM tuning.',
    permutations = 'Number of permutations used by the significance test.',
    perplexity = 't-SNE perplexity; it must satisfy the sample-size constraint.',
    trymax = 'Maximum number of random starts attempted by NMDS.',
    max_iter = 'Maximum number of iterations allowed during model fitting.',
    predI = 'Number of predictive components fitted by ropls.',
    orthoI = 'Number of orthogonal components fitted by OPLS-DA.',
    formula = 'Model formula defining the response and grouping variables.',
    comparison = 'Two outcome levels ordered as case and control.',
    measure = 'Effect-size measure passed to `metafor::escalc()`.',
    alternative = 'Alternative hypothesis used to calculate the empirical P value.',
    p_adjust_method = 'Multiple-testing correction method passed to `stats::p.adjust()`.',
    cor_method = 'Correlation method: Pearson, Spearman, or Kendall.',
    adjust = 'Multiple-testing correction method for correlation P values.',
    direction = 'Direction used for prediction, ROC calculation, or network extraction.',
    conf_type = 'Confidence-region geometry: ellipse, encircle, or none.',
    display_type = 'How samples are displayed, such as points alone or points joined to centroids.',
    display = 'Whether composition is summarized by sample or by group.',
    theme = 'Plot theme preset; supported values are shown in Usage.',
    theme_style = 'ROC theme preset; supported values are shown in Usage.',
    aspect_ratio = 'Panel aspect ratio passed to `ggplot2::theme()`.',
    legend_title = 'Legend title; `NULL` uses a context-dependent default.',
    label_pos = 'Length-two numeric vector giving the x and y position of a label.',
    constraint_cols = 'Metadata columns included as constrained variables in dbRDA.',
    covariate_cols = 'Metadata columns whose effects should be retained during batch correction.',
    sub_sample = 'Optional sample identifiers retained before analysis.',
    sub_group = 'Optional group values retained before analysis.',
    sample_metadata = 'Sample metadata required by normalization methods that model library size.',
    tree = 'Phylogenetic tree required for UniFrac distances.',
    weighted = 'Whether to calculate weighted rather than unweighted UniFrac.',
    pseudocount = 'Positive value used to replace or offset zeros before logarithmic operations.',
    use_half_minimum = 'Whether non-positive values are replaced by half the smallest positive value.',
    min_abundance = 'Minimum abundance above which a feature is considered present.',
    min_prevalence = 'Minimum fraction of samples in which a feature must be present.',
    min_n = 'Minimum number of samples in which a feature must be present.',
    n_group = 'Minimum number of groups that must pass the filtering criterion.',
    all_group = 'Whether every group, rather than at least `n_group`, must pass filtering.',
    by_group = 'Whether calculations are performed separately within each group.',
    top_n = 'Number of highest-ranking features or categories retained.',
    top_frac = 'Fraction of highest-ranking features retained.',
    top_list = 'Optional feature or category names retained regardless of abundance rank.',
    sort_method = 'Summary statistic or function used to rank features.',
    sort_value = 'Optional ascending or descending ordering of plotted values.',
    sort_decreasing = 'Whether sorting is performed in decreasing order.',
    na_fill = 'Value used to replace missing observations before analysis.',
    fill_value = 'Value assigned to observations that fail the replacement threshold.',
    fill_missing = 'Distance assigned to feature pairs absent from a pairwise table.',
    unknown_pattern = 'Regular expression identifying unknown or unclassified annotations.',
    other_name = 'Label assigned to features combined into the residual category.',
    library_size = 'External library sizes as a table, file path, or named numeric vector.',
    gene_length = 'Feature-length table used for TPM or FPKM/RPKM normalization.',
    adjacency = 'Adjacency matrix used to create or summarize a network.',
    graph = 'An igraph object to summarize.',
    output_dir = 'Existing directory used for optional file export.',
    output = 'Requested output representation or output path.',
    export = 'Optional export format; `NULL` keeps the operation in memory.',
    file = 'Path to an input file.',
    filename = 'Path of the workbook or output file.',
    path = 'Path to the required input file.',
    path_file = 'Path to the Reactome parent-child relation file.',
    path_info_file = 'Path to the Reactome pathway-information file.',
    xml_file = 'Path to the HMDB XML file.',
    sheets = 'Worksheet names to read; `NULL` reads every worksheet.',
    sheet = 'Worksheet name used for Excel input or output.',
    sep = 'Field separator used when reading or writing a text file.',
    interaction = 'Interaction class or direction retained in network results.',
    type = 'Analysis or value type; supported values are shown in Usage.',
    mode = 'Processing mode; supported values are shown in Usage.',
    style = 'Visual style preset; supported values are shown in Usage.',
    scale_method = 'Normalization method applied before network construction.',
    scale_env = 'Whether environmental variables are standardized before correlation analysis.',
    p_filter = 'Maximum P value retained for plotting or downstream analysis.',
    p_breaks = 'Numeric breakpoints used to categorize P values.',
    p_labels = 'Labels corresponding to intervals defined by `p_breaks`.',
    r_breaks = 'Numeric breakpoints used to categorize correlation or Mantel r.',
    r_labels = 'Labels corresponding to intervals defined by `r_breaks`.',
    r_cut_abs = 'Whether absolute r values are used when assigning line-width categories.',
    curve_bend = 'Signed curvature used for Mantel connection paths.',
    curve_n = 'Number of interpolation points used for each connection path.',
    grid_breaks = 'Reference values used to draw circular or radar grid lines.',
    grid_n = 'Number of grid intervals.',
    grid_max = 'Upper limit of the radar grid; `NULL` derives it from the data.',
    palette = 'Color palette name or vector supplied to the plot.',
    taxonomy = 'Feature taxonomy table used for annotation or aggregation.',
    taxa = 'Taxon names or taxonomy data used by the operation.',
    group_by = 'Metadata column or grouping definition used for aggregation.',
    ids = 'Identifiers to query or retain.',
    set_ID = 'Gene-set identifier selected from a GSEA result.',
    gseaResult = 'A GSEA result object containing ranked genes and enrichment results.',
    cpd_list = 'Metabolite or compound identifiers to annotate.',
    ncm_result = 'Optional result returned by `calcu_NCM()` to avoid refitting the model.',
    adonis_object = 'Object returned by a PERMANOVA/adonis calculation.',
    test = 'Statistical-test result table used for plotting.',
    result = 'Result object or table to summarize.',
    roc = 'A single object returned by `pROC::roc()`.',
    roc_list = 'Named list of objects returned by `pROC::roc()`.',
    x = 'Primary vector or object supplied to the utility.',
    y = 'Secondary vector or object supplied to the utility.',
    abs_value = 'Whether to cluster using absolute pairwise values.',
    adonis2 = 'Whether to run PERMANOVA with `vegan::adonis2()` and attach its result.',
    anosim = 'Whether to run ANOSIM and annotate its statistic and P value.',
    bend = 'Signed curvature of the generated connection path.',
    bias_adjust = 'Whether to apply the small-sample bias correction in `vegan::betadisper()`.',
    by = 'Spacing between successive specificity values used to calculate the ROC confidence band.',
    category_keep = 'Optional environmental categories retained for Mantel tests.',
    ci = 'Whether to calculate confidence intervals for AUC estimates.',
    circular_label = 'Whether category labels should follow the circular plotting direction.',
    cluster_cols = 'Logical value or clustering object controlling column clustering in the heatmap.',
    col_annotation = 'Optional data frame supplying annotations for heatmap columns.',
    collapse_group = 'Whether samples within each group are aggregated before plotting.',
    comment = 'Character marker written before comment rows in the worksheet.',
    coord_flip = 'Whether to exchange the x and y axes with `ggplot2::coord_flip()`.',
    core_enrichment = 'Gene identifiers in the leading-edge or core-enrichment subset.',
    count = 'Whether prevalence is returned as sample counts instead of proportions.',
    cumulative_eig = 'Cumulative explained-variance threshold used to retain ordination axes.',
    dark = 'Text color returned for a sufficiently light background.',
    decreasing = 'Whether categories are sorted in decreasing abundance order.',
    dim = 'Number of ordination dimensions retained in the returned coordinate table.',
    distinct = 'Whether duplicated KEGG mapping records are removed.',
    duplicate_fun = 'Function used to combine duplicated feature-pair values.',
    empty_bar = 'Angular width reserved as an empty separator in the circular bar plot.',
    enriched_by = 'Column used to identify which comparison group is enriched.',
    env_cols = 'Environmental-variable columns selected from `envs`.',
    envs = 'Environmental or host-variable table used in Mantel analyses.',
    exact = 'Whether an exact test is requested when supported by the selected method.',
    expand = 'Numeric expansion applied around groups drawn with `ggforce::geom_mark_ellipse()`.',
    facet = 'Whether enrichment results are separated into facets.',
    fc_method = 'Method used to calculate fold change; supported values are shown in Usage.',
    fc_pseudo = 'Pseudocount added before calculating fold change.',
    filter = 'Optional expression or criterion used to retain converted records.',
    frac = 'Fraction of the highest-ranking features retained.',
    from = 'Source taxonomy rank or identifier type to be converted.',
    hemisphere = 'Portion of the circle used for the pie layout; supported values are shown in Usage.',
    keep_desc = 'Whether KEGG descriptions are retained in the returned table.',
    label = 'Outcome labels corresponding to the profile samples.',
    label_pad = 'Radial padding between circular bars and their labels.',
    label_wrap = 'Maximum label width before circular labels are wrapped.',
    legend.position = 'Legend position passed to `ggplot2::theme()`.',
    light = 'Text color returned for a sufficiently dark background.',
    linetype = 'Line type used for the accumulation or rarefaction curve.',
    linkage_method = 'Hierarchical-clustering linkage method passed to `stats::hclust()`.',
    log_pseudo_factor = 'Multiplier used to derive a log-scale pseudocount from the minimum positive value.',
    logical = 'Whether the adjacency matrix contains logical presence/absence values.',
    mean_only = 'Whether ComBat adjusts batch-specific means without adjusting variances.',
    nrow = 'Number of rows used when arranging multiple plots.',
    obs = 'Observed statistic compared with the randomization distribution.',
    one_minus = 'Whether pairwise values are converted to `1 - value` before clustering.',
    only_signif = 'Whether only statistically significant correlations are retained.',
    organism = 'Reactome organism name used to restrict pathway records.',
    other_last = 'Whether the residual `Other` category is placed last.',
    out_all = 'Whether all hierarchy levels or intermediate results are returned.',
    out_other = 'Whether discarded features are combined into an `Other` category.',
    par_prior = 'Whether parametric empirical-Bayes priors are used by ComBat.',
    paste = 'Whether palette colors are returned as copy-ready quoted text.',
    path_ID = 'Reactome pathway identifier used as the starting node.',
    pca = 'Whether an initial PCA reduction is performed before t-SNE.',
    percent_digits = 'Number of decimal places shown in percentage labels.',
    plab_fmt = 'Format string or function used to display significance labels.',
    prefix = 'Prefix used when naming derived coordinates or labels.',
    prefix_length = 'Number of leading characters used to derive node prefixes.',
    prefix_pattern = 'Regular expression used to extract prefixes from node identifiers.',
    prior_plots = 'Whether ComBat diagnostic prior plots are produced.',
    progress = 'Whether progress information is printed during repeated analyses.',
    progresults = 'Whether XML parsing progress is reported.',
    pval = 'P-value column or matrix paired with the correlation coefficients.',
    pvalue = 'Maximum raw P value retained in the result.',
    qvalue = 'Maximum adjusted P value retained in the result.',
    r = 'Correlation-coefficient column or matrix to be converted to tidy form.',
    random = 'Randomized statistics forming the empirical null distribution.',
    ref_batch = 'Optional reference-batch level passed to ComBat.',
    relation = 'Relationship table or relation type used to connect database identifiers.',
    replace_sheet = 'Whether an existing worksheet with the same name is replaced.',
    rm_suffix = 'Whether suffix text is removed from split taxonomy labels.',
    rownames_fmt = 'Format used to construct row names after KEGG conversion.',
    s_shape = 'Shape or smoothness parameter used for group-enclosure geometry.',
    sort_out = 'Whether the returned set-operation result is sorted.',
    spec_hjust = 'Horizontal justification applied to species-side labels.',
    spec_offset = 'Offset between species labels and the heatmap boundary.',
    spec_range = 'Relative range reserved for the species-label region.',
    split_multi = 'Whether records containing multiple identifiers are expanded into separate rows.',
    star_radius = 'Inner radius used to position the radar-chart star polygons.',
    start = 'Starting angle, in radians, used for the circular coordinate system.',
    step = 'Increment in sequencing depth between rarefaction points.',
    subtitle_keywords = 'Keywords whose occurrences are highlighted in the plot subtitle.',
    suffix = 'Suffix removed from or appended to derived identifiers.',
    test_trans = 'Transformation applied to feature values before statistical testing.',
    theta = 't-SNE speed-accuracy trade-off; zero requests exact t-SNE.',
    to = 'Target taxonomy rank or identifier type produced by the conversion.',
    trans = 'Transformation applied to abundance values before plotting.',
    trans_ra = 'Whether abundances are converted to relative abundance before analysis.',
    unique_out = 'Whether duplicate values are removed from the set-operation result.',
    unknown = 'Label assigned to unknown or unclassified taxonomy entries.',
    var_equal = 'Whether two-sample t tests assume equal group variances.',
    x_limit = 'Optional numeric limits for the x axis.',
    x_text_angle = 'Rotation angle, in degrees, for x-axis text.',
    xend = 'Ending x coordinate of the generated connection path.',
    yend = 'Ending y coordinate of the generated connection path.'
  )
  if (param %in% names(exact)) return(unname(exact[[param]]))

  label <- gsub('_', ' ', param, fixed = TRUE)
  if (grepl('^(add|show|remove|drop|write|rotate|flip|normalize|symmetric|parallel)', param)) {
    return(paste0('Whether to enable the ', label, ' behavior.'))
  }
  if (grepl('(_file|file)$', param)) return(paste0('Path used for ', label, '.'))
  if (grepl('(_mat|matrix)$', param)) return(paste0('Matrix supplying ', label, '.'))
  if (grepl('(_breaks)$', param)) return(paste0('Numeric breakpoints for ', label, '.'))
  if (grepl('(_labels)$', param)) return(paste0('Labels corresponding to ', label, '.'))
  if (grepl('(_name|name)$', param)) return(paste0('Display or identifier name for ', label, '.'))
  if (grepl('^(min|max|top|n_)', param)) return(paste0('Numeric limit controlling ', label, '.'))
  if (grepl('(color|fill)', param)) return(paste0('Color specification for ', label, '.'))
  if (grepl('(title|xlab|ylab)$', param)) return(paste0('Optional label used for ', label, '.'))
  paste0('Value used to configure ', label, ' in this analysis.')
}

patch_lines <- '*** Begin Patch'
for (source_file in source_files) {
  source_lines <- readLines(source_file, warn = FALSE, encoding = 'UTF-8')
  matches <- grep(
    paste0(
      "^#' @param [^ ]+ (Input controlling `[^`]+`|",
      "Value used to configure .+ in this analysis)[.]$"
    ),
    source_lines
  )
  if (!length(matches)) next

  file_hunks <- character()
  for (line_no in matches) {
    param <- sub("^#' @param ([^ ]+) .*$", '\\1', source_lines[line_no])
    replacement <- paste0("#' @param ", param, ' ', describe_param(param))
    file_hunks <- c(
      file_hunks, '@@', paste0('-', source_lines[line_no]),
      paste0('+', replacement)
    )
  }
  patch_lines <- c(
    patch_lines,
    paste0('*** Update File: ', normalizePath(source_file, winslash = '\\')),
    file_hunks
  )
}
writeLines(c(patch_lines, '*** End Patch'), useBytes = TRUE)
