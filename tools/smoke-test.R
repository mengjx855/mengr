## mengR source-level smoke tests. This script does not build or install the package.

package_dir <- normalizePath(file.path(getwd(), 'mengR'), winslash = '/')
source_files <- list.files(
  file.path(package_dir, 'R'), pattern = '[.]R$', full.names = TRUE
)
test_env <- new.env(parent = globalenv())
for (source_file in source_files) sys.source(source_file, envir = test_env)

with(test_env, {
  profile <- data.frame(
    S1 = c(10, 20, 30), S2 = c(20, 40, 40),
    row.names = c('g1', 'g2', 'g3'), check.names = FALSE
  )

  ## 1. read-count conversion and external library-size alignment
  rpm_default <- rc2rpm(profile)
  stopifnot(all.equal(colSums(rpm_default), c(S1 = 1e6, S2 = 1e6)))
  library_size_df <- data.frame(
    sample = c('S2', 'S1'), library_size = c(200, 100)
  )
  rpm_table <- rc2rpm(profile, library_size_df)
  stopifnot(all.equal(unname(rpm_table[1, ]), c(1e5, 1e5)))

  library_size_file <- tempfile(fileext = '.tsv')
  utils::write.table(
    library_size_df, library_size_file, sep = '\t', row.names = FALSE,
    quote = FALSE
  )
  rpm_file <- rc2rpm(profile, library_size_file)
  stopifnot(all.equal(rpm_file, rpm_table))

  ## 2. core transformations and chemistry helpers
  stopifnot(all.equal(unname(calcu_DBE(c('C6H6', 'C2H5Cl'))), c(4, 0)))
  stopifnot(all.equal(
    unname(colSums(profile_trans_ra(profile, base = 1))), c(1, 1)
  ))
  stopifnot(nrow(profile_top_n(profile, n = 2)) == 2)

  ## 2b. difference-analysis non-positive replacement and main workflow
  raw_mat <- matrix(
    c(0, 2, 4, 8, 1, 3, 5, 7), nrow = 2,
    dimnames = list(c('g1', 'g2'), paste0('S', 1:4))
  )
  log2_mat <- log2(.replace_nonpositive(raw_mat, pseudo_factor = 0.5))
  stopifnot(all(is.finite(log2_mat)), identical(dim(log2_mat), dim(raw_mat)))

  difference_df <- suppressMessages(difference_analysis(
    raw_mat,
    data.frame(
      sample = paste0('S', 1:4), group = c('case', 'case', 'control', 'control')
    ),
    comparison = c('case', 'control'), method = 't', test_trans = 'none',
    progress = FALSE
  ))
  stopifnot(
    nrow(difference_df) == 2,
    all(c('name', 'FC', 'log2FC', 'pval', 'padj') %in% colnames(difference_df))
  )

  ## 3. ecological distance and enterotype helpers
  distance <- calcu_distance(profile, dist_method = 'bray')
  stopifnot(inherits(distance, 'dist'))
  jsd <- calcu_jsd_dist(profile)
  stopifnot(inherits(jsd, 'dist'), all(as.matrix(jsd) >= 0))

  ## 4. rarefaction and categorical tests
  spec_df <- calcu_specaccum(profile, method = 'collector')
  stopifnot(all(c('sample_n', 'richness', 'sd') %in% colnames(spec_df)))
  fisher_df <- calcu_fisher(
    data.frame(
      feature = c('a', 'b'), x_pos = c(3, 1), y_pos = c(2, 4),
      x_neg = c(7, 9), y_neg = c(8, 6)
    ),
    feature_col = 'feature'
  )
  stopifnot(nrow(fisher_df) == 2, all(is.finite(fisher_df$pval)))

  ## 5. annotation order should be aligned by name
  row_annotation <- data.frame(
    name = c('g3', 'g1', 'g2'), class = c('B', 'A', 'A')
  )
  col_annotation <- data.frame(
    name = c('S2', 'S1'), group = c('case', 'control')
  )
  heatmap <- plot_heatmap(
    profile, scale = 'none', row_annotation = row_annotation,
    col_annotation = col_annotation, cellwidth = NA, cellheight = NA,
    run_draw = FALSE
  )
  stopifnot(inherits(heatmap, 'Heatmap'))

  ## 6. CellMarker database build and match
  marker_dir <- tempfile('cellmarker-')
  marker_df <- data.frame(
    species = c('Mouse', 'Mouse', 'Human'),
    tissue_class = c('Intestine', 'Intestine', 'Intestine'),
    cell_name = c('Stem cell', 'Stem cell', 'Stem cell'),
    marker = c('Lgr5', 'Smoc2', 'LGR5'),
    Symbol = c('Lgr5', 'Smoc2', 'LGR5'),
    GeneID = c('123', '456', '789')
  )
  marker_file <- scRNA_marker_db_build(marker_df, database = marker_dir)
  stopifnot(file.exists(marker_file))
  marker_match <- scRNA_marker_match(
    c('Lgr5', 'Mki67'), database = marker_dir
  )
  stopifnot(marker_match$matched_n[1] == 1)

  ## 7. Seurat expression summary/dotplot and CellChat radar extraction
  count_mat <- Matrix::Matrix(
    matrix(
      c(3, 0, 1, 0, 2, 0, 1, 1), nrow = 2,
      dimnames = list(c('G1', 'G2'), paste0('C', 1:4))
    ),
    sparse = TRUE
  )
  seurat <- Seurat::CreateSeuratObject(counts = count_mat)
  seurat$celltype <- c('A', 'A', 'B', 'B')
  seurat <- Seurat::NormalizeData(seurat, verbose = FALSE)
  expression_df <- scRNA_summarise_genes(
    seurat, group_cols = 'celltype', genes = c('G1', 'G2')
  )
  stopifnot(nrow(expression_df) == 4)
  expression_plot <- scRNA_expression_dotplot(
    expression_df, x_col = 'celltype', y_col = 'gene'
  )
  stopifnot(inherits(expression_plot, 'ggplot'))

  if (!methods::isClass('mengR_mock_cellchat')) {
    methods::setClass('mengR_mock_cellchat', slots = c(net = 'list'))
  }
  network_mat <- matrix(
    c(1, 2, 3, 4), nrow = 2,
    dimnames = list(c('A', 'B'), c('A', 'B'))
  )
  mock_cellchat <- methods::new(
    'mengR_mock_cellchat', net = list(count = network_mat, weight = network_mat)
  )
  radar_plot <- cellchat_radar(mock_cellchat, cell = 'A')
  stopifnot(inherits(radar_plot, 'ggplot'))
})

message('All source-level smoke tests passed.')
