#### Jin-Xin Meng, jinxmeng@zju.edu.cn, 20251229, 20260923 ####

# time-stamp:
# 20251229: create script, add function 'plot_gsea_barcode()'.
# 20260527: add functions `plot_go_bar()` and `plot_go_circular_bar()`.
# 20260822: add function 'run_limma_diff()', 'run_voom_diff()', 'run_deseq2_diff()'
#           for RNA-seq differential analysis in different scenarios.
# 20260916: rename GO plotting functions to lowercase, standardize documentation, and adopt `*_df` data-frame names.
# 20260923: clarify metadata argument names and remove Chinese text from Roxygen documentation.


#### run_limma_diff ####
# 使用 limma 对连续型 profile 数据进行组间差异分析
# 使用场景：
#   适用于已经标准化并转换到适合线性模型分析尺度的连续型 profile，例如：
#   - log2(TPM + 1)
#   - log2(FPKM + 1)
#   - log2(CPM + 1)
#   - 标准化后的 microarray expression matrix
#   - 其他经过适当转换的连续型 abundance/profile 数据
#   如果有 RNA-seq raw count，优先考虑：
#   - run_voom_diff()
#   - run_deseq2_diff()
# 输入：
#   profile:
#     数值型 matrix 或 data.frame。
#     行为 feature，列为 sample；
#     rownames 为 feature 名称，colnames 为样本名称。
#   metadata:
#     样本信息表，至少包含 sample_col 和 group_col 指定的两列。
#   comparisons:
#     pairwise comparison list。
#     例如：
#       list(
#         c("IBD", "HC"),
#         c("CD", "HC")
#       )
#     c("IBD", "HC") 表示 IBD vs HC：
#       log2FC > 0 -> IBD 更高
#       log2FC < 0 -> HC 更高
#   sample_col:
#     metadata 中样本名称所在列。
#   group_col:
#     metadata 中分组信息所在列。
#   fc_cutoff:
#     判断差异 feature 的 |log2FC| 阈值，默认 1。
#   padj_cutoff:
#     判断差异 feature 的 BH-adjusted P value 阈值，默认 0.05。
# 返回：
#   list：
#   - difference:
#       每个 comparison 的完整 limma 差异分析结果。
#   - summary:
#       每个 comparison 中两个差异方向及非差异 feature 的数量。
#
# 注意事项：
#   1. 不要将 RNA-seq raw count 直接输入本函数。
#   2. TPM/FPKM 等通常应先进行 log2 转换，例如：
#        profile <- log2(profile + 1)
#   3. log2FC 实际表示输入 log2 scale 上的组间均值差。
#      对 log2(TPM + 1) 等数据，不严格等同于两个组原始 TPM 均值之比的 log2。
#   4. comparisons 中第一个组为 numerator，第二个组为 denominator。
#   5. limma::makeContrasts() 要求 group 名称适合作为 R model coefficient，
#      建议使用 WT、KO、Control、Treat 等简单名称，避免空格和 "-"。

run_limma_diff <- function(
    profile, sample_meta, comparisons, sample_col = "sample", group_col = "group",
    fc_cutoff = 1, padj_cutoff = .05
) {
  metadata <- sample_meta
  
  # 匹配 metadata 与 profile 中的样本，并统一样本顺序
  metadata <- metadata |>
    dplyr::filter(.data[[sample_col]] %in% colnames(profile)) |>
    dplyr::arrange(match(.data[[sample_col]], colnames(profile)))
  
  profile <- profile[, metadata[[sample_col]], drop = FALSE]
  
  # 构建分组因子，保持 metadata 中首次出现的分组顺序
  group <- factor(metadata[[group_col]], unique(metadata[[group_col]]))
  
  # 无截距 design，每个 group 对应一个 coefficient
  design <- stats::model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  
  # 根据 comparisons 构建 contrasts
  # 例如 c("IBD", "HC") -> IBD - HC
  contrasts <- limma::makeContrasts(
    contrasts = purrr::map_vec(
      comparisons, ~ paste(.x, collapse = "-")
    ),
    levels = design
  )
  
  # 拟合线性模型并进行 empirical Bayes moderation
  fit <- limma::lmFit(profile, design) |>
    limma::contrasts.fit(contrasts) |>
    limma::eBayes()
  
  # 提取每个 comparison 的差异分析结果
  difference <- purrr::map(
    comparisons, ~
      limma::topTable(
        fit, coef = paste(.x, collapse = "-"),
        adjust.method = "BH", number = Inf
      ) |>
      tibble::rownames_to_column("name") |>
      dplyr::rename(
        log2FC  = logFC,
        pvalue  = P.Value,
        padjust = adj.P.Val
      ) |>
      dplyr::mutate(
        enriched = dplyr::case_when(
          .data[["log2FC"]] >  fc_cutoff & .data[["padjust"]] < padj_cutoff ~ .x[1],
          .data[["log2FC"]] < -fc_cutoff & .data[["padjust"]] < padj_cutoff ~ .x[2],
          TRUE ~ "none"
        )
      )
  ) |>
    purrr::set_names(
      purrr::map_vec(
        comparisons, ~ paste(.x, collapse = "_vs_")
      )
    )
  
  # 汇总各 comparison 中的差异 feature 数量
  summary <- purrr::map2_dfr(
    names(difference), comparisons, \(x, y)
    data.frame(
      comparison = x,
      up   = sum(difference[[x]]$enriched == y[1]),
      down = sum(difference[[x]]$enriched == y[2]),
      none = sum(difference[[x]]$enriched == "none")
    )
  )
  
  list(
    difference = difference,
    summary = summary
  )
}

#### run_voom_diff ####
# 使用 voom/voomLmFit + limma 对 raw count 数据进行差异分析
# 使用场景：
#   适用于 RNA-seq 等 sequencing count 数据，例如：
#   - featureCounts gene count
#   - HTSeq count
#   - 其他 gene/feature-level sequencing count matrix
#   与其他两个函数的区别：
#   run_limma_diff()
#     输入已经转换后的连续型 profile。
#   run_voom_diff()
#     输入 count matrix，
#     使用 voom 建模 mean-variance relationship 后进入 limma。
#   run_deseq2_diff()
#     输入 raw integer count，
#     使用 negative binomial GLM。
# 输入：
#   rc:
#     非负 count matrix 或 data.frame。
#     行为 feature，列为 sample；
#     rownames 为 feature 名称，colnames 为样本名称。
#   metadata:
#     样本信息表，至少包含 sample_col 和 group_col 指定的两列。
#   comparisons:
#     pairwise comparison list。
#     例如：
#       list(
#         c("KO", "WT"),
#         c("Treat", "Control")
#       )
#     c("KO", "WT") 表示 KO vs WT：
#       log2FC > 0 -> KO 更高
#       log2FC < 0 -> WT 更高
#   sample_col:
#     metadata 中样本名称所在列。
#   group_col:
#     metadata 中分组信息所在列。
#   method:
#     voom 分析方法：
#     "voomLmFit":
#       使用 edgeR::voomLmFit()。
#       默认方法，对含较多 exact zero count 的数据处理更完善。
#     "voom":
#       使用经典 limma::voom() + limma::lmFit()。
#   fc_cutoff:
#     判断差异 feature 的 |log2FC| 阈值，默认 1。
#   padj_cutoff:
#     判断差异 feature 的 BH-adjusted P value 阈值，默认 0.05。
#   plot:
#     是否绘制 voom mean-variance trend，默认 FALSE。
# 返回：
#   list：
#   - difference:
#       每个 comparison 的 limma 差异分析结果。
#   - logCPM:
#       经过低表达过滤、TMM normalization 和 voom 转换后的
#       normalized log2-CPM matrix。
#       可用于 PCA、heatmap、聚类和样本相关性分析。
#   - summary:
#       每个 comparison 的差异 feature 数量。
# 注意事项：
#   1. 输入应为 sequencing count，不要输入 TPM、FPKM 或 log2 expression。
#   2. 函数内部使用 edgeR::filterByExpr() 去除低表达 feature。
#   3. 使用 edgeR::calcNormFactors() 计算 TMM normalization factor。
#      calcNormFactors() 不会直接修改原始 count，而是修正 effective library size。
#   4. voom 根据 count 的 mean-variance relationship 计算 precision weights。
#   5. logCPM 只包含通过 filterByExpr() 的 feature。
#   6. comparisons 中第一个组为 numerator，第二个组为 denominator。
#   7. group 名称建议使用简单且合法的 R coefficient 名称。

run_voom_diff <- function(
    rc, sample_meta, comparisons, sample_col = "sample", group_col = "group",
    method = c("voomLmFit", "voom"), fc_cutoff = 1, padj_cutoff = .05,
    plot = FALSE
) {
  metadata <- sample_meta
  
  method <- match.arg(method)
  
  # 匹配 metadata 与 count matrix，并统一样本顺序
  metadata <- metadata |>
    dplyr::filter(.data[[sample_col]] %in% colnames(rc)) |>
    dplyr::arrange(match(.data[[sample_col]], colnames(rc)))
  
  rc <- rc[, metadata[[sample_col]], drop = FALSE] |>
    as.matrix()
  
  # 构建分组因子
  group <- factor(metadata[[group_col]], unique(metadata[[group_col]]))
  
  # 无截距 design，每个 group 对应一个 coefficient
  design <- stats::model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  
  # 构建 edgeR count object
  dge <- edgeR::DGEList(counts = rc, group = group)
  
  # 根据表达水平、library size 和组内样本量过滤低表达 feature
  keep <- edgeR::filterByExpr(dge, group = group)
  
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  # TMM normalization
  # normalization factor 会参与 effective library size 的计算
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  
  if (method == "voomLmFit") {
    # 同时进行 voom mean-variance 建模和 linear model fitting
    fit <- edgeR::voomLmFit(dge, design = design, plot = plot, keep.EList = TRUE)
    
    # 提取 normalized log2-CPM
    logCPM <- fit$EList$E
    
  } else {
    
    # 经典 voom：生成 log2-CPM 和 observation-level weights
    voom <- limma::voom(dge, design = design, plot = plot)
    
    # 使用 voom weights 拟合 linear model
    fit <- limma::lmFit(voom, design)
    
    logCPM <- voom$E
  }
  
  # 根据 comparisons 构建 contrast matrix
  contrasts <- limma::makeContrasts(
    contrasts = purrr::map_vec(
      comparisons, ~ paste(.x, collapse = "-")
    ),
    levels = design
  )
  
  # 应用 contrasts 并进行 empirical Bayes moderation
  fit <- fit |>
    limma::contrasts.fit(contrasts) |>
    limma::eBayes()
  
  # 提取每个 comparison 的差异分析结果
  difference <- purrr::map(
    comparisons, ~
      limma::topTable(
        fit, coef = paste(.x, collapse = "-"), 
        adjust.method = "BH", number = Inf
      ) |>
      tibble::rownames_to_column("name") |>
      dplyr::rename(
        log2FC  = logFC,
        pvalue  = P.Value,
        padjust = adj.P.Val
      ) |>
      dplyr::mutate(
        enriched = dplyr::case_when(
          .data[["log2FC"]] >  fc_cutoff & .data[["padjust"]] < padj_cutoff ~ .x[1],
          .data[["log2FC"]] < -fc_cutoff & .data[["padjust"]] < padj_cutoff ~ .x[2],
          TRUE ~ "none"
        )
      )
  ) |>
    purrr::set_names(
      purrr::map_vec(
        comparisons, ~ paste(.x, collapse = "_vs_")
      )
    )
  
  # 汇总差异 feature 数量
  summary <- purrr::map2_dfr(
    names(difference), comparisons, \(x, y)
    data.frame(
      comparison = x,
      up   = sum(difference[[x]]$enriched == y[1]),
      down = sum(difference[[x]]$enriched == y[2]),
      none = sum(difference[[x]]$enriched == "none")
    )
  )
  
  list(
    difference = difference,
    logCPM = logCPM,
    summary = summary
  )
}

#### run_deseq2_diff ####
# 使用 DESeq2 对 raw count 数据进行组间差异分析
# 使用场景：
#   适用于 RNA-seq 等基于 sequencing count 的差异分析，例如：
#   - featureCounts gene count
#   - HTSeq count
#   - 其他非负整数型 gene/feature count matrix
#   输入必须是 raw integer count，不能是：
#   - TPM
#   - FPKM
#   - CPM
#   - log2 transformed expression
#   - VST/rlog transformed expression
# 输入：
#   rc:
#     非负整数型 raw count matrix 或 data.frame。
#     行为 feature，列为 sample；
#     rownames 为 feature 名称，colnames 为样本名称。
#   metadata:
#     样本信息表，至少包含 sample_col 和 group_col 指定的两列。
#   comparisons:
#     pairwise comparison list。
#     例如：
#       list(
#         c("KO", "WT"),
#         c("Treat", "Control")
#       )
#     c("KO", "WT") 表示 KO vs WT：
#       log2FC > 0 -> KO 更高
#       log2FC < 0 -> WT 更高
#   sample_col:
#     metadata 中样本名称所在列。
#   group_col:
#     metadata 中分组信息所在列。
#   pairwise:
#     FALSE：
#       所有组共同建立一个 DESeq2 模型，
#       然后从同一模型中提取不同 contrasts。
#       默认并推荐这种方式。
#     TRUE：
#       每个 comparison 单独提取两个组，
#       并分别进行 filtering、normalization、
#       dispersion estimation 和模型拟合。
#   fc_cutoff:
#     判断差异 feature 的 |log2FC| 阈值，默认 1。
#   padj_cutoff:
#     判断差异 feature 的 BH-adjusted P value 阈值，默认 0.05。
#     同时作为 results() independent filtering 的 alpha。
# 返回：
#   list：
#   - difference:
#       每个 comparison 的完整 DESeq2 差异分析结果。
#       independent filtering 或 outlier detection 可能产生 padjust = NA。
#   - norm_count:
#       经过低表达预过滤后，
#       基于 DESeq2 size factor normalization 的 count matrix。
#   - vst:
#       经过低表达预过滤后的 variance stabilized expression matrix。
#       可用于 PCA、heatmap、聚类和样本相关性分析。
#   - summary:
#       每个 comparison 的差异 feature 数量。
# 注意事项：
#   1. rc 必须是真实的非负整数 count。
#      本函数不会通过 round() 将其他类型的数据伪装成 count。
#   2. 函数使用 edgeR::filterByExpr() 进行低表达 pre-filtering，
#      这里只使用 edgeR 的过滤策略，不进行 TMM normalization。
#   3. DESeq2 使用自己的 size factor normalization。
#      不要在输入 DESeq2 前使用 calcNormFactors() 后的 normalized values。
#   4. DESeq2 的差异检验直接基于 raw count model，
#      不使用 norm_count 或 vst 进行差异检验。
#   5. results() 仍会执行 DESeq2 自身的 independent filtering。
#   6. padjust = NA 的 feature 会保留在 difference 中，
#      enriched 自动记为 "none"。
#   7. vst 使用 blind = FALSE，
#      适用于已知 experimental design 后的下游可视化和探索分析。
#   8. pairwise = TRUE 时，各 comparison 独立拟合；
#      但返回的 norm_count 和 vst 仍来自全部样本建立的整体 DESeq2 对象。
#   9. comparisons 中第一个组为 numerator，第二个组为 denominator。

run_deseq2_diff <- function(
    rc, sample_meta, comparisons, sample_col = "sample", group_col = "group",
    pairwise = FALSE, fc_cutoff = 1, padj_cutoff = .05
) {
  metadata <- sample_meta
  
  # 匹配 metadata 与 raw count matrix，并统一样本顺序
  metadata <- metadata |>
    dplyr::filter(.data[[sample_col]] %in% colnames(rc)) |>
    dplyr::arrange(match(.data[[sample_col]], colnames(rc)))
  
  rc <- rc[, metadata[[sample_col]], drop = FALSE] |>
    as.matrix()
  
  # 构建全样本分组因子
  group <- factor(metadata[[group_col]], unique(metadata[[group_col]]))
  
  # 低表达 pre-filter
  # 这里只使用 edgeR 的 filtering strategy，不进行 TMM normalization
  keep <- edgeR::filterByExpr(rc, group = group)
  
  rc_all <- rc[keep, , drop = FALSE]
  
  # 构建 DESeq2 sample metadata
  metadata_all <- data.frame(row.names = metadata[[sample_col]], group = group)
  
  # 创建 DESeqDataSet
  # countData 必须是未经 normalization 的 raw integer count
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = rc_all, colData = metadata_all, design = ~ group
  )
  
  # 完成 size factor、dispersion 和 negative binomial GLM estimation
  des <- DESeq2::DESeq(dds)
  
  if (!pairwise) {
    
    # 所有 comparison 共用同一个 DESeq2 model
    difference <- purrr::map(
      comparisons, ~
        DESeq2::results(
          des, contrast = c("group", .x[1], .x[2]),
          alpha = padj_cutoff # 用实际目标 FDR 优化 independent filtering
        ) |>
        data.frame(check.names = FALSE) |>
        tibble::rownames_to_column("name") |>
        dplyr::rename(
          log2FC  = log2FoldChange,
          padjust = padj
        ) |>
        dplyr::mutate(
          enriched = dplyr::case_when(
            .data[["log2FC"]] >  fc_cutoff & .data[["padjust"]] < padj_cutoff ~ .x[1],
            .data[["log2FC"]] < -fc_cutoff & .data[["padjust"]] < padj_cutoff ~ .x[2],
            TRUE ~ "none"
          )
        )
    )
    
  } else {
    
    # 每个 comparison 独立重新拟合 DESeq2 model
    difference <- purrr::map(
      comparisons, \(x) {
        # 提取当前 comparison 的样本
        metadata_x <- metadata |> dplyr::filter(.data[[group_col]] %in% x)
        
        rc_x <- rc[, metadata_x[[sample_col]], drop = FALSE]
        
        # 重新建立当前 pairwise comparison 的 group factor，
        # 避免 metadata 原 factor 中残留其他 unused levels
        group_x <- factor(metadata_x[[group_col]], levels = x)
        
        # 根据当前两个组重新进行低表达过滤
        keep_x <- edgeR::filterByExpr(rc_x, group = group_x)
        
        rc_x <- rc_x[keep_x, , drop = FALSE]
        
        # 第二个组作为 reference level
        metadata_x <- data.frame(
          row.names = metadata_x[[sample_col]],
          group = factor(metadata_x[[group_col]], levels = c(x[2], x[1]))
        )
        
        # 当前 comparison 独立建立 DESeq2 model
        dds_x <- DESeq2::DESeqDataSetFromMatrix(
          countData = rc_x, colData = metadata_x, design = ~ group
        )
        
        des_x <- DESeq2::DESeq(dds_x)
        
        # 提取 x[1] vs x[2]
        DESeq2::results(
          des_x, contrast = c("group", x[1], x[2]),
          alpha = padj_cutoff
        ) |>
          data.frame(check.names = FALSE) |>
          tibble::rownames_to_column("name") |>
          dplyr::rename(
            log2FC  = .data[["log2FoldChange"]],
            padjust = .data[["padj"]]
          ) |>
          dplyr::mutate(
            enriched = dplyr::case_when(
              .data[["log2FC"]] >  fc_cutoff & .data[["padjust"]] < padj_cutoff ~ x[1],
              .data[["log2FC"]] < -fc_cutoff & .data[["padjust"]] < padj_cutoff ~ x[2],
              TRUE ~ "none"
            )
          )
      }
    )
  }
  
  # 为 difference list 添加 comparison 名称
  difference <- difference |>
    purrr::set_names(
      purrr::map_vec(
        comparisons, ~ paste(.x, collapse = "_vs_")
      )
    )
  
  # DESeq2 size-factor normalized count
  norm_count <- DESeq2::counts(des, normalized = TRUE)
  
  # variance stabilizing transformation
  # 用于 PCA、heatmap、聚类等，不用于差异检验
  vst <- DESeq2::varianceStabilizingTransformation(des, blind = FALSE) |>
    SummarizedExperiment::assay()
  
  # 汇总各 comparison 中的差异 feature 数量
  summary <- purrr::map2_dfr(
    names(difference), comparisons, \(x, y)
    data.frame(
      comparison = x,
      up   = sum(difference[[x]]$enriched == y[1]),
      down = sum(difference[[x]]$enriched == y[2]),
      none = sum(difference[[x]]$enriched == "none")
    )
  )
  
  list(
    difference = difference,
    norm_count = norm_count,
    vst = vst,
    summary = summary
  )
}

#### plot_gsea_barcode ####
#' Plot Gsea Barcode utility
#'
#'
#'
#' @param gseaResult A GSEA result object containing ranked genes and enrichment results.
#' @param set_ID Gene-set identifier selected from a GSEA result.
#' @param bar_color Color specification for `bar_color`.
#' @param bar_width Width of plotted bars.
#' @param core_enrichment Gene identifiers in the leading-edge or core-enrichment subset.
#' @param add_table Whether to append a compact enrichment-statistics table.
#' @return A ggplot-compatible plot object; computed data or fitted objects are retained as attributes when applicable.
#' @export
plot_gsea_barcode <- function(
  gseaResult, set_ID, bar_color = c(up = "#de77ae", down = "#35978f"),
  bar_width = 0.8, core_enrichment = FALSE, add_table = TRUE
) {
  if (!inherits(gseaResult, "gseaResult")) {
    stop("Need a gseaResult object")
  }

  if (!all(set_ID %in% gseaResult@result$ID)) {
    stop("IDs not all in gseaResult object")
  }

  if (length(bar_color) != 2) {
    bar_color <- c(up = "#de77ae", down = "#35978f")
  }

  if (is.null(names(bar_color))) {
    bar_color <- structure(bar_color, names = c("up", "down"))
  }

  if (isFALSE(core_enrichment)) {
    plot_bar <- purrr::map_dfr(
      set_ID, ~
        data.frame(
          TERM = .x,
          GENE = gseaResult@geneSets[[.x]]
        )
    )
  } else {
    plot_bar <- dplyr::filter(gseaResult@result, ID %in% set_ID) |>
      dplyr::group_by(ID) |>
      dplyr::group_modify(~ data.frame(GENE = stringr::str_split_1(.x$core_enrichment, "/"))) |>
      dplyr::ungroup() |>
      dplyr::rename(TERM = ID)
  }

  plot_enriched <- dplyr::filter(gseaResult@result, ID %in% set_ID) |>
    dplyr::select(ID, NES) |>
    dplyr::mutate(enriched = ifelse(NES > 0, "up", "down"))

  rank_df <- data.frame(
    GENE = names(gseaResult@geneList),
    RANK = seq_along(gseaResult@geneList)
  ) |>
    dplyr::left_join(plot_bar, by = "GENE") |>
    dplyr::filter(!is.na(TERM)) |>
    tibble::add_column(value = 1) |>
    dplyr::left_join(dplyr::select(plot_enriched, TERM = ID, enriched), by = "TERM")

  p <- ggplot2::ggplot(
    rank_df, ggplot2::aes(x = RANK, y = value)
  ) +
    ggplot2::geom_bar(
      ggplot2::aes(fill = enriched),
      color = NA, stat = "identity",
      width = bar_width, position = ggplot2::position_identity()
    ) +
    ggplot2::facet_grid(rows = ggplot2::vars(TERM)) +
    ggplot2::scale_fill_manual(values = bar_color) +
    ggplot2::scale_x_continuous(
      limits = c(1, length(gseaResult@geneList)), expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(
      x = "Rank in Ordered Dataset", y = "",
      fill = "Enriched in"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.ticks.length = grid::unit(1.8, "mm"),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = "black", size = 10),
      axis.text.y = ggplot2::element_blank(),
      legend.position = "none",
      panel.spacing.y = grid::unit(0, "cm"),
      strip.text.y = ggplot2::element_text(angle = 0, size = 12, hjust = 0, face = "italic"),
      strip.background.y = ggplot2::element_blank(),
      aspect.ratio = 1 / 10
    )

  if (isTRUE(add_table)) {
    label_df <- dplyr::filter(gseaResult@result, ID %in% set_ID) |>
      dplyr::select(ID, Description, NES, FDR = p.adjust) |>
      dplyr::mutate(
        NES = signif(NES, 4),
        FDR = signif(FDR, 4)
      ) |>
      data.frame(row.names = NULL)

    p_add <- ggpubr::ggtexttable(
      label_df,
      theme = ggpubr::ttheme("blank"), rows = NULL
    ) |>
      ggpubr::tab_add_hline(
        at.row = c(1, 2), row.side = "top",
        linewidth = 3, linetype = 1
      ) |>
      ggpubr::tab_add_hline(
        at.row = nrow(label_df) + 1, row.side = "bottom",
        linewidth = 3, linetype = 1
      )

    p <- cowplot::plot_grid(p, p_add, ncol = 1, align = "v")
  }
  return(p)
}

#### plot_go_bar ####
# GO 富集分析普通横向柱状图
# data: GO enrichment 结果表
# term_col: GO term 名称列，默认 Description
# group_col: GO 分类列，默认 ONTOLOGY
# value_col: 柱子高度列，默认 FoldEnrichment
# sort_col: top_n 筛选前的排序列；NULL 表示按输入顺序
# top_n: 每个 group 保留前 top_n 个 term；NULL 表示不筛选
# sort_decreasing: sort_col 是否降序排序
# group_level: 指定 group 顺序
# group_name: 是否将 BP/CC/MF 改成完整名称
# facet: 是否按 group 分面

#' Plot GO Bar utility
#'
#'
#'
#' @param data An input data frame or compatible object.
#' @param term_col Name of the `term_col` input column.
#' @param group_col Name of the grouping column.
#' @param value_col Name of the `value_col` input column.
#' @param sort_col Name of the `sort_col` input column.
#' @param top_n Number of highest-ranking features or categories retained.
#' @param sort_decreasing Whether sorting is performed in decreasing order.
#' @param group_level Optional order of group levels.
#' @param group_name Named character vector mapping group codes to display labels.
#' @param palette Color palette name or vector supplied to the plot.
#' @param title Optional plot or result title.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param bar_width Width of plotted bars.
#' @param bar_alpha Opacity of plotted bars.
#' @param show_legend Whether to display the plot legend.
#' @param facet Whether enrichment results are separated into facets.
#' @return A ggplot-compatible plot object; computed data or fitted objects are retained as attributes when applicable.
#' @export
plot_go_bar <- function(
  data,
  term_col = "Description", group_col = "ONTOLOGY",
  value_col = "FoldEnrichment", sort_col = NULL,
  top_n = 10, sort_decreasing = FALSE,
  group_level = NULL,
  group_name = c(
    BP = "Biological Process",
    CC = "Cellular Component",
    MF = "Molecular Function"
  ),
  palette = c("#66c2a5", "#fc8d62", "#8da0cb"),
  title = "GO Enrichment Analysis",
  xlab = "FoldEnrichment", ylab = "",
  bar_width = .75, bar_alpha = .85,
  show_legend = FALSE,
  facet = TRUE
) {
  data <- data.frame(data, check.names = FALSE)

  need_col <- c(term_col, group_col, value_col)

  if (!is.null(sort_col)) {
    need_col <- c(need_col, sort_col)
  }

  if (!all(need_col %in% colnames(data))) {
    stop(
      "data should contain columns: ",
      paste(need_col, collapse = " | ")
    )
  }

  ## 1. 整理数据
  data$.input_order <- seq_len(nrow(data))

  if (is.null(sort_col)) {
    data$.sort <- data$.input_order
  } else {
    data$.sort <- data[[sort_col]]
  }

  data <- data |>
    dplyr::transmute(
      term = .data[[term_col]],
      group = .data[[group_col]],
      value = .data[[value_col]],
      .sort = .data[[".sort"]]
    ) |>
    dplyr::filter(!is.na(term), !is.na(group), !is.na(value)) |>
    dplyr::mutate(
      term = as.character(term),
      group = as.character(group),
      value = suppressWarnings(as.numeric(value)),
      .sort = suppressWarnings(as.numeric(.sort))
    ) |>
    dplyr::filter(is.finite(value))

  if (nrow(data) == 0) {
    stop("No valid enrichment terms remained for plotting.")
  }

  ## 2. GO ontology 名称转换
  if (!is.null(group_name)) {
    data <- data |>
      dplyr::mutate(
        group = dplyr::case_when(
          group %in% names(group_name) ~ unname(group_name[group]),
          TRUE ~ group
        )
      )

    if (!is.null(group_level)) {
      group_level <- ifelse(
        group_level %in% names(group_name),
        unname(group_name[group_level]),
        group_level
      )
    }
  }

  if (is.null(group_level)) {
    group_level <- unique(data$group)
  }

  group_level <- group_level[group_level %in% unique(as.character(data$group))]

  if (length(group_level) == 0) {
    stop("No valid group remained for plotting.")
  }

  data <- data |>
    dplyr::mutate(group = factor(group, levels = group_level)) |>
    dplyr::filter(!is.na(group))

  ## 3. 每个 ontology 取前 top_n 个 term
  if (!is.null(top_n)) {
    if (isTRUE(sort_decreasing)) {
      data <- data |>
        dplyr::arrange(group, dplyr::desc(.sort))
    } else {
      data <- data |>
        dplyr::arrange(group, .sort)
    }

    data <- data |>
      dplyr::group_by(group) |>
      dplyr::slice_head(n = top_n) |>
      dplyr::ungroup()
  }

  ## 4. 按 value 排序，用唯一 ID 避免不同 group 中 term 重名
  data <- data |>
    dplyr::arrange(group, value) |>
    dplyr::mutate(
      .term_id = paste0(dplyr::row_number(), "__", term),
      .term_id = factor(.term_id, levels = .term_id)
    )

  term_label <- stats::setNames(data$term, data$.term_id)

  ## 5. 颜色
  palette <- rep(palette, length.out = length(group_level))
  names(palette) <- group_level

  ## 6. 作图
  p <- ggplot2::ggplot(data, ggplot2::aes(x = value, y = .term_id, fill = group)) +
    ggplot2::geom_col(width = bar_width, alpha = bar_alpha, color = NA) +
    ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
    ggplot2::scale_y_discrete(labels = term_label) +
    ggplot2::labs(x = xlab, y = ylab, title = title, fill = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 10, color = "black"),
      axis.title = ggplot2::element_text(size = 12, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = .5, color = "black"),
      axis.ticks.length = grid::unit(2, "mm"),
      axis.line = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(linewidth = .5, color = "black", fill = NA),
      panel.grid.major = ggplot2::element_line(color = "grey88", linewidth = .5),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = .5, size = 13, face = "bold"),
      strip.text = ggplot2::element_text(face = "bold", size = 11),
      strip.background = ggplot2::element_blank(),
      legend.position = ifelse(show_legend, "right", "none")
    )

  if (isTRUE(facet)) {
    p <- p +
      ggplot2::facet_grid(
        rows = ggplot2::vars(group),
        scales = "free_y",
        space = "free_y"
      ) +
      ggplot2::theme(
        panel.spacing.y = grid::unit(1, "mm")
      )
  }

  return(p)
}

#### plot_go_circular_bar ####
# GO 富集分析环状柱状图
# data: GO enrichment 结果表
# term_col: GO term 名称列，默认 Description
# group_col: GO 分类列，默认 ONTOLOGY
# value_col: 柱子高度列，默认 FoldEnrichment
# sort_col: top_n 筛选前的排序列；NULL 表示按输入顺序
# top_n: 每个 group 保留前 top_n 个 term；NULL 表示不筛选
# sort_decreasing: sort_col 是否降序排序
# group_level: 指定 group 顺序
# group_name: 是否将 BP/CC/MF 改成完整名称
# empty_bar: 每个 group 后面添加的空白 bar 数量
# grid_breaks: 环形辅助线刻度；NULL 自动生成
# grid_n: 自动生成 grid_breaks 时的参考数量
# show_grid: 是否显示 grid line
# inner_size: 中间空白区域大小，数值越大中间越空
# label_pad: term 标签距离柱子顶端的距离
# label_wrap: term 标签换行宽度；NULL 表示不换行

#' Plot GO Circular Bar utility
#'
#'
#'
#' @param data An input data frame or compatible object.
#' @param term_col Name of the `term_col` input column.
#' @param group_col Name of the grouping column.
#' @param value_col Name of the `value_col` input column.
#' @param sort_col Name of the `sort_col` input column.
#' @param top_n Number of highest-ranking features or categories retained.
#' @param sort_decreasing Whether sorting is performed in decreasing order.
#' @param group_level Optional order of group levels.
#' @param group_name Named character vector mapping group codes to display labels.
#' @param empty_bar Angular width reserved as an empty separator in the circular bar plot.
#' @param palette Color palette name or vector supplied to the plot.
#' @param grid_breaks Reference values used to draw circular or radar grid lines.
#' @param grid_n Number of grid intervals.
#' @param show_grid Whether to draw panel grid lines.
#' @param inner_size Radius of the empty inner circle in bar-height units.
#' @param label_pad Radial padding between circular bars and their labels.
#' @param label_wrap Maximum label width before circular labels are wrapped.
#' @param bar_width Width of plotted bars.
#' @param bar_alpha Opacity of plotted bars.
#' @param label_size Text size for sample or group labels.
#' @param group_label_size Text size for group labels.
#' @param title Optional plot or result title.
#' @param show_legend Whether to display the plot legend.
#' @return A ggplot-compatible plot object; computed data or fitted objects are retained as attributes when applicable.
#' @export
plot_go_circular_bar <- function(
  data,
  term_col = "Description", group_col = "ONTOLOGY",
  value_col = "FoldEnrichment", sort_col = NULL,
  top_n = 10, sort_decreasing = FALSE,
  group_level = NULL,
  group_name = c(
    BP = "Biological Process",
    CC = "Cellular Component",
    MF = "Molecular Function"
  ),
  empty_bar = 2, palette = c("#66c2a5", "#fc8d62", "#8da0cb"),
  grid_breaks = NULL, grid_n = 4, show_grid = TRUE,
  inner_size = 1, label_pad = NULL, label_wrap = NULL,
  bar_width = .75, bar_alpha = .65,
  label_size = 3, group_label_size = 3,
  title = "GO Enrichment Analysis",
  show_legend = FALSE
) {
  data <- data.frame(data, check.names = FALSE)

  need_col <- c(term_col, group_col, value_col)

  if (!is.null(sort_col)) {
    need_col <- c(need_col, sort_col)
  }

  if (!all(need_col %in% colnames(data))) {
    stop(
      "data should contain columns: ",
      paste(need_col, collapse = " | ")
    )
  }

  ## 1. 整理数据
  data$.input_order <- seq_len(nrow(data))

  if (is.null(sort_col)) {
    data$.sort <- data$.input_order
  } else {
    data$.sort <- data[[sort_col]]
  }

  data <- data |>
    dplyr::transmute(
      term = .data[[term_col]],
      group = .data[[group_col]],
      value = .data[[value_col]],
      .sort = .data[[".sort"]]
    ) |>
    dplyr::filter(!is.na(term), !is.na(group), !is.na(value)) |>
    dplyr::mutate(
      term = as.character(term),
      group = as.character(group),
      value = suppressWarnings(as.numeric(value)),
      .sort = suppressWarnings(as.numeric(.sort))
    ) |>
    dplyr::filter(is.finite(value))

  if (nrow(data) == 0) {
    stop("No valid enrichment terms remained for plotting.")
  }

  ## 2. GO ontology 名称转换
  if (!is.null(group_name)) {
    data <- data |>
      dplyr::mutate(
        group = dplyr::case_when(
          group %in% names(group_name) ~ unname(group_name[group]),
          TRUE ~ group
        )
      )

    ## 如果 group_level 输入的是 BP/CC/MF，也同步转换
    if (!is.null(group_level)) {
      group_level <- ifelse(
        group_level %in% names(group_name),
        unname(group_name[group_level]),
        group_level
      )
    }
  }

  if (is.null(group_level)) {
    group_level <- unique(data$group)
  }

  group_level <- group_level[group_level %in% unique(as.character(data$group))]

  if (length(group_level) == 0) {
    stop("No valid group remained for plotting.")
  }

  data <- data |>
    dplyr::mutate(group = factor(group, levels = group_level)) |>
    dplyr::filter(!is.na(group))

  ## 3. 每个 ontology 取前 top_n 个 term
  if (!is.null(top_n)) {
    if (isTRUE(sort_decreasing)) {
      data <- data |>
        dplyr::arrange(group, dplyr::desc(.sort))
    } else {
      data <- data |>
        dplyr::arrange(group, .sort)
    }

    data <- data |>
      dplyr::group_by(group) |>
      dplyr::slice_head(n = top_n) |>
      dplyr::ungroup()
  }

  ## 4. 组内按 value 排序，并记录顺序
  data <- data |>
    dplyr::arrange(group, value) |>
    dplyr::group_by(group) |>
    dplyr::mutate(
      .order = dplyr::row_number(),
      .empty = FALSE
    ) |>
    dplyr::ungroup()

  ## 5. 添加 empty bar
  if (empty_bar > 0) {
    to_add <- data.frame(
      term = NA_character_,
      group = rep(group_level, each = empty_bar),
      value = NA_real_,
      .sort = NA_real_,
      .order = rep(seq_len(empty_bar), times = length(group_level)) + 1e6,
      .empty = TRUE,
      check.names = FALSE
    )

    to_add$group <- factor(to_add$group, levels = group_level)

    data <- dplyr::bind_rows(data, to_add)
  }

  data <- data |>
    dplyr::arrange(group, .empty, .order) |>
    dplyr::mutate(id = dplyr::row_number())

  ## 6. label 角度
  label_df <- data |>
    dplyr::filter(!.empty, !is.na(term), is.finite(value))

  number_of_bar <- nrow(data)

  label_df <- label_df |>
    dplyr::mutate(
      angle = 90 - 360 * (id - 0.5) / number_of_bar,
      hjust = ifelse(angle < -90, 1, 0),
      angle = ifelse(angle < -90, angle + 180, angle)
    )

  if (!is.null(label_wrap)) {
    label_df$term <- stringr::str_wrap(label_df$term, width = label_wrap)
  }

  ## 7. group baseline 数据
  base_df <- data |>
    dplyr::filter(!.empty) |>
    dplyr::group_by(group) |>
    dplyr::summarise(
      start = min(id),
      end = max(id),
      .groups = "drop"
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(title = mean(c(start, end))) |>
    dplyr::ungroup()

  ## 8. 坐标范围
  max_value <- max(data$value, na.rm = TRUE)

  if (!is.finite(max_value) || max_value <= 0) {
    stop("value_col should contain positive numeric values.")
  }

  if (is.null(label_pad)) {
    label_pad <- max_value * .05
  }

  if (is.null(grid_breaks)) {
    grid_breaks <- pretty(c(0, max_value), n = grid_n)
    grid_breaks <- grid_breaks[grid_breaks > 0 & grid_breaks < max_value]
  }

  y_min <- -max_value * inner_size
  y_max <- max_value + label_pad * 4

  ## 9. grid line 数据
  ## 只在每个 group 后面的 empty_bar 区域画 grid
  ## 最后一个 group 后面的 empty_bar 也保留

  grid_df <- base_df |>
    dplyr::mutate(
      gap_start = end + 1,
      gap_end = dplyr::lead(start) - 1
    )

  grid_df$gap_end[nrow(grid_df)] <- nrow(data)

  grid_df <- grid_df |>
    dplyr::filter(gap_end >= gap_start)

  ## 10. 颜色
  palette <- rep(palette, length.out = length(group_level))
  names(palette) <- group_level

  ## 11. 作图
  p <- ggplot2::ggplot(data, ggplot2::aes(x = id, y = value, fill = group))

  ## grid lines 放在前面，让柱子盖在上面
  ## 这里只在 empty_bar 空白区域画 grid

  if (isTRUE(show_grid) && nrow(grid_df) > 0 && length(grid_breaks) > 0) {
    for (g in grid_breaks) {
      p <- p +
        ggplot2::geom_segment(
          data = grid_df, ggplot2::aes(x = gap_start, xend = gap_end),
          y = g, yend = g, colour = "grey80", linewidth = .3,
          linetype = "longdash", inherit.aes = FALSE
        )
    }

    p <- p +
      ggplot2::annotate(
        "text",
        x = max(data$id), y = grid_breaks, label = grid_breaks,
        color = "grey50", size = 3, fontface = "bold", hjust = 1
      )
  }

  p <- p +
    ggplot2::geom_col(width = bar_width, alpha = bar_alpha, color = NA, na.rm = TRUE) +
    ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
    ggplot2::scale_y_continuous(limits = c(y_min, y_max), expand = c(0, 0)) +
    ggplot2::coord_polar() +
    ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(
        x = id, y = value + label_pad,
        label = term, hjust = hjust, angle = angle
      ),
      color = "black", size = label_size, alpha = .75,
      fontface = "plain", inherit.aes = FALSE,
      na.rm = TRUE
    ) +
    ggplot2::geom_segment(
      data = base_df,
      ggplot2::aes(x = start, y = y_min * .08, xend = end, yend = y_min * .08),
      colour = "black", alpha = .8,
      linewidth = .6, inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = base_df,
      ggplot2::aes(x = title, y = y_min * .22, label = group),
      colour = "black", alpha = .8, size = group_label_size,
      fontface = "bold", inherit.aes = FALSE
    ) +
    ggplot2::labs(title = title, fill = "") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = ifelse(show_legend, "right", "none"),
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = .5, size = 13, face = "bold"),
      plot.margin = grid::unit(c(5, 5, 5, 5), "mm")
    )

  return(p)
}
