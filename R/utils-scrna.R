#' Build a local CellMarker reference database
#'
#' Build a named-list reference database from a CellMarker-style table and save
#' it as an RDS file. The input is filtered to the requested organism and tissue
#' after names are normalized to lower snake case.
#'
#' 从 CellMarker 格式的表格构建命名基因集列表，并保存为 RDS 文件。物种、
#' 组织和细胞名会先统一为小写 snake case，再按 `org` 和 `tissue` 筛选。
#'
#' @param data A data frame containing `species`, `tissue_class`, `cell_name`,
#'   `marker`, `Symbol`, and `GeneID` columns.
#' @param tissue Tissue name used to filter the input and name the RDS file.
#' @param org Organism name used to filter the input and name the RDS file.
#' @param database Directory in which the database file is saved.
#'
#' @return Invisibly returns the path of the saved RDS file.
#' @export
scRNA_marker_db_build <- function(
  data, tissue = "intestine", org = "mouse",
  database = .mengR_db_file("CellMarker-2.0", "cell_marker_rds")
) {
  required_cols <- c(
    "species", "tissue_class", "cell_name", "marker", "Symbol", "GeneID"
  )
  data_df <- .as_df(data)
  .check_columns(data_df, required_cols, "data")

  normalize_name <- function(x) {
    x |>
      as.character() |>
      stringr::str_to_lower() |>
      stringr::str_replace_all("[^[:alnum:]]+", "_") |>
      stringr::str_replace_all("^_|_$", "")
  }

  tissue <- normalize_name(tissue)
  org <- normalize_name(org)

  # 1. 统一分类名称，并只保留目标物种和组织
  data_df <- data_df |>
    dplyr::transmute(
      species = normalize_name(.data$species),
      tissue_class = normalize_name(.data$tissue_class),
      cell_name = normalize_name(.data$cell_name),
      marker = as.character(.data$marker),
      symbol = as.character(.data$Symbol),
      gene_id = as.character(.data$GeneID)
    ) |>
    dplyr::filter(.data$species == org, .data$tissue_class == tissue) |>
    dplyr::distinct()

  if (nrow(data_df) == 0) {
    stop('No records match org = "', org, '" and tissue = "', tissue, '".')
  }

  # 2. 以 species|tissue|cell 为名称，分别保存三种基因标识
  key_vec <- paste(
    data_df$species, data_df$tissue_class, data_df$cell_name,
    sep = "|"
  )
  split_df <- split(data_df, key_vec)
  clean_values <- function(x) unique(x[!is.na(x) & nzchar(x)])

  database_list <- list(
    marker = lapply(split_df, \(x) clean_values(x$marker)),
    symbol = lapply(split_df, \(x) clean_values(x$symbol)),
    gene_id = lapply(split_df, \(x) clean_values(x$gene_id))
  )
  database_list <- lapply(database_list, \(x) x[lengths(x) > 0])

  # 3. 创建目录并写出数据库
  if (!dir.exists(database)) {
    dir.create(database, recursive = TRUE, showWarnings = FALSE)
  }
  database_file <- file.path(database, paste0(org, ".", tissue, ".rds"))
  saveRDS(database_list, database_file)
  invisible(database_file)
}

#' Match marker genes against a local CellMarker database
#'
#' Compare an input marker vector with every cell-type gene set in a database
#' created by [scRNA_marker_db_build()]. Results include overlap counts, input
#' coverage, reference coverage, and matched genes.
#'
#' 将输入 marker 与 [scRNA_marker_db_build()] 构建的各细胞类型基因集比对，
#' 返回交集数量、输入覆盖率、参考集覆盖率和匹配基因。
#'
#' @param markers Character vector of input marker genes.
#' @param tissue Tissue name used in the database file name.
#' @param org Organism name used in the database file name.
#' @param gene_type Gene identifier type: `"symbol"`, `"marker"`, or
#'   `"gene_id"`.
#' @param database Directory containing CellMarker RDS files.
#'
#' @return A data frame sorted by decreasing input-marker coverage.
#' @export
scRNA_marker_match <- function(
  markers, tissue = "intestine", org = "mouse",
  gene_type = c("symbol", "marker", "gene_id"),
  database = .mengR_db_file("CellMarker-2.0", "cell_marker_rds")
) {
  gene_type <- match.arg(gene_type)
  markers <- unique(as.character(markers))
  markers <- markers[!is.na(markers) & nzchar(markers)]
  if (length(markers) == 0) stop("markers should contain at least one gene.")

  database_file <- file.path(database, paste0(org, ".", tissue, ".rds"))
  if (!file.exists(database_file)) {
    stop("CellMarker database file not found: ", database_file)
  }

  database_list <- readRDS(database_file)
  if (!gene_type %in% names(database_list)) {
    stop("gene_type is not available in the database: ", gene_type)
  }
  reference_list <- database_list[[gene_type]]

  # 逐个参考基因集计算交集及双向覆盖率
  result_df <- purrr::imap_dfr(reference_list, \(reference_vec, item) {
    matched_vec <- intersect(markers, reference_vec)
    item_vec <- strsplit(item, "\\|")[[1]]
    data.frame(
      item = item,
      species = item_vec[1],
      tissue = item_vec[2],
      cell = paste(item_vec[-c(1, 2)], collapse = "|"),
      input_n = length(markers),
      reference_n = length(reference_vec),
      matched_n = length(matched_vec),
      input_coverage = length(matched_vec) / length(markers) * 100,
      reference_coverage = length(matched_vec) / length(reference_vec) * 100,
      matched_genes = paste(matched_vec, collapse = ","),
      check.names = FALSE
    )
  })

  result_df |>
    dplyr::arrange(
      dplyr::desc(.data$input_coverage),
      dplyr::desc(.data$matched_n)
    )
}

#' Plot grouped marker expression on a Seurat embedding
#'
#' Create one row of Seurat feature plots for each element of a named marker
#' list. Missing genes are removed, and rows are padded to equal widths.
#'
#' 对命名 marker 列表中的每一组基因绘制一行 Seurat FeaturePlot。函数会
#' 删除不存在的基因，并用空白图保持各行宽度一致。
#'
#' @param seurat A Seurat object.
#' @param markers A named list of marker vectors.
#' @param reduction Dimensional reduction passed to `Seurat::FeaturePlot()`.
#'
#' @return A combined ggplot/cowplot object.
#' @export
scRNA_feature_plots <- function(seurat, markers, reduction = "umap") {
  if (!is.list(markers) || is.null(names(markers))) {
    stop("markers should be a named list of character vectors.")
  }

  marker_list <- lapply(markers, \(x) intersect(as.character(x), rownames(seurat)))
  marker_list <- marker_list[lengths(marker_list) > 0]
  if (length(marker_list) == 0) {
    stop("None of the marker genes were found in the Seurat object.")
  }
  max_n <- max(lengths(marker_list))

  # 每一组 marker 构建一行：组名、基因图和补齐用空白图
  row_plot_list <- purrr::imap(marker_list, \(gene_vec, group_name) {
    title_plot <- ggpubr::ggtext(
      data.frame(x = 1, y = 1), "x", "y",
      size = 14, face = "bold",
      label = stringr::str_to_title(gsub("_", " ", group_name))
    ) + ggplot2::theme_void()
    feature_plot_list <- lapply(gene_vec, \(gene) {
      Seurat::FeaturePlot(seurat, reduction = reduction, features = gene) +
        ggplot2::theme(aspect.ratio = 1)
    })
    blank_plot <- ggplot2::ggplot() +
      ggplot2::theme_void()
    cowplot::plot_grid(
      plotlist = c(
        list(title_plot), feature_plot_list,
        rep(list(blank_plot), max_n - length(gene_vec))
      ),
      nrow = 1
    )
  })

  cowplot::plot_grid(plotlist = row_plot_list, align = "v", ncol = 1)
}

#' Plot marker-gene violin plots from a Seurat object
#'
#' Reshape selected marker expression and draw vertically faceted violin plots.
#' Markers may be supplied as a vector or a named list; list names are included
#' in facet labels.
#'
#' 整理 Seurat 对象中的 marker 表达量，绘制纵向分面小提琴图。`markers` 可以
#' 是字符向量或命名列表；命名列表的组名会写入分面标签。
#'
#' @param seurat A Seurat object.
#' @param markers A character vector or a named list of marker vectors.
#' @param group_col Metadata column used for the x-axis. If `NULL`, active
#'   identities are used.
#' @param group_level Optional order of cell groups.
#' @param group_colors Optional named color vector for cell groups.
#'
#' @return A ggplot object.
#' @export
scRNA_violin_plots <- function(
  seurat, markers, group_col = NULL, group_level = NULL,
  group_colors = NULL
) {
  is_grouped <- is.list(markers) && !is.null(names(markers))
  marker_list <- if (is.list(markers)) markers else list(markers)
  marker_list <- lapply(
    marker_list, \(x) intersect(as.character(x), rownames(seurat))
  )
  marker_list <- marker_list[lengths(marker_list) > 0]
  if (length(marker_list) == 0) {
    stop("None of the marker genes were found in the Seurat object.")
  }
  gene_vec <- unique(unlist(marker_list, use.names = FALSE))

  # 1. 提取表达量，并按 FetchData 的细胞顺序对齐分组
  expression_df <- Seurat::FetchData(seurat, vars = gene_vec)
  if (is.null(group_col)) {
    cell_group <- as.character(Seurat::Idents(seurat)[rownames(expression_df)])
  } else {
    if (!group_col %in% colnames(seurat@meta.data)) {
      stop("group_col was not found in seurat@meta.data: ", group_col)
    }
    cell_group <- as.character(
      seurat@meta.data[rownames(expression_df), group_col, drop = TRUE]
    )
  }
  if (is.null(group_level)) group_level <- unique(cell_group)

  plot_df <- expression_df |>
    tibble::rownames_to_column("cell_id") |>
    dplyr::mutate(group = factor(cell_group, levels = group_level)) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(gene_vec), names_to = "gene", values_to = "value"
    )

  # 2. 如果 marker 按组输入，则在分面标签中保留组名
  gene_label <- gene_vec
  if (is_grouped) {
    gene_group <- rep(names(marker_list), lengths(marker_list))
    names(gene_group) <- unlist(marker_list, use.names = FALSE)
    gene_label <- paste0("(", gene_group[gene_vec], ") ", gene_vec)
  }
  label_map <- stats::setNames(gene_label, gene_vec)
  plot_df$gene_label <- factor(label_map[plot_df$gene], levels = gene_label)

  group_colors <- .resolve_group_colors(group_level, group_colors)

  ggplot2::ggplot(plot_df, ggplot2::aes(.data$group, .data$value)) +
    ggplot2::geom_violin(ggplot2::aes(fill = .data$group), scale = "width") +
    ggplot2::scale_fill_manual(values = group_colors) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::facet_grid(.data$gene_label ~ ., scales = "free_y") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 45, hjust = 1, color = "black", size = 14
      ),
      axis.text.y = ggplot2::element_blank(),
      legend.position = "none",
      panel.spacing.y = grid::unit(0, "cm"),
      strip.text.y = ggplot2::element_text(
        angle = 0, size = 14, hjust = 0, face = "italic"
      ),
      strip.background.y = ggplot2::element_blank()
    )
}

#' Run a standard CellChat analysis from a Seurat object
#'
#' Prepare an RNA assay, create a CellChat object, infer ligand-receptor
#' communication probabilities, aggregate pathway networks, and calculate
#' sender/receiver centrality. The input Seurat object is copied by R semantics;
#' the returned value is the completed CellChat object.
#'
#' 从 Seurat 对象整理 RNA assay，创建 CellChat 对象，依次完成配体-受体
#' 通讯概率、通路聚合网络和 sender/receiver 中心性分析。
#'
#' @param seurat A Seurat object.
#' @param sample_name A single sample or dataset label stored in CellChat
#'   metadata.
#' @param celltype_col Metadata column containing cell-type labels.
#' @param species Ligand-receptor database species: `"mouse"` or `"human"`.
#' @param assay Assay used for CellChat analysis.
#' @param min_cells Minimum number of cells required by
#'   `CellChat::filterCommunication()`.
#' @param population_size Whether communication probabilities should account
#'   for the size of each cell population.
#' @param prob_method Average-expression method passed as `type` to
#'   `CellChat::computeCommunProb()`.
#' @param trim Trim fraction used when `prob_method = "truncatedMean"`.
#'
#' @return A fully processed CellChat object.
#' @export
cellchat_run <- function(
  seurat, sample_name, celltype_col = "celltype",
  species = c("mouse", "human"), assay = "RNA", min_cells = 10,
  population_size = FALSE,
  prob_method = c("truncatedMean", "triMean", "thresholdedMean", "median"),
  trim = 0.1
) {
  species <- match.arg(species)
  prob_method <- match.arg(prob_method)
  if (length(sample_name) != 1 || is.na(sample_name) || !nzchar(sample_name)) {
    stop("sample_name should be one non-empty value.")
  }
  if (!is.numeric(min_cells) || length(min_cells) != 1 || min_cells < 1) {
    stop("min_cells should be a positive integer.")
  }
  if (!assay %in% names(seurat@assays)) {
    stop("assay was not found in the Seurat object: ", assay)
  }

  # 1. 合并 Seurat v5 分层，并在需要时创建 log-normalized data layer
  Seurat::DefaultAssay(seurat) <- assay
  layer_vec <- Seurat::Layers(seurat[[assay]])
  if (any(grepl("^(counts|data)[.]", layer_vec))) {
    seurat <- Seurat::JoinLayers(seurat, assay = assay)
    layer_vec <- Seurat::Layers(seurat[[assay]])
  }
  if (!"data" %in% layer_vec) {
    seurat <- Seurat::NormalizeData(
      seurat,
      assay = assay, normalization.method = "LogNormalize",
      scale.factor = 10000, verbose = FALSE
    )
  }

  expression_mat <- Seurat::GetAssayData(
    seurat,
    assay = assay, layer = "data"
  )
  metadata_df <- seurat@meta.data[colnames(expression_mat), , drop = FALSE]
  .check_columns(metadata_df, celltype_col, "seurat@meta.data")
  label_vec <- as.character(metadata_df[[celltype_col]])
  keep_cell <- !is.na(label_vec) & nzchar(label_vec)
  if (!any(keep_cell)) stop("No cells have valid cell-type labels.")

  expression_mat <- expression_mat[, keep_cell, drop = FALSE]
  metadata_df <- metadata_df[keep_cell, , drop = FALSE]
  metadata_df$labels <- droplevels(factor(label_vec[keep_cell]))
  metadata_df$samples <- factor(sample_name)

  # 2. 创建 CellChat 对象并选择物种数据库
  cellchat <- CellChat::createCellChat(
    object = expression_mat, meta = metadata_df, group.by = "labels"
  )
  cellchat@DB <- if (species == "mouse") {
    CellChat::CellChatDB.mouse
  } else {
    CellChat::CellChatDB.human
  }

  # 3. 识别高表达交互，计算通讯概率并按通路聚合
  cellchat <- CellChat::subsetData(cellchat)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::computeCommunProb(
    cellchat,
    type = prob_method, trim = trim,
    population.size = population_size
  )
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = min_cells)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)
  CellChat::netAnalysis_computeCentrality(cellchat, slot.name = "netP")
}

#' Plot a CellChat interaction profile as a radar chart
#'
#' Extract outgoing or incoming interaction counts/weights for one cell group
#' and display the values across partner groups as a single-series radar chart.
#'
#' 提取某一细胞群的 outgoing 或 incoming 交互数量/强度，并将各互作
#' 细胞群的数值绘制为单序列雷达图。
#'
#' @param cellchat A processed CellChat object.
#' @param cell Cell-group name present in the selected network matrix.
#' @param type Network value type: `"count"` or `"weight"`.
#' @param direction Use a matrix row for `"outgoing"` or a column for
#'   `"incoming"` interactions.
#' @param title Optional plot title.
#' @param color Line and point color.
#' @param grid_min,grid_max Radar-grid limits. `grid_max = NULL` derives a limit
#'   from the data.
#' @param axis_label_size,grid_label_size Text sizes passed to `ggradar`.
#'
#' @return A ggplot object returned by `ggradar::ggradar()`.
#' @export
#' @param grid_max Upper limit of the radar grid; `NULL` derives it from the data.
#' @param grid_label_size Numeric setting for `grid_label_size`.
cellchat_radar <- function(
  cellchat, cell, type = c("count", "weight"),
  direction = c("outgoing", "incoming"), title = NULL,
  color = "#4DBBD5", grid_min = 0, grid_max = NULL,
  axis_label_size = 3, grid_label_size = 3
) {
  type <- match.arg(type)
  direction <- match.arg(direction)
  network_mat <- cellchat@net[[type]]
  if (is.null(network_mat)) stop("The CellChat network does not contain: ", type)

  # 行表示 outgoing，列表示 incoming
  if (direction == "outgoing") {
    if (!cell %in% rownames(network_mat)) stop("Unknown sender cell group: ", cell)
    value_vec <- network_mat[cell, ]
  } else {
    if (!cell %in% colnames(network_mat)) stop("Unknown receiver cell group: ", cell)
    value_vec <- network_mat[, cell]
  }
  value_vec <- sort(as.numeric(value_vec), decreasing = TRUE)
  names(value_vec) <- if (direction == "outgoing") {
    colnames(network_mat)[order(network_mat[cell, ], decreasing = TRUE)]
  } else {
    rownames(network_mat)[order(network_mat[, cell], decreasing = TRUE)]
  }
  radar_df <- data.frame(
    group = type, as.list(value_vec), check.names = FALSE
  )

  if (is.null(grid_max)) {
    grid_max <- signif(max(value_vec, na.rm = TRUE) * 1.05, 2)
    if (!is.finite(grid_max) || grid_max <= grid_min) grid_max <- grid_min + 1
  }
  if (grid_max <= grid_min) stop("grid_max should be larger than grid_min.")
  grid_mid <- (grid_min + grid_max) / 2
  if (is.null(title)) {
    metric <- if (type == "count") "interaction number" else "interaction strength"
    title <- paste(cell, direction, metric, sep = ": ")
  }

  ggradar::ggradar(
    radar_df,
    grid.min = grid_min,
    grid.mid = grid_mid,
    grid.max = grid_max,
    values.radar = signif(c(grid_min, grid_mid, grid_max), 2),
    group.colours = color,
    group.line.width = 0.7,
    group.point.size = 2.2,
    background.circle.colour = "white",
    gridline.min.colour = "grey88",
    gridline.mid.colour = "grey82",
    gridline.max.colour = "grey70",
    grid.line.width = 0.4,
    axis.label.offset = 1.1,
    axis.line.colour = "grey82",
    axis.label.size = axis_label_size,
    grid.label.size = grid_label_size,
    legend.position = "none"
  ) +
    ggplot2::labs(title = title) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5)
    )
}

#' Summarize selected genes across single-cell groups
#'
#' For each selected gene and metadata group, calculate the number and
#' percentage of expressing cells plus mean count and normalized expression in
#' all cells and positive cells.
#'
#' 针对每个基因和 metadata 分组，计算阳性细胞数及比例，以及全部细胞和
#' 阳性细胞中的平均 count 与标准化表达量。
#'
#' @param seurat A Seurat object.
#' @param group_cols Metadata columns defining groups.
#' @param genes Character vector of genes. Genes absent from the object are
#'   omitted.
#' @param assay Assay containing count and normalized-expression layers.
#' @param count_layer Layer containing raw counts.
#' @param expression_layer Layer containing normalized expression.
#'
#' @return A long-format tibble with one row per group-gene combination.
#' @export
scRNA_summarise_genes <- function(
  seurat, group_cols, genes, assay = "RNA", count_layer = "counts",
  expression_layer = "data"
) {
  metadata_df <- seurat@meta.data |>
    tibble::rownames_to_column("cell")
  .check_columns(metadata_df, group_cols, "seurat@meta.data")
  gene_vec <- intersect(unique(as.character(genes)), rownames(seurat[[assay]]))
  if (length(gene_vec) == 0) {
    stop("None of genes were found in the selected assay.")
  }

  count_mat <- Seurat::GetAssayData(
    seurat,
    assay = assay, layer = count_layer
  )[gene_vec, , drop = FALSE]
  expression_mat <- Seurat::GetAssayData(
    seurat,
    assay = assay, layer = expression_layer
  )[gene_vec, , drop = FALSE]
  metadata_df <- metadata_df |>
    dplyr::select("cell", dplyr::all_of(group_cols))

  # 逐基因整理单细胞数值，再按用户指定的 metadata 字段汇总
  purrr::map_dfr(gene_vec, \(gene) {
    gene_df <- data.frame(
      cell = colnames(count_mat),
      gene = gene,
      count = as.numeric(Matrix::colSums(count_mat[gene, , drop = FALSE])),
      expression = as.numeric(
        Matrix::colSums(expression_mat[gene, , drop = FALSE])
      ),
      check.names = FALSE
    ) |>
      dplyr::mutate(positive = .data$count > 0) |>
      dplyr::left_join(metadata_df, by = "cell")

    gene_df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols)), .data$gene) |>
      dplyr::summarise(
        n_cells = dplyr::n(),
        positive_cells = sum(.data$positive),
        percent = mean(.data$positive) * 100,
        mean_count = mean(.data$count),
        mean_count_positive = ifelse(
          any(.data$positive), mean(.data$count[.data$positive]), 0
        ),
        mean_expression = mean(.data$expression),
        mean_expression_positive = ifelse(
          any(.data$positive), mean(.data$expression[.data$positive]), 0
        ),
        .groups = "drop"
      )
  })
}

#' Plot summarized single-cell gene expression as dot plots
#'
#' Draw one dot-plot facet per gene from the output of
#' [scRNA_summarise_genes()] or a compatible table. Dot size represents the
#' positive-cell percentage and color represents mean expression by default.
#'
#' 使用 [scRNA_summarise_genes()] 的输出或兼容数据绘图。默认点大小表示
#' 阳性细胞比例，颜色表示平均标准化表达量，每个基因独立分面。
#'
#' @param data A data frame containing grouping, gene, percentage, and
#'   expression columns.
#' @param x_col,y_col Columns mapped to the x and y axes.
#' @param gene_col Column used for facets.
#' @param percent_col Column mapped to point size.
#' @param expression_col Column mapped to point color.
#' @param title Optional plot title.
#' @param facet_nrow Number of facet rows.
#' @param color_option Viridis palette option from `"A"` through `"H"`.
#'
#' @return A ggplot object.
#' @export
#' @param y_col Name of the `y_col` input column.
scRNA_expression_dotplot <- function(
  data, x_col, y_col, gene_col = "gene", percent_col = "percent",
  expression_col = "mean_expression", title = NULL, facet_nrow = 1,
  color_option = LETTERS[1:8]
) {
  color_option <- match.arg(color_option)
  plot_df <- .as_df(data)
  .check_columns(
    plot_df, c(x_col, y_col, gene_col, percent_col, expression_col), "data"
  )

  ggplot2::ggplot(
    plot_df, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]])
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = .data[[percent_col]], color = .data[[expression_col]]
      )
    ) +
    ggplot2::facet_wrap(
      stats::as.formula(paste("~", gene_col)),
      nrow = facet_nrow
    ) +
    ggplot2::scale_color_viridis_c(option = color_option) +
    ggplot2::labs(
      x = NULL, y = NULL, size = "Positive (%)",
      color = "Mean expression", title = title
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}
