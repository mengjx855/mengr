#### Jin-Xin Meng, jinxmeng@zju.edu.cn, 20260911, 20260923 ####

# 20260911: create new file for single-cell spatial RNA-seq analysis functions.
# 20260911: add 'strna_rank_roi()' to rank spatial ROIs by target-cell enrichment.
# 20260911: add 'strna_roi_plot()' to plot spatial ROIs with cell polygons.
# 20260911: add 'strna_overview_plot()' to plot spatial overview with cell centroids.
# 20260916: correct ROI indexing and boundary handling, add validation, and document all spatial functions.
# 20260923: correct Roxygen export tags and the window-height variable reference.


#### strna_rank_roi ####
#' Rank spatial ROIs by target-cell enrichmen
#'
#' Scan rectangular windows across a spatial cell-centroid table and rank
#' candidate ROIs according to enrichment of one or more target cell types.
#'
#' Enrichment is calculated relative to the cell-type fraction in the whole
#' spatial image:
#'
#'   expected = n_cells_in_roi * background_fraction
#'
#'   log2_enrichment =
#'     log2((observed + pseudocount) / (expected + pseudocount))
#'
#' For multiple target cell types, enrichment scores can be combined using
#' either their mean or minimum value.
#'
#' @param data A data frame containing exactly one row per cell.
#' @param cell_col Column containing unique cell IDs.
#' @param x_col,y_col Numeric columns containing cell-centroid coordinates.
#' @param celltype_col Column containing cell-type annotations.
#' @param target_celltype Character vector specifying one or more targe
#'   cell types.
#' @param window_width,window_height Width and height of each sliding ROI.
#' @param step_x,step_y Sliding step along the x and y axes.
#' @param min_cells Minimum total number of cells required in an ROI.
#' @param min_target_cells Minimum number of each target cell type required
#'   in an ROI.
#' @param combine_method Method used to combine enrichment scores when multiple
#'   target cell types are supplied. "mean" uses the mean enrichment, while
#'   "min" uses the weakest enrichment and is stricter for co-enrichment.
#' @param pseudocount Pseudocount added to observed and expected cell counts.
#' @param show_top_roi Non-negative integer giving the number of top-ranked ROIs
#'   for which comma-separated cell IDs are retained in the output.
#'
#' @details Candidate windows start at the minimum observed coordinate and move
#'   by `step_x` and `step_y`. Only full windows that fit inside the observed
#'   coordinate range are evaluated. Cells on a window boundary are included
#'   in every candidate window whose closed rectangle contains that cell.
#'
#' @return A data frame containing ROI coordinates, target-cell counts,
#'   fractions, expected counts, enrichment scores, and ROI rank.
#'
#' @export
strna_rank_roi <- function(
  data, cell_col = "cell_id", x_col = "X", y_col = "Y", celltype_col = "celltype",
  target_celltype, window_width = 700, window_height = 700, step_x = 200,
  step_y = 200, min_cells = 100, min_target_cells = 1,
  combine_method = c("mean", "min"), pseudocount = 0.5, show_top_roi = 20
) {
  combine_method <- match.arg(combine_method)

  if (!is.data.frame(data)) {
    stop("data must be a data frame.")
  }

  # 1. Check inpu
  required_cols <- c(cell_col, x_col, y_col, celltype_col)

  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!is.numeric(data[[x_col]]) || !is.numeric(data[[y_col]])) {
    stop("x_col and y_col must be numeric.")
  }

  size_vec <- c(window_width, window_height, step_x, step_y)
  if (length(size_vec) != 4L || any(!is.finite(size_vec)) || any(size_vec <= 0)) {
    stop("window_width, window_height, step_x and step_y must be > 0.")
  }

  count_vec <- c(min_cells, min_target_cells, show_top_roi)
  if (any(!is.finite(count_vec)) || any(count_vec < 0) || any(count_vec %% 1 != 0)) {
    stop("min_cells, min_target_cells and show_top_roi must be non-negative integers.")
  }
  if (min_cells < 1L) {
    stop("min_cells must be at least 1.")
  }
  if (length(pseudocount) != 1L || !is.finite(pseudocount) || pseudocount <= 0) {
    stop("pseudocount must be one positive finite number.")
  }

  target_celltype <- unique(as.character(target_celltype))
  target_celltype <- target_celltype[!is.na(target_celltype) & nzchar(target_celltype)]
  if (!length(target_celltype)) {
    stop("target_celltype must contain at least one non-missing cell type.")
  }

  all_celltype <- unique(as.character(data[[celltype_col]]))

  missing_target <- setdiff(target_celltype, all_celltype)

  if (length(missing_target) > 0) {
    stop("Target cell types not found: ", paste(missing_target, collapse = ", "))
  }

  if (anyNA(data[[cell_col]]) || anyDuplicated(data[[cell_col]])) {
    stop("cell_col must contain unique, non-missing cell IDs.")
  }

  # 2. Prepare centroid table
  cell_df <- data.table::data.table(
    cell_id = data[[cell_col]], x = data[[x_col]], y = data[[y_col]],
    celltype = as.character(data[[celltype_col]])
  )

  cell_df <- cell_df[is.finite(x) & is.finite(y) & !is.na(celltype)]

  if (nrow(cell_df) == 0L) {
    stop("No valid cells remain after removing missing coordinates.")
  }

  # 3. Background abundance in the whole image
  # Example:
  #   Eosinophil = 300 / 10000 = 0.03
  # Thus, the expected Eosinophil count in an ROI containing 500 cells is:
  #   500 * 0.03 = 15
  #
  # 这里首先计算“整张图”的 cell type 组成。
  # .N 是 data.table 中非常常用的特殊符号，表示当前分组中有多少行
  background_df <- cell_df[, .(background_n = .N), by = celltype]
  background_df[, background_fraction := background_n / sum(background_n)]
  background_df <- background_df[celltype %in% target_celltype]

  # 4. Build sliding-window grid
  # x_start 保存“所有候选 ROI 的左边界坐标”
  # 假设：
  # min(x) = 100，window_width = 700，step_x = 200
  # 那么可能得到：# 100, 300, 500, 700, 900, ...
  # max(dt$x) - window_width
  # 是为了避免最后一个窗口超过整张图的右边界
  x_range <- range(cell_df$x)
  y_range <- range(cell_df$y)
  if (diff(x_range) < window_width || diff(y_range) < window_height) {
    stop("The requested window is larger than the spatial image.")
  }
  x_start <- seq(x_range[[1L]], x_range[[2L]] - window_width, by = step_x)

  # y_start 同理，保存每个候选 ROI 的上边界/起始 y 坐标
  # 对图像坐标来说：y 越大通常越靠下，但这里只是在数值坐标上进行运算，不影响计算
  y_start <- seq(y_range[[1L]], y_range[[2L]] - window_height, by = step_y)

  # 如果连一个完整窗口都放不进去，
  # 说明用户设置的 window_width / height 比整张图还大。
  if (length(x_start) == 0 || length(y_start) == 0) {
    stop("The requested window is larger than the spatial image.")
  }

  # CJ = Cross Join，即笛卡尔积。
  # 假设 x 有 3 个候选位置：ix = 1, 2, 3
  # y 有 2 个候选位置：iy = 1, 2
  # 那么二维 ROI 一共有：
  # 3 × 2 = 6
  # 即：(1,1), (1,2), (2,1), (2,2), (3,1), (3,2)
  grid_df <- data.table::CJ(ix = seq_along(x_start), iy = seq_along(y_start))

  # 根据 ix / iy，把“窗口索引”转换成实际空间坐标。
  # 一个 ROI 最终由四个边界定义：
  # x_min -------- x_max
  #   |              |
  #   |     ROI      |
  #   |              |
  # y_min -------- y_max
  # 注意对于图像坐标：
  # y_min 在视觉上通常更靠上，
  # y_max 在视觉上通常更靠下。
  grid_df[, `:=`(
    x_min = x_start[ix],
    x_max = x_start[ix] + window_width,
    y_min = y_start[iy],
    y_max = y_start[iy] + window_height
  )
  ]

  # .I 表示当前行号，因此每个候选窗口获得一个唯一 roi_id：ROI 1, ROI 2, ROI 3, ...
  grid_df[, roi_id := .I]

  # nx = x 方向一共有多少个候选窗口位置
  nx <- length(x_start)
  # ny = y 方向一共有多少个候选窗口位置
  ny <- length(y_start)

  # 4. Assign cells to candidate windows
  # A cell does not need to be compared with every ROI. Because the windows
  # lie on a regular grid, a cell can only belong to a small number of nearby
  # windows.
  #
  # Example: window = 700 um, step = 200 um
  #
  # A cell can belong to at most approximately: ceiling(700 / 200)^2 = 16
  # windows, which is substantially faster than scanning all cells separately
  # for every ROI.
  #
  # 因为窗口的位置是规则的，所以我们可以直接算出“哪些窗口可能覆盖这个细胞”，
  # 没必要把这个细胞和所有窗口逐一比较。
  # 原来：一个 cell 需要 ROI1, ROI2, ... ROI1000？全部检查
  # 优化：一个 cell 根据坐标直接知道 “只有附近这十几个 ROI 有可能包含我” 只检查这些

  # overlap_x 表示：
  # 在 x 方向上，一个 cell 最多可能同时被多少个重叠窗口覆盖。
  # 例如：
  # window_width = 700, step_x = 200
  # 700 / 200 = 3.5; ceiling(3.5) = 4
  # 所以一个 cell 在 x 方向最多可能属于 4 个候选窗口。
  overlap_x <- floor(window_width / step_x) + 1L
  overlap_y <- floor(window_height / step_y) + 1L

  # ceiling(W / step) → 一个 cell 最多可能被几个窗口覆盖？
  # floor((x - 起点) / step) + 1 → 离这个 cell 最近、且不超过它的窗口是第几个？
  # 0:(overlap - 1) → 从这个最近窗口向前找几个？
  # CJ(dx, dy) → 把 x 的所有可能 × y 的所有可能组合起来
  #              ↓
  # 得到这个 cell 可能所属的所有二维 ROI

  # 这一行是在计算：
  # 对每一个 cell，找到“起点最靠近该 cell，但又没有超过该 cell”的那个 x 窗口索引。
  # 假设：
  # x_start = 100, 300, 500, 700, 900, ...
  # cell x = 950
  # 那么最靠近且 <= 950 的起点是 900，
  # 它是第 5 个窗口：ix_latest = 5
  # 数学上：
  # floor((950 - 100) / 200) + 1
  # = floor(4.25) + 1
  # = 5

  cell_df[, ix_latest := floor((x - x_start[1]) / step_x) + 1L]
  cell_df[, iy_latest := floor((y - y_start[1]) / step_y) + 1L]

  # 假设 overlap_x = 4，overlap_y = 4，
  # 那么：dx = 0, 1, 2, 3; dy = 0, 1, 2, 3
  # CJ() 会生成它们的所有组合，一共：4 × 4 = 16
  # dx/dy 不是“多少微米”，而是“从最近那个窗口向前退几个 step”。
  # 例如：
  # ix_latest = 5
  # dx = 0 → ix = 5
  # dx = 1 → ix = 4
  # dx = 2 → ix = 3
  # dx = 3 → ix = 2
  offsets <- data.table::CJ(
    dx = seq.int(0L, overlap_x - 1L),
    dy = seq.int(0L, overlap_y - 1L)
  )

  # Replicate each cell only for windows that it could potentially belong to.
  # 这里开始构造 cell × candidate ROI 的候选关系。
  # 假设：
  # dt 有 10000 个细胞，offsets 有 16 行，
  # 那么每个细胞暂时复制 16 次
  # 复制的目的不是把数据“重复分析 16 次”，
  # 而是为这个 cell 创建最多 16 个“可能所属 ROI”的候选记录。
  membership_df <- cell_df[rep(seq_len(.N), each = nrow(offsets))]

  # 把 offsets 中的 dx / dy 分配给每个 cell 的候选记录。
  # 所以对于每个 cell，它都会得到完整的一套：(dx, dy) 候选组合。
  membership_df[, `:=`(
    dx = rep(offsets$dx, times = nrow(cell_df)),
    dy = rep(offsets$dy, times = nrow(cell_df))
  )
  ]

  # 由 最近窗口 ix_latest / iy_lates
  # 向前退 dx / dy 个窗口，得到这个 cell 可能所属窗口的真正索引。
  # 例如：ix_latest = 5，dx = 2 则：ix = 5 - 2 = 3
  membership_df[, `:=`(
    ix = ix_latest - dx,
    iy = iy_latest - dy
  )
  ]

  # Remove impossible grid indices.
  # 某些位于图像边缘的 cell，往前退几步后可能得到：
  # ix = 0 或 ix = -1 或超过最大窗口编号
  # 这些索引在实际 grid 中不存在，所以先删除
  membership_df <- membership_df[ix >= 1L & ix <= nx & iy >= 1L & iy <= ny]

  # Exact geometric check.
  # Closed intervals retain cells on the outer image boundary. A cell on a
  # shared boundary is intentionally counted in every overlapping window.
  # 前面得到的只是“数学上可能包含这个 cell 的候选 ROI”
  # 这里再做一次精确坐标判断：
  # x_start[ix] <= x <= x_start[ix] + window_width
  # y_start[iy] <= y <= y_start[iy] + window_heigh
  # 只有真正落在窗口内部的 cell × ROI 关系才会留下
  # 也就是说：
  # 前面 = 快速缩小候选范围
  # 这里 = 最终精确确认
  membership_df <- membership_df[
    x >= x_start[ix] & x <= x_start[ix] + window_width &
      y >= y_start[iy] & y <= y_start[iy] + window_height
  ]

  # 现在每一行 membership 已经确定：“某个 cell 属于某个 ix, iy 对应的 ROI”
  # 下一步需要把二维索引 (ix, iy) 转换成一个一维唯一 roi_id。
  # CJ() 先按 ix、再按 iy 排序，因此 roi_id 使用相同的展开顺序。
  membership_df[, roi_id := (ix - 1L) * ny + iy]

  if (nrow(membership_df) == 0L) {
    stop("No cells fall inside the candidate ROIs.")
  }

  # 5. Total number of cells in each ROI
  # 此时 membership 的每一行代表：“一个 cell 位于一个 ROI 中”
  # 所以按照 roi_id 分组后计算 .N，就得到每个 ROI 中总共有多少细胞。
  roi_total_df <- membership_df[, .(n_cells = .N), by = roi_id]

  # 6. Target-cell counts in each ROI
  # 先筛出目标细胞，再按照 ROI + celltype 分组计数。
  # 如果 target_celltype = "Eosinophil"：
  # ROI1  Eosinophil  20
  # ROI2  Eosinophil   5
  # ROI3  Eosinophil  12
  # 如果 target_celltype 有多个，则每个 ROI 会有多行。
  target_count_df <- membership_df[
    celltype %in% target_celltype, .(n_target = .N),
    by = .(roi_id, celltype)
  ]

  # Create every ROI x target-cell combination, including zero counts.
  # 这一点非常重要，target_count 只包含“真的出现过 target cell”的组合。
  # 例如 ROI3 中没有 Tuft，那么 target_count 中可能根本没有：ROI3 + Tuf
  # 但统计时我们需要明确知道：ROI3 中 Tuft = 0
  # 所以这里主动构造：所有 ROI × 所有 target celltype 的完整笛卡尔积。
  target_stat_df <- data.table::CJ(
    roi_id = grid_df$roi_id, celltype = target_celltype
  )

  # 把每个 ROI 的总细胞数加进来。
  target_stat_df <- merge(target_stat_df, roi_total_df, by = "roi_id", all.x = TRUE)

  # 再加入每个 ROI 中每种 target cell 的实际数量。
  # all.x = TRUE 表示：即使某个 ROI 中 target cell 数量为 0，
  # 这个 ROI × celltype 组合仍然保留。
  target_stat_df <- merge(
    target_stat_df, target_count_df, by = c("roi_id", "celltype"), all.x = TRUE
  )

  # 再加入目标细胞在整张图中的背景数量和背景比例。
  # 这样每一行就同时拥有：ROI 局部信息 + 全图背景信息
  target_stat_df <- merge(
    target_stat_df, background_df, by = "celltype", all.x = TRUE
  )

  # 如果某个 ROI 没有任何 cell，merge 后 n_cells 会是 NA
  # 从定义上说它应该是 0，因此改为 0
  target_stat_df[is.na(n_cells), n_cells := 0L]

  # 如果某个 ROI 没有某种 target cell，n_target 会是 NA。
  # 实际上应该理解为：“观察到 0 个”，因此改成 0。
  target_stat_df[is.na(n_target), n_target := 0L]

  # 7. Calculate ROI composition and enrichmen
  # ROI fraction:
  #   fraction = target cells in ROI / all cells in ROI
  #
  # Expected count under random mixing:
  #   expected = n_cells * background_fraction
  #
  # Enrichment:
  #   log2 enrichment = log2((observed + pseudocount) / (expected + pseudocount))
  #
  # Interpretation:
  #    0 = same as image background
  #   +1 = approximately 2-fold enrichmen
  #   +2 = approximately 4-fold enrichmen
  #   -1 = approximately 2-fold depletion

  # target_fraction 是最直观的局部组成比例：
  # target_fraction = ROI 内 target cell 数 / ROI 内所有 cell 数
  # 例如：ROI 中有 500 cells，其中 Eosinophil = 40
  # target_fraction = 40 / 500 = 0.08，即 8%。
  target_stat_df[
    ,
    target_fraction := data.table::fifelse(
      n_cells > 0,
      n_target / n_cells,
      NA_real_
    )
  ]

  # expected_target 是“如果 target cell 在空间中没有局部富集，
  # 只是按照整张图的平均比例随机分布”，那么这个 ROI 理论上应该有多少 target cell。
  # 例如：ROI 总细胞 = 500，全图 Eosinophil 比例 = 3%
  # expected_target = 500 × 0.03 = 15
  target_stat_df[, expected_target := n_cells * background_fraction]

  # observed = n_target; expected = expected_targe
  # 最核心的富集指标：
  # log2_enrichment = log2((observed + 0.5) / (expected + 0.5))
  # 如果 observed = expected,     O/E ≈ 1, log2(1) = 0
  # 如果 observed ≈ expected × 2：O/E ≈ 2, log2(2) = 1
  # 如果 observed ≈ expected × 4：O/E ≈ 4, log2(4) = 2
  # pseudocount 的作用主要是避免 observed = 0 时产生：log2(0) = -Inf
  target_stat_df[,
    log2_enrichment := log2(
      (n_target + pseudocount) / (expected_target + pseudocount)
    )
  ]

  # 8. Hard filtering
  # High enrichment from only one or two rare cells can be unstable and is
  # usually not useful when selecting representative spatial fields.
  # 富集倍数很高，并不一定意味着这个 ROI 很可靠。例如：
  # expected = 0.1; observed = 1
  # O/E 会非常高，但其实 ROI 中只有 1 个 target cell。
  # 所以增加最低 target cell 数量要求。
  # target_pass = TRUE 表示这个 ROI 中该目标细胞数量达到最低要求。
  target_stat_df[, target_pass := n_target >= min_target_cells]

  # 9. Combine multiple target-cell enrichment scores
  # mean:
  #   favors overall enrichment across the requested populations.
  # min:
  #   the weakest target determines the score and is therefore preferable when
  #   all requested cell types should be enriched simultaneously.

  # 如果只有一种 target cell：score 就等于这个 cell type 的 log2_enrichment。
  # 如果有多种 target cell：
  # ROI1:
  # Tuft        = 2.0
  # Eosinophil  = 1.0
  # mean: score = (2 + 1) / 2 = 1.5
  # min: score = min(2, 1) = 1
  # min 更严格：只要有一种 target cell 富集得不好，整个 ROI 的 score 就会被它限制。
  roi_score_df <- target_stat_df[
    ,
    .(
      # 同一个 roi_id 下的 n_cells 是一样的，所以取第一个即可。
      n_cells = data.table::first(n_cells),
      # 如果有多个 target cell，all(target_pass) 要求它们全部达到 min_target_cells
      target_pass = all(target_pass),
      score = if (combine_method == "mean") {
        mean(log2_enrichment)
      } else {
        min(log2_enrichment)
      }
    ),
    by = roi_id
  ]

  # 一个 ROI 最终是否有资格进入排序，需要同时满足：
  # 1. ROI 总细胞数 >= min_cells
  # 2. 所有目标细胞都通过 min_target_cells 要求
  roi_score_df[, eligible := n_cells >= min_cells & target_pass]

  # 10. Convert target-specific statistics to wide forma
  # target_stat 当前是 long format，例如：
  # roi_id   celltype      n_targe
  # 1        Tuft          10
  # 1        Eosinophil    20
  # 2        Tuft           5
  # 2        Eosinophil    12
  # dcast() 后变成 wide format：
  # roi_id  n_target_Tuft  n_target_Eosinophil
  # 1       10             20
  # 2        5             12
  # 这样最终 ranking table 更方便阅读和筛选。

  target_wide_df <- data.table::dcast(
    target_stat_df,
    roi_id ~ celltype,
    value.var = c(
      "n_target",
      "target_fraction",
      "background_fraction",
      "expected_target",
      "log2_enrichment"
    ),
    fill = 0
  )

  # 把 ROI 空间坐标 + 总体评分合并。
  result_df <- merge(grid_df, roi_score_df, by = "roi_id", all.x = TRUE)

  # 再加入每种 target cell 的详细统计信息。
  result_df <- merge(result_df, target_wide_df, by = "roi_id", all.x = TRUE)

  # Only eligible ROIs enter the ranking.
  # 不符合最低细胞数量要求的 ROI，不让它参与排名
  # 与其给它一个很低的 score，直接设为 NA 更明确：“这个 ROI 没有资格参与评分排序”
  result_df[eligible == FALSE, score := NA_real_]

  # 先按照 eligible 排：TRUE 在前，FALSE 在后。#
  # 再按照 score 从高到低排序。前面的 "-" 表示降序。
  data.table::setorder(result_df, -eligible, -score)

  # 只给 eligible == TRUE 的 ROI 编排名
  # 因为前一步已经按 score 降序排好了，所以 rank = 1 就是 score 最高的合格 ROI。
  result_df[eligible == TRUE, rank := seq_len(.N)]

  # 11. Add cell IDs for top-ranked ROIs
  # membership 已经保存了每个 cell 实际属于哪个 ROI。
  # 为避免结果表过大，只为排名前 show_top_roi 的 ROI
  # 保存 cell ID；其余 ROI 的 cell_id 保持 NA。
  result_df[, cell_id := NA_character_]
  top_roi <- result_df[eligible == TRUE & rank <= show_top_roi, roi_id]

  # 将同一个 ROI 中的所有 cell ID 用 "," 连接。
  # 例如：cell_1,cell_2,cell_3,cell_4
  roi_cell_df <- membership_df[
    roi_id %in% top_roi,
    .(cell_id = paste(unique(cell_id), collapse = ",")),
    by = roi_id
  ]

  # 根据 roi_id 将 cell ID 字符串写回结果表。
  result_df[roi_cell_df, on = "roi_id", cell_id := i.cell_id]

  # Put important columns first.
  # 把最常用的信息放在结果表最前面，方便人工查看：
  # ROI 编号、排名、是否合格、enrichment score、总细胞数、ROI 坐标边界
  first_cols <- c(
    "roi_id", "rank", "eligible", "score", "n_cells",
    "x_min", "x_max", "y_min", "y_max", "cell_id"
  )

  # setcolorder() 只改变列顺序，不改变数据内容。
  data.table::setcolorder(
    result_df, c(first_cols, setdiff(names(result_df), first_cols))
  )

  # 函数内部用 data.table 提高计算效率，输出为df
  as.data.frame(result_df)
}

#### strna_roi_plot ####
#' Plot cell polygons inside a spatial region of interes
#'
#' Select every cell polygon whose bounding box overlaps a rectangular region
#' and draw the complete polygons, clipped to the requested plotting limits.
#' Polygon vertices for each cell must already be stored in drawing order.
#'
#' @param boundary Data frame containing polygon vertices and cell annotations.
#' @param xmin,xmax,ymin,ymax Numeric limits of the rectangular ROI.
#' @param cell_col Column containing cell identifiers.
#' @param x_col,y_col Numeric columns containing polygon-vertex coordinates.
#' @param celltype_col Column containing the cell type used for polygon fill.
#' @param palette Optional named character vector mapping cell types to colors.
#' @param reverse_y Whether to reverse the y axis for image-coordinate data.
#' @param title,subtitle Optional plot title and subtitle. When `subtitle` is
#'   `NULL`, ROI dimensions and the number of overlapping cells are shown.
#' @param show_axes Whether to display coordinate axes.
#' @param panel_fill Panel background color.
#' @param polygon_color Polygon-border color; use `NA` for no border.
#' @param polygon_linewidth Width of polygon borders.
#' @param scalebar_position Corner in which to draw the scale bar.
#' @param scalebar_length Scale-bar length in the same units as the coordinates.
#' @param scalebar_color Scale-bar and label color.
#' @param scalebar_textsize Scale-bar label size.
#' @param scalebar_linewidth Scale-bar line width.
#'
#' @return A ggplot object containing the clipped cell polygons.
#' @export
strna_roi_plot <- function(
  boundary, xmin, xmax, ymin, ymax,
  cell_col = "cell_id", x_col = "vertex_x", y_col = "vertex_y",
  celltype_col = "celltype", palette = NULL, reverse_y = TRUE,
  title = NULL, subtitle = NULL,
  show_axes = FALSE, panel_fill = "black",
  polygon_color = NA, polygon_linewidth = 0.08,
  scalebar_position = c("bottomleft", "bottomright", "topleft", "topright"),
  scalebar_length = 100, scalebar_color = "white",
  scalebar_textsize = 3, scalebar_linewidth = 1
) {
  scalebar_position <- match.arg(scalebar_position)

  if (!is.data.frame(boundary)) {
    stop("boundary must be a data frame.")
  }

  # 1. Check inpu
  required_cols <- c(cell_col, x_col, y_col, celltype_col)
  missing_cols <- setdiff(required_cols, colnames(boundary))

  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!is.numeric(boundary[[x_col]]) || !is.numeric(boundary[[y_col]])) {
    stop("x_col and y_col must be numeric.")
  }

  limit_vec <- c(xmin, xmax, ymin, ymax)
  if (length(limit_vec) != 4L || any(!is.finite(limit_vec))) {
    stop("xmin, xmax, ymin and ymax must be finite numeric scalars.")
  }

  if (xmin >= xmax) {
    stop("xmin must be smaller than xmax.")
  }

  if (ymin >= ymax) {
    stop("ymin must be smaller than ymax.")
  }

  if (length(scalebar_length) != 1L || !is.finite(scalebar_length) ||
      scalebar_length <= 0) {
    stop("scalebar_length must be > 0.")
  }

  if (scalebar_length >= (xmax - xmin)) {
    stop("scalebar_length must be smaller than the ROI width.")
  }

  # 2. Prepare boundary table
  # 统一内部列名，避免后面反复处理用户指定的列名
  boundary_df <- data.table::data.table(
    cell_id = boundary[[cell_col]],
    vertex_x = boundary[[x_col]],
    vertex_y = boundary[[y_col]],
    celltype = as.character(boundary[[celltype_col]])
  )

  boundary_df <- boundary_df[
    is.finite(vertex_x) & is.finite(vertex_y) & !is.na(cell_id)
  ]

  if (nrow(boundary_df) == 0L) {
    stop("No valid polygon vertices remain.")
  }

  if (!is.null(palette)) {
    if (!is.character(palette) || is.null(names(palette)) ||
        any(!nzchar(names(palette))) || anyDuplicated(names(palette))) {
      stop("palette must be a character vector with unique, non-empty names.")
    }
    missing_color <- setdiff(
      unique(stats::na.omit(boundary_df$celltype)), names(palette)
    )
    if (length(missing_color) > 0L) {
      warning(
        "Cell types missing from palette: ",
        paste(missing_color, collapse = ", ")
      )
    }
    boundary_df[, celltype := factor(celltype, levels = names(palette))]
  }

  # 3. Find cells overlapping the ROI
  # 找到与指定 ROI 有空间重叠的 cell
  # 这里先计算每个 cell polygon 的外接矩形：
  # poly_xmin -------- poly_xmax
  #     |                  |
  #     |      cell        |
  #     |                  |
  # poly_ymin -------- poly_ymax
  # 只要 cell 的外接矩形和 ROI 有重叠，就保留整个 cell polygon。
  # 最终超出 ROI 的部分会由 coord_fixed() 自动裁掉
  bbox_df <- boundary_df[
    ,
    .(
      poly_xmin = min(vertex_x, na.rm = TRUE),
      poly_xmax = max(vertex_x, na.rm = TRUE),
      poly_ymin = min(vertex_y, na.rm = TRUE),
      poly_ymax = max(vertex_y, na.rm = TRUE)
    ),
    by = cell_id
  ]

  keep_cell <- bbox_df[
    poly_xmax >= xmin & poly_xmin <= xmax &
      poly_ymax >= ymin & poly_ymin <= ymax,
    cell_id
  ]

  roi_df <- boundary_df[cell_id %in% keep_cell]

  if (nrow(roi_df) == 0L) {
    stop("No cells were found in the specified ROI.")
  }

  n_cell <- data.table::uniqueN(roi_df$cell_id)

  # 4. Default subtitle
  if (is.null(subtitle)) {
    subtitle <- paste0(
      round(xmax - xmin), " \u00d7 ", round(ymax - ymin),
      " \u00b5m | ", format(n_cell, big.mark = ","), " cells"
    )
  }

  # 5. Scale bar position
  roi_width <- xmax - xmin
  roi_height <- ymax - ymin

  x_pad <- roi_width * 0.05
  y_pad <- roi_height * 0.06
  text_gap <- roi_height * 0.04

  bottom_y <- if (reverse_y) ymax - y_pad else ymin + y_pad
  bottom_text_y <- if (reverse_y) bottom_y - text_gap else bottom_y + text_gap
  top_y <- if (reverse_y) ymin + y_pad else ymax - y_pad
  top_text_y <- if (reverse_y) top_y + text_gap else top_y - text_gap

  if (scalebar_position == "bottomleft") {
    bar_x1 <- xmin + x_pad
    bar_x2 <- bar_x1 + scalebar_length
    bar_y <- bottom_y
    text_y <- bottom_text_y
  }

  if (scalebar_position == "bottomright") {
    bar_x2 <- xmax - x_pad
    bar_x1 <- bar_x2 - scalebar_length
    bar_y <- bottom_y
    text_y <- bottom_text_y
  }

  if (scalebar_position == "topleft") {
    bar_x1 <- xmin + x_pad
    bar_x2 <- bar_x1 + scalebar_length
    bar_y <- top_y
    text_y <- top_text_y
  }

  if (scalebar_position == "topright") {
    bar_x2 <- xmax - x_pad
    bar_x1 <- bar_x2 - scalebar_length
    bar_y <- top_y
    text_y <- top_text_y
  }

  # 6. Plot polygons
  p <- ggplot2::ggplot(
    roi_df, ggplot2::aes(
      x = vertex_x, y = vertex_y, group = cell_id, fill = celltype
    )
  ) +
    ggplot2::geom_polygon(color = polygon_color, linewidth = polygon_linewidth)

  if (!is.null(palette)) {
    p <- p +
      ggplot2::scale_fill_manual(
        values = palette, drop = FALSE, na.value = "grey80", name = NULL
      )
  }

  if (reverse_y) {
    p <- p +
      ggplot2::scale_y_reverse(expand = c(0, 0))
  } else {
    p <- p +
      ggplot2::scale_y_continuous(expand = c(0, 0))
  }

  p <- p +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::coord_fixed(
      xlim = c(xmin, xmax), ylim = c(ymin, ymax),
      expand = FALSE, clip = "on"
    ) +
    ggplot2::annotate(
      "segment", x = bar_x1, xend = bar_x2, y = bar_y, yend = bar_y,
      linewidth = scalebar_linewidth, color = scalebar_color
    ) +
    ggplot2::annotate(
      "text", x = (bar_x1 + bar_x2) / 2, y = text_y,
      label = paste0(scalebar_length, " \u00b5m"), size = scalebar_textsize,
      color = scalebar_color, fontface = "bold"
    ) +
    ggplot2::labs(
      title = title, subtitle = subtitle,
      x = if (show_axes) "X (\u00b5m)" else NULL,
      y = if (show_axes) "Y (\u00b5m)" else NULL
    )

  # 7. Theme
  if (show_axes) {
    p <- p +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.text = ggplot2::element_text(size = 11, color = "black"),
        axis.title = ggplot2::element_text(size = 11, color = "black"),
        axis.line = ggplot2::element_blank(),
        axis.ticks.length = grid::unit(1, "mm"),
        axis.ticks = ggplot2::element_line(linewidth = 0.5, color = "black"),
        panel.background = ggplot2::element_rect(fill = panel_fill, color = NA),
        panel.border = ggplot2::element_rect(
          color = "black", linewidth = 0.5, fill = NA
        ),
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(hjust = 0.5),
        legend.position = "right",
        legend.key.spacing.y = grid::unit(1, "mm"),
        legend.background = ggplot2::element_rect(fill = "white", color = NA)
      )
  } else {
    p <- p +
      ggplot2::theme_void() +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = panel_fill, color = NA),
        panel.border = ggplot2::element_rect(
          color = "black", linewidth = 0.5, fill = NA
        ),
        plot.title = ggplot2::element_text(
          face = "bold", hjust = 0.5, color = "black"
        ),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "black"),
        legend.position = "right",
        legend.key.spacing.y = grid::unit(1, "mm"),
        legend.background = ggplot2::element_rect(fill = "white", color = NA)
      )
  }

  p
}

#### strna_overview_plot ####
#' Plot a centroid overview of a spatial cell map
#'
#' Collapse each cell polygon to the mean of its vertex coordinates and draw a
#' fast overview colored by cell type. The computed centers are intended for
#' navigation and visualization, not geometric centroid analysis.
#'
#' @param boundary Data frame containing polygon vertices and cell annotations.
#' @param cell_col Column containing cell identifiers.
#' @param x_col,y_col Numeric columns containing polygon-vertex coordinates.
#' @param celltype_col Column containing cell types.
#' @param palette Optional named character vector mapping cell types to colors.
#' @param reverse_y Whether to reverse the y axis for image-coordinate data.
#' @param point_size,point_alpha Point size and opacity.
#' @param raster Whether to use `ggrastr::geom_point_rast()` when the optional
#'   package is installed. Otherwise regular ggplot2 points are used.
#' @param title Optional plot title.
#'
#' @return A ggplot object with one point per cell.
#' @export
strna_overview_plot <- function(
  boundary, cell_col = "cell_id", x_col = "vertex_x", y_col = "vertex_y",
  celltype_col = "celltype", palette = NULL, reverse_y = TRUE,
  point_size = 0.5, point_alpha = 0.9, raster = TRUE, title = NULL
) {
  if (!is.data.frame(boundary)) {
    stop("boundary must be a data frame.")
  }

  # 1. Check inpu
  required_cols <- c(cell_col, x_col, y_col, celltype_col)

  missing_cols <- setdiff(required_cols, colnames(boundary))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # 2. Prepare boundary table
  boundary_df <- data.table::data.table(
    cell_id = boundary[[cell_col]],
    vertex_x = boundary[[x_col]],
    vertex_y = boundary[[y_col]],
    celltype = as.character(boundary[[celltype_col]])
  )

  boundary_df <- boundary_df[
    is.finite(vertex_x) & is.finite(vertex_y) & !is.na(cell_id)
  ]
  if (nrow(boundary_df) == 0L) {
    stop("No valid polygon vertices remain.")
  }

  # 3. Calculate a representative center for each cell
  #
  # 这里的 centroid 只用于整张组织的快速空间导航，
  # 不用于正式的几何/形态分析。
  #
  # 因此直接使用 polygon 顶点坐标的平均值即可，
  # 计算简单、速度快，而且足够反映 cell 的空间位置。
  centroid_df <- boundary_df[
    ,
    .(
      X = mean(vertex_x),
      Y = mean(vertex_y),
      celltype = data.table::first(celltype)
    ),
    by = cell_id
  ]

  if (!is.null(palette)) {
    if (!is.character(palette) || is.null(names(palette)) ||
        any(!nzchar(names(palette))) || anyDuplicated(names(palette))) {
      stop("palette must be a character vector with unique, non-empty names.")
    }
    centroid_df[
      ,
      celltype := factor(
        celltype,
        levels = names(palette)
      )
    ]
  }

  # 4. Plo
  p <- ggplot2::ggplot(
    centroid_df, ggplot2::aes(x = X, y = Y, color = celltype)
  )
  if (raster && requireNamespace("ggrastr", quietly = TRUE)) {
    p <- p +
      ggrastr::geom_point_rast(
        size = point_size, alpha = point_alpha, stroke = 0, raster.dpi = 300
      )
  } else {
    p <- p +
      ggplot2::geom_point(size = point_size, alpha = point_alpha, stroke = 0)
  }

  if (!is.null(palette)) {
    p <- p +
      ggplot2::scale_color_manual(
        values = palette, drop = FALSE, na.value = "grey80", name = NULL
      )
  }

  if (reverse_y) {
    p <- p + ggplot2::scale_y_reverse()
  }

  p +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = title, x = "X (\u00b5m)", y = "Y (\u00b5m)") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 10, color = "black"),
      axis.title = ggplot2::element_text(size = 11, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.4),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.key.spacing.y = grid::unit(1.5, "mm")
    )
}
