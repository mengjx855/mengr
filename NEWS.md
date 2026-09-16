# mengR 0.0.2

## 函数体与行为修改

- `format_stat_label()`：增加 P 值范围、`digits` 和 `scientific` 的输入检查。
- `add_strip()`：条带背景颜色改为常量映射，不再添加
  `scale_fill_identity()`，避免覆盖原图已有的 fill scale；同时增加参数检查。
- `strna_rank_roi()`：修正二维 ROI 网格到 `roi_id` 的索引顺序；边界细胞按闭区间
  纳入；修正每个细胞可能覆盖的滑动窗口数量，并补充输入与空结果检查。
- `strna_roi_plot()`：修正 `palette = NULL` 时访问调色板名称的问题，以及非反转
  y 轴时比例尺的垂直位置；补充输入检查。
- `strna_overview_plot()`：补充输入、空数据和调色板检查。
- `scrna_reduction_plot()`：将错误的 `ggplot2::unit()` 改为 `grid::unit()`，并补充
  Seurat、降维结果、抽样比例和绘图参数检查。
- `scrna_seurat_subsample()`：补充输入与抽样检查，新增 `verbose` 控制消息输出。
- `cellchat_run()`：Seurat v5 layer 操作改用实际导出这些函数的
  `SeuratObject::Layers()` 和 `SeuratObject::JoinLayers()`。
- `palc()`：把新增的 32 组连续渐变锚点整合进内部 palette 列表，按需生成颜色；
  移除脚本加载时的 `library(grDevices)` 和 32 个顶层生成器函数。
- `extract_hmdb_xrefs()`：将 `progresults` 参数改为 `progress`，并修正不存在的
  `txtProgresultsBar()` / `setTxtProgresultsBar()` 调用。
- `plot_compos_multiple()`：移除 taxonomy rank 列表末尾的缺失元素，并先捕获
  `...` 再安全转发给 `plot_compos()`。

## API 与命名修改

- 所有公开和内部函数名统一为小写 snake_case；原有大写缩写函数名不再导出。
- 主要改名包括 `mengr_config()`、`calcu_pca()`、`plot_pca()`、`calcu_ncm()`、
  `plot_ncm()`、`plot_dbrda()`、`plot_pls()`、`plot_opls()`、`calcu_men()`、
  `calcu_dbe()`、`extract_hmdb_xrefs()`、`parse_mchromatograms()`、
  `mbt_ekegg()`、`load_kegg_info()`、`profile_kegg_trans()`、`tidy_lefse()`、
  `tidy_cazyme_profile()`、`plot_go_bar()` 和 `plot_go_circular_bar()`。
- 向量转换函数统一为 `transform_log2()`、`transform_log10()` 和
  `transform_clr()`。
- data.frame、matrix、vector、list 和通用对象的局部变量分别采用
  `*_df`、`*_mat`、`*_vec`、`*_list` 和 `*_obj` 风格；保留既有公开返回组件名和
  plot attributes，避免仅因内部变量整理造成额外 API 破坏。

## 文档与包结构

- 版本升级到 0.0.2。
- README 改为个人使用手册，记录本地环境、常用入口、安装更新方式和长期未使用后的
  查找顺序。
- 49 个 R 脚本统一文件头、按时间排列的修改记录和
  `#### function_name ####` 函数章节标题。
- 逐函数补写参数与真实返回结构，并重新生成 Rd 和 NAMESPACE。
- 补齐 `data.table`、`ggrastr`、`DESeq2`、`edgeR`、`SeuratObject` 和
  `SummarizedExperiment` 等实际使用的建议依赖；`rlang` 改列为 Imports。
