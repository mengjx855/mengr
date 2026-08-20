# mengR 0.0.1 迁移与修改记录

本文记录从 `R_func` 整理到 `mengR` 最初版本时做出的取舍。原始 `R_func`
目录没有修改；需要比较旧行为时，仍以原文件为准。

## 1. 未纳入的外部代码

以下文件的来源标记、作者信息或代码结构表明它们不是 Jin-Xin Meng 独立编写，
因此没有复制到包中：

- `calcu_GSEA.R`：原文件注明 “by Prof Li”。
- `calcu_REI.R`：来自 MicrobeDS。
- `trans_diff.R`：外部 R6 实现。
- `public/info.centrality.R`：公共/外部脚本。
- `others/zy_corr.R`：他人函数。
- `plot_flower.R`：主体改写自外部博客且直接操作 base graphics，不作为包 API。

## 2. 删除的重复或不成形函数

| 原函数或文件 | 处理方式 | 推荐替代 |
| --- | --- | --- |
| `get.nwk.attr()` | 删除旧点号版本 | `get_nwk_attr()` |
| `calcu_specaccum_manual()` | 删除重复实现 | `calcu_specaccum(method = "random")` |
| `get_adj()` | 删除重复整理函数 | `tidy_correlation(only_signif = TRUE)` |
| `corr_process()` / `corr_merge()` | 删除旧流程 | `calcu_correlation()` + `tidy_correlation()` |
| `rc2rpkm()` | 删除同公式重复函数 | `rc2fpkm()`；read count 输入时结果即 RPKM |
| `calu_adjusted_r2()` | 删除拼写错误副本 | `calcu_adjusted_r2()` |
| `add_plab()` | 删除公开重复函数 | 包内辅助函数 `.add_plab()` |
| `plot_venn.R` | 空文件，不纳入 | 需要时直接使用成熟 Venn/UpSet 包 |
| `plot_boxplots()` | 与 taxa boxplot 重复 | `plot_taxa_boxplot()` |
| 旧版 `calcu_metafor.1` | 删除旧版本 | `calcu_metafor()` |
| `mjx_Sparcc_stat.R` 旧包装 | 删除点号风格包装 | correlation/network 系列函数 |

## 3. 主要函数改名

### 拼写、大小写与统一命名

| 原名 | 新名 |
| --- | --- |
| `lasso_next_vaildate()` | `lasso_next_validate()` |
| `rf_cross_dataset_vaildate()` | `rf_cross_dataset_validate()` |
| `rf_next_vaildate()` | `rf_next_validate()` |
| `svm_next_vaildate()` | `svm_next_validate()` |
| `PAM_clustering()` | `pam_clustering()` |
| `rf_loom()` | `rf_leave_one_out()` |
| `lasso_Kfold()` | `lasso_kfold()` |
| `rf_Kfold()` | `rf_kfold()` |
| `rf_Kfold_rep()` | `rf_repeated_kfold()` |
| `svm_Kfold()` | `svm_kfold()` |
| `plot_Procrustes()` | `plot_procrustes()` |
| `removeBatch_by_combat()` | `remove_batch_combat()` |
| `removeBatch_by_limma()` | `remove_batch_limma()` |
| `profile_transLOG2()` | `profile_trans_log2()` |
| `profile_transLOG10()` | `profile_trans_log10()` |
| `profile_transCLR()` | `profile_trans_clr()` |
| `profile_transSqrt()` | `profile_trans_sqrt()` |
| `profile_transRA()` | `profile_trans_ra()` |
| `profile_transHell()` | `profile_trans_hellinger()` |
| `profile_rmZeroVar()` | `profile_remove_zero_var()` |
| `scRNA_featurePlots()` | `scRNA_feature_plots()` |
| `scRNA_violinPlots()` | `scRNA_violin_plots()` |

### 新增的四个 scRNA/CellChat 函数

| 原始新增名 | mengR 名称 | 命名理由 |
| --- | --- | --- |
| `run_cellchat()` | `cellchat_run()` | 以分析领域作前缀，便于自动补全检索 |
| `cc_interactionRadarPlot()` | `cellchat_radar()` | 保留 CellChat 与 radar 两个核心信息 |
| `summarise_key_genes()` | `scRNA_summarise_genes()` | 明确是单细胞分组表达汇总 |
| `plot_sc_keygene_dot()` | `scRNA_expression_dotplot()` | 明确图中表达量与阳性比例的语义 |

## 4. 参数和行为变化

- `rc2rpm()` 的 `library_size` 可为 `NULL`、带名称数值向量、两列数据框，或
  CSV/TSV 文件路径。表格通过 `sample_col` 和 `library_size_col` 按名称对齐；
  `NULL` 时使用 `colSums(profile)`。
- `dis_method` 统一为 `dist_method`；距离和转换方法通过 `match.arg()` 检查。
- 分组接口尽量统一为 `sample_col`、`group_col`、`group_level` 和
  `group_color`。
- `transRA` 改为 `trans_ra`；`convert2group` 改为 `collapse_group`。
- t-SNE 修正了旧代码忽略输入 `perplexity`、始终使用 1 的问题。
- `calcu_specaccum(method = "collector")` 没有随机标准差时，`sd` 返回 `NA`
  而不是生成长度不一致的数据框。
- heatmap 注释始终根据 `name` 对齐，不依赖输入行顺序。
- `get_nwk_attr()` 默认只返回对象；只有明确提供导出位置时才写文件。
- 绘图函数不自动保存文件，也不再打印固定的 `ggsave()` 命令。

## 5. 代码组织与依赖

- 文件名统一为小写 kebab-case，并使用 `plot-`、`stats-`、`profile-`、
  `model-` 和 `utils-` 前缀。
- 外部调用使用 `package::function()`；管道统一为 base R `|>`。
- 点前缀仅用于包级内部 helper 和常量，局部变量使用带语义的普通名称。
- 公共数据库根目录可通过 `mengR_config(database = ...)`、R option
  `mengR.database` 或环境变量 `MENGR_DATABASE` 配置。
- CellChat 是新增 `cellchat_run()` 的可选依赖。本地审计时其余依赖已安装，
  `CellChat` 尚未安装。

## 6. 当前验证边界

- 已逐文件执行 `parse()` 和隔离环境 `sys.source()`。
- 已对 read-count 转换、profile 转换、距离、JSD、rarefaction、Fisher、
  heatmap 注释对齐和 CellMarker 数据库流程运行小数据 smoke test。
- 按要求尚未运行 package build、install 或完整 `R CMD check`。
- CellChat 主流程需要安装 CellChat 并使用真实 Seurat 对象进一步验证。
