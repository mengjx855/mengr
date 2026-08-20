# mengR 0.0.1 函数总览

本文汇总 `mengR` 当前源码中的公开函数及主要功能，便于长期未使用后快速定位。
函数的完整参数、输入格式和边界条件以各函数的 roxygen2 帮助文档为准；函数改名、
删除和行为变化见 [MIGRATION_NOTES.md](MIGRATION_NOTES.md)，统一代码习惯见
[STYLE_GUIDE.md](STYLE_GUIDE.md)。

## 1. 统一输入与返回约定

- `profile` 通常表示 feature × sample 数值表，行为 feature，列为 sample。
- `group` 通常包含 `sample_col` 和 `group_col`，并按样本名与 profile 对齐。
- `plot_df`、`test_df` 和 `result_df` 分别表示作图、检验和最终结果数据。
- 绘图函数返回 ggplot、patchwork、ComplexHeatmap 或相应绘图对象，不自动保存。
- 复杂分析通常返回命名 list，或把模型和作图数据保存在返回对象的属性中。
- 外部数据库默认根目录可通过 `mengR_config()` 配置。

## 2. 配置与通用工具

| 函数 | 功能描述 |
| --- | --- |
| `mengR_config()` | 设置或读取 mengR 配置；目前主要管理公共数据库根目录。 |
| `floor_n()` | 按指定的十进制位数向下取整，例如按百位或千位取整。 |
| `ceiling_n()` | 按指定的十进制位数向上取整。 |
| `set_calcu()` | 对多个向量或一个向量列表计算交集、并集或集合差。 |
| `get_text_color()` | 根据背景色亮度选择深色或浅色文字，提高文字对比度。 |
| `write_xlsx_with_comment()` | 将顶部说明和数据写入 Excel；支持追加、替换 worksheet。 |
| `read_xlsx_multiple()` | 一次读取 Excel 文件中的多个 worksheet，返回命名列表。 |
| `pairwise_cluster()` | 将 feature1-feature2-value 长表转换为距离矩阵并进行层次聚类。 |
| `theme_bw_clean()` | 返回简洁的黑白 ggplot2 theme，作为其他绘图的统一基础主题。 |

## 3. 调色板工具

| 函数 | 功能描述 |
| --- | --- |
| `pald()` | 按名称取得离散配色，可指定数量、反转顺序或输出可复制文本。 |
| `pald_show()` | 预览一个离散调色板，或列出可用离散调色板名称。 |
| `palc()` | 根据预设名称生成指定数量的连续渐变颜色。 |
| `palc_show()` | 预览指定的连续调色板。 |

## 4. Profile 转换与基础处理

### 4.1 数值转换

| 函数 | 功能描述 |
| --- | --- |
| `LOG2()` | 对一个数值向量进行 log2 转换，并处理零值和缺失值。 |
| `LOG10()` | 对一个数值向量进行 log10 转换，并处理零值和缺失值。 |
| `CLR()` | 对组成型数值向量执行 centered log-ratio 转换。 |
| `profile_trans_log2()` | 对 feature × sample profile 按 feature 执行 log2 转换。 |
| `profile_trans_log10()` | 对 feature × sample profile 按 feature 执行 log10 转换。 |
| `profile_trans_clr()` | 对每个样本执行 CLR 转换，保持 feature × sample 方向。 |
| `profile_trans_sqrt()` | 对 profile 执行平方根转换。 |
| `profile_trans_ra()` | 将每个样本转换为相对丰度或百分比，并处理空样本。 |
| `profile_trans_hellinger()` | 对 profile 执行 Hellinger 转换，适合部分生态距离分析。 |

### 4.2 筛选、汇总与统计

| 函数 | 功能描述 |
| --- | --- |
| `profile_collapse()` | 按样本分组对 profile 求 mean、median、sum 或自定义统计量。 |
| `profile_filter()` | 按总体或组内 prevalence、出现样本数和最低丰度筛选 feature。 |
| `profile_top_n()` | 保留总体丰度最高的 n 个 feature，可把其余 feature 合并为 Other。 |
| `profile_top_frac()` | 按比例保留丰度最高的 feature，可把其余部分合并为 Other。 |
| `profile_replace()` | 按阈值替换低丰度值，并可在替换前转换为相对丰度。 |
| `profile_adjacency()` | 将丰度表转换为 feature × sample 的存在/缺失矩阵。 |
| `profile_prevalence()` | 计算 feature 在总体或各组中的出现样本数或 prevalence。 |
| `profile_statistics()` | 计算 feature 的 mean、SD、median 和 prevalence，可按组汇总。 |
| `profile_aggregate()` | 根据 feature metadata 聚合 profile，并处理未知分类。 |
| `profile_remove_zero_var()` | 删除跨样本零方差的 feature，避免降维或模型拟合失败。 |

## 5. Read-count 标准化

| 函数 | 功能描述 |
| --- | --- |
| `rc2tpm()` | 根据 feature 长度把 read counts 转换为 TPM。 |
| `rc2fpkm()` | 根据 feature 长度和样本文库大小计算 FPKM/RPKM。 |
| `rc2rpm()` | 将 counts 转换为 RPM；文库大小可来自列和、命名向量、表格或 CSV/TSV 文件。 |

## 6. Taxonomy 与组成图

| 函数 | 功能描述 |
| --- | --- |
| `profile_mpa()` | 处理 MetaPhlAn 风格 profile，拆分层级并输出指定分类水平。 |
| `taxa_split()` | 将完整 taxonomy 字符串拆分为 domain 至 species 等分类列。 |
| `taxa_trans()` | 按 taxonomy 层级聚合 profile，选择 top taxa，并可按组或相对丰度转换。 |
| `plot_compos()` | 绘制单个分类层级的样本或分组组成堆叠柱状图。 |
| `plot_compos_multiple()` | 对多个 taxonomy 层级批量绘制组成图并组合输出。 |
| `plot_compos_manual()` | 对已经整理好的组成数据绘制可精细控制的堆叠柱状图。 |
| `plot_taxa_boxplot()` | 对每个 taxa 绘制分组箱线图，并添加显著性比较。 |

## 7. 批次效应处理

| 函数 | 功能描述 |
| --- | --- |
| `remove_batch_combat()` | 使用 `sva::ComBat()` 校正一个或两个批次变量，可保留协变量效应。 |
| `remove_batch_limma()` | 使用 `limma::removeBatchEffect()` 校正批次并保留设计矩阵中的协变量。 |

## 8. Alpha/Beta diversity 与距离

| 函数 | 功能描述 |
| --- | --- |
| `calcu_alpha()` | 计算 richness、Shannon、Simpson、Chao1、ACE 等 alpha diversity 指标。 |
| `plot_alpha()` | 绘制 alpha diversity 分组图，并支持排序、显著性及参考线。 |
| `calcu_distance()` | 计算 Bray、Jaccard、Euclidean、UniFrac 等样本距离，可先转换 profile。 |
| `calcu_beta()` | 基于距离矩阵执行整体 beta-diversity 组间检验。 |
| `calcu_adjusted_r2()` | 从 adonis 类结果中计算或提取 adjusted R-squared。 |
| `calcu_adonis_r2()` | 根据距离对象和分组标签计算 PERMANOVA R-squared。 |
| `calcu_pairwise_adonis()` | 对各组组合执行 pairwise PERMANOVA，并校正 P 值。 |
| `plot_pairwise_adonis()` | 将 pairwise PERMANOVA 结果绘制为矩阵或热图式结果图。 |
| `calcu_betadisper()` | 检验各组到中心或中位数中心的 multivariate dispersion。 |
| `plot_betadisper()` | 绘制组内离散度及组间比较结果。 |

## 9. 稀释曲线与丰度排序

| 函数 | 功能描述 |
| --- | --- |
| `calcu_specaccum()` | 计算总体 species/feature accumulation curve。 |
| `plot_specaccum()` | 绘制总体累积曲线，可添加误差线或置信 ribbon。 |
| `calcu_specaccum_by_group()` | 分组计算 species/feature accumulation curve。 |
| `plot_specaccum_by_group()` | 绘制多组累积曲线并使用统一分组配色。 |
| `calcu_specaccum_by_depth()` | 逐步增加测序深度，计算随机稀释后的 feature 数量。 |
| `plot_specaccum_by_depth()` | 绘制测序深度与观测 feature 数量的关系。 |
| `calcu_rankabund()` | 计算每个样本的 rank-abundance 数据及 log abundance。 |
| `plot_rankabund()` | 绘制 rank-abundance 曲线。 |

## 10. 差异检验与富集统计

| 函数 | 功能描述 |
| --- | --- |
| `calcu_empirical_p()` | 根据置换或背景分布计算双侧、左侧或右侧 empirical P 值。 |
| `calcu_diff()` | 对公式指定的数值与分组执行 Wilcoxon、ANOVA 或 t-test 两两比较。 |
| `calcu_diff_profile()` | 对 profile 中的每个 feature 批量执行分组差异检验。 |
| `difference_analysis()` | 汇总差异检验、fold change、均值和 prevalence 的完整 feature 分析流程。 |
| `calcu_fisher()` | 对每个 feature 的 2×2 计数表执行 Fisher exact test。 |
| `calcu_feature_auc()` | 逐 feature 计算二分类 ROC AUC、置信区间、方向和区分强度。 |
| `calcu_stamp()` | 计算两组 STAMP 风格均值差、置信区间和显著性结果。 |
| `plot_stamp()` | 将 `calcu_stamp()` 结果绘制为 STAMP 风格多面板图。 |
| `calcu_stamp_multiple()` | 对多组数据计算 STAMP 风格的整体差异统计。 |

## 11. 相关性、Mantel 与网络

| 函数 | 功能描述 |
| --- | --- |
| `calcu_correlation()` | 计算 feature 间 Pearson、Spearman 或 Kendall 相关及校正 P 值。 |
| `tidy_correlation()` | 将相关矩阵整理为 edge long table，并按相关性和显著性筛选。 |
| `get_nwk_attr()` | 从相关 edge table 构建 igraph 网络，整理 edge/node 属性并可选导出。 |
| `get_nwk_stat()` | 计算网络节点的 degree、strength、centrality 等拓扑统计量。 |
| `calcu_mantel()` | 在多个 profile 或环境变量之间批量执行 Mantel test。 |
| `make_curve_path()` | 为 Mantel 图中的节点连线生成平滑曲线路径坐标。 |
| `plot_mantel_lower()` | 在环境相关矩阵下三角区域叠加 Mantel 关联曲线。 |
| `plot_mantel_upper()` | 在环境相关矩阵上三角区域叠加 Mantel 关联曲线。 |
| `calcu_MEN()` | 使用 ggClusterNet 流程构建 microbial ecological network。 |

## 12. Meta-analysis 与 enterotype

| 函数 | 功能描述 |
| --- | --- |
| `calcu_metafor()` | 按项目计算 feature 效应量，并使用 metafor 合并多队列结果。 |
| `calcu_jsd_dist()` | 对组成型 profile 计算 Jensen-Shannon divergence 距离。 |
| `pam_clustering()` | 对距离对象执行 partitioning around medoids 聚类。 |

## 13. 降维计算与通用绘图

| 函数 | 功能描述 |
| --- | --- |
| `plot_dim()` | 统一绘制二维降维结果，支持 point/line、置信区域、标签和主题。 |
| `calcu_PCA()` | 对 feature × sample profile 执行 PCA 并返回坐标和解释方差。 |
| `plot_PCA()` | 对齐 profile 与分组后执行并绘制 PCA。 |
| `calcu_PCoA()` | 从 profile 或距离对象计算 PCoA 坐标和特征值。 |
| `plot_PCoA()` | 计算并绘制 PCoA，可附加 PERMANOVA 统计结果。 |
| `plot_PCoA_box()` | 同时展示 PCoA 散点图及主要坐标轴的边际箱线图。 |
| `plot_NMDS()` | 从 profile 或距离对象执行 NMDS，可附加 ANOSIM 结果。 |
| `calcu_pairwise_anosim()` | 对各分组组合执行 pairwise ANOSIM。 |
| `plot_pairwise_anosim()` | 将 pairwise ANOSIM 结果绘制为矩阵式结果图。 |
| `plot_dbRDA()` | 从 profile 或距离对象执行 constrained dbRDA，并绘制样本和变量箭头。 |
| `plot_tsne()` | 对齐样本分组后执行 t-SNE，并返回统一风格降维图。 |
| `plot_umap()` | 对齐样本分组后执行 UMAP，并返回统一风格降维图。 |
| `plot_PLS()` | 使用 ropls 执行 PLS-DA，并绘制前两个 predictive components。 |
| `plot_OPLS()` | 对二分类数据执行 OPLS-DA 并绘制 predictive/orthogonal components。 |
| `plot_procrustes()` | 比较两个 profile 或距离空间，执行 Procrustes 和 permutation test。 |

## 14. 常用统计图和生物信息图

| 函数 | 功能描述 |
| --- | --- |
| `plot_heatmap()` | 绘制 ComplexHeatmap 热图，并按 `name` 对齐行列注释。 |
| `plot_volcano()` | 根据 effect、P 值和富集方向绘制 volcano plot。 |
| `plot_pie()` | 从指定名称列和值列绘制饼图或环形图。 |
| `plot_NCM()` | 绘制 neutral community model 的拟合曲线和置信边界。 |
| `calcu_NCM()` | 拟合 neutral community model，返回摘要、拟合数据和模型。 |
| `plot_msa()` | 读取 multiple sequence alignment 并按氨基酸类别分块绘图。 |
| `plot_ggradar_pair()` | 将成对比较数据整理为 ggradar 所需格式并绘制雷达图。 |
| `plot_pair_radar()` | 使用 ggplot2 绘制成对比较的极坐标雷达图。 |

## 15. ROC 工具

| 函数 | 功能描述 |
| --- | --- |
| `theme_roc()` | 提供 default、grid、classic、minimal、nature 和 lancet ROC 主题。 |
| `roc_auc_label()` | 从 pROC 对象生成包含 AUC 和置信区间的标签。 |
| `roc_se_data()` | 提取 ROC sensitivity 置信区间并整理为 ribbon 数据。 |
| `plot_roc()` | 绘制单个 pROC ROC 曲线，可添加 sensitivity 置信 ribbon。 |
| `plot_roc_multiple()` | 在同一图中绘制多个 ROC 对象及各自 AUC 标签。 |

## 16. 机器学习

### 16.1 LASSO

| 函数 | 功能描述 |
| --- | --- |
| `lasso_kfold()` | 使用 glmnet 和分层 folds 执行 LASSO k-fold cross-validation。 |
| `lasso_next_validate()` | 在独立数据集上应用已训练的 LASSO 模型并评估预测。 |

### 16.2 Random forest

| 函数 | 功能描述 |
| --- | --- |
| `rf_kfold()` | 执行分层 random-forest k-fold cross-validation。 |
| `rf_repeated_kfold()` | 多次重复 random-forest k-fold，并标记每次重复和 fold。 |
| `rf_cross_dataset_validate()` | 按 dataset 留一验证，评估跨队列 random-forest 泛化能力。 |
| `rf_importance()` | 提取并整理 random-forest feature importance。 |
| `rf_leave_one_out()` | 执行 leave-one-out random-forest 预测。 |
| `rf_next_validate()` | 用训练数据拟合 random forest，并在独立数据中验证。 |

### 16.3 Support vector machine

| 函数 | 功能描述 |
| --- | --- |
| `svm_base()` | 拟合带参数搜索的基础 SVM，并返回样本预测和性能。 |
| `svm_kfold()` | 执行分层 SVM k-fold cross-validation。 |
| `svm_next_validate()` | 在训练集调参拟合 SVM，并在独立数据中验证。 |

## 17. 单细胞与 CellChat

| 函数 | 功能描述 |
| --- | --- |
| `scRNA_marker_db_build()` | 从 CellMarker 风格表格构建按物种和组织划分的本地 RDS marker 数据库。 |
| `scRNA_marker_match()` | 将输入 marker 与本地细胞类型基因集匹配并计算双向覆盖率。 |
| `scRNA_feature_plots()` | 对命名 marker 列表批量绘制 Seurat FeaturePlot 并按行组合。 |
| `scRNA_violin_plots()` | 从 Seurat 提取 marker 表达并绘制纵向分面 violin plots。 |
| `scRNA_summarise_genes()` | 按 metadata 分组汇总基因阳性率、count 和标准化表达量。 |
| `scRNA_expression_dotplot()` | 将基因表达汇总表绘制为大小代表阳性率、颜色代表表达量的 dotplot。 |
| `cellchat_run()` | 从 Seurat RNA assay 完成 CellChat 通讯概率、通路网络和中心性分析。 |
| `cellchat_radar()` | 提取一个细胞群的 outgoing/incoming CellChat count 或 weight 并绘制雷达图。 |

## 18. 宏基因组工具

| 函数 | 功能描述 |
| --- | --- |
| `plot_maaslin3_multi()` | 绘制 MaAsLin3 多分类层级的正负效应及 taxonomy 色带。 |
| `plot_maaslin3_abundance()` | 绘制 MaAsLin3 abundance model 的 effect 结果。 |
| `plot_maaslin3_both()` | 联合展示 MaAsLin3 abundance 和 prevalence 模型结果。 |
| `load_KEGG_info()` | 读取并整理本地 KEGG 注释表。 |
| `profile_KEGG_trans()` | 将 KO profile 映射并聚合到 KEGG A/B/C 层级。 |
| `tidy_LEfSe()` | 统一 LEfSe 结果列和 taxonomy 字段，便于后续筛选绘图。 |
| `tidy_CAZyme_profile()` | 整理 CAZyme profile 名称并按目标分类层级聚合。 |

## 19. 代谢组、化学式与数据库

| 函数 | 功能描述 |
| --- | --- |
| `calcu_DBE()` | 从分子式计算 double-bond equivalents，并正确处理卤素元素。 |
| `parse_MChromatograms()` | 汇总每条 chromatogram 的保留时间、最高峰位置和强度分位数。 |
| `MBT_eKEGG()` | 将代谢物列表与本地 eKEGG 注释关联并整理 pathway 结果。 |
| `extract_HMDB_xrefs()` | 流式解析 HMDB XML，提取代谢物名称、标识符和外部数据库交叉引用。 |

## 20. Reactome 与 RNA-seq 绘图

| 函数 | 功能描述 |
| --- | --- |
| `reactome_longest_path()` | 根据本地 Reactome 父子关系表提取目标通路的最长层级路径。 |
| `plot_gsea_barcode()` | 绘制 GSEA barcode/running-score 风格结果图。 |
| `plot_GO_bar()` | 绘制 GO enrichment 柱状图，并支持分组、排序和标签。 |
| `plot_GO_circular_bar()` | 将 GO enrichment 结果绘制为 circular barplot。 |

## 21. 当前包内隐藏 helper

以下对象以点开头，只服务于包内函数，不应作为稳定的用户接口调用。

| Helper | 内部作用 |
| --- | --- |
| `.as_df()` | 在保留原始列名的前提下统一转换为 data frame。 |
| `.as_profile_df()` | 检查并整理 feature × sample profile，可强制转换为数值。 |
| `.check_columns()` | 统一检查输入对象是否缺少必要列。 |
| `.match_distance_method()` | 集中维护可用 distance method。 |
| `.match_transform_method()` | 集中维护 profile transformation method。 |
| `.resolve_group_colors()` | 根据 group level 生成或对齐命名颜色。 |
| `.align_profile_group()` | 按样本名对齐 profile 与 group table。 |
| `.profile_long_df()` | 将 feature × sample profile 转换为带分组的长表。 |
| `.mengR_db_file()` | 根据配置的数据库根目录生成数据库文件路径。 |
| `.align_gene_length()` | 按 feature 名对齐 gene-length table。 |
| `.prepare_batch_data()` | 为 ComBat/limma 批次校正整理表达矩阵和 metadata。 |
| `.prepare_rf_data()` | 为 random forest 对齐 profile、分组和 outcome。 |
| `.rf_probability_df()` | 从 random-forest 模型提取类别概率。 |
| `.prepare_svm_data()` | 为 SVM 对齐 profile、分组和 outcome。 |
| `.fit_tuned_svm()` | 执行 SVM 参数搜索并返回最佳模型。 |
| `.svm_prediction_df()` | 将 SVM 类别和概率预测整理为标准数据框。 |
| `.binary_metrics()` | 从真实类别和预测类别计算二分类指标。 |
| `.lasso_outcome()` | 根据 glmnet family 准备 LASSO outcome。 |
| `.heatmap_annotation_palette()` | 为 heatmap 分类注释生成配色。 |
| `.prepare_heatmap_annotation()` | 检查并按 name 对齐 heatmap annotation。 |
| `.check_group()` | 检查 STAMP 分析所需的分组信息。 |
| `.add_p_q_label()` | 为 STAMP 结果生成 P/Q 值标签。 |
| `.add_plab()` | 将 P 值转换为星号或格式化显著性标签。 |

## 22. 可选依赖提示

- `cellchat_run()` 需要 `CellChat`；当前本地依赖审计显示尚未安装。
- Seurat 绘图和汇总函数需要 `Seurat` 与 `Matrix`。
- 生态统计主要使用 `vegan`、`phyloseq`、`picante` 和 `BiodiversityR`。
- 机器学习分别使用 `glmnet`、`randomForest`、`e1071` 和 `caret`。
- 大部分绘图依赖 `ggplot2`，组合图可能使用 `patchwork` 或 `cowplot`。
- 本文只描述当前源码接口；包尚未 build 或 install。
