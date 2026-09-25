# mengr

`mengr` 是我的个人 R 函数包，用来集中保存经过整理、以后会重复使用的
生物信息学、组学分析、统计、机器学习和科研绘图函数。

这份 README 主要写给以后的自己：隔一段时间没有使用后，可以快速确认包在哪里、
如何安装、输入数据采用什么方向，以及常见任务应该从哪个函数开始。

## 当前信息

| 项目 | 当前本地设置 | 当前工作站设置 |
| --- | --- | --- |
| 版本 | `0.0.3` | `0.0.3` |
| GitHub | `mengjx855/mengr` | `mengjx855/mengr` |
| 源码目录 | `F:/R_proj/mengr` | `/data/mengjx/R_proj/mengr` |
| R 版本 | R 4.4.3 | R 4.4.3 |
| 个人包目录 | `E:/SoftwareData/R/win-library/4.4` | `/home/mengjx/.R/library/4.4` |
| 默认数据库目录 | `F:/database` | `/home/mengjx/database` |

## 使用前先记住

- `profile` 一般是 **feature × sample**：行为 feature，列为 sample。
- `sample_meta` 是样本元数据表，通常至少包含 `sample` 和 `group` 两列。
- `feature_meta` 是与 profile 行名对应的 feature 注释表。
- 函数名使用小写 `snake_case`，例如 `calcu_pca()` 和 `plot_pcoa()`。
- `sample_col`、`group_col` 用来指定元数据中的样本列和分组列。
- 绘图函数通常返回 ggplot、ggraph、patchwork 或 ComplexHeatmap 对象，不自动保存图片。
- 专业依赖按需安装。例如生态统计常用 `vegan`，机器学习可能需要 `caret`、
  `randomForest`、`glmnet` 或 `e1071`，单细胞函数通常需要 `Seurat` 或 `CellChat`。

## 安装或更新

### 从本地源码包安装

```r
lib_dir <- "E:/SoftwareData/R/win-library/4.4"


install.packages(
  "F:/R_proj/mengr/mengr_0.0.3.tar.gz",
  repos = NULL,
  type = "source",
  lib = lib_dir
)

.libPaths(c(lib_dir, .libPaths()))
library(mengr)
packageVersion("mengr")
```

### 从 GitHub 安装

```r
remotes::install_github(
  "mengjx855/mengr",
  lib = "E:/SoftwareData/R/win-library/4.4",
  dependencies = FALSE,
  upgrade = "never"
)
```

包名已经由 `mengR` 改为 `mengr`。旧代码中的 `library(mengR)` 需要相应改成
`library(mengr)`。只有明确需要更新已安装版本时才执行安装；日常修改源码不必反复
覆盖本地包。

## 常用函数入口

| 任务 | 建议先看 |
| --- | --- |
| profile 转换 | `profile_trans_ra()`、`profile_trans_log2()`、`profile_trans_clr()`、`profile_trans_hellinger()` |
| profile 筛选与统计 | `profile_filter()`、`profile_top_n()`、`profile_prevalence()`、`profile_statistics()` |
| 样本或 feature 汇总 | `profile_collapse()`、`profile_aggregate()`、`aggregate_df()` |
| taxonomy 与组成图 | `taxa_split()`、`taxa_trans()`、`plot_compos()`、`plot_taxa_boxplot()` |
| 饼图与圆形打包图 | `plot_pie()`、`plot_circlepack()` |
| alpha/beta diversity | `calcu_alpha()`、`plot_alpha()`、`calcu_distance()`、`calcu_beta()` |
| 稀释与 rank abundance | `calcu_specaccum()`、`plot_specaccum()`、`calcu_rankabund()`、`plot_rankabund()` |
| PCA/PCoA/NMDS | `calcu_pca()`、`plot_pca()`、`plot_pcoa()`、`plot_nmds()` |
| 其他降维与约束排序 | `plot_dbrda()`、`plot_tsne()`、`plot_umap()`、`plot_pls()`、`plot_opls()` |
| 差异分析 | `difference_analysis()`、`calcu_diff_profile()`、`calcu_stamp()`、`plot_stamp()` |
| 相关性与网络 | `calcu_correlation()`、`tidy_correlation()`、`get_nwk_attr()`、`calcu_men()` |
| Mantel 分析 | `calcu_mantel()`、`plot_mantel_lower()`、`plot_mantel_upper()` |
| ROC | `calcu_feature_auc()`、`plot_roc()`、`plot_roc_multiple()` |
| LASSO | `lasso_kfold()`、`lasso_next_validate()` |
| Random forest | `rf_kfold()`、`rf_repeated_kfold()`、`rf_importance()`、`rf_next_validate()` |
| SVM | `svm_base()`、`svm_kfold()`、`svm_next_validate()` |
| 批次效应 | `remove_batch_combat()`、`remove_batch_limma()` |
| RNA-seq count | `rc2tpm()`、`rc2fpkm()`、`rc2rpm()` |
| 单细胞与空间转录组 | `scrna_marker_match()`、`scrna_expression_dotplot()`、`strna_rank_roi()`、`cellchat_run()` |
| 宏基因组与代谢组 | `profile_kegg_trans()`、`tidy_lefse()`、`extract_hmdb_xrefs()` |
| 调色板 | `pald()`、`palc()`、`pald_show()`、`palc_show()` |

完整的领域分类见 [FUNCTION_REFERENCE.md](FUNCTION_REFERENCE.md)。准确的参数、默认值、
返回结构和边界条件以函数帮助为准：

```r
help(package = "mengr")
?profile_filter
?plot_circlepack
?rf_kfold
```

## 常见工作流

### 1. profile 转换与 PCA

```r
library(mengr)

profile_df <- data.frame(
  sample_a = c(10, 30, 60),
  sample_b = c(20, 20, 60),
  sample_c = c(40, 10, 50),
  sample_d = c(30, 30, 40),
  row.names = c("feature_1", "feature_2", "feature_3"),
  check.names = FALSE
)

profile_ra_df <- profile_trans_ra(profile_df, base = 1)
pca_result <- calcu_pca(profile_ra_df)
pca_result$points
```

### 2. 对齐样本元数据并绘制 PCoA

```r
sample_meta <- data.frame(
  sample = c("sample_a", "sample_b", "sample_c", "sample_d"),
  group = c("control", "control", "case", "case")
)

pcoa_plot <- plot_pcoa(
  profile = profile_ra_df,
  sample_meta = sample_meta,
  sample_col = "sample",
  group_col = "group"
)
```

### 3. 汇总类别并绘制组成图

```r
composition_df <- data.frame(
  name = c("A", "B", "C", "D"),
  n = c(40, 30, 20, 10)
)

plot_pie(composition_df)

# 需要 ggraph、tidygraph 和 withr
plot_circlepack(composition_df, add_percent = TRUE)
```

### 4. 取得和预览调色板

```r
pald("Set2", n = 5)
pald_show("Set2")
palc("BlWhRe", n = 20)
```

## 数据库路径

部分注释和数据库函数从公共数据库根目录读取文件：

```r
mengr_config(database = "F:/database")
mengr_config()
```

也可以通过 R option `mengr.database` 或环境变量 `MENGR_DATABASE` 设置。当前默认值为
`F:/database`。

## 隔久了以后从哪里看

1. 先看 [FUNCTION_REFERENCE.md](FUNCTION_REFERENCE.md)，按分析领域找函数。
2. 再运行 `?function_name`，确认参数、输入方向和返回对象。
3. 旧代码出现函数名或参数错误时，看 [MIGRATION_NOTES.md](MIGRATION_NOTES.md)。
4. 想知道当前版本修改了什么，看 [NEWS.md](NEWS.md)。
5. 新写或整理函数时，看 [STYLE_GUIDE.md](STYLE_GUIDE.md)。

## 我自己维护时

帮助文档写在 `R/` 源码的 Roxygen 注释中，不直接修改 `man/*.Rd`。重新生成帮助：

```r
roxygen2::roxygenise(
  ".",
  roclets = c("collate", "namespace", "rd")
)
```

运行源码级检查：

```powershell
$env:R_LIBS_USER = "E:/SoftwareData/R/win-library/4.4"
& "D:\R\R-4.4.3\bin\Rscript.exe" "tools/smoke-test.R"
```

构建并检查正式源码包：

```powershell
$env:R_LIBS_USER = "E:/SoftwareData/R/win-library/4.4"
& "D:\R\R-4.4.3\bin\R.exe" CMD build --no-manual --no-build-vignettes .
& "D:\R\R-4.4.3\bin\R.exe" CMD check --no-manual mengr_0.0.3.tar.gz
```

代码风格和修改记录按照个人 skill `$mengr-package-style` 维护。

## License

MIT License，见 [LICENSE](LICENSE)。
