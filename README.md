# mengR

这是我自己的 R 函数包，用来集中保存已经整理过、以后还会重复使用的
生物信息学、组学分析、统计、机器学习和科研绘图函数。

这份 README 的用途是：隔一段时间没有使用后，能够快速想起包放在哪里、怎样
安装、数据方向如何约定，以及去哪里查函数。

## 当前信息

| 项目 | 当前设置 |
| --- | --- |
| 版本 | `0.0.2` |
| 源码目录 | `F:/R_proj/mengR` |
| R 版本 | R 4.4.3 |
| 个人包目录 | `E:/SoftwareData/R/win-library/4.4` |
| 默认数据库目录 | `F:/database` |

## 使用前先记住

- `profile` 一般是 **feature × sample**：行为 feature，列为 sample。
- `group` 一般至少包含样本列和分组列，默认名称为 `sample` 和 `group`。
- 函数名全部使用小写 `snake_case`，例如 `calcu_pca()`、`plot_pcoa()`。
- 绘图函数通常返回 ggplot、patchwork 或 ComplexHeatmap 对象，不自动保存图片。
- 专业依赖按需安装；例如 `cellchat_run()` 需要 `CellChat`，单细胞函数通常需要
  `Seurat`，生态统计函数通常需要 `vegan`、`phyloseq` 或 `picante`。

## 安装或更新本地包

已经构建好的源码包位于包目录时：

```r
lib_dir <- "E:/SoftwareData/R/win-library/4.4"

install.packages(
  "F:/R_proj/mengR/mengR_0.0.2.tar.gz",
  repos = NULL,
  type = "source",
  lib = lib_dir
)

.libPaths(c(lib_dir, .libPaths()))
library(mengR)
packageVersion("mengR")
```

只有明确需要更新已安装版本时才执行安装；平时修改源码不必反复覆盖本地包。

## 常用函数入口

| 任务 | 常用函数 |
| --- | --- |
| profile 转换 | `profile_trans_ra()`、`profile_trans_log2()`、`profile_trans_clr()` |
| profile 筛选与汇总 | `profile_filter()`、`profile_top_n()`、`profile_aggregate()` |
| 分类组成 | `taxa_trans()`、`plot_compos()`、`plot_taxa_boxplot()` |
| alpha/beta diversity | `calcu_alpha()`、`calcu_distance()`、`calcu_beta()` |
| 降维 | `calcu_pca()`、`plot_pcoa()`、`plot_nmds()`、`plot_dbrda()` |
| 差异与相关性 | `difference_analysis()`、`calcu_correlation()`、`calcu_mantel()` |
| ROC 与机器学习 | `plot_roc()`、`lasso_kfold()`、`rf_kfold()`、`svm_kfold()` |
| 单细胞 | `scrna_marker_match()`、`scrna_expression_dotplot()`、`cellchat_run()` |
| 宏基因组与代谢组 | `profile_kegg_trans()`、`tidy_lefse()`、`extract_hmdb_xrefs()` |
| 调色板 | `pald()`、`palc()`、`pald_show()`、`palc_show()` |

全部函数按领域整理在 [FUNCTION_REFERENCE.md](FUNCTION_REFERENCE.md) 中。具体参数和
返回值直接查看帮助：

```r
help(package = "mengR")
?calcu_pca
?profile_trans_ra
```

## 最小示例

```r
library(mengR)

profile_df <- data.frame(
  sample_a = c(10, 30, 60),
  sample_b = c(20, 20, 60),
  row.names = c("feature_1", "feature_2", "feature_3"),
  check.names = FALSE
)

profile_ra_df <- profile_trans_ra(profile_df, base = 1)
pca_result <- calcu_pca(profile_ra_df)

pca_result$points
pald("Set2", n = 5)
```

## 数据库路径

部分注释和数据库函数从公共数据库根目录读取文件：

```r
mengr_config(database = "F:/database")
mengr_config()
```

还可以通过 R option `mengR.database` 或环境变量 `MENGR_DATABASE` 设置。当前默认值
是 `F:/database`。

## 隔久了以后从哪里看

1. 先看 [FUNCTION_REFERENCE.md](FUNCTION_REFERENCE.md)，找到需要的函数。
2. 再运行 `?function_name`，确认参数、输入方向和返回对象。
3. 如果旧代码报函数不存在，看 [MIGRATION_NOTES.md](MIGRATION_NOTES.md)。
4. 想知道 0.0.2 改了什么，看 [NEWS.md](NEWS.md)。
5. 新写或整理函数时，看 [STYLE_GUIDE.md](STYLE_GUIDE.md)。

## 我自己维护时

帮助文档写在 `R/` 源码的 Roxygen 注释中，不直接修改 `man/*.Rd`。需要重新生成时，
在包根目录运行：

```r
roxygen2::roxygenise(
  ".",
  roclets = c("rd", "namespace"),
  load_code = roxygen2::load_source
)
```

快速运行源码测试：

```powershell
$env:R_LIBS_USER = "E:/SoftwareData/R/win-library/4.4"
& "D:\R\R-4.4.3\bin\Rscript.exe" "tools/smoke-test.R"
```

准备安装或做较大更新时，再构建和检查：

```powershell
$env:R_LIBS_USER = "E:/SoftwareData/R/win-library/4.4"
& "D:\R\R-4.4.3\bin\R.exe" CMD build --no-manual --no-build-vignettes .
& "D:\R\R-4.4.3\bin\R.exe" CMD check --no-manual mengR_0.0.2.tar.gz
```

代码风格和函数记录按照个人 skill `$mengr-package-style` 维护。

## License

MIT License，见 [LICENSE](LICENSE)。
