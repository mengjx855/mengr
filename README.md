# mengR

`mengR` is Jin-Xin Meng's personal R package for reusable bioinformatics,
omics analysis, statistics, machine-learning, and visualization functions.

`mengR` 是 Jin-Xin Meng 的个人 R 工具包，主要整理可复用的生物信息学、组学
数据处理、统计分析、机器学习和绘图函数。

Contact: `jinxmeng@zju.edu.cn` / `mengjx855@163.com`.

Version: `0.0.1` (initial development version).

The package is currently prepared as source files for review and testing. It has
not been built or installed. See [STYLE_GUIDE.md](STYLE_GUIDE.md) for coding
conventions and [MIGRATION_NOTES.md](MIGRATION_NOTES.md) for renamed, removed,
and behavior-changed functions. A categorized overview of every current
function is available in [FUNCTION_REFERENCE.md](FUNCTION_REFERENCE.md).

当前只准备源码和测试文件，尚未构建或安装。代码规范见
[STYLE_GUIDE.md](STYLE_GUIDE.md)，函数改名、删除和行为变化见
[MIGRATION_NOTES.md](MIGRATION_NOTES.md)，全部函数的分类索引见
[FUNCTION_REFERENCE.md](FUNCTION_REFERENCE.md)。

## Source-level test / 源码级测试

Run from the parent directory of `mengR`:

```powershell
$env:R_LIBS_USER='E:/SoftwareData/R/win-library/4.4'
& 'D:\R\R-4.4.3\bin\Rscript.exe' 'mengR/tools/smoke-test.R'
```

The optional `cellchat_run()` workflow additionally requires `CellChat`.

## Help documentation / 帮助文档

Function documentation is written next to the source code with roxygen2.
English usage details are followed by a Chinese functional summary. Generated
help pages are stored in `man/`; edit the roxygen comments in `R/` instead of
editing `.Rd` files directly.

函数帮助使用 roxygen2 写在 `R/` 源码中，英文用法说明后附中文功能摘要。生成的
帮助页位于 `man/`；后续修改应编辑源码中的 roxygen 注释，不要直接编辑 `.Rd`。

Regenerate and audit the documentation without building or installing the
package:

```powershell
$env:R_LIBS_USER='E:/SoftwareData/R/win-library/4.4'
& 'D:\R\R-4.4.3\bin\Rscript.exe' -e `
  "roxygen2::roxygenise('mengR', roclets = c('rd', 'namespace'), clean = TRUE)"
& 'D:\R\R-4.4.3\bin\Rscript.exe' 'mengR/tools/audit-roxygen.R' 'mengR'
```
