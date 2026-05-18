# Mass Analyzer

Mass Analyzer 是一个基于 R Shiny 的蛋白候选分析工具。项目当前包含一个主应用文件 `protein_candidate_explorer_shiny_app_v29.R`，用于从已筛选的基因/蛋白列表或原始丰度表出发，完成候选筛选、差异分析、GO/KEGG 富集、火山图和金属离子相关注释。

## 功能概览

- Gene-list mode：面向已经筛好的基因/蛋白表，自动提取 GeneName，进行离子注释和 BP/MF/CC/KEGG 富集分析。
- Raw-data mode：面向原始 abundance 表，完成质控、污染物过滤、重复基因处理、差异分析、Up/Down 候选筛选、火山图、富集分析和离子注释。
- 支持 `.xlsx`、`.xls`、`.csv`、`.tsv`、`.txt` 输入。
- 支持导出表格、富集结果、火山图和 dotplot，图形可选 PNG、JPG、EMF、PDF。
- 使用 `renv` 锁定 R 包环境，便于在不同电脑上复现运行环境。

## 项目结构

```text
.
├── .Rprofile
├── protein_candidate_explorer_shiny_app_v29.R
├── renv/
└── renv.lock
```

主要文件说明：

- `protein_candidate_explorer_shiny_app_v29.R`：Shiny 应用主文件，包含界面、分析流程和导出逻辑。
- `renv.lock`：项目依赖锁定文件。
- `.Rprofile`：进入项目时自动启用 `renv` 环境。
- `renv/`：`renv` 项目环境相关文件。

## 环境要求

- R 4.6.0 或兼容版本
- Bioconductor 3.23 或兼容版本
- RStudio 或可运行 R 脚本的终端环境

主要依赖包括：

- Shiny/UI：`shiny`、`DT`
- 数据读取与处理：`readxl`、`readr`、`dplyr`、`stringr`
- 绘图与导出：`ggplot2`、`openxlsx`、`devEMF`
- 注释与富集：`AnnotationDbi`、`org.Hs.eg.db`、`GO.db`、`clusterProfiler`、`enrichplot`

## 安装依赖

首次打开项目后，建议先恢复 `renv` 环境：

```r
install.packages("renv")
renv::restore()
```

如果不使用 `renv`，也可以手动安装依赖：

```r
install.packages(c(
  "shiny", "DT", "readxl", "readr", "dplyr", "stringr",
  "ggplot2", "openxlsx", "devEMF"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "AnnotationDbi", "org.Hs.eg.db", "GO.db",
  "clusterProfiler", "enrichplot"
))
```

## 运行应用

在项目根目录运行：

```r
shiny::runApp("protein_candidate_explorer_shiny_app_v29.R")
```

也可以在 RStudio 中打开 `protein_candidate_explorer_shiny_app_v29.R`，点击 `Run App`。

## 输入数据要求

### Gene-list mode

适用于已经筛选好的基因或蛋白列表。

基本要求：

- 输入文件格式支持 `.xlsx`、`.xls`、`.csv`、`.tsv`、`.txt`。
- 表格中需要能识别到蛋白描述列，默认会尝试查找 `Description`、`description`、`Protein Description`、`protein_description`。
- 如果描述列里包含类似 `GN=TP53` 的字段，应用会自动提取 `GeneName`。
- 如果自动识别失败，可以在界面中手动填写 Description 列名。

主要输出：

- 输入预览
- Ion Protein summary
- Ion GO hits
- Annotated 原表
- GO BP/MF/CC 富集表和 dotplot
- KEGG 富集表和 dotplot

### Raw-data mode

适用于原始蛋白丰度表。

基本要求：

- 输入文件格式支持 `.xlsx`、`.xls`、`.csv`、`.tsv`、`.txt`。
- 丰度列需要命名为 `Test1`、`Test2`、... 和 `Ctrl1`、`Ctrl2`、...。
- 推荐包含可识别的 Description 列，用于提取 `GeneName`。
- 如果存在 `Protein FDR Confidence: Combined` 或类似 confidence 列，可以在界面中选择过滤范围。

Raw-data mode 包含三种分析模式：

- 模式 A：宽松候选筛选。空白先补 0，Test 组均值只使用有效正值，Ctrl 组按当前规则将 0 计入均值。
- 模式 B：平衡模式。对 Test 和 Ctrl 两组做更均衡的清洗和比较。
- 模式 C：严格模式。在模式 B 基础上结合 pvalue 和 fold change 阈值筛选差异候选。

主要处理步骤：

- 自动识别 Test/Ctrl 丰度列
- Protein FDR Confidence 过滤
- 常见污染物或 Keratin 污染物过滤
- 可选按 `GeneName` 去重，保留整体丰度最高的一行
- 重复值清洗和异常值处理
- Fold change、log2FC、pvalue 计算
- Up/Down 候选筛选
- 火山图绘制
- GO/KEGG 富集分析
- calcium、magnesium、manganese、zinc 离子注释

主要输出：

- 差异结果总表
- Filtered out 清单
- Contaminants removed 清单
- Duplicate 处理记录
- Annotated 差异结果
- 当前模式 Up/Down 列表
- 火山图
- Up/Down 的 GO BP/MF/CC 和 KEGG 富集结果
- `GO enrichment.xlsx`
- `KEGG enrichment.xlsx`
- `Ion annotation.xlsx`

## 常用参数说明

- `Fold change cutoff`：差异倍数阈值，默认 2。
- `Differential pvalue cutoff`：差异分析 pvalue 阈值，默认 0.05。
- `Enrichment pvalue cutoff`：富集分析 pvalue 阈值，默认 0.05。
- `GO qvalue cutoff`：GO 富集 qvalue 阈值，默认 0.2。
- `dotplot 显示条目数`：富集 dotplot 中展示的条目数量。
- `Contaminant filter`：污染物过滤方式，可选择 293T 常见污染物、Keratin only 或不过滤。
- `Enrichment 策略`：`Exploratory enrichment` 更适合候选解释，`Context-aware enrichment` 会引入本次分析背景。
- `Volcano display style`：控制火山图展示方式，不改变候选分类逻辑。

## 注意事项

- KEGG 富集可能需要联网访问 KEGG 服务；如果网络不可用，GO 分析仍可正常尝试运行。
- Raw-data mode 默认只计算当前分析模式。若需要在火山图中跨模式切换，请勾选“预计算全部模式”。
- GeneName 的提取依赖 Description 列中的 `GN=` 字段；如果源数据没有该字段，需要先补充基因名信息或调整输入表。
- EMF 导出依赖 `devEMF`，在不同系统上的可用性可能不同；论文或汇报后期编辑优先尝试 PDF 或 EMF。

## 更新依赖

如果修改了依赖，运行：

```r
renv::snapshot()
```

如果需要在新环境恢复依赖，运行：

```r
renv::restore()
```

