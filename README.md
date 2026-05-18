# Mass Analyzer

Mass Analyzer 是一个基于 R Shiny 的蛋白候选分析工具。项目当前的主应用文件是 `app.R`，用于从已筛选的基因/蛋白列表或原始丰度表出发，完成候选筛选、差异分析、GO/KEGG 富集、火山图和金属离子相关注释。

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
├── app.R
├── README.md
├── renv/
├── renv.lock
└── rsconnect/
```

主要文件说明：

- `app.R`：Shiny 应用主文件，包含界面、分析流程和导出逻辑。
- `renv.lock`：项目依赖锁定文件。
- `.Rprofile`：进入项目时自动启用 `renv` 环境。
- `renv/`：`renv` 项目环境相关文件。
- `rsconnect/`：Shiny 部署相关配置目录。

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
shiny::runApp()
```

也可以在 RStudio 中打开 `app.R`，点击 `Run App`。

## 程序运行逻辑

应用启动后会先检查必需 R 包是否可用；缺包时会停止运行并提示安装命令。随后加载 UI 与 server，最后通过 Shiny 启动两个工作流页面：`Gene-list mode` 和 `Raw-data mode`。

通用读表逻辑：

- 自动读取 Excel、CSV、TSV 或 TXT 文件。
- 自动修复重复列名。
- 自动识别 Description 列，优先查找 `Description`、`description`、`Protein Description`、`protein_description`。
- 从 Description 文本里的 `GN=xxx` 字段提取 `GeneName`。

## Gene-list mode

适用于已经筛选好的基因或蛋白列表。

基本要求：

- 输入文件格式支持 `.xlsx`、`.xls`、`.csv`、`.tsv`、`.txt`。
- 表格中需要能识别到蛋白描述列，或者手动填写 Description 列名。
- Description 列中建议包含类似 `GN=TP53` 的字段，用于提取 `GeneName`。

运行步骤：

1. 上传已筛选好的基因/蛋白表。
2. 自动或手动指定 Description 列。
3. 选择要注释的离子：calcium、magnesium、manganese、zinc。
4. 设置 enrichment pvalue cutoff、GO qvalue cutoff 和 dotplot 展示条目数。
5. 点击运行，程序会完成离子注释和 BP/MF/CC/KEGG 富集。

主要输出：

- 输入预览
- Ion Protein summary
- Ion GO hits
- Annotated 原表
- GO BP/MF/CC 富集表和 dotplot
- KEGG 富集表和 dotplot

## Raw-data mode

适用于原始蛋白丰度表。

基本要求：

- 输入文件格式支持 `.xlsx`、`.xls`、`.csv`、`.tsv`、`.txt`。
- 丰度列必须命名为 `Test1`、`Test2`、... 和 `Ctrl1`、`Ctrl2`、...。
- Test 和 Ctrl 至少各 2 列。
- 推荐包含可识别的 Description 列，用于提取 `GeneName`。
- 如果存在 `Protein FDR Confidence: Combined` 或类似 confidence 列，可以在界面中选择过滤范围。

Raw-data mode 的处理顺序：

1. 读取原始丰度表。
2. 从 Description 中提取 `GeneName`。
3. 按 Protein FDR Confidence 过滤，可选 High only、High + Medium 或 All。
4. 按污染物规则过滤，可选 293T 常见污染物、Keratin only 或不过滤。
5. 自动识别 `Test1...TestN` 和 `Ctrl1...CtrlN` 丰度列。
6. 将丰度列安全转换为数字。
7. 检查重复 GeneName，可选择保留整体丰度最高的一行。
8. 对每个 GeneName 逐行计算均值、FC、log2FC、pvalue 和重复质量标记。
9. 根据当前分析模式 A/B/C 标记 up、down、presence_up、presence_down、ns 或 filtered_out。
10. 生成火山图坐标、Up/Down 候选列表、GO/KEGG 富集和离子注释。

### 分析模式

- 模式 A：宽松候选筛选。缺失值先补 0；Test 组均值只使用有效正值；Ctrl 组均值会把 0 计入。适合发现只在 Test 中出现的候选。
- 模式 B：平衡模式。保留 NA，不直接补 0；Test 和 Ctrl 都做重复清洗；主要根据 fold change 判断 Up/Down。
- 模式 C：严格模式。沿用模式 B 的清洗和计算，但 Up/Down 还必须满足 pvalue 阈值。

### 状态分类

- `up`：Test 和 Ctrl 都有信号，且 `log2FC >= log2(Fold change cutoff)`。
- `down`：Test 和 Ctrl 都有信号，且 `log2FC <= -log2(Fold change cutoff)`。
- `presence_up`：Test 有信号，Ctrl 没有信号，且重复数足够。
- `presence_down`：Ctrl 有信号，Test 没有信号，且重复数足够。
- `ns`：不满足 Up/Down/Presence 条件。
- `filtered_out`：重复不足、重复差异过大、清洗后不符合规则等。

### 污染物过滤

默认使用 `Common contaminants (293T default)`，会过滤常见 keratin、trypsin、Lys-C、albumin/BSA/casein、hemoglobin、immunoglobulin 等污染物。也可以选择只过滤 Keratin，或完全不过滤。被过滤的行会进入 `contaminants removed` 清单。

### 重复基因处理

如果勾选“按 GeneName 去重”，同一个 GeneName 多行时，程序会根据 Test/Ctrl 丰度列计算整体丰度得分，并保留整体丰度最高的一行。处理记录会输出到 duplicate 处理表。

### 火山图规则

火山图默认使用当前分析模式的 `active_log2FC` 和 `active_neglog10p`。Presence 类型点由于一组没有信号，程序提供固定最小值补值或左尾随机补值用于可视化；这些补值只影响画图位置，不改变候选分类。

### 富集和离子注释

富集分析会把 Gene Symbol 转换为 Entrez ID，然后分别运行 GO BP、GO MF、GO CC 和 KEGG enrichment。Raw-data mode 会对 Up 和 Down 分开富集。

离子注释基于 GO term 关键词匹配，当前支持 calcium、magnesium、manganese、zinc。输出包括每个离子的 GO term summary、Protein summary 和 GO hits。

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
- `Enrichment 策略`：`Exploratory enrichment` 更适合候选解释；`Context-aware enrichment` 会把当前分析中未 filtered out 的基因作为背景。
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
