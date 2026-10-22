# Raw RCC + DESeq2 完整分析流程

这个版本保留原 R 脚本的样本分组和后续分析顺序，只替换了错误的表达数据入口。原脚本没有被覆盖。

## 已连接的流程

1. 从 `GSE224450_RAW.tar` 解压并读取原始 `.RCC.gz` 文件。
2. 检查计数必须为非负整数，不进行四舍五入。
3. 保留原脚本的 5 个 Low_CP 和 17 个 High_CP 样本。
4. 使用原脚本过滤规则：至少两个样本中的计数不低于 10。
5. 仅使用 housekeeping probes 估计 DESeq2 size factors。
6. 用 DESeq2 比较 `High_CP - Low_CP`。
7. 差异分析完成后才进行 rlog，rlog 只用于后续作图和表达分析。
8. 继续执行 PCA、火山图、差异基因热图、GO、KEGG、免疫浸润、hub genes、机器学习、WGCNA、通路相关和生存分析。

所有 R 代码注释均为简单英文。

## 本次运行结果

- 样本：5 Low_CP，17 High_CP
- 通过过滤的 endogenous probes：730
- 严格阈值：`padj < 0.01` 且 `|log2FC| > 2`
- 显著探针：78
- 上调：11
- 下调：67
- GO：上调305条，下调754条
- KEGG：上调24条，下调39条
- WGCNA：7个模块，soft power为7
- 生存分析：17个基因
- 机器学习：6个模型完成；XGBoost训练失败并被单独记录
- 免疫浸润：CIBERSORT、MCP-counter、EPIC、quanTIseq和ESTIMATE共5种方法完成

结构测试和完整输出测试均已通过。九个模块均成功完成，所有结果表和PNG图片均可读取且非空。

## 一键运行

从本目录（`sequencing/rna`）运行；脚本通过自身路径自动定位项目根：

```bash
Rscript run_full_pipeline.R
```

运行后可执行：

```bash
Rscript tests/validate_output.R
```

也可以单独运行模块：

```bash
Rscript run_analysis.R --modules=data,figures,enrichment
```

## 数据目录

输入数据（GSE224450 原始 RCC、GSE102349 表型与表达文件等）位于仓库之外，通过环境变量指定；脚本不会回退到任何个人电脑路径：

- `RNA_INPUT_DIR`：输入数据目录（**必设**；未设置时脚本尝试仓库内 `data/inputs/`，缺失则报错并提示）
- `RNA_OUTPUT_DIR`：输出目录（默认 `<项目根>/output`）
- `RNA_SERIES_MATRIX`：敏感性分析使用的 GEO 系列矩阵文件（默认 `<RNA_INPUT_DIR>/GSE224450_series_matrix.txt.gz`）
- `RNA_EXTERNAL_LIBRARY`（可选）：额外的 R 包目录，设置后加入 `.libPaths()`

示例：

```bash
RNA_INPUT_DIR=/path/to/inputs Rscript run_full_pipeline.R
```

## 需要谨慎解释的部分

### CIBERSORT

如本地 R 包未满足依赖，可设置 `RNA_EXTERNAL_LIBRARY` 指向自有 R 包目录，脚本会在 `.libPaths()` 中加入该目录。CIBERSORT使用本次DESeq2产生的非对数标准化计数、LM22 signature、1000次置换和`QN = FALSE`重新计算，没有复用旧免疫结果。详细参数、signature覆盖率和警告写入`output/tables/immune_method_audit.csv`。

同时重新运行MCP-counter、EPIC、quanTIseq、xCell和ESTIMATE。CIBERSORT、MCP-counter、EPIC、quanTIseq和ESTIMATE成功；xCell因IO 360靶向panel的signature覆盖不足而失败，并保留在审计表中，未被静默删除。

在多重检验校正后，CIBERSORT显示High_CP组中性粒细胞和活化肥大细胞比例较低。其他方法仅提供未通过FDR校正的方向性证据，例如MCP-counter细胞毒性评分升高、EPIC CD4 T细胞升高及quanTIseq M1巨噬细胞降低。由于IO 360仅覆盖部分去卷积signature，这些结果应作为探索性分析。

### 机器学习

`xgboost` 和 `kernlab` 已安装在交付文件夹内的 `R_library`，运行程序会自动加载这个本地软件包目录。Random Forest、Elastic Net、Lasso、KNN、Ridge Logistic和SVM均已完成。XGBoost在当前 `caret` 交叉验证中产生全部缺失的ROC，已排除，其失败状态写入 `output/tables/model_skips.csv`，没有混入模型比较图。

模型ROC非常高，但只有22个样本，而且 Low_CP 只有5个。交叉验证中的训练样本仍很少，因此这些ROC只适合作为探索性结果，不能当作已完成独立验证的模型性能。

原Python notebook还存在两个问题：它读取旧的硬编码输入文件；并且在每个仅含一个测试样本的LOOCV fold内计算balanced accuracy再取平均。交付的`python_ml/run_ml_analysis.py`改为先收集全部22个留一预测，再计算accuracy、balanced accuracy和ROC AUC，并把标准化限制在每个训练fold内。

Python重跑中，Elastic Net的LOOCV balanced accuracy为0.871、ROC AUC为0.953；LASSO分别为0.841和0.965。LASSO仅保留中性粒细胞和活化肥大细胞两个特征，方向均指向Low_CP。该特征选择是在完整数据上拟合得到，仍需要独立队列验证，不能解释为已验证的生物标志物。

### 原分组

为了保持与原 R 脚本一致，本次继续使用原来的样本ID。其中包含一个 EBV DNA 缺失样本，以及一个位于1045边界的样本。确定最终分组规则后，应再次运行本流程。

### WGCNA

WGCNA可以运行并生成结果，但22个样本仍属于较小队列。模块与分组的相关性应作为探索性结果，而不是稳定网络验证。

## 主要文件

- `run_full_pipeline.R`：一键运行全部模块。
- `run_analysis.R`：模块调度程序。
- `R/01_data.R`：raw RCC读取、housekeeping normalization、DESeq2和rlog。
- `R/02_figures.R`至`R/09_survival.R`：接回的后续分析模块。
- `output/tables/differential_expression.csv`：完整DESeq2结果。
- `output/tables/filtered_raw_rcc_counts.csv`：过滤后的原始计数。
- `output/tables/rlog_expression.csv`：后续分析使用的rlog矩阵。
- `output/tables/deseq2_normalized_counts.csv`：免疫去卷积使用的非对数DESeq2标准化计数。
- `output/tables/raw_rcc_normalization_audit.csv`：housekeeping normalization审计。
- `output/tables/immune_method_audit.csv`：多方法免疫分析状态、覆盖率和警告。
- `output/tables/immune_multimethod_group_comparison.csv`：各方法的组间比较和FDR。
- `output/tables/pipeline_summary.csv`：结果数量汇总。
- `output/tables/pipeline_validation_checks.csv`：完整测试结果。
- `python_ml/output/ML_Model_Performance.csv`：Python LOOCV模型结果。
- `IMMUNE_AND_ML_FINDINGS.md`：免疫和Python机器学习结果摘要。
- `output/figures/`：重新生成的全部图形。
