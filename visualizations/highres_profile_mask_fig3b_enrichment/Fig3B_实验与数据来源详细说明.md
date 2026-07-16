# Fig. 3B 实验与数据来源详细说明

## 1. 面板身份与替换关系

本目录中的 Fig. 3B 是针对高分辨率、细胞级 profile-masking 实验重新计算并生成的正式替代面板：

- PNG：`fig3b_highres_profile_mask_enrichment.png`
- 可编辑 SVG：`fig3b_highres_profile_mask_enrichment.svg`
- PDF：`fig3b_highres_profile_mask_enrichment.pdf`
- 生成脚本：`scripts/build_fig3b_highres_profile_mask_enrichment.py`
- 冻结配置：`configs/fig3b_highres_profile_mask_enrichment.json`

原论文 PDF 中可见的 Peak ES `0.068` 和 `0.080` 缺少完整 source-value 链，不能继续作为仓库验证值。本次替代面板从当前表达矩阵和当前映射结果重新计算，最终值为：

- CytoSPACE Peak ES：`0.07404515018728003`
- SVTuner Peak ES：`0.09147893720835686`
- SVTuner - CytoSPACE：`0.017433787021076827`

图中按三位小数显示为 `0.074` 和 `0.091`。

## 2. 生物学问题

该面板检验的问题是：在目标细胞类型的空间表达 profile 被强烈削弱后，映射结果中较低的 target-like signal 是否更集中于残余目标 marker support 较低的空间单元。

这里的 target-like signal 不是外部病理真值，也不是细胞类型分类准确率。它是根据映射到每个空间单元的参考细胞组成以及这些参考细胞类型的目标 marker 表达构建的映射后代理分数。

因此，图中指标必须称为：

```text
Peak enrichment score (Peak ES)
```

不应称为 NES、accuracy、AUROC 或 ground-truth recovery accuracy。

## 3. 数据来源

### 3.1 公开数据来源

五个高分辨率数据集来自 Vizgen 2022 年 5 月发布的 FFPE Human Immuno-oncology data release，平台为 MERSCOPE/MERFISH，项目审计将其记录为 500-gene FFPE MERSCOPE 数据。公开来源入口为：

- Vizgen FFPE Human Immuno-oncology release：<https://vizgen.com/human-ffpe-immunooncology-release-roadmap/>
- 项目审计辅助引用：<https://doi.org/10.1038/s41587-025-02811-9>

`HumanBreastCancerPatient1`、`HumanColonCancerPatient1`、`HumanLungCancerPatient1`、`HumanMelanomaPatient1` 和 `HumanMelanomaPatient2` 是该公开 release 中使用的正式数据发布标识。它们不应被解释成项目获得的临床患者编号或独立 accession。

仓库内公开来源审计记录位于：

```text
visualizations/manuscript_public_source_resolution/public_source_resolution.csv
```

### 3.2 五个 profile-mask 场景

| 数据发布标识 | profile-mask 目标 | 空间单元数 | 原空间数据中的目标细胞数 | Stage3A 检测缺失类型 | 移除参考细胞数 |
|---|---:|---:|---:|---|---:|
| HumanBreastCancerPatient1 | Monocytes and Macrophages | 2,561 | 411 | Monocytes and Macrophages | 820 |
| HumanColonCancerPatient1 | Fibroblasts | 2,561 | 459 | Fibroblasts | 746 |
| HumanLungCancerPatient1 | Plasma cells | 2,561 | 141 | Plasma cells | 232 |
| HumanMelanomaPatient1 | Fibroblasts | 2,561 | 928 | Fibroblasts | 1,876 |
| HumanMelanomaPatient2 | B cells | 2,561 | 150 | B cells; NK cells | 294 |

场景清单和上述数字来自：

```text
visualizations/highres_profile_mask_mapping/highres_profile_mask_mapping_manifest.csv
visualizations/highres_profile_mask_mapping/highres_profile_mask_mapping_manifest.json
```

### 3.3 Fig. 3B 代表案例

主图显示 `HumanMelanomaPatient2` 的 B-cell profile-mask 场景。该案例来自历史 Fig. 3B 脚本所使用的代表场景，并在本次重算配置中明确冻结，不是依据本次最终 Peak ES 再从五个结果中临时挑选。

代表案例的输入构成为：

- 2,561 个细胞级空间单元，仓库字段名为 `spot_id`；
- 5,200 个同平台、同 assay、与空间单元原始细胞不重叠的参考细胞；
- `spatial_sc_raw_overlap = 0`；
- 8 个参考细胞类型；
- 原参考中 B cells 为 228 个，NK cells 为 66 个；
- 原空间选择中 B cells 为 150 个；
- profile-mask 缩放系数为 `0.02`。

参考细胞构成为：

| 参考细胞类型 | 细胞数 |
|---|---:|
| B cells | 228 |
| Endothelial cells | 282 |
| Epithelial cells | 1,260 |
| Fibroblasts | 1,262 |
| Monocytes and Macrophages | 970 |
| NK cells | 66 |
| Plasma cells | 288 |
| T cells | 844 |

代表案例的 profile-mask 元数据位于：

```text
data/processed/highres_humanmelanomapatient2_profile_mask_b_cells/
  stage1_preprocess/fig2d_profile_mask_info.json
```

## 4. Profile masking

profile masking 针对目标细胞类型的 marker profile 进行表达缩放。代表案例将 B-cell target profile 按 `mask_scale = 0.02` 处理，即保留原尺度的 2%。

代表案例存储的 14 个 target-panel 通道为：

```text
MS4A1, CD79B, TNFRSF13C, CD22, BLK, CD19, CD79A,
PAX5, TNFRSF13B, FCER2, CXCR5, STAP1, CR2, Blank-34
```

`Blank-34` 是数据 panel 中保留的 blank-control 通道。当前重算忠实使用仓库冻结的 panel，没有在看过结果后删除该通道。其他数据集最多使用其 panel 中前 30 个、同时存在于表达矩阵的通道；肺癌场景实际使用 18 个，melanoma-2 场景实际使用 14 个。

## 5. 两条映射路线

### 5.1 CytoSPACE baseline

baseline 使用完整参考，不执行缺失类型过滤：

```text
filter_mode = none
cell_type_column = sc_meta
missing_type = __NO_MISSING__
stage4_suffix = _baseline_highres
```

### 5.2 SVTuner route

SVTuner route 先运行 Stage3A 类型支持检测，再使用插件输出过滤检测到的缺失参考类型：

```text
filter_mode = plugin_unknown
cell_type_column = plugin_type
missing_type = __AUTO__
filter_scope = missing_detected_only
stage4_suffix = _route2_highres
```

两条路线共同使用：

```text
mapping_cells_per_spot = 1
sc_expr_source = normalized
no_sampling_sub_spots = true
n_processors = 1
```

在代表案例中，Stage3A 自动检测到 B cells 和 NK cells 为 missing types。5,200 个参考细胞中保留 4,906 个，移除 294 个，其中包括 228 个 B cells 和 66 个 NK cells。

需要特别说明：该 Fig. 3B 比较的是 Stage3A 参考过滤后的 SVTuner route 与 CytoSPACE baseline，不包含 Stage3B 空间 unsupported-region blanking。

Stage3A 审计文件位于：

```text
result/highres_humanmelanomapatient2_profile_mask_b_cells/
  stage3_typematch/stage3_summary.json
```

映射输入位于：

```text
result/<profile_mask_sample>/stage4_cytospace_baseline_highres/
  cytospace_output/cell_type_assignments_by_spot.csv

result/<profile_mask_sample>/stage4_cytospace_route2_highres/
  cytospace_output/cell_type_assignments_by_spot.csv
```

## 6. Target-like score 的计算

对每个场景，从冻结的 target marker panel 中读取最多 30 个可用通道，记为集合 \(G\)。

### 6.1 残余空间 support

对空间单元 \(i\)，残余 support 为归一化空间表达的 marker 均值：

\[
S_i = \frac{1}{|G|}\sum_{g \in G} X^{ST}_{ig}.
\]

随后在该数据集内部做 min-max 标准化：

\[
\widetilde S_i = \frac{S_i-\min(S)}{\max(S)-\min(S)+10^{-12}}.
\]

所有空间单元按 \(\widetilde S_i\) 从低到高稳定排序，构成横轴 `Relative support rank`。

### 6.2 参考细胞类型 marker 分数

对每个参考细胞先计算 marker 均值，再按参考细胞类型 \(t\) 求平均，得到：

\[
M_t = \operatorname{mean}_{c \in t}
\left(\frac{1}{|G|}\sum_{g \in G}X^{SC}_{cg}\right).
\]

### 6.3 映射后 target-like score

设 CytoSPACE 在空间单元 \(i\) 上分配的类型 \(t\) 细胞数为 \(A_{it}\)，总分配数为 \(N_i\)。映射后 target-like score 为：

\[
R_i = \frac{\sum_t A_{it}M_t}{\max(N_i,1)}.
\]

baseline 和 SVTuner 分别从各自的 `cell_type_assignments_by_spot.csv` 计算 \(R_i\)。

## 7. Bottom-10% hits 与 Peak ES

### 7.1 精确 10% 有效 hit mass

每种方法在本数据集内部取 target-like score 最低的 10%。由于映射分数是离散值，10% 分位点存在大量并列。直接使用 `score <= quantile(0.1)` 会错误地纳入远多于 10% 的空间单元。

当前实现采用 fractional hit mass：

- 低于阈值的单元 hit weight 为 1；
- 位于阈值并列组的单元共享剩余 hit mass；
- 每种方法的总有效 hit mass 严格为 `0.1 x 2561 = 256.1`。

设单元 \(i\) 的 hit weight 为 \(h_i\)，总 hit mass 为 \(H=256.1\)，空间单元总数为 \(n=2561\)。运行增量为：

\[
\Delta_i = \frac{h_i}{H} - \frac{1-h_i}{n-H}.
\]

运行富集曲线为：

\[
ES(k)=\sum_{i=1}^{k}\Delta_i.
\]

最终指标为：

\[
Peak\ ES=\max_k ES(k).
\]

图形绘制时对运行曲线使用约 15 个空间单元的居中 rolling mean，仅用于视觉展示；CSV 和图中文字中的 Peak ES 均来自未平滑的原始运行曲线。

## 8. 图形元素解释

- 上方曲线：CytoSPACE baseline；
- 下方曲线：SVTuner Stage3A route；
- 橙色曲线：沿 residual-support rank 的 running ES；
- 橙色竖线：低 target-like score hit，浅色竖线表示阈值并列组的 fractional weight；
- 灰色底部轨道：min-max 标准化后的 residual target-marker support；
- 左端：较低 masked support；
- 右端：较高 masked support。

Peak ES 越高，表示低 target-like score 越集中于 residual support 较低的一端。它只能说明排序一致性增强，不能证明绝对定位正确，也不能替代外部组织学真值。

## 9. 五个数据集的完整结果

| 数据集 | mask 目标 | marker 数 | CytoSPACE Peak ES | SVTuner Peak ES | 差值 |
|---|---|---:|---:|---:|---:|
| HumanBreastCancerPatient1 | Monocytes and Macrophages | 30 | 0.112862 | 0.085565 | -0.027297 |
| HumanColonCancerPatient1 | Fibroblasts | 30 | 0.069460 | 0.045717 | -0.023743 |
| HumanLungCancerPatient1 | Plasma cells | 18 | 0.066178 | 0.064902 | -0.001275 |
| HumanMelanomaPatient1 | Fibroblasts | 30 | 0.250535 | 0.157611 | -0.092925 |
| HumanMelanomaPatient2 | B cells | 14 | 0.074045 | 0.091479 | +0.017434 |

只有用于主图的历史代表案例为正向差值，其余四个数据集均为负向差值。因此论文中必须将 Fig. 3B 描述为 `representative case` 或 `representative high-resolution profile-mask example`，不能写成五个数据集一致改善，也不能将该面板单独作为总体优越性的证据。

完整数值位于：

```text
fig3b_highres_profile_mask_enrichment_metrics.csv
fig3b_highres_profile_mask_enrichment_source_values.csv
fig3b_highres_profile_mask_enrichment_selected_source_values.csv
```

## 10. 推荐图例文字

```text
b, Representative high-resolution profile-mask case in the Vizgen
HumanMelanomaPatient2 FFPE MERSCOPE dataset. Spatial units were ranked by
residual B-cell marker support after profile masking. Low mapped target-like
scores were evaluated using an exact 10% fractional hit mass to account for
tied mapping scores. Peak ES increased from 0.074 for CytoSPACE to 0.091 for
the Stage3A-filtered SVTuner route. This panel is a representative case and
does not imply consistent improvement across all five high-resolution datasets.
```

## 11. 推荐 Methods 表述

```text
For the cell-resolution profile-mask enrichment analysis, spatial units were
ranked by the mean residual expression of up to 30 masked-target panel genes.
A mapped target-like score was calculated as the assigned-cell-count-weighted
mean reference marker score at each spatial unit. The lowest-scoring 10% was
represented by an exact fractional hit mass, with equal fractional weights for
scores tied at the boundary. Peak ES was defined as the maximum unweighted
running enrichment score along the residual-support rank. The displayed case
was the historically selected HumanMelanomaPatient2 B-cell profile-mask
scenario and was fixed before the current rerun.
```

## 12. SVG 可编辑性

SVG 使用 Matplotlib 的：

```text
svg.fonttype = none
```

因此：

- 标题、方法名、轴标签和数值保留为 `<text>`；
- 曲线、坐标轴、hit ticks 和 support track 保留为矢量路径；
- 不含 `<image>` 嵌入位图；
- 不嵌入 base64 字体；
- 可在 Adobe Illustrator 中分别选择和编辑文本、曲线与线段。

当前自动检查记录为 14 个 text 元素、1,320 个 path 元素、0 个 raster image 和 0 个 embedded font。

## 13. 复现与审计

在项目根目录运行：

```bash
python scripts/build_fig3b_highres_profile_mask_enrichment.py --project_root .
```

完整 manifest 位于：

```text
fig3b_highres_profile_mask_enrichment_manifest.json
```

manifest 记录：

- 生成脚本和冻结配置的 SHA256；
- 五个数据集的 marker panel、ST/SC 表达矩阵和 SC metadata 的 SHA256；
- baseline 与 route2 assignment 文件的 SHA256；
- PNG、SVG、PDF、metrics 和 source-values 的 SHA256；
- 代表案例、选择规则、指标定义和精确数值。

这套记录用于保证 Fig. 3B 的图片、数值、输入和实验定义可以相互追溯。
