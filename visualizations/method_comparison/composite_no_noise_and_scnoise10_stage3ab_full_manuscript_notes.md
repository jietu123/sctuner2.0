# Composite no-noise 与 10% SC-reference noise 实验及可视化论文材料

## 1. 文档用途与证据范围

本文档对应以下两个保留结果目录：

```text
visualizations/method_comparison/composite_no_noise/
visualizations/method_comparison/composite_scnoise10_stage3ab_full/
```

它们分别记录：

1. 无额外 SC 表达噪声时，7 种映射方法在 9 个联合 Stage3A/Stage3B 模拟场景中的全空间组成恢复结果；
2. 在 SC reference 中加入 10% entry-wise within-gene replacement noise 后，同一组 9 个场景和 7 种方法的结果，其中 SVTuner 使用完整 Stage3A + Stage3B + CytoSPACE 路线。

本文所有数值来自上述目录中的 CSV/JSON，而不是从图片人工读取。两个箱线图均为描述性汇总，不包含显著性检验。`n=9` 指 9 个工程化模拟场景，不是 9 个独立患者或生物学重复。

## 2. 核心实验问题

实验检验两个互补问题：

1. 当 SC reference 中含有 ST 不支持的细胞类型时，Stage3A 是否能识别并过滤这些 unsupported-reference cells，从而减少错误空间分配；
2. 当 ST 中存在 SC reference 未覆盖的空间区域时，Stage3B 是否能在映射前识别并留白这些 unsupported spatial regions，从而避免 CytoSPACE 被迫分配错误细胞。

SVTuner 路线不是一种新的基础 mapping solver。本实验中的 SVTuner 指：

```text
Stage3A unsupported-reference filtering
+ Stage3B unsupported-region detection and pre-mapping blanking
+ CytoSPACE mapping on the retained cells and spots
```

CytoSPACE 组指不使用 Stage3A filtering 和 Stage3B abstention 的 baseline mapping。其余比较方法为 Tangram（all genes）、Tangram（marker genes）、novoSpaRc、SpaOTsc 和 CellTrek。

## 3. 联合模拟场景设计

### 3.1 三个数据背景和九个场景

每个数据背景包含 control、single missing 和 double missing 三种 Stage3A 条件，同时固定一个 Stage3B SC-reference dropout target。因此共有 `3 datasets x 3 conditions = 9 scenarios`。

| 数据背景 | Stage3A control | Stage3A single missing | Stage3A double missing | Stage3B target（从 SC reference 删除） |
|---|---|---|---|---|
| Real BRCA | 无额外 ST-missing type | Epithelial cells | Epithelial cells + PCs | Endothelial cells |
| Human lung 5-location | 无额外 ST-missing type | AT2 | AT2 + Fibroblast | B_cell |
| Mouse brain refined | 无额外 ST-missing type | Micro | Micro + Oligo_2 | Ext_L56 |

这里的 Stage3A missing type 表示该类型仍存在于 SC reference，但在对应模拟 ST 中缺乏支持；Stage3B target 则相反，目标类型保留在 ST simulation truth 中，但从 SC reference 删除。

九个无噪声 sample ID 为：

```text
real_brca7_endothelial_marker_control_sc_missing_endothelial_cells
real_brca7_endothelial_marker_missing_epithelial_cells_sc_missing_endothelial_cells
real_brca7_endothelial_marker_missing_epithelial_cells_pcs_sc_missing_endothelial_cells

human_lung_5loc_fine9_clustered_sim_sc_missing_b_cell
human_lung_5loc_fine9_clustered_sim_missing_at2_sc_missing_b_cell
human_lung_5loc_fine9_clustered_sim_missing_at2_fibroblast_sc_missing_b_cell

mouse_brain_refined7_balanced_clustered_sim_sc_missing_ext_l56
mouse_brain_refined7_balanced_clustered_sim_missing_micro_fill_inh_pvalb_sc_missing_ext_l56
mouse_brain_refined7_balanced_clustered_sim_missing_micro_oligo_2_fill_inh_pvalb_sc_missing_ext_l56
```

10% noise 版本在相同 ID 后增加 `_scnoise10`。每个场景的 simulation truth 位于：

```text
data/sim/<group>/<sample>/sim_truth_spot_type_fraction.csv
```

场景构造信息位于同目录的 `sim_info.json`。其中 `missing_types` 记录 Stage3A 工程化缺失类型，`sc_reference_drop_types` 记录 Stage3B target。评分程序只使用这些字段定义事后 simulation-truth evaluation，不将 truth label 输入 Stage3A 或 Stage3B 的检测步骤。

### 3.2 Stage3A 和 Stage3B 的运行方式

无噪声 SVTuner 的 Stage4 目录为：

```text
result/<sample>/stage4_cytospace_stage3b_blank/
```

10% noise 完整 SVTuner 的 Stage3A、Stage3B 和 Stage4 使用隔离后缀：

```text
stage3_typematch_scnoise10_stage3ab_full
stage3b_st_unsupported_scnoise10_stage3ab_full
stage4_cytospace_scnoise10_stage3ab_full
```

对三个 control 场景，Stage3A 使用：

```text
filter_mode = none
cell_type_column = sc_meta
n_filtered = 0
```

对六个 single/double missing 场景，Stage3A 使用：

```text
filter_mode = plugin_unknown
filter_scope = unsupported_all
cell_type_column = plugin_type
truth_filter_enabled = False
```

Stage3B blank mask 在 CytoSPACE mapping 前应用，留白 spot 不参与 mapping capacity；映射结束后恢复完整坐标集合，并将留白 spot 写为全零组成行。两个噪声条件均使用 `n_subspots=800` 和 `mapping_cells_per_spot=5`。

### 3.3 其他映射方法及运行参数

外部方法由 `scripts/run_nine_scenario_method_benchmark.py` 调度：

| 方法 | 主要运行设置 |
|---|---|
| Tangram, all genes | `gene_mode=all`, `top_n_marker=50`, `num_epochs=200`, CPU |
| Tangram, marker genes | `gene_mode=marker`, `top_n_marker=50`, `num_epochs=200`, CPU |
| novoSpaRc | `max_genes=500`, `n_pcs=30` |
| SpaOTsc | `max_genes=500`, `n_pcs=30` |
| CellTrek | `max_genes=2000`, `n_pcs=30`, `ntree=500` |

这些方法的标准输出位置为：

```text
result/<sample>/stage4_mapping/<method>/spot_type_fraction.csv
```

CytoSPACE baseline 输出为：

```text
result/<sample>/stage4_cytospace_baseline/cytospace_output/fractional_abundances_by_spot.csv
```

## 4. 10% noise 的准确实现

噪声由 `scripts/generate_sc_noise_from_processed_sim.py` 生成，参数为：

```text
noise_fraction = 0.10
seed = 42
noise_type = sc_reference_entry_permutation
```

实现不是对所有表达值加入 10% 幅度的高斯噪声，也不是扰动 10% 的细胞或 10% 的基因。实际过程为：

1. 对 SC expression 矩阵的每个 entry 独立生成 Bernoulli mask，选择概率为 `0.10`；
2. 对每个基因列，所有被选中的 entry 均由该基因在随机 donor cell 中的表达值替换；
3. 因而保留每个基因的经验边际取值范围，但破坏约 10% 的 cell-by-gene 配对结构；
4. `sc_expression_normalized.csv`、`sc_expression_data.csv` 和存在时的 `sc_expression_counts.csv` 分别被扰动；
5. 当前 Stage3A/Stage4 route 明确使用 `sc_expression_normalized.csv`；
6. ST expression、ST coordinates 和 simulation truth 保持不变。

九个正式场景中 normalized SC matrix 的实际扰动比例为：

| 数据背景 | observed fraction |
|---|---:|
| Real BRCA（三个场景相同） | 0.0999122 |
| Human lung（三个场景相同） | 0.0999693 |
| Mouse brain（三个场景相同） | 0.0999333 |

`composite_scnoise10_stage3ab_full_input_audit.csv` 含 `36` 行输入审计：9 个 noisy SC normalized matrices 均与无噪声源文件不同，另外 `27` 个 ST expression、coordinates 和 truth 文件均与对应无噪声源文件哈希一致。

10% noise 中，CytoSPACE 和五种外部方法的既有 noisy-input mapping 被复用；完整 SVTuner route 重新运行。`composite_scnoise10_stage3ab_full_reused_method_audit.csv` 含 `54` 行，即 `9 scenarios x 6 reused methods`，全部通过 output existence、spot universe 和 noisy-input provenance 检查，`rerun_required=0`。

## 5. Abstention-aware whole-space recovery score

### 5.1 普通组成恢复

令第 (i) 个 spot 的预测和真值组成分别为 (p_{ik}) 和 (t_{ik})。评分前对每一行归一化，并对预测与真值细胞类型列取并集。普通 spot 的组成重叠为：

\[
O_i = \sum_k \min(p_{ik}, t_{ik}).
\]

当预测和真值均为归一化组成时，(O_i\in[0,1])。完全一致为 1，完全无重叠为 0。

评分使用 simulation truth 中的完整 spot universe。预测文件中缺失的 truth spot 会在 reindex 后成为全零行，不会从评价集合中删除。因此该指标不是 predicted-supported-region score。

### 5.2 留白的判定和计分

对每个场景，从 `sim_info.json` 读取 `sc_reference_drop_types`。若该类型是 spot 真值组成中占比最高的类型，则该 spot 被定义为 truth-unsupported spot：

\[
U_i = 1\{\arg\max_k t_{ik} \in K_{drop}\}.
\]

SVTuner 的最终 spot score 为：

\[
S_i =
\begin{cases}
1, & \text{SVTuner 留白且 } U_i=1;\\
0, & \text{SVTuner 留白且 } U_i=0;\\
O_i, & \text{SVTuner 未留白}.
\end{cases}
\]

场景总分为所有 truth spots 的算术均值：

\[
S = \frac{1}{N}\sum_i S_i.
\]

其他六种方法没有 abstention 输出，全部按 (O_i) 评分。该规则不会把正确留白误判为零，也不会无条件奖励留白：只有 target-dominant truth spot 中的留白得 1，非目标区域的错误留白明确得 0。

本评分中的 simulation truth 只用于最终评价和审计。Stage4 summaries 中 `truth_filter_enabled=False`，Stage3A 使用 `unsupported_all` 而非 `missing_only`，因此没有使用 oracle missing labels 执行过滤。

### 5.3 留白辅助指标

```text
blank precision = correct_abstention_spots / predicted_blank_spots
blank recall    = correct_abstention_spots / truth_unsupported_spots
```

这些指标以 target-dominant truth rule 为参照。它们衡量 engineered simulation truth 下的 abstention correctness，不应直接称为真实组织中的 biological accuracy。

## 6. 可视化设计

### 6.1 无噪声图

```text
visualizations/method_comparison/composite_no_noise/
composition_recovery_7mapping_methods_composite_no_noise_abstention_aware_boxplot.png
```

图题为：

```text
Spatial cell-type composition recovery, abstention-aware whole-space evaluation
```

### 6.2 10% noise 图

```text
visualizations/method_comparison/composite_scnoise10_stage3ab_full/
fig1d_scnoise10_stage3ab_full.png
```

图题为：

```text
Spatial cell-type composition recovery, 10% sc-reference noise
```

### 6.3 两图共同编码规则

1. x 轴顺序固定为 CytoSPACE、SVTuner、Tangram all genes、Tangram marker genes、novoSpaRc、SpaOTsc、CellTrek；
2. y 轴为 `Abstention-aware recovery score`，范围 `0-1.02`；
3. 每种方法有 9 个场景级观测点；灰色散点显示所有 9 个值；
4. 箱体显示四分位范围，中线为中位数，whisker 使用 seaborn boxplot 默认规则；
5. boxplot 自身隐藏离群点符号，但所有原始场景值仍通过 stripplot 显示；
6. 两图使用相同方法顺序和颜色，因此可以并排比较；
7. 图中没有 paired line，也没有按数据集或 missing condition 编码。逐场景配对结论必须来自 source CSV，不能只凭箱体位置判断。

## 7. 七种方法的总体结果

| 方法 | 0% mean | 0% median | 0% SD | 10% mean | 10% median | 10% SD | 10%-0% mean change |
|---|---:|---:|---:|---:|---:|---:|---:|
| CytoSPACE | 0.5154 | 0.5492 | 0.1078 | 0.5046 | 0.5171 | 0.1031 | -0.0109 |
| SVTuner | 0.7104 | 0.7119 | 0.1064 | 0.6553 | 0.6788 | 0.1632 | -0.0551 |
| Tangram, all genes | 0.4377 | 0.4639 | 0.1106 | 0.4513 | 0.4573 | 0.1168 | +0.0136 |
| Tangram, marker genes | 0.4628 | 0.4517 | 0.1204 | 0.4687 | 0.4554 | 0.1259 | +0.0059 |
| novoSpaRc | 0.4914 | 0.4766 | 0.1020 | 0.4889 | 0.4708 | 0.1032 | -0.0026 |
| SpaOTsc | 0.3950 | 0.4087 | 0.1178 | 0.3836 | 0.3852 | 0.1217 | -0.0115 |
| CellTrek | 0.5309 | 0.5403 | 0.1166 | 0.4352 | 0.4810 | 0.1196 | -0.0957 |

这些跨噪声变化是同一场景集合上的描述性变化，但 0% 和 10% SVTuner 不只是输入噪声不同：10% route 的 Stage3A 在若干场景发生漏检。因此 `SVTuner mean change=-0.0551` 不能被解释为纯粹的 mapping solver noise sensitivity。

## 8. 无噪声结果

### 8.1 SVTuner 与 CytoSPACE

无噪声时：

```text
SVTuner mean              = 0.710372
CytoSPACE mean            = 0.515413
absolute mean improvement = 0.194960
relative improvement      = 37.83%
paired wins               = 9/9
```

SVTuner 的平均 raw composition overlap 为 `0.514982`，平均 abstention credit 为 `0.195391`，两者之和对应平均总分 `0.710372`。这说明总体优势同时来自 supported-region composition 和正确 Stage3B abstention，不能将全部差异归因于 Stage3A。

### 8.2 与每种方法的配对比较

| 比较方法 | SVTuner mean difference | wins | losses |
|---|---:|---:|---:|
| CytoSPACE | +0.1950 | 9 | 0 |
| Tangram, all genes | +0.2727 | 9 | 0 |
| Tangram, marker genes | +0.2476 | 9 | 0 |
| novoSpaRc | +0.2189 | 9 | 0 |
| SpaOTsc | +0.3153 | 9 | 0 |
| CellTrek | +0.1795 | 8 | 1 |

因此可以写“SVTuner 在全部 9 个场景中高于 CytoSPACE”，但不能写“SVTuner 在每个场景都高于所有六种方法”，因为相对 CellTrek 为 `8/9`。

### 8.3 逐场景 SVTuner 与 CytoSPACE

| 数据背景 | Stage3A 条件 | CytoSPACE | SVTuner | difference |
|---|---|---:|---:|---:|
| Real BRCA | control | 0.6847 | 0.7845 | +0.0998 |
| Real BRCA | single missing | 0.5492 | 0.7119 | +0.1627 |
| Real BRCA | double missing | 0.3771 | 0.5735 | +0.1964 |
| Human lung | control | 0.5935 | 0.6728 | +0.0793 |
| Human lung | single missing | 0.4902 | 0.6305 | +0.1403 |
| Human lung | double missing | 0.3673 | 0.5666 | +0.1993 |
| Mouse brain | control | 0.5787 | 0.8215 | +0.2427 |
| Mouse brain | single missing | 0.5727 | 0.8562 | +0.2835 |
| Mouse brain | double missing | 0.4254 | 0.7759 | +0.3505 |

### 8.4 无噪声 Stage3A 过滤审计

从对应 Stage4 summaries 读取：

| 数据背景 | control n_filtered | single missing n_filtered | double missing n_filtered |
|---|---:|---:|---:|
| Real BRCA | 0 | 500 | 991 |
| Human lung | 0 | 2,779 | 5,870 |
| Mouse brain | 0 | 1,000 | 2,000 |

六个 missing 场景均使用 `plugin_unknown`，三个 control 均使用 `none`，所有场景 `truth_filter_enabled=False`。这些数量与工程化 missing cell counts 一致，支持无噪声 Stage3A route 已实际执行，而不是仅运行 Stage3B blanking。

### 8.5 无噪声留白审计

九个场景合并后：

```text
predicted blank spots = 10,445
correct abstentions   = 10,389
incorrect abstentions = 56
truth-unsupported     = 10,890
blank precision       = 0.9946
blank recall          = 0.9540
```

逐场景 precision/recall：

| 数据背景 | 条件 | precision | recall |
|---|---|---:|---:|
| Real BRCA | control | 1.0000 | 0.9948 |
| Real BRCA | single | 1.0000 | 0.9948 |
| Real BRCA | double | 1.0000 | 0.9948 |
| Human lung | control | 0.9744 | 0.9896 |
| Human lung | single | 0.9897 | 0.9723 |
| Human lung | double | 0.9873 | 0.9798 |
| Mouse brain | control | 0.9941 | 0.9411 |
| Mouse brain | single | 0.9960 | 0.9544 |
| Mouse brain | double | 0.9963 | 0.9390 |

## 9. 10% noise 完整 Stage3A+Stage3B 结果

### 9.1 SVTuner 与 CytoSPACE

10% noise 时：

```text
SVTuner mean              = 0.655260
CytoSPACE mean            = 0.504556
absolute mean improvement = 0.150704
relative improvement      = 29.87%
paired wins               = 8/9
```

SVTuner 的平均 raw composition overlap 为 `0.487245`，平均 abstention credit 为 `0.168015`。总体结果仍高于 CytoSPACE，但不是所有场景均改善。

### 9.2 与每种方法的配对比较

| 比较方法 | SVTuner mean difference | wins | losses |
|---|---:|---:|---:|
| CytoSPACE | +0.1507 | 8 | 1 |
| Tangram, all genes | +0.2040 | 8 | 1 |
| Tangram, marker genes | +0.1866 | 8 | 1 |
| novoSpaRc | +0.1664 | 8 | 1 |
| SpaOTsc | +0.2717 | 9 | 0 |
| CellTrek | +0.2201 | 6 | 3 |

### 9.3 逐场景 SVTuner 与 CytoSPACE

| 数据背景 | Stage3A 条件 | CytoSPACE | SVTuner | difference |
|---|---|---:|---:|---:|
| Real BRCA | control | 0.6501 | 0.7687 | +0.1186 |
| Real BRCA | single missing | 0.5171 | 0.6788 | +0.1616 |
| Real BRCA | double missing | 0.3636 | 0.5186 | +0.1550 |
| Human lung | control | 0.5874 | 0.5894 | +0.0020 |
| Human lung | single missing | 0.4831 | 0.4802 | -0.0030 |
| Human lung | double missing | 0.3634 | 0.4057 | +0.0424 |
| Mouse brain | control | 0.5799 | 0.8267 | +0.2468 |
| Mouse brain | single missing | 0.5719 | 0.8570 | +0.2851 |
| Mouse brain | double missing | 0.4245 | 0.7723 | +0.3478 |

唯一相对 CytoSPACE 的 loss 是 human-lung single-missing 场景，差值为 `-0.002959`。因此正文不能写“10% noise 下全部九个场景均改善”。

### 9.4 数据背景分层均值

| 数据背景 | 0% CytoSPACE | 0% SVTuner | 10% CytoSPACE | 10% SVTuner |
|---|---:|---:|---:|---:|
| Real BRCA | 0.5370 | 0.6900 | 0.5103 | 0.6554 |
| Human lung | 0.4836 | 0.6233 | 0.4780 | 0.4917 |
| Mouse brain | 0.5256 | 0.8179 | 0.5254 | 0.8187 |

主要噪声敏感性集中在 human lung，而 mouse brain 的三个 SVTuner 分数基本稳定。该异质性与下面的 Stage3A 漏检和 Stage3B lung recall 下降一致。

### 9.5 10% noise Stage3A 审计

| 数据背景 | 条件 | expected types（仅审计） | detected types | n_filtered | status |
|---|---|---|---|---:|---|
| Real BRCA | control | 无 | 无 | 0 | not applicable |
| Real BRCA | single | Epithelial cells | Epithelial cells | 500 | detected |
| Real BRCA | double | Epithelial cells + PCs | Epithelial cells | 500 | partial detection |
| Human lung | control | 无 | 无 | 0 | not applicable |
| Human lung | single | AT2 | 无 | 0 | detection failure |
| Human lung | double | AT2 + Fibroblast | AT2 | 2,779 | partial detection |
| Mouse brain | control | 无 | 无 | 0 | not applicable |
| Mouse brain | single | Micro | Micro | 1,000 | detected |
| Mouse brain | double | Micro + Oligo_2 | Micro + Oligo_2 | 2,000 | detected |

所有 missing 场景均按配置请求 `plugin_unknown`，但 human-lung single 场景实际 `n_filtered=0`。因此 guardrail 中：

```text
missing_n_filtered_positive_quality_check = false
all_expected_stage3a_types_detected_quality_check = false
decision = PARTIAL
```

`PARTIAL` 表示核心执行、输入 provenance、truth-filter 禁用、9x7 结果完整性和 Stage3B integration 均通过，但 Stage3A detection quality 未全部通过。它不能被改写成完整的 noise-robustness PASS。

### 9.6 10% noise Stage3B 审计

| 数据背景 | 条件 | blank spots | dominant precision | dominant recall | target-mass recall | zero-target FPR |
|---|---|---:|---:|---:|---:|---:|
| Real BRCA | control | 382 | 1.0000 | 0.9948 | 0.9948 | 0.0000 |
| Real BRCA | single | 382 | 1.0000 | 0.9948 | 0.9948 | 0.0000 |
| Real BRCA | double | 382 | 1.0000 | 0.9948 | 0.9948 | 0.0000 |
| Human lung | control | 205 | 0.6829 | 0.3733 | 0.3662 | 0.02334 |
| Human lung | single | 22 | 0.8636 | 0.0491 | 0.0491 | 0.00110 |
| Human lung | double | 8 | 0.6250 | 0.0129 | 0.0131 | 0.00110 |
| Mouse brain | control | 2,914 | 0.9749 | 0.9958 | 0.9821 | 0.00131 |
| Mouse brain | single | 2,902 | 0.9793 | 0.9961 | 0.9804 | 0.00164 |
| Mouse brain | double | 3,042 | 0.9372 | 0.9993 | 0.9944 | 0.00098 |

九场景合并后的 target-dominant blank 指标为：

```text
predicted blank spots = 10,239
correct abstentions   = 9,844
incorrect abstentions = 395
truth-unsupported     = 10,890
blank precision       = 0.9614
blank recall          = 0.9039
```

聚合 precision/recall 较高，但受 BRCA 和 mouse brain 的大量高质量留白主导。Human lung single/double 场景的 dominant recall 分别仅为 `0.0491` 和 `0.0129`，必须单独披露，不能只报告合并值。

## 10. 统计与解释边界

### 10.1 可以支持的表述

1. 在无噪声联合模拟中，SVTuner 的平均 abstention-aware whole-space recovery 高于 CytoSPACE，且 9/9 场景配对提高；
2. 在 10% entry-wise within-gene replacement noise 下，完整 SVTuner route 的平均分仍高于 CytoSPACE，配对结果为 8/9；
3. 正确 Stage3B 留白不会被计为零，错误留白仍被惩罚为零；
4. 无噪声 Stage3A 过滤数量与工程化 missing counts 一致；
5. 10% noise 下性能具有明显数据背景异质性，Stage3A 和 Stage3B 在 human lung 场景中较脆弱。

### 10.2 不应使用的表述

1. “10% noise 下所有场景均改善”；实际为 8/9 对 CytoSPACE；
2. “SVTuner 在所有场景都优于所有方法”；无噪声对 CellTrek 为 8/9，10% noise 为 6/9；
3. “10% noise 完整审计 PASS”；当前正式 decision 为 `PARTIAL`；
4. “所有 Stage3A missing types 在 10% noise 下均被识别”；BRCA double 和 lung double 为 partial，lung single 为 failure；
5. “总体 blank recall=0.9039 说明所有组织都稳定”；lung 分层结果不支持这一结论；
6. “箱线图显示统计显著性”；当前没有 p-value、置信区间或预先规定的 inferential test；
7. “9 个点是生物学重复”；它们是三个模拟背景下的九个工程化场景；
8. “结果证明新的 mapping engine 优于其他方法”；SVTuner 在这里是 filtering + abstention + CytoSPACE route。

## 11. 可复现命令

### 11.1 运行或续跑五种外部方法

```powershell
# No noise
python scripts/run_nine_scenario_method_benchmark.py `
  --project_root . `
  --scenario_preset composite

# 10% SC-reference noise
python scripts/run_nine_scenario_method_benchmark.py `
  --project_root . `
  --scenario_preset composite `
  --sample_suffix _scnoise10
```

### 11.2 重建无噪声汇总和图

```powershell
python scripts/plot_method_comparison_composition_recovery.py `
  --project_root . `
  --scenario_preset composite `
  --route2_stage4_dir stage4_cytospace_stage3b_blank `
  --reward_correct_stage3b_abstention `
  --abstention_truth_rule target_dominant `
  --out_dir visualizations/method_comparison/composite_no_noise `
  --output_prefix composition_recovery_7mapping_methods_composite_no_noise_abstention_aware `
  --title "Spatial cell-type composition recovery, abstention-aware whole-space evaluation" `
  --ylabel "Abstention-aware recovery score"
```

### 11.3 审计并运行 10% noise 完整路线

```powershell
# 仅检查输入和可复用 mapping provenance
python scripts/run_composite_scnoise10_stage3ab_full.py `
  --project_root .

# 运行 Stage3A、Stage3B、Stage4，随后汇总并绘图
python scripts/run_composite_scnoise10_stage3ab_full.py `
  --project_root . `
  --execute `
  --resume `
  --n_subspots 800 `
  --mapping_cells_per_spot 5
```

若质量审计仍为 `PARTIAL`，runner 会在写出完整结果后返回非零状态码 `2`。这表示质量门未完全通过，不表示结果文件没有生成。

## 12. 关键证据文件

### 12.1 无噪声

```text
visualizations/method_comparison/composite_no_noise/
  composition_recovery_7mapping_methods_composite_no_noise_abstention_aware_boxplot.png
  composition_recovery_7mapping_methods_composite_no_noise_abstention_aware_boxplot.pdf
  composition_recovery_7mapping_methods_composite_no_noise_abstention_aware_summary.csv
  composition_recovery_7mapping_methods_composite_no_noise_abstention_aware_scenario.csv
  composition_recovery_7mapping_methods_composite_no_noise_abstention_aware_cell_type.csv
```

`scenario.csv` 是主图每个灰点和所有配对比较的直接数据源；`cell_type.csv` 是更细粒度的 per-cell-type recovery 辅助表。

### 12.2 10% noise 完整路线

```text
visualizations/method_comparison/composite_scnoise10_stage3ab_full/
  fig1d_scnoise10_stage3ab_full.png
  fig1d_scnoise10_stage3ab_full.pdf
  composite_scnoise10_stage3ab_full_source_values.csv
  composite_scnoise10_stage3ab_full_method_summary.csv
  composite_scnoise10_stage3ab_full_paired_comparison.csv
  composite_scnoise10_stage3ab_full_svtuner_method_comparisons.csv
  composite_scnoise10_stage3ab_full_blank_summary.csv
  composite_scnoise10_stage3ab_full_stage3a_audit.csv
  composite_scnoise10_stage3ab_full_stage3b_audit.csv
  composite_scnoise10_stage3ab_full_input_audit.csv
  composite_scnoise10_stage3ab_full_reused_method_audit.csv
  composite_scnoise10_stage3ab_full_scenario_manifest.csv
  composite_scnoise10_stage3ab_full_guardrails.json
  composite_scnoise10_stage3ab_full_summary.json
```

当前 `guardrails.json` 是 2026-07-11 正式运行时写出的历史审计快照，其中 `old_output_hashes_*` 字段记录了当时对已退休 Stage3B-only 结果的只读保护检查。旧结果目录后来已删除；这些历史字段不参与当前 source values、method summary、Stage3A audit 或 Stage3B audit 的数值计算。现行 runner 已解除对退休目录的运行依赖。

## 13. 建议论文结果段落

### 13.1 中文结果表述

> 在三个模拟数据背景构成的九个联合 reference-mismatch 场景中，我们采用全空间 abstention-aware composition overlap 同时评价支持区域的组成恢复与不支持区域的正确留白。无额外 SC-reference 噪声时，SVTuner 的平均恢复分数为 0.7104，高于 CytoSPACE 的 0.5154，绝对提高 0.1950，并在全部九个配对场景中取得更高分。加入 10% entry-wise within-gene replacement noise 后，完整 Stage3A+Stage3B SVTuner 的平均分数为 0.6553，而 CytoSPACE 为 0.5046，绝对提高 0.1507，配对结果为 8/9。该噪声效应具有明显数据背景差异：mouse-brain 场景基本稳定，而 human-lung 场景中 Stage3A 出现漏检且 Stage3B 留白召回下降。因此，10% noise 结果支持总体优势仍然存在，但不支持跨场景一致鲁棒性的结论。

### 13.2 英文结果表述

> Across nine joint reference-mismatch scenarios spanning three simulation backgrounds, we evaluated supported-region composition recovery and unsupported-region abstention using a whole-space, target-dominant abstention-aware overlap score. Without additional single-cell-reference noise, SVTuner achieved a mean score of 0.7104 compared with 0.5154 for CytoSPACE, an absolute improvement of 0.1950, and improved all nine paired scenarios. Under 10% entry-wise within-gene replacement noise, the full Stage3A+Stage3B SVTuner route achieved a mean score of 0.6553 compared with 0.5046 for CytoSPACE, with improvements in eight of nine scenarios. Robustness was dataset-dependent: the mouse-brain scenarios remained stable, whereas Stage3A detection and Stage3B abstention recall deteriorated in the human-lung scenarios.

## 14. 建议图例

> **Spatial cell-type composition recovery under reference mismatch and SC-reference noise.** Boxplots summarize nine engineered simulation scenarios across Real BRCA, human lung and mouse brain backgrounds. Each point represents one scenario. Predictions and simulation truth were normalized within spot and evaluated over the complete truth spot universe. For methods without abstention, the score is the mean spot-wise composition overlap. For SVTuner, a blanked spot receives a score of one only when the cell type removed from the SC reference is dominant in simulation truth, while an incorrect blank receives zero; all nonblank spots retain their composition-overlap score. The no-noise panel uses the Stage3A+Stage3B SVTuner route without added expression perturbation. The 10% noise panel uses independently selected SC expression entries replaced by values from donor cells for the same gene (seed 42), while ST expression, coordinates and simulation truth remain unchanged. Boxes show the interquartile range and median; whiskers follow the default 1.5-IQR convention, and all scenario-level observations are overlaid.
