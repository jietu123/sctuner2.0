# 论文材料准备B：Stage3B 与 ST-only unsupported-region

> 状态：已按项目现存脚本、CSV、JSON 和论文候选图重新核对；2026-07-11 补充无噪声与 10% scRNA noise 的七方法 abstention-aware whole-space benchmark。本文所有结果值均来自项目内可追溯文件，不使用已删除旧结果。相对路径均以 `sctuner2.0/` 为根。

## 0. 结论摘要

Stage3B 处理的不是 Stage3A 的“SC reference 中存在、ST 中缺失”的问题，而是相反方向的 reference mismatch：**ST 中存在真实表达结构，但 scRNA-seq reference 缺少相应类型或状态**。如果仍进行强制映射，缺失类型区域会被解释成 reference 中剩余的相似类型；Stage3B 的任务是在映射前识别这类 reference-unsupported spots，并允许后端在这些位置 abstain/blank。

项目现有证据支持以下结论：

1. 在 9 个 Stage3A+Stage3B 联合模拟场景中，Stage3A 和 Stage3B 检查均为 `9/9` 通过；Stage3B dominant-region recall 为 `0.9948`，dominant-region precision 为 `0.9736`，target mass recall 为 `0.9908`，zero-target false-positive rate 为 `0.000593`。
2. 在无噪声联合 benchmark 中，SVTuner 的 abstention-aware whole-space score 为 `0.7104`，高于 CytoSPACE 的 `0.5154`，并在 `9/9` 场景中取胜；正确留白 precision 为 `0.9946`、recall 为 `0.9540`。
3. 在 10% scRNA noise 下，SVTuner mean score 为 `0.6220`，高于 CytoSPACE 的 `0.5046`，但只在 `7/9` 场景中取胜；正确留白 precision/recall 降至 `0.9614/0.9039`，说明收益保留但噪声鲁棒性并非无损。
4. 模拟结果不能写成“真实数据中可完美留白”。在 6 个真实 reference-dropout 场景、top15 marker core 定义下，blank spots 的平均 target-associated fraction 为 `0.9379`，平均 core recall 仅为 `0.5778`，说明当前方法偏向高精度、保守覆盖。
5. 在真实 mouse brain `ST8059051 / Thalamic excitatory` 主案例中，Stage3B 对 top15 marker region 的 recall 为 `0.9586`、precision 为 `0.7199`；abstention-aware accuracy 从 CytoSPACE 的 `0.3295` 提升到 `0.5721`。
6. 当 `Thalamic excitatory` 从 reference 中移除后，CytoSPACE 将 362 个目标 spots 中的 `52.76%` 强制分配为 inhibitory neuron、`24.31%` 分配为 astrocyte、`18.78%` 分配为 excitatory neuron，说明错误不是自然留白，而是系统性替代。
7. 在 BRCA HER2 FFPE plasma-cell dropout 下，Stage3B 将 reference-relative coupling absolute error 从 `1.0000` 降到 `0.3501`，material spurious pairs 从 `10` 降到 `0`；代表性 `KITLG->KIT` local Moran 均值下降 `84.28%`，固定阈值 hotspots 从 `110` 降到 `84`。

论文中最稳健的定位是：

> Stage3B is a reference-adequacy diagnostic and abstention layer that detects spatial expression regions insufficiently explained by the available single-cell reference, thereby reducing forced surrogate assignments and their downstream spatial artifacts.

## 1. Stage3A 与 Stage3B 的边界

| 模块 | mismatch 方向 | 诊断对象 | 主要动作 | 主要风险 |
|---|---|---|---|---|
| Stage3A | SC-only | reference 中不被 ST 支持的类型 | 从 mapping pool 过滤或重标记 SC cells | 误删相似但真实存在的类型 |
| Stage3B | ST-only | ST 中不被 reference 解释的 spots/regions | 对相应空间位置 abstain/blank | 漏掉边界或弱信号区域，或在 domain shift 下产生额外 blank |

Stage3B 的核心实现：

```text
src/stages/stage3b_st_unsupported.py
```

配套脚本：

```text
scripts/generate_st_only_reference_dropout_from_sim.py
scripts/visualize_stage3b_st_only_triptych.py
scripts/evaluate_stage3b_simulation.py
scripts/prepare_stage3b_realdata_reference_dropout_scenarios.py
scripts/plot_stage3b_reference_dropout_spatial_stack.py
scripts/plot_stage3b_reference_dropout_panel_a_blank_composition.py
scripts/plot_stage3b_threshold_robustness.py
```

## 2. Stage3B 算法逻辑

### 2.1 输入和输出

Stage3B 使用：

- SC expression 与 cell-type labels；
- ST expression；
- SC/ST 共同基因；
- ST 坐标，仅用于空间连通和区域检验；
- 预先设定的统计参数，例如默认 FDR `0.05` 和空间置换次数 `200`。

主要输出包括 spot-level unsupported score、candidate/region flags、blank fraction、区域表和 summary。下游 Stage4 可以据此阻止在 unsupported spots 上继续强制分配。

### 2.2 Reference reconstruction 与异常特征

1. 按 SC label 构建 cell-type expression profiles。
2. 对每个 ST spot 用非负最小二乘拟合 reference type profiles。
3. 从 observed ST 与 reference reconstruction 的差异中计算四类特征：relative reconstruction error、cosine deficit、positive residual fraction 和 residual concentration。
4. 从真实 SC cells 生成 supported pseudo-ST calibration，使异常定义相对于“reference 能够解释的混合”建立，而不是相对于人工目标类型建立。
5. 对异常方向做 one-sided Stouffer 聚合和 Benjamini-Hochberg 校正。

### 2.3 空间区域与 residual-program 分支

逐 spot 异常只作为候选。Stage3B 进一步基于坐标构建空间边，并对连通区域进行置换检验。只有通过区域显著性检验的结构才进入 unsupported region。

对于整体 profile reconstruction 尚可、但存在局部正残差程序的情况，代码还包含 residual-program branch：

- 对 novelty-weighted positive residual 做低秩分解；
- 检验 residual components 的 signed tails；
- 要求区域与 whole-profile anomaly 有足够重叠；
- 要求 reference-orthogonal score 达到门槛；
- 再经过空间区域 gate。

该分支用于处理“部分未被 reference 表达程序解释”的情况，而不是依赖某个已知 marker 白名单。

### 2.4 白名单与评估真值的严格区分

Stage3B 主算法不接收“待留白目标类型”作为检测白名单。模拟场景中的 target label、真实数据中的 marker-defined region 用于**事后评估和画图**，不用于 Stage3B 的 score 或阈值选择。

需要同时承认两个限制：

1. 模拟数据的生成过程知道被删除类型，因而真值边界比真实数据清晰；高分不等于真实数据可完美恢复。
2. 真实数据的 top10/top15/top20/top25 marker core 是评估口径，不是单细胞级 ground truth；因此应写作 marker-defined target-associated region，而不是 definitive cell-type truth。

## 3. Stage3A+Stage3B 联合模拟

### 3.1 设计

每个数据集包含三行：Stage3A 无缺失、单缺失、双缺失；同时从 SC reference 删除一个在 ST 中仍存在的类型，形成 Stage3B 目标。

| 数据组 | Stage3B 目标 | Stage3A 单缺失 | Stage3A 双缺失 |
|---|---|---|---|
| Real BRCA | Endothelial cells | Epithelial cells | Epithelial cells + PCs |
| Human lung 5-location | B cell | AT2 | AT2 + Fibroblast |
| Mouse brain refined | Ext_L56 | Micro | Micro + Oligo_2 |

论文总览图：

```text
visualizations/simulations/simulation_stage3ab_joint_triptych_overview_stack_3datasets.svg
visualizations/simulations/simulation_stage3ab_joint_triptych_overview_stack_3datasets.png
```

审计表：

```text
visualizations/simulations/composite_stage3a_stage3b_recheck.csv
```

### 3.2 九场景汇总

| 指标 | 9 场景均值/计数 |
|---|---:|
| Stage3A expected vs detected 完全一致 | 9/9 |
| Stage3B scenario check 通过 | 9/9 |
| dominant-region recall | 0.9948 |
| dominant-region precision | 0.9736 |
| target-positive recall | 0.8060 |
| target mass recall | 0.9908 |
| zero-target false-positive rate | 0.000593 |

分数据集结果：

| 数据组 | dominant recall | dominant precision | positive recall | zero-target FP rate |
|---|---:|---:|---:|---:|
| Real BRCA | 0.9948 | 1.0000 | 0.9948 | 0.000000 |
| Human lung | 0.9924 | 0.9555 | 0.8451 | 0.000593 |
| Mouse brain refined | 0.9971 | 0.9654 | 0.5780 | 0.001184 |

Mouse brain 的 dominant-region 指标很高，但 positive recall 只有 `0.5780`。这解释了为什么空间图中主要核心区域可以被正确留白，而混合 spots、边界 spots 和低占比目标 spots 仍然保留。论文中应强调“dominant unsupported region recovery”，不能将其扩展为所有 target-positive spots 的完整恢复。

### 3.3 Abstention-aware 全空间多方法 comparison

新版评分始终把全部 truth spots 保留在分母中。非留白 spot 使用逐行归一化组成的 overlap score `sum_k min(pred_ik, truth_ik)`；SVTuner 留白仅在独立 simulation truth 表明 SC-reference-dropped 类型为该 spot 的 dominant type 时计 `1`，其余错误留白计 `0`。因此，正确留白不会被误判为零分，错误留白也不能通过删除评估 spot 而逃避惩罚。

无噪声结果：

| 方法 | n | mean | median | SD |
|---|---:|---:|---:|---:|
| CytoSPACE | 9 | 0.5154 | 0.5492 | 0.1078 |
| SVTuner | 9 | 0.7104 | 0.7119 | 0.1064 |
| Tangram, all genes | 9 | 0.4377 | 0.4639 | 0.1106 |
| Tangram, marker genes | 9 | 0.4628 | 0.4517 | 0.1204 |
| novoSpaRc | 9 | 0.4914 | 0.4766 | 0.1020 |
| SpaOTsc | 9 | 0.3950 | 0.4087 | 0.1178 |
| CellTrek | 9 | 0.5309 | 0.5403 | 0.1166 |

SVTuner 相对 CytoSPACE 的 mean absolute improvement 为 `+0.1950`，relative improvement 为 `+37.83%`，逐场景为 `9/9` 胜；相对 CellTrek 为 `8/9` 胜，相对其余方法均为 `9/9` 胜。共预测 `10,445` 个 blank spots，其中 `10,389` 个正确、`56` 个错误；相对于 `10,890` 个 truth unsupported spots，abstention precision 为 `0.9946`、recall 为 `0.9540`。

10% scRNA noise 结果：

| 方法 | n | mean | median | SD |
|---|---:|---:|---:|---:|
| CytoSPACE | 9 | 0.5046 | 0.5171 | 0.1031 |
| SVTuner | 9 | 0.6220 | 0.6346 | 0.1576 |
| Tangram, all genes | 9 | 0.4513 | 0.4573 | 0.1168 |
| Tangram, marker genes | 9 | 0.4687 | 0.4554 | 0.1259 |
| novoSpaRc | 9 | 0.4889 | 0.4708 | 0.1032 |
| SpaOTsc | 9 | 0.3836 | 0.3852 | 0.1217 |
| CellTrek | 9 | 0.4352 | 0.4810 | 0.1196 |

10% noise 下，SVTuner 相对 CytoSPACE 的 mean absolute improvement 为 `+0.1174`，relative improvement 为 `+23.27%`，逐场景为 `7/9` 胜。共预测 `10,239` 个 blank spots，其中 `9,844` 个正确、`395` 个错误；abstention precision 为 `0.9614`、recall 为 `0.9039`。两项 human-lung missing 场景出现轻微负 delta：single missing `-0.0030`、double missing `-0.0035`。

按数据组比较 CytoSPACE -> SVTuner，无噪声为 BRCA `0.5370 -> 0.6900`、lung `0.4837 -> 0.6233`、mouse brain `0.5256 -> 0.8179`；10% noise 为 BRCA `0.5103 -> 0.6327`、lung `0.4780 -> 0.4765`、mouse brain `0.5254 -> 0.7568`。与无噪声相比，SVTuner mean 降低 `0.0884`，abstention precision/recall 分别降低 `0.0332/0.0500`，表明 reference noise 主要削弱 Stage3B abstention 的稳定性。

旧版 naive whole-space 口径得到的 SVTuner mean `0.5166` 会把正确留白当作 `0`；旧版 predicted-supported-region 口径得到的 `0.6501` 则使用方法依赖的评估子集。两套旧文件已删除，以上数值只用于解释口径变更，不应进入论文结果表。新版 benchmark 支持 SVTuner 在预定义 selective-utility 任务上的总体优势，但仍不能外推为传统 composition recovery 或任意噪声条件下的普遍最优。

## 4. 真实主案例：ST8059051 Thalamic excitatory

主目录：

```text
visualizations/cell2location_stage3b_case/thalamic_top15/
```

### 4.1 Target region 定义

- ST spots：`2,409`
- target SC cells：`1,868`
- marker genes：`30`
- marker-region quantile：`0.85`，即 top15%
- marker-region spots：`362`，占全部 spots 的 `15.03%`
- marker-score threshold：`0.5541`

目标 broad label 为 `Thalamic excitatory`；reference dropout 对应 `Ext_Thal_1` 和 `Ext_Thal_2`。

### 4.2 Stage3B 空间检测

| 指标 | 数值 |
|---|---:|
| Stage3B blank spots | 482 |
| target-region overlap | 347 |
| target-region missed | 15 |
| blank outside marker region | 135 |
| precision | 0.7199 |
| recall | 0.9586 |
| target core spots | 302 |
| core hits | 301 |
| boundary spots | 60 |
| boundary hits | 46 |

命中 spots 的 median marker score 为 `1.0836`，漏检 spots 为 `0.5993`；命中 spots 的 median residual z 为 `4.6145`，漏检 spots 为 `1.4276`。因此误差主要集中在弱信号和边界，而不是核心区域。

### 4.3 Baseline forced assignment

362 个 target-region spots 在 CytoSPACE baseline 中被解释为：

| 被强制分配类型 | spots | 比例 |
|---|---:|---:|
| Inhibitory neuron | 191 | 52.76% |
| Astrocyte | 88 | 24.31% |
| Excitatory neuron | 68 | 18.78% |
| Oligodendrocyte/OPC | 13 | 3.59% |
| Neuroblast | 1 | 0.28% |
| Microglia | 1 | 0.28% |

该结果说明 reference 缺失会把目标区域转换成具有解释诱惑力的相似类型，而不是产生显式失败信号。

### 4.4 Abstention-aware accuracy

评估集包含 `1,290` spots：362 个 unsupported core spots，加上 928 个 supported-region spots。

| 方法 | unsupported correct | supported correct | accuracy |
|---|---:|---:|---:|
| CytoSPACE | 0 | 425 | 0.3295 |
| SVTuner | 347 | 391 | 0.5721 |

绝对提升为 `0.2426`。SVTuner 在 supported-region correct 上从 425 降至 391，因此结果不是“无代价改善”；收益来自正确 abstention 大于 supported-region 损失。

## 5. 六个真实 reference-dropout 场景

空间图：

```text
visualizations/stage3b_realdata_candidate_scan/spatial_9x2/stage3b_reference_dropout_spatial_stack_recommended_6x2.svg
```

定量图：

```text
visualizations/stage3b_realdata_candidate_scan/panel_a_blank_composition/stage3b_reference_dropout_panel_a_blank_region_composition.png
```

top15 结果：

| 场景 | blank spots | inside core | near core | outside associated | inside fraction | target-associated fraction |
|---|---:|---:|---:|---:|---:|---:|
| Mouse embryo Endoderm/Gut | 832 | 702 | 5 | 125 | 0.8438 | 0.8498 |
| Mouse embryo Erythroid | 812 | 686 | 21 | 105 | 0.8448 | 0.8707 |
| BRCA TNBC Plasma cells | 97 | 91 | 6 | 0 | 0.9381 | 1.0000 |
| CRC B cells | 298 | 234 | 43 | 21 | 0.7852 | 0.9295 |
| BRCA HER2 FFPE Plasma cells | 337 | 231 | 100 | 6 | 0.6855 | 0.9822 |
| BRCA HER2 FFPE Epithelial cells | 204 | 159 | 44 | 1 | 0.7794 | 0.9951 |

跨场景均值：inside-core fraction `0.8128`，near-core fraction `0.1251`，outside fraction `0.0621`，target-associated fraction `0.9379`，median target-associated fraction `0.9559`。

这组真实结果应写成“blank spots predominantly localized to marker-defined target-associated regions”，而不是“blank 完全正确”。

## 6. Marker-core 阈值鲁棒性

| 定义 | mean target-associated | median target-associated | mean inside | mean outside | mean core recall |
|---|---:|---:|---:|---:|---:|
| top10 | 0.8972 | 0.8961 | 0.6817 | 0.1028 | 0.7260 |
| top15 | 0.9379 | 0.9559 | 0.8128 | 0.0621 | 0.5778 |
| top20 | 0.9479 | 0.9635 | 0.8732 | 0.0521 | 0.4644 |
| top25 | 0.9567 | 0.9733 | 0.8977 | 0.0433 | 0.3824 |

随着 marker core 扩大，blank composition 看起来更“纯”，但 core recall 下降。原因是评估 denominator 增大，而 Stage3B blank set 相对保守。top15 是 purity 与 coverage 的折中，并非通过最大化某个结果指标后选择。

## 7. Reference-missing stress test

主图：

```text
visualizations/stage3b_reference_missing_stress/cell2location_reference_missing_multitarget_panels_bc.png
```

7 个 broad categories 均在同一 362-spot target-region 口径下分别从 reference 删除。forced-assignment heatmap 和 expression-similarity heatmap 共同说明 surrogate assignment 与类型相似性有关。例如：

- Excitatory neuron 缺失后，`79.01%` 被分配为 inhibitory neuron；两类 expression cosine similarity 为 `0.9688`。
- Inhibitory neuron 缺失后，`50.55%` 被分配为 excitatory neuron、`34.81%` 为 astrocyte。
- Oligodendrocyte/OPC 缺失后，`64.64%` 被分配为 astrocyte。
- Neuroblast 缺失后，`54.70%` 被分配为 excitatory neuron；两类 similarity 为 `0.8666`。

这些例子支持“错误具有表达相似性驱动的结构”，但当前图没有给出跨矩阵单元的正式相关检验，因此不应写成已证明的 causal correlation。

## 8. False spatial niche / coupling

场景：BRCA HER2 FFPE，Plasma cells reference dropout。

主图：

```text
visualizations/stage3b_false_spatial_niche/candidate_validation/brca_her2_plasma/cytospace_fig2d_tme_brca_her2_ffpe_sc_missing_plasma_cells_ncem_style.png
```

关键值：

| 指标 | 数值 |
|---|---:|
| total spots | 2,518 |
| marker target spots | 378 |
| Stage3B blank spots | 337 |
| target/blank overlap | 231 |
| precision | 0.6855 |
| recall | 0.6111 |
| CytoSPACE absolute coupling error | 1.0000 |
| SVTuner absolute coupling error | 0.3501 |
| material-effect threshold | 0.02530 |
| CytoSPACE material spurious pairs | 10 |
| SVTuner material spurious pairs | 0 |

`uses_target_markers_in_stage3b=False`，因此 target markers 仅用于评估 target region，没有进入 Stage3B 检测。

建议表述：

> Reference dropout induced broad reference-relative neighborhood coupling changes, whereas Stage3B abstention reduced the aggregate coupling error and eliminated material spurious pairs under the predefined effect threshold.

这里的 coupling 是 project-specific reference-relative diagnostic，不应等同于经过独立生物学验证的真实细胞互作。

## 9. Communication validation

主图：

```text
visualizations/stage3b_communication_validation/brca_her2_ffpe_plasma/strict_reference_style/stage3b_communication_strict_reference_style.png
```

### 9.1 Pathway-level false-pair summary

| Pathway | selected pairs | false fraction | Fisher p |
|---|---:|---:|---:|
| NFkB | 8 | 0.625 | 0.0900 |
| TGFb | 3 | 0.333 | 0.7850 |
| TNFa | 3 | 0.333 | 0.7850 |
| Hypoxia | 1 | 1.000 | 0.3810 |

这些 pathway p-values 未达到常用显著性门槛，因此 Panel A 只能作为描述性汇总，不应用于显著通路发现声明。

### 9.2 KITLG->KIT local hotspot

| 指标 | CytoSPACE dropout | SVTuner | 变化 |
|---|---:|---:|---:|
| top10% spots | 507 | 474 | -33 |
| mean local Moran | 0.02765 | 0.00435 | -84.28% |
| fixed-threshold hotspots | 110 | 84 | -23.64% |

Panel C 的 ligand-target matrix 和 ligand scores 用于提供表达响应背景。最高 ligand scores 包括 `ADM=1.0000`、`VWF=0.8752`、`TGFB1=0.6916`、`VEGFA=0.5785`，但这些是当前验证框架中的归一化 scores，不应作为新 ligand discovery 独立报告。

## 10. 论文图组织

正文建议分成三层：

1. **算法和联合模拟**：Stage3A+Stage3B joint overview，加 9 场景检测指标以及无噪声/10% noise abstention-aware 全空间比较。
2. **真实 reference-dropout 主证据**：6x2 spatial stack + blank composition；Thalamic top15 深入案例作为定量放大。
3. **下游后果**：forced-assignment/similarity heatmap 和 false niche；communication 放补充材料。

推荐正文图候选：

```text
visualizations/simulations/simulation_stage3ab_joint_triptych_overview_stack_3datasets.svg
visualizations/method_comparison/composite_no_noise/composition_recovery_7mapping_methods_composite_no_noise_abstention_aware_boxplot.png
visualizations/method_comparison/composite_scnoise10/composition_recovery_7mapping_methods_composite_scnoise10_abstention_aware_boxplot.png
visualizations/stage3b_realdata_candidate_scan/spatial_9x2/stage3b_reference_dropout_spatial_stack_recommended_6x2.svg
visualizations/stage3b_realdata_candidate_scan/panel_a_blank_composition/stage3b_reference_dropout_panel_a_blank_region_composition.png
visualizations/cell2location_stage3b_case/thalamic_top15/cell2location_ST8059051_thalamic_excitatory_abstention_aware_region_accuracy.png
visualizations/stage3b_false_spatial_niche/candidate_validation/brca_her2_plasma/cytospace_fig2d_tme_brca_her2_ffpe_sc_missing_plasma_cells_ncem_style.png
```

补充材料候选：top20 sensitivity、threshold robustness、multitarget forced-assignment heatmaps 和 communication panels。若正文只保留无噪声 comparison，则 10% noise 箱线图和逐场景审计表放补充材料。

## 11. 写作边界

可以写：

- Stage3B detects spatial regions insufficiently explained by the available reference.
- Blank calls were predominantly associated with independently defined marker-rich regions across six reference-dropout settings.
- Stage3B reduced forced surrogate assignment and reference-relative downstream artifacts.
- Detection was conservative in mixed and boundary spots.
- The mean abstention-aware whole-space score exceeded CytoSPACE under both no-noise and 10% reference-noise conditions, with 9/9 and 7/9 paired scenario wins, respectively.

不能写：

- Stage3B perfectly recovers every missing type or every target-positive spot.
- Real-data blank regions are definitive single-cell ground truth.
- Stage3B used no truth information anywhere in the experiment；正确说法是算法未使用 target labels，但模拟生成和事后评价使用了真值。
- Communication analysis discovered new biological signaling pathways.
- Composite composition recovery proves broad superiority over all mapping methods.
- SVTuner improves every scenario under 10% scRNA-seq noise.

## 12. 可直接用于 Results 的英文段落

> In nine joint mismatch simulations, Stage3A recovered all prespecified SC-only missing-type sets and Stage3B recovered dominant ST-only regions with a mean recall of 0.995 and precision of 0.974. The mean target-mass recall was 0.991, while the zero-target false-positive rate was 5.93e-4. Recovery was more conservative for mixed target-positive spots in the mouse-brain simulations, where positive-spot recall was 0.578 despite a dominant-region recall of 0.997.

> In the same nine simulations, SVTuner increased the mean abstention-aware whole-space recovery score from 0.515 for CytoSPACE to 0.710 under no added noise and achieved a higher score in all nine paired scenarios. Correct-abstention precision and recall were 0.995 and 0.954, respectively. Under 10% scRNA-seq expression perturbation, the mean score increased from 0.505 to 0.622, with improvements in seven of nine scenarios; correct-abstention precision and recall decreased to 0.961 and 0.904.

> Across six real reference-dropout settings, 93.8% of Stage3B blank spots were located within or near independently defined top-15% marker cores, whereas mean core recall was 57.8%. This precision-coverage pattern indicates that Stage3B preferentially abstains in strongly unsupported spatial compartments rather than exhaustively masking all spots associated with the removed type.

> In the ST8059051 thalamic case, Stage3B overlapped 347 of 362 marker-defined target spots, yielding a recall of 0.959 and precision of 0.720. An abstention-aware evaluation increased accuracy from 0.329 for CytoSPACE to 0.572 for SVTuner, while revealing a modest reduction in supported-region correct calls (425 to 391). Without abstention, 52.8% of target-region spots were forced to inhibitory neurons and 24.3% to astrocytes.

## 13. 关键证据文件

```text
visualizations/simulations/composite_stage3a_stage3b_recheck.csv
visualizations/method_comparison/composite_no_noise/*abstention_aware_summary.csv
visualizations/method_comparison/composite_no_noise/*abstention_aware_scenario.csv
visualizations/method_comparison/composite_scnoise10/*abstention_aware_summary.csv
visualizations/method_comparison/composite_scnoise10/*abstention_aware_scenario.csv
visualizations/method_comparison/composite_scnoise10/composite_scnoise10_execution_audit.csv
visualizations/cell2location_stage3b_case/thalamic_top15/*summary.json
visualizations/stage3b_realdata_candidate_scan/panel_a_blank_composition/*summary.csv
visualizations/stage3b_realdata_candidate_scan/stage3b_threshold_robustness/*summary.csv
visualizations/stage3b_reference_missing_stress/*forced_assignment.csv
visualizations/stage3b_reference_missing_stress/*expression_similarity.csv
visualizations/stage3b_false_spatial_niche/candidate_validation/brca_her2_plasma/*summary.csv
visualizations/stage3b_communication_validation/brca_her2_ffpe_plasma/strict_reference_style/*.csv
```
