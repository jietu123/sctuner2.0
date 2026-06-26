# 论文材料准备B：Stage3B 与 ST-only unsupported-region 实验材料

本文档记录项目中围绕 Stage3B 开发、验证和保留下来的实验材料。这里的 Stage3B 指的是：当某些细胞类型在 ST 数据中存在，但在 scRNA-seq reference 中不存在时，SVTuner 应识别这些不被 reference 支持的 ST 空间区域，并在后续映射中将对应位置保留为空白/abstention，而不是强行分配给 reference 中剩余的细胞类型。

与 Stage3A 的区别：

- Stage3A 处理的是 `SC-only unsupported type`：scRNA-seq reference 中存在，但 ST 数据中不存在或缺乏支持的类型，应在 SC reference 侧过滤。
- Stage3B 处理的是 `ST-only unsupported region`：ST 数据中存在，但 scRNA-seq reference 中不存在的类型或区域，应在 ST 空间侧留白。

核心目标：

1. 避免 CytoSPACE 等强制映射方法把 ST-only 区域错误分配给相似或邻近的 reference 细胞类型。
2. 在不使用白名单的前提下，由 Stage3B 自动检测 unsupported ST spots。
3. 在模拟数据和真实数据上证明：Stage3B 能识别需要留白的空间区域，并减少由强制映射导致的虚假解释。

## 1. Stage3B 核心实现

### 1.1 核心脚本

主实现：

```text
src/stages/stage3b_st_unsupported.py
```

相关 CLI / pipeline 调用：

```text
src/svtuner/cli.py
scripts/generate_st_only_reference_dropout_from_sim.py
scripts/visualize_stage3b_st_only_triptych.py
scripts/evaluate_stage3b_simulation.py
```

### 1.2 设计逻辑

Stage3B 的基本思想是：当 reference 中缺少某一类 ST 中实际存在的细胞类型时，该区域不能被简单地强制映射到现有 reference 类型上。Stage3B 通过 ST 表达结构、reference 支持程度、局部空间一致性和相似类型干扰处理来判断哪些 ST spots 属于 unsupported region。

核心输出形式：

- 对每个 spot 给出是否属于 unsupported / blank region 的判定。
- 后续 mapping 中，这些 spot 不再被强行映射，而是保留为空白。
- 可视化中通常用灰色或 blank 标记这些区域。

### 1.3 关键改进

在真实 cell2location mouse brain case 中，最初 Stage3B 对 `Oligodendrocyte/OPC` 区域识别不足。原因不是简单阈值问题，而是：

1. 目标类型和 reference 中相似类型之间存在表达相似性，容易被误判为仍有支持。
2. 原始逐 spot 判定容易受边界 spot 和相似类型互相解释的影响。
3. 对真实数据中目标区域定义的 top marker threshold 会影响评估口径。

后续优化方向：

- 增加相似类型诊断与抑制逻辑。
- 使用空间一致性约束减少零散误判。
- 在真实数据评估中使用 marker-defined region，例如 top15% 或 top20%，区分算法检测和评估定义。

当前主分析中以 top15% 作为主要真实目标区域定义，top20% 保留为敏感性/备选分析。

## 2. 复合模拟数据实验：Stage3A + Stage3B 同时存在

### 2.1 实验目的

原始模拟实验只覆盖了 Stage3A：即 scRNA-seq reference 中存在但 ST 中不存在的类型。Stage3B 新增后，需要构造复合型模拟场景，同时满足：

- ST 中存在但 SC reference 中不存在的细胞类型，需要 Stage3B 识别并留白。
- SC reference 中存在但 ST 中不存在的目标类型，需要 Stage3A 正确过滤。
- 无缺失、单缺失、双缺失三类 Stage3A 场景都能与 Stage3B 留白目标共存。

因此每个数据集保留三组复合型模拟场景：

1. Stage3A 无缺失 + Stage3B ST-only 留白。
2. Stage3A 单类型缺失 + Stage3B ST-only 留白。
3. Stage3A 双类型缺失 + Stage3B ST-only 留白。

### 2.2 BRCA 复合模拟场景

保留场景：

```text
real_brca7_endothelial_marker_control_sc_missing_endothelial_cells
real_brca7_endothelial_marker_missing_epithelial_cells_sc_missing_endothelial_cells
real_brca7_endothelial_marker_missing_epithelial_cells_pcs_sc_missing_endothelial_cells
```

实验含义：

- Stage3B 留白目标：`Endothelial cells`。
- Stage3A 单缺失目标：`Epithelial cells`。
- Stage3A 双缺失目标：`Epithelial cells + PCs`。

保留可视化：

```text
visualizations/simulations/real_brca/real_brca7_endothelial_marker_control_sc_missing_endothelial_cells/st_only_stage3b_triptych.png
visualizations/simulations/real_brca/real_brca7_endothelial_marker_missing_epithelial_cells_sc_missing_endothelial_cells/st_only_stage3b_triptych.png
visualizations/simulations/real_brca/real_brca7_endothelial_marker_missing_epithelial_cells_pcs_sc_missing_endothelial_cells/st_only_stage3b_triptych.png
```

配套配置：

```text
configs/datasets/real_brca7_endothelial_marker_control_sc_missing_endothelial_cells.yaml
configs/datasets/real_brca7_endothelial_marker_missing_epithelial_cells_sc_missing_endothelial_cells.yaml
configs/datasets/real_brca7_endothelial_marker_missing_epithelial_cells_pcs_sc_missing_endothelial_cells.yaml
```

可用于论文说明：

- 同一个 ST-only 留白目标可以嵌入 Stage3A 无缺失、单缺失、双缺失场景。
- Stage3B 留白不依赖人工白名单，而由算法输出。
- Stage3A 与 Stage3B 可以同时运行，互不破坏。

### 2.3 Human Lung 5-location 复合模拟场景

保留场景：

```text
human_lung_5loc_fine9_clustered_sim_sc_missing_b_cell
human_lung_5loc_fine9_clustered_sim_missing_at2_sc_missing_b_cell
human_lung_5loc_fine9_clustered_sim_missing_at2_fibroblast_sc_missing_b_cell
```

实验含义：

- Stage3B 留白目标：`B cell`。
- Stage3A 单缺失目标：`AT2`。
- Stage3A 双缺失目标：`AT2 + Fibroblast`。

保留可视化：

```text
visualizations/simulations/human_lung_5loc/human_lung_5loc_fine9_clustered_sim_sc_missing_b_cell/st_only_stage3b_triptych.png
visualizations/simulations/human_lung_5loc/human_lung_5loc_fine9_clustered_sim_missing_at2_sc_missing_b_cell/st_only_stage3b_triptych.png
visualizations/simulations/human_lung_5loc/human_lung_5loc_fine9_clustered_sim_missing_at2_fibroblast_sc_missing_b_cell/st_only_stage3b_triptych.png
```

配套配置：

```text
configs/datasets/human_lung_5loc_fine9_clustered_sim_sc_missing_b_cell.yaml
configs/datasets/human_lung_5loc_fine9_clustered_sim_missing_at2_sc_missing_b_cell.yaml
configs/datasets/human_lung_5loc_fine9_clustered_sim_missing_at2_fibroblast_sc_missing_b_cell.yaml
```

可用于论文说明：

- 在 Lung 这类空间结构明显的数据中，Stage3B 可以在保留 Stage3A 过滤结果的同时识别 ST-only 留白区域。
- 该组可以和原 Stage3A Lung 无缺失/单缺失/双缺失场景一起讲，形成完整的 mismatch 类型扩展。

### 2.4 Mouse Brain refined 复合模拟场景

保留场景：

```text
mouse_brain_refined7_balanced_clustered_sim_sc_missing_ext_l56
mouse_brain_refined7_balanced_clustered_sim_missing_micro_fill_inh_pvalb_sc_missing_ext_l56
mouse_brain_refined7_balanced_clustered_sim_missing_micro_oligo_2_fill_inh_pvalb_sc_missing_ext_l56
```

实验含义：

- Stage3B 留白目标：`Ext_L56`。
- Stage3A 单缺失目标：`Micro`，并使用 `Inh_Pvalb` 作为补充细胞类型。
- Stage3A 双缺失目标：`Micro + Oligo_2`，并使用 `Inh_Pvalb` 作为补充细胞类型。

保留可视化：

```text
visualizations/simulations/mouse_brain_refined/mouse_brain_refined7_balanced_clustered_sim_sc_missing_ext_l56/st_only_stage3b_triptych.png
visualizations/simulations/mouse_brain_refined/mouse_brain_refined7_balanced_clustered_sim_missing_micro_fill_inh_pvalb_sc_missing_ext_l56/st_only_stage3b_triptych.png
visualizations/simulations/mouse_brain_refined/mouse_brain_refined7_balanced_clustered_sim_missing_micro_oligo_2_fill_inh_pvalb_sc_missing_ext_l56/st_only_stage3b_triptych.png
```

配套配置：

```text
configs/datasets/mouse_brain_refined7_balanced_clustered_sim_sc_missing_ext_l56.yaml
configs/datasets/mouse_brain_refined7_balanced_clustered_sim_missing_micro_fill_inh_pvalb_sc_missing_ext_l56.yaml
configs/datasets/mouse_brain_refined7_balanced_clustered_sim_missing_micro_oligo_2_fill_inh_pvalb_sc_missing_ext_l56.yaml
```

可用于论文说明：

- Mouse brain 中存在更强的相似神经细胞类型干扰，适合说明 Stage3B 的相似类型处理必要性。
- 该组也用于证明 Stage3B 不是只在肿瘤组织或 Lung 结构中有效。

### 2.5 复合模拟总览与检查

保留汇总文件：

```text
visualizations/simulations/composite_stage3a_stage3b_recheck.csv
```

含义：

- 用于记录 9 个复合模拟场景中 Stage3A 与 Stage3B 的联合检查结果。
- 应作为论文结果整理时的内部证据表，不一定直接作为主图。

保留总览图：

```text
visualizations/simulations/simulation_triptych_overview_stack_3datasets.png
```

含义：

- 三个数据集的空间映射结果总览。
- 可作为补充材料或内部核验图。

## 3. 复合模拟多方法比较

### 3.1 实验目的

在复合型 Stage3A + Stage3B 模拟场景上，比较 SVTuner 与多种 mapping 方法的 composition recovery 表现。

需要注意：

- 该实验不应该作为 Stage3B 最强核心证据，因为传统 composition recovery 指标可能把 SVTuner 的正确留白视作未映射，从而低估优势。
- 该实验更适合作为补充分析，说明在复合 mismatch 条件下 SVTuner 整体表现没有崩溃，并可在 supported region 上进行较公平比较。

### 3.2 保留可视化

```text
visualizations/method_comparison/composite_no_noise/composition_recovery_7mapping_methods_composite_no_noise_boxplot.png
visualizations/method_comparison/composite_no_noise/composition_recovery_7mapping_methods_composite_no_noise_supported_regions_boxplot.png
```

配套 CSV：

```text
visualizations/method_comparison/composite_no_noise/composition_recovery_7mapping_methods_composite_no_noise_cell_type.csv
visualizations/method_comparison/composite_no_noise/composition_recovery_7mapping_methods_composite_no_noise_scenario.csv
visualizations/method_comparison/composite_no_noise/composition_recovery_7mapping_methods_composite_no_noise_summary.csv
visualizations/method_comparison/composite_no_noise/composition_recovery_7mapping_methods_composite_no_noise_supported_regions_cell_type.csv
visualizations/method_comparison/composite_no_noise/composition_recovery_7mapping_methods_composite_no_noise_supported_regions_scenario.csv
visualizations/method_comparison/composite_no_noise/composition_recovery_7mapping_methods_composite_no_noise_supported_regions_summary.csv
```

论文建议：

- 如果放正文，建议只放 supported-region-aware 版本。
- 如果 SVTuner 在全局 composition recovery 中不是第一，不应过度强调全局指标。
- 需要在图注中说明：blank region 的正确 abstention 不适合用传统强制映射准确率直接评价。

## 4. 真实 cell2location mouse brain case：Thalamic excitatory

### 4.1 实验目的

该实验是 Stage3B 在真实 ST 数据上的核心 case study。构造方式不是模拟 ST，而是：

1. 使用真实 ST 数据。
2. 从 scRNA-seq reference 中剔除目标细胞类型。
3. 检查 Stage3B 是否能识别目标类型在 ST 中富集、但 reference 不支持的空间区域。
4. 比较 baseline CytoSPACE 在该区域中的强制分配，以及 SVTuner 的留白结果。

最终定版目标：

```text
ST8059051 + Thalamic excitatory
```

主要阈值：

```text
top15% marker-defined target region
```

top20% 作为敏感性和备选结果保留。

### 4.2 top15 主结果

主目录：

```text
visualizations/cell2location_stage3b_case/thalamic_top15/
```

#### 面板 A：真实目标类型空间信号

```text
visualizations/cell2location_stage3b_case/thalamic_top15/cell2location_ST8059051_thalamic_excitatory_marker_region.png
```

含义：

- 展示真实 ST 数据中 `Thalamic excitatory` marker-rich 区域。
- 目标区域由 marker score 的 top15% 定义。
- 该图不是 Stage3B 输出，而是真实数据目标区域定义。

#### 面板 B：Stage3B 留白识别区域

```text
visualizations/cell2location_stage3b_case/thalamic_top15/cell2location_ST8059051_thalamic_excitatory_stage3b_miss_diagnostics_spatial.png
```

含义：

- 展示 Stage3B 识别出的 blank / unsupported region。
- 可与面板 A 对比，判断 Stage3B 是否主要落在真实目标类型富集区域。

#### 面板 C：baseline 强制映射空间分布

```text
visualizations/cell2location_stage3b_case/thalamic_top15/cell2location_ST8059051_thalamic_excitatory_baseline_forced_assignment_region.png
```

含义：

- 只关注真实目标区域内，CytoSPACE baseline 被迫映射成哪些 reference 类型。
- 这个图说明：如果没有 Stage3B 留白，目标区域会被解释为 reference 中剩余类型。

#### 面板 D：baseline 强制映射组成比例

```text
visualizations/cell2location_stage3b_case/thalamic_top15/cell2location_ST8059051_thalamic_excitatory_forced_assignment_composition.png
```

含义：

- 对面板 C 的量化汇总。
- 展示目标区域中 baseline 强制分配到各 reference 类型的比例。

#### 面板 E：abstention-aware accuracy

```text
visualizations/cell2location_stage3b_case/thalamic_top15/cell2location_ST8059051_thalamic_excitatory_abstention_aware_region_accuracy.png
```

含义：

- 用区域感知方式评价 baseline 与 SVTuner。
- supported region 中正常评价映射是否合理。
- unsupported target region 中，如果应该留白且模型正确留白，则计为正确。
- 避免把正确 blank 误算为 mapping failure。

### 4.3 top20 备选结果

目录：

```text
visualizations/cell2location_stage3b_case/thalamic_top20/
```

用途：

- 用作阈值敏感性分析。
- top20 比 top15 更严格地要求 Stage3B 识别更多目标区域，因此通常更难。
- 不建议作为主图，但可以作为补充材料或内部审查。

保留主图类型与 top15 相同：

```text
marker_region.png
stage3b_miss_diagnostics_spatial.png
baseline_forced_assignment_region.png
forced_assignment_composition.png
abstention_aware_region_accuracy.png
```

## 5. 真实多场景 Stage3B reference-dropout 验证

### 5.1 实验目的

单一真实 case 容易被认为是特例。因此构建多个真实数据 reference-dropout 场景，检验 Stage3B blank region 是否与真实目标 marker-defined region 一致。

最终保留 6 个推荐场景：

1. Mouse embryo：`Endoderm/Gut`
2. Mouse embryo：`Erythroid`
3. BRCA TNBC：`Plasma cells`
4. CRC fresh frozen：`B cells`
5. BRCA HER2 FFPE：`Plasma cells`
6. BRCA HER2 FFPE：`Epithelial cells`

### 5.2 空间对比图

保留主图：

```text
visualizations/stage3b_realdata_candidate_scan/spatial_9x2/stage3b_reference_dropout_spatial_stack_recommended_6x2.png
```

含义：

- 每个场景一行。
- 左列：真实 ST target marker signal，并用轮廓标出 marker top15 region。
- 右列：Stage3B blank result，并用轮廓标出 Stage3B 留白区域。
- 该图主要证明 Stage3B 识别的 blank region 与真实 target marker-rich region 空间上高度对应。

配套文件：

```text
visualizations/stage3b_realdata_candidate_scan/spatial_9x2/stage3b_reference_dropout_spatial_stack_recommended_6x2_manifest.csv
visualizations/stage3b_realdata_candidate_scan/spatial_9x2/stage3b_reference_dropout_spatial_stack_recommended_6x2_metadata.csv
visualizations/stage3b_realdata_candidate_scan/spatial_9x2/stage3b_reference_dropout_spatial_stack_recommended_6x2_scenes/
```

### 5.3 新 Panel A：blank-region composition accuracy

保留主图：

```text
visualizations/stage3b_realdata_candidate_scan/panel_a_blank_composition/stage3b_reference_dropout_panel_a_blank_region_composition.png
```

含义：

- 替代早期 violin 设计。
- 不再强调误差大小，而是直接展示 Stage3B blank region 与 marker-defined target region 的组成关系。
- 用来回答：Stage3B 留白区域到底有多少来自真实目标区域。

配套表：

```text
visualizations/stage3b_realdata_candidate_scan/panel_a_blank_composition/stage3b_reference_dropout_panel_a_blank_region_composition_summary.csv
```

论文建议：

- 该图适合与 6x2 空间图配套使用。
- 6x2 空间图提供直观证据，Panel A 提供跨场景定量总结。

### 5.4 阈值鲁棒性分析

保留目录：

```text
visualizations/stage3b_realdata_candidate_scan/stage3b_threshold_robustness/
```

主图：

```text
visualizations/stage3b_realdata_candidate_scan/stage3b_threshold_robustness/stage3b_threshold_robustness_summary.png
```

含义：

- 比较 top10%、top15%、top20%、top25% 等 target-region 定义下，Stage3B blank region 与真实目标区域的关系是否稳定。
- 该实验来自对 top15/top20 选择问题的延伸。

论文建议：

- 可作为补充材料。
- 如果正文空间有限，可以不放，但可以在方法或补充说明中引用。

## 6. Reference-missing stress test：多目标强制映射与表达相似性

### 6.1 实验目的

该实验不是空间图，而是解释机制：当 reference 缺失某些 ST-only 类型时，baseline CytoSPACE 会把这些区域强制分配给哪些剩余 reference 类型，以及这种错误分配是否与表达相似性有关。

保留主图：

```text
visualizations/stage3b_reference_missing_stress/cell2location_reference_missing_multitarget_panels_bc.png
```

### 6.2 Panel B：forced assignment heatmap

配套文件：

```text
visualizations/stage3b_reference_missing_stress/cell2location_reference_missing_multitarget_panels_bc_cytospace_forced_assignment.csv
```

含义：

- 行：被从 reference 中移除的 ST-only 目标类型。
- 列：baseline CytoSPACE 在目标区域中实际强制分配的 reference 类型。
- 值：分配比例或加权比例。

该图用于说明：没有 Stage3B 时，目标区域不会自然留白，而是会被 baseline 分配给某些相似或常见 reference 类型。

### 6.3 Panel C：expression similarity heatmap

配套文件：

```text
visualizations/stage3b_reference_missing_stress/cell2location_reference_missing_multitarget_panels_bc_expression_similarity.csv
```

含义：

- 行：ST-only 目标类型。
- 列：reference 中剩余类型。
- 值：表达相似性。

该图解释 Panel B：baseline 的错误分配往往不是随机的，而是受表达相似性影响。因此 Stage3B 需要显式处理相似类型干扰，而不是只靠简单阈值。

### 6.4 thalamic top15 B/C 备选图

保留但不是主线：

```text
visualizations/stage3b_reference_missing_stress/cell2location_thalamic_top15_reference_missing_panels_bc.png
```

用途：

- 这是单目标 `Thalamic excitatory` case 的 B/C 热图。
- 当前主线更推荐使用 multitarget 版本。

## 7. False spatial niche / neighborhood coupling 实验

### 7.1 实验目的

该实验属于更深一层的下游验证：不仅说明 baseline 把 unsupported region 映射错，还说明这种错误映射可能制造虚假的空间邻域结构或细胞类型 coupling 关系。

选择场景：

```text
BRCA HER2 FFPE + Plasma cells reference dropout
```

保留主图：

```text
visualizations/stage3b_false_spatial_niche/candidate_validation/brca_her2_plasma/cytospace_fig2d_tme_brca_her2_ffpe_sc_missing_plasma_cells_ncem_style.png
```

配套文件：

```text
visualizations/stage3b_false_spatial_niche/candidate_validation/brca_her2_plasma/cytospace_fig2d_tme_brca_her2_ffpe_sc_missing_plasma_cells_coupling_matrices.csv
visualizations/stage3b_false_spatial_niche/candidate_validation/brca_her2_plasma/cytospace_fig2d_tme_brca_her2_ffpe_sc_missing_plasma_cells_pair_statistics.csv
visualizations/stage3b_false_spatial_niche/candidate_validation/brca_her2_plasma/cytospace_fig2d_tme_brca_her2_ffpe_sc_missing_plasma_cells_summary.csv
```

### 7.2 图的核心解释

该图借鉴 NCEM / spatial graph 类型论文的可视化风格，重点不是复现原论文生物结论，而是把其空间邻域分析逻辑迁移到 Stage3B 问题上。

核心比较：

- Full reference：作为参考基线。
- CytoSPACE dropout：Plasma cells 从 reference 中删除后，baseline 被迫将目标区域分配给其他类型。
- SVTuner Stage3B：unsupported region 被留白，减少由强制分配造成的虚假邻域组成变化。

论文可用结论：

- Reference dropout 会诱导 baseline 出现 reference-relative neighborhood composition changes。
- Stage3B abstention 可以减少这类由 reference 缺失造成的虚假空间邻域解释。

### 7.3 compact 版本

保留但不是主线：

```text
visualizations/stage3b_false_spatial_niche/brca_her2_plasma_compact/
```

用途：

- 这是精简版尝试。
- 当前更建议使用 `candidate_validation/brca_her2_plasma` 下的 NCEM-style 图。

## 8. 通信下游实验：false communication validation

### 8.1 实验目的

该实验是 Stage3B 的下游升级实验：验证 reference 缺失导致的强制映射是否会进一步诱导虚假的 spatial communication / ligand-receptor 解释。

选择场景：

```text
BRCA HER2 FFPE + Plasma cells reference dropout
```

保留主图：

```text
visualizations/stage3b_communication_validation/brca_her2_ffpe_plasma/strict_reference_style/stage3b_communication_strict_reference_style.png
```

配套脚本：

```text
scripts/run_stage3b_spatial_communication_inference.py
scripts/validate_stage3b_communication_with_observed_st.py
scripts/plot_stage3b_communication_strict_reference_style.py
```

配套数据：

```text
data/processed/stage3b_spatial_communication/brca_her2_ffpe_plasma/
data/processed/stage3b_communication_downstream_validation/brca_her2_ffpe_plasma/
```

### 8.2 Panel A：dropout-induced false communication pathway enrichment

配套文件：

```text
visualizations/stage3b_communication_validation/brca_her2_ffpe_plasma/strict_reference_style/panel_a_spatialdm_fig2c_pathway_dotplot_values.csv
```

含义：

- 统计 reference dropout 后被诱导出来的 false communication LR pairs。
- 按 pathway 汇总，展示哪些通路更容易受到 reference 缺失和强制映射影响。

注意：

- 当前 Panel A 点数不多，因为候选 LR pairs 被汇总到有限 pathway 中。
- 它更适合作为通信实验的入口说明，而不是唯一核心证据。

### 8.3 Panel B：local LR hotspot comparison

配套文件：

```text
visualizations/stage3b_communication_validation/brca_her2_ffpe_plasma/strict_reference_style/panel_b_spatialdm_fig3f_hotspot_summary.csv
```

含义：

- 借鉴 SpatialDM 局部热点图风格。
- 展示特定 LR pair 在 CytoSPACE dropout 与 SVTuner Stage3B 下的 local Moran / hotspot 分布。
- 证明 Stage3B 留白后，某些由强制映射产生的局部通信热点被削弱。

当前代表 pair：

```text
KITLG -> KIT
```

### 8.4 Panel C：ligand-target evidence heatmap

配套文件：

```text
visualizations/stage3b_communication_validation/brca_her2_ffpe_plasma/strict_reference_style/panel_c_renoir_fig4gh_ligand_target_matrix.csv
visualizations/stage3b_communication_validation/brca_her2_ffpe_plasma/strict_reference_style/panel_c_renoir_fig4gh_ligand_score.csv
visualizations/stage3b_communication_validation/brca_her2_ffpe_plasma/strict_reference_style/panel_c_ligand_annotations.csv
```

含义：

- 借鉴 Renoir Fig.4g/h 的 ligand-target heatmap 视觉逻辑。
- 展示 dropout-induced false communication 中 ligand 与 downstream target evidence 的关系。
- 该图用于把空间通信候选和下游表达响应联系起来。

注意：

- 当前图是风格借鉴，不是严格复现 Renoir 数据结构。
- ligand score 与 target regulatory potential 来自本项目通信验证结果。

## 9. 论文图组织建议

### 9.1 推荐正文主线

建议将 Stage3B 部分组织为四层证据：

1. 复合模拟数据：证明 Stage3A 与 Stage3B 可以同时运行。
2. 真实多场景 reference-dropout：证明 Stage3B blank region 与真实 marker-defined target region 一致。
3. 单真实 case 深入：展示 baseline 强制分配、SVTuner 留白、abstention-aware accuracy。
4. 下游验证：证明错误强制映射会影响 neighborhood / communication 解释，而 Stage3B 可以缓解。

### 9.2 推荐主图候选

主图候选 1：Stage3B 多场景空间验证

```text
visualizations/stage3b_realdata_candidate_scan/spatial_9x2/stage3b_reference_dropout_spatial_stack_recommended_6x2.png
visualizations/stage3b_realdata_candidate_scan/panel_a_blank_composition/stage3b_reference_dropout_panel_a_blank_region_composition.png
```

主图候选 2：Thalamic excitatory case study

```text
visualizations/cell2location_stage3b_case/thalamic_top15/cell2location_ST8059051_thalamic_excitatory_marker_region.png
visualizations/cell2location_stage3b_case/thalamic_top15/cell2location_ST8059051_thalamic_excitatory_stage3b_miss_diagnostics_spatial.png
visualizations/cell2location_stage3b_case/thalamic_top15/cell2location_ST8059051_thalamic_excitatory_baseline_forced_assignment_region.png
visualizations/cell2location_stage3b_case/thalamic_top15/cell2location_ST8059051_thalamic_excitatory_forced_assignment_composition.png
visualizations/cell2location_stage3b_case/thalamic_top15/cell2location_ST8059051_thalamic_excitatory_abstention_aware_region_accuracy.png
```

主图候选 3：下游通信验证

```text
visualizations/stage3b_communication_validation/brca_her2_ffpe_plasma/strict_reference_style/stage3b_communication_strict_reference_style.png
```

### 9.3 推荐补充材料

复合模拟 9 场景：

```text
visualizations/simulations/real_brca/*_sc_missing_endothelial_cells/st_only_stage3b_triptych.png
visualizations/simulations/human_lung_5loc/*_sc_missing_b_cell/st_only_stage3b_triptych.png
visualizations/simulations/mouse_brain_refined/*_sc_missing_ext_l56/st_only_stage3b_triptych.png
```

阈值鲁棒性：

```text
visualizations/stage3b_realdata_candidate_scan/stage3b_threshold_robustness/stage3b_threshold_robustness_summary.png
```

多方法复合模拟：

```text
visualizations/method_comparison/composite_no_noise/composition_recovery_7mapping_methods_composite_no_noise_supported_regions_boxplot.png
```

多目标强制映射解释：

```text
visualizations/stage3b_reference_missing_stress/cell2location_reference_missing_multitarget_panels_bc.png
```

## 10. 当前保留但需谨慎使用的材料

### 10.1 top20 thalamic case

```text
visualizations/cell2location_stage3b_case/thalamic_top20/
```

说明：

- top20 是更严格的 target-region 定义。
- 可以证明结论对阈值并非完全依赖，但主图建议使用 top15。

### 10.2 Panel A violin 旧版

```text
visualizations/stage3b_realdata_candidate_scan/panel_a_violin/
```

说明：

- 这是早期设计，用 marker percentile violin 比较真实目标区域与 Stage3B blank region。
- 视觉上 Stage3B 的低分尾部较明显，容易削弱结论。
- 当前已被 blank-region composition Panel A 替代。

### 10.3 thalamic-only B/C heatmap

```text
visualizations/stage3b_reference_missing_stress/cell2location_thalamic_top15_reference_missing_panels_bc.png
```

说明：

- 单目标版本解释力不如 multitarget。
- 如果版面有限，优先使用 multitarget panels B/C。

### 10.4 composite method comparison

```text
visualizations/method_comparison/composite_no_noise/
```

说明：

- 可以作为补充，但不建议作为 Stage3B 的最核心证据。
- 原因是传统 mapping 指标天然偏向强制映射方法，可能低估正确留白。

## 11. 需要避免的表述

避免说：

- Stage3B 在所有真实数据中完全准确识别所有 ST-only spots。
- Stage3B 的优势主要体现在传统全局 composition recovery 上。
- 所有下游 communication 分析都是严格 ligand-receptor 生物实验验证。

建议说：

- Stage3B identifies reference-unsupported ST regions and abstains from forced assignment.
- Correct abstention should be credited in unsupported regions rather than treated as mapping failure.
- In reference-dropout experiments, Stage3B blank regions are enriched for marker-defined target regions across multiple real ST datasets.
- Forced assignment can create misleading neighborhood or communication patterns, while Stage3B reduces this failure mode by leaving unsupported regions unassigned.

## 12. 当前阶段结论

目前 Stage3B 相关实验已经形成完整证据链：

1. 模拟数据证明：Stage3B 可以与 Stage3A 同时处理双向 mismatch。
2. 真实多场景证明：Stage3B blank region 与 marker-defined ST-only target region 空间一致。
3. 单 case study 证明：baseline 在 unsupported region 中会强制分配，而 SVTuner 可以留白。
4. 机制热图证明：baseline 的强制分配与表达相似性有关，相似类型会造成混淆。
5. 下游实验证明：错误强制分配可能进一步影响 spatial niche 和 communication 解释。

这些材料可以支持论文中一个独立的 Stage3B 小节：`Detection and abstention of reference-unsupported ST regions`。
