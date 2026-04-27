# HLN 双缺失场景中 CD4 Treg 自动缺失识别失败问题说明

## 背景

本项目的 Stage3 目标是在不通过命令行硬指定缺失类型的前提下，基于 scRNA 参考与 ST 表达之间的证据，自动识别模拟场景中应当被过滤的缺失细胞类型。随后 Stage4 在 `missing_detected_only` 范围内只过滤 Stage3 标记为缺失的类型，从而避免把非缺失类型误过滤。

在 real_brca 和 adult_mouse_kidney 的单缺失、双缺失验证中，Stage3 已经能够正确识别目标缺失类型，并且 Stage5 审计显示非缺失类型误过滤为 0。HLN 数据集的 NK cell 单缺失场景也已经可以被当前 Stage3 机制正确识别。

但是，在 Human Lymph Node 的双缺失场景中，仍然存在一个稳定缺陷：Stage3 能识别 NK cell，但无法稳定识别 CD4 Treg。

## 受影响场景

样本：

`human_lymph_node_sim_mt_nk_cd4treg_fill_cd4t`

该场景的设计目标是：

- 缺失类型：`NK cell`, `CD4 Treg`
- 补充类型：`CD4 T cell`
- 运行方式：`--missing_type __NO_MISSING__`
- 过滤范围：`--route2_filter_scope missing_detected_only`

因此，正确行为应当是：Stage3 自主识别出 `NK cell` 和 `CD4 Treg`，Stage4 route2 只过滤这两个被 Stage3 识别出的缺失类型。

## 实际观察

当前结果显示：

- Stage3 自动识别到：`NK cell`
- Stage3 未识别到：`CD4 Treg`
- Stage4 route2 实际过滤：302 个 `NK cell`
- Stage4 route2 未过滤：`CD4 Treg`
- Stage5 route2 仍检测到 `CD4 Treg` 泄漏
- Stage5 过滤审计显示非缺失类型误过滤为 0

这说明当前问题不是“误过滤”，而是“漏检”。

## 关键证据

### 单缺失 HLN NK cell 场景

样本：

`human_lymph_node_clustered_sim_mt_nk_cell_d100_fill_cd4_t_cell`

结果：

- Stage3 识别缺失类型：`NK cell`
- Stage4 route2 过滤细胞数：302
- Stage5 审计：`non_missing_filtered = 0`

该结果说明当前 Stage3 机制已经能够在 HLN 中识别 NK cell，不是 HLN 数据集整体不可用。

### 双缺失 HLN NK cell + CD4 Treg 场景

样本：

`human_lymph_node_sim_mt_nk_cd4treg_fill_cd4t`

结果：

- Stage3 识别缺失类型：`NK cell`
- Stage3 未识别：`CD4 Treg`
- Stage4 route2 过滤细胞数：302
- Stage5 route2 中仍存在 `CD4 Treg` assignment
- Stage5 审计：`non_missing_filtered = 0`

该结果说明当前机制在双缺失场景中仍然是保守的：它没有误删非缺失类型，但也没有完全识别出所有真实缺失类型。

## 原因分析

CD4 Treg 漏检的核心原因是：`CD4 Treg` 与替代类型 `CD4 T cell` 以及其他 T cell 相关亚群在表达空间中高度相似，导致缺失后仍然能从 ST 中获得较高的支持度。

具体表现为：

1. `CD4 Treg` 的整体表达 profile 与 `T cell` / `CD4 T cell` 非常接近。
2. 当 `CD4 Treg` 在模拟 ST 中被 `CD4 T cell` 补充后，许多 T cell 相关 marker 仍然保留。
3. Stage3 的全局 support score 因此不会把 `CD4 Treg` 判定为低支持异常值。
4. 新增的 identity-marker depletion 证据能够识别 NK cell，但对 CD4 Treg 仍不够强，因为 CD4 Treg 的许多可检测 marker 与 T cell 背景高度重叠。
5. 如果为了识别 CD4 Treg 而继续放宽阈值，会引入 Myofibroblast、Vascular smooth muscle cell 等非缺失类型被误过滤的风险。

因此，当前选择是保守处理：宁可保留 CD4 Treg 漏检，也不通过过度激进的阈值造成非缺失类型误过滤。

## 与已成功场景的区别

real_brca 和 adult_mouse_kidney 中的缺失类型与补充类型之间通常存在更明显的表达差异，缺失后支持度会明显下降，因此 Stage3 的自适应低支持机制能够识别缺失类型。

HLN 的 NK cell 也能被识别，是因为其 identity marker 与 T cell 背景之间仍有足够区分度。

但 CD4 Treg 属于更细粒度的 T cell 亚型，且补充类型正是 `CD4 T cell`。这种“高度相似亚型被相近主类替代”的场景，是当前自动缺失识别机制的困难边界。

## 对当前项目结论的影响

这个问题不影响以下结论：

- Stage3 可以在多个数据集和多个模拟场景中实现自动缺失类型识别。
- 当前机制能够避免明显的非缺失类型误过滤。
- real_brca、adult_mouse_kidney 以及 HLN NK cell 单缺失场景的过滤结论仍然成立。

但这个问题限制了以下结论：

- 不能声称当前 Stage3 能稳定识别所有高度相似免疫亚型的缺失。
- 不能把 HLN 双缺失场景作为“完全自动识别所有缺失类型”的正例。
- 对 CD4 Treg 这类细粒度亚型，应当在论文中作为当前方法的已知局限进行说明。

## 建议在论文中的表述口径

可以将该问题表述为：

> In highly similar immune subtypes, such as CD4 Treg replaced by CD4 T cells in the Human Lymph Node simulation, the current unsupervised missing-type detector may fail to distinguish true subtype loss from retained pan-T-cell expression. The detector is intentionally conservative to avoid false removal of non-missing cell types, resulting in a false-negative missing-type detection for CD4 Treg in this scenario.

中文口径可以写作：

> 在高度相似的免疫亚型中，例如 Human Lymph Node 数据集中由 CD4 T cell 补充 CD4 Treg 缺失区域时，当前无监督缺失类型识别机制可能无法区分真实亚型缺失与泛 T cell 表达信号保留。为了避免误过滤非缺失类型，当前策略保持保守，因此在该场景中出现了 CD4 Treg 的漏检。

## 后续可能改进方向

后续如果需要进一步解决该问题，可以考虑：

1. 引入人工确认的 lineage-specific marker panel，例如 FOXP3、IL2RA、CTLA4、IKZF2、TNFRSF18 等用于 CD4 Treg 的专门证据。
2. 在 T cell 亚型内部单独建立 subtype-level contrast，而不是只使用全局 support score。
3. 将缺失检测拆分为两层：粗粒度 lineage detection 和细粒度 subtype detection。
4. 对高度相似亚型采用弱监督或 marker-prior 约束，而不是完全无监督判断。

这些方向不属于当前版本的核心实现范围，当前版本将该问题记录为固定局限。
