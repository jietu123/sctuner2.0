# Stage6 真实数据审计设计（含项目实际路径）

## 0. Stage6 一句话定义

Stage6 = 在真实数据（无真值）上，对 baseline 与 route2(V5.3.1) 的映射输出做“工程一致性 + 分布合理性 + 弱监督 marker 一致性”的审计验收，并产出结构化 JSON（可直接被论文/Stage7 复用）。

---

## 1. 目标与非目标

### 1.1 目标（必须完成）

1. 工程正确性（IO 三件套一致、ID 对齐、行和合理、meta 可追溯、容量/ledger 对账不穿帮）
2. 分布合理性（不崩盘：不是全单一类型、不是全零、不是极端硬分配；类型质量/出现范围合理）
3. 弱监督一致性（marker 与预测类型比例不“反着来”；至少多数 marker 呈正向趋势）

### 1.2 非目标（明确不做）

- 不计算 L1/JS/corr（无真值）
- 不对结果做任何再优化/重分配（Stage6 只读审计）
- 不画最终论文图（Stage7 画），但 Stage6 必须落盘“可画图的数据指标”

---

## 2. 输入与约定（数据契约）

Stage6 对每个 run（baseline / route2）读取同一套“映射产物”。marker 检验需要 ST 表达矩阵（来自 Stage1 的产物）。

### 2.1 Run 必需输入（每个 run 一套）

**A. spot_type_fraction（必需）**

- 形态：宽表（每行一个 spot，每列一个 type 的 fraction）
- 必须列：`spot_id` + `type_*`（或任意 type 列集合）
- 约束：fraction >= 0；每行 sum≈1（允许浮点误差）
- 实际路径（优先使用 run 归档目录）：
  - baseline: `D:\Experiment\SVTuner_Project\result\real_brca\stage4_cytospace\runs\baseline\cell_type_assignments_by_spot.csv`
  - route2: `D:\Experiment\SVTuner_Project\result\real_brca\stage4_cytospace\runs\route2_v5_3_1\cell_type_assignments_by_spot.csv`
- 若未归档，使用默认输出目录（仅作为回退）：
  - `D:\Experiment\SVTuner_Project\result\real_brca\stage4_cytospace\cytospace_output\cell_type_assignments_by_spot.csv`

**B. cell_assignment（必需）**

- 必须列：`cell_id`, `spot_id`
- 推荐列（若存在可做更强一致性检查）：`cell_type` / `sc_type` / `assigned_type`
- 实际路径（优先 run 归档目录）：
  - baseline: `D:\Experiment\SVTuner_Project\result\real_brca\stage4_cytospace\runs\baseline\cell_assignment.csv`
  - route2: `D:\Experiment\SVTuner_Project\result\real_brca\stage4_cytospace\runs\route2_v5_3_1\cell_assignment.csv`
- 回退路径：
  - `D:\Experiment\SVTuner_Project\result\real_brca\stage4_cytospace\cytospace_output\cell_assignment.csv`

**C. meta.json（必需）**

- 用于追溯：`module_sha1` / `runner_sha1` / `config_effective_subset` / `capacity_audit` 等
- 实际路径（优先 run 归档目录）：
  - baseline: `D:\Experiment\SVTuner_Project\result\real_brca\stage4_cytospace\runs\baseline\stage4_summary.json`
  - route2: `D:\Experiment\SVTuner_Project\result\real_brca\stage4_cytospace\runs\route2_v5_3_1\stage4_summary.json`
- 回退路径：
  - `D:\Experiment\SVTuner_Project\result\real_brca\stage4_cytospace\cytospace_output\stage4_summary.json`

**D. st_coords（推荐）**

- 列：`spot_id`, `row/col` 或 `x/y`
- 实际路径：
  - `D:\Experiment\SVTuner_Project\data\processed\real_brca\stage1_preprocess\exported\st_coordinates.csv`

**E. stage3_summary.json（route2 推荐 / baseline 可无）**

- 用于 ledger（filter/rescue 对账）
- 实际路径：
  - `D:\Experiment\SVTuner_Project\result\real_brca\stage3_typematch\stage3_summary.json`

> 备注：若有 `cell_spot_matrix`，Stage6 会把它作为可选增强一致性检查（不强依赖）。

### 2.2 Marker 检验额外输入（至少一种）

Stage6 需要 spot-level gene expression（用于相关/分组差）。

- 推荐输入：
  - `D:\Experiment\SVTuner_Project\data\processed\real_brca\stage1_preprocess\exported\st_expression_normalized.csv`

Stage6 只要求：能得到 `expr[spot_id, marker_gene]` 的数值矩阵即可。

---

## 3. 输出（文件与 JSON Schema）

Stage6 输出分“单 run 审计”与“两 run 对比”。

### 3.1 单 run 审计输出

- `D:\Experiment\SVTuner_Project\result\real_brca\stage6_real_audit\<run_id>\stage6_audit_<run_id>.json`
- `D:\Experiment\SVTuner_Project\result\real_brca\stage6_real_audit\<run_id>\stage6_plotdata_<run_id>.csv`（可选但强烈推荐）
  - 拼表：`spot_id + coords + top1_type + top1_fraction + entropy + (marker_expr可选)`

### 3.2 对比输出

- `D:\Experiment\SVTuner_Project\result\real_brca\stage6_real_audit\stage6_compare_baseline_vs_route2.json`
- `D:\Experiment\SVTuner_Project\result\real_brca\stage6_real_audit\stage6_plotdata_compare.csv`（可选）

---

## 4. 单 run 审计 JSON 结构（建议定稿版）

```json
{
  "stage": "stage6_real_audit",
  "sample": "real_brca",
  "run_id": "baseline|route2_v5_3_1",
  "created_at": "ISO_TIME",
  "input_fingerprint": {
    "spot_type_fraction_sha1": "...",
    "cell_assignment_sha1": "...",
    "meta_sha1": "...",
    "st_expr_sha1": "optional",
    "stage3_summary_sha1": "optional"
  },
  "meta_trace": {
    "module_sha1": "...",
    "runner_sha1": "...",
    "config_effective_subset": { "...": "..." },
    "knn_mode": "optional",
    "harden_method": "optional",
    "cells_per_spot": "optional"
  },
  "io_integrity": { "...": "..." },
  "ledger_integrity": { "...": "..." },
  "distribution_sanity": { "...": "..." },
  "weak_supervision_marker_check": { "...": "..." },
  "overall_status": "PASS|WARN|FAIL",
  "reasons": [ "..." ]
}
```

---

## 5. 具体审计模块设计

### 5.1 io_integrity（工程一致性）

**目的**：证明“输出三件套一致 + 没有明显表结构/ID/数值错误”。

建议字段：

- `n_spots_fraction`
- `n_types_fraction`
- `n_cells_assignment`
- `cell_id_unique_ok`
- `spot_id_intersection_ok`（fraction 与 coords 的 spot_id 交集是否完整）
- `assignment_spot_id_subset_fraction_ok`
- `fraction_nonneg_ok`
- `fraction_row_sum_stats`：`min/median/max` + `n_outside_tol`
- `fraction_has_nan_or_inf`
- `cells_per_spot_stats`：`min/median/max`（基于 cell_assignment 统计）

**可选增强：cell_spot_matrix 一致性**

- `cell_spot_matrix_present`
- 行集合/列集合一致性 + argmax spot 一致率

**门禁建议（FAIL 条件）**

- cell_id 非唯一
- assignment spot_id 不在 fraction spot_id 中（非空）
- fraction 行和大量不≈1（例如 `n_outside_tol > 0.5% * n_spots`）
- fraction 存在 NaN/Inf 或明显负值

---

### 5.2 ledger_integrity（账本一致性，route2 专属）

**目的**：证明 route2 的 filter/rescue 没有写了但没生效、没有对账失败。

输入依赖：`stage3_summary.json`（或 meta 里落的 ledger 字段）。

建议输出字段：

- `status`: `OK|SKIPPED|WARN|FAIL`
- `filter_ledger`:
  - `expected_n_filtered`
  - `observed_n_filtered`
  - `ledger_check_ok`
- `rescue_ledger`:
  - `rescued_types`
  - `rescued_cells_count`
  - `rescued_cells_marked_unknown`
  - `rescue_check_ok`
- `warnings`

baseline：没有 route2 ledger 就输出 `status="SKIPPED"` + reasons。

**门禁建议（route2 FAIL 条件）**

- `ledger_check_ok=false`
- `rescued_cells_marked_unknown>0`

---

### 5.3 distribution_sanity（分布合理性）

**目的**：无真值下证明输出不崩盘，并给 Stage7 出图数据。

#### 5.3.1 Spot 级统计（从 spot_type_fraction 得到）

对每个 spot 的概率向量 `p`：

- `entropy`：`H(p) = -sum(p_i log p_i)`
- 推荐归一化：`H_norm = H / log(K)`（K=type数）
- `top1_fraction = max(p)`
- `top1_type = argmax type`
- `n_active_types = count(p_i > eps)`

输出汇总：

- `spot_entropy_stats`: min/median/max/p05/p95
- `spot_top1_fraction_stats`: min/median/max/p05/p95
- `spot_active_type_count_stats`
- `spot_top1_type_frequency_top10`

#### 5.3.2 Type 级统计

对每个 type 列 `t`：

- `type_total_mass = sum(fraction_t)`
- `type_mean_mass = mean(fraction_t)`
- `type_presence_spots_{thr}`（thr=0.01/0.05/0.1）
- `type_max_fraction`

输出：

- `type_mass_table`（可截断 TopN，完整表写入 plotdata）

#### 5.3.3 分布“崩盘”判定（建议门禁）

FAIL：

- `median(top1_fraction) >= 0.98`
- `median(entropy_norm) <= 0.01`
- 或 `top1_type` 某一类型占比 > 95%

WARN：

- `median(top1_fraction) in [0.95,0.98)`
- `median(entropy_norm) in (0.01,0.05]`
- 大量 spot `n_active_types` 很低（p95 <= 2）

---

### 5.4 weak_supervision_marker_check（弱监督 marker 一致性）

**目的**：用 marker 表达作为弱监督信号，证明预测类型比例至少方向一致。

#### 5.4.1 marker 配置（YAML）

配置文件：

- `D:\Experiment\SVTuner_Project\configs\stage6_markers.yaml`

示例结构：

```yaml
version: 1
expr_source:
  kind: "csv"
  path: "D:/Experiment/SVTuner_Project/data/processed/real_brca/stage1_preprocess/exported/st_expression_normalized.csv"
  transform: "none"
marker_sets:
  - type: "T cells CD8"
    markers: ["CD3D", "CD3E", "CD8A"]
  - type: "B cells"
    markers: ["MS4A1", "CD79A"]
  - type: "Myeloid"
    markers: ["LYZ", "S100A8", "S100A9"]
  - type: "Epithelial"
    markers: ["EPCAM", "KRT19"]
test:
  method: "spearman_corr + top_bottom"
  top_quantile: 0.2
  bottom_quantile: 0.2
  min_spots: 200
  min_var: 1e-8
```

#### 5.4.2 统计方法

对每个 `{type, marker}`：

**A. Spearman 相关**

- `x = predicted_fraction(type)`
- `y = expr(marker)`
- 输出 `rho`, `n_spots_used`, `valid_spots_fraction`

**B. Top-vs-Bottom 分组差**

- 按 `x` 排序，取 top20% 与 bottom20%
- `delta_mean = mean(y_top) - mean(y_bottom)`
- 可选：`effect_size`、`pvalue`（Mann–Whitney U）

#### 5.4.3 marker 模块输出结构

- `status`: `OK|SKIPPED|WARN|FAIL`
- `tests[]`: `type`, `marker`, `rho_spearman`, `delta_mean_top_bottom`, `pvalue_optional`, `direction_ok`
- `summary`: `n_tests`, `n_direction_ok`, `direction_ok_fraction`, `pass_rule`

**门禁建议**

- marker 输入缺失：`SKIPPED`（overall 最高 WARN）
- `direction_ok_fraction < 0.5`：FAIL
- `direction_ok_fraction in [0.5,0.6)`：WARN

---

## 6. 两 run 对比（baseline vs route2）设计

输出：`stage6_compare_baseline_vs_route2.json`

### 6.1 delta 指标（最少集）

- `delta_spot_entropy_median`
- `delta_spot_top1_fraction_median`
- `global_type_mass_baseline` / `global_type_mass_route2`
- `delta_type_total_mass`（每个 type）
- `top1_type_shift`

### 6.2 漂移/异常检测（可选增强）

- `global_distribution_js`
- `spotwise_kl_stats`

### 6.3 对比结论字段

- `interpretation`（可解释）
- `warnings`

---

## 7. 总体通过标准（Stage6 Gate）

**PASS**

- io_integrity 全 OK
- route2 ledger OK（baseline 可 SKIPPED）
- distribution 无崩盘（不触发 FAIL）
- marker OK 或多数正向

**WARN**

- IO OK，但 marker 缺失/不足，或分布接近阈值

**FAIL**

- 任何 io_integrity 强失败
- route2 ledger 对账失败
- 分布崩盘
- marker 多数反向

---

## 8. CLI 与配置（带实际路径）

入口脚本（待实现）：`D:\Experiment\SVTuner_Project\src\stages\stage6_real_audit.py`

示例命令：

```powershell
python D:\Experiment\SVTuner_Project\src\stages\stage6_real_audit.py `
  --sample real_brca `
  --baseline_id baseline `
  --route2_id route2_v5_3_1 `
  --fraction_baseline D:\Experiment\SVTuner_Project\result\real_brca\stage4_cytospace\runs\baseline\cell_type_assignments_by_spot.csv `
  --assignment_baseline D:\Experiment\SVTuner_Project\result\real_brca\stage4_cytospace\runs\baseline\cell_assignment.csv `
  --meta_baseline D:\Experiment\SVTuner_Project\result\real_brca\stage4_cytospace\runs\baseline\stage4_summary.json `
  --fraction_route2 D:\Experiment\SVTuner_Project\result\real_brca\stage4_cytospace\runs\route2_v5_3_1\cell_type_assignments_by_spot.csv `
  --assignment_route2 D:\Experiment\SVTuner_Project\result\real_brca\stage4_cytospace\runs\route2_v5_3_1\cell_assignment.csv `
  --meta_route2 D:\Experiment\SVTuner_Project\result\real_brca\stage4_cytospace\runs\route2_v5_3_1\stage4_summary.json `
  --stage3_summary D:\Experiment\SVTuner_Project\result\real_brca\stage3_typematch\stage3_summary.json `
  --st_expr D:\Experiment\SVTuner_Project\data\processed\real_brca\stage1_preprocess\exported\st_expression_normalized.csv `
  --st_coords D:\Experiment\SVTuner_Project\data\processed\real_brca\stage1_preprocess\exported\st_coordinates.csv `
  --marker_config D:\Experiment\SVTuner_Project\configs\stage6_markers.yaml
```

---

## 9. 代码组织建议（最小可维护拆分）

```
src/stages/stage6_real_audit.py           # 入口：加载两 run、跑审计、写 json、写 plotdata
src/stages/stage6/
  io_checks.py                             # io_integrity
  ledger_checks.py                         # ledger_integrity
  distribution_checks.py                   # distribution_sanity + per-spot/per-type 表
  marker_checks.py                         # marker 一致性（相关 + 分组）
  compare.py                               # baseline vs route2 delta
  schema.py                                # 字段名/默认阈值/小工具
```

---

## 10. 最小验收用例（可做 CI）

1. Happy path：两套 run 产出 audit JSON，overall 至少 PASS/WARN
2. 缺列：spot_type_fraction 缺 `spot_id` 或 fraction 列 -> FAIL
3. 重复 cell_id -> FAIL
4. row sum 不为 1（构造一批 spot）-> WARN/FAIL
5. marker 缺失 -> SKIPPED（overall 最高 WARN）
6. route2 ledger 对账失败 -> FAIL

---

## 11. 运行注意事项（项目约定）

1. Stage4 默认输出目录会被覆盖，建议在运行 Stage4 后将产物复制到 `result\real_brca\stage4_cytospace\runs\<run_id>\` 归档，再执行 Stage6。
2. Stage6 是只读审计，不修改任何上游结果。
3. 如果 marker 配置缺失，Stage6 仍可输出 io/distribution/ledger，但 overall 最多 WARN。
