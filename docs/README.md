# SVTuner 2.0 项目说明

SVTuner 是一个面向空间单细胞映射的类型感知协调层。它不是为了替代 CytoSPACE、Tangram、novoSpaRc、SpaOTsc 或 CellTrek 这类映射器，而是在正式映射之前增加一层 `type-aware / mismatch-aware` 判断：先判断单细胞参考中的每一种细胞类型是否应该被允许进入空间映射池，再把经过协调的输入交给下游映射方法。

项目要解决的核心问题是：单细胞参考中存在某类细胞，但空间转录组样本中并不支持或并不存在该类型时，传统表达匹配型映射器仍可能把这类细胞“合法地”分配到 ST spot 上。这个问题不是普通数值噪声，而是会在空间图上形成生物学上不应出现的伪分布。SVTuner 的目标就是在映射发生之前识别并处理这种类型错配。

当前代码库已经包含以下完整内容：

- Stage3 类型错配检测与 `plugin_type` 标签生成。
- CytoSPACE baseline 映射和 SVTuner + CytoSPACE route2 映射。
- 真实数据 profile-mask 掩盖实验。
- 可控缺失类型的模拟数据实验。
- 10% scRNA reference noise 噪声实验。
- Tangram、novoSpaRc、SpaOTsc、CellTrek 等对照方法接入。
- 真实数据掩盖可视化、模拟数据三联图和多方法箱线图。

## 1. 核心创新逻辑

SVTuner 的核心创新不是重新设计一个空间分配求解器，而是把“哪些细胞类型应该进入映射池”这个问题从下游映射器中独立出来。

当前逻辑可以分为两层：

1. Stage3 先根据 ST 表达证据判断每个单细胞类型是否被空间样本支持。
2. Stage4 再根据 Stage3 结果重建或标注单细胞池，然后运行实际映射。

在项目中：

- `Baseline` 指直接使用原始单细胞类型标签运行 CytoSPACE。
- `Route2` 指使用 Stage3 协调后的 `plugin_type` 和过滤结果运行 CytoSPACE。
- `SVTuner + CytoSPACE` 是当前主结果图中使用的方法名称。

理想情况下：

- 如果目标细胞类型在 ST 样本中确实缺失，Stage3 应该把它识别为 unsupported 或 missing。
- Route2 应该阻止该类型进入最终映射结果。
- 非缺失类型即使信号偏弱，也不应该被粗暴误过滤。
- 在同一个下游求解器 CytoSPACE 上比较时，baseline 和 route2 的差异应主要来自 SVTuner 的前置协调层。

## 2. 主流程结构

当前项目主流程可以概括为 Stage1、Stage3、Stage4 和评估可视化四个部分。

### 2.1 Stage1：标准化预处理

Stage1 把真实数据或模拟数据整理为统一格式，供后续 Python 阶段读取。核心导出文件通常位于：

```text
stage1_preprocess/exported/
```

常用文件包括：

```text
sc_expression_normalized.csv
sc_expression_data.csv
sc_expression_counts.csv
st_expression_normalized.csv
st_coordinates.csv
sc_metadata.csv
```

模拟数据通常从已有模拟源重建 Stage1：

```powershell
python scripts\prepare_stage1_from_sim_source.py --project_root . --sample <sample>
```

真实低分辨率数据通常使用 R 版 Stage1：

```powershell
Rscript r_scripts\stage1_preprocess.R --sample <sample> --project_root .
```

### 2.2 Stage3：类型错配检测

Stage3 是 SVTuner 的核心模块。它不做最终映射，而是判断单细胞参考中每个细胞类型在 ST 样本中是否有足够支持。

常用运行命令：

```powershell
python -m src.stages.stage3_type_plugin --sample <sample> --sc_expr_source normalized
```

主要输出：

```text
stage3_typematch/type_support.csv
stage3_typematch/sc_metadata_with_plugin_type.csv 或等价 metadata 输出
```

`type_support.csv` 中最关键的列包括：

- `orig_type`：原始单细胞类型。
- `n_cells`：该类型在参考中的细胞数量。
- `support_score`：该类型在 ST 中的支持分数。
- `support_category`：强支持、弱支持或不支持等类别。
- `auto_missing`：是否被自动判定为疑似缺失。
- `auto_missing_confirmed`：是否最终确认缺失。
- `auto_missing_confirmation_reason`：确认或拒绝缺失的原因。
- `masked_missing_candidate`：是否像被掩盖的缺失类型。
- `marker_identity_candidate`：marker 证据指向的身份候选。
- `marker_identity_z`：marker 证据 z 分数。
- `Action`：最终操作，例如 keep、drop、mark unknown 或 relabel。

Stage3 的主要判断逻辑包括：

- 根据单细胞表达选择 marker genes。
- 过滤非特异 marker，避免弱特异性基因误导类型判断。
- 比较 marker 在 ST 中的表达支持情况。
- 结合 `support_score`、marker identity 和缺失候选状态判断是否 missing。
- 对强支持类型提供保护，避免把真实存在的类型误判为缺失。
- 对弱支持类型保留拯救空间，而不是一律硬过滤。
- 为下游生成 `plugin_type`，让 Stage4 可以使用协调后的细胞池。

因此 Stage3 的定位是“映射前的类型准入控制”，不是最终空间分配器。

### 2.3 Stage4：Baseline 与 Route2 映射

Stage4 负责实际运行映射。当前主线使用 CytoSPACE，并固定生成两套结果。

CytoSPACE baseline：

```powershell
python -m src.stages.stage4_cytospace `
  --sample <sample> `
  --project_root . `
  --missing_type "__AUTO__" `
  --n_processors 1 `
  --n_subspots 800 `
  --mapping_cells_per_spot <2-or-5> `
  --sc_expr_source normalized `
  --filter_mode none `
  --cell_type_column sc_meta `
  --filter_scope unsupported_all `
  --stage4_suffix _baseline
```

SVTuner route2：

```powershell
python -m src.stages.stage4_cytospace `
  --sample <sample> `
  --project_root . `
  --missing_type "__AUTO__" `
  --n_processors 1 `
  --n_subspots 800 `
  --mapping_cells_per_spot <2-or-5> `
  --sc_expr_source normalized `
  --filter_mode plugin_unknown `
  --cell_type_column plugin_type `
  --filter_scope missing_only `
  --stage4_suffix _route2
```

主要输出：

```text
result/<sample>/stage4_cytospace_baseline/cytospace_output/cell_assignment.csv
result/<sample>/stage4_cytospace_route2/cytospace_output/cell_assignment.csv
result/<sample>/stage4_cytospace_*/cytospace_output/fractional_abundances_by_spot.csv
result/<sample>/stage4_cytospace_*/stage4_summary.json
```

`stage4_summary.json` 用来记录：

- 映射前后细胞数量。
- 过滤掉的细胞数量。
- 缺失类型是否仍出现在最终 assignment 中。
- Stage3 是否检测到缺失类型。
- 是否存在非缺失类型被误标 unknown。
- route2 是否恢复了非缺失 unknown 细胞。

### 2.4 评估与可视化

模拟数据有真值，因此可以计算映射质量指标。真实数据没有显式真值，因此主要依靠 marker signature、空间表达模式和映射叠加图进行解释。

当前主方法对比指标是 `composition_recovery`，它比较每个 spot 的预测细胞类型组成和模拟真值组成：

```text
composition_recovery = sum(min(predicted_fraction, truth_fraction))，再按 truth mass 加权
```

无噪声方法对比图：

```powershell
python scripts\plot_method_comparison_composition_recovery.py `
  --project_root . `
  --out_dir visualizations\method_comparison\no_noise `
  --output_prefix composition_recovery_7mapping_methods_no_noise
```

10% scRNA reference noise 方法对比图：

```powershell
python scripts\plot_method_comparison_composition_recovery.py `
  --project_root . `
  --out_dir visualizations\method_comparison\scnoise10 `
  --output_prefix composition_recovery_7mapping_methods_scnoise10 `
  --sample_suffix _scnoise10 `
  --title "Spatial cell-type composition recovery, 10% sc noise"
```

## 3. 模拟数据构造

模拟数据不是随机生成的点，而是在真实 ST 空间结构基础上构造可控缺失类型场景。核心原则是：保留真实空间坐标、尽量保留细胞类型聚集结构，并生成可用于评估的真值文件。

当前模拟数据组包括：

```text
real_brca
human_lung_5loc
mouse_brain_refined
```

### 3.1 构造原则

模拟数据构造遵循以下逻辑：

1. 使用真实或真实风格的 ST 坐标作为空间骨架。
2. 构造 spot 级细胞类型组成真值。
3. 为每个 query cell 记录真实 spot 或类型来源。
4. 在指定场景中移除目标缺失类型。
5. 使用替代类型进行填补或重平衡，保证 ST 表达矩阵仍可运行。
6. 保存 truth 文件，用于评价 baseline、route2 和其他映射方法是否把缺失类型错误映射回来。

关键真值文件：

```text
data/sim/<group>/<sample>/sim_truth_spot_type_fraction.csv
data/sim/<group>/<sample>/sim_truth_query_cell_spot.csv
data/sim/<group>/<sample>/sim_info.json
```

### 3.2 当前 12 个无噪声模拟场景

```text
real_brca_clustered_sim
real_brca_clustered_sim_missing_epithelial_cells
real_brca_clustered_sim_missing_epithelial_monocytes_macrophages
real_brca_clustered_sim_missing_epithelial_monocytes_endothelial
real_brca_clustered_sim_missing_epithelial_monocytes_endothelial_fibroblasts

human_lung_5loc_fine9_clustered_sim
human_lung_5loc_fine9_clustered_sim_missing_ciliated
human_lung_5loc_fine9_clustered_sim_missing_ciliated_endothelia_vascular

mouse_brain_refined8_balanced_clustered_sim
mouse_brain_refined8_balanced_clustered_sim_missing_micro_fill_ext_l56
mouse_brain_refined8_balanced_clustered_sim_missing_micro_astro_ctx_fill_ext_l56
mouse_brain_refined8_balanced_clustered_sim_missing_micro_astro_ctx_oligo_2_fill_ext_l56
```

`real_brca` 主要检验乳腺癌场景中的 epithelial、immune、endothelial 和 fibroblast 相关缺失类型。`human_lung_5loc` 主要检验肺组织细粒度类型，例如 `Ciliated` 和 `Endothelia_vascular`。`mouse_brain_refined` 主要检验脑组织细粒度类型，例如 `Micro`、`Astro_CTX` 和 `Oligo_2`。

## 4. 真实数据 profile-mask 掩盖实验

真实数据没有显式真值。为了在真实组织结构中构造可解释的缺失类型压力测试，项目建立了 profile-mask 场景：保留真实 ST 坐标和表达背景，同时对指定目标类型进行掩盖，使 Stage3 和 route2 可以检验该类型是否被识别为不应参与映射。

当前真实 profile-mask 场景包括：

```text
adult_mouse_kidney_real_profile_mask_endo
ffpe_mouse_brain_sagittal_real_profile_mask_astrocyte
human_breast_cancer_real_profile_mask_basal_cell
human_breast_cancer_visium_ff_wta_real_profile_mask_macrophage
human_breast_cancer_wta_120_real_profile_mask_endothelial_cell
human_cervical_cancer_real_profile_mask_epithelial_cell
human_heart_ff_real_profile_mask_endothelial_cell
human_intestine_cancer_real_profile_mask_endothelial_cell
human_lymph_node_real_profile_mask_b_cell
mouse_embryo_real_profile_mask_erythroid
```

真实掩盖场景可视化脚本：

```powershell
python scripts\render_real_profile_mask_visualization.py --project_root . --sample <profile_mask_sample>
```

输出目录：

```text
visualizations/masked_scenarios/real/<sample>/
```

当前真实数据掩盖图采用四面板标准：

- A：原始真实 ST 的目标 marker signature score。
- B：profile-mask 后真实 ST 的目标 marker signature score。
- C：baseline 映射后目标类型叠加到 ST 表达信号上的空间图。
- D：route2 映射后目标类型叠加到 ST 表达信号上的空间图。

可视化标准已经统一为深紫色表达背景和青色空心轮廓叠加。每个面板底部统计标注和 D 面板额外文字标注已经删除，以保证图面简洁。

## 5. 10% scRNA Reference Noise 场景

噪声实验用于检验 Stage3 和各映射方法在参考表达被扰动时是否仍稳定。

噪声构造方式：

- 只扰动 scRNA reference 表达矩阵。
- ST 表达矩阵、空间坐标和 truth 文件保持不变。
- 对 sc expression 中约 10% 的元素进行扰动。
- 每个表达矩阵独立生成 mask，避免不同矩阵维度不一致导致错误。

生成脚本：

```powershell
python scripts\generate_sc_noise_from_processed_sim.py `
  --project_root . `
  --sim_group <group> `
  --source_sample <source_sample> `
  --target_sample <source_sample>_scnoise10 `
  --noise_fraction 0.10 `
  --seed 42 `
  --overwrite
```

生成文件包括：

```text
data/sim/<group>/<sample>_scnoise10
data/processed/simulation_experiments/<group>/<sample>_scnoise10
configs/datasets/<sample>_scnoise10.yaml
result/<sample>_scnoise10/sc_noise_generation_summary.json
```

当前 12 个 `_scnoise10` 场景已经完成 7 方法运行。需要如实记录的 Stage3 诊断现象是：

- 最终 route2 映射层面，缺失类型没有残留在最终 assignment 中。
- 如果严格按 Stage3 自动检测逻辑，10% noise 下有两个诊断问题：一个 human lung 双缺失场景漏检 `Endothelia_vascular`，一个 no-missing mouse brain 场景误标 `Inh_Sst`。
- route2 最终通过恢复或 truth-aware 过滤避免了实际错误过滤，但 Stage3 诊断结果需要在后续论文表述中区分清楚。

## 6. 多方法对照实验

当前正式对照实验包含 7 种方法：

1. `CytoSPACE`
2. `SVTuner + CytoSPACE`
3. `Tangram (all genes)`
4. `Tangram (marker genes)`
5. `novoSpaRc`
6. `SpaOTsc`
7. `CellTrek`

早期使用过的 Pearson correlation 和 Euclidean distance 已经从正式对照中移除，因为它们更像简单相似度 baseline，而不是完整空间映射方法。对应脚本、日志和结果目录已经删除。

### 6.1 Tangram

Marker genes 版本：

```powershell
python scripts\run_tangram_marker_mapping.py `
  --project_root . `
  --group <group> `
  --sample <sample> `
  --top_n_marker 50 `
  --num_epochs 200 `
  --device cpu
```

All genes 版本：

```powershell
python scripts\run_tangram_marker_mapping.py `
  --project_root . `
  --group <group> `
  --sample <sample> `
  --gene_mode all `
  --num_epochs 200 `
  --device cpu
```

输出目录：

```text
result/<sample>/stage4_mapping/tangram_marker/
result/<sample>/stage4_mapping/tangram_all/
```

### 6.2 CellTrek

当前项目使用 Python runner 生成标准化 CellTrek-style 对照结果：

```powershell
python scripts\run_celltrek_mapping.py `
  --project_root . `
  --group <group> `
  --sample <sample> `
  --max_genes 2000 `
  --n_pcs 30 `
  --ntree 500
```

输出目录：

```text
result/<sample>/stage4_mapping/celltrek/
```

### 6.3 novoSpaRc

```powershell
python scripts\run_novosparc_mapping.py `
  --project_root . `
  --group <group> `
  --sample <sample> `
  --max_genes 500 `
  --n_pcs 30 `
  --alpha_linear 0.5 `
  --epsilon 0.005
```

输出目录：

```text
result/<sample>/stage4_mapping/novosparc/
```

### 6.4 SpaOTsc

```powershell
python scripts\run_spaotsc_mapping.py `
  --project_root . `
  --group <group> `
  --sample <sample> `
  --max_genes 500 `
  --n_pcs 30 `
  --alpha 0.1 `
  --epsilon 0.1 `
  --rho inf `
  --niter 10
```

输出目录：

```text
result/<sample>/stage4_mapping/spaotsc/
```

### 6.5 标准化输出

所有外部方法 runner 都尽量输出统一文件：

```text
cell_assignment.csv
spot_type_fraction.csv
metrics_simulation.json
run.log
```

这样后续评估脚本可以用同一套逻辑读取不同方法的结果。

## 7. 当前方法对比结果

无噪声正式图：

```text
visualizations/method_comparison/no_noise/composition_recovery_7mapping_methods_no_noise_boxplot.png
```

10% sc-noise 正式图：

```text
visualizations/method_comparison/scnoise10/composition_recovery_7mapping_methods_scnoise10_boxplot.png
```

无噪声均值排序：

1. `SVTuner + CytoSPACE`: 0.6536
2. `novoSpaRc`: 0.6226
3. `Tangram (marker genes)`: 0.5935
4. `CytoSPACE`: 0.5848
5. `Tangram (all genes)`: 0.5674
6. `CellTrek`: 0.4988
7. `SpaOTsc`: 0.4843

10% sc-noise 均值排序：

1. `SVTuner + CytoSPACE`: 0.6554
2. `novoSpaRc`: 0.6212
3. `Tangram (marker genes)`: 0.6000
4. `CytoSPACE`: 0.5882
5. `Tangram (all genes)`: 0.5671
6. `CellTrek`: 0.4861
7. `SpaOTsc`: 0.4807

这些数值来自 12 个模拟场景的 scenario-level composition recovery。

## 8. 当前可视化资产

主要可视化目录：

```text
visualizations/masked_scenarios/real/
visualizations/simulations/real_brca/
visualizations/simulations/human_lung_5loc/
visualizations/simulations/mouse_brain_refined/
visualizations/method_comparison/no_noise/
visualizations/method_comparison/scnoise10/
```

模拟数据可视化通常是三联图：

- truth 空间分布。
- CytoSPACE baseline 映射。
- SVTuner + CytoSPACE route2 映射。

真实 profile-mask 可视化是四面板图。方法对比可视化是基于 12 个场景的箱线图。

## 9. 环境和运行注意事项

主环境：

```powershell
conda activate cytospace_v1.1.0_py310
```

常用起始命令：

```powershell
Set-Location "E:\AAA文件\Experiment\SVTuner\sctuner2.0"
$py = "E:\ANACONDA\envs\cytospace_v1.1.0_py310\python.exe"
```

运行 CytoSPACE 主线时，通常建议禁用 user site：

```powershell
$env:PYTHONNOUSERSITE = "1"
$env:CYTOSPACE_SKIP_ASSIGNED_EXPRESSION = "1"
```

运行 Tangram、novoSpaRc 和 SpaOTsc 时，不能启用 `PYTHONNOUSERSITE=1`，因为这些包安装在 user site：

```powershell
Remove-Item Env:PYTHONNOUSERSITE -ErrorAction SilentlyContinue
```

减少线程争用：

```powershell
$env:OMP_NUM_THREADS = "1"
$env:MKL_NUM_THREADS = "1"
$env:NUMBA_CACHE_DIR = Join-Path $PWD ".numba_cache"
```

## 10. 仓库和存储策略

仓库追踪代码、配置、文档和精选可视化。大型运行数据和本地中间产物不进入 Git。

`.gitignore` 当前忽略：

```text
data/
logs/
result/
.numba_cache/
tmp_*.pdf
tmp_*.txt
external/CellTrek/
```

这意味着仓库保留可复现流程和正式图表，而原始大矩阵、每次运行的详细结果和本地缓存留在本机。

## 11. 当前项目状态

当前检查点已经完成：

- 真实 profile-mask 掩盖场景可视化。
- `real_brca`、`human_lung_5loc`、`mouse_brain_refined` 模拟场景可视化。
- 无噪声 7 方法对照实验。
- 10% scRNA reference noise 7 方法对照实验。
- Pearson correlation 和 Euclidean distance 从正式对照中移除。
- 当前正式方法集固定为 CytoSPACE、SVTuner + CytoSPACE、Tangram all genes、Tangram marker genes、novoSpaRc、SpaOTsc 和 CellTrek。

下一步可以进入论文级结果组织：整理主图、补充图、方法表、噪声鲁棒性表述和 Stage3 诊断边界。