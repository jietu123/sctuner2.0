# SVTuner 场景运行指南

本文档说明当前仓库中 Stage 1、Stage 3A、Stage 3B、Stage 4，以及 C1–C8 审稿补充实验的实际入口和文件约定。默认原则是：**不覆盖历史输出、不隐式清理目录、不把本地大数据提交到 Git。**

## 1. 运行前准备

在项目根目录执行：

~~~powershell
conda env create -f configs/environment.yml
conda activate cytospace_v1.1.0_py310
pip install -e .
svtuner version
svtuner envcheck
~~~

建议显式指定项目根目录：

~~~powershell
$ProjectRoot = "E:\AAA文件\Experiment\SVTuner\sctuner2.0"
Set-Location -LiteralPath $ProjectRoot
~~~

本文档不提供默认递归删除命令。需要重跑时，应先确认精确的 sample、stage 和 suffix，再单独处理目标目录。

## 2. 配置与目录解析

每个 sample 对应：

~~~text
configs/datasets/<sample>.yaml
~~~

如果 YAML 定义了 storage.group：

~~~text
data/raw/<group>/<sample>/
data/processed/<group>/<sample>/
result/<group>/<sample>/
~~~

否则：

~~~text
data/raw/<sample>/
data/processed/<sample>/
result/<sample>/
~~~

Stage 1、Stage 3A、Stage 3B 使用上述 group-aware 路径。当前 Stage 4 的输出例外地写入：

~~~text
result/<sample>/stage4_cytospace<suffix>/
~~~

运行前至少核对：

- dataset YAML 中的 sample 与 storage.group；
- SC/ST 表达和 metadata 文件名；
- cell-type 列名；
- Stage3A/Stage3B 参数；
- 输出 suffix 是否会覆盖历史目录。

## 3. Stage 1：预处理

正式 R 入口：

~~~powershell
Rscript r_scripts/stage1_preprocess.R --sample <sample> --project_root $ProjectRoot --export_csv
~~~

典型输出：

~~~text
stage1_preprocess/exported/sc_expression_normalized.csv
stage1_preprocess/exported/sc_expression_data.csv
stage1_preprocess/exported/sc_expression_counts.csv
stage1_preprocess/exported/sc_metadata.csv
stage1_preprocess/exported/st_expression_normalized.csv
stage1_preprocess/exported/st_coordinates.csv
stage1_preprocess/hvg_genes.txt
stage1_preprocess/stage1_summary.json
~~~

这些矩阵通常体积很大，位于 Git 忽略的 data/processed/ 下。不要为了提交 GitHub 而复制或强制加入。

## 4. Stage 3A：reference-level diagnostics

推荐显式指定 SC expression source，并用 suffix 隔离测试：

~~~powershell
python -m src.stages.stage3_type_plugin --project_root $ProjectRoot --sample <sample> --sc_expr_source normalized --output_suffix _trial
~~~

主要输出目录：

~~~text
stage3_typematch<suffix>/
~~~

关键文件：

- stage3_summary.json：最终 missing/unsupported types、参数与诊断摘要；
- type_support.csv：每个 cell type 的 support score、category、n_cells 等；
- cell_type_relabel.csv：原始与调整后标签；
- stage3_adjusted_annotations.csv：供后续 mapping 使用的 annotations。

正式历史目录 stage3_typematch 没有 suffix。参数敏感性或诊断运行必须使用独立 suffix，不能覆盖正式目录。

## 5. Stage 3B：spatial-unit diagnostics

若希望沿用 dataset YAML 与程序默认值，使用模块入口且不要额外覆盖数值参数：

~~~powershell
python -m src.stages.stage3b_st_unsupported --project_root $ProjectRoot --sample <sample> --output_suffix _trial
~~~

说明：

- 模块入口会解析 dataset YAML 和内部默认值；
- svtuner stage3b 会传递 CLI 中定义的显式默认值；
- 精确复现已有结果时，应先确认历史运行采用哪种入口；
- standalone 与 sequential 两种 Stage3B 必须使用独立输出目录和独立 mask。

关键输出通常包括：

- spot_unsupported_scores.csv：spot-level observed statistics、p/q values、score 和最终 withheld decision；
- Stage3B summary/region 文件：候选与最终空间区域；
- 用于 Stage 4 blanking 的 final mask。

## 6. Stage 4：四条 decomposition 路线

四条路线的科学定义如下：

| 路线 | SC reference | Stage3B mask |
|---|---|---|
| Baseline | 完整 reference | 无 |
| Stage3A-only | Stage3A retained reference | 无 |
| Stage3B-only | 完整 reference | standalone Stage3B mask |
| Full | Stage3A retained reference | sequential Stage3B mask |

### Baseline

~~~powershell
python -m src.stages.stage4_cytospace --project_root $ProjectRoot --sample <sample> --filter_scope none --sc_expr_source normalized --stage4_suffix _baseline
~~~

### Stage3A-only

~~~powershell
python -m src.stages.stage4_cytospace --project_root $ProjectRoot --sample <sample> --filter_scope unsupported_all --sc_expr_source normalized --stage4_suffix _stage3a_only
~~~

如果历史路线使用 plugin_unknown、plugin_type 或特殊 missing_type 约定，应严格沿用对应 frozen 命令，而不是根据名称猜测。

### Stage3B-only

使用完整 reference，并通过以下参数传入 **standalone** Stage3B mask：

~~~text
--stage3b_blank_regions
--stage3b_scores_path <standalone spot_unsupported_scores.csv>
~~~

### Full

使用 Stage3A retained reference，并传入基于该 retained reference 重新计算的 **sequential** Stage3B mask。不能把 standalone mask 复用于 Full。

Stage 4 还支持 --stage3_suffix、--stage4_suffix、--filter_scope、--keep_redundant 和 --sc_expr_source。

## 7. 便捷命令的边界

~~~powershell
svtuner run --sample <sample> --project-root $ProjectRoot
~~~

当前 svtuner run 会运行：

1. Stage 1；
2. Stage 3A；
3. baseline Stage 4；
4. Stage3A-controlled Stage 4。

它**不会运行 Stage 3B**。此外，configs/pipeline_presets.yaml 当前为空，因此 --from-scratch 不能在没有补充 preset 的情况下作为通用入口。

## 8. 九个 composite no-noise 场景

| ID | config / sample | ST missing | SC dropout |
|---|---|---|---|
| BC-1 | real_brca7_endothelial_marker_control_sc_missing_endothelial_cells | None | Endothelial cells |
| BC-2 | real_brca7_endothelial_marker_missing_epithelial_cells_sc_missing_endothelial_cells | Epithelial cells | Endothelial cells |
| BC-3 | real_brca7_endothelial_marker_missing_epithelial_cells_pcs_sc_missing_endothelial_cells | Epithelial cells + PCs | Endothelial cells |
| Lung-1 | human_lung_5loc_fine9_clustered_sim_sc_missing_b_cell | None | B cell |
| Lung-2 | human_lung_5loc_fine9_clustered_sim_missing_at2_sc_missing_b_cell | AT2 | B cell |
| Lung-3 | human_lung_5loc_fine9_clustered_sim_missing_at2_fibroblast_sc_missing_b_cell | AT2 + Fibroblast | B cell |
| Brain-1 | mouse_brain_refined7_balanced_clustered_sim_sc_missing_ext_l56 | None | Ext_L56 |
| Brain-2 | mouse_brain_refined7_balanced_clustered_sim_missing_micro_fill_inh_pvalb_sc_missing_ext_l56 | Micro | Ext_L56 |
| Brain-3 | mouse_brain_refined7_balanced_clustered_sim_missing_micro_oligo_2_fill_inh_pvalb_sc_missing_ext_l56 | Micro + Oligo_2 | Ext_L56 |

## 9. C1：multi-seed simulation

C1 使用 simulation seeds 42–51；Stage3A、Stage3B 和 CytoSPACE 的 algorithm seed 固定为 42。runner 每次处理一个 scenario × seed：

~~~powershell
python scripts/run_c1_seed_repeat.py --scenario BC-2 --simulation-seed 43
~~~

实际参数名以脚本 --help 为准；不要用同一 sample/output path 覆盖不同 seed。正式汇总与图位于：

~~~text
visualizations/method_comparison/c1_multiseed/
~~~

本地逐次运行结果位于 result/，不会上传 GitHub。

## 10. C2–C8 入口索引

| 实验 | 当前入口 / 资产 | 注意事项 |
|---|---|---|
| C2 | **CLOSED / CANCELLED** | 不发布正式 sensitivity experiment；临时诊断运行不属于论文证据，也不应作为 C2 结果重建或报告 |
| C3 | scripts/plot_stage3b_threshold_robustness.py | 基于已有 reference-dropout 输出；不重跑 Stage3B |
| C4.1 | scripts/export_c4_1_stage3b_raw_calibration_null.py；scripts/run_c4_1_observed_vs_null_comparison.py | 导出/比较 Stage3B raw calibration null |
| C4.2 | scripts/run_c4_2_technical_perturbation_pilot.py | 研究脚本；只运行 Stage3B |
| C4.3 | scripts/run_c4_3_stage3b_statistical_robustness_pilot.py | BH/BY 与 spatial permutation robustness |
| C5 | scripts/run_c5_decomposition.py | 研究 runner，含本机历史路径假设；迁移机器前先检查，不要盲目批跑 |
| C6 | scripts/generate_real_brca_clustered_sim.py；scripts/generate_st_only_reference_dropout_from_sim.py；scripts/evaluate_c6_independent_source.py | 独立 reference/profile source；生成数据与 result 保持本地 |
| C7 | scripts/repair_c7_metric_provenance.py | 只读冻结输入，确定性修复 CTA 指标 provenance |
| C8 | C8 paired analysis 与 plotting scripts | 推断单位为 independent experiment；不要把内部 readout 当独立重复 |

部分 C4/C5/C8 脚本是冻结研究脚本而非通用产品 CLI，可能包含固定数据集、路径或输出约定。运行前只核对脚本顶部常量和输入路径，不应重构算法来适配本机。

## 11. C6 independent-source 场景

冻结样本：

~~~text
c6_breast_independent_control
c6_breast_independent_fibroblast_dropout
~~~

配置：

~~~text
configs/datasets/c6_breast_independent_control.yaml
configs/datasets/c6_breast_independent_fibroblast_dropout.yaml
~~~

设计边界：

- Dataset A 提供 SC reference、composition truth、坐标与 library size 框架；
- Dataset B 仅提供独立的 type profile；
- generator 的 --profile_source_sample 未提供时保持历史行为；
- dropout utility 只删除 SC reference 中的 Fibroblasts，ST 与 truth 必须保持不变。

## 12. 输出完成性检查

每次正式运行后只做一次最终检查：

1. 命令退出状态为 0；
2. 目标 stage 目录存在；
3. summary JSON/CSV 可以读取；
4. barcode 或 spatial-unit ID 与 coordinates 一一对应；
5. sample、suffix、seed 与 provenance 一致；
6. 没有覆盖无 suffix 的历史结果；
7. 没有把 standalone mask 用于 Full；
8. Git 跟踪边界符合预期。

只读检查示例：

~~~powershell
git check-ignore -v result
git ls-files -- "result/**"
git status --short
~~~

## 13. GitHub 收录边界

.gitignore 明确忽略 data/、result/ 和 logs/。因此：

- data/raw、data/sim、data/processed 不上传；
- result/ 中的本地正式 CSV、JSON、mask、mapping 结果也不上传；
- GitHub 上公开的是代码、configs、reproducibility_release，以及 visualizations/ 中经过筛选的图和 source table；
- 不要使用 git add -f result/ 或 git clean 处理实验数据；
- 需要公开的结果应复制为小型、可再分发、带 provenance 的 curated artifact，再提交到明确的 tracked 目录。

## 14. 安全规则

- 不对项目根目录执行递归删除；
- 不把目录名相似视为 provenance 相同；
- 不修改原始 dataset YAML 来做临时 sweep；
- 不覆盖无 suffix 的 Stage3A/Stage3B/Stage4 正式目录；
- 不因单个 warning 自动调参或改变实验定义；
- 运行长批次时优先本地终端执行，结束后一次性汇总；
- commit 前使用 git status 和 git diff 精确确认要上传的文件。
