# SVTuner 项目说明（2025）

SVTuner：面向空间变异性（SVG）与类型感知的细胞–空间位点映射调谐框架。当前仓库已接入 GitHub，Stage0 环境自检、Stage1 预处理、Stage3/4/5/6 脚本已可运行（Stage2 已移除）；后续规划见 SVTuner.md。

## 目录结构（关键路径）
- `configs/`：全局 & 数据集配置
  - `project_config.yaml`：全局默认
  - `datasets/<sample>.yaml`：样本级配置（推荐在此覆盖 QC/路径）
  - `env_archive/`：历史环境 yml 备份
- `data/`
  - `raw/<sample>/`：原始输入（scRNA/ST）
  - `processed/<sample>/stage1_preprocess/`：Stage1 输出（sc/st processed、common_genes、HVG）
- `result/<sample>/stage1_preprocess/`：Stage1 摘要与 QC 图
- `r_scripts/`：Stage1 预处理（Seurat）
- `src/`
  - `config.py`：默认参数与路径
  - `sim/`：SimGen 场景构建
  - `stages/`：阶段脚本（Stage0/3/4/5/6）
  - `utils/`：通用工具
- `external/`：第三方源码（如 cytospace 副本）

## 环境要求
- Python 3.10（已在 `cytospace_v1.1.0_py310` 环境验证）
- R 4.5.x + Seurat 5.3.x（推荐通过 conda 安装，见下）
- Git LFS（已用于跟踪 data/raw 大文件）

### R/Seurat 安装（推荐）
```pwsh
conda activate cytospace_v1.1.0_py310
conda install -c conda-forge r-base r-seurat r-data.table r-matrix r-jsonlite -y
```

## 数据输入（real_brca 示例）
放置于 `data/raw/real_brca/`，默认文件名（可在 YAML 配置中覆盖）：
- `real_brca_scRNA_GEP.txt`：scRNA 表达，首列 Gene，后续列为细胞
- `real_brca_scRNA_celllabels.txt`：两列 `SpotID / CellType`
- `real_brca_STdata_GEP.txt`：ST 表达，首列 Gene，后续列为 spot
- `real_brca_STdata_coordinates.txt`：三列 `SpotID / row / col`

## 配置说明
- `configs/project_config.yaml`：全局默认（可为空）
- `configs/datasets/<sample>.yaml`：推荐在此覆盖 QC/路径，例如：
```yaml
paths:
  sc_expr: real_brca_scRNA_GEP.txt
  sc_meta: real_brca_scRNA_celllabels.txt
  st_expr: real_brca_STdata_GEP.txt
  st_meta: real_brca_STdata_coordinates.txt
  svg_marker_whitelist: configs/svg_whitelist.txt  # 可选
qc:
  sc_min_genes: 200
  sc_max_genes: 6000
  sc_max_mt: 10
  st_min_genes: 100
  st_max_genes: Inf
  st_max_mt: 20
  hvg_nfeatures: 2000
  mt_pattern: "^MT-"
gene_filter:
  min_cells_sc: 0
  min_cells_st: 0
```
注：未指定时，脚本会尝试用 `<sample>_...` 前缀自动匹配；否则报错。

## 运行方式（当前以分阶段脚本为主）
- Stage0 环境自检（可选）：
```pwsh
python src/stages/stage0_envcheck.py
```
- Stage1 预处理：使用 `r_scripts/stage1_preprocess.R`（见下方示例）
- Stage3/4/5/6：分别运行对应脚本（见下方“脚本运行方式”清单）

## Stage1 预处理要点（r_scripts/stage1_preprocess.R）
- 读取配置：`configs/project_config.yaml`、`configs/datasets/<sample>.yaml`，支持覆盖 QC/路径。
- 温和过滤（可调）：`qc` 与 `gene_filter` 参数，支持 `mt_pattern`，支持 SVG/marker 白名单强制保留。
- 输出：
  - `data/processed/<sample>/stage1_preprocess/`：
    - `common_genes.txt`, `hvg_genes.txt`
    - 若加 `--export_csv`：`exported/` 下生成 `sc_expression_normalized.csv`、`st_expression_normalized.csv`、`sc_metadata.csv`、`st_coordinates.csv`（Stage3/4/5/6 直接复用）
  - `result/<sample>/stage1_preprocess/`：
    - `stage1_summary.json`
- Summary 包含：过滤前后细胞/基因数、whitelist 保留数、公共基因数与路径、HVG 请求/实际数、mt% 统计等。
- 运行（独立调用示例，避免 MinGW 冲突，先激活 env 再清理 PATH）：  
```pwsh
$env:PATH="E:\ANACONDA\envs\cytospace_v1.1.0_py310\Library\bin;E:\ANACONDA\envs\cytospace_v1.1.0_py310\bin;E:\ANACONDA\envs\cytospace_v1.1.0_py310\Scripts;C:\Windows\System32"
& "E:\ANACONDA\envs\cytospace_v1.1.0_py310\Scripts\Rscript.exe" r_scripts/stage1_preprocess.R --sample real_brca --project_root D:/Experiment/SVTuner_Project --export_csv
```
  也可用 `conda run -n cytospace_v1.1.0_py310 Rscript ...`（若无 DLL 冲突）。

## Stage3 类型不匹配 / Unknown-aware 插件（src/stages/stage3_type_plugin.py）
- 配置/参数（CLI 可覆盖）：`strong_th`、`weak_th`、`st_cluster_k`、`unknown_floor`、`min_cells_rare_type`、`eps`，相似度度量固定为 `weighted_cosine`；support_score 定义固定为 top-3 相似度均值。
- 输入：
  - Stage1 导出：`data/processed/<sample>/stage1_preprocess/exported/` 下的 `sc_expression_normalized.csv`、`st_expression_normalized.csv`、`sc_metadata.csv`（类型列为 `cell_type` 或 `celltype`）、`st_coordinates.csv`
  - 插件基因：默认使用 `data/processed/<sample>/stage1_preprocess/hvg_genes.txt`，或在 `configs/datasets/<sample>.yaml` 的 `stage3.plugin_genes_path` 中显式指定
  - 可选基因权重：在 `stage3.gene_weights_path` 指向权重文件（若存在会用于加权相似度）
- 输出：
  - `data/processed/<sample>/stage3_typematch/`：`type_support.csv`、`cell_type_relabel.csv`（Unknown 统一为 `Unknown_sc_only`）、`type_prior_matrix.csv`
  - `result/<sample>/stage3_typematch/`：`stage3_summary.json`（记录参数、支持度分布、Unknown 占比、稀有类型、plugin_type 列表）
- 运行示例（已激活 env，复用导出的 CSV）：
```pwsh
cd D:\Experiment\SVTuner_Project
python src/stages/stage3_type_plugin.py --sample real_brca
# 可调阈值/聚类等：--strong_th --weak_th --st_cluster_k --unknown_floor --min_cells_rare_type
```

## SimGen 场景构建（src/sim/*.py）
- 作用：生成模拟场景（S0/S1/M0），并同步生成 Stage1 export 所需 CSV。
- 配置：`configs/simgen/*.yaml`
- 输出：
  - `data/sim/<sample>/<scenario>/<seed_*>/`：sc/st 模拟数据与 truth 文件
  - `data/processed/<scenario>/stage1_preprocess/exported/`：Stage1 export CSV（Stage3/4/5 直接使用）

## Stage4 映射（src/stages/stage4_cytospace.py）
- Route2：基于 Stage3 的 `plugin_type` 过滤映射
- Baseline：不经过过滤，等价官方 CytoSPACE

## Stage5 评估（src/stages/stage5_route2_s0.py）
- 产出 composition/coverage/accuracy 等指标，并支持 baseline 对比

## Stage6 真实数据审计（src/stages/stage6_real_audit.py）
- 对 baseline 与 route2_v5_3_1 的映射输出进行一致性/分布/marker 审计，产出结构化 JSON
- 运行示例：
```pwsh
python src/stages/stage6_real_audit.py --sample real_brca --baseline_id baseline --route2_id route2_v5_3_1
```

## Git/LFS 注意
- 大文件（raw 数据）已由 Git LFS 跟踪；请确保本地安装并启用 Git LFS。
- 默认 `.gitignore` 忽略 `result/`，`data/processed/` 目前未忽略，如需忽略可自行调整。

## 后续规划
- Stage7 可视化与报告：详见 `dosc/Stage7 可视化与报告模块设计.md` 与 `dosc/SVTuner.md`。

## 常见问题
- Rscript 找不到：确认 `--r-cmd` 指向有效 Rscript（建议 conda run）。
- DLL 相关报错（Windows）：使用默认 `--r-cmd "conda run -n cytospace_v1.1.0_py310 Rscript"`，或手动设置 `CONDA_DLL_SEARCH_MODIFICATION_ENABLE=1`。
- 输入文件名不匹配：在 `configs/datasets/<sample>.yaml` 显式配置 `paths`。

## 脚本运行方式（逐文件清单）
说明：以下列出了当前仓库内所有可执行脚本的运行方式。
- 有 CLI 参数的脚本：先运行 `python <script> --help` 查看完整参数。
- 无 CLI 参数的脚本：需在脚本顶部修改 sample/path 常量后运行。

### 主入口 / Stage 脚本
- `src/stages/stage0_envcheck.py`：环境自检
```pwsh
python src/stages/stage0_envcheck.py
```
- `r_scripts/stage1_preprocess.R`：Stage1 R 预处理
```pwsh
$env:PATH="E:\ANACONDA\envs\cytospace_v1.1.0_py310\Library\bin;E:\ANACONDA\envs\cytospace_v1.1.0_py310\bin;E:\ANACONDA\envs\cytospace_v1.1.0_py310\Scripts;C:\Windows\System32"
& "E:\ANACONDA\envs\cytospace_v1.1.0_py310\Scripts\Rscript.exe" r_scripts/stage1_preprocess.R --sample real_brca_simS0_seed0 --project_root D:/Experiment/SVTuner_Project --export_csv
```
- `src/stages/stage3_type_plugin.py`：Stage3 V4/V5 unknown-aware
```pwsh
python src/stages/stage3_type_plugin.py --sample real_brca_simS0_seed0
# 指定数据配置：--dataset_config configs/datasets/real_brca_simS0_seed0.yaml
```
- `src/stages/stage4_cytospace.py`：Stage4 Route2/Baseline
```pwsh
# Route2 (plugin_unknown)
python src/stages/stage4_cytospace.py --sample real_brca_simS0_seed0 --filter_mode plugin_unknown --cell_type_column plugin_type --missing_type "T cells CD8" --mean_cell_numbers 5 --solver_method lap_CSPR --seed 0 --sampling_sub_spots --n_subspots 800 --n_processors 1 --filter_scope unsupported_all

# Baseline (official CytoSPACE)
python src/stages/stage4_cytospace.py --sample real_brca_simS0_seed0 --filter_mode none --cell_type_column sc_meta --seed 0 --sampling_sub_spots --n_subspots 800 --n_processors 1
```
- `src/stages/stage4_route2.py`：Stage4 MVP（仅过滤验收）
```pwsh
python src/stages/stage4_route2.py --sample real_brca_simS0_seed0 --missing_type "T cells CD8" --filter_mode plugin_unknown
```
- `src/stages/stage5_route2_s0.py`：Stage5 评估 / 对比
```pwsh
python src/stages/stage5_route2_s0.py --sample real_brca_simS0_seed0 --run_tag route2_v5_3_1 --stage4_dir result/real_brca_simS0_seed0/stage4_cytospace/cytospace_output --sim_dir data/sim/real_brca/S0/t_cells_cd8/seed_0 --out_dir result/real_brca_simS0_seed0/stage5_route2_s0/route2_v5_3_1

# 对比 json（baseline vs route2）
python src/stages/stage5_route2_s0.py --compare_baseline <baseline.json> --compare_route2 <route2.json> --compare_out <compare.json>
```
- `src/stages/stage6_real_audit.py`：Stage6 真实数据审计
```pwsh
python src/stages/stage6_real_audit.py --sample real_brca --baseline_id baseline --route2_id route2_v5_3_1
```
- `src/stages/stage7_report.py`：Stage7 报告生成
```pwsh
python src/stages/stage7_report.py --sample real_brca --baseline_id baseline --route2_id route2_v5_3_1
```

### SimGen 脚本（src/sim/*.py）
- `src/sim/simgen_s0.py`
```pwsh
python src/sim/simgen_s0.py --sample real_brca --sim_config configs/simgen/s0.yaml --seed 0 --missing_type "T cells CD8"
```
- `src/sim/simgen_s1.py`
```pwsh
python src/sim/simgen_s1.py --sample real_brca --sim_config configs/simgen/s1.yaml
```
- `src/sim/simgen_m0.py`
```pwsh
python src/sim/simgen_m0.py --sample real_brca --sim_config configs/simgen/m0.yaml --seed 42 --missing_type "T cells CD8" --dispersed_type "T cells CD4" --scenario_tag cd4_10pct_cd8_3pct_clustered
```

### scripts/ 辅助工具
- `scripts/auto_tune_support_threshold.py`
```pwsh
python scripts/auto_tune_support_threshold.py --sample real_brca_simS0_seed42 --missing_type "T cells CD8" --seed 42 --sim_dir data/sim/real_brca/S0/t_cells_cd8/seed_42
```
- `scripts/auto_tune_weak_th.py`
```pwsh
python scripts/auto_tune_weak_th.py --sample real_brca_simS0_mt_b_cells_seed42 --missing_type "B cells" --seed 42 --sim_dir data/sim/real_brca/S0/b_cells/seed_42
```
- `scripts/prepare_sample_for_tuning.py`
```pwsh
python scripts/prepare_sample_for_tuning.py --sample real_brca_simS0_seed42 --missing_type "T cells CD8" --seed 42
```
- `scripts/run_single_weak_th.py`
```pwsh
python scripts/run_single_weak_th.py --sample real_brca_simS0_seed42 --missing_type "T cells CD8" --weak_th 0.55 --seed 42
```
- `scripts/run_until_cleared.py`
```pwsh
python scripts/run_until_cleared.py --sample real_brca_simS0_seed42 --missing_type "T cells CD8" --seed 42 --start_weak_th 0.50 --max_weak_th 0.70 --step 0.02
```
- `scripts/transpose_sample.py`
```pwsh
python scripts/transpose_sample.py --sample real_brca_simS0_seed42
```
- `scripts/transpose_seed.py`
```pwsh
python scripts/transpose_seed.py --seed 42
```
- `scripts/check_st_expression.py`（传入 3 个路径）
```pwsh
python scripts/check_st_expression.py <st_expression_normalized.csv> <sc_expression_normalized.csv> <sc_metadata.csv>
```
- `scripts/check_stage3_cd4.py`
```pwsh
python scripts/check_stage3_cd4.py real_brca_simS0_seed42
```
- `scripts/check_truth_types.py`
```pwsh
python scripts/check_truth_types.py data/sim/real_brca/S0/t_cells_cd8/seed_42/sim_truth_spot_type_fraction.csv
```
- `scripts/filter_audit.py`
```pwsh
python scripts/filter_audit.py real_brca_simS0_seed42 "T cells CD8"
```
- `scripts/verify_filter_audit.py`
```pwsh
python scripts/verify_filter_audit.py result/real_brca_simS0_seed42/stage5_route2_s0/stage5_route2_s0__route2_v5_3_1.json
```
- `scripts/aggregate_seed_summary.py`（脚本内 SEEDS/base_tag 可修改）
```pwsh
python scripts/aggregate_seed_summary.py
```
- `scripts/analyze_false_positives.py`（脚本内 sample/missing_type/threshold 可修改）
```pwsh
python scripts/analyze_false_positives.py
```
- `scripts/audit_sha1.py`（脚本内路径可修改）
```pwsh
python scripts/audit_sha1.py
```
- `scripts/check_false_positives.py`（脚本内样本可修改）
```pwsh
python scripts/check_false_positives.py
```
- `scripts/compare_metrics_v41.py`（脚本内路径可修改）
```pwsh
python scripts/compare_metrics_v41.py
```
- `scripts/inspect_sc.py`（脚本内路径可修改）
```pwsh
python scripts/inspect_sc.py
```
- `scripts/print_metrics.py`（脚本内路径可修改）
```pwsh
python scripts/print_metrics.py
```
- `scripts/test_v42_pvalue.py`（脚本内参数可修改）
```pwsh
python scripts/test_v42_pvalue.py
```
- `scripts/transpose_seed0.py`（脚本内路径可修改）
```pwsh
python scripts/transpose_seed0.py
```
- `scripts/type_filter_audit_report.py`（脚本内路径可修改）
```pwsh
python scripts/type_filter_audit_report.py
```

### tools/ 工具
- `tools/check_stage4_stats.py`（脚本内 sample 可修改）
```pwsh
python tools/check_stage4_stats.py
```
