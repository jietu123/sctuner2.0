# SVTuner 数据场景运行说明

本文档说明两类数据场景的基本运行方式：

- 模拟数据场景：数据通常位于 `data/sim/<group>/<sample>`，配置位于 `configs/datasets/<sample>.yaml`。
- 真实数据场景：数据通常位于 `data/raw/low_resolution_experiments/<sample>`，配置位于 `configs/datasets/<sample>.yaml`。

以下命令默认在 PowerShell 中运行。

## 0. 通用环境

先进入项目根目录，并避免 Python 用户环境污染：

```powershell
Set-Location "E:\AAA文件\Experiment\SVTuner\sctuner2.0"
$env:PYTHONNOUSERSITE = "1"
```

如果需要显式使用指定 conda 环境中的 Python：

```powershell
$envPrefix = "E:\ANACONDA\envs\cytospace_v1.1.0_py310"
```

## 1. 查看可运行样本

所有可运行样本都应有对应配置文件：

```powershell
Get-ChildItem configs\datasets -Filter *.yaml |
  Select-Object -ExpandProperty BaseName |
  Sort-Object
```

当前重点模拟场景包括：

```text
human_lung_5loc_fine9_clustered_sim
human_lung_5loc_fine9_clustered_sim_missing_ciliated
human_lung_5loc_fine9_clustered_sim_missing_ciliated_endothelia_vascular

mouse_brain_refined8_balanced_clustered_sim
mouse_brain_refined8_balanced_clustered_sim_missing_micro_fill_ext_l56
mouse_brain_refined8_balanced_clustered_sim_missing_micro_astro_ctx_fill_ext_l56
mouse_brain_refined8_balanced_clustered_sim_missing_micro_astro_ctx_oligo_2_fill_ext_l56

real_brca_clustered_sim
real_brca_clustered_sim_missing_epithelial_cells
real_brca_clustered_sim_missing_epithelial_monocytes_endothelial
real_brca_clustered_sim_missing_epithelial_monocytes_endothelial_fibroblasts
real_brca_clustered_sim_missing_epithelial_monocytes_macrophages
```

真实数据场景示例包括：

```text
adult_mouse_kidney_real
human_breast_cancer_real
human_cervical_cancer_real
human_heart_ff_real
human_intestine_cancer_real
human_lymph_node_real
```

## 2. 模拟数据场景运行方式

模拟数据通常已经生成在 `data/sim` 下。运行时建议先从模拟源数据重建 Stage1 中间文件，再运行 Stage3。

### 2.1 Stage1 + Stage3 最小命令

只需要替换 `$sample`：

```powershell
Set-Location "E:\AAA文件\Experiment\SVTuner\sctuner2.0"
$env:PYTHONNOUSERSITE = "1"
$envPrefix = "E:\ANACONDA\envs\cytospace_v1.1.0_py310"

$sample = "<simulation_sample_name>"

Remove-Item "data\processed\simulation_experiments\*\$sample" -Recurse -Force -ErrorAction SilentlyContinue
Remove-Item "result\simulation_experiments\*\$sample" -Recurse -Force -ErrorAction SilentlyContinue

python scripts\prepare_stage1_from_sim_source.py `
  --project_root . `
  --sample $sample

if ($LASTEXITCODE -ne 0) { throw "Stage1 sim export failed: $sample" }

& "$envPrefix\python.exe" -m src.stages.stage3_type_plugin `
  --sample $sample `
  --sc_expr_source normalized

if ($LASTEXITCODE -ne 0) { throw "Stage3 failed: $sample" }
```

### 2.2 查看 Stage3 结果

根据模拟组选择对应路径。下面示例使用 `human_lung_5loc`：

```powershell
$group = "human_lung_5loc"
$sample = "human_lung_5loc_fine9_clustered_sim_missing_ciliated"

Import-Csv "data\processed\simulation_experiments\$group\$sample\stage3_typematch\type_support.csv" |
  Select-Object orig_type,n_cells,support_score,support_category,auto_missing,auto_missing_confirmed,auto_missing_confirmation_reason,masked_missing_candidate,marker_identity_candidate,marker_identity_z,Action |
  Sort-Object {[double]$_.support_score} |
  Format-Table -AutoSize
```

### 2.3 模拟场景完整主流程

如果要继续运行 Stage4 和 Stage5，可以使用主流程脚本：

```powershell
$sample = "<simulation_sample_name>"

python scripts\run_project_mainline.py `
  --sample $sample `
  --project_root . `
  --stage3_sc_expr_source normalized `
  --stage4_sc_expr_source normalized `
  --route2_filter_scope missing_only
```

输出通常位于：

```text
data/processed/simulation_experiments/<group>/<sample>/
result/simulation_experiments/<group>/<sample>/
```

## 3. 真实数据场景运行方式

真实数据场景通常从 `data/raw/low_resolution_experiments/<sample>` 读取。运行时可以使用 R 版 Stage1，或直接使用主流程脚本。

### 3.1 Stage1 + Stage3 最小命令

只需要替换 `$sample`：

```powershell
Set-Location "E:\AAA文件\Experiment\SVTuner\sctuner2.0"
$env:PYTHONNOUSERSITE = "1"
$envPrefix = "E:\ANACONDA\envs\cytospace_v1.1.0_py310"

$sample = "<real_sample_name>"

Remove-Item "data\processed\low_resolution_experiments\$sample\stage3_typematch" -Recurse -Force -ErrorAction SilentlyContinue
Remove-Item "result\low_resolution_experiments\$sample\stage3_typematch" -Recurse -Force -ErrorAction SilentlyContinue

Rscript r_scripts\stage1_preprocess.R `
  --sample $sample `
  --project_root .

if ($LASTEXITCODE -ne 0) { throw "Stage1 failed: $sample" }

& "$envPrefix\python.exe" -m src.stages.stage3_type_plugin `
  --sample $sample `
  --sc_expr_source normalized

if ($LASTEXITCODE -ne 0) { throw "Stage3 failed: $sample" }
```

### 3.2 查看 Stage3 结果

```powershell
$sample = "<real_sample_name>"

Import-Csv "data\processed\low_resolution_experiments\$sample\stage3_typematch\type_support.csv" |
  Select-Object orig_type,n_cells,support_score,support_category,auto_missing,auto_missing_confirmed,auto_missing_confirmation_reason,masked_missing_candidate,marker_identity_candidate,marker_identity_z,Action |
  Sort-Object {[double]$_.support_score} |
  Format-Table -AutoSize
```

### 3.3 真实数据完整主流程

如果要继续运行 Stage4 和 Stage5：

```powershell
$sample = "<real_sample_name>"

python scripts\run_project_mainline.py `
  --sample $sample `
  --project_root . `
  --stage3_sc_expr_source normalized `
  --stage4_sc_expr_source normalized `
  --route2_filter_scope missing_only
```

输出通常位于：

```text
data/processed/low_resolution_experiments/<sample>/
result/low_resolution_experiments/<sample>/
```

## 4. 简单命令样式

最常用的命令样式可以简化成：

```powershell
# 模拟数据
$sample = "<simulation_sample_name>"
python scripts\prepare_stage1_from_sim_source.py --project_root . --sample $sample
python -m src.stages.stage3_type_plugin --sample $sample --sc_expr_source normalized
```

```powershell
# 真实数据
$sample = "<real_sample_name>"
Rscript r_scripts\stage1_preprocess.R --sample $sample --project_root .
python -m src.stages.stage3_type_plugin --sample $sample --sc_expr_source normalized
```

```powershell
# 主流程：Stage1 -> Stage3 -> Stage4 -> Stage5
$sample = "<sample_name>"
python scripts\run_project_mainline.py --sample $sample --project_root .
```

## 5. 注意事项

- 样本名必须和 `configs/datasets/<sample>.yaml` 的文件名一致。
- 模拟数据场景推荐使用 `scripts/prepare_stage1_from_sim_source.py` 准备 Stage1，而不是直接跑 R 版 Stage1。
- 真实数据场景默认使用 `Rscript r_scripts/stage1_preprocess.R` 准备 Stage1。
- 如果只是验证缺失类型检测机制，通常跑到 Stage3 并检查 `type_support.csv` 即可。
- 如果要比较 baseline 与 route2 映射效果，需要继续运行 Stage4 和 Stage5。

## 6. 最近一次 Stage3 验证命令记录

以下命令是最近一次实际用于验证的命令。该验证只运行到 Stage3，用于检查类型缺失检测是否正常；没有继续运行 Stage4、Stage5 或后续映射/评估流程。

### 6.1 模拟数据单缺失场景

验证样本：

```text
human_lung_5loc_fine9_clustered_sim_missing_ciliated
```

预期结果：`Ciliated` 被识别为确认缺失类型，其他类型不应被确认缺失。

```powershell
$ErrorActionPreference = 'Stop'
Set-Location "E:\AAA文件\Experiment\SVTuner\sctuner2.0"
$env:PYTHONNOUSERSITE = "1"
$envPrefix = "E:\ANACONDA\envs\cytospace_v1.1.0_py310"
$py = "$envPrefix\python.exe"

$simGroup = "human_lung_5loc"
$simSample = "human_lung_5loc_fine9_clustered_sim_missing_ciliated"

Remove-Item "data\processed\simulation_experiments\$simGroup\$simSample" -Recurse -Force -ErrorAction SilentlyContinue
Remove-Item "result\simulation_experiments\$simGroup\$simSample" -Recurse -Force -ErrorAction SilentlyContinue

& $py scripts\prepare_stage1_from_sim_source.py `
  --project_root . `
  --sample $simSample

if ($LASTEXITCODE -ne 0) { throw "Stage1 sim export failed: $simSample" }

& $py -m src.stages.stage3_type_plugin `
  --sample $simSample `
  --sc_expr_source normalized

if ($LASTEXITCODE -ne 0) { throw "Stage3 failed: $simSample" }

Import-Csv "data\processed\simulation_experiments\$simGroup\$simSample\stage3_typematch\type_support.csv" |
  Select-Object orig_type,n_cells,support_score,support_category,auto_missing,auto_missing_confirmed,auto_missing_confirmation_reason,masked_missing_candidate,marker_identity_candidate,marker_identity_z,Action |
  Sort-Object {[double]$_.support_score} |
  Format-Table -AutoSize
```

### 6.2 真实数据无类型掩盖场景

验证样本：

```text
adult_mouse_kidney_real
```

预期结果：不应出现确认缺失类型，即 `auto_missing_confirmed=Yes` 的数量应为 0。

```powershell
$ErrorActionPreference = 'Stop'
Set-Location "E:\AAA文件\Experiment\SVTuner\sctuner2.0"
$env:PYTHONNOUSERSITE = "1"
$envPrefix = "E:\ANACONDA\envs\cytospace_v1.1.0_py310"
$py = "$envPrefix\python.exe"

$realSample = "adult_mouse_kidney_real"

Remove-Item "data\processed\low_resolution_experiments\$realSample\stage3_typematch" -Recurse -Force -ErrorAction SilentlyContinue
Remove-Item "result\low_resolution_experiments\$realSample\stage3_typematch" -Recurse -Force -ErrorAction SilentlyContinue

& $py -m src.stages.stage3_type_plugin `
  --sample $realSample `
  --sc_expr_source normalized

if ($LASTEXITCODE -ne 0) { throw "Stage3 failed: $realSample" }

Import-Csv "data\processed\low_resolution_experiments\$realSample\stage3_typematch\type_support.csv" |
  Select-Object orig_type,n_cells,support_score,support_category,auto_missing,auto_missing_confirmed,auto_missing_confirmation_reason,masked_missing_candidate,marker_identity_candidate,marker_identity_z,Action |
  Sort-Object {[double]$_.support_score} |
  Format-Table -AutoSize
```

### 6.3 真实数据单类型掩盖场景

验证样本：

```text
adult_mouse_kidney_real_profile_mask_endo
```

预期结果：`Endo` 被识别为确认缺失类型，其他类型不应被确认缺失。

```powershell
$ErrorActionPreference = 'Stop'
Set-Location "E:\AAA文件\Experiment\SVTuner\sctuner2.0"
$env:PYTHONNOUSERSITE = "1"
$envPrefix = "E:\ANACONDA\envs\cytospace_v1.1.0_py310"
$py = "$envPrefix\python.exe"

$maskSample = "adult_mouse_kidney_real_profile_mask_endo"

Remove-Item "data\processed\low_resolution_experiments\$maskSample\stage3_typematch" -Recurse -Force -ErrorAction SilentlyContinue
Remove-Item "result\low_resolution_experiments\$maskSample\stage3_typematch" -Recurse -Force -ErrorAction SilentlyContinue

& $py -m src.stages.stage3_type_plugin `
  --sample $maskSample `
  --sc_expr_source normalized

if ($LASTEXITCODE -ne 0) { throw "Stage3 failed: $maskSample" }

Import-Csv "data\processed\low_resolution_experiments\$maskSample\stage3_typematch\type_support.csv" |
  Select-Object orig_type,n_cells,support_score,support_category,auto_missing,auto_missing_confirmed,auto_missing_confirmation_reason,masked_missing_candidate,marker_identity_candidate,marker_identity_z,Action |
  Sort-Object {[double]$_.support_score} |
  Format-Table -AutoSize
```
