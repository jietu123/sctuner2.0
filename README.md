# SVTuner 项目说明（2025）

SVTuner：面向空间变异性（SVG）与类型感知的细胞–空间位点映射调谐框架。当前仓库已接入 GitHub，Stage0 占位，Stage1 预处理已实现并可运行，后续阶段按 SVTuner.md 规划逐步落地。

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
  - `main.py`：管线中控（当前实现 Stage0 占位、Stage1）
  - `config.py`：默认参数与路径
  - `stages/`：各阶段入口（目前 Stage0/Stage1 占位/空）
  - `stages/backends/`：后端占位（cytospace 等）
- `logs/`：运行日志（Stage0 envcheck 占位）
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

## 运行管线（当前实现 Stage0 占位 + Stage1）
默认仅跑 Stage1：
```pwsh
cd D:\Experiment\SVTuner_Project
python src/main.py --sample real_brca
```
如需包含占位 Stage0（仅打印“环境配置成功”）：
```pwsh
python src/main.py --sample real_brca --stages 0,1
```
如需自定义 R 命令：
```pwsh
python src/main.py --sample real_brca --r-cmd "Rscript" --stages 1
```
当前主程序仅调度 Stage0（占位）与 Stage1（R 预处理），其余阶段会提示“未实现”。

## Stage1 预处理要点（r_scripts/stage1_preprocess.R）
- 读取配置：`configs/project_config.yaml`、`configs/datasets/<sample>.yaml`，支持覆盖 QC/路径。
- 温和过滤（可调）：`qc` 与 `gene_filter` 参数，支持 mt_pattern；支持 SVG/marker 白名单强制保留。
- 输出：
  - `data/processed/<sample>/stage1_preprocess/`：
    - `sc_processed.rds`, `st_processed.rds`
    - `common_genes.txt`, `hvg_genes.txt`
  - `result/<sample>/stage1_preprocess/`：
    - `stage1_summary.json`
    - `qc_plots/`（基础直方图：counts/genes/percent.mt）
- Summary 包含：过滤前后细胞/基因数、whitelist 保留数、公共基因数与路径、HVG 请求/实际数、mt% 统计等。
- 运行（独立调用）：  
```pwsh
conda run -n cytospace_v1.1.0_py310 Rscript r_scripts/stage1_preprocess.R --sample real_brca
Rscript r_scripts/stage1_preprocess.R --sample <sample>  # 如已将 Rscript 加入 PATH
```
（Windows 建议 `CONDA_DLL_SEARCH_MODIFICATION_ENABLE=1`，main.py 已默认设置）

## Git/LFS 注意
- 大文件（raw 数据）已由 Git LFS 跟踪；请确保本地安装并启用 Git LFS。
- 默认 `.gitignore` 忽略 `result/`，`data/processed/` 目前未忽略，如需忽略可自行调整。

## 后续阶段（占位说明）
- Stage0：当前为占位提示，未来可接入真实环境检查。
- Stage2/3/4/5/6/7：按 `dosc/SVTuner.md` 规划依次接入（SVG/HVG 加权、类型对齐、映射 orchestrator、模拟/真实评估、报告）。

## 常见问题
- Rscript 找不到：确认 `--r-cmd` 指向有效 Rscript（建议 conda run）。
- DLL 相关报错（Windows）：使用默认 `--r-cmd "conda run -n cytospace_v1.1.0_py310 Rscript"`，或手动设置 `CONDA_DLL_SEARCH_MODIFICATION_ENABLE=1`。
- 输入文件名不匹配：在 `configs/datasets/<sample>.yaml` 显式配置 `paths`。
