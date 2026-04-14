# SVTuner Project

Current scope has been simplified to the core real-data pipeline (`real_brca`) and the Stage3 plugin workflow.
Historical simulation batches and large visualization bundles were removed.

## Current repository scope

- One maintained dataset config: `configs/datasets/real_brca.yaml`
- Core stages:
  - Stage0: environment check
  - Stage1: preprocess (R)
  - Stage3: plugin type matching / filtering
  - Stage4: CytoSPACE mapping (baseline + route2)
  - Stage5: evaluation
  - Stage6: real-data audit
  - Stage7: report assembly (outputs under `result/<sample>/stage7_report/`)
- Utility scripts still in use:
  - `scripts/extract_xenium_panel_genes.py`
  - `scripts/build_xenium_pseudospots.py`
  - `scripts/cleanup_redundant_storage.py`

## Quick start

### 1) Environment check

```bash
python -m src.stages.stage0_envcheck
```

### 2) Run full pipeline via CLI

```bash
python -m svtuner run --sample real_brca
```

### 2.2) Run the packaged mainline script (Stage1/3/4/5)

```bash
python scripts/run_project_mainline.py --sample real_brca --project_root .
```

This script executes:
- Stage1 preprocessing
- Stage3 plugin filtering
- Stage4 baseline + route2 mapping
- Stage5 metrics summary

It intentionally does **not** run Stage7 visualization/reporting.
Stage5 behavior is auto:
- if sim truth files exist (`sim_info.json`, `sim_truth_spot_type_fraction.csv`), it runs simulation Stage5 evaluator;
- otherwise it runs real-data Stage5 summary (`src.stages.stage5_real_summary`).

### 2.1) Build a distributable bundle

```bash
# default: code/config/docs only (no raw data, no results)
python -m svtuner bundle

# include raw data and results if needed
python -m svtuner bundle --include-raw-data --include-results --name svtuner_full_delivery
```

Bundle output is written under `dist/` as `<bundle_name>.zip`.

### 3) Run key stages manually

```bash
# Stage1 (R preprocess)
Rscript r_scripts/stage1_preprocess.R --sample real_brca --project_root .

# Stage3
python -m src.stages.stage3_type_plugin --sample real_brca

# Stage4 baseline
python -m src.stages.stage4_cytospace --sample real_brca --filter_mode none --cell_type_column sc_meta --stage4_suffix _baseline

# Stage4 route2
python -m src.stages.stage4_cytospace --sample real_brca --filter_mode plugin_unknown --cell_type_column plugin_type --filter_scope unsupported_all --stage4_suffix _route2

# Stage5 baseline / route2 (example)
python -m src.stages.stage5_route2_s0 --sample real_brca --run_tag baseline --stage4_dir result/real_brca/stage4_cytospace_baseline/cytospace_output --out_dir result/real_brca/stage5_eval_baseline
python -m src.stages.stage5_route2_s0 --sample real_brca --run_tag route2 --stage4_dir result/real_brca/stage4_cytospace_route2/cytospace_output --out_dir result/real_brca/stage5_eval_route2

# Stage6
python -m src.stages.stage6_real_audit --sample real_brca --stage4_id baseline
python -m src.stages.stage6_real_audit --sample real_brca --stage4_id route2

# Stage7
python -m src.stages.stage7_report --sample real_brca --baseline_id baseline --route2_id route2
```

## Notes

- `visualizations/`, simulation outputs, and deprecated helper scripts were intentionally removed.
- Keep new additions aligned with this reduced scope; avoid reintroducing obsolete simulation-only flows.
