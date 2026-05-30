# SVTuner 2.0

SVTuner is a pre-mapping type-coordination workflow for spatial transcriptomics cell-to-spot assignment. The project focuses on a specific failure mode: a single-cell reference may contain cell types or states that are not supported by the spatial transcriptomics sample, yet a downstream mapper can still assign those unsupported cells into spatial locations. SVTuner addresses this by detecting unsupported or profile-masked cell types before mapping, then running the same mapping backend with a cleaned cell pool.

In this repository, the maintained primary comparison is:

- `CytoSPACE baseline`: CytoSPACE mapping with the original single-cell labels.
- `SVTuner + CytoSPACE`: Stage3 type diagnostics followed by CytoSPACE route2 mapping using `plugin_type` labels and filtered unsupported cells.

The repository also contains benchmarking utilities and curated visualizations for simulated missing-type settings, real profile-mask settings, CytoSPACE Fig.2-style reproductions, and cell-level high-resolution profile-mask experiments.

## Repository Layout

```text
configs/datasets/      Dataset YAML files for real, simulated, profile-mask, Fig.2-style, and high-resolution scenarios.
r_scripts/             R preprocessing entry point for standard real-data Stage1 export.
scripts/               Data construction, mapping wrappers, metric calculation, and visualization scripts.
src/stages/            Core pipeline stages: environment check, Stage1 IO helpers, Stage3 type diagnostics, Stage4 CytoSPACE mapping.
src/svtuner/           Command-line wrapper and bundle helper.
src/utils/             Shared path and cell-type-name utilities.
visualizations/        Curated figures retained for paper-level inspection.
data/                  Local raw/processed data; ignored by Git.
result/                Local mapping outputs and metric intermediates; ignored by Git.
logs/                  Local run logs; ignored by Git.
external/              External method code or local method resources; large folders are ignored where appropriate.
```

Large data and run outputs are intentionally not versioned. The repository tracks code, configs, documentation, and selected figure outputs.

## Core Pipeline

The current maintained pipeline stops after producing two Stage4 mapping outputs: CytoSPACE baseline and SVTuner + CytoSPACE route2.

### Stage0: Environment Check

```powershell
python -m src.stages.stage0_envcheck
```

### Stage1: Preprocessing

Stage1 converts raw or simulated inputs into a common exported format used by Stage3 and Stage4.

Common exported files:

```text
data/processed/<sample>/stage1_preprocess/exported/sc_expression_normalized.csv
data/processed/<sample>/stage1_preprocess/exported/st_expression_normalized.csv
data/processed/<sample>/stage1_preprocess/exported/sc_metadata.csv
data/processed/<sample>/stage1_preprocess/exported/st_coordinates.csv
```

For real datasets using the R preprocessing path:

```powershell
Rscript r_scripts\stage1_preprocess.R --sample <sample> --project_root . --export_csv
```

For simulation-derived datasets, Stage1 exports can be prepared from an existing simulation source:

```powershell
python scripts\prepare_stage1_from_sim_source.py --project_root . --sample <sample>
```

### Stage3: Unsupported Type Diagnostics

Stage3 is the core SVTuner layer. It does not assign cells to spots. It evaluates whether each single-cell type is supported by the spatial sample and writes adjusted metadata for Stage4.

```powershell
python -m src.stages.stage3_type_plugin --sample <sample> --sc_expr_source normalized
```

Main Stage3 outputs:

```text
data/processed/<sample>/stage3_typematch/type_support.csv
data/processed/<sample>/stage3_typematch/cell_type_relabel.csv
data/processed/<sample>/stage3_typematch/type_prior_matrix.csv
data/processed/<sample>/stage3_typematch/stage3_adjusted_annotations.csv
result/<sample>/stage3_typematch/stage3_summary.json
```

Current Stage3 logic is intentionally simplified to the mechanisms used in the retained experiments:

- Build marker evidence for each reference cell type.
- Score support of each type in the ST expression matrix.
- Detect low-support or profile-masked candidates.
- Use marker-identity diagnostics to avoid assigning unsupported evidence to the wrong type.
- Apply T/NK similarity protection where related immune states are difficult to separate from weak support alone.
- Generate `plugin_type` labels and a type prior matrix for downstream route2 mapping.

Unused legacy rescue/protection branches were removed from the maintained workflow. The retained Stage3 output is therefore easier to audit: unsupported-type decisions are represented directly in `type_support.csv`, `stage3_summary.json`, and `stage3_adjusted_annotations.csv`.

### Stage4: CytoSPACE Baseline and Route2

Baseline mapping uses the original `sc_meta` cell-type labels:

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

Route2 mapping uses Stage3-adjusted `plugin_type` labels and filters Stage3-detected missing/unsupported cells:

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

Main Stage4 outputs:

```text
result/<sample>/stage4_cytospace_baseline/cytospace_output/cell_assignment.csv
result/<sample>/stage4_cytospace_baseline/cytospace_output/fractional_abundances_by_spot.csv
result/<sample>/stage4_cytospace_route2/cytospace_output/cell_assignment.csv
result/<sample>/stage4_cytospace_route2/cytospace_output/fractional_abundances_by_spot.csv
result/<sample>/stage4_cytospace_*/stage4_summary.json
```

To avoid verbose CytoSPACE output and large assigned-expression exports during repeated runs:

```powershell
$env:CYTOSPACE_SKIP_ASSIGNED_EXPRESSION = "1"
```

## Maintained CLI Entry Points

Run the maintained Stage1/Stage3/Stage4 pipeline:

```powershell
python -m svtuner run --sample <sample>
```

Run the project mainline wrapper:

```powershell
python scripts\run_project_mainline.py --sample <sample> --project_root .
```

Build a distributable code/config/docs bundle:

```powershell
python -m svtuner bundle
python -m svtuner bundle --include-raw-data --include-results --name svtuner_full_delivery
```

## Data Scenario Groups

### 1. Low-resolution real profile-mask scenarios

These are real ST datasets with one target type masked or suppressed to test whether Stage3 and route2 can avoid reconstructing unsupported target-like signal.

Current representative config names include:

```text
adult_mouse_kidney_real_profile_mask_endo
ffpe_mouse_brain_sagittal_real_profile_mask_microglia
human_breast_cancer_real_profile_mask_basal_cell
human_breast_cancer_visium_ff_wta_real_profile_mask_macrophage
human_breast_cancer_wta_120_real_profile_mask_endothelial_cell
human_cervical_cancer_real_profile_mask_epithelial_cell
human_heart_ff_real_profile_mask_endothelial_cell
human_intestine_cancer_real_profile_mask_endothelial_cell
human_lymph_node_real_profile_mask_b_cell
mouse_embryo_real_profile_mask_erythroid
```

The key visualization family is:

```text
visualizations/masked_scenarios/
visualizations/simulations/real_profile_mask_fig2c_only/
visualizations/simulations/real_profile_mask_fig2d_only/
visualizations/simulations/real_profile_mask_expression_recovery/
```

### 2. Simulated missing-type scenarios

Simulations are built on real or real-like spatial scaffolds and contain explicit truth files. They are used for direct mapping-quality comparisons.

Maintained simulation groups:

```text
real_brca
human_lung_5loc
mouse_brain_refined
```

Representative samples:

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

Simulation truth files are expected under `data/sim/<group>/<sample>/` and copied into Stage1 exports when needed:

```text
sim_info.json
sim_truth_query_cell_spot.csv
sim_truth_spot_type_fraction.csv
```

### 3. 10% scRNA reference-noise scenarios

The `_scnoise10` scenarios perturb the single-cell reference expression while keeping ST expression, coordinates, and truth unchanged. They are used for robustness checks across the same method set.

Generator:

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

### 4. CytoSPACE Fig.2-style profile-mask experiments

The repository keeps several CytoSPACE Fig.2-inspired experiments. These are not direct claims that the datasets are identical to the paper in every detail; they are controlled experiments using the same style of biological readout and visualization logic.

Retained figure families include:

```text
visualizations/cytospace_fig2c_melanoma_stage3_profile_mask/
visualizations/cytospace_fig2d_profile_mask_mapping/
visualizations/cytospace_fig2d_profile_mask_benchmark/
visualizations/cytospace_fig2e_stage3_profile_mask/
visualizations/cytospace_fig2i_mouse_kidney_stage3_unsupported_decoy_sensitivity/
visualizations/cytospace_fig2k_tcell_states_stage3_decoy/
```

These experiments are designed to distinguish two cases:

- Forced/white-list filtering, which is not sufficient as evidence for Stage3.
- Stage3-detected profile-mask or decoy scenarios, where unsupported targets are identified from expression evidence before route2 mapping.

Only the latter are retained as the main SVTuner-supporting experiments.

### 5. Cell-level high-resolution profile-mask scenarios

The high-resolution cell-level scenarios use cell-level spatial data from selected Vizgen-style datasets. These are treated as high-resolution/cell-level experiments, distinct from spot-level low-resolution ST scenarios.

Current retained high-resolution profile-mask configs:

```text
highres_humanbreastcancerpatient1_profile_mask_monocytes_and_macrophages
highres_humancoloncancerpatient1_profile_mask_fibroblasts
highres_humanlungcancerpatient1_profile_mask_plasma_cells
highres_humanmelanomapatient1_profile_mask_fibroblasts
highres_humanmelanomapatient2_profile_mask_b_cells
```

Main high-resolution visualization outputs:

```text
visualizations/highres_profile_mask_mapping/highres_profile_mask_mapping_stack_5x4.png
visualizations/highres_profile_mask_fig2c_only/fig2_panel_c_highres_cell_profile_mask.png
visualizations/highres_profile_mask_fig2d_benchmark/fig2d_highres_profile_mask_benchmark.png
visualizations/highres_profile_mask_fig2c_expression_enrichment/
```

## External Mapping Methods

The retained method-comparison runners are:

### Tangram

Marker-gene mode:

```powershell
python scripts\run_tangram_marker_mapping.py `
  --project_root . `
  --group <group> `
  --sample <sample> `
  --top_n_marker 50 `
  --num_epochs 200 `
  --device cpu
```

All-gene mode:

```powershell
python scripts\run_tangram_marker_mapping.py `
  --project_root . `
  --group <group> `
  --sample <sample> `
  --gene_mode all `
  --num_epochs 200 `
  --device cpu
```

### CellTrek-style runner

```powershell
python scripts\run_celltrek_mapping.py `
  --project_root . `
  --group <group> `
  --sample <sample> `
  --max_genes 2000 `
  --n_pcs 30 `
  --ntree 500
```

Standardized method outputs are written under:

```text
result/<sample>/stage4_mapping/<method>/
```

Expected files include:

```text
cell_assignment.csv
spot_type_fraction.csv
metrics_simulation.json
run.log
```

Pearson correlation, Euclidean distance, novoSpaRc, and SpaOTsc were used during exploration but are no longer part of the maintained runnable script set.

## Key Visualization Outputs

The following retained outputs are currently the most relevant for paper-level review:

```text
visualizations/simulations/simulation_triptych_overview_stack_3datasets.png
visualizations/masked_scenarios/masked_scenarios_real_stack_10x4.png

visualizations/simulations/real_profile_mask_fig2c_only/fig2_panel_c_real_profile_mask.png
visualizations/simulations/real_profile_mask_fig2d_only/fig2_panel_d_real_profile_mask.png
visualizations/simulations/real_profile_mask_expression_recovery/expression_recovery_cosine_summary_bar.png

visualizations/cytospace_fig2c_melanoma_stage3_profile_mask/fig2c_cd4_ce9_stage3_detected_macrophage_mask_baseline_vs_route2.png
visualizations/cytospace_fig2d_profile_mask_mapping/cytospace_fig2d_profile_mask_mapping_stack_6x4.png
visualizations/cytospace_fig2d_profile_mask_benchmark/fig2d_profile_mask_benchmark.png
visualizations/cytospace_fig2e_stage3_profile_mask/fig2e_stage3_profile_mask_route2_ce9_ce10.png
visualizations/cytospace_fig2i_mouse_kidney_stage3_unsupported_decoy_sensitivity/fig2i_stage3_detected_state32like_n1000_baseline_vs_route2_four_panel.png
visualizations/cytospace_fig2k_tcell_states_stage3_decoy/fig2k_stage3_detected_cd4_state_decoy_baseline_vs_route2.png

visualizations/highres_profile_mask_mapping/highres_profile_mask_mapping_stack_5x4.png
visualizations/highres_profile_mask_fig2c_only/fig2_panel_c_highres_cell_profile_mask.png
visualizations/highres_profile_mask_fig2d_benchmark/fig2d_highres_profile_mask_benchmark.png
```

Simulation overview figures compare truth, CytoSPACE baseline mapping, and SVTuner + CytoSPACE route2 mapping. Real profile-mask and Fig.2-style figures focus on whether route2 suppresses unsupported target-like signal or improves downstream biological readouts after Stage3-detected profile masking.

## Environment Notes

Primary environment:

```powershell
conda activate cytospace_v1.1.0_py310
Set-Location "E:\AAA文件\Experiment\SVTuner\sctuner2.0"
$py = "E:\ANACONDA\envs\cytospace_v1.1.0_py310\python.exe"
```

For CytoSPACE-heavy runs:

```powershell
$env:PYTHONNOUSERSITE = "1"
$env:CYTOSPACE_SKIP_ASSIGNED_EXPRESSION = "1"
```

For Tangram, CellTrek, and other user-site installed tools:

```powershell
Remove-Item Env:PYTHONNOUSERSITE -ErrorAction SilentlyContinue
$env:OMP_NUM_THREADS = "1"
$env:MKL_NUM_THREADS = "1"
$env:NUMBA_CACHE_DIR = Join-Path $PWD ".numba_cache"
```

Some CytoSPACE Fig.2 reproduction work used a separate R/Seurat environment during exploration. The maintained primary workflow does not require merging that R environment into the Python/CytoSPACE environment.

## Storage and Git Policy

The project intentionally ignores large local data, run outputs, logs, and cache directories:

```text
data/
result/
logs/
.numba_cache/
external/CellTrek/
external/cytospace/data/
external/cytospace/images/
```

This keeps the Git repository focused on reproducible code, configs, documentation, and curated visualization outputs.

## Current Project Positioning

SVTuner should be described as a pre-mapping unsupported-type diagnostic and coordination layer rather than a replacement for CytoSPACE. The cleanest comparison is to run the same CytoSPACE backend twice:

1. Baseline with original reference labels.
2. Route2 after Stage3 detects and filters unsupported or profile-masked cell types.

The retained experiments are organized to support this claim across:

- Simulated missing-type datasets with explicit truth.
- Low-resolution real profile-mask scenarios.
- CytoSPACE Fig.2-style biological readouts.
- Cell-level high-resolution profile-mask scenarios.

When interpreting results, avoid framing route2 as a generic improvement for all mapping tasks. Its intended advantage is strongest when the single-cell reference contains unsupported, missing, or profile-masked cell types/states that would otherwise be assigned into spatial locations by a downstream mapper.

## Known Boundary Cases

Stage3 is intentionally conservative. This is important for avoiding false removal of valid cell types, but it also means that some highly similar immune subtypes can be difficult to call as missing without stronger lineage-specific evidence.

One historical example is a Human Lymph Node simulation where `CD4 Treg` was replaced by `CD4 T cell`. Because these two labels share broad pan-T-cell expression, the detector can retain enough support signal for `CD4 Treg` and fail to mark it as missing. This is a false-negative missing-type diagnosis rather than a false-positive filtering problem. The current interpretation is:

- SVTuner is suitable for unsupported types or states with detectable expression-profile loss.
- It should not be claimed to perfectly resolve every closely related subtype replacement.
- For very similar immune subtypes, additional lineage-specific markers or subtype-level priors may be needed.

Obsolete design ideas that are no longer part of the maintained workflow include Stage2 SVG-aware weighting, Stage5/6/7 reporting layers, V5 rescue branches, and the exploratory generic backend API. The maintained project is the simpler Stage1 -> Stage3 -> Stage4 workflow documented above.
