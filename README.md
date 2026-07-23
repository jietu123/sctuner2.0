# SVTuner 2.0

SVTuner is a reference-adequacy and abstention layer for spatial transcriptomics mapping. It detects two reciprocal forms of single-cell/spatial-reference mismatch before or alongside a downstream mapper:

- **Stage3A, SC-only mismatch:** a cell type or state exists in the single-cell reference but is unsupported by the spatial sample. SVTuner filters or relabels those reference cells before mapping.
- **Stage3B, ST-only mismatch:** a spatial expression region exists in ST but is unsupported by the available single-cell reference. SVTuner identifies those spots and allows the mapping workflow to abstain instead of forcing a surrogate assignment.

The maintained backend is [CytoSPACE](external/cytospace). SVTuner is not a replacement for CytoSPACE and should not be described as a generally superior mapping algorithm. Its intended use is to make a mapper mismatch-aware when the reference and spatial sample do not contain the same biological types or states.

This is a research repository. Code, dataset configurations, audit tables, reproducibility records, and selected paper figures are versioned. Most raw data and large run outputs are local and intentionally excluded from Git.

## Release Candidate

This repository is associated with the manuscript reproducibility release candidate `v1.0.0-rc1`. The frozen reproducibility package, manifests, source-value records, environment specifications, configuration records, and audit summaries are indexed in [`reproducibility_release/`](reproducibility_release/README.md).

GitHub releases are available from the [SVTuner Releases page](https://github.com/jietu123/sctuner2.0/releases). Large raw datasets and complete mapping intermediates are not stored in Git. Third-party inputs must be reacquired from the accessions and official provider routes recorded in the frozen dataset manifest.

## Current Evidence Snapshot

The current simulation endpoint is an **abstention-aware whole-space recovery score**. Every truth spot remains in the denominator. A non-withheld spot receives its normalized composition-overlap score; an SVTuner-withheld spatial unit receives `1` only when independent simulation truth confirms that the SC-reference-dropped type is dominant at that spot, and an incorrect withholding decision receives `0`.

| Method | No noise, mean (n=9) | 10% SC noise, mean (n=9) |
|---|---:|---:|
| CytoSPACE | 0.5154 | 0.5046 |
| **SVTuner** | **0.7104** | **0.6553** |
| Tangram, all genes | 0.4377 | 0.4513 |
| Tangram, marker genes | 0.4628 | 0.4687 |
| novoSpaRc | 0.4914 | 0.4889 |
| SpaOTsc | 0.3950 | 0.3836 |
| CellTrek | 0.5309 | 0.4352 |

Under no added noise, SVTuner exceeded CytoSPACE in `9/9` scenarios. Of `10,445` predicted withheld spatial units, `10,389` were correct, giving abstention precision `0.9946` and recall `0.9540` against `10,890` truth-unsupported spots.

Under 10% scRNA-seq expression perturbation, the complete Stage3A plus Stage3B route achieved a mean score of `0.6553`, compared with `0.5046` for CytoSPACE, and exceeded CytoSPACE in `8/9` scenarios. Abstention precision and recall were `0.9614` and `0.9039`. The result supports an aggregate benefit, not uniform robustness across every dataset background.

The final biological application uses an independent computational-pathology annotation (CTA) Immune endpoint in breast-cancer Visium, not the retired Xenium/B-lineage candidate. Under an immune-all-dropout reference perturbation, the frozen analysis reported:

| Biological-application metric | Value |
|---|---:|
| Analysis spots, positive / negative | 134 / 1,754 |
| Mean withheld score, positive / negative | 0.5206 / 0.3752 |
| AUROC / average precision (AP) | 0.8305 / 0.2837 |
| Binary withheld rate, positive / negative | 44.78% / 2.17% |
| Risk ratio / odds ratio | 20.67 / 36.61 |
| Fisher exact P | 1.95e-49 |
| CTA-positive vs withheld Jaccard | 0.2691 |
| Binary / continuous forced-burden prevention | 0.4478 / 0.5206 |

These are withheld-aware burden estimates under a frozen reference-dropout experiment, not post-remapping cell-composition improvements or new biological discoveries.

## Repository Layout

```text
configs/                 Environment, pipeline, and dataset YAML files
data/                    Local raw, simulated, and processed data (Git-ignored)
docs/                    Development audits (currently Git-ignored)
external/cytospace/      Vendored CytoSPACE backend
r_scripts/               R Stage1 preprocessing
result/                  Local mappings and intermediates (Git-ignored)
scripts/                 Experiment, benchmark, audit, and figure entry points
src/stages/              Stage0, Stage1 IO, Stage3A, Stage3B, and Stage4
src/svtuner/             CLI and pipeline orchestration
src/utils/               Shared path and cell-type utilities
tests/                   Stage3B and Stage4 withholding-restoration tests
visualizations/          Curated figures and source-value tables
reproducibility_release/ Frozen release candidate, manifests, and public audit indexes
```

## Installation

The validated environment is Windows with Python 3.10.19. The complete pinned environment, including R/Seurat and the local CytoSPACE install, is in [`configs/environment.yml`](configs/environment.yml).

```powershell
conda env create -f configs/environment.yml
conda activate cytospace_v1.1.0_py310
python -m pip install -e .
```

If the environment already exists:

```powershell
conda activate cytospace_v1.1.0_py310
python -m pip install -e external/cytospace
python -m pip install -e .
```

The package metadata intentionally does not install the scientific stack by itself; use the pinned environment rather than relying on `pip install -e .` alone.

When `Rscript` is not on `PATH`, copy the local configuration template and set the machine-specific path:

```powershell
Copy-Item configs/project_config.local.yaml.example configs/project_config.local.yaml
```

`configs/project_config.local.yaml` is Git-ignored. Do not commit absolute machine paths.

Verify the environment and CLI:

```powershell
svtuner envcheck
svtuner version
```

For repeated CytoSPACE runs, assigned-expression export can be disabled to reduce storage:

```powershell
$env:PYTHONNOUSERSITE = "1"
$env:CYTOSPACE_SKIP_ASSIGNED_EXPRESSION = "1"
```

The external Tangram, novoSpaRc, SpaOTsc, and CellTrek benchmark methods may use a separate Python environment. Their runners write a common output contract and are not required for the core SVTuner/CytoSPACE workflow.

## Input Contract

Each sample requires `configs/datasets/<sample>.yaml` and Stage1 exports under:

```text
data/processed/<sample>/stage1_preprocess/exported/
```

Core files are:

```text
sc_expression_normalized.csv   cells x genes
st_expression_normalized.csv   spots x genes
sc_metadata.csv                cell annotations
st_coordinates.csv             spot coordinates
```

Simulation evaluation additionally uses:

```text
sim_info.json
sim_truth_query_cell_spot.csv
sim_truth_spot_type_fraction.csv
```

Real-data preprocessing:

```powershell
Rscript r_scripts/stage1_preprocess.R `
  --sample <sample> `
  --project_root . `
  --export_csv
```

Simulation-derived Stage1 preparation:

```powershell
python scripts/prepare_stage1_from_sim_source.py `
  --project_root . `
  --sample <sample>
```

## Core Workflows

### Stage3A: unsupported SC reference types

Stage3A builds type-level marker evidence, evaluates support in ST, applies the maintained identity/similarity protections, and produces `plugin_type` annotations for mapping.

```powershell
python -m src.stages.stage3_type_plugin `
  --sample <sample> `
  --sc_expr_source normalized
```

Main outputs:

```text
data/processed/<sample>/stage3_typematch/type_support.csv
data/processed/<sample>/stage3_typematch/cell_type_relabel.csv
data/processed/<sample>/stage3_typematch/type_prior_matrix.csv
data/processed/<sample>/stage3_typematch/stage3_adjusted_annotations.csv
result/<sample>/stage3_typematch/stage3_summary.json
```

Run the same CytoSPACE backend as an unfiltered baseline and as the Stage3A-adjusted route:

```powershell
# Baseline
python -m src.stages.stage4_cytospace `
  --sample <sample> `
  --filter_mode none `
  --cell_type_column sc_meta `
  --stage4_suffix _baseline

# SVTuner Stage3A route
python -m src.stages.stage4_cytospace `
  --sample <sample> `
  --filter_mode plugin_unknown `
  --filter_scope unsupported_all `
  --cell_type_column plugin_type `
  --stage4_suffix _route2
```

The convenience wrapper executes Stage1, Stage3A, baseline Stage4, and route2 Stage4:

```powershell
svtuner run --sample <sample>
```

For simulation inputs that do not need the R preprocessing path:

```powershell
svtuner run --sample <sample> --use-python-stage1
```

### Stage3B: unsupported ST regions

Stage3B fits supported SC-reference mixtures to ST and combines reconstruction-error and reference-orthogonal residual evidence. It uses self-calibration, spot-level FDR, spatial coherence, and permutation-based region validation. The algorithm does not read simulation truth, missing-type labels, CTA endpoint labels, or a target whitelist; truth and external endpoints are used only after execution for evaluation.

```powershell
svtuner stage3b `
  --sample <sample> `
  --fdr 0.05 `
  --n-spatial-permutations 200
```

Equivalent module entry point:

```powershell
python -m src.stages.stage3b_st_unsupported `
  --sample <sample> `
  --n_spatial_permutations 200
```

Main outputs:

```text
data/processed/<sample>/stage3b_st_unsupported/spot_unsupported_scores.csv
data/processed/<sample>/stage3b_st_unsupported/unsupported_regions.csv
data/processed/<sample>/stage3b_st_unsupported/unsupported_region_residual_genes.csv
data/processed/<sample>/stage3b_st_unsupported/supported_mixture_weights.csv
result/<sample>/stage3b_st_unsupported/stage3b_summary.json
```

To preserve Stage3B regions as withheld spatial units, pass the mask to Stage4. The `blank` token in the following option and output suffix is retained as a legacy compatibility identifier:

```powershell
python -m src.stages.stage4_cytospace `
  --sample <sample> `
  --filter_mode none `
  --cell_type_column sc_meta `
  --stage3b_blank_regions `
  --stage4_suffix _stage3b_blank
```

Masked spots are excluded from mapping capacity. Final spot-level outputs restore the full coordinate universe and represent masked spots as all-zero rows. `stage4_summary.json` records mask and zero-row audits.

Main Stage4 outputs are:

```text
result/<sample>/stage4_cytospace_<suffix>/cytospace_output/cell_assignment.csv
result/<sample>/stage4_cytospace_<suffix>/cytospace_output/fractional_abundances_by_spot.csv
result/<sample>/stage4_cytospace_<suffix>/stage4_summary.json
```

## Joint Simulation Benchmark

The retained composite benchmark crosses Stage3A missing-type conditions with one Stage3B reference-dropout target:

| Group | Stage3B target | Stage3A rows |
|---|---|---|
| Real BRCA | Endothelial cells | control; Epithelial cells; Epithelial cells + PCs |
| Human lung 5-location | B cell | control; AT2; AT2 + Fibroblast |
| Mouse brain refined | Ext_L56 | control; Micro; Micro + Oligo_2 |

The same nine scenarios are evaluated with no added noise and with 10% SC-expression perturbation. In the noise experiment, the measured perturbation fraction is `0.099912-0.099969`; ST expression, coordinates, and simulation truth remain unchanged. Seven methods produce `63` method-scenario rows per noise level.

External-method runs can be resumed with:

```powershell
# No noise
python scripts/run_nine_scenario_method_benchmark.py `
  --project_root . `
  --scenario_preset composite

# 5% SC noise
python scripts/run_nine_scenario_method_benchmark.py `
  --project_root . `
  --scenario_preset composite `
  --sample_suffix _scnoise05

# 10% SC noise
python scripts/run_nine_scenario_method_benchmark.py `
  --project_root . `
  --scenario_preset composite `
  --sample_suffix _scnoise10
```

These commands run Tangram, novoSpaRc, SpaOTsc, and CellTrek. The no-noise summary additionally requires the CytoSPACE baseline and SVTuner Stage3B-withholding outputs. The historical output directory still contains `blank` in its compatibility name. The full 5% and 10% noise SVTuner routes are run separately below.

Rebuild the abstention-aware summaries and boxplots:

```powershell
# No noise
python scripts/plot_method_comparison_composition_recovery.py `
  --project_root . `
  --scenario_preset composite `
  --route2_stage4_dir stage4_cytospace_stage3b_blank `
  --reward_correct_stage3b_abstention `
  --abstention_truth_rule target_dominant `
  --out_dir visualizations/method_comparison/composite_no_noise `
  --output_prefix composition_recovery_7mapping_methods_composite_no_noise_abstention_aware `
  --title "Spatial cell-type composition recovery, abstention-aware whole-space evaluation" `
  --ylabel "Abstention-aware recovery score"

# 5% SC noise, full Stage3A + Stage3B SVTuner route and all baselines.
# Set these to the local CytoSPACE and multi-method Python executables.
$cytoPython = "path\to\cytospace\python.exe"
$multiMethodPython = "path\to\multimethod\python.exe"
& $cytoPython scripts/run_composite_scnoise05_stage3ab_full.py `
  --project_root . `
  --python $cytoPython `
  --external_python $multiMethodPython `
  --prepare_inputs `
  --run_baselines `
  --run_external `
  --execute `
  --resume

# 10% SC noise, full Stage3A + Stage3B SVTuner route
python scripts/run_composite_scnoise10_stage3ab_full.py `
  --project_root . `
  --execute `
  --resume
```

Authoritative result tables:

```text
visualizations/method_comparison/composite_no_noise/
  composition_recovery_7mapping_methods_composite_no_noise_abstention_aware_summary.csv
  composition_recovery_7mapping_methods_composite_no_noise_abstention_aware_scenario.csv

visualizations/method_comparison/composite_scnoise05_stage3ab_full/
  composite_scnoise05_stage3ab_full_method_summary.csv
  composite_scnoise05_stage3ab_full_source_values.csv
  composite_scnoise05_stage3ab_full_stage3a_audit.csv
  composite_scnoise05_stage3ab_full_stage3b_audit.csv
  composite_scnoise05_stage3ab_full_summary.json

visualizations/method_comparison/composite_scnoise10_stage3ab_full/
  composite_scnoise10_stage3ab_full_method_summary.csv
  composite_scnoise10_stage3ab_full_source_values.csv
  composite_scnoise10_stage3ab_full_stage3a_audit.csv
  composite_scnoise10_stage3ab_full_stage3b_audit.csv
  composite_scnoise10_stage3ab_full_summary.json
```

The previous naive whole-space score incorrectly treated correct SVTuner withholding decisions as zero. The previous predicted-supported-region score used a method-dependent evaluation subset. Both are retired and must not be mixed with the current benchmark.

## Biological Application

The completed application is a breast-cancer Visium experiment with a frozen, independent CTA-defined `Immune cells` endpoint:

1. CTA objects were registered to `2,248` Visium spots before mapping evaluation.
2. The frozen main set contains `134` endpoint-positive and `1,754` endpoint-negative spots; `290` ambiguous and `70` unmatched spots are excluded from the primary comparison.
3. All `2,414` immune cells were removed from a `4,014`-cell SC reference, leaving `1,600` non-immune cells.
4. CytoSPACE, which has no abstention mechanism, necessarily assigned the immune-positive spots to the remaining non-immune reference.
5. Stage3B ran at its predefined FDR `0.05` without CTA labels or endpoint-tuned thresholds.
6. Phase 8 evaluated the frozen Stage3B output; Phase 9 locked the allowed interpretation.

The Xenium-defined B-lineage/humoral candidate was screened in Phase 0 but retired after the formal CytoSPACE feasibility mainline stopped. It is not the completed biological application and should not be presented as such.

Phase scripts and frozen outputs are under:

```text
scripts/run_bioapp_phase*.py
visualizations/bioapp_experiment/
```

The final audit can be rerun when its local prerequisites are present:

```powershell
python scripts/run_bioapp_phase9_final_biological_application_audit.py
```

The recommended paper figure is:

```text
visualizations/bioapp_experiment/
  bioapp_main_figure_v3_12_panel_A_evidence_chain_redesign/
  fig_bioapp_main_composite_v3_12.svg
```

The exploratory downstream D0-D3 analyses test morphology, interface, and local microenvironment context. They support nonrandom spatial context but do not establish a new pathological niche or causal biology.

## Other Retained Evidence

Beyond the joint simulations and CTA application, the repository retains:

- ten low-resolution real profile-mask scenarios;
- six real Stage3B reference-dropout settings;
- a detailed mouse-brain `ST8059051 / Thalamic excitatory` Stage3B case;
- CytoSPACE Fig.2-style state-enrichment and ordering experiments;
- five cell-level high-resolution profile-mask datasets;
- forced-assignment, false-spatial-niche, and communication stress tests.

Representative figures:

```text
visualizations/simulations/simulation_stage3ab_joint_triptych_overview_stack_3datasets.svg
visualizations/stage3b_realdata_candidate_scan/spatial_9x2/
  stage3b_reference_dropout_spatial_stack_recommended_6x2.svg
visualizations/masked_scenarios/masked_scenarios_real_stack_10x4.png
visualizations/highres_profile_mask_fig2d/targeted_validation/
  targeted_fig2d_highres_fixed_panel.png
```

## Testing

The maintained unit tests cover Stage3B calibration/spatial logic and Stage4 withholding-mask restoration:

```powershell
python -m pytest -q tests
```

The explicit `tests` path is required because local ignored raw-data folders may contain vendored upstream repositories with their own optional-dependency test suites.

Useful lightweight checks before committing documentation or analysis changes:

```powershell
git diff --check
svtuner envcheck
```

## Data and Git Policy

The repository ignores large or machine-specific content, including:

```text
data/
result/
logs/
.numba_cache/
configs/project_config.local.yaml
external/CellTrek/
external/cytospace/data/
external/cytospace/images/
```

Selected visualizations and source-value CSVs are committed as the auditable paper evidence layer. Reproducing mapping from scratch requires obtaining the corresponding local raw/reference data and satisfying each dataset's license and access conditions.

The `v1.0.0-rc1` public reproducibility layer is indexed in [`reproducibility_release/README.md`](reproducibility_release/README.md). It excludes restricted provider data, complete GEO/10x copies, complete Space Ranger outputs, raw MERSCOPE provider files, large expression matrices, and complete CytoSPACE intermediate outputs.

Create a distributable code/configuration bundle with:

```powershell
svtuner bundle
svtuner bundle --include-raw-data --include-results --name svtuner_full_delivery
```

## Interpretation Boundaries

Supported statements:

- Stage3A can reduce unsupported reference-type reconstruction under tested mismatch conditions.
- Stage3B can identify spatial regions insufficiently explained by the available SC reference and withhold them from forced mapping.
- SVTuner improves the mean abstention-aware whole-space score in the tested no-noise and 10% SC-noise composite benchmarks.
- The frozen CTA application shows concordance between SVTuner withholding and an independent immune-associated pathology endpoint under immune-all-dropout.

Unsupported statements:

- SVTuner perfectly detects every missing type or every target-positive spot.
- SVTuner is uniformly superior to all mapping methods or under all noise conditions.
- A correct withholding decision is equivalent to a conventional composition prediction.
- The CTA experiment proves a new breast-cancer immune mechanism or reports post-remapping composition improvement.
- Simulation truth, target labels, or CTA endpoint labels are used by the Stage3B algorithm. They are used only for post hoc evaluation.

Closely related cell states remain a difficult boundary case because broad lineage markers may preserve apparent support after subtype replacement. Real-data marker cores are evaluation proxies, not definitive single-cell ground truth.

## License and Status

`pyproject.toml` currently declares a proprietary license, but this repository does not yet contain a project-level `LICENSE` file. No release license is inferred here. External components retain their own licenses; see [`external/cytospace/LICENSE`](external/cytospace/LICENSE).

The manuscript-associated reproducibility version is `v1.0.0-rc1`; the Python package metadata remains `0.1.0`. The maintained command surface is `svtuner run`, `svtuner envcheck`, `svtuner stage3b`, `svtuner bundle`, and `svtuner version`. The many scripts under `scripts/` are experiment-specific research entry points and are not a stable public API.
