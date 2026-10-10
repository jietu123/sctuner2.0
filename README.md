# SVTuner

SVTuner is a spatial-transcriptomics workflow that detects reference populations unsupported by the observed tissue and prevents unsupported mappings from being treated as biological signal. This repository contains the maintained pipeline, the primary simulation and real-data evaluations, the biological application, frozen reproducibility assets, and targeted robustness analyses added during manuscript revision.

The core workflow has four stages:

1. **Stage 1 — preprocessing:** harmonize single-cell and spatial expression inputs.
2. **Stage 3A — reference-level diagnostics:** identify unsupported cell types and filter the single-cell reference.
3. **Stage 3B — spatial-unit diagnostics:** score spatial units, calibrate against pseudo-ST null data, and create a final withheld mask.
4. **Stage 4 — mapping:** run CytoSPACE with the selected reference and, optionally, the Stage3B mask.

> **Repository boundary:** raw data, processed matrices, and the working result/ tree are intentionally excluded from Git. GitHub contains code, configurations, frozen manifests, and curated tables/figures under visualizations/ and reproducibility_release/.

## Quick links

- [Scenario-running guide](RUNNING_SCENARIOS.md)
- [Frozen reproducibility release](reproducibility_release/README.md)
- [C1–C8 provenance audit](C_revision_provenance_audit.md)
- [Environment specification](configs/environment.yml)
- [Dataset configurations](configs/datasets)
- [Pipeline source](src)
- [Research scripts](scripts)
- [Tests](tests)

## Installation

The reproducibility environment is defined in `configs/environment.yml`
and pins Python 3.10.19, R 4.5.2, and the principal scientific
dependencies used by the maintained workflow.

~~~bash
conda env create -f configs/environment.yml
conda activate cytospace_v1.1.0_py310
pip install -e .
svtuner version
svtuner envcheck
~~~

External tools used by specific routes, notably CytoSPACE, must also be available in the active environment. Dataset files are not redistributed automatically.

## Data and storage contract

Each dataset has a YAML file in configs/datasets/. Important fields include the sample name, input filenames, optional storage.group, cell-type columns, and Stage3 parameters.

When storage.group is set, the main paths are:

~~~text
data/raw/<group>/<sample>/
data/processed/<group>/<sample>/
result/<group>/<sample>/
~~~

Without a group, they are:

~~~text
data/raw/<sample>/
data/processed/<sample>/
result/<sample>/
~~~

Stage 1, Stage 3A, and Stage 3B use this group-aware layout. Stage 4 currently writes to:

~~~text
result/<sample>/stage4_cytospace<suffix>/
~~~

This Stage 4 path is documented explicitly because it differs from the group-aware Stage 1–3 layout.

Expected raw inputs are normally declared in the dataset YAML and commonly include:

~~~text
brca_scRNA_GEP.txt
brca_scRNA_celllabels.txt
brca_STdata_GEP.txt
brca_STdata_coordinates.txt
~~~

Do not infer missing inputs from similarly named samples. Use the dataset configuration and provenance files.

## Running the core stages

The detailed and safer command reference is in [RUNNING_SCENARIOS.md](RUNNING_SCENARIOS.md). The minimal pattern is:

~~~bash
# Stage 1
Rscript r_scripts/stage1_preprocess.R \
  --sample <sample> \
  --project_root <project_root> \
  --export_csv

# Stage 3A
python -m src.stages.stage3_type_plugin \
  --project_root <project_root> \
  --sample <sample> \
  --sc_expr_source normalized \
  --output_suffix _trial

# Stage 3B
python -m src.stages.stage3b_st_unsupported \
  --project_root <project_root> \
  --sample <sample> \
  --output_suffix _trial
~~~

For Stage 3B, invoking the module without numerical overrides lets the dataset YAML and program defaults resolve the parameters. The convenience svtuner stage3b command exposes explicit command-line defaults; use it only when those values are intended.

Stage 4 supports four decomposition routes:

| Route | Reference | Spatial mask |
|---|---|---|
| CytoSPACE baseline | Experiment-specific Stage1 reference | None |
| Stage3A-only | Stage3A-admissible reference | None |
| Stage3B-only | Unfiltered Stage1 reference supplied to that route | Standalone Stage3B mask |
| Full | Stage3A-admissible reference | Sequential Stage3B mask |

Standalone and sequential Stage3B masks are different analyses and must not be shared.

Example Stage 4 pattern:

~~~bash
python -m src.stages.stage4_cytospace \
  --project_root <project_root> \
  --sample <sample> \
  --filter_mode none \
  --cell_type_column sc_meta \
  --sc_expr_source normalized \
  --stage4_suffix _baseline
~~~

For the generic Stage3A-only route, use the Stage3A-adjusted
`plugin_type` annotations and exclude `Unknown_sc_only` cells:

~~~bash
python -m src.stages.stage4_cytospace \
  --project_root <project_root> \
  --sample <sample> \
  --filter_mode plugin_unknown \
  --cell_type_column plugin_type \
  --filter_scope unsupported_all \
  --sc_expr_source normalized \
  --stage4_suffix _stage3a_only
~~~

Use --stage3_suffix, --stage3b_blank_regions, and --stage3b_scores_path as appropriate for the Stage3A-only, Stage3B-only, and Full routes. Always use unique suffixes for exploratory work.

The svtuner run convenience command runs Stage 1, Stage 3A, baseline Stage 4, and the Stage3A-controlled Stage 4 route. It does **not** run Stage 3B. In addition, configs/pipeline_presets.yaml is currently empty, so --from-scratch requires a maintained preset before use.

## Evaluation framework

The scientific evaluation is organized around three primary experiment families. The later C1–C8 labels are internal revision-work-package identifiers; they are not biological cohorts, dataset names, or the main experiment hierarchy.

### Primary experiment 1: controlled simulation benchmark

The main simulation benchmark asks whether SVTuner improves mapping when the spatial data and the single-cell reference are incomplete in different ways. Each scenario combines two deliberately distinct challenges:

- **ST-side absence:** one or more cell types remain in the SC reference but are absent from the simulated ST expression. This tests whether Stage3A prevents unsupported reference populations from being mapped.
- **SC-reference dropout:** a cell type contributes to the simulated ST truth but is removed from the reference supplied to mapping. This tests whether Stage3B identifies spatial units that the available reference cannot adequately explain.

The legacy scenario IDs are compact file labels:

- **BC** means the breast-cancer simulation family; **Lung** means human lung; **Brain** means mouse brain.
- Suffix **1** is the ST-side control, suffix **2** removes one ST population, and suffix **3** removes two ST populations.
- All nine scenarios also contain the tissue-specific SC-reference dropout shown in the final column.

| Readable scenario | Legacy ID | ST-side population(s) absent | SC-reference dropout |
|---|---|---|---|
| Breast cancer, control | BC-1 | None | Endothelial cells |
| Breast cancer, single ST absence | BC-2 | Epithelial cells | Endothelial cells |
| Breast cancer, double ST absence | BC-3 | Epithelial cells + PCs | Endothelial cells |
| Human lung, control | Lung-1 | None | B cell |
| Human lung, single ST absence | Lung-2 | AT2 | B cell |
| Human lung, double ST absence | Lung-3 | AT2 + Fibroblast | B cell |
| Mouse brain, control | Brain-1 | None | Ext_L56 |
| Mouse brain, single ST absence | Brain-2 | Micro | Ext_L56 |
| Mouse brain, double ST absence | Brain-3 | Micro + Oligo_2 | Ext_L56 |

Cell-type names follow the corresponding dataset annotations: AT2 denotes alveolar type 2 epithelial cells, PCs denotes plasma cells, and Micro denotes microglia; Ext_L56 and Oligo_2 are the exact labels used by the mouse-brain reference. The primary comparison evaluates SVTuner and mapping baselines with abstention-aware composition recovery. Frozen source values are under [visualizations/method_comparison/composite_no_noise/](visualizations/method_comparison/composite_no_noise/); exact config names and runnable scenario commands are in [RUNNING_SCENARIOS.md](RUNNING_SCENARIOS.md).

### Primary experiment 2: real-data profile-masking validation

The second main experiment family uses profile masking in real spatial datasets. Target-associated profile genes are deliberately attenuated in the spatial expression data while the scRNA-seq reference is kept unchanged. This creates a controlled unsupported-reference challenge without claiming that the edited data are a new biological sample. It tests whether:

1. Stage3A detects and excludes the deliberately masked target population from the admissible reference;
2. standalone and sequential Stage3B masks quantify the incremental spatial abstention before and after the Stage3A reference decision; and
3. the combined route suppresses unsupported target-like assignments under a fixed low-support evaluation set.

The maintained collection contains 15 profile-masking experiments: 10 spot-resolution datasets and 5 MERSCOPE datasets evaluated at capacity 1. Curated benchmark and decomposition assets are under [visualizations/cytospace_fig2d_profile_mask_benchmark/](visualizations/cytospace_fig2d_profile_mask_benchmark/), [visualizations/highres_profile_mask_mapping/](visualizations/highres_profile_mask_mapping/), and [visualizations/c5_decomposition/](visualizations/c5_decomposition/).

### Primary experiment 3: biological application

The biological application evaluates SVTuner against a frozen CTA-derived immune endpoint after validated CTA-to-spot matching. CytoSPACE and SVTuner outputs are compared on the same formal barcode-matched analysis universe. The maintained evidence chain, endpoint freeze, mapping outputs, and final evaluation assets are under [visualizations/bioapp_experiment/](visualizations/bioapp_experiment/).

## Targeted robustness and revision checks

The C1–C8 labels were introduced only to track focused questions raised during revision. They extend or audit the primary experiments above; they do not replace them. Their purpose is to show whether the main conclusions depend on a particular seed, cutoff, null calibration, stage combination, expression source, or reporting choice.

| Internal ID | Why this check was added | Public entry points / artifacts |
|---|---|---|
| C1 | Test whether the nine-scenario simulation conclusion persists across ten simulation seeds rather than one realization | scripts/run_c1_seed_repeat.py; visualizations/method_comparison/c1_multiseed/ |
| C2 | A parameter-sensitivity attempt was considered, but no formal result was released | **CLOSED / CANCELLED.** Temporary diagnostic runs are not manuscript evidence; parameter transparency is handled in the Methods and reproducibility records |
| C3 | Test whether the spatial marker-core conclusion depends on choosing a 15% cutoff by repeating the readout at 10%, 15%, 20%, and 25% | scripts/plot_stage3b_threshold_robustness.py; visualizations/stage3b_realdata_candidate_scan/stage3b_threshold_robustness/ |
| C4 | Audit Stage3B calibration against pseudo-null data, technical perturbations, BH versus BY correction, and 200 versus 1,000 spatial permutations | scripts/export_c4_1_stage3b_raw_calibration_null.py; scripts/run_c4_*; visualizations/c4_stage3b_calibration/ |
| C5 | Decompose the profile-masking results into Stage3A-only, standalone Stage3B, and sequential Stage3A + Stage3B contributions | scripts/run_c5_decomposition.py; visualizations/c5_decomposition/ |
| C6 | Remove a potential shared-source advantage by deriving simulated-ST expression profiles from an independent breast-cancer study while retaining a separate SC reference | C6 dataset configs and generator/evaluation scripts; visualizations/c6_independent_source_evaluation/ |
| C7 | Add the reviewer-requested precision–recall and frozen-binary metrics for the existing CTA application, then repair their deterministic provenance | scripts/repair_c7_metric_provenance.py; visualizations/bioapp_experiment/C7_cta_extended_performance/ |
| C8 | Summarize direct paired method comparisons using the independent experiment, rather than internal readouts, as the inferential unit | C8 paired-analysis/plotting scripts; visualizations/c8_direct_paired_comparison/ |

Examples of the resulting secondary checks include seed-level stability across the nine simulation scenarios, marker-core robustness across four cutoffs, full Stage3A exclusion of the masked target in all 15 profile-masking experiments, and deterministic recovery of the CTA AUROC/AP and binary-decision provenance. Exact values, configurations, input paths, and known limitations are recorded in [C_revision_provenance_audit.md](C_revision_provenance_audit.md).

## Reproducibility release

reproducibility_release/ contains the frozen release-candidate provenance package, including manifests, source-value records, configurations, environment records, audit records, redistribution metadata, file-integrity records, and reconstruction instructions.

The revised-manuscript analyses, including reviewer-response experiments C1–C8, are additionally represented in the maintained repository through the corresponding scripts, configurations, provenance records, and curated outputs described above. The repository root and reproducibility_release/ therefore provide complementary layers of the public reproducibility record.

Its manifest records 497 release files and 482 hashed payload records. Read its [README](reproducibility_release/README.md) before reproducing a frozen analysis.

## What is actually on GitHub?

Tracked and uploaded:

- source code, scripts, tests, and YAML configurations;
- frozen provenance and reproducibility manifests;
- curated publication figures and source tables under visualizations/;
- the contents of reproducibility_release/.

Not tracked or uploaded:

- data/ (raw, simulated, and processed matrices);
- result/ (working and formal local result trees);
- logs/, local environments, caches, and temporary workspaces.

The exclusion is enforced by .gitignore. A local result such as result/c6_independent_source_evaluation/c6_supplementary_table.csv can be scientifically frozen locally while still being absent from GitHub. If an artifact must be published, copy only the curated, redistributable file to an approved tracked location such as visualizations/<analysis>/ or a versioned reproducibility package and document its provenance. Do not force-add the entire result/ tree.

## Repository layout

~~~text
configs/                 project, environment, and dataset YAML files
data/                    local raw/simulated/processed data (Git-ignored)
external/                external integration code retained in the repository
r_scripts/               Stage 1 and supporting R scripts
reproducibility_release/ frozen reproducibility package
result/                  local working/formal outputs (Git-ignored)
scripts/                 experiment, audit, evaluation, and plotting scripts
src/                     maintained SVTuner implementation
tests/                   focused Stage3B/Stage4 tests
visualizations/          curated figures and source tables tracked by Git
~~~

## Validation

Run the focused tests with:

~~~bash
python -m pytest -q tests
~~~

Before committing a result, verify its input provenance, dataset YAML, output suffix, and tracking status:

~~~bash
git check-ignore -v <path>
git ls-files -- <path>
~~~

## Interpretation

SVTuner supports abstention-aware spatial mapping; it does not replace biological validation. “Unsupported” is evidence relative to a reference, expression-processing path, and spatial context. Stage3A diagnoses reference populations, while Stage3B diagnoses spatial units; their outputs are not interchangeable.

## License

The package metadata currently declares a proprietary license. Third-party datasets and software remain subject to their original licenses and redistribution terms.
