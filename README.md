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
- [Revision provenance audit](C_revision_provenance_audit.md)
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

## Experiments and evidence

The evaluation follows the two failure modes addressed by SVTuner. Stage3A tests whether a cell type present in the single-cell (SC) reference is supported by spatial transcriptomics (ST). Stage3B tests whether a spatial unit can be explained by the reference that is available for mapping. The experiments examine these decisions in controlled simulations, real tissues, multiple spatial resolutions, downstream analyses, and a biological application.

The frozen [experiment manifest](reproducibility_release/manifests/experiment_manifest_v1.0.tsv) indexes 40 experiment records in the following families. It describes the release-candidate evidence; later analyses are linked below alongside the relevant family.

| Experiment family | Records in frozen manifest | Main question | Starting point |
|---|---:|---|---|
| Joint Stage3A–Stage3B simulations | 9 | Can the workflow handle unsupported reference types and missing reference types in the same spatial sample? | [Simulation comparison](visualizations/method_comparison/composite_no_noise/) |
| Spot-resolution profile masking | 10 | Does Stage3A reject a target profile deliberately attenuated in real ST expression? | [Profile-masking experiments](visualizations/masked_scenarios/real/) |
| MERSCOPE profile masking | 5 | Do the reference and spatial-unit decisions remain interpretable at single-cell spatial resolution and mapping capacity 1? | [High-resolution mapping](visualizations/highres_profile_mask_mapping/) |
| External tumour comparisons | 6 | Do profile-masking effects extend to independent melanoma, breast-cancer and colorectal-cancer samples and downstream state analyses? | [External profile-mask benchmark](visualizations/cytospace_fig2d_profile_mask_benchmark/) |
| Real-data reference dropout | 6 | Does Stage3B flag spatial regions when a population is removed from the SC reference while ST remains? | [Reference-missing analyses](visualizations/stage3b_realdata_candidate_scan/) |
| Detailed thalamic-excitatory case | 1 | Where does reference absence cause forced assignment in a mouse-brain spatial dataset? | [Thalamic case](visualizations/cell2location_stage3b_case/) |
| State-decoy tests | 2 | How does Stage3A respond to unsupported or misleading reference cell states? | [Kidney decoy](visualizations/cytospace_fig2i_mouse_kidney_stage3_unsupported_decoy_sensitivity/) and [breast-cancer state decoy](visualizations/cytospace_fig2k_tcell_states_stage3_decoy/) |
| CTA biological application | 1 | Does the mapped output agree with an independently frozen CTA-derived immune endpoint? | [Biological application](visualizations/bioapp_experiment/) |

### Controlled simulation and method comparison

The joint simulation crosses two kinds of mismatch. In the Stage3A condition, a cell type remains in the SC reference but its profile is absent from simulated ST. In the Stage3B condition, a cell type remains in ST simulation truth but is removed from the SC reference supplied to mapping. Breast-cancer, human-lung and mouse-brain backgrounds each have a control, a single ST-side missing type and a double ST-side missing type, yielding nine scenarios. The short IDs in data files are **BC** (breast cancer), **Lung**, and **Brain**, with suffixes **1**, **2**, and **3** for those three conditions. The SC-reference dropout is present across all three conditions within each tissue.

The primary benchmark compares abstention-aware composition recovery across SVTuner, CytoSPACE and other mapping methods. It also examines added SC-reference noise. Source values and experiment notes are in [method_comparison/](visualizations/method_comparison/); [RUNNING_SCENARIOS.md](RUNNING_SCENARIOS.md) maps all nine IDs to exact dataset configs. Repeated simulation seeds and an independent breast-cancer source for simulated-ST expression profiles test how much these results depend on a single realization or a shared expression source. Their curated data are in [multiseed](visualizations/method_comparison/c1_multiseed/) and [independent-source evaluation](visualizations/c6_independent_source_evaluation/).

### Real-tissue validation and mapping consequences

Profile masking attenuates target-associated genes in real ST expression while retaining the SC reference. This controlled perturbation tests Stage3A's exclusion decision and the effect on mapping. The spot-resolution experiments span ten datasets; the five MERSCOPE experiments examine the same question at capacity 1. The later decomposition compares baseline, Stage3A-only, standalone Stage3B and sequential Stage3A + Stage3B routes with the appropriate reference and mask for each route. See the [profile-mask benchmark](visualizations/cytospace_fig2d_profile_mask_benchmark/), [MERSCOPE analyses](visualizations/highres_profile_mask_mapping/) and [stage-decomposition figures and source tables](visualizations/c5_decomposition/).

Complementary reference-dropout experiments remove a population from the SC reference while keeping the spatial measurements. They probe Stage3B's marker-core localization, spatial withholding and failure modes, including a detailed thalamic-excitatory example. The [real-data candidate scan](visualizations/stage3b_realdata_candidate_scan/), [reference-missing stress tests](visualizations/stage3b_reference_missing_stress/) and [thalamic case](visualizations/cell2location_stage3b_case/) provide the associated evidence. Calibration against pseudo-null data, technical perturbations, alternative multiple-testing correction and spatial-permutation depth were examined in [Stage3B robustness analyses](visualizations/c4_stage3b_calibration/). External tumour comparisons and state-decoy tests assess how these decisions affect assignments and downstream readouts in other biological settings.

### Biological application and cross-experiment comparison

The biological application uses a frozen CTA-derived immune endpoint after validated CTA-to-spot matching. CytoSPACE and SVTuner are compared on the same barcode-matched analysis set. The [application archive](visualizations/bioapp_experiment/) contains the endpoint, mapping and evaluation trail; the [precision–recall supplement](visualizations/bioapp_experiment/C7_cta_extended_performance/) adds operating-point and continuous-score readouts. A [direct paired comparison](visualizations/c8_direct_paired_comparison/) summarizes selected mapping and downstream outcomes at the independent-experiment level.

For exact sample names, configurations, source-value files and known provenance limits, use the frozen [experiment manifest](reproducibility_release/manifests/experiment_manifest_v1.0.tsv) together with the [revision provenance audit](C_revision_provenance_audit.md). Internal B/C task labels appear in some filenames and audit records; they are work-package identifiers, not the structure of the scientific evaluation. The attempted Stage3A parameter-sensitivity work package was closed without a formal result and is not part of the reported evidence.

## Reproducibility release

reproducibility_release/ contains the frozen release-candidate provenance package, including manifests, source-value records, configurations, environment records, audit records, redistribution metadata, file-integrity records, and reconstruction instructions.

The later analyses described above are represented in the maintained repository through their scripts, configurations, provenance records, and curated outputs. The repository root and reproducibility_release/ provide complementary layers of the public reproducibility record.

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
