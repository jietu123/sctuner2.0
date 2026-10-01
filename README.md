# SVTuner

SVTuner is a spatial-transcriptomics workflow that detects reference populations unsupported by the observed tissue and prevents unsupported mappings from being treated as biological signal. This repository contains the maintained pipeline, dataset configurations, frozen reproducibility assets, and reviewer-response experiments C1–C8.

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
  --filter_scope none \
  --sc_expr_source normalized \
  --stage4_suffix _baseline
~~~

Use --stage3_suffix, --stage3b_blank_regions, and --stage3b_scores_path as appropriate for the Stage3A-only, Stage3B-only, and Full routes. Always use unique suffixes for exploratory work.

The svtuner run convenience command runs Stage 1, Stage 3A, baseline Stage 4, and the Stage3A-controlled Stage 4 route. It does **not** run Stage 3B. In addition, configs/pipeline_presets.yaml is currently empty, so --from-scratch requires a maintained preset before use.

## Composite simulation scenarios

The frozen no-noise composite suite contains nine scenarios:

| ID | Tissue | ST missing population(s) | SC reference dropout |
|---|---|---|---|
| BC-1 | Breast cancer | None | Endothelial cells |
| BC-2 | Breast cancer | Epithelial cells | Endothelial cells |
| BC-3 | Breast cancer | Epithelial cells + PCs | Endothelial cells |
| Lung-1 | Human lung | None | B cell |
| Lung-2 | Human lung | AT2 | B cell |
| Lung-3 | Human lung | AT2 + Fibroblast | B cell |
| Brain-1 | Mouse brain | None | Ext_L56 |
| Brain-2 | Mouse brain | Micro | Ext_L56 |
| Brain-3 | Mouse brain | Micro + Oligo_2 | Ext_L56 |

The exact config names and the C1 command are listed in [RUNNING_SCENARIOS.md](RUNNING_SCENARIOS.md).

## Revision experiments C1–C8

| Experiment | Purpose | Main entry points / public artifacts |
|---|---|---|
| C1 | Ten-seed robustness across nine composite scenarios | scripts/run_c1_seed_repeat.py; visualizations/method_comparison/c1_multiseed/ |
| C2 | Stage3A parameter-sensitivity attempt | **CLOSED / CANCELLED.** No formal sensitivity result is released; temporary diagnostic runs are not manuscript evidence. Parameter transparency is handled in the Methods and reproducibility records. |
| C3 | Marker-core robustness at 10%, 15%, 20%, and 25% | scripts/plot_stage3b_threshold_robustness.py; visualizations/stage3b_realdata_candidate_scan/stage3b_threshold_robustness/ |
| C4 | Stage3B null calibration, technical perturbation, and statistical robustness | scripts/export_c4_1_stage3b_raw_calibration_null.py and the run_c4_* scripts; visualizations/c4_stage3b_calibration/ |
| C5 | Stage3A/Stage3B decomposition in 15 profile-masking experiments | scripts/run_c5_decomposition.py; visualizations/c5_decomposition/ |
| C6 | Independent-source simulated-ST profiles | C6 configs plus generator, dropout, and evaluation scripts |
| C7 | CTA performance extension and deterministic provenance repair | scripts/repair_c7_metric_provenance.py; visualizations/bioapp_experiment/C7_cta_extended_performance/ |
| C8 | Direct paired comparison at the independent experiment level | C8 paired-analysis/plotting scripts; visualizations/c8_direct_paired_comparison/ |

Selected frozen results:

- **C1:** nine scenarios × ten simulation seeds; the overall seed-level mean paired difference was approximately +0.0975, with all ten overall seed means positive.
- **C3:** at the original 15% marker-core threshold, mean marker-associated fraction was 93.79% and mean core withheld rate was 57.78% across six scenarios.
- **C5:** all 15 masked targets were detected and fully excluded by Stage3A; melanoma patient 2 also had the recorded NK-cell off-target exclusion.
- **C7:** formal matched universe n=1888 (134 CTA-positive, 1754 CTA-negative); AUROC 0.8305238346466074 and AP 0.28365511700029472. The frozen binary decision is not exactly reproducible by a one-dimensional withheld_score threshold; the minimum disagreement is 93.

For exact configurations, input paths, known limitations, and frozen-vs-local distinctions, see [C_revision_provenance_audit.md](C_revision_provenance_audit.md).

## Independent-source simulation (C6)

C6 separates the single-cell reference source from the source used to derive simulated-ST cell-type profiles:

~~~text
SC reference / composition source: human_breast_cancer_real
Independent profile source:       Wu et al. GSE176078 HER2+ breast cancer
Harmonized broad types:           5
Three-way shared genes:           16,401
Seed / mix_alpha / depth_scale:   42 / 1.0 / 1.0
~~~

The generator adds the optional --profile_source_sample interface. Without it, the historical generator path is preserved. C6-generated data and result files remain local because data/ and result/ are ignored.

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
