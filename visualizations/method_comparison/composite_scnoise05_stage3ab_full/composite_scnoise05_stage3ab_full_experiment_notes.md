# 5% SC-reference noise: full Stage3A + Stage3B benchmark

## 1. Experiment status

This directory contains the complete 5% single-cell-reference noise benchmark that mirrors the existing 10% experiment. The benchmark includes nine simulation scenarios, seven mapping routes, the full SVTuner Stage3A + Stage3B route, source-value tables, provenance audits, and the final boxplot.

The frozen quality decision is **PARTIAL**, not PASS. All execution and provenance guardrails passed, but Stage3A did not detect every expected unsupported reference type in the human-lung scenarios. This is an observed method limitation and was not corrected with truth labels or a whitelist.

## 2. Noise implementation

The noise is implemented as independent entry-wise within-gene replacement with:

- requested probability: `p = 0.05`;
- random seed: `42`;
- perturbed object: the SC-reference expression matrices only;
- replacement rule: each selected matrix entry is replaced by the value of the same gene from a randomly selected donor cell;
- unchanged inputs: ST expression, spatial coordinates, and simulation truth.

This is not 5% additive-amplitude noise, cell-level corruption, gene-level corruption, or dropout. Separate Bernoulli masks are generated for the normalized, data, and counts matrices.

Observed perturbed-entry fractions were:

| Dataset group | Normalized | Data | Counts |
|---|---:|---:|---:|
| Real BRCA | 0.049995 | 0.049967 | 0.049940 |
| Human lung 5-location | 0.049981 | 0.050112 | 0.050004 |
| Mouse brain refined | 0.049976 | 0.049924 | 0.050040 |

The input audit contains `36` rows. All `9` noisy SC normalized matrices differ from their no-noise sources, while all `27` ST-expression, coordinate, and truth inputs are byte-identical to their corresponding no-noise sources.

## 3. Experimental design

The benchmark uses three dataset groups and three conditions per group:

| Group | Control | Single Stage3A missing type | Double Stage3A missing types | Stage3B unsupported ST target |
|---|---|---|---|---|
| Real BRCA | none | Epithelial cells | Epithelial cells + PCs | Endothelial cells |
| Human lung 5-location | none | AT2 | AT2 + Fibroblast | B_cell |
| Mouse brain refined | none | Micro | Micro + Oligo_2 | Ext_L56 |

Control scenarios use `filter_mode=none` and remove zero reference cells. Missing-reference scenarios use:

```text
filter_mode = plugin_unknown
filter_scope = unsupported_all
cell_type_column = plugin_type
truth_filter = false
```

The complete SVTuner route is:

```text
5% noisy SC reference
-> Stage3A unsupported-reference detection
-> filter plugin_unknown reference cells
-> Stage3B unsupported-ST-region detection
-> CytoSPACE mapping
-> Stage3B blank-mask integration
```

The seven compared routes are CytoSPACE, SVTuner, Tangram (all genes), Tangram (marker genes), novoSpaRc, SpaOTsc, and CellTrek. Each method contributes one result for each of the nine scenarios, yielding `63` scenario-level values. The cell-type-level source table contains `420` rows.

CytoSPACE and the five external-method outputs were generated from the same 5% noisy inputs. The `54 = 9 scenarios x 6 non-SVTuner routes` output-provenance rows all passed file-existence, spot-universe, timestamp, deterministic-noise-recreation, and input-hash checks; `rerun_required=False` for all rows.

## 4. Evaluation metric

The figure uses whole-space target-dominant abstention-aware composition overlap.

For a nonblank spot, predicted and truth cell-type compositions are normalized over the complete cell-type union and scored as:

```text
spot overlap = sum_c min(predicted_fraction_c, truth_fraction_c)
```

For an SVTuner-blanked spot:

- score `1` when the Stage3B target is the truth-dominant type;
- score `0` otherwise.

The scenario score is the mean over all evaluated spots. Thus a correct SVTuner abstention is not incorrectly counted as a zero, while an incorrect abstention receives no reward. Other methods have no abstention state and are scored by ordinary whole-space composition overlap.

## 5. Main results

| Method | n | Mean | Median | SD |
|---|---:|---:|---:|---:|
| CytoSPACE | 9 | 0.505557 | 0.519654 | 0.103455 |
| SVTuner | 9 | **0.659712** | **0.681771** | 0.159112 |
| Tangram (all genes) | 9 | 0.448806 | 0.457472 | 0.115926 |
| Tangram (marker genes) | 9 | 0.466261 | 0.455603 | 0.124594 |
| novoSpaRc | 9 | 0.489779 | 0.472817 | 0.102491 |
| SpaOTsc | 9 | 0.385585 | 0.389270 | 0.121520 |
| CellTrek | 9 | 0.459861 | 0.486228 | 0.109603 |

Relative to CytoSPACE, SVTuner showed:

- absolute mean improvement: `0.154155`;
- relative mean improvement: `30.4922%`;
- paired wins: `7/9`;
- paired losses: `2/9`;
- pooled blank precision: `0.957936`;
- pooled blank recall: `0.903398`.

SVTuner comparisons with all other methods were:

| Comparison | Mean difference | Wins | Losses | Ties |
|---|---:|---:|---:|---:|
| vs CytoSPACE | 0.154155 | 7 | 2 | 0 |
| vs Tangram (all genes) | 0.210906 | 8 | 1 | 0 |
| vs Tangram (marker genes) | 0.193451 | 8 | 1 | 0 |
| vs novoSpaRc | 0.169933 | 8 | 1 | 0 |
| vs SpaOTsc | 0.274128 | 9 | 0 | 0 |
| vs CellTrek | 0.199851 | 6 | 3 | 0 |

## 6. SVTuner versus CytoSPACE by scenario

| Group | Condition | CytoSPACE | SVTuner | SVTuner - CytoSPACE |
|---|---|---:|---:|---:|
| Real BRCA | control | 0.652655 | 0.769384 | 0.116729 |
| Real BRCA | single missing | 0.519654 | 0.681771 | 0.162117 |
| Real BRCA | double missing | 0.363320 | 0.554192 | 0.190872 |
| Human lung | control | 0.588419 | 0.587772 | -0.000647 |
| Human lung | single missing | 0.484215 | 0.482535 | -0.001680 |
| Human lung | double missing | 0.365083 | 0.407648 | 0.042566 |
| Mouse brain | control | 0.579690 | 0.826086 | 0.246396 |
| Mouse brain | single missing | 0.572331 | 0.855748 | 0.283416 |
| Mouse brain | double missing | 0.424647 | 0.772276 | 0.347629 |

The two paired losses are both small human-lung differences. They must remain in the reported `7/9` result and must not be removed as outliers.

## 7. Stage3A audit and PARTIAL decision

| Group | Condition | Expected types | Detected types | Unknown/filtered cells | Detection status |
|---|---|---|---|---:|---|
| Real BRCA | control | none | none | 0 | not applicable |
| Real BRCA | single missing | Epithelial cells | Epithelial cells | 500 | detected |
| Real BRCA | double missing | Epithelial cells; PCs | Epithelial cells; PCs | 991 | detected |
| Human lung | control | none | none | 0 | not applicable |
| Human lung | single missing | AT2 | none | 0 | detection failure |
| Human lung | double missing | AT2; Fibroblast | AT2 | 2779 | partial detection |
| Mouse brain | control | none | none | 0 | not applicable |
| Mouse brain | single missing | Micro | Micro | 1000 | detected |
| Mouse brain | double missing | Micro; Oligo_2 | Micro; Oligo_2 | 2000 | detected |

All critical execution guardrails passed: controls used no filtering, missing scenarios used `plugin_unknown`, Stage3A outputs and logs were present, `unsupported_all` was used, truth filtering was disabled, all 63 method-scenario values were present, and Stage3B masks were integrated in all nine Stage4 results.

Two quality-only guardrails failed because human-lung single missing filtered `0` cells and not all expected human-lung types were detected. Therefore the final decision is `PARTIAL`. The result supports a qualified robustness statement but not a claim that Stage3A recovered every 5%-noise missing-reference perturbation.

## 8. Stage3B abstention audit

| Group/condition | Predicted blanks | Correct blanks | Incorrect blanks | Truth unsupported | Precision | Recall |
|---|---:|---:|---:|---:|---:|---:|
| BRCA control | 382 | 382 | 0 | 384 | 1.000000 | 0.994792 |
| BRCA single missing | 382 | 382 | 0 | 384 | 1.000000 | 0.994792 |
| BRCA double missing | 382 | 382 | 0 | 384 | 1.000000 | 0.994792 |
| Lung control | 194 | 128 | 66 | 385 | 0.659794 | 0.332468 |
| Lung single missing | 28 | 22 | 6 | 397 | 0.785714 | 0.055416 |
| Lung double missing | 8 | 5 | 3 | 397 | 0.625000 | 0.012594 |
| Brain control | 2933 | 2843 | 90 | 2853 | 0.969315 | 0.996495 |
| Brain single missing | 2918 | 2843 | 75 | 2853 | 0.974297 | 0.996495 |
| Brain double missing | 3043 | 2851 | 192 | 2853 | 0.936904 | 0.999299 |

The large cross-dataset variation in abstention recall is visible in the audit and must not be hidden by the pooled `0.903398` value.

## 9. Figure and authoritative files

Main figure:

```text
fig1d_scnoise05_stage3ab_full.png
fig1d_scnoise05_stage3ab_full.pdf
```

Authoritative data and audit files:

```text
composite_scnoise05_stage3ab_full_source_values.csv
composite_scnoise05_stage3ab_full_cell_type_values.csv
composite_scnoise05_stage3ab_full_method_summary.csv
composite_scnoise05_stage3ab_full_paired_comparison.csv
composite_scnoise05_stage3ab_full_svtuner_method_comparisons.csv
composite_scnoise05_stage3ab_full_blank_summary.csv
composite_scnoise05_stage3ab_full_stage3a_audit.csv
composite_scnoise05_stage3ab_full_stage3b_audit.csv
composite_scnoise05_stage3ab_full_input_audit.csv
composite_scnoise05_stage3ab_full_reused_method_audit.csv
composite_scnoise05_stage3ab_full_scenario_manifest.csv
composite_scnoise05_stage3ab_full_guardrails.json
composite_scnoise05_stage3ab_full_summary.json
```

## 10. Reproduction

The 5% wrapper calls the shared percentage-configurable runner. On the machine used for this run:

```powershell
E:\ANACONDA\envs\cytospace_v1.1.0_py310\python.exe `
  scripts\run_composite_scnoise05_stage3ab_full.py `
  --project_root . `
  --python E:\ANACONDA\envs\cytospace_v1.1.0_py310\python.exe `
  --external_python E:\tmp\conda_envs\svtuner_multimethod_py310\python.exe `
  --prepare_inputs `
  --run_baselines `
  --run_external `
  --execute `
  --resume `
  --n_processors 1 `
  --n_subspots 800 `
  --mapping_cells_per_spot 5
```

On Windows paths containing non-ASCII characters, a temporary ASCII `subst` path may be required by isolated Conda/R environments. This does not alter data identity; the hashes in the audit tables are authoritative.
