# C revision provenance and configuration audit

## 1. Audit scope

- Read-only audit of completed C1/C3/C4/C5/C6/C7/C8 revision artifacts. C2 is excluded by instruction.
- No Stage1, Stage3A, Stage3B, CytoSPACE, simulation, evaluator, or statistical analysis pipeline was rerun.
- Evidence hierarchy used: persisted resolved output/log > persisted run manifest/provenance > persisted explicit config with direct run linkage > source-value/result table. Current defaults and merely similar configs were not used to fill historical fields.
- Missing values are written exactly as `NOT DOCUMENTED / NOT RECOVERABLE`.
- File timestamps below are filesystem LastWriteTime evidence, not embedded run timestamps unless explicitly stated.

## 2. C5 exact Stage3A configuration

- C5 reused the pre-existing unsuffixed `stage3_typematch` outputs for the 15 profile-masking experiments; C5 did not rerun Stage3A.
- Authoritative per-experiment resolved parameters are the `params` objects in each formal `stage3_summary.json`.
- Fourteen current dataset YAMLs match representative resolved values, but the historical command/config linkage is not recorded. Therefore `eps`, `plugin_genes_path`, and `gene_weights_path`, which are absent from the resolved summary, remain unrecovered in the TSV.
- `mouse_embryo_real_profile_mask_erythroid` has a Stage3 summary but no current dataset YAML.
- Exact Stage3A command/launcher, working directory, git commit, and historical code hash are not persisted for these 15 historical runs.
- `src/stages/stage3_type_plugin.py` is the identifiable module implementation, but the exact historical checkout used by each pre-existing output is not recorded.

## 3. C5 per-experiment parameter table

| ID | dataset | target | resolution | capacity | C5 role | matching config | current config SHA-256 | Stage3 summary LastWriteTime | resolved result | status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| C5-01 | adult_mouse_kidney_real_profile_mask_endo | Endo | spot-resolution | 5 | Stage3A exclusion; Stage3B standalone/sequential | configs/datasets/adult_mouse_kidney_real_profile_mask_endo.yaml | 4d4295b1515af070a5f1e01a14ee521c4a8cbdf1b4775c649cdcdb970e670ed0 | 2026-05-29T16:22:49.100799+08:00 | result/low_resolution_experiments/adult_mouse_kidney_real_profile_mask_endo/stage3_typematch/stage3_summary.json | PARTIAL |
| C5-02 | ffpe_mouse_brain_sagittal_real_profile_mask_microglia | Microglia | spot-resolution | 1 | Stage3A exclusion; Stage3B standalone/sequential | configs/datasets/ffpe_mouse_brain_sagittal_real_profile_mask_microglia.yaml | 840a3f427373700b55a6c4962ec9684e3554d64bdcb51269191fee6325c2b894 | 2026-05-29T16:26:06.278985+08:00 | result/low_resolution_experiments/ffpe_mouse_brain_sagittal_real_profile_mask_microglia/stage3_typematch/stage3_summary.json | PARTIAL |
| C5-03 | human_breast_cancer_real_profile_mask_basal_cell | Basal cell | spot-resolution | 5 | Stage3A exclusion; Stage3B standalone/sequential | configs/datasets/human_breast_cancer_real_profile_mask_basal_cell.yaml | 4da98bf69ba47f72a136272bb46f193f3c50b24f6981aefa46c902f9077700ba | 2026-05-07T13:28:58.870325+08:00 | result/low_resolution_experiments/human_breast_cancer_real_profile_mask_basal_cell/stage3_typematch/stage3_summary.json | PARTIAL |
| C5-04 | human_breast_cancer_visium_ff_wta_real_profile_mask_macrophage | Macrophage | spot-resolution | 1 | Stage3A exclusion; Stage3B standalone/sequential | configs/datasets/human_breast_cancer_visium_ff_wta_real_profile_mask_macrophage.yaml | 468d429005ce3d3b7933765c5404a5b56c971501843c4b488011be60a3e6b532 | 2026-05-02T14:53:54.676223+08:00 | result/low_resolution_experiments/human_breast_cancer_visium_ff_wta_real_profile_mask_macrophage/stage3_typematch/stage3_summary.json | PARTIAL |
| C5-05 | human_breast_cancer_wta_120_real_profile_mask_endothelial_cell | Endothelial cell | spot-resolution | 1 | Stage3A exclusion; Stage3B standalone/sequential | configs/datasets/human_breast_cancer_wta_120_real_profile_mask_endothelial_cell.yaml | a202c8fd3c6847d480e4d7c8976a48059c01f04aff045d16861807f82d61dc02 | 2026-05-02T14:48:38.248650+08:00 | result/low_resolution_experiments/human_breast_cancer_wta_120_real_profile_mask_endothelial_cell/stage3_typematch/stage3_summary.json | PARTIAL |
| C5-06 | human_cervical_cancer_real_profile_mask_epithelial_cell | Epithelial cell | spot-resolution | 1 | Stage3A exclusion; Stage3B standalone/sequential | configs/datasets/human_cervical_cancer_real_profile_mask_epithelial_cell.yaml | aa8a6183de7d8739707428ba5235d87a7b3dc873364c88a64b183a50e488e757 | 2026-05-06T21:22:37.832436+08:00 | result/low_resolution_experiments/human_cervical_cancer_real_profile_mask_epithelial_cell/stage3_typematch/stage3_summary.json | PARTIAL |
| C5-07 | human_heart_ff_real_profile_mask_endothelial_cell | Endothelial cell | spot-resolution | 1 | Stage3A exclusion; Stage3B standalone/sequential | configs/datasets/human_heart_ff_real_profile_mask_endothelial_cell.yaml | 39b77470d07ea1fd2d7aec90219995a488eceb72b4ed36d8a693cd1e0f57ad2b | 2026-05-06T21:26:45.452467+08:00 | result/low_resolution_experiments/human_heart_ff_real_profile_mask_endothelial_cell/stage3_typematch/stage3_summary.json | PARTIAL |
| C5-08 | human_intestine_cancer_real_profile_mask_endothelial_cell | Endothelial cell | spot-resolution | 1 | Stage3A exclusion; Stage3B standalone/sequential | configs/datasets/human_intestine_cancer_real_profile_mask_endothelial_cell.yaml | ee4f475f9fc9121d040c5939fe855cb81838ec60fbca9c66b942a9c3594ea8ba | 2026-05-06T21:30:28.824230+08:00 | result/low_resolution_experiments/human_intestine_cancer_real_profile_mask_endothelial_cell/stage3_typematch/stage3_summary.json | PARTIAL |
| C5-09 | human_lymph_node_real_profile_mask_b_cell | B cell | spot-resolution | 1 | Stage3A exclusion; Stage3B standalone/sequential | configs/datasets/human_lymph_node_real_profile_mask_b_cell.yaml | 61ea1706de57252b5ec3bd59ac085f4f773f5858952492690eaa6f5449be35eb | 2026-04-30T20:12:27.316407+08:00 | result/low_resolution_experiments/human_lymph_node_real_profile_mask_b_cell/stage3_typematch/stage3_summary.json | PARTIAL |
| C5-10 | mouse_embryo_real_profile_mask_erythroid | Erythroid | spot-resolution | 1 | Stage3A exclusion; Stage3B standalone/sequential | NOT DOCUMENTED / NOT RECOVERABLE | NOT DOCUMENTED / NOT RECOVERABLE | 2026-05-04T19:53:00.249511+08:00 | result/low_resolution_experiments/mouse_embryo_real_profile_mask_erythroid/stage3_typematch/stage3_summary.json | PARTIAL |
| C5-11 | highres_humanbreastcancerpatient1_profile_mask_monocytes_and_macrophages | Monocytes and Macrophages | MERSCOPE | 1 | Stage3A exclusion; Stage3B standalone/sequential; MERSCOPE four-route ablation | configs/datasets/highres_humanbreastcancerpatient1_profile_mask_monocytes_and_macrophages.yaml | b424c6db3545f25b4a72cb201b8fb47531064f144918f9b6096f49e23a1c98c3 | 2026-05-29T16:26:13.588971+08:00 | result/highres_humanbreastcancerpatient1_profile_mask_monocytes_and_macrophages/stage3_typematch/stage3_summary.json | PARTIAL |
| C5-12 | highres_humancoloncancerpatient1_profile_mask_fibroblasts | Fibroblasts | MERSCOPE | 1 | Stage3A exclusion; Stage3B standalone/sequential; MERSCOPE four-route ablation | configs/datasets/highres_humancoloncancerpatient1_profile_mask_fibroblasts.yaml | 4a7d5518c9cdf68d4e804fbe38ec964ea35c453b8382c4a892b5ddca1b44b477 | 2026-05-27T16:15:43.573645+08:00 | result/highres_humancoloncancerpatient1_profile_mask_fibroblasts/stage3_typematch/stage3_summary.json | PARTIAL |
| C5-13 | highres_humanlungcancerpatient1_profile_mask_plasma_cells | Plasma cells | MERSCOPE | 1 | Stage3A exclusion; Stage3B standalone/sequential; MERSCOPE four-route ablation | configs/datasets/highres_humanlungcancerpatient1_profile_mask_plasma_cells.yaml | e27e1ab8ca37dca8c5df206cbd0987b04e134befef3e1e9f57dacdc94e7973a0 | 2026-05-27T16:16:25.841165+08:00 | result/highres_humanlungcancerpatient1_profile_mask_plasma_cells/stage3_typematch/stage3_summary.json | PARTIAL |
| C5-14 | highres_humanmelanomapatient1_profile_mask_fibroblasts | Fibroblasts | MERSCOPE | 1 | Stage3A exclusion; Stage3B standalone/sequential; MERSCOPE four-route ablation | configs/datasets/highres_humanmelanomapatient1_profile_mask_fibroblasts.yaml | 445100b9d31b5baad3381dee02220a89d9a4a97f6464e6b8204c2afa2526e800 | 2026-05-27T16:17:11.629009+08:00 | result/highres_humanmelanomapatient1_profile_mask_fibroblasts/stage3_typematch/stage3_summary.json | PARTIAL |
| C5-15 | highres_humanmelanomapatient2_profile_mask_b_cells | B cells | MERSCOPE | 1 | Stage3A exclusion; Stage3B standalone/sequential; MERSCOPE four-route ablation | configs/datasets/highres_humanmelanomapatient2_profile_mask_b_cells.yaml | f2018a429c020cdf1eaeed2aa400967bccf9a824e0482fa5bacac7f3e8c3132e | 2026-05-27T16:17:44.204580+08:00 | result/highres_humanmelanomapatient2_profile_mask_b_cells/stage3_typematch/stage3_summary.json | PARTIAL |

Full resolved values are in `C5_stage3a_resolved_parameters.tsv`. `PARTIAL` means the resolved output is available but the exact historical invocation/commit/config binding is not.

## 4. C5 parameter-family summary

| Parameter | Unique resolved values and experiment IDs | Direct source |
| --- | --- | --- |
| strong_th | `0.7`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10, C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| weak_th | `0.4`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10, C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| st_cluster_k | `30`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10<br>`0`: C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| unknown_floor | `0.3`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10, C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| min_cells_rare_type | `20`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10, C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| auto_detection_enable | `true`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10, C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| auto_method | `adaptive_low_support`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10, C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| auto_min_cells | `50`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10, C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| robust_z_th | `-2.0`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10<br>`-1.8`: C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| soft_z_th | `-0.9`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10<br>`-0.75`: C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| require_masked_for_soft | `false`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10, C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| require_masked_for_hard | `true`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10<br>`false`: C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| max_fraction_types | `0.3`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10<br>`0.35`: C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| auto_max_types | `2`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10, C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| auto_action | `mark_unknown`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10, C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| require_confirmation | `true`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10, C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| confirmation_max_support_score | `null (resolved output)`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10<br>`0.75`: C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| confirmation_use_masked_missing | `true`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10, C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| confirmation_use_marker_identity | `true`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10, C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| confirmation_marker_identity_z_th | `-1.5`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10<br>`0.1`: C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| confirmation_marker_support_score_th | `0.5`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10<br>`0.65`: C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| confirmation_marker_max_support_score | `null (resolved output)`: C5-01, C5-02, C5-03, C5-06, C5-07, C5-08, C5-11, C5-12, C5-13, C5-14, C5-15<br>`NOT DOCUMENTED / NOT RECOVERABLE`: C5-04, C5-05, C5-09, C5-10 | stage3_summary.json::params (resolved output) |
| masked_detection_enable | `true`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10, C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| apply_to_auto_missing | `true`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10, C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| neighbor_cosine_th | `0.9`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10<br>`0.86`: C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| neighbor_cell_ratio_min | `1.0`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10<br>`0.35`: C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| marker_top_n | `20`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10, C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| min_identity_markers | `3`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10, C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| min_marker_specificity | `1.2`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10<br>`1.05`: C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| min_marker_type_mean | `0.001`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10, C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| min_marker_st_detect_frac | `0.005`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10<br>`0.001`: C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| st_presence_quantile | `0.9`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10, C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| identity_z_th | `-0.8`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10<br>`-0.55`: C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| pressure_z_th | `1.0`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10<br>`0.0`: C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| min_support_score_for_apply | `0.8`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10<br>`0.7`: C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| masked_max_types | `2`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10, C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| marker_diag_enable | `true`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10, C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| marker_diag_marker_top_n | `80`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10, C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| marker_diag_min_identity_markers | `5`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10<br>`3`: C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| marker_diag_min_all_specificity | `1.3`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10<br>`1.05`: C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| marker_diag_min_neighbor_specificity | `1.1`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10<br>`1.02`: C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| marker_diag_min_marker_type_mean | `0.001`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10, C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| marker_diag_min_marker_st_detect_frac | `0.005`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10<br>`0.001`: C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| marker_diag_st_presence_quantile | `0.9`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10, C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |
| marker_diag_depleted_z_th | `-0.9`: C5-01, C5-02, C5-03, C5-04, C5-05, C5-06, C5-07, C5-08, C5-09, C5-10<br>`-0.55`: C5-11, C5-12, C5-13, C5-14, C5-15 | stage3_summary.json::params (resolved output) |

### Family classification

- Common to all 15 resolved outputs: `strong_th=0.7`, `weak_th=0.4`, `unknown_floor=0.3`, `min_cells_rare_type=20`, auto-missing enabled with `adaptive_low_support`, `auto_min_cells=50`, confirmation enabled, masked-missing enabled, marker-identity diagnostics enabled, and the shared values enumerated in the table above.
- Spot-resolution family (C5-01–C5-10): `st_cluster_k=30`, `robust_z_th=-2.0`, `soft_z_th=-0.9`, `require_masked_for_hard=true`, `max_fraction_types=0.3`, `neighbor_cosine_th=0.9`, `neighbor_cell_ratio_min=1.0`, `identity_z_th=-0.8`, `pressure_z_th=1.0`, `min_support_score_for_apply=0.8`, and the spot-family specificity/detection values shown above.
- MERSCOPE family (C5-11–C5-15): `st_cluster_k=0`, `robust_z_th=-1.8`, `soft_z_th=-0.75`, `require_masked_for_hard=false`, `max_fraction_types=0.35`, `confirmation_max_support_score=0.75`, `neighbor_cosine_th=0.86`, `neighbor_cell_ratio_min=0.35`, `identity_z_th=-0.55`, `pressure_z_th=0.0`, `min_support_score_for_apply=0.7`, and the MERSCOPE-family specificity/detection values shown above.
- No dataset-specific resolved numeric variation is evidenced within either resolution family. The only apparent extra split is a schema-recording difference: `confirmation_marker_max_support_score` is absent in four older spot summaries and recorded as null in the others; this is not evidence of a different numerical setting.
- Unrecoverable across historical Stage3A outputs: `eps`, `plugin_genes_path`, `gene_weights_path`, exact CLI overrides, exact working directory, git commit, and historical config hash. Matching current YAML hashes are listed only as present-day evidence, not as historical-run proof.
- Protection/guards: recorded `bt_neighbor_guard_*`, population-protection, weak-mismatch, support-margin, and presence-gate fields are null or absent in the resolved outputs. No active T/NK-like, B/T-neighbour, population-protection, keep-weak-mismatch, immune-similarity, or other guard can be claimed from C5 run evidence. Configurability in source code is not treated as use.

## 5. C5 source-value provenance

| file | evidence | relevant fields | rows/items | status |
| --- | --- | --- | --- | --- |
| result/c5_stage3a_decomposition/c5_stage3a_exclusion_by_dataset.csv | 15/15 target detection/exclusion; off-target flag/details; `NK cells:66` | dataset,target,target_total_cells,target_dropped_cells,target_exclusion_fraction,target_detected_missing,off_target_exclusion,off_target_details | 15 | canonical C5 dataset-level result; no separate immutable freeze manifest |
| result/c5_decomposition/c5_stage3b_decomposition_by_dataset.csv | standalone versus sequential Stage3B masks/counts/fractions/Jaccard | standalone_withheld,standalone_fraction,sequential_withheld,sequential_fraction,mask_intersection,mask_union,jaccard,standalone_mask_path,sequential_mask_path | 15 | canonical C5 decomposition result |
| result/c5_decomposition/c5_merscope_four_route_suppression_long.csv | 5 MERSCOPE × 4 routes with fixed low-support denominator | dataset,route,n_total_units,n_low_support_units,n_withheld_total,n_withheld_in_low_support,low_support_suppression_score,assignment_path,stage3b_mask_path,fixed_support_source | 20 | canonical C5 four-route long result |
| result/c5_decomposition/c5_merscope_four_route_suppression_wide.csv | five-dataset route/effect matrix | baseline,stage3a_only,stage3b_only,full,effect_stage3a,effect_stage3b,effect_full,increment_stage3b_after_stage3a,increment_stage3a_after_stage3b | 5 | canonical C5 four-route wide result |
| visualizations/c5_decomposition/c5_fig3_spatial_transition_map_source.csv | representative Breast P1 transition states and coordinates | spot_id,x,y,standalone_withheld,sequential_withheld,transition_state | 2561 | canonical Figure 3 source table |
| result/c5_decomposition/c5_run_runtime.json | per-dataset Stage3B paths, retained/dropped reference counts, reuse status | JSON fields: retained_cells,dropped_cells,dropped_populations,standalone_mask,sequential_mask | 15 | runtime manifest |

The Stage3A exclusion CSV is traceable to the corresponding unsuffixed `stage3_summary.json` and `cell_type_relabel.csv` outputs. The Stage3B runtime manifest and suppression tables trace the sequential route to Stage3A-retained references. Exact historical Stage3A command/commit provenance remains partial as stated above.

## 6. C6 configuration and provenance

### Formal Fibroblast reference-dropout run

| Item | Recovered value | Direct source |
| --- | --- | --- |
| Control config | `configs/datasets/c6_breast_independent_control.yaml`; SHA-256 `9782f85af11294cd1aeb6c7e310600f67e79b2188476908e5cb167718080b2a2` | configs/datasets/c6_breast_independent_control.yaml |
| Dropout config | `configs/datasets/c6_breast_independent_fibroblast_dropout.yaml`; SHA-256 `f7a427b82162116ca8be44707a60aa652ae95764b9f0dc26f1188c4861b489ec` | configs/datasets/c6_breast_independent_fibroblast_dropout.yaml |
| Stage3A normalized config identity | {"c6_breast_independent_control": "09bfe93cc5bf6a03af36c7c5d6fc8901c008e05485e12c333d577d950fe439bc", "c6_breast_independent_fibroblast_dropout": "09bfe93cc5bf6a03af36c7c5d6fc8901c008e05485e12c333d577d950fe439bc"} | result/c6_independent_source_evaluation/c6_provenance.json |
| Stage3B normalized config identity | {"c6_breast_independent_control": "8c03f5ae38944fc3ee9d2f6ffe4682dda23146dc7410251b743da5dd6420e7db", "c6_breast_independent_fibroblast_dropout": "8c03f5ae38944fc3ee9d2f6ffe4682dda23146dc7410251b743da5dd6420e7db"} | result/c6_independent_source_evaluation/c6_provenance.json |
| Seed | 42 | control/dropout sim_info.json; c6_provenance.json |
| Reference source | human_breast_cancer_real; hECA / Tabula Sapiens breast, TSP4 | data/sim/c6_independent_source/*/sim_info.json; c6_reporting_notes.txt |
| Simulated-ST profile source | real_brca; Wu et al. GSE176078 HER2+ breast cancer | data/sim/c6_independent_source/*/sim_info.json; c6_reporting_notes.txt |
| Three-way shared genes | 16,401 | data/sim/c6_independent_source/*/sim_info.json; c6_provenance.json |
| Dropout | Fibroblasts; 2,044 SC cells removed; ST/truth unchanged | dropout sim_info.json |
| Exact successful command lines | NOT DOCUMENTED / NOT RECOVERABLE | execution logs preserve modules/inputs/results, but not complete successful command lines consistently |
| Working directory | E:/AAA文件/Experiment/SVTuner/sctuner2.0 (task/run context; logs contain project paths) | C6 execution context and logs |
| Git commit/code version | NOT DOCUMENTED / NOT RECOVERABLE | no run commit field in C6 provenance/logs |

### C6 scripts and execution evidence

| Stage | Implementation/evidence |
| --- | --- |
| Control simulation | scripts/generate_real_brca_clustered_sim.py; control sim_info.json |
| Reference dropout | scripts/generate_st_only_reference_dropout_from_sim.py; dropout sim_info.json |
| Stage1 | r_scripts/stage1_preprocess.R; result/c6_independent_source_evaluation/execution_logs/*_stage1.log |
| Stage3A | src/stages/stage3_type_plugin.py; execution_logs/*_stage3a.log; formal stage3_summary.json |
| Stage3B | src/stages/stage3b_st_unsupported.py; execution_logs/*_stage3b.log; formal stage3b_summary.json |
| CytoSPACE | src/stages/stage4_cytospace.py; execution_logs/*_stage4_*.log |
| Evaluation | scripts/evaluate_c6_independent_source.py |

### C6 resolved settings

- Stage3A uses an independent C6 dataset config, not a C5 config file. Both control and dropout explicitly contain the spot-resolution parameter family: `strong_th=0.7`, `weak_th=0.4`, `st_cluster_k=30`, `unknown_floor=0.3`, `min_cells_rare_type=20`, `eps=1e-8`, `robust_z_th=-2.0`, `soft_z_th=-0.9`, and the same confirmation/masked/marker settings. Formal Stage3 summaries confirm the resolved values; `c6_provenance.json` records identical normalized Stage3A identities.
- Stage3B resolved config: `{"enable_compensatory_diagnostics": true, "enable_residual_program_branch": true, "expression_scale": "log1p", "fdr": 0.05, "max_genes": 0, "n_calibration": 0, "n_spatial_permutations": 200, "random_seed": 42, "residual_program_candidate_fdr_factor": 2.0, "residual_program_components": 8, "residual_program_min_reference_orthogonal_score": 0.3, "residual_program_min_whole_profile_overlap": 0.5, "sc_expr_source": "normalized", "sc_profile_scale": null, "sc_profile_source": null}`.
- CytoSPACE settings: `{"baseline_filter_mode": "none", "mapping_cells_per_spot": 5, "n_processors": 1, "n_subspots": 800, "sc_expr_source": "normalized", "seed": 42, "svtuner_cell_type_column": "plugin_type", "svtuner_filter_mode": "plugin_unknown", "svtuner_filter_scope": "unsupported_all", "svtuner_stage3b_blank_regions": true}`.

### C6 formal values and direct sources

| Value | Formal value | Direct source |
| --- | --- | --- |
| Precision | 0.8796680497925311 (reported 0.8797) | result/c6_independent_source_evaluation/c6_stage3b_dropout_metrics.csv::precision |
| Recall | 0.4742729306487696 (reported 0.4743) | result/c6_independent_source_evaluation/c6_stage3b_dropout_metrics.csv::recall |
| Truth-supported FPR | 0.014002897151134718 (reported 0.0140) | result/c6_independent_source_evaluation/c6_stage3b_dropout_metrics.csv::supported_spot_FPR |
| AUROC | 0.6663253170176844 (reported 0.6663) | result/c6_independent_source_evaluation/c6_stage3b_dropout_metrics.csv::AUROC |
| AP | 0.5202356989910382 (reported 0.5202) | result/c6_independent_source_evaluation/c6_stage3b_dropout_metrics.csv::AP |
| Recovery | 0.5922201547901481 → 0.6401778198077774; Δ=0.047957665017629325 | result/c6_independent_source_evaluation/c6_summary_metrics.csv; c6_supplementary_table.csv |

## 7. C1/C3/C4/C7/C8 provenance

### C1

| Audit item | Recovered evidence |
| --- | --- |
| Design | 90 rows = 9 scenarios × 10 seeds; seeds 42–51 |
| Seed 42 source | visualizations/method_comparison/composite_no_noise/composition_recovery_7mapping_methods_composite_no_noise_abstention_aware_scenario.csv |
| Seeds 43–51 source | result/c1_seed_scores.csv |
| Analysis/reconstruction script | scripts/plot_c1_multiseed_paired_dotbox.py |
| Scenario/overall summary | visualizations/method_comparison/c1_multiseed/c1_paired_difference_summary.csv |
| Seed-level nine-scenario mean | first average within each simulation_seed across 9 scenarios; then summarize 10 seed-level values |
| CytoSPACE mean | 0.5753428212 (reported 0.5753) |
| SVTuner mean | 0.6728616075 (reported 0.6729) |
| Mean paired delta | 0.0975187864 (reported +0.0975) |
| 95% CI | 0.0673056253–0.1277319475 |
| Positive overall seeds | 10/10 |
| Scenario CIs > 0 | 7/9 |

The script that originally wrote `c1_paired_difference_summary.csv` is not separately persisted; the current plotting script reconstructs the same 90-row frozen dataset.

### C3

| Marker-core threshold | Mean target-associated fraction (%) | Mean core withheld rate (%) |
| --- | --- | --- |
| 10% | 89.72 | 72.60 |
| 15% | 93.79 | 57.78 |
| 20% | 94.79 | 46.44 |
| 25% | 95.67 | 38.24 |

- Direct source: `visualizations/stage3b_realdata_candidate_scan/stage3b_threshold_robustness/stage3b_threshold_robustness_summary.csv`; six-scene source: adjacent `stage3b_threshold_robustness_detail.csv`.
- Analysis script: `scripts/plot_stage3b_threshold_robustness.py`.
- The 15% row is the original threshold and records 93.787889% target-associated fraction and 57.781721% core withheld rate, matching the stated 93.79% and 57.78%.

### C4

- Technical perturbation: `result/c4_stage3b_technical_perturbation/c4_2_technical_perturbation_pilot_summary.csv`; provenance JSON in the same directory; script `scripts/run_c4_2_technical_perturbation_pilot.py`.
- BH→BY and 200→1000 permutation results: `result/c4_stage3b_statistical_robustness/` per-dataset and combined CSV/JSON outputs; script `scripts/run_c4_3_stage3b_statistical_robustness_pilot.py`.
- Compact final values: `result/c4_stage3b_calibration/c4_supplementary_table.csv`. It records final-mask Jaccard=1.0 for all three datasets and seed=42.
- Pseudo-null raw distributions exist: one `c4_1_raw_null_statistics.csv` per dataset under `result/c4_stage3b_calibration/<sample>/`, with `pseudo_null_id,relative_reconstruction_error,cosine_deficit,positive_residual_fraction,residual_concentration`.
- Observed raw diagnostic distributions are not copied into a dedicated C4 QQ/ECDF source-value file. They remain in each existing Stage3B `spot_unsupported_scores.csv`, using the same four statistic columns. The summary file is `result/c4_stage3b_calibration/c4_1_observed_vs_null_summary.csv`.
- QQ/ECDF plotting/comparison script: `scripts/run_c4_1_observed_vs_null_comparison.py`; raw-null export script: `scripts/export_c4_1_stage3b_raw_calibration_null.py`.

### C7

| Item | Value/source |
| --- | --- |
| Frozen endpoint | visualizations/bioapp_experiment/bioapp_phase2c_endpoint_freeze_with_validated_cta_to_spot_mapping/spot_level_endpoint_freeze.csv::primary_endpoint_status |
| Frozen SVTuner output | visualizations/bioapp_experiment/bioapp_phase7_svtuner_aware_execution_against_frozen_cta_immune_endpoint/svtuner_immune_all_dropout/svtuner_immune_all_dropout_spot_level_raw_output_contract.csv::withheld_score,withheld_binary |
| Analysis set | n=1888; positive=134; negative=1754 |
| Binary metrics | TP=60, FP=38, TN=1716, FN=74, precision=0.6122448979591837, recall=0.44776119402985076, F1=0.5172413793103449, specificity=0.9783352337514253; source c7_metrics_summary.csv |
| AP | 0.2836551170002947; source c7_metrics_summary.csv |
| AUROC | 0.8305; persisted only in final_supplementary/C7_Supp_metrics_table.csv and hard-coded in scripts/plot_c7_final_supplementary.py; the metric-computation record is NOT DOCUMENTED / NOT RECOVERABLE |
| P≥0.25 | threshold=0.5016863844421177, precision=0.2500, recall=0.6044776119402985; recoverable from c7_pr_curve_points.csv |
| P≥0.30 | selected frozen point threshold=0.5280528671035196, precision=0.3008474576271186, recall=0.5298507462686567; curve source c7_pr_curve_points.csv; plotting script rounds label to 0.3008/0.5299 |
| P≥0.50 | not achieved by any finite threshold; c7_fixed_precision_summary.csv and final supplementary table |
| Closest simple threshold differs by 93 units | Value was reported in the interactive C7-1.5 step but is not persisted in a C7 result file or computation script: NOT DOCUMENTED / NOT RECOVERABLE as file provenance |
| Figure/table script | scripts/plot_c7_final_supplementary.py |

### C8

| Metric ID | Metric | n | mean paired Δ | bootstrap 95% CI | Cohen dz | favorable | two-sided Wilcoxon P | bootstrap | seed |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A1_reciprocal_suppression | Reciprocal suppression | 10 | 0.010648844982 | [0.00481044352448, 0.0169992540615] | 1.01602231882 | 9/10 | 0.009765625 | 10000 | 20260927 |
| A2_cosine_similarity | Reconstructed-expression cosine similarity | 10 | 0.00644454851649 | [0.00204277983152, 0.0119332636671] | 0.750090763526 | 8/10 | 0.02734375 | 10000 | 20260927 |
| A3_ecotyper_experiment_mean | EcoTyper normalized enrichment (experiment mean) | 6 | 0.186323161234 | [0.0791912714464, 0.297021782765] | 1.22871094203 | 5/6 | 0.0625 | 10000 | 20260927 |
| B1_merscope_suppression | MERSCOPE low-support suppression score | 5 | 0.0405351454458 | [0.0129583652898, 0.0838427163825] | 0.851179641975 | 5/5 | 0.0625 | 10000 | 20260927 |
| B2_merscope_peak_es | MERSCOPE peak ES / low-support enrichment | 5 | 0.00380841949706 | [-0.0877467428721, 0.108795387243] | 0.0287031628505 | 2/5 | 1 | 10000 | 20260927 |

- Formal statistics source: `result/c8_direct_paired_comparison/c8_paired_statistics.csv`.
- Formal paired source values: `result/c8_direct_paired_comparison/c8_canonical_pairs.csv` and `visualizations/c8_direct_paired_comparison/C8_direct_paired_comparison_source_values.csv`.
- Publication table: `visualizations/c8_direct_paired_comparison/C8_direct_paired_comparison_stats_table.csv`.
- Analysis script: `scripts/run_c8_direct_paired_comparison.py`; figure scripts: `scripts/plot_c8_direct_paired_comparison.py` and `_v2.py`.
- The statistics file explicitly records two-sided Wilcoxon tests, Cohen dz, 10,000 paired-bootstrap iterations, and seed 20260927 for all five analyses.

## 8. Missing or unrecoverable information

- C5: exact historical Stage3A command/launcher for all 15 experiments.
- C5: exact historical working directory and git commit/code checkout.
- C5: historical config binding/hash. Fourteen current YAMLs match resolved values but are not direct proof of the config file used at run time; mouse embryo YAML is absent.
- C5: resolved `eps`, `plugin_genes_path`, and `gene_weights_path` are not stored in `stage3_summary.json`.
- C5: no active protection/guard can be established; relevant resolved fields are null or absent.
- C6: complete successful command lines and git commit are not consistently persisted, although configs, logs, formal outputs, and provenance are present.
- C7: the AUROC computation artifact/script is not persisted; 0.8305 survives only in the final table and plotting-script constant.
- C7: the 93-discordant-unit score-threshold check is not persisted in a result file or script.
- C1: the original script that wrote `c1_paired_difference_summary.csv` is not separately persisted; the current plot script reconstructs the same frozen 90-row input.

## 9. Final authoritative-file list

- `C5_stage3a_resolved_parameters.tsv`
- `result/c5_stage3a_decomposition/c5_stage3a_exclusion_by_dataset.csv`
- `result/c5_decomposition/c5_stage3b_decomposition_by_dataset.csv`
- `result/c5_decomposition/c5_merscope_four_route_suppression_long.csv`
- `result/c5_decomposition/c5_merscope_four_route_suppression_wide.csv`
- `result/c5_decomposition/c5_run_runtime.json`
- `visualizations/c5_decomposition/c5_fig3_spatial_transition_map_source.csv`
- `configs/datasets/c6_breast_independent_control.yaml`
- `configs/datasets/c6_breast_independent_fibroblast_dropout.yaml`
- `data/sim/c6_independent_source/c6_breast_independent_control/sim_info.json`
- `data/sim/c6_independent_source/c6_breast_independent_fibroblast_dropout/sim_info.json`
- `result/c6_independent_source_evaluation/c6_provenance.json`
- `result/c6_independent_source_evaluation/c6_summary_metrics.csv`
- `result/c6_independent_source_evaluation/c6_stage3a_audit.csv`
- `result/c6_independent_source_evaluation/c6_stage3b_dropout_metrics.csv`
- `result/c6_independent_source_evaluation/c6_supplementary_table.csv`
- `result/c1_seed_scores.csv`
- `visualizations/method_comparison/composite_no_noise/composition_recovery_7mapping_methods_composite_no_noise_abstention_aware_scenario.csv`
- `visualizations/method_comparison/c1_multiseed/c1_paired_difference_summary.csv`
- `visualizations/stage3b_realdata_candidate_scan/stage3b_threshold_robustness/stage3b_threshold_robustness_detail.csv`
- `visualizations/stage3b_realdata_candidate_scan/stage3b_threshold_robustness/stage3b_threshold_robustness_summary.csv`
- `result/c4_stage3b_calibration/c4_supplementary_table.csv`
- `result/c4_stage3b_calibration/c4_1_observed_vs_null_summary.csv`
- `result/c4_stage3b_statistical_robustness/`
- `result/c4_stage3b_technical_perturbation/`
- `visualizations/bioapp_experiment/C7_cta_extended_performance/c7_metrics_summary.csv`
- `visualizations/bioapp_experiment/C7_cta_extended_performance/c7_pr_curve_points.csv`
- `visualizations/bioapp_experiment/C7_cta_extended_performance/c7_fixed_precision_summary.csv`
- `visualizations/bioapp_experiment/C7_cta_extended_performance/final_supplementary/C7_Supp_metrics_table.csv`
- `result/c8_direct_paired_comparison/c8_canonical_pairs.csv`
- `result/c8_direct_paired_comparison/c8_paired_statistics.csv`
- `visualizations/c8_direct_paired_comparison/C8_direct_paired_comparison_source_values.csv`
- `visualizations/c8_direct_paired_comparison/C8_direct_paired_comparison_stats_table.csv`
