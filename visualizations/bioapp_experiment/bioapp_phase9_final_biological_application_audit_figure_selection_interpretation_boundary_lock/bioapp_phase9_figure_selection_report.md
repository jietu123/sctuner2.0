# BioApp Phase 9 Figure Selection Report

## Decision

BioApp Phase 9 decision: **PASS**

Biological application allowed: **true**

## Main Evidence

The Phase 8 evidence chain supports a formal biological-application evidence claim under the locked boundary:

- Endpoint: CTA-defined Immune cells
- Baseline anchor: baseline_immune_all_dropout
- SVTuner condition: svtuner_immune_all_dropout
- Main analysis spots: 1888
- Effective gene count: 2000

Key Phase 8 metrics:

- withheld_AUROC = 0.8305
- withheld_AUPRC = 0.2837
- withheld_enrichment_delta = 0.1453
- positive withheld rate = 0.4478
- negative withheld rate = 0.0217
- operational contradiction prevention rate = 0.4478

## Figure Triage

| figure_file                                                                                                                                            | classification                         | recommended_use                         | layout_review   |
|:-------------------------------------------------------------------------------------------------------------------------------------------------------|:---------------------------------------|:----------------------------------------|:----------------|
| visualizations/bioapp_phase8_svtuner_vs_endpoint_evaluation_and_baseline_svtuner_comparison/endpoint_baseline_svtuner_three_way_spatial_comparison.pdf | main figure candidate                  | Primary spatial evidence module         | required        |
| visualizations/bioapp_phase8_svtuner_vs_endpoint_evaluation_and_baseline_svtuner_comparison/endpoint_positive_forced_nonimmune_burden_reduction.pdf    | main or supplementary figure candidate | Quantitative burden-reduction support   | minor           |
| visualizations/bioapp_phase8_svtuner_vs_endpoint_evaluation_and_baseline_svtuner_comparison/phase8_svtuner_endpoint_metric_summary.pdf                 | supplementary figure candidate         | Metric summary support                  | minor           |
| visualizations/bioapp_phase8_svtuner_vs_endpoint_evaluation_and_baseline_svtuner_comparison/withheld_score_positive_vs_negative_boxplot.pdf            | supplementary figure candidate         | Distribution-level evidence             | minor           |
| visualizations/bioapp_phase8_svtuner_vs_endpoint_evaluation_and_baseline_svtuner_comparison/withheld_score_roc_curve.pdf                               | supplementary figure candidate         | Metric evidence                         | minor           |
| visualizations/bioapp_phase8_svtuner_vs_endpoint_evaluation_and_baseline_svtuner_comparison/withheld_score_pr_curve.pdf                                | supplementary figure candidate         | Metric evidence under endpoint sparsity | minor           |
| visualizations/bioapp_phase8_svtuner_vs_endpoint_evaluation_and_baseline_svtuner_comparison/binary_withheld_endpoint_contingency_plot.pdf              | QC/audit-only figure                   | Internal or supplementary audit         | not required    |

## Recommended Main Figure Candidate

`endpoint_baseline_svtuner_three_way_spatial_comparison.pdf`

This is the strongest main figure candidate because it connects the frozen endpoint, baseline forced nonimmune burden, SVTuner withheld score, binary withheld output, and endpoint/withheld spatial overlay.

Layout review is required before manuscript use. Only layout-level changes are allowed.

## Interpretation Boundary

Allowed claim:

> SVTuner provides biological-application evidence that its unsupported/withheld output is concordant with an independent CTA-defined Immune endpoint under an immune-all-dropout condition, reducing forced non-immune assignment burden relative to the baseline mapping.

Disallowed claims remain locked in `bioapp_phase9_interpretation_boundary_lock.json`.
