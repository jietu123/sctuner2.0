# Fig. 5b precision-recall formal audit

This directory is an independent, read-only audit of the precision-recall metric shown in Fig. 5b. It does not modify or rerun Stage3, Stage3B, Stage4, CytoSPACE, the CTA endpoint, the shared-gene set, the manuscript, the formal figure, or existing source-value files.

## Reproduction

Run from the repository root:

```powershell
python scripts/audit_fig5b_auprc_provenance.py --formal
```

Formal inputs:

- endpoint: `visualizations/bioapp_experiment/bioapp_phase2c_endpoint_freeze_with_validated_cta_to_spot_mapping/spot_level_endpoint_freeze.csv`
- score: `visualizations/bioapp_experiment/bioapp_phase7_svtuner_aware_execution_against_frozen_cta_immune_endpoint/svtuner_immune_all_dropout/svtuner_immune_all_dropout_spot_level_raw_output_contract.csv`
- endpoint column: `primary_endpoint_status`
- score column: `withheld_score`
- binary column: `withheld_binary`

The formal analysis is fixed to CTA endpoint-positive versus endpoint-negative spots. Ambiguous and excluded spots are retained in `excluded_endpoint_360_spots.csv` and are not used for formal discrimination metrics.

`average precision` and `trapezoidal PR-AUC` are reported separately because they are not interchangeable.
