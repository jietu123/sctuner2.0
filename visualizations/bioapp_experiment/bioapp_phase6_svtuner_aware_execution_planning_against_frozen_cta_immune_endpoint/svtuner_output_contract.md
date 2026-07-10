# SVTuner Output Contract

Phase 7 must produce a spot-level table with these fields:

- `barcode`
- `primary_endpoint_status`
- `baseline_anchor_run`
- `svtuner_condition`
- `reference_supported_status`
- `reference_unrepresented_score`
- `withheld_score`
- `withheld_binary`
- `supported_assignment_score`
- `svtuner_immune_score`
- `svtuner_nonimmune_score`
- `svtuner_dominant_label`
- `svtuner_dominant_fraction`
- `baseline_immune_score_from_phase5`
- `baseline_nonimmune_score_from_phase5`
- `baseline_dominant_label_from_phase5`
- `baseline_dominant_fraction_from_phase5`
- `score_parse_status`
- `svtuner_output_available`

Baseline fields must be imported from Phase 5 without recomputation or redefinition.
