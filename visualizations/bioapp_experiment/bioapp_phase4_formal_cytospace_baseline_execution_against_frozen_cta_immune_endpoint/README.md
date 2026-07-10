# BioApp Phase 4

Purpose: execute formal CytoSPACE baseline runs against the frozen CTA-defined Immune cells endpoint input universe.

Inputs come from Phase 2C and Phase 3R. The endpoint is only an anchor for later evaluation and was not used for mapping.

ST expression uses the full 2248 frozen endpoint spot universe. The CytoSPACE ST data matrix is `log1p_CPM_from_counts_in_Phase3R`; it was not used to define or select the endpoint.

Exact-overlap genes used: 2000.

Baseline runs:
- baseline_full_reference
- baseline_immune_all_dropout
- baseline_nonimmune_all_dropout_control

Control limitation: true non-immune size-matched control is impossible because non-immune pool size is 1600 and immune dropout size is 2414. The primary available label-level clean control is nonimmune-all-dropout.

This phase does not compute endpoint metrics and does not run SVTuner, Stage3, or Stage4.

Decision: `PASS`
