# BioApp Phase 3 Reference Audit

Purpose: audit whether the frozen Phase 2C CTA-defined Immune cells endpoint has enough reference and input support for later formal baseline-vs-endpoint comparison.

Phase 2C input: `visualizations/bioapp_phase2c_endpoint_freeze_with_validated_cta_to_spot_mapping`.

Reference audit method: scan candidate files, read recommended scRNA metadata, audit label vocabulary, and match immune-related labels using predefined keywords.

Immune label keywords: immune, T cell, T_cells, T-cell, CD4, CD8, Treg, NK, B cell, B_cells, B-cell, Plasma, Macrophage, Monocyte, Myeloid, Dendritic, DC, Mast, Neutrophil, Lymphocyte, Leukocyte.

Dropout design: primary design removes all immune-related labels; size-matched non-immune control is only designed, not executed.

Formal baseline feasibility: checks ST expression availability, scRNA expression, metadata, gene overlap, barcode compatibility, and CytoSPACE entry availability without running CytoSPACE.

Endpoint-specific metrics are only planned in `endpoint_specific_metric_plan.json`; no metric values are computed in Phase 3.

Decision: REVIEW_REQUIRED
Reason: formal baseline inputs are incomplete or need manual preparation before Phase 4

Allowed claims:
- Reference label vocabulary was audited against the frozen CTA Immune cells endpoint.
- A formal baseline feasibility plan was established without running CytoSPACE or SVTuner.
- Endpoint-specific metrics were predefined for later baseline/SVTuner comparison.

Disallowed claims:
- Baseline creates false niche calls.
- SVTuner improves endpoint recovery.
- SVTuner prevents contradicted interpretation.
- Biological application completed.
- Biological discovery made.
