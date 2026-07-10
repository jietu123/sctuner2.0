# BioApp Phase 2C Endpoint Freeze

Purpose: freeze CTA-derived external endpoint status at Visium spot level using the validated Phase 2B-R2 CTA-to-spot mapping.

Input source: `visualizations/bioapp_phase2br2_explicit_rscript_cta_align_locale_safe_dry_run`.

Endpoint freeze rule:
- positive: CTA_mapped is true, total_CTA_objects >= 2, class_count >= 2, and class_fraction >= 0.60
- negative: CTA_mapped is true, total_CTA_objects >= 2, and class_fraction <= 0.20
- ambiguous: CTA_mapped is true and neither positive nor negative
- excluded: CTA_mapped is false

Decision: PASS

Primary endpoint: Immune cells

Allowed interpretation:
- CTA-derived endpoint was frozen at Visium spot level.
- CTA-to-spot mapping enabled endpoint-positive, endpoint-negative, ambiguous, and excluded spot definitions.

Disallowed interpretation:
- SVTuner improves endpoint recovery.
- SVTuner prevents contradicted interpretation.
- Baseline creates false niche calls.
- Biological application completed.
- Biological discovery made.

Next condition:
- PASS: BioApp Phase 3 - reference audit and formal baseline feasibility for CTA-defined endpoint.
- REVIEW_REQUIRED: review endpoint thresholds, CTA spot composition, and spatial continuity before Phase 3.
- FAIL: fix input mapping or boundary violation before retrying Phase 2C.
