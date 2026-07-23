# Task 4E: formal Fig. 2 enrichment-backend audit

Audit date: 2026-07-22

## Scope

This is a static, read-only provenance audit of the retained Fig. 2b, 2c, 2d and 2f enrichment assets. It did not access the formal manuscript or bibliography, rerun any experimental stage, rerun CytoSPACE, or alter a formal figure/source-value table.

## Decision

**CONDITIONAL PASS.** Formal source values and figure inputs are located and all requested headline values reproduce exactly. Fig. 2b, 2d and 2f have direct custom-Python chains. Fig. 2c is also classified as `CUSTOM_PYTHON` with high confidence: the formal runner defaults to that branch with 1,000 permutations, and every retained P value lies on the exact `1/1001` lattice. Both retained R candidates instead fix `nperm=10000`. The missing historical command/log and pair-level intermediate summaries prevent an unconditional PASS.

No audit-only reconstruction is required.

## Panel conclusions

| Panel | Backend | Evidence | Status |
|---|---|---|---|
| Fig. 2b | `CUSTOM_PYTHON` independent running-enrichment statistic | Direct formal plot source and direct input table | Confirmed |
| Fig. 2c | `CUSTOM_PYTHON` 1,000-permutation enrichment | Runner default plus formal P-value lattice and retained aggregate source | High-confidence conditional |
| Fig. 2d | `CUSTOM_PYTHON` (`python_permutation_gsea`) | Backend and `nperm=1000` recorded in formal/upstream rows | Confirmed |
| Fig. 2f | `CUSTOM_PYTHON` 1,000-permutation enrichment | Direct metrics generator and retained 120-row table | Confirmed |

## Key source checks

- Fig. 2b Peak ES: CytoSPACE `0.064410337138`, SVTuner `0.074773711137`.
- Fig. 2c mean NES: CytoSPACE `1.3458403003` -> `1.3458`; SVTuner `1.5321634615` -> `1.5322`; wins `11/12`.
- Fig. 2d NES: `2.171340694744` and `2.523755126438`; both P values round to `0.0010`.
- Fig. 2f: `50` paired valid readouts; SVTuner higher in `28/50`.

## Important boundary

The existence of `compute_cytospace_fig2c_official_enrichment_seurat.R` does not establish formal use. Its fixed 10,000-permutation fgsea path is inconsistent with the complete 1,000-permutation P-value lattice in the retained Fig. 2c source table. No Seurat enrichment-backend claim should be added to the manuscript.

See `backend_decision_by_panel.tsv`, `formal_source_value_consistency.tsv`, `audit_summary.json`, and `final_decision.md` for the formal audit record.
