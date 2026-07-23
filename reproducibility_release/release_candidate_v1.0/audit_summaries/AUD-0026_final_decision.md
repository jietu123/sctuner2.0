# Final decision

**Task:** Task 4E formal Fig. 2 enrichment-backend provenance audit  
**Date:** 2026-07-22  
**Decision:** **CONDITIONAL PASS**

## Findings

- **Fig. 2b:** `CUSTOM_PYTHON`. This is an independent unweighted running-enrichment/Peak ES display calculation and has no NES, permutation, Seurat, or fgsea step.
- **Fig. 2c:** `CUSTOM_PYTHON`, high-confidence conditional. The runner defaults to the Python branch at 1,000 permutations and all 48 retained P values are exact multiples of `1/1001`. Both retained R alternatives use 10,000 permutations. The historical command/log and pair-level summaries are absent.
- **Fig. 2d:** `CUSTOM_PYTHON`, confirmed by `backend=python_permutation_gsea` and `nperm=1000` in both formal rows and both upstream summaries.
- **Fig. 2f:** `CUSTOM_PYTHON`, confirmed by the direct metrics generator; the retained table reproduces 50 valid paired readouts and 28/50 SVTuner wins.

## Consistency

All requested formal source-value checks pass: Fig. 2c `1.3458`, `1.5322`, and `11/12`; Fig. 2d `2.17`, `2.52`, and both `P=0.0010`; Fig. 2f `50` and `28/50`.

## Manuscript implication

Do not add a Seurat or fgsea enrichment-backend claim. A manual Methods update may describe the project-specific Python permutation implementation. Exact formal runtime package versions should remain unspecified unless a historical run-local environment record is recovered.

## Reconstruction

Not required. The remaining gaps concern historical command/runtime documentation, not numerical identity or algorithm definition.

## Guardrails

- Formal manuscript accessed: false
- Formal manuscript modified: false
- Formal bibliography accessed: false
- Formal bibliography modified: false
- Formal figures modified: false
- Formal source-value tables modified: false
- Experimental code modified: false
- Experimental outputs modified: false
- Stage0--Stage5 rerun: false
- CytoSPACE rerun: false
- GitHub access required: false
