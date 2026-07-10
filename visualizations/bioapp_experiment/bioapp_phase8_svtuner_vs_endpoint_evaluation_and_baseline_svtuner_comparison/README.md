# BioApp Phase 8

Stage: SVTuner-vs-endpoint evaluation and baseline/SVTuner comparison.

This phase reads frozen outputs from Phase 2C, Phase 5, Phase 6, and Phase 7. It does not rerun SVTuner, Stage3, Stage3B, Stage4, CytoSPACE, or any Phase 4 baseline.

Frozen endpoint:
- CTA-defined Immune cells.
- The endpoint is used only for evaluation.
- The endpoint is not redefined in Phase 8.

Baseline anchor:
- baseline_immune_all_dropout

SVTuner condition:
- svtuner_immune_all_dropout

Gene count audit:
- final_effective_gene_count_for_phase8 = 2000
- material_mismatch = False
- explanation = formal_gene_intersection.txt contains 2000 non-empty unique gene lines; Stage3B also reports 2000 genes. The previous 1999 count is attributable to a line-count/read-mode artifact, not a material mismatch.

Scores:
- withheld_score and withheld_binary are imported from Phase 7 raw SVTuner contract.
- binary withheld status uses the predefined Phase 7 Stage3B rule; Phase 8 does not tune thresholds.

Endpoint analysis set:
- Main metrics use CTA Immune-positive and CTA Immune-negative spots.
- Ambiguous and excluded spots remain in spot-level tables but are excluded from main metrics.

Operational definitions:
- Burden reduction is a withheld-aware forced-call burden estimate.
- It is not a post-Stage4 remapping result.

Allowed claims:
- SVTuner raw unsupported / withheld outputs were evaluated against the frozen CTA Immune endpoint.
- Baseline and SVTuner-aware outputs were compared under the immune-all-dropout condition.
- Endpoint-specific comparison metrics were generated for final biological-application audit.

Disallowed claims:
- Biological application completed.
- Biological discovery made.
- SVTuner definitively improves biological interpretation.

Decision:
- PASS

Next:
- BioApp Phase 9 final biological-application audit, figure selection, and interpretation boundary lock.
