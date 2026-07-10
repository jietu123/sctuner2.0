# BioApp Phase 7

This phase performs SVTuner-aware raw execution for the pre-authorized
immune-all-dropout condition only.

Manual authorization was granted only for raw SVTuner-aware execution of the
immune-all-dropout condition. It does not authorize final prevention analysis
or biological application claims.

Inputs come from Phase 2C, Phase 3R, Phase 4, Phase 5, and Phase 6. The frozen
CTA Immune endpoint is retained only for spot-universe alignment and the raw
contract status column; it is not used for mapping, threshold selection, or
Phase 7 evaluation.

Primary SVTuner condition: `svtuner_immune_all_dropout`.
Baseline anchor: `baseline_immune_all_dropout`.

This phase does not compute prevention, SVTuner-vs-endpoint metrics, or
baseline/SVTuner endpoint-specific comparisons.
