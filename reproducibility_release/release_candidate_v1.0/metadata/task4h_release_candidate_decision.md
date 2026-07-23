# Task 4H release-candidate decision

## Decision

`PASS_WITH_NONBLOCKING_GAPS`

The candidate contains 482 hashed payload assets. All included paths,
SHA-1 values and byte sizes passed validation.
Formal source-value coverage is `COMPLETE` across the required Fig. 1-5
groups. No prohibited third-party raw dataset is included.

## Non-blocking gaps

The six human-approved provenance gaps remain unresolved. Historical run-local
environment snapshots are also unavailable for mainline Stage1, legacy Fig. 2
enrichment, Tangram, novoSpaRc, SpaOTsc and the CellTrek-labelled Python fallback.
These versions were not inferred from the current machine.

## Git release readiness

The repository is on `main` at `47c644715976853717bb3d3093b407a7f54e8006`. The working tree was not clean
before Task 4H, so this commit is not suitable as the manuscript release commit.
An intentional scoped commit and review are required before tag `v1.0.0-rc1`.

## Guardrails

No SVTuner/CytoSPACE stage, metric, mapping or figure was rerun. No manuscript,
BibTeX, figure, endpoint, threshold, result number, frozen manifest or existing
Task 4E/4F inventory record was modified.
