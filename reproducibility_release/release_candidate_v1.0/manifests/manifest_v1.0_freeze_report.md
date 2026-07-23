# Task 4F v1.0 Manifest Freeze Report

## Decision

`PASS_MANIFESTS_FROZEN`

Human approval was recorded on 2026-07-22. Candidate hashes, record counts,
referential integrity, redistribution guardrails, and the Task 4E release
inventory baseline were verified before the formal v1.0 artifacts were created.

## Frozen manifests

| Artifact | Records | SHA-1 | Candidate byte identity |
|---|---:|---|---|
| `reproducibility_release/manifests/dataset_manifest_v1.0.tsv` | 37 | `64223e7846fc3cc8228a076b0f3e069742f58afe` | true |
| `reproducibility_release/manifests/experiment_manifest_v1.0.tsv` | 40 | `347539ed5677dc89318772333552d5944821a355` | true |

Formal metadata SHA-1: `0141fc28d9a106a8c7995a5a5658bf60fe298ef5`

Human approval record SHA-1: `711ae2b8a6d973725ab28e0423d5925fe90ca495`

## Accepted unresolved non-blocking fields

- human-lung read-run accession
- CID4465 slide identifier
- CID4465 capture area
- CTA exact chemistry or kit version
- Vizgen provider-supplied checksum
- Vizgen explicit redistribution permission

These fields remain unresolved. Approval does not authorize inferred values or
redistribution of third-party provider files whose permission is unresolved.

## Release constraints

- Vizgen original provider files are not included in the public release while explicit redistribution permission remains unresolved.
- CID4465 slide and capture-area fields remain unresolved.
- CTA exact chemistry or kit version remains unresolved.
- The human-lung ENA sample accession ERS20065156 must not be relabelled as a sequencing-run accession.

## Immutability rule

`dataset_manifest_v1.0.tsv` and `experiment_manifest_v1.0.tsv` are immutable.
Any later substantive field change requires v1.0.1, v1.1, or a later version.
Silent modification of v1.0 is prohibited.

## Guardrails

- Existing Task 4E release records are preserved byte-for-byte.
- The release inventory receives append-only Task 4F records.
- No provider-hosted raw data are included.
- No experiment stage, CytoSPACE run, figure generation, or data download was performed.
- No formal manuscript or bibliography was accessed or modified.
