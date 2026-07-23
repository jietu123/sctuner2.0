# Task 4F manifest reconciliation and v1.0 freeze candidate

Generated: 2026-07-22

## Decision

**CONDITIONAL_PASS_TO_FREEZE_REVIEW**

The v0.2 source manifests remain unchanged at 37 dataset records and 40 experiment records. Their status anchors also remain unchanged:

- datasets: {'CONFIRMED': 28, 'PARTIAL': 7, 'PROJECT_GENERATED': 2}
- experiments: {'CONFIRMED': 29, 'PARTIAL': 9, 'PROJECT_GENERATED': 2}

The v1.0 files are candidates only and are marked `PENDING_HUMAN_APPROVAL`. They are not `FROZEN`, `FINAL` or `RELEASED`.

## Reconciliation

The Task 4F user-supplied manuscript fact ledger was converted into one check per fact without reading the formal manuscript or bibliography. Reconciliation counts are:

{'MATCH': 174, 'MATCH_WITH_QUALIFIED_WORDING': 17}

No blocking manuscript-manifest conflict was found. Six explicitly bounded non-blocking fields remain unresolved: the human-lung run accession, CID4465 slide, CID4465 capture area, CTA exact chemistry, Vizgen provider checksum and Vizgen redistribution permission.

## Classification

All 37 dataset records have an allowed source-layer classification. Dataset and experiment redistribution classes are conservative; combined counts are:

{'REDISTRIBUTABLE_DERIVATIVE_ONLY': 53, 'LINK_AND_RECONSTRUCT_ONLY': 18, 'RESTRICTED': 1, 'PERMISSION_UNRESOLVED': 5}

No provider file with unresolved permission is classified as `REDISTRIBUTABLE`. Project-generated derivatives are not represented as independent public datasets.

## Hash validation

All retained input and formal configuration paths referenced by v0.2 exist and have SHA-1 values. Current bytes were checked against those hashes; large-file hashes are cached only when path, size and nanosecond mtime remain unchanged.

## Boundary

This task did not access or modify the formal manuscript, bibliography, figures, source-value tables, experimental code, configurations or experimental outputs. It did not run Stage0--Stage5, CytoSPACE or any figure generator, and it did not access a GitHub remote. The existing Task 4E release manifest was protected and not edited.
