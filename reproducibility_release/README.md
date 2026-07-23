# SVTuner reproducibility release

This directory is the public entry point for the manuscript-associated
SVTuner `v1.0.0-rc1` release candidate. It packages the frozen provenance,
source-value, code, configuration, environment, reconstruction, and audit
records needed to inspect the reported analyses without redistributing the
complete third-party inputs or full mapping intermediates.

## Release status

- Task 4H decision: `PASS_WITH_NONBLOCKING_GAPS`
- Candidate files: `497`
- Hashed payload records: `482`
- Candidate size: `46,364,491` bytes
- Hash validation: `PASS`
- Fig. 1 through Fig. 5 source-value dependency groups: `19`, complete
- Accepted nonblocking gaps: `13`
- Blocking gaps: `0`

The frozen candidate is in [`release_candidate_v1.0/`](release_candidate_v1.0/).
Its existing files were not regenerated or modified during release
preparation.

## Directory index

| Entry | Purpose |
|---|---|
| [`manifests/`](manifests/) | Frozen dataset and experiment manifests and freeze records |
| [`source_values/`](source_values/) | Index to the 161 hashed source-value records |
| [`audit_records/`](audit_records/) | Index to the 51 selected public audit summaries |
| [`environments/`](environments/) | Environment specifications and runtime evidence |
| [`configurations/`](configurations/) | Dataset and pipeline configuration inventory |
| [`inventories/`](inventories/) | Task 4H asset, hash, source-value, and script inventories |
| [`redistribution/`](redistribution/) | Redistribution decisions and release-package checksum |
| [`release_candidate_v1.0/reconstruction/`](release_candidate_v1.0/reconstruction/) | Data acquisition and reconstruction instructions |

## Frozen manifests

| File | Records | SHA-1 |
|---|---:|---|
| [`dataset_manifest_v1.0.tsv`](manifests/dataset_manifest_v1.0.tsv) | 37 | `64223e7846fc3cc8228a076b0f3e069742f58afe` |
| [`experiment_manifest_v1.0.tsv`](manifests/experiment_manifest_v1.0.tsv) | 40 | `347539ed5677dc89318772333552d5944821a355` |
| [`manifest_v1.0_metadata.json`](manifests/manifest_v1.0_metadata.json) | n/a | `0141fc28d9a106a8c7995a5a5658bf60fe298ef5` |

The manuscript-to-manifest reconciliation contains `174 MATCH` and `17
MATCH_WITH_QUALIFIED_WORDING` records, with no missing records and no
conflicts.

## Public data boundary

The 482 candidate payload records comprise:

- `264` `REDISTRIBUTABLE` records;
- `161` `PROJECT_GENERATED_REDISTRIBUTABLE` records;
- `57` `METADATA_ONLY` records.

Complete raw or retained third-party data are not redistributed. This includes
raw Vizgen MERSCOPE provider files, complete GEO and 10x Genomics copies,
complete Space Ranger outputs, large scRNA-seq expression matrices, complete
CytoSPACE intermediate outputs, local caches, and failed experiment
directories.

The redistribution matrix classifies `18` datasets or resources as
`LINK_AND_RECONSTRUCT`, `5` Vizgen MERSCOPE datasets as
`UNRESOLVED_PERMISSION`, and `1` dataset as `INTERNAL_NOT_FOR_RELEASE`. None
of those 24 records is included as a data file. Their accession, provider,
sample, checksum, acquisition, and reconstruction metadata are retained where
available.

Use the frozen dataset manifest and
[`reconstruction/README.md`](release_candidate_v1.0/reconstruction/README.md)
to reacquire third-party inputs from official sources. Validate downloaded
files against retained checksums where a provider checksum is available.

## Integrity

The Task 4H hash report contains 482 passing payload checks:

```powershell
Import-Csv `
  reproducibility_release/release_candidate_v1.0/inventory/task4h_hash_validation_report.tsv `
  -Delimiter "`t" |
  Group-Object status
```

The release asset was created as a distribution copy of the existing 497-file
candidate directory. It contains 497 file entries and has SHA-256:

```text
03ecadb41f002fd96b2949a0fe5ad80e2204cd78137bbcafed2d3ed3f86dde80
```

See [`redistribution/`](redistribution/) for the checksum record.

## Recorded local paths

Some frozen provenance records and scripts retain `E:\...` paths from the
original execution environment. They are preserved because changing them
would invalidate frozen hashes. The public candidate contains no personal
user-profile paths, username strings, credentials, tokens, private keys, or
cloud access secrets.

## Known gaps

The 13 accepted nonblocking gaps are recorded in
[`task4h_release_candidate_summary.json`](release_candidate_v1.0/metadata/task4h_release_candidate_summary.json).
They concern bounded accession, sample-detail, provider-checksum,
redistribution-permission, and historical runtime-snapshot gaps. There are no
blocking gaps.

This is a release candidate. Final `v1.0.0` publication remains subject to
figure harmonization and the final manuscript-figure-source-value audit.
