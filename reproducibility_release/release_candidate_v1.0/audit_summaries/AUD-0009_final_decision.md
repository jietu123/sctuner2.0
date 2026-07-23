# Task 4B-2 Final Decision

## 1. Decision

`CONDITIONAL PASS`

The major evidence search was completed without guessing. CTA slide/capture, exact retained hECA lineage, five Vizgen sample identities, and CRC slide metadata were recovered. CTA chemistry and Vizgen checksum/redistribution evidence remain partial.

## 2. CTA result

- Sample identity: `BCSA2TumB1` (**CONFIRMED**)
- Platform: 10x Visium Spatial Gene Expression (**CONFIRMED**)
- Chemistry/version: exact version **UNRESOLVED**
- Slide: `V10F24-112` (**CONFIRMED**)
- Capture area: `C1` (**CONFIRMED**)
- Space Ranger: `1.0.0` (**CONFIRMED**)
- Status: `PARTIAL`

## 3. Vizgen result

| Sample ID | Raw input | Local checksum | Download route | Redistribution |
|---|---|---|---|---|
| HumanBreastCancerPatient1 | confirmed | confirmed | confirmed | terms ambiguous |
| HumanColonCancerPatient1 | confirmed | confirmed | confirmed | terms ambiguous |
| HumanLungCancerPatient1 | confirmed | confirmed | confirmed | terms ambiguous |
| HumanMelanomaPatient1 | confirmed | confirmed | confirmed | terms ambiguous |
| HumanMelanomaPatient2 | confirmed | confirmed | confirmed | terms ambiguous |

Provider-issued checksums and historical download receipts/dates were not preserved.

## 4. hECA result

| Reference | Actual input | Lineage | Cell-level source mapping | Local checksum | Safe claim |
|---|---|---|---|---|---|
| LR03_BRCA_FFPE | confirmed | CONFIRMED_EXACT_LINEAGE | confirmed | confirmed | retained hECA breast input |
| LR04_BRCA_FF | confirmed | CONFIRMED_EXACT_LINEAGE | confirmed | confirmed | retained hECA breast input |
| LR05_BRCA_ILC | confirmed | CONFIRMED_EXACT_LINEAGE | confirmed | confirmed | retained hECA breast input |
| LR06_CERVICAL | confirmed | CONFIRMED_EXACT_LINEAGE | confirmed | confirmed | retained hECA uterus input |
| LR08_INTESTINE | confirmed | CONFIRMED_EXACT_LINEAGE | confirmed | confirmed | retained hECA intestine input |

This classification describes exact lineage within the retained input. It does not assert byte identity with the current compressed Zenodo archive.

## 5. Minor metadata

- `CID4465`: sample identity confirmed; slide and capture area remain blank (`PARTIAL`).
- Colorectal section: slide `V10A13-206`, capture `C1`, chemistry `Spatial 3' v1`, Space Ranger `1.2.0` (`CONFIRMED`).

## 6. Manifest v0.2

- Dataset records: `37`
- Experiment records: `40`
- Link integrity: `PASS`
- Status: `DRAFT_NOT_FROZEN`

## 7. Remaining blockers

### CRITICAL

- Formal manuscript/bibliography reconciliation and user-approved manifest freeze were outside this task and remain outstanding.

### MAJOR

- CTA exact Visium chemistry/kit version and original local Space Ranger bundle.
- Vizgen provider checksums and historical download receipts/dates.
- Vizgen explicit redistribution permission.

### MINOR

- CID4465 slide serial and capture area.
- Human-lung verified read-run accession.

## 8. Task status

```text
Task 4A: COMPLETED
Task 4B-1: CONDITIONAL PASS
Task 4B-2: CONDITIONAL PASS
Overall Task 4: NOT YET RESOLVED
Data availability ready: false
Dataset manifest frozen: false
```

## 9. Guardrails

```text
GitHub repository access required: false
Stage0 rerun: false
Stage1 rerun: false
Stage3A rerun: false
Stage3B rerun: false
Stage4 rerun: false
Stage5 rerun: false
CytoSPACE rerun: false
Experimental outputs modified: false
Figures modified: false
Formal manuscript accessed: false
Formal manuscript modified: false
Formal bibliography accessed: false
Formal bibliography modified: false
Task 4A outputs overwritten: false
Task 4B-1 outputs overwritten: false
```
