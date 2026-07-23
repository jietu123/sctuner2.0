# Task 4B-2 Major Provenance Evidence Recovery

This directory is an audit-only, reproducible update of the Task 4B-1 draft manifest.

- Decision: **CONDITIONAL PASS**
- Dataset manifest: **37 records**, version **0.2**, `DRAFT_NOT_FROZEN`
- Experiment manifest: **40 records**
- Link integrity: **PASS**
- Formal manuscript/BibTeX accessed: **false**
- Experimental outputs or figures modified: **false**

## Recovered evidence

- CTA `BCSA2TumB1`: official metadata maps the sample to slide `V10F24-112`, capture area `C1`; exact chemistry version remains unresolved.
- Vizgen: all five exact provider sample routes and local input cell totals are confirmed; provider checksums, historical receipts, and unambiguous redistribution permission remain unresolved.
- hECA: all five retained inputs have exact local SHA-1, cell counts, and cell-level `study_id`/`donor_ID` lineage.
- CRC: slide `V10A13-206`, capture `C1`, Spatial 3' v1, and Space Ranger 1.2.0 are confirmed.
- CID4465: sample identity is confirmed; slide/capture remain blank.

The official values were verified on 2026-07-22. The generator records URLs but does not require network access during deterministic regeneration.
