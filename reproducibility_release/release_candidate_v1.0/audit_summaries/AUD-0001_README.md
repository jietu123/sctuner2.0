# Task 4A dataset provenance formal audit

This directory is an independent, read-only provenance audit generated on 2026-07-22. It does not replace a formal dataset manifest.

## Decision

**HOLD**

The repository contains extensive experiment-level provenance, but the formal manuscript `.tex`, compiled `sn-bibliography.bib`, and unified formal dataset manifest were not present. A submission-level manuscript/BibTeX/manifest consistency check is therefore impossible. The human-lung record also contains a material accession-type conflict: `ERS20065156` is an ENA sample accession, not a run accession.

## Scope

- 37 provenance master records.
- 40 experiment crosswalk rows.
- Required coverage: 9 joint simulations, 10 low-resolution scenarios, 5 MERSCOPE samples, 2 state-decoy experiments, 6 real reference-dropout settings, one detailed thalamic case, and one CTA biological application.
- Additional external tumour pairings and EcoTyper, LIANA/OmniPath, and CytoSPACE supplementary resources are classified separately.

## Evidence hierarchy

Local formal inputs and frozen experiment manifests were read first, followed by preprocessing scripts/configs, official provider/database pages, source publications, and finally directory names. Existing audits were treated as evidence, not as formal manifests.

## Important terminology locks

- `ERS20065156`: ENA **sample** accession.
- `SAMEA115633909`: BioSample linked to that sample.
- MERSCOPE: platform; MERFISH: assay.
- `Human...Patient1`: Vizgen release sample identifier, not a verified clinical patient ID or repository accession.
- Simulation truth, disjoint same-assay references, decoys, CTA-to-spot endpoint, and the 4,014-cell Wu subset are project-generated derivatives.
- ST8059051/C05717-021/B1 is distinct from the V52B25-081/B FFPE CytAssist scaffold.

## Guardrails

No experimental stage was rerun. No manuscript, bibliography, formal manifest, source-value table, result, or figure was modified.
