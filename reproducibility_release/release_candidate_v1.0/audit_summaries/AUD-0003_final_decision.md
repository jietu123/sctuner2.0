# Final decision

## 1. Decision

**HOLD**

The project is on HOLD for submission-level dataset provenance. Core source identities are mostly recoverable, but the formal manuscript, formal bibliography and unified formal dataset manifest are absent, and one accession-type conflict remains.

## 2. Audit coverage

| Family | Expected | Audited | Confirmed | Partial | Project-generated | Unresolved | Conflict |
|---|---:|---:|---:|---:|---:|---:|---:|
| Joint simulations | 9 | 9 | 6 | 0 | 0 | 0 | 3 |
| Low-resolution profile masking | 10 | 10 | 5 | 5 | 0 | 0 | 0 |
| MERSCOPE | 5 | 5 | 0 | 5 | 0 | 0 | 0 |
| State-decoy experiments | 2 | 2 | 0 | 0 | 2 | 0 | 0 |
| Real reference dropout | 6 | 6 | 6 | 0 | 0 | 0 | 0 |
| Detailed Thalamic case | 1 | 1 | 1 | 0 | 0 | 0 | 0 |
| CTA biological application | 1 | 1 | 0 | 1 | 0 | 0 | 0 |

## 3. Confirmed canonical records

Confirmed canonical identities include the 10x breast/kidney/brain/embryo releases and slide/capture-area identifiers listed in the master table; Wu GSE176078 including CID45171/GSM5354535; cell2location E-MTAB-11114/11115; Ransick GSE129798/SRP192559/PRJNA532850; Lee GSE132465; the six real-dropout experiment identities; ST8059051 = Visium-29B = C05717-021/B1; and the CTA article plus raw/processed repository records.

## 4. Unresolved records

- **CRITICAL U001** (FORMAL_ARTIFACTS): Formal LaTeX source was not found. Required resolution: Provide actual compilation directory and rerun occurrence scan.
- **CRITICAL U002** (FORMAL_ARTIFACTS): Only candidate audit BibTeX files exist. Required resolution: Provide the compiled bibliography and run key-level conflict audit.
- **CRITICAL U003** (FORMAL_ARTIFACTS): No unified formal dataset manifest is present; dispersed experiment manifests are not equivalent. Required resolution: Create and freeze a formal manifest only after this audit is reviewed.
- **MAJOR U004** (LR03_BRCA_FFPE;LR04_BRCA_FF;LR05_BRCA_ILC): The current public hECA project export and source study are identified, but the historical local file lacks a preserved byte-for-byte merge/checksum manifest. Required resolution: Recover the original export/merge log or describe the reference as a curated hECA export without overclaiming exact component provenance.
- **MAJOR U005** (LR06_CERVICAL): Source studies are documented, but exact contributing-cell/accession membership of the local historical export is incomplete. Required resolution: Recover a cell-level source manifest.
- **MAJOR U006** (LR08_INTESTINE): Source studies are documented, but exact contributing-cell/accession membership of the local historical export is incomplete. Required resolution: Recover a cell-level source manifest.
- **MAJOR U007** (MERS01_BREAST;MERS02_COLON;MERS03_LUNG;MERS04_MELANOMA1;MERS05_MELANOMA2): Provider release and release sample IDs are confirmed, but repository lacks preserved license text and stable per-sample download URLs/checksums. Required resolution: Archive provider terms and original download receipt/URL/checksum before Data availability is finalized.
- **MAJOR U008** (CTA01_BCSA2TUMB1): CTA paper confirms frozen Visium sections but does not expose the exact chemistry version or slide serial/capture area for BCSA2TumB1 in the evidence currently present. Required resolution: Recover sample sheet or raw Space Ranger metadata.
- **MINOR U009** (RD03_TNBC_PLASMA): CID4465 identity and Zenodo source are confirmed, but slide serial/capture area are not exposed. Required resolution: Leave fields unresolved unless primary metadata is recovered.
- **MINOR U010** (RD04_CRC_B): Formal 10x release and capture area C1 are confirmed, but slide serial is not. Required resolution: Leave slide serial unresolved unless provider metadata is recovered.
- **MAJOR U011** (SIM03_HUMAN_LUNG): ERS20065156 is a sample accession and no linked run accession was recovered from the ENA read-run query. Required resolution: Do not invent a run accession; report project, study, sample and BioSample identifiers only.

## 5. Conflicts

- **C001** (SIM03_HUMAN_LUNG): ERS20065156 accession type. Canonical recommendation: ERS20065156 = ENA sample accession; SAMEA115633909 = BioSample; run accession unresolved.
- **C002** (FORMAL_ARTIFACTS): submission readiness. Canonical recommendation: Task 4A decision HOLD until formal artifacts are supplied and audited.

## 6. Highest-priority manuscript risks

- **CRITICAL:** The formal LaTeX and compiled BibTeX are absent, so Results/Methods/legends/key consistency cannot be certified.
- **CRITICAL:** No unified formal dataset manifest exists; future-tense manifest language is not yet supportable.
- **CRITICAL:** `ERS20065156` is mis-typed as a run in the current claim; ENA identifies it as a sample.
- **MAJOR:** Historical hECA composite exports lack preserved cell-level merge/checksum manifests.
- **MAJOR:** Vizgen per-sample download receipts/license terms and CTA exact Visium chemistry/slide metadata remain incomplete.

## 7. Required future manuscript changes

Changes are listed only in `manuscript_required_updates.tsv`. None were applied in this task.

## 8. Data availability readiness

**Not ready.** The inventory is sufficiently structured for drafting, but formal manuscript linkage, Vizgen terms, CTA exact platform metadata, composite-reference lineage, and processed-derivative deposition/location must be resolved first.

## 9. Dataset manifest readiness

**Not ready to freeze.** Review this audit, resolve CRITICAL/MAJOR records, supply the formal manuscript/BibTeX, then generate a separate formal manifest with approval.

## 10. Guardrails

```text
Stage1 rerun: false
Stage3A rerun: false
Stage3B rerun: false
Stage4 rerun: false
CytoSPACE rerun: false
Experimental outputs modified: false
Formal manuscript modified: false
Formal bibliography modified: false
Formal manifest overwritten: false
Figures modified: false
```
