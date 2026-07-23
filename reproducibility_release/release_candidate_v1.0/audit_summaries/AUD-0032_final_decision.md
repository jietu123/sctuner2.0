# Task 4F final decision

**Decision:** `CONDITIONAL_PASS_TO_FREEZE_REVIEW`  
**Candidate status:** `PENDING_HUMAN_APPROVAL`  
**Ready for final manifest freeze:** `false`

## Basis

- Current v0.2 manifests reconcile to exactly 37 dataset and 40 experiment records; no record was added, deleted, split or merged.
- All source-layer and redistribution classifications use the frozen Task 4F vocabularies.
- All manuscript fact checks avoid overclaiming unresolved metadata.
- Reconciliation counts: `{'MATCH': 174, 'MATCH_WITH_QUALIFIED_WORDING': 17}`.
- Blocking conflicts: `0`.
- Major unexplained manuscript-manifest conflicts: `0`.
- Explicit unresolved non-blocking records: `6`.
- Candidate dataset SHA-1: `64223e7846fc3cc8228a076b0f3e069742f58afe`.
- Candidate experiment SHA-1: `347539ed5677dc89318772333552d5944821a355`.

## Why conditional

Major provenance gaps remain for Vizgen provider checksums/redistribution permission and the exact CTA chemistry. The manuscript fact ledger already uses qualified wording, the candidate records preserve `PARTIAL` where appropriate, and no raw provider file is declared redistributable. These gaps do not block freeze review but require human acceptance before final freeze.

## Required next action

Review the two candidate manifests, `manuscript_manifest_reconciliation.tsv`, `candidate_change_log.tsv`, `audit_summary.json` and this decision. Do not execute the final freeze until human approval is recorded.

## Guardrails

- Formal manuscript accessed: false
- Formal manuscript modified: false
- Formal bibliography accessed: false
- Formal bibliography modified: false
- Experimental code modified: false
- Experimental outputs modified: false
- Formal stages rerun: false
- Formal figures modified: false
- Formal source-value tables modified: false
- GitHub remote access required: false
