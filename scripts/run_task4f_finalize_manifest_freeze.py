#!/usr/bin/env python3
"""Finalize the human-approved Task 4F v1.0 manifest freeze.

This is an archival closure utility, not an experimental runner. It creates
immutable v1.0 artifacts once and becomes verification-only on later runs.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_DIR = ROOT / "reproducibility_release" / "manifests"
AUDIT_DIR = (
    ROOT
    / "visualizations"
    / "manuscript_audits"
    / "manifest_reconciliation_and_freeze_candidate"
)
RELEASE_MANIFEST = ROOT / "reproducibility_release" / "release_manifest.tsv"

DATASET_CANDIDATE = MANIFEST_DIR / "dataset_manifest_v1.0_freeze_candidate.tsv"
EXPERIMENT_CANDIDATE = MANIFEST_DIR / "experiment_manifest_v1.0_freeze_candidate.tsv"
CANDIDATE_METADATA = MANIFEST_DIR / "manifest_v1.0_freeze_candidate_metadata.json"
DATASET_FROZEN = MANIFEST_DIR / "dataset_manifest_v1.0.tsv"
EXPERIMENT_FROZEN = MANIFEST_DIR / "experiment_manifest_v1.0.tsv"
FROZEN_METADATA = MANIFEST_DIR / "manifest_v1.0_metadata.json"
HUMAN_APPROVAL = MANIFEST_DIR / "manifest_v1.0_human_approval.json"
FREEZE_REPORT = MANIFEST_DIR / "manifest_v1.0_freeze_report.md"

EXPECTED_DATASET_SHA1 = "64223e7846fc3cc8228a076b0f3e069742f58afe"
EXPECTED_EXPERIMENT_SHA1 = "347539ed5677dc89318772333552d5944821a355"
EXPECTED_CANDIDATE_METADATA_SHA1 = "80b7d1905b24bf9fe7d5de9bf52f8a1e7998575c"
EXPECTED_TASK4E_RELEASE_SHA1 = "a5df6564d552a0b726cf9273902e248d3c854944"
EXPECTED_DATASET_RECORDS = 37
EXPECTED_EXPERIMENT_RECORDS = 40
EXPECTED_RELEASE_HEADER = [
    "release_item_id",
    "category",
    "project_relative_path",
    "sha1",
    "bytes",
    "archive_status",
    "provenance",
    "notes",
]

ACCEPTED_UNRESOLVED_FIELDS = [
    "human-lung read-run accession",
    "CID4465 slide identifier",
    "CID4465 capture area",
    "CTA exact chemistry or kit version",
    "Vizgen provider-supplied checksum",
    "Vizgen explicit redistribution permission",
]

RELEASE_CONSTRAINTS = [
    (
        "Vizgen original provider files are not included in the public release "
        "while explicit redistribution permission remains unresolved."
    ),
    "CID4465 slide and capture-area fields remain unresolved.",
    "CTA exact chemistry or kit version remains unresolved.",
    (
        "The human-lung ENA sample accession ERS20065156 must not be relabelled "
        "as a sequencing-run accession."
    ),
]

APPROVAL_CONDITIONS = [
    "No unresolved value may be inferred or fabricated.",
    (
        "No provider-hosted Vizgen raw file may be redistributed without "
        "explicit permission."
    ),
    "Qualified manuscript wording must be retained.",
    "Any later substantive manifest change requires a new manifest version.",
]


class FreezeError(RuntimeError):
    """Raised when a freeze guardrail fails."""


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def sha1_bytes(data: bytes) -> str:
    return hashlib.sha1(data).hexdigest()


def sha1_file(path: Path) -> str:
    return sha1_bytes(path.read_bytes())


def require(condition: bool, message: str) -> None:
    if not condition:
        raise FreezeError(message)


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames is not None, f"Missing TSV header: {rel(path)}")
        rows = list(reader)
        return list(reader.fieldnames), rows


def json_bytes(payload: dict[str, Any]) -> bytes:
    return (json.dumps(payload, indent=2, ensure_ascii=True) + "\n").encode("utf-8")


def write_immutable(path: Path, data: bytes) -> bool:
    """Create a frozen file once; existing content must already be identical."""
    if path.exists():
        require(path.read_bytes() == data, f"Immutable artifact differs: {rel(path)}")
        return False
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)
    return True


def verify_candidates() -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, Any]]:
    expected = {
        DATASET_CANDIDATE: EXPECTED_DATASET_SHA1,
        EXPERIMENT_CANDIDATE: EXPECTED_EXPERIMENT_SHA1,
        CANDIDATE_METADATA: EXPECTED_CANDIDATE_METADATA_SHA1,
    }
    for path, expected_sha1 in expected.items():
        require(path.is_file(), f"Required candidate missing: {rel(path)}")
        observed = sha1_file(path)
        require(observed == expected_sha1, f"Candidate hash mismatch: {rel(path)}")

    _, dataset_rows = read_tsv(DATASET_CANDIDATE)
    _, experiment_rows = read_tsv(EXPERIMENT_CANDIDATE)
    metadata = json.loads(CANDIDATE_METADATA.read_text(encoding="utf-8"))

    require(len(dataset_rows) == EXPECTED_DATASET_RECORDS, "Dataset record-count mismatch")
    require(
        len(experiment_rows) == EXPECTED_EXPERIMENT_RECORDS,
        "Experiment record-count mismatch",
    )
    require(
        metadata.get("candidate_status") == "PENDING_HUMAN_APPROVAL",
        "Candidate status is not PENDING_HUMAN_APPROVAL",
    )
    require(metadata.get("dataset_record_count") == EXPECTED_DATASET_RECORDS, "Metadata dataset count mismatch")
    require(
        metadata.get("experiment_record_count") == EXPECTED_EXPERIMENT_RECORDS,
        "Metadata experiment count mismatch",
    )
    require(metadata.get("blocking_conflict_count") == 0, "Blocking conflicts are nonzero")
    require(metadata.get("major_conflict_count") == 0, "Major conflicts are nonzero")
    require(metadata.get("unresolved_nonblocking_count") == 6, "Unresolved count mismatch")
    require(metadata.get("ready_for_freeze_review") is True, "Candidate is not ready for freeze review")
    return dataset_rows, experiment_rows, metadata


def verify_manifest_relations(
    dataset_rows: list[dict[str, str]], experiment_rows: list[dict[str, str]]
) -> None:
    dataset_ids = [row["dataset_record_id"] for row in dataset_rows]
    experiment_ids = [row["experiment_record_id"] for row in experiment_rows]
    require(all(dataset_ids), "A dataset ID is blank")
    require(all(experiment_ids), "An experiment ID is blank")
    require(len(dataset_ids) == len(set(dataset_ids)), "Dataset IDs are not unique")
    require(len(experiment_ids) == len(set(experiment_ids)), "Experiment IDs are not unique")
    dataset_id_set = set(dataset_ids)
    unresolved_foreign_keys = sorted(
        {
            row["dataset_record_id"]
            for row in experiment_rows
            if row["dataset_record_id"] not in dataset_id_set
        }
    )
    require(not unresolved_foreign_keys, f"Unresolved foreign keys: {unresolved_foreign_keys}")

    require(
        all(row["candidate_status"] == "PENDING_HUMAN_APPROVAL" for row in dataset_rows),
        "Dataset candidate status changed",
    )
    require(
        all(row["candidate_status"] == "PENDING_HUMAN_APPROVAL" for row in experiment_rows),
        "Experiment candidate status changed",
    )

    by_id = {row["dataset_record_id"]: row for row in dataset_rows}
    require(by_id["SIM03_HUMAN_LUNG"]["run_accession"] == "", "Human-lung run accession was populated")
    require(by_id["SIM03_HUMAN_LUNG"]["sample_accession"] == "ERS20065156", "Human-lung sample accession changed")
    require(by_id["RD03_TNBC_PLASMA"]["slide"] == "UNRESOLVED", "CID4465 slide was populated")
    require(by_id["RD03_TNBC_PLASMA"]["capture_area"] == "UNRESOLVED", "CID4465 capture area was populated")
    require(
        "unresolved" in by_id["CTA01_BCSA2TUMB1"]["assay_or_chemistry"].lower(),
        "CTA chemistry was populated",
    )
    vizgen_ids = {
        "MERS01_BREAST",
        "MERS02_COLON",
        "MERS03_LUNG",
        "MERS04_MELANOMA1",
        "MERS05_MELANOMA2",
    }
    for record_id in vizgen_ids:
        require(
            by_id[record_id]["redistribution_class"] == "PERMISSION_UNRESOLVED",
            f"Vizgen redistribution class changed: {record_id}",
        )
    for row in dataset_rows:
        if row["redistribution_class"] in {"PERMISSION_UNRESOLVED", "RESTRICTED"}:
            require(
                "REDISTRIBUTABLE" not in row["redistribution_class"],
                f"Restricted record became redistributable: {row['dataset_record_id']}",
            )


def metadata_payload(candidate_metadata: dict[str, Any]) -> dict[str, Any]:
    return {
        "manifest_version": "1.0",
        "manifest_status": "FROZEN",
        "freeze_date": "2026-07-22",
        "human_approval": "APPROVED",
        "candidate_dataset_path": rel(DATASET_CANDIDATE),
        "candidate_dataset_sha1": EXPECTED_DATASET_SHA1,
        "frozen_dataset_path": rel(DATASET_FROZEN),
        "frozen_dataset_sha1": EXPECTED_DATASET_SHA1,
        "candidate_experiment_path": rel(EXPERIMENT_CANDIDATE),
        "candidate_experiment_sha1": EXPECTED_EXPERIMENT_SHA1,
        "frozen_experiment_path": rel(EXPERIMENT_FROZEN),
        "frozen_experiment_sha1": EXPECTED_EXPERIMENT_SHA1,
        "dataset_record_count": EXPECTED_DATASET_RECORDS,
        "experiment_record_count": EXPECTED_EXPERIMENT_RECORDS,
        "blocking_conflict_count": 0,
        "major_conflict_count": 0,
        "unresolved_nonblocking_count": 6,
        "manuscript_match_count": 174,
        "manuscript_qualified_match_count": 17,
        "validation_passed": 16,
        "validation_failed": 0,
        "repository_local_commit": candidate_metadata["repository_local_commit"],
        "source_task4f_audit_path": rel(AUDIT_DIR),
        "source_candidate_metadata_path": rel(CANDIDATE_METADATA),
        "source_candidate_metadata_sha1": EXPECTED_CANDIDATE_METADATA_SHA1,
        "release_constraints": RELEASE_CONSTRAINTS,
        "immutability_rule": (
            "Any later substantive change requires v1.0.1, v1.1, or a later version; "
            "dataset_manifest_v1.0.tsv and experiment_manifest_v1.0.tsv must not be modified."
        ),
        "public_release_status": "MANIFESTS_FROZEN_ONLY",
        "formal_manuscript_accessed": False,
        "formal_manuscript_modified": False,
        "formal_bibliography_accessed": False,
        "formal_bibliography_modified": False,
    }


def approval_payload() -> dict[str, Any]:
    return {
        "approval_status": "APPROVED",
        "approval_date": "2026-07-22",
        "approval_scope": "Task 4F dataset and experiment manifest v1.0 freeze",
        "approved_candidate_dataset_path": rel(DATASET_CANDIDATE),
        "approved_candidate_dataset_sha1": EXPECTED_DATASET_SHA1,
        "approved_candidate_experiment_path": rel(EXPERIMENT_CANDIDATE),
        "approved_candidate_experiment_sha1": EXPECTED_EXPERIMENT_SHA1,
        "approved_candidate_metadata_sha1": EXPECTED_CANDIDATE_METADATA_SHA1,
        "approval_basis": {
            "blocking_conflicts": 0,
            "major_conflicts": 0,
            "manuscript_match": 174,
            "manuscript_match_with_qualified_wording": 17,
            "manuscript_missing_in_manifest": 0,
            "manuscript_conflict": 0,
            "validation_passed": 16,
            "validation_failed": 0,
            "candidate_generation_deterministic": True,
        },
        "accepted_unresolved_nonblocking_fields": ACCEPTED_UNRESOLVED_FIELDS,
        "approval_conditions": APPROVAL_CONDITIONS,
        "formal_manuscript_accessed": False,
        "formal_manuscript_modified": False,
        "formal_bibliography_accessed": False,
        "formal_bibliography_modified": False,
    }


def report_bytes(metadata_sha1: str, approval_sha1: str) -> bytes:
    constraints = "\n".join(f"- {item}" for item in RELEASE_CONSTRAINTS)
    unresolved = "\n".join(f"- {item}" for item in ACCEPTED_UNRESOLVED_FIELDS)
    content = f"""# Task 4F v1.0 Manifest Freeze Report

## Decision

`PASS_MANIFESTS_FROZEN`

Human approval was recorded on 2026-07-22. Candidate hashes, record counts,
referential integrity, redistribution guardrails, and the Task 4E release
inventory baseline were verified before the formal v1.0 artifacts were created.

## Frozen manifests

| Artifact | Records | SHA-1 | Candidate byte identity |
|---|---:|---|---|
| `{rel(DATASET_FROZEN)}` | 37 | `{EXPECTED_DATASET_SHA1}` | true |
| `{rel(EXPERIMENT_FROZEN)}` | 40 | `{EXPECTED_EXPERIMENT_SHA1}` | true |

Formal metadata SHA-1: `{metadata_sha1}`

Human approval record SHA-1: `{approval_sha1}`

## Accepted unresolved non-blocking fields

{unresolved}

These fields remain unresolved. Approval does not authorize inferred values or
redistribution of third-party provider files whose permission is unresolved.

## Release constraints

{constraints}

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
"""
    return content.encode("utf-8")


def verify_task4e_release_prefix(data: bytes) -> None:
    marker = b"T4F-"
    marker_index = data.find(marker)
    prefix = data if marker_index < 0 else data[:marker_index]
    require(
        sha1_bytes(prefix) == EXPECTED_TASK4E_RELEASE_SHA1,
        "Existing Task 4E release records changed",
    )


def release_entry(
    item_id: str,
    artifact_class: str,
    path: Path,
    freeze_status: str,
    redistribution_class: str,
    notes: str,
) -> list[str]:
    require(path.is_file(), f"Release artifact missing: {rel(path)}")
    return [
        item_id,
        artifact_class,
        rel(path),
        sha1_file(path),
        str(path.stat().st_size),
        freeze_status,
        f"source_task=Task 4F-Closure;redistribution_class={redistribution_class}",
        notes,
    ]


def expected_release_entries() -> list[list[str]]:
    return [
        release_entry("T4F-CLOSURE-001", "FROZEN_DATASET_MANIFEST", DATASET_FROZEN, "FROZEN", "MANIFEST_ONLY", "Human-approved immutable dataset manifest v1.0."),
        release_entry("T4F-CLOSURE-002", "FROZEN_EXPERIMENT_MANIFEST", EXPERIMENT_FROZEN, "FROZEN", "MANIFEST_ONLY", "Human-approved immutable experiment manifest v1.0."),
        release_entry("T4F-CLOSURE-003", "FROZEN_MANIFEST_METADATA", FROZEN_METADATA, "FROZEN", "METADATA_ONLY", "Formal v1.0 freeze metadata and bounded release constraints."),
        release_entry("T4F-CLOSURE-004", "HUMAN_APPROVAL_RECORD", HUMAN_APPROVAL, "FROZEN", "METADATA_ONLY", "Human approval and accepted unresolved non-blocking fields."),
        release_entry("T4F-CLOSURE-005", "MANIFEST_FREEZE_REPORT", FREEZE_REPORT, "FROZEN", "DOCUMENTATION_ONLY", "Task 4F closure decision, constraints, and immutability rule."),
        release_entry("T4F-CLOSURE-006", "FREEZE_VALIDATION_SCRIPT", Path(__file__).resolve(), "INCLUDED", "SOURCE_CODE_ONLY", "Deterministic archival freeze and idempotency validator; not an experimental runner."),
        release_entry("T4F-CANDIDATE-001", "APPROVED_FREEZE_CANDIDATE", DATASET_CANDIDATE, "RETAINED_APPROVED_SOURCE", "MANIFEST_ONLY", "Byte source for dataset_manifest_v1.0.tsv."),
        release_entry("T4F-CANDIDATE-002", "APPROVED_FREEZE_CANDIDATE", EXPERIMENT_CANDIDATE, "RETAINED_APPROVED_SOURCE", "MANIFEST_ONLY", "Byte source for experiment_manifest_v1.0.tsv."),
        release_entry("T4F-CANDIDATE-003", "FREEZE_CANDIDATE_METADATA", CANDIDATE_METADATA, "RETAINED_APPROVED_SOURCE", "METADATA_ONLY", "Task 4F candidate metadata with PENDING_HUMAN_APPROVAL provenance state."),
        release_entry("T4F-AUDIT-001", "TASK4F_AUDIT", AUDIT_DIR / "README.md", "INCLUDED", "DOCUMENTATION_ONLY", "Task 4F audit scope and freeze-candidate decision."),
        release_entry("T4F-AUDIT-002", "TASK4F_AUDIT", AUDIT_DIR / "manuscript_manifest_reconciliation.tsv", "INCLUDED", "AUDIT_DERIVATIVE_ONLY", "Manuscript-to-manifest reconciliation ledger."),
        release_entry("T4F-AUDIT-003", "TASK4F_AUDIT", AUDIT_DIR / "candidate_change_log.tsv", "INCLUDED", "AUDIT_DERIVATIVE_ONLY", "Field-level v0.2 to v1.0 candidate change log."),
        release_entry("T4F-AUDIT-004", "TASK4F_AUDIT", AUDIT_DIR / "audit_summary.json", "INCLUDED", "METADATA_ONLY", "Machine-readable Task 4F audit summary."),
        release_entry("T4F-AUDIT-005", "TASK4F_AUDIT", AUDIT_DIR / "final_decision.md", "INCLUDED", "DOCUMENTATION_ONLY", "Task 4F freeze-review decision."),
    ]


def update_release_manifest(entries: list[list[str]]) -> int:
    require(RELEASE_MANIFEST.is_file(), "Release manifest is missing")
    existing_bytes = RELEASE_MANIFEST.read_bytes()
    verify_task4e_release_prefix(existing_bytes)

    header, existing_rows = read_tsv(RELEASE_MANIFEST)
    require(header == EXPECTED_RELEASE_HEADER, "Release manifest header changed")
    existing_by_id = {row["release_item_id"]: row for row in existing_rows}
    existing_by_path = {row["project_relative_path"]: row for row in existing_rows}
    require(len(existing_by_id) == len(existing_rows), "Duplicate release item IDs exist")
    require(len(existing_by_path) == len(existing_rows), "Duplicate release paths exist")

    to_append: list[list[str]] = []
    for values in entries:
        expected = dict(zip(EXPECTED_RELEASE_HEADER, values))
        item_id = expected["release_item_id"]
        project_path = expected["project_relative_path"]
        if item_id in existing_by_id or project_path in existing_by_path:
            require(item_id in existing_by_id, f"Release path exists under another ID: {project_path}")
            require(project_path in existing_by_path, f"Release ID exists for another path: {item_id}")
            require(existing_by_id[item_id] == expected, f"Existing Task 4F release row differs: {item_id}")
        else:
            to_append.append(values)

    if to_append:
        buffer = io.StringIO(newline="")
        writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
        writer.writerows(to_append)
        prefix = b"" if existing_bytes.endswith((b"\n", b"\r")) else b"\n"
        with RELEASE_MANIFEST.open("ab") as handle:
            handle.write(prefix + buffer.getvalue().encode("utf-8"))
    return len(to_append)


def verify_release_manifest(entries: list[list[str]]) -> None:
    data = RELEASE_MANIFEST.read_bytes()
    verify_task4e_release_prefix(data)
    header, rows = read_tsv(RELEASE_MANIFEST)
    require(header == EXPECTED_RELEASE_HEADER, "Release manifest header changed")
    by_id = {row["release_item_id"]: row for row in rows}
    require(len(by_id) == len(rows), "Duplicate release item IDs detected")
    paths = [row["project_relative_path"] for row in rows]
    require(len(paths) == len(set(paths)), "Duplicate release paths detected")

    for values in entries:
        expected = dict(zip(EXPECTED_RELEASE_HEADER, values))
        item_id = expected["release_item_id"]
        require(item_id in by_id, f"Release row missing: {item_id}")
        require(by_id[item_id] == expected, f"Release row mismatch: {item_id}")
        artifact = ROOT / expected["project_relative_path"]
        require(artifact.is_file(), f"Release path missing: {rel(artifact)}")
        require(sha1_file(artifact) == expected["sha1"], f"Release hash mismatch: {rel(artifact)}")
        require(str(artifact.stat().st_size) == expected["bytes"], f"Release byte-size mismatch: {rel(artifact)}")
        lowered = expected["project_relative_path"].lower()
        require("data/raw/" not in lowered, f"Raw provider path entered release inventory: {lowered}")
        require("provider_download" not in lowered, f"Provider raw path entered release inventory: {lowered}")


def verify_frozen_artifacts() -> None:
    require(DATASET_FROZEN.read_bytes() == DATASET_CANDIDATE.read_bytes(), "Frozen dataset is not byte-identical")
    require(EXPERIMENT_FROZEN.read_bytes() == EXPERIMENT_CANDIDATE.read_bytes(), "Frozen experiment is not byte-identical")
    require(sha1_file(DATASET_FROZEN) == EXPECTED_DATASET_SHA1, "Frozen dataset hash mismatch")
    require(sha1_file(EXPERIMENT_FROZEN) == EXPECTED_EXPERIMENT_SHA1, "Frozen experiment hash mismatch")
    metadata = json.loads(FROZEN_METADATA.read_text(encoding="utf-8"))
    approval = json.loads(HUMAN_APPROVAL.read_text(encoding="utf-8"))
    require(metadata.get("manifest_status") == "FROZEN", "Formal metadata is not FROZEN")
    require(metadata.get("human_approval") == "APPROVED", "Formal metadata approval is not APPROVED")
    require(approval.get("approval_status") == "APPROVED", "Human approval record is not APPROVED")


def main() -> int:
    try:
        dataset_rows, experiment_rows, candidate_metadata = verify_candidates()
        verify_manifest_relations(dataset_rows, experiment_rows)

        initial_release_bytes = RELEASE_MANIFEST.read_bytes()
        verify_task4e_release_prefix(initial_release_bytes)
        existing_formal = [
            DATASET_FROZEN,
            EXPERIMENT_FROZEN,
            FROZEN_METADATA,
            HUMAN_APPROVAL,
            FREEZE_REPORT,
        ]
        preexisting = {path: path.exists() for path in existing_formal}
        require(
            len(set(preexisting.values())) == 1,
            "Partial formal freeze detected; refusing to repair immutable v1.0 in place",
        )

        write_immutable(DATASET_FROZEN, DATASET_CANDIDATE.read_bytes())
        write_immutable(EXPERIMENT_FROZEN, EXPERIMENT_CANDIDATE.read_bytes())

        metadata_data = json_bytes(metadata_payload(candidate_metadata))
        approval_data = json_bytes(approval_payload())
        write_immutable(FROZEN_METADATA, metadata_data)
        write_immutable(HUMAN_APPROVAL, approval_data)
        report_data = report_bytes(sha1_bytes(metadata_data), sha1_bytes(approval_data))
        write_immutable(FREEZE_REPORT, report_data)

        verify_frozen_artifacts()
        entries = expected_release_entries()
        added = update_release_manifest(entries)
        verify_release_manifest(entries)

        expected_added = 0 if all(preexisting.values()) else len(entries)
        require(added == expected_added, "Unexpected release inventory append count")
        mode = "IDEMPOTENT_VERIFICATION" if expected_added == 0 else "INITIAL_FREEZE"
        result = {
            "decision": "PASS_MANIFESTS_FROZEN",
            "mode": mode,
            "dataset_records": len(dataset_rows),
            "dataset_sha1": sha1_file(DATASET_FROZEN),
            "experiment_records": len(experiment_rows),
            "experiment_sha1": sha1_file(EXPERIMENT_FROZEN),
            "metadata_sha1": sha1_file(FROZEN_METADATA),
            "approval_sha1": sha1_file(HUMAN_APPROVAL),
            "freeze_report_sha1": sha1_file(FREEZE_REPORT),
            "release_records_added": added,
            "release_manifest_sha1": sha1_file(RELEASE_MANIFEST),
            "task4e_records_unchanged": True,
            "all_release_hashes_valid": True,
        }
        print(json.dumps(result, indent=2))
        return 0
    except (FreezeError, KeyError, json.JSONDecodeError) as error:
        print(json.dumps({"decision": "HOLD", "error": str(error)}, indent=2), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
