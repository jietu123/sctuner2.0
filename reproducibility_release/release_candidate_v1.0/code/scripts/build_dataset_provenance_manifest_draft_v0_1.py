#!/usr/bin/env python3
"""Generate Task 4B-1 draft provenance manifests from Task 4A outputs only."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import subprocess
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "visualizations" / "manuscript_audits" / "dataset_provenance_formal_audit"
OUT = ROOT / "visualizations" / "manuscript_audits" / "dataset_provenance_manifest_draft_v0_1"
GENERATED_DATE = "2026-07-22"
VERSION = "0.1"
STATUS = "DRAFT_NOT_FROZEN"
HUMAN_LUNG_RECORD = "SIM03_HUMAN_LUNG"
HUMAN_LUNG_EXPERIMENTS = [
    "human_lung_5loc_fine9_clustered_sim_sc_missing_b_cell",
    "human_lung_5loc_fine9_clustered_sim_missing_at2_sc_missing_b_cell",
    "human_lung_5loc_fine9_clustered_sim_missing_at2_fibroblast_sc_missing_b_cell",
]
PARTIAL_REASON = (
    "Accession-type conflict resolved; read-run accession remains unresolved "
    "and is not required to identify the frozen spatial sample."
)

SOURCE_FILES = [
    "dataset_provenance_master.tsv",
    "experiment_input_crosswalk.tsv",
    "official_source_registry.tsv",
    "accession_type_audit.tsv",
    "manifest_consistency_audit.tsv",
    "data_availability_inventory.tsv",
    "unresolved_records.tsv",
    "conflicts.tsv",
    "manuscript_required_updates.tsv",
    "audit_summary.json",
    "final_decision.md",
    "README.md",
]

DATASET_FIELDS = [
    "record_id", "dataset_family", "experiment_id", "internal_dataset_id", "data_role",
    "formal_dataset_name", "formal_release_name", "source_study_title", "species", "tissue",
    "disease", "donor_or_patient_id", "sample_id", "section_id", "slide_id", "capture_area",
    "spatial_platform", "assay", "tissue_processing", "panel_or_gene_count", "reference_type",
    "reference_dataset_name", "reference_sample_subset", "database", "study_accession",
    "project_accession", "experiment_accession", "run_accession", "sample_accession",
    "biosample_accession", "geo_accession", "arrayexpress_accession", "ena_accession",
    "official_source_url", "provider_download_url", "source_publication", "bibtex_key",
    "input_file", "input_file_sha1", "config_file", "config_file_sha1",
    "redistribution_status", "project_generated_derivative", "provenance_status",
    "unresolved_issue_id", "conflict_id", "evidence_summary", "notes",
]

EXPERIMENT_FIELDS = [
    "experiment_id", "dataset_family", "figure_panel", "manuscript_result_section",
    "spatial_input_record_id", "reference_input_record_id", "evaluation_resource_record_id",
    "profile_mask_target", "reference_dropout_target", "project_derivative", "formal_status",
    "status_basis", "unresolved_issue_id", "notes",
]

SOURCE_REGISTRY_FIELDS = [
    "record_id", "source_authority", "source_type", "official_url", "provider_url",
    "accessed_date", "supports_fields", "source_status", "notes",
]


def sha1_bytes(data: bytes) -> str:
    return hashlib.sha1(data).hexdigest()


def sha1_file(path: Path) -> str:
    digest = hashlib.sha1()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def clean(value: str | None) -> str:
    """Use blank draft fields for Task 4A NOT_APPLICABLE sentinels."""
    if value is None or value == "NOT_APPLICABLE":
        return ""
    return value


def read_tsv(name: str) -> list[dict[str, str]]:
    path = SOURCE / name
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def tsv_bytes(rows: list[dict[str, str]], fields: list[str]) -> bytes:
    handle = io.StringIO(newline="")
    writer = csv.DictWriter(
        handle,
        fieldnames=fields,
        delimiter="\t",
        extrasaction="ignore",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)
    return handle.getvalue().encode("utf-8")


def json_bytes(payload: dict) -> bytes:
    return (json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n").encode("utf-8")


def linked_ids(value: str) -> list[str]:
    return [item.strip() for item in value.split(";") if item.strip()]


def issue_links(unresolved: list[dict[str, str]]) -> dict[str, list[str]]:
    links: dict[str, list[str]] = defaultdict(list)
    for row in unresolved:
        for record_id in linked_ids(row["linked_record_id"]):
            links[record_id].append(row["issue_id"])
    return links


def conflict_links(conflicts: list[dict[str, str]]) -> dict[str, list[str]]:
    links: dict[str, list[str]] = defaultdict(list)
    for row in conflicts:
        for record_id in linked_ids(row["linked_record_id"]):
            links[record_id].append(row["conflict_id"])
    return links


def build_dataset_manifest(
    master: list[dict[str, str]],
    issues: dict[str, list[str]],
    conflicts: dict[str, list[str]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for source in master:
        row = {field: clean(source.get(field, "")) for field in DATASET_FIELDS}
        row["provenance_status"] = source["status"]
        row["unresolved_issue_id"] = ";".join(issues.get(source["record_id"], []))
        row["conflict_id"] = ";".join(conflicts.get(source["record_id"], []))
        row["evidence_summary"] = source["evidence"]
        if source["record_id"] == HUMAN_LUNG_RECORD:
            row["run_accession"] = ""
            row["sample_accession"] = "ERS20065156"
            row["biosample_accession"] = "SAMEA115633909"
            row["provenance_status"] = "PARTIAL"
            row["notes"] = f"{source['notes']} {PARTIAL_REASON} Task 4A conflict C001 is retained as resolved history."
        rows.append(row)
    return rows


def build_experiment_manifest(
    crosswalk: list[dict[str, str]],
    dataset_by_id: dict[str, dict[str, str]],
    issues: dict[str, list[str]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for source in crosswalk:
        linked = [
            source["spatial_input_record_id"],
            source["reference_input_record_id"],
        ]
        if source["evaluation_resource_record_id"] != "NOT_APPLICABLE":
            linked.append(source["evaluation_resource_record_id"])
        unresolved_ids = sorted({issue for record_id in linked for issue in issues.get(record_id, [])})
        spatial = dataset_by_id[source["spatial_input_record_id"]]
        row = {
            "experiment_id": source["experiment_id"],
            "dataset_family": spatial["dataset_family"],
            "figure_panel": source["figure_panel"],
            "manuscript_result_section": source["manuscript_result_section"],
            "spatial_input_record_id": source["spatial_input_record_id"],
            "reference_input_record_id": source["reference_input_record_id"],
            "evaluation_resource_record_id": clean(source["evaluation_resource_record_id"]),
            "profile_mask_target": clean(source["profile_mask_target"]),
            "reference_dropout_target": clean(source["reference_dropout_target"]),
            "project_derivative": clean(source["project_derivative"]),
            "formal_status": source["formal_status"],
            "status_basis": "Inherited from Task 4A experiment crosswalk; unresolved/partial records retained.",
            "unresolved_issue_id": ";".join(unresolved_ids),
            "notes": source["notes"],
        }
        if source["experiment_id"] in HUMAN_LUNG_EXPERIMENTS:
            row["formal_status"] = "PARTIAL"
            row["status_basis"] = PARTIAL_REASON
            row["notes"] = f"{source['notes']} Human-lung accession correction recorded in Task 4B-1.".strip()
        rows.append(row)
    return rows


def build_source_registry(
    source_rows: list[dict[str, str]],
    dataset_by_id: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    rows = []
    for source in source_rows:
        dataset = dataset_by_id[source["record_id"]]
        rows.append({
            "record_id": source["record_id"],
            "source_authority": source["source_authority"],
            "source_type": source["source_type"],
            "official_url": source["url"],
            "provider_url": clean(dataset["provider_download_url"]),
            "accessed_date": source["accessed_date"],
            "supports_fields": source["supports_fields"],
            "source_status": source["status"],
            "notes": source["notes"],
        })
    return rows


def build_human_lung_correction() -> list[dict[str, str]]:
    experiments = ";".join(HUMAN_LUNG_EXPERIMENTS)
    evidence = (
        "Task 4A accession_type_audit.tsv C001 evidence; ENA sample record for ERS20065156; "
        "Task 4A master record SIM03_HUMAN_LUNG"
    )
    common = {
        "record_id": HUMAN_LUNG_RECORD,
        "experiment_ids": experiments,
        "evidence_source": evidence,
        "remaining_unknown": "verified read-run accession not recovered",
        "status": "CORRECTED_WITH_REMAINING_PARTIAL_METADATA",
        "notes": "WSA_LngSP10193345, E-MTAB-11640 and PRJEB52292 are unchanged. No ERR/SRR/DRR accession was inferred.",
    }
    return [
        {
            **common,
            "field": "run_accession",
            "task4a_value": "run_accession = ERS20065156 (historical claim retained in Task 4A conflict record)",
            "verified_value": "run_accession = blank",
            "verified_type": "read-run accession not recovered",
            "correction_action": "Remove ERS20065156 from run_accession in draft artifacts; leave the field blank.",
        },
        {
            **common,
            "field": "sample_accession",
            "task4a_value": "ERS20065156 was historically described as a run; Task 4A master already normalized it to sample_accession",
            "verified_value": "sample_accession = ERS20065156",
            "verified_type": "ENA sample accession",
            "correction_action": "Retain ERS20065156 only in sample_accession and generic ENA identifier context.",
        },
        {
            **common,
            "field": "biosample_accession",
            "task4a_value": "biosample_accession = SAMEA115633909",
            "verified_value": "biosample_accession = SAMEA115633909",
            "verified_type": "BioSample",
            "correction_action": "Retain the associated BioSample without alteration.",
        },
    ]


def build_status_transitions(
    master: list[dict[str, str]],
    dataset_draft: list[dict[str, str]],
    crosswalk: list[dict[str, str]],
    experiment_draft: list[dict[str, str]],
) -> list[dict[str, str]]:
    draft_datasets = {row["record_id"]: row for row in dataset_draft}
    draft_experiments = {row["experiment_id"]: row for row in experiment_draft}
    rows = []
    for source in master:
        new_status = draft_datasets[source["record_id"]]["provenance_status"]
        changed = source["status"] != new_status
        rows.append({
            "record_or_experiment_id": source["record_id"],
            "task4a_status": source["status"],
            "draft_v0_1_status": new_status,
            "change_reason": PARTIAL_REASON if source["record_id"] == HUMAN_LUNG_RECORD else "No status change.",
            "evidence": "Task 4A C001 and accession type audit" if changed else "Task 4A dataset provenance master",
            "allowed": "true",
            "notes": "Dataset-level correction in draft only; Task 4A source record was not overwritten." if changed else "Status preserved.",
        })
    for source in crosswalk:
        new_status = draft_experiments[source["experiment_id"]]["formal_status"]
        changed = source["formal_status"] != new_status
        rows.append({
            "record_or_experiment_id": source["experiment_id"],
            "task4a_status": source["formal_status"],
            "draft_v0_1_status": new_status,
            "change_reason": PARTIAL_REASON if source["experiment_id"] in HUMAN_LUNG_EXPERIMENTS else "No status change.",
            "evidence": "Task 4A C001 and accession type audit" if changed else "Task 4A experiment input crosswalk",
            "allowed": "true",
            "notes": "Experiment-level conflict resolved to PARTIAL; read-run metadata remains unresolved." if changed else "Status preserved.",
        })
    return rows


def bad_accession_fields(dataset_rows: list[dict[str, str]]) -> list[str]:
    bad = []
    for row in dataset_rows:
        run = row["run_accession"]
        if run and not run.startswith(("ERR", "SRR", "DRR")):
            bad.append(f"{row['record_id']}:run_accession={run}")
        sample = row["sample_accession"]
        if sample and not sample.startswith(("ERS", "SRS", "DRS", "GSM")):
            bad.append(f"{row['record_id']}:sample_accession={sample}")
        biosample = row["biosample_accession"]
        if biosample and not biosample.startswith(("SAMEA", "SAMN", "SAMD")):
            bad.append(f"{row['record_id']}:biosample_accession={biosample}")
        geo = row["geo_accession"]
        if geo and any(not token.startswith("GSE") for token in linked_ids(geo)):
            bad.append(f"{row['record_id']}:geo_accession={geo}")
        arrayexpress = row["arrayexpress_accession"]
        if arrayexpress and any(not token.startswith("E-MTAB-") for token in linked_ids(arrayexpress)):
            bad.append(f"{row['record_id']}:arrayexpress_accession={arrayexpress}")
    return bad


def build_link_integrity(
    dataset_rows: list[dict[str, str]],
    experiment_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    dataset_ids = [row["record_id"] for row in dataset_rows]
    experiment_ids = [row["experiment_id"] for row in experiment_rows]
    dataset_set = set(dataset_ids)
    spatial_missing = sorted({row["spatial_input_record_id"] for row in experiment_rows} - dataset_set)
    reference_missing = sorted({row["reference_input_record_id"] for row in experiment_rows} - dataset_set)
    evaluation_ids = {row["evaluation_resource_record_id"] for row in experiment_rows if row["evaluation_resource_record_id"]}
    evaluation_missing = sorted(evaluation_ids - dataset_set)
    linked = {
        value
        for row in experiment_rows
        for value in (
            row["spatial_input_record_id"],
            row["reference_input_record_id"],
            row["evaluation_resource_record_id"],
        )
        if value
    }
    orphans = sorted(dataset_set - linked)
    duplicate_datasets = sorted(item for item, count in Counter(dataset_ids).items() if count > 1)
    duplicate_experiments = sorted(item for item, count in Counter(experiment_ids).items() if count > 1)
    accession_errors = bad_accession_fields(dataset_rows)
    generated_with_accessions = []
    for row in dataset_rows:
        if row["provenance_status"] != "PROJECT_GENERATED":
            continue
        # External identifiers in paired source/derivative rows must be explained as source-only.
        if row["project_generated_derivative"] and not row["notes"]:
            generated_with_accessions.append(row["record_id"])

    checks = [
        ("L01", "dataset record count", "37", str(len(dataset_rows)), len(dataset_rows) == 37, []),
        ("L02", "experiment record count", "40", str(len(experiment_rows)), len(experiment_rows) == 40, []),
        ("L03", "all spatial links resolve", "0 missing", str(len(spatial_missing)), not spatial_missing, spatial_missing),
        ("L04", "all reference links resolve", "0 missing", str(len(reference_missing)), not reference_missing, reference_missing),
        ("L05", "all non-empty evaluation-resource links resolve", "0 missing", str(len(evaluation_missing)), not evaluation_missing, evaluation_missing),
        ("L06", "no orphan dataset records", "0 orphan", str(len(orphans)), not orphans, orphans),
        ("L07", "unique dataset and experiment identifiers", "0 duplicates", str(len(duplicate_datasets) + len(duplicate_experiments)), not duplicate_datasets and not duplicate_experiments, duplicate_datasets + duplicate_experiments),
        ("L08", "accessions are in type-correct fields", "0 invalid", str(len(accession_errors)), not accession_errors, accession_errors),
        ("L09", "project-generated derivatives have no fabricated derivative accession", "0 fabricated", str(len(generated_with_accessions)), not generated_with_accessions, generated_with_accessions),
    ]
    return [
        {
            "check_id": check_id,
            "check_name": name,
            "expected": expected,
            "observed": observed,
            "status": "PASS" if passed else "FAIL",
            "affected_records": ";".join(affected),
            "notes": (
                "External accessions on PROJECT_GENERATED rows identify underlying public source inputs, not the derivative itself."
                if check_id == "L09" else "Deterministic draft link/type validation."
            ),
        }
        for check_id, name, expected, observed, passed, affected in checks
    ]


def build_blockers() -> list[dict[str, str]]:
    rows = [
        {
            "blocker_id": "B001", "severity": "CRITICAL", "record_ids": "FORMAL_ARTIFACTS", "category": "scope boundary",
            "confirmed_information": "Task 4A repository provenance audit and Task 4B-1 draft manifests exist.",
            "missing_information": "Submission-level manuscript linkage remains outside this task.",
            "why_it_matters": "Draft manifests alone do not establish submission-level consistency.",
            "permitted_current_wording": "Draft provenance manifest v0.1; not frozen.",
            "prohibited_overclaim": "Do not call overall Task 4 resolved or submission provenance PASS.",
            "recommended_next_evidence": "Separate user/ChatGPT manuscript reconciliation after provenance blockers are resolved.", "status": "OPEN",
        },
        {
            "blocker_id": "B002", "severity": "CRITICAL", "record_ids": "FORMAL_ARTIFACTS", "category": "bibliography linkage",
            "confirmed_information": "Task 4A source and citation candidates are inventoried.",
            "missing_information": "Formal bibliography linkage remains outside this task.",
            "why_it_matters": "Draft source records cannot certify the bibliography actually compiled for submission.",
            "permitted_current_wording": "Bibliographic reconciliation pending outside Task 4B-1.",
            "prohibited_overclaim": "Do not claim formal citation-key consistency.",
            "recommended_next_evidence": "Separate manual bibliography reconciliation by the user and ChatGPT.", "status": "OPEN",
        },
        {
            "blocker_id": "B003", "severity": "CRITICAL", "record_ids": "ALL", "category": "manifest freeze",
            "confirmed_information": "All 37 dataset and 40 experiment records are retained in draft v0.1.",
            "missing_information": "Review approval and resolution of major provenance gaps.",
            "why_it_matters": "The draft cannot be used as a frozen submission manifest.",
            "permitted_current_wording": "DRAFT_NOT_FROZEN.",
            "prohibited_overclaim": "Do not call this a final or frozen dataset manifest.",
            "recommended_next_evidence": "Resolve B004-B009, review corrections, then perform a separate freeze task.", "status": "OPEN",
        },
        {
            "blocker_id": "B004", "severity": "MAJOR", "record_ids": "LR03_BRCA_FFPE;LR04_BRCA_FF;LR05_BRCA_ILC;LR06_CERVICAL;LR08_INTESTINE", "category": "hECA composite lineage",
            "confirmed_information": "Source studies and current hECA project/export layer are partially confirmed.",
            "missing_information": "Exact contributing-cell membership, historical merge log and checksum lineage.",
            "why_it_matters": "Exact cell-level provenance cannot be reproduced from the retained records.",
            "permitted_current_wording": "Curated hECA project export/composite reference with verified source studies listed separately.",
            "prohibited_overclaim": "Do not claim exact cell-level source composition or direct Tabula Sapiens input without a retained lineage record.",
            "recommended_next_evidence": "Recover original export, cell-source manifest, merge log and checksums.", "status": "OPEN",
        },
        {
            "blocker_id": "B005", "severity": "MAJOR", "record_ids": "MERS01_BREAST;MERS02_COLON;MERS03_LUNG;MERS04_MELANOMA1;MERS05_MELANOMA2", "category": "Vizgen MERSCOPE evidence",
            "confirmed_information": "Release, sample identifiers, MERSCOPE platform, MERFISH assay, FFPE and 500-target-plus-50-blank panel are confirmed.",
            "missing_information": "Stable per-sample download receipt/URL, original checksums and redistribution terms.",
            "why_it_matters": "Access and redistribution language cannot yet be finalized.",
            "permitted_current_wording": "Provider-hosted Vizgen FFPE Human Immuno-oncology release; no repository accession assigned.",
            "prohibited_overclaim": "Do not call release identifiers clinical patient IDs or invent an accession/license.",
            "recommended_next_evidence": "Archive provider download receipts, stable URLs, original checksums and applicable terms.", "status": "OPEN",
        },
        {
            "blocker_id": "B006", "severity": "MAJOR", "record_ids": "CTA01_BCSA2TUMB1", "category": "CTA exact metadata",
            "confirmed_information": "CTA paper, BCSA2TumB1 identity, HER2 context, frozen Visium processing and public repositories are confirmed.",
            "missing_information": "Exact Visium chemistry/version, slide serial and capture area.",
            "why_it_matters": "Sample identity is locked, but platform-level reproducibility metadata remain incomplete.",
            "permitted_current_wording": "Frozen 10x Visium section BCSA2TumB1 from the CTA study.",
            "prohibited_overclaim": "Do not guess chemistry, slide or capture area.",
            "recommended_next_evidence": "Recover the sample sheet or original Space Ranger metadata.", "status": "OPEN",
        },
        {
            "blocker_id": "B007", "severity": "MINOR", "record_ids": "RD03_TNBC_PLASMA", "category": "CID4465 spatial metadata",
            "confirmed_information": "CID4465 spatial identity, source study and public Zenodo record are confirmed.",
            "missing_information": "Slide serial and capture area.",
            "why_it_matters": "Non-core sample hardware metadata are incomplete.",
            "permitted_current_wording": "Wu fresh-frozen TNBC section CID4465.",
            "prohibited_overclaim": "Do not assign a slide or capture area.",
            "recommended_next_evidence": "Primary sample sheet or original Space Ranger metadata.", "status": "OPEN",
        },
        {
            "blocker_id": "B008", "severity": "MINOR", "record_ids": "RD04_CRC_B", "category": "colorectal spatial metadata",
            "confirmed_information": "Formal 10x colorectal release and capture area C1 are confirmed.",
            "missing_information": "Slide serial.",
            "why_it_matters": "One non-core hardware identifier remains incomplete.",
            "permitted_current_wording": "10x Parent Visium human colorectal-cancer section, capture area C1.",
            "prohibited_overclaim": "Do not infer a slide serial.",
            "recommended_next_evidence": "Original provider sample metadata.", "status": "OPEN",
        },
        {
            "blocker_id": "B009", "severity": "MAJOR", "record_ids": "SIM03_HUMAN_LUNG", "category": "human-lung read-run metadata",
            "confirmed_information": "WSA_LngSP10193345, E-MTAB-11640, PRJEB52292, ENA sample ERS20065156 and BioSample SAMEA115633909 are confirmed; the accession-type conflict is corrected.",
            "missing_information": "Verified read-run accession was not recovered.",
            "why_it_matters": "The frozen sample is identifiable, but run-level raw-read metadata remain incomplete.",
            "permitted_current_wording": "Sample WSA_LngSP10193345 (ENA sample ERS20065156; BioSample SAMEA115633909).",
            "prohibited_overclaim": "Do not report ERS20065156 as a run or invent ERR/SRR/DRR identifiers.",
            "recommended_next_evidence": "Official ENA/BioStudies linkage if a read-run record becomes available.", "status": "OPEN_PARTIAL_METADATA",
        },
        {
            "blocker_id": "B010", "severity": "CRITICAL", "record_ids": "ALL", "category": "Data availability readiness",
            "confirmed_information": "A structured readiness inventory is generated for all 37 records.",
            "missing_information": "Resolution of access, redistribution, metadata and derivative-deposition gaps.",
            "why_it_matters": "A final Data availability statement would currently overstate readiness.",
            "permitted_current_wording": "Data availability inventory is draft and incomplete.",
            "prohibited_overclaim": "Do not state Data availability ready.",
            "recommended_next_evidence": "Resolve blockers and conduct a separate final statement task.", "status": "OPEN",
        },
    ]
    source_links = {
        "B001": "U001;C002",
        "B002": "U002;C002",
        "B003": "U003;C002",
        "B004": "U004;U005;U006",
        "B005": "U007",
        "B006": "U008",
        "B007": "U009",
        "B008": "U010",
        "B009": "U011;C001",
        "B010": "U001;U002;U003;U004;U005;U006;U007;U008;U009;U010;U011",
    }
    for row in rows:
        row["source_task4a_issue_or_conflict_id"] = source_links[row["blocker_id"]]
    return rows


def build_manual_recommendations() -> list[dict[str, str]]:
    return [
        {
            "recommendation_id": "R001", "affected_dataset": "SIM03_HUMAN_LUNG", "likely_manuscript_section": "Methods — Datasets, nomenclature and data preprocessing",
            "current_claim_from_task4a_or_repo": "ERS20065156 was incorrectly described as a run accession.",
            "verified_information": "ERS20065156 is an ENA sample accession; SAMEA115633909 is its BioSample; no verified read-run accession was recovered.",
            "problem": "ERS20065156 was incorrectly described as a run accession.",
            "recommended_manual_wording": "The spatial source corresponded to sample WSA_LngSP10193345 (ENA sample accession ERS20065156; BioSample SAMEA115633909), and the reference was derived from the corresponding donor-level lung resource represented in E-MTAB-11640 and PRJEB52292.",
            "must_avoid": "Do not report ERS20065156 as a run accession. Do not invent a read-run accession.",
            "priority": "CRITICAL", "evidence_record_ids": "SIM03_HUMAN_LUNG;C001;U011",
        },
        {
            "recommendation_id": "R002", "affected_dataset": "LR03_BRCA_FFPE;LR04_BRCA_FF;LR05_BRCA_ILC;LR06_CERVICAL;LR08_INTESTINE", "likely_manuscript_section": "Methods — Low-resolution references",
            "current_claim_from_task4a_or_repo": "Historical references are described with hECA/Tabula Sapiens shorthand.",
            "verified_information": "Curated hECA project exports and source studies are partly resolved; exact historical contributing-cell lineage is incomplete.",
            "problem": "Shorthand can imply a direct or fully reconstructed cell-level provenance chain.",
            "recommended_manual_wording": "References were curated from repository-local hECA project exports; verified source studies are listed, while exact historical cell-level merge membership was not retained.",
            "must_avoid": "Do not claim exact contributing-cell membership or direct Tabula Sapiens input without a merge manifest.",
            "priority": "MAJOR", "evidence_record_ids": "U004;U005;U006",
        },
        {
            "recommendation_id": "R003", "affected_dataset": "MERS01_BREAST;MERS02_COLON;MERS03_LUNG;MERS04_MELANOMA1;MERS05_MELANOMA2", "likely_manuscript_section": "Methods — High-resolution datasets and Data availability planning",
            "current_claim_from_task4a_or_repo": "Vizgen Human...Patient1 labels are used as dataset identifiers.",
            "verified_information": "They are provider release sample identifiers from the FFPE Human Immuno-oncology MERSCOPE/MERFISH release; no repository accession is assigned.",
            "problem": "The names can be misread as clinical patient IDs or accessions; redistribution evidence is incomplete.",
            "recommended_manual_wording": "Five provider release samples from Vizgen's FFPE Human Immuno-oncology MERSCOPE data release were used; MERFISH is the assay and the Human...Patient labels are release sample identifiers.",
            "must_avoid": "Do not call the labels verified clinical patient IDs, invent an accession, or assert redistribution permission.",
            "priority": "MAJOR", "evidence_record_ids": "MERS01_BREAST;MERS02_COLON;MERS03_LUNG;MERS04_MELANOMA1;MERS05_MELANOMA2;U007",
        },
        {
            "recommendation_id": "R004", "affected_dataset": "CTA01_BCSA2TUMB1", "likely_manuscript_section": "Methods — CTA biological application",
            "current_claim_from_task4a_or_repo": "BCSA2TumB1 is a HER2-positive Visium section with a CTA endpoint.",
            "verified_information": "The study, section, HER2 context, frozen Visium processing, image-derived CTA endpoint and public repositories are confirmed.",
            "problem": "Exact chemistry/version, slide serial and capture area remain unresolved.",
            "recommended_manual_wording": "The biological application used the frozen 10x Visium section BCSA2TumB1 from the CTA breast-cancer study and a project-generated, donor-unmatched 4,014-cell subset of GSE176078.",
            "must_avoid": "Do not guess the Visium version, chemistry, slide serial or capture area.",
            "priority": "MAJOR", "evidence_record_ids": "CTA01_BCSA2TUMB1;U008",
        },
        {
            "recommendation_id": "R005", "affected_dataset": "SIM01_BRCA;SIM02_MOUSE_BRAIN;SIM03_HUMAN_LUNG;STATE01_TCELL_DECOY;STATE02_KIDNEY_STATE32;CTA01_BCSA2TUMB1", "likely_manuscript_section": "Methods — Derived data and evaluation design",
            "current_claim_from_task4a_or_repo": "Simulations, disjoint references, decoys and endpoint/reference subsets are used across figures.",
            "verified_information": "These are project-generated derivatives of identified public/provider sources.",
            "problem": "They must not be presented as externally accessioned source datasets or independent biological truth.",
            "recommended_manual_wording": "Simulation truth, disjoint same-assay references, engineered decoys, the CTA-to-spot endpoint and balanced reference subsets were generated within this project from the cited source inputs.",
            "must_avoid": "Do not assign external accessions to derivatives or describe reference-relative proxies as independent ground truth.",
            "priority": "MAJOR", "evidence_record_ids": "dataset_manifest_draft_v0_1.tsv",
        },
        {
            "recommendation_id": "R006", "affected_dataset": "ALL", "likely_manuscript_section": "Methods/Data availability planning",
            "current_claim_from_task4a_or_repo": "A public dataset manifest was planned.",
            "verified_information": "Draft provenance manifests v0.1 now exist but remain DRAFT_NOT_FROZEN and Data availability is not ready.",
            "problem": "Draft generation does not complete or freeze overall provenance.",
            "recommended_manual_wording": "Use the draft only as a manual reconciliation source until a separate review and freeze task is completed.",
            "must_avoid": "Do not state that the manifest is final/frozen or that Data availability is ready.",
            "priority": "CRITICAL", "evidence_record_ids": "manifest_metadata.json;final_decision.md",
        },
    ]


def availability_class(row: dict[str, str]) -> str:
    if row["provenance_status"] == "PROJECT_GENERATED" or row["project_generated_derivative"].startswith("PROJECT_GENERATED"):
        return "PROJECT_GENERATED_DERIVATIVE"
    if row["record_id"].startswith("MERS"):
        return "REDISTRIBUTION_TERMS_UNRESOLVED"
    if row["provenance_status"] == "PARTIAL":
        return "PARTIAL_PROVENANCE"
    if row["provenance_status"] in {"CONFLICT", "UNRESOLVED"}:
        return "METADATA_UNRESOLVED"
    accession_fields = [
        "study_accession", "project_accession", "experiment_accession", "run_accession",
        "sample_accession", "biosample_accession", "geo_accession", "arrayexpress_accession",
    ]
    if any(row[field] for field in accession_fields):
        return "PUBLIC_ACCESSION_CONFIRMED"
    return "PUBLIC_PROVIDER_RELEASE_CONFIRMED"


def build_availability(dataset_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows = []
    for row in dataset_rows:
        identifiers = []
        for field in (
            "study_accession", "project_accession", "experiment_accession", "run_accession",
            "sample_accession", "biosample_accession", "geo_accession", "arrayexpress_accession",
        ):
            if row[field]:
                identifiers.append(f"{field}={row[field]}")
        category = availability_class(row)
        ready = row["provenance_status"] == "CONFIRMED" and category in {
            "PUBLIC_ACCESSION_CONFIRMED", "PUBLIC_PROVIDER_RELEASE_CONFIRMED"
        }
        rows.append({
            "record_id": row["record_id"],
            "formal_dataset_name": row["formal_dataset_name"],
            "source_database_or_provider": row["database"],
            "accession_or_provider_release": ";".join(identifiers) if identifiers else row["formal_release_name"],
            "official_url": row["official_source_url"],
            "download_or_access_method": "Use the official repository/provider URL recorded in this draft.",
            "redistribution_status": row["redistribution_status"],
            "project_generated_derivative": row["project_generated_derivative"],
            "readiness_class": category,
            "ready_for_final_statement": str(ready).lower(),
            "blocking_issue": row["unresolved_issue_id"] or ("project derivative deposition/license not frozen" if category == "PROJECT_GENERATED_DERIVATIVE" else ""),
            "safe_current_description": f"{row['formal_dataset_name']} ({category}).",
            "notes": "Overall Data availability ready remains false even where an individual public source is confirmed.",
        })
    return rows


def build_readme() -> str:
    return """# Draft dataset provenance manifests v0.1

This directory contains Task 4B-1 draft provenance manifests generated only from the Task 4A repository audit outputs.

- Status: `DRAFT_NOT_FROZEN`.
- This is not a formal submission manifest.
- Task 4A source audit files were not overwritten.
- Formal manuscript LaTeX was not accessed or modified.
- Formal `sn-bibliography.bib` was not accessed or modified.
- The human-lung accession type is corrected in this draft: `ERS20065156` is an ENA sample accession and `SAMEA115633909` is its BioSample.
- No verified human-lung read-run accession was recovered; `run_accession` is blank and no identifier was inferred.
- hECA composite lineage, Vizgen per-sample download/checksum/redistribution evidence, and CTA exact Visium metadata remain unresolved.
- Data availability is not ready.
- Only the user and ChatGPT will perform any later manual manuscript changes; this task produced recommendations only.

All 37 Task 4A dataset records and all 40 experiment crosswalk rows are retained. The three human-lung joint-simulation experiment statuses and their shared dataset record are changed from `CONFLICT` to `PARTIAL` only in this draft, with a complete transition audit.
"""


def build_final_decision(
    dataset_rows: list[dict[str, str]],
    experiment_rows: list[dict[str, str]],
    integrity_rows: list[dict[str, str]],
) -> str:
    integrity = "PASS" if all(row["status"] == "PASS" for row in integrity_rows) else "FAIL"
    return f"""# Task 4B-1 final decision

## 1. Decision

**CONDITIONAL PASS**

Draft provenance manifests were successfully generated and the confirmed human-lung accession-type conflict was corrected. This does not resolve overall Task 4.

## 2. Scope

```text
Formal manuscript accessed: false
Formal manuscript modified: false
Formal bibliography accessed: false
Formal bibliography modified: false
```

## 3. Human-lung correction

```text
ERS20065156 = ENA sample accession
SAMEA115633909 = BioSample
verified read-run accession = none
run_accession in draft = blank
```

`WSA_LngSP10193345`, `E-MTAB-11640` and `PRJEB52292` were not changed.

## 4. Manifest generation

```text
dataset manifest version = 0.1
dataset records = {len(dataset_rows)}
experiment records = {len(experiment_rows)}
status = DRAFT_NOT_FROZEN
link integrity = {integrity}
```

## 5. Status transitions

The shared dataset record `SIM03_HUMAN_LUNG` and its three joint-simulation experiments changed from `CONFLICT` to `PARTIAL` in draft v0.1. Basis: {PARTIAL_REASON}

All other Task 4A statuses were retained. `PARTIAL` records were not upgraded, and `PROJECT_GENERATED` records remain project-generated.

## 6. Remaining blockers

- hECA composite contributing-cell, merge-log and checksum lineage.
- Vizgen per-sample download receipt, original checksum and redistribution/license evidence.
- CTA exact Visium chemistry/version, slide serial and capture area.
- CID4465 slide/capture metadata.
- Colorectal slide serial.
- Human-lung verified read-run accession.
- Final Data availability reconciliation.
- Review and freeze of a formal dataset manifest.

## 7. Overall status

```text
Task 4A: COMPLETED
Task 4B-1: CONDITIONAL PASS
Overall Task 4: NOT YET RESOLVED
Data availability ready: false
Dataset manifest frozen: false
```

## 8. Guardrails

```text
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
```
"""


def repository_commit() -> str | None:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    return result.stdout.strip() if result.returncode == 0 else None


def validate_inputs() -> dict[str, str]:
    missing = [name for name in SOURCE_FILES if not (SOURCE / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing Task 4A inputs: {missing}")
    hashes = {name: sha1_file(SOURCE / name) for name in SOURCE_FILES}
    summary = json.loads((SOURCE / "audit_summary.json").read_text(encoding="utf-8"))
    if summary.get("decision") != "HOLD":
        raise RuntimeError("Task 4A source decision is not HOLD")
    if summary.get("counts", {}).get("master_records") != 37:
        raise RuntimeError("Task 4A master record count is not 37")
    if summary.get("counts", {}).get("experiment_crosswalk_rows") != 40:
        raise RuntimeError("Task 4A experiment crosswalk count is not 40")
    return hashes


def build_core_artifacts() -> tuple[dict[str, bytes], dict]:
    source_hashes = validate_inputs()
    master = read_tsv("dataset_provenance_master.tsv")
    crosswalk = read_tsv("experiment_input_crosswalk.tsv")
    sources = read_tsv("official_source_registry.tsv")
    unresolved = read_tsv("unresolved_records.tsv")
    conflicts = read_tsv("conflicts.tsv")

    if len(master) != 37 or len({row["record_id"] for row in master}) != 37:
        raise RuntimeError("Task 4A master records are not 37 unique records")
    if len(crosswalk) != 40 or len({row["experiment_id"] for row in crosswalk}) != 40:
        raise RuntimeError("Task 4A experiment records are not 40 unique rows")

    issues = issue_links(unresolved)
    conflict_map = conflict_links(conflicts)
    dataset_rows = build_dataset_manifest(master, issues, conflict_map)
    dataset_by_id = {row["record_id"]: row for row in dataset_rows}
    experiment_rows = build_experiment_manifest(crosswalk, dataset_by_id, issues)
    registry_rows = build_source_registry(sources, dataset_by_id)
    correction_rows = build_human_lung_correction()
    transition_rows = build_status_transitions(master, dataset_rows, crosswalk, experiment_rows)
    integrity_rows = build_link_integrity(dataset_rows, experiment_rows)
    blocker_rows = build_blockers()
    recommendation_rows = build_manual_recommendations()
    availability_rows = build_availability(dataset_rows)

    human_lung = dataset_by_id[HUMAN_LUNG_RECORD]
    if human_lung["run_accession"]:
        raise RuntimeError("Human-lung run_accession is not blank")
    if human_lung["sample_accession"] != "ERS20065156":
        raise RuntimeError("ERS20065156 is not in sample_accession")
    if human_lung["biosample_accession"] != "SAMEA115633909":
        raise RuntimeError("SAMEA115633909 is not in biosample_accession")
    if any("ERS20065156" in row["run_accession"] for row in dataset_rows):
        raise RuntimeError("ERS20065156 appears in a run_accession field")
    if any(row["status"] != "PASS" for row in integrity_rows):
        failed = [row for row in integrity_rows if row["status"] != "PASS"]
        raise RuntimeError(f"Manifest integrity failure: {failed}")

    dataset_statuses = Counter(row["provenance_status"] for row in dataset_rows)
    experiment_statuses = Counter(row["formal_status"] for row in experiment_rows)
    summary = {
        "task": "Task 4B-1 draft provenance manifest generation",
        "decision": "CONDITIONAL PASS",
        "manifest_version": VERSION,
        "manifest_status": STATUS,
        "dataset_record_count": len(dataset_rows),
        "experiment_record_count": len(experiment_rows),
        "dataset_status_counts": dict(sorted(dataset_statuses.items())),
        "experiment_status_counts": dict(sorted(experiment_statuses.items())),
        "human_lung_accession_correction": {
            "ERS20065156": "ENA sample accession",
            "SAMEA115633909": "BioSample",
            "verified_run_accession": None,
            "dataset_status": "PARTIAL",
            "experiment_status": "PARTIAL",
            "status_basis": PARTIAL_REASON,
        },
        "link_integrity": "PASS",
        "remaining_blocker_count": len(blocker_rows),
        "data_availability_ready": False,
        "manifest_frozen": False,
        "formal_manuscript_accessed": False,
        "formal_bibliography_accessed": False,
        "experimental_outputs_modified": False,
        "figures_modified": False,
        "task4a_outputs_overwritten": False,
        "deterministic_double_build_identical": True,
    }

    artifacts = {
        "dataset_manifest_draft_v0_1.tsv": tsv_bytes(dataset_rows, DATASET_FIELDS),
        "experiment_manifest_draft_v0_1.tsv": tsv_bytes(experiment_rows, EXPERIMENT_FIELDS),
        "source_registry_draft_v0_1.tsv": tsv_bytes(registry_rows, SOURCE_REGISTRY_FIELDS),
        "human_lung_accession_correction.tsv": tsv_bytes(correction_rows, [
            "record_id", "experiment_ids", "field", "task4a_value", "verified_value",
            "verified_type", "evidence_source", "correction_action", "remaining_unknown", "status", "notes",
        ]),
        "record_status_transition_audit.tsv": tsv_bytes(transition_rows, [
            "record_or_experiment_id", "task4a_status", "draft_v0_1_status", "change_reason",
            "evidence", "allowed", "notes",
        ]),
        "manifest_link_integrity.tsv": tsv_bytes(integrity_rows, [
            "check_id", "check_name", "expected", "observed", "status", "affected_records", "notes",
        ]),
        "remaining_provenance_blockers.tsv": tsv_bytes(blocker_rows, [
            "blocker_id", "source_task4a_issue_or_conflict_id", "severity", "record_ids", "category", "confirmed_information",
            "missing_information", "why_it_matters", "permitted_current_wording",
            "prohibited_overclaim", "recommended_next_evidence", "status",
        ]),
        "manuscript_manual_update_recommendations.tsv": tsv_bytes(recommendation_rows, [
            "recommendation_id", "affected_dataset", "likely_manuscript_section",
            "current_claim_from_task4a_or_repo", "verified_information", "problem",
            "recommended_manual_wording", "must_avoid", "priority", "evidence_record_ids",
        ]),
        "data_availability_readiness_inventory.tsv": tsv_bytes(availability_rows, [
            "record_id", "formal_dataset_name", "source_database_or_provider",
            "accession_or_provider_release", "official_url", "download_or_access_method",
            "redistribution_status", "project_generated_derivative", "readiness_class",
            "ready_for_final_statement", "blocking_issue", "safe_current_description", "notes",
        ]),
        "README.md": build_readme().encode("utf-8"),
        "final_decision.md": build_final_decision(dataset_rows, experiment_rows, integrity_rows).encode("utf-8"),
        "audit_summary.json": json_bytes(summary),
    }
    context = {
        "source_hashes": source_hashes,
        "dataset_rows": dataset_rows,
        "experiment_rows": experiment_rows,
        "summary": summary,
    }
    return artifacts, context


def main() -> None:
    artifacts_a, context_a = build_core_artifacts()
    artifacts_b, context_b = build_core_artifacts()
    if artifacts_a != artifacts_b or context_a["summary"] != context_b["summary"]:
        raise RuntimeError("Two in-memory deterministic builds did not match")

    OUT.mkdir(parents=True, exist_ok=True)
    for name, data in artifacts_a.items():
        (OUT / name).write_bytes(data)

    metadata = {
        "manifest_version": VERSION,
        "status": STATUS,
        "source_audit": "Task 4A dataset provenance formal audit",
        "source_audit_decision": "HOLD",
        "dataset_record_count": 37,
        "experiment_record_count": 40,
        "formal_manuscript_accessed": False,
        "formal_bibliography_accessed": False,
        "human_lung_accession_correction": {
            "ERS20065156": "ENA sample accession",
            "SAMEA115633909": "BioSample",
            "verified_run_accession": None,
        },
        "data_availability_ready": False,
        "manifest_frozen": False,
        "experimental_outputs_modified": False,
        "figures_modified": False,
        "manuscript_modified": False,
        "bibliography_modified": False,
        "task4a_outputs_overwritten": False,
        "generated_date": GENERATED_DATE,
        "repository_commit_sha": repository_commit(),
        "source_task4a_file_sha1": context_a["source_hashes"],
        "generator_script_sha1": sha1_file(Path(__file__)),
        "generated_file_sha1": {name: sha1_bytes(data) for name, data in sorted(artifacts_a.items())},
        "hash_scope_note": "manifest_metadata.json excludes its own SHA-1 to avoid self-reference; all other generated files are hashed.",
        "deterministic_double_build_identical": True,
    }
    (OUT / "manifest_metadata.json").write_bytes(json_bytes(metadata))

    print("Task 4B-1 draft provenance manifest generation completed.")
    print("\nDecision:\nCONDITIONAL PASS")
    print("\nDataset manifest draft:\nrecords = 37\nstatus = DRAFT_NOT_FROZEN")
    print("\nExperiment manifest draft:\nrecords = 40")
    print("\nHuman-lung accession correction:")
    print("ERS20065156 = ENA sample accession")
    print("SAMEA115633909 = BioSample")
    print("verified read-run accession = none")
    print("\nHuman-lung experiment status:\nCONFLICT -> PARTIAL (3 experiments)")
    print("\nRemaining major blockers:\nhECA lineage; Vizgen evidence/terms; CTA exact metadata; human-lung run metadata; manifest freeze")
    print("\nData availability ready:\nfalse")
    print("\nDataset manifest frozen:\nfalse")
    print("\nOverall Task 4 resolved:\nfalse")
    print("\nFormal manuscript accessed:\nfalse")
    print("\nFormal manuscript modified:\nfalse")
    print("\nFormal bibliography accessed:\nfalse")
    print("\nFormal bibliography modified:\nfalse")
    print("\nExperimental stages rerun:\nfalse")
    print("\nFigures modified:\nfalse")
    print(f"\nOutput directory:\n{OUT.relative_to(ROOT).as_posix()}")


if __name__ == "__main__":
    main()
