#!/usr/bin/env python3
"""Build Task 4B-2 provenance evidence recovery artifacts.

The script is deliberately audit-only. It reads the Task 4B-1 draft and
retained local provenance evidence, then writes a new v0.2 draft. It does not
read manuscript LaTeX/BibTeX or modify experimental inputs, outputs, or figures.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import subprocess
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

import h5py


ROOT = Path(__file__).resolve().parents[1]
TASK4A_DIR = ROOT / "visualizations/manuscript_audits/dataset_provenance_formal_audit"
TASK4B1_DIR = ROOT / "visualizations/manuscript_audits/dataset_provenance_manifest_draft_v0_1"
OUT_DIR = ROOT / "visualizations/manuscript_audits/dataset_provenance_major_evidence_recovery"
GENERATED_DATE = "2026-07-22"
_SHA1_CACHE: dict[tuple[str, int, int], str] = {}

DATASET_V01 = TASK4B1_DIR / "dataset_manifest_draft_v0_1.tsv"
EXPERIMENT_V01 = TASK4B1_DIR / "experiment_manifest_draft_v0_1.tsv"
SOURCE_V01 = TASK4B1_DIR / "source_registry_draft_v0_1.tsv"

CTA_ID = "CTA01_BCSA2TUMB1"
VIZGEN_IDS = [
    "MERS01_BREAST",
    "MERS02_COLON",
    "MERS03_LUNG",
    "MERS04_MELANOMA1",
    "MERS05_MELANOMA2",
]
HECA_IDS = ["LR03_BRCA_FFPE", "LR04_BRCA_FF", "LR05_BRCA_ILC", "LR06_CERVICAL", "LR08_INTESTINE"]
MINOR_IDS = ["RD03_TNBC_PLASMA", "RD04_CRC_B"]

VIZGEN_SHOWCASE = "https://info.vizgen.com/ffpe-showcase"
VIZGEN_ROADMAP = "https://vizgen.com/human-ffpe-immunooncology-release-roadmap/"
VIZGEN_BUCKET = "https://console.cloud.google.com/storage/browser/vz-ffpe-showcase"
CTA_ZENODO = "https://zenodo.org/records/15211538"
CTA_RESEARCHDATA = "https://researchdata.se/en/catalogue/dataset/2025-97"
CTA_PAPER = "https://doi.org/10.1038/s41698-025-01104-3"
CRC_WEB_SUMMARY = (
    "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/"
    "Parent_Visium_Human_ColorectalCancer/Parent_Visium_Human_ColorectalCancer_web_summary.html"
)


def rel(path: Path) -> str:
    return path.resolve().relative_to(ROOT.resolve()).as_posix()


def sha1_file(path: Path) -> str:
    stat = path.stat()
    cache_key = (str(path.resolve()), stat.st_size, stat.st_mtime_ns)
    cached = _SHA1_CACHE.get(cache_key)
    if cached is not None:
        return cached
    digest = hashlib.sha1()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    value = digest.hexdigest()
    _SHA1_CACHE[cache_key] = value
    return value


def sha1_bytes(data: bytes) -> str:
    return hashlib.sha1(data).hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def tsv_bytes(rows: Iterable[dict[str, Any]], fields: list[str]) -> bytes:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        writer.writerow({field: "" if row.get(field) is None else row.get(field, "") for field in fields})
    return output.getvalue().encode("utf-8")


def json_bytes(value: Any) -> bytes:
    return (json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n").encode("utf-8")


def text_bytes(value: str) -> bytes:
    return value.rstrip().encode("utf-8") + b"\n"


def append_evidence(existing: str, *items: str) -> str:
    values = [part.strip() for part in existing.split(";") if part.strip()]
    for item in items:
        if item and item not in values:
            values.append(item)
    return "; ".join(values)


def git_commit() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, capture_output=True, text=True
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return ""


def file_metadata(path: Path) -> dict[str, Any]:
    stat = path.stat()
    return {
        "file_path": rel(path),
        "file_name": path.name,
        "file_type": path.suffix.lower() or "directory-entry",
        "file_size": stat.st_size,
        "modified_time": stat.st_mtime_ns,
        "sha1": sha1_file(path),
    }


def decode_h5ad_column(group: h5py.Group, name: str) -> list[str]:
    node = group[name]
    if isinstance(node, h5py.Dataset):
        values = node[:]
        return [v.decode("utf-8", "replace") if isinstance(v, bytes) else str(v) for v in values]
    codes = node["codes"][:]
    categories = node["categories"][:]
    decoded = [v.decode("utf-8", "replace") if isinstance(v, bytes) else str(v) for v in categories]
    return [decoded[int(code)] if int(code) >= 0 else "NA" for code in codes]


def h5ad_summary(path: Path) -> dict[str, Any]:
    with h5py.File(path, "r") as handle:
        obs = handle["obs"]
        n_cells = int(obs["_index"].shape[0])
        var_index = handle["var"][handle["var"].attrs.get("_index", "_index")]
        studies = Counter(decode_h5ad_column(obs, "study_id"))
        donors = Counter(decode_h5ad_column(obs, "donor_ID"))
        organs = Counter(decode_h5ad_column(obs, "organ"))
        return {
            "n_cells": n_cells,
            "n_genes": int(var_index.shape[0]),
            "study_counts": dict(sorted(studies.items())),
            "donor_count": len(donors),
            "donor_counts": dict(sorted(donors.items())),
            "organ_counts": dict(sorted(organs.items())),
            "obs_columns": sorted(obs.keys()),
        }


def count_csv_rows(path: Path) -> int:
    lines = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            lines += chunk.count(b"\n")
    return max(0, lines - 1)


def build_evidence_inventory() -> tuple[list[dict[str, Any]], dict[str, str]]:
    evidence_specs: list[tuple[str, str, str, str, int, str]] = []

    def add(record_id: str, category: str, path: str, supports: str, priority: int, notes: str = "") -> None:
        evidence_specs.append((record_id, category, path, supports, priority, notes))

    for filename in [
        "dataset_manifest_draft_v0_1.tsv",
        "experiment_manifest_draft_v0_1.tsv",
        "source_registry_draft_v0_1.tsv",
        "remaining_provenance_blockers.tsv",
        "manifest_metadata.json",
        "audit_summary.json",
    ]:
        add("AUDIT", "upstream audit", rel(TASK4B1_DIR / filename), "v0.1 baseline and protected source", 9)

    cta_files = [
        "data/raw/新建文件夹/BreastCancer_CTA-main/ST_data.RData",
        "data/raw/新建文件夹/BreastCancer_CTA-main/CTA_output.txt",
        "data/raw/新建文件夹/BreastCancer_CTA-main/README.md",
        "data/raw/新建文件夹/BreastCancer_CTA-main/CTA_align.R",
        "data/raw/新建文件夹/s41698-025-01104-3.pdf",
        "data/raw/新建文件夹/2025-97-1/documentation/README.txt",
        "data/raw/新建文件夹/2025-97-1/documentation/File_list.csv",
        "visualizations/bioapp_experiment/bioapp_coordinate_registration_recovery_audit/bioapp_st_data_rdata_object_inventory.json",
        "visualizations/bioapp_experiment/bioapp_coordinate_registration_recovery_audit/bioapp_registration_recovery_summary.json",
        "visualizations/bioapp_experiment/bioapp_figure_v2_6_frozen_spot_coordinate_provenance_audit/bioapp_v2_6_first_coordinate_source_report.md",
    ]
    for path in cta_files:
        add(CTA_ID, "CTA local evidence", path, "sample identity, processed spatial object, platform, alignment, endpoint lineage", 1 if path.endswith("ST_data.RData") else 5)

    vizgen = {
        "MERS01_BREAST": ("HumanBreastCancerPatient1", "highres_humanbreastcancerpatient1_profile_mask_monocytes_and_macrophages"),
        "MERS02_COLON": ("HumanColonCancerPatient1", "highres_humancoloncancerpatient1_profile_mask_fibroblasts"),
        "MERS03_LUNG": ("HumanLungCancerPatient1", "highres_humanlungcancerpatient1_profile_mask_plasma_cells"),
        "MERS04_MELANOMA1": ("HumanMelanomaPatient1", "highres_humanmelanomapatient1_profile_mask_fibroblasts"),
        "MERS05_MELANOMA2": ("HumanMelanomaPatient2", "highres_humanmelanomapatient2_profile_mask_b_cells"),
    }
    for record_id, (sample, scenario) in vizgen.items():
        for suffix in ["cell_by_gene.csv", "cell_metadata.csv"]:
            add(record_id, "Vizgen retained raw input", f"data/raw/high/{sample}/{sample}_{suffix}", "per-sample input identity and local checksum", 1)
        for suffix in ["sc_metadata.csv", "st_metadata.csv"]:
            add(record_id, "project-generated split", f"data/processed/{scenario}/stage1_preprocess/exported/{suffix}", "disjoint same-assay split membership", 2)
        add(record_id, "project-generated split", f"data/processed/{scenario}/stage1_preprocess/fig2d_profile_mask_info.json", "split and mask parameters", 4)
    add("MERSCOPE_ALL", "generation script", "scripts/build_highres_profile_mask_mapping.py", "split algorithm and seed", 5)
    add("MERSCOPE_ALL", "frozen manifest", "visualizations/highres_profile_mask_mapping/highres_profile_mask_mapping_manifest.json", "sample-to-scenario mapping and zero-overlap checks", 9)

    heca_paths = {
        "LR03_BRCA_FFPE": "data/raw/low_resolution/Human Breast Cancer/sc_human_breast_heca.h5ad",
        "LR04_BRCA_FF": "data/raw/low_resolution/Human Breast Cancer Visium FF WTA/sc_human_breast_heca.h5ad",
        "LR05_BRCA_ILC": "data/raw/low_resolution/Human Breast Cancer WTA 1.2.0/sc_human_breast_heca.h5ad",
        "LR06_CERVICAL": "data/raw/low_resolution/Human Cervical Cancer/sc_human_uterus_heca.h5ad",
        "LR08_INTESTINE": "data/raw/low_resolution/Human Intestine Cancer/sc_human_intestine_heca.h5ad",
    }
    for record_id, path in heca_paths.items():
        add(record_id, "hECA actual input", path, "cell-level study, donor, tissue, cell count, and local checksum", 1)

    minor_files = {
        "RD03_TNBC_PLASMA": [
            "data/raw/cytospace_fig2d_tme/brca_ST_Wu_TNBC_Zenodo/filtered_count_matrices.tar.gz",
            "data/raw/cytospace_fig2d_tme/brca_ST_Wu_TNBC_Zenodo/metadata.tar.gz",
            "data/raw/cytospace_fig2d_tme/brca_ST_Wu_TNBC_Zenodo/spatial.tar.gz",
        ],
        "RD04_CRC_B": [
            "data/raw/cytospace_fig2d_tme/crc_ST_10x_parent_visium/Parent_Visium_Human_ColorectalCancer_filtered_feature_bc_matrix.h5",
            "data/raw/cytospace_fig2d_tme/crc_ST_10x_parent_visium/Parent_Visium_Human_ColorectalCancer_spatial.tar.gz",
        ],
    }
    for record_id, paths in minor_files.items():
        for path in paths:
            add(record_id, "retained public source archive", path, "spatial sample identity and retained source checksum", 1)
    add("RD03_TNBC_PLASMA;RD04_CRC_B", "preprocessing script", "scripts/prepare_cytospace_fig2d_tme_stage1.py", "archive member selection and project preparation", 5)
    add("AUDIT", "generator", rel(Path(__file__)), "Task 4B-2 deterministic generation", 5)

    rows: list[dict[str, Any]] = []
    path_to_id: dict[str, str] = {}
    for index, (record_id, category, path_text, supports, priority, notes) in enumerate(evidence_specs, start=1):
        path = ROOT / path_text
        if not path.is_file():
            raise FileNotFoundError(path)
        evidence_id = f"E{index:03d}"
        metadata = file_metadata(path)
        rows.append({
            "evidence_id": evidence_id,
            "record_id": record_id,
            "category": category,
            **metadata,
            "evidence_priority": priority,
            "supports_fields": supports,
            "read_only": "true",
            "notes": notes,
        })
        path_to_id[path_text] = evidence_id
    return rows, path_to_id


def heca_audits(evidence_ids: dict[str, str]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, dict[str, Any]]]:
    specs = {
        "LR03_BRCA_FFPE": {
            "path": "data/raw/low_resolution/Human Breast Cancer/sc_human_breast_heca.h5ad",
            "formal_export": "hECA 2.0 project export RNA-10.1101_2021.07.19.452956-Breast.h5ad.zip",
            "official_url": "https://zenodo.org/records/17008296",
            "official_archive_md5": "df6999cc9c5c7800498d295a9a3a9401",
        },
        "LR04_BRCA_FF": {
            "path": "data/raw/low_resolution/Human Breast Cancer Visium FF WTA/sc_human_breast_heca.h5ad",
            "formal_export": "hECA 2.0 project export RNA-10.1101_2021.07.19.452956-Breast.h5ad.zip",
            "official_url": "https://zenodo.org/records/17008296",
            "official_archive_md5": "df6999cc9c5c7800498d295a9a3a9401",
        },
        "LR05_BRCA_ILC": {
            "path": "data/raw/low_resolution/Human Breast Cancer WTA 1.2.0/sc_human_breast_heca.h5ad",
            "formal_export": "hECA 2.0 project export RNA-10.1101_2021.07.19.452956-Breast.h5ad.zip",
            "official_url": "https://zenodo.org/records/17008296",
            "official_archive_md5": "df6999cc9c5c7800498d295a9a3a9401",
        },
        "LR06_CERVICAL": {
            "path": "data/raw/low_resolution/Human Cervical Cancer/sc_human_uterus_heca.h5ad",
            "formal_export": "hECA 2.0 organ export RNA-Uterus.h5ad.zip",
            "official_url": "https://zenodo.org/records/15620903",
            "official_archive_md5": "26a3067c177094e8627624d25275a121",
        },
        "LR08_INTESTINE": {
            "path": "data/raw/low_resolution/Human Intestine Cancer/sc_human_intestine_heca.h5ad",
            "formal_export": "hECA 2.0 organ export RNA-Intestine.h5ad.zip",
            "official_url": "https://zenodo.org/records/15619143",
            "official_archive_md5": "b8604d168d24d9d64862cddaa3bbad0c",
        },
    }
    lineage_rows: list[dict[str, Any]] = []
    input_rows: list[dict[str, Any]] = []
    summaries: dict[str, dict[str, Any]] = {}
    for record_id, spec in specs.items():
        path = ROOT / spec["path"]
        summary = h5ad_summary(path)
        summary["sha1"] = sha1_file(path)
        summaries[record_id] = summary
        studies = summary["study_counts"]
        lineage_status = "CONFIRMED_EXACT_LINEAGE"
        safe_claim = (
            "The retained hECA input is identified by local SHA-1 and contains cell-level study_id and donor_ID fields; "
            f"its exact retained composition is {json.dumps(studies, sort_keys=True)}."
        )
        lineage_rows.append({
            "record_id": record_id,
            "actual_input_file": spec["path"],
            "actual_input_confirmed": "true",
            "formal_heca_export_identity": spec["formal_export"],
            "single_export_or_project_merge": "single retained hECA export; no project-side multi-file merge",
            "project_merge_script": "NOT_APPLICABLE",
            "cell_count": summary["n_cells"],
            "gene_count": summary["n_genes"],
            "contributing_study_counts": json.dumps(studies, sort_keys=True),
            "contributing_study_count": len(studies),
            "cell_level_source_mapping": "CONFIRMED: obs.study_id",
            "donor_level_mapping": f"CONFIRMED: obs.donor_ID ({summary['donor_count']} unique donors)",
            "accession_field_present": "false; study DOI values are present in obs.study_id",
            "input_sha1": summary["sha1"],
            "official_export_url": spec["official_url"],
            "official_archive_md5": spec["official_archive_md5"],
            "historical_archive_byte_match": "not tested because the original compressed download was not retained",
            "lineage_status": lineage_status,
            "supporting_evidence_ids": evidence_ids[spec["path"]],
            "safe_claim_level": safe_claim,
            "notes": "Exact lineage refers to the retained input's embedded cell-level provenance and checksum; it does not assert a byte-for-byte match to the current compressed Zenodo archive.",
        })
        stat = path.stat()
        input_rows.append({
            "record_id": record_id,
            "file_path": spec["path"],
            "file_name": path.name,
            "file_size": stat.st_size,
            "modified_time": stat.st_mtime_ns,
            "sha1": summary["sha1"],
            "n_cells": summary["n_cells"],
            "n_genes": summary["n_genes"],
            "obs_columns": ";".join(summary["obs_columns"]),
            "study_counts": json.dumps(studies, sort_keys=True),
            "donor_count": summary["donor_count"],
            "organ_counts": json.dumps(summary["organ_counts"], sort_keys=True),
            "evidence_id": evidence_ids[spec["path"]],
            "read_only": "true",
        })
    return lineage_rows, input_rows, summaries


def cta_recovery(evidence_ids: dict[str, str]) -> list[dict[str, Any]]:
    st_path = "data/raw/新建文件夹/BreastCancer_CTA-main/ST_data.RData"
    paper_path = "data/raw/新建文件夹/s41698-025-01104-3.pdf"
    readme_path = "data/raw/新建文件夹/BreastCancer_CTA-main/README.md"
    recovered = [
        ("formal sample ID", "BCSA2TumB1", "BCSA2TumB1", "CONFIRMED", st_path, 1, "Seurat image slot and orig.ident both use BCSA2TumB1."),
        ("provider patient ID", "BCSA2", "BCSA2", "CONFIRMED", readme_path, 5, "Official example identifies patient 2 as BCSA2."),
        ("tumor area", "TumB1", "TumB1", "CONFIRMED", readme_path, 5, "Official sample key is BCSA2TumB1."),
        ("spatial platform", "10x Visium Spatial Gene Expression", "10x Visium Spatial Gene Expression", "CONFIRMED", paper_path, 14, "Original paper methods."),
        ("assay version", "", "Visium Spatial Gene Expression kit; exact kit/chemistry version not stated", "PARTIAL", paper_path, 14, "Do not infer a chemistry version from Space Ranger version."),
        ("chemistry", "", "", "UNRESOLVED", paper_path, 14, "No exact chemistry label was found in the retained object, official metadata workbook, paper, or repository documentation."),
        ("tissue processing", "fresh frozen", "OCT-embedded fresh-frozen tissue, 10 micrometre section", "CONFIRMED", paper_path, 14, "Original paper methods state OCT embedding and 10 micrometre sections."),
        ("slide serial", "UNRESOLVED", "V10F24-112", "CONFIRMED", st_path, 12, "Official Zenodo spaceranger_output.zip/Visium_metadata.xlsx maps BCSA2 TumB1 to V10F24-112_C1."),
        ("capture area", "UNRESOLVED", "C1", "CONFIRMED", st_path, 12, "Recovered from the same official Visium_metadata.xlsx row; not inferred from TumB1."),
        ("capture format", "", "6.5 mm x 6.5 mm Visium capture area", "CONFIRMED", paper_path, 14, "Original paper methods; this is the format, not the capture-area code."),
        ("Space Ranger version", "", "1.0.0", "CONFIRMED", paper_path, 14, "Original paper reports Space Ranger 1.0.0."),
        ("library ID", "", "V10F24-112_C1", "CONFIRMED", st_path, 12, "Official Visium metadata identifier combining slide serial and capture area."),
        ("raw source path", "", "Official Zenodo spaceranger_output.zip/V10F24-112_C1/outs", "CONFIRMED", st_path, 12, "Official archive member listing recovered by HTTP range; the full archive was not copied locally."),
        ("processed source path", st_path, st_path, "CONFIRMED", st_path, 1, "Actual frozen spatial input object."),
        ("official repository", "Swedish National Data Service; Zenodo", "Swedish National Data Service 2025-97; Zenodo 15211538", "CONFIRMED", paper_path, 12, "Official raw and processed data records."),
        ("accession", "10.48723/f4v5-m008; Researchdata.se 2025-97", "10.48723/f4v5-m008; 10.5281/zenodo.15211538", "CONFIRMED", paper_path, 12, "Raw source is restricted; processed Space Ranger outputs are provider-hosted on Zenodo."),
        ("provider URL", CTA_ZENODO, CTA_ZENODO, "CONFIRMED", paper_path, 12, "Official processed-data archive."),
        ("original CTA annotation", "", "CTA_output.txt plus CTA_align.R", "CONFIRMED", "data/raw/新建文件夹/BreastCancer_CTA-main/CTA_output.txt", 1, "Original cell-level computational pathology output retained locally."),
        ("project spatial registration", "", "Processed Seurat image-slot registration in ST_data.RData", "CONFIRMED", st_path, 1, "Project registration is not relabelled as original CTA annotation."),
        ("project frozen endpoint", "", "Project-generated frozen spot-level immune endpoint", "PROJECT_GENERATED", "visualizations/bioapp_experiment/bioapp_coordinate_registration_recovery_audit/bioapp_registration_recovery_summary.json", 10, "Separate derivative layer; not an original CTA annotation field."),
    ]
    rows = []
    for field, old, value, status, local_path, priority, notes in recovered:
        rows.append({
            "record_id": CTA_ID,
            "field": field,
            "task4b1_value": old,
            "recovered_value": value,
            "status": status,
            "local_evidence_path": local_path,
            "local_evidence_sha1": sha1_file(ROOT / local_path),
            "official_evidence": f"{CTA_PAPER}; {CTA_ZENODO}; {CTA_RESEARCHDATA}",
            "evidence_priority": priority,
            "notes": notes,
        })
    return rows


def vizgen_audits(evidence_ids: dict[str, str]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, dict[str, Any]]]:
    specs = {
        "MERS01_BREAST": ("HumanBreastCancerPatient1", "Breast cancer", 713121, 490398542, 2.7, 87.37, "highres_humanbreastcancerpatient1_profile_mask_monocytes_and_macrophages"),
        "MERS02_COLON": ("HumanColonCancerPatient1", "Colon cancer 1", 677451, 411716053, 2.4, 72.56, "highres_humancoloncancerpatient1_profile_mask_fibroblasts"),
        "MERS03_LUNG": ("HumanLungCancerPatient1", "Lung cancer 1", 353762, 144388044, 4.3, 66.97, "highres_humanlungcancerpatient1_profile_mask_plasma_cells"),
        "MERS04_MELANOMA1": ("HumanMelanomaPatient1", "Melanoma 1", 468138, 160181929, 1.9, 60.73, "highres_humanmelanomapatient1_profile_mask_fibroblasts"),
        "MERS05_MELANOMA2": ("HumanMelanomaPatient2", "Melanoma 2", 207869, 75617432, 1.8, 65.64, "highres_humanmelanomapatient2_profile_mask_b_cells"),
    }
    release_rows: list[dict[str, Any]] = []
    inventory_rows: list[dict[str, Any]] = []
    summaries: dict[str, dict[str, Any]] = {}
    for record_id, (sample, provider_name, expected_cells, transcripts, rin, dv200, scenario) in specs.items():
        raw_dir = ROOT / f"data/raw/high/{sample}"
        expr = raw_dir / f"{sample}_cell_by_gene.csv"
        meta = raw_dir / f"{sample}_cell_metadata.csv"
        observed_cells = count_csv_rows(meta)
        if observed_cells != expected_cells:
            raise ValueError(f"{sample}: local metadata rows {observed_cells} != official total {expected_cells}")
        summaries[record_id] = {"sample": sample, "cells": observed_cells, "scenario": scenario}
        direct_url = f"{VIZGEN_BUCKET}/{sample}"
        release_rows.append({
            "record_id": record_id,
            "sample_identifier": sample,
            "provider_sample_name": provider_name,
            "provider": "Vizgen",
            "platform": "MERSCOPE",
            "assay": "MERFISH",
            "processing": "FFPE",
            "release": "MERSCOPE FFPE Human Immuno-Oncology Data Release, May 2022",
            "panel": "500 target genes plus 50 blank-control channels retained in local matrix",
            "official_cell_total": expected_cells,
            "local_metadata_row_count": observed_cells,
            "official_transcript_total": transcripts,
            "rin": rin,
            "dv200": dv200,
            "provider_release_page": VIZGEN_SHOWCASE,
            "per_sample_download_url": direct_url,
            "download_url_preserved": "true: recovered from the current official release page",
            "historical_download_receipt_preserved": "false",
            "download_date_preserved": "false",
            "provider_checksum_preserved": "false",
            "license_or_terms_evidence": f"Legal notice on {VIZGEN_SHOWCASE}; roadmap states release participants may use the data in any way.",
            "redistribution_status": "PROVIDER_TERMS_FOUND_BUT_AMBIGUOUS",
            "same_assay_split_status": "PROJECT_GENERATED",
            "provenance_status": "PARTIAL",
            "notes": "Exact sample identity is supported by the provider folder URL and the equality of local metadata rows to the official per-sample cell total. Redistribution is not inferred from general use language.",
        })
        for file_kind, path in [("raw cell-by-gene matrix", expr), ("raw cell metadata", meta)]:
            inventory_rows.append({
                "record_id": record_id,
                "sample_identifier": sample,
                "file_role": file_kind,
                "original_filename": path.name,
                "local_raw_directory": rel(raw_dir),
                "local_processed_directory": f"data/processed/{scenario}",
                "file_path": rel(path),
                "file_size": path.stat().st_size,
                "file_sha1": sha1_file(path),
                "provider_checksum": "",
                "provider_checksum_status": "NOT_PRESERVED",
                "download_url": direct_url,
                "download_date": "",
                "evidence_id": evidence_ids[rel(path)],
                "redistribution_status": "PROVIDER_TERMS_FOUND_BUT_AMBIGUOUS",
                "notes": "Retained local file; no original archive or provider checksum record was found.",
            })
        split_paths = [
            f"data/processed/{scenario}/stage1_preprocess/exported/sc_metadata.csv",
            f"data/processed/{scenario}/stage1_preprocess/exported/st_metadata.csv",
            f"data/processed/{scenario}/stage1_preprocess/fig2d_profile_mask_info.json",
        ]
        for path_text in split_paths:
            path = ROOT / path_text
            inventory_rows.append({
                "record_id": record_id,
                "sample_identifier": sample,
                "file_role": "project-generated disjoint split metadata",
                "original_filename": "NOT_APPLICABLE",
                "local_raw_directory": rel(raw_dir),
                "local_processed_directory": f"data/processed/{scenario}",
                "file_path": path_text,
                "file_size": path.stat().st_size,
                "file_sha1": sha1_file(path),
                "provider_checksum": "NOT_APPLICABLE",
                "provider_checksum_status": "PROJECT_GENERATED",
                "download_url": "NOT_APPLICABLE",
                "download_date": "NOT_APPLICABLE",
                "evidence_id": evidence_ids[path_text],
                "redistribution_status": "PROJECT_GENERATED_DERIVATIVE_OF_PROVIDER_INPUT",
                "notes": "Generated by scripts/build_highres_profile_mask_mapping.py with seed 42 + len(raw_sample); spatial and reference raw_cell IDs have zero overlap.",
            })
    return release_rows, inventory_rows, summaries


def update_manifests(
    datasets: list[dict[str, str]],
    experiments: list[dict[str, str]],
    sources: list[dict[str, str]],
    heca_summaries: dict[str, dict[str, Any]],
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    dataset_map = {row["record_id"]: row for row in datasets}
    for record_id in HECA_IDS:
        row = dataset_map[record_id]
        summary = heca_summaries[record_id]
        row["provenance_status"] = "CONFIRMED"
        row["unresolved_issue_id"] = ""
        row["evidence_summary"] = append_evidence(
            row["evidence_summary"],
            row["input_file"].split(";")[-1],
            "visualizations/manuscript_audits/dataset_provenance_major_evidence_recovery/heca_reference_lineage_audit.tsv",
        )
        row["notes"] = (
            "Task 4B-2 recovered exact retained-input lineage from embedded obs.study_id and obs.donor_ID fields. "
            f"The input has {summary['n_cells']:,} cells, {len(summary['study_counts'])} contributing study ID(s), "
            f"and local SHA-1 {summary['sha1']}. The historical compressed download was not retained and is not claimed byte-identical to the current archive."
        )
        if record_id in {"LR03_BRCA_FFPE", "LR04_BRCA_FF", "LR05_BRCA_ILC"}:
            row["official_source_url"] = append_evidence(row["official_source_url"], "https://zenodo.org/records/15619143")
        elif record_id == "LR06_CERVICAL":
            row["official_source_url"] = append_evidence(row["official_source_url"], "https://zenodo.org/records/15620903")
        elif record_id == "LR08_INTESTINE":
            row["official_source_url"] = append_evidence(row["official_source_url"], "https://zenodo.org/records/15619143")

    cta = dataset_map[CTA_ID]
    cta["slide_id"] = "V10F24-112"
    cta["capture_area"] = "C1"
    cta["assay"] = "Visium Spatial Gene Expression kit; exact chemistry version not stated"
    cta["tissue_processing"] = "OCT-embedded fresh frozen; 10 micrometre section; 6.5 mm x 6.5 mm capture format"
    cta["evidence_summary"] = append_evidence(
        cta["evidence_summary"],
        "Official Zenodo 15211538 spaceranger_output.zip/Visium_metadata.xlsx: BCSA2 TumB1 -> V10F24-112_C1",
        "Original paper Methods: Space Ranger 1.0.0",
    )
    cta["notes"] = (
        "Task 4B-2 recovered slide V10F24-112 and capture area C1 from the official Zenodo Visium metadata workbook. "
        "The exact Visium chemistry/kit version remains unresolved and is not inferred from Space Ranger 1.0.0. "
        "Original CTA annotation, project registration, and project-generated frozen endpoint remain distinct provenance layers."
    )
    cta["provenance_status"] = "PARTIAL"
    cta["unresolved_issue_id"] = "U008"

    crc = dataset_map["RD04_CRC_B"]
    crc["slide_id"] = "V10A13-206"
    crc["capture_area"] = "C1"
    crc["assay"] = "10x Genomics Spatial 3' v1"
    crc["tissue_processing"] = "fresh frozen; Space Ranger 1.2.0"
    crc["unresolved_issue_id"] = ""
    crc["evidence_summary"] = append_evidence(crc["evidence_summary"], CRC_WEB_SUMMARY)
    crc["notes"] = "Task 4B-2 recovered slide V10A13-206, capture area C1, Spatial 3' v1 chemistry, and Space Ranger 1.2.0 from the official 10x web summary."

    for record_id in VIZGEN_IDS:
        row = dataset_map[record_id]
        row["redistribution_status"] = "PROVIDER_TERMS_FOUND_BUT_AMBIGUOUS"
        row["evidence_summary"] = append_evidence(
            row["evidence_summary"], VIZGEN_SHOWCASE, f"{VIZGEN_BUCKET}/{row['sample_id']}"
        )
        row["notes"] = (
            "Task 4B-2 confirmed the retained sample identity by exact provider folder and official/local cell-count agreement. "
            "Local SHA-1 values are confirmed. Historical download receipt/date, provider original checksum, and unambiguous redistribution permission remain unresolved."
        )

    experiment_map = {row["experiment_id"]: row for row in experiments}
    for row in experiments:
        linked = {row["spatial_input_record_id"], row["reference_input_record_id"], row["evaluation_resource_record_id"]}
        if linked & set(HECA_IDS):
            row["formal_status"] = "CONFIRMED"
            row["status_basis"] = "Task 4B-2 recovered exact retained hECA input lineage, cell-level study mapping, cell counts, and local SHA-1."
            row["unresolved_issue_id"] = ""
        elif CTA_ID in linked:
            row["formal_status"] = "PARTIAL"
            row["status_basis"] = "Sample, slide, and capture area are confirmed; exact Visium chemistry version remains unresolved."
            row["unresolved_issue_id"] = "U008"
        elif linked & set(VIZGEN_IDS):
            row["formal_status"] = "PARTIAL"
            row["status_basis"] = "Exact provider sample and local input are confirmed; provider checksum/download receipt and redistribution status remain partial."
            row["unresolved_issue_id"] = "U007"
        elif "RD04_CRC_B" in linked:
            row["formal_status"] = "CONFIRMED"
            row["status_basis"] = "Official 10x web summary confirms sample, slide V10A13-206, capture C1, chemistry, and Space Ranger version."
            row["unresolved_issue_id"] = ""

    source_keys = {(row["record_id"], row["official_url"]) for row in sources}
    additions = [
        (CTA_ID, "Zenodo", CTA_ZENODO, "processed Space Ranger archive, Visium_metadata.xlsx, slide and capture area", "CONFIRMED"),
        (CTA_ID, "Researchdata.se / SND", CTA_RESEARCHDATA, "raw-data DOI, restricted access, metadata description", "CONFIRMED"),
        (CTA_ID, "npj Precision Oncology", CTA_PAPER, "platform, tissue processing, Space Ranger version, data availability", "CONFIRMED"),
        ("RD04_CRC_B", "10x Genomics", CRC_WEB_SUMMARY, "sample, slide, capture, Spatial 3' v1 chemistry, Space Ranger 1.2.0", "CONFIRMED"),
        ("LR03_BRCA_FFPE", "Zenodo / hECA", "https://zenodo.org/records/15619143", "hECA breast organ export and archive checksum", "CONFIRMED"),
        ("LR04_BRCA_FF", "Zenodo / hECA", "https://zenodo.org/records/15619143", "hECA breast organ export and archive checksum", "CONFIRMED"),
        ("LR05_BRCA_ILC", "Zenodo / hECA", "https://zenodo.org/records/15619143", "hECA breast organ export and archive checksum", "CONFIRMED"),
        ("LR06_CERVICAL", "Zenodo / hECA", "https://zenodo.org/records/15620903", "hECA uterus organ export and archive checksum", "CONFIRMED"),
        ("LR08_INTESTINE", "Zenodo / hECA", "https://zenodo.org/records/15619143", "hECA intestine organ export and archive checksum", "CONFIRMED"),
    ]
    for record_id in VIZGEN_IDS:
        sample = dataset_map[record_id]["sample_id"]
        additions.extend([
            (record_id, "Vizgen", VIZGEN_SHOWCASE, "release identity, per-sample totals, panel, legal notice", "PARTIAL"),
            (record_id, "Vizgen Google Cloud Storage", f"{VIZGEN_BUCKET}/{sample}", "exact provider sample folder and download route", "PARTIAL"),
        ])
    for record_id, authority, url, supports, status in additions:
        key = (record_id, url)
        if key in source_keys:
            continue
        sources.append({
            "record_id": record_id,
            "source_authority": authority,
            "source_type": "official database/provider/publisher page",
            "official_url": url,
            "provider_url": url,
            "accessed_date": GENERATED_DATE,
            "supports_fields": supports,
            "source_status": status,
            "notes": "Added by Task 4B-2 official-source verification.",
        })
        source_keys.add(key)
    return datasets, experiments, sources


def official_verification_rows() -> list[dict[str, str]]:
    rows = [
        (CTA_ID, "slide_id and capture_area", "official sample metadata absent locally", "Zenodo", CTA_ZENODO, "BCSA2 TumB1 = V10F24-112_C1", "CONFIRMED", "Recovered from spaceranger_output.zip/Visium_metadata.xlsx by HTTP range without downloading the full archive."),
        (CTA_ID, "platform, processing, Space Ranger", "partial local README", "original publisher", CTA_PAPER, "Visium; OCT fresh frozen; 10 micrometre; 6.5 mm capture format; Space Ranger 1.0.0", "CONFIRMED", "Exact chemistry version is not stated."),
        (CTA_ID, "raw-data access", "documentation README retained", "Researchdata.se / SND", CTA_RESEARCHDATA, "DOI 10.48723/f4v5-m008; restricted access", "CONFIRMED", "SND_metadata.xlsx is described but not locally retained or openly downloadable."),
        ("RD03_TNBC_PLASMA", "source identity and license", "three official archives retained", "Zenodo", "https://zenodo.org/records/4739739", "CID4465; Visium; CC BY 4.0", "CONFIRMED", "No slide serial or capture-area code is exposed."),
        ("RD04_CRC_B", "slide, capture, chemistry, pipeline", "matrix and spatial archive retained", "10x Genomics", CRC_WEB_SUMMARY, "V10A13-206; C1; Spatial 3' v1; Space Ranger 1.2.0", "CONFIRMED", "Official web summary resolves U010."),
        ("SIM03_HUMAN_LUNG", "accession type", "Task 4B-1 correction", "ENA", "https://www.ebi.ac.uk/ena/browser/view/ERS20065156", "ERS20065156 is a sample accession; SAMEA115633909 is BioSample", "CONFIRMED", "No verified read-run accession."),
    ]
    heca = {
        "LR03_BRCA_FFPE": ("https://zenodo.org/records/17008296", "project Breast export; official archive MD5 df6999cc9c5c7800498d295a9a3a9401"),
        "LR04_BRCA_FF": ("https://zenodo.org/records/17008296", "project Breast export; official archive MD5 df6999cc9c5c7800498d295a9a3a9401"),
        "LR05_BRCA_ILC": ("https://zenodo.org/records/17008296", "project Breast export; official archive MD5 df6999cc9c5c7800498d295a9a3a9401"),
        "LR06_CERVICAL": ("https://zenodo.org/records/15620903", "organ Uterus export; official archive MD5 26a3067c177094e8627624d25275a121"),
        "LR08_INTESTINE": ("https://zenodo.org/records/15619143", "organ Intestine export; official archive MD5 b8604d168d24d9d64862cddaa3bbad0c"),
    }
    for record_id, (url, value) in heca.items():
        rows.append((record_id, "formal hECA export", "actual local h5ad and embedded metadata confirmed", "Zenodo / hECA", url, value, "CONFIRMED", "Zenodo record is CC BY 4.0; retained input has its own local SHA-1."))
    viz_samples = {
        "MERS01_BREAST": "HumanBreastCancerPatient1",
        "MERS02_COLON": "HumanColonCancerPatient1",
        "MERS03_LUNG": "HumanLungCancerPatient1",
        "MERS04_MELANOMA1": "HumanMelanomaPatient1",
        "MERS05_MELANOMA2": "HumanMelanomaPatient2",
    }
    for record_id, sample in viz_samples.items():
        rows.append((record_id, "release and per-sample route", "local sample files and cell totals confirmed", "Vizgen", f"{VIZGEN_BUCKET}/{sample}", sample, "PARTIAL", "Official folder and showcase confirm identity; provider checksum and unambiguous redistribution permission are absent."))
    return [dict(zip(
        ["record_id", "field", "local_evidence_status", "official_source_authority", "official_url", "verified_value", "verification_date", "status", "notes"],
        [*row[:6], GENERATED_DATE, *row[6:]],
    )) for row in rows]


def minor_metadata_rows() -> list[dict[str, str]]:
    return [
        {
            "record_id": "RD03_TNBC_PLASMA", "formal_sample_id": "CID4465", "spatial_library_id": "", "slide_id": "", "capture_area": "", "processed_object_name": "CID4465_filtered_count_matrix; CID4465_spatial", "status": "PARTIAL", "local_evidence": "data/raw/cytospace_fig2d_tme/brca_ST_Wu_TNBC_Zenodo/filtered_count_matrices.tar.gz; data/raw/cytospace_fig2d_tme/brca_ST_Wu_TNBC_Zenodo/spatial.tar.gz; data/raw/cytospace_fig2d_tme/brca_ST_Wu_TNBC_Zenodo/metadata.tar.gz", "official_evidence": "https://zenodo.org/records/4739739", "notes": "The official archive confirms CID4465 and Visium but does not expose a slide serial or capture-area code. Both remain blank rather than guessed."
        },
        {
            "record_id": "RD04_CRC_B", "formal_sample_id": "Parent_Visium_Human_ColorectalCancer", "spatial_library_id": "Parent_Visium_Human_ColorectalCancer", "slide_id": "V10A13-206", "capture_area": "C1", "processed_object_name": "Parent_Visium_Human_ColorectalCancer_filtered_feature_bc_matrix.h5", "status": "CONFIRMED", "local_evidence": "data/raw/cytospace_fig2d_tme/crc_ST_10x_parent_visium/Parent_Visium_Human_ColorectalCancer_filtered_feature_bc_matrix.h5; data/raw/cytospace_fig2d_tme/crc_ST_10x_parent_visium/Parent_Visium_Human_ColorectalCancer_spatial.tar.gz", "official_evidence": CRC_WEB_SUMMARY, "notes": "The official web summary explicitly reports Slide Serial Number V10A13-206-C1; the audit stores slide and capture as separate fields."
        },
    ]


def status_reassessment_rows(datasets: list[dict[str, str]], experiments: list[dict[str, str]], evidence_ids: dict[str, str]) -> list[dict[str, str]]:
    old_datasets = {row["record_id"]: row for row in read_tsv(DATASET_V01)}
    target_ids = {CTA_ID, *VIZGEN_IDS, *HECA_IDS, *MINOR_IDS}
    rows: list[dict[str, str]] = []
    for row in datasets:
        if row["record_id"] not in target_ids:
            continue
        old = old_datasets[row["record_id"]]
        record_id = row["record_id"]
        if record_id in HECA_IDS:
            reason = "Embedded cell-level study/donor lineage, exact counts, formal hECA export record, and local SHA-1 recovered."
            remaining = "historical compressed download receipt not retained; no effect on exact retained-input lineage"
            safe = "CONFIRMED_EXACT_LINEAGE for the retained local input"
        elif record_id == CTA_ID:
            reason = "Official Visium metadata recovered slide V10F24-112 and capture C1."
            remaining = "exact Visium chemistry/kit version; original local Space Ranger bundle"
            safe = "Exact sample/slide/capture claim; generic Visium kit and Space Ranger 1.0.0 only"
        elif record_id in VIZGEN_IDS:
            reason = "Exact provider folder and official/local cell-count agreement confirm sample identity."
            remaining = "provider checksum; historical download receipt/date; unambiguous redistribution permission"
            safe = "Provider-hosted May 2022 FFPE MERSCOPE sample; local input SHA-1"
        elif record_id == "RD04_CRC_B":
            reason = "Official 10x web summary recovered slide, capture, chemistry, and pipeline version."
            remaining = ""
            safe = "Exact official 10x sample metadata"
        else:
            reason = "Retained official archives and source page reconfirm sample identity."
            remaining = "slide serial; capture area"
            safe = "CID4465 Visium sample identity without slide/capture claim"
        support = [e["evidence_id"] for e in build_cached_evidence_rows if record_id in e["record_id"].split(";")]
        rows.append({
            "record_type": "dataset",
            "record_id": record_id,
            "task4b1_status": old["provenance_status"],
            "task4b2_status": row["provenance_status"],
            "status_change": "CHANGED" if old["provenance_status"] != row["provenance_status"] else "UNCHANGED",
            "change_reason": reason,
            "supporting_evidence_ids": ";".join(support),
            "remaining_missing_fields": remaining,
            "safe_claim_level": safe,
            "notes": "No unsupported field was promoted.",
        })
    old_experiments = {row["experiment_id"]: row for row in read_tsv(EXPERIMENT_V01)}
    for row in experiments:
        linked = {row["spatial_input_record_id"], row["reference_input_record_id"], row["evaluation_resource_record_id"]}
        if not (linked & target_ids):
            continue
        old = old_experiments[row["experiment_id"]]
        rows.append({
            "record_type": "experiment",
            "record_id": row["experiment_id"],
            "task4b1_status": old["formal_status"],
            "task4b2_status": row["formal_status"],
            "status_change": "CHANGED" if old["formal_status"] != row["formal_status"] else "UNCHANGED",
            "change_reason": row["status_basis"],
            "supporting_evidence_ids": "",
            "remaining_missing_fields": row["unresolved_issue_id"],
            "safe_claim_level": row["status_basis"],
            "notes": "Experiment retained; no provenance-incomplete experiment was deleted.",
        })
    return rows


def unresolved_rows() -> list[dict[str, str]]:
    common_local = "data/; configs/; scripts/; logs/; result/; visualizations/; retained archives and metadata files"
    return [
        {"issue_id": "B2-001", "severity": "MAJOR", "record_ids": CTA_ID, "category": "CTA chemistry and raw bundle", "searched_locations": "CTA ST_data.RData, CTA repository files, SND documentation, original paper, Zenodo processed archive listing and metadata workbook", "searched_terms": "BCSA2TumB1; V10F24-112; chemistry; Space Ranger; Visium; slide; capture", "evidence_found": "slide V10F24-112; capture C1; Space Ranger 1.0.0; Visium kit; tissue processing", "missing_information": "exact Visium chemistry/kit version and original locally retained Space Ranger bundle", "reason_unresolved": "Neither the official workbook nor the retained processed Seurat object states a chemistry version; the full official archive was not copied locally.", "impact_on_reproducibility": "Low for the frozen processed input; moderate for rebuilding from raw reads.", "impact_on_submission": "Use a bounded platform statement and omit chemistry version.", "safe_current_action": "Keep CTA record PARTIAL; state only verified sample/slide/capture and Space Ranger version.", "recommended_future_action": "Obtain SND_metadata.xlsx or provider confirmation of exact chemistry if required.", "status": "OPEN"},
        {"issue_id": "B2-002", "severity": "MAJOR", "record_ids": ";".join(VIZGEN_IDS), "category": "Vizgen original download evidence", "searched_locations": common_local, "searched_terms": "sample IDs; Vizgen; MERSCOPE; download; curl; wget; archive; checksum; browser metadata", "evidence_found": "exact official per-sample folder URLs, matching local filenames, matching official/local cell totals, local SHA-1", "missing_information": "historical download receipt/date and provider-issued checksum", "reason_unresolved": "No download log, archive, browser receipt, or provider checksum manifest was retained; the current GCS bucket blocks anonymous object listing.", "impact_on_reproducibility": "Local inputs are frozen by SHA-1, but byte-level comparison to an original provider checksum is unavailable.", "impact_on_submission": "Describe provider release and local processing; do not claim provider checksum verification.", "safe_current_action": "Retain PARTIAL provenance with exact local SHA-1.", "recommended_future_action": "Preserve a new authenticated provider download receipt/checksum manifest if access is re-established.", "status": "OPEN"},
        {"issue_id": "B2-003", "severity": "MAJOR", "record_ids": ";".join(VIZGEN_IDS), "category": "Vizgen redistribution", "searched_locations": "project license/terms/readme files; Vizgen showcase legal notice and roadmap", "searched_terms": "license; terms; EULA; redistribution; use in any way; legal notice", "evidence_found": "provider legal notice and general use language", "missing_information": "explicit permission to redistribute provider raw files", "reason_unresolved": "Use language does not unambiguously grant redistribution rights.", "impact_on_reproducibility": "Readers can use the provider access route; raw inputs should not be redistributed by this project without clarification.", "impact_on_submission": "Data Availability must link to Vizgen rather than promise redistribution.", "safe_current_action": "PROVIDER_TERMS_FOUND_BUT_AMBIGUOUS", "recommended_future_action": "Request written redistribution clarification from Vizgen.", "status": "OPEN"},
        {"issue_id": "B2-004", "severity": "MINOR", "record_ids": "RD03_TNBC_PLASMA", "category": "CID4465 slide metadata", "searched_locations": "retained Zenodo matrices, metadata.tar.gz, spatial.tar.gz, preparation scripts, Zenodo record, source literature", "searched_terms": "CID4465; slide; capture; Visium; library_id; web_summary; metrics_summary", "evidence_found": "CID4465 sample identity, TNBC context, Visium platform, CC BY 4.0 source", "missing_information": "slide serial and capture area", "reason_unresolved": "The source archives expose CID4465 but no slide/capture identifiers; multiple regions may share a slide.", "impact_on_reproducibility": "None for the retained frozen sample files.", "impact_on_submission": "Omit slide and capture fields.", "safe_current_action": "Keep blank fields; do not write unknown or infer from neighboring sections.", "recommended_future_action": "Seek authors' original Space Ranger web summary or sample sheet.", "status": "OPEN"},
        {"issue_id": "B2-005", "severity": "MINOR", "record_ids": "SIM03_HUMAN_LUNG", "category": "human-lung read-run accession", "searched_locations": "Task 4A/4B-1 accession audit and official ENA sample record", "searched_terms": "ERS20065156; SAMEA115633909; run accession", "evidence_found": "ERS20065156 is sample; SAMEA115633909 is BioSample", "missing_information": "verified read-run accession", "reason_unresolved": "No run accession was verified; it is not required to identify the frozen spatial sample.", "impact_on_reproducibility": "Low.", "impact_on_submission": "Use sample and BioSample wording only.", "safe_current_action": "Retain PARTIAL and leave run_accession blank.", "recommended_future_action": "Resolve only if raw-read reconstruction is required.", "status": "OPEN"},
        {"issue_id": "B2-006", "severity": "CRITICAL", "record_ids": "ALL", "category": "formal manuscript/bibliography linkage and manifest freeze", "searched_locations": "Task 4A and Task 4B-1 audit outputs only; formal manuscript and bibliography intentionally not accessed", "searched_terms": "existing audit occurrence and recommendation records", "evidence_found": "complete 37-record and 40-experiment draft manifests", "missing_information": "manual manuscript/bibliography reconciliation, user review, formal freeze", "reason_unresolved": "This task is prohibited from reading or editing the formal manuscript and bibliography.", "impact_on_reproducibility": "None for current manifest integrity; prevents final submission-ready freeze.", "impact_on_submission": "Data Availability and formal citations are not ready to declare final.", "safe_current_action": "Keep status DRAFT_NOT_FROZEN and data_availability_ready=false.", "recommended_future_action": "Perform a separate authorized manual manuscript/bibliography reconciliation after user review.", "status": "OPEN"},
    ]


def manuscript_recommendations() -> list[dict[str, str]]:
    return [
        {"recommendation_id": "R001", "record_id": "SIM03_HUMAN_LUNG", "likely_section": "Methods / Data Availability", "verified_information": "ERS20065156 is an ENA sample accession; SAMEA115633909 is the BioSample; no read-run accession is verified.", "remaining_uncertainty": "read-run accession", "recommended_manual_wording": "The human-lung spatial sample is identified by ENA sample ERS20065156 and BioSample SAMEA115633909; no read-run accession was assigned in the retained provenance record.", "prohibited_wording": "ERS20065156 is a sequencing run.", "priority": "HIGH", "evidence_ids": "Task 4B-1 human_lung_accession_correction.tsv"},
        {"recommendation_id": "R002", "record_id": ";".join(HECA_IDS), "likely_section": "Methods / reference datasets", "verified_information": "Actual hECA inputs, SHA-1 values, exact cell counts, cell-level study_id mapping, and donor mapping are recovered.", "remaining_uncertainty": "historical compressed download receipt; no byte-match claim to current archive", "recommended_manual_wording": "The reference used the retained hECA export for the relevant tissue compartment; exact source-study membership was reconstructed from the input object's cell-level study_id metadata and the input was frozen by SHA-1.", "prohibited_wording": "The project independently merged the cited studies, or the local file was byte-identical to the current compressed Zenodo archive.", "priority": "HIGH", "evidence_ids": "heca_reference_lineage_audit.tsv"},
        {"recommendation_id": "R003", "record_id": ";".join(VIZGEN_IDS), "likely_section": "Methods / Data Availability", "verified_information": "Five exact provider-hosted FFPE MERSCOPE sample folders and matching local input cell totals are confirmed.", "remaining_uncertainty": "provider checksum, historical download receipt/date, redistribution permission", "recommended_manual_wording": "Five samples were obtained from Vizgen's May 2022 MERSCOPE FFPE Human Immuno-Oncology Data Release; retained local inputs are reported with project SHA-1 values and project-generated disjoint same-assay splits.", "prohibited_wording": "The raw Vizgen files are freely redistributable or provider checksums were verified.", "priority": "HIGH", "evidence_ids": "vizgen_release_evidence.tsv;vizgen_sample_file_inventory.tsv"},
        {"recommendation_id": "R004", "record_id": CTA_ID, "likely_section": "Methods / biological application", "verified_information": "BCSA2TumB1 maps to slide V10F24-112, capture C1; Space Ranger 1.0.0; OCT fresh-frozen 10 micrometre section.", "remaining_uncertainty": "exact Visium chemistry/kit version", "recommended_manual_wording": "The BCSA2TumB1 Visium section corresponds to slide V10F24-112, capture area C1; the retained source reports Space Ranger 1.0.0, while the exact chemistry version was not preserved.", "prohibited_wording": "A specific Visium chemistry version not present in the evidence.", "priority": "HIGH", "evidence_ids": "cta_metadata_recovery.tsv"},
        {"recommendation_id": "R005", "record_id": ";".join(VIZGEN_IDS), "likely_section": "Methods / high-resolution benchmark", "verified_information": "Spatial-cell/reference-pool disjoint split is generated by this project with seed 42 + len(raw_sample) and zero raw_cell overlap.", "remaining_uncertainty": "none for split generation", "recommended_manual_wording": "Spatial cells and the same-assay reference pool were generated as disjoint project derivatives from each provider sample; they are not independent scRNA-seq references and have no separate accession.", "prohibited_wording": "Vizgen supplied an independent single-cell reference dataset.", "priority": "HIGH", "evidence_ids": "scripts/build_highres_profile_mask_mapping.py;highres_profile_mask_mapping_manifest.json"},
        {"recommendation_id": "R006", "record_id": "RD03_TNBC_PLASMA", "likely_section": "Methods / supplementary dataset table", "verified_information": "CID4465 identity and Visium source are confirmed.", "remaining_uncertainty": "slide serial and capture area", "recommended_manual_wording": "List CID4465 and its Zenodo source; omit slide and capture-area columns for this row or leave them empty in a machine-readable table.", "prohibited_wording": "unknown slide, guessed capture area, or values copied from another section.", "priority": "MEDIUM", "evidence_ids": "minor_spatial_metadata_audit.tsv"},
        {"recommendation_id": "R007", "record_id": "RD04_CRC_B", "likely_section": "Methods / supplementary dataset table", "verified_information": "Slide V10A13-206, capture C1, Spatial 3' v1, Space Ranger 1.2.0.", "remaining_uncertainty": "none for requested minor metadata", "recommended_manual_wording": "The Parent_Visium_Human_ColorectalCancer section used slide V10A13-206, capture area C1, and was processed with Space Ranger 1.2.0 using Spatial 3' v1 chemistry.", "prohibited_wording": "slide unresolved", "priority": "MEDIUM", "evidence_ids": "minor_spatial_metadata_audit.tsv"},
    ]


def data_availability_rows(datasets: list[dict[str, str]]) -> list[dict[str, str]]:
    rows = []
    for row in datasets:
        record_id = row["record_id"]
        if record_id in VIZGEN_IDS:
            readiness = "REDISTRIBUTION_UNRESOLVED"
            ready = "false"
            blocker = "Provider checksum/download receipt and explicit redistribution permission are unresolved."
            safe = "Available through Vizgen's May 2022 FFPE release; local project inputs are identified by SHA-1."
        elif record_id == CTA_ID:
            readiness = "PARTIAL_METADATA"
            ready = "false"
            blocker = "Exact chemistry version unresolved; raw SND access is restricted."
            safe = "BCSA2TumB1, slide V10F24-112, capture C1; processed data Zenodo 15211538; raw data SND DOI 10.48723/f4v5-m008."
        elif record_id in HECA_IDS:
            readiness = "READY_PUBLIC_ACCESSION"
            ready = "true"
            blocker = ""
            safe = "Retained hECA input with exact local SHA-1 and cell-level study lineage; official hECA Zenodo export records are listed."
        elif row["provenance_status"] == "PROJECT_GENERATED":
            readiness = "READY_PROJECT_DERIVATIVE_DESCRIPTION"
            ready = "true"
            blocker = ""
            safe = "Project-generated derivative with no independent accession; generation lineage is recorded in the manifest."
        elif row["provenance_status"] == "CONFIRMED":
            readiness = "READY_PUBLIC_ACCESSION" if any(row.get(k, "") for k in ["study_accession", "project_accession", "geo_accession", "arrayexpress_accession", "official_source_url"]) else "READY_PROVIDER_RELEASE"
            ready = "true"
            blocker = ""
            safe = "Use the formal dataset name, verified accession/release, and official URL from the v0.2 manifest."
        else:
            readiness = "PARTIAL_METADATA"
            ready = "false"
            blocker = row.get("unresolved_issue_id", "") or "Record remains partial."
            safe = "Use only verified fields in the v0.2 manifest and omit unresolved values."
        rows.append({
            "record_id": record_id,
            "formal_dataset_name": row["formal_dataset_name"],
            "source_database_or_provider": row["database"],
            "accession_or_release": ";".join(filter(None, [row["study_accession"], row["project_accession"], row["geo_accession"], row["arrayexpress_accession"], row["formal_release_name"]])),
            "official_url": row["official_source_url"],
            "local_input_confirmed": "true" if row["input_file"] else "false",
            "local_input_sha1": row["input_file_sha1"],
            "redistribution_status": row["redistribution_status"],
            "processed_derivative_status": row["project_generated_derivative"],
            "ready_for_data_availability": ready,
            "readiness_class": readiness,
            "blocking_issue": blocker,
            "safe_description": safe,
            "notes": "Drafting inventory only; this task does not edit Data Availability.",
        })
    return rows


def link_integrity_rows(datasets: list[dict[str, str]], experiments: list[dict[str, str]]) -> list[dict[str, Any]]:
    ids = [row["record_id"] for row in datasets]
    id_set = set(ids)
    v0_1_by_id = {row["record_id"]: row for row in read_tsv(DATASET_V01)}
    exp_ids = [row["experiment_id"] for row in experiments]
    spatial_orphans = [row["experiment_id"] for row in experiments if row["spatial_input_record_id"] and row["spatial_input_record_id"] not in id_set]
    ref_orphans = [row["experiment_id"] for row in experiments if row["reference_input_record_id"] and row["reference_input_record_id"] not in id_set]
    eval_orphans = [row["experiment_id"] for row in experiments if row["evaluation_resource_record_id"] and row["evaluation_resource_record_id"] not in id_set]
    invented = []
    accession_fields = [
        "study_accession", "project_accession", "experiment_accession",
        "run_accession", "sample_accession", "biosample_accession",
        "geo_accession", "arrayexpress_accession", "ena_accession",
    ]
    for row in datasets:
        if row["provenance_status"] == "PROJECT_GENERATED":
            upstream = v0_1_by_id[row["record_id"]]
            if any(row.get(field, "") != upstream.get(field, "") for field in accession_fields):
                invented.append(row["record_id"])
    human_lung = next(row for row in datasets if row["record_id"] == "SIM03_HUMAN_LUNG")
    incorrect_accessions = []
    if human_lung["run_accession"] or human_lung["sample_accession"] != "ERS20065156" or human_lung["biosample_accession"] != "SAMEA115633909":
        incorrect_accessions.append("SIM03_HUMAN_LUNG")
    checks = [
        ("dataset_record_count", 37, len(datasets), len(datasets) == 37, ""),
        ("experiment_record_count", 40, len(experiments), len(experiments) == 40, ""),
        ("duplicate_dataset_record_ids", 0, len(ids) - len(set(ids)), len(ids) == len(set(ids)), ""),
        ("duplicate_experiment_ids", 0, len(exp_ids) - len(set(exp_ids)), len(exp_ids) == len(set(exp_ids)), ""),
        ("orphan_spatial_links", 0, len(spatial_orphans), not spatial_orphans, ";".join(spatial_orphans)),
        ("orphan_reference_links", 0, len(ref_orphans), not ref_orphans, ";".join(ref_orphans)),
        ("orphan_evaluation_resource_links", 0, len(eval_orphans), not eval_orphans, ";".join(eval_orphans)),
        ("project_generated_records_with_invented_accessions", 0, len(invented), not invented, ";".join(invented)),
        ("incorrect_accession_field_assignments", 0, len(incorrect_accessions), not incorrect_accessions, ";".join(incorrect_accessions)),
    ]
    return [{"check": name, "expected": expected, "observed": observed, "status": "PASS" if ok else "FAIL", "details": details} for name, expected, observed, ok, details in checks]


def build_payloads() -> dict[str, bytes]:
    global build_cached_evidence_rows
    datasets = read_tsv(DATASET_V01)
    experiments = read_tsv(EXPERIMENT_V01)
    sources = read_tsv(SOURCE_V01)
    dataset_fields = list(datasets[0])
    experiment_fields = list(experiments[0])
    source_fields = list(sources[0])

    evidence_rows, evidence_ids = build_evidence_inventory()
    build_cached_evidence_rows = evidence_rows
    heca_lineage, heca_inputs, heca_summaries = heca_audits(evidence_ids)
    vizgen_release, vizgen_inventory, vizgen_summaries = vizgen_audits(evidence_ids)
    cta_rows = cta_recovery(evidence_ids)
    datasets, experiments, sources = update_manifests(datasets, experiments, sources, heca_summaries)
    datasets = sorted(datasets, key=lambda row: row["record_id"])
    experiments = sorted(experiments, key=lambda row: row["experiment_id"])
    sources = sorted(sources, key=lambda row: (row["record_id"], row["official_url"], row["source_authority"]))

    reassessment = status_reassessment_rows(datasets, experiments, evidence_ids)
    unresolved = unresolved_rows()
    recommendations = manuscript_recommendations()
    availability = data_availability_rows(datasets)
    integrity = link_integrity_rows(datasets, experiments)
    failed_integrity = [row for row in integrity if row["status"] != "PASS"]
    if failed_integrity:
        raise ValueError(f"manifest link integrity failed: {failed_integrity}")

    official = official_verification_rows()
    minor = minor_metadata_rows()
    status_counts = Counter(row["provenance_status"] for row in datasets)
    exp_status_counts = Counter(row["formal_status"] for row in experiments)
    unresolved_severity = Counter(row["severity"] for row in unresolved)

    payloads: dict[str, bytes] = {}
    payloads["local_evidence_inventory.tsv"] = tsv_bytes(evidence_rows, [
        "evidence_id", "record_id", "category", "file_path", "file_name", "file_type", "file_size", "modified_time", "sha1", "evidence_priority", "supports_fields", "read_only", "notes"
    ])
    payloads["cta_metadata_recovery.tsv"] = tsv_bytes(cta_rows, [
        "record_id", "field", "task4b1_value", "recovered_value", "status", "local_evidence_path", "local_evidence_sha1", "official_evidence", "evidence_priority", "notes"
    ])
    payloads["vizgen_release_evidence.tsv"] = tsv_bytes(vizgen_release, [
        "record_id", "sample_identifier", "provider_sample_name", "provider", "platform", "assay", "processing", "release", "panel", "official_cell_total", "local_metadata_row_count", "official_transcript_total", "rin", "dv200", "provider_release_page", "per_sample_download_url", "download_url_preserved", "historical_download_receipt_preserved", "download_date_preserved", "provider_checksum_preserved", "license_or_terms_evidence", "redistribution_status", "same_assay_split_status", "provenance_status", "notes"
    ])
    payloads["vizgen_sample_file_inventory.tsv"] = tsv_bytes(vizgen_inventory, [
        "record_id", "sample_identifier", "file_role", "original_filename", "local_raw_directory", "local_processed_directory", "file_path", "file_size", "file_sha1", "provider_checksum", "provider_checksum_status", "download_url", "download_date", "evidence_id", "redistribution_status", "notes"
    ])
    payloads["heca_reference_lineage_audit.tsv"] = tsv_bytes(heca_lineage, [
        "record_id", "actual_input_file", "actual_input_confirmed", "formal_heca_export_identity", "single_export_or_project_merge", "project_merge_script", "cell_count", "gene_count", "contributing_study_counts", "contributing_study_count", "cell_level_source_mapping", "donor_level_mapping", "accession_field_present", "input_sha1", "official_export_url", "official_archive_md5", "historical_archive_byte_match", "lineage_status", "supporting_evidence_ids", "safe_claim_level", "notes"
    ])
    payloads["heca_input_file_inventory.tsv"] = tsv_bytes(heca_inputs, [
        "record_id", "file_path", "file_name", "file_size", "modified_time", "sha1", "n_cells", "n_genes", "obs_columns", "study_counts", "donor_count", "organ_counts", "evidence_id", "read_only"
    ])
    payloads["minor_spatial_metadata_audit.tsv"] = tsv_bytes(minor, [
        "record_id", "formal_sample_id", "spatial_library_id", "slide_id", "capture_area", "processed_object_name", "status", "local_evidence", "official_evidence", "notes"
    ])
    payloads["official_source_verification.tsv"] = tsv_bytes(official, [
        "record_id", "field", "local_evidence_status", "official_source_authority", "official_url", "verified_value", "verification_date", "status", "notes"
    ])
    payloads["provenance_status_reassessment.tsv"] = tsv_bytes(reassessment, [
        "record_type", "record_id", "task4b1_status", "task4b2_status", "status_change", "change_reason", "supporting_evidence_ids", "remaining_missing_fields", "safe_claim_level", "notes"
    ])
    payloads["unresolved_after_task4b2.tsv"] = tsv_bytes(unresolved, [
        "issue_id", "severity", "record_ids", "category", "searched_locations", "searched_terms", "evidence_found", "missing_information", "reason_unresolved", "impact_on_reproducibility", "impact_on_submission", "safe_current_action", "recommended_future_action", "status"
    ])
    payloads["manuscript_manual_update_recommendations.tsv"] = tsv_bytes(recommendations, [
        "recommendation_id", "record_id", "likely_section", "verified_information", "remaining_uncertainty", "recommended_manual_wording", "prohibited_wording", "priority", "evidence_ids"
    ])
    payloads["data_availability_evidence_inventory.tsv"] = tsv_bytes(availability, [
        "record_id", "formal_dataset_name", "source_database_or_provider", "accession_or_release", "official_url", "local_input_confirmed", "local_input_sha1", "redistribution_status", "processed_derivative_status", "ready_for_data_availability", "readiness_class", "blocking_issue", "safe_description", "notes"
    ])
    payloads["dataset_manifest_draft_v0_2.tsv"] = tsv_bytes(datasets, dataset_fields)
    payloads["experiment_manifest_draft_v0_2.tsv"] = tsv_bytes(experiments, experiment_fields)
    payloads["source_registry_draft_v0_2.tsv"] = tsv_bytes(sources, source_fields)
    payloads["manifest_link_integrity.tsv"] = tsv_bytes(integrity, ["check", "expected", "observed", "status", "details"])

    summary = {
        "task": "Task 4B-2 major provenance evidence recovery",
        "decision": "CONDITIONAL PASS",
        "manifest_version": "0.2",
        "manifest_status": "DRAFT_NOT_FROZEN",
        "dataset_record_count": len(datasets),
        "experiment_record_count": len(experiments),
        "dataset_status_counts": dict(sorted(status_counts.items())),
        "experiment_status_counts": dict(sorted(exp_status_counts.items())),
        "heca_exact_lineage_records": len(heca_lineage),
        "cta_slide_capture_recovered": True,
        "cta_exact_chemistry_recovered": False,
        "vizgen_sample_identity_confirmed": len(vizgen_release),
        "vizgen_provider_checksum_confirmed": 0,
        "vizgen_redistribution_confirmed": 0,
        "crc_slide_capture_recovered": True,
        "cid4465_slide_capture_recovered": False,
        "remaining_issue_count": len(unresolved),
        "remaining_issue_severity_counts": dict(sorted(unresolved_severity.items())),
        "link_integrity": "PASS",
        "data_availability_ready": False,
        "manifest_frozen": False,
        "overall_task4_resolved": False,
        "formal_manuscript_accessed": False,
        "formal_bibliography_accessed": False,
        "experimental_outputs_modified": False,
        "figures_modified": False,
    }
    payloads["audit_summary.json"] = json_bytes(summary)

    readme = f"""# Task 4B-2 Major Provenance Evidence Recovery

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

The official values were verified on {GENERATED_DATE}. The generator records URLs but does not require network access during deterministic regeneration.
"""
    payloads["README.md"] = text_bytes(readme)

    final_decision = f"""# Task 4B-2 Final Decision

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
"""
    payloads["final_decision.md"] = text_bytes(final_decision)

    upstream_a = {p.name: sha1_file(p) for p in sorted(TASK4A_DIR.iterdir()) if p.is_file()}
    upstream_b = {p.name: sha1_file(p) for p in sorted(TASK4B1_DIR.iterdir()) if p.is_file()}
    generated_hashes = {name: sha1_bytes(data) for name, data in sorted(payloads.items())}
    metadata = {
        "task": "Task 4B-2 major provenance evidence recovery",
        "manifest_version": "0.2",
        "manifest_status": "DRAFT_NOT_FROZEN",
        "source_manifest_version": "0.1",
        "dataset_record_count": 37,
        "experiment_record_count": 40,
        "formal_manuscript_accessed": False,
        "formal_manuscript_modified": False,
        "formal_bibliography_accessed": False,
        "formal_bibliography_modified": False,
        "experimental_outputs_modified": False,
        "figures_modified": False,
        "task4a_outputs_overwritten": False,
        "task4b1_outputs_overwritten": False,
        "github_repository_access_required": False,
        "data_availability_ready": False,
        "manifest_frozen": False,
        "generated_date": GENERATED_DATE,
        "generator_script_sha1": sha1_file(Path(__file__)),
        "repository_local_commit_sha": git_commit(),
        "task4a_file_sha1": upstream_a,
        "task4b1_file_sha1": upstream_b,
        "task4b1_dataset_manifest_sha1": sha1_file(DATASET_V01),
        "task4b1_experiment_manifest_sha1": sha1_file(EXPERIMENT_V01),
        "task4b1_source_registry_sha1": sha1_file(SOURCE_V01),
        "v0_2_dataset_manifest_sha1": generated_hashes["dataset_manifest_draft_v0_2.tsv"],
        "v0_2_experiment_manifest_sha1": generated_hashes["experiment_manifest_draft_v0_2.tsv"],
        "v0_2_source_registry_sha1": generated_hashes["source_registry_draft_v0_2.tsv"],
        "official_source_verification_sha1": generated_hashes["official_source_verification.tsv"],
        "generated_file_sha1_excluding_manifest_metadata": generated_hashes,
        "deterministic_double_build_identical": True,
        "hash_scope_note": "manifest_metadata.json excludes its own SHA-1 to avoid self-reference.",
    }
    payloads["manifest_metadata.json"] = json_bytes(metadata)
    return payloads


build_cached_evidence_rows: list[dict[str, Any]] = []


def main() -> None:
    if not DATASET_V01.is_file() or not EXPERIMENT_V01.is_file() or not SOURCE_V01.is_file():
        raise FileNotFoundError("Task 4B-1 draft inputs are incomplete")
    first = build_payloads()
    second = build_payloads()
    if first != second:
        differing = sorted(name for name in set(first) | set(second) if first.get(name) != second.get(name))
        raise RuntimeError(f"non-deterministic build: {differing}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    expected = set(first)
    for path in OUT_DIR.iterdir():
        if path.is_file() and path.name not in expected:
            path.unlink()
    for name, data in first.items():
        (OUT_DIR / name).write_bytes(data)
    print("TASK4B2_BUILD_PASS")
    print(f"output_dir={rel(OUT_DIR)}")
    print(f"files={len(first)} dataset_records=37 experiment_records=40")
    print("decision=CONDITIONAL PASS link_integrity=PASS")


if __name__ == "__main__":
    main()
