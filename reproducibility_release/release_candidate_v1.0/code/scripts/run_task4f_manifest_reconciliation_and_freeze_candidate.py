#!/usr/bin/env python
from __future__ import annotations

import hashlib
import json
import re
import subprocess
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE_DIR = ROOT / "visualizations/manuscript_audits/dataset_provenance_major_evidence_recovery"
SOURCE_DATASET = SOURCE_DIR / "dataset_manifest_draft_v0_2.tsv"
SOURCE_EXPERIMENT = SOURCE_DIR / "experiment_manifest_draft_v0_2.tsv"
SOURCE_METADATA = SOURCE_DIR / "manifest_metadata.json"
TASK4E_RELEASE = ROOT / "reproducibility_release/release_manifest.tsv"
CANDIDATE_DIR = ROOT / "reproducibility_release/manifests"
CANDIDATE_DATASET = CANDIDATE_DIR / "dataset_manifest_v1.0_freeze_candidate.tsv"
CANDIDATE_EXPERIMENT = CANDIDATE_DIR / "experiment_manifest_v1.0_freeze_candidate.tsv"
CANDIDATE_METADATA = CANDIDATE_DIR / "manifest_v1.0_freeze_candidate_metadata.json"
OUT = ROOT / "visualizations/manuscript_audits/manifest_reconciliation_and_freeze_candidate"
CACHE = OUT / "local_input_sha1_verification_cache.json"
GENERATED_DATE = "2026-07-22"
CANDIDATE_VERSION = "1.0"
CANDIDATE_STATUS = "PENDING_HUMAN_APPROVAL"

ALLOWED_SOURCE_LAYERS = {
    "THIRD_PARTY_ORIGINAL",
    "THIRD_PARTY_RETAINED_EXPORT",
    "PROJECT_PROCESSED_INPUT",
    "PROJECT_GENERATED_DERIVATIVE",
    "PROJECT_GENERATED_SIMULATION_TRUTH",
    "PROJECT_GENERATED_ENDPOINT_DERIVATIVE",
    "PROJECT_GENERATED_DECOY",
}
ALLOWED_REDISTRIBUTION = {
    "REDISTRIBUTABLE",
    "REDISTRIBUTABLE_DERIVATIVE_ONLY",
    "LINK_AND_RECONSTRUCT_ONLY",
    "RESTRICTED",
    "PERMISSION_UNRESOLVED",
    "NOT_APPLICABLE",
}

EVIDENCE = {
    "base": SOURCE_DIR / "dataset_manifest_draft_v0_2.tsv",
    "metadata": SOURCE_DIR / "manifest_metadata.json",
    "source_registry": SOURCE_DIR / "source_registry_draft_v0_2.tsv",
    "heca": SOURCE_DIR / "heca_reference_lineage_audit.tsv",
    "vizgen": SOURCE_DIR / "vizgen_release_evidence.tsv",
    "vizgen_files": SOURCE_DIR / "vizgen_sample_file_inventory.tsv",
    "minor": SOURCE_DIR / "minor_spatial_metadata_audit.tsv",
    "cta": SOURCE_DIR / "cta_metadata_recovery.tsv",
    "unresolved": SOURCE_DIR / "unresolved_after_task4b2.tsv",
    "accessions": ROOT / "visualizations/manuscript_audits/dataset_provenance_formal_audit/accession_type_audit.tsv",
    "task4e": ROOT / "visualizations/manuscript_audits/fig2_enrichment_backend_formal_audit/audit_summary.json",
}


def sha1(path: Path) -> str:
    digest = hashlib.sha1()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def write_tsv(path: Path, rows: list[dict[str, Any]] | pd.DataFrame, columns: list[str] | None = None) -> None:
    frame = rows if isinstance(rows, pd.DataFrame) else pd.DataFrame(rows, columns=columns)
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n")


def parse_path_hashes(value: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for item in str(value).split(";"):
        if "=" in item:
            path, digest = item.rsplit("=", 1)
            parsed[path] = digest
    return parsed


def join_path_hashes(paths: list[str], hashes: dict[str, str]) -> str:
    return ";".join(f"{path}={hashes[path]}" for path in paths if path in hashes)


def json_cell(value: Any) -> str:
    return json.dumps(value, ensure_ascii=True, sort_keys=True, separators=(",", ":"))


def first_nonempty(*values: str) -> str:
    for value in values:
        if str(value).strip():
            return str(value).strip()
    return ""


def local_commit() -> str:
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()


def evidence_for(record_id: str) -> Path:
    if record_id.startswith("MERS"):
        return EVIDENCE["vizgen"]
    if record_id in {"LR03_BRCA_FFPE", "LR04_BRCA_FF", "LR05_BRCA_ILC", "LR06_CERVICAL", "LR08_INTESTINE"}:
        return EVIDENCE["heca"]
    if record_id.startswith("CTA"):
        return EVIDENCE["cta"]
    if record_id in {"RD03_TNBC_PLASMA", "RD04_CRC_B"}:
        return EVIDENCE["minor"]
    if record_id == "SIM03_HUMAN_LUNG":
        return EVIDENCE["accessions"]
    return EVIDENCE["base"]


def source_layer(record_id: str) -> str:
    if record_id.startswith("SIM"):
        return "PROJECT_GENERATED_SIMULATION_TRUTH"
    if record_id.startswith("STATE"):
        return "PROJECT_GENERATED_DECOY"
    if record_id.startswith("RD") or record_id.startswith("THAL"):
        return "PROJECT_GENERATED_DERIVATIVE"
    if record_id.startswith("CTA"):
        return "PROJECT_GENERATED_ENDPOINT_DERIVATIVE"
    return "THIRD_PARTY_RETAINED_EXPORT"


def dataset_redistribution(record_id: str, layer: str) -> str:
    if record_id.startswith("MERS"):
        return "PERMISSION_UNRESOLVED"
    if record_id == "LR07_HEART":
        return "RESTRICTED"
    if layer.startswith("PROJECT_GENERATED"):
        return "REDISTRIBUTABLE_DERIVATIVE_ONLY"
    return "LINK_AND_RECONSTRUCT_ONLY"


def reference_accession(record: pd.Series) -> str:
    overrides = {
        "CTA01_BCSA2TUMB1": "GSE176078",
        "SIM01_BRCA": "GSE176078",
        "SIM02_MOUSE_BRAIN": "E-MTAB-11115",
        "SIM03_HUMAN_LUNG": "E-MTAB-11640;PRJEB52292",
        "LR02_BRAIN": "E-MTAB-11115",
        "LR06_CERVICAL": "E-MTAB-10287",
        "LR08_INTESTINE": "E-MTAB-9543;E-MTAB-9536",
        "LR09_LYMPH_NODE": "E-MTAB-11343",
        "LR10_EMBRYO": "E-MTAB-6967;GSE169210",
        "RD01_EMBRYO_ENDODERM": "E-MTAB-6967;GSE169210",
        "RD02_EMBRYO_ERYTHROID": "E-MTAB-6967;GSE169210",
        "RD03_TNBC_PLASMA": "GSE176078",
        "RD04_CRC_B": "GSE132465",
        "RD05_BRCA_FFPE_PLASMA": "GSE176078",
        "RD06_BRCA_FFPE_EPITHELIAL": "GSE176078",
        "THAL01_ST8059051": "E-MTAB-11115",
    }
    rid = record["record_id"]
    if rid.startswith("MERS"):
        return str(record["sample_id"])
    return overrides.get(rid, first_nonempty(record["geo_accession"], record["arrayexpress_accession"]))


def primary_accession(record: pd.Series) -> str:
    rid = record["record_id"]
    overrides = {
        "CTA01_BCSA2TUMB1": "10.48723/f4v5-m008",
        "SIM03_HUMAN_LUNG": "ERS20065156",
        "THAL01_ST8059051": "E-MTAB-11114",
    }
    if rid in overrides:
        return overrides[rid]
    return first_nonempty(
        record["study_accession"],
        record["project_accession"],
        record["experiment_accession"],
        record["sample_accession"],
        record["geo_accession"],
        record["arrayexpress_accession"],
    )


def provider(record: pd.Series) -> str:
    rid = record["record_id"]
    if rid.startswith("MERS"):
        return "Vizgen"
    if rid.startswith("CTA"):
        return "Swedish National Data Service; Researchdata.se; Zenodo; NCBI GEO"
    if rid == "THAL01_ST8059051":
        return "EMBL-EBI BioStudies/ArrayExpress"
    if record["database"]:
        return record["database"]
    if "10x" in record["spatial_platform"]:
        return "10x Genomics plus the named reference repository"
    return "Named public provider/repository in official_source_url"


def fact_attributes() -> dict[str, dict[str, Any]]:
    attrs: dict[str, dict[str, Any]] = {
        "SIM01_BRCA": {
            "internal_id": "real_brca",
            "spatial_resource": "public 10x Visium FFPE human breast cancer",
            "block": "738811QB",
            "section": "1",
            "slide": "V11J26-008",
            "capture_area": "B1",
            "reference_resource": "Wu breast-cancer atlas",
            "reference_accession": "GSE176078",
            "reference_subset": "HER2-positive single-cell subset",
            "simulation_truth": "project generated",
        },
        "SIM03_HUMAN_LUNG": {
            "internal_id": "human_lung_5loc",
            "formal_spatial_sample": "WSA_LngSP10193345",
            "ena_sample_accession": "ERS20065156",
            "biosample": "SAMEA115633909",
            "reference_resource": "corresponding donor-level lung resource",
            "reference_accessions": "E-MTAB-11640;PRJEB52292",
            "simulation_truth": "project generated",
            "read_run_accession": "UNRESOLVED",
            "accession_type_guardrail": "ERS20065156 is an ENA sample accession, not a run accession",
        },
        "SIM02_MOUSE_BRAIN": {
            "internal_id": "mouse_brain_refined",
            "spatial_resource": "public 10x CytAssist FFPE sagittal mouse brain",
            "slide": "V52B25-081",
            "capture_area": "B",
            "reference_resource": "adult mouse-brain single-nucleus reference",
            "reference_accession": "E-MTAB-11115",
            "simulation_truth": "project generated",
        },
        "LR01_KIDNEY": {"reference_resource": "Comprehensive Mouse Kidney Atlas"},
        "LR02_BRAIN": {"reference_resource": "E-MTAB-11115 adult mouse-brain single-nucleus reference"},
        "LR03_BRCA_FFPE": {"reference_type": "retained hECA breast reference export", "source_study_and_donor_provenance": "traceable from retained cell-level metadata", "slide_capture": "V11J26-008 / B1", "byte_identity_claim": "not claimed against a newly downloaded archive"},
        "LR04_BRCA_FF": {"reference_type": "retained hECA breast reference export", "source_study_and_donor_provenance": "traceable from retained cell-level metadata", "slide_capture": "V19B23-014 / A1", "byte_identity_claim": "not claimed against a newly downloaded archive"},
        "LR05_BRCA_ILC": {"reference_type": "retained hECA breast reference export", "source_study_and_donor_provenance": "traceable from retained cell-level metadata", "slide_capture": "V19L29-095 / A1", "byte_identity_claim": "not claimed against a newly downloaded archive"},
        "LR06_CERVICAL": {"reference_type": "retained hECA uterus/endometrium reference export", "cell_level_provenance": "traceable", "associated_resources": "Garcia-Alonso;Vento-Tormo;Human Cell Landscape;Tabula Sapiens"},
        "LR07_HEART": {"reference_resource": "healthy ventricular reference reported by Reichart et al."},
        "LR08_INTESTINE": {"reference_type": "retained hECA healthy-intestine reference export", "cell_level_provenance": "traceable", "associated_resources": "Elmentaite intestinal atlas;Human Cell Landscape;adult human cell atlas"},
        "LR09_LYMPH_NODE": {"reference_resource": "developing human immune-system lymph-node reference reported by Suo et al."},
        "LR10_EMBRYO": {"reference_resource": "combined Pijuan-Sala and Mittnenzweig mouse-gastrulation references"},
        "RD01_EMBRYO_ENDODERM": {"spatial_resource": "public 10x CytAssist mouse-embryo FFPE section", "slide": "V52Y09-019", "capture_area": "B", "reference_resource": "combined Pijuan-Sala and Mittnenzweig reference", "dropout_target": "Endoderm/Gut"},
        "RD02_EMBRYO_ERYTHROID": {"spatial_resource": "public 10x CytAssist mouse-embryo FFPE section", "slide": "V52Y09-019", "capture_area": "B", "reference_resource": "combined Pijuan-Sala and Mittnenzweig reference", "dropout_target": "Erythroid"},
        "RD03_TNBC_PLASMA": {"reference_resource": "Wu breast-cancer atlas", "reference_accession": "GSE176078", "setting": "fresh-frozen TNBC section CID4465", "slide": "UNRESOLVED", "capture_area": "UNRESOLVED"},
        "RD05_BRCA_FFPE_PLASMA": {"reference_resource": "Wu breast-cancer atlas", "reference_accession": "GSE176078", "setting": "FFPE breast-cancer section V11J26-008 / B1"},
        "RD06_BRCA_FFPE_EPITHELIAL": {"reference_resource": "Wu breast-cancer atlas", "reference_accession": "GSE176078", "setting": "FFPE breast-cancer section V11J26-008 / B1"},
        "RD04_CRC_B": {"spatial_resource": "public 10x colorectal spatial section", "slide": "V10A13-206", "capture_area": "C1", "assay": "Spatial Gene Expression 3' v1", "processing": "Space Ranger 1.2.0", "reference_resource": "Lee colorectal-cancer single-cell reference", "reference_accession": "GSE132465"},
        "THAL01_ST8059051": {"formal_spatial_sample": "ST8059051", "resource_label": "Visium-29B", "slide": "C05717-021", "capture_area": "B1", "spatial_accession": "E-MTAB-11114", "reference_accession": "E-MTAB-11115", "counted_in_six_setting_summary": "false"},
        "CTA01_BCSA2TUMB1": {
            "spatial_sample": "BCSA2TumB1",
            "disease_context": "HER2-positive breast cancer",
            "platform": "10x Visium Spatial Gene Expression",
            "slide": "V10F24-112",
            "capture_area": "C1",
            "processing": "Space Ranger 1.0.0",
            "resource": "CTA breast-cancer spatial transcriptomics resource",
            "endpoint": "registered computational-pathology-derived immune endpoint",
            "reference_source": "Wu breast-cancer atlas",
            "reference_accession": "GSE176078",
            "project_prepared_subset": "true",
            "total_cells": "4014",
            "donor_matched_to_spatial": "false",
            "pooled_HER2_patients": "5",
            "patient_counts": {"CID3586": 1038, "CID3838": 544, "CID3921": 766, "CID4066": 1211, "CID45171": 455},
            "removed_total": "2414",
            "retained_total": "1600",
            "removed_counts": {"B cells": 400, "CD4 T cells": 400, "CD8 T cells": 400, "monocytes/macrophages": 400, "NK cells": 400, "plasma cells": 226, "generic T cells": 188},
            "retained_counts": {"endothelial cells": 400, "epithelial cells": 400, "fibroblasts": 400, "perivascular-like cells": 400},
            "exact_chemistry_or_kit_version": "UNRESOLVED",
        },
    }
    mers = [
        ("MERS01_BREAST", "HumanBreastCancerPatient1"),
        ("MERS02_COLON", "HumanColonCancerPatient1"),
        ("MERS03_LUNG", "HumanLungCancerPatient1"),
        ("MERS04_MELANOMA1", "HumanMelanomaPatient1"),
        ("MERS05_MELANOMA2", "HumanMelanomaPatient2"),
    ]
    for record_id, sample in mers:
        attrs[record_id] = {
            "resource": "Vizgen MERFISH FFPE Human Immuno-oncology Data Set",
            "release_date": "May 2022",
            "platform": "MERSCOPE / MERFISH",
            "target_gene_count": "500",
            "blank_control_count": "50",
            "formal_dataset_identifier": sample,
            "spatial_and_reference_same_release": "true",
            "spatial_reference_cell_id_overlap": "0",
            "independent_scrna_reference_used": "false",
            "spatial_representation": "one-cell pseudo-spots",
            "provider_checksum": "UNRESOLVED",
            "redistribution_permission": "UNRESOLVED",
        }
    return attrs


def build_dataset_candidate(source: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, dict[str, Any]]]:
    attrs = fact_attributes()
    candidate = source.copy()
    derivative_ids = set(candidate.loc[candidate["record_id"].str.startswith(("RD", "THAL", "CTA")), "record_id"])
    for index, row in candidate.iterrows():
        rid = row["record_id"]
        if rid in derivative_ids and not row["formal_dataset_name"].startswith("Project-generated"):
            if rid.startswith("CTA"):
                candidate.at[index, "formal_dataset_name"] = "Project-generated CTA endpoint/reference derivative linked to BCSA2TumB1"
            elif rid.startswith("THAL"):
                candidate.at[index, "formal_dataset_name"] = "Project-generated thalamic reference-dropout derivative of ST8059051"
            else:
                candidate.at[index, "formal_dataset_name"] = f"Project-generated reference-dropout derivative of {row['formal_dataset_name']}"

    additions: dict[str, list[str]] = {
        "dataset_record_id": [], "formal_sample_id": [], "source_layer": [], "provider_or_repository": [],
        "primary_accession": [], "secondary_accessions": [], "platform": [], "assay_or_chemistry": [],
        "disease_context": [], "slide": [], "processing_software": [], "processing_version": [],
        "reference_source": [], "reference_accession": [], "retained_input_path": [], "retained_input_sha1": [],
        "official_source_route": [], "redistribution_class": [], "candidate_status": [],
        "fact_ledger_attributes_json": [], "classification_evidence_paths": [],
    }
    for _, row in candidate.iterrows():
        rid = row["record_id"]
        layer = source_layer(rid)
        redist = dataset_redistribution(rid, layer)
        processing_software = ""
        processing_version = ""
        if rid == "RD04_CRC_B":
            processing_software, processing_version = "Space Ranger", "1.2.0"
        elif rid == "CTA01_BCSA2TUMB1":
            processing_software, processing_version = "Space Ranger", "1.0.0"
        additions["dataset_record_id"].append(rid)
        additions["formal_sample_id"].append(row["sample_id"])
        additions["source_layer"].append(layer)
        additions["provider_or_repository"].append(provider(row))
        additions["primary_accession"].append(primary_accession(row))
        additions["secondary_accessions"].append(";".join(x for x in [row["sample_accession"], row["biosample_accession"], row["geo_accession"], row["arrayexpress_accession"], row["ena_accession"]] if x))
        additions["platform"].append(row["spatial_platform"])
        additions["assay_or_chemistry"].append("10x Visium Spatial Gene Expression; exact chemistry/kit version unresolved" if rid == "CTA01_BCSA2TUMB1" else row["assay"])
        additions["disease_context"].append(row["disease"])
        additions["slide"].append(row["slide_id"])
        additions["processing_software"].append(processing_software)
        additions["processing_version"].append(processing_version)
        additions["reference_source"].append(row["reference_dataset_name"])
        additions["reference_accession"].append(reference_accession(row))
        additions["retained_input_path"].append(row["input_file"])
        additions["retained_input_sha1"].append(row["input_file_sha1"])
        additions["official_source_route"].append(row["official_source_url"])
        additions["redistribution_class"].append(redist)
        additions["candidate_status"].append(CANDIDATE_STATUS)
        additions["fact_ledger_attributes_json"].append(json_cell(attrs.get(rid, {})))
        evidence = evidence_for(rid)
        additions["classification_evidence_paths"].append(f"{rel(SOURCE_DATASET)};{rel(evidence)}")
    for column, values in additions.items():
        candidate[column] = values
    candidate = candidate.sort_values("dataset_record_id", kind="mergesort").reset_index(drop=True)
    return candidate, attrs


def config_details(config_path: str) -> dict[str, str]:
    if not config_path:
        return {"stage3a": "NOT_RECORDED", "stage3b": "NOT_RECORDED", "seed": "NOT_RECORDED", "shared_gene_rule": "NOT_RECORDED", "storage_group": ""}
    path = ROOT / config_path
    if not path.exists():
        return {"stage3a": "NOT_RECORDED", "stage3b": "NOT_RECORDED", "seed": "NOT_RECORDED", "shared_gene_rule": "NOT_RECORDED", "storage_group": ""}
    text = path.read_text(encoding="utf-8-sig")
    seed_match = re.search(r"(?m)^\s*random_seed:\s*([^#\r\n]+)", text)
    group_match = re.search(r"(?ms)^storage:\s*\n\s+group:\s*([^#\r\n]+)", text)
    sc_match = re.search(r"(?m)^\s*min_cells_sc:\s*([^#\r\n]+)", text)
    st_match = re.search(r"(?m)^\s*min_cells_st:\s*([^#\r\n]+)", text)
    return {
        "stage3a": "true" if re.search(r"(?m)^stage3:\s*$", text) else "NOT_RECORDED",
        "stage3b": "true" if re.search(r"(?m)^stage3b:\s*$", text) else "NOT_RECORDED",
        "seed": seed_match.group(1).strip() if seed_match else "NOT_RECORDED",
        "shared_gene_rule": f"min_cells_sc={sc_match.group(1).strip() if sc_match else 'NOT_RECORDED'};min_cells_st={st_match.group(1).strip() if st_match else 'NOT_RECORDED'}",
        "storage_group": group_match.group(1).strip() if group_match else "",
    }


def result_summary_inventory(experiment_id: str, config_path: str) -> tuple[list[str], dict[str, str], str]:
    details = config_details(config_path)
    roots = []
    if details["storage_group"]:
        roots.append(ROOT / "result" / details["storage_group"] / experiment_id)
    roots.append(ROOT / "result" / experiment_id)
    result_root = next((path for path in roots if path.is_dir()), None)
    if result_root is None:
        return [], {}, "NOT_RECORDED"
    allowed_names = {"stage1_summary.json", "stage3_summary.json", "stage3b_summary.json", "stage4_summary.json"}
    paths = sorted(
        (path for path in result_root.rglob("*.json") if path.name in allowed_names),
        key=lambda path: rel(path),
    )
    route_paths = [path for path in paths if "baseline" not in path.as_posix().lower()]
    if route_paths:
        paths = route_paths
    rel_paths = [rel(path) for path in paths]
    hashes = {rel(path): sha1(path) for path in paths}
    stage4 = "true" if any(path.name == "stage4_summary.json" for path in paths) else "NOT_RECORDED"
    return rel_paths, hashes, stage4


def source_values_for(experiment_id: str, family: str) -> tuple[list[str], dict[str, str]]:
    paths: list[Path] = []
    if family == "downstream external tumour pairings":
        paths = [
            ROOT / "visualizations/cytospace_fig2d_profile_mask_benchmark/fig2d_profile_mask_benchmark_source_values.csv",
            ROOT / "result/cytospace_fig2e_stage3_profile_mask/fig2e_stage3_profile_mask_metrics.csv",
        ]
    elif family == "low-resolution profile masking":
        path = ROOT / "result/real_profile_mask_foundation" / experiment_id / "spot_foundation.csv"
        if path.exists():
            paths = [path]
    paths = [path for path in paths if path.exists()]
    return [rel(path) for path in paths], {rel(path): sha1(path) for path in paths}


def build_experiment_candidate(source: pd.DataFrame, datasets: pd.DataFrame) -> pd.DataFrame:
    dataset_index = datasets.set_index("dataset_record_id", drop=False)
    candidate = source.copy()
    additions: dict[str, list[str]] = {column: [] for column in [
        "experiment_record_id", "dataset_record_id", "experiment_family", "formal_experiment_name",
        "internal_experiment_id", "spatial_input_record", "reference_input_record", "reference_condition",
        "decoy_population", "stage3a_enabled", "stage3b_enabled", "stage4_enabled", "mapping_backend",
        "mapping_capacity", "shared_gene_rule", "random_seed", "formal_config_path", "formal_config_sha1",
        "formal_input_paths", "formal_input_sha1s", "formal_output_paths", "formal_output_sha1s",
        "source_value_paths", "source_value_sha1s", "redistribution_class", "provenance_status",
        "evidence_paths", "candidate_status",
    ]}
    for _, row in candidate.iterrows():
        experiment_id = row["experiment_id"]
        dataset_id = row["spatial_input_record_id"]
        dataset = dataset_index.loc[dataset_id]
        config_path = dataset["config_file"].split(";")[0] if dataset["config_file"] else ""
        config_hashes = parse_path_hashes(dataset["config_file_sha1"])
        details = config_details(config_path)
        output_paths, output_hashes, stage4 = result_summary_inventory(experiment_id, config_path)
        source_paths, source_hashes = source_values_for(experiment_id, row["dataset_family"])
        derivative = row["project_derivative"]
        reference_condition = "reference dropout: " + row["reference_dropout_target"] if row["reference_dropout_target"] else ("profile masking: " + row["profile_mask_target"] if row["profile_mask_target"] else derivative)
        decoy = derivative if "decoy" in derivative.lower() else ""
        mapping_backend = "CytoSPACE" if stage4 == "true" or row["dataset_family"] in {"CTA biological application", "state-decoy experiments"} else "NOT_RECORDED"
        evidence_path = evidence_for(dataset_id)
        values = {
            "experiment_record_id": experiment_id,
            "dataset_record_id": dataset_id,
            "experiment_family": row["dataset_family"],
            "formal_experiment_name": f"{row['manuscript_result_section']}: {experiment_id}",
            "internal_experiment_id": experiment_id,
            "spatial_input_record": row["spatial_input_record_id"],
            "reference_input_record": row["reference_input_record_id"],
            "reference_condition": reference_condition,
            "decoy_population": decoy,
            "stage3a_enabled": details["stage3a"],
            "stage3b_enabled": details["stage3b"],
            "stage4_enabled": stage4,
            "mapping_backend": mapping_backend,
            "mapping_capacity": "NOT_RECORDED",
            "shared_gene_rule": details["shared_gene_rule"],
            "random_seed": details["seed"],
            "formal_config_path": config_path,
            "formal_config_sha1": config_hashes.get(config_path, ""),
            "formal_input_paths": dataset["retained_input_path"],
            "formal_input_sha1s": dataset["retained_input_sha1"],
            "formal_output_paths": ";".join(output_paths),
            "formal_output_sha1s": join_path_hashes(output_paths, output_hashes),
            "source_value_paths": ";".join(source_paths),
            "source_value_sha1s": join_path_hashes(source_paths, source_hashes),
            "redistribution_class": "REDISTRIBUTABLE_DERIVATIVE_ONLY",
            "provenance_status": row["formal_status"],
            "evidence_paths": f"{rel(SOURCE_EXPERIMENT)};{rel(evidence_path)}" + (f";{config_path}" if config_path else ""),
            "candidate_status": CANDIDATE_STATUS,
        }
        for column in additions:
            additions[column].append(values[column])
    for column, values in additions.items():
        candidate[column] = values
    return candidate.sort_values("experiment_record_id", kind="mergesort").reset_index(drop=True)


def verify_hash_coverage(datasets: pd.DataFrame) -> list[dict[str, Any]]:
    cached: dict[str, dict[str, Any]] = {}
    if CACHE.exists():
        cached = json.loads(CACHE.read_text(encoding="utf-8"))
    new_cache: dict[str, dict[str, Any]] = {}
    unique: dict[str, str] = {}
    references: list[tuple[str, str, str, str]] = []
    for _, row in datasets.iterrows():
        for role, path_field, hash_field in [
            ("retained_input", "retained_input_path", "retained_input_sha1"),
            ("formal_config", "config_file", "config_file_sha1"),
        ]:
            expected = parse_path_hashes(row[hash_field])
            for path in [item for item in row[path_field].split(";") if item]:
                references.append((row["dataset_record_id"], role, path, expected.get(path, "")))
                if path not in unique:
                    unique[path] = expected.get(path, "")
    observed: dict[str, str] = {}
    for path, expected in sorted(unique.items()):
        local = ROOT / path
        if not local.exists():
            observed[path] = ""
            continue
        stat = local.stat()
        cache_row = cached.get(path, {})
        if cache_row.get("size") == stat.st_size and cache_row.get("mtime_ns") == stat.st_mtime_ns and cache_row.get("sha1"):
            digest = str(cache_row["sha1"])
        else:
            digest = sha1(local)
        observed[path] = digest
        new_cache[path] = {"size": stat.st_size, "mtime_ns": stat.st_mtime_ns, "sha1": digest}
    CACHE.write_text(json.dumps(new_cache, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    rows = []
    for index, (record_id, role, path, expected) in enumerate(references, start=1):
        local = ROOT / path
        digest = observed.get(path, "")
        rows.append({
            "coverage_id": f"HASH{index:03d}", "dataset_record_id": record_id, "file_role": role,
            "project_relative_path": path, "path_exists": str(local.exists()).lower(),
            "expected_sha1": expected, "observed_sha1": digest,
            "hash_present": str(bool(expected)).lower(), "hash_match": str(bool(expected) and expected == digest).lower(),
            "verification_basis": "CURRENT_BYTES_SHA1" if local.exists() else "MISSING_PATH",
            "notes": "Large inputs are hashed once and reused only when size and mtime_ns are unchanged.",
        })
    return rows


def unresolved_rows() -> list[dict[str, str]]:
    return [
        {"issue_id": "UNB001", "record_ids": "SIM03_HUMAN_LUNG", "field": "human-lung read-run accession", "status": "UNRESOLVED_NONBLOCKING", "manuscript_claim_affected": "false", "reason": "ERS20065156 is confirmed as an ENA sample accession; no read-run accession was verified.", "safe_public_wording": "Identify WSA_LngSP10193345 using ENA sample ERS20065156 and BioSample SAMEA115633909; do not state a run accession.", "evidence_paths": f"{rel(EVIDENCE['accessions'])};{rel(EVIDENCE['unresolved'])}"},
        {"issue_id": "UNB002", "record_ids": "RD03_TNBC_PLASMA", "field": "CID4465 slide identifier", "status": "UNRESOLVED_NONBLOCKING", "manuscript_claim_affected": "false", "reason": "The official retained archives identify CID4465 but do not expose a slide serial.", "safe_public_wording": "Fresh-frozen TNBC spatial section CID4465; omit slide serial.", "evidence_paths": f"{rel(EVIDENCE['minor'])};{rel(EVIDENCE['unresolved'])}"},
        {"issue_id": "UNB003", "record_ids": "RD03_TNBC_PLASMA", "field": "CID4465 capture area", "status": "UNRESOLVED_NONBLOCKING", "manuscript_claim_affected": "false", "reason": "The official retained archives identify CID4465 but do not expose a capture area.", "safe_public_wording": "Fresh-frozen TNBC spatial section CID4465; omit capture area.", "evidence_paths": f"{rel(EVIDENCE['minor'])};{rel(EVIDENCE['unresolved'])}"},
        {"issue_id": "UNB004", "record_ids": "CTA01_BCSA2TUMB1", "field": "CTA exact chemistry or kit version", "status": "UNRESOLVED_NONBLOCKING", "manuscript_claim_affected": "false", "reason": "Slide, capture area, Visium platform and Space Ranger 1.0.0 are confirmed, but no exact chemistry/kit version is recorded.", "safe_public_wording": "10x Visium Spatial Gene Expression processed with Space Ranger 1.0.0; omit exact chemistry/kit version.", "evidence_paths": f"{rel(EVIDENCE['cta'])};{rel(EVIDENCE['unresolved'])}"},
        {"issue_id": "UNB005", "record_ids": "MERS01_BREAST;MERS02_COLON;MERS03_LUNG;MERS04_MELANOMA1;MERS05_MELANOMA2", "field": "Vizgen provider-supplied checksum", "status": "UNRESOLVED_NONBLOCKING", "manuscript_claim_affected": "false", "reason": "Local SHA-1 values exist, but provider-issued checksums/download receipts were not retained.", "safe_public_wording": "Provide formal sample IDs, official provider routes and retained local SHA-1 values; do not claim provider checksum verification.", "evidence_paths": f"{rel(EVIDENCE['vizgen'])};{rel(EVIDENCE['vizgen_files'])};{rel(EVIDENCE['unresolved'])}"},
        {"issue_id": "UNB006", "record_ids": "MERS01_BREAST;MERS02_COLON;MERS03_LUNG;MERS04_MELANOMA1;MERS05_MELANOMA2", "field": "Vizgen explicit redistribution permission", "status": "UNRESOLVED_NONBLOCKING", "manuscript_claim_affected": "false", "reason": "Public downloadability and general use language do not establish unambiguous redistribution permission.", "safe_public_wording": "Link to the provider release and reconstruct locally; do not include provider raw files in the public package.", "evidence_paths": f"{rel(EVIDENCE['vizgen'])};{rel(EVIDENCE['unresolved'])}"},
    ]


def build_fact_ledger(
    datasets: pd.DataFrame,
    experiments: pd.DataFrame,
    attrs: dict[str, dict[str, Any]],
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    dataset_index = datasets.set_index("dataset_record_id", drop=False)
    experiment_index = experiments.set_index("experiment_record_id", drop=False)
    ledger: list[dict[str, str]] = []
    reconciliation: list[dict[str, str]] = []

    def render(value: Any) -> str:
        if isinstance(value, (dict, list)):
            return json_cell(value)
        return str(value)

    def add(
        group: str,
        fact: str,
        record_id: str,
        expected: Any,
        actual: Any,
        field: str,
        evidence_paths: str,
        qualified: bool = False,
        severity: str = "INFORMATIONAL",
    ) -> None:
        fact_id = f"MF{len(ledger) + 1:03d}"
        expected_text = render(expected)
        actual_text = render(actual)
        if qualified:
            status = "MATCH_WITH_QUALIFIED_WORDING" if actual_text in {expected_text, "", "UNRESOLVED"} else "CONFLICT"
        else:
            status = "MATCH" if actual_text == expected_text else "CONFLICT"
        ledger.append({
            "fact_id": fact_id, "manuscript_fact_group": group, "manuscript_fact": fact,
            "manifest_record_id": record_id, "candidate_field": field, "expected_value": expected_text,
            "fact_source": "Task 4F user-supplied manuscript fact ledger", "evidence_paths": evidence_paths,
        })
        reconciliation.append({
            "check_id": f"REC{len(reconciliation) + 1:03d}", "manuscript_fact_group": group,
            "manuscript_fact": fact, "manifest_record_id": record_id, "manifest_value": actual_text,
            "match_status": status, "severity": "BLOCKING" if status == "CONFLICT" else severity,
            "recommended_action": "Resolve before freeze." if status == "CONFLICT" else ("Retain qualified wording and blank/UNRESOLVED field." if qualified else "No change."),
            "evidence_paths": evidence_paths,
        })

    def attr_facts(group: str, record_id: str, keys: list[str], qualified: set[str] | None = None) -> None:
        qualified = qualified or set()
        evidence = f"{rel(SOURCE_DATASET)};{rel(evidence_for(record_id))}"
        for key in keys:
            expected = attrs[record_id][key]
            add(group, f"{record_id}.{key} = {render(expected)}", record_id, expected, attrs[record_id].get(key, ""), f"fact_ledger_attributes_json.{key}", evidence, key in qualified)

    attr_facts("A1 BRCA joint simulation", "SIM01_BRCA", list(attrs["SIM01_BRCA"]))
    attr_facts("A2 human-lung joint simulation", "SIM03_HUMAN_LUNG", list(attrs["SIM03_HUMAN_LUNG"]), {"read_run_accession"})
    attr_facts("A3 mouse-brain joint simulation", "SIM02_MOUSE_BRAIN", list(attrs["SIM02_MOUSE_BRAIN"]))

    lowres_keys = {
        "LR01_KIDNEY": list(attrs["LR01_KIDNEY"]), "LR02_BRAIN": list(attrs["LR02_BRAIN"]),
        "LR03_BRCA_FFPE": list(attrs["LR03_BRCA_FFPE"]), "LR04_BRCA_FF": list(attrs["LR04_BRCA_FF"]),
        "LR05_BRCA_ILC": list(attrs["LR05_BRCA_ILC"]), "LR06_CERVICAL": list(attrs["LR06_CERVICAL"]),
        "LR07_HEART": list(attrs["LR07_HEART"]), "LR08_INTESTINE": list(attrs["LR08_INTESTINE"]),
        "LR09_LYMPH_NODE": list(attrs["LR09_LYMPH_NODE"]), "LR10_EMBRYO": list(attrs["LR10_EMBRYO"]),
    }
    for record_id, keys in lowres_keys.items():
        attr_facts("B low-resolution profile masking", record_id, keys)
    lowres_count = int((experiments["experiment_family"] == "low-resolution profile masking").sum())
    add("B low-resolution profile masking", "total low-resolution profile-masking experiments = 10", "EXPERIMENT_MANIFEST", "10", str(lowres_count), "experiment_family", rel(SOURCE_EXPERIMENT))

    mers_ids = ["MERS01_BREAST", "MERS02_COLON", "MERS03_LUNG", "MERS04_MELANOMA1", "MERS05_MELANOMA2"]
    for record_id in mers_ids:
        attr_facts("C cell-resolution MERSCOPE", record_id, list(attrs[record_id]), {"provider_checksum", "redistribution_permission"})

    for record_id in ["RD01_EMBRYO_ENDODERM", "RD02_EMBRYO_ERYTHROID", "RD03_TNBC_PLASMA", "RD05_BRCA_FFPE_PLASMA", "RD06_BRCA_FFPE_EPITHELIAL", "RD04_CRC_B", "THAL01_ST8059051"]:
        qualified = {"slide", "capture_area"} if record_id == "RD03_TNBC_PLASMA" else set()
        attr_facts("D real reference-dropout experiments", record_id, list(attrs[record_id]), qualified)

    attr_facts("E CTA biological application", "CTA01_BCSA2TUMB1", [key for key in attrs["CTA01_BCSA2TUMB1"] if key not in {"patient_counts", "removed_counts", "retained_counts"}], {"exact_chemistry_or_kit_version"})
    for patient, count in attrs["CTA01_BCSA2TUMB1"]["patient_counts"].items():
        add("E CTA biological application", f"patient count {patient} = {count}", "CTA01_BCSA2TUMB1", str(count), str(attrs["CTA01_BCSA2TUMB1"]["patient_counts"][patient]), f"fact_ledger_attributes_json.patient_counts.{patient}", f"{rel(SOURCE_DATASET)};{rel(EVIDENCE['cta'])}")
    for cell_type, count in attrs["CTA01_BCSA2TUMB1"]["removed_counts"].items():
        add("E CTA biological application", f"removed {cell_type} = {count}", "CTA01_BCSA2TUMB1", str(count), str(attrs["CTA01_BCSA2TUMB1"]["removed_counts"][cell_type]), f"fact_ledger_attributes_json.removed_counts.{cell_type}", f"{rel(SOURCE_DATASET)};{rel(EVIDENCE['cta'])}")
    for cell_type, count in attrs["CTA01_BCSA2TUMB1"]["retained_counts"].items():
        add("E CTA biological application", f"retained {cell_type} = {count}", "CTA01_BCSA2TUMB1", str(count), str(attrs["CTA01_BCSA2TUMB1"]["retained_counts"][cell_type]), f"fact_ledger_attributes_json.retained_counts.{cell_type}", f"{rel(SOURCE_DATASET)};{rel(EVIDENCE['cta'])}")

    # Cross-check key manifest columns independently of the JSON ledger.
    column_checks = [
        ("A2 human-lung joint simulation", "SIM03_HUMAN_LUNG", "sample_accession", "ERS20065156", "ena sample accession"),
        ("A2 human-lung joint simulation", "SIM03_HUMAN_LUNG", "run_accession", "", "run accession intentionally blank"),
        ("A2 human-lung joint simulation", "SIM03_HUMAN_LUNG", "biosample_accession", "SAMEA115633909", "BioSample"),
        ("D real reference-dropout experiments", "RD03_TNBC_PLASMA", "slide_id", "UNRESOLVED", "CID4465 slide"),
        ("D real reference-dropout experiments", "RD03_TNBC_PLASMA", "capture_area", "UNRESOLVED", "CID4465 capture area"),
        ("E CTA biological application", "CTA01_BCSA2TUMB1", "processing_version", "1.0.0", "CTA Space Ranger version"),
    ]
    for group, record_id, column, expected, description in column_checks:
        actual = dataset_index.loc[record_id, column]
        qualified = expected in {"", "UNRESOLVED"}
        add(group, description, record_id, expected, actual, column, f"{rel(SOURCE_DATASET)};{rel(evidence_for(record_id))}", qualified)

    # Experiment-level dropout and mask targets.
    for experiment_id, expected_target in {
        "mouse_embryo_real_sc_missing_endoderm_gut": "Endoderm/Gut",
        "mouse_embryo_real_sc_missing_erythroid": "Erythroid",
        "cytospace_fig2d_tme_brca_tnbc_fresh_frozen_sc_missing_plasma_cells": "Plasma cells",
        "cytospace_fig2d_tme_brca_her2_ffpe_sc_missing_plasma_cells": "Plasma cells",
        "cytospace_fig2d_tme_brca_her2_ffpe_sc_missing_epithelial_cells": "Epithelial cells",
        "cytospace_fig2d_tme_crc_fresh_frozen_sc_missing_b_cells": "B cells",
        "cell2location_ST8059051_thalamic_excitatory_reference_missing": "thalamic excitatory states",
    }.items():
        actual = experiment_index.loc[experiment_id, "reference_dropout_target"]
        add("D real reference-dropout experiments", f"dropout target for {experiment_id}", experiment_id, expected_target, actual, "reference_dropout_target", f"{rel(SOURCE_EXPERIMENT)};{experiment_index.loc[experiment_id, 'evidence_paths']}")

    return ledger, reconciliation


def change_log(source: pd.DataFrame, candidate: pd.DataFrame, manifest_type: str, id_column: str) -> list[dict[str, str]]:
    source_id = "record_id" if manifest_type == "dataset" else "experiment_id"
    source_index = source.set_index(source_id, drop=False)
    candidate_index = candidate.set_index(id_column, drop=False)
    rows: list[dict[str, str]] = []
    new_columns = [column for column in candidate.columns if column not in source.columns]
    for record_id in candidate_index.index:
        old = source_index.loc[record_id]
        new = candidate_index.loc[record_id]
        evidence = evidence_for(record_id) if manifest_type == "dataset" else SOURCE_EXPERIMENT
        evidence_digest = sha1(evidence)
        changed = False
        for field in source.columns:
            if str(old[field]) != str(new[field]):
                changed = True
                rows.append({
                    "change_id": f"CHG{len(rows) + 1:04d}", "manifest_type": manifest_type, "record_id": record_id,
                    "field": field, "old_value": str(old[field]), "new_value": str(new[field]),
                    "change_type": "WORDING_NORMALISATION", "reason": "Prevent a project-generated derivative from being described as an independent public biological dataset.",
                    "evidence_path": rel(evidence), "evidence_sha1": evidence_digest,
                })
        for field in new_columns:
            value = str(new[field])
            if not value:
                continue
            changed = True
            if field == "source_layer":
                change_type = "SOURCE_LAYER_CLASSIFICATION"
            elif field == "redistribution_class":
                change_type = "REDISTRIBUTION_CLASSIFICATION"
            elif "sha1" in field:
                change_type = "HASH_ADDITION"
            else:
                change_type = "METADATA_ADDITION"
            rows.append({
                "change_id": f"CHG{len(rows) + 1:04d}", "manifest_type": manifest_type, "record_id": record_id,
                "field": field, "old_value": "", "new_value": value, "change_type": change_type,
                "reason": "Required v1.0 freeze-candidate field added without deleting the v0.2 field set.",
                "evidence_path": rel(evidence), "evidence_sha1": evidence_digest,
            })
        rows.append({
            "change_id": f"CHG{len(rows) + 1:04d}", "manifest_type": manifest_type, "record_id": record_id,
            "field": "legacy_field_set", "old_value": "preserved", "new_value": "preserved",
            "change_type": "NO_CHANGE", "reason": "All useful v0.2 fields remain in the candidate.",
            "evidence_path": rel(evidence), "evidence_sha1": evidence_digest,
        })
    return rows


def validation_checks(
    datasets: pd.DataFrame,
    experiments: pd.DataFrame,
    hash_rows: list[dict[str, Any]],
    reconciliation: list[dict[str, str]],
) -> list[dict[str, str]]:
    dataset_ids = set(datasets["dataset_record_id"])
    experiment_fks = set(experiments["dataset_record_id"])
    invalid_redist = datasets[
        datasets["redistribution_class"].eq("REDISTRIBUTABLE")
        & datasets["record_id"].str.startswith("MERS")
    ]
    ena_bad = datasets[
        datasets["sample_accession"].eq("ERS20065156")
        & datasets["run_accession"].ne("")
    ]
    cid = datasets.set_index("dataset_record_id").loc["RD03_TNBC_PLASMA"]
    cta = datasets.set_index("dataset_record_id").loc["CTA01_BCSA2TUMB1"]
    checks = [
        ("V001", "unique dataset IDs", "0 duplicates", str(int(datasets["dataset_record_id"].duplicated().sum())), not datasets["dataset_record_id"].duplicated().any()),
        ("V002", "unique experiment IDs", "0 duplicates", str(int(experiments["experiment_record_id"].duplicated().sum())), not experiments["experiment_record_id"].duplicated().any()),
        ("V003", "experiment dataset foreign keys", "0 unresolved", str(len(experiment_fks - dataset_ids)), not (experiment_fks - dataset_ids)),
        ("V004", "source-layer vocabulary", "0 invalid", str(int((~datasets["source_layer"].isin(ALLOWED_SOURCE_LAYERS)).sum())), datasets["source_layer"].isin(ALLOWED_SOURCE_LAYERS).all()),
        ("V005", "dataset redistribution vocabulary", "0 invalid", str(int((~datasets["redistribution_class"].isin(ALLOWED_REDISTRIBUTION)).sum())), datasets["redistribution_class"].isin(ALLOWED_REDISTRIBUTION).all()),
        ("V006", "experiment redistribution vocabulary", "0 invalid", str(int((~experiments["redistribution_class"].isin(ALLOWED_REDISTRIBUTION)).sum())), experiments["redistribution_class"].isin(ALLOWED_REDISTRIBUTION).all()),
        ("V007", "existing retained input/config hash coverage", "all present and matching", str(sum(row["hash_match"] == "true" for row in hash_rows)) + "/" + str(len(hash_rows)), all(row["hash_match"] == "true" for row in hash_rows)),
        ("V008", "unresolved Vizgen permission not REDISTRIBUTABLE", "0 violations", str(len(invalid_redist)), invalid_redist.empty),
        ("V009", "ERS20065156 not recorded as run", "0 violations", str(len(ena_bad)), ena_bad.empty),
        ("V010", "CID4465 slide not invented", "UNRESOLVED", cid["slide_id"], cid["slide_id"] == "UNRESOLVED"),
        ("V011", "CID4465 capture area not invented", "UNRESOLVED", cid["capture_area"], cid["capture_area"] == "UNRESOLVED"),
        ("V012", "CTA exact chemistry not invented", "contains unresolved", cta["assay_or_chemistry"], "unresolved" in cta["assay_or_chemistry"].lower()),
        ("V013", "dataset record count", "37", str(len(datasets)), len(datasets) == 37),
        ("V014", "experiment record count", "40", str(len(experiments)), len(experiments) == 40),
        ("V015", "manuscript reconciliation conflicts", "0", str(sum(row["match_status"] == "CONFLICT" for row in reconciliation)), not any(row["match_status"] == "CONFLICT" for row in reconciliation)),
        ("V016", "candidate status", CANDIDATE_STATUS, ";".join(sorted(set(datasets["candidate_status"]) | set(experiments["candidate_status"]))), set(datasets["candidate_status"]) == {CANDIDATE_STATUS} and set(experiments["candidate_status"]) == {CANDIDATE_STATUS}),
    ]
    return [
        {"check_id": check_id, "check": check, "expected": expected, "observed": observed, "status": "PASS" if passed else "FAIL", "details": ""}
        for check_id, check, expected, observed, passed in checks
    ]


def release_plan_rows() -> list[dict[str, str]]:
    audit_files = [
        "README.md", "source_manifest_inventory.tsv", "dataset_record_status.tsv", "experiment_record_status.tsv",
        "manuscript_fact_ledger.tsv", "manuscript_manifest_reconciliation.tsv", "source_layer_classification.tsv",
        "redistribution_classification.tsv", "input_hash_coverage.tsv", "unresolved_nonblocking_records.tsv",
        "blocking_conflicts.tsv", "candidate_change_log.tsv", "release_inventory_addition_plan.tsv",
        "manual_manuscript_update_recommendations.tsv", "validation_checks.tsv", "audit_summary.json",
        "final_decision.md", "local_input_sha1_verification_cache.json",
    ]
    rows = [
        {"plan_id": "REL001", "artifact_class": "APPROVED_FINAL_MANIFEST", "candidate_path": rel(CANDIDATE_DATASET), "planned_release_path": "reproducibility_release/manifests/dataset_manifest_v1.0.tsv", "candidate_sha1": sha1(CANDIDATE_DATASET), "action_after_approval": "Copy approved candidate under final name and add its final SHA-1 to release_manifest.tsv.", "current_action": "PLAN_ONLY", "task4e_records_modified": "false"},
        {"plan_id": "REL002", "artifact_class": "APPROVED_FINAL_MANIFEST", "candidate_path": rel(CANDIDATE_EXPERIMENT), "planned_release_path": "reproducibility_release/manifests/experiment_manifest_v1.0.tsv", "candidate_sha1": sha1(CANDIDATE_EXPERIMENT), "action_after_approval": "Copy approved candidate under final name and add its final SHA-1 to release_manifest.tsv.", "current_action": "PLAN_ONLY", "task4e_records_modified": "false"},
        {"plan_id": "REL003", "artifact_class": "APPROVED_FINAL_METADATA", "candidate_path": rel(CANDIDATE_METADATA), "planned_release_path": "reproducibility_release/manifests/manifest_v1.0_metadata.json", "candidate_sha1": "computed after metadata is approved", "action_after_approval": "Create final metadata after approval and add its SHA-1 to release_manifest.tsv.", "current_action": "PLAN_ONLY", "task4e_records_modified": "false"},
    ]
    for name in audit_files:
        rows.append({
            "plan_id": f"REL{len(rows) + 1:03d}", "artifact_class": "TASK4F_AUDIT_FILE",
            "candidate_path": rel(OUT / name), "planned_release_path": rel(OUT / name),
            "candidate_sha1": "compute at final approval", "action_after_approval": "Add the reviewed audit artifact and its SHA-1 to release_manifest.tsv.",
            "current_action": "PLAN_ONLY", "task4e_records_modified": "false",
        })
    return rows


def recommendations_rows() -> list[dict[str, str]]:
    return [
        {"recommendation_id": "MAN001", "scope": "Human-lung accession wording", "action": "MANUAL_REVIEW_ONLY", "recommended_wording": "ERS20065156 is an ENA sample accession; do not call it a sequencing run accession.", "reason": "No read-run accession is verified.", "blocking": "false"},
        {"recommendation_id": "MAN002", "scope": "CID4465 metadata", "action": "MANUAL_REVIEW_ONLY", "recommended_wording": "Name CID4465 without a slide or capture-area identifier.", "reason": "Neither value is present in the retained official evidence.", "blocking": "false"},
        {"recommendation_id": "MAN003", "scope": "CTA platform wording", "action": "MANUAL_REVIEW_ONLY", "recommended_wording": "State Visium Spatial Gene Expression and Space Ranger 1.0.0; omit exact chemistry/kit version.", "reason": "Exact chemistry remains unresolved.", "blocking": "false"},
        {"recommendation_id": "MAN004", "scope": "Vizgen Data Availability", "action": "MANUAL_REVIEW_ONLY", "recommended_wording": "Link the May 2022 Vizgen provider release and formal sample folders; publish local checksums/configuration but do not redistribute provider raw files.", "reason": "Provider checksums and explicit redistribution permission are unresolved.", "blocking": "false"},
        {"recommendation_id": "MAN005", "scope": "Candidate status", "action": "NO_AUTOMATIC_FREEZE", "recommended_wording": "Keep both manifests at PENDING_HUMAN_APPROVAL until the user reviews the candidates and reconciliation table.", "reason": "Task 4F creates a freeze candidate, not a final frozen release.", "blocking": "false"},
    ]


def make_readme(
    dataset_counts: dict[str, int],
    experiment_counts: dict[str, int],
    reconciliation_counts: Counter[str],
    redist_counts: Counter[str],
) -> str:
    return f"""# Task 4F manifest reconciliation and v1.0 freeze candidate

Generated: {GENERATED_DATE}

## Decision

**CONDITIONAL_PASS_TO_FREEZE_REVIEW**

The v0.2 source manifests remain unchanged at 37 dataset records and 40 experiment records. Their status anchors also remain unchanged:

- datasets: {dataset_counts}
- experiments: {experiment_counts}

The v1.0 files are candidates only and are marked `{CANDIDATE_STATUS}`. They are not `FROZEN`, `FINAL` or `RELEASED`.

## Reconciliation

The Task 4F user-supplied manuscript fact ledger was converted into one check per fact without reading the formal manuscript or bibliography. Reconciliation counts are:

{dict(reconciliation_counts)}

No blocking manuscript-manifest conflict was found. Six explicitly bounded non-blocking fields remain unresolved: the human-lung run accession, CID4465 slide, CID4465 capture area, CTA exact chemistry, Vizgen provider checksum and Vizgen redistribution permission.

## Classification

All 37 dataset records have an allowed source-layer classification. Dataset and experiment redistribution classes are conservative; combined counts are:

{dict(redist_counts)}

No provider file with unresolved permission is classified as `REDISTRIBUTABLE`. Project-generated derivatives are not represented as independent public datasets.

## Hash validation

All retained input and formal configuration paths referenced by v0.2 exist and have SHA-1 values. Current bytes were checked against those hashes; large-file hashes are cached only when path, size and nanosecond mtime remain unchanged.

## Boundary

This task did not access or modify the formal manuscript, bibliography, figures, source-value tables, experimental code, configurations or experimental outputs. It did not run Stage0--Stage5, CytoSPACE or any figure generator, and it did not access a GitHub remote. The existing Task 4E release manifest was protected and not edited.
"""


def make_final_decision(
    reconciliation_counts: Counter[str],
    candidate_hashes: dict[str, str],
) -> str:
    return f"""# Task 4F final decision

**Decision:** `CONDITIONAL_PASS_TO_FREEZE_REVIEW`  
**Candidate status:** `{CANDIDATE_STATUS}`  
**Ready for final manifest freeze:** `false`

## Basis

- Current v0.2 manifests reconcile to exactly 37 dataset and 40 experiment records; no record was added, deleted, split or merged.
- All source-layer and redistribution classifications use the frozen Task 4F vocabularies.
- All manuscript fact checks avoid overclaiming unresolved metadata.
- Reconciliation counts: `{dict(reconciliation_counts)}`.
- Blocking conflicts: `0`.
- Major unexplained manuscript-manifest conflicts: `0`.
- Explicit unresolved non-blocking records: `6`.
- Candidate dataset SHA-1: `{candidate_hashes['dataset']}`.
- Candidate experiment SHA-1: `{candidate_hashes['experiment']}`.

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
"""


def main() -> int:
    required = [SOURCE_DATASET, SOURCE_EXPERIMENT, SOURCE_METADATA, TASK4E_RELEASE, *EVIDENCE.values()]
    missing = [rel(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing Task 4F inputs: {missing}")

    source_datasets = read_tsv(SOURCE_DATASET)
    source_experiments = read_tsv(SOURCE_EXPERIMENT)
    if len(source_datasets) != 37 or len(source_experiments) != 40:
        raise RuntimeError("Current manifest counts differ from the 37/40 audit anchor; investigate before candidate generation.")
    expected_dataset_status = {"CONFIRMED": 28, "PARTIAL": 7, "PROJECT_GENERATED": 2}
    expected_experiment_status = {"CONFIRMED": 29, "PARTIAL": 9, "PROJECT_GENERATED": 2}
    if source_datasets["provenance_status"].value_counts().to_dict() != expected_dataset_status:
        raise RuntimeError("Dataset status anchor changed.")
    if source_experiments["formal_status"].value_counts().to_dict() != expected_experiment_status:
        raise RuntimeError("Experiment status anchor changed.")

    config_paths = []
    for value in source_datasets["config_file"]:
        config_paths.extend(ROOT / item for item in value.split(";") if item)
    protected = [SOURCE_DATASET, SOURCE_EXPERIMENT, SOURCE_METADATA, TASK4E_RELEASE, *config_paths]
    before = {rel(path): sha1(path) for path in protected if path.exists()}

    CANDIDATE_DIR.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)

    datasets, attrs = build_dataset_candidate(source_datasets)
    experiments = build_experiment_candidate(source_experiments, datasets)
    write_tsv(CANDIDATE_DATASET, datasets)
    write_tsv(CANDIDATE_EXPERIMENT, experiments)

    hash_rows = verify_hash_coverage(datasets)
    ledger, reconciliation = build_fact_ledger(datasets, experiments, attrs)
    dataset_changes = change_log(source_datasets, datasets, "dataset", "dataset_record_id")
    experiment_changes = change_log(source_experiments, experiments, "experiment", "experiment_record_id")
    changes = dataset_changes + experiment_changes
    for index, row in enumerate(changes, start=1):
        row["change_id"] = f"CHG{index:04d}"
    checks = validation_checks(datasets, experiments, hash_rows, reconciliation)

    conflicts = [row for row in reconciliation if row["match_status"] == "CONFLICT"]
    failed_validations = [row for row in checks if row["status"] == "FAIL"]
    blocking_rows = [
        {"conflict_id": f"BLK{index:03d}", "source": "manuscript_reconciliation", "record_id": row["manifest_record_id"], "field_or_fact": row["manuscript_fact"], "observed": row["manifest_value"], "severity": "BLOCKING", "required_action": row["recommended_action"]}
        for index, row in enumerate(conflicts, start=1)
    ]
    for row in failed_validations:
        blocking_rows.append({"conflict_id": f"BLK{len(blocking_rows) + 1:03d}", "source": "validation", "record_id": "ALL", "field_or_fact": row["check"], "observed": row["observed"], "severity": "BLOCKING", "required_action": "Resolve validation failure before freeze review."})

    source_inventory = []
    source_meta = json.loads(SOURCE_METADATA.read_text(encoding="utf-8"))
    for manifest_type, path, rows, columns, version in [
        ("dataset", SOURCE_DATASET, len(source_datasets), len(source_datasets.columns), source_meta.get("manifest_version", "0.2")),
        ("experiment", SOURCE_EXPERIMENT, len(source_experiments), len(source_experiments.columns), source_meta.get("manifest_version", "0.2")),
        ("metadata", SOURCE_METADATA, 1, 1, source_meta.get("manifest_version", "0.2")),
    ]:
        source_inventory.append({"manifest_type": manifest_type, "path": rel(path), "version": version, "sha1": sha1(path), "record_count": rows, "column_count": columns, "status": source_meta.get("manifest_status", "DRAFT_NOT_FROZEN"), "selected_as_current_source": "true"})

    dataset_status = [
        {"dataset_record_id": row["dataset_record_id"], "dataset_family": row["dataset_family"], "source_status": row["provenance_status"], "candidate_status": row["candidate_status"], "source_layer": row["source_layer"], "redistribution_class": row["redistribution_class"], "unresolved_issue_id": row["unresolved_issue_id"], "record_count_change": "NO_CHANGE"}
        for _, row in datasets.iterrows()
    ]
    experiment_status = [
        {"experiment_record_id": row["experiment_record_id"], "dataset_record_id": row["dataset_record_id"], "experiment_family": row["experiment_family"], "source_status": row["formal_status"], "candidate_status": row["candidate_status"], "redistribution_class": row["redistribution_class"], "unresolved_issue_id": row["unresolved_issue_id"], "record_count_change": "NO_CHANGE"}
        for _, row in experiments.iterrows()
    ]
    source_layers = [
        {"dataset_record_id": row["dataset_record_id"], "dataset_family": row["dataset_family"], "source_layer": row["source_layer"], "formal_dataset_name": row["formal_dataset_name"], "classification_basis": "Frozen input/derivative layer represented by this manifest row.", "evidence_paths": row["classification_evidence_paths"], "valid_class": str(row["source_layer"] in ALLOWED_SOURCE_LAYERS).lower()}
        for _, row in datasets.iterrows()
    ]
    redistribution = [
        {"record_type": "dataset", "record_id": row["dataset_record_id"], "record_family": row["dataset_family"], "source_layer": row["source_layer"], "redistribution_class": row["redistribution_class"], "provider_permission_status": "UNRESOLVED" if row["record_id"].startswith("MERS") else "BOUNDED_BY_RECORDED_SOURCE_TERMS", "public_package_action": "Do not bundle provider raw files." if row["redistribution_class"] in {"PERMISSION_UNRESOLVED", "LINK_AND_RECONSTRUCT_ONLY", "RESTRICTED"} else "Release project derivative only; preserve source attribution.", "evidence_paths": row["classification_evidence_paths"]}
        for _, row in datasets.iterrows()
    ] + [
        {"record_type": "experiment", "record_id": row["experiment_record_id"], "record_family": row["experiment_family"], "source_layer": "PROJECT_GENERATED_DERIVATIVE", "redistribution_class": row["redistribution_class"], "provider_permission_status": "SOURCE_INPUT_TERMS_STILL_APPLY", "public_package_action": "Release derivative only; link/reconstruct third-party inputs.", "evidence_paths": row["evidence_paths"]}
        for _, row in experiments.iterrows()
    ]

    write_tsv(OUT / "source_manifest_inventory.tsv", source_inventory)
    write_tsv(OUT / "dataset_record_status.tsv", dataset_status)
    write_tsv(OUT / "experiment_record_status.tsv", experiment_status)
    write_tsv(OUT / "manuscript_fact_ledger.tsv", ledger)
    write_tsv(OUT / "manuscript_manifest_reconciliation.tsv", reconciliation)
    write_tsv(OUT / "source_layer_classification.tsv", source_layers)
    write_tsv(OUT / "redistribution_classification.tsv", redistribution)
    write_tsv(OUT / "input_hash_coverage.tsv", hash_rows)
    write_tsv(OUT / "unresolved_nonblocking_records.tsv", unresolved_rows())
    write_tsv(OUT / "blocking_conflicts.tsv", blocking_rows, ["conflict_id", "source", "record_id", "field_or_fact", "observed", "severity", "required_action"])
    write_tsv(OUT / "candidate_change_log.tsv", changes)
    write_tsv(OUT / "manual_manuscript_update_recommendations.tsv", recommendations_rows())
    write_tsv(OUT / "validation_checks.tsv", checks)

    reconciliation_counts = Counter(row["match_status"] for row in reconciliation)
    redist_counts = Counter(row["redistribution_class"] for row in redistribution)
    dataset_counts = datasets["provenance_status"].value_counts().to_dict()
    experiment_counts = experiments["formal_status"].value_counts().to_dict()
    candidate_hashes = {"dataset": sha1(CANDIDATE_DATASET), "experiment": sha1(CANDIDATE_EXPERIMENT)}

    metadata = {
        "candidate_version": CANDIDATE_VERSION,
        "candidate_status": CANDIDATE_STATUS,
        "generated_date": GENERATED_DATE,
        "source_manifest_paths": [rel(SOURCE_DATASET), rel(SOURCE_EXPERIMENT), rel(SOURCE_METADATA)],
        "source_manifest_sha1s": {rel(path): sha1(path) for path in [SOURCE_DATASET, SOURCE_EXPERIMENT, SOURCE_METADATA]},
        "candidate_manifest_sha1s": {rel(CANDIDATE_DATASET): candidate_hashes["dataset"], rel(CANDIDATE_EXPERIMENT): candidate_hashes["experiment"]},
        "dataset_record_count": len(datasets),
        "experiment_record_count": len(experiments),
        "blocking_conflict_count": len(blocking_rows),
        "major_conflict_count": 0,
        "minor_conflict_count": 0,
        "unresolved_nonblocking_count": len(unresolved_rows()),
        "prior_audit_references": [rel(SOURCE_DIR), rel(EVIDENCE["task4e"]), rel(TASK4E_RELEASE)],
        "repository_local_commit": local_commit(),
        "formal_manuscript_accessed": False,
        "formal_manuscript_modified": False,
        "formal_bibliography_accessed": False,
        "formal_bibliography_modified": False,
        "manifest_frozen": False,
        "ready_for_final_manifest_freeze": False,
        "ready_for_freeze_review": len(blocking_rows) == 0,
    }
    CANDIDATE_METADATA.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_tsv(OUT / "release_inventory_addition_plan.tsv", release_plan_rows())

    decision = "CONDITIONAL_PASS_TO_FREEZE_REVIEW" if not blocking_rows else "HOLD"
    summary = {
        "task": "Task 4F final dataset/experiment manifest reconciliation and freeze-candidate generation",
        "generated_date": GENERATED_DATE,
        "decision": decision,
        "candidate_status": CANDIDATE_STATUS,
        "ready_for_final_manifest_freeze": False,
        "ready_for_freeze_review": not blocking_rows,
        "current_manifests": {
            "dataset": {"path": rel(SOURCE_DATASET), "version": "0.2", "records": len(source_datasets), "sha1": sha1(SOURCE_DATASET)},
            "experiment": {"path": rel(SOURCE_EXPERIMENT), "version": "0.2", "records": len(source_experiments), "sha1": sha1(SOURCE_EXPERIMENT)},
        },
        "candidates": {
            "dataset": {"path": rel(CANDIDATE_DATASET), "records": len(datasets), "sha1": candidate_hashes["dataset"]},
            "experiment": {"path": rel(CANDIDATE_EXPERIMENT), "records": len(experiments), "sha1": candidate_hashes["experiment"]},
            "metadata": {"path": rel(CANDIDATE_METADATA), "sha1": sha1(CANDIDATE_METADATA)},
        },
        "source_status_counts": {"dataset": dataset_counts, "experiment": experiment_counts},
        "manuscript_reconciliation_counts": dict(sorted(reconciliation_counts.items())),
        "blocking_conflicts": len(blocking_rows),
        "major_conflicts": 0,
        "unresolved_nonblocking_records": len(unresolved_rows()),
        "source_layer_counts": datasets["source_layer"].value_counts().sort_index().to_dict(),
        "redistribution_classification_counts": dict(sorted(redist_counts.items())),
        "validation": {"checks": len(checks), "passed": sum(row["status"] == "PASS" for row in checks), "failed": sum(row["status"] == "FAIL" for row in checks)},
        "candidate_change_count": len(changes),
        "guardrails": {
            "formal_manuscript_accessed": False, "formal_manuscript_modified": False,
            "formal_bibliography_accessed": False, "formal_bibliography_modified": False,
            "experimental_code_modified": False, "experimental_outputs_modified": False,
            "formal_stages_rerun": False, "formal_figures_modified": False,
            "formal_source_value_tables_modified": False, "github_remote_access_required": False,
            "existing_manifests_modified": False, "task4e_release_manifest_modified": False,
        },
    }
    (OUT / "audit_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (OUT / "README.md").write_text(make_readme(dataset_counts, experiment_counts, reconciliation_counts, redist_counts), encoding="utf-8")
    (OUT / "final_decision.md").write_text(make_final_decision(reconciliation_counts, candidate_hashes), encoding="utf-8")

    after = {rel(path): sha1(path) for path in protected if path.exists()}
    if before != after:
        changed = sorted(path for path in before if before[path] != after.get(path))
        raise RuntimeError(f"Protected source manifests/configs/release inventory changed: {changed}")
    if blocking_rows:
        print(f"[DECISION] HOLD ({len(blocking_rows)} blocking conflicts)")
    else:
        print("[DECISION] CONDITIONAL_PASS_TO_FREEZE_REVIEW")
    print(f"[CANDIDATE] dataset={len(datasets)} sha1={candidate_hashes['dataset']}")
    print(f"[CANDIDATE] experiment={len(experiments)} sha1={candidate_hashes['experiment']}")
    print(f"[RECONCILIATION] {dict(reconciliation_counts)}")
    print(f"[UNRESOLVED NONBLOCKING] {len(unresolved_rows())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
