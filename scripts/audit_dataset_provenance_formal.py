#!/usr/bin/env python3
"""Build the submission-facing dataset provenance audit without changing project outputs."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "visualizations" / "manuscript_audits" / "dataset_provenance_formal_audit"
OLD_AUDIT = ROOT / "manuscript" / "reference_audit" / "lowres_and_reference_dropout_dataset_source_audit.csv"
PUBLIC_AUDIT = ROOT / "visualizations" / "manuscript_public_source_resolution" / "public_source_resolution.csv"
ACCESSED = "2026-07-22"

STATUSES = {
    "CONFIRMED",
    "PARTIAL",
    "UNRESOLVED",
    "CONFLICT",
    "NOT_APPLICABLE",
    "PROJECT_GENERATED",
}

MASTER_FIELDS = [
    "record_id", "dataset_family", "experiment_id", "internal_dataset_id", "data_role",
    "formal_dataset_name", "formal_release_name", "source_study_title", "species", "tissue",
    "disease", "donor_or_patient_id", "sample_id", "section_id", "slide_id", "capture_area",
    "spatial_platform", "assay", "tissue_processing", "panel_or_gene_count", "reference_type",
    "reference_dataset_name", "reference_sample_subset", "database", "study_accession",
    "project_accession", "experiment_accession", "run_accession", "sample_accession",
    "biosample_accession", "geo_accession", "arrayexpress_accession", "ena_accession",
    "official_source_url", "provider_download_url", "source_publication", "bibtex_key",
    "input_file", "input_file_sha1", "config_file", "config_file_sha1", "manifest_file",
    "manifest_row", "manuscript_location", "redistribution_status", "project_generated_derivative",
    "status", "evidence", "notes",
]

CROSSWALK_FIELDS = [
    "experiment_id", "figure_panel", "manuscript_result_section", "spatial_input_record_id",
    "reference_input_record_id", "evaluation_resource_record_id", "profile_mask_target",
    "reference_dropout_target", "project_derivative", "formal_status", "notes",
]


def rel(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return str(path)


def sha1(path_text: str) -> str:
    if not path_text or path_text in STATUSES:
        return path_text or "NOT_APPLICABLE"
    path = ROOT / path_text
    if not path.is_file():
        return "UNRESOLVED"
    digest = hashlib.sha1()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha1_list(value: str) -> str:
    paths = [item.strip() for item in value.split(";") if item.strip()]
    if not paths:
        return "NOT_APPLICABLE"
    return ";".join(f"{item}={sha1(item)}" for item in paths)


def base_record(**values: str) -> dict[str, str]:
    row = {field: "NOT_APPLICABLE" for field in MASTER_FIELDS}
    row.update({key: str(value) for key, value in values.items()})
    row.setdefault("manuscript_location", "UNRESOLVED")
    return row


def assign_accessions(row: dict[str, str], accessions: str) -> None:
    buckets: dict[str, list[str]] = defaultdict(list)
    for token in re.split(r"[;,]", accessions or ""):
        token = token.strip()
        if not token or token in STATUSES or token.lower() in {"none", "multiple"}:
            continue
        if token.startswith("GSE"):
            buckets["geo_accession"].append(token)
        elif token.startswith("GSM"):
            buckets["sample_accession"].append(token)
        elif token.startswith("E-MTAB-"):
            buckets["arrayexpress_accession"].append(token)
        elif token.startswith(("PRJEB", "PRJNA")):
            buckets["project_accession"].append(token)
            buckets["ena_accession"].append(token)
        elif token.startswith("ERR"):
            buckets["run_accession"].append(token)
            buckets["ena_accession"].append(token)
        elif token.startswith("ERS"):
            buckets["sample_accession"].append(token)
            buckets["ena_accession"].append(token)
        elif token.startswith("SAMEA"):
            buckets["biosample_accession"].append(token)
        elif token.startswith("SRP"):
            buckets["study_accession"].append(token)
        elif token.startswith("EGAS"):
            buckets["study_accession"].append(token)
    for field, items in buckets.items():
        existing = [] if row[field] == "NOT_APPLICABLE" else row[field].split(";")
        row[field] = ";".join(dict.fromkeys(existing + items))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


OLD_ROWS = {row["internal_scenario_id"]: row for row in read_csv(OLD_AUDIT)}
PUBLIC_ROWS = read_csv(PUBLIC_AUDIT)


LOWRES_INPUTS = {
    "adult_mouse_kidney_real": ("Adult Mouse Kidney", "Visium_FFPE_Mouse_Kidney_filtered_feature_bc_matrix.h5", "adata.h5ad"),
    "ffpe_mouse_brain_sagittal_real": ("FFPE Mouse Brain Sagittal", "CytAssist_FFPE_Sagittal_Mouse_Brain_filtered_feature_bc_matrix.h5", "sc.h5ad"),
    "human_breast_cancer_real": ("Human Breast Cancer", "Visium_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix.h5", "sc_human_breast_heca.h5ad"),
    "human_breast_cancer_visium_ff_wta_real": ("Human Breast Cancer Visium FF WTA", "Visium_Human_Breast_Cancer_filtered_feature_bc_matrix.h5", "sc_human_breast_heca.h5ad"),
    "human_breast_cancer_wta_120_real": ("Human Breast Cancer WTA 1.2.0", "Parent_Visium_Human_BreastCancer_filtered_feature_bc_matrix.h5", "sc_human_breast_heca.h5ad"),
    "human_cervical_cancer_real": ("Human Cervical Cancer", "Visium_FFPE_Human_Cervical_Cancer_filtered_feature_bc_matrix.h5", "sc_human_uterus_heca.h5ad"),
    "human_heart_ff_real": ("Human Heart (FF)", "V1_Human_Heart_filtered_feature_bc_matrix.h5", "sc_human_heart_abo1984.h5ad"),
    "human_intestine_cancer_real": ("Human Intestine Cancer", "Visium_FFPE_Human_Intestinal_Cancer_filtered_feature_bc_matrix.h5", "sc_human_intestine_heca.h5ad"),
    "human_lymph_node_real": ("Human Lymph Node", "V1_Human_Lymph_Node_filtered_feature_bc_matrix.h5", "sc_lymph_node_abo0510.h5ad"),
    "mouse_embryo_real": ("Mouse Embryo", "dataset.h5ad", "mouse_atlas_gottgens_stelzer.h5ad"),
}


def lowres_input(base_id: str) -> str:
    folder, st_file, ref_file = LOWRES_INPUTS[base_id]
    root = Path("data/raw/low_resolution") / folder
    return f"{(root / st_file).as_posix()};{(root / ref_file).as_posix()}"


def pair_from_old(record_id: str, scenario_id: str, family: str, input_file: str, status: str | None = None) -> dict[str, str]:
    source = OLD_ROWS[scenario_id]
    species = "Mus musculus" if "mouse" in source["spatial_tissue"].lower() else "Homo sapiens"
    result = base_record(
        record_id=record_id,
        dataset_family=family,
        experiment_id=scenario_id,
        internal_dataset_id=scenario_id,
        data_role="paired spatial input and single-cell/single-nucleus reference input",
        formal_dataset_name=source["spatial_formal_dataset_name"],
        formal_release_name=source["spatial_formal_dataset_name"],
        source_study_title=source["spatial_original_paper_title"],
        species=species,
        tissue=source["spatial_tissue"],
        disease=source["spatial_disease"],
        donor_or_patient_id=source["spatial_sample_id"],
        sample_id=source["spatial_sample_id"],
        section_id=source["spatial_sample_id"],
        slide_id=source["spatial_slide_id"],
        capture_area=source["spatial_capture_area"],
        spatial_platform=source["spatial_platform"],
        assay=source["spatial_assay"],
        tissue_processing=source["spatial_processing_type"],
        panel_or_gene_count="whole transcriptome or probe-based WTA as specified by release",
        reference_type=source["reference_assay"],
        reference_dataset_name=source["reference_formal_dataset_name"],
        reference_sample_subset=source["reference_sample_or_cohort"],
        database=source["spatial_accession_database"],
        official_source_url=";".join(dict.fromkeys(filter(None, [source["spatial_public_url"], source["reference_public_url"]]))),
        provider_download_url=source["spatial_public_url"],
        source_publication=";".join(dict.fromkeys(filter(None, [source["spatial_original_paper_title"], source["reference_original_paper_title"]]))),
        bibtex_key=";".join(filter(None, [source["recommended_spatial_bibtex_key"], source["recommended_reference_bibtex_key"]])),
        input_file=input_file,
        config_file=source["config_path"],
        manifest_file=source["final_source_table_path"],
        manifest_row=f"internal_scenario_id={scenario_id}",
        manuscript_location=f"{source['figure']} {source['figure_panel']}; formal LaTeX file absent",
        redistribution_status="See provider/repository terms; project does not establish redistribution permission",
        project_generated_derivative="profile-masked or reference-dropout copy generated by this project",
        status=status or ("CONFIRMED" if source["confidence"] == "HIGH" else "PARTIAL"),
        evidence=f"{source['repository_evidence']}; {source['public_source_evidence']}",
        notes=source["notes"] + (f" Unresolved: {source['unresolved_fields']}." if source["unresolved_fields"] != "none" else ""),
    )
    assign_accessions(result, source["spatial_accession"])
    assign_accessions(result, source["reference_accession"])
    return result


def public_row(internal_id: str) -> dict[str, str]:
    return next(row for row in PUBLIC_ROWS if row["internal_dataset_id"] == internal_id)


def build_records() -> list[dict[str, str]]:
    records: list[dict[str, str]] = []

    lowres_specs = [
        ("LR01_KIDNEY", "adult_mouse_kidney_real_profile_mask_endo", "adult_mouse_kidney_real"),
        ("LR02_BRAIN", "ffpe_mouse_brain_sagittal_real_profile_mask_microglia", "ffpe_mouse_brain_sagittal_real"),
        ("LR03_BRCA_FFPE", "human_breast_cancer_real_profile_mask_basal_cell", "human_breast_cancer_real"),
        ("LR04_BRCA_FF", "human_breast_cancer_visium_ff_wta_real_profile_mask_macrophage", "human_breast_cancer_visium_ff_wta_real"),
        ("LR05_BRCA_ILC", "human_breast_cancer_wta_120_real_profile_mask_endothelial_cell", "human_breast_cancer_wta_120_real"),
        ("LR06_CERVICAL", "human_cervical_cancer_real_profile_mask_epithelial_cell", "human_cervical_cancer_real"),
        ("LR07_HEART", "human_heart_ff_real_profile_mask_endothelial_cell", "human_heart_ff_real"),
        ("LR08_INTESTINE", "human_intestine_cancer_real_profile_mask_endothelial_cell", "human_intestine_cancer_real"),
        ("LR09_LYMPH_NODE", "human_lymph_node_real_profile_mask_b_cell", "human_lymph_node_real"),
        ("LR10_EMBRYO", "mouse_embryo_real_profile_mask_erythroid", "mouse_embryo_real"),
    ]
    for record_id, scenario_id, base_id in lowres_specs:
        records.append(pair_from_old(record_id, scenario_id, "low-resolution profile masking", lowres_input(base_id)))

    dropout_specs = [
        ("RD01_EMBRYO_ENDODERM", "mouse_embryo_real_sc_missing_endoderm_gut", lowres_input("mouse_embryo_real")),
        ("RD02_EMBRYO_ERYTHROID", "mouse_embryo_real_sc_missing_erythroid", lowres_input("mouse_embryo_real")),
        ("RD03_TNBC_PLASMA", "cytospace_fig2d_tme_brca_tnbc_fresh_frozen_sc_missing_plasma_cells", "data/processed/cytospace_fig2d_tme/cytospace_fig2d_tme_brca_tnbc_fresh_frozen/stage1_preprocess/exported/st_expression_normalized.csv;data/processed/cytospace_fig2d_tme/cytospace_fig2d_tme_brca_tnbc_fresh_frozen/stage1_preprocess/exported/sc_expression_normalized.csv"),
        ("RD04_CRC_B", "cytospace_fig2d_tme_crc_fresh_frozen_sc_missing_b_cells", "data/processed/cytospace_fig2d_tme/cytospace_fig2d_tme_crc_fresh_frozen/stage1_preprocess/exported/st_expression_normalized.csv;data/processed/cytospace_fig2d_tme/cytospace_fig2d_tme_crc_fresh_frozen/stage1_preprocess/exported/sc_expression_normalized.csv"),
        ("RD05_BRCA_FFPE_PLASMA", "cytospace_fig2d_tme_brca_her2_ffpe_sc_missing_plasma_cells", "data/processed/cytospace_fig2d_tme/cytospace_fig2d_tme_brca_her2_ffpe/stage1_preprocess/exported/st_expression_normalized.csv;data/processed/cytospace_fig2d_tme/cytospace_fig2d_tme_brca_her2_ffpe/stage1_preprocess/exported/sc_expression_normalized.csv"),
        ("RD06_BRCA_FFPE_EPITHELIAL", "cytospace_fig2d_tme_brca_her2_ffpe_sc_missing_epithelial_cells", "data/processed/cytospace_fig2d_tme/cytospace_fig2d_tme_brca_her2_ffpe/stage1_preprocess/exported/st_expression_normalized.csv;data/processed/cytospace_fig2d_tme/cytospace_fig2d_tme_brca_her2_ffpe/stage1_preprocess/exported/sc_expression_normalized.csv"),
    ]
    for record_id, scenario_id, input_file in dropout_specs:
        records.append(pair_from_old(record_id, scenario_id, "real reference dropout", input_file))

    sim_specs = [
        {
            "record_id": "SIM01_BRCA",
            "experiment_id": "real_brca_joint_stage3ab",
            "internal_dataset_id": "real_brca",
            "formal_dataset_name": "Project-generated BRCA clustered simulation derived from 10x FFPE breast scaffold and Wu atlas",
            "formal_release_name": "Human Breast Cancer: Ductal Carcinoma In Situ, Invasive Carcinoma (FFPE)",
            "source_study_title": "A single-cell and spatially resolved atlas of human breast cancers",
            "species": "Homo sapiens", "tissue": "breast", "disease": "breast cancer",
            "donor_or_patient_id": "738811QB spatial block; Wu HER2+ reference patients",
            "sample_id": "Visium_FFPE_Human_Breast_Cancer", "section_id": "Block 738811QB, Section 1",
            "slide_id": "V11J26-008", "capture_area": "B1", "spatial_platform": "10x Visium Spatial Gene Expression",
            "assay": "FFPE whole-transcriptome probe-based spatial expression", "tissue_processing": "FFPE direct placement",
            "reference_type": "scRNA-seq", "reference_dataset_name": "Wu breast-cancer atlas",
            "reference_sample_subset": "3,977 HER2+ cells: CID3921 1,597; CID45171 1,138; CID3838 1,242; differs from the 4,014-cell CTA subset",
            "database": "10x Genomics; NCBI GEO", "geo_accession": "GSE176078",
            "official_source_url": "https://www.10xgenomics.com/datasets/human-breast-cancer-ductal-carcinoma-in-situ-invasive-carcinoma-ffpe-1-standard-1-3-0;https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078",
            "provider_download_url": "https://www.10xgenomics.com/datasets/human-breast-cancer-ductal-carcinoma-in-situ-invasive-carcinoma-ffpe-1-standard-1-3-0",
            "source_publication": "Wu et al., Nature Genetics (2021), doi:10.1038/s41588-021-00911-1", "bibtex_key": "Wu2021BreastCancerAtlas;TenXFFPEBreastCancer",
            "input_file": "data/sim/real_brca/real_brca7_endothelial_marker_control_sc_missing_endothelial_cells/brca_STdata_GEP.txt;data/processed/simulation_experiments/real_brca/real_brca7_endothelial_marker_control_sc_missing_endothelial_cells/stage1_preprocess/exported/sc_expression_normalized.csv",
            "config_file": "configs/datasets/real_brca7_endothelial_marker_control_sc_missing_endothelial_cells.yaml",
            "status": "CONFIRMED",
        },
        {
            "record_id": "SIM02_MOUSE_BRAIN", "experiment_id": "mouse_brain_refined_joint_stage3ab", "internal_dataset_id": "mouse_brain_refined",
            "formal_dataset_name": "Project-generated refined mouse-brain clustered simulation derived from a 10x CytAssist FFPE sagittal-brain scaffold",
            "formal_release_name": "Preservation Method Comparison on Visium CytAssist: FFPE Mouse Brain, Sagittal, 11 mm",
            "source_study_title": "Cell2location maps fine-grained cell types in spatial transcriptomics", "species": "Mus musculus", "tissue": "adult sagittal brain", "disease": "healthy",
            "sample_id": "CytAssist_FFPE_Sagittal_Mouse_Brain", "section_id": "10x provider demonstration section", "slide_id": "V52B25-081", "capture_area": "B",
            "spatial_platform": "10x Visium CytAssist", "assay": "FFPE whole-transcriptome probe-based spatial expression", "tissue_processing": "FFPE CytAssist, 11 mm",
            "reference_type": "single-nucleus RNA-seq", "reference_dataset_name": "Cell2location adult mouse-brain reference", "reference_sample_subset": "6,112 locally refined nuclei/cells across eight broad types",
            "database": "10x Genomics; BioStudies/ArrayExpress", "arrayexpress_accession": "E-MTAB-11115",
            "official_source_url": "https://www.10xgenomics.com/datasets/preservation-method-comparison-on-visium-cytassist-ffpe-mouse-brain-sagittal-11-mm-capture-area-2-standard;https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-11115",
            "provider_download_url": "https://www.10xgenomics.com/datasets/preservation-method-comparison-on-visium-cytassist-ffpe-mouse-brain-sagittal-11-mm-capture-area-2-standard",
            "source_publication": "Kleshchevnikov et al., Nature Biotechnology (2022), doi:10.1038/s41587-021-01139-4", "bibtex_key": "Kleshchevnikov2022Cell2location;TenXFFPEMouseBrainSagittal",
            "input_file": "data/sim/mouse_brain_refined/mouse_brain_refined7_balanced_clustered_sim_sc_missing_ext_l56/brca_STdata_GEP.txt;data/sim/mouse_brain_refined/mouse_brain_refined7_balanced_clustered_sim_sc_missing_ext_l56/brca_scRNA_GEP.txt",
            "config_file": "configs/datasets/mouse_brain_refined7_balanced_clustered_sim_sc_missing_ext_l56.yaml", "status": "CONFIRMED",
        },
        {
            "record_id": "SIM03_HUMAN_LUNG", "experiment_id": "human_lung_5loc_joint_stage3ab", "internal_dataset_id": "human_lung_5loc",
            "formal_dataset_name": "Project-generated clustered simulation from Madissoon human lung atlas sample WSA_LngSP10193345",
            "formal_release_name": "A spatially resolved atlas of the human lung characterizes a gland-associated immune niche",
            "source_study_title": "A spatially resolved atlas of the human lung characterizes a gland-associated immune niche", "species": "Homo sapiens", "tissue": "upper-left-lobe lung parenchyma", "disease": "normal",
            "donor_or_patient_id": "A48", "sample_id": "WSA_LngSP10193345", "section_id": "e_TopLeftPar", "slide_id": "UNRESOLVED", "capture_area": "UNRESOLVED",
            "spatial_platform": "10x Visium Spatial Gene Expression", "assay": "Spatial 3' v1 whole-transcriptome expression", "tissue_processing": "OCT embedded, flash frozen; 10 micrometre section",
            "reference_type": "single-cell RNA-seq", "reference_dataset_name": "Madissoon lung atlas reference", "reference_sample_subset": "15,000 selected cells, donor A48/e_TopLeftPar, nine labels",
            "database": "BioStudies/ArrayExpress; ENA", "arrayexpress_accession": "E-MTAB-11640", "project_accession": "PRJEB52292", "sample_accession": "ERS20065156", "biosample_accession": "SAMEA115633909", "ena_accession": "PRJEB52292;ERS20065156",
            "official_source_url": "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-11640;https://www.ebi.ac.uk/ena/browser/view/ERS20065156;https://doi.org/10.1038/s41588-022-01243-4",
            "source_publication": "Madissoon et al., Nature Genetics (2023), doi:10.1038/s41588-022-01243-4", "bibtex_key": "Madissoon2022LungAtlas",
            "input_file": "data/sim/human_lung_5loc/human_lung_5loc_fine9_clustered_sim_sc_missing_b_cell/brca_STdata_GEP.txt;data/processed/simulation_experiments/human_lung_5loc/human_lung_5loc_fine9_clustered_sim_sc_missing_b_cell/stage1_preprocess/exported/sc_expression_normalized.csv",
            "config_file": "configs/datasets/human_lung_5loc_fine9_clustered_sim_sc_missing_b_cell.yaml", "status": "CONFLICT",
            "notes": "ERS20065156 is an ENA sample accession, not a run accession. The five-location study sampled trachea, 2nd/3rd-generation bronchi, 4th-generation bronchi, upper-left-lobe parenchyma and lower-left-lobe parenchyma; this benchmark uses one A48 upper-lobe sample.",
        },
    ]
    for spec in sim_specs:
        records.append(base_record(
            dataset_family="joint Stage3A-Stage3B simulations", data_role="project-generated simulation input and truth",
            panel_or_gene_count="project-generated clustered spot mixtures", project_generated_derivative="PROJECT_GENERATED",
            manifest_file="visualizations/method_comparison/composite_no_noise/composite_no_noise_scenario_manifest.csv",
            manifest_row=f"dataset_group={spec['internal_dataset_id']}", manuscript_location="Fig. 1C-D; formal LaTeX file absent",
            redistribution_status="Project derivative; source-data terms still apply", evidence="simulation config and local generated matrices; official source pages; manuscript public-source audit", **spec,
        ))

    merscope = [
        ("MERS01_BREAST", "HumanBreastCancerPatient1", "breast cancer", "713121", "Monocytes and Macrophages"),
        ("MERS02_COLON", "HumanColonCancerPatient1", "colon cancer", "677451", "Fibroblasts"),
        ("MERS03_LUNG", "HumanLungCancerPatient1", "lung cancer", "353762", "Plasma cells"),
        ("MERS04_MELANOMA1", "HumanMelanomaPatient1", "melanoma", "468138", "Fibroblasts"),
        ("MERS05_MELANOMA2", "HumanMelanomaPatient2", "melanoma", "207869", "B cells"),
    ]
    for record_id, sample, disease, folder_id, target in merscope:
        internal = sample.lower()
        input_file = f"data/raw/high/{sample}/{sample}_cell_by_gene.csv;data/raw/high/{sample}/{sample}_cell_metadata.csv"
        records.append(base_record(
            record_id=record_id, dataset_family="high-resolution MERSCOPE profile masking", experiment_id=f"highres_{internal}_profile_mask", internal_dataset_id=f"highres_{internal}",
            data_role="MERSCOPE spatial input plus project-generated disjoint same-assay reference split", formal_dataset_name=f"Vizgen release sample {sample}",
            formal_release_name="Vizgen MERFISH FFPE Human Immuno-oncology Data Set, May 2022", source_study_title="Vizgen Human FFPE Immuno-oncology Data Release",
            species="Homo sapiens", tissue=disease.replace(" cancer", ""), disease=disease, donor_or_patient_id="UNRESOLVED; release identifier is not a clinical patient ID",
            sample_id=sample, section_id=f"provider release folder token {folder_id}", slide_id="NOT_APPLICABLE", capture_area="NOT_APPLICABLE",
            spatial_platform="MERSCOPE", assay="MERFISH", tissue_processing="FFPE", panel_or_gene_count="500 target genes + 50 blank controls (551 CSV columns including cell ID)",
            reference_type="project-generated disjoint same-assay reference split", reference_dataset_name=f"Disjoint reference pool from {sample}",
            reference_sample_subset=f"Spatial and reference cell IDs split with zero overlap; masking target: {target}", database="Vizgen provider-hosted data release",
            official_source_url="https://vizgen.com/human-ffpe-immunooncology-release-roadmap/", provider_download_url="https://info.vizgen.com/merscope-ffpe-solution",
            source_publication="Vizgen provider data release; no external repository accession assigned", bibtex_key="Vizgen2022FFPEImmunoOncology",
            input_file=input_file, config_file=f"configs/datasets/highres_{internal}_profile_mask_{target.lower().replace(' ', '_')}.yaml" if target != "Monocytes and Macrophages" else "configs/datasets/highres_humanbreastcancerpatient1_profile_mask_monocytes_and_macrophages.yaml",
            manifest_file="visualizations/highres_profile_mask_mapping/highres_profile_mask_mapping_manifest.csv", manifest_row=f"raw_sample={sample}",
            manuscript_location="Fig. 3A-D; formal LaTeX file absent", redistribution_status="Public provider download; exact downstream redistribution terms not preserved in repository",
            project_generated_derivative="PROJECT_GENERATED disjoint same-assay reference split", status="PARTIAL",
            evidence=f"raw cell_by_gene and metadata files; highres mapping manifest; Vizgen provider page; release folder token {folder_id} retained from prior source audit",
            notes="The sample string is a provider release identifier, not a repository accession or verified clinical patient identifier.",
        ))

    records.extend([
        base_record(
            record_id="STATE01_TCELL_DECOY", dataset_family="state-decoy experiments", experiment_id="cytospace_fig2k_breast_state_decoy_stage3_detected",
            internal_dataset_id="cytospace_fig2k_breast_state_decoy_stage3_detected", data_role="project-generated engineered decoy plus literature-derived ordering",
            formal_dataset_name="HumanBreastCancerPatient1 T-cell-state decoy experiment", formal_release_name="Vizgen MERFISH FFPE Human Immuno-oncology Data Set, May 2022",
            source_study_title="Pan-cancer single-cell landscape of tumor-infiltrating T cells", species="Homo sapiens", tissue="breast", disease="breast cancer",
            sample_id="HumanBreastCancerPatient1", spatial_platform="MERSCOPE", assay="MERFISH", tissue_processing="FFPE", panel_or_gene_count="500 target genes + 50 blanks",
            reference_type="engineered decoy", reference_dataset_name="23 CD4 T-cell state ordering", reference_sample_subset="Zheng et al. Supplementary Table S3; numerical ordering transcribed via CytoSPACE Supplementary Table S9",
            database="Vizgen provider release; literature supplement", official_source_url="https://vizgen.com/human-ffpe-immunooncology-release-roadmap/;https://doi.org/10.1126/science.abe6474;https://doi.org/10.1038/s41587-023-01697-9",
            source_publication="Zheng et al. Science (2021); Vahid et al. Nature Biotechnology (2023)", bibtex_key="Zheng2021PanCancerTCells;Vahid2023CytoSPACE",
            input_file="data/raw/high/HumanBreastCancerPatient1/HumanBreastCancerPatient1_cell_by_gene.csv", config_file="configs/datasets/cytospace_fig2k_breast_state_decoy_stage3_detected.yaml",
            manifest_file="visualizations/cytospace_fig2k_tcell_states_stage3_decoy/fig2k_stage3_detected_state_decoy_baseline_vs_route2_manifest.json", manifest_row="entire manifest",
            manuscript_location="Fig. 3E; formal LaTeX file absent", redistribution_status="Project decoy; source terms apply", project_generated_derivative="PROJECT_GENERATED",
            status="PROJECT_GENERATED", evidence="local config/manifest; Zheng paper supplement; CytoSPACE supplement used only as transcription source",
            notes="CytoSPACE Supplementary Table S9 is an evaluation/transcription resource, not a biological dataset or independent truth.",
        ),
        base_record(
            record_id="STATE02_KIDNEY_STATE32", dataset_family="state-decoy experiments", experiment_id="cytospace_fig2i_mouse_kidney_stage3_unsupported_decoy_n1000",
            internal_dataset_id="cytospace_fig2i_mouse_kidney_stage3_unsupported_decoy_n1000", data_role="external spatial/reference input plus project-generated decoys",
            formal_dataset_name="Adult Mouse Kidney (FFPE) with Ransick/KidneyCellExplorer reference and State32-like decoys", formal_release_name="Adult Mouse Kidney (FFPE)",
            source_study_title="Single-cell transcriptomic profiling of the adult mouse kidney", species="Mus musculus", tissue="adult kidney", disease="healthy",
            sample_id="Visium_FFPE_Mouse_Kidney", slide_id="V11A13-021", capture_area="C1", spatial_platform="10x Visium Spatial Gene Expression",
            assay="FFPE probe-based spatial expression", tissue_processing="FFPE direct placement", panel_or_gene_count="whole-transcriptome probe panel",
            reference_type="scRNA-seq epithelial-state reference", reference_dataset_name="Ransick/KidneyCellExplorer", reference_sample_subset="13,195 authentic cells plus 1,000 project-generated State32-like decoys",
            database="10x Genomics; NCBI GEO/SRA", project_accession="PRJNA532850", study_accession="SRP192559", geo_accession="GSE129798",
            official_source_url="https://www.10xgenomics.com/datasets/adult-mouse-kidney-ffpe-1-standard-1-3-0;https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129798;https://doi.org/10.1016/j.devcel.2019.10.005",
            source_publication="Ransick et al., Developmental Cell (2019)", bibtex_key="Ransick2019KidneyStates;Vahid2023CytoSPACE",
            input_file="data/raw/cytospace_fig2i_mouse_kidney_strict/kidneycellexplorer-master/data/expr_data.rds;data/processed/cytospace_fig2i_mouse_kidney_stage3_unsupported_decoy_n1000/stage1_preprocess/exported/st_expression_normalized.csv", config_file="configs/datasets/cytospace_fig2i_mouse_kidney_stage3_unsupported_decoy_n1000.yaml",
            manuscript_location="Fig. 3F; formal LaTeX file absent", redistribution_status="Project decoy; source terms apply", project_generated_derivative="PROJECT_GENERATED 1,000 State32-like decoys",
            status="PROJECT_GENERATED", evidence="10x official page; GEO GSE129798; local strict-source directory and config",
            notes="KidneyCellExplorer is a portal, not an accession. State 32 is 'Deep medullary epithelium of pelvis'. Original CytoSPACE Fig. 2i used Ferreira GSE171406, which is not this local reference.",
        ),
    ])

    records.append(base_record(
        record_id="THAL01_ST8059051", dataset_family="detailed Thalamic case", experiment_id="cell2location_ST8059051_thalamic_excitatory_reference_missing",
        internal_dataset_id="cell2location_ST8059051_thalamic_excitatory_reference_missing", data_role="paired spatial and snRNA reference input",
        formal_dataset_name="Cell2location mouse-brain section ST8059051", formal_release_name="E-MTAB-11114 spatial study paired with E-MTAB-11115 reference",
        source_study_title="Cell2location maps fine-grained cell types in spatial transcriptomics", species="Mus musculus", tissue="adult mouse brain", disease="healthy",
        sample_id="ST8059051", section_id="Visium-29B", slide_id="C05717-021", capture_area="B1", spatial_platform="10x Visium Spatial Gene Expression",
        assay="fresh-frozen whole-transcriptome spatial expression", tissue_processing="fresh frozen cryosection", reference_type="single-nucleus RNA-seq",
        reference_dataset_name="Cell2location adult mouse-brain reference", reference_sample_subset="thalamic excitatory states removed in project perturbation",
        database="BioStudies/ArrayExpress", arrayexpress_accession="E-MTAB-11114;E-MTAB-11115", official_source_url="https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-11114;https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-11115;https://doi.org/10.1038/s41587-021-01139-4",
        source_publication="Kleshchevnikov et al., Nature Biotechnology (2022)", bibtex_key="Kleshchevnikov2022Cell2location",
        input_file="data/processed/cell2location_mouse_brain/cell2loc_scan_st8059051_sc_missing_thalamic_excitatory/stage1_preprocess/exported/st_expression_normalized.csv;data/processed/cell2location_mouse_brain/cell2loc_scan_st8059051_sc_missing_thalamic_excitatory/stage1_preprocess/exported/sc_expression_normalized.csv",
        config_file="configs/datasets/cell2loc_scan_st8059051_sc_missing_thalamic_excitatory.yaml", manuscript_location="Fig. 4A-E/H; formal LaTeX file absent",
        redistribution_status="Public study; repository terms apply", project_generated_derivative="reference-dropout copy and marker-region proxy",
        status="CONFIRMED", evidence="config and Stage1 exported matrices; E-MTAB-11114/11115; source publication",
        notes="This fresh-frozen C05717-021/B1 section is not the FFPE CytAssist V52B25-081/B scaffold used in the joint simulation. Accuracy-like readouts are reference-relative proxies.",
    ))

    records.append(base_record(
        record_id="CTA01_BCSA2TUMB1", dataset_family="CTA biological application", experiment_id="bioapp_BCSA2TumB1_CTA_immune_endpoint",
        internal_dataset_id="bioapp_BCSA2TumB1_CTA_immune_endpoint", data_role="external spatial input, external image-derived endpoint, and project-balanced external reference subset",
        formal_dataset_name="Breast Cancer CTA section BCSA2TumB1", formal_release_name="Computational pathology annotation enhances the resolution and interpretation of breast cancer spatial transcriptomics data",
        source_study_title="Computational pathology annotation enhances the resolution and interpretation of breast cancer spatial transcriptomics data", species="Homo sapiens", tissue="breast", disease="HER2-positive breast cancer",
        donor_or_patient_id="BCSA2", sample_id="BCSA2TumB1", section_id="TumB1", slide_id="UNRESOLVED", capture_area="UNRESOLVED", spatial_platform="10x Visium Spatial Gene Expression",
        assay="whole-transcriptome spatial expression", tissue_processing="fresh frozen", panel_or_gene_count="2,248 repository-frozen spots; 2,000-gene mapping intersection",
        reference_type="scRNA-seq project-balanced donor-unmatched subset", reference_dataset_name="Wu breast-cancer atlas", reference_sample_subset="4,014 cells: CID4066 1,211; CID3586 1,038; CID3921 766; CID3838 544; CID45171 455; most cell types capped at 400, plasma 226, generic T cells 188",
        database="Swedish National Data Service; Zenodo; Researchdata.se; NCBI GEO", study_accession="10.48723/f4v5-m008", geo_accession="GSE176078", project_accession="Researchdata.se 2025-97",
        official_source_url="https://doi.org/10.1038/s41698-025-01104-3;https://doi.org/10.48723/f4v5-m008;https://doi.org/10.5281/zenodo.15211538;https://www.researchdata.se/en/catalogue/dataset/2025-97;https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078",
        provider_download_url="https://doi.org/10.5281/zenodo.15211538", source_publication="Li et al., npj Precision Oncology 9:310 (2025), doi:10.1038/s41698-025-01104-3; Wu et al. (2021)",
        bibtex_key="Li2025CTA;Wu2021BreastCancerAtlas", input_file="visualizations/bioapp_experiment/bioapp_phase3r_review_resolution_for_cta_endpoint_baseline_feasibility/st_expression_full_2248_counts.csv.gz;data/processed/cytospace_fig2d_tme/cytospace_fig2d_tme_brca_her2_ffpe/stage1_preprocess/exported/sc_metadata.csv",
        manifest_file="visualizations/bioapp_experiment/bioapp_phase2c_endpoint_freeze_with_validated_cta_to_spot_mapping/endpoint_freeze_manifest.json", manifest_row="entire manifest",
        manuscript_location="Fig. 5A-I; formal LaTeX file absent", redistribution_status="Raw data via SND DOI; processed data openly deposited in Zenodo; project subset/endpoint derivative planned separately",
        project_generated_derivative="PROJECT_GENERATED 4,014-cell balanced subset and CTA-to-spot frozen endpoint", status="PARTIAL",
        evidence="CTA paper and Data availability; local frozen endpoint manifest; CytoSPACE input manifest; Wu GEO record and cell IDs",
        notes="CID45171 is confirmed as Wu sample GSM5354535. BCSA2 is HER2-positive. The exact Visium chemistry/version and slide/capture-area identifiers remain unresolved. Spatial and scRNA donors are unmatched.",
    ))

    external_specs = [
        ("EXT01_MEL1", "cytospace_fig2c_melanoma_mel1_rep2_profile_mask_endothelial_cells", "ST_mel1_rep2", "Endothelial cells"),
        ("EXT02_MEL2", "cytospace_fig2c_melanoma_mel2_rep1_profile_mask_nk_cells", "ST_mel2_rep1", "NK cells"),
        ("EXT03_BRCA_V1", "cytospace_fig2d_tme_brca_er_her2_fresh_frozen_profile_mask_t_cells", "V1_Breast_Cancer_Block_A_Section_1", "T cells"),
        ("EXT04_BRCA_FFPE", "cytospace_fig2d_tme_brca_her2_ffpe_profile_mask_t_cells", "Visium_FFPE_Human_Breast_Cancer", "T cells"),
        ("EXT05_TNBC", "cytospace_fig2d_tme_brca_tnbc_fresh_frozen_profile_mask_fibroblasts", "CID4465", "Fibroblasts"),
        ("EXT06_CRC", "cytospace_fig2d_tme_crc_fresh_frozen_profile_mask_monocytes_and_macrophages", "Parent_Visium_Human_ColorectalCancer", "Monocytes and Macrophages"),
    ]
    for record_id, internal_id, sample, target in external_specs:
        source = public_row(internal_id)
        is_melanoma = "melanoma" in internal_id
        if sample == "ST_mel2_rep1":
            melanoma_input = "data/raw/cytospace_fig2d_tme/melanoma_ST_Thrane_legacyST/berglund2018spatial_ST_mel2_rep1.h5ad"
        else:
            melanoma_input = f"data/raw/cytospace_fig2c_melanoma/berglund2018spatial_{sample}.h5ad"
        records.append(base_record(
            record_id=record_id, dataset_family="downstream external tumour pairings", experiment_id=internal_id, internal_dataset_id=internal_id,
            data_role="paired spatial/reference input plus EcoTyper evaluation resource", formal_dataset_name=source["formal_dataset_name"], formal_release_name=source["formal_dataset_name"],
            source_study_title=source["spatial_source_title"], species=source["species"], tissue=source["tissue"], disease=source["disease"], sample_id=sample,
            spatial_platform=source["spatial_platform"], assay="spatial transcriptomics", tissue_processing=source["sample_preparation"], reference_type="scRNA-seq",
            reference_dataset_name=source["scrna_reference_name"], reference_sample_subset=source["scrna_sample_id"], database="official provider/repository records",
            geo_accession=source["scrna_accession"] if source["scrna_accession"].startswith("GSE") else "NOT_APPLICABLE",
            official_source_url=source["official_source_urls"], source_publication=f"{source['spatial_source_title']}; {source['scrna_source_title']}", bibtex_key=source["verified_citation_keys"],
            input_file=f"data/processed/cytospace_fig2d_tme/{internal_id.replace('_profile_mask_' + target.lower().replace(' ', '_'), '')}/stage1_preprocess/exported/st_expression_normalized.csv" if not is_melanoma else f"{melanoma_input};data/raw/cytospace_fig2c_melanoma/GSE72056_melanoma_single_cell_revised_v2.txt",
            config_file=f"configs/datasets/{internal_id}.yaml", manuscript_location="Fig. 2C/D/F; formal LaTeX file absent", redistribution_status="See source provider terms",
            project_generated_derivative=f"profile masking target={target}", status="CONFIRMED", evidence=source["repository_evidence_paths"] + ";" + source["official_source_urls"], notes=source["unresolved_fields"],
        ))

    records.extend([
        base_record(
            record_id="RES01_ECOTYPER", dataset_family="downstream resources", experiment_id="multiple Fig. 2/3 enrichment experiments", internal_dataset_id="EcoTyper_state_programs",
            data_role="gene-program/evaluation resource", formal_dataset_name="EcoTyper carcinoma ecotype and cell-state programs", formal_release_name="EcoTyper",
            source_study_title="Archetypes of human tumor microenvironment cell states and ecosystems", species="Homo sapiens", tissue="multiple carcinomas", disease="cancer",
            reference_type="gene-program resource", reference_dataset_name="EcoTyper CE9/CE10 and state programs", database="publication supplement/software resource",
            official_source_url="https://doi.org/10.1016/j.cell.2021.09.014", source_publication="Luca et al., Cell (2021)", bibtex_key="Luca2021EcoTyper",
            manuscript_location="Fig. 2C/F and Fig. 3C; formal LaTeX file absent", redistribution_status="See EcoTyper terms", status="CONFIRMED",
            evidence="targeted/all-candidate manifests and original paper", notes="Evaluation resource, not a spatial or single-cell input dataset.",
        ),
        base_record(
            record_id="RES02_LIANA_OMNIPATH", dataset_family="downstream resources", experiment_id="stage3b_communication_validation", internal_dataset_id="liana_1.7.3_omni_resource",
            data_role="ligand-receptor interaction resource", formal_dataset_name="LIANA consensus resource derived from OmniPath", formal_release_name="LIANA 1.7.3 omni_resource export",
            source_study_title="Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data", species="Homo sapiens", tissue="multiple", disease="multiple",
            reference_type="ligand-receptor prior", reference_dataset_name="LIANA consensus/OmniPath resource", database="LIANA/OmniPath software resources",
            official_source_url="https://doi.org/10.1038/s41467-022-30755-0;https://doi.org/10.1038/s41592-019-0409-6", source_publication="Dimitrov et al. (2022); Turei et al. (2021)", bibtex_key="Dimitrov2022LIANA;TureiOmniPath",
            input_file="data/raw/spatial_communication_reference_resources/liana_1.7.3_omni_resource.csv", manuscript_location="Fig. 4I communication validation; formal LaTeX file absent",
            redistribution_status="Derived resource; original resource terms apply", status="CONFIRMED", evidence="local frozen CSV and LIANA/OmniPath publications",
            notes="This resource materially contributes ligand-receptor pairs; it is not biological ground truth.",
        ),
        base_record(
            record_id="RES03_CYTOSPACE_SUPP_ORDERING", dataset_family="downstream resources", experiment_id="state-decoy experiments", internal_dataset_id="CytoSPACE_supplementary_tables_S8_S9",
            data_role="secondary numerical transcription/evaluation resource", formal_dataset_name="CytoSPACE Supplementary Tables S8/S9", formal_release_name="CytoSPACE supplementary information",
            source_study_title="CytoSPACE: a method for mapping single-cell transcriptomic data to spatial expression profiles", species="NOT_APPLICABLE", tissue="NOT_APPLICABLE", disease="NOT_APPLICABLE",
            reference_type="secondary transcription resource", reference_dataset_name="State-ordering values transcribed from source literature", database="journal supplementary material",
            official_source_url="https://doi.org/10.1038/s41587-023-01697-9", source_publication="Vahid et al., Nature Biotechnology (2023)", bibtex_key="Vahid2023CytoSPACE",
            manuscript_location="Fig. 3E-F; formal LaTeX file absent", redistribution_status="Publication supplement terms apply", status="CONFIRMED",
            evidence="state-decoy configs/reports and CytoSPACE publication supplement", notes="Not an original dataset and must not replace Zheng or Ransick citations.",
        ),
    ])

    for row in records:
        if row["input_file"] != "NOT_APPLICABLE":
            row["input_file_sha1"] = sha1_list(row["input_file"])
        if row["config_file"] != "NOT_APPLICABLE":
            row["config_file_sha1"] = sha1_list(row["config_file"])
        if row["status"] not in STATUSES:
            raise ValueError(f"Invalid status for {row['record_id']}: {row['status']}")
    return records


def build_crosswalk() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []

    def add(exp: str, panel: str, section: str, record: str, *, evaluation: str = "NOT_APPLICABLE", profile: str = "NOT_APPLICABLE", dropout: str = "NOT_APPLICABLE", derivative: str = "NOT_APPLICABLE", status: str = "CONFIRMED", notes: str = "") -> None:
        rows.append({
            "experiment_id": exp, "figure_panel": panel, "manuscript_result_section": section,
            "spatial_input_record_id": record, "reference_input_record_id": record,
            "evaluation_resource_record_id": evaluation, "profile_mask_target": profile,
            "reference_dropout_target": dropout, "project_derivative": derivative,
            "formal_status": status, "notes": notes,
        })

    sim = [
        ("real_brca7_endothelial_marker_control_sc_missing_endothelial_cells", "SIM01_BRCA", "no ST-profile masking; SC endothelial missing", "CONFIRMED"),
        ("real_brca7_endothelial_marker_missing_epithelial_cells_sc_missing_endothelial_cells", "SIM01_BRCA", "epithelial profile masking; SC endothelial missing", "CONFIRMED"),
        ("real_brca7_endothelial_marker_missing_epithelial_cells_pcs_sc_missing_endothelial_cells", "SIM01_BRCA", "epithelial+PC profile masking; SC endothelial missing", "CONFIRMED"),
        ("mouse_brain_refined7_balanced_clustered_sim_sc_missing_ext_l56", "SIM02_MOUSE_BRAIN", "no ST-profile masking; SC Ext_L56 missing", "CONFIRMED"),
        ("mouse_brain_refined7_balanced_clustered_sim_missing_micro_fill_inh_pvalb_sc_missing_ext_l56", "SIM02_MOUSE_BRAIN", "Microglia profile masking; SC Ext_L56 missing", "CONFIRMED"),
        ("mouse_brain_refined7_balanced_clustered_sim_missing_micro_oligo_2_fill_inh_pvalb_sc_missing_ext_l56", "SIM02_MOUSE_BRAIN", "Microglia+Oligo_2 profile masking; SC Ext_L56 missing", "CONFIRMED"),
        ("human_lung_5loc_fine9_clustered_sim_sc_missing_b_cell", "SIM03_HUMAN_LUNG", "no ST-profile masking; SC B cell missing", "CONFLICT"),
        ("human_lung_5loc_fine9_clustered_sim_missing_at2_sc_missing_b_cell", "SIM03_HUMAN_LUNG", "AT2 profile masking; SC B cell missing", "CONFLICT"),
        ("human_lung_5loc_fine9_clustered_sim_missing_at2_fibroblast_sc_missing_b_cell", "SIM03_HUMAN_LUNG", "AT2+Fibroblast profile masking; SC B cell missing", "CONFLICT"),
    ]
    for exp, record, derivative, status in sim:
        add(exp, "Fig. 1C-D", "Joint Stage3A-Stage3B simulations", record, derivative=derivative, status=status)

    lowres = [
        ("adult_mouse_kidney_real_profile_mask_endo", "LR01_KIDNEY", "Endo", "CONFIRMED"),
        ("ffpe_mouse_brain_sagittal_real_profile_mask_microglia", "LR02_BRAIN", "Microglia", "CONFIRMED"),
        ("human_breast_cancer_real_profile_mask_basal_cell", "LR03_BRCA_FFPE", "Basal cell", "PARTIAL"),
        ("human_breast_cancer_visium_ff_wta_real_profile_mask_macrophage", "LR04_BRCA_FF", "Macrophage", "PARTIAL"),
        ("human_breast_cancer_wta_120_real_profile_mask_endothelial_cell", "LR05_BRCA_ILC", "Endothelial cell", "PARTIAL"),
        ("human_cervical_cancer_real_profile_mask_epithelial_cell", "LR06_CERVICAL", "Epithelial cell", "PARTIAL"),
        ("human_heart_ff_real_profile_mask_endothelial_cell", "LR07_HEART", "Endothelial cell", "CONFIRMED"),
        ("human_intestine_cancer_real_profile_mask_endothelial_cell", "LR08_INTESTINE", "Endothelial cell", "PARTIAL"),
        ("human_lymph_node_real_profile_mask_b_cell", "LR09_LYMPH_NODE", "B cell", "CONFIRMED"),
        ("mouse_embryo_real_profile_mask_erythroid", "LR10_EMBRYO", "Erythroid", "CONFIRMED"),
    ]
    for exp, record, target, status in lowres:
        add(exp, "Fig. 2A-B/E", "Low-resolution profile masking", record, profile=target, derivative="project-generated masked profile", status=status)

    merscope = [
        ("highres_humanbreastcancerpatient1_profile_mask_monocytes_and_macrophages", "MERS01_BREAST", "Monocytes and Macrophages"),
        ("highres_humancoloncancerpatient1_profile_mask_fibroblasts", "MERS02_COLON", "Fibroblasts"),
        ("highres_humanlungcancerpatient1_profile_mask_plasma_cells", "MERS03_LUNG", "Plasma cells"),
        ("highres_humanmelanomapatient1_profile_mask_fibroblasts", "MERS04_MELANOMA1", "Fibroblasts"),
        ("highres_humanmelanomapatient2_profile_mask_b_cells", "MERS05_MELANOMA2", "B cells"),
    ]
    for exp, record, target in merscope:
        add(exp, "Fig. 3A-D", "High-resolution profile masking", record, evaluation="RES01_ECOTYPER", profile=target, derivative="project-generated disjoint same-assay reference split", status="PARTIAL")

    add("cytospace_fig2k_breast_state_decoy_stage3_detected", "Fig. 3E", "T-cell state-decoy experiment", "STATE01_TCELL_DECOY", evaluation="RES03_CYTOSPACE_SUPP_ORDERING", derivative="engineered T-cell-state decoy", status="PROJECT_GENERATED")
    add("cytospace_fig2i_mouse_kidney_stage3_unsupported_decoy_n1000", "Fig. 3F", "Kidney State32 decoy experiment", "STATE02_KIDNEY_STATE32", evaluation="RES03_CYTOSPACE_SUPP_ORDERING", derivative="1,000 State32-like decoys", status="PROJECT_GENERATED")

    dropout = [
        ("mouse_embryo_real_sc_missing_endoderm_gut", "RD01_EMBRYO_ENDODERM", "Endoderm/Gut"),
        ("mouse_embryo_real_sc_missing_erythroid", "RD02_EMBRYO_ERYTHROID", "Erythroid"),
        ("cytospace_fig2d_tme_brca_tnbc_fresh_frozen_sc_missing_plasma_cells", "RD03_TNBC_PLASMA", "Plasma cells"),
        ("cytospace_fig2d_tme_crc_fresh_frozen_sc_missing_b_cells", "RD04_CRC_B", "B cells"),
        ("cytospace_fig2d_tme_brca_her2_ffpe_sc_missing_plasma_cells", "RD05_BRCA_FFPE_PLASMA", "Plasma cells"),
        ("cytospace_fig2d_tme_brca_her2_ffpe_sc_missing_epithelial_cells", "RD06_BRCA_FFPE_EPITHELIAL", "Epithelial cells"),
    ]
    for exp, record, target in dropout:
        evaluation = "RES02_LIANA_OMNIPATH" if record == "RD05_BRCA_FFPE_PLASMA" else "NOT_APPLICABLE"
        add(exp, "Fig. 4F-G/I", "Real reference dropout", record, evaluation=evaluation, dropout=target, derivative="project-generated reference-dropout copy")

    add("cell2location_ST8059051_thalamic_excitatory_reference_missing", "Fig. 4A-E/H", "Detailed thalamic-excitatory case", "THAL01_ST8059051", dropout="thalamic excitatory states", derivative="reference-dropout copy and marker-region proxy")
    add("bioapp_BCSA2TumB1_CTA_immune_endpoint", "Fig. 5A-I", "CTA biological application", "CTA01_BCSA2TUMB1", derivative="4,014-cell balanced reference subset and frozen CTA immune endpoint", status="PARTIAL")

    external = [
        ("cytospace_fig2c_melanoma_mel1_rep2_profile_mask_endothelial_cells", "EXT01_MEL1", "Endothelial cells"),
        ("cytospace_fig2c_melanoma_mel2_rep1_profile_mask_nk_cells", "EXT02_MEL2", "NK cells"),
        ("cytospace_fig2d_tme_brca_er_her2_fresh_frozen_profile_mask_t_cells", "EXT03_BRCA_V1", "T cells"),
        ("cytospace_fig2d_tme_brca_her2_ffpe_profile_mask_t_cells", "EXT04_BRCA_FFPE", "T cells"),
        ("cytospace_fig2d_tme_brca_tnbc_fresh_frozen_profile_mask_fibroblasts", "EXT05_TNBC", "Fibroblasts"),
        ("cytospace_fig2d_tme_crc_fresh_frozen_profile_mask_monocytes_and_macrophages", "EXT06_CRC", "Monocytes and Macrophages"),
    ]
    for exp, record, target in external:
        add(exp, "Fig. 2C/F", "External tumour state-enrichment pairings", record, evaluation="RES01_ECOTYPER", profile=target)
    return rows


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_bib() -> dict[str, dict[str, str]]:
    files = [
        ROOT / "visualizations/manuscript_public_source_resolution/proposed_verified_bibliography_entries.bib",
        ROOT / "manuscript/reference_audit/dataset_reference_candidates.bib",
    ]
    records: dict[str, dict[str, str]] = {}
    for path in files:
        text = path.read_text(encoding="utf-8")
        for match in re.finditer(r"@\w+\{([^,]+),(.*?)(?=\n\})", text, flags=re.S):
            key, body = match.groups()
            fields = {}
            for name in ("title", "author", "year", "doi"):
                hit = re.search(rf"\b{name}\s*=\s*\{{(.*?)\}}", body, flags=re.S | re.I)
                fields[name] = re.sub(r"\s+", " ", hit.group(1)).strip() if hit else "UNRESOLVED"
            records.setdefault(key.strip(), fields)
    return records


def build_bibliography(records: list[dict[str, str]]) -> list[dict[str, str]]:
    parsed = parse_bib()
    links: dict[str, list[str]] = defaultdict(list)
    for row in records:
        for key in row["bibtex_key"].split(";"):
            key = key.strip()
            if key and key not in STATUSES:
                links[key].append(row["record_id"])
    output = []
    for key, record_ids in sorted(links.items()):
        meta = parsed.get(key, {})
        output.append({
            "bibtex_key": key,
            "publication_title": meta.get("title", "UNRESOLVED"),
            "authors": meta.get("author", "UNRESOLVED"),
            "year": meta.get("year", "UNRESOLVED"),
            "doi": meta.get("doi", "UNRESOLVED"),
            "dataset_role": "source publication, provider release, or evaluation-resource citation",
            "linked_record_ids": ";".join(sorted(set(record_ids))),
            "supports_dataset_claim": "yes, with role distinctions in notes",
            "status": "PARTIAL",
            "notes": "Candidate entry exists only in audit bibliography. Formal sn-bibliography.bib was not found and was not modified.",
        })
    return output


def build_occurrences() -> list[dict[str, str]]:
    rows = [{
        "section": "repository-wide search", "subsection": "formal manuscript", "figure_panel": "Fig. 1-5",
        "quoted_context": "No .tex file found in repository/workspace; no formal manuscript could be scanned.",
        "dataset_or_identifier": "FORMAL_LATEX_NOT_FOUND", "identifier_type": "formal manuscript source",
        "linked_record_id": "UNRESOLVED", "status": "UNRESOLVED",
        "required_action": "Provide the actual manuscript compilation directory and repeat this occurrence scan before submission.",
        "notes": "Markdown preparation files and PDFs are not substitutes for the formal LaTeX source.",
    }]
    support_files = ["论文材料准备.md", "论文材料准备B.md", "生物学应用实验论文材料准备.md", "README.md"]
    tokens = {
        "GSE176078": "CTA01_BCSA2TUMB1", "BCSA2TumB1": "CTA01_BCSA2TUMB1",
        "E-MTAB-11114": "THAL01_ST8059051", "E-MTAB-11115": "THAL01_ST8059051",
        "ERS20065156": "SIM03_HUMAN_LUNG", "HumanBreastCancerPatient1": "MERS01_BREAST",
        "V11J26-008": "LR03_BRCA_FFPE", "V52B25-081": "SIM02_MOUSE_BRAIN",
    }
    for file_name in support_files:
        path = ROOT / file_name
        if not path.is_file():
            continue
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            for token, record_id in tokens.items():
                if token in line:
                    rows.append({
                        "section": file_name, "subsection": f"line {number}", "figure_panel": "supporting material only",
                        "quoted_context": line.strip()[:400], "dataset_or_identifier": token, "identifier_type": "dataset/sample/accession token",
                        "linked_record_id": record_id, "status": "PARTIAL",
                        "required_action": "Recheck against formal LaTeX when supplied.",
                        "notes": "Occurrence is in a preparation Markdown file, not the formal manuscript.",
                    })
    return rows


def build_accession_audit() -> list[dict[str, str]]:
    specs = [
        ("GSE176078", "dataset/reference accession", "NCBI GEO Series accession", "SIM01_BRCA;CTA01_BCSA2TUMB1;RD03_TNBC_PLASMA;RD05_BRCA_FFPE_PLASMA;RD06_BRCA_FFPE_EPITHELIAL", "CONFIRMED", "GEO Series GSE176078"),
        ("GSM5354535", "Wu patient identifier CID45171", "NCBI GEO Sample accession for CID45171", "CTA01_BCSA2TUMB1", "CONFIRMED", "CID45171 (GEO Sample GSM5354535)"),
        ("E-MTAB-11114", "spatial accession", "BioStudies/ArrayExpress study accession", "THAL01_ST8059051", "CONFIRMED", "spatial study E-MTAB-11114"),
        ("E-MTAB-11115", "reference accession", "BioStudies/ArrayExpress study accession", "SIM02_MOUSE_BRAIN;LR02_BRAIN;THAL01_ST8059051", "CONFIRMED", "snRNA-seq reference study E-MTAB-11115"),
        ("E-MTAB-11640", "study/reference record", "BioStudies/ArrayExpress study accession", "SIM03_HUMAN_LUNG", "CONFIRMED", "study E-MTAB-11640"),
        ("PRJEB52292", "study record", "ENA project accession", "SIM03_HUMAN_LUNG", "CONFIRMED", "ENA project PRJEB52292"),
        ("ERS20065156", "run", "ENA sample accession", "SIM03_HUMAN_LUNG", "CONFLICT", "sample accession ERS20065156; do not call it a run"),
        ("SAMEA115633909", "BioSample", "BioSample accession linked to ERS20065156", "SIM03_HUMAN_LUNG", "CONFIRMED", "BioSample SAMEA115633909"),
        ("GSE129798", "kidney reference", "NCBI GEO Series accession", "STATE02_KIDNEY_STATE32", "CONFIRMED", "GEO Series GSE129798"),
        ("SRP192559", "kidney reference", "SRA study accession", "STATE02_KIDNEY_STATE32", "CONFIRMED", "SRA study SRP192559"),
        ("PRJNA532850", "kidney reference", "NCBI BioProject accession", "STATE02_KIDNEY_STATE32", "CONFIRMED", "BioProject PRJNA532850"),
        ("GSE132465", "CRC reference", "NCBI GEO Series accession", "RD04_CRC_B;EXT06_CRC", "CONFIRMED", "GEO Series GSE132465"),
        ("10.5281/zenodo.15211538", "CTA accession", "Zenodo dataset DOI", "CTA01_BCSA2TUMB1", "CONFIRMED", "processed-data DOI 10.5281/zenodo.15211538"),
        ("10.48723/f4v5-m008", "CTA accession", "Swedish National Data Service raw-data DOI", "CTA01_BCSA2TUMB1", "CONFIRMED", "raw-data DOI 10.48723/f4v5-m008"),
        ("Researchdata.se 2025-97", "CTA accession", "Researchdata.se catalogue identifier", "CTA01_BCSA2TUMB1", "CONFIRMED", "Researchdata.se dataset 2025-97"),
        ("V11J26-008", "dataset accession", "10x slide serial", "SIM01_BRCA;LR03_BRCA_FFPE;RD05_BRCA_FFPE_PLASMA;RD06_BRCA_FFPE_EPITHELIAL", "CONFIRMED", "slide V11J26-008, capture area B1; not a repository accession"),
        ("HumanBreastCancerPatient1", "patient ID", "Vizgen provider release sample identifier", "MERS01_BREAST;STATE01_TCELL_DECOY", "PARTIAL", "release sample identifier HumanBreastCancerPatient1; do not imply verified clinical identity"),
        ("ST8059051", "accession", "cell2location spatial section identifier", "THAL01_ST8059051", "CONFIRMED", "section ST8059051 (Visium-29B; C05717-021/B1)"),
    ]
    return [dict(zip(["accession", "claimed_type_in_manuscript", "verified_type", "linked_record_id", "status", "recommended_wording"], spec)) for spec in specs]


def build_unresolved() -> list[dict[str, str]]:
    issues = [
        ("U001", "FORMAL_ARTIFACTS", "formal manuscript .tex", "Formal LaTeX source was not found.", "repository, parent workspace, and prior bibliography inventory", "Provide actual compilation directory and rerun occurrence scan.", "CRITICAL"),
        ("U002", "FORMAL_ARTIFACTS", "formal sn-bibliography.bib", "Only candidate audit BibTeX files exist.", "repository and prior full-user-profile inventory", "Provide the compiled bibliography and run key-level conflict audit.", "CRITICAL"),
        ("U003", "FORMAL_ARTIFACTS", "formal dataset manifest", "No unified formal dataset manifest is present; dispersed experiment manifests are not equivalent.", "repository manifest/provenance/accession search", "Create and freeze a formal manifest only after this audit is reviewed.", "CRITICAL"),
        ("U004", "LR03_BRCA_FFPE;LR04_BRCA_FF;LR05_BRCA_ILC", "historical hECA breast export lineage", "The current public hECA project export and source study are identified, but the historical local file lacks a preserved byte-for-byte merge/checksum manifest.", "local h5ad metadata, current hECA/Zenodo records, prior source audit", "Recover the original export/merge log or describe the reference as a curated hECA export without overclaiming exact component provenance.", "MAJOR"),
        ("U005", "LR06_CERVICAL", "healthy uterus composite contributing cells", "Source studies are documented, but exact contributing-cell/accession membership of the local historical export is incomplete.", "local h5ad, hECA metadata, prior source audit", "Recover a cell-level source manifest.", "MAJOR"),
        ("U006", "LR08_INTESTINE", "healthy intestine composite contributing cells", "Source studies are documented, but exact contributing-cell/accession membership of the local historical export is incomplete.", "local h5ad, hECA metadata, prior source audit", "Recover a cell-level source manifest.", "MAJOR"),
        ("U007", "MERS01_BREAST;MERS02_COLON;MERS03_LUNG;MERS04_MELANOMA1;MERS05_MELANOMA2", "Vizgen redistribution and stable per-sample URL", "Provider release and release sample IDs are confirmed, but repository lacks preserved license text and stable per-sample download URLs/checksums.", "Vizgen provider page, raw local files, high-resolution manifest", "Archive provider terms and original download receipt/URL/checksum before Data availability is finalized.", "MAJOR"),
        ("U008", "CTA01_BCSA2TUMB1", "Visium chemistry/version and slide/capture area", "CTA paper confirms frozen Visium sections but does not expose the exact chemistry version or slide serial/capture area for BCSA2TumB1 in the evidence currently present.", "CTA paper, processed data, local endpoint manifests", "Recover sample sheet or raw Space Ranger metadata.", "MAJOR"),
        ("U009", "RD03_TNBC_PLASMA", "CID4465 slide and capture area", "CID4465 identity and Zenodo source are confirmed, but slide serial/capture area are not exposed.", "Zenodo 4739739, local data, prior source audit", "Leave fields unresolved unless primary metadata is recovered.", "MINOR"),
        ("U010", "RD04_CRC_B", "colorectal slide serial", "Formal 10x release and capture area C1 are confirmed, but slide serial is not.", "10x dataset page, local input, prior source audit", "Leave slide serial unresolved unless provider metadata is recovered.", "MINOR"),
        ("U011", "SIM03_HUMAN_LUNG", "run accession", "ERS20065156 is a sample accession and no linked run accession was recovered from the ENA read-run query.", "ENA sample XML and read-run query", "Do not invent a run accession; report project, study, sample and BioSample identifiers only.", "MAJOR"),
    ]
    fields = ["issue_id", "linked_record_id", "field", "issue", "searched_evidence", "required_resolution", "priority"]
    return [dict(zip(fields, item)) for item in issues]


def build_conflicts() -> list[dict[str, str]]:
    return [
        {
            "conflict_id": "C001", "linked_record_id": "SIM03_HUMAN_LUNG", "field": "ERS20065156 accession type",
            "evidence_a": "Task/manuscript claim describes ERS20065156 as a run.",
            "evidence_b": "ENA XML returns SAMPLE accession=ERS20065156 and external BioSample SAMEA115633909; no read run was returned.",
            "recommended_canonical_value": "ERS20065156 = ENA sample accession; SAMEA115633909 = BioSample; run accession unresolved.",
            "status": "CONFLICT", "required_action": "Correct future manuscript/manifest wording after user review.",
        },
        {
            "conflict_id": "C002", "linked_record_id": "FORMAL_ARTIFACTS", "field": "submission readiness",
            "evidence_a": "Earlier repository-only audit issued PASS for a narrower internal inventory.",
            "evidence_b": "No formal LaTeX, formal bibliography, or unified formal dataset manifest is present for cross-file consistency testing.",
            "recommended_canonical_value": "Task 4A decision HOLD until formal artifacts are supplied and audited.",
            "status": "CONFLICT", "required_action": "Do not reuse the earlier PASS as a submission-level provenance decision.",
        },
    ]


def build_updates() -> list[dict[str, str]]:
    specs = [
        ("formal manuscript .tex", "Methods/Data availability", "future-tense public manifest language", "No formal manifest exists and formal LaTeX is absent.", "FORMAL_ARTIFACTS", "Replace only after a reviewed manifest is frozen; provide compilation source first.", "CRITICAL"),
        ("formal manuscript .tex", "Human-lung dataset description", "ERS20065156 as run", "Wrong accession type.", "SIM03_HUMAN_LUNG", "Call ERS20065156 an ENA sample accession and SAMEA115633909 its BioSample; omit run accession.", "CRITICAL"),
        ("formal manuscript .tex", "Fig. 3 datasets", "Human...Patient1 as patient IDs", "These are provider release sample identifiers, not verified clinical IDs or accessions.", "MERS01_BREAST", "Use 'Vizgen release sample identifier' and separate MERSCOPE platform from MERFISH assay.", "MAJOR"),
        ("formal manuscript .tex", "Low-resolution references", "hECA/Tabula Sapiens as a single exact source", "Historical local export lineage is incomplete.", "LR03_BRCA_FFPE", "Describe curated hECA project exports and list only verified contributing source studies.", "MAJOR"),
        ("formal manuscript .tex", "Fig. 1 mouse-brain simulation", "mouse-brain spatial study phrasing", "10x FFPE CytAssist scaffold is not ST8059051/E-MTAB-11114.", "SIM02_MOUSE_BRAIN", "Name V52B25-081/B provider scaffold and E-MTAB-11115 reference; distinguish from detailed thalamic case.", "MAJOR"),
        ("formal manuscript .tex", "Fig. 4 detailed case", "accuracy", "Evaluation regions are reference-derived marker proxies.", "THAL01_ST8059051", "Use reference-relative recovery/agreement terminology unless an independent truth is supplied.", "MAJOR"),
        ("formal manuscript .tex", "Fig. 4 communication", "coupling absolute error", "LIANA/OmniPath coupling is relative to a reference mapping, not biological truth.", "RES02_LIANA_OMNIPATH", "Use reference-relative coupling deviation and cite LIANA/OmniPath.", "MAJOR"),
        ("formal manuscript .tex", "Fig. 5 data and endpoint", "CTA endpoint/reference description", "Project derivatives and donor mismatch need explicit disclosure.", "CTA01_BCSA2TUMB1", "State that the endpoint is image-derived CTA converted/frozen at spot level; 4,014-cell Wu subset is project-generated and donor-unmatched.", "MAJOR"),
        ("formal manuscript .tex", "Data availability", "MERSCOPE accession", "No repository accession was found.", "MERS01_BREAST", "Cite the official Vizgen provider release URL; do not invent an accession; disclose access/redistribution terms once archived.", "MAJOR"),
        ("sn-bibliography.bib", "dataset-related keys", "candidate keys", "Formal bibliography absent; key conflicts cannot be tested.", "FORMAL_ARTIFACTS", "Merge only after the compiled BibTeX is supplied and each key is checked.", "CRITICAL"),
    ]
    fields = ["file", "section", "current_text", "problem", "canonical_record_id", "recommended_replacement", "priority"]
    return [dict(zip(fields, item)) for item in specs]


def build_manifest_audit(records: list[dict[str, str]]) -> list[dict[str, str]]:
    rows = []
    for record in records:
        rows.append({
            "record_id": record["record_id"],
            "manuscript_value": "UNRESOLVED: formal LaTeX absent",
            "manifest_value": record["manifest_file"] if record["manifest_file"] != "NOT_APPLICABLE" else "formal dataset manifest not yet present",
            "config_value": record["config_file"],
            "input_metadata_value": record["input_file"],
            "official_source_value": record["formal_dataset_name"],
            "consistency_status": "CONFLICT" if record["status"] == "CONFLICT" else ("PARTIAL" if record["status"] in {"PARTIAL", "UNRESOLVED"} else "CONFIRMED"),
            "recommended_canonical_value": f"{record['formal_dataset_name']} | {record['sample_id']} | {record['slide_id']}/{record['capture_area']}",
            "notes": "Dispersed experiment manifests were audited as evidence but are not treated as the formal dataset manifest.",
        })
    return rows


def build_source_registry(records: list[dict[str, str]]) -> list[dict[str, str]]:
    rows = []
    for record in records:
        urls = [url for url in record["official_source_url"].split(";") if url and url not in STATUSES]
        if not urls:
            continue
        for index, url in enumerate(urls, 1):
            authority = "Official DOI/publisher"
            if "10xgenomics" in url:
                authority = "10x Genomics"
            elif "ncbi.nlm.nih.gov" in url:
                authority = "NCBI GEO/SRA"
            elif "ebi.ac.uk" in url:
                authority = "EMBL-EBI"
            elif "vizgen" in url:
                authority = "Vizgen"
            elif "zenodo" in url:
                authority = "Zenodo"
            rows.append({
                "record_id": record["record_id"], "source_authority": authority,
                "source_type": "official database/provider/publisher page", "url": url,
                "accessed_date": ACCESSED, "supports_fields": "dataset identity, platform, sample/accession, publication, or access route",
                "status": "PARTIAL" if record["status"] == "PARTIAL" else "CONFIRMED",
                "notes": f"Source {index} of {len(urls)} for this record.",
            })
    return rows


def build_availability(records: list[dict[str, str]]) -> list[dict[str, str]]:
    rows = []
    for record in records:
        if record["status"] == "PROJECT_GENERATED" or record["project_generated_derivative"].startswith("PROJECT_GENERATED"):
            category = "project-generated derivative; source-data terms apply"
        elif "provider" in record["database"].lower() or "10x" in record["database"].lower() or "vizgen" in record["database"].lower():
            category = "provider-hosted release"
        else:
            category = "public repository/publisher source"
        accession = ";".join(value for field in ("study_accession", "project_accession", "run_accession", "sample_accession", "biosample_accession", "geo_accession", "arrayexpress_accession") if (value := record[field]) not in STATUSES)
        rows.append({
            "record_id": record["record_id"], "formal_dataset_name": record["formal_dataset_name"],
            "source_database": record["database"], "accession": accession or "NOT_APPLICABLE",
            "official_url": record["official_source_url"], "download_or_access_method": category,
            "redistributable": "UNRESOLVED" if "not preserved" in record["redistribution_status"].lower() else record["redistribution_status"],
            "project_processed_data_planned": "yes for project-generated derivatives; location/license not yet frozen",
            "required_statement": f"Name the formal source and exact identifier; classify as {category}.",
            "unresolved_issue": record["notes"] if record["status"] in {"PARTIAL", "CONFLICT", "UNRESOLVED"} else "none",
        })
    return rows


def coverage(crosswalk: list[dict[str, str]]) -> dict[str, dict[str, int]]:
    groups = {
        "Joint simulations": [row for row in crosswalk if row["manuscript_result_section"] == "Joint Stage3A-Stage3B simulations"],
        "Low-resolution profile masking": [row for row in crosswalk if row["manuscript_result_section"] == "Low-resolution profile masking"],
        "MERSCOPE": [row for row in crosswalk if row["manuscript_result_section"] == "High-resolution profile masking"],
        "State-decoy experiments": [row for row in crosswalk if "decoy experiment" in row["manuscript_result_section"]],
        "Real reference dropout": [row for row in crosswalk if row["manuscript_result_section"] == "Real reference dropout"],
        "Detailed Thalamic case": [row for row in crosswalk if row["manuscript_result_section"] == "Detailed thalamic-excitatory case"],
        "CTA biological application": [row for row in crosswalk if row["manuscript_result_section"] == "CTA biological application"],
    }
    result = {}
    for name, rows in groups.items():
        counts = Counter(row["formal_status"] for row in rows)
        result[name] = {
            "expected": {"Joint simulations": 9, "Low-resolution profile masking": 10, "MERSCOPE": 5, "State-decoy experiments": 2, "Real reference dropout": 6, "Detailed Thalamic case": 1, "CTA biological application": 1}[name],
            "audited": len(rows), "confirmed": counts["CONFIRMED"], "partial": counts["PARTIAL"],
            "unresolved": counts["UNRESOLVED"], "conflict": counts["CONFLICT"], "project_generated": counts["PROJECT_GENERATED"],
        }
    return result


def build_readme(summary: dict) -> str:
    return f"""# Task 4A dataset provenance formal audit

This directory is an independent, read-only provenance audit generated on {ACCESSED}. It does not replace a formal dataset manifest.

## Decision

**{summary['decision']}**

The repository contains extensive experiment-level provenance, but the formal manuscript `.tex`, compiled `sn-bibliography.bib`, and unified formal dataset manifest were not present. A submission-level manuscript/BibTeX/manifest consistency check is therefore impossible. The human-lung record also contains a material accession-type conflict: `ERS20065156` is an ENA sample accession, not a run accession.

## Scope

- {summary['counts']['master_records']} provenance master records.
- {summary['counts']['experiment_crosswalk_rows']} experiment crosswalk rows.
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
"""


def build_final(summary: dict, unresolved: list[dict[str, str]], conflicts: list[dict[str, str]]) -> str:
    lines = [
        "# Final decision", "", "## 1. Decision", "", f"**{summary['decision']}**", "",
        "The project is on HOLD for submission-level dataset provenance. Core source identities are mostly recoverable, but the formal manuscript, formal bibliography and unified formal dataset manifest are absent, and one accession-type conflict remains.", "",
        "## 2. Audit coverage", "",
        "| Family | Expected | Audited | Confirmed | Partial | Project-generated | Unresolved | Conflict |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for name, values in summary["coverage"].items():
        lines.append(f"| {name} | {values['expected']} | {values['audited']} | {values['confirmed']} | {values['partial']} | {values['project_generated']} | {values['unresolved']} | {values['conflict']} |")
    lines += [
        "", "## 3. Confirmed canonical records", "",
        "Confirmed canonical identities include the 10x breast/kidney/brain/embryo releases and slide/capture-area identifiers listed in the master table; Wu GSE176078 including CID45171/GSM5354535; cell2location E-MTAB-11114/11115; Ransick GSE129798/SRP192559/PRJNA532850; Lee GSE132465; the six real-dropout experiment identities; ST8059051 = Visium-29B = C05717-021/B1; and the CTA article plus raw/processed repository records.", "",
        "## 4. Unresolved records", "",
    ]
    for issue in unresolved:
        lines.append(f"- **{issue['priority']} {issue['issue_id']}** ({issue['linked_record_id']}): {issue['issue']} Required resolution: {issue['required_resolution']}")
    lines += ["", "## 5. Conflicts", ""]
    for conflict in conflicts:
        lines.append(f"- **{conflict['conflict_id']}** ({conflict['linked_record_id']}): {conflict['field']}. Canonical recommendation: {conflict['recommended_canonical_value']}")
    lines += [
        "", "## 6. Highest-priority manuscript risks", "",
        "- **CRITICAL:** The formal LaTeX and compiled BibTeX are absent, so Results/Methods/legends/key consistency cannot be certified.",
        "- **CRITICAL:** No unified formal dataset manifest exists; future-tense manifest language is not yet supportable.",
        "- **CRITICAL:** `ERS20065156` is mis-typed as a run in the current claim; ENA identifies it as a sample.",
        "- **MAJOR:** Historical hECA composite exports lack preserved cell-level merge/checksum manifests.",
        "- **MAJOR:** Vizgen per-sample download receipts/license terms and CTA exact Visium chemistry/slide metadata remain incomplete.",
        "", "## 7. Required future manuscript changes", "",
        "Changes are listed only in `manuscript_required_updates.tsv`. None were applied in this task.",
        "", "## 8. Data availability readiness", "",
        "**Not ready.** The inventory is sufficiently structured for drafting, but formal manuscript linkage, Vizgen terms, CTA exact platform metadata, composite-reference lineage, and processed-derivative deposition/location must be resolved first.",
        "", "## 9. Dataset manifest readiness", "",
        "**Not ready to freeze.** Review this audit, resolve CRITICAL/MAJOR records, supply the formal manuscript/BibTeX, then generate a separate formal manifest with approval.",
        "", "## 10. Guardrails", "",
        "```text",
        "Stage1 rerun: false", "Stage3A rerun: false", "Stage3B rerun: false", "Stage4 rerun: false", "CytoSPACE rerun: false",
        "Experimental outputs modified: false", "Formal manuscript modified: false", "Formal bibliography modified: false",
        "Formal manifest overwritten: false", "Figures modified: false", "```", "",
    ]
    return "\n".join(lines)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    records = build_records()
    crosswalk = build_crosswalk()
    unresolved = build_unresolved()
    conflicts = build_conflicts()
    coverage_summary = coverage(crosswalk)
    status_counts = Counter(row["status"] for row in records)
    summary = {
        "audit": "Task 4A dataset provenance, formal identifier, platform and accession audit",
        "generated_date": ACCESSED,
        "decision": "HOLD",
        "data_availability_ready": False,
        "dataset_manifest_ready_to_freeze": False,
        "formal_manuscript_found": False,
        "formal_bibliography_found": False,
        "formal_dataset_manifest_found": False,
        "counts": {
            "master_records": len(records), "experiment_crosswalk_rows": len(crosswalk),
            "unresolved_issues": len(unresolved), "conflicts": len(conflicts), "record_status": dict(status_counts),
        },
        "coverage": coverage_summary,
        "guardrails": {
            "stage1_rerun": False, "stage3a_rerun": False, "stage3b_rerun": False,
            "stage4_rerun": False, "cytospace_rerun": False, "experimental_outputs_modified": False,
            "formal_manuscript_modified": False, "formal_bibliography_modified": False,
            "formal_manifest_overwritten": False, "figures_modified": False,
        },
    }

    write_tsv(OUT / "dataset_provenance_master.tsv", records, MASTER_FIELDS)
    write_tsv(OUT / "experiment_input_crosswalk.tsv", crosswalk, CROSSWALK_FIELDS)
    write_tsv(OUT / "manuscript_dataset_occurrences.tsv", build_occurrences(), ["section", "subsection", "figure_panel", "quoted_context", "dataset_or_identifier", "identifier_type", "linked_record_id", "status", "required_action", "notes"])
    write_tsv(OUT / "bibliography_dataset_crosswalk.tsv", build_bibliography(records), ["bibtex_key", "publication_title", "authors", "year", "doi", "dataset_role", "linked_record_ids", "supports_dataset_claim", "status", "notes"])
    write_tsv(OUT / "manifest_consistency_audit.tsv", build_manifest_audit(records), ["record_id", "manuscript_value", "manifest_value", "config_value", "input_metadata_value", "official_source_value", "consistency_status", "recommended_canonical_value", "notes"])
    write_tsv(OUT / "accession_type_audit.tsv", build_accession_audit(), ["accession", "claimed_type_in_manuscript", "verified_type", "linked_record_id", "status", "recommended_wording"])
    write_tsv(OUT / "official_source_registry.tsv", build_source_registry(records), ["record_id", "source_authority", "source_type", "url", "accessed_date", "supports_fields", "status", "notes"])
    write_tsv(OUT / "unresolved_records.tsv", unresolved, ["issue_id", "linked_record_id", "field", "issue", "searched_evidence", "required_resolution", "priority"])
    write_tsv(OUT / "conflicts.tsv", conflicts, ["conflict_id", "linked_record_id", "field", "evidence_a", "evidence_b", "recommended_canonical_value", "status", "required_action"])
    write_tsv(OUT / "manuscript_required_updates.tsv", build_updates(), ["file", "section", "current_text", "problem", "canonical_record_id", "recommended_replacement", "priority"])
    write_tsv(OUT / "data_availability_inventory.tsv", build_availability(records), ["record_id", "formal_dataset_name", "source_database", "accession", "official_url", "download_or_access_method", "redistributable", "project_processed_data_planned", "required_statement", "unresolved_issue"])
    (OUT / "audit_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (OUT / "README.md").write_text(build_readme(summary), encoding="utf-8")
    (OUT / "final_decision.md").write_text(build_final(summary, unresolved, conflicts), encoding="utf-8")

    print("Task 4A dataset provenance audit completed.")
    print("\nDecision:\nHOLD")
    print("\nExperiment coverage:")
    for name, values in coverage_summary.items():
        print(f"{name}: {values['audited']} / {values['expected']}")
    print(f"\nConfirmed records:\n{status_counts['CONFIRMED']}")
    print(f"\nPartial records:\n{status_counts['PARTIAL']}")
    print(f"\nUnresolved records:\n{len(unresolved)} issues")
    print(f"\nConflicts:\n{len(conflicts)}")
    print("\nCritical manuscript risks:\nformal LaTeX absent; formal bibliography absent; formal dataset manifest absent; ERS20065156 accession type conflict")
    print("\nData availability ready:\nfalse")
    print("\nDataset manifest ready to freeze:\nfalse")
    print("\nFormal manuscript modified: false")
    print("Formal bibliography modified: false")
    print("Formal manifest overwritten: false")
    print("Experimental stages rerun: false")
    print(f"\nAudit script:\n{rel(Path(__file__))}")
    print(f"Audit output directory:\n{rel(OUT)}")


if __name__ == "__main__":
    main()
