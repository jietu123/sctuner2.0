from __future__ import annotations

import csv
import json
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PHASE0R_DIR = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase0r_candidate_reranking_after_janesick_mapping_stop"
OUT_DIR = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase1_pathology_endpoint_and_mapping_feasibility_gate"

CANDIDATE = "Computational pathology annotation / breast cancer Visium CTA"
ARTICLE_URL = "https://www.nature.com/articles/s41698-025-01104-3"
DATA_URL = "https://researchdata.se/en/catalogue/dataset/2025-97"
CODE_URL = "https://github.com/JohanHartmanGroupBioteam/BreastCancer_CTA"
PROCESSED_DOI = "https://doi.org/10.5281/zenodo.15211538"
WU_REF = "Wu et al. 2021 human breast cancer scRNA-seq reference"


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def inventory_rows() -> list[dict[str, object]]:
    return [
        {
            "file_id": "CTA_OUTPUT",
            "file_name_or_dataset_name": "CTA_output.txt / QuPath detection measurements",
            "modality": "Computational pathology annotation",
            "platform": "QuPath object classifier on H&E image",
            "sample_id": "BCSA patient tumor-region samples",
            "sample_description": "Cell-level tumor / immune / stroma annotations from H&E-associated breast cancer sections",
            "source_url_or_accession": CODE_URL,
            "file_type": "tabular annotation / QuPath measurement export",
            "estimated_size_mb": "small to moderate; exact file size requires Phase 2 inventory",
            "required_for_endpoint": True,
            "required_for_visium_mapping": False,
            "required_for_sc_reference": False,
            "required_for_phase2": True,
            "required_for_phase3": True,
            "download_priority": "Priority 1",
            "risk": "Must verify actual file availability and sample IDs before Phase 2",
            "notes": "GitHub vignette states CTA_output can be obtained from QuPath detection measurements and aligned to Visium spots.",
        },
        {
            "file_id": "HE_IMAGE_OR_CTA_IMAGE_INPUT",
            "file_name_or_dataset_name": "Paired H&E image / image-derived CTA input",
            "modality": "Pathology image",
            "platform": "H&E / microscopy image",
            "sample_id": "BCSA1/BCSA2/BCSA3/BCSA4 tumor regions",
            "sample_description": "Matched H&E images from Visium spatial transcriptomics analysis",
            "source_url_or_accession": ARTICLE_URL,
            "file_type": "image or image metadata",
            "estimated_size_mb": "potentially large if full-resolution image is needed",
            "required_for_endpoint": True,
            "required_for_visium_mapping": False,
            "required_for_sc_reference": False,
            "required_for_phase2": "Only if spot-level CTA output is not already available",
            "required_for_phase3": False,
            "download_priority": "Priority 3",
            "risk": "Full-resolution images may be large; avoid unless endpoint registration cannot be verified from processed output",
            "notes": "Endpoint is acceptable only because it is derived from H&E/pathology image, not SVTuner/CytoSPACE output.",
        },
        {
            "file_id": "CTA_SPOT_ALIGNMENT",
            "file_name_or_dataset_name": "CTA_align output / spot barcode cell-composition matrix",
            "modality": "Endpoint-to-spot mapping",
            "platform": "CTA_align R function",
            "sample_id": "Same Visium sample as CTA output",
            "sample_description": "Spot-level tumor / immune / stroma percentages or counts",
            "source_url_or_accession": CODE_URL,
            "file_type": "CSV/RData/table with spot barcode IDs",
            "estimated_size_mb": "small",
            "required_for_endpoint": True,
            "required_for_visium_mapping": True,
            "required_for_sc_reference": False,
            "required_for_phase2": True,
            "required_for_phase3": True,
            "download_priority": "Priority 1",
            "risk": "If only cell-level CTA is available, Phase 2 must reproduce alignment before mapping",
            "notes": "GitHub README states alignment result contains barcode ID and can quantify cell type percentages.",
        },
        {
            "file_id": "VISIUM_CTA_ST_OBJECT",
            "file_name_or_dataset_name": "Processed Visium Seurat object / ST_data.RData / expression matrix",
            "modality": "Spatial transcriptomics",
            "platform": "10x Visium",
            "sample_id": "BCSA tumor-region sample matching CTA output",
            "sample_description": "Visium breast cancer spatial expression object used for CTA/deconvolution comparison",
            "source_url_or_accession": f"{DATA_URL}; {PROCESSED_DOI}",
            "file_type": "RData / expression matrix / spot metadata",
            "estimated_size_mb": "moderate; processed preferred over raw FASTQ",
            "required_for_endpoint": False,
            "required_for_visium_mapping": True,
            "required_for_sc_reference": False,
            "required_for_phase2": True,
            "required_for_phase3": True,
            "download_priority": "Priority 2",
            "risk": "Researchdata raw data are restricted and about 300 GB; processed data should be preferred",
            "notes": "Paper states raw data are in SND/DORIS and processed data are in Zenodo.",
        },
        {
            "file_id": "VISIUM_COORDINATES_METADATA",
            "file_name_or_dataset_name": "Visium spot coordinates / tissue positions / image scaling metadata",
            "modality": "Spatial metadata",
            "platform": "10x Visium SpaceRanger / Seurat image slot",
            "sample_id": "Same as selected ST object",
            "sample_description": "Spot barcode, tissue position, image coordinate, scale factor metadata",
            "source_url_or_accession": f"{CODE_URL}; {PROCESSED_DOI}",
            "file_type": "CSV/RData metadata",
            "estimated_size_mb": "small",
            "required_for_endpoint": True,
            "required_for_visium_mapping": True,
            "required_for_sc_reference": False,
            "required_for_phase2": True,
            "required_for_phase3": True,
            "download_priority": "Priority 1",
            "risk": "Coordinate convention must be checked because the paper mentions y-axis reversal and scale factors",
            "notes": "Required to define endpoint-positive/negative spots before mapping.",
        },
        {
            "file_id": "SCRNA_REFERENCE_EXPRESSION",
            "file_name_or_dataset_name": WU_REF,
            "modality": "scRNA-seq reference",
            "platform": "Public human breast cancer scRNA-seq",
            "sample_id": "External breast tumor reference",
            "sample_description": "Reference used by the paper to deconvolute Visium data",
            "source_url_or_accession": "Wu et al. 2021; cited by CTA paper",
            "file_type": "expression matrix",
            "estimated_size_mb": "large but already familiar to project; use processed reference if possible",
            "required_for_endpoint": False,
            "required_for_visium_mapping": False,
            "required_for_sc_reference": True,
            "required_for_phase2": True,
            "required_for_phase3": True,
            "download_priority": "Priority 2",
            "risk": "Need collapse from 29 cell types to tumor/immune/stroma or SVTuner target labels consistently",
            "notes": "Paper used public scRNA reference for seven deconvolution methods including CytoSPACE.",
        },
        {
            "file_id": "SCRNA_CELL_METADATA",
            "file_name_or_dataset_name": "Wu et al. cell metadata / cell-type annotations",
            "modality": "scRNA-seq metadata",
            "platform": "Public human breast cancer scRNA-seq",
            "sample_id": "External breast tumor reference",
            "sample_description": "Cell type minor labels and collapsed labels",
            "source_url_or_accession": "Wu et al. 2021; cited by CTA paper",
            "file_type": "metadata table",
            "estimated_size_mb": "small to moderate",
            "required_for_endpoint": False,
            "required_for_visium_mapping": False,
            "required_for_sc_reference": True,
            "required_for_phase2": True,
            "required_for_phase3": True,
            "download_priority": "Priority 2",
            "risk": "Must preserve identical labels for baseline and SVTuner branches",
            "notes": "Required for CytoSPACE and SVTuner traceability.",
        },
        {
            "file_id": "DOCUMENTATION_README_METADATA",
            "file_name_or_dataset_name": "README, File_list.csv, SND_metadata.xlsx, md5sum.txt",
            "modality": "Documentation / sample metadata",
            "platform": "SND/Researchdata/Zenodo/GitHub",
            "sample_id": "All candidate samples",
            "sample_description": "Sample-to-file and patient-to-region metadata",
            "source_url_or_accession": DATA_URL,
            "file_type": "README/CSV/XLSX/checksum",
            "estimated_size_mb": "small",
            "required_for_endpoint": True,
            "required_for_visium_mapping": True,
            "required_for_sc_reference": True,
            "required_for_phase2": True,
            "required_for_phase3": False,
            "download_priority": "Priority 1",
            "risk": "Needed before downloading any large files",
            "notes": "Researchdata page exposes documentation files even though raw data are restricted.",
        },
    ]


def endpoint_independence_rows() -> list[dict[str, object]]:
    return [
        {
            "question": "What generated the computational pathology annotation?",
            "answer": "H&E image-derived computational tissue annotation using QuPath object classifier and pathology review",
            "status": "PASS",
            "evidence": "Article describes H&E stain correction, nucleus segmentation, object classifier using random trees, and pathologist review.",
            "source": ARTICLE_URL,
        },
        {
            "question": "Is it derived from H&E / pathology image / expert annotation / pathology model?",
            "answer": "Yes. It is derived from paired H&E images and pathology image analysis.",
            "status": "PASS",
            "evidence": "CTA was performed on matched H&E images from Visium spatial transcriptomics analysis.",
            "source": ARTICLE_URL,
        },
        {
            "question": "Independent from SVTuner output?",
            "answer": "Yes. SVTuner was not part of the source pipeline.",
            "status": "PASS",
            "evidence": "CTA predates our analysis and is generated from image morphology.",
            "source": ARTICLE_URL,
        },
        {
            "question": "Independent from CytoSPACE mapping output?",
            "answer": "Yes. CTA was used as an external reference to assess deconvolution methods including CytoSPACE, not derived from them.",
            "status": "PASS",
            "evidence": "Paper used CTA results as a reference to compare deconvolution outputs.",
            "source": ARTICLE_URL,
        },
        {
            "question": "Independent from Visium expression matrix?",
            "answer": "Yes for endpoint generation. CTA is generated from H&E morphology, then aligned to Visium spots.",
            "status": "PASS_WITH_CAVEAT",
            "evidence": "CTA is image-based; alignment to Visium uses coordinates. Phase 2 must ensure endpoint labels are not reconstructed from marker expression.",
            "source": ARTICLE_URL,
        },
        {
            "question": "Does it exist before SVTuner analysis?",
            "answer": "Yes. It is a published external endpoint and code/data resource.",
            "status": "PASS",
            "evidence": "Publication and GitHub pipeline exist before the current SVTuner analysis.",
            "source": f"{ARTICLE_URL}; {CODE_URL}",
        },
        {
            "question": "Can it define endpoint-positive / endpoint-negative spots?",
            "answer": "Yes in principle. CTA_align can generate spot-level cell composition via spot barcode IDs.",
            "status": "PASS",
            "evidence": "GitHub README states aligned result contains barcode ID and can calculate percentage of cell types.",
            "source": CODE_URL,
        },
        {
            "question": "Is it pre-defined rather than post hoc?",
            "answer": "Yes, if Phase 2 uses published CTA tumor/immune/stroma endpoint labels and freezes thresholds before mapping.",
            "status": "PASS_WITH_CAVEAT",
            "evidence": "Endpoint class family is predefined; exact positive/negative threshold must be frozen in Phase 2 before mapping.",
            "source": ARTICLE_URL,
        },
    ]


def spot_definability_rows() -> list[dict[str, object]]:
    return [
        {
            "sample_or_object": "CTA output",
            "annotation_resolution": "cell-level QuPath object classification",
            "spot_barcode_available_or_derivable": "derivable by CTA_align",
            "roi_polygon_or_image_coordinate_available": "image coordinates and cell centroid coordinates required",
            "additional_registration_needed": "No additional cross-sample registration if using paired H&E/Visium sample; coordinate convention must be checked",
            "endpoint_positive_spots_definable_before_mapping": True,
            "endpoint_negative_spots_definable_before_mapping": True,
            "endpoint_spatial_map_independently_plottable": True,
            "risk": "Phase 2 must verify real files include spot barcodes or enough image/coordinate metadata to reproduce CTA_align.",
        }
    ]


def sample_matching_rows() -> list[dict[str, object]]:
    return [
        {
            "component": "CTA endpoint vs Visium",
            "matching_status": "same or paired section / same tumor region according to paper workflow",
            "sample_ids": "BCSA tumor-region samples; exact sample must be selected in Phase 2",
            "mismatch_risk": "One sample was excluded in the paper due to image metadata mismatch; Phase 2 must avoid mismatched sample",
            "sufficient_for_endpoint_specific_evaluation": True,
            "notes": "Use only samples with verified CTA-to-spot alignment.",
        },
        {
            "component": "scRNA-seq reference vs Visium",
            "matching_status": "external public breast tumor scRNA reference, not necessarily same donor",
            "sample_ids": "Wu et al. 2021 reference",
            "mismatch_risk": "Disease/cohort reference mismatch possible but acceptable for CytoSPACE-style reference mapping if labels/gene overlap pass",
            "sufficient_for_endpoint_specific_evaluation": True,
            "notes": "Phase 2 must validate reference labels and gene intersection.",
        },
    ]


def formal_mapping_rows() -> list[dict[str, object]]:
    return [
        {
            "check": "identifiable ST mapping object",
            "status": "PASS",
            "evidence": "Processed Visium spatial data and ST_data.RData workflow are described.",
            "risk": "Need exact processed object for selected sample.",
        },
        {
            "check": "usable scRNA-seq reference",
            "status": "PASS",
            "evidence": "Paper used public Wu et al. breast tumor scRNA-seq reference for seven deconvolution methods.",
            "risk": "Need local processed expression and metadata import.",
        },
        {
            "check": "sufficient gene intersection likely",
            "status": "PASS_WITH_CAVEAT",
            "evidence": "Same paper performed deconvolution of Visium data with scRNA reference.",
            "risk": "Phase 2 must compute actual shared genes.",
        },
        {
            "check": "compatible preprocessing",
            "status": "PASS_WITH_CAVEAT",
            "evidence": "Paper provides method-specific folders and code, including CytoSPACE.",
            "risk": "Our SVTuner input contract may require conversion from R/Seurat to CSV.",
        },
        {
            "check": "baseline CytoSPACE feasibility",
            "status": "PASS_WITH_CAVEAT",
            "evidence": "Paper evaluated CytoSPACE among seven deconvolution methods.",
            "risk": "Need reproduce our own CytoSPACE input without relying on paper's final outputs.",
        },
        {
            "check": "SVTuner-enhanced mapping feasibility",
            "status": "PASS_WITH_CAVEAT",
            "evidence": "If CytoSPACE input tables can be built, SVTuner Stage3/Stage4 can share identical basis.",
            "risk": "Requires Phase 2 preflight and no use of post hoc endpoint labels inside mapping.",
        },
        {
            "check": "identical input basis for baseline and SVTuner",
            "status": "PASS_WITH_CAVEAT",
            "evidence": "Can be designed by freezing same ST matrix, scRNA reference, genes, and cell labels.",
            "risk": "Must be audited before Phase 3.",
        },
        {
            "check": "endpoint-specific comparison after mapping",
            "status": "PASS",
            "evidence": "CTA provides spot-level tumor/immune/stroma endpoint percentages for comparison.",
            "risk": "Thresholds for endpoint-positive/negative spots must be frozen before mapping.",
        },
        {
            "check": "expected compute feasibility on RTX 4060 / i5-12600 / 16GB RAM",
            "status": "PASS_WITH_CAVEAT",
            "evidence": "Minimal sample/downsampled mode should be possible if using processed matrices rather than raw FASTQ.",
            "risk": "Full raw dataset is about 300GB and not appropriate for Phase 1.",
        },
    ]


def minimal_download_rows() -> list[dict[str, object]]:
    return [
        {
            "priority": "Priority 1",
            "item": "README / File_list.csv / sample metadata / SND_metadata.xlsx / md5sum",
            "source": DATA_URL,
            "purpose": "Confirm sample IDs, file naming, access mode, and minimal processed targets.",
            "download_now": False,
            "forbidden": False,
            "notes": "Small documentation only; do not download raw FASTQ.",
        },
        {
            "priority": "Priority 1",
            "item": "CTA_output.txt or spot-level CTA annotation output for one selected sample",
            "source": CODE_URL,
            "purpose": "Define endpoint-positive/negative spots before mapping.",
            "download_now": False,
            "forbidden": False,
            "notes": "Prefer already aligned spot-level output if available.",
        },
        {
            "priority": "Priority 1",
            "item": "Visium spot coordinates / tissue positions / image metadata for same selected sample",
            "source": f"{CODE_URL}; {PROCESSED_DOI}",
            "purpose": "Verify endpoint-to-spot definability.",
            "download_now": False,
            "forbidden": False,
            "notes": "Needed before expression matrix download.",
        },
        {
            "priority": "Priority 2",
            "item": "Processed Visium expression object for selected sample",
            "source": PROCESSED_DOI,
            "purpose": "Build baseline/SVTuner ST mapping object.",
            "download_now": False,
            "forbidden": False,
            "notes": "Download only selected sample, not all samples blindly.",
        },
        {
            "priority": "Priority 2",
            "item": "Processed Wu et al. scRNA expression and metadata",
            "source": "Wu et al. 2021 source used by CTA paper",
            "purpose": "Build scRNA reference for CytoSPACE/SVTuner.",
            "download_now": False,
            "forbidden": False,
            "notes": "Use existing local reference if available and not deleted; otherwise download processed reference only.",
        },
        {
            "priority": "Priority 3",
            "item": "H&E/pathology full-resolution images",
            "source": DATA_URL,
            "purpose": "Only needed if CTA-to-spot output cannot be obtained directly.",
            "download_now": False,
            "forbidden": False,
            "notes": "Potentially large; avoid unless registration validation requires it.",
        },
        {
            "priority": "Forbidden",
            "item": "Raw FASTQ / all raw sequencing files / all replicates blindly",
            "source": DATA_URL,
            "purpose": "Not needed for Phase 1 feasibility gate.",
            "download_now": False,
            "forbidden": True,
            "notes": "Researchdata page reports raw dataset total size about 300GB.",
        },
    ]


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    phase0r_summary = read_json(PHASE0R_DIR / "bioapp_phase0r_summary.json")
    recommended = read_csv(PHASE0R_DIR / "bioapp_phase0r_recommended_candidates.csv")

    inventory_fields = [
        "file_id",
        "file_name_or_dataset_name",
        "modality",
        "platform",
        "sample_id",
        "sample_description",
        "source_url_or_accession",
        "file_type",
        "estimated_size_mb",
        "required_for_endpoint",
        "required_for_visium_mapping",
        "required_for_sc_reference",
        "required_for_phase2",
        "required_for_phase3",
        "download_priority",
        "risk",
        "notes",
    ]
    write_csv(OUT_DIR / "bioapp_phase1_pathology_data_inventory.csv", inventory_rows(), inventory_fields)

    independence_fields = ["question", "answer", "status", "evidence", "source"]
    independence = endpoint_independence_rows()
    write_csv(OUT_DIR / "bioapp_phase1_endpoint_independence_audit.csv", independence, independence_fields)
    write_text(
        OUT_DIR / "bioapp_phase1_endpoint_independence_audit_report.txt",
        "\n".join(
            [
                "Endpoint independence audit",
                "",
                "Conclusion: PASS_WITH_CAVEAT.",
                "CTA is acceptable as an independent external endpoint if Phase 2 uses published H&E/QuPath-derived CTA labels or CTA-to-spot outputs, not marker-expression-derived labels.",
                "It is independent from SVTuner and CytoSPACE because those methods do not generate the CTA endpoint.",
                "It is independent from Visium expression for endpoint generation because CTA is image-derived, although it is aligned onto Visium coordinates.",
                "Phase 2 must freeze endpoint-positive and endpoint-negative definitions before any mapping.",
                "",
                f"Primary sources: {ARTICLE_URL}; {CODE_URL}; {DATA_URL}",
            ]
        )
        + "\n",
    )

    spot_rows = spot_definability_rows()
    spot_fields = [
        "sample_or_object",
        "annotation_resolution",
        "spot_barcode_available_or_derivable",
        "roi_polygon_or_image_coordinate_available",
        "additional_registration_needed",
        "endpoint_positive_spots_definable_before_mapping",
        "endpoint_negative_spots_definable_before_mapping",
        "endpoint_spatial_map_independently_plottable",
        "risk",
    ]
    write_csv(OUT_DIR / "bioapp_phase1_endpoint_spot_definability_table.csv", spot_rows, spot_fields)
    write_text(
        OUT_DIR / "bioapp_phase1_endpoint_spatial_registration_report.txt",
        "\n".join(
            [
                "Endpoint spatial registration gate",
                "",
                "Conclusion: PASS_WITH_CAVEAT.",
                "The endpoint can be mapped to Visium spots because the CTA workflow aligns QuPath cell/object annotations to Visium spot barcodes.",
                "The GitHub vignette documents CTA_align with ST_data, image name, scaling factor, pixel size, raw image dimensions, and CTA output.",
                "The publication states cells whose centroid coordinates fall within Visium spot areas are assigned spot IDs and cell type percentages are computed.",
                "Phase 2 must select only samples with verified image/metadata matching; the paper notes one sample was excluded due to image metadata mismatch.",
            ]
        )
        + "\n",
    )

    sample_fields = [
        "component",
        "matching_status",
        "sample_ids",
        "mismatch_risk",
        "sufficient_for_endpoint_specific_evaluation",
        "notes",
    ]
    write_csv(OUT_DIR / "bioapp_phase1_sample_matching_table.csv", sample_matching_rows(), sample_fields)
    write_text(
        OUT_DIR / "bioapp_phase1_sample_matching_report.txt",
        "\n".join(
            [
                "Sample matching gate",
                "",
                "Conclusion: PASS_WITH_CAVEAT.",
                "CTA and Visium are paired by the study design for breast cancer tumor-region samples.",
                "The scRNA reference is public breast tumor scRNA-seq, not necessarily same donor; this is acceptable for mapping feasibility but must be labeled as external reference.",
                "Phase 2 must identify one selected sample where CTA output, Visium expression, spot coordinates, and image metadata agree.",
            ]
        )
        + "\n",
    )

    mapping_fields = ["check", "status", "evidence", "risk"]
    mapping = formal_mapping_rows()
    write_csv(OUT_DIR / "bioapp_phase1_formal_mapping_feasibility_table.csv", mapping, mapping_fields)
    write_text(
        OUT_DIR / "bioapp_phase1_formal_mapping_feasibility_report.txt",
        "\n".join(
            [
                "Formal mapping feasibility gate",
                "",
                "Conclusion: PASS_WITH_CAVEAT.",
                "The candidate satisfies Golden Rule 18 at Phase 1 evidence level: a Visium mapping object, scRNA reference, and endpoint-specific comparison route are identifiable.",
                "The paper evaluated several deconvolution methods including CytoSPACE using the public scRNA reference and CTA endpoint, so baseline mapping is plausible.",
                "This is not a formal mapping result. Phase 2 must perform file-level import, shared-gene preflight, and identical-input-basis design before any CytoSPACE/SVTuner run.",
            ]
        )
        + "\n",
    )

    download_fields = ["priority", "item", "source", "purpose", "download_now", "forbidden", "notes"]
    write_csv(OUT_DIR / "bioapp_phase1_minimal_download_plan.csv", minimal_download_rows(), download_fields)

    checks = {
        "stage_type": "computational pathology endpoint and mapping feasibility gate",
        "biological_question_defined": True,
        "external_endpoint_predefined": True,
        "endpoint_independent_from_SVTuner": True,
        "endpoint_independent_from_CytoSPACE": True,
        "endpoint_independent_from_Visium_expression": True,
        "endpoint_spatial_registration_feasible": True,
        "endpoint_positive_negative_spots_definable_before_mapping": True,
        "sample_matching_feasible": True,
        "visium_CTA_data_identifiable": True,
        "scRNAseq_reference_identifiable": True,
        "formal_baseline_CytoSPACE_feasible": True,
        "SVTuner_mapping_feasible": True,
        "identical_input_basis_for_baseline_and_SVTuner_feasible": True,
        "endpoint_specific_quantitative_metric_feasible": True,
        "endpoint_baseline_SVTuner_spatial_comparison_feasible": True,
        "compute_feasible_on_16GB_RAM_or_downsample_mode": True,
        "biological_application_allowed_to_phase2": True,
        "decision": "PASS",
        "mandatory_phase2_caveats": [
            "Use processed files only; do not download raw FASTQ.",
            "Freeze endpoint-positive/negative spot rule before mapping.",
            "Verify CTA is from H&E/QuPath output and not reconstructed from expression.",
            "Verify exact selected sample has matching CTA, Visium expression, coordinates, and metadata.",
            "Run shared-gene/input-basis preflight before CytoSPACE/SVTuner.",
        ],
    }
    write_json(OUT_DIR / "bioapp_phase1_golden_rules_v2_1_check.json", checks)

    risk_report = "\n".join(
        [
            "BioApp Phase 1 risk report",
            "",
            "Main residual risks:",
            "1. Raw SND/DORIS data are restricted and approximately 300GB; Phase 2 must avoid raw FASTQ.",
            "2. Processed data availability and exact selected-sample files must be verified before mapping.",
            "3. One sample in the paper was excluded because of image metadata mismatch; Phase 2 must avoid mismatched samples.",
            "4. Endpoint thresholds must be frozen before any CytoSPACE/SVTuner mapping.",
            "5. CytoSPACE feasibility is literature-supported, not yet locally executed for this candidate.",
            "",
            "Decision impact: these are Phase 2 inventory/preflight risks, not Phase 1 blockers.",
        ]
    )
    write_text(OUT_DIR / "bioapp_phase1_risk_report.txt", risk_report + "\n")

    decision = "PASS"
    summary = {
        "phase": "BioApp Phase 1 - computational pathology endpoint and mapping feasibility gate",
        "decision": decision,
        "candidate": CANDIDATE,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "phase0r_input_available": bool(phase0r_summary),
        "phase0r_recommended_candidate_available": bool(recommended),
        "endpoint_independent_from_SVTuner": True,
        "endpoint_independent_from_CytoSPACE": True,
        "endpoint_independent_from_Visium_expression": True,
        "endpoint_positive_negative_spots_definable_before_mapping": True,
        "sample_matching_feasible": True,
        "formal_baseline_CytoSPACE_feasible": True,
        "SVTuner_mapping_feasible": True,
        "identical_input_basis_feasible": True,
        "minimal_download_plan_generated": True,
        "large_raw_data_downloaded": False,
        "raw_FASTQ_downloaded": False,
        "cytospace_run": False,
        "svtuner_run": False,
        "stage4_run": False,
        "formal_metrics_recomputed": False,
        "biological_application_allowed_to_phase2": True,
        "output_directory": str(OUT_DIR),
        "next": "Manual review before BioApp Phase 2 endpoint spatial registration",
    }
    write_json(OUT_DIR / "bioapp_phase1_summary.json", summary)
    write_text(
        OUT_DIR / "bioapp_phase1_readme.txt",
        "BioApp Phase 1 checked the computational pathology endpoint and mapping feasibility gate for the breast cancer Visium CTA candidate. "
        "No data download, CytoSPACE, SVTuner, Stage4, mapping, or formal metric recomputation was performed.\n",
    )
    write_text(OUT_DIR / "bioapp_phase1_decision.txt", decision + "\n")

    print("BioApp Phase 1 completed.")
    print("\nDecision:")
    print(decision)
    print("\nCandidate:")
    print(CANDIDATE)
    print("\nEndpoint independent from SVTuner:")
    print("true")
    print("\nEndpoint independent from CytoSPACE:")
    print("true")
    print("\nEndpoint independent from Visium expression:")
    print("true")
    print("\nEndpoint positive/negative spots definable before mapping:")
    print("true")
    print("\nSample matching feasible:")
    print("true")
    print("\nFormal baseline CytoSPACE feasible:")
    print("true")
    print("\nSVTuner mapping feasible:")
    print("true")
    print("\nIdentical input basis feasible:")
    print("true")
    print("\nMinimal download plan generated:")
    print("true")
    print("\nLarge raw data downloaded:")
    print("false")
    print("\nRaw FASTQ downloaded:")
    print("false")
    print("\nCytoSPACE run:")
    print("false")
    print("\nSVTuner run:")
    print("false")
    print("\nStage4 run:")
    print("false")
    print("\nFormal metrics recomputed:")
    print("false")
    print("\nOutput directory:")
    print(OUT_DIR)
    print("\nNext:")
    print("Manual review before BioApp Phase 2 endpoint spatial registration")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
