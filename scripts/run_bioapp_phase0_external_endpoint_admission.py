from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = PROJECT_ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase0_external_endpoint_admission"


@dataclass
class Candidate:
    candidate_id: str
    dataset_name: str
    publication_or_source: str
    organ_or_disease: str
    biological_question_candidate: str
    external_endpoint_type: str
    external_endpoint_description: str
    endpoint_predefined_before_analysis: str
    endpoint_independent_from_SVTuner: str
    endpoint_spatially_registered_or_mappable: str
    endpoint_coordinate_or_roi_available: str
    visium_or_ST_data_available: str
    scRNAseq_reference_available: str
    same_sample_or_matched_sample_status: str
    baseline_mapping_possible: str
    endpoint_specific_metric_possible: str
    spatial_comparison_possible: str
    estimated_data_size: str
    download_required_for_phase1: str
    main_risk: str
    golden_rules_preliminary_status: str
    recommendation: str
    source_link_or_accession: str
    notes: str


def classify(row: Candidate) -> str:
    positives = {
        "endpoint": row.endpoint_independent_from_SVTuner.lower().startswith("yes"),
        "registered": row.endpoint_spatially_registered_or_mappable.lower().startswith("yes"),
        "coords": row.endpoint_coordinate_or_roi_available.lower().startswith("yes"),
        "st": row.visium_or_ST_data_available.lower().startswith("yes"),
        "ref": row.scRNAseq_reference_available.lower().startswith("yes"),
        "metric": row.endpoint_specific_metric_possible.lower().startswith("yes"),
        "spatial": row.spatial_comparison_possible.lower().startswith("yes"),
    }
    if all(positives.values()):
        return "Tier 1"
    if positives["endpoint"] and (positives["st"] or positives["ref"] or positives["registered"]):
        return "Tier 2"
    return "Tier 3"


def golden_rules(row: Candidate) -> dict[str, bool | str]:
    tier = classify(row)
    allowed = tier == "Tier 1"
    return {
        "candidate_id": row.candidate_id,
        "tier": tier,
        "biological_question_defined": row.biological_question_candidate.strip().lower() not in {"", "none", "unknown"},
        "external_endpoint_predefined": row.endpoint_predefined_before_analysis.lower().startswith("yes"),
        "endpoint_independent_from_SVTuner": row.endpoint_independent_from_SVTuner.lower().startswith("yes"),
        "endpoint_spatially_registered": row.endpoint_spatially_registered_or_mappable.lower().startswith("yes"),
        "baseline_comparison_available": row.baseline_mapping_possible.lower().startswith("yes"),
        "endpoint_specific_improvement_defined": row.endpoint_specific_metric_possible.lower().startswith("yes"),
        "endpoint_specific_quantitative_metric_available": row.endpoint_specific_metric_possible.lower().startswith("yes"),
        "endpoint_baseline_SVTuner_spatial_comparison_available": row.spatial_comparison_possible.lower().startswith("yes"),
        "interpretation_boundary_defined": "no manuscript" in row.notes.lower()
        or "endpoint admission only" in row.notes.lower()
        or "not a success claim" in row.notes.lower(),
        "biological_application_allowed_preliminary": allowed,
    }


def candidates() -> list[Candidate]:
    return [
        Candidate(
            candidate_id="T1_XENIUM_JANESICK_BREAST_2023",
            dataset_name="Janesick Xenium-Visium breast cancer",
            publication_or_source=(
                "Janesick et al., Nature Communications 2023; 10x Genomics companion "
                "data and GitHub registration notebooks; GEO GSE243280"
            ),
            organ_or_disease="Human breast cancer FFPE",
            biological_question_candidate=(
                "Does SVTuner better recover or avoid forced contradiction in a "
                "Xenium-defined humoral/B-lineage spatial region in serial breast cancer sections?"
            ),
            external_endpoint_type="Xenium cell annotation / ROI endpoint",
            external_endpoint_description=(
                "Xenium-derived single-cell spatial annotations can define B-cell, plasma-cell, "
                "or humoral immune endpoint-positive regions independently of SVTuner."
            ),
            endpoint_predefined_before_analysis="Yes - endpoint family is specified before SVTuner analysis",
            endpoint_independent_from_SVTuner="Yes - Xenium annotations are generated outside SVTuner",
            endpoint_spatially_registered_or_mappable=(
                "Yes - public companion code describes Xenium/Visium serial-section keypoint registration"
            ),
            endpoint_coordinate_or_roi_available="Yes - Xenium cell coordinates and registration workflow are public",
            visium_or_ST_data_available="Yes - Visium CytAssist whole-transcriptome data are part of the study",
            scRNAseq_reference_available="Yes - scFFPE-seq reference data are part of the study",
            same_sample_or_matched_sample_status="Matched serial FFPE sections from one breast cancer tissue block",
            baseline_mapping_possible=(
                "Yes in principle - CytoSPACE baseline can be attempted after Phase 1 inventory"
            ),
            endpoint_specific_metric_possible=(
                "Yes - ROI enrichment, AUROC/AUPRC, effect size, and spatial overlap/correlation are definable"
            ),
            spatial_comparison_possible=(
                "Yes - endpoint ROI can be compared with baseline and SVTuner spatial outputs after mapping"
            ),
            estimated_data_size="Large; exact download size not assessed in Phase 0",
            download_required_for_phase1="Yes - local raw endpoint package was not found in this workspace",
            main_risk=(
                "The old local Janesick formal-mapping line was stopped after CytoSPACE timeout/fallback collapse; "
                "Phase 1 must re-inventory files and feasibility from scratch."
            ),
            golden_rules_preliminary_status="PASS_PRELIMINARY_TIER1",
            recommendation="Tier 1 - recommended for manual review before Phase 1",
            source_link_or_accession=(
                "https://www.nature.com/articles/s41467-023-43458-x ; "
                "https://github.com/10XGenomics/janesick_nature_comms_2023_companion ; "
                "GSE243280"
            ),
            notes=(
                "Endpoint admission only; not a success claim. Do not reuse deleted Janesick results. "
                "No manuscript text, no mapping, no formal metrics in Phase 0."
            ),
        ),
        Candidate(
            candidate_id="T2_XENIUM_COSMX_SPATIAL_TOUCHSTONE",
            dataset_name="Spatial Touchstone Xenium/CosMx multi-tissue benchmark",
            publication_or_source="Spatial Touchstone / Nature Biotechnology 2025 / GSE277080",
            organ_or_disease="Six tissue types; benchmark samples",
            biological_question_candidate=(
                "Can an externally annotated imaging endpoint be used to test endpoint-concordant "
                "SVTuner behavior across standardized tissues?"
            ),
            external_endpoint_type="Xenium/CosMx imaging-based cell annotation endpoint",
            external_endpoint_description="Large public imaging-based spatial repository across Xenium and CosMx platforms",
            endpoint_predefined_before_analysis="Yes - platform-derived imaging annotations are external",
            endpoint_independent_from_SVTuner="Yes - generated outside SVTuner",
            endpoint_spatially_registered_or_mappable="Unclear - Visium/ST mapping object is not established in Phase 0",
            endpoint_coordinate_or_roi_available="Yes - imaging-based spatial coordinates are expected",
            visium_or_ST_data_available="Unclear - primary resource is imaging-based benchmark, not a confirmed Visium mapping pair",
            scRNAseq_reference_available="Unclear - annotations/reference relationship needs manual confirmation",
            same_sample_or_matched_sample_status="Unclear",
            baseline_mapping_possible="No - not established without a paired ST mapping object",
            endpoint_specific_metric_possible="Yes - possible if a paired ST object is identified",
            spatial_comparison_possible="Unclear - requires paired ST or aggregation strategy",
            estimated_data_size="Large",
            download_required_for_phase1="Yes, if selected",
            main_risk="Strong endpoint resource but not yet a proven SVTuner/CytoSPACE mapping scenario",
            golden_rules_preliminary_status="REVIEW_REQUIRED_TIER2",
            recommendation="Tier 2 - manual confirmation required",
            source_link_or_accession=(
                "https://www.nature.com/articles/s41587-025-02811-9 ; GSE277080 ; "
                "https://www.spatialtouchstone.org"
            ),
            notes="Endpoint admission only; no manuscript text and no formal metrics in Phase 0.",
        ),
        Candidate(
            candidate_id="T2_MULTI_PLATFORM_TUMOR_BENCHMARK",
            dataset_name="Matched FFPE tumor multi-platform ST benchmark",
            publication_or_source="Systematic benchmarking / technical comparison of Visium, Xenium, CosMx and related platforms",
            organ_or_disease="Human tumor FFPE samples",
            biological_question_candidate=(
                "Can a platform-defined imaging endpoint be compared against Visium/Visium HD mapping in matched tumors?"
            ),
            external_endpoint_type="Xenium/CosMx imaging endpoint",
            external_endpoint_description="Matched-platform benchmark includes imaging-based spatial modalities and Visium-family data",
            endpoint_predefined_before_analysis="Yes in principle",
            endpoint_independent_from_SVTuner="Yes in principle",
            endpoint_spatially_registered_or_mappable="Unclear - sample-level matching does not prove coordinate registration",
            endpoint_coordinate_or_roi_available="Unclear",
            visium_or_ST_data_available="Yes in principle - Visium-family platforms are described",
            scRNAseq_reference_available="Yes in principle - matched scRNA-seq is described in benchmark reporting",
            same_sample_or_matched_sample_status="Matched FFPE tumors, exact per-sample pairing requires review",
            baseline_mapping_possible="Unclear - needs data access and compatible spot-level matrix",
            endpoint_specific_metric_possible="Yes in principle",
            spatial_comparison_possible="Unclear - coordinate/ROI pairing must be verified",
            estimated_data_size="Large",
            download_required_for_phase1="Yes, if selected",
            main_risk="May be a methods benchmark rather than a clean biological endpoint application",
            golden_rules_preliminary_status="REVIEW_REQUIRED_TIER2",
            recommendation="Tier 2 - reserve backup candidate",
            source_link_or_accession=(
                "https://pmc.ncbi.nlm.nih.gov/articles/PMC12888464/ ; "
                "https://www.nature.com/articles/s41467-025-64292-3"
            ),
            notes="Endpoint admission only; not a success claim and no manuscript text in Phase 0.",
        ),
        Candidate(
            candidate_id="T3_COSMX_NSCLC_PUBLIC",
            dataset_name="CosMx SMI FFPE human NSCLC public dataset",
            publication_or_source="Bruker/NanoString public CosMx SMI datasets",
            organ_or_disease="Human NSCLC FFPE",
            biological_question_candidate="Potential immune niche endpoint in lung cancer",
            external_endpoint_type="CosMx cell annotation endpoint",
            external_endpoint_description="Public CosMx single-cell spatial dataset with cell-level imaging coordinates",
            endpoint_predefined_before_analysis="Yes in principle",
            endpoint_independent_from_SVTuner="Yes",
            endpoint_spatially_registered_or_mappable="No - no paired Visium/ST mapping object verified",
            endpoint_coordinate_or_roi_available="Yes",
            visium_or_ST_data_available="No paired Visium/ST object verified",
            scRNAseq_reference_available="No matched scRNA-seq reference verified",
            same_sample_or_matched_sample_status="Not established",
            baseline_mapping_possible="No",
            endpoint_specific_metric_possible="No - endpoint cannot be compared to mapping without ST object",
            spatial_comparison_possible="No",
            estimated_data_size="Large",
            download_required_for_phase1="Not recommended",
            main_risk="Independent endpoint exists, but not an admissible SVTuner biological application scenario",
            golden_rules_preliminary_status="REJECT_TIER3",
            recommendation="Tier 3 - reject for Phase 1",
            source_link_or_accession="https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/",
            notes="Endpoint admission only; no manuscript text and no mapping in Phase 0.",
        ),
        Candidate(
            candidate_id="T3_LOCAL_BREAST_MERSCOPE",
            dataset_name="Local breast MERSCOPE/high-resolution branch",
            publication_or_source="Local configs/data under cytospace_fig2k_breast_merscope and highres branches",
            organ_or_disease="Human breast cancer",
            biological_question_candidate="Possible high-resolution immune endpoint, not Xenium/CosMx/MERFISH+Visium admitted",
            external_endpoint_type="MERSCOPE/Vizgen cell metadata",
            external_endpoint_description="Local high-resolution cell metadata exists, but not the requested admitted endpoint-plus-ST scenario",
            endpoint_predefined_before_analysis="No - not specified for the fourth experiment",
            endpoint_independent_from_SVTuner="Yes in principle",
            endpoint_spatially_registered_or_mappable="No - no paired Visium/ST mapping endpoint registration established",
            endpoint_coordinate_or_roi_available="Yes in local high-resolution branch",
            visium_or_ST_data_available="No paired mapping object established for this endpoint",
            scRNAseq_reference_available="Yes in related breast scRNA resources, not admitted here",
            same_sample_or_matched_sample_status="Unclear",
            baseline_mapping_possible="No for this endpoint design",
            endpoint_specific_metric_possible="No",
            spatial_comparison_possible="No",
            estimated_data_size="Already local, but not admissible",
            download_required_for_phase1="No",
            main_risk="Would reuse old/different branch instead of the new endpoint admission",
            golden_rules_preliminary_status="REJECT_TIER3",
            recommendation="Tier 3 - reject for this experiment",
            source_link_or_accession="Local workspace only",
            notes="Endpoint admission only; do not reuse old results and no manuscript text in Phase 0.",
        ),
        Candidate(
            candidate_id="T3_LOCAL_BRCA_PATHOLOGY_VISIUM",
            dataset_name="Local BRCA Visium/plasma/pathology legacy branch",
            publication_or_source="Local BRCA Visium/scRNA configs and previous plasma-cell outputs",
            organ_or_disease="Human breast cancer",
            biological_question_candidate="Potential pathology or histology ROI endpoint, but not independently admitted",
            external_endpoint_type="Pathology/histology ROI candidate",
            external_endpoint_description="Visium and scRNA resources exist, but no predefined independent ROI endpoint was verified",
            endpoint_predefined_before_analysis="No",
            endpoint_independent_from_SVTuner="No - current endpoint would be reconstructed post hoc",
            endpoint_spatially_registered_or_mappable="No independent endpoint registration verified",
            endpoint_coordinate_or_roi_available="No admitted ROI file verified",
            visium_or_ST_data_available="Yes - local Visium breast branches exist",
            scRNAseq_reference_available="Yes - local BRCA scRNA branches exist",
            same_sample_or_matched_sample_status="Variable legacy branches",
            baseline_mapping_possible="Yes technically, but not biologically admissible",
            endpoint_specific_metric_possible="No admitted endpoint metric",
            spatial_comparison_possible="No admitted endpoint comparison",
            estimated_data_size="Already local",
            download_required_for_phase1="No",
            main_risk="Real data plus marker pattern would be mistaken for biological application",
            golden_rules_preliminary_status="REJECT_TIER3",
            recommendation="Tier 3 - reject unless an independent ROI endpoint is provided",
            source_link_or_accession="Local workspace only",
            notes="Endpoint admission only; do not reuse BCSA/starfysh/deleted or legacy results.",
        ),
    ]


def write_csv(path: Path, rows: list[Candidate]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()) + ["tier"])
        writer.writeheader()
        for row in rows:
            data = asdict(row)
            data["tier"] = classify(row)
            writer.writerow(data)


def write_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = candidates()
    tier1 = [row for row in rows if classify(row) == "Tier 1"]
    tier2 = [row for row in rows if classify(row) == "Tier 2"]
    rejected = [row for row in rows if classify(row) == "Tier 3"]
    recommended = tier1[0] if tier1 else None
    decision = "PASS" if tier1 else ("REVIEW_REQUIRED" if tier2 else "FAIL")

    write_csv(OUT_DIR / "bioapp_phase0_candidate_dataset_table.csv", rows)
    write_csv(OUT_DIR / "bioapp_phase0_tier1_candidates.csv", tier1)
    write_csv(OUT_DIR / "bioapp_phase0_tier2_candidates.csv", tier2)
    write_csv(OUT_DIR / "bioapp_phase0_rejected_candidates.csv", rejected)

    checks = {
        row.candidate_id: golden_rules(row)
        for row in rows
    }
    write_json(OUT_DIR / "bioapp_phase0_golden_rules_check.json", checks)

    summary = {
        "phase": "BioApp Phase 0 鈥?external endpoint admission",
        "decision": decision,
        "candidate_datasets_screened": len(rows),
        "tier1_candidates": len(tier1),
        "tier2_candidates": len(tier2),
        "rejected_candidates": len(rejected),
        "recommended_candidate": recommended.dataset_name if recommended else None,
        "external_endpoint": recommended.external_endpoint_description if recommended else None,
        "st_mapping_object": recommended.visium_or_ST_data_available if recommended else None,
        "scRNAseq_reference": recommended.scRNAseq_reference_available if recommended else None,
        "biological_application_allowed_preliminary": bool(recommended),
        "large_raw_data_downloaded": False,
        "cytospace_run": False,
        "svtuner_run": False,
        "stage4_run": False,
        "formal_metrics_recomputed": False,
        "output_directory": str(OUT_DIR).replace("\\", "/"),
        "next": "Manual review before BioApp Phase 1 data inventory and feasibility",
    }
    write_json(OUT_DIR / "bioapp_phase0_summary.json", summary)

    report_lines = [
        "BioApp Phase 0 鈥?external endpoint admission",
        "",
        "Stage objective:",
        "Screen public and local candidate data scenarios for an independent external spatial endpoint paired with a Visium/ST mapping object and scRNA-seq reference.",
        "",
        f"Decision: {decision}",
        f"Candidate datasets screened: {len(rows)}",
        f"Tier 1 candidates: {len(tier1)}",
        f"Tier 2 candidates: {len(tier2)}",
        f"Rejected candidates: {len(rejected)}",
        "",
        "Does at least one Tier 1 candidate exist?",
        "Yes" if tier1 else "No",
        "",
    ]
    if recommended:
        report_lines += [
            f"Most recommended candidate: {recommended.dataset_name}",
            f"External endpoint: {recommended.external_endpoint_description}",
            f"Biological question: {recommended.biological_question_candidate}",
            f"ST mapping object: {recommended.visium_or_ST_data_available}",
            f"scRNA-seq reference: {recommended.scRNAseq_reference_available}",
            "",
            "Why it preliminarily satisfies Golden Rules v2.0:",
            "- The biological question is endpoint-specific and defined before analysis.",
            "- The endpoint is a Xenium-derived cell/ROI annotation independent from SVTuner.",
            "- Public companion material describes Xenium/Visium serial-section registration.",
            "- Visium CytAssist data and scFFPE-seq reference data are part of the same study.",
            "- Endpoint-specific metrics and spatial comparison are definable after Phase 1 inventory and later mapping.",
            "- Interpretation is bounded to endpoint concordance / prevention of forced contradictory calls, not discovery of new biology from real data alone.",
            "",
            f"Main risk: {recommended.main_risk}",
            "",
            "Is Phase 1 worth starting?",
            "Yes, after manual review, because Phase 0 found one Tier 1 public candidate. Phase 1 must still verify files, local paths, registration, and run feasibility before any mapping.",
        ]
    else:
        report_lines += [
            "Most recommended candidate: null",
            "Phase 1 is not justified without manual intervention.",
        ]
    report_lines += [
        "",
        "Hard boundaries honored:",
        "- No large raw datasets downloaded.",
        "- CytoSPACE was not run.",
        "- SVTuner was not run.",
        "- Stage4 was not run.",
        "- No formal metrics were recomputed.",
        "- No dropout experiment was created.",
        "- No manuscript text was written.",
    ]
    (OUT_DIR / "bioapp_phase0_endpoint_admission_report.txt").write_text(
        "\n".join(report_lines) + "\n",
        encoding="utf-8",
    )

    readme = [
        "BioApp Phase 0 output directory",
        "",
        "This directory contains endpoint-admission outputs only.",
        "The PASS decision means at least one candidate is eligible for manual review before Phase 1.",
        "It does not mean that SVTuner succeeded, that mapping is feasible, or that the biological application is complete.",
        "",
        "Files:",
        "- bioapp_phase0_summary.json",
        "- bioapp_phase0_readme.txt",
        "- bioapp_phase0_candidate_dataset_table.csv",
        "- bioapp_phase0_endpoint_admission_report.txt",
        "- bioapp_phase0_golden_rules_check.json",
        "- bioapp_phase0_tier1_candidates.csv",
        "- bioapp_phase0_tier2_candidates.csv",
        "- bioapp_phase0_rejected_candidates.csv",
        "- bioapp_phase0_decision.txt",
    ]
    (OUT_DIR / "bioapp_phase0_readme.txt").write_text("\n".join(readme) + "\n", encoding="utf-8")

    decision_text = [
        decision,
        "",
        "Reason:",
        "At least one Tier 1 candidate was found." if tier1 else "No Tier 1 candidate was found.",
        "",
        "Manual review is required before BioApp Phase 1 data inventory and feasibility.",
        "No mapping may be started directly from Phase 0.",
    ]
    (OUT_DIR / "bioapp_phase0_decision.txt").write_text("\n".join(decision_text) + "\n", encoding="utf-8")

    print("BioApp Phase 0 completed.")
    print()
    print("Decision:")
    print(decision)
    print()
    print("Candidate datasets screened:")
    print(len(rows))
    print()
    print("Tier 1 candidates:")
    print(len(tier1))
    print()
    print("Tier 2 candidates:")
    print(len(tier2))
    print()
    print("Rejected candidates:")
    print(len(rejected))
    print()
    print("Recommended candidate:")
    print(recommended.dataset_name if recommended else "null")
    print()
    print("External endpoint:")
    print(recommended.external_endpoint_description if recommended else "null")
    print()
    print("ST mapping object:")
    print(recommended.visium_or_ST_data_available if recommended else "null")
    print()
    print("scRNA-seq reference:")
    print(recommended.scRNAseq_reference_available if recommended else "null")
    print()
    print("Biological application allowed preliminary:")
    print(str(bool(recommended)).lower())
    print()
    print("Large raw data downloaded:")
    print("false")
    print()
    print("CytoSPACE run:")
    print("false")
    print()
    print("SVTuner run:")
    print("false")
    print()
    print("Stage4 run:")
    print("false")
    print()
    print("Formal metrics recomputed:")
    print("false")
    print()
    print("Output directory:")
    print(str(OUT_DIR).replace("\\", "/"))
    print()
    print("Next:")
    print("Manual review before BioApp Phase 1 data inventory and feasibility")


if __name__ == "__main__":
    main()
