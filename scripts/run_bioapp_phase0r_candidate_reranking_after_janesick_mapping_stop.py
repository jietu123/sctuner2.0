from pathlib import Path
import csv
import json
from datetime import datetime

ROOT = Path(__file__).resolve().parents[1]
PHASE0_DIR = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase0_external_endpoint_admission"
PHASE1_DIR = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase1_janesick_feasibility_and_prior_risk_gate"
OUT_DIR = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase0r_candidate_reranking_after_janesick_mapping_stop"


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fields})


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def load_json(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def yn(value: str) -> bool:
    return str(value).strip().lower().startswith("yes") or str(value).strip().lower() == "true"


def rerank_previous_candidates() -> list[dict[str, str]]:
    rows = []
    for source in [
        PHASE0_DIR / "bioapp_phase0_candidate_dataset_table.csv",
        PHASE0_DIR / "bioapp_phase0_tier2_candidates.csv",
        PHASE0_DIR / "bioapp_phase0_rejected_candidates.csv",
    ]:
        rows.extend(read_csv(source))

    seen = set()
    output = []
    for row in rows:
        cid = row.get("candidate_id", "")
        if not cid or cid in seen:
            continue
        seen.add(cid)
        dataset_name = row.get("dataset_name", "")
        previous_tier = row.get("tier", "") or row.get("golden_rules_preliminary_status", "")

        is_janesick = "janesick" in cid.lower() or "janesick" in dataset_name.lower()
        endpoint_available = yn(row.get("endpoint_coordinate_or_roi_available", "")) or yn(row.get("external_endpoint_available", "")) or yn(row.get("endpoint_independent_from_SVTuner", ""))
        endpoint_independent = yn(row.get("endpoint_independent_from_SVTuner", ""))
        endpoint_mappable = yn(row.get("endpoint_spatially_registered_or_mappable", ""))
        st_available = yn(row.get("visium_or_ST_data_available", ""))
        sc_available = yn(row.get("scRNAseq_reference_available", ""))
        baseline_feasible = yn(row.get("baseline_mapping_possible", ""))
        endpoint_metric = yn(row.get("endpoint_specific_metric_possible", ""))

        if is_janesick:
            new_tier = "RETIRED"
            status = "RETIRED"
            recommendation = "Retire from current main BioApp candidate; retain only as endpoint feasibility/cautionary reference."
            blocking = "Prior formal CytoSPACE mapping mainline stopped; endpoint feasibility alone is insufficient."
            formal_baseline = "No - prior formal mapping line stopped"
            svtuner_feasible = "No - blocked by formal baseline mapping stop"
        elif endpoint_available and endpoint_independent and endpoint_mappable and st_available and sc_available and baseline_feasible and endpoint_metric:
            new_tier = "ACTIVE_TIER1"
            status = "ACTIVE_TIER1"
            recommendation = "Can enter manual Phase 1 review if data inventory confirms mapping object/reference sizes."
            blocking = "Requires fresh Phase 1 inventory; no large raw download in Phase 0R."
            formal_baseline = "Plausible from Phase 0 metadata, not executed"
            svtuner_feasible = "Plausible if baseline input contract is confirmed"
        elif endpoint_available and endpoint_independent and st_available and sc_available:
            new_tier = "TIER2_REVIEW_REQUIRED"
            status = "TIER2_REVIEW_REQUIRED"
            recommendation = "Needs manual registration/mapping feasibility review before Phase 1."
            blocking = "Endpoint or formal mapping feasibility remains unclear."
            formal_baseline = "Unclear"
            svtuner_feasible = "Unclear"
        else:
            new_tier = "REJECTED"
            status = "REJECTED"
            recommendation = "Do not advance without independent endpoint plus mapping feasibility."
            blocking = "Missing endpoint, ST mapping object, scRNA-seq reference, or endpoint-specific comparison path."
            formal_baseline = "No or not established"
            svtuner_feasible = "No or not established"

        output.append(
            {
                "candidate_id": cid,
                "dataset_name": dataset_name,
                "previous_tier": previous_tier,
                "new_tier": new_tier,
                "external_endpoint_available": str(endpoint_available),
                "endpoint_independent": str(endpoint_independent),
                "endpoint_spatially_registered_or_mappable": str(endpoint_mappable),
                "ST_mapping_object_available": str(st_available),
                "scRNAseq_reference_available": str(sc_available),
                "formal_baseline_mapping_feasible": formal_baseline,
                "SVTuner_mapping_feasible": svtuner_feasible,
                "endpoint_specific_metric_feasible": str(endpoint_metric),
                "expected_compute_feasibility_on_16GB_RAM": "Unknown until Phase 1 inventory" if status != "REJECTED" else "Not relevant",
                "main_blocking_risk": blocking,
                "candidate_status": status,
                "recommendation": recommendation,
            }
        )
    return output


def expanded_candidates() -> list[dict[str, str]]:
    # Lightweight web-screened candidate table. No large data downloaded; links are for manual Phase 1 review.
    return [
        {
            "candidate_id": "P0R_CTA_BREAST_PATHOLOGY_VISIUM_2025",
            "dataset_name": "Computational pathology annotation / breast cancer Visium CTA",
            "source_or_publication": "npj Precision Oncology 2025; dataset DOI landing page",
            "organ_or_disease": "Human breast cancer",
            "endpoint_type": "Pathology ROI / computational tissue annotation endpoint",
            "endpoint_description": "CTA-defined tumor, immune, and stromal tissue compartments from H&E-associated spatial transcriptomics samples.",
            "endpoint_source": "Independent computational pathology annotation; external to SVTuner mapping.",
            "ST_data_type": "10x Visium spatial transcriptomics",
            "sc_reference_type": "Public breast tumor scRNA-seq reference used for deconvolution in the paper",
            "same_sample_or_matched_status": "Same Visium samples with pathology annotation; scRNA-seq reference is public/matched by disease rather than necessarily same patient",
            "formal_mapping_feasibility": "Plausible - spot-level Visium object and scRNA reference are described; must verify file-level import and gene overlap in Phase 1",
            "estimated_data_size": "Moderate; dataset page reports 14 samples from 3 patients, but Phase 1 should use minimal endpoint/ST/reference files only",
            "compute_risk": "Medium; likely feasible with minimal sample/downsampled Phase 1 on 16GB RAM",
            "tier": "ACTIVE_TIER1",
            "recommendation": "Recommended first Phase 1 candidate, but do not reuse deleted legacy BCSA outputs; perform fresh inventory only.",
            "source_link_or_accession": "https://www.nature.com/articles/s41698-025-01104-3 ; https://researchdata.se/en/catalogue/dataset/2025-97",
            "notes": "Best fit because it has an independent external endpoint plus a plausible Visium/scRNA mapping route; still requires Phase 1 formal mapping feasibility gate.",
        },
        {
            "candidate_id": "P0R_SPATIAL_TOUCHSTONE_XENIUM_COSMX",
            "dataset_name": "Spatial Touchstone Xenium/CosMx benchmark",
            "source_or_publication": "Nature Biotechnology 2025; GSE277080 / Spatial Touchstone Portal",
            "organ_or_disease": "Six tissue types / imaging-based spatial benchmark",
            "endpoint_type": "Xenium/CosMx imaging endpoint",
            "endpoint_description": "Platform-derived imaging-based cell annotations and spatial coordinates across many profiles.",
            "endpoint_source": "Spatial Touchstone imaging-based ST repository",
            "ST_data_type": "Xenium and CosMx; Visium-like mapping object not confirmed in Phase 0R",
            "sc_reference_type": "Annotations can be imputed from single-cell datasets according to source description; exact reference needs review",
            "same_sample_or_matched_status": "Standardized multi-site samples; exact paired Visium/ST mapping object not confirmed",
            "formal_mapping_feasibility": "Unclear - strong endpoint resource but missing confirmed spot-level ST mapping object for CytoSPACE",
            "estimated_data_size": "Large; source describes 203 spatial profiles / large imaging repository",
            "compute_risk": "High unless a minimal subset is identified",
            "tier": "TIER2_REVIEW_REQUIRED",
            "recommendation": "Manual review only; do not advance unless a spot-level mapping object and scRNA reference are confirmed.",
            "source_link_or_accession": "https://www.nature.com/articles/s41587-025-02811-9 ; https://www.omicsdi.org/dataset/geo/GSE277080",
            "notes": "Useful backup for endpoint-rich screening, but not yet a clean CytoSPACE/SVTuner biological application candidate.",
        },
        {
            "candidate_id": "P0R_COSMX_NSCLC_PUBLIC",
            "dataset_name": "CosMx SMI FFPE NSCLC public dataset",
            "source_or_publication": "Bruker/NanoString public CosMx SMI NSCLC dataset",
            "organ_or_disease": "Human non-small-cell lung cancer FFPE",
            "endpoint_type": "CosMx cell annotation endpoint",
            "endpoint_description": "Single-cell spatial atlas / cell typing endpoint from CosMx SMI FFPE NSCLC samples.",
            "endpoint_source": "CosMx public dataset",
            "ST_data_type": "CosMx imaging data; no paired Visium/ST mapping object verified",
            "sc_reference_type": "Not confirmed for CytoSPACE mapping in Phase 0R",
            "same_sample_or_matched_status": "CosMx samples public; paired spot-level ST reference not established",
            "formal_mapping_feasibility": "No - endpoint exists but formal baseline CytoSPACE mapping route is not established",
            "estimated_data_size": "Large public imaging dataset",
            "compute_risk": "High; missing mapping object is the main blocker",
            "tier": "REJECTED",
            "recommendation": "Reject for current BioApp design unless paired ST + scRNA reference are separately identified.",
            "source_link_or_accession": "https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/",
            "notes": "Endpoint-rich but violates Golden Rule 18 unless a formal baseline mapping object is identified.",
        },
        {
            "candidate_id": "P0R_VISIUM_HD_COLON_10X",
            "dataset_name": "10x Human Colon Cancer Visium HD",
            "source_or_publication": "10x Genomics HumanColonCancer_VisiumHD GitHub",
            "organ_or_disease": "Human colorectal cancer",
            "endpoint_type": "Potential pathology/region endpoint, not yet independently frozen",
            "endpoint_description": "High-definition Visium spatial profile may support tissue-region analysis, but independent external endpoint is not established in Phase 0R.",
            "endpoint_source": "No independent endpoint confirmed; possible H&E/pathology endpoint would need separate source",
            "ST_data_type": "Visium HD",
            "sc_reference_type": "Needs matched or public CRC scRNA reference identification",
            "same_sample_or_matched_status": "Not established",
            "formal_mapping_feasibility": "Unclear - mapping object likely available, independent endpoint missing",
            "estimated_data_size": "Potentially large; downsampled feasibility mode likely required",
            "compute_risk": "Medium-high",
            "tier": "TIER2_REVIEW_REQUIRED",
            "recommendation": "Reserve only if independent pathology ROI endpoint can be obtained before mapping.",
            "source_link_or_accession": "https://github.com/10XGenomics/HumanColonCancer_VisiumHD",
            "notes": "Mapping data alone is insufficient; independent endpoint must be added before Phase 1.",
        },
    ]


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    phase0_summary = load_json(PHASE0_DIR / "bioapp_phase0_summary.json")
    phase1_summary = load_json(PHASE1_DIR / "bioapp_phase1_summary.json")

    retirement = {
        "candidate_id": "T1_XENIUM_JANESICK_BREAST_2023",
        "candidate_name": "Janesick Xenium-Visium breast cancer",
        "candidate_status": "RETIRED_FROM_MAIN_BIOAPP_CANDIDATE",
        "reason": "Prior formal CytoSPACE mapping mainline stopped.",
        "not_deleted_as_scientific_reference": True,
        "retired_only_as_current_main_svtuner_biological_application_candidate": True,
        "endpoint_roi_feasibility_informative": True,
        "formal_mapping_branch_reliable_enough_for_current_experiment": False,
        "allowed_future_roles": ["endpoint feasibility reference", "ROI validation precedent", "negative candidate / cautionary case"],
        "disallowed_future_roles": ["main BioApp experiment", "formal SVTuner mapping biological application", "endpoint-specific mapping comparison"],
        "historical_status": ["Janesick Phase 3.1: formal stop", "Janesick Phase 3.2: formal mapping mainline stopped"],
    }
    write_text(
        OUT_DIR / "bioapp_phase0r_janesick_retirement_record.txt",
        "Janesick is not deleted as a scientific reference.\n"
        "Janesick is retired only as the current main SVTuner biological application candidate.\n"
        "Reason: prior formal CytoSPACE mapping mainline stopped.\n"
        "Endpoint/ROI feasibility alone is insufficient because the current experiment requires formal baseline mapping + SVTuner mapping + endpoint-specific comparison.\n\n"
        "Allowed future roles:\n- endpoint feasibility reference\n- ROI validation precedent\n- negative candidate / cautionary case\n\n"
        "Disallowed roles:\n- main BioApp experiment\n- formal SVTuner mapping biological application\n- endpoint-specific mapping comparison\n",
    )
    write_text(OUT_DIR / "bioapp_phase0r_janesick_status_update.json", json.dumps(retirement, indent=2),)

    golden_rules_text = """Golden Rules v2.1 update\n\nGolden Rule 18:\nNo reliable formal baseline mapping, no SVTuner biological application.\n\nA candidate cannot enter BioApp Phase 2 unless formal baseline mapping feasibility is supported.\nThe candidate must support:\n1. identifiable ST mapping object\n2. usable scRNA-seq reference\n3. sufficient gene intersection\n4. compatible preprocessing\n5. baseline CytoSPACE feasibility\n6. SVTuner-enhanced mapping feasibility\n7. identical input basis for baseline and SVTuner\n8. endpoint-specific comparison after mapping\n\nIf external endpoint is feasible but mapping mainline is not feasible, the candidate must be downgraded to:\nendpoint feasibility only, not main biological application.\n"""
    write_text(OUT_DIR / "bioapp_phase0r_golden_rules_v2_1_update.txt", golden_rules_text)
    schema = {
        "version": "Golden Rules v2.1",
        "new_rule_id": 18,
        "new_rule": "No reliable formal baseline mapping, no SVTuner biological application.",
        "phase2_entry_required_checks": [
            "identifiable ST mapping object",
            "usable scRNA-seq reference",
            "sufficient gene intersection",
            "compatible preprocessing",
            "baseline CytoSPACE feasibility",
            "SVTuner-enhanced mapping feasibility",
            "identical input basis for baseline and SVTuner",
            "endpoint-specific comparison after mapping",
        ],
        "downgrade_if_endpoint_only": "endpoint feasibility only; not main biological application",
    }
    write_text(OUT_DIR / "bioapp_phase0r_golden_rules_v2_1_check_schema.json", json.dumps(schema, indent=2))

    reranked = rerank_previous_candidates()
    rerank_fields = [
        "candidate_id",
        "dataset_name",
        "previous_tier",
        "new_tier",
        "external_endpoint_available",
        "endpoint_independent",
        "endpoint_spatially_registered_or_mappable",
        "ST_mapping_object_available",
        "scRNAseq_reference_available",
        "formal_baseline_mapping_feasible",
        "SVTuner_mapping_feasible",
        "endpoint_specific_metric_feasible",
        "expected_compute_feasibility_on_16GB_RAM",
        "main_blocking_risk",
        "candidate_status",
        "recommendation",
    ]
    write_csv(OUT_DIR / "bioapp_phase0r_candidate_reranking_table.csv", reranked, rerank_fields)

    expanded = expanded_candidates()
    expanded_fields = [
        "candidate_id",
        "dataset_name",
        "source_or_publication",
        "organ_or_disease",
        "endpoint_type",
        "endpoint_description",
        "endpoint_source",
        "ST_data_type",
        "sc_reference_type",
        "same_sample_or_matched_status",
        "formal_mapping_feasibility",
        "estimated_data_size",
        "compute_risk",
        "tier",
        "recommendation",
        "source_link_or_accession",
        "notes",
    ]
    write_csv(OUT_DIR / "bioapp_phase0r_expanded_candidate_search_table.csv", expanded, expanded_fields)

    recommended = [row for row in expanded if row["tier"] == "ACTIVE_TIER1"]
    tier2 = [row for row in expanded if row["tier"] == "TIER2_REVIEW_REQUIRED"]
    rejected = [row for row in expanded if row["tier"] in {"REJECTED", "RETIRED"}]
    rec_fields = expanded_fields + ["priority"]
    rec_rows = []
    for i, row in enumerate(recommended, start=1):
        rec_rows.append({**row, "priority": f"Priority {i}"})
    write_csv(OUT_DIR / "bioapp_phase0r_recommended_candidates.csv", rec_rows, rec_fields)

    if recommended:
        decision = "PASS"
        recommended_name = recommended[0]["dataset_name"]
        rec_report = (
            "BioApp Phase 0R recommendation\n\n"
            f"Recommended candidate: {recommended_name}\n"
            "Reason: it currently best satisfies external endpoint + ST + sc reference + plausible formal mapping feasibility + manageable minimal Phase 1 path.\n\n"
            "Important constraint: this is not a mapping success claim. Phase 1 must still verify file inventory, endpoint/ST alignment, scRNA reference import, gene intersection, and baseline/SVTuner input compatibility.\n\n"
            "Janesick is retired from the main BioApp candidate role because the prior formal CytoSPACE mapping mainline stopped.\n"
        )
    elif tier2:
        decision = "REVIEW_REQUIRED"
        recommended_name = None
        rec_report = "No active Tier 1 candidate found. Tier 2 candidates require manual mapping/registration review.\n"
    else:
        decision = "FAIL"
        recommended_name = None
        rec_report = "No active Tier 1 candidate found. Manual dataset search required.\n"
    write_text(OUT_DIR / "bioapp_phase0r_recommendation_report.txt", rec_report)

    summary = {
        "phase": "BioApp Phase 0R - candidate reranking after Janesick formal-mapping stop",
        "decision": decision,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "janesick_retired_as_main_bioapp_candidate": True,
        "janesick_retirement_reason": "prior formal CytoSPACE mapping mainline stopped",
        "golden_rules_v2_1_mapping_feasibility_gate_added": True,
        "previous_candidates_re_ranked": len(reranked),
        "expanded_candidates_screened": len(expanded),
        "active_tier1_candidates": len(recommended),
        "tier2_candidates": len(tier2),
        "rejected_or_retired_candidates": len(rejected) + len([row for row in reranked if row["candidate_status"] in {"RETIRED", "REJECTED"}]),
        "recommended_candidate": recommended_name,
        "large_raw_data_downloaded": False,
        "cytospace_run": False,
        "svtuner_run": False,
        "stage4_run": False,
        "formal_metrics_recomputed": False,
        "phase0_input_available": bool(phase0_summary),
        "phase1_input_available": bool(phase1_summary),
        "output_directory": str(OUT_DIR),
        "next": "Manual review before BioApp Phase 1 for the newly recommended candidate",
    }
    write_text(OUT_DIR / "bioapp_phase0r_summary.json", json.dumps(summary, indent=2))
    readme = """BioApp Phase 0R\n\nPurpose: retire Janesick as the active main BioApp candidate after prior formal mapping stop, add Golden Rule 18, rerank previous candidates, and lightly screen expanded candidates.\n\nNo CytoSPACE, SVTuner, Stage4, formal metric recomputation, large raw download, FASTQ download, or manuscript writing was performed.\n"""
    write_text(OUT_DIR / "bioapp_phase0r_readme.txt", readme)
    write_text(OUT_DIR / "bioapp_phase0r_decision.txt", decision + "\n")

    print("BioApp Phase 0R completed.")
    print("\nDecision:")
    print(decision)
    print("\nJanesick retired as main BioApp candidate:")
    print("true")
    print("\nReason:")
    print("prior formal CytoSPACE mapping mainline stopped")
    print("\nGolden Rules v2.1 mapping-feasibility gate added:")
    print("true")
    print("\nPrevious candidates re-ranked:")
    print(len(reranked))
    print("\nExpanded candidates screened:")
    print(len(expanded))
    print("\nActive Tier 1 candidates:")
    print(len(recommended))
    print("\nTier 2 candidates:")
    print(len(tier2))
    print("\nRejected / retired candidates:")
    print(summary["rejected_or_retired_candidates"])
    print("\nRecommended candidate:")
    print(recommended_name if recommended_name else "null")
    print("\nLarge raw data downloaded:")
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
    print("Manual review before BioApp Phase 1 for the newly recommended candidate")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
