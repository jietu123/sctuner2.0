#!/usr/bin/env python
"""BioApp Phase 2C endpoint freeze with validated CTA-to-spot mapping.

This stage freezes CTA-derived external endpoints at Visium spot level using
only the validated Phase 2B-R2 CTA-to-spot mapping outputs. It must not run
CytoSPACE, SVTuner, Stage4, or any downstream formal metric.
"""

from __future__ import annotations

import json
from collections import deque
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PHASE2BR2_DIR = (
    ROOT
    / "visualizations"
    / "bioapp_phase2br2_explicit_rscript_cta_align_locale_safe_dry_run"
)
OUT_DIR = (
    ROOT
    / "visualizations"
    / "bioapp_phase2c_endpoint_freeze_with_validated_cta_to_spot_mapping"
)

PHASE_NAME = "BioApp Phase 2C endpoint freeze with validated CTA-to-spot mapping"
STAGE_TYPE = "external endpoint freeze / biological application candidate preparation"

CLASSES = ["Tumor", "Immune cells", "Stroma"]
CLASS_SAFE = {
    "Tumor": "Tumor",
    "Immune cells": "Immune_cells",
    "Stroma": "Stroma",
}
COLORS_CLASS = {
    "Tumor": "#C44E52",
    "Immune cells": "#4C72B0",
    "Stroma": "#55A868",
    "unmapped": "#BDBDBD",
}
COLORS_STATUS = {
    "positive": "#D55E00",
    "negative": "#C7DCEF",
    "ambiguous": "#FDE0C5",
    "excluded": "#D9D9D9",
}

FREEZE_RULE = {
    "min_total_CTA_objects": 2,
    "min_class_count_positive": 2,
    "positive_fraction_threshold": 0.60,
    "negative_fraction_threshold": 0.20,
}

ALLOWED_CLAIMS = [
    "CTA-derived endpoint was frozen at Visium spot level.",
    "CTA-to-spot mapping enabled endpoint-positive, endpoint-negative, ambiguous, and excluded spot definitions.",
]
DISALLOWED_CLAIMS = [
    "SVTuner improves endpoint recovery.",
    "SVTuner prevents contradicted interpretation.",
    "Baseline creates false niche calls.",
    "Biological application completed.",
    "Biological discovery made.",
]


def rel(path: Path) -> str:
    try:
        return path.relative_to(ROOT).as_posix()
    except ValueError:
        return path.as_posix()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def fail_outputs(reason: str, inputs: dict[str, Path] | None = None) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary = {
        "phase": PHASE_NAME,
        "decision": "FAIL",
        "reason": reason,
        "stage_type": STAGE_TYPE,
        "input_files": {k: rel(v) for k, v in (inputs or {}).items()},
        "endpoint_frozen": False,
        "ready_for_phase3": False,
        "validated_CTA_to_spot_mapping_used": False,
        "CTA_mapping_recomputed": False,
        "CTA_align_R_modified": False,
        "expression_markers_used_to_define_endpoint": False,
        "mapping_outputs_used_to_define_endpoint": False,
        "CytoSPACE_run": False,
        "SVTuner_run": False,
        "Stage4_run": False,
        "formal_metrics_recomputed": False,
        "allowed_claims": ALLOWED_CLAIMS,
        "disallowed_claims": DISALLOWED_CLAIMS,
    }
    write_json(OUT_DIR / "bioapp_phase2c_endpoint_freeze_summary.json", summary)
    write_json(
        OUT_DIR / "bioapp_phase2c_golden_rules_v2_1_check.json",
        {
            "stage_type": STAGE_TYPE,
            "biological_question_defined": True,
            "external_endpoint_predefined": True,
            "endpoint_independent_from_SVTuner": True,
            "endpoint_spatially_registered": False,
            "baseline_comparison_available": False,
            "endpoint_specific_improvement_defined": False,
            "endpoint_specific_quantitative_metric_available": False,
            "endpoint_baseline_SVTuner_spatial_comparison_available": False,
            "interpretation_boundary_defined": True,
            "biological_application_allowed": False,
            "allowed_claim_level": "endpoint freeze only",
            "decision": "FAIL",
            "reason": reason,
        },
    )
    write_text(
        OUT_DIR / "decision.txt",
        "\n".join(
            [
                "BioApp Phase 2C - endpoint freeze with validated CTA-to-spot mapping",
                "",
                "Decision: FAIL",
                "",
                "Endpoint frozen: false",
                "Primary endpoint: null",
                "Ready for Phase 3: false",
                "",
                f"Reason: {reason}",
                "",
                "Next:",
                "Stop. Fix input mapping or boundary violation before retrying Phase 2C.",
                "",
            ]
        ),
    )
    print("BioApp Phase 2C completed.")
    print("\nDecision:\nFAIL")
    print("\nEndpoint frozen:\nfalse")
    print("\nPrimary endpoint:\nnull")
    print("\nNext:\nStop. Fix input mapping or boundary violation before retrying Phase 2C.")


def find_inputs() -> dict[str, Path]:
    return {
        "summary": PHASE2BR2_DIR / "bioapp_phase2br2_summary.json",
        "coordinates": PHASE2BR2_DIR / "bioapp_phase2br2_st_data_coordinates.csv",
        "mapping": PHASE2BR2_DIR / "bioapp_phase2br2_cta_to_spot_mapping_example.csv",
        "invocation_record": PHASE2BR2_DIR / "bioapp_phase2br2_cta_align_invocation_record.json",
        "st_inventory": PHASE2BR2_DIR / "bioapp_phase2br2_st_data_object_inventory.json",
        "golden_rules": PHASE2BR2_DIR / "bioapp_phase2br2_golden_rules_v2_1_check.json",
    }


def validate_inputs(paths: dict[str, Path]) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame]:
    missing = [name for name, path in paths.items() if name in {"summary", "coordinates", "mapping"} and not path.exists()]
    if missing:
        raise RuntimeError(f"validated CTA-to-spot mapping not found; missing: {', '.join(missing)}")

    summary = json.loads(paths["summary"].read_text(encoding="utf-8"))
    if summary.get("decision") != "PASS":
        raise RuntimeError("Phase 2B-R2 decision not PASS")
    if not summary.get("CTA_to_spot_mapping_generated"):
        raise RuntimeError("Phase 2B-R2 did not generate validated CTA-to-spot mapping")
    if summary.get("endpoint_frozen"):
        raise RuntimeError("Phase 2B-R2 unexpectedly reports endpoint already frozen")
    if summary.get("CytoSPACE_run") or summary.get("SVTuner_run") or summary.get("Stage4_run"):
        raise RuntimeError("boundary violation detected in Phase 2B-R2 summary")
    if summary.get("expression_markers_used_to_define_endpoint") or summary.get("mapping_outputs_used_to_define_endpoint"):
        raise RuntimeError("endpoint definition boundary violation detected in Phase 2B-R2 summary")

    coords = pd.read_csv(paths["coordinates"])
    mapping = pd.read_csv(paths["mapping"])
    if "barcode" not in coords.columns or "barcode" not in mapping.columns:
        raise RuntimeError("input files inconsistent: barcode column missing")
    required_mapping_cols = {"CTA_class_or_label", "CTA_score_or_count"}
    if not required_mapping_cols.issubset(mapping.columns):
        raise RuntimeError("input files inconsistent: mapping class/count columns missing")
    return summary, coords, mapping


def build_spot_level_composition(coords: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    coords = coords.copy()
    coords["barcode"] = coords["barcode"].astype(str)
    mapping = mapping.copy()
    mapping["barcode"] = mapping["barcode"].astype(str)
    mapping["CTA_score_or_count"] = pd.to_numeric(mapping["CTA_score_or_count"], errors="coerce").fillna(0)
    mapping = mapping[mapping["CTA_class_or_label"].isin(CLASSES)]

    pivot = (
        mapping.pivot_table(
            index="barcode",
            columns="CTA_class_or_label",
            values="CTA_score_or_count",
            aggfunc="sum",
            fill_value=0,
        )
        .reset_index()
        .rename_axis(None, axis=1)
    )
    for cls in CLASSES:
        if cls not in pivot.columns:
            pivot[cls] = 0

    keep_cols = [
        "barcode",
        "tissue",
        "row",
        "col",
        "imagerow",
        "imagecol",
        "pxl_col_in_fullres",
        "pxl_row_in_fullres",
    ]
    spot = coords[keep_cols].merge(pivot[["barcode", *CLASSES]], on="barcode", how="left")

    for cls in CLASSES:
        safe = CLASS_SAFE[cls]
        spot[cls] = pd.to_numeric(spot[cls], errors="coerce").fillna(0).astype(int)
        spot[f"{safe}_count"] = spot[cls]
    spot["total_CTA_objects"] = spot[[CLASS_SAFE[cls] + "_count" for cls in CLASSES]].sum(axis=1).astype(int)
    spot["CTA_mapped"] = spot["total_CTA_objects"] > 0

    for cls in CLASSES:
        safe = CLASS_SAFE[cls]
        spot[f"{safe}_fraction"] = np.where(
            spot["total_CTA_objects"] > 0,
            spot[f"{safe}_count"] / spot["total_CTA_objects"],
            0.0,
        )

    fraction_cols = [CLASS_SAFE[cls] + "_fraction" for cls in CLASSES]
    count_cols = [CLASS_SAFE[cls] + "_count" for cls in CLASSES]
    dominant_index = spot[fraction_cols].to_numpy().argmax(axis=1)
    dominant = np.array(CLASSES, dtype=object)[dominant_index]
    dominant_fraction = spot[fraction_cols].to_numpy().max(axis=1)
    dominant[~spot["CTA_mapped"].to_numpy()] = "unmapped"
    dominant_fraction[~spot["CTA_mapped"].to_numpy()] = 0.0
    spot["dominant_CTA_class"] = dominant
    spot["dominant_CTA_fraction"] = dominant_fraction
    spot["unmatched_reason"] = np.where(spot["CTA_mapped"], "", "no_CTA_object_captured")

    output_cols = [
        "barcode",
        "tissue",
        "row",
        "col",
        "imagerow",
        "imagecol",
        "pxl_col_in_fullres",
        "pxl_row_in_fullres",
        "total_CTA_objects",
        *count_cols,
        *fraction_cols,
        "dominant_CTA_class",
        "dominant_CTA_fraction",
        "CTA_mapped",
        "unmatched_reason",
    ]
    return spot[output_cols].copy()


def assign_endpoint_status(spot: pd.DataFrame) -> pd.DataFrame:
    out = spot.copy()
    for cls in CLASSES:
        safe = CLASS_SAFE[cls]
        count_col = f"{safe}_count"
        frac_col = f"{safe}_fraction"
        status = np.full(len(out), "ambiguous", dtype=object)
        mapped = out["CTA_mapped"].to_numpy()
        total_ok = out["total_CTA_objects"].to_numpy() >= FREEZE_RULE["min_total_CTA_objects"]
        positive = (
            mapped
            & total_ok
            & (out[count_col].to_numpy() >= FREEZE_RULE["min_class_count_positive"])
            & (out[frac_col].to_numpy() >= FREEZE_RULE["positive_fraction_threshold"])
        )
        negative = mapped & total_ok & (out[frac_col].to_numpy() <= FREEZE_RULE["negative_fraction_threshold"])
        status[positive] = "positive"
        status[negative & ~positive] = "negative"
        status[~mapped] = "excluded"
        out[f"{safe}_endpoint_status"] = status
    return out


def build_knn_adjacency(df: pd.DataFrame, k: int = 6) -> list[set[int]]:
    coords = df[["imagecol", "imagerow"]].to_numpy(dtype=float)
    n = coords.shape[0]
    adjacency: list[set[int]] = [set() for _ in range(n)]
    try:
        from scipy.spatial import cKDTree  # type: ignore

        tree = cKDTree(coords)
        _, idx = tree.query(coords, k=min(k + 1, n))
        for i, neighbors in enumerate(np.atleast_2d(idx)):
            for j in neighbors:
                if i == int(j):
                    continue
                adjacency[i].add(int(j))
                adjacency[int(j)].add(i)
    except Exception:
        diff = coords[:, None, :] - coords[None, :, :]
        dist = np.sqrt(np.sum(diff * diff, axis=2))
        np.fill_diagonal(dist, np.inf)
        idx = np.argsort(dist, axis=1)[:, : min(k, n - 1)]
        for i in range(n):
            for j in idx[i]:
                adjacency[i].add(int(j))
                adjacency[int(j)].add(i)
    return adjacency


def positive_component_stats(status: pd.Series, adjacency: list[set[int]]) -> dict[str, Any]:
    positive = set(np.flatnonzero(status.to_numpy() == "positive").tolist())
    if not positive:
        return {
            "number_of_positive_components": 0,
            "largest_positive_component_size": 0,
            "largest_positive_component_fraction": 0.0,
        }
    seen: set[int] = set()
    sizes: list[int] = []
    for start in positive:
        if start in seen:
            continue
        queue: deque[int] = deque([start])
        seen.add(start)
        size = 0
        while queue:
            node = queue.popleft()
            size += 1
            for nxt in adjacency[node]:
                if nxt in positive and nxt not in seen:
                    seen.add(nxt)
                    queue.append(nxt)
        sizes.append(size)
    largest = max(sizes) if sizes else 0
    return {
        "number_of_positive_components": len(sizes),
        "largest_positive_component_size": int(largest),
        "largest_positive_component_fraction": float(largest / len(positive)) if positive else 0.0,
    }


def endpoint_qc(df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, dict[str, Any]]]:
    adjacency = build_knn_adjacency(df)
    mapped_count = int(df["CTA_mapped"].sum())
    rows: list[dict[str, Any]] = []
    qc: dict[str, dict[str, Any]] = {}
    for cls in CLASSES:
        safe = CLASS_SAFE[cls]
        status_col = f"{safe}_endpoint_status"
        frac_col = f"{safe}_fraction"
        positive_mask = df[status_col] == "positive"
        negative_mask = df[status_col] == "negative"
        ambiguous_mask = df[status_col] == "ambiguous"
        excluded_mask = df[status_col] == "excluded"
        comp = positive_component_stats(df[status_col], adjacency)
        n_positive = int(positive_mask.sum())
        n_negative = int(negative_mask.sum())
        n_ambiguous = int(ambiguous_mask.sum())
        n_excluded = int(excluded_mask.sum())
        positive_fracs = df.loc[positive_mask, frac_col]
        dominant_purity = (
            float((df.loc[positive_mask, "dominant_CTA_class"] == cls).mean()) if n_positive > 0 else 0.0
        )
        eligible = n_positive >= 30 and n_negative >= 30 and comp["largest_positive_component_size"] >= 10
        row = {
            "endpoint_class": cls,
            "n_positive": n_positive,
            "n_negative": n_negative,
            "n_ambiguous": n_ambiguous,
            "n_excluded": n_excluded,
            "positive_fraction_among_mapped": float(n_positive / mapped_count) if mapped_count else 0.0,
            "negative_fraction_among_mapped": float(n_negative / mapped_count) if mapped_count else 0.0,
            "median_C_fraction_positive": float(positive_fracs.median()) if n_positive > 0 else 0.0,
            "mean_C_fraction_positive": float(positive_fracs.mean()) if n_positive > 0 else 0.0,
            "dominant_class_purity_positive": dominant_purity,
            **comp,
            "eligible_endpoint": bool(eligible),
        }
        rows.append(row)
        qc[cls] = {k: (bool(v) if isinstance(v, np.bool_) else v) for k, v in row.items() if k != "endpoint_class"}
    return pd.DataFrame(rows), qc


def select_primary_endpoint(qc: dict[str, dict[str, Any]]) -> str | None:
    for cls in ["Immune cells", "Tumor", "Stroma"]:
        if qc.get(cls, {}).get("eligible_endpoint"):
            return cls
    return None


def add_primary_status(df: pd.DataFrame, primary: str | None) -> pd.DataFrame:
    out = df.copy()
    if primary is None:
        out["primary_endpoint_status"] = "not_selected"
        out["primary_endpoint_class"] = ""
        return out
    safe = CLASS_SAFE[primary]
    out["primary_endpoint_status"] = out[f"{safe}_endpoint_status"]
    out["primary_endpoint_class"] = primary
    return out


def write_spot_lists(df: pd.DataFrame, primary: str | None) -> None:
    for cls in CLASSES:
        safe = CLASS_SAFE[cls]
        status_col = f"{safe}_endpoint_status"
        for status in ["positive", "negative", "ambiguous", "excluded"]:
            path = OUT_DIR / f"endpoint_{safe}_{status}_spots.txt"
            values = df.loc[df[status_col] == status, "barcode"].astype(str).tolist()
            write_text(path, "\n".join(values) + ("\n" if values else ""))
    if primary is not None:
        safe = CLASS_SAFE[primary]
        status_col = f"{safe}_endpoint_status"
        for status in ["positive", "negative", "ambiguous", "excluded"]:
            values = df.loc[df[status_col] == status, "barcode"].astype(str).tolist()
            write_text(OUT_DIR / f"primary_endpoint_{status}_spots.txt", "\n".join(values) + ("\n" if values else ""))


def setup_spatial_axis(ax: plt.Axes) -> None:
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.invert_yaxis()
    for spine in ax.spines.values():
        spine.set_visible(False)


def save_spatial(fig: plt.Figure, basename: str) -> None:
    for ext in ["pdf", "svg"]:
        fig.savefig(OUT_DIR / f"{basename}.{ext}", bbox_inches="tight")
    plt.close(fig)


def plot_numeric_map(df: pd.DataFrame, value_col: str, title: str, basename: str, cmap: str = "magma") -> None:
    fig, ax = plt.subplots(figsize=(6.2, 6.0))
    sc = ax.scatter(
        df["imagecol"],
        df["imagerow"],
        c=df[value_col],
        s=13,
        cmap=cmap,
        linewidths=0,
    )
    setup_spatial_axis(ax)
    ax.set_title(title, fontsize=11, weight="bold")
    cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(value_col, fontsize=8)
    save_spatial(fig, basename)


def plot_dominant_class_map(df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 6.0))
    for label in ["unmapped", *CLASSES]:
        mask = df["dominant_CTA_class"] == label
        if not mask.any():
            continue
        ax.scatter(
            df.loc[mask, "imagecol"],
            df.loc[mask, "imagerow"],
            s=13,
            color=COLORS_CLASS[label],
            label=label,
            linewidths=0,
        )
    setup_spatial_axis(ax)
    ax.set_title("CTA dominant class", fontsize=11, weight="bold")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=8)
    save_spatial(fig, "cta_dominant_class_map")


def plot_status_map(df: pd.DataFrame, status_col: str, title: str, basename: str) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 6.0))
    for status in ["excluded", "negative", "ambiguous", "positive"]:
        mask = df[status_col] == status
        if not mask.any():
            continue
        ax.scatter(
            df.loc[mask, "imagecol"],
            df.loc[mask, "imagerow"],
            s=13,
            color=COLORS_STATUS[status],
            label=status,
            linewidths=0,
        )
    setup_spatial_axis(ax)
    ax.set_title(title, fontsize=11, weight="bold")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=8)
    save_spatial(fig, basename)


def generate_maps(df: pd.DataFrame, primary: str | None) -> list[Path]:
    before = set(OUT_DIR.glob("*"))
    plot_numeric_map(df, "total_CTA_objects", "CTA total object count", "cta_total_object_count_map", cmap="viridis")
    plot_dominant_class_map(df)
    for cls in CLASSES:
        safe = CLASS_SAFE[cls]
        plot_numeric_map(
            df,
            f"{safe}_fraction",
            f"CTA {cls} fraction",
            f"cta_{safe}_fraction_map",
            cmap="magma",
        )
        plot_status_map(
            df,
            f"{safe}_endpoint_status",
            f"Endpoint status: {cls}",
            f"endpoint_{safe}_status_map",
        )
    if primary is not None:
        safe = CLASS_SAFE[primary]
        plot_status_map(
            df,
            f"{safe}_endpoint_status",
            f"Primary endpoint status: {primary}",
            "primary_endpoint_status_map",
        )
    after = set(OUT_DIR.glob("*"))
    return sorted(after - before)


def write_manifest(paths: list[Path], extra_rows: list[dict[str, str]]) -> None:
    rows = []
    for path in sorted(paths):
        if path.is_file():
            rows.append(
                {
                    "file": rel(path),
                    "type": path.suffix.lstrip(".") or "file",
                    "description": "BioApp Phase 2C generated output",
                    "created_by_phase": "BioApp Phase 2C",
                }
            )
    rows.extend(extra_rows)
    pd.DataFrame(rows).drop_duplicates(subset=["file"], keep="last").to_csv(OUT_DIR / "manifest.csv", index=False)


def write_readme(primary: str | None, decision: str) -> None:
    readme = f"""# BioApp Phase 2C Endpoint Freeze

Purpose: freeze CTA-derived external endpoint status at Visium spot level using the validated Phase 2B-R2 CTA-to-spot mapping.

Input source: `{rel(PHASE2BR2_DIR)}`.

Endpoint freeze rule:
- positive: CTA_mapped is true, total_CTA_objects >= 2, class_count >= 2, and class_fraction >= 0.60
- negative: CTA_mapped is true, total_CTA_objects >= 2, and class_fraction <= 0.20
- ambiguous: CTA_mapped is true and neither positive nor negative
- excluded: CTA_mapped is false

Decision: {decision}

Primary endpoint: {primary if primary is not None else "null"}

Allowed interpretation:
- CTA-derived endpoint was frozen at Visium spot level.
- CTA-to-spot mapping enabled endpoint-positive, endpoint-negative, ambiguous, and excluded spot definitions.

Disallowed interpretation:
- SVTuner improves endpoint recovery.
- SVTuner prevents contradicted interpretation.
- Baseline creates false niche calls.
- Biological application completed.
- Biological discovery made.

Next condition:
- PASS: BioApp Phase 3 - reference audit and formal baseline feasibility for CTA-defined endpoint.
- REVIEW_REQUIRED: review endpoint thresholds, CTA spot composition, and spatial continuity before Phase 3.
- FAIL: fix input mapping or boundary violation before retrying Phase 2C.
"""
    write_text(OUT_DIR / "README.md", readme)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    inputs = find_inputs()
    try:
        phase2_summary, coords, mapping = validate_inputs(inputs)
    except Exception as exc:  # noqa: BLE001
        fail_outputs(str(exc), inputs)
        return 1

    composition = build_spot_level_composition(coords, mapping)
    if len(composition) != int(phase2_summary.get("number_of_spots_in_ST_data", len(composition))):
        fail_outputs("spot count mismatch cannot be resolved", inputs)
        return 1

    endpoint_df = assign_endpoint_status(composition)
    qc_table, qc = endpoint_qc(endpoint_df)
    primary = select_primary_endpoint(qc)
    endpoint_df = add_primary_status(endpoint_df, primary)

    decision = "PASS" if primary is not None else "REVIEW_REQUIRED"
    endpoint_frozen = decision == "PASS"
    ready_for_phase3 = decision == "PASS"
    reason = "" if decision == "PASS" else "no CTA compartment endpoint passed minimum freeze QC"

    composition.to_csv(OUT_DIR / "spot_level_CTA_composition.csv", index=False)
    endpoint_df.to_csv(OUT_DIR / "spot_level_endpoint_freeze.csv", index=False)
    qc_table.to_csv(OUT_DIR / "endpoint_eligibility_qc.csv", index=False)

    write_spot_lists(endpoint_df, primary)
    map_paths = generate_maps(endpoint_df, primary)

    st_spots_total = int(len(composition))
    cta_mapped_spots = int(composition["CTA_mapped"].sum())
    unmatched_spots = int((~composition["CTA_mapped"]).sum())
    cta_objects_captured = int(composition["total_CTA_objects"].sum())

    endpoint_definition = {
        "endpoint_type": "CTA-defined compartment endpoint",
        "endpoint_classes": CLASSES,
        "primary_endpoint": primary,
        "freeze_rule": FREEZE_RULE,
        "primary_endpoint_recommendation_rule": [
            "Immune cells if eligible",
            "else Tumor if eligible",
            "else Stroma if eligible",
            "else null",
        ],
        "selection_independent_from_mapping_or_SVTuner_outputs": True,
    }
    write_json(OUT_DIR / "endpoint_definition.json", endpoint_definition)

    freeze_manifest = {
        "phase": PHASE_NAME,
        "decision": decision,
        "endpoint_frozen": endpoint_frozen,
        "primary_endpoint": primary,
        "spot_level_endpoint_freeze": rel(OUT_DIR / "spot_level_endpoint_freeze.csv"),
        "spot_level_CTA_composition": rel(OUT_DIR / "spot_level_CTA_composition.csv"),
        "endpoint_eligibility_qc": rel(OUT_DIR / "endpoint_eligibility_qc.csv"),
        "endpoint_definition": rel(OUT_DIR / "endpoint_definition.json"),
    }
    write_json(OUT_DIR / "endpoint_freeze_manifest.json", freeze_manifest)

    summary = {
        "phase": PHASE_NAME,
        "decision": decision,
        "reason": reason,
        "stage_type": STAGE_TYPE,
        "input_phase2br2_decision": phase2_summary.get("decision"),
        "validated_CTA_to_spot_mapping_used": True,
        "CTA_mapping_recomputed": False,
        "CTA_align_R_modified": False,
        "ST_spots_total": st_spots_total,
        "CTA_mapped_spots": cta_mapped_spots,
        "unmatched_spots": unmatched_spots,
        "CTA_objects_total": int(phase2_summary.get("CTA_objects_total", 0)),
        "CTA_objects_captured": cta_objects_captured,
        "CTA_objects_not_captured": int(phase2_summary.get("CTA_objects_not_captured", 0)),
        "endpoint_frozen": endpoint_frozen,
        "primary_endpoint": primary,
        "ready_for_phase3": ready_for_phase3,
        "endpoint_freeze_rule": FREEZE_RULE,
        "endpoint_QC": qc,
        "expression_markers_used_to_define_endpoint": False,
        "mapping_outputs_used_to_define_endpoint": False,
        "CytoSPACE_run": False,
        "SVTuner_run": False,
        "Stage4_run": False,
        "formal_metrics_recomputed": False,
        "allowed_claims": ALLOWED_CLAIMS,
        "disallowed_claims": DISALLOWED_CLAIMS,
        "output_dir": rel(OUT_DIR),
    }
    write_json(OUT_DIR / "bioapp_phase2c_endpoint_freeze_summary.json", summary)

    golden = {
        "stage_type": STAGE_TYPE,
        "biological_question_defined": True,
        "external_endpoint_predefined": True,
        "endpoint_independent_from_SVTuner": True,
        "endpoint_spatially_registered": bool(cta_mapped_spots > 0),
        "baseline_comparison_available": False,
        "endpoint_specific_improvement_defined": False,
        "endpoint_specific_quantitative_metric_available": False,
        "endpoint_baseline_SVTuner_spatial_comparison_available": False,
        "interpretation_boundary_defined": True,
        "biological_application_allowed": False,
        "allowed_claim_level": "endpoint freeze only",
        "decision": decision,
    }
    write_json(OUT_DIR / "bioapp_phase2c_golden_rules_v2_1_check.json", golden)

    if decision == "PASS":
        next_text = "BioApp Phase 3 - reference audit and formal baseline feasibility for CTA-defined endpoint"
    elif decision == "REVIEW_REQUIRED":
        next_text = "Review endpoint threshold, CTA spot composition, and spatial continuity before Phase 3.\nDo not run CytoSPACE/SVTuner/Stage4."
    else:
        next_text = "Stop. Fix input mapping or boundary violation before retrying Phase 2C."
    decision_txt = [
        "BioApp Phase 2C - endpoint freeze with validated CTA-to-spot mapping",
        "",
        f"Decision: {decision}",
        "",
        f"Endpoint frozen: {str(endpoint_frozen).lower()}",
        f"Primary endpoint: {primary if primary is not None else 'null'}",
        f"Ready for Phase 3: {str(ready_for_phase3).lower()}",
        "",
        "Allowed claim:",
        *[f"- {claim}" for claim in ALLOWED_CLAIMS],
        "",
        "Disallowed claims:",
        *[f"- {claim}" for claim in DISALLOWED_CLAIMS],
        "",
        "Next:",
        next_text,
        "",
    ]
    write_text(OUT_DIR / "decision.txt", "\n".join(decision_txt))
    write_readme(primary, decision)

    generated = list(OUT_DIR.glob("*"))
    extra_manifest_rows = [
        {
            "file": rel(path),
            "type": "input",
            "description": f"Phase 2B-R2 input used by Phase 2C: {name}",
            "created_by_phase": "BioApp Phase 2B-R2",
        }
        for name, path in inputs.items()
        if path.exists()
    ]
    write_manifest(generated + map_paths, extra_manifest_rows)

    print("BioApp Phase 2C completed.")
    print("\nDecision:")
    print(decision)
    print("\nEndpoint frozen:")
    print(str(endpoint_frozen).lower())
    print("\nPrimary endpoint:")
    print(primary if primary is not None else "null")
    print("\nST spots:")
    print(st_spots_total)
    print("\nCTA-mapped spots:")
    print(cta_mapped_spots)
    print("\nUnmatched spots:")
    print(unmatched_spots)
    print("\nEndpoint QC:")
    for cls in CLASSES:
        c = qc[cls]
        print(
            f"{cls}: positive={c['n_positive']}, negative={c['n_negative']}, "
            f"ambiguous={c['n_ambiguous']}, excluded={c['n_excluded']}, "
            f"eligible={str(c['eligible_endpoint']).lower()}"
        )
    print("\nBoundary checks:")
    print("CytoSPACE run: false")
    print("SVTuner run: false")
    print("Stage4 run: false")
    print("Expression markers used to define endpoint: false")
    print("Mapping outputs used to define endpoint: false")
    print("\nNext:")
    print(next_text)
    return 0 if decision in {"PASS", "REVIEW_REQUIRED"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
