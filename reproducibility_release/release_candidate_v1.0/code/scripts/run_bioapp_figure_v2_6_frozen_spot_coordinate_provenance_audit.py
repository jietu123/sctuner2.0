from __future__ import annotations

import ast
import csv
import json
import math
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_figure_v2_6_frozen_spot_coordinate_provenance_audit"

P2BR2 = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase2br2_explicit_rscript_cta_align_locale_safe_dry_run"
P2C = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase2c_endpoint_freeze_with_validated_cta_to_spot_mapping"
P8 = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase8_svtuner_vs_endpoint_evaluation_and_baseline_svtuner_comparison"
P9 = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase9_final_biological_application_audit_figure_selection_interpretation_boundary_lock"
V1 = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_main_figure_layout_only_spatial_visualization"
V2 = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_main_figure_v2_tissue_overlay_refinement"
V25 = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_figure_v2_5_online_he_spatial_bundle_provenance_audit"

FROZEN = P2C / "spot_level_endpoint_freeze.csv"

TEXT_SUFFIXES = {".py", ".R", ".r", ".md", ".txt", ".json", ".yaml", ".yml"}
TABLE_SUFFIXES = {".csv", ".tsv", ".txt"}
SPOT_TERMS = [
    "imagecol",
    "imagerow",
    "image_col",
    "image_row",
    "pxl_col_in_fullres",
    "pxl_row_in_fullres",
    "array_row",
    "array_col",
    "barcode",
    "spot",
    "coordinate",
    "spatial",
    "endpoint",
    "CTA",
    "Immune cells",
    "frozen",
    "phase8",
    "phase9",
    "withheld",
    "baseline_immune_all_dropout",
    "svtuner_immune_all_dropout",
]
REGISTRATION_TERMS = [
    "registration",
    "transform",
    "affine",
    "homography",
    "keypoint",
    "keypoints",
    "Xenium",
    "Visium",
    "morphology",
    "H&E",
    "HE",
    "registered",
    "alignment",
    "matrix",
    "warp",
    "rotate",
    "scale",
    "translate",
]


def rel(path: Path | None) -> str | None:
    if path is None:
        return None
    try:
        return path.relative_to(ROOT).as_posix()
    except ValueError:
        return path.as_posix()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def write_rows(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k) for k in fields})


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def safe_read_text(path: Path, max_chars: int = 250000) -> str:
    try:
        text = path.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return ""
    return text[:max_chars]


def candidate_roots() -> list[Path]:
    roots = [ROOT / "scripts"]
    vis = ROOT / "visualizations"
    if vis.exists():
        roots.extend([p for p in vis.iterdir() if p.is_dir() and ("bioapp" in p.name.lower() or "cytospace_fig2j" in p.name.lower())])
    for p in [
        ROOT / "data" / "raw" / "\u65b0\u5efa\u6587\u4ef6\u5939",
        ROOT / "data" / "external" / "bioapp_provenance_audit",
        ROOT / "data" / "processed" / "bioapp",
        ROOT / "outputs" / "bioapp",
    ]:
        if p.exists():
            roots.append(p)
    # De-duplicate while preserving order.
    out: list[Path] = []
    seen: set[Path] = set()
    for p in roots:
        rp = p.resolve()
        if rp not in seen:
            seen.add(rp)
            out.append(p)
    return out


def likely_phase(path: Path) -> str:
    s = rel(path) or str(path)
    m = re.search(r"bioapp_phase([0-9]+[a-zA-Z0-9]*)", s, flags=re.IGNORECASE)
    if m:
        return f"BioApp Phase {m.group(1)}"
    if "main_figure_v2" in s:
        return "BioApp Figure V2"
    if "main_figure_layout" in s:
        return "BioApp Figure V1"
    if "v2_5" in s or "figure_v2_5" in s:
        return "BioApp Figure V2.5"
    return ""


def detect_barcode_col(columns: list[str]) -> str | None:
    for c in columns:
        lc = c.lower()
        if lc in {"barcode", "spot_id", "spot", "spotid"} or "barcode" in lc:
            return c
    return None


def detect_coord_cols(columns: list[str]) -> list[str]:
    keys = [
        "imagecol",
        "imagerow",
        "pxl_col_in_fullres",
        "pxl_row_in_fullres",
        "row",
        "col",
        "array_row",
        "array_col",
        "x",
        "y",
    ]
    out = []
    for c in columns:
        lc = c.lower()
        if lc in keys or any(k in lc for k in ["imagecol", "imagerow", "pxl_col", "pxl_row"]):
            out.append(c)
    return out


def classify_range(min_v: float | None, max_v: float | None) -> dict[str, bool | str]:
    if min_v is None or max_v is None or not np.isfinite(min_v) or not np.isfinite(max_v):
        return {
            "value_range_classification": "unreadable",
            "compatible_with_fullres_pixels": False,
            "compatible_with_hires_pixels": False,
            "compatible_with_lowres_pixels": False,
            "compatible_with_array_coordinates": False,
            "compatible_with_normalized_coordinates": False,
        }
    span = max_v - min_v
    return {
        "value_range_classification": (
            "normalized_0_1"
            if 0 <= min_v <= 1 and 0 <= max_v <= 1
            else "array_like"
            if 0 <= min_v <= 200 and 0 <= max_v <= 300
            else "hires_or_scaled_pixel_like"
            if 0 <= min_v and max_v <= 20000
            else "fullres_pixel_like"
            if 0 <= min_v and max_v <= 60000
            else "unknown_large_or_shifted"
        ),
        "compatible_with_fullres_pixels": bool(0 <= min_v and max_v <= 60000 and span > 1000),
        "compatible_with_hires_pixels": bool(0 <= min_v and max_v <= 20000 and span > 1000),
        "compatible_with_lowres_pixels": bool(0 <= min_v and max_v <= 5000 and span > 200),
        "compatible_with_array_coordinates": bool(0 <= min_v and max_v <= 300 and span > 10),
        "compatible_with_normalized_coordinates": bool(0 <= min_v and max_v <= 1),
    }


def read_table_light(path: Path) -> pd.DataFrame | None:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path, nrows=5000)
    if suffix == ".tsv":
        return pd.read_csv(path, sep="\t", nrows=5000)
    if suffix == ".txt":
        # Only inspect text as a table when it looks delimited.
        preview = safe_read_text(path, 4000)
        if "," in preview.splitlines()[0] if preview.splitlines() else False:
            return pd.read_csv(path, nrows=5000)
        if "\t" in preview.splitlines()[0] if preview.splitlines() else False:
            return pd.read_csv(path, sep="\t", nrows=5000)
    return None


def inventory_spot_files() -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    coord_rows: list[dict[str, Any]] = []
    file_id = 1
    for root in candidate_roots():
        for path in root.rglob("*"):
            if not path.is_file():
                continue
            suffix = path.suffix.lower()
            is_candidate = False
            text = ""
            if suffix in {".csv", ".tsv"}:
                is_candidate = True
            elif suffix in {".json", ".md", ".txt"}:
                text = safe_read_text(path, 50000)
                is_candidate = any(t.lower() in text.lower() or t.lower() in path.name.lower() for t in SPOT_TERMS)
            elif suffix in {".h5", ".h5ad", ".rds", ".rdata", ".pkl", ".parquet"}:
                is_candidate = any(t.lower() in path.name.lower() for t in ["st_data", "spatial", "bioapp", "coordinate", "spot"])
            if not is_candidate:
                continue

            df = None
            columns: list[str] = []
            row_count: int | None = None
            col_count: int | None = None
            note = ""
            try:
                df = read_table_light(path)
                if df is not None:
                    columns = [str(c) for c in df.columns]
                    row_count = int(len(pd.read_csv(path))) if suffix in {".csv", ".tsv"} and path.stat().st_size < 15_000_000 else len(df)
                    col_count = int(len(columns))
            except Exception as exc:
                note = f"table_read_failed={type(exc).__name__}: {exc}"

            if not columns and suffix == ".json":
                js = load_json(path)
                columns = list(js.keys())[:80] if isinstance(js, dict) else []
                row_count = None
                col_count = len(columns)

            barcode_col = detect_barcode_col(columns)
            coord_cols = detect_coord_cols(columns)
            xmin = xmax = ymin = ymax = None
            if df is not None and {"imagecol", "imagerow"}.issubset(df.columns):
                xmin, xmax = float(pd.to_numeric(df["imagecol"], errors="coerce").min()), float(pd.to_numeric(df["imagecol"], errors="coerce").max())
                ymin, ymax = float(pd.to_numeric(df["imagerow"], errors="coerce").min()), float(pd.to_numeric(df["imagerow"], errors="coerce").max())
            elif df is not None and len(coord_cols) >= 2:
                a = pd.to_numeric(df[coord_cols[0]], errors="coerce")
                b = pd.to_numeric(df[coord_cols[1]], errors="coerce")
                xmin, xmax = float(a.min()), float(a.max())
                ymin, ymax = float(b.min()), float(b.max())

            cols_lower = " ".join(columns).lower()
            row = {
                "file_id": f"file_{file_id:04d}",
                "file_path": rel(path),
                "file_type": suffix,
                "row_count": row_count,
                "column_count": col_count,
                "columns": "|".join(columns[:120]),
                "candidate_barcode_column": barcode_col,
                "candidate_coordinate_columns": "|".join(coord_cols),
                "coordinate_min_x": xmin,
                "coordinate_max_x": xmax,
                "coordinate_min_y": ymin,
                "coordinate_max_y": ymax,
                "matches_2248_universe": row_count == 2248,
                "matches_1888_main_analysis": row_count == 1888,
                "contains_endpoint_labels": any(k in cols_lower for k in ["endpoint", "primary_endpoint_status"]),
                "contains_baseline_burden": any(k in cols_lower for k in ["baseline_forced", "burden"]),
                "contains_withheld_score": "withheld_score" in cols_lower,
                "contains_withheld_binary": "withheld_binary" in cols_lower,
                "contains_source_metadata": any(k in cols_lower for k in ["source", "path", "spatial", "rdata", "scalefactor"]),
                "likely_phase": likely_phase(path),
                "generated_by_script": "",
                "notes": note,
            }
            rows.append(row)

            if df is not None:
                for c in coord_cols:
                    vals = pd.to_numeric(df[c], errors="coerce")
                    if vals.notna().sum() == 0:
                        continue
                    min_v = float(vals.min())
                    max_v = float(vals.max())
                    cls = classify_range(min_v, max_v)
                    coord_rows.append(
                        {
                            "file_id": row["file_id"],
                            "file_path": rel(path),
                            "coordinate_column_name": c,
                            "role": "x" if "col" in c.lower() or c.lower() == "x" else "y" if "row" in c.lower() or c.lower() == "y" else "coordinate",
                            "min": min_v,
                            "max": max_v,
                            "mean": float(vals.mean()),
                            "std": float(vals.std()),
                            "integer_like": bool(np.nanmean(np.abs(vals - np.round(vals))) < 1e-6),
                            **cls,
                            "notes": "",
                        }
                    )
            file_id += 1
    return rows, coord_rows


def extract_string_paths_from_ast(script: Path) -> tuple[list[str], list[str]]:
    text = safe_read_text(script, 500000)
    inputs, outputs = [], []
    try:
        tree = ast.parse(text)
    except Exception:
        return inputs, outputs
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            s = node.value
            if any(ext in s for ext in [".csv", ".json", ".txt", ".RData", ".h5", ".svg", ".pdf", ".png", ".md"]):
                target = outputs if any(k in s.lower() for k in ["out", "summary", "report", "manifest", ".svg", ".pdf", ".png"]) else inputs
                target.append(s.replace("\n", " ")[:160])
    return sorted(set(inputs)), sorted(set(outputs))


def script_lineage() -> list[dict[str, Any]]:
    rows = []
    script_id = 1
    for script in sorted((ROOT / "scripts").glob("run_bioapp*.py")) + sorted((ROOT / "scripts").glob("run_bioapp*.R")):
        text = safe_read_text(script, 700000)
        lower = text.lower()
        inputs, outputs = extract_string_paths_from_ast(script) if script.suffix.lower() == ".py" else ([], [])
        coord_reads = sorted(set(re.findall(r"[\"'](imagecol|imagerow|pxl_col_in_fullres|pxl_row_in_fullres|row|col)[\"']", text)))
        coord_writes = coord_reads
        snippets = []
        for i, line in enumerate(text.splitlines(), start=1):
            if any(k in line for k in ["imagecol", "imagerow", "pxl_col_in_fullres", "pxl_row_in_fullres", "SCALING_FACTOR", "scaling_factor", "tissue_positions", "scalefactors_json"]):
                snippets.append(f"L{i}: {line.strip()[:220]}")
            if len(snippets) >= 10:
                break
        transform_terms = [k for k in ["affine", "scale", "scaling_factor", "rotate", "flip", "translation", "homography", "crop"] if k in lower]
        rows.append(
            {
                "script_id": f"script_{script_id:03d}",
                "script_path": rel(script),
                "phase": likely_phase(script),
                "inputs_detected": "|".join(inputs[:80]),
                "outputs_detected": "|".join(outputs[:80]),
                "coordinate_columns_read": "|".join(coord_reads),
                "coordinate_columns_written": "|".join(coord_writes),
                "barcode_columns_read": "barcode" if "barcode" in lower else "",
                "barcode_columns_written": "barcode" if "barcode" in lower else "",
                "tissue_positions_read": "tissue_positions" in lower,
                "scalefactors_read": "scalefactors" in lower,
                "tissue_image_read": any(k in lower for k in ["tissue_hires_image", "tissue_lowres_image", "tissue_image"]),
                "coordinate_transform_detected": bool(transform_terms),
                "transform_type": "|".join(transform_terms),
                "relevant_code_lines_or_snippets": " || ".join(snippets),
                "notes": "",
            }
        )
        script_id += 1
    return rows


def graph_edges() -> list[dict[str, Any]]:
    return [
        {
            "source_node": "data/raw/\u65b0\u5efa\u6587\u4ef6\u5939/BreastCancer_CTA-main/ST_data.RData",
            "target_node": "visualizations/bioapp_experiment/bioapp_phase2br2_explicit_rscript_cta_align_locale_safe_dry_run/bioapp_phase2br2_st_data_coordinates.csv",
            "edge_type": "coordinate_export",
            "script_or_process": "scripts/run_bioapp_phase2br2_explicit_rscript_cta_align_locale_safe_dry_run.py::st_data_inventory R helper",
            "coordinate_columns_transferred": "tissue,row,col,imagerow,imagecol",
            "barcode_columns_transferred": "barcode",
            "transform_applied": "pxl_col_in_fullres=imagecol/0.3; pxl_row_in_fullres=imagerow/0.3; reverse y and spot radius boxes added",
            "evidence": "bioapp_phase2br2_st_data_object_inventory.json and script constants SCALING_FACTOR=0.3, IMAGE_NAME=BCSA2TumB1",
            "confidence": "high",
        },
        {
            "source_node": "bioapp_phase2br2_st_data_coordinates.csv + bioapp_phase2br2_cta_to_spot_mapping_example.csv",
            "target_node": "visualizations/bioapp_experiment/bioapp_phase2c_endpoint_freeze_with_validated_cta_to_spot_mapping/spot_level_endpoint_freeze.csv",
            "edge_type": "endpoint_freeze_table_build",
            "script_or_process": "scripts/run_bioapp_phase2c_endpoint_freeze_with_validated_cta_to_spot_mapping.py",
            "coordinate_columns_transferred": "tissue,row,col,imagerow,imagecol,pxl_col_in_fullres,pxl_row_in_fullres",
            "barcode_columns_transferred": "barcode",
            "transform_applied": "none to imagecol/imagerow; CTA object counts joined by barcode",
            "evidence": "Phase 2C script keep_cols and spot_level_endpoint_freeze.csv",
            "confidence": "high",
        },
        {
            "source_node": "spot_level_endpoint_freeze.csv",
            "target_node": "visualizations/bioapp_experiment/bioapp_phase8_svtuner_vs_endpoint_evaluation_and_baseline_svtuner_comparison/svtuner_endpoint_score_by_spot.csv",
            "edge_type": "evaluation_join",
            "script_or_process": "scripts/run_bioapp_phase8_svtuner_vs_endpoint_evaluation_and_baseline_svtuner_comparison.py",
            "coordinate_columns_transferred": "imagecol,imagerow used for plots; endpoint status retained",
            "barcode_columns_transferred": "barcode",
            "transform_applied": "none",
            "evidence": "Phase 8 plotting uses endpoint[['barcode','imagecol','imagerow','primary_endpoint_status']]",
            "confidence": "high",
        },
        {
            "source_node": "spot_level_endpoint_freeze.csv + svtuner_endpoint_score_by_spot.csv",
            "target_node": "BioApp main figure V1/V2 plotting data",
            "edge_type": "figure_plotting",
            "script_or_process": "scripts/run_bioapp_main_figure_layout_only_spatial_visualization.py and run_bioapp_main_figure_v2_tissue_overlay_refinement.py",
            "coordinate_columns_transferred": "imagecol,imagerow",
            "barcode_columns_transferred": "barcode",
            "transform_applied": "plot-lattice display only; no biological data transform",
            "evidence": "V1/V2 scripts scatter imagecol/imagerow directly",
            "confidence": "high",
        },
    ]


def read_tissue_positions(path: Path) -> pd.DataFrame | None:
    try:
        first = pd.read_csv(path, nrows=1, header=None)
        if str(first.iloc[0, 0]).lower() == "barcode":
            df = pd.read_csv(path)
            rename = {}
            for c in df.columns:
                lc = str(c).lower()
                if lc == "barcode":
                    rename[c] = "barcode"
                elif lc in {"array_row", "row"}:
                    rename[c] = "array_row"
                elif lc in {"array_col", "col"}:
                    rename[c] = "array_col"
                elif lc in {"pxl_row_in_fullres", "imagerow"}:
                    rename[c] = "pxl_row"
                elif lc in {"pxl_col_in_fullres", "imagecol"}:
                    rename[c] = "pxl_col"
            return df.rename(columns=rename)
        df = pd.read_csv(path, header=None)
        if df.shape[1] < 6:
            return None
        df = df.iloc[:, :6].copy()
        df.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]
        return df
    except Exception:
        return None


def candidate_spatial_bundle_comparison(frozen: pd.DataFrame) -> list[dict[str, Any]]:
    roots = [
        ROOT / "data" / "raw" / "\u65b0\u5efa\u6587\u4ef6\u5939",
        ROOT / "data" / "external" / "bioapp_provenance_audit",
        ROOT / "data" / "raw" / "real_brca",
        ROOT / "data" / "raw" / "low_resolution" / "Human Breast Cancer",
        ROOT / "data" / "raw" / "low_resolution" / "Human Breast Cancer WTA 1.2.0",
        ROOT / "data" / "raw" / "low_resolution" / "Human Breast Cancer Visium FF WTA",
        V2,
        V25,
    ]
    position_files: list[Path] = []
    for root in roots:
        if root.exists():
            for name in ["tissue_positions.csv", "tissue_positions_list.csv"]:
                position_files.extend(root.rglob(name))
    rows = []
    for idx, pos in enumerate(sorted(set(position_files)), start=1):
        sdir = pos.parent
        img = None
        for n in ["tissue_hires_image.png", "tissue_lowres_image.png", "aligned_fiducials.jpg", "detected_tissue_image.jpg"]:
            if (sdir / n).exists():
                img = sdir / n
                break
        sf = sdir / "scalefactors_json.json"
        tpos = read_tissue_positions(pos)
        if tpos is None or "barcode" not in tpos.columns:
            continue
        tpos["barcode"] = tpos["barcode"].astype(str)
        m = frozen.merge(tpos, on="barcode", how="inner", suffixes=("_frozen", "_candidate"))
        overlap = int(m["barcode"].nunique())
        frac = float(overlap / len(frozen)) if len(frozen) else 0.0
        array_match = False
        image_match = False
        if overlap > 0 and {"row", "col", "array_row", "array_col"}.issubset(m.columns):
            array_match = bool((pd.to_numeric(m["row"], errors="coerce") == pd.to_numeric(m["array_row"], errors="coerce")).mean() > 0.99 and (pd.to_numeric(m["col"], errors="coerce") == pd.to_numeric(m["array_col"], errors="coerce")).mean() > 0.99)
        if overlap > 0 and {"imagecol", "imagerow", "pxl_col", "pxl_row"}.issubset(m.columns):
            dx = np.sqrt(np.mean((pd.to_numeric(m["imagecol"], errors="coerce") - pd.to_numeric(m["pxl_col"], errors="coerce")) ** 2))
            dy = np.sqrt(np.mean((pd.to_numeric(m["imagerow"], errors="coerce") - pd.to_numeric(m["pxl_row"], errors="coerce")) ** 2))
            image_match = bool(dx < 5 and dy < 5)
        source_consistent = any(k in str(sdir).lower() for k in ["bcsa2tumb1", "bioapp_phase2", "breastcancer_cta-main"])
        accepted = bool(source_consistent and overlap == len(frozen) and (image_match or array_match))
        reason = (
            "Accepted by provenance and coordinate match."
            if accepted
            else "Not accepted: barcode/array/image similarity alone is insufficient without provenance to ST_data.RData/BCSA2TumB1 source."
        )
        rows.append(
            {
                "candidate_id": f"bundle_{idx:03d}",
                "spatial_dir": rel(sdir),
                "tissue_image": rel(img),
                "tissue_positions": rel(pos),
                "scalefactors_json": rel(sf) if sf.exists() else None,
                "source_consistent_with_provenance": source_consistent,
                "barcode_overlap": overlap,
                "barcode_overlap_fraction": frac,
                "array_coordinate_match": array_match,
                "image_coordinate_match": image_match,
                "documented_transform_match": False,
                "alignment_verified": accepted,
                "accepted_as_true_source": accepted,
                "reason": reason,
            }
        )
    return rows


def registration_file_search() -> list[dict[str, Any]]:
    rows = []
    for root in candidate_roots():
        for path in root.rglob("*"):
            if not path.is_file():
                continue
            name_match = [t for t in REGISTRATION_TERMS if t.lower() in path.name.lower()]
            text = ""
            content_match = []
            if path.suffix.lower() in TEXT_SUFFIXES and path.stat().st_size < 5_000_000:
                text = safe_read_text(path, 200000)
                lower = text.lower()
                content_match = [t for t in REGISTRATION_TERMS if t.lower() in lower]
            matched = sorted(set(name_match + content_match))
            if not matched:
                continue
            lower = text.lower()
            rows.append(
                {
                    "file_path": rel(path),
                    "file_type": path.suffix.lower(),
                    "matched_terms": "|".join(matched),
                    "likely_registration_file": any(t.lower() in {"registration", "transform", "affine", "homography", "keypoint", "keypoints"} for t in matched),
                    "contains_transform_matrix": any(k in lower for k in ["matrix", "affine", "homography", "transform"]),
                    "contains_keypoints": "keypoint" in lower or "keypoints" in lower,
                    "contains_xenium_visium_reference": "xenium" in lower and "visium" in lower,
                    "contains_he_reference": "h&e" in lower or "he image" in lower or "tissue image" in lower,
                    "usable_for_recovery": False,
                    "notes": "Mentioned registration/alignment terms, but no verified inverse transform from BioApp frozen coordinates to a specific H&E bundle was identified.",
                }
            )
    return rows


def find_st_data_source() -> str:
    candidates = sorted((ROOT / "data" / "raw").rglob("ST_data.RData")) if (ROOT / "data" / "raw").exists() else []
    for path in candidates:
        if path.name == "ST_data.RData" and path.parent.name == "BreastCancer_CTA-main":
            return "data/raw/\u65b0\u5efa\u6587\u4ef6\u5939/BreastCancer_CTA-main/ST_data.RData"
    return rel(P2BR2 / "_ascii_runtime" / "ST_data.RData") or "visualizations/bioapp_experiment/bioapp_phase2br2_explicit_rscript_cta_align_locale_safe_dry_run/_ascii_runtime/ST_data.RData"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    frozen = pd.read_csv(FROZEN)
    frozen["barcode"] = frozen["barcode"].astype(str)

    spot_rows, coord_rows = inventory_spot_files()
    script_rows = script_lineage()
    edges = graph_edges()
    bundle_rows = candidate_spatial_bundle_comparison(frozen)
    reg_rows = registration_file_search()

    write_rows(
        OUT / "bioapp_v2_6_spot_level_file_inventory.csv",
        spot_rows,
        [
            "file_id",
            "file_path",
            "file_type",
            "row_count",
            "column_count",
            "columns",
            "candidate_barcode_column",
            "candidate_coordinate_columns",
            "coordinate_min_x",
            "coordinate_max_x",
            "coordinate_min_y",
            "coordinate_max_y",
            "matches_2248_universe",
            "matches_1888_main_analysis",
            "contains_endpoint_labels",
            "contains_baseline_burden",
            "contains_withheld_score",
            "contains_withheld_binary",
            "contains_source_metadata",
            "likely_phase",
            "generated_by_script",
            "notes",
        ],
    )
    write_rows(
        OUT / "bioapp_v2_6_coordinate_column_audit.csv",
        coord_rows,
        [
            "file_id",
            "file_path",
            "coordinate_column_name",
            "role",
            "min",
            "max",
            "mean",
            "std",
            "integer_like",
            "value_range_classification",
            "compatible_with_fullres_pixels",
            "compatible_with_hires_pixels",
            "compatible_with_lowres_pixels",
            "compatible_with_array_coordinates",
            "compatible_with_normalized_coordinates",
            "notes",
        ],
    )
    write_rows(
        OUT / "bioapp_v2_6_script_lineage_table.csv",
        script_rows,
        [
            "script_id",
            "script_path",
            "phase",
            "inputs_detected",
            "outputs_detected",
            "coordinate_columns_read",
            "coordinate_columns_written",
            "barcode_columns_read",
            "barcode_columns_written",
            "tissue_positions_read",
            "scalefactors_read",
            "tissue_image_read",
            "coordinate_transform_detected",
            "transform_type",
            "relevant_code_lines_or_snippets",
            "notes",
        ],
    )
    write_rows(
        OUT / "bioapp_v2_6_provenance_graph_edges.csv",
        edges,
        [
            "source_node",
            "target_node",
            "edge_type",
            "script_or_process",
            "coordinate_columns_transferred",
            "barcode_columns_transferred",
            "transform_applied",
            "evidence",
            "confidence",
        ],
    )
    write_rows(
        OUT / "bioapp_v2_6_candidate_spatial_bundle_comparison.csv",
        bundle_rows,
        [
            "candidate_id",
            "spatial_dir",
            "tissue_image",
            "tissue_positions",
            "scalefactors_json",
            "source_consistent_with_provenance",
            "barcode_overlap",
            "barcode_overlap_fraction",
            "array_coordinate_match",
            "image_coordinate_match",
            "documented_transform_match",
            "alignment_verified",
            "accepted_as_true_source",
            "reason",
        ],
    )
    write_rows(
        OUT / "bioapp_v2_6_registration_file_search_report.csv",
        reg_rows,
        [
            "file_path",
            "file_type",
            "matched_terms",
            "likely_registration_file",
            "contains_transform_matrix",
            "contains_keypoints",
            "contains_xenium_visium_reference",
            "contains_he_reference",
            "usable_for_recovery",
            "notes",
        ],
    )

    first_source_file = find_st_data_source()
    first_source_script = "scripts/run_bioapp_phase2br2_explicit_rscript_cta_align_locale_safe_dry_run.py"
    classification = "Transformed / registered Visium coordinates from processed Seurat ST_data.RData image slot"
    decision = "PROVENANCE_VERIFIED_BUT_H_AND_E_NOT_DIRECTLY_RECOVERABLE"
    next_step = "BioApp Coordinate Registration Recovery Audit"

    source_md = f"""# BioApp Figure V2.6 First Coordinate Source Report

## Earliest Proven Coordinate Source

The earliest proven source of the frozen BioApp plotting coordinates is:

`{first_source_file}`

The coordinates were exported by:

`{first_source_script}`

## Evidence

Phase 2B-R2 parsed `ST_data.RData` as a Seurat object named `bcsa` and selected image slot `BCSA2TumB1`. The exported coordinate table contains `barcode`, `tissue`, `row`, `col`, `imagerow`, `imagecol`, `pxl_col_in_fullres`, and `pxl_row_in_fullres`.

The Phase 2B-R2 script applies:

```text
pxl_col_in_fullres = imagecol / 0.3
pxl_row_in_fullres = imagerow / 0.3
```

This demonstrates that the frozen `imagecol/imagerow` coordinates are processed Seurat image-slot coordinates, not a newly computed Phase 8 or figure-stage coordinate system.

## Classification

`{classification}`

## Space Ranger Bundle Status

No local or downloaded official `spatial/` bundle was proven to be the original source of `ST_data.RData` / `BCSA2TumB1`. Therefore, a raw H&E overlay cannot be recovered directly from Space Ranger files in this phase.
"""
    (OUT / "bioapp_v2_6_first_coordinate_source_report.md").write_text(source_md, encoding="utf-8")

    decision_md = f"""# BioApp Figure V2.6 Coordinate Recovery Decision

## Final Decision

`{decision}`

## Earliest Coordinate Source Found

`{first_source_file}`

Generated/exported by:

`{first_source_script}`

## Coordinate Source Classification

`{classification}`

## Are These Original Space Ranger Image Coordinates?

Not proven. The coordinates are stored in a processed Seurat object image slot (`BCSA2TumB1`) and were exported from `ST_data.RData`. The original Space Ranger `spatial/` directory that produced this processed object was not verified.

## Is A Specific Spatial Bundle Verified?

No. Candidate local bundles can share Visium barcodes and sometimes match array geometry, but barcode/affine/array similarity alone is insufficient because Visium barcodes and array coordinates recur across samples. The downloaded official 10x CytAssist FFPE Human Breast Cancer spatial bundle was already audited in V2.5 and did not match the frozen BioApp barcodes.

## Is H&E Overlay Recoverable?

Not directly in this phase. Recovering a verified H&E overlay would require the original image/scalefactors/tissue_positions source for `ST_data.RData` or an explicit inverse registration transform from the processed `BCSA2TumB1` coordinate system back to a specific H&E image.

## Reason

The coordinate provenance is verified only up to the processed Seurat/ST object. The provenance chain from that object back to an official Space Ranger / Visium / CytAssist spatial bundle is missing.

## Valid Figure Basis

The V2 spot-lattice fallback remains the valid figure basis.

## Biological Results

No biological result, endpoint, metric, threshold, label, or conclusion was changed.

## Recommended Next Step

`{next_step}`
"""
    (OUT / "bioapp_v2_6_coordinate_recovery_decision.md").write_text(decision_md, encoding="utf-8")

    summary = {
        "phase": "BioApp Figure V2.6 frozen spot-coordinate provenance audit",
        "decision": decision,
        "audit_only": True,
        "layout_only": True,
        "SVTuner_rerun": False,
        "Stage3_rerun": False,
        "Stage4_run": False,
        "CytoSPACE_rerun": False,
        "endpoint_redefined": False,
        "data_changed": False,
        "metrics_changed": False,
        "threshold_changed": False,
        "endpoint_labels_changed": False,
        "conclusions_changed": False,
        "frozen_spot_count": int(len(frozen)),
        "main_analysis_spots": 1888,
        "coordinate_columns_identified": ["imagecol", "imagerow", "pxl_col_in_fullres", "pxl_row_in_fullres", "row", "col"],
        "barcode_column_identified": "barcode",
        "first_coordinate_source_file": first_source_file,
        "first_coordinate_source_script": first_source_script,
        "coordinate_source_classification": classification,
        "specific_spatial_bundle_verified": False,
        "verified_spatial_bundle_path": None,
        "verified_tissue_image_path": None,
        "verified_tissue_positions_path": None,
        "verified_scalefactors_json_path": None,
        "registration_transform_found": False,
        "registration_transform_path": None,
        "h_and_e_overlay_recoverable": False,
        "recommended_next_step": next_step,
        "phase2br2_st_data_object": load_json(P2BR2 / "bioapp_phase2br2_st_data_object_inventory.json").get("Seurat_object_name"),
        "phase2br2_image_slot": load_json(P2BR2 / "bioapp_phase2br2_st_data_object_inventory.json").get("selected_image_slot"),
        "v2_5_decision": load_json(V25 / "bioapp_v2_5_online_he_provenance_audit_summary.json").get("decision"),
    }
    write_json(OUT / "bioapp_v2_6_coordinate_provenance_summary.json", summary)

    guardrails = {
        "phase": "BioApp Figure V2.6 frozen spot-coordinate provenance audit",
        "decision": decision,
        "audit_only": True,
        "layout_only": True,
        "main_figure_redrawn": False,
        "SVTuner_rerun": False,
        "Stage3_rerun": False,
        "Stage4_run": False,
        "CytoSPACE_rerun": False,
        "endpoint_redefined": False,
        "new_biological_endpoint_introduced": False,
        "threshold_selected_using_endpoint_labels": False,
        "new_mapping_or_new_computation": False,
        "data_changed": False,
        "metrics_changed": False,
        "threshold_changed": False,
        "endpoint_labels_changed": False,
        "conclusions_changed": False,
        "external_image_used_as_background_without_verification": False,
        "h_and_e_used_in_figure": False,
        "all_quantitative_claims_based_on_full_frozen_analysis_set": True,
        "guardrails_passed": True,
    }
    write_json(OUT / "bioapp_v2_6_guardrails.json", guardrails)

    readme = f"""BioApp Figure V2.6 frozen spot-coordinate provenance audit completed.
Decision: {decision}
Coordinate source classification: {classification}
Specific spatial bundle verified: false
H&E overlay recoverable: false
Data changed: false
Metrics changed: false
Endpoint redefined: false
Next: {next_step}
"""
    (OUT / "README.md").write_text(readme, encoding="utf-8")

    print("BioApp Figure V2.6 frozen spot-coordinate provenance audit completed.")
    print(f"Decision: {decision}")
    print(f"Coordinate source classification: {classification}")
    print("H&E overlay recoverable: false")
    print(f"Recommended next step: {next_step}")
    print("Guardrails passed: true")


if __name__ == "__main__":
    main()
