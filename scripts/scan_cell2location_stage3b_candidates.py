#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

_HERE = Path(__file__).resolve()
_ROOT = _HERE.parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from scripts.prepare_cell2location_mouse_brain_stage3b_case import (  # noqa: E402
    read_10x_h5,
    target_markers_from_counts,
    unique_gene_positions,
)


SECTIONS = ["ST8059048", "ST8059049", "ST8059050", "ST8059051", "ST8059052"]

TARGETS = {
    "oligodendrocyte_opc": ["Oligo_1", "Oligo_2", "OPC_1", "OPC_2"],
    "oligodendrocyte": ["Oligo_1", "Oligo_2"],
    "opc": ["OPC_1", "OPC_2"],
    "thalamic_excitatory": ["Ext_Thal_1", "Ext_Thal_2"],
    "ext_thal_1": ["Ext_Thal_1"],
    "ext_thal_2": ["Ext_Thal_2"],
    "hippocampal_excitatory": [
        "Ext_Hpc_CA1",
        "Ext_Hpc_CA2",
        "Ext_Hpc_CA3",
        "Ext_Hpc_DG1",
        "Ext_Hpc_DG2",
    ],
    "astrocyte": [
        "Astro_AMY",
        "Astro_AMY_CTX",
        "Astro_CTX",
        "Astro_HPC",
        "Astro_HYPO",
        "Astro_STR",
        "Astro_THAL_hab",
        "Astro_THAL_lat",
        "Astro_THAL_med",
        "Astro_WM",
    ],
    "microglia": ["Micro"],
}

DEFAULT_CANDIDATES = [
    ("ST8059048", "oligodendrocyte_opc"),
    ("ST8059049", "oligodendrocyte_opc"),
    ("ST8059051", "oligodendrocyte_opc"),
    ("ST8059052", "oligodendrocyte_opc"),
    ("ST8059048", "oligodendrocyte"),
    ("ST8059049", "oligodendrocyte"),
    ("ST8059050", "oligodendrocyte"),
    ("ST8059051", "thalamic_excitatory"),
    ("ST8059048", "thalamic_excitatory"),
    ("ST8059049", "thalamic_excitatory"),
    ("ST8059052", "hippocampal_excitatory"),
    ("ST8059048", "hippocampal_excitatory"),
    ("ST8059049", "astrocyte"),
    ("ST8059050", "astrocyte"),
    ("ST8059052", "astrocyte"),
    ("ST8059050", "microglia"),
    ("ST8059052", "microglia"),
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Scan cell2location mouse brain section/target candidates for Stage3B."
    )
    p.add_argument("--project_root", default=".")
    p.add_argument("--out_dir", default="visualizations/cell2location_stage3b_case")
    p.add_argument("--quantile", type=float, default=0.80)
    p.add_argument("--n_markers", type=int, default=30)
    p.add_argument("--n_spatial_permutations", type=int, default=200)
    p.add_argument("--overwrite", action="store_true")
    p.add_argument(
        "--candidate",
        action="append",
        default=[],
        help="Candidate as SECTION:target_key. May be repeated.",
    )
    return p.parse_args()


def log_normalize_rows(values: sparse.spmatrix) -> sparse.csr_matrix:
    values = values.astype(np.float32).tocsr()
    totals = np.asarray(values.sum(axis=1)).ravel()
    totals[totals <= 0] = 1.0
    return values.multiply(10000.0 / totals[:, None]).log1p().tocsr()


def sample_name(section: str, target: str) -> str:
    clean = re.sub(r"[^a-z0-9]+", "_", target.lower()).strip("_")
    return f"cell2loc_scan_{section.lower()}_sc_missing_{clean}"


def run_command(cmd: list[str], cwd: Path) -> None:
    env = os.environ.copy()
    env["PYTHONPATH"] = str(cwd / "src") + os.pathsep + env.get("PYTHONPATH", "")
    ret = subprocess.run(cmd, cwd=cwd, env=env)
    if ret.returncode != 0:
        raise RuntimeError(f"Command failed ({ret.returncode}): {' '.join(cmd)}")


def load_marker_region(
    root: Path,
    section: str,
    drop_types: list[str],
    quantile: float,
    n_markers: int,
) -> tuple[pd.DataFrame, list[str]]:
    import anndata as ad

    source = root / "data" / "raw" / "cell2location_mouse_brain"
    sc_data = ad.read_h5ad(source / "sc.h5ad")
    if sc_data.raw is None:
        raise ValueError("sc.h5ad lacks raw counts")
    labels = sc_data.obs["annotation_1"].astype(str)
    raw_counts = sc_data.raw.X
    if not sparse.issparse(raw_counts):
        raw_counts = sparse.csr_matrix(raw_counts)

    section_dir = source / "mouse_brain_visium_wo_cloupe_data" / "rawdata" / section
    counts, barcodes, st_genes = read_10x_h5(section_dir / "filtered_feature_bc_matrix.h5")
    positions = pd.read_csv(
        section_dir / "spatial" / "tissue_positions_list.csv",
        header=None,
        names=["spot_id", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"],
    )
    positions["spot_id"] = positions["spot_id"].astype(str)
    barcode_pos = {barcode: i for i, barcode in enumerate(barcodes.astype(str))}
    tissue_spots = [
        spot
        for spot in positions.loc[positions["in_tissue"] == 1, "spot_id"].tolist()
        if spot in barcode_pos
    ]
    row_idx = np.array([barcode_pos[spot] for spot in tissue_spots], dtype=int)
    counts = counts[row_idx, :].tocsr()

    st_pos = unique_gene_positions(st_genes)
    markers = target_markers_from_counts(raw_counts, sc_data.raw.var, labels, drop_types, st_pos)
    markers = [gene for gene in markers[:n_markers] if gene in st_pos]
    if len(markers) < 8:
        raise ValueError(f"Too few marker genes for {drop_types}: {markers}")
    marker_idx = np.array([st_pos[gene] for gene in markers], dtype=int)
    marker_score = np.asarray(log_normalize_rows(counts)[:, marker_idx].mean(axis=1)).ravel()
    threshold = float(np.quantile(marker_score, quantile))

    marker = positions.set_index("spot_id").loc[tissue_spots].copy()
    marker["marker_score"] = marker_score
    marker["marker_region"] = marker_score >= threshold
    return marker, markers


def top_residual_genes(path: Path, selected_model: str, n: int = 8) -> str:
    if not path.exists() or not selected_model:
        return ""
    match = re.search(r"c(\d+)", selected_model)
    if not match:
        return ""
    comp_id = int(match.group(1))
    comp = f"component_{comp_id}"
    frame = pd.read_csv(path)
    if {"component", "gene", "abs_loading"}.issubset(frame.columns):
        order = frame.loc[frame["component"].astype(int) == comp_id].sort_values(
            "abs_loading", ascending=False
        )
        return ";".join(order["gene"].astype(str).head(n).tolist())
    if comp not in frame.columns or "gene" not in frame.columns:
        return ""
    order = frame.assign(abs_loading=frame[comp].abs()).sort_values("abs_loading", ascending=False)
    return ";".join(order["gene"].astype(str).head(n).tolist())


def evaluate_candidate(root: Path, section: str, target: str, args: argparse.Namespace) -> dict[str, object]:
    drop_types = TARGETS[target]
    sample = sample_name(section, target)

    stage3b_scores = (
        root
        / "data"
        / "processed"
        / "cell2location_mouse_brain"
        / sample
        / "stage3b_st_unsupported"
        / "spot_unsupported_scores.csv"
    )
    if args.overwrite or not stage3b_scores.exists():
        run_command(
            [
                sys.executable,
                "scripts/prepare_cell2location_mouse_brain_stage3b_case.py",
                "--project_root",
                ".",
                "--section",
                section,
                "--target_type",
                drop_types[0],
                "--drop_types",
                ",".join(drop_types),
                "--target_label",
                target,
                "--sample",
                sample,
                "--overwrite",
            ],
            root,
        )
        run_command(
            [
                sys.executable,
                "-m",
                "svtuner.cli",
                "stage3b",
                "--project-root",
                ".",
                "--sample",
                sample,
                "--fdr",
                "0.05",
                "--n-calibration",
                "0",
                "--n-spatial-permutations",
                str(args.n_spatial_permutations),
                "--random-seed",
                "42",
                "--sc-expr-source",
                "counts",
                "--sc-profile-source",
                "counts",
                "--expression-scale",
                "linear",
                "--sc-profile-scale",
                "linear",
                "--max-genes",
                "0",
            ],
            root,
        )

    marker, markers = load_marker_region(root, section, drop_types, args.quantile, args.n_markers)
    scores = pd.read_csv(stage3b_scores, index_col=0)
    scores.index = scores.index.astype(str)
    common = marker.index.intersection(scores.index)
    marker_region = marker.loc[common, "marker_region"].astype(bool)
    final_blank = scores.loc[common, "is_unsupported_region"].astype(bool)
    residual_blank = scores.loc[common, "is_residual_program_unsupported_region"].astype(bool)
    whole_blank = scores.loc[common, "is_whole_profile_unsupported_region"].astype(bool)

    marker_n = int(marker_region.sum())
    final_n = int(final_blank.sum())
    residual_n = int(residual_blank.sum())
    final_overlap = int((marker_region & final_blank).sum())
    residual_overlap = int((marker_region & residual_blank).sum())
    whole_overlap = int((marker_region & whole_blank).sum())
    selected_model = (
        scores["residual_program_selected_model"].dropna().astype(str).iloc[0]
        if "residual_program_selected_model" in scores.columns
        else ""
    )
    substitute = ""
    if "compensatory_substitute_type" in scores.columns and residual_overlap:
        sub = scores.loc[common[marker_region & residual_blank], "compensatory_substitute_type"]
        substitute = ";".join(f"{k}:{v}" for k, v in sub.value_counts().head(5).items())

    loading_path = stage3b_scores.parent / "residual_program_gene_loadings.csv"
    info_path = stage3b_scores.parent.parent / "stage1_preprocess" / "stage3b_case_info.json"
    info = json.loads(info_path.read_text(encoding="utf-8")) if info_path.exists() else {}

    result = {
        "section": section,
        "target": target,
        "drop_types": ",".join(drop_types),
        "sample": sample,
        "st_spots": int(len(common)),
        "target_reference_cells_removed": int(info.get("target_reference_cells_removed", -1)),
        "marker_region_spots": marker_n,
        "whole_blank_spots": int(whole_blank.sum()),
        "whole_marker_overlap": whole_overlap,
        "residual_blank_spots": residual_n,
        "residual_marker_overlap": residual_overlap,
        "residual_precision": residual_overlap / residual_n if residual_n else 0.0,
        "residual_recall": residual_overlap / marker_n if marker_n else 0.0,
        "final_blank_spots": final_n,
        "final_marker_overlap": final_overlap,
        "final_precision": final_overlap / final_n if final_n else 0.0,
        "final_recall": final_overlap / marker_n if marker_n else 0.0,
        "selected_residual_model": selected_model,
        "top_residual_genes": top_residual_genes(loading_path, selected_model),
        "top_compensatory_substitutes": substitute,
        "markers_used": ";".join(markers[:10]),
    }
    return result


def requested_candidates(items: list[str]) -> list[tuple[str, str]]:
    if not items:
        return DEFAULT_CANDIDATES
    parsed: list[tuple[str, str]] = []
    for item in items:
        if ":" not in item:
            raise ValueError(f"Candidate must be SECTION:target_key, got {item}")
        section, target = item.split(":", 1)
        section = section.strip()
        target = target.strip()
        if section not in SECTIONS:
            raise ValueError(f"Unknown section: {section}")
        if target not in TARGETS:
            raise ValueError(f"Unknown target: {target}")
        parsed.append((section, target))
    return parsed


def main() -> int:
    args = parse_args()
    root = Path(args.project_root).resolve()
    out_dir = root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for section, target in requested_candidates(args.candidate):
        print(f"[scan] {section} {target}", flush=True)
        try:
            row = evaluate_candidate(root, section, target, args)
            row["status"] = "ok"
        except Exception as exc:  # noqa: BLE001
            row = {
                "section": section,
                "target": target,
                "drop_types": ",".join(TARGETS.get(target, [])),
                "sample": sample_name(section, target),
                "status": f"failed: {exc}",
            }
        rows.append(row)
        print(json.dumps(row, ensure_ascii=False), flush=True)

    frame = pd.DataFrame(rows)
    out_path = out_dir / "cell2location_stage3b_candidate_scan.csv"
    frame.to_csv(out_path, index=False)
    print(f"[scan] wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
