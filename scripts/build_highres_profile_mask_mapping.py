from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from matplotlib import colormaps
from matplotlib.colors import LinearSegmentedColormap, Normalize
from PIL import Image, ImageDraw, ImageFont


DATASETS = {
    "HumanBreastCancerPatient1": "Monocytes and Macrophages",
    "HumanColonCancerPatient1": "Fibroblasts",
    "HumanLungCancerPatient1": "Plasma cells",
    "HumanMelanomaPatient1": "Fibroblasts",
    "HumanMelanomaPatient2": "B cells",
}

MARKERS = {
    "T cells": ["CD3D", "CD3E", "CD3G", "TRAC", "CD2", "CD4", "CD8A", "CD8B"],
    "B cells": ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1"],
    "Plasma cells": ["MZB1", "JCHAIN", "IGKC", "IGHG1", "TNFRSF17", "SDC1"],
    "Monocytes and Macrophages": ["CD68", "CD163", "LYZ", "LST1", "FCGR3A", "C1QA", "C1QB", "C1QC", "CSF1R"],
    "Fibroblasts": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "ACTA2"],
    "Endothelial cells": ["PECAM1", "VWF", "KDR", "CLDN5", "CDH5", "ENG"],
    "Epithelial cells": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CDH1", "MUC1"],
    "NK cells": ["NKG7", "GNLY", "PRF1", "KLRD1", "GZMB", "KLRF1"],
}

TYPE_COLORS = {
    "T cells": "#7B3FB2",
    "B cells": "#2B83BA",
    "Plasma cells": "#35B7B4",
    "Monocytes and Macrophages": "#B8860B",
    "Fibroblasts": "#F2B31A",
    "Endothelial cells": "#4DA3D9",
    "Epithelial cells": "#E76F51",
    "NK cells": "#D94A8C",
    "Unassigned": "#BDBDBD",
}

PANEL_BG = "#E8E8E8"
CYAN = "#00CFE8"


def _slug(value: str) -> str:
    return "".join(ch.lower() if ch.isalnum() else "_" for ch in str(value)).strip("_")


def _deep_purple_magma() -> LinearSegmentedColormap:
    base = colormaps["magma"]
    colors = base(np.linspace(0.18, 1.0, 256))
    return LinearSegmentedColormap.from_list("magma_deep_purple_highres", colors)


SIGNATURE_CMAP = _deep_purple_magma()


def _font(size: int, bold: bool = False) -> ImageFont.ImageFont:
    candidates = []
    if bold:
        candidates.extend([r"C:\Windows\Fonts\arialbd.ttf", r"C:\Windows\Fonts\calibrib.ttf"])
    candidates.extend([r"C:\Windows\Fonts\arial.ttf", r"C:\Windows\Fonts\calibri.ttf"])
    for path in candidates:
        try:
            return ImageFont.truetype(path, size=size)
        except Exception:
            pass
    return ImageFont.load_default()


def _header(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return next(csv.reader(handle))


def _score_table(expr_path: Path, meta_path: Path, chunk_size: int) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    header = _header(expr_path)
    genes = set(header)
    available = {k: [g for g in v if g in genes] for k, v in MARKERS.items()}
    usecols = ["cell", *sorted({g for vals in available.values() for g in vals})]
    chunks: list[pd.DataFrame] = []
    for chunk in pd.read_csv(expr_path, usecols=usecols, chunksize=chunk_size):
        cell_id = chunk["cell"].astype(str)
        scores = {}
        for cell_type, marker_genes in available.items():
            if marker_genes:
                scores[cell_type] = np.log1p(chunk[marker_genes].astype("float32")).mean(axis=1).to_numpy()
            else:
                scores[cell_type] = np.zeros(len(chunk), dtype="float32")
        keys = list(MARKERS)
        mat = np.vstack([scores[k] for k in keys]).T
        order = np.argsort(mat, axis=1)
        best = order[:, -1]
        second = order[:, -2]
        bestv = mat[np.arange(len(chunk)), best]
        secondv = mat[np.arange(len(chunk)), second]
        labels = np.array(keys, dtype=object)[best]
        labels[(bestv <= 0.08) | ((bestv - secondv) <= 0.025)] = "Unassigned"
        chunks.append(
            pd.DataFrame(
                {
                    "raw_cell": cell_id,
                    "best_score": bestv,
                    "score_margin": bestv - secondv,
                    "cell_type": labels,
                }
            )
        )
    scored = pd.concat(chunks, ignore_index=True)
    meta = pd.read_csv(meta_path, usecols=["center_x", "center_y"])
    if len(meta) != len(scored):
        raise ValueError(f"metadata/expression row mismatch: {meta_path}")
    scored["center_x"] = pd.to_numeric(meta["center_x"], errors="coerce")
    scored["center_y"] = pd.to_numeric(meta["center_y"], errors="coerce")
    scored = scored.dropna(subset=["center_x", "center_y"]).copy()
    return scored, available


def _select_spatial_cells(scored: pd.DataFrame, target_type: str, n_spots: int, min_target: int) -> pd.DataFrame:
    target = scored[scored["cell_type"].eq(target_type)].copy()
    if len(target) < min_target:
        raise ValueError(f"target {target_type} has only {len(target)} cells")
    high = target.nlargest(min(2000, len(target)), "best_score")
    anchor_x = float(high["center_x"].median())
    anchor_y = float(high["center_y"].median())
    dist = (scored["center_x"] - anchor_x) ** 2 + (scored["center_y"] - anchor_y) ** 2
    selected = scored.assign(_dist=dist).nsmallest(n_spots, "_dist").copy()
    if int(selected["cell_type"].eq(target_type).sum()) < min_target:
        selected_ids = set(selected["raw_cell"].astype(str))
        add = target[~target["raw_cell"].astype(str).isin(selected_ids)].assign(
            _dist=(target["center_x"] - anchor_x) ** 2 + (target["center_y"] - anchor_y) ** 2
        ).nsmallest(min_target - int(selected["cell_type"].eq(target_type).sum()), "_dist")
        removable = selected[~selected["cell_type"].eq(target_type)].nlargest(len(add), "_dist").index
        selected = pd.concat([selected.drop(index=removable), add], ignore_index=True)
    return _trim_spatial_outliers(selected.drop(columns=["_dist"], errors="ignore"), keep_quantile=0.985)


def _trim_spatial_outliers(df: pd.DataFrame, keep_quantile: float) -> pd.DataFrame:
    """Remove edge outliers that make cell-level maps look like a core plus scattered satellites."""
    x = df["center_x"].astype(float)
    y = df["center_y"].astype(float)
    cx = float(x.median())
    cy = float(y.median())
    dist = np.sqrt((x - cx) ** 2 + (y - cy) ** 2)
    cutoff = float(dist.quantile(keep_quantile))
    trimmed = df.loc[dist <= cutoff].copy()
    return trimmed.reset_index(drop=True)


def _read_selected_expr(expr_path: Path, selected: pd.DataFrame, chunk_size: int) -> pd.DataFrame:
    wanted = set(selected["raw_cell"].astype(str))
    pieces = []
    for chunk in pd.read_csv(expr_path, chunksize=chunk_size):
        mask = chunk["cell"].astype(str).isin(wanted)
        if mask.any():
            pieces.append(chunk.loc[mask].copy())
    if not pieces:
        raise ValueError(f"no selected cells found in {expr_path}")
    expr = pd.concat(pieces, ignore_index=True)
    expr["raw_cell"] = expr["cell"].astype(str)
    expr = expr.set_index("raw_cell")
    expr = expr.drop(columns=["cell"], errors="ignore")
    expr = expr.apply(pd.to_numeric, errors="coerce").fillna(0.0).astype("float32")
    expr = expr.loc[selected["raw_cell"].astype(str)]
    expr.index = selected["cell_id"].tolist()
    return expr


def _read_selected_expr_with_id(expr_path: Path, selected: pd.DataFrame, id_col: str, out_col: str, chunk_size: int) -> pd.DataFrame:
    wanted = set(selected[id_col].astype(str))
    pieces = []
    for chunk in pd.read_csv(expr_path, chunksize=chunk_size):
        mask = chunk["cell"].astype(str).isin(wanted)
        if mask.any():
            pieces.append(chunk.loc[mask].copy())
    if not pieces:
        raise ValueError(f"no selected cells found in {expr_path}")
    expr = pd.concat(pieces, ignore_index=True)
    expr["raw_cell"] = expr["cell"].astype(str)
    expr = expr.set_index("raw_cell")
    expr = expr.drop(columns=["cell"], errors="ignore")
    expr = expr.apply(pd.to_numeric, errors="coerce").fillna(0.0).astype("float32")
    expr = expr.loc[selected[id_col].astype(str)]
    expr.index = selected[out_col].tolist()
    return expr


def _select_sc_reference(scored: pd.DataFrame, spatial: pd.DataFrame, n_cells: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    used = set(spatial["raw_cell"].astype(str))
    pool = scored[~scored["raw_cell"].astype(str).isin(used)].copy()
    pool = pool[pool["cell_type"].astype(str).ne("Unassigned")].copy()
    spatial_counts = spatial["cell_type"].value_counts()
    pieces = []
    for cell_type, count in spatial_counts.items():
        type_pool = pool[pool["cell_type"].eq(cell_type)]
        if type_pool.empty:
            continue
        take = min(len(type_pool), max(20, int(count)))
        idx = rng.choice(type_pool.index.to_numpy(), size=take, replace=False)
        pieces.append(type_pool.loc[idx])
    selected = pd.concat(pieces, ignore_index=True) if pieces else pool.sample(min(n_cells, len(pool)), random_state=seed)
    if len(selected) < n_cells:
        used2 = set(selected["raw_cell"].astype(str))
        rest = pool[~pool["raw_cell"].astype(str).isin(used2)]
        add_n = min(len(rest), n_cells - len(selected))
        if add_n > 0:
            selected = pd.concat([selected, rest.sample(add_n, random_state=seed + 1)], ignore_index=True)
    elif len(selected) > n_cells:
        selected = selected.sample(n_cells, random_state=seed).reset_index(drop=True)
    selected = selected.reset_index(drop=True)
    selected["cell_id"] = [f"sc_{i}_{x}" for i, x in enumerate(selected["raw_cell"].astype(str))]
    return selected


def _build_panel(expr: pd.DataFrame, cell_types: pd.Series, target_type: str, max_genes: int) -> pd.DataFrame:
    type_means = expr.groupby(cell_types).mean(numeric_only=True)
    target = type_means.loc[target_type]
    others = type_means.drop(index=target_type, errors="ignore")
    other_mean = others.mean(axis=0)
    other_max = others.max(axis=0)
    specificity = (target + 1e-6) / (other_max + 1e-6)
    detect_frac = (expr > 0).mean(axis=0)
    df = pd.DataFrame(
        {
            "gene": expr.columns,
            "target_type": target_type,
            "target_mean": target.to_numpy(),
            "other_mean": other_mean.to_numpy(),
            "other_max": other_max.to_numpy(),
            "specificity_max": specificity.to_numpy(),
            "detect_frac": detect_frac.to_numpy(),
        }
    )
    marker_bonus = df["gene"].isin(MARKERS.get(target_type, []))
    keep = (df["target_mean"] > 0.01) & (df["detect_frac"] > 0.005) & ((df["specificity_max"] > 1.15) | marker_bonus)
    df = df.loc[keep].copy()
    df["score"] = df["target_mean"] * np.log1p(df["specificity_max"].clip(0, 20)) + marker_bonus.loc[df.index].astype(float) * 0.25
    df = df.sort_values(["score", "specificity_max", "target_mean"], ascending=False).head(max_genes)
    if len(df) < 8:
        raise ValueError(f"too few target panel genes for {target_type}: {len(df)}")
    return df


def _write_config(root: Path, sample: str, hvg_path: Path) -> None:
    cfg = {
        "paths": {"sc_expr": "unused", "sc_meta": None, "st_expr": None, "st_meta": None, "svg_marker_whitelist": None},
        "qc": {"sc_min_genes": 0, "sc_max_genes": 100000, "sc_max_mt": 100, "st_min_genes": 0, "st_max_genes": "Inf", "st_max_mt": 100, "hvg_nfeatures": 2000, "mt_pattern": "^(MT-|mt-)"},
        "gene_filter": {"min_cells_sc": 0, "min_cells_st": 0},
        "stage3": {
            "strong_th": 0.7,
            "weak_th": 0.4,
            "st_cluster_k": 0,
            "unknown_floor": 0.3,
            "min_cells_rare_type": 20,
            "eps": 1.0e-8,
            "plugin_genes_path": str(hvg_path.relative_to(root)).replace("\\", "/"),
            "gene_weights_path": None,
            "auto_missing_detection": {
                "enable": True,
                "method": "adaptive_low_support",
                "min_cells": 50,
                "robust_z_th": -1.8,
                "soft_z_th": -0.75,
                "require_masked_for_soft": False,
                "require_masked_for_hard": False,
                "max_fraction_types": 0.35,
                "max_types": 2,
                "action": "mark_unknown",
                "require_confirmation": True,
                "confirmation_max_support_score": 0.75,
                "confirmation_use_masked_missing": True,
                "confirmation_use_marker_identity": True,
                "confirmation_marker_identity_z_th": 0.1,
                "confirmation_marker_support_score_th": 0.65,
            },
            "masked_missing_detection": {
                "enable": True,
                "apply_to_auto_missing": True,
                "neighbor_cosine_th": 0.86,
                "neighbor_cell_ratio_min": 0.35,
                "marker_top_n": 20,
                "min_identity_markers": 3,
                "min_marker_specificity": 1.05,
                "min_marker_type_mean": 0.001,
                "min_marker_st_detect_frac": 0.001,
                "st_presence_quantile": 0.90,
                "identity_z_th": -0.55,
                "pressure_z_th": 0.0,
                "min_support_score_for_apply": 0.70,
                "max_types": 2,
            },
            "marker_identity_diagnostics": {
                "enable": True,
                "marker_top_n": 80,
                "min_identity_markers": 3,
                "min_all_specificity": 1.05,
                "min_neighbor_specificity": 1.02,
                "min_marker_type_mean": 0.001,
                "min_marker_st_detect_frac": 0.001,
                "st_presence_quantile": 0.9,
                "depleted_z_th": -0.55,
            },
        },
    }
    path = root / "configs" / "datasets" / f"{sample}.yaml"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(cfg, sort_keys=False, allow_unicode=True), encoding="utf-8")


def _prepare_dataset(root: Path, raw_sample: str, target_type: str, args: argparse.Namespace) -> dict[str, object]:
    raw_dir = root / "data" / "raw" / "high" / raw_sample
    expr_path = raw_dir / f"{raw_sample}_cell_by_gene.csv"
    meta_path = raw_dir / f"{raw_sample}_cell_metadata.csv"
    if not expr_path.exists() or not meta_path.exists():
        raise FileNotFoundError(raw_dir)

    source_sample = f"highres_{_slug(raw_sample)}"
    masked_sample = f"{source_sample}_profile_mask_{_slug(target_type)}"
    source_root = root / "data" / "processed" / source_sample / "stage1_preprocess"
    masked_root = root / "data" / "processed" / masked_sample / "stage1_preprocess"
    if args.overwrite:
        shutil.rmtree(source_root.parent, ignore_errors=True)
        shutil.rmtree(masked_root.parent, ignore_errors=True)
    if (masked_root / "exported" / "st_expression_normalized.csv").exists():
        info_path = masked_root / "fig2d_profile_mask_info.json"
        info = json.loads(info_path.read_text(encoding="utf-8")) if info_path.exists() else {}
        return {
            "raw_sample": raw_sample,
            "source_sample": source_sample,
            "profile_mask_sample": masked_sample,
            "masked_target_type": target_type,
            "n_selected_spots": info.get("n_selected_spots"),
            "target_cells_in_selected": info.get("target_cells_in_selected"),
            "target_panel_genes": info.get("target_panel_genes"),
            "prepared": "existing",
        }

    scored, available_markers = _score_table(expr_path, meta_path, args.chunk_size)
    spatial = _select_spatial_cells(scored, target_type, args.n_spots, args.min_target_cells)
    spatial = spatial.reset_index(drop=True)
    spatial["cell_id"] = [f"spatial_{i}_{x}" for i, x in enumerate(spatial["raw_cell"].astype(str))]
    spatial["spot_id"] = [f"spot_{i}_{x}" for i, x in enumerate(spatial["raw_cell"].astype(str))]
    sc_ref = _select_sc_reference(scored, spatial, args.n_sc_cells, seed=42 + len(raw_sample))

    sc_expr = _read_selected_expr_with_id(expr_path, sc_ref, "raw_cell", "cell_id", args.chunk_size)
    st_expr = _read_selected_expr_with_id(expr_path, spatial, "raw_cell", "spot_id", args.chunk_size)
    sc_types = sc_ref.set_index("cell_id")["cell_type"].reindex(sc_expr.index)
    panel = _build_panel(sc_expr, sc_types, target_type, args.max_panel_genes)
    masked_expr = st_expr.copy()
    marker_genes = panel["gene"].tolist()
    masked_expr.loc[:, marker_genes] = masked_expr.loc[:, marker_genes] * args.mask_scale

    for stage1, spatial_expr in [(source_root, st_expr), (masked_root, masked_expr)]:
        export = stage1 / "exported"
        export.mkdir(parents=True, exist_ok=True)
        sc_expr.to_csv(export / "sc_expression_normalized.csv", index_label="cell_id")
        spatial_expr.to_csv(export / "st_expression_normalized.csv", index_label="spot_id")
        coords = spatial[["spot_id", "center_y", "center_x"]].rename(columns={"center_y": "row", "center_x": "col"}).set_index("spot_id")
        coords.to_csv(export / "st_coordinates.csv", index_label="spot_id")
        meta = sc_ref[["cell_id", "cell_type", "raw_cell", "best_score", "score_margin", "center_x", "center_y"]].copy()
        meta["sc_meta"] = meta["cell_type"]
        meta.to_csv(export / "sc_metadata.csv", index=False)
        st_meta = spatial[["spot_id", "cell_type", "raw_cell", "best_score", "score_margin", "center_x", "center_y"]].copy()
        st_meta.to_csv(export / "st_metadata.csv", index=False)
        (export / "sim_info.json").write_text(json.dumps({"cells_per_spot": 1}, indent=2), encoding="utf-8")
        (stage1 / "hvg_genes.txt").write_text("\n".join(sc_expr.columns.astype(str)) + "\n", encoding="utf-8")
        _write_config(root, stage1.parent.name, stage1 / "hvg_genes.txt")

    panel.to_csv(masked_root / "fig2d_profile_mask_gene_panel.csv", index=False)
    (masked_root / "fig2d_profile_mask_info.json").write_text(
        json.dumps(
            {
                "raw_sample": raw_sample,
                "source_sample": source_sample,
                "profile_mask_sample": masked_sample,
                "target_type": target_type,
                "n_selected_spots": int(len(spatial)),
                "n_sc_reference_cells": int(len(sc_ref)),
                "target_cells_in_selected": int(spatial["cell_type"].eq(target_type).sum()),
                "target_cells_in_sc_reference": int(sc_ref["cell_type"].eq(target_type).sum()),
                "spatial_sc_raw_overlap": int(len(set(spatial["raw_cell"].astype(str)) & set(sc_ref["raw_cell"].astype(str)))),
                "mask_scale": float(args.mask_scale),
                "available_markers": available_markers,
                "target_panel_genes": marker_genes,
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    return {
        "raw_sample": raw_sample,
        "source_sample": source_sample,
        "profile_mask_sample": masked_sample,
        "masked_target_type": target_type,
        "n_selected_spots": int(len(spatial)),
        "n_sc_reference_cells": int(len(sc_ref)),
        "target_cells_in_selected": int(spatial["cell_type"].eq(target_type).sum()),
        "target_cells_in_sc_reference": int(sc_ref["cell_type"].eq(target_type).sum()),
        "spatial_sc_raw_overlap": int(len(set(spatial["raw_cell"].astype(str)) & set(sc_ref["raw_cell"].astype(str)))),
        "target_panel_genes": marker_genes,
        "prepared": "new",
    }


def _run(cmd: list[str], root: Path, log_path: Path, env: dict[str, str] | None = None) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as handle:
        proc = subprocess.run(cmd, cwd=root, stdout=handle, stderr=subprocess.STDOUT, env=env)
    if proc.returncode != 0:
        tail = "\n".join(log_path.read_text(encoding="utf-8", errors="replace").splitlines()[-80:])
        raise RuntimeError(f"command failed: {' '.join(cmd)}\n--- log tail ---\n{tail}")


def _stage3_stage4(root: Path, sample: str, args: argparse.Namespace) -> dict[str, object]:
    py = sys.executable
    log_dir = root / "logs" / "highres_profile_mask_mapping"
    _run([py, "-m", "src.stages.stage3_type_plugin", "--sample", sample, "--project_root", str(root), "--sc_expr_source", "normalized"], root, log_dir / f"{sample}.stage3.log")
    env = os.environ.copy()
    env["CYTOSPACE_SKIP_ASSIGNED_EXPRESSION"] = "1"
    common = [
        py,
        "-m",
        "src.stages.stage4_cytospace",
        "--sample",
        sample,
        "--project_root",
        str(root),
        "--n_processors",
        "1",
        "--mapping_cells_per_spot",
        "1",
        "--sc_expr_source",
        "normalized",
        "--no_sampling_sub_spots",
    ]
    _run(
        common
        + [
            "--missing_type",
            "__NO_MISSING__",
            "--filter_mode",
            "none",
            "--cell_type_column",
            "sc_meta",
            "--filter_scope",
            "unsupported_all",
            "--stage4_suffix",
            "_baseline_highres",
        ],
        root,
        log_dir / f"{sample}.stage4_baseline.log",
        env=env,
    )
    _run(
        common
        + [
            "--missing_type",
            "__AUTO__",
            "--filter_mode",
            "plugin_unknown",
            "--cell_type_column",
            "plugin_type",
            "--filter_scope",
            "missing_detected_only",
            "--stage4_suffix",
            "_route2_highres",
        ],
        root,
        log_dir / f"{sample}.stage4_route2.log",
        env=env,
    )
    summary_path = root / "result" / sample / "stage3_typematch" / "stage3_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8-sig")) if summary_path.exists() else {}
    auto = summary.get("auto_missing_detection", {}) or {}
    action = summary.get("action_overview", {}) or {}
    return {
        "stage3_missing_types": auto.get("auto_missing_types") or action.get("missing_types") or [],
        "stage3_dropped_cells": int((action.get("by_cell_count", {}) or {}).get("Dropped", 0)),
    }


def _stage3_existing_summary(root: Path, sample: str) -> dict[str, object]:
    summary_path = root / "result" / sample / "stage3_typematch" / "stage3_summary.json"
    if not summary_path.exists():
        return {}
    summary = json.loads(summary_path.read_text(encoding="utf-8-sig"))
    auto = summary.get("auto_missing_detection", {}) or {}
    action = summary.get("action_overview", {}) or {}
    return {
        "stage3_missing_types": auto.get("auto_missing_types") or action.get("missing_types") or [],
        "stage3_dropped_cells": int((action.get("by_cell_count", {}) or {}).get("Dropped", 0)),
    }


def _load_target_counts(root: Path, sample: str, target_type: str, suffix: str) -> pd.DataFrame:
    path = root / "result" / sample / f"stage4_cytospace{suffix}" / "cytospace_output" / "cell_type_assignments_by_spot.csv"
    df = pd.read_csv(path).rename(columns={pd.read_csv(path, nrows=0).columns[0]: "spot_id"})
    col = next((c for c in df.columns if str(c).casefold() == target_type.casefold()), None)
    df["target_present"] = pd.to_numeric(df[col], errors="coerce").fillna(0).gt(0) if col else False
    return df[["spot_id", "target_present"]]


def _clip_plot_body(df: pd.DataFrame, keep_quantile: float) -> pd.DataFrame:
    """Clip only the rendered view, not the underlying mapping results."""
    if len(df) == 0 or keep_quantile >= 1:
        return df
    x = df["col"].astype(float)
    y = df["row"].astype(float)
    cx = float(x.median())
    cy = float(y.median())
    dist = np.sqrt((x - cx) ** 2 + (y - cy) ** 2)
    return df.loc[dist <= float(dist.quantile(keep_quantile))].copy()


def _render_one(root: Path, row: dict[str, object], out_root: Path) -> Path:
    source = str(row["source_sample"])
    sample = str(row["profile_mask_sample"])
    target = str(row["masked_target_type"])
    src_export = root / "data" / "processed" / source / "stage1_preprocess" / "exported"
    mask_export = root / "data" / "processed" / sample / "stage1_preprocess" / "exported"
    panel = pd.read_csv(root / "data" / "processed" / sample / "stage1_preprocess" / "fig2d_profile_mask_gene_panel.csv")
    genes = panel["gene"].astype(str).head(30).tolist()

    coords = pd.read_csv(src_export / "st_coordinates.csv")
    source_expr = pd.read_csv(src_export / "st_expression_normalized.csv", usecols=["spot_id", *genes])
    mask_expr = pd.read_csv(mask_export / "st_expression_normalized.csv", usecols=["spot_id", *genes])
    meta = pd.read_csv(src_export / "st_metadata.csv", usecols=["spot_id", "cell_type"])

    def build(expr: pd.DataFrame) -> pd.DataFrame:
        df = coords.merge(expr, on="spot_id", how="inner")
        df["target_signature"] = np.log1p(df[genes].astype(float)).mean(axis=1)
        return df

    source_df = build(source_expr).merge(meta, on="spot_id", how="left")
    masked_df = build(mask_expr).merge(meta, on="spot_id", how="left")
    baseline = masked_df.merge(_load_target_counts(root, sample, target, "_baseline_highres"), on="spot_id", how="left")
    route2 = masked_df.merge(_load_target_counts(root, sample, target, "_route2_highres"), on="spot_id", how="left")
    baseline["target_present"] = baseline["target_present"].fillna(False).astype(bool)
    route2["target_present"] = route2["target_present"].fillna(False).astype(bool)

    plot_keep_quantile = 0.96 if str(row["raw_sample"]) in {"HumanLungCancerPatient1", "HumanMelanomaPatient2"} else 0.985
    source_df = _clip_plot_body(source_df, plot_keep_quantile)
    masked_df = _clip_plot_body(masked_df, plot_keep_quantile)
    baseline = _clip_plot_body(baseline, plot_keep_quantile)
    route2 = _clip_plot_body(route2, plot_keep_quantile)

    combined = pd.concat([source_df["target_signature"], masked_df["target_signature"]], ignore_index=True)
    norm = Normalize(float(np.percentile(combined, 1)), float(np.percentile(combined, 99)))

    plt.rcParams.update({"font.family": "Arial", "svg.fonttype": "none"})
    fig, axes = plt.subplots(1, 4, figsize=(15.8, 4.6), dpi=220)
    titles = [
        ("A  Original cell-level map", "target marker signature", source_df, "truth"),
        ("B  Profile-mask map", "target signature after masking", masked_df, None),
        ("C  CytoSPACE baseline", "mapped target outline", baseline, "target_present"),
        ("D  SVTuner + CytoSPACE", "mapped target outline", route2, "target_present"),
    ]
    sc = None
    for ax, (title, subtitle, df, outline_col) in zip(axes, titles):
        ax.set_facecolor(PANEL_BG)
        sc = ax.scatter(
            df["col"],
            -df["row"],
            c=df["target_signature"],
            s=15.0,
            cmap=SIGNATURE_CMAP,
            norm=norm,
            linewidths=0,
            alpha=0.96,
            rasterized=True,
        )
        if outline_col is None:
            hi = df.iloc[0:0]
        elif outline_col == "truth":
            hi = df[df["cell_type"].astype(str).eq(target)]
        else:
            hi = df[df[outline_col].astype(bool)]
        if not hi.empty:
            # Keep cyan outlines inside the main tissue body. Some cell-level
            # MERSCOPE target calls are sparse satellites at the edge; large
            # hollow outlines make them look detached from the tissue map.
            x = df["col"].astype(float)
            y = df["row"].astype(float)
            cx = float(x.median())
            cy = float(y.median())
            dist = np.sqrt((x - cx) ** 2 + (y - cy) ** 2)
            outline_cutoff = float(dist.quantile(0.92))
            hi_dist = np.sqrt((hi["col"].astype(float) - cx) ** 2 + (hi["row"].astype(float) - cy) ** 2)
            hi = hi.loc[hi_dist <= outline_cutoff]
        if not hi.empty:
            ax.scatter(hi["col"], -hi["row"], s=34, facecolors="none", edgecolors=CYAN, linewidths=0.9, alpha=0.95)
        ax.set_title(f"{title}\n{subtitle}", loc="left", fontsize=9.2, fontweight="bold")
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([])
        ax.set_yticks([])
    axes[0].set_ylabel("spatial row", fontsize=10)
    fig.suptitle(f"Cell-level profile-mask mapping: {row['raw_sample']} ({target})", x=0.03, ha="left", fontsize=13.5, fontweight="bold")
    fig.subplots_adjust(left=0.03, right=0.92, top=0.78, bottom=0.12, wspace=0.03)
    cax = fig.add_axes([0.935, 0.20, 0.012, 0.55])
    cb = fig.colorbar(sc, cax=cax)
    cb.set_label(f"{target} marker signature", fontsize=7.4)
    cb.ax.tick_params(labelsize=7)

    out_dir = out_root / sample
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / f"{_slug(target)}_mask_ABCD_mapping.png"
    fig.savefig(out_png, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[OK] wrote: {out_png}")
    return out_png


def _stack(pngs: list[Path], out_root: Path) -> Path:
    images = [Image.open(p).convert("RGB") for p in pngs]
    target_w = 2300
    resized = []
    for img in images:
        scale = target_w / img.size[0]
        resized.append(img.resize((target_w, int(round(img.size[1] * scale))), Image.Resampling.LANCZOS))
    outer, title_h, gap = 28, 70, 16
    w = target_w + outer * 2
    h = outer * 2 + title_h + sum(i.size[1] for i in resized) + gap * (len(resized) - 1)
    canvas = Image.new("RGB", (w, h), "white")
    draw = ImageDraw.Draw(canvas)
    title = "Cell-level profile-mask mapping across five high-resolution datasets"
    font = _font(34, bold=True)
    bb = draw.textbbox((0, 0), title, font=font)
    draw.text(((w - (bb[2] - bb[0])) // 2, outer), title, fill=(20, 20, 20), font=font)
    y = outer + title_h
    for img in resized:
        canvas.paste(img, (outer, y))
        y += img.size[1] + gap
    out = out_root / "highres_profile_mask_mapping_stack_5x4.png"
    canvas.save(out, optimize=True)
    print(f"[OK] wrote: {out}")
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project_root", default=".")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--prepare_only", action="store_true")
    parser.add_argument("--skip_mapping", action="store_true")
    parser.add_argument("--chunk_size", type=int, default=50000)
    parser.add_argument("--n_spots", type=int, default=2600)
    parser.add_argument("--n_sc_cells", type=int, default=5200)
    parser.add_argument("--min_target_cells", type=int, default=180)
    parser.add_argument("--max_panel_genes", type=int, default=80)
    parser.add_argument("--mask_scale", type=float, default=0.02)
    args = parser.parse_args()

    root = Path(args.project_root).resolve()
    rows: list[dict[str, object]] = []
    for raw_sample, target_type in DATASETS.items():
        print(f"[PREP] {raw_sample}: target={target_type}")
        row = _prepare_dataset(root, raw_sample, target_type, args)
        if not args.prepare_only and not args.skip_mapping:
            row.update(_stage3_stage4(root, str(row["profile_mask_sample"]), args))
        elif args.skip_mapping:
            row.update(_stage3_existing_summary(root, str(row["profile_mask_sample"])))
        rows.append(row)

    out_root = root / "visualizations" / "highres_profile_mask_mapping"
    out_root.mkdir(parents=True, exist_ok=True)
    manifest = pd.DataFrame(rows)
    manifest_path = out_root / "highres_profile_mask_mapping_manifest.csv"
    manifest.to_csv(manifest_path, index=False)
    (out_root / "highres_profile_mask_mapping_manifest.json").write_text(
        json.dumps(rows, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(f"[OK] wrote: {manifest_path}")

    if not args.prepare_only:
        pngs = [_render_one(root, row, out_root) for row in rows]
        _stack(pngs, out_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
