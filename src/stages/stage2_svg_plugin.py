"""
Stage 2: SVG + HVG dynamic weighting plugin.

Responsibilities:
- Ensure Stage1 RDS exports are available (R helper for RDS->CSV/NPZ).
- Compute HVG and SVG scores, normalize, and fuse into final_weight.
- Output gene_weights.csv / plugin_genes.txt to processed dir.
- Output svg_filter_sensitivity.json, stage2_summary.json to result dir.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import yaml
from scipy.spatial import cKDTree

# Optional plotting removed (no PNG outputs).


# ----------------------------
# Config helpers
# ----------------------------


def load_yaml(path: Path) -> dict:
    if not path.exists() or path.stat().st_size == 0:
        return {}
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def merge_configs(project_cfg: dict, dataset_cfg: dict, defaults: dict) -> dict:
    cfg = defaults.copy()
    for src in (project_cfg, dataset_cfg):
        if not isinstance(src, dict):
            continue
        for k, v in src.get("stage2", {}).items():
            cfg[k] = v
    return cfg


# ----------------------------
# Data I/O
# ----------------------------


def detect_project_root() -> Path:
    here = Path(__file__).resolve()
    # .../src/stages/stage2_svg_plugin.py -> root at parents[2]
    try:
        return here.parents[2]
    except IndexError:
        return here.parent.parent


def run_r_export(project_root: Path, sample: str, r_helper: Path, skip: bool = False):
    if skip:
        return
    if not r_helper.exists():
        raise FileNotFoundError(f"R export helper not found: {r_helper}")
    cmd = [
        "Rscript",
        str(r_helper),
        "--sample",
        sample,
        "--project_root",
        str(project_root),
    ]
    print(f"[Stage2] Running R export helper: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr)
        raise RuntimeError("R export helper failed; see logs above.")


def read_expression_table(path: Path, index_col: str) -> pd.DataFrame:
    """
    Expect CSV with first column as index.
    """
    df = pd.read_csv(path, index_col=0)
    df.index.name = index_col
    return df


def ensure_exports(project_root: Path, sample: str, r_helper: Path, skip_export: bool):
    stage1_dir = project_root / "data" / "processed" / sample / "stage1_preprocess"
    export_dir = stage1_dir / "exported"
    export_dir.mkdir(parents=True, exist_ok=True)

    sc_expr = export_dir / "sc_expression_normalized.csv"
    st_expr = export_dir / "st_expression_normalized.csv"
    sc_meta = export_dir / "sc_metadata.csv"
    st_coords = export_dir / "st_coordinates.csv"

    required = (sc_expr, st_expr, st_coords)
    if all(p.exists() for p in required):
        if sc_meta.exists():
            print("[Stage2] Found cached Stage1 exports.")
        else:
            print("[Stage2] Cached exports found (sc_meta missing; proceeding without it).")
    else:
        print("[Stage2] Stage1 exports missing; invoking R helper for RDS->CSV.")
        run_r_export(project_root, sample, r_helper, skip=skip_export)

    return {
        "export_dir": export_dir,
        "sc_expr": sc_expr,
        "st_expr": st_expr,
        "sc_meta": sc_meta,
        "st_coords": st_coords,
        "common_genes": stage1_dir / "common_genes.txt",
        "hvg_genes": stage1_dir / "hvg_genes.txt",
    }


def _resolve_truth_spot_type_path(project_root: Path, sample: str, cfg: dict) -> Optional[Path]:
    hint = cfg.get("svg_cost_truth_path")
    if isinstance(hint, str) and hint.strip():
        p = Path(hint)
        return p if p.exists() else None
    sim_dir = project_root / "data" / "sim" / sample
    for name in ("sim_truth_spot_type_fraction.csv", "spot_true_type_fraction.csv"):
        p = sim_dir / name
        if p.exists():
            return p
    return None


def _compute_truth_corr(
    st_expr: pd.DataFrame,
    truth_df: pd.DataFrame,
    metric: str,
    eps: float,
) -> Tuple[pd.Series, pd.Series, Dict[str, object]]:
    if "spot_id" in truth_df.columns:
        truth_df = truth_df.set_index("spot_id")
    common = st_expr.index.intersection(truth_df.index)
    info = {"n_spots_common": int(len(common)), "n_types": int(truth_df.shape[1])}
    if len(common) < 3:
        return (
            pd.Series(np.nan, index=st_expr.columns),
            pd.Series([None] * st_expr.shape[1], index=st_expr.columns),
            {**info, "status": "insufficient_spots"},
        )
    X = st_expr.loc[common].to_numpy(dtype=float)
    Y = truth_df.loc[common].to_numpy(dtype=float)
    type_names = list(truth_df.columns)

    std_y = Y.std(axis=0, ddof=0)
    keep = std_y > 0
    if not keep.any():
        return (
            pd.Series(np.nan, index=st_expr.columns),
            pd.Series([None] * st_expr.shape[1], index=st_expr.columns),
            {**info, "status": "no_type_variance"},
        )
    Y = Y[:, keep]
    type_names = [t for t, ok in zip(type_names, keep) if ok]

    if metric == "spearman":
        X = pd.DataFrame(X).rank(axis=0).to_numpy(dtype=float)
        Y = pd.DataFrame(Y).rank(axis=0).to_numpy(dtype=float)

    Xc = X - X.mean(axis=0, keepdims=True)
    Yc = Y - Y.mean(axis=0, keepdims=True)
    std_x = Xc.std(axis=0, ddof=0)
    std_y = Yc.std(axis=0, ddof=0)
    denom = (std_x[:, None] * std_y[None, :]) + eps
    corr = (Xc.T @ Yc) / (len(common) * denom)
    corr = np.nan_to_num(corr, nan=-np.inf, neginf=-np.inf, posinf=np.inf)
    best_idx = np.argmax(corr, axis=1)
    corr_max = corr[np.arange(corr.shape[0]), best_idx]
    best_type = [type_names[i] if np.isfinite(corr_max[j]) else None for j, i in enumerate(best_idx)]
    return (
        pd.Series(corr_max, index=st_expr.columns),
        pd.Series(best_type, index=st_expr.columns),
        {**info, "status": "ok"},
    )


# ----------------------------
# SVG / HVG calculations
# ----------------------------


def normalize_scores(series: pd.Series, method: str = "rank") -> pd.Series:
    if method == "rank":
        return series.rank(method="average") / len(series)
    if method == "zscore":
        mean = series.mean()
        std = series.std(ddof=0)
        if std == 0:
            return pd.Series(np.zeros_like(series), index=series.index)
        return (series - mean) / std
    if method == "minmax":
        min_v, max_v = series.min(), series.max()
        if max_v == min_v:
            return pd.Series(np.zeros_like(series), index=series.index)
        return (series - min_v) / (max_v - min_v)
    raise ValueError(f"Unknown normalization method: {method}")


def compute_hvg_scores(sc_expr: pd.DataFrame, norm_method: str) -> Tuple[pd.Series, pd.Series]:
    # Variance across cells per gene
    var_sc = sc_expr.var(axis=0)
    hvg_score = np.log1p(var_sc)
    hvg_norm = normalize_scores(pd.Series(hvg_score, index=sc_expr.columns), method=norm_method)
    return hvg_score, hvg_norm


@dataclass
class NeighborGraph:
    indices: np.ndarray  # shape (n, k)
    weights: np.ndarray  # shape (n, k)
    weight_sum: float


def _select_coord_cols(coords: pd.DataFrame) -> pd.DataFrame:
    if {"row", "col"}.issubset(coords.columns):
        return coords[["row", "col"]]
    if {"x", "y"}.issubset(coords.columns):
        return coords[["x", "y"]]
    raise KeyError("st_coordinates.csv 需要包含 row/col 或 x/y 列用于构建邻接。")


def build_knn_graph(coords: pd.DataFrame, k_neighbors: int = 6) -> NeighborGraph:
    coord_df = _select_coord_cols(coords)
    points = coord_df.to_numpy()
    tree = cKDTree(points)
    dists, idx = tree.query(points, k=k_neighbors + 1)  # include self at position 0
    idx = idx[:, 1:]  # drop self
    dists = dists[:, 1:]
    # Avoid div by zero; use uniform weights if distance == 0
    weights = np.where(dists > 0, 1.0 / dists, 1.0)
    weight_sum = float(weights.sum())
    return NeighborGraph(indices=idx, weights=weights, weight_sum=weight_sum)


def morans_I(x: np.ndarray, graph: NeighborGraph) -> float:
    x_mean = x.mean()
    diff = x - x_mean
    denom = np.sum(diff ** 2)
    if denom <= 0 or graph.weight_sum <= 0:
        return 0.0
    num = np.sum(graph.weights * diff[:, None] * diff[graph.indices])
    n = len(x)
    return float((n / graph.weight_sum) * (num / denom))


def compute_svg_scores(
    st_expr: pd.DataFrame,
    graph: NeighborGraph,
    detect_thresh: float,
    min_var: float,
    method: str,
    norm_method: str,
    n_jobs: int = 1,
) -> Tuple[pd.Series, pd.Series]:
    genes = st_expr.columns
    values = st_expr.to_numpy(dtype=float)

    def worker(col_idx: int) -> float:
        x = values[:, col_idx]
        if np.count_nonzero(x > 0) / len(x) < detect_thresh:
            return 0.0
        if np.var(x) < min_var:
            return 0.0
        if method == "morans":
            return morans_I(x, graph)
        raise ValueError(f"Unsupported svg_method: {method}")

    if n_jobs <= 1:
        svg_raw = [worker(i) for i in range(len(genes))]
    else:
        with ThreadPoolExecutor(max_workers=n_jobs) as ex:
            svg_raw = list(ex.map(worker, range(len(genes))))

    svg_series = pd.Series(svg_raw, index=genes)
    svg_norm = normalize_scores(svg_series, method=norm_method)
    return svg_series, svg_norm


# ----------------------------
# Sensitivity analysis
# ----------------------------


def read_raw_st_expression(raw_path: Path, st_coords_index: pd.Index) -> Tuple[Optional[pd.DataFrame], Dict]:
    if not raw_path.exists():
        return None, {"status": "skipped", "reason": "raw_st_missing", "path": str(raw_path)}
    try:
        dt = pd.read_csv(raw_path, sep=None, engine="python")
    except Exception as e:
        return None, {"status": "skipped", "reason": "read_failed", "error": str(e), "path": str(raw_path)}
    if dt.shape[1] < 2:
        return None, {"status": "skipped", "reason": "raw_shape_invalid", "path": str(raw_path)}

    st_set = set(st_coords_index.astype(str))
    first_col = dt.iloc[:, 0].astype(str)
    other_cols = [str(c) for c in dt.columns[1:]]
    overlap_cols = [c for c in other_cols if c in st_set]
    overlap_rows = [r for r in first_col if r in st_set]

    # Case A: gene x spot (columns are spot ids)
    if len(overlap_cols) > 0:
        expr = pd.DataFrame(dt.iloc[:, 1:].to_numpy(dtype=float), index=first_col, columns=other_cols)
        # 按 st_coords 顺序对齐列
        cols = [c for c in st_coords_index.astype(str) if c in expr.columns]
        expr = expr.loc[:, cols]
        info = {
            "status": "ok",
            "orientation": "gene_spot",
            "n_genes": expr.shape[0],
            "n_spots": expr.shape[1],
        }
        return expr, info

    # Case B: spot x gene（行是 spot）
    if len(overlap_rows) > 0:
        expr = pd.DataFrame(dt.iloc[:, 1:].to_numpy(dtype=float), index=first_col, columns=other_cols)
        expr = expr.loc[expr.index.intersection(st_coords_index.astype(str))]
        expr = expr.T  # 现在变为 gene x spot
        cols = [c for c in st_coords_index.astype(str) if c in expr.columns]
        expr = expr.loc[:, cols]
        info = {
            "status": "ok",
            "orientation": "spot_gene",
            "n_genes": expr.shape[0],
            "n_spots": expr.shape[1],
        }
        return expr, info

    return None, {
        "status": "skipped",
        "reason": "orientation_not_detected",
        "path": str(raw_path),
    }


def svg_filter_sensitivity(
    raw_st_df: pd.DataFrame,
    graph: NeighborGraph,
    filtered_genes: List[str],
    st_coords_index: pd.Index,
    detect_thresh: float,
    min_var: float,
    method: str,
    norm_method: str,
    top_n: int = 500,
    n_jobs: int = 1,
) -> Dict:
    if raw_st_df is None:
        return {
            "status": "skipped",
            "reason": "raw_st_missing",
        }

    if raw_st_df.shape[1] != len(st_coords_index):
        return {
            "status": "skipped",
            "reason": "spot_mismatch",
            "n_spots_raw": raw_st_df.shape[1],
            "n_spots_coords": len(st_coords_index),
        }

    # 确保列顺序与 st_coords 一致
    raw_st_df = raw_st_df.loc[:, st_coords_index.astype(str)]

    # Align to graph spots if possible
    raw_values = raw_st_df.to_numpy(dtype=float)
    # Fast approximate: reuse existing graph size
    def worker(row_idx: int) -> float:
        x = raw_values[row_idx, :]
        if np.count_nonzero(x > 0) / len(x) < detect_thresh:
            return 0.0
        if np.var(x) < min_var:
            return 0.0
        if method == "morans":
            return morans_I(x, graph)
        return 0.0

    genes = raw_st_df.index.to_list()
    if n_jobs <= 1:
        scores = [worker(i) for i in range(len(genes))]
    else:
        with ThreadPoolExecutor(max_workers=n_jobs) as ex:
            scores = list(ex.map(worker, range(len(genes))))

    series = pd.Series(scores, index=genes)
    series_norm = normalize_scores(series, method=norm_method)
    top_genes = series_norm.sort_values(ascending=False).head(top_n).index

    filtered_set = set(filtered_genes)
    lost_genes = [g for g in top_genes if g not in filtered_set]
    lost_ratio = len(lost_genes) / len(top_genes) if len(top_genes) else 0.0

    return {
        "status": "ok",
        "top_n": top_n,
        "n_lost": len(lost_genes),
        "lost_ratio": lost_ratio,
        "lost_genes_sample": lost_genes[:50],
    }




def classify_gene_category(gene: str) -> str:
    """
    简单的基因类别标签，便于 summary 查看。非严格生物注释，仅作展示。
    """
    g = gene.upper()
    epithelial_prefix = ("KRT", "EPCAM", "MUC", "MUC1", "KRTAP")
    immune_prefix = ("CD", "HLA", "IGH", "IGK", "IGL", "MS4A", "LILR")
    stromal_prefix = ("COL", "SPARC", "FN1", "VIM", "ACTA", "TAGLN", "DCN", "THBS", "LUM", "IGFBP", "CAV", "VWF", "KDR")
    if g.startswith(epithelial_prefix):
        return "epithelial_like"
    if g.startswith(immune_prefix):
        return "immune_like"
    if g.startswith(stromal_prefix):
        return "stromal_like"
    return "other"


# ----------------------------
# Main
# ----------------------------


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Stage 2: SVG+HVG dynamic weighting plugin")
    p.add_argument("--sample", default="real_brca", help="sample id")
    p.add_argument("--project_root", default=None, help="override project root")
    p.add_argument("--alpha", type=float, default=None, help="weight for HVG score")
    p.add_argument("--beta", type=float, default=None, help="weight for SVG score")
    p.add_argument("--svg_method", default=None, choices=["morans"], help="SVG metric")
    p.add_argument("--k_neighbors", type=int, default=None, help="k for kNN graph")
    p.add_argument("--svg_min_detect", type=float, default=None, help="min detect fraction")
    p.add_argument("--svg_min_var", type=float, default=None, help="min variance for SVG")
    p.add_argument("--svg_topk", type=int, default=None, help="top-K genes for plugin")
    p.add_argument(
        "--plugin_genes_mode",
        default=None,
        choices=["topk_final_weight", "union_hvg_svg"],
        help="plugin genes selection mode",
    )
    p.add_argument("--hvg_nfeatures", type=int, default=None, help="HVG count used in union_hvg_svg mode")
    p.add_argument("--svg_cost_topk", type=int, default=None, help="top-K SVG genes to upweight in cost")
    p.add_argument("--svg_cost_kappa", type=float, default=None, help="cost weight strength (scaled by beta)")
    p.add_argument("--svg_cost_beta", type=float, default=None, help="extra multiplier for cost weight strength")
    p.add_argument("--svg_cost_weight_max", type=float, default=None, help="max cost_weight for SVG genes")
    p.add_argument(
        "--svg_cost_mode",
        default=None,
        choices=["linear", "exp"],
        help="cost weight mode for SVG genes",
    )
    p.add_argument("--svg_cost_only_up", action="store_true", default=None, help="force cost_weight >= 1")
    p.add_argument("--svg_cost_truth_filter", action="store_true", default=None, help="enable truth-guided SVG cost filter (sim only)")
    p.add_argument("--svg_cost_truth_min_corr", type=float, default=None, help="min positive corr with truth spot×type")
    p.add_argument(
        "--svg_cost_truth_metric",
        default=None,
        choices=["pearson", "spearman"],
        help="correlation metric for truth-guided filter",
    )
    p.add_argument("--svg_cost_truth_path", default=None, help="override sim truth spot×type CSV path")
    p.add_argument("--weight_norm", default=None, choices=["rank", "zscore", "minmax"])
    p.add_argument("--weight_clip", type=float, nargs=2, default=None, help="clip final_weight to [min max]")
    p.add_argument("--n_jobs", type=int, default=None, help="parallel threads")
    p.add_argument("--r_helper", default=None, help="path to R export helper")
    p.add_argument("--skip_export", action="store_true", help="skip calling R helper")
    p.add_argument("--config", default="configs/project_config.yaml", help="project config YAML")
    p.add_argument("--dataset_config", default=None, help="dataset-specific YAML")
    return p.parse_args()


def _parse_weight_clip(value: Optional[object]) -> Optional[List[float]]:
    if value is None:
        return None
    if isinstance(value, str):
        parts = [p for p in value.replace(",", " ").split() if p]
        if len(parts) != 2:
            raise ValueError(f"weight_clip must have 2 values, got: {value!r}")
        return [float(parts[0]), float(parts[1])]
    if isinstance(value, (list, tuple)) and len(value) == 2:
        return [float(value[0]), float(value[1])]
    raise ValueError(f"weight_clip must be a 2-item list or string, got: {type(value)}")


def main():
    args = parse_args()
    project_root = Path(args.project_root) if args.project_root else detect_project_root()
    project_cfg_path = project_root / args.config
    project_cfg = load_yaml(project_cfg_path)
    dataset_cfg_map = (project_cfg or {}).get("dataset_config_map", {}) or {}
    if args.dataset_config:
        dataset_cfg_path = Path(args.dataset_config)
    else:
        mapped_name = dataset_cfg_map.get(args.sample)
        dataset_cfg_path = project_root / "configs" / "datasets" / (mapped_name or f"{args.sample}.yaml")

    defaults = {
        "alpha": 0.5,
        "beta": 0.5,
        "svg_method": "morans",
        "k_neighbors": 6,
        "svg_min_detect": 0.05,
        "svg_min_var": 1e-6,
        "svg_topk": 1000,
        "plugin_genes_mode": "topk_final_weight",
        "hvg_nfeatures": 2000,
        "svg_cost_topk": 30,
        "svg_cost_kappa": 0.4,
        "svg_cost_beta": 1.0,
        "svg_cost_weight_max": 1.5,
        "svg_cost_mode": "linear",
        "svg_cost_only_up": True,
        "svg_cost_truth_filter": False,
        "svg_cost_truth_min_corr": 0.1,
        "svg_cost_truth_metric": "pearson",
        "svg_cost_truth_path": None,
        "weight_norm": "rank",
        "weight_clip": None,
        "n_jobs": 1,
        "raw_st_expr": None,
    }

    cfg = merge_configs(project_cfg, load_yaml(dataset_cfg_path), defaults)

    # CLI overrides
    for key in [
        "alpha",
        "beta",
        "svg_method",
        "k_neighbors",
        "svg_min_detect",
        "svg_min_var",
        "svg_topk",
        "plugin_genes_mode",
        "hvg_nfeatures",
        "svg_cost_topk",
        "svg_cost_kappa",
        "svg_cost_beta",
        "svg_cost_weight_max",
        "svg_cost_mode",
        "svg_cost_only_up",
        "svg_cost_truth_filter",
        "svg_cost_truth_min_corr",
        "svg_cost_truth_metric",
        "svg_cost_truth_path",
        "weight_norm",
        "weight_clip",
        "n_jobs",
    ]:
        v = getattr(args, key)
        if v is not None:
            cfg[key] = v

    r_helper = Path(args.r_helper) if args.r_helper else project_root / "r_scripts" / "export_stage1_to_csv.R"

    paths = ensure_exports(project_root, args.sample, r_helper=r_helper, skip_export=args.skip_export)
    export_dir = paths["export_dir"]
    sc_expr = read_expression_table(paths["sc_expr"], index_col="cell_id")
    st_expr = read_expression_table(paths["st_expr"], index_col="spot_id")
    st_coords = read_expression_table(paths["st_coords"], index_col="spot_id")

    # Align spot ordering
    st_expr = st_expr.loc[st_coords.index]

    # Common genes alignment
    common_genes_path = paths["common_genes"]
    if common_genes_path.exists():
        common_genes = [g.strip() for g in common_genes_path.read_text().splitlines() if g.strip()]
        sc_expr = sc_expr[common_genes]
        st_expr = st_expr[common_genes]
    else:
        common_genes = list(sorted(set(sc_expr.columns) & set(st_expr.columns)))
        sc_expr = sc_expr[common_genes]
        st_expr = st_expr[common_genes]

    # HVG scores
    hvg_score, hvg_norm = compute_hvg_scores(sc_expr, norm_method=cfg["weight_norm"])
    is_sc_hvg = pd.Series(False, index=common_genes)
    hvg_file = paths["hvg_genes"]
    if hvg_file.exists():
        hvg_list = [g.strip() for g in hvg_file.read_text().splitlines() if g.strip()]
        is_sc_hvg.loc[is_sc_hvg.index.isin(hvg_list)] = True

    # SVG scores
    graph = build_knn_graph(st_coords, k_neighbors=cfg["k_neighbors"])
    svg_score, svg_norm = compute_svg_scores(
        st_expr,
        graph,
        detect_thresh=cfg["svg_min_detect"],
        min_var=cfg["svg_min_var"],
        method=cfg["svg_method"],
        norm_method=cfg["weight_norm"],
        n_jobs=int(cfg["n_jobs"]),
    )
    svg_cost_norm = normalize_scores(svg_score, method="minmax")

    # Optional truth-guided filter (sim only)
    truth_filter_enabled = bool(cfg.get("svg_cost_truth_filter", False))
    truth_min_corr = float(cfg.get("svg_cost_truth_min_corr") or 0.0)
    truth_metric = str(cfg.get("svg_cost_truth_metric") or "pearson").lower()
    truth_info: Dict[str, object] = {"status": "disabled"}
    truth_corr_max = pd.Series(np.nan, index=common_genes)
    truth_best_type = pd.Series([None] * len(common_genes), index=common_genes)
    if truth_filter_enabled:
        truth_path = _resolve_truth_spot_type_path(project_root, args.sample, cfg)
        if truth_path and truth_path.exists():
            truth_df = pd.read_csv(truth_path)
            truth_corr_max, truth_best_type, truth_info = _compute_truth_corr(
                st_expr=st_expr,
                truth_df=truth_df,
                metric=truth_metric,
                eps=1e-8,
            )
            truth_info["path"] = str(truth_path)
            truth_info["metric"] = truth_metric
            truth_info["min_corr"] = truth_min_corr
        else:
            truth_info = {"status": "missing_truth", "path": str(truth_path) if truth_path else None}
            truth_filter_enabled = False

    # Final weights
    final_weight = cfg["alpha"] * hvg_norm + cfg["beta"] * svg_norm
    weight_clip = _parse_weight_clip(cfg.get("weight_clip"))
    if weight_clip is not None:
        low, high = weight_clip
        if low > high:
            raise ValueError(f"weight_clip must be [min max], got: {weight_clip}")
        final_weight = final_weight.clip(lower=low, upper=high)
        cfg["weight_clip"] = [float(low), float(high)]
    genes_df = pd.DataFrame(
        {
            "gene": common_genes,
            "hvg_score": hvg_score.values,
            "svg_score": svg_score.values,
            "hvg_norm": hvg_norm.values,
            "svg_norm": svg_norm.values,
            "svg_cost_norm": svg_cost_norm.values,
            "final_weight": final_weight.values,
            "is_sc_hvg": is_sc_hvg.values,
        }
    )
    genes_df["truth_corr_max"] = truth_corr_max.values
    genes_df["truth_corr_best_type"] = truth_best_type.values
    if truth_filter_enabled and truth_info.get("status") == "ok":
        truth_pass = truth_corr_max >= truth_min_corr
    else:
        truth_pass = pd.Series(False, index=truth_corr_max.index)
    genes_df["truth_pass"] = truth_pass.values.astype(int)
    genes_df["selected_in_plugin"] = 0

    # Plugin genes
    mode = str(cfg.get("plugin_genes_mode") or "topk_final_weight")
    if mode == "topk_final_weight":
        topk = min(cfg["svg_topk"], len(genes_df))
        plugin_genes = genes_df.sort_values("final_weight", ascending=False).head(topk)["gene"].tolist()
    elif mode == "union_hvg_svg":
        svg_topk = min(cfg["svg_topk"], len(genes_df))
        svg_genes = genes_df.sort_values("svg_norm", ascending=False).head(svg_topk)["gene"].tolist()
        hvg_n = int(cfg.get("hvg_nfeatures") or 0)
        hvg_candidates = genes_df[genes_df["is_sc_hvg"] == True]
        if hvg_candidates.empty:
            hvg_candidates = genes_df
        if hvg_n <= 0:
            hvg_n = len(hvg_candidates)
        hvg_n = min(hvg_n, len(genes_df))
        if len(hvg_candidates) >= hvg_n:
            hvg_genes = (
                hvg_candidates.sort_values("hvg_norm", ascending=False).head(hvg_n)["gene"].tolist()
            )
        else:
            hvg_genes = hvg_candidates.sort_values("hvg_norm", ascending=False)["gene"].tolist()
            needed = hvg_n - len(hvg_genes)
            if needed > 0:
                extra = (
                    genes_df[~genes_df["gene"].isin(hvg_genes)]
                    .sort_values("hvg_norm", ascending=False)
                    .head(needed)["gene"]
                    .tolist()
                )
                hvg_genes += extra
        union_set = set(hvg_genes) | set(svg_genes)
        plugin_genes = (
            genes_df[genes_df["gene"].isin(union_set)]
            .sort_values("final_weight", ascending=False)["gene"]
            .tolist()
        )
    else:
        raise ValueError(f"Unknown plugin_genes_mode: {mode}")
    genes_df.loc[genes_df["gene"].isin(plugin_genes), "selected_in_plugin"] = 1

    # Cost weights: default 1, upweight only top SVG genes.
    svg_cost_topk = int(cfg.get("svg_cost_topk") or 0)
    svg_cost_topk = min(max(svg_cost_topk, 0), len(genes_df))
    svg_cost_candidates = svg_cost_norm
    truth_filter_applied = False
    truth_filter_empty = False
    if truth_filter_enabled and truth_info.get("status") == "ok":
        truth_filter_applied = True
        svg_cost_candidates = svg_cost_candidates[truth_pass]
        if svg_cost_candidates.empty:
            truth_filter_empty = True
    svg_cost_genes = (
        svg_cost_candidates.sort_values(ascending=False).head(svg_cost_topk).index.tolist()
        if svg_cost_topk > 0 and len(svg_cost_candidates) > 0
        else []
    )
    kappa = float(cfg.get("svg_cost_kappa") or 0.0)
    cost_beta = float(cfg.get("svg_cost_beta") or 1.0)
    kappa_eff = kappa * cost_beta
    if bool(cfg.get("svg_cost_only_up", True)):
        kappa_eff = max(kappa_eff, 0.0)
    if str(cfg.get("svg_cost_mode") or "linear") == "exp":
        weight_vals = np.exp(kappa_eff * svg_cost_norm)
    else:
        weight_vals = 1.0 + kappa_eff * svg_cost_norm
    weight_vals = weight_vals.clip(lower=1.0)
    weight_max = float(cfg.get("svg_cost_weight_max") or 1.0)
    if weight_max < 1.0:
        raise ValueError(f"svg_cost_weight_max must be >= 1, got {weight_max}")
    weight_vals = weight_vals.clip(upper=weight_max)
    cost_weight = pd.Series(1.0, index=svg_cost_norm.index, dtype=float)
    if svg_cost_genes:
        cost_weight.loc[svg_cost_genes] = weight_vals.loc[svg_cost_genes]
    genes_df["cost_weight"] = cost_weight.values
    genes_df["selected_svg_cost"] = genes_df["gene"].isin(svg_cost_genes).astype(int)

    # Output dirs
    processed_out = project_root / "data" / "processed" / args.sample / "stage2_svg_plugin"
    result_out = project_root / "result" / args.sample / "stage2_svg_plugin"
    processed_out.mkdir(parents=True, exist_ok=True)
    result_out.mkdir(parents=True, exist_ok=True)

    # Save gene weights and plugin list
    genes_df.to_csv(processed_out / "gene_weights.csv", index=False)
    (processed_out / "plugin_genes.txt").write_text("\n".join(plugin_genes), encoding="utf-8")

    # SVG filter sensitivity (optional for simulated samples)
    raw_st_df = None
    raw_info: Dict = {"status": "skipped", "reason": "raw_st_expr_not_set"}
    raw_st_path = cfg.get("raw_st_expr")
    if raw_st_path is None:
        default_raw = project_root / "data" / "raw" / args.sample / f"{args.sample}_STdata_GEP.txt"
        if default_raw.exists():
            raw_st_path = default_raw
        else:
            try:
                rel_cfg = dataset_cfg_path.relative_to(project_root)
                cfg_hint_path = rel_cfg.as_posix()
            except Exception:
                cfg_hint_path = str(dataset_cfg_path)
            raw_info = {
                "status": "skipped",
                "reason": "raw_st_missing",
                "default_path": str(default_raw),
                "hint": f"若需要敏感性分析，请在 {cfg_hint_path} 的 stage2.raw_st_expr 指定原始 ST 表达矩阵路径",
            }
            print(f"[Stage2] Raw ST sensitivity skipped: {raw_info}")
            raw_st_path = None
    if raw_st_path is not None:
        raw_st_path = Path(raw_st_path)
        if not raw_st_path.is_absolute():
            raw_st_path = project_root / raw_st_path
        raw_st_df, raw_info = read_raw_st_expression(raw_st_path, st_coords_index=st_coords.index)
        if raw_info.get("status") != "ok":
            print(f"[Stage2] Raw ST sensitivity skipped: {raw_info}")
            raw_st_df = None

    filtered_genes = list(st_expr.columns)
    if raw_st_df is None:
        sensitivity = {"status": "skipped", "raw_info": raw_info}
    else:
        sensitivity = svg_filter_sensitivity(
            raw_st_df,
            graph,
            filtered_genes,
            st_coords_index=st_coords.index.astype(str),
            detect_thresh=cfg["svg_min_detect"],
            min_var=cfg["svg_min_var"],
            method=cfg["svg_method"],
            norm_method=cfg["weight_norm"],
            top_n=500,
            n_jobs=int(cfg["n_jobs"]),
        )
    with (result_out / "svg_filter_sensitivity.json").open("w", encoding="utf-8") as f:
        json.dump(sensitivity, f, ensure_ascii=False, indent=2)

    # Summary / params
    top_plugin_info = [
        {
            "gene": g,
            "final_weight": float(genes_df.set_index("gene").loc[g, "final_weight"]),
            "category": classify_gene_category(g),
        }
        for g in plugin_genes[:20]
    ]
    cw = genes_df["cost_weight"].to_numpy(dtype=float)
    truth_keep = int(truth_pass.sum()) if truth_filter_enabled and truth_info.get("status") == "ok" else None
    truth_keep_ratio = (truth_keep / len(genes_df)) if truth_keep is not None and len(genes_df) else None
    cost_weight_stats = {
        "min": float(np.min(cw)) if len(cw) else None,
        "p50": float(np.percentile(cw, 50)) if len(cw) else None,
        "p95": float(np.percentile(cw, 95)) if len(cw) else None,
        "max": float(np.max(cw)) if len(cw) else None,
        "n_weighted": int((cw > 1.0).sum()) if len(cw) else 0,
        "weighted_ratio": float((cw > 1.0).mean()) if len(cw) else 0.0,
        "svg_cost_topk": int(svg_cost_topk),
        "truth_filter_enabled": bool(truth_filter_enabled),
        "truth_filter_applied": bool(truth_filter_applied),
        "truth_filter_empty": bool(truth_filter_empty),
        "truth_filter_keep": truth_keep,
        "truth_filter_keep_ratio": truth_keep_ratio,
        "truth_filter_status": truth_info.get("status"),
        "truth_filter_path": truth_info.get("path"),
        "truth_filter_min_corr": truth_min_corr,
    }
    summary = {
        "sample": args.sample,
        "n_genes_total": len(genes_df),
        "n_genes_common": len(filtered_genes),
        "n_genes_plugin": len(plugin_genes),
        "top_plugin_genes": top_plugin_info,
        "sensitivity": sensitivity,
        "cost_weight_stats": cost_weight_stats,
        "truth_filter_info": truth_info,
        "params": cfg,
    }
    with (result_out / "stage2_summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    with (result_out / "stage2_params.json").open("w", encoding="utf-8") as f:
        json.dump(cfg, f, ensure_ascii=False, indent=2)

    print(f"[Stage2] gene_weights.csv -> {processed_out / 'gene_weights.csv'}")
    print(f"[Stage2] plugin_genes.txt -> {processed_out / 'plugin_genes.txt'}")
    print(f"[Stage2] summary -> {result_out / 'stage2_summary.json'}")


if __name__ == "__main__":
    main()
