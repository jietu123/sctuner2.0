"""CytoSPACE backend: baseline vs plus (SVG+type plugins) with optional refine."""
from __future__ import annotations

import hashlib
import json
import os
import sys
import tempfile
import time
from pathlib import Path
from typing import Dict, Any, List, Tuple, Optional

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse import coo_matrix, csr_matrix
from scipy.spatial.distance import cdist
from scipy.spatial import cKDTree

# Stage4-plus（CytoSPACEBackend.run_plus）关键路径唯一入口（防回归）：
# - soft matrix：_parse_assigned_locations -> _assignments_to_matrix
# - SVG refine：_refine_spot_cell_matrix_svg
# - type prior refine：_type_prior_refine
# - harden（容量可行化）：_harden_assignment_quota_matching
# - outputs：_build_outputs

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "src"
for p in (ROOT, SRC):
    if str(p) not in sys.path:
        sys.path.append(str(p))

from stages.backends.base_backend import MappingBackend
try:
    from cytospace.cytospace import main_cytospace
except ModuleNotFoundError:
    # Fallback to local editable source (external/cytospace) if not installed in current env.
    ext = ROOT / "external" / "cytospace"
    if ext.exists() and str(ext) not in sys.path:
        sys.path.append(str(ext))
    from cytospace.cytospace import main_cytospace


def _assert_unique_index(df: pd.DataFrame, name: str):
    if not df.index.is_unique:
        dup = df.index[df.index.duplicated()].unique()[:5]
        raise ValueError(f"{name} index 存在重复: {list(dup)}")


def _ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)


def _module_fingerprint() -> Dict[str, Optional[str]]:
    p = Path(__file__).resolve()
    try:
        sha1 = hashlib.sha1(p.read_bytes()).hexdigest().lower()
    except Exception:
        sha1 = None
    return {"module_file": str(p), "module_sha1": sha1}


def _sha1_file(path: Path) -> Optional[str]:
    try:
        return hashlib.sha1(path.read_bytes()).hexdigest().lower()
    except Exception:
        return None


_DEPRECATED_KEYS: Dict[str, str] = {
    "capacity_normalize_mode": "废弃/从未完整实现：请改用 cells_per_spot_* 族配置（cells_per_spot_source/clip/rounding）",
    "capacity_normalize_factor": "废弃：请改用 umi_to_cell_norm 或 default_cells_per_spot（配合 cells_per_spot_source）",
    "capacity_norm_mode": "废弃：同 capacity_normalize_mode",
}

_TRACE_KEYS: set[str] = {
    "sample",
    "backend",
    "config_id",
    "project_root",
    "project_config_path",
    "dataset_config_path",
    "runner_file",
    "runner_sha1",
}

_USED_KEYS: set[str] = {
    # run identity / outputs
    "mode",
    "run_id",
    "variant",
    # cyto/cost/refine
    "seed",
    "solver_method",
    "distance_metric",
    "eps",
    "svg_refine_lambda",
    "svg_refine_k",
    "svg_refine_batch_size",
    "lambda_prior",
    "prior_candidate_topk",
    "prior_candidate_weight",
    "abstain_unknown_sc_only",
    "type_prior_apply_refine",
    "type_prior_apply_harden",
    "knn_metric",
    "knn_block_size",
    "knn_max_dense_n",
    "harden_topk",
    "prior_ablation_enabled",
    "min_gene_overlap_ratio",
    "max_cells_missing_type_prior_ratio",
    "min_prior_row_nonzero_ratio",
    "post_assign_adjust",
    "cost_expr_norm",
    "cost_expr_norm_eps",
    "cost_scale_source",
    "cost_scale_metric",
    "cost_scale_eps",
    "cost_scale_floor",
    "cost_scale_clip",
    "cost_scale_report_quantiles",
    "cost_scale_beta",
    "cost_scale_clip_ratio_soft",
    "cost_scale_clip_ratio_max",
    "cost_scale_downweight_only",
    "spot_weight_mode",
    "spot_weight_topk",
    "spot_weight_k",
    "spot_weight_kappa",
    "spot_weight_weight_max",
    "spot_weight_scale",
    "spot_weight_only_up",
    "spot_weight_truth_filter",
    "spot_specific_basin",
    "spot_specific_threshold_kappa",
    "default_config_id",
    # cells_per_spot / capacity
    "cells_per_spot_source",
    "cells_per_spot_clip_min",
    "cells_per_spot_clip_max",
    "cells_per_spot_rounding",
    "umi_to_cell_norm",
    "default_cells_per_spot",
}

_CONFIG_EFFECTIVE_KEYS: list[str] = [
    "mode",
    "run_id",
    "variant",
    "seed",
    "solver_method",
    "distance_metric",
    "eps",
    "svg_refine_lambda",
    "svg_refine_k",
    "svg_refine_batch_size",
    "lambda_prior",
    "prior_candidate_topk",
    "prior_candidate_weight",
    "abstain_unknown_sc_only",
    "type_prior_apply_refine",
    "type_prior_apply_harden",
    "prior_ablation_enabled",
    "knn_metric",
    "knn_block_size",
    "knn_max_dense_n",
    "harden_topk",
    "min_gene_overlap_ratio",
    "max_cells_missing_type_prior_ratio",
    "min_prior_row_nonzero_ratio",
    "post_assign_adjust",
    "cost_expr_norm",
    "cost_expr_norm_eps",
    "cost_scale_source",
    "cost_scale_metric",
    "cost_scale_eps",
    "cost_scale_floor",
    "cost_scale_clip",
    "cost_scale_report_quantiles",
    "cost_scale_beta",
    "cost_scale_clip_ratio_soft",
    "cost_scale_clip_ratio_max",
    "cost_scale_downweight_only",
    "spot_weight_mode",
    "spot_weight_topk",
    "spot_weight_k",
    "spot_weight_kappa",
    "spot_weight_weight_max",
    "spot_weight_scale",
    "spot_weight_only_up",
    "spot_weight_truth_filter",
    "spot_specific_basin",
    "spot_specific_threshold_kappa",
    "default_config_id",
    "cells_per_spot_source",
    "umi_to_cell_norm",
    "default_cells_per_spot",
    "cells_per_spot_rounding",
    "cells_per_spot_clip_min",
    "cells_per_spot_clip_max",
    "strict_config",
]


def _validate_and_resolve_config(config: Dict[str, Any], *, context: str) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    Risk#7: 配置严格校验/归一化，避免“配了但没生效/静默忽略”。
    返回 (cfg_effective, config_validation)。
    """
    if config is None:
        config = {}
    if not isinstance(config, dict):
        raise TypeError(f"[{context}] config must be a dict, got {type(config)!r}")

    cfg = dict(config)

    def _as_bool(v: Any, key: str) -> bool:
        if isinstance(v, bool):
            return v
        if isinstance(v, (int, np.integer)):
            return bool(int(v))
        if isinstance(v, str):
            s = v.strip().lower()
            if s in ("1", "true", "yes", "y", "on"):
                return True
            if s in ("0", "false", "no", "n", "off"):
                return False
        raise ValueError(f"[{context}] invalid bool {key}={v!r}")

    def _as_int(v: Any, key: str, *, min_v: Optional[int] = None) -> int:
        if isinstance(v, bool):
            raise ValueError(f"[{context}] invalid int {key}={v!r}")
        try:
            iv = int(v)
        except Exception as e:
            raise ValueError(f"[{context}] invalid int {key}={v!r}: {e}") from e
        if min_v is not None and iv < min_v:
            raise ValueError(f"[{context}] invalid int {key}={iv} < {min_v}")
        return iv

    def _as_float(v: Any, key: str, *, min_v: Optional[float] = None, max_v: Optional[float] = None) -> float:
        if isinstance(v, bool):
            raise ValueError(f"[{context}] invalid float {key}={v!r}")
        try:
            fv = float(v)
        except Exception as e:
            raise ValueError(f"[{context}] invalid float {key}={v!r}: {e}") from e
        if min_v is not None and fv < min_v:
            raise ValueError(f"[{context}] invalid float {key}={fv} < {min_v}")
        if max_v is not None and fv > max_v:
            raise ValueError(f"[{context}] invalid float {key}={fv} > {max_v}")
        if not np.isfinite(fv):
            raise ValueError(f"[{context}] invalid float {key}={fv} (non-finite)")
        return fv

    strict = _as_bool(cfg.get("strict_config", True), "strict_config")
    cfg["strict_config"] = strict

    deprecated_found = sorted(set(cfg.keys()) & set(_DEPRECATED_KEYS.keys()))
    if deprecated_found:
        details = "; ".join([f"{k}: {_DEPRECATED_KEYS[k]}" for k in deprecated_found])
        raise ValueError(f"[{context}] deprecated config keys found: {deprecated_found}. {details}")

    allowed_keys = set(_TRACE_KEYS) | set(_USED_KEYS) | {"strict_config"}
    unknown_keys = sorted([k for k in cfg.keys() if k not in allowed_keys])
    if unknown_keys and strict:
        raise ValueError(f"[{context}] unknown config keys (strict_config=True): {unknown_keys}")

    # Enums / ranges
    if "solver_method" in cfg and cfg["solver_method"] is not None:
        solver = str(cfg["solver_method"])
        allowed_solver = {"lap_CSPR", "lapjv", "lapjv_compat"}
        if solver not in allowed_solver:
            raise ValueError(f"[{context}] invalid solver_method={solver!r}; allowed={sorted(allowed_solver)}")
        cfg["solver_method"] = solver

    if "distance_metric" in cfg and cfg["distance_metric"] is not None:
        dm = str(cfg["distance_metric"])
        allowed_dm = {"Pearson_correlation", "Spearman_correlation", "Euclidean"}
        if dm not in allowed_dm:
            raise ValueError(f"[{context}] invalid distance_metric={dm!r}; allowed={sorted(allowed_dm)}")
        cfg["distance_metric"] = dm

    if "cost_expr_norm" in cfg and cfg["cost_expr_norm"] is not None:
        norm = str(cfg["cost_expr_norm"]).strip().lower()
        if norm in ("", "none"):
            norm = "none"
        allowed_norm = {"none", "st_center", "st_zscore", "gene_scale"}
        if norm not in allowed_norm:
            raise ValueError(f"[{context}] invalid cost_expr_norm={norm!r}; allowed={sorted(allowed_norm)}")
        cfg["cost_expr_norm"] = norm
    if "cost_expr_norm_eps" in cfg and cfg["cost_expr_norm_eps"] is not None:
        cfg["cost_expr_norm_eps"] = _as_float(cfg["cost_expr_norm_eps"], "cost_expr_norm_eps", min_v=1e-12)
    if "cost_scale_source" in cfg and cfg["cost_scale_source"] is not None:
        src = str(cfg["cost_scale_source"]).strip().lower()
        allowed_src = {"st", "joint"}
        if src not in allowed_src:
            raise ValueError(f"[{context}] invalid cost_scale_source={src!r}; allowed={sorted(allowed_src)}")
        cfg["cost_scale_source"] = src
    if "cost_scale_metric" in cfg and cfg["cost_scale_metric"] is not None:
        metric = str(cfg["cost_scale_metric"]).strip().lower()
        allowed_metric = {"std", "mean", "p95"}
        if metric not in allowed_metric:
            raise ValueError(f"[{context}] invalid cost_scale_metric={metric!r}; allowed={sorted(allowed_metric)}")
        cfg["cost_scale_metric"] = metric
    if "cost_scale_eps" in cfg and cfg["cost_scale_eps"] is not None:
        cfg["cost_scale_eps"] = _as_float(cfg["cost_scale_eps"], "cost_scale_eps", min_v=1e-12)
    if "cost_scale_floor" in cfg and cfg["cost_scale_floor"] is not None:
        cfg["cost_scale_floor"] = _as_float(cfg["cost_scale_floor"], "cost_scale_floor", min_v=0.0)
    if "cost_scale_clip" in cfg and cfg["cost_scale_clip"] is not None:
        clip = cfg["cost_scale_clip"]
        if isinstance(clip, str):
            parts = [p for p in clip.replace(",", " ").split() if p]
            if len(parts) != 2:
                raise ValueError(f"[{context}] cost_scale_clip must have 2 values, got {clip!r}")
            clip = [float(parts[0]), float(parts[1])]
        if not isinstance(clip, (list, tuple)) or len(clip) != 2:
            raise ValueError(f"[{context}] cost_scale_clip must be 2-item list, got {clip!r}")
        clip_low, clip_high = float(clip[0]), float(clip[1])
        if clip_low <= 0 or clip_high <= 0 or clip_low > clip_high:
            raise ValueError(f"[{context}] invalid cost_scale_clip={clip!r}")
        cfg["cost_scale_clip"] = [clip_low, clip_high]
    if "cost_scale_report_quantiles" in cfg and cfg["cost_scale_report_quantiles"] is not None:
        q = cfg["cost_scale_report_quantiles"]
        if isinstance(q, str):
            parts = [p for p in q.replace(",", " ").split() if p]
            q = [float(x) for x in parts]
        if not isinstance(q, (list, tuple)) or len(q) == 0:
            raise ValueError(f"[{context}] cost_scale_report_quantiles must be list, got {q!r}")
        cfg["cost_scale_report_quantiles"] = [float(x) for x in q]
    if "cost_scale_beta" in cfg and cfg["cost_scale_beta"] is not None:
        cfg["cost_scale_beta"] = _as_float(cfg["cost_scale_beta"], "cost_scale_beta", min_v=0.0, max_v=1.0)
    if "cost_scale_clip_ratio_soft" in cfg and cfg["cost_scale_clip_ratio_soft"] is not None:
        cfg["cost_scale_clip_ratio_soft"] = _as_float(
            cfg["cost_scale_clip_ratio_soft"], "cost_scale_clip_ratio_soft", min_v=0.0, max_v=1.0
        )
    if "cost_scale_clip_ratio_max" in cfg and cfg["cost_scale_clip_ratio_max"] is not None:
        cfg["cost_scale_clip_ratio_max"] = _as_float(
            cfg["cost_scale_clip_ratio_max"], "cost_scale_clip_ratio_max", min_v=0.0, max_v=1.0
        )
    if "cost_scale_downweight_only" in cfg and cfg["cost_scale_downweight_only"] is not None:
        cfg["cost_scale_downweight_only"] = _as_bool(cfg["cost_scale_downweight_only"], "cost_scale_downweight_only")

    if "knn_metric" in cfg and cfg["knn_metric"] is not None:
        km = str(cfg["knn_metric"]).lower()
        allowed_km = {"euclidean", "cosine"}
        if km not in allowed_km:
            raise ValueError(f"[{context}] invalid knn_metric={km!r}; allowed={sorted(allowed_km)}")
        cfg["knn_metric"] = km

    if "spot_weight_mode" in cfg and cfg["spot_weight_mode"] is not None:
        mode = str(cfg["spot_weight_mode"]).strip().lower()
        if mode in ("", "none"):
            mode = "none"
        allowed_mode = {"none", "local_zscore"}
        if mode not in allowed_mode:
            raise ValueError(f"[{context}] invalid spot_weight_mode={mode!r}; allowed={sorted(allowed_mode)}")
        cfg["spot_weight_mode"] = mode
    if "spot_weight_scale" in cfg and cfg["spot_weight_scale"] is not None:
        scale = str(cfg["spot_weight_scale"]).strip().lower()
        allowed_scale = {"linear", "exp"}
        if scale not in allowed_scale:
            raise ValueError(f"[{context}] invalid spot_weight_scale={scale!r}; allowed={sorted(allowed_scale)}")
        cfg["spot_weight_scale"] = scale

    if "cells_per_spot_source" in cfg and cfg["cells_per_spot_source"] is not None:
        src0 = str(cfg["cells_per_spot_source"]).strip()
        src_map = {
            "auto": "auto",
            "spot_cell_counts": "spot_cell_counts",
            "umi_total": "UMI_total",
            "umi_total ": "UMI_total",
            "umi_total\n": "UMI_total",
            "UMI_total": "UMI_total",
            "uniform": "uniform",
        }
        src_norm = src0 if src0 in src_map else src0.lower()
        if src_norm not in src_map:
            raise ValueError(
                f"[{context}] invalid cells_per_spot_source={src0!r}; allowed={['auto','spot_cell_counts','UMI_total','uniform']}"
            )
        cfg["cells_per_spot_source"] = src_map[src_norm]

    if "cells_per_spot_rounding" in cfg and cfg["cells_per_spot_rounding"] is not None:
        r = str(cfg["cells_per_spot_rounding"]).strip().lower()
        allowed_r = {"round", "floor", "ceil"}
        if r not in allowed_r:
            raise ValueError(f"[{context}] invalid cells_per_spot_rounding={r!r}; allowed={sorted(allowed_r)}")
        cfg["cells_per_spot_rounding"] = r

    if "seed" in cfg and cfg["seed"] is not None:
        cfg["seed"] = _as_int(cfg["seed"], "seed")
    if "eps" in cfg and cfg["eps"] is not None:
        cfg["eps"] = _as_float(cfg["eps"], "eps", min_v=1e-16)
    if "svg_refine_lambda" in cfg and cfg["svg_refine_lambda"] is not None:
        cfg["svg_refine_lambda"] = _as_float(cfg["svg_refine_lambda"], "svg_refine_lambda", min_v=0.0)
    if "svg_refine_k" in cfg and cfg["svg_refine_k"] is not None:
        cfg["svg_refine_k"] = _as_int(cfg["svg_refine_k"], "svg_refine_k", min_v=1)
    if "svg_refine_batch_size" in cfg:
        if cfg["svg_refine_batch_size"] in (None, ""):
            cfg["svg_refine_batch_size"] = None
        else:
            cfg["svg_refine_batch_size"] = _as_int(cfg["svg_refine_batch_size"], "svg_refine_batch_size", min_v=1)

    if "lambda_prior" in cfg and cfg["lambda_prior"] is not None:
        cfg["lambda_prior"] = _as_float(cfg["lambda_prior"], "lambda_prior", min_v=0.0)
    if "prior_candidate_topk" in cfg and cfg["prior_candidate_topk"] is not None:
        cfg["prior_candidate_topk"] = _as_int(cfg["prior_candidate_topk"], "prior_candidate_topk", min_v=0)
    if "prior_candidate_weight" in cfg and cfg["prior_candidate_weight"] is not None:
        cfg["prior_candidate_weight"] = _as_float(cfg["prior_candidate_weight"], "prior_candidate_weight", min_v=0.0)
    if "abstain_unknown_sc_only" in cfg and cfg["abstain_unknown_sc_only"] is not None:
        cfg["abstain_unknown_sc_only"] = _as_bool(cfg["abstain_unknown_sc_only"], "abstain_unknown_sc_only")
    if "harden_topk" in cfg and cfg["harden_topk"] is not None:
        cfg["harden_topk"] = _as_int(cfg["harden_topk"], "harden_topk", min_v=1)
    if "spot_weight_topk" in cfg and cfg["spot_weight_topk"] is not None:
        cfg["spot_weight_topk"] = _as_int(cfg["spot_weight_topk"], "spot_weight_topk", min_v=0)
    if "spot_weight_k" in cfg and cfg["spot_weight_k"] is not None:
        cfg["spot_weight_k"] = _as_int(cfg["spot_weight_k"], "spot_weight_k", min_v=1)
    if "spot_weight_kappa" in cfg and cfg["spot_weight_kappa"] is not None:
        cfg["spot_weight_kappa"] = _as_float(cfg["spot_weight_kappa"], "spot_weight_kappa", min_v=0.0)
    if "spot_weight_weight_max" in cfg and cfg["spot_weight_weight_max"] is not None:
        cfg["spot_weight_weight_max"] = _as_float(cfg["spot_weight_weight_max"], "spot_weight_weight_max", min_v=1.0)
    if "spot_weight_only_up" in cfg and cfg["spot_weight_only_up"] is not None:
        cfg["spot_weight_only_up"] = _as_bool(cfg["spot_weight_only_up"], "spot_weight_only_up")
    if "spot_weight_truth_filter" in cfg and cfg["spot_weight_truth_filter"] is not None:
        cfg["spot_weight_truth_filter"] = _as_bool(cfg["spot_weight_truth_filter"], "spot_weight_truth_filter")

    for k in ("knn_block_size", "knn_max_dense_n"):
        if k in cfg and cfg[k] is not None:
            cfg[k] = _as_int(cfg[k], k, min_v=1)

    for k in ("min_gene_overlap_ratio", "max_cells_missing_type_prior_ratio", "min_prior_row_nonzero_ratio"):
        if k in cfg and cfg[k] is not None:
            cfg[k] = _as_float(cfg[k], k, min_v=0.0, max_v=1.0)

    for k in ("type_prior_apply_refine", "type_prior_apply_harden", "prior_ablation_enabled"):
        if k in cfg and cfg[k] is not None:
            cfg[k] = _as_bool(cfg[k], k)
    if "post_assign_adjust" in cfg and cfg["post_assign_adjust"] is not None:
        cfg["post_assign_adjust"] = _as_bool(cfg["post_assign_adjust"], "post_assign_adjust")

    if "cells_per_spot_clip_min" in cfg and cfg["cells_per_spot_clip_min"] is not None:
        cfg["cells_per_spot_clip_min"] = _as_int(cfg["cells_per_spot_clip_min"], "cells_per_spot_clip_min", min_v=0)
    if "cells_per_spot_clip_max" in cfg:
        if cfg["cells_per_spot_clip_max"] in (None, ""):
            cfg["cells_per_spot_clip_max"] = None
        else:
            cfg["cells_per_spot_clip_max"] = _as_int(cfg["cells_per_spot_clip_max"], "cells_per_spot_clip_max", min_v=0)
    if cfg.get("cells_per_spot_clip_max") is not None and cfg.get("cells_per_spot_clip_min") is not None:
        if int(cfg["cells_per_spot_clip_max"]) < int(cfg["cells_per_spot_clip_min"]):
            raise ValueError(
                f"[{context}] cells_per_spot_clip_max {cfg['cells_per_spot_clip_max']} < cells_per_spot_clip_min {cfg['cells_per_spot_clip_min']}"
            )

    if "default_cells_per_spot" in cfg and cfg["default_cells_per_spot"] is not None:
        cfg["default_cells_per_spot"] = _as_float(cfg["default_cells_per_spot"], "default_cells_per_spot", min_v=0.0)
    if "umi_to_cell_norm" in cfg and cfg["umi_to_cell_norm"] is not None:
        norm = cfg["umi_to_cell_norm"]
        if isinstance(norm, str) and norm.strip().lower() == "median":
            cfg["umi_to_cell_norm"] = "median"
        else:
            cfg["umi_to_cell_norm"] = _as_float(norm, "umi_to_cell_norm", min_v=1e-12)

    config_validation = {
        "schema_version": "cytospace_backend_config_v1",
        "strict_config": strict,
        "unknown_keys": unknown_keys,
        "deprecated_keys_found": deprecated_found,
    }
    return cfg, config_validation


def _get_coord_xy(df: pd.DataFrame) -> pd.DataFrame:
    if {"row", "col"}.issubset(df.columns):
        return df[["row", "col"]]
    if {"x", "y"}.issubset(df.columns):
        return df[["x", "y"]].rename(columns={"x": "row", "y": "col"})
    raise KeyError("st_coordinates 需要包含 row/col 或 x/y 坐标列")


def _load_stage1(stage1_dir: Path):
    exp_dir = stage1_dir / "exported"
    sc_expr = pd.read_csv(exp_dir / "sc_expression_normalized.csv", index_col=0)
    st_expr = pd.read_csv(exp_dir / "st_expression_normalized.csv", index_col=0)
    sc_meta = pd.read_csv(exp_dir / "sc_metadata.csv")
    st_coords = pd.read_csv(exp_dir / "st_coordinates.csv", index_col=0)
    # 确保 cell_id 对齐
    if "cell_id" in sc_meta.columns:
        sc_meta = sc_meta.set_index("cell_id")
    else:
        raise KeyError("sc_metadata.csv 需要包含 cell_id 列用于与 sc_expression_normalized 对齐")
    missing = set(sc_expr.index) - set(sc_meta.index)
    if missing:
        raise ValueError(f"sc_metadata 缺少以下 cell_id: {list(missing)[:5]} ...")
    sc_meta = sc_meta.loc[sc_expr.index]
    # 类型列兼容 + 护栏（避免 cell_type/celltype 同时存在但不一致）
    if "cell_type" in sc_meta.columns and "celltype" in sc_meta.columns:
        a = sc_meta["cell_type"].astype(str).to_numpy()
        b = sc_meta["celltype"].astype(str).to_numpy()
        if a.shape == b.shape and (a != b).any():
            raise ValueError("sc_metadata.csv 同时存在 cell_type 与 celltype，但两列不一致（禁止 silent bug）")
    if "cell_type" in sc_meta.columns:
        type_col = "cell_type"
    elif "celltype" in sc_meta.columns:
        type_col = "celltype"
    elif "type" in sc_meta.columns:
        sc_meta = sc_meta.copy()
        sc_meta["celltype"] = sc_meta["type"].astype(str)
        type_col = "celltype"
    else:
        raise KeyError("sc_metadata.csv 需要包含 cell_type/celltype/type 之一作为类型列")
    sc_meta["type_col"] = sc_meta[type_col]
    # 坐标列容错
    coords_xy = _get_coord_xy(st_coords)
    st_coords = pd.concat([st_coords.drop(columns=[c for c in ["row", "col", "x", "y"] if c in st_coords.columns]), coords_xy], axis=1)
    return sc_expr, st_expr, sc_meta, st_coords, type_col



def _load_stage2(stage2_dir: Path):
    plugin_path = stage2_dir / "plugin_genes.txt"
    plugin_genes = [g.strip() for g in plugin_path.read_text(encoding="utf-8").splitlines() if g.strip()]
    gw_path = stage2_dir / "gene_weights.csv"
    gene_weights = pd.read_csv(gw_path)
    if "gene" not in gene_weights.columns:
        raise ValueError("gene_weights.csv ??????????????? gene")
    if "cost_weight" in gene_weights.columns:
        weight_field = "cost_weight"
        weights = pd.to_numeric(gene_weights["cost_weight"], errors="coerce").fillna(1.0).astype(float)
    elif "final_weight" in gene_weights.columns:
        weight_field = "final_weight_compat"
        weights = pd.to_numeric(gene_weights["final_weight"], errors="coerce").fillna(1.0).astype(float)
        weights = np.maximum(weights, 1.0)
    else:
        raise ValueError("gene_weights.csv ???????????? cost_weight ???final_weight")
    weight_map = dict(zip(gene_weights["gene"], weights))
    return plugin_genes, weight_map, weight_field, gene_weights


def _load_stage3(stage3_dir: Path):
    relabel_path = stage3_dir / "cell_type_relabel.csv"
    type_prior_path = stage3_dir / "type_prior_matrix.csv"
    relabel = pd.read_csv(relabel_path)
    if "cell_id" not in relabel.columns or "plugin_type" not in relabel.columns:
        raise ValueError("cell_type_relabel.csv 需要包含列 cell_id 与 plugin_type")
    type_prior = pd.read_csv(type_prior_path, index_col=0)
    return relabel, type_prior


def _compute_cells_per_spot(st_coords: pd.DataFrame, cfg: Dict[str, Any]) -> Tuple[pd.Series, str, Optional[float]]:
    """
    Compute raw (float) cells-per-spot series using the configured source.
    Risk#7: the *final* rounding/clipping to int must be done by _resolve_cells_per_spot.
    Returns (capacity_raw, source_resolved, norm_val_if_any).
    """
    src = cfg.get("cells_per_spot_source", "auto") or "auto"
    src = str(src)

    def _has_col(name: str) -> bool:
        return name in st_coords.columns and st_coords[name].notna().any()

    if src.lower() == "auto":
        if _has_col("spot_cell_counts"):
            src = "spot_cell_counts"
        elif _has_col("UMI_total") or _has_col("nCount_RNA") or _has_col("nCount_Spatial"):
            src = "UMI_total"
        else:
            src = "uniform"

    if src == "spot_cell_counts":
        if "spot_cell_counts" not in st_coords.columns:
            raise KeyError("cells_per_spot_source=spot_cell_counts 但 st_coords 缺少 spot_cell_counts 列")
        cps = st_coords["spot_cell_counts"].fillna(0).astype(float)
        cps.index = st_coords.index
        return cps, "spot_cell_counts", None

    if src == "UMI_total":
        umi_col = None
        for c in ("UMI_total", "nCount_RNA", "nCount_Spatial"):
            if c in st_coords.columns and st_coords[c].notna().any():
                umi_col = c
                break
        if umi_col is None:
            raise KeyError("cells_per_spot_source=UMI_total 但 st_coords 缺少 UMI_total/nCount_RNA/nCount_Spatial 列")
        umi = pd.to_numeric(st_coords[umi_col], errors="coerce").fillna(0).astype(float)
        umi.index = st_coords.index
        norm = cfg.get("umi_to_cell_norm", 1000)
        if isinstance(norm, str) and norm.strip().lower() == "median":
            norm_val = float(np.median(umi[umi > 0])) if (umi > 0).any() else 1000.0
        else:
            norm_val = float(norm)
        cps = umi / max(norm_val, 1e-8)
        cps = cps.astype(float)
        cps.index = st_coords.index
        src_resolved = "UMI_total" if umi_col == "UMI_total" else f"UMI_total[{umi_col}]"
        return cps, src_resolved, norm_val

    if src == "uniform":
        default_c = float(cfg.get("default_cells_per_spot", 1.0))
        cps = pd.Series(float(default_c), index=st_coords.index, dtype=float)
        return cps, "uniform", None

    raise ValueError(f"invalid cells_per_spot_source={src!r}; expected auto/spot_cell_counts/UMI_total/uniform")


def _resolve_cells_per_spot(capacity_raw: pd.Series, cfg: Dict[str, Any], *, context: str) -> Tuple[pd.Series, Dict[str, Any]]:
    """
    Risk#7: produce final int capacity series for cells_per_spot.csv with audit.
    - rounding: round|floor|ceil
    - clip: [clip_min, clip_max]
    """
    if not isinstance(capacity_raw, pd.Series):
        raise TypeError(f"[{context}] capacity_raw must be pd.Series, got {type(capacity_raw)!r}")

    rounding = str(cfg.get("cells_per_spot_rounding", "round") or "round").strip().lower()
    if rounding not in ("round", "floor", "ceil"):
        raise ValueError(f"[{context}] invalid cells_per_spot_rounding={rounding!r}")

    clip_min_v = cfg.get("cells_per_spot_clip_min", 1)
    if clip_min_v in (None, ""):
        clip_min_v = 1
    clip_min = int(clip_min_v)
    clip_max_v = cfg.get("cells_per_spot_clip_max", None)
    clip_max = None if clip_max_v in (None, "") else int(clip_max_v)
    if clip_max is not None and clip_max < clip_min:
        raise ValueError(f"[{context}] cells_per_spot_clip_max {clip_max} < clip_min {clip_min}")

    x0 = np.asarray(capacity_raw, dtype=float)
    n_nan = int(np.isnan(x0).sum())
    x0 = np.nan_to_num(x0, nan=0.0, posinf=0.0, neginf=0.0)

    if rounding == "round":
        xr = np.rint(x0)
    elif rounding == "floor":
        xr = np.floor(x0)
    else:
        xr = np.ceil(x0)

    xr_int = xr.astype(int)
    n_nonpositive_before_clip = int((xr_int <= 0).sum())
    n_clipped_min = int((xr_int < clip_min).sum())
    n_clipped_max = int((xr_int > clip_max).sum()) if clip_max is not None else 0

    x = xr_int
    if clip_min is not None:
        x = np.maximum(x, clip_min)
    if clip_max is not None:
        x = np.minimum(x, clip_max)
    x = x.astype(int)

    cap_int = pd.Series(x, index=capacity_raw.index, dtype=int)
    audit = {
        "context": context,
        "rounding": rounding,
        "clip_min": int(clip_min),
        "clip_max": int(clip_max) if clip_max is not None else None,
        "n_spots": int(len(cap_int)),
        "n_nan_before": int(n_nan),
        "n_nonpositive_before_clip": int(n_nonpositive_before_clip),
        "n_clipped_min": int(n_clipped_min),
        "n_clipped_max": int(n_clipped_max),
        "min": int(cap_int.min()) if len(cap_int) else None,
        "mean": float(cap_int.mean()) if len(cap_int) else None,
        "max": int(cap_int.max()) if len(cap_int) else None,
        "sum": int(cap_int.sum()) if len(cap_int) else None,
    }
    return cap_int, audit


def _build_feature_mats(sc_expr: pd.DataFrame, st_expr: pd.DataFrame, genes: List[str], weights: Dict[str, float]):
    genes_use = [g for g in genes if g in sc_expr.columns and g in st_expr.columns]
    if len(genes_use) == 0:
        raise ValueError("可用基因为空，请检查插件基因与表达矩阵交集")
    w = np.array([weights.get(g, 1.0) for g in genes_use], dtype=float)
    sc_mat = sc_expr[genes_use].to_numpy(dtype=float)
    st_mat = st_expr[genes_use].to_numpy(dtype=float)
    return sc_mat, st_mat, w, genes_use


def _apply_cost_expr_norm(
    sc_expr: pd.DataFrame,
    st_expr: pd.DataFrame,
    genes_use: List[str],
    cfg: Dict[str, Any],
    eps: float,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, Any]]:
    mode = str(cfg.get("cost_expr_norm") or "none").lower()
    if mode in ("", "none"):
        return sc_expr, st_expr, {"mode": "none"}

    st_vals = st_expr[genes_use].to_numpy(dtype=float)
    st_mean = st_vals.mean(axis=0)
    audit: Dict[str, Any] = {"mode": mode}

    if mode in ("st_center", "st_zscore"):
        if mode == "st_center":
            sc_vals = sc_expr[genes_use].to_numpy(dtype=float) - st_mean
            st_vals = st_vals - st_mean
        else:
            st_std = st_vals.std(axis=0, ddof=0)
            st_std = np.where(st_std > eps, st_std, 1.0)
            sc_vals = (sc_expr[genes_use].to_numpy(dtype=float) - st_mean) / st_std
            st_vals = (st_vals - st_mean) / st_std
            audit["std_min"] = float(np.min(st_std))
            audit["std_max"] = float(np.max(st_std))

        min_val = float(min(np.min(sc_vals), np.min(st_vals)))
        if min_val < 0:
            shift = -min_val
            sc_vals = sc_vals + shift
            st_vals = st_vals + shift
            audit["shift_nonneg"] = float(shift)
        else:
            audit["shift_nonneg"] = 0.0

        audit["mean_min"] = float(np.min(st_mean))
        audit["mean_max"] = float(np.max(st_mean))
        sc_out = sc_expr.copy()
        st_out = st_expr.copy()
        sc_out.loc[:, genes_use] = sc_vals
        st_out.loc[:, genes_use] = st_vals
        return sc_out, st_out, audit

    if mode != "gene_scale":
        raise ValueError(f"invalid cost_expr_norm={mode!r}")

    beta_raw = cfg.get("cost_scale_beta", None)
    scale_beta = 1.0 if beta_raw is None else float(beta_raw)
    scale_beta = float(np.clip(scale_beta, 0.0, 1.0))
    if scale_beta <= 0.0:
        q_list = cfg.get("cost_scale_report_quantiles") or [0.01, 0.05, 0.5, 0.95, 0.99]
        q_list = [float(q) for q in q_list]
        q_stats = {f"q{int(round(q * 100)):02d}": 1.0 for q in q_list}
        clip_ratio_soft = float(cfg.get("cost_scale_clip_ratio_soft") or 0.85)
        clip_ratio_max = float(cfg.get("cost_scale_clip_ratio_max") or 0.95)
        audit = {
            "mode": "gene_scale",
            "scale_source": str(cfg.get("cost_scale_source") or "st").lower(),
            "scale_metric": str(cfg.get("cost_scale_metric") or "std").lower(),
            "scale_beta": scale_beta,
            "scale_downweight_only": bool(cfg.get("cost_scale_downweight_only", False)),
            "scale_eps": float(cfg.get("cost_scale_eps") or eps),
            "scale_floor": float(cfg.get("cost_scale_floor") or 0.0),
            "scale_clip": cfg.get("cost_scale_clip"),
            "clip_ratio": 0.0,
            "clip_ratio_soft": float(clip_ratio_soft),
            "clip_ratio_max": float(clip_ratio_max),
            "scale_stats": {
                "scale_min": None,
                "scale_max": None,
                "n_floor_applied": 0,
                "n_clip_low": 0,
                "n_clip_high": 0,
                **{k: None for k in q_stats.keys()},
            },
            "multiplier_stats_planned": {
                "min": 1.0,
                "max": 1.0,
                "n_clip_low": 0,
                "n_clip_high": 0,
                **q_stats,
            },
            "multiplier_stats_applied": {
                "min": 1.0,
                "max": 1.0,
                "n_clip_low": 0,
                "n_clip_high": 0,
                **q_stats,
            },
            "genes_use_count": int(len(genes_use)),
            "norm_applied": False,
            "norm_fallback": False,
            "norm_reason": "beta_is_zero_skip",
            "norm_fallback_reason": [],
            "scale_min_threshold": 1e-4,
        }
        return sc_expr, st_expr, audit

    source = str(cfg.get("cost_scale_source") or "st").lower()
    if source == "joint":
        base = pd.concat([sc_expr[genes_use], st_expr[genes_use]], axis=0)
    else:
        source = "st"
        base = st_expr[genes_use]

    scale_metric = str(cfg.get("cost_scale_metric") or "std").lower()
    base_vals = base.to_numpy(dtype=float)
    if scale_metric == "std":
        scale = base_vals.std(axis=0, ddof=0)
    elif scale_metric == "mean":
        scale = base_vals.mean(axis=0)
    elif scale_metric == "p95":
        scale = np.quantile(base_vals, 0.95, axis=0)
    else:
        raise ValueError(f"invalid cost_scale_metric={scale_metric!r}")
    scale = np.where(np.isfinite(scale), scale, 0.0)
    scale_eps = float(cfg.get("cost_scale_eps") or eps)
    scale_floor = float(cfg.get("cost_scale_floor") or 0.0)
    scale = np.maximum(scale, scale_eps)
    scale = np.maximum(scale, scale_floor)

    mult_raw = 1.0 / scale
    mult = 1.0 + scale_beta * (mult_raw - 1.0)
    downweight_only = bool(cfg.get("cost_scale_downweight_only", False))
    if downweight_only:
        mult = np.minimum(mult, 1.0)
    clip = cfg.get("cost_scale_clip")
    clip_low = clip_high = None
    n_clip_low = n_clip_high = 0
    if clip:
        clip_low, clip_high = float(clip[0]), float(clip[1])
        n_clip_low = int((mult < clip_low).sum())
        n_clip_high = int((mult > clip_high).sum())
        mult = np.clip(mult, clip_low, clip_high)

    sc_vals = sc_expr[genes_use].to_numpy(dtype=float) * mult
    st_vals = st_expr[genes_use].to_numpy(dtype=float) * mult

    q_list = cfg.get("cost_scale_report_quantiles") or [0.01, 0.05, 0.5, 0.95, 0.99]
    q_list = [float(q) for q in q_list]
    scale_stats = {}
    mult_stats_before = {}
    mult_stats_after = {}
    for q in q_list:
        label = f"q{int(round(q * 100)):02d}"
        scale_stats[label] = float(np.quantile(scale, q))
        mult_stats_before[label] = float(np.quantile(mult_raw, q))
        mult_stats_after[label] = float(np.quantile(mult, q))

    n_floor_applied = int((scale <= scale_floor + 1e-12).sum())
    n_genes = int(len(scale))
    floor_ratio = n_floor_applied / max(n_genes, 1)
    clip_ratio = (n_clip_low + n_clip_high) / max(n_genes, 1)
    scale_min = float(np.min(scale)) if n_genes else None
    scale_q01 = scale_stats.get("q01")
    min_threshold = 1e-4
    clip_ratio_soft = float(cfg.get("cost_scale_clip_ratio_soft") or 0.85)
    clip_ratio_max = float(cfg.get("cost_scale_clip_ratio_max") or 0.95)
    if clip_ratio > clip_ratio_soft:
        audit["clip_ratio_warning"] = float(clip_ratio)

    fallback_reasons = []
    if floor_ratio > 0.3:
        fallback_reasons.append(f"floor_ratio_too_high:{floor_ratio:.3f}")
    if scale_min is not None and scale_min < min_threshold:
        fallback_reasons.append(f"scale_too_small:{scale_min:.6g}")
    if scale_q01 is not None and scale_q01 < min_threshold:
        fallback_reasons.append(f"scale_q01_too_small:{scale_q01:.6g}")
    if clip and clip_ratio > clip_ratio_max:
        fallback_reasons.append(f"clip_ratio_too_high:{clip_ratio:.3f}")

    audit.update(
        {
            "mode": "gene_scale",
            "scale_source": source,
            "scale_metric": scale_metric,
            "scale_beta": scale_beta,
            "scale_downweight_only": downweight_only,
            "scale_eps": scale_eps,
            "scale_floor": scale_floor,
            "scale_clip": clip if clip else None,
            "clip_ratio": float(clip_ratio),
            "clip_ratio_soft": float(clip_ratio_soft),
            "clip_ratio_max": float(clip_ratio_max),
            "scale_stats": {
                "scale_min": scale_min,
                "scale_max": float(np.max(scale)) if n_genes else None,
                "n_floor_applied": n_floor_applied,
                "n_clip_low": n_clip_low,
                "n_clip_high": n_clip_high,
                **scale_stats,
            },
            "multiplier_stats_before_beta": {
                "min": float(np.min(mult_raw)) if n_genes else None,
                "max": float(np.max(mult_raw)) if n_genes else None,
                **mult_stats_before,
            },
            "genes_use_count": n_genes,
            "norm_fallback": len(fallback_reasons) > 0,
            "norm_fallback_reason": fallback_reasons,
            "scale_min_threshold": min_threshold,
        }
    )

    planned_stats = {
        "min": float(np.min(mult)) if n_genes else None,
        "max": float(np.max(mult)) if n_genes else None,
        "n_clip_low": n_clip_low,
        "n_clip_high": n_clip_high,
        **mult_stats_after,
    }
    audit["multiplier_stats_planned"] = planned_stats
    if fallback_reasons:
        applied_stats = {
            "min": 1.0,
            "max": 1.0,
            "n_clip_low": 0,
            "n_clip_high": 0,
            **{k: 1.0 for k in mult_stats_after.keys()},
        }
        audit["multiplier_stats_applied"] = applied_stats
        audit["norm_applied"] = False
        return sc_expr, st_expr, audit

    audit["multiplier_stats_applied"] = planned_stats
    audit["norm_applied"] = True

    sc_out = sc_expr.copy()
    st_out = st_expr.copy()
    sc_out.loc[:, genes_use] = sc_vals
    st_out.loc[:, genes_use] = st_vals
    return sc_out, st_out, audit


def _compute_spot_weight_matrix(
    st_expr: pd.DataFrame,
    st_coords: pd.DataFrame,
    genes_use: List[str],
    cfg: Dict[str, Any],
    *,
    eps: float,
    candidate_mask: Optional[np.ndarray] = None,
) -> Tuple[Optional[np.ndarray], Dict[str, Any]]:
    mode = str(cfg.get("spot_weight_mode") or "none").lower()
    if mode in ("", "none"):
        return None, {"status": "disabled", "mode": "none"}
    if mode != "local_zscore":
        raise ValueError(f"invalid spot_weight_mode={mode!r}")

    topk = int(cfg.get("spot_weight_topk") or 0)
    if topk <= 0:
        return None, {"status": "disabled", "mode": mode, "reason": "topk<=0"}

    k = int(cfg.get("spot_weight_k") or 8)
    kappa = float(cfg.get("spot_weight_kappa") or 0.0)
    weight_max = float(cfg.get("spot_weight_weight_max") or 1.0)
    scale = str(cfg.get("spot_weight_scale") or "linear").lower()
    only_up = bool(cfg.get("spot_weight_only_up", True))

    coords = _get_coord_xy(st_coords).to_numpy(dtype=float)
    n_spots = coords.shape[0]
    if n_spots == 0:
        return None, {"status": "disabled", "mode": mode, "reason": "no_spots"}
    k_use = min(max(k, 1), n_spots)
    try:
        tree = cKDTree(coords)
        _, idxs = tree.query(coords, k=k_use)
    except Exception:
        dists = cdist(coords, coords)
        idxs = np.argsort(dists, axis=1)[:, :k_use]

    if idxs.ndim == 1:
        idxs = idxs[:, None]

    st_vals = st_expr[genes_use].to_numpy(dtype=float)
    n_genes = st_vals.shape[1]
    topk = min(topk, n_genes)

    if candidate_mask is not None:
        candidate_mask = np.asarray(candidate_mask).astype(bool)
        if candidate_mask.shape[0] != n_genes:
            raise ValueError("spot_weight candidate_mask shape mismatch")
        cand_idx = np.where(candidate_mask)[0]
    else:
        cand_idx = None

    weights = np.ones((n_spots, n_genes), dtype=float)
    weighted_per_spot = []
    for i in range(n_spots):
        neigh_idx = idxs[i]
        local = st_vals[neigh_idx]
        mean = local.mean(axis=0)
        std = local.std(axis=0)
        score = np.abs(st_vals[i] - mean) / (std + eps)
        if cand_idx is not None:
            if len(cand_idx) == 0:
                continue
            score_cand = score[cand_idx]
            topk_i = min(topk, len(cand_idx))
            pick = cand_idx[np.argpartition(-score_cand, topk_i - 1)[:topk_i]]
        else:
            topk_i = topk
            pick = np.argpartition(-score, topk_i - 1)[:topk_i]
        if topk_i <= 0:
            continue
        smax = float(np.max(score[pick]))
        if smax <= 0:
            continue
        score_norm = score[pick] / (smax + eps)
        if scale == "exp":
            w = np.exp(kappa * score_norm)
        else:
            w = 1.0 + kappa * score_norm
        if only_up:
            w = np.maximum(w, 1.0)
        w = np.clip(w, 1.0, weight_max)
        weights[i, pick] = w
        weighted_per_spot.append(int(len(pick)))

    stats = {
        "status": "ok",
        "mode": mode,
        "scale": scale,
        "topk": int(topk),
        "k": int(k_use),
        "kappa": float(kappa),
        "weight_max": float(weight_max),
        "n_spots": int(n_spots),
        "n_genes": int(n_genes),
        "weighted_per_spot_mean": float(np.mean(weighted_per_spot)) if weighted_per_spot else 0.0,
        "weights_min": float(np.min(weights)),
        "weights_p50": float(np.percentile(weights, 50)),
        "weights_p95": float(np.percentile(weights, 95)),
        "weights_max": float(np.max(weights)),
        "weighted_ratio": float((weights > 1.0).mean()),
        "candidate_n": int(len(cand_idx)) if cand_idx is not None else None,
    }
    return weights, stats


def _spot_fraction_from_mat(
    mat: csr_matrix,
    type_labels: List[str],
    st_index: List[str],
    type_order: List[str],
    *,
    context: str = "",
    min_type_intersection: float = 0.2,
    max_unique_ratio: float = 0.8,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """
    从 spot×cell 矩阵聚合得到 spot×type 计数矩阵（未归一化），并返回校验/审计信息。
    目的：防止把 cell_id 等错误输入当成 type_labels 传入，导致 prior_effect_* 指标“沉默失真”。
    """
    if not sparse.isspmatrix_csr(mat):
        mat = mat.tocsr()

    n_spots, n_cells = mat.shape
    if len(type_labels) != n_cells:
        raise ValueError(
            f"[{context}] type_labels 长度 {len(type_labels)} != mat.shape[1] {n_cells}，疑似 cell 维度未对齐"
        )

    type_order = list(type_order)
    if len(type_order) == 0:
        raise ValueError(f"[{context}] type_order 为空，无法计算 spot×type")

    uniq = pd.Series(type_labels).astype(str).unique().tolist()
    inter = len(set(uniq) & set(type_order))
    inter_ratio = inter / max(len(set(uniq)), 1)
    unique_ratio = len(set(type_labels)) / max(n_cells, 1)

    if inter == 0:
        raise ValueError(
            f"[{context}] type_labels 与 type_order 无交集（inter=0），高度疑似把 cell_id 当 cell_types 传入"
        )
    if inter_ratio < float(min_type_intersection):
        raise ValueError(
            f"[{context}] type_labels 与 type_order 交集比例过低 inter_ratio={inter_ratio:.3f}，疑似类型体系不一致/传参错误"
        )
    if unique_ratio > float(max_unique_ratio):
        raise ValueError(
            f"[{context}] type_labels 唯一值比例过高 unique_ratio={unique_ratio:.3f}，高度疑似传入近似唯一的 cell_id 列表"
        )

    df = pd.DataFrame(0.0, index=st_index, columns=type_order)
    type_arr = np.array([str(x) for x in type_labels], dtype=str)

    missing_mask = ~np.isin(type_arr, np.array(type_order, dtype=str))
    n_missing = int(missing_mask.sum())
    missing_examples = sorted(set(type_arr[missing_mask]))[:5] if n_missing > 0 else []

    for t in type_order:
        mask = type_arr == t
        if mask.any():
            df[t] = mat[:, mask].sum(axis=1).A1

    audit = {
        "context": context,
        "n_spots": int(n_spots),
        "n_cells": int(n_cells),
        "n_types_unique": int(len(set(type_labels))),
        "unique_ratio": float(unique_ratio),
        "type_intersection": int(inter),
        "type_intersection_ratio": float(inter_ratio),
        "n_missing_type_labels": int(n_missing),
        "missing_type_label_examples": missing_examples,
        "min_type_intersection": float(min_type_intersection),
        "max_unique_ratio": float(max_unique_ratio),
    }
    return df, audit


def _align_spot_inputs(
    *,
    spot_ids: List[str],
    st_coords: pd.DataFrame,
    capacity: pd.Series,
    type_prior_raw: Optional[pd.DataFrame] = None,
    context: str = "",
) -> Tuple[pd.DataFrame, pd.Series, Optional[pd.DataFrame], Dict[str, Any]]:
    """
    Stage4 Risk#5 防护：以 spot_ids 为唯一锚点，对齐所有 spot 相关输入并生成可审计信息。

    约束：
    - 不允许静默丢 spot：任何 missing/extra 都直接 raise
    - 对齐仅允许在此函数完成（统一 reindex 到 spot_ids 顺序）
    """

    def _examples(xs) -> List[str]:
        return sorted([str(x) for x in xs])[:5]

    if len(spot_ids) != len(set(spot_ids)):
        raise ValueError(f"[{context}] spot_ids 存在重复（len={len(spot_ids)}, unique={len(set(spot_ids))}）")

    if not st_coords.index.is_unique:
        dup = st_coords.index[st_coords.index.duplicated()].unique()[:5]
        raise ValueError(f"[{context}] st_coords index 存在重复: {list(dup)}")
    if not capacity.index.is_unique:
        dup = capacity.index[capacity.index.duplicated()].unique()[:5]
        raise ValueError(f"[{context}] capacity index 存在重复: {list(dup)}")
    if type_prior_raw is not None and not type_prior_raw.index.is_unique:
        dup = type_prior_raw.index[type_prior_raw.index.duplicated()].unique()[:5]
        raise ValueError(f"[{context}] type_prior_raw index 存在重复: {list(dup)}")

    spot_set = set(spot_ids)
    coords_set = set(st_coords.index)
    capacity_set = set(capacity.index)

    missing_in_coords = spot_set - coords_set
    extra_in_coords = coords_set - spot_set
    missing_in_capacity = spot_set - capacity_set
    extra_in_capacity = capacity_set - spot_set

    if missing_in_coords or extra_in_coords:
        raise ValueError(
            f"[{context}] st_coords spot_id mismatch: missing={len(missing_in_coords)}, extra={len(extra_in_coords)}, "
            f"missing_examples={_examples(missing_in_coords)}, extra_examples={_examples(extra_in_coords)}"
        )
    if missing_in_capacity or extra_in_capacity:
        raise ValueError(
            f"[{context}] capacity spot_id mismatch: missing={len(missing_in_capacity)}, extra={len(extra_in_capacity)}, "
            f"missing_examples={_examples(missing_in_capacity)}, extra_examples={_examples(extra_in_capacity)}"
        )

    type_prior_set = None
    missing_in_type_prior: set = set()
    extra_in_type_prior: set = set()
    if type_prior_raw is not None:
        type_prior_set = set(type_prior_raw.index)
        missing_in_type_prior = spot_set - type_prior_set
        extra_in_type_prior = type_prior_set - spot_set
        if missing_in_type_prior or extra_in_type_prior:
            raise ValueError(
                f"[{context}] type_prior_raw spot_id mismatch: missing={len(missing_in_type_prior)}, extra={len(extra_in_type_prior)}, "
                f"missing_examples={_examples(missing_in_type_prior)}, extra_examples={_examples(extra_in_type_prior)}"
            )

    anchor = pd.Index(spot_ids)
    order_coords = bool(st_coords.index.equals(anchor))
    order_capacity = bool(capacity.index.equals(anchor))
    order_type_prior = bool(type_prior_raw.index.equals(anchor)) if type_prior_raw is not None else None
    order_all = order_coords and order_capacity and (order_type_prior if order_type_prior is not None else True)

    coords_cols_used = "unknown"
    if {"row", "col"}.issubset(st_coords.columns):
        coords_cols_used = "row/col"
    elif {"x", "y"}.issubset(st_coords.columns):
        coords_cols_used = "x/y"

    spot_ids_sha1 = hashlib.sha1(",".join([str(x) for x in spot_ids]).encode("utf-8")).hexdigest().lower()

    audit = {
        "context": context,
        "spot_ids_sha1": spot_ids_sha1,
        "order_equal_before": bool(order_all),
        "order_equal_before_coords": bool(order_coords),
        "order_equal_before_capacity": bool(order_capacity),
        "order_equal_before_type_prior": order_type_prior,
        "n_spots_expr": int(len(spot_ids)),
        "n_spots_coords": int(len(st_coords.index)),
        "n_spots_capacity": int(len(capacity.index)),
        "n_spots_type_prior": int(len(type_prior_raw.index)) if type_prior_raw is not None else None,
        "n_missing_in_coords": int(len(missing_in_coords)),
        "example_missing_in_coords": _examples(missing_in_coords),
        "n_missing_in_capacity": int(len(missing_in_capacity)),
        "example_missing_in_capacity": _examples(missing_in_capacity),
        "n_missing_in_type_prior": int(len(missing_in_type_prior)) if type_prior_raw is not None else None,
        "example_missing_in_type_prior": _examples(missing_in_type_prior) if type_prior_raw is not None else None,
        "coords_cols_used": coords_cols_used,
        "head_spot_ids": [str(x) for x in spot_ids[:5]],
        "tail_spot_ids": [str(x) for x in spot_ids[-5:]],
    }

    st_coords_aligned = st_coords.reindex(spot_ids)
    capacity_aligned = capacity.reindex(spot_ids)
    if capacity_aligned.isna().any():
        raise ValueError(f"[{context}] capacity reindex 后存在 NA（疑似 spot_ids 锚点不一致）")
    type_prior_aligned = type_prior_raw.reindex(spot_ids) if type_prior_raw is not None else None
    return st_coords_aligned, capacity_aligned, type_prior_aligned, audit


def _write_cytospace_inputs(
    sc_expr: pd.DataFrame,
    sc_meta: pd.DataFrame,
    st_expr: pd.DataFrame,
    st_coords: pd.DataFrame,
    cell_types: List[str],
    genes_use: List[str],
    capacity: pd.Series,
    spot_ids: List[str],
    tmp_dir: Path,
    prefix: str,
):
    sc_sub = sc_expr[genes_use]
    st_sub = st_expr[genes_use]

    anchor = pd.Index(spot_ids)
    if not st_expr.index.equals(anchor):
        raise ValueError(f"[{prefix}] st_expr.index 与 spot_ids 锚点不一致（可能存在隐形错位）")
    if not st_coords.index.equals(anchor):
        raise ValueError(f"[{prefix}] st_coords.index 与 spot_ids 锚点不一致（可能存在隐形错位）")
    if not capacity.index.equals(anchor):
        raise ValueError(f"[{prefix}] capacity.index 与 spot_ids 锚点不一致（可能存在隐形错位）")

    coords_xy = _get_coord_xy(st_coords)

    sc_path = tmp_dir / f"{prefix}scRNA.csv"
    st_path = tmp_dir / f"{prefix}st.csv"
    ct_path = tmp_dir / f"{prefix}cell_type.csv"
    coord_path = tmp_dir / f"{prefix}coords.csv"
    cps_path = tmp_dir / f"{prefix}cells_per_spot.csv"

    sc_sub.T.to_csv(sc_path)
    st_sub.T.to_csv(st_path)
    # CytoSPACE 的 read_file 约定：CSV 第一列为 rownames，且必须有 header
    ct_df = pd.DataFrame({"cell_type": cell_types}, index=sc_expr.index)
    ct_df.index.name = "cell_id"
    ct_df.to_csv(ct_path)

    coord_df = pd.DataFrame(
        {"row": coords_xy.iloc[:, 0].values, "col": coords_xy.iloc[:, 1].values},
        index=anchor,
    )
    coord_df.index.name = "spot_id"
    coord_df.to_csv(coord_path)

    cps_df = pd.DataFrame({"cells": np.asarray(capacity, dtype=int)}, index=anchor)
    cps_df.index.name = "spot_id"
    cps_df.to_csv(cps_path)
    return sc_path, ct_path, st_path, coord_path, cps_path


def _parse_assigned_locations(
    assigned: pd.DataFrame,
    st_ids: List[str],
    require_original_cid: bool = False,
) -> Tuple[List[str], List[str], Optional[List[str]], List[Tuple[int, int, float]]]:
    """
    解析 CytoSPACE 的 assigned_locations.csv：
    - 使用 UniqueCID 作为 cell_id（每个“cell instance”唯一）
    - 使用 CellType 作为 type（baseline 为 orig_type；plus 为 plugin_type）
    - SpotID 必须能在 Stage1 的 spot_id 集合中找到
    """
    required_cols = {"UniqueCID", "SpotID", "CellType"}
    missing = sorted(required_cols - set(assigned.columns))
    if missing:
        raise ValueError(f"assigned_locations.csv 缺少必需列: {missing}")
    has_original = "OriginalCID" in assigned.columns
    if require_original_cid and not has_original:
        raise ValueError("assigned_locations.csv 缺少列 OriginalCID（用于从 Stage1 sc_expr 构建 cell pool 表达）")

    spot_pos = {sid: i for i, sid in enumerate(st_ids)}
    cell_ids = assigned["UniqueCID"].astype(str).tolist()
    cell_types = assigned["CellType"].astype(str).tolist()
    orig_cell_ids = assigned["OriginalCID"].astype(str).tolist() if has_original else None

    assignments: List[Tuple[int, int, float]] = []
    bad_spots: List[str] = []
    for cell_idx, sid in enumerate(assigned["SpotID"].astype(str).tolist()):
        if sid not in spot_pos:
            if len(bad_spots) < 5:
                bad_spots.append(sid)
            continue
        assignments.append((cell_idx, spot_pos[sid], 1.0))
    if len(assignments) != len(cell_ids):
        raise ValueError(f"assigned_locations 中有 spot_id 未匹配 Stage1：{len(cell_ids)-len(assignments)} 行，示例: {bad_spots}")
    return cell_ids, cell_types, orig_cell_ids, assignments


def _assignments_to_matrix(assignments: List[Tuple[int, int, float]], n_spots: int, n_cells: int) -> sparse.csr_matrix:
    if len(assignments) == 0:
        return sparse.csr_matrix((n_spots, n_cells), dtype=float)
    data = np.array([score for _, _, score in assignments], dtype=float)
    row_ind = [spot_idx for _, spot_idx, _ in assignments]
    col_ind = [cell_idx for cell_idx, _, _ in assignments]
    return sparse.csr_matrix((data, (row_ind, col_ind)), shape=(n_spots, n_cells))


def _build_sparse_knn_safe(
    coords: pd.DataFrame,
    k: int,
    metric: str = "euclidean",
    block_size: int = 1024,
    max_dense_n: int = 5000,
) -> Tuple[csr_matrix, str]:
    """
    构建稀疏KNN 邻接矩阵（避免大规模 SxS dense）：
    - S <= max_dense_n: dense 计算后转 CSR，mode="dense"
    - S  大: 分块 brute-force，内存占用约 block_size x S，mode="block"
    返回 (W, mode)
    """
    xy = coords[["row", "col"]].to_numpy(dtype=float)
    n = xy.shape[0]
    if n <= max_dense_n:
        dist = cdist(xy, xy, metric=metric)
        np.fill_diagonal(dist, np.inf)
        idx_knn = np.argpartition(dist, kth=k, axis=1)[:, :k]
        rows = np.repeat(np.arange(n), k)
        cols = idx_knn.reshape(-1)
        data = np.ones_like(cols, dtype=float)
        W = sparse.csr_matrix((data, (rows, cols)), shape=(n, n))
        row_sum = np.asarray(W.sum(axis=1)).ravel()
        row_sum[row_sum == 0] = 1.0
        W = W.multiply(1.0 / row_sum[:, None])
        return W, "dense"

    rows_list, cols_list, data_list = [], [], []
    all_xy = xy.astype(np.float32, copy=False)
    for start in range(0, n, block_size):
        end = min(start + block_size, n)
        block_xy = all_xy[start:end]
        dist_block = cdist(block_xy, all_xy, metric=metric)
        for i in range(end - start):
            dist_block[i, start + i] = np.inf
        idx_knn = np.argpartition(dist_block, kth=k, axis=1)[:, :k]
        block_rows = np.repeat(np.arange(start, end), k)
        block_cols = idx_knn.reshape(-1)
        block_data = np.ones_like(block_cols, dtype=float)
        rows_list.append(block_rows)
        cols_list.append(block_cols)
        data_list.append(block_data)
    rows = np.concatenate(rows_list)
    cols = np.concatenate(cols_list)
    data = np.concatenate(data_list)
    W = sparse.csr_matrix((data, (rows, cols)), shape=(n, n))
    row_sum = np.asarray(W.sum(axis=1)).ravel()
    row_sum[row_sum == 0] = 1.0
    W = W.multiply(1.0 / row_sum[:, None])
    return W, "block"


def _row_norm(x: np.ndarray, eps: float) -> np.ndarray:
    norm = np.linalg.norm(x, axis=1, keepdims=True)
    return x / np.maximum(norm, eps)


def _refine_spot_cell_matrix_svg(
    mat: csr_matrix,
    st_coords: pd.DataFrame,
    st_svg: np.ndarray,
    sc_svg: np.ndarray,
    lambda_base: float,
    k: int,
    eps: float = 1e-8,
    batch_size: Optional[int] = None,
    knn_metric: str = "euclidean",
    knn_block_size: int = 1024,
    knn_max_dense_n: int = 5000,
) -> Tuple[csr_matrix, str]:
    """
    基于 SVG 重构误差的加权平滑：
    - 使用当前分配 mat 预测 ST 表达（plugin_genes）
    - 计算观测 vs 预测的余弦误差 per-spot
    - 将误差归一化到 [0,1] 并乘 lambda_base 得到每个 spot 的局部平滑系数
    - 使用稀疏 KNN 邻接矩阵在 spot 维度做加权平滑，保持稀疏
    """
    if lambda_base <= 0:
        return mat if sparse.isspmatrix_csr(mat) else mat.tocsr(), "skipped_lambda0"
    if not sparse.isspmatrix_csr(mat):
        mat = mat.tocsr()

    if st_svg.size == 0 or sc_svg.size == 0:
        return mat, "skipped_empty_svg"
    if st_svg.shape[0] != mat.shape[0]:
        raise ValueError(f"st_svg 行数 {st_svg.shape[0]} 与 mat spot 数 {mat.shape[0]} 不一致")
    if sc_svg.shape[0] != mat.shape[1]:
        raise ValueError(f"sc_svg 行数 {sc_svg.shape[0]} 与 mat cell 数 {mat.shape[1]} 不一致")

    # 预测表达：SxC · CxG -> SxG
    pred_svg = mat.dot(sc_svg)
    pred_n = _row_norm(pred_svg, eps)
    obs_n = _row_norm(st_svg, eps)
    cos = np.sum(pred_n * obs_n, axis=1)
    cos = np.clip(cos, -1.0, 1.0)
    err = 1.0 - cos
    err = np.nan_to_num(err, nan=0.0, posinf=1.0, neginf=0.0)
    err = err - err.min()
    if err.max() > 0:
        err = err / err.max()
    lambda_vec = lambda_base * err  # (S,)

    W, knn_mode = _build_sparse_knn_safe(
        st_coords,
        k,
        metric=knn_metric,
        block_size=knn_block_size,
        max_dense_n=knn_max_dense_n,
    )
    n_spots, n_cells = mat.shape
    if batch_size is None:
        smoothed = W.dot(mat)
    else:
        blocks = []
        for start in range(0, n_cells, batch_size):
            end = min(start + batch_size, n_cells)
            smoothed_block = W.dot(mat[:, start:end])
            blocks.append(smoothed_block)
        smoothed = sparse.hstack(blocks).tocsr()

    lambda_vec = lambda_vec.astype(float)
    part1 = mat.multiply((1.0 - lambda_vec)[:, None])
    part2 = smoothed.multiply(lambda_vec[:, None])
    mat_ref = (part1 + part2).tocsr()
    return mat_ref, knn_mode


def _write_meta(out_dir: Path, meta: Dict[str, Any]):
    with (out_dir / "meta.json").open("w", encoding="utf-8") as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)


def _type_prior_refine(mat: csr_matrix, cell_types: List[str], type_cols: List[str], type_prior_spot: pd.DataFrame, lambda_prior: float, eps: float = 1e-6) -> csr_matrix:
    """
    使用 per-spot type_prior 进行软约束微调：
    - 计算当前 spot×type 分布 Curr
    - 使用 type_prior_spot（行归一化）作为 Prior
    - 权重 w = (Prior+eps) / (Curr+eps) 的 lambda_prior 次方
    - 按 type 对应的列乘以 w，再按行保持原总量
    """
    if lambda_prior <= 0:
        return mat
    mat_csr = mat.tocsr()
    n_spots, n_cells = mat_csr.shape
    type_cols = list(type_cols)
    # map cell -> type index
    type_map = [type_cols.index(t) if t in type_cols else None for t in cell_types]
    type_map = np.array(type_map)
    # 当前分布
    row_sum = np.asarray(mat_csr.sum(axis=1)).ravel()
    curr = np.zeros((n_spots, len(type_cols)), dtype=float)
    for i, t in enumerate(type_cols):
        mask = type_map == i
        if mask.any():
            curr[:, i] = mat_csr[:, mask].sum(axis=1).A1
    curr_frac = curr / np.maximum(row_sum[:, None], eps)
    # 先验分布（行归一化）
    prior = type_prior_spot[type_cols]
    prior = prior.div(prior.sum(axis=1).replace(0, 1.0), axis=0)
    prior_np = prior.to_numpy(dtype=float)
    w = ((prior_np + eps) / (curr_frac + eps)) ** float(lambda_prior)
    # 应用权重：对每个非零条目，乘以对应 type 的权重
    coo = mat_csr.tocoo()
    w_lookup = w[coo.row, type_map[coo.col]]
    data_new = coo.data * w_lookup
    mat_new = coo_matrix((data_new, (coo.row, coo.col)), shape=mat_csr.shape).tocsr()
    # 行重新缩放保持总量
    new_row_sum = np.asarray(mat_new.sum(axis=1)).ravel()
    scale = row_sum / np.maximum(new_row_sum, eps)
    mat_new = mat_new.multiply(scale[:, None])
    return mat_new


def _compute_ablation_change_rate_checked(
    base_assign: List[Tuple[int, int, float]],
    ablate_assign: List[Tuple[int, int, float]],
    *,
    context: str = "prior_ablation",
) -> Tuple[float, Dict[str, Any]]:
    """
    计算 prior ablation 的 change_rate：同一 cell 在 base vs ablate 是否换 spot。
    - 强校验：禁止长度不一致/对齐失败时写“假值”
    - 以 cell_idx 为键做对齐，避免仅靠列表顺序导致误判
    返回 (change_rate, audit)
    """

    def build_map(assign: List[Tuple[int, int, float]]) -> Dict[int, int]:
        m: Dict[int, int] = {}
        dup = 0
        for cid, sid, _ in assign:
            if cid in m:
                dup += 1
            m[cid] = sid
        if dup > 0:
            raise ValueError(f"[{context}] assignments 内存在重复 cell_idx，dup={dup}")
        return m

    base_map = build_map(base_assign)
    abl_map = build_map(ablate_assign)

    base_keys = set(base_map.keys())
    abl_keys = set(abl_map.keys())
    common = base_keys & abl_keys
    only_base = base_keys - abl_keys
    only_abl = abl_keys - base_keys

    audit = {
        "context": context,
        "n_base": len(base_assign),
        "n_ablate": len(ablate_assign),
        "n_base_unique_cells": len(base_map),
        "n_ablate_unique_cells": len(abl_map),
        "n_common_cells": len(common),
        "n_only_in_base": len(only_base),
        "n_only_in_ablate": len(only_abl),
        "example_only_in_base": sorted(list(only_base))[:5],
        "example_only_in_ablate": sorted(list(only_abl))[:5],
    }

    if len(common) == 0:
        raise ValueError(f"[{context}] base 与 ablate 无共同 cell，无法计算 change_rate")
    if len(only_base) > 0 or len(only_abl) > 0:
        raise ValueError(
            f"[{context}] base/ablate cell 集合不一致：only_base={len(only_base)}, only_ablate={len(only_abl)}"
        )

    changed = 0
    for cid in common:
        if base_map[cid] != abl_map[cid]:
            changed += 1
    audit["n_changed"] = int(changed)
    change_rate = changed / max(len(common), 1)
    return float(change_rate), audit


def _harden_assignment_quota_matching(
    mat: csr_matrix,
    capacity: np.ndarray,
    cell_types: List[str],
    type_prior: pd.DataFrame,
    spot_ids: List[str],
    lambda_prior: float,
    eps: float = 1e-8,
    topk: int = 5,
    prior_candidate_topk: int = 0,
    prior_candidate_weight: float = 1.0,
    fallback_assignments: Optional[List[Tuple[int, int, float]]] = None,
    ablate_lambda: bool = False,
    min_prior_row_nonzero_ratio: float = 0.0,
) -> Tuple[List[Tuple[int, int, float]], csr_matrix, Dict[str, Any]]:
    """
    将 soft 矩阵投影回满足容量的 hard assignment：
    - 对每个 cell 取 topK 候选 + fallback 的原始 CytoSPACE 结果
    - 效用 U = log(score) + lambda_prior*log(type_prior)
    - 用带容量的延迟接受（deferred acceptance）匹配
    """
    if not sparse.isspmatrix_csr(mat):
        mat = mat.tocsr()
    cap = np.asarray(capacity, dtype=int)
    cap[cap < 0] = 0
    n_spots, n_cells = mat.shape
    mat_csc = mat.tocsc()
    fallback_map = {}
    if fallback_assignments:
        for cidx, sidx, _ in fallback_assignments:
            fallback_map[cidx] = sidx
    # 对齐 type_prior：按 spot 顺序、type 顺序重排
    type_cols = list(type_prior.columns)
    type_prior_aligned = type_prior.reindex(index=spot_ids).reindex(columns=type_cols)
    if type_prior_aligned.isna().all(axis=None):
        raise ValueError("type_prior 对齐后全部为 NaN，spot_id 或列名未匹配")
    prior_row_nonzero = (~(type_prior_aligned.fillna(0.0) == 0).all(axis=1)).mean()
    prior_row_entropy = []
    prior_filled = type_prior_aligned.fillna(0.0)
    for _, row in prior_filled.iterrows():
        p = row.to_numpy(dtype=float)
        p_sum = p.sum()
        if p_sum <= 0:
            prior_row_entropy.append(0.0)
        else:
            p = p / p_sum
            prior_row_entropy.append(-float(np.sum(p * np.log(p + eps))))
    prior_entropy_mean = float(np.mean(prior_row_entropy)) if prior_row_entropy else 0.0
    prior_entropy_max = float(np.max(prior_row_entropy)) if prior_row_entropy else 0.0
    prior_np = type_prior_aligned.fillna(0.0).to_numpy(dtype=float, copy=False)
    type_to_col = {t: i for i, t in enumerate(type_cols)}
    prior_candidate_topk = int(prior_candidate_topk or 0)
    prior_candidate_weight = float(1.0 if prior_candidate_weight is None else prior_candidate_weight)
    prior_topk_by_type: Optional[List[List[int]]] = None
    if lambda_prior and prior_candidate_topk > 0:
        k = min(prior_candidate_topk, n_spots)
        prior_topk_by_type = []
        for t_idx in range(len(type_cols)):
            col = prior_np[:, t_idx]
            if k >= n_spots:
                idx = np.argsort(-col)
            else:
                idx = np.argpartition(-col, k - 1)[:k]
                idx = idx[np.argsort(-col[idx])]
            prior_topk_by_type.append(idx.tolist())
    if prior_row_nonzero < min_prior_row_nonzero_ratio:
        raise ValueError(f"type_prior 非零行比例过低: {prior_row_nonzero:.3f} < {min_prior_row_nonzero_ratio}")

    candidates: List[List[Tuple[int, float]]] = []
    local_change = 0
    for c in range(n_cells):
        col = mat_csc.getcol(c)
        rows = col.indices
        data = col.data
        if len(rows) > 0:
            order = np.argsort(data)[::-1]
            if topk and topk > 0:
                order = order[:topk]
            rows = rows[order]
            data = data[order]
        cand_dict: Dict[int, float] = {}
        for r, v in zip(rows, data):
            cand_dict[r] = max(cand_dict.get(r, 0.0), v)
        if c in fallback_map:
            cand_dict.setdefault(fallback_map[c], 0.0)
        cand_list: List[Tuple[int, float]] = []
        type_col = type_to_col.get(cell_types[c], None)
        best_score_only = None
        if lambda_prior and prior_topk_by_type is not None and type_col is not None:
            for r in prior_topk_by_type[type_col]:
                if r not in cand_dict:
                    cand_dict[r] = prior_candidate_weight * float(prior_np[r, type_col])
        for r, v in cand_dict.items():
            prior_val = prior_np[r, type_col] if (type_col is not None and r < prior_np.shape[0]) else 0.0
            u = np.log(v + eps)
            if lambda_prior:
                u += float(lambda_prior) * np.log(prior_val + eps)
            cand_list.append((r, u))
            if best_score_only is None or v > best_score_only[1]:
                best_score_only = (r, v)
        # 候选为空时，尝试 fallback spot
        if not cand_list and c in fallback_map:
            r = fallback_map[c]
            prior_val = prior_np[r, type_col] if (type_col is not None and r < prior_np.shape[0]) else 0.0
            u = float(lambda_prior) * np.log(prior_val + eps) if lambda_prior else 0.0
            cand_list.append((r, u))
        cand_list.sort(key=lambda x: x[1], reverse=True)
        candidates.append(cand_list)
        if best_score_only and cand_list:
            if best_score_only[0] != cand_list[0][0]:
                local_change += 1

    assignment = [-1] * n_cells
    spot_pool: Dict[int, List[Tuple[int, float]]] = {i: [] for i in range(n_spots)}
    next_pos = [0] * n_cells
    free = [i for i in range(n_cells) if candidates[i]]
    while free:
        c = free.pop()
        cand = candidates[c]
        if next_pos[c] >= len(cand):
            continue
        s, util = cand[next_pos[c]]
        next_pos[c] += 1
        spot_pool[s].append((c, util))
        if len(spot_pool[s]) > cap[s]:
            spot_pool[s].sort(key=lambda x: x[1], reverse=True)
            dropped = spot_pool[s][cap[s]:]
            spot_pool[s] = spot_pool[s][: cap[s]]
            for dc, _ in dropped:
                free.append(dc)

    for s, pool in spot_pool.items():
        for c, util in pool:
            assignment[c] = s

    data = []
    rows = []
    cols = []
    for c, s in enumerate(assignment):
        if s >= 0:
            rows.append(s)
            cols.append(c)
            data.append(1.0)
    hard_mat = sparse.csr_matrix((data, (rows, cols)), shape=mat.shape)
    assigned_per_spot = np.asarray(hard_mat.sum(axis=1)).ravel()

    fallback_used = 0
    if fallback_map:
        for c, a in enumerate(assignment):
            if a == -1 and c in fallback_map:
                s = fallback_map[c]
                if assigned_per_spot[s] < cap[s]:
                    assignment[c] = s
                    assigned_per_spot[s] += 1
                    fallback_used += 1
    # 重新生成 hard 矩阵，保证包含 fallback
    data = []
    rows = []
    cols = []
    for c, s in enumerate(assignment):
        if s >= 0:
            rows.append(s)
            cols.append(c)
            data.append(1.0)
    hard_mat = sparse.csr_matrix((data, (rows, cols)), shape=mat.shape)
    assigned_per_spot = np.asarray(hard_mat.sum(axis=1)).ravel()
    assignment_arr = np.asarray(assignment, dtype=int)
    n_unassigned_before_repair = int((assignment_arr < 0).sum())

    # 修复：在 topK 偏好过窄时，延迟接受可能留下未分配 cell；用“按剩余容量+type_prior”填满剩余槽位
    repair_used = 0
    if n_unassigned_before_repair > 0:
        rem = cap - assigned_per_spot
        unassigned_cells = np.where(assignment_arr < 0)[0]
        for c in unassigned_cells:
            avail = np.where(rem > 0)[0]
            if avail.size == 0:
                break
            type_col = type_to_col.get(cell_types[c], None)
            if type_col is None:
                s = int(avail[0])
            else:
                col = prior_np[avail, type_col]
                s = int(avail[int(np.argmax(col))])
            assignment[c] = s
            rem[s] -= 1
            assigned_per_spot[s] += 1
            repair_used += 1
        # 重建 hard 矩阵（包含 repair）
        data = []
        rows = []
        cols = []
        for c, s in enumerate(assignment):
            if s >= 0:
                rows.append(s)
                cols.append(c)
                data.append(1.0)
        hard_mat = sparse.csr_matrix((data, (rows, cols)), shape=mat.shape)
        assigned_per_spot = np.asarray(hard_mat.sum(axis=1)).ravel()
        assignment_arr = np.asarray(assignment, dtype=int)

    n_unassigned = int((assignment_arr < 0).sum())
    cap_diff = assigned_per_spot - cap
    overflow = np.maximum(cap_diff, 0)
    underfill = np.maximum(-cap_diff, 0)
    meta = {
        "harden_method": "quota_matching",
        "harden_topk": topk,
        "harden_fallback_n_cells": fallback_used,
        "harden_repair_n_cells": repair_used,
        "n_unassigned_before_repair": n_unassigned_before_repair,
        "capacity_check": {
            "max_overflow": float(overflow.max()) if len(overflow) else 0.0,
            "max_underfill": float(underfill.max()) if len(underfill) else 0.0,
            "n_spots_overflow": int((overflow > 0).sum()) if len(overflow) else 0,
            "n_spots_underfill": int((underfill > 0).sum()) if len(underfill) else 0,
        },
        "n_unassigned_cells": n_unassigned,
        "n_cells_total": int(n_cells),
        "type_prior_row_index_type": str(type_prior.index.__class__.__name__),
        "type_prior_row_nonzero_ratio": float(prior_row_nonzero),
        "type_prior_row_entropy_mean": prior_entropy_mean,
        "type_prior_row_entropy_max": prior_entropy_max,
        "prior_local_top1_change_rate": float(local_change / n_cells) if n_cells > 0 else 0.0,
        "prior_candidate_topk": int(prior_candidate_topk),
        "prior_candidate_weight": float(prior_candidate_weight),
    }
    if n_unassigned > 0:
        meta["status"] = "failed_unassigned"
        meta["error_msg"] = f"harden left {n_unassigned}/{n_cells} cells unassigned"
        raise RuntimeError(meta["error_msg"])
    if overflow.max() > 0:
        meta["status"] = "failed_overflow"
        meta["error_msg"] = f"harden overflow detected max_overflow={overflow.max()}"
        raise RuntimeError(meta["error_msg"])
    if ablate_lambda:
        meta["ablation_mode"] = "lambda0"
    assignments_hard = [(c, s, 1.0) for c, s in enumerate(assignment) if s >= 0]
    return assignments_hard, hard_mat, meta


def _write_global_fraction(cell_types: List[str], type_order: List[str], out_path: Path):
    ser = pd.Series(cell_types, name="cell_type").value_counts(normalize=True)
    frac_row = ser.reindex(type_order).fillna(0).to_frame().T
    frac_row.index = ["global"]
    frac_row.to_csv(out_path)
    return out_path


def _ensure_type_prior_columns(type_prior: pd.DataFrame, type_cols: List[str]) -> pd.DataFrame:
    for col in type_cols:
        if col not in type_prior.columns:
            type_prior[col] = 0.0
    return type_prior[type_cols]


def _build_outputs(
    assignments: List[Tuple[int, int, float]],
    sc_index: List[str],
    st_index: List[str],
    cell_types: List[str],
    out_dir: Path,
    mode: str,
    type_order: Optional[List[str]] = None,
    hard_matrix: Optional[sparse.csr_matrix] = None,
    unique_cell_ids: Optional[List[str]] = None,
):
    if unique_cell_ids is not None and len(unique_cell_ids) != len(sc_index):
        raise ValueError(
            f"[outputs] unique_cell_ids length {len(unique_cell_ids)} != sc_index length {len(sc_index)}"
        )
    rows = []
    for cell_idx, spot_idx, score in assignments:
        r = {
            "cell_id": sc_index[cell_idx],
            "spot_id": st_index[spot_idx],
            "type": cell_types[cell_idx],
            "backend": "cytospace",
            "mode": mode,
            "assign_score": score,
        }
        if unique_cell_ids is not None:
            r["unique_cid"] = unique_cell_ids[cell_idx]
        rows.append(r)
    ca = pd.DataFrame(rows)
    ca.to_csv(out_dir / f"cell_assignment_{mode}.csv", index=False)

    if hard_matrix is not None:
        mat = hard_matrix if sparse.isspmatrix_csr(hard_matrix) else hard_matrix.tocsr()
    else:
        data = np.ones(len(assignments), dtype=float)
        row_ind = [spot_idx for _, spot_idx, _ in assignments]
        col_ind = [cell_idx for cell_idx, _, _ in assignments]
        mat = sparse.csr_matrix((data, (row_ind, col_ind)), shape=(len(st_index), len(sc_index)))
    # API 约定：cell×spot（行=cell，列=spot）
    sparse.save_npz(out_dir / f"cell_spot_matrix_{mode}.npz", mat.T.tocsr())

    if type_order is None:
        type_order = sorted(list(dict.fromkeys(cell_types)))
    df = pd.DataFrame(0.0, index=st_index, columns=type_order)
    mat_csr = mat if sparse.isspmatrix_csr(mat) else mat.tocsr()
    cell_type_arr = np.array([cell_types[i] for i in range(len(sc_index))])
    for t in type_order:
        mask = cell_type_arr == t
        if mask.any():
            df[t] = mat_csr[:, mask].sum(axis=1).A1
    df = df.div(df.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    df.index.name = "spot_id"
    df.reset_index().to_csv(out_dir / f"spot_type_fraction_{mode}.csv", index=False)


class CytoSPACEBackend(MappingBackend):
    def __init__(self):
        super().__init__(name="cytospace")

    def run_baseline(self, stage1_dir: Path, out_dir: Path, config: Dict[str, Any]) -> None:
        t0 = time.time()
        _ensure_dir(out_dir)
        status = "success"
        error_msg = None
        type_order = None
        config_validation = None
        capacity_audit = None
        config_effective_subset = None
        try:
            config, config_validation = _validate_and_resolve_config(config, context="baseline")
            config_effective_subset = {k: config.get(k) for k in _CONFIG_EFFECTIVE_KEYS if k in config}

            sc_expr, st_expr, sc_meta, st_coords, type_col = _load_stage1(stage1_dir)
            _assert_unique_index(sc_expr, "sc_expr")
            _assert_unique_index(st_expr, "st_expr")
            _assert_unique_index(sc_meta, "sc_meta")
            _assert_unique_index(st_coords, "st_coords")

            mode = str(config.get("mode") or "")
            run_id = config.get("run_id")
            variant = config.get("variant")
            if mode not in ("baseline", "plus"):
                raise ValueError(f"[baseline] invalid mode={mode!r}; expected 'baseline' or 'plus'")
            if mode != "baseline":
                raise ValueError(f"[baseline] invalid mode={mode!r}; expected 'baseline'")
            if run_id is None:
                raise ValueError("[baseline] missing config.run_id")
            run_id = str(run_id)
            if not run_id.startswith("baseline"):
                raise ValueError(f"[baseline] invalid run_id={run_id!r}; expected startswith 'baseline'")
            if out_dir.name != run_id:
                raise ValueError(f"[baseline] out_dir.name={out_dir.name!r} != run_id={run_id!r}")

            spot_ids = list(st_expr.index)
            capacity_raw, cps_source, norm_val = _compute_cells_per_spot(st_coords, config)
            st_coords, capacity_raw, _, spot_alignment_audit = _align_spot_inputs(
                spot_ids=spot_ids,
                st_coords=st_coords,
                capacity=capacity_raw,
                type_prior_raw=None,
                context="baseline",
            )
            capacity, cap_audit = _resolve_cells_per_spot(capacity_raw, config, context="baseline")
            capacity_audit = {
                **cap_audit,
                "source_config": config.get("cells_per_spot_source"),
                "source_resolved": cps_source,
                "umi_to_cell_norm": norm_val,
                "default_cells_per_spot": config.get("default_cells_per_spot"),
            }
            genes_use = sorted(set(sc_expr.columns) & set(st_expr.columns))
            cell_types = sc_meta["type_col"].tolist()
            type_order = sorted(list(dict.fromkeys(cell_types)))

            with tempfile.TemporaryDirectory() as tmp:
                tmp_dir = Path(tmp)
                os.environ.setdefault("PYTHON", sys.executable)
                os.environ.setdefault("PYTHON3", sys.executable)
                ct_frac_path = _write_global_fraction(cell_types, type_order, tmp_dir / "ct_fraction.csv")
                sc_path, ct_path, st_path, coord_path, cps_path = _write_cytospace_inputs(
                    sc_expr,
                    sc_meta,
                    st_expr,
                    st_coords,
                    cell_types,
                    genes_use,
                    capacity,
                    spot_ids,
                    tmp_dir,
                    prefix="base_",
                )
                out_tmp = tmp_dir / "cyto_out"
                out_tmp.mkdir(exist_ok=True)
                solver_method = config.get("solver_method", "lap_CSPR")
                main_cytospace(
                    scRNA_path=str(sc_path),
                    cell_type_path=str(ct_path),
                    n_cells_per_spot_path=str(cps_path),
                    st_cell_type_path=None,
                    cell_type_fraction_estimation_path=str(ct_frac_path),
                    spaceranger_path=None,
                    st_path=str(st_path),
                    coordinates_path=str(coord_path),
                    output_folder=str(out_tmp),
                    output_prefix="",
                    mean_cell_numbers=int(max(1, np.mean(capacity))),
                    downsample_off=True,
                    scRNA_max_transcripts_per_cell=1500,
                    solver_method=solver_method,
                    distance_metric="Pearson_correlation",
                    sampling_method="duplicates",
                    single_cell=False,
                    number_of_selected_spots=10000,
                    sampling_sub_spots=False,
                    number_of_selected_sub_spots=10000,
                    number_of_processors=1,
                    seed=config.get("seed", 1),
                    plot_off=True,
                    geometry="honeycomb",
                    max_num_cells_plot=50000,
                    num_column=3,
                )
                assigned_path = out_tmp / "assigned_locations.csv"
                if not assigned_path.exists():
                    raise FileNotFoundError("CytoSPACE 输出缺少 assigned_locations.csv")
                assigned = pd.read_csv(assigned_path)
                cell_ids, pool_types, pool_orig_ids, assignments = _parse_assigned_locations(
                    assigned, list(st_expr.index), require_original_cid=True
                )
                try:
                    assigned[["UniqueCID", "OriginalCID", "SpotID", "CellType"]].to_csv(
                        out_dir / "cell_pool_map_baseline.csv", index=False
                    )
                except Exception:
                    pass

                if not pool_orig_ids:
                    raise ValueError("[baseline] assigned_locations 缺少 OriginalCID，无法与 Stage1/SimGen Query 真值对齐")
                mat = _assignments_to_matrix(assignments, len(st_expr.index), len(cell_ids)).tocsr()
                _build_outputs(
                    assignments,
                    [str(x) for x in pool_orig_ids],
                    list(st_expr.index),
                    pool_types,
                    out_dir,
                    mode=mode,
                    type_order=type_order,
                    hard_matrix=mat,
                    unique_cell_ids=[str(x) for x in cell_ids],
                )
        except Exception as e:
            import traceback
            status = "failed"
            error_msg = f"{e}\n{traceback.format_exc()}"
            for fname in ["cell_assignment_baseline.csv", "cell_spot_matrix_baseline.npz", "spot_type_fraction_baseline.csv"]:
                try:
                    (out_dir / fname).unlink()
                except Exception:
                    pass
        meta = {
            "backend": "cytospace",
            "mode": config.get("mode"),
            "run_id": config.get("run_id"),
            "variant": config.get("variant"),
            "out_dir_name": out_dir.name,
            "sample": config.get("sample"),
            "seed": config.get("seed"),
            **_module_fingerprint(),
            "runner_file": config.get("runner_file"),
            "runner_sha1": (str(config.get("runner_sha1")).lower() if config.get("runner_sha1") else None),
            "svg_refine_lambda": None,
            "feature_set": "all_genes",
            "type_set": "orig_type",
            "weight_field": None,
            "stage1_dir": str(stage1_dir),
            "stage2_dir": None,
            "stage3_dir": None,
            "config_id": config.get("config_id"),
            "project_config_path": config.get("project_config_path"),
            "dataset_config_path": config.get("dataset_config_path"),
            "config_validation": locals().get("config_validation", None),
            "config_effective_subset": locals().get("config_effective_subset", None),
            "capacity_audit": locals().get("capacity_audit", None),
            "cell_id_space": "original_cid",
            "cell_instance_id_column": "unique_cid",
            "cells_per_spot_source": locals().get("cps_source", None),
            "umi_to_cell_norm": locals().get("norm_val", None),
            "has_capacity_constraint": status == "success",
            "cost_cfg": {"solver_method": config.get("solver_method", "lap_CSPR")},
            "refine_cfg": {},
            "n_spots_stage1": len(st_expr.index),
            "n_spots_coords": len(st_coords.index),
            "n_spots_intersection": len(set(st_expr.index) & set(st_coords.index)),
            "n_cells_sc_expr": len(sc_expr.index),
            "n_cells_sc_meta": len(sc_meta.index),
            "n_cells_intersection": len(set(sc_expr.index) & set(sc_meta.index)),
            "n_cells_pool": len(locals().get("cell_ids", [])) if status == "success" else None,
            "type_columns": type_order,
            "spot_alignment_audit": locals().get("spot_alignment_audit", None),
            "runtime_sec": time.time() - t0,
            "status": status,
            "error_msg": error_msg,
            "resolved_mapping_config": config,
        }
        _write_meta(out_dir, meta)

    def run_plus(
        self,
        stage1_dir: Path,
        stage2_dir: Path,
        stage3_dir: Path,
        out_dir: Path,
        config: Dict[str, Any],
    ) -> None:
        t0 = time.time()
        _ensure_dir(out_dir)
        status = "success"
        error_msg = None
        config_validation = None
        capacity_audit = None
        config_effective_subset = None
        try:
            config, config_validation = _validate_and_resolve_config(config, context="plus")
            config_effective_subset = {k: config.get(k) for k in _CONFIG_EFFECTIVE_KEYS if k in config}

            sc_expr, st_expr, sc_meta, st_coords, type_col = _load_stage1(stage1_dir)
            _assert_unique_index(sc_expr, "sc_expr")
            _assert_unique_index(st_expr, "st_expr")
            _assert_unique_index(sc_meta, "sc_meta")
            _assert_unique_index(st_coords, "st_coords")

            mode = str(config.get("mode") or "")
            run_id = config.get("run_id")
            variant = config.get("variant")
            if mode not in ("baseline", "plus"):
                raise ValueError(f"[plus] invalid mode={mode!r}; expected 'baseline' or 'plus'")
            if mode != "plus":
                raise ValueError(f"[plus] invalid mode={mode!r}; expected 'plus'")
            if run_id is None:
                raise ValueError("[plus] missing config.run_id")
            run_id = str(run_id)
            if not run_id.startswith("plus"):
                raise ValueError(f"[plus] invalid run_id={run_id!r}; expected startswith 'plus'")
            if out_dir.name != run_id:
                raise ValueError(f"[plus] out_dir.name={out_dir.name!r} != run_id={run_id!r}")

            spot_ids = list(st_expr.index)
            plugin_genes, weight_map, weight_field, gene_weights = _load_stage2(stage2_dir)
            relabel, type_prior = _load_stage3(stage3_dir)
            stage2_gene_weights_sha1 = _sha1_file(stage2_dir / "gene_weights.csv")
            plugin_genes_sha1 = _sha1_file(stage2_dir / "plugin_genes.txt")
            if len(plugin_genes) == 0:
                raise ValueError("plugin_genes 为空")

            relabel_map = dict(zip(relabel["cell_id"], relabel["plugin_type"]))
            eps = float(config.get("eps", 1e-8) or 1e-8)
            type_cols = sorted(list(dict.fromkeys(relabel_map.get(cid, "Unknown_sc_only") for cid in sc_expr.index)))
            cell_types = [relabel_map.get(cid, "Unknown_sc_only") for cid in sc_expr.index]
            prior_cols_raw = set(type_prior.columns)
            sc_types = set(type_cols)
            missing_type_cols = sorted(sc_types - prior_cols_raw)
            n_cells_missing_type = sum(1 for t in cell_types if t in missing_type_cols)
            ratio_cells_missing_type = n_cells_missing_type / max(1, len(cell_types))
            if ratio_cells_missing_type > float(config.get("max_cells_missing_type_prior_ratio", 0.0)):
                raise ValueError(f"{ratio_cells_missing_type:.3f} 比例的细胞类型缺失在 type_prior 列中，超过阈值")

            
            genes_use = sorted(set(sc_expr.columns) & set(st_expr.columns))
            if not genes_use:
                raise ValueError("sc/st ??????")
            gene_to_idx = {g: i for i, g in enumerate(genes_use)}
            plugin_n = max(1, len(plugin_genes))
            plugin_in_sc = [g for g in plugin_genes if g in sc_expr.columns]
            plugin_in_st = [g for g in plugin_genes if g in st_expr.columns]
            sc_cov = len(plugin_in_sc) / plugin_n
            st_cov = len(plugin_in_st) / plugin_n
            plugin_overlap_genes = [g for g in plugin_genes if g in gene_to_idx]
            plugin_overlap_ratio = len(plugin_overlap_genes) / plugin_n

            truth_pass_map = None
            if "truth_pass" in gene_weights.columns:
                truth_pass_map = dict(zip(gene_weights["gene"], gene_weights["truth_pass"]))
            spot_weight_candidate_mask = None
            if config.get("spot_weight_truth_filter", False) and truth_pass_map is not None:
                spot_weight_candidate_mask = np.array(
                    [bool(truth_pass_map.get(g, 0)) for g in genes_use],
                    dtype=bool,
                )

            w_full = np.ones(len(genes_use), dtype=float)
            weighted_genes = []
            for g in plugin_overlap_genes:
                idx = gene_to_idx[g]
                w_val = float(weight_map.get(g, 1.0))
                if not np.isfinite(w_val):
                    w_val = 1.0
                w_full[idx] = w_val
                if abs(w_val - 1.0) > 1e-12:
                    weighted_genes.append(g)

            w_min = float(np.min(w_full))
            w_p50 = float(np.percentile(w_full, 50))
            w_p95 = float(np.percentile(w_full, 95))
            w_max = float(np.max(w_full))
            weighted_ratio = len(weighted_genes) / max(len(genes_use), 1)

            svg_fallback_reasons = []
            if plugin_overlap_ratio < float(config.get("min_gene_overlap_ratio", 0.0)):
                svg_fallback_reasons.append(f"plugin_overlap_low:{plugin_overlap_ratio:.3f}")
            if weighted_ratio > 0.2:
                svg_fallback_reasons.append(f"weighted_ratio>0.2:{weighted_ratio:.3f}")
            if abs(w_p50 - 1.0) > 1e-8:
                svg_fallback_reasons.append(f"p50!=1:{w_p50:.6f}")
            if w_min < 1.0 - 1e-8:
                svg_fallback_reasons.append(f"weight_below_1:{w_min:.6f}")

            svg_fallback = len(svg_fallback_reasons) > 0
            w_stats_raw = {"min": w_min, "p50": w_p50, "p95": w_p95, "max": w_max}
            if svg_fallback:
                w_full = np.ones(len(genes_use), dtype=float)
                weighted_genes_applied = []
            else:
                weighted_genes_applied = weighted_genes
            weighted_ratio_applied = len(weighted_genes_applied) / max(len(genes_use), 1)
            w_stats = {
                "min": float(np.min(w_full)),
                "p50": float(np.percentile(w_full, 50)),
                "p95": float(np.percentile(w_full, 95)),
                "max": float(np.max(w_full)),
            }

            capacity_raw, cps_source, norm_val = _compute_cells_per_spot(st_coords, config)
            st_coords, capacity_raw, type_prior_raw, spot_alignment_audit = _align_spot_inputs(
                spot_ids=spot_ids,
                st_coords=st_coords,
                capacity=capacity_raw,
                type_prior_raw=type_prior,
                context="plus_svg_type",
            )
            capacity, cap_audit = _resolve_cells_per_spot(capacity_raw, config, context="plus")
            capacity_audit = {
                **cap_audit,
                "source_config": config.get("cells_per_spot_source"),
                "source_resolved": cps_source,
                "umi_to_cell_norm": norm_val,
                "default_cells_per_spot": config.get("default_cells_per_spot"),
            }
            if type_prior_raw is None:
                raise ValueError("type_prior_raw 缺失（align_spot_inputs 返回 None）")
            type_prior_raw = _ensure_type_prior_columns(type_prior_raw, type_cols)
            prior_intersection = len(set(spot_ids) & set(type_prior_raw.index))
            missing_prior_spots = list(set(spot_ids) - set(type_prior_raw.index))
            missing_prior_examples = missing_prior_spots[:5]
            cost_expr_guard_reason = None
            norm_mode = str(config.get("cost_expr_norm") or "none").lower()
            sample_id = str(config.get("sample") or "")
            is_sim = bool((ROOT / "data" / "sim" / sample_id).exists())
            norm_cfg = config
            if is_sim and norm_mode == "st_zscore":
                norm_cfg = dict(config)
                norm_cfg["cost_expr_norm"] = "none"
                cost_expr_guard_reason = "sim_zscore_forbidden"
            sc_cost, st_cost, cost_expr_audit = _apply_cost_expr_norm(
                sc_expr,
                st_expr,
                genes_use,
                norm_cfg,
                eps=float(norm_cfg.get("cost_expr_norm_eps") or 1e-8),
            )
            if cost_expr_guard_reason:
                cost_expr_audit["guard_reason"] = cost_expr_guard_reason
            spot_weight_matrix, spot_weight_stats = _compute_spot_weight_matrix(
                st_cost,
                st_coords,
                genes_use,
                config,
                eps=float(config.get("cost_expr_norm_eps") or 1e-8),
                candidate_mask=spot_weight_candidate_mask,
            )
            refine_lambda = float(config.get("svg_refine_lambda", 0.0) or 0.0)
            refine_k = int(config.get("svg_refine_k", 8))
            distance_metric = config.get("distance_metric", "Pearson_correlation")
            lambda_prior = float(config.get("lambda_prior", 1.0))
            effective_lambda_refine = lambda_prior if config.get("type_prior_apply_refine", True) else 0.0
            effective_lambda_harden = lambda_prior if config.get("type_prior_apply_harden", True) else 0.0
            type_prior_eps = eps

            with tempfile.TemporaryDirectory() as tmp:
                tmp_dir = Path(tmp)
                os.environ.setdefault("PYTHON", sys.executable)
                os.environ.setdefault("PYTHON3", sys.executable)
                sc_weighted = sc_cost.copy()
                st_weighted = st_cost.copy()
                sc_vals = sc_cost[genes_use].to_numpy(dtype=float) * w_full
                st_vals = st_cost[genes_use].to_numpy(dtype=float) * w_full
                if spot_weight_matrix is not None:
                    st_vals = st_vals * spot_weight_matrix
                sc_weighted.loc[:, genes_use] = sc_vals
                st_weighted.loc[:, genes_use] = st_vals
                sc_path, ct_path, st_path, coord_path, cps_path = _write_cytospace_inputs(
                    sc_weighted,
                    sc_meta,
                    st_weighted,
                    st_coords,
                    cell_types,
                    genes_use,
                    capacity,
                    spot_ids,
                    tmp_dir,
                    prefix="plus_",
                )
                # NOTE (Stage5 truth alignment):
                # CytoSPACE samples a "cell pool" by the provided cell-type fractions. If we derive those
                # fractions from ST priors, it can over-demand rare types and trigger duplicate sampling,
                # breaking one-to-one Query cell_id alignment with SimGen truth. We therefore sample by the
                # Query sc composition (plugin_type) and apply type_prior only in our own refine/harden.
                st_ct_frac_path = _write_global_fraction(cell_types, type_cols, tmp_dir / "ct_fraction.csv")
                # spot×type -> 全局 1×type（喂 CytoSPACE）
                row_sum = type_prior_raw.sum(axis=1).replace(0, eps)
                type_prior_norm = type_prior_raw.div(row_sum, axis=0)

                out_tmp = tmp_dir / "cyto_out"
                out_tmp.mkdir(exist_ok=True)
                solver_method = config.get("solver_method", "lap_CSPR")
                main_cytospace(
                    scRNA_path=str(sc_path),
                    cell_type_path=str(ct_path),
                    n_cells_per_spot_path=str(cps_path),
                    st_cell_type_path=None,
                    cell_type_fraction_estimation_path=str(st_ct_frac_path),
                    spaceranger_path=None,
                    st_path=str(st_path),
                    coordinates_path=str(coord_path),
                    output_folder=str(out_tmp),
                    output_prefix="",
                    mean_cell_numbers=int(max(1, np.mean(capacity))),
                    downsample_off=True,
                    scRNA_max_transcripts_per_cell=1500,
                    solver_method=solver_method,
                    distance_metric=distance_metric,
                    sampling_method="duplicates",
                    single_cell=False,
                    number_of_selected_spots=10000,
                    sampling_sub_spots=False,
                    number_of_selected_sub_spots=10000,
                    number_of_processors=1,
                    seed=config.get("seed", 1),
                    plot_off=True,
                    geometry="honeycomb",
                    max_num_cells_plot=50000,
                    num_column=3,
                )
                assigned_path = out_tmp / "assigned_locations.csv"
                if not assigned_path.exists():
                    raise FileNotFoundError("CytoSPACE 输出缺少 assigned_locations.csv")
                assigned = pd.read_csv(assigned_path)
                cell_ids, pool_types, pool_orig_ids, assignments = _parse_assigned_locations(
                    assigned, list(st_expr.index), require_original_cid=True
                )
                if pool_orig_ids is None:
                    raise ValueError("pool_orig_ids 缺失")
                mat0 = _assignments_to_matrix(assignments, len(st_expr.index), len(cell_ids)).tocsr()
                mat = mat0
                capacity_arr = np.asarray(capacity, dtype=int)
                cap_sum = int(capacity_arr.sum())
                if cap_sum != len(cell_ids):
                    raise ValueError(f"cells_per_spot 总和 {cap_sum} 与 CytoSPACE cell pool 大小 {len(cell_ids)} 不一致")
                post_adjust = bool(config.get("post_assign_adjust", True))
                knn_mode = None
                if post_adjust and refine_lambda > 0:
                    missing_pool = [cid for cid in pool_orig_ids if cid not in sc_expr.index]
                    if missing_pool:
                        raise ValueError(f"CytoSPACE cell pool ?????????sc_expr ???????????? OriginalCID????????? {missing_pool[:5]}")
                    sc_svg = sc_expr.loc[pool_orig_ids, genes_use].to_numpy(dtype=float) * w_full
                    st_svg = st_expr[genes_use].to_numpy(dtype=float) * w_full
                    mat, knn_mode = _refine_spot_cell_matrix_svg(
                        mat=mat,
                        st_coords=st_coords,
                        st_svg=st_svg,
                        sc_svg=sc_svg,
                        lambda_base=refine_lambda,
                        k=refine_k,
                        eps=eps,
                        batch_size=config.get("svg_refine_batch_size"),
                        knn_metric=config.get("knn_metric", "euclidean"),
                        knn_block_size=int(config.get("knn_block_size", 1024) or 1024),
                        knn_max_dense_n=int(config.get("knn_max_dense_n", 5000) or 5000),
                    )
                    refine_used = True
                else:
                    refine_used = False
                # ???????????? refine???per-spot soft adjust???
                if post_adjust and effective_lambda_refine > 0 and config.get("type_prior_apply_refine", True):
                    missing_pool_types = sorted(set(pool_types) - set(type_cols))
                    if missing_pool_types:
                        raise ValueError(f"cell pool ?????? type_cols ??????????????????: {missing_pool_types[:5]}")
                    mat = _type_prior_refine(mat, pool_types, type_cols, type_prior_norm, effective_lambda_refine, eps=type_prior_eps)
                hard_topk = int(config.get("harden_topk", 5) or 5)
                prior_candidate_topk = int(config.get("prior_candidate_topk", 0) or 0)
                prior_candidate_weight = float(config.get("prior_candidate_weight", 1.0) or 1.0)
                if post_adjust:
                    assignments_hard, hard_mat, hard_meta = _harden_assignment_quota_matching(
                        mat=mat,
                        capacity=capacity_arr,
                        cell_types=pool_types,
                        type_prior=type_prior_norm,
                        spot_ids=list(st_expr.index),
                        lambda_prior=effective_lambda_harden,
                        eps=eps,
                        topk=hard_topk,
                        prior_candidate_topk=prior_candidate_topk,
                        prior_candidate_weight=prior_candidate_weight,
                        fallback_assignments=assignments,
                        ablate_lambda=False,
                        min_prior_row_nonzero_ratio=float(config.get("min_prior_row_nonzero_ratio", 0.0)),
                    )
                else:
                    assignments_hard = assignments
                    hard_mat = mat0
                    hard_meta = {
                        "harden_method": "skipped_post_adjust",
                        "harden_topk": None,
                        "capacity_check": None,
                        "harden_fallback_n_cells": None,
                        "harden_repair_n_cells": None,
                        "n_unassigned_before_repair": None,
                        "n_unassigned_cells": None,
                        "type_prior_row_nonzero_ratio": None,
                        "type_prior_row_entropy_mean": None,
                        "type_prior_row_entropy_max": None,
                    }
                if post_adjust and config.get("prior_ablation_enabled", False):
                    # ablation 仅用于证据链；即使失败也不应影响 plus 主输出
                    hard_meta["prior_ablation_enabled"] = True
                    try:
                        ablate_assign, ablate_mat, ablate_meta = _harden_assignment_quota_matching(
                            mat=mat,
                            capacity=capacity_arr,
                            cell_types=pool_types,
                            type_prior=type_prior_norm,
                            spot_ids=list(st_expr.index),
                            lambda_prior=0.0,
                            eps=eps,
                            topk=hard_topk,
                            prior_candidate_topk=prior_candidate_topk,
                            prior_candidate_weight=prior_candidate_weight,
                            fallback_assignments=assignments,
                            ablate_lambda=True,
                            min_prior_row_nonzero_ratio=float(config.get("min_prior_row_nonzero_ratio", 0.0)),
                        )
                        hard_meta["prior_ablation_meta"] = ablate_meta
                        try:
                            change_rate, audit = _compute_ablation_change_rate_checked(
                                assignments_hard,
                                ablate_assign,
                                context="prior_ablation_change_rate",
                            )
                            hard_meta["prior_ablation_status"] = "ok"
                            hard_meta["prior_ablation_change_rate"] = change_rate
                            hard_meta["prior_ablation_audit"] = audit
                        except Exception as e:
                            hard_meta["prior_ablation_status"] = "invalid"
                            hard_meta["prior_ablation_change_rate"] = None
                            hard_meta["prior_ablation_error_reason"] = str(e)
                            hard_meta["prior_ablation_audit"] = {
                                "context": "prior_ablation_change_rate",
                                "n_base": len(assignments_hard),
                                "n_ablate": len(ablate_assign),
                            }
                    except Exception as e:
                        hard_meta["prior_ablation_status"] = "invalid"
                        hard_meta["prior_ablation_change_rate"] = None
                        hard_meta["prior_ablation_error_reason"] = f"ablation_harden_failed: {e}"
                        hard_meta["prior_ablation_audit"] = {
                            "context": "prior_ablation_change_rate",
                            "n_base": len(assignments_hard),
                            "n_ablate": None,
                        }
                frac_before, audit_b = _spot_fraction_from_mat(
                    mat0,
                    pool_types,
                    list(st_expr.index),
                    type_cols,
                    context="prior_effect_fraction_delta/before",
                )
                frac_after, audit_a = _spot_fraction_from_mat(
                    hard_mat,
                    pool_types,
                    list(st_expr.index),
                    type_cols,
                    context="prior_effect_fraction_delta/after",
                )
                frac_before_norm = frac_before.div(frac_before.sum(axis=1).replace(0, 1.0), axis=0)
                frac_after_norm = frac_after.div(frac_after.sum(axis=1).replace(0, 1.0), axis=0)
                per_spot_l1 = np.abs(frac_before_norm.values - frac_after_norm.values).sum(axis=1)
                hard_meta["prior_effect_fraction_delta"] = float(per_spot_l1.sum())
                hard_meta["prior_effect_fraction_delta_mean_per_spot"] = float(per_spot_l1.mean())
                hard_meta["prior_effect_fraction_delta_audit_before"] = audit_b
                hard_meta["prior_effect_fraction_delta_audit_after"] = audit_a
                try:
                    assigned[["UniqueCID", "OriginalCID", "SpotID", "CellType"]].to_csv(
                        out_dir / "cell_pool_map_plus.csv", index=False
                    )
                except Exception:
                    pass

                if not pool_orig_ids:
                    raise ValueError("[plus] assigned_locations 缺少 OriginalCID，无法与 Stage1/SimGen Query 真值对齐")
                # ⚠️cell_id 统一用 OriginalCID；UniqueCID 作为额外列 unique_cid 供审计
                abstain_unknown = bool(config.get("abstain_unknown_sc_only", False)) if post_adjust else False
                unknown_label = "Unknown_sc_only"
                unknown_mask = [t == unknown_label for t in pool_types]
                n_unknown_pool = int(sum(unknown_mask))
                n_unknown_dropped = 0
                if abstain_unknown and n_unknown_pool:
                    keep_indices = [i for i, is_unknown in enumerate(unknown_mask) if not is_unknown]
                    n_unknown_dropped = int(len(pool_types) - len(keep_indices))
                    if keep_indices:
                        keep_set = set(keep_indices)
                        index_map = {old_i: new_i for new_i, old_i in enumerate(keep_indices)}
                        assignments_hard = [
                            (index_map[cell_i], spot_i, score)
                            for cell_i, spot_i, score in assignments_hard
                            if cell_i in keep_set
                        ]
                        hard_mat = hard_mat[:, keep_indices].tocsr()
                        pool_types = [pool_types[i] for i in keep_indices]
                        pool_orig_ids = [pool_orig_ids[i] for i in keep_indices]
                        cell_ids = [cell_ids[i] for i in keep_indices]
                    else:
                        assignments_hard = []
                        hard_mat = hard_mat[:, :0].tocsr()
                        pool_types = []
                        pool_orig_ids = []
                        cell_ids = []
                hard_meta["abstain_unknown_sc_only"] = abstain_unknown
                hard_meta["abstain_unknown_label"] = unknown_label
                hard_meta["abstain_unknown_pool_count"] = n_unknown_pool
                hard_meta["abstain_unknown_dropped"] = n_unknown_dropped
                hard_meta["abstain_unknown_remaining"] = len(pool_orig_ids)
                _build_outputs(
                    assignments_hard,
                    [str(x) for x in pool_orig_ids],
                    list(st_expr.index),
                    pool_types,
                    out_dir,
                    mode=mode,
                    type_order=type_cols,
                    hard_matrix=hard_mat,
                    unique_cell_ids=[str(x) for x in cell_ids],
                )
        except Exception as e:
            import traceback
            status = "failed"
            error_msg = f"{e}\n{traceback.format_exc()}"
            for fname in ["cell_assignment_plus.csv", "cell_spot_matrix_plus.npz", "spot_type_fraction_plus.csv"]:
                try:
                    (out_dir / fname).unlink()
                except Exception:
                    pass
        meta = {
            "backend": "cytospace",
            "mode": config.get("mode"),
            "run_id": config.get("run_id"),
            "variant": config.get("variant"),
            "out_dir_name": out_dir.name,
            "sample": config.get("sample"),
            "seed": config.get("seed"),
            **_module_fingerprint(),
            "runner_file": config.get("runner_file"),
            "runner_sha1": (str(config.get("runner_sha1")).lower() if config.get("runner_sha1") else None),
            "svg_refine_lambda": config.get("svg_refine_lambda"),
            "feature_set": "all_common_genes",
            "type_set": "plugin_type",
            "weight_field": locals().get("weight_field"),
            "stage1_dir": str(stage1_dir),
            "stage2_dir": str(stage2_dir),
            "stage3_dir": str(stage3_dir),
            "stage2_gene_weights_sha1": locals().get("stage2_gene_weights_sha1"),
            "plugin_genes_sha1": locals().get("plugin_genes_sha1"),
            "n_genes_use": len(genes_use) if "genes_use" in locals() else None,
            "n_weighted_genes": len(weighted_genes_applied) if "weighted_genes_applied" in locals() else None,
            "weighted_gene_examples": (locals().get("weighted_genes_applied") or [])[:10],
            "weighted_gene_ratio": locals().get("weighted_ratio_applied"),
            "w_stats": locals().get("w_stats"),
            "w_stats_raw": locals().get("w_stats_raw"),
            "svg_fallback": locals().get("svg_fallback"),
            "svg_fallback_reason": locals().get("svg_fallback_reasons"),
            "plugin_overlap_ratio": locals().get("plugin_overlap_ratio"),
            "cost_expr_norm": locals().get("cost_expr_audit", None),
            "spot_weight_stats": locals().get("spot_weight_stats", None),
            "spot_specific_basin": config.get("spot_specific_basin"),
            "spot_specific_threshold_kappa": config.get("spot_specific_threshold_kappa"),
            "norm_applied": (locals().get("cost_expr_audit") or {}).get("norm_applied"),
            "norm_fallback": (locals().get("cost_expr_audit") or {}).get("norm_fallback"),
            "norm_fallback_reason": (locals().get("cost_expr_audit") or {}).get("norm_fallback_reason"),
            "norm_reason": (locals().get("cost_expr_audit") or {}).get("norm_reason"),
            "config_id": config.get("config_id"),
            "project_config_path": config.get("project_config_path"),
            "dataset_config_path": config.get("dataset_config_path"),
            "config_validation": locals().get("config_validation", None),
            "config_effective_subset": locals().get("config_effective_subset", None),
            "capacity_audit": locals().get("capacity_audit", None),
            "cell_id_space": "original_cid",
            "cell_instance_id_column": "unique_cid",
            "cells_per_spot_source": locals().get("cps_source", None),
            "umi_to_cell_norm": locals().get("norm_val", None),
            "has_capacity_constraint": status == "success",
            "cost_cfg": {
                "lambda_prior": config.get("lambda_prior", 1.0),
                "solver_method": config.get("solver_method", "lap_CSPR"),
                "distance_metric": locals().get("distance_metric", None),
            },
            "refine_cfg": {
                "refine_used": locals().get("refine_used", False),
                "svg_refine_lambda": locals().get("refine_lambda", 0.0),
                "svg_refine_k": locals().get("refine_k", None),
                "eps": config.get("eps", None),
                "svg_refine_mode": "error_weighted_svg_loss",
                "lambda_prior": locals().get("lambda_prior", None),
                "knn_block_size": config.get("knn_block_size", None),
                "knn_max_dense_n": config.get("knn_max_dense_n", None),
                "knn_metric": config.get("knn_metric", None),
                "knn_mode": locals().get("knn_mode", None),
                "harden_method": locals().get("hard_meta", {}).get("harden_method"),
                "harden_topk": locals().get("hard_meta", {}).get("harden_topk"),
                "type_prior_apply_refine": config.get("type_prior_apply_refine", True),
                "type_prior_apply_harden": config.get("type_prior_apply_harden", True),
                "svg_refine_batch_size": config.get("svg_refine_batch_size", None),
                "effective_lambda_prior_refine": locals().get("effective_lambda_refine", None),
                "effective_lambda_prior_harden": locals().get("effective_lambda_harden", None),
                "prior_ablation_enabled": config.get("prior_ablation_enabled", False),
                "prior_ablation_change_rate": locals().get("hard_meta", {}).get("prior_ablation_change_rate"),
                "prior_ablation_status": locals().get("hard_meta", {}).get("prior_ablation_status"),
                "prior_ablation_error_reason": locals().get("hard_meta", {}).get("prior_ablation_error_reason"),
            },
            "type_prior_mode": "global_mean_from_spot_matrix",
            "type_prior_mode_for_cytospace": "global_mean_1xT",
            "type_prior_mode_for_refine": "per_spot_SxT",
            "type_columns": locals().get("type_cols", None),
            "spot_alignment_audit": locals().get("spot_alignment_audit", None),
            "n_spots_stage1": len(st_expr.index),
            "n_spots_coords": len(st_coords.index),
            "n_spots_type_prior": len(type_prior.index) if "type_prior" in locals() else None,
            "n_spots_intersection": len(set(st_expr.index) & set(st_coords.index)),
            "n_spots_expr_type_prior_intersection": prior_intersection,
            "example_missing_in_type_prior": missing_prior_examples,
            "missing_type_cols": locals().get("missing_type_cols", []),
            "n_cells_missing_type_col": locals().get("n_cells_missing_type", 0),
            "ratio_cells_missing_type_col": locals().get("ratio_cells_missing_type", 0.0),
            "n_cells_sc_expr": len(sc_expr.index),
            "n_cells_sc_meta": len(sc_meta.index),
            "n_cells_intersection": len(set(sc_expr.index) & set(sc_meta.index)),
            "n_cells_pool": len(locals().get("cell_ids", [])) if status == "success" else None,
            "capacity_check": locals().get("hard_meta", {}).get("capacity_check"),
            "harden_fallback_n_cells": locals().get("hard_meta", {}).get("harden_fallback_n_cells"),
            "harden_repair_n_cells": locals().get("hard_meta", {}).get("harden_repair_n_cells"),
            "n_unassigned_before_repair": locals().get("hard_meta", {}).get("n_unassigned_before_repair"),
            "n_unassigned_cells": locals().get("hard_meta", {}).get("n_unassigned_cells"),
            "prior_effect_fraction_delta": locals().get("hard_meta", {}).get("prior_effect_fraction_delta"),
            "prior_effect_fraction_delta_mean_per_spot": locals().get("hard_meta", {}).get(
                "prior_effect_fraction_delta_mean_per_spot"
            ),
            "prior_effect_fraction_delta_audit_before": locals().get("hard_meta", {}).get(
                "prior_effect_fraction_delta_audit_before"
            ),
            "prior_effect_fraction_delta_audit_after": locals().get("hard_meta", {}).get(
                "prior_effect_fraction_delta_audit_after"
            ),
            "prior_local_top1_change_rate": locals().get("hard_meta", {}).get("prior_local_top1_change_rate"),
            "prior_ablation_status": locals().get("hard_meta", {}).get("prior_ablation_status"),
            "prior_ablation_change_rate": locals().get("hard_meta", {}).get("prior_ablation_change_rate"),
            "prior_ablation_error_reason": locals().get("hard_meta", {}).get("prior_ablation_error_reason"),
            "prior_ablation_audit": locals().get("hard_meta", {}).get("prior_ablation_audit"),
            "type_prior_row_nonzero_ratio": locals().get("hard_meta", {}).get("type_prior_row_nonzero_ratio"),
            "type_prior_row_entropy_mean": locals().get("hard_meta", {}).get("type_prior_row_entropy_mean"),
            "type_prior_row_entropy_max": locals().get("hard_meta", {}).get("type_prior_row_entropy_max"),
            "gene_usage_stats": {
                "plugin_genes_total": len(plugin_genes),
                "genes_use": len(genes_use) if "genes_use" in locals() else None,
                "genes_missing_in_sc": int(len(set(plugin_genes) - set(sc_expr.columns))),
                "genes_missing_in_st": int(len(set(plugin_genes) - set(st_expr.columns))),
                "weight_min": (locals().get("w_stats") or {}).get("min"),
                "weight_p50": (locals().get("w_stats") or {}).get("p50"),
                "weight_p95": (locals().get("w_stats") or {}).get("p95"),
                "weight_max": (locals().get("w_stats") or {}).get("max"),
                "weighted_gene_count": len(weighted_genes_applied) if "weighted_genes_applied" in locals() else None,
                "weighted_gene_ratio": locals().get("weighted_ratio_applied"),
            },
            "runtime_sec": time.time() - t0,
            "status": status,
            "error_msg": error_msg,
            "resolved_mapping_config": config,
        }
        _write_meta(out_dir, meta)
