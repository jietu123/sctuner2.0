"""
SimGen: simulate ST scenarios with ground-truth for Stage5 evaluation.

Key contracts (⚠️):
- sim_root is unified (default: data/sim/<scenario_id>/...)
- Stage1-exported sc MUST be Query (not Ref)
- Truth naming is stable with alias compatibility:
  - world truth: sim_truth_cell_spot.csv (World cell -> spot)
  - query truth: sim_truth_query_cell_spot.csv (Query cell -> spot) + alias cell_true_spot.csv
  - spot truth:  sim_truth_spot_type_fraction.csv + alias spot_true_type_fraction.csv
- spot_id set/order must be identical between st_expression_normalized.csv and st_coordinates.csv
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import yaml
from scipy.spatial import cKDTree


def detect_project_root() -> Path:
    here = Path(__file__).resolve()
    try:
        return here.parents[2]
    except IndexError:
        return here.parent.parent


def sha1_bytes(data: bytes) -> str:
    return hashlib.sha1(data).hexdigest().lower()


def sha1_file(path: Path) -> Optional[str]:
    try:
        return sha1_bytes(path.read_bytes())
    except Exception:
        return None


def load_yaml(path: Path) -> dict:
    if not path.exists() or path.stat().st_size == 0:
        return {}
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def _now_iso() -> str:
    return datetime.now().isoformat(timespec="seconds")


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _write_csv_with_index(df: pd.DataFrame, path: Path, index_name: str) -> None:
    df2 = df.copy()
    df2.index.name = index_name
    df2.to_csv(path)


def _write_csv_no_index(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, index=False)


def _read_lines(path: Path) -> List[str]:
    return [ln.strip() for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]


def _select_type_col(sc_meta: pd.DataFrame) -> str:
    if "celltype" in sc_meta.columns:
        return "celltype"
    if "cell_type" in sc_meta.columns:
        return "cell_type"
    if "type" in sc_meta.columns:
        return "type"
    raise KeyError("sc_metadata.csv 需要包含 celltype/cell_type/type 之一作为类型列")


def _normalize_seurat_log1p(counts: np.ndarray, scale_factor: float = 1e4) -> np.ndarray:
    counts = np.asarray(counts, dtype=float)
    lib = counts.sum(axis=1, keepdims=True)
    lib_safe = np.where(lib > 0, lib, 1.0)
    x = counts / lib_safe * float(scale_factor)
    return np.log1p(x)


def _pseudo_counts_from_norm(
    expr_norm: np.ndarray,
    rng: np.random.Generator,
    *,
    target_depth: int,
) -> np.ndarray:
    x = np.asarray(expr_norm, dtype=float)
    w = np.expm1(x)
    w[w < 0] = 0.0
    row_sum = w.sum(axis=1, keepdims=True)
    out = np.zeros_like(w, dtype=np.int32)
    for i in range(w.shape[0]):
        s = float(row_sum[i, 0])
        if s <= 0:
            continue
        p = w[i] / s
        out[i] = rng.multinomial(int(target_depth), p)
    return out


@dataclass(frozen=True)
class SpotGrid:
    spot_ids: List[str]
    centers_xy: np.ndarray  # shape (S,2) in continuous space
    grid_rc: np.ndarray  # shape (S,2) integer row/col (grid indices)


def _make_spot_grid(space_w: float, space_h: float, nx: int, ny: int) -> SpotGrid:
    xs = np.linspace(0, space_w, nx, endpoint=False) + (space_w / nx) / 2.0
    ys = np.linspace(0, space_h, ny, endpoint=False) + (space_h / ny) / 2.0
    centers = []
    grid_rc = []
    spot_ids = []
    k = 0
    for i, x in enumerate(xs):
        for j, y in enumerate(ys):
            spot_ids.append(f"spot_{k:04d}")
            centers.append((x, y))
            grid_rc.append((i, j))
            k += 1
    return SpotGrid(
        spot_ids=spot_ids,
        centers_xy=np.asarray(centers, dtype=float),
        grid_rc=np.asarray(grid_rc, dtype=int),
    )


def _block_bounds(space_w: float, space_h: float, block_nx: int, block_ny: int) -> Dict[str, Tuple[float, float, float, float]]:
    bw = space_w / block_nx
    bh = space_h / block_ny
    out: Dict[str, Tuple[float, float, float, float]] = {}
    for i in range(block_nx):
        for j in range(block_ny):
            x0 = i * bw
            x1 = (i + 1) * bw
            y0 = j * bh
            y1 = (j + 1) * bh
            out[f"block_{i}_{j}"] = (x0, x1, y0, y1)
    return out


def _neighbors_block(i: int, j: int, block_nx: int, block_ny: int) -> List[Tuple[int, int]]:
    cand = []
    for di, dj in ((0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)):
        ni, nj = i + di, j + dj
        if 0 <= ni < block_nx and 0 <= nj < block_ny:
            cand.append((ni, nj))
    return cand


def _build_type_block_mixture(
    types: Sequence[str],
    block_nx: int,
    block_ny: int,
    rng: np.random.Generator,
    *,
    primary_weight: float = 0.75,
    neighbor_weight: float = 0.23,
    background_weight: float = 0.02,
) -> Dict[str, Dict[str, float]]:
    if not (0 < primary_weight < 1):
        raise ValueError("primary_weight must be in (0,1)")
    if primary_weight + neighbor_weight + background_weight <= 0:
        raise ValueError("invalid weights")
    block_ids = [f"block_{i}_{j}" for i in range(block_nx) for j in range(block_ny)]
    out: Dict[str, Dict[str, float]] = {}
    for t in types:
        pi = int(rng.integers(0, block_nx))
        pj = int(rng.integers(0, block_ny))
        neigh = _neighbors_block(pi, pj, block_nx, block_ny)
        neigh_ids = [f"block_{i}_{j}" for i, j in neigh]
        p: Dict[str, float] = {bid: 0.0 for bid in block_ids}
        for bid in neigh_ids:
            p[bid] += neighbor_weight / max(1, len(neigh_ids))
        p[f"block_{pi}_{pj}"] += primary_weight
        for bid in block_ids:
            p[bid] += background_weight / len(block_ids)
        s = sum(p.values())
        for k in list(p.keys()):
            p[k] = float(p[k] / s)
        out[str(t)] = p
    return out


def _sample_xy_for_cells(
    cell_types: Sequence[str],
    block_mix: Dict[str, Dict[str, float]],
    block_bounds_map: Dict[str, Tuple[float, float, float, float]],
    rng: np.random.Generator,
) -> Tuple[np.ndarray, List[str]]:
    block_ids = list(block_bounds_map.keys())
    coords = np.zeros((len(cell_types), 2), dtype=float)
    chosen_blocks: List[str] = []
    for i, t in enumerate(cell_types):
        mix = block_mix.get(str(t))
        if not mix:
            probs = np.ones(len(block_ids), dtype=float) / len(block_ids)
        else:
            probs = np.array([mix[bid] for bid in block_ids], dtype=float)
            probs = probs / probs.sum()
        b = str(rng.choice(block_ids, p=probs))
        x0, x1, y0, y1 = block_bounds_map[b]
        coords[i, 0] = float(rng.uniform(x0, x1))
        coords[i, 1] = float(rng.uniform(y0, y1))
        chosen_blocks.append(b)
    return coords, chosen_blocks


def _assign_to_spots_nearest_with_radius(
    cell_xy: np.ndarray,
    grid: SpotGrid,
    *,
    radius: float,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    tree = cKDTree(grid.centers_xy)
    dist, idx = tree.query(cell_xy, k=1)
    idx = np.asarray(idx, dtype=int)
    dist = np.asarray(dist, dtype=float)
    assigned = np.where(dist <= float(radius), idx, -1)
    audit = {
        "radius": float(radius),
        "n_cells": int(cell_xy.shape[0]),
        "n_uncovered_cells": int((assigned < 0).sum()),
        "uncovered_fraction": float((assigned < 0).mean()) if cell_xy.shape[0] > 0 else 0.0,
        "dist_min": float(dist.min()) if dist.size else None,
        "dist_median": float(np.median(dist)) if dist.size else None,
        "dist_p90": float(np.quantile(dist, 0.9)) if dist.size else None,
        "dist_max": float(dist.max()) if dist.size else None,
    }
    return assigned, audit


def _ensure_gene_subset(df: pd.DataFrame, genes: List[str], context: str) -> pd.DataFrame:
    missing = [g for g in genes if g not in df.columns]
    if missing:
        raise KeyError(f"[{context}] 缺少基因列: {missing[:5]} (total_missing={len(missing)})")
    return df.loc[:, genes]


def _load_real_scrna_inputs(
    stage1_export_dir: Path,
    genes: List[str],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    sc_meta = pd.read_csv(stage1_export_dir / "sc_metadata.csv")
    if "cell_id" not in sc_meta.columns:
        raise KeyError("Stage1 sc_metadata.csv 缺少 cell_id 列")
    type_col_src = _select_type_col(sc_meta)
    sc_meta = sc_meta.copy()
    sc_meta["celltype"] = sc_meta[type_col_src].astype(str)
    sc_meta = sc_meta[["cell_id", "celltype"]].copy()

    sc_expr = pd.read_csv(stage1_export_dir / "sc_expression_normalized.csv", index_col=0, usecols=["cell_id"] + genes)
    sc_expr.index = sc_expr.index.astype(str)
    sc_expr = _ensure_gene_subset(sc_expr, genes, context="sc_expression_normalized")
    return sc_expr, sc_meta


def _apply_type_cleaning(
    sc_meta: pd.DataFrame,
    *,
    min_cells_per_type: int,
    low_count_policy: str,
    other_name: str,
    merge_map_path: Optional[Path],
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    df = sc_meta.copy()
    if merge_map_path and merge_map_path.exists() and merge_map_path.stat().st_size > 0:
        mp = pd.read_csv(merge_map_path)
        if not {"orig_type", "new_type"}.issubset(mp.columns):
            raise ValueError("merge_map_path CSV 需要包含列 orig_type,new_type")
        m = {str(a): str(b) for a, b in zip(mp["orig_type"], mp["new_type"])}
        df["celltype"] = df["celltype"].astype(str).map(lambda x: m.get(str(x), str(x)))

    vc = df["celltype"].value_counts()
    small = set(vc[vc < int(min_cells_per_type)].index.astype(str).tolist())
    policy = str(low_count_policy).strip().lower()
    if policy not in ("drop", "merge_to_other"):
        raise ValueError("type_cleaning.low_count_policy 必须是 drop 或 merge_to_other")
    if policy == "drop":
        df = df[~df["celltype"].astype(str).isin(small)].copy()
    else:
        df.loc[df["celltype"].astype(str).isin(small), "celltype"] = str(other_name)

    audit = {
        "min_cells_per_type": int(min_cells_per_type),
        "low_count_policy": policy,
        "other_name": str(other_name),
        "merge_map_path": str(merge_map_path) if merge_map_path else None,
        "n_types_before": int(vc.shape[0]),
        "n_types_after": int(df["celltype"].nunique()),
        "small_types": sorted(list(small))[:50],
        "n_small_types": int(len(small)),
    }
    return df.reset_index(drop=True), audit


def _build_type_plan(
    types: Sequence[str],
    counts: Dict[str, int],
    cfg: Dict[str, Any],
) -> pd.DataFrame:
    sc_missing = set([str(x) for x in (cfg.get("sc_missing_types") or [])])
    st_missing = set([str(x) for x in (cfg.get("st_missing_types") or [])])
    wf = float(cfg.get("world_fraction", 0.6))
    qf = float(cfg.get("query_fraction", 0.4))
    rf = float(cfg.get("ref_fraction", 0.0))

    rows = []
    for t in types:
        n = int(counts.get(t, 0))
        world_n = int(np.floor(n * wf))
        query_n = int(np.floor(n * qf))
        ref_n = int(np.floor(n * rf))

        world_action = "keep"
        query_action = "keep"
        ref_action = "keep"

        if t in st_missing:
            world_action = "remove"
            world_n = 0
        if t in sc_missing:
            query_action = "remove"
            query_n = 0

        rows.append(
            {
                "type": t,
                "n_available": n,
                "world_action": world_action,
                "world_target_n": world_n,
                "query_action": query_action,
                "query_target_n": query_n,
                "ref_action": ref_action,
                "ref_target_n": ref_n,
            }
        )
    return pd.DataFrame(rows)


def _sample_sets_by_type(
    sc_meta: pd.DataFrame,
    type_plan: pd.DataFrame,
    rng: np.random.Generator,
    *,
    allow_world_ref_overlap: bool,
    allow_world_query_overlap: bool,
    allow_ref_query_overlap: bool,
) -> Tuple[List[str], List[str], List[str], Dict[str, Any], pd.DataFrame]:
    by_type: Dict[str, List[str]] = {}
    for t, sub in sc_meta.groupby("celltype"):
        by_type[str(t)] = sub["cell_id"].astype(str).tolist()

    world_ids: List[str] = []
    query_ids: List[str] = []
    ref_ids: List[str] = []

    per_type_rows = []
    for _, r in type_plan.iterrows():
        t = str(r["type"])
        pool_all = list(by_type.get(t, []))
        rng.shuffle(pool_all)
        n_avail = int(r["n_available"])
        if n_avail != len(pool_all):
            n_avail = len(pool_all)

        w_n = int(r["world_target_n"])
        q_n = int(r["query_target_n"])
        ref_n = int(r["ref_target_n"])

        if w_n > n_avail:
            raise ValueError(f"type={t}: world_target_n={w_n} > n_available={n_avail}")
        world_t = pool_all[:w_n]

        pool_for_query = pool_all if allow_world_query_overlap else [x for x in pool_all if x not in set(world_t)]
        if q_n > len(pool_for_query):
            raise ValueError(f"type={t}: query_target_n={q_n} > available_for_query={len(pool_for_query)}")
        query_t = pool_for_query[:q_n]

        pool_for_ref = pool_all
        if not allow_world_ref_overlap:
            pool_for_ref = [x for x in pool_for_ref if x not in set(world_t)]
        if not allow_ref_query_overlap:
            pool_for_ref = [x for x in pool_for_ref if x not in set(query_t)]
        if ref_n > len(pool_for_ref):
            raise ValueError(f"type={t}: ref_target_n={ref_n} > available_for_ref={len(pool_for_ref)}")
        ref_t = pool_for_ref[:ref_n]

        world_ids.extend(world_t)
        query_ids.extend(query_t)
        ref_ids.extend(ref_t)

        per_type_rows.append(
            {
                "type": t,
                "n_available": int(n_avail),
                "n_world": int(len(world_t)),
                "n_query": int(len(query_t)),
                "n_ref": int(len(ref_t)),
            }
        )

    world_set = set(world_ids)
    query_set = set(query_ids)
    ref_set = set(ref_ids)

    audit = {
        "allow_world_ref_overlap": bool(allow_world_ref_overlap),
        "allow_world_query_overlap": bool(allow_world_query_overlap),
        "allow_ref_query_overlap": bool(allow_ref_query_overlap),
        "n_world": int(len(world_ids)),
        "n_query": int(len(query_ids)),
        "n_ref": int(len(ref_ids)),
        "overlap_world_query": int(len(world_set & query_set)),
        "overlap_world_ref": int(len(world_set & ref_set)),
        "overlap_ref_query": int(len(ref_set & query_set)),
    }
    per_type_counts = pd.DataFrame(per_type_rows)
    return world_ids, query_ids, ref_ids, audit, per_type_counts


def _allocate_capacity_from_world_density(
    world_spot_counts: np.ndarray,
    rng: np.random.Generator,
    *,
    n_query_cells: int,
    clip_min: int = 1,
) -> np.ndarray:
    n_spots = int(world_spot_counts.shape[0])
    if n_query_cells < n_spots * int(clip_min):
        raise ValueError(f"n_query_cells={n_query_cells} < n_spots*clip_min={n_spots*clip_min}; 请减少 spot 数或增大 query 规模")
    cap = np.full(n_spots, int(clip_min), dtype=int)
    remaining = int(n_query_cells - int(cap.sum()))
    if remaining <= 0:
        return cap
    w = np.asarray(world_spot_counts, dtype=float)
    w = np.where(w > 0, w, 0.0)
    if w.sum() <= 0:
        w = np.ones_like(w, dtype=float)
    p = w / w.sum()
    extra = rng.choice(np.arange(n_spots), size=remaining, replace=True, p=p)
    for idx in extra:
        cap[int(idx)] += 1
    return cap


def validate_scenario(
    *,
    sim_out_dir: Path,
    stage1_export_dir: Path,
    require_query_truth: bool,
    require_stage1_gene_lists: bool = False,
    scenario_meta_path: Optional[Path] = None,
) -> None:
    required_sim = [
        "sim_sc_expression.csv",
        "sim_sc_metadata.csv",
        "sim_st_expression.csv",
        "sim_st_coordinates.csv",
        "sim_truth_cell_spot.csv",
        "sim_truth_spot_type_fraction.csv",
        "scenario_meta.json",
    ]
    if require_query_truth:
        required_sim.append("sim_truth_query_cell_spot.csv")
    missing = [f for f in required_sim if not (sim_out_dir / f).exists()]
    if missing:
        raise FileNotFoundError(f"[validate] sim_out_dir 缺少文件: {missing}")

    required_stage1 = [
        "sc_expression_normalized.csv",
        "st_expression_normalized.csv",
        "sc_metadata.csv",
        "st_coordinates.csv",
    ]
    missing2 = [f for f in required_stage1 if not (stage1_export_dir / f).exists()]
    if missing2:
        raise FileNotFoundError(f"[validate] stage1_export_dir 缺少文件: {missing2}")

    st_expr = pd.read_csv(stage1_export_dir / "st_expression_normalized.csv", index_col=0)
    st_coords = pd.read_csv(stage1_export_dir / "st_coordinates.csv", index_col=0)
    if not st_expr.index.is_unique:
        raise ValueError("[validate] st_expression_normalized.csv spot_id 存在重复")
    if not st_coords.index.is_unique:
        raise ValueError("[validate] st_coordinates.csv spot_id 存在重复")
    if set(st_expr.index.astype(str)) != set(st_coords.index.astype(str)):
        raise ValueError("[validate] st_expression 与 st_coordinates 的 spot_id 集合不一致")
    if not st_expr.index.astype(str).equals(st_coords.index.astype(str)):
        raise ValueError("[validate] st_expression 与 st_coordinates 的 spot_id 顺序不一致（禁止隐性错配）")

    sc_meta = pd.read_csv(stage1_export_dir / "sc_metadata.csv")
    if "cell_id" not in sc_meta.columns or "celltype" not in sc_meta.columns:
        raise ValueError("[validate] sc_metadata.csv 必须包含 cell_id 与 celltype 列")

    if require_query_truth:
        dt = pd.read_csv(sim_out_dir / "sim_truth_query_cell_spot.csv")
        if "cell_id" not in dt.columns or "true_spot_id" not in dt.columns:
            raise ValueError("[validate] sim_truth_query_cell_spot.csv 必须包含 cell_id,true_spot_id")

    if require_stage1_gene_lists:
        stage1_dir = stage1_export_dir.parent
        for name in ("common_genes.txt", "hvg_genes.txt"):
            if not (stage1_dir / name).exists():
                raise FileNotFoundError(f"[validate] stage1_preprocess 缺少基因列表文件: {name}")
        try:
            common = set(_read_lines(stage1_dir / "common_genes.txt"))
            hvg = _read_lines(stage1_dir / "hvg_genes.txt")
            if not set(hvg).issubset(common):
                raise ValueError("[validate] hvg_genes.txt 必须是 common_genes.txt 的子集")
            if len(common) > 0 and len(hvg) == len(common):
                print(f"[validate][WARNING] HVG equals common (n={len(common)}). Consider lowering stage1_compat_gene_lists.hvg_cap_n_top.")
        except Exception as e:
            raise ValueError(f"[validate] 读取/校验基因列表失败: {e}")

    # Optional feasibility warning + write back to scenario_meta.json (no extra files)
    if scenario_meta_path and scenario_meta_path.exists():
        try:
            meta = json.loads(scenario_meta_path.read_text(encoding="utf-8"))
        except Exception:
            meta = {}
        try:
            n_spots = int(pd.read_csv(stage1_export_dir / "st_expression_normalized.csv", index_col=0).shape[0])
        except Exception:
            n_spots = None
        try:
            n_world = int((meta.get("overlap_audit") or {}).get("n_world")) if (meta.get("overlap_audit") or {}).get("n_world") is not None else None
        except Exception:
            n_world = None
        try:
            n_query = int((meta.get("overlap_audit") or {}).get("n_query")) if (meta.get("overlap_audit") or {}).get("n_query") is not None else None
        except Exception:
            n_query = None
        try:
            min_cells = (meta.get("coverage_audit") or {}).get("min_cells_per_spot")
            min_cells = int(min_cells) if min_cells is not None else None
        except Exception:
            min_cells = None

        # Coverage explicit check/warn (no spot drop; spot_id must stay unchanged)
        cov = meta.get("coverage_audit") or {}
        try:
            max_empty = cov.get("max_empty_spot_fraction")
            max_empty = float(max_empty) if max_empty is not None else None
        except Exception:
            max_empty = None
        try:
            empty_frac = cov.get("empty_spot_fraction")
            empty_frac = float(empty_frac) if empty_frac is not None else None
        except Exception:
            empty_frac = None
        try:
            satisfied_frac = cov.get("min_cells_per_spot_satisfied_fraction")
            satisfied_frac = float(satisfied_frac) if satisfied_frac is not None else None
        except Exception:
            satisfied_frac = None
        try:
            warn_thr = cov.get("min_cells_per_spot_satisfied_fraction_warn")
            warn_thr = float(warn_thr) if warn_thr is not None else 0.50
        except Exception:
            warn_thr = 0.50

        if empty_frac is None or satisfied_frac is None:
            try:
                steps = list(cov.get("fallback_steps") or [])
                for s in reversed(steps):
                    if str(s.get("action")) == "coverage_check":
                        if empty_frac is None:
                            empty_frac = float(s.get("empty_spot_fraction"))
                        if satisfied_frac is None:
                            satisfied_frac = float(s.get("min_cells_per_spot_satisfied_fraction"))
                        break
            except Exception:
                pass

        # Fallback: derive from exported st_coordinates (world_cells_in_spot)
        try:
            if "world_cells_in_spot" in st_coords.columns:
                world_counts = pd.to_numeric(st_coords["world_cells_in_spot"], errors="coerce").fillna(0.0)
                if empty_frac is None:
                    empty_frac = float((world_counts <= 0).mean())
                if satisfied_frac is None and min_cells is not None:
                    satisfied_frac = float((world_counts >= float(min_cells)).mean())
        except Exception:
            pass

        if max_empty is not None and empty_frac is not None and float(empty_frac) > float(max_empty):
            raise ValueError(f"[validate] coverage: empty_spot_fraction={empty_frac} > max_empty_spot_fraction={max_empty}")
        if satisfied_frac is not None and float(satisfied_frac) < float(warn_thr):
            print(
                f"[validate][WARNING] coverage: min_cells_per_spot_satisfied_fraction={satisfied_frac} < warn_threshold={warn_thr}; "
                f"consider lowering min_cells_per_spot or reducing spot grid."
            )

        avg_world = (float(n_world) / float(n_spots)) if (n_world is not None and n_spots) else None
        avg_query = (float(n_query) / float(n_spots)) if (n_query is not None and n_spots) else None

        # Heuristic: when min is much higher than average, satisfied_fraction can be misleadingly low.
        feasibility_warn = False
        rule = "min_cells_per_spot > avg_world_cells_per_spot + 0.5"
        if min_cells is not None and avg_world is not None:
            feasibility_warn = bool(float(min_cells) > float(avg_world) + 0.5)
            if feasibility_warn:
                print(
                    f"[validate][WARNING] min_cells_per_spot={min_cells} > avg_world_cells_per_spot={avg_world:.3f} + 0.5 "
                    f"(n_world={n_world}, n_spots={n_spots}); expected low satisfied_fraction. Consider lowering min or reducing spot grid."
                )

        meta["feasibility_audit"] = {
            "rule": rule,
            "feasibility_warn": bool(feasibility_warn),
            "min_cells_per_spot": min_cells,
            "n_spots": n_spots,
            "n_world_cells": n_world,
            "n_query_cells": n_query,
            "avg_world_cells_per_spot": avg_world,
            "avg_query_cells_per_spot": avg_query,
        }
        scenario_meta_path.write_text(json.dumps(meta, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_stage1_compat_gene_lists(
    *,
    project_root: Path,
    scenario_id: str,
    cfg: Dict[str, Any],
    input_sample: str,
    stage1_dir: Path,
    query_expr: pd.DataFrame,
    st_expr: pd.DataFrame,
) -> Dict[str, Any]:
    """
    Fill minimal gene list files under data/processed/<scenario_id>/stage1_preprocess/ (non-exported):
      - common_genes.txt: intersection of sc/st exported genes (always real for this scenario)
      - hvg_genes.txt: prefer copy+filter from source_sample; fallback to variance top-N (reproducible)
      - genes_manifest.json: small provenance (optional, default on)
    """
    compat = dict(cfg.get("stage1_compat_gene_lists") or {})
    enabled = bool(compat.get("enabled", True))
    if not enabled:
        return {"enabled": False}

    source_sample = str(compat.get("source_sample") or input_sample or "real_brca")
    filter_to_sim_common = bool(compat.get("filter_to_sim_common", True))
    generate_if_missing = bool(compat.get("generate_if_missing", True))
    copy_files = compat.get("copy_files") or ["common_genes.txt", "hvg_genes.txt"]
    write_manifest = bool(compat.get("write_manifest", True))
    hvg_cap_n_top = compat.get("hvg_cap_n_top", None)
    if hvg_cap_n_top in (None, "", 0):
        hvg_cap_n_top_int: Optional[int] = None
    else:
        hvg_cap_n_top_int = int(hvg_cap_n_top)
    min_hvg_after_filter = compat.get("min_hvg_after_filter", None)
    hvg_fallback = dict(compat.get("hvg_fallback") or {})
    hvg_method = str(hvg_fallback.get("method", "variance")).strip().lower()
    hvg_n_top = int(hvg_fallback.get("n_top", hvg_cap_n_top_int or 2000))
    if min_hvg_after_filter in (None, "", 0):
        min_after_filter = int(hvg_fallback.get("min_after_filter", 200))
    else:
        min_after_filter = int(min_hvg_after_filter)

    stage1_dir.mkdir(parents=True, exist_ok=True)

    sc_genes = [str(c) for c in list(query_expr.columns)]
    st_genes = [str(c) for c in list(st_expr.columns)]
    sim_common = sorted(set(sc_genes) & set(st_genes))
    common_path = stage1_dir / "common_genes.txt"
    common_text = "\n".join(sim_common) + ("\n" if sim_common else "")
    with common_path.open("w", encoding="utf-8", newline="\n") as f:
        f.write(common_text)

    common_set = set(sim_common)
    hvg_genes: List[str] = []
    hvg_source: str = "missing"
    hvg_cap_applied = False
    src_stage1_dir = project_root / "data" / "processed" / source_sample / "stage1_preprocess"
    src_hvg = src_stage1_dir / "hvg_genes.txt"
    if "hvg_genes.txt" in set(map(str, copy_files)) and src_hvg.exists():
        raw = _read_lines(src_hvg)
        seen: set = set()
        dedup = []
        for g in raw:
            if g in seen:
                continue
            seen.add(g)
            dedup.append(g)
        if filter_to_sim_common:
            dedup = [g for g in dedup if g in common_set]
        if hvg_cap_n_top_int is not None and len(dedup) > hvg_cap_n_top_int:
            dedup = dedup[:hvg_cap_n_top_int]
            hvg_cap_applied = True
        if len(dedup) >= max(1, min_after_filter) or not generate_if_missing:
            hvg_genes = dedup
            hvg_source = "copied_filtered"

    if (not hvg_genes) and generate_if_missing:
        if hvg_method != "variance":
            raise ValueError(f"stage1_compat_gene_lists.hvg_fallback.method 仅支持 variance，当前: {hvg_method}")
        if not sim_common:
            raise ValueError("sim_common genes 为空，无法生成 HVG 列表")
        df_for_var = query_expr[sim_common].copy()
        var = df_for_var.var(axis=0)
        ranked = [str(x) for x in list(var.sort_values(ascending=False).index)]
        if filter_to_sim_common:
            ranked = [g for g in ranked if g in common_set]
        hvg_genes = ranked[: max(0, min(int(hvg_n_top), len(ranked)))]
        if hvg_cap_n_top_int is not None and len(hvg_genes) > hvg_cap_n_top_int:
            hvg_genes = hvg_genes[:hvg_cap_n_top_int]
            hvg_cap_applied = True
        hvg_source = "generated_variance"

    hvg_path = stage1_dir / "hvg_genes.txt"
    hvg_text = "\n".join(hvg_genes) + ("\n" if hvg_genes else "")
    with hvg_path.open("w", encoding="utf-8", newline="\n") as f:
        f.write(hvg_text)

    audit = {
        "enabled": True,
        "source_sample": source_sample,
        "filter_to_sim_common": filter_to_sim_common,
        "generate_if_missing": generate_if_missing,
        "hvg_cap_n_top": hvg_cap_n_top_int,
        "hvg_cap_applied": bool(hvg_cap_applied),
        "min_hvg_after_filter": int(min_after_filter),
        "common_genes_count": int(len(sim_common)),
        "hvg_genes_count": int(len(hvg_genes)),
        "hvg_source": hvg_source,
        "hvg_fallback_method": hvg_method,
        "hvg_fallback_n_top": int(hvg_n_top),
    }
    if write_manifest:
        sha1_common_file = sha1_file(common_path) if common_path.exists() else None
        sha1_hvg_file = sha1_file(hvg_path) if hvg_path.exists() else None
        manifest = {
            "scenario_id": scenario_id,
            "generated_at": _now_iso(),
            "source_sample": source_sample,
            "common_genes_count": int(len(sim_common)),
            "hvg_genes_count": int(len(hvg_genes)),
            "copy_used": bool(hvg_source.startswith("copied")),
            "filter_to_sim_common": filter_to_sim_common,
            "generate_if_missing": generate_if_missing,
            "hvg_cap_n_top": hvg_cap_n_top_int,
            "hvg_cap_applied": bool(hvg_cap_applied),
            "min_hvg_after_filter": int(min_after_filter),
            "hvg_source": hvg_source,
            "hvg_fallback": {"method": hvg_method, "n_top": int(hvg_n_top)},
            "sha1_common_genes": sha1_common_file,
            "sha1_hvg_genes": sha1_hvg_file,
        }
        (stage1_dir / "genes_manifest.json").write_text(
            json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        audit["manifest_written"] = True
    else:
        audit["manifest_written"] = False
    return audit


def run_simgen_for_scenario(
    *,
    project_root: Path,
    scenario_cfg: Dict[str, Any],
    scenario_id: str,
    sim_root: Path,
    overwrite: bool,
) -> Path:
    cfg = dict(scenario_cfg or {})
    cfg.setdefault("scenario_id", scenario_id)
    cfg.setdefault("seed", 42)
    cfg.setdefault("sim_root", str(sim_root.as_posix()))

    out_dir = project_root / sim_root / scenario_id
    stage1_export_dir = project_root / "data" / "processed" / scenario_id / "stage1_preprocess" / "exported"
    stage1_dir = stage1_export_dir.parent

    if out_dir.exists() and not overwrite:
        raise FileExistsError(f"SimGen 输出目录已存在: {out_dir}；如需覆盖请加 --overwrite")
    if stage1_dir.exists() and not overwrite:
        raise FileExistsError(f"Stage1-exported 目录已存在: {stage1_dir}；如需覆盖请加 --overwrite")
    if overwrite:
        # Keep output directories clean: avoid stale/legacy files accumulating.
        for p in [out_dir, stage1_dir]:
            if p.exists():
                shutil.rmtree(p)

    _ensure_dir(out_dir)
    _ensure_dir(stage1_export_dir)

    rng = np.random.default_rng(int(cfg.get("seed", 42)))

    # ---- Input: real scRNA (normalized) as template
    input_sample = str(cfg.get("input_sample") or "real_brca")
    real_stage1_dir = project_root / "data" / "processed" / input_sample / "stage1_preprocess"
    real_stage1_export = real_stage1_dir / "exported"
    if not real_stage1_export.exists():
        raise FileNotFoundError(f"找不到真实 Stage1 exported: {real_stage1_export}")

    # ---- Gene panel
    gene_panel = str(cfg.get("gene_panel") or "hvg").strip().lower()
    gene_panel_path = cfg.get("gene_panel_path")
    if gene_panel_path:
        gene_path = (project_root / gene_panel_path) if not Path(gene_panel_path).is_absolute() else Path(gene_panel_path)
        genes = _read_lines(gene_path)
    else:
        if gene_panel == "hvg":
            genes = _read_lines(real_stage1_dir / "hvg_genes.txt")
        elif gene_panel == "common":
            genes = _read_lines(real_stage1_dir / "common_genes.txt")
        else:
            raise ValueError("gene_panel 仅支持 hvg/common 或提供 gene_panel_path")
    n_genes_max = cfg.get("n_genes_max", None)
    if n_genes_max not in (None, "", 0):
        genes = genes[: int(n_genes_max)]

    # ---- Load sc inputs (normalized)
    sc_expr_all, sc_meta_all = _load_real_scrna_inputs(real_stage1_export, genes=genes)

    # ---- Type cleaning
    tc = cfg.get("type_cleaning") or {}
    merge_map_path = tc.get("merge_map_path")
    merge_map_path = (project_root / merge_map_path) if merge_map_path else None
    sc_meta_all, type_clean_audit = _apply_type_cleaning(
        sc_meta_all,
        min_cells_per_type=int(tc.get("min_cells_per_type", 30)),
        low_count_policy=str(tc.get("low_count_policy", "merge_to_other")),
        other_name=str(tc.get("other_name", "Other")),
        merge_map_path=merge_map_path,
    )

    keep_cells = set(sc_meta_all["cell_id"].astype(str).tolist())
    sc_expr_all = sc_expr_all.loc[[cid for cid in sc_expr_all.index.astype(str) if cid in keep_cells], :]
    sc_expr_all.index = sc_expr_all.index.astype(str)

    type_counts = sc_meta_all["celltype"].value_counts().to_dict()
    type_vocab = sorted([str(x) for x in sc_meta_all["celltype"].unique().tolist()])
    type_plan = _build_type_plan(type_vocab, {str(k): int(v) for k, v in type_counts.items()}, cfg)

    # ---- Overlap constraints (P5)
    allow_world_ref_overlap = bool(cfg.get("allow_world_ref_overlap", False))
    allow_world_query_overlap = bool(cfg.get("allow_world_query_overlap", False))
    allow_ref_query_overlap = bool(cfg.get("allow_ref_query_overlap", False))
    world_ids, query_ids, ref_ids, overlap_audit, per_type_counts = _sample_sets_by_type(
        sc_meta_all,
        type_plan,
        rng,
        allow_world_ref_overlap=allow_world_ref_overlap,
        allow_world_query_overlap=allow_world_query_overlap,
        allow_ref_query_overlap=allow_ref_query_overlap,
    )

    role_map = {"stage1_sc": "query", "world": "world", "ref": "ref(optional)"}

    sc_meta_by_id = sc_meta_all.set_index("cell_id", drop=False)
    query_ids = [str(x) for x in query_ids]
    world_ids = [str(x) for x in world_ids]
    ref_ids = [str(x) for x in ref_ids]

    query_expr = sc_expr_all.loc[query_ids, :]
    world_expr = sc_expr_all.loc[world_ids, :]
    ref_expr = sc_expr_all.loc[ref_ids, :] if ref_ids else None

    query_meta = pd.DataFrame(
        {
            "cell_id": query_ids,
            "type": [str(sc_meta_by_id.loc[cid, "celltype"]) for cid in query_ids],
        }
    )
    query_meta["celltype"] = query_meta["type"]
    query_meta["cell_type"] = query_meta["type"]

    ref_meta = None
    if ref_ids:
        ref_meta = pd.DataFrame(
            {
                "cell_id": ref_ids,
                "type": [str(sc_meta_by_id.loc[cid, "celltype"]) for cid in ref_ids],
            }
        )
        ref_meta["celltype"] = ref_meta["type"]
        ref_meta["cell_type"] = ref_meta["type"]

    world_types = [str(sc_meta_by_id.loc[cid, "celltype"]) for cid in world_ids]

    # ---- Spatial layout
    space_w = float(cfg.get("space_width", 1000))
    space_h = float(cfg.get("space_height", 1000))
    block_nx = int(cfg.get("block_nx", 5))
    block_ny = int(cfg.get("block_ny", 5))
    spot_grid_nx = int(cfg.get("spot_grid_nx", 25))
    spot_grid_ny = int(cfg.get("spot_grid_ny", 25))
    spot_radius = float(cfg.get("spot_radius", 80))
    min_cells_per_spot = int(cfg.get("min_cells_per_spot", 3))
    min_cells_per_spot_satisfied_fraction_warn = float(cfg.get("min_cells_per_spot_satisfied_fraction_warn", 0.50))

    grid = _make_spot_grid(space_w, space_h, spot_grid_nx, spot_grid_ny)
    n_spots = len(grid.spot_ids)
    if len(world_ids) < n_spots * min_cells_per_spot:
        raise ValueError(
            f"world_cells={len(world_ids)} < n_spots*min_cells_per_spot={n_spots*min_cells_per_spot}; "
            f"请减少 spot_grid 或降低 min_cells_per_spot 或增大 world_fraction/允许 overlap"
        )

    bounds = _block_bounds(space_w, space_h, block_nx, block_ny)
    mix_cfg = cfg.get("block_mixture") or {}
    block_mix = _build_type_block_mixture(
        type_vocab,
        block_nx,
        block_ny,
        rng,
        primary_weight=float(mix_cfg.get("primary_weight", 0.50)),
        neighbor_weight=float(mix_cfg.get("neighbor_weight", 0.30)),
        background_weight=float(mix_cfg.get("background_weight", 0.20)),
    )

    # ---- Coverage fallback (P4: spot_id 不变，禁止 drop)
    max_resample_rounds = int(cfg.get("max_resample_rounds", 5))
    max_radius_multiplier = float(cfg.get("max_radius_multiplier", 2.0))
    max_empty_spot_fraction = float(cfg.get("max_empty_spot_fraction", 0.02))
    fallback_steps = list(cfg.get("coverage_fallback") or ["expand_radius", "resample_world"])
    if any(str(x).strip().lower() == "drop_empty_spots" for x in fallback_steps):
        raise ValueError("coverage_fallback 不允许使用 drop_empty_spots（spot_id 必须保持不变）")

    allow_expand_radius = any(str(x).strip().lower() == "expand_radius" for x in fallback_steps)
    allow_resample_world = any(str(x).strip().lower() == "resample_world" for x in fallback_steps)

    radius = float(spot_radius)
    coverage_audit_steps: List[Dict[str, Any]] = []
    world_xy = None
    world_blocks = None
    world_assigned = None
    last_cov: Optional[Dict[str, float]] = None
    resample_reason = "initial"

    for attempt in range(max_resample_rounds + 1):
        radius_start = float(radius)
        resample_step: Dict[str, Any] = {
            "action": "resample_world",
            "attempt": int(attempt),
            "reason": str(resample_reason),
            "radius_start": radius_start,
        }
        if last_cov is not None:
            resample_step.update(
                {
                    "prev_empty_spot_fraction": float(last_cov.get("empty_spot_fraction", float("nan"))),
                    "prev_min_cells_per_spot_satisfied_fraction": float(
                        last_cov.get("min_cells_per_spot_satisfied_fraction", float("nan"))
                    ),
                }
            )
        coverage_audit_steps.append(resample_step)

        world_xy, world_blocks = _sample_xy_for_cells(world_types, block_mix, bounds, rng)
        world_assigned, assign_audit = _assign_to_spots_nearest_with_radius(world_xy, grid, radius=radius_start)

        # uncovered -> expand_radius (if allowed), else resample_world / fail
        if (world_assigned < 0).any():
            resample_step.update(
                {
                    "n_uncovered_cells_before": int(assign_audit.get("n_uncovered_cells", 0)),
                    "uncovered_fraction_before": float(assign_audit.get("uncovered_fraction", 0.0)),
                }
            )
            if not allow_expand_radius:
                resample_step.update(
                    {
                        "radius_end": radius_start,
                        "n_uncovered_cells_after": int(assign_audit.get("n_uncovered_cells", 0)),
                        "uncovered_fraction_after": float(assign_audit.get("uncovered_fraction", 0.0)),
                        "status_after": "fail_uncovered",
                    }
                )
                if (attempt >= max_resample_rounds) or (not allow_resample_world):
                    raise RuntimeError(
                        f"覆盖率兜底失败：存在 uncovered_cells 且未允许 expand_radius/resample_world（attempt={attempt}, radius={radius_start}）"
                    )
                resample_reason = "uncovered_cells"
                continue

            for expand_iter in range(8):
                if not (world_assigned < 0).any():
                    break
                if float(radius) > float(spot_radius) * float(max_radius_multiplier):
                    break
                radius_before = float(radius)
                audit_before = dict(assign_audit)
                radius = float(radius) * 1.15
                world_assigned, assign_audit = _assign_to_spots_nearest_with_radius(world_xy, grid, radius=radius)
                expand_step: Dict[str, Any] = {
                    "action": "expand_radius",
                    "attempt": int(attempt),
                    "iter": int(expand_iter),
                    "radius_before": radius_before,
                    "radius_after": float(radius),
                    "n_uncovered_cells_before": int(audit_before.get("n_uncovered_cells", 0)),
                    "n_uncovered_cells_after": int(assign_audit.get("n_uncovered_cells", 0)),
                    "uncovered_fraction_before": float(audit_before.get("uncovered_fraction", 0.0)),
                    "uncovered_fraction_after": float(assign_audit.get("uncovered_fraction", 0.0)),
                }
                if not (world_assigned < 0).any():
                    tmp_counts = np.bincount(world_assigned, minlength=n_spots)
                    expand_step.update(
                        {
                            "empty_spot_fraction_after": float((tmp_counts == 0).mean()),
                            "min_cells_per_spot_satisfied_fraction_after": float((tmp_counts >= min_cells_per_spot).mean()),
                        }
                    )
                coverage_audit_steps.append(expand_step)

            if (world_assigned < 0).any():
                resample_step.update(
                    {
                        "radius_end": float(radius),
                        "n_uncovered_cells_after": int(assign_audit.get("n_uncovered_cells", 0)),
                        "uncovered_fraction_after": float(assign_audit.get("uncovered_fraction", 0.0)),
                        "status_after": "fail_uncovered",
                    }
                )
                if (attempt >= max_resample_rounds) or (not allow_resample_world):
                    raise RuntimeError(
                        f"覆盖率兜底失败：uncovered_cells 无法通过 expand_radius 解决（attempt={attempt}, radius={radius}）"
                    )
                resample_reason = "uncovered_cells"
                continue

        world_spot_counts = np.bincount(world_assigned, minlength=n_spots)
        empty_frac = float((world_spot_counts == 0).mean())
        satisfied_frac = float((world_spot_counts >= min_cells_per_spot).mean())

        resample_step.update(
            {
                "radius_end": float(radius),
                "empty_spot_fraction_after": empty_frac,
                "min_cells_per_spot_satisfied_fraction_after": satisfied_frac,
                "max_empty_spot_fraction": float(max_empty_spot_fraction),
            }
        )
        if last_cov is not None:
            resample_step.update(
                {
                    "delta_empty_spot_fraction": float(empty_frac - float(last_cov.get("empty_spot_fraction", empty_frac))),
                    "delta_min_cells_per_spot_satisfied_fraction": float(
                        satisfied_frac - float(last_cov.get("min_cells_per_spot_satisfied_fraction", satisfied_frac))
                    ),
                }
            )

        check_step: Dict[str, Any] = {
            "action": "coverage_check",
            "attempt": int(attempt),
            "radius_start": radius_start,
            "radius_end": float(radius),
            "min_cells_per_spot": int(min_cells_per_spot),
            "empty_spot_fraction": empty_frac,
            "min_cells_per_spot_satisfied_fraction": satisfied_frac,
            "max_empty_spot_fraction": float(max_empty_spot_fraction),
        }
        if last_cov is not None:
            check_step.update(
                {
                    "prev_empty_spot_fraction": float(last_cov.get("empty_spot_fraction", empty_frac)),
                    "prev_min_cells_per_spot_satisfied_fraction": float(
                        last_cov.get("min_cells_per_spot_satisfied_fraction", satisfied_frac)
                    ),
                    "delta_empty_spot_fraction": float(empty_frac - float(last_cov.get("empty_spot_fraction", empty_frac))),
                    "delta_min_cells_per_spot_satisfied_fraction": float(
                        satisfied_frac - float(last_cov.get("min_cells_per_spot_satisfied_fraction", satisfied_frac))
                    ),
                }
            )
        coverage_audit_steps.append(check_step)
        last_cov = {"empty_spot_fraction": empty_frac, "min_cells_per_spot_satisfied_fraction": satisfied_frac}

        if empty_frac <= max_empty_spot_fraction:
            resample_step["status_after"] = "pass"
            check_step["status"] = "pass"
            break

        resample_step["status_after"] = "fail_empty"
        check_step["status"] = "fail_empty"
        if attempt >= max_resample_rounds:
            raise RuntimeError(f"覆盖率兜底失败：empty_spot_fraction={empty_frac} > {max_empty_spot_fraction}")
        if not allow_resample_world:
            raise RuntimeError(
                f"覆盖率兜底失败：empty_spot_fraction={empty_frac} > {max_empty_spot_fraction} 且 coverage_fallback 未允许 resample_world"
            )
        resample_reason = "empty_spot_fraction>max"

    assert world_xy is not None and world_blocks is not None and world_assigned is not None
    world_spot_counts = np.bincount(world_assigned, minlength=n_spots)
    empty_spot_fraction_final = float((world_spot_counts == 0).mean())
    min_cells_per_spot_satisfied_fraction_final = float((world_spot_counts >= min_cells_per_spot).mean())

    # ---- Counts-space simulation (R2/P3)
    counts_cfg = cfg.get("counts") or {}
    counts_source = str(counts_cfg.get("counts_source", "pseudo_from_norm"))
    target_depth = int(counts_cfg.get("target_depth", 5000))
    scale_factor = float(counts_cfg.get("scale_factor", 1e4))
    if counts_source != "pseudo_from_norm":
        raise ValueError("当前实现仅支持 counts.counts_source=pseudo_from_norm（Stage1 export 为 log-normalized）")
    world_counts_cell = _pseudo_counts_from_norm(world_expr.to_numpy(dtype=float), rng, target_depth=target_depth)

    spot_counts = np.zeros((n_spots, len(genes)), dtype=np.int32)
    for i_cell, s_idx in enumerate(world_assigned):
        spot_counts[int(s_idx), :] += world_counts_cell[i_cell, :]

    # ---- Spot noise in counts space
    noise_cfg = cfg.get("noise") or {}
    libsize_mode = str(noise_cfg.get("libsize_mode", "lognormal"))
    libsize_sigma = float(noise_cfg.get("libsize_sigma", 0.3))
    noise_level = float(noise_cfg.get("noise_level", 1.0))
    if libsize_mode == "lognormal":
        lib_factor = np.exp(rng.normal(loc=0.0, scale=libsize_sigma, size=(n_spots, 1)))
    elif libsize_mode in ("none", "uniform"):
        lib_factor = np.ones((n_spots, 1), dtype=float)
    else:
        raise ValueError("noise.libsize_mode 仅支持 lognormal/none/uniform")
    lam = np.asarray(spot_counts, dtype=float) * lib_factor
    lam[lam < 0] = 0.0
    if noise_level <= 0:
        spot_counts_noised = np.rint(lam).astype(np.int32)
        noise_model = "deterministic_round"
    else:
        spot_counts_noised = rng.poisson(lam=lam).astype(np.int32)
        noise_model = "poisson"

    st_norm = _normalize_seurat_log1p(spot_counts_noised, scale_factor=scale_factor)
    st_expr = pd.DataFrame(st_norm, index=pd.Index(grid.spot_ids, name="spot_id"), columns=genes)

    # ---- World truth (World cell -> spot)
    world_df = pd.DataFrame(
        {
            "cell_id": world_ids,
            "true_type": world_types,
            "block_id": world_blocks,
            "sim_x": world_xy[:, 0],
            "sim_y": world_xy[:, 1],
            "true_spot_id": [grid.spot_ids[int(s)] for s in world_assigned],
        }
    )
    truth_world = world_df[["cell_id", "true_spot_id", "true_type"]].copy()
    # spot×type truth (rows: spot_id, cols: type_vocab) — include missing-type columns as all-zeros (M2)
    spot_type_counts = pd.crosstab(world_df["true_spot_id"], world_df["true_type"])
    spot_type_counts = spot_type_counts.reindex(index=grid.spot_ids, fill_value=0).reindex(columns=type_vocab, fill_value=0)
    spot_type_fraction = spot_type_counts.div(spot_type_counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    spot_type_fraction.index = pd.Index(grid.spot_ids, name="spot_id")

    # ---- Query truth (P2): sample latent xy then nearest spot; st_missing_types -> NA
    query_truth_policy = str(cfg.get("query_truth_policy", "nearest_spot_after_sampling_xy"))
    st_missing_types = set([str(x) for x in (cfg.get("st_missing_types") or [])])
    query_types = query_meta["celltype"].astype(str).tolist()
    query_xy, query_blocks = _sample_xy_for_cells(query_types, block_mix, bounds, rng)
    q_assigned, q_assign_audit = _assign_to_spots_nearest_with_radius(query_xy, grid, radius=radius)
    query_true_spot = []
    truth_status = []
    for i, t in enumerate(query_types):
        if t in st_missing_types:
            query_true_spot.append("")
            truth_status.append("st_missing_type")
        else:
            si = int(q_assigned[i])
            if si < 0:
                query_true_spot.append("")
                truth_status.append("uncovered")
            else:
                query_true_spot.append(grid.spot_ids[si])
                truth_status.append("ok")
    truth_query = pd.DataFrame(
        {
            "cell_id": query_ids,
            "true_spot_id": query_true_spot,
            "celltype": query_types,
            "truth_status": truth_status,
            "block_id": query_blocks,
            "sim_x": query_xy[:, 0],
            "sim_y": query_xy[:, 1],
        }
    )

    # ---- Capacity: ensure sum == n_query (helps Stage4 outputs align with Query truth)
    cap_policy = str(cfg.get("capacity_policy", "world_density"))
    if cap_policy not in ("world_density", "uniform"):
        raise ValueError("capacity_policy 仅支持 world_density/uniform")
    weight_for_cap = world_spot_counts.astype(float) if cap_policy == "world_density" else np.ones(n_spots, dtype=float)
    capacity = _allocate_capacity_from_world_density(
        weight_for_cap,
        rng,
        n_query_cells=len(query_ids),
        clip_min=int(cfg.get("capacity_clip_min", 1)),
    )

    # ---- st_coordinates (Stage1-compatible)
    nCount_RNA = spot_counts_noised.sum(axis=1).astype(int)
    nFeature_RNA = (spot_counts_noised > 0).sum(axis=1).astype(int)
    st_coords = pd.DataFrame(
        {
            "row": grid.grid_rc[:, 0],
            "col": grid.grid_rc[:, 1],
            "x": grid.centers_xy[:, 0],
            "y": grid.centers_xy[:, 1],
            "world_cells_in_spot": world_spot_counts.astype(int),
            "in_tissue": (world_spot_counts > 0).astype(int),
            "spot_cell_counts": capacity.astype(int),
            "nCount_RNA": nCount_RNA,
            "nFeature_RNA": nFeature_RNA,
            "UMI_total": nCount_RNA,
        },
        index=pd.Index(grid.spot_ids, name="spot_id"),
    )

    # ---- Output policy: optionally drop non-tissue spots from all ST outputs
    output_in_tissue_only = bool(cfg.get("output_in_tissue_only", True))
    spot_output_audit = {
        "output_in_tissue_only": output_in_tissue_only,
        "n_spots_total": int(len(st_coords)),
        "n_spots_output": int(len(st_coords)),
        "n_spots_dropped": 0,
    }
    if output_in_tissue_only:
        tissue_mask = st_coords["in_tissue"].astype(int).to_numpy() == 1
        n_drop = int((~tissue_mask).sum())
        spot_output_audit["n_spots_output"] = int(tissue_mask.sum())
        spot_output_audit["n_spots_dropped"] = n_drop
        if n_drop > 0:
            keep_spot_ids = st_coords.index[tissue_mask]
            drop_spot_ids = set(st_coords.index[~tissue_mask].astype(str))
            st_expr = st_expr.loc[keep_spot_ids]
            st_coords = st_coords.loc[keep_spot_ids].copy()
            spot_type_fraction = spot_type_fraction.loc[keep_spot_ids]

            # Recompute capacity for kept spots so sum == n_query (Stage4 expects full pool)
            capacity = _allocate_capacity_from_world_density(
                world_spot_counts[tissue_mask],
                rng,
                n_query_cells=len(query_ids),
                clip_min=int(cfg.get("capacity_clip_min", 1)),
            )
            st_coords["spot_cell_counts"] = capacity.astype(int)

            # Query truth: if assigned to a dropped spot, reassign to nearest tissue spot
            if "true_spot_id" in truth_query.columns:
                truth_policy = str(cfg.get("non_tissue_truth_policy", "reassign")).strip().lower()
                if truth_policy not in ("reassign", "missing"):
                    raise ValueError("non_tissue_truth_policy must be reassign or missing")
                mask_drop = truth_query["true_spot_id"].astype(str).isin(drop_spot_ids)
                if mask_drop.any():
                    if truth_policy == "missing":
                        truth_query.loc[mask_drop, "true_spot_id"] = ""
                        if "truth_status" in truth_query.columns:
                            truth_query.loc[mask_drop, "truth_status"] = "non_tissue_spot_filtered"
                    else:
                        spot_xy = st_coords[["x", "y"]].to_numpy(dtype=float)
                        tree = cKDTree(spot_xy)
                        q_xy = truth_query.loc[mask_drop, ["sim_x", "sim_y"]].to_numpy(dtype=float)
                        _, nn_idx = tree.query(q_xy, k=1)
                        new_spot_ids = st_coords.index.to_numpy()[nn_idx]
                        truth_query.loc[mask_drop, "true_spot_id"] = new_spot_ids
                        if "truth_status" in truth_query.columns:
                            truth_query.loc[mask_drop, "truth_status"] = "non_tissue_spot_reassigned"
                    spot_output_audit["non_tissue_truth_policy"] = truth_policy
                    spot_output_audit["n_truth_reassigned"] = int(mask_drop.sum()) if truth_policy == "reassign" else 0
                    spot_output_audit["n_truth_missing_due_to_non_tissue"] = int(mask_drop.sum()) if truth_policy == "missing" else 0

    # ---- Write outputs
    _write_csv_with_index(query_expr, out_dir / "sim_sc_expression.csv", index_name="cell_id")
    _write_csv_no_index(query_meta, out_dir / "sim_sc_metadata.csv")
    _write_csv_with_index(st_expr, out_dir / "sim_st_expression.csv", index_name="spot_id")
    _write_csv_with_index(st_coords, out_dir / "sim_st_coordinates.csv", index_name="spot_id")
    _write_csv_no_index(truth_world, out_dir / "sim_truth_cell_spot.csv")
    _write_csv_with_index(spot_type_fraction, out_dir / "sim_truth_spot_type_fraction.csv", index_name="spot_id")
    _write_csv_no_index(truth_query, out_dir / "sim_truth_query_cell_spot.csv")
    if ref_expr is not None and ref_meta is not None:
        _write_csv_with_index(ref_expr, out_dir / "sim_ref_sc_expression.csv", index_name="cell_id")
        _write_csv_no_index(ref_meta, out_dir / "sim_ref_sc_metadata.csv")
    # Optional debug outputs (off by default; keep SimGen outputs minimal per SVTuner.md)
    if bool(cfg.get("write_debug_outputs", False)):
        _write_csv_no_index(type_plan, out_dir / "type_plan.csv")
        _write_csv_no_index(per_type_counts, out_dir / "per_type_counts_world_query_ref.csv")
        _write_csv_no_index(truth_query, out_dir / "query_debug.csv")

    # Aliases for compatibility
    shutil.copy2(out_dir / "sim_truth_spot_type_fraction.csv", out_dir / "spot_true_type_fraction.csv")
    shutil.copy2(out_dir / "sim_truth_query_cell_spot.csv", out_dir / "cell_true_spot.csv")

    # ---- Stage1-exported: sc must be Query (P0)
    _write_csv_with_index(query_expr, stage1_export_dir / "sc_expression_normalized.csv", index_name="cell_id")
    _write_csv_with_index(st_expr, stage1_export_dir / "st_expression_normalized.csv", index_name="spot_id")
    sc_meta_export = query_meta[["cell_id", "celltype"]].copy()
    sc_meta_export["cell_type"] = sc_meta_export["celltype"]
    _write_csv_no_index(sc_meta_export, stage1_export_dir / "sc_metadata.csv")
    _write_csv_with_index(st_coords, stage1_export_dir / "st_coordinates.csv", index_name="spot_id")
    # Stage1 compat gene lists (non-exported): keep small, avoid future missing-file rework.
    stage1_gene_audit = write_stage1_compat_gene_lists(
        project_root=project_root,
        scenario_id=scenario_id,
        cfg=cfg,
        input_sample=input_sample,
        stage1_dir=stage1_dir,
        query_expr=query_expr,
        st_expr=st_expr,
    )

    # ---- scenario_meta.json
    scenario_meta = {
        "scenario_id": scenario_id,
        "description": cfg.get("description"),
        "generated_at": _now_iso(),
        "sim_root": str(sim_root.as_posix()),
        "output_dir": str(out_dir.as_posix()),
        "role_map": role_map,
        "cell_level_mode": "query",
        "query_truth_policy": query_truth_policy,
        "truth_files": {
            "spot_type_fraction": {"canonical": "sim_truth_spot_type_fraction.csv", "aliases": ["spot_true_type_fraction.csv"]},
            "cell_spot_query": {"canonical": "sim_truth_query_cell_spot.csv", "aliases": ["cell_true_spot.csv"]},
            "cell_spot_world": {"canonical": "sim_truth_cell_spot.csv", "aliases": []},
        },
        "overlap_audit": overlap_audit,
        "type_cleaning_audit": type_clean_audit,
        "type_vocab": type_vocab,
        "coverage_audit": {
            "min_cells_per_spot": min_cells_per_spot,
            "max_empty_spot_fraction": max_empty_spot_fraction,
            "max_resample_rounds": max_resample_rounds,
            "max_radius_multiplier": max_radius_multiplier,
            "coverage_fallback": fallback_steps,
            "radius_final": float(radius),
            "empty_spot_fraction": empty_spot_fraction_final,
            "min_cells_per_spot_satisfied_fraction": min_cells_per_spot_satisfied_fraction_final,
            "min_cells_per_spot_satisfied_fraction_warn": float(min_cells_per_spot_satisfied_fraction_warn),
            "fallback_steps": coverage_audit_steps,
        },
        "spot_output_audit": spot_output_audit,
        "expression_audit": {
            "counts_source": counts_source,
            "target_depth": target_depth,
            "scale_factor": scale_factor,
            "noise_model": noise_model,
            "libsize_mode": libsize_mode,
            "libsize_sigma": libsize_sigma,
            "noise_level": noise_level,
        },
        "capacity_audit": {
            "capacity_policy": cap_policy,
            "capacity_clip_min": int(cfg.get("capacity_clip_min", 1)),
            "capacity_sum": int(capacity.sum()),
            "n_query_cells": int(len(query_ids)),
        },
        "input_files": [
            {"path": str((real_stage1_export / "sc_expression_normalized.csv").as_posix()), "sha1": sha1_file(real_stage1_export / "sc_expression_normalized.csv")},
            {"path": str((real_stage1_export / "sc_metadata.csv").as_posix()), "sha1": sha1_file(real_stage1_export / "sc_metadata.csv")},
        ],
        "code": {"file": str(Path(__file__).resolve().as_posix()), "sha1": sha1_file(Path(__file__).resolve())},
        "runtime": {"python": sys.version, "numpy": np.__version__, "pandas": pd.__version__},
        "stage1_compat_gene_lists_audit": stage1_gene_audit,
        "config_effective_subset": cfg,
    }
    (out_dir / "scenario_meta.json").write_text(json.dumps(scenario_meta, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    validate_scenario(
        sim_out_dir=out_dir,
        stage1_export_dir=stage1_export_dir,
        require_query_truth=True,
        require_stage1_gene_lists=bool((cfg.get("stage1_compat_gene_lists") or {}).get("enabled", True)),
        scenario_meta_path=(out_dir / "scenario_meta.json"),
    )
    return out_dir


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="SimGen: generate simulation scenario under data/sim/<scenario_id>/")
    p.add_argument("--scenario_id", required=True, help="Scenario ID, e.g. S0_matched")
    p.add_argument("--config", required=True, help="Scenario YAML path, e.g. configs/simgen/S0_matched.yaml")
    p.add_argument("--sim_root", default="data/sim", help='Sim root dir (default: "data/sim")')
    p.add_argument("--overwrite", action="store_true", default=False, help="Overwrite existing outputs")
    p.add_argument("--migrate_old_root", action="store_true", default=False, help="Migrate from data/simulations/<sid>/ if present")
    return p.parse_args()


def maybe_migrate_old_root(project_root: Path, scenario_id: str, sim_root: Path) -> None:
    old_root = project_root / "data" / "simulations" / scenario_id
    new_root = project_root / sim_root / scenario_id
    if old_root.exists() and not new_root.exists():
        print(f"[SimGen] Found old dir {old_root} -> copying to {new_root}")
        shutil.copytree(old_root, new_root)


def main() -> None:
    args = parse_args()
    project_root = detect_project_root()
    cfg_path = Path(args.config)
    if not cfg_path.is_absolute():
        cfg_path = project_root / cfg_path
    scenario_cfg = load_yaml(cfg_path)

    sim_root = Path(str(args.sim_root))
    if args.migrate_old_root:
        maybe_migrate_old_root(project_root, args.scenario_id, sim_root)

    out_dir = run_simgen_for_scenario(
        project_root=project_root,
        scenario_cfg=scenario_cfg,
        scenario_id=args.scenario_id,
        sim_root=sim_root,
        overwrite=bool(args.overwrite),
    )
    print(f"[SimGen] Done: {out_dir}")


if __name__ == "__main__":
    main()
