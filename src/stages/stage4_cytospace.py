"""
Stage4 (Route2, full): wrap external CytoSPACE mapping with Stage3 filtering.

Pipeline:
- Load Stage1 exports (sc_expression_normalized.csv, st_expression_normalized.csv, st_coordinates.csv)
- Load Stage3 outputs (cell_type_relabel.csv, type_support.csv)
- Filter cells: either drop plugin_type == Unknown_sc_only (default) or drop support_category == unsupported
- Prepare CytoSPACE input files (gene x cell, gene x spot, cell type table, type fractions)
- Run CytoSPACE main_cytospace
- Convert outputs to Stage4 format (cell_assignment.csv, stage4_summary.json)

Outputs are written to result/<sample>/stage4_cytospace/
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

# Redirect numba JIT cache to a writable temp directory.
# This avoids PermissionError when scanpy/__pycache__ is not writable
# (common in conda envs on Windows without admin rights).
_numba_cache = Path(tempfile.gettempdir()) / "numba_cache_svtuner"
_numba_cache.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("NUMBA_CACHE_DIR", str(_numba_cache))
from typing import Literal, Sequence

import hashlib
import numpy as np
import pandas as pd
import sys as _sys
import yaml

# 确保项目根目录在 sys.path 中，便于在多种运行方式下都能导入 src.*
_HERE = Path(__file__).resolve()
_ROOT = _HERE.parents[2]
if str(_ROOT) not in _sys.path:
    _sys.path.insert(0, str(_ROOT))

from src.stages.stage1_io import load_stage1
from src.utils.type_name import load_alias_map, canonicalize_type_name, normalize_type_name
from src.utils.sample_paths import sample_dir_candidates


FilterMode = Literal["plugin_unknown", "unsupported", "none"]
FilterScope = Literal["unsupported_all", "missing_only", "missing_detected_only"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Route2 Stage4 full (CytoSPACE wrapper)")
    p.add_argument("--sample", default="real_brca", help="sample id")
    p.add_argument("--project_root", default=None, help="override project root")
    p.add_argument("--stage4_suffix", default="", help="suffix appended to stage4_cytospace dir (e.g. '_baseline', '_route2')")
    p.add_argument("--filter_mode", choices=["plugin_unknown", "unsupported", "none"], default="none",
                   help="Filter mode: 'none' for baseline (official CytoSPACE), 'plugin_unknown' or 'unsupported' for Route2")
    p.add_argument("--cell_type_column", choices=["plugin_type", "orig_type", "sc_meta"], default="sc_meta",
                   help="which column to use as cell_type: 'sc_meta' for baseline (original types), 'plugin_type' for Route2")
    p.add_argument("--missing_type", default="T cells CD8", help="target missing type for reporting")
    p.add_argument("--mean_cell_numbers", type=int, default=5, help="CytoSPACE mean_cell_numbers")
    p.add_argument("--solver_method", default="lap_CSPR", help="CytoSPACE solver_method (lap_CSPR avoids lapjv dependency)")
    p.add_argument("--distance_metric", default="Pearson_correlation", help="CytoSPACE distance_metric")
    p.add_argument("--seed", type=int, default=42, help="random seed")
    p.add_argument("--sampling_sub_spots", action="store_true", default=True, help="enable sub-spot sampling to reduce memory")
    p.add_argument("--n_subspots", type=int, default=800, help="number of selected sub-spots when sampling_sub_spots is on")
    p.add_argument(
        "--n_processors",
        type=int,
        default=1,
        help="number of processors (reduce to 1 to avoid memory issues)",
    )
    p.add_argument(
        "--mapping_cells_per_spot",
        type=int,
        default=None,
        help="override cells_per_spot for CytoSPACE mapping (optional)",
    )
    p.add_argument(
        "--filter_scope",
        choices=["unsupported_all", "missing_only", "missing_detected_only"],
        default="unsupported_all",
        help=(
            "scope of Route2 filtering when filter_mode uses Stage3 plugin/unsupported info; "
            "'unsupported_all' (default) drops all unsupported/Unknown_sc_only cells; "
            "'missing_only' drops all cells whose orig_type matches the specified missing_type "
            "(oracle-style fallback, regardless of Stage3 detection); "
            "'missing_detected_only' only drops cells that are both Stage3-marked and match the specified missing_type"
        ),
    )
    p.add_argument(
        "--filter_to_sim_truth",
        action="store_true",
        help="filter cell pool to sim_truth_query_cell_spot.csv (simulation-only evaluation alignment)",
    )
    p.add_argument(
        "--hotspot_rescue_post",
        action="store_true",
        help="enable post-mapping hotspot rescue (stage4-level)",
    )
    p.add_argument(
        "--hotspot_rescue_qvalue_min",
        default=None,
        help="minimum qvalue for rescue window (use 'final_fdr' to follow Stage3 alpha)",
    )
    p.add_argument(
        "--hotspot_rescue_qvalue_max",
        type=float,
        default=None,
        help="maximum qvalue for rescue window (default 0.1)",
    )
    p.add_argument(
        "--hotspot_rescue_min_genes",
        type=int,
        default=None,
        help="minimum marker genes required to flag hotspot (default 3)",
    )
    p.add_argument(
        "--hotspot_rescue_expr_quantile",
        type=float,
        default=None,
        help="expression quantile for hotspot mask (default 0.9)",
    )
    p.add_argument(
        "--hotspot_rescue_min_cluster_fraction",
        type=float,
        default=None,
        help="minimum hotspot spot fraction to accept rescue (default 0.01)",
    )
    p.add_argument(
        "--hotspot_rescue_dedupe",
        action="store_true",
        help="dedupe cell_assignment by cell_id before filling empty slots",
    )
    p.add_argument(
        "--hotspot_rescue_cells_per_spot",
        type=int,
        default=None,
        help="cells_per_spot used when filling rescue empty slots (optional)",
    )
    p.add_argument(
        "--hotspot_rescue_allow_nonsignificant",
        action="store_true",
        help="allow Significant=No types to enter rescue candidate window (still filtered by qvalue)",
    )
    p.add_argument(
        "--hotspot_rescue_max_replace_fraction",
        type=float,
        default=None,
        help="max fraction of total assignments to replace during rescue (0-1, optional)",
    )
    p.add_argument(
        "--hotspot_rescue_replace_low_confidence",
        action="store_true",
        help="replace lowest marker-support cells within hotspot instead of random picks",
    )
    p.add_argument(
        "--hotspot_rescue_replace_support_max",
        type=float,
        default=None,
        help="only replace spots with marker support <= threshold (optional)",
    )
    p.add_argument(
        "--keep_redundant",
        action="store_true",
        help="keep assigned_expression and sc/st_expression_for_cytospace (default: delete to save ~6GB per run)",
    )
    p.add_argument(
        "--sc_expr_source",
        type=str,
        choices=["normalized", "data", "counts", "auto"],
        default="normalized",
        help="SC expression source from stage1 export used for CytoSPACE input.",
    )
    return p.parse_args()


def detect_root() -> Path:
    here = Path(__file__).resolve()
    try:
        return here.parents[2]
    except IndexError:
        return here.parent.parent


def _clean_spot_index(idx: pd.Index) -> pd.Index:
    return pd.Index([str(x).split("\t")[0] for x in idx])


def load_stage3(root: Path, sample: str) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    base = root / "data" / "processed" / sample / "stage3_typematch"
    relabel = pd.read_csv(base / "cell_type_relabel.csv")
    ts_path = base / "type_support.csv"
    ts = pd.read_csv(ts_path) if ts_path.exists() else None
    return relabel, ts


def load_stage3_thresholds(root: Path, sample: str) -> tuple[float | None, float | None]:
    cfg_path = root / "configs" / "datasets" / f"{sample}.yaml"
    if not cfg_path.exists():
        return None, None
    try:
        import yaml  # type: ignore
    except Exception:
        return None, None
    with cfg_path.open("r", encoding="utf-8") as f:
        data = yaml.safe_load(f)
    stage3 = (data or {}).get("stage3") or {}
    return stage3.get("weak_th"), stage3.get("strong_th")


def load_stage4_truth_filter(root: Path, sample: str) -> bool | None:
    cfg_path = root / "configs" / "datasets" / f"{sample}.yaml"
    if not cfg_path.exists():
        return None
    try:
        import yaml  # type: ignore
    except Exception:
        return None
    with cfg_path.open("r", encoding="utf-8") as f:
        data = yaml.safe_load(f)
    stage4 = (data or {}).get("stage4") or {}
    val = stage4.get("filter_to_sim_truth")
    if val is None:
        return None
    return bool(val)


def load_dataset_cfg(root: Path, sample: str) -> dict:
    cfg_path = root / "configs" / "datasets" / f"{sample}.yaml"
    if not cfg_path.exists():
        return {}
    with cfg_path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def parse_overrides(value, cast_type):
    overrides = {}
    if isinstance(value, str):
        for part in value.split(","):
            part = part.strip()
            if not part:
                continue
            if "=" in part:
                key, val = part.split("=", 1)
            elif ":" in part:
                key, val = part.split(":", 1)
            else:
                continue
            key = key.strip().lower()
            val = val.strip()
            try:
                overrides[key] = cast_type(val)
            except (TypeError, ValueError):
                continue
    elif isinstance(value, dict):
        for key, val in value.items():
            key = str(key).strip().lower()
            try:
                overrides[key] = cast_type(val)
            except (TypeError, ValueError):
                continue
    return overrides


def resolve_final_fdr(stage3_cfg: dict, stage3_summary_path: Path) -> float:
    alpha = stage3_cfg.get("alpha")
    if stage3_summary_path.exists():
        try:
            summary = json.loads(stage3_summary_path.read_text(encoding="utf-8"))
            params = summary.get("params", {})
            alpha = params.get("alpha", alpha)
        except Exception:
            pass
    try:
        return float(alpha) if alpha is not None else 0.05
    except (TypeError, ValueError):
        return 0.05


def resolve_qvalue_min(value, final_fdr: float) -> float:
    if value is None:
        return final_fdr
    if isinstance(value, str) and value.strip().lower() == "final_fdr":
        return final_fdr
    try:
        return float(value)
    except (TypeError, ValueError):
        return final_fdr


def build_hotspot_mask(
    st_expr: pd.DataFrame,
    marker_genes: list[str],
    expr_quantile: float,
    min_genes: int,
) -> tuple[pd.Series, pd.Series]:
    genes = [g for g in marker_genes if g in st_expr.columns]
    if len(genes) < min_genes:
        empty = pd.Series(False, index=st_expr.index)
        counts = pd.Series(0, index=st_expr.index)
        return empty, counts
    expr = st_expr[genes]
    high_flags = []
    for g in genes:
        vals = expr[g].to_numpy(dtype=float)
        thr = float(np.quantile(vals, expr_quantile))
        high_flags.append(vals > thr)
    high_flags = np.stack(high_flags, axis=1)
    counts = pd.Series(high_flags.sum(axis=1), index=st_expr.index)
    mask = counts >= min_genes
    return mask, counts


def prepare_post_rescue_assignment(
    assignment_path: Path,
    relabel: pd.DataFrame,
    sc_expr: pd.DataFrame,
    sc_meta: pd.DataFrame,
    st_expr: pd.DataFrame,
    stage1_export: Path,
    stage3_cfg: dict,
    q_min: float,
    q_max: float,
    expr_quantile: float,
    min_genes: int,
    min_cluster_fraction: float,
    cells_per_spot: int,
    seed: int,
    dedupe: bool,
    allow_nonsignificant: bool,
    max_replace_fraction: float | None,
    replace_low_confidence: bool,
    replace_support_max: float | None,
) -> tuple[pd.DataFrame, dict]:
    from src.stages.stage3_type_plugin import select_marker_genes_with_specificity, select_marker_genes

    assignment = pd.read_csv(assignment_path)
    base_assignment = assignment
    slot_assignment = assignment
    stats = {
        "enabled": True,
        "deduped": False,
        "dedup_removed": 0,
        "allow_nonsignificant": bool(allow_nonsignificant),
        "candidate_types": [],
        "rescued_types": [],
        "rescued_cells_total": 0,
        "rescued_cells_by_type": {},
        "replaced_cells_total": 0,
        "replaced_cells_by_type": {},
        "hotspot_fraction_by_type": {},
        "max_replace_fraction": float(max_replace_fraction) if max_replace_fraction is not None else 0.0,
        "replace_low_confidence": bool(replace_low_confidence),
        "replace_support_max": float(replace_support_max) if replace_support_max is not None else None,
        "replace_support_filtered": 0,
    }

    base_assignment = base_assignment.copy()
    base_assignment["assigned_spot"] = base_assignment["assigned_spot"].astype(str).str.split(r"\s|\t").str[0]

    if dedupe:
        before = len(slot_assignment)
        slot_assignment = slot_assignment.drop_duplicates(subset=["cell_id"], keep="first")
        stats["deduped"] = True
        stats["dedup_removed"] = int(before - len(slot_assignment))

    type_support_path = stage1_export.parent.parent / "stage3_typematch" / "type_support.csv"
    if not type_support_path.exists():
        return base_assignment, stats
    type_support = pd.read_csv(type_support_path)
    if "QValue" not in type_support.columns or "Significant" not in type_support.columns:
        return base_assignment, stats

    if allow_nonsignificant:
        sig_mask = pd.Series(True, index=type_support.index)
    else:
        sig_mask = type_support["Significant"].astype(str).str.strip().str.lower() == "yes"
    qvals = pd.to_numeric(type_support["QValue"], errors="coerce")
    rescue_mask = sig_mask & qvals.ge(q_min) & qvals.le(q_max)
    missing_types = type_support.loc[rescue_mask, "orig_type"].astype(str).tolist()
    if not missing_types:
        return base_assignment, stats

    stats["candidate_types"] = missing_types

    plugin_genes_path = stage3_cfg.get("plugin_genes_path") or (stage1_export.parent / "hvg_genes.txt")
    if not Path(plugin_genes_path).exists():
        return base_assignment, stats
    plugin_genes = [g.strip() for g in Path(plugin_genes_path).read_text(encoding="utf-8").splitlines() if g.strip()]
    plugin_genes = [g for g in plugin_genes if g in sc_expr.columns]
    if not plugin_genes:
        return base_assignment, stats

    if "cell_type" in sc_meta.columns:
        type_col = "cell_type"
    elif "celltype" in sc_meta.columns:
        type_col = "celltype"
    else:
        type_col = sc_meta.columns[1]

    marker_gene_count = int(stage3_cfg.get("marker_gene_count", 30) or 30)
    min_marker_specificity = float(stage3_cfg.get("min_marker_specificity", 0.0) or 0.0)
    marker_specificity_mode = stage3_cfg.get("marker_specificity_mode", "mean")
    min_marker_genes = int(stage3_cfg.get("min_marker_genes", 2) or 2)
    fallback_specificity_mode = stage3_cfg.get("fallback_specificity_mode")
    fallback_specificity_types = stage3_cfg.get("fallback_specificity_types") or []
    if isinstance(fallback_specificity_types, str):
        fallback_specificity_types = [t.strip().lower() for t in fallback_specificity_types.split(",") if t.strip()]
    else:
        fallback_specificity_types = [str(t).strip().lower() for t in fallback_specificity_types]

    spec_overrides = parse_overrides(stage3_cfg.get("min_marker_specificity_overrides"), float)
    min_genes_overrides = parse_overrides(stage3_cfg.get("min_marker_genes_overrides"), int)
    mode_overrides = parse_overrides(stage3_cfg.get("marker_specificity_mode_overrides"), str)

    rng = np.random.RandomState(seed)

    assignment_spot = slot_assignment["assigned_spot"].astype(str).str.split(r"\s|\t").str[0]
    slot_assignment = slot_assignment.copy()
    slot_assignment["assigned_spot"] = assignment_spot

    spot_counts = slot_assignment.groupby("assigned_spot").size()
    spot_counts = spot_counts.reindex(st_expr.index, fill_value=0)
    empty_slots = (cells_per_spot - spot_counts).clip(lower=0)

    rescued_rows = []
    used_cells = set(base_assignment["cell_id"].astype(str))
    max_replace = 0
    if max_replace_fraction is not None:
        try:
            max_replace = int(len(base_assignment) * float(max_replace_fraction))
        except (TypeError, ValueError):
            max_replace = 0
    max_replace = max(0, max_replace)

    for ct in missing_types:
        key = str(ct).strip().lower()
        eff_min_spec = spec_overrides.get(key, min_marker_specificity)
        eff_min_genes = min_genes_overrides.get(key, min_marker_genes)
        eff_mode = mode_overrides.get(key, marker_specificity_mode)
        fallback_mode = (
            fallback_specificity_mode
            if (key in fallback_specificity_types or not fallback_specificity_types)
            else None
        )
        if eff_min_spec > 0:
            marker_genes = select_marker_genes_with_specificity(
                sc_expr,
                sc_meta,
                ct,
                type_col,
                plugin_genes,
                marker_gene_count,
                eff_min_spec,
                eff_mode,
                eff_min_genes,
                fallback_mode,
            )
        else:
            marker_genes = select_marker_genes(sc_expr, sc_meta, ct, type_col, plugin_genes, marker_gene_count)

        mask, counts = build_hotspot_mask(st_expr, marker_genes, expr_quantile, min_genes)
        hotspot_fraction = float(mask.mean()) if len(mask) else 0.0
        stats["hotspot_fraction_by_type"][ct] = hotspot_fraction
        if hotspot_fraction < min_cluster_fraction:
            continue

        candidate_spots = mask[mask].index
        if len(candidate_spots) == 0:
            continue

        available_cells = (
            relabel.loc[relabel["orig_type"].astype(str) == str(ct), "cell_id"].astype(str).tolist()
        )
        available_cells = [cid for cid in available_cells if cid not in used_cells]
        if not available_cells:
            continue

        weights = counts.reindex(candidate_spots).fillna(0).to_numpy(dtype=float)
        if weights.sum() <= 0:
            weights = np.ones(len(candidate_spots), dtype=float)
        weights = weights / weights.sum()

        rescued = 0
        for cid in available_cells:
            available_spots = candidate_spots[empty_slots.reindex(candidate_spots).to_numpy() > 0]
            if len(available_spots) == 0:
                break
            avail_weights = weights[: len(available_spots)]
            avail_weights = avail_weights / avail_weights.sum()
            spot = rng.choice(available_spots, p=avail_weights)
            rescued_rows.append({"cell_id": cid, "assigned_spot": spot, "cell_type": ct})
            used_cells.add(cid)
            empty_slots.loc[spot] -= 1
            rescued += 1
        if rescued > 0:
            stats["rescued_types"].append(ct)
            stats["rescued_cells_by_type"][ct] = rescued
            stats["rescued_cells_total"] += rescued

        # Optional replacement within hotspot spots (keeps total assignments constant)
        if max_replace > 0:
            available_cells = [cid for cid in available_cells if cid not in used_cells]
            if available_cells:
                candidate_rows = base_assignment.loc[
                    base_assignment["assigned_spot"].astype(str).isin(candidate_spots)
                    & (base_assignment["cell_type"].astype(str) != str(ct))
                ]
                if replace_low_confidence or replace_support_max is not None:
                    def spot_support_for_type(cell_type: str) -> pd.Series:
                        key2 = str(cell_type).strip().lower()
                        eff_min_spec2 = spec_overrides.get(key2, min_marker_specificity)
                        eff_min_genes2 = min_genes_overrides.get(key2, min_marker_genes)
                        eff_mode2 = mode_overrides.get(key2, marker_specificity_mode)
                        fallback_mode2 = (
                            fallback_specificity_mode
                            if (key2 in fallback_specificity_types or not fallback_specificity_types)
                            else None
                        )
                        if eff_min_spec2 > 0:
                            genes2 = select_marker_genes_with_specificity(
                                sc_expr,
                                sc_meta,
                                cell_type,
                                type_col,
                                plugin_genes,
                                marker_gene_count,
                                eff_min_spec2,
                                eff_mode2,
                                eff_min_genes2,
                                fallback_mode2,
                            )
                        else:
                            genes2 = select_marker_genes(sc_expr, sc_meta, cell_type, type_col, plugin_genes, marker_gene_count)
                        available = [g for g in genes2 if g in st_expr.columns]
                        if len(available) < 2:
                            return pd.Series(0.0, index=st_expr.index)
                        type_cells2 = sc_meta[sc_meta[type_col] == cell_type]["cell_id"]
                        type_cells2 = [cid for cid in type_cells2 if cid in sc_expr.index]
                        if len(type_cells2) == 0:
                            return pd.Series(0.0, index=st_expr.index)
                        type_expr2 = sc_expr.loc[type_cells2, available]
                        type_profile2 = type_expr2.mean(axis=0)
                        type_vec = type_profile2.to_numpy(dtype=float)
                        st_subset = st_expr[available].to_numpy(dtype=float)
                        type_center = type_vec - type_vec.mean()
                        st_centered = st_subset - st_subset.mean(axis=1, keepdims=True)
                        denom = (np.linalg.norm(type_center) * np.linalg.norm(st_centered, axis=1))
                        denom = np.where(denom <= 1e-10, np.inf, denom)
                        r_vals = (st_centered @ type_center) / denom
                        r_vals = np.clip(r_vals, -1.0 + 1e-10, 1.0 - 1e-10)
                        z_vals = np.arctanh(r_vals)
                        return pd.Series(z_vals, index=st_expr.index)

                    support_cache: dict[str, pd.Series] = {}
                    scores = []
                    for _, row in candidate_rows.iterrows():
                        ctype = str(row["cell_type"])
                        if ctype not in support_cache:
                            support_cache[ctype] = spot_support_for_type(ctype)
                        spot_id = str(row["assigned_spot"])
                        scores.append(float(support_cache[ctype].get(spot_id, 0.0)))
                    candidate_rows = candidate_rows.copy()
                    candidate_rows["support_score"] = scores
                    if replace_support_max is not None:
                        before_filter = len(candidate_rows)
                        candidate_rows = candidate_rows[candidate_rows["support_score"] <= float(replace_support_max)]
                        stats["replace_support_filtered"] += int(before_filter - len(candidate_rows))
                    if replace_low_confidence:
                        candidate_rows = candidate_rows.sort_values("support_score", ascending=True)

                replace_cap = min(max_replace, len(available_cells), len(candidate_rows))
                if replace_cap > 0:
                    if replace_low_confidence:
                        pick_idx = candidate_rows.index[:replace_cap]
                    else:
                        pick_idx = rng.choice(candidate_rows.index.to_numpy(), size=replace_cap, replace=False)
                    replace_cells = available_cells[:replace_cap]
                    base_assignment.loc[pick_idx, "cell_id"] = replace_cells
                    base_assignment.loc[pick_idx, "cell_type"] = ct
                    used_cells.update(replace_cells)
                    max_replace -= replace_cap
                    stats["replaced_cells_total"] += replace_cap
                    stats["replaced_cells_by_type"][ct] = stats["replaced_cells_by_type"].get(ct, 0) + replace_cap

    if rescued_rows:
        rescue_df = pd.DataFrame(rescued_rows)
        base_assignment = pd.concat([base_assignment, rescue_df], ignore_index=True)

    return base_assignment, stats


def build_by_spot_counts_from_assignment(assignment: pd.DataFrame, st_index: pd.Index) -> pd.DataFrame:
    df = assignment.copy()
    df["assigned_spot"] = df["assigned_spot"].astype(str).str.split(r"\s|\t").str[0]
    counts = (
        df.groupby(["assigned_spot", "cell_type"]).size().unstack(fill_value=0)
        if not df.empty
        else pd.DataFrame(index=[], columns=[])
    )
    counts = counts.reindex(st_index, fill_value=0)
    counts.index.name = "SpotID"
    counts["Total cells"] = counts.sum(axis=1)
    return counts


def _sim_info_candidates(root: Path, sample: str) -> list[Path]:
    out = [root / "data" / "processed" / sample / "stage1_preprocess" / "exported" / "sim_info.json"]
    out.extend([p / "sim_info.json" for p in sample_dir_candidates(root, sample, sim_group="real_brca")])
    return out


def _load_sim_info(root: Path, sample: str) -> dict | None:
    for p in _sim_info_candidates(root, sample):
        if not p.exists():
            continue
        try:
            return json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            continue
    return None


def _normalize_cli_missing_type(value: str | None) -> str | None:
    if value is None:
        return None
    v = str(value).strip()
    if not v:
        return None
    if v.upper() == "__NO_MISSING__":
        return None
    return v


def collect_missing_types_from_sim_chain(
    root: Path,
    sample: str,
    alias_map: dict[str, str],
    max_depth: int = 12,
) -> list[str]:
    """Collect all missing types from sim_info and source_sample chain."""
    cur = str(sample)
    visited: set[str] = set()
    out: list[str] = []

    for _ in range(max_depth):
        if not cur or cur in visited:
            break
        visited.add(cur)
        info = _load_sim_info(root, cur)
        if not info:
            break

        raw_list: list[str] = []
        mts = info.get("missing_types")
        if isinstance(mts, list):
            raw_list.extend([str(x) for x in mts if str(x).strip()])
        mt = info.get("missing_type")
        if mt is not None and str(mt).strip():
            raw_list.append(str(mt))

        for item in raw_list:
            canon = canonicalize_type_name(item, alias_map)
            if canon and canon not in out:
                out.append(canon)

        parent = info.get("source_sample")
        if parent is None or not str(parent).strip():
            break
        cur = str(parent).strip()

    return out


def resolve_effective_missing_types(
    root: Path,
    sample: str,
    cli_missing_type: str | None,
    alias_map: dict[str, str],
) -> list[str]:
    """Union of CLI missing_type and sim chain missing_types (canonicalized)."""
    out: list[str] = []
    cli_norm = _normalize_cli_missing_type(cli_missing_type)
    if cli_norm:
        for item in str(cli_norm).replace(";", ",").split(","):
            item = item.strip()
            if not item:
                continue
            canon = canonicalize_type_name(item, alias_map)
            if canon and canon not in out:
                out.append(canon)
    for mt in collect_missing_types_from_sim_chain(root, sample, alias_map=alias_map):
        if mt not in out:
            out.append(mt)
    return out


def collect_stage3_auto_missing_types(stage3_summary_path: Path, alias_map: dict[str, str]) -> list[str]:
    if not stage3_summary_path.exists():
        return []
    try:
        summary = json.loads(stage3_summary_path.read_text(encoding="utf-8"))
    except Exception:
        return []
    auto_obj = summary.get("auto_missing_detection") or {}
    raw_types = auto_obj.get("auto_missing_types") or []
    out: list[str] = []
    if isinstance(raw_types, list):
        for item in raw_types:
            canon = canonicalize_type_name(str(item), alias_map)
            if canon and canon not in out:
                out.append(canon)
    return out


def load_missing_truth(root: Path, sample: str) -> str | None:
    data = _load_sim_info(root, sample)
    if not data:
        return None
    mts = data.get("missing_types")
    if isinstance(mts, list) and mts:
        head = str(mts[0]).strip()
        if head:
            return head
    mt = data.get("missing_type")
    return str(mt).strip() if mt is not None and str(mt).strip() else None


def _infer_sim_scenario(sample: str) -> str | None:
    s = str(sample or "").lower()
    if "sims0" in s:
        return "S0"
    if "sims1" in s:
        return "S1"
    return None


def load_sim_truth_cell_ids(root: Path, sample: str) -> set[str] | None:
    data = _load_sim_info(root, sample)
    if not data:
        return None
    output_dir = data.get("output_dir")
    seed = data.get("seed")
    missing_type = data.get("missing_type")
    sim_sample = data.get("sample")
    if not output_dir or seed is None:
        return None
    out_dir = Path(output_dir)
    candidates = [
        out_dir / f"seed_{seed}" / "sim_truth_query_cell_spot.csv",
        out_dir / "sim_truth_query_cell_spot.csv",
    ]
    truth_path = None
    for cand in candidates:
        if cand.exists():
            truth_path = cand
            break
    if truth_path is None and sim_sample and missing_type:
        scenario = _infer_sim_scenario(sample)
        if scenario:
            mt_slug = str(missing_type).strip().lower().replace(" ", "_")
            fallback = (
                root
                / "data"
                / "sim"
                / str(sim_sample)
                / scenario
                / mt_slug
                / f"seed_{seed}"
                / "sim_truth_query_cell_spot.csv"
            )
            if fallback.exists():
                truth_path = fallback
    if truth_path is None:
        return None
    try:
        truth = pd.read_csv(truth_path, usecols=["cell_id"])
    except Exception:
        return None
    return set(truth["cell_id"].astype(str))


def load_cells_per_spot_override(root: Path, sample: str) -> int | None:
    data = _load_sim_info(root, sample)
    if not data:
        return None
    try:
        return data.get("cells_per_spot")
    except Exception:
        return None


def sha1(path: Path) -> str | None:
    if not path.exists():
        return None
    h = hashlib.sha1()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def apply_filter(
    relabel: pd.DataFrame,
    type_support: pd.DataFrame | None,
    mode: FilterMode,
    missing_types: Sequence[str] | None = None,
    filter_scope: FilterScope = "unsupported_all",
) -> tuple[pd.DataFrame, int, dict]:
    """
    根据Stage3输出决定在Stage4中要丢弃哪些细胞。

    - mode == "plugin_unknown": 使用 plugin_type == Unknown_sc_only
    - mode == "unsupported": 使用 type_support.support_category == unsupported
    - mode == "none": 不做过滤

    filter_scope 进一步约束过滤范围：
    - "unsupported_all": 按上述规则丢弃所有命中的细胞（当前默认行为）
    - "missing_only": 只丢弃 orig_type 属于 missing_type 的细胞（oracle fallback）
    - "missing_detected_only": 只丢弃同时满足“被 Stage3 标记”且 orig_type==missing_type 的细胞
    """
    missing_list = [str(x).strip() for x in (missing_types or []) if str(x).strip()]
    missing_set = set(missing_list)
    if filter_scope == "missing_only" and not missing_set:
        raise ValueError(f"filter_scope={filter_scope} requires non-empty missing_types.")
    if mode == "plugin_unknown":
        base_mask = relabel["plugin_type"] == "Unknown_sc_only"
    elif mode == "unsupported":
        col_name = "orig_type_canonical" if "orig_type_canonical" in relabel.columns else "orig_type"
        ts_col = "orig_type_canonical" if type_support is not None and "orig_type_canonical" in type_support.columns else "orig_type"
        unsupported_types = (
            set(type_support.loc[type_support["support_category"] == "unsupported", ts_col])
            if type_support is not None
            else set()
        )
        base_mask = relabel[col_name].isin(unsupported_types)
    else:
        base_mask = pd.Series(False, index=relabel.index)

    # 统计所有被Stage3标记为"低支持/Unknown"的规模（与filter_scope无关）
    col_canonical = "orig_type_canonical" if "orig_type_canonical" in relabel.columns else "orig_type"
    missing_mask = relabel[col_canonical].isin(missing_set) if missing_set else pd.Series(False, index=relabel.index)
    missing_totals_by_type = (
        relabel.loc[missing_mask, col_canonical].value_counts().to_dict() if missing_set else {}
    )
    marked_total = int(base_mask.sum())
    marked_missing = int((base_mask & missing_mask).sum())
    marked_non_missing = marked_total - marked_missing

    # 进一步按scope约束：
    # - missing_only：直接过滤指定 missing_type（oracle fallback）
    # - missing_detected_only：仅过滤被 Stage3 标记且命中 missing_type 的细胞
    if filter_scope == "missing_only" and missing_set:
        mask_drop = missing_mask
    elif filter_scope == "missing_detected_only":
        mask_drop = base_mask & missing_mask
    else:
        # 默认行为：所有 unsupported/Unknown_sc_only 都被丢弃
        mask_drop = base_mask

    filtered = relabel.loc[~mask_drop].copy()
    n_filtered = int(mask_drop.sum())
    missing_total = int(missing_mask.sum())
    oracle_only_filtered = int((missing_mask & ~base_mask).sum()) if missing_set else 0

    mark_stats = {
        "marked_unknown_total": marked_total,
        "marked_unknown_missing": marked_missing,
        "marked_unknown_non_missing": marked_non_missing,
        "missing_type_list": sorted(missing_set),
        "missing_type_total": missing_total,
        "missing_type_totals_by_type": {str(k): int(v) for k, v in missing_totals_by_type.items()},
        "missing_type_stage3_detected": bool(marked_missing > 0),
        "oracle_only_filtered_missing": oracle_only_filtered,
    }
    return filtered, n_filtered, mark_stats


def apply_effective_route2_labels(
    relabel_filtered: pd.DataFrame,
    cell_type_column: str,
    filter_scope: FilterScope,
    missing_types: Sequence[str] | None,
) -> tuple[pd.DataFrame, str, dict]:
    """
    For Route2 in missing-focused scopes, restore non-missing Unknown labels back to orig_type.

    Rationale:
    - `missing_only` / `missing_detected_only` semantics are "focus on missing type filtering".
    - Keeping non-missing cells as `Unknown_sc_only` pollutes downstream mapping and acts like
      implicit over-filtering in practice.
    """
    stats = {
        "restored_unknown_non_missing": 0,
        "kept_unknown_missing": 0,
        "effective_column": cell_type_column,
    }
    if cell_type_column != "plugin_type":
        return relabel_filtered, cell_type_column, stats
    if filter_scope not in {"missing_only", "missing_detected_only"}:
        return relabel_filtered, cell_type_column, stats
    if "plugin_type" not in relabel_filtered.columns:
        return relabel_filtered, cell_type_column, stats

    # Prefer canonicalized orig type to keep naming consistent with alias map.
    orig_col = None
    if "orig_type_canonical" in relabel_filtered.columns:
        orig_col = "orig_type_canonical"
    elif "orig_type" in relabel_filtered.columns:
        orig_col = "orig_type"
    if orig_col is None:
        return relabel_filtered, cell_type_column, stats

    unknown_mask = relabel_filtered["plugin_type"].astype(str).eq("Unknown_sc_only")
    missing_set = {str(x).strip() for x in (missing_types or []) if str(x).strip()}
    missing_mask = relabel_filtered[orig_col].astype(str).isin(missing_set) if missing_set else pd.Series(False, index=relabel_filtered.index)
    restore_mask = unknown_mask & ~missing_mask
    keep_mask = unknown_mask & missing_mask

    if not bool(restore_mask.any()):
        stats["kept_unknown_missing"] = int(keep_mask.sum())
        return relabel_filtered, cell_type_column, stats

    out = relabel_filtered.copy()
    effective_col = "plugin_type_effective"
    out[effective_col] = out["plugin_type"]
    out.loc[restore_mask, effective_col] = out.loc[restore_mask, orig_col]
    if "orig_type" in out.columns:
        out[effective_col] = out[effective_col].fillna(out["orig_type"])
    out[effective_col] = out[effective_col].fillna("Unknown_sc_only")
    out[effective_col] = out[effective_col].astype(str)

    stats["restored_unknown_non_missing"] = int(restore_mask.sum())
    stats["kept_unknown_missing"] = int(keep_mask.sum())
    stats["effective_column"] = effective_col
    return out, effective_col, stats


def prepare_cytospace_inputs(
    sc_expr: pd.DataFrame,
    st_expr: pd.DataFrame,
    st_coords: pd.DataFrame,
    sc_meta_stage1: pd.DataFrame,
    relabel_filtered: pd.DataFrame,
    work_dir: Path,
    cells_per_spot_override: int | None = None,
    cell_type_column: str = "plugin_type",
    rescue_types: list[str] | None = None,
    rescue_weight: float | None = None,
) -> dict[str, Path]:
    work_dir.mkdir(parents=True, exist_ok=True)

    # align SC expression to filtered cells
    # 如果使用sc_meta（baseline），直接使用所有SC细胞（不过滤）
    if cell_type_column == "sc_meta":
        keep_ids = [cid for cid in sc_meta_stage1["cell_id"] if cid in sc_expr.index]
    else:
        keep_ids = [cid for cid in relabel_filtered["cell_id"] if cid in sc_expr.index]
    sc_expr_filtered = sc_expr.loc[keep_ids]

    # CytoSPACE expects genes as rows, samples as columns
    sc_expr_gene_by_cell = sc_expr_filtered.T
    st_expr_gene_by_spot = st_expr.T

    sc_expr_path = work_dir / "sc_expression_for_cytospace.csv"
    st_expr_path = work_dir / "st_expression_for_cytospace.csv"
    coords_path = work_dir / "st_coordinates_for_cytospace.csv"
    cell_type_path = work_dir / "cell_types_for_cytospace.csv"
    fractions_path = work_dir / "cell_type_fractions.csv"
    n_cells_path = work_dir / "n_cells_per_spot.csv"

    # 始终覆盖写入，避免使用旧格式的缓存文件导致 CytoSPACE 报错
    print(f"[stage4] writing sc expression ({sc_expr_gene_by_cell.shape}) → {sc_expr_path.name}")
    sc_expr_gene_by_cell.to_csv(sc_expr_path, sep=",")
    print(f"[stage4] writing st expression ({st_expr_gene_by_spot.shape}) → {st_expr_path.name}")
    st_expr_gene_by_spot.to_csv(st_expr_path, sep=",")
    if "row" in st_coords.columns and "col" in st_coords.columns:
        coords_for_cyto = st_coords[["row", "col"]].copy()
    else:
        coords_for_cyto = st_coords.iloc[:, :2].copy()
    coords_for_cyto.columns = ["row", "col"]
    coords_for_cyto.to_csv(coords_path, sep=",", index_label="SpotID")

    if cell_type_column == "sc_meta":
        if "cell_type" not in sc_meta_stage1.columns:
            raise ValueError("sc_metadata.csv missing cell_type column")
        cell_type_df = sc_meta_stage1[["cell_id", "cell_type"]].copy()
    else:
        if cell_type_column not in relabel_filtered.columns:
            raise ValueError(f"cell_type_column={cell_type_column} not found in relabel file")
        cell_type_df = relabel_filtered[["cell_id", cell_type_column]].copy()
        cell_type_df.rename(columns={cell_type_column: "cell_type"}, inplace=True)
    cell_type_df.set_index("cell_id", inplace=True)
    cell_type_df.to_csv(cell_type_path, header=["cell_type"], sep=",")

    counts = cell_type_df["cell_type"].value_counts()
    if rescue_types and rescue_weight is not None and rescue_weight < 1.0:
        counts = counts.astype(float)
        for ct in rescue_types:
            if ct in counts:
                counts[ct] *= rescue_weight
    frac_series = counts / counts.sum()
    frac_df = pd.DataFrame([frac_series.values], columns=frac_series.index)
    frac_df.index = ["Fraction"]
    frac_df.to_csv(fractions_path, sep=",")

    inputs = {
        "sc_expr": sc_expr_path,
        "st_expr": st_expr_path,
        "coords": coords_path,
        "cell_types": cell_type_path,
        "fractions": fractions_path,
    }
    if cells_per_spot_override is not None:
        df_n = pd.DataFrame({"n_cells": [cells_per_spot_override] * len(st_coords)}, index=st_coords.index)
        df_n.to_csv(n_cells_path)
        inputs["n_cells_per_spot"] = n_cells_path
    return inputs


def _patch_datatable_fread() -> None:
    """Monkeypatch datatable.fread to handle non-ASCII file paths on Windows.

    datatable's C backend cannot parse paths containing Chinese/Unicode characters
    because it misinterprets the UTF-8 multi-byte sequences as escape chars.
    The fix: intercept fread calls, read file bytes via Python (which handles
    Unicode paths correctly), then pass raw bytes to the original fread.
    """
    try:
        import datatable as dt
    except ImportError:
        return

    _orig = dt.fread

    def _safe_fread(file=None, *args, **kwargs):
        if isinstance(file, (str, Path)):
            p = Path(file)
            try:
                if p.is_file() and not all(ord(c) < 128 for c in str(p)):
                    # Read bytes via Python (handles Unicode paths), pass content to fread
                    return _orig(p.read_bytes(), *args, **kwargs)
            except Exception:
                pass
        return _orig(file, *args, **kwargs)

    dt.fread = _safe_fread

    # Also patch on any already-imported cytospace submodule
    for mod_name, mod in list(sys.modules.items()):
        if "cytospace" in mod_name and hasattr(mod, "dt"):
            try:
                mod.dt.fread = _safe_fread
            except Exception:
                pass


def run_cytospace(
    inputs: dict[str, Path],
    output_dir: Path,
    project_root: Path,
    args: argparse.Namespace,
):
    # Patch datatable BEFORE importing cytospace so non-ASCII paths work
    _patch_datatable_fread()
    if not args.keep_redundant:
        os.environ["CYTOSPACE_SKIP_ASSIGNED_EXPRESSION"] = "1"
    else:
        os.environ.pop("CYTOSPACE_SKIP_ASSIGNED_EXPRESSION", None)

    sys.path.insert(0, str(project_root / "external" / "cytospace"))
    from cytospace.cytospace import main_cytospace  # type: ignore

    # Re-apply patch on freshly imported cytospace modules
    _patch_datatable_fread()

    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"[stage4] cytospace output dir: {output_dir}")
    main_cytospace(
        scRNA_path=str(inputs["sc_expr"]),
        cell_type_path=str(inputs["cell_types"]),
        n_cells_per_spot_path=str(inputs["n_cells_per_spot"]) if "n_cells_per_spot" in inputs else None,
        st_cell_type_path=None,
        cell_type_fraction_estimation_path=str(inputs["fractions"]),
        spaceranger_path=None,
        st_path=str(inputs["st_expr"]),
        coordinates_path=str(inputs["coords"]),
        output_folder=str(output_dir),
        output_prefix="",
        mean_cell_numbers=args.mean_cell_numbers,
        downsample_off=False,
        scRNA_max_transcripts_per_cell=1500,
        solver_method=args.solver_method,
        distance_metric=args.distance_metric,
        sampling_method="duplicates",
        single_cell=False,
        number_of_selected_spots=10000,
        sampling_sub_spots=args.sampling_sub_spots,
        number_of_selected_sub_spots=args.n_subspots,
        number_of_processors=args.n_processors,
        seed=args.seed,
        plot_off=True,
        geometry="honeycomb",
        max_num_cells_plot=50000,
        num_column=3,
    )


def _looks_like_non_type(cols: list[str]) -> bool:
    if not cols:
        return False
    for col in cols:
        name = str(col)
        if name.startswith("Unnamed"):
            continue
        if name.isdigit():
            continue
        return False
    return True


def _canonical_type(name: str, alias_map: dict[str, str]) -> str:
    return canonicalize_type_name(name, alias_map)


def _merge_type_columns(
    df: pd.DataFrame,
    type_cols: list[str],
    alias_map: dict[str, str],
) -> tuple[pd.DataFrame, list[str]]:
    groups: dict[str, list[str]] = {}
    for col in type_cols:
        canon = _canonical_type(col, alias_map)
        if not canon:
            continue
        groups.setdefault(canon, []).append(col)
    merged = pd.DataFrame(index=df.index)
    for canon, cols in groups.items():
        merged[canon] = df[cols].sum(axis=1)
    return merged, list(groups.keys())


def _load_fraction_table(
    path: Path,
    active_type_list: list[str],
    alias_map: dict[str, str],
) -> tuple[pd.DataFrame, list[str], str | None]:
    df = pd.read_csv(path)
    spot_col = "spot_id"
    if "spot_id" not in df.columns:
        if "SpotID" in df.columns:
            spot_col = "SpotID"
        else:
            spot_col = df.columns[0]
        df = df.rename(columns={spot_col: "spot_id"})
    df["spot_id"] = df["spot_id"].astype(str).map(lambda x: str(x).split("	")[0])

    total_col = "Total cells" if "Total cells" in df.columns else None
    type_cols = [c for c in df.columns if c not in {"spot_id", total_col}]
    if not type_cols:
        raise ValueError(f"{path} has no type columns")

    if _looks_like_non_type(type_cols):
        active = [_canonical_type(t, alias_map) for t in active_type_list if t]
        if len(active) != len(type_cols):
            raise ValueError(
                f"{path} columns look unnamed but length mismatch: "
                f"{len(type_cols)} vs active_type_list {len(active)}"
            )
        rename_map = dict(zip(type_cols, active))
        df = df.rename(columns=rename_map)
        raw_output_types = list(active)
        merged_df = df[["spot_id"] + raw_output_types].copy()
    else:
        merged_types_df, raw_output_types = _merge_type_columns(df, type_cols, alias_map)
        merged_df = pd.concat([df[["spot_id"]], merged_types_df], axis=1)

    if total_col is not None:
        merged_df["Total cells"] = pd.to_numeric(df[total_col], errors="coerce").fillna(0.0)

    return merged_df, raw_output_types, total_col


def _align_to_full_types(
    df: pd.DataFrame,
    full_type_list: list[str],
) -> pd.DataFrame:
    type_cols = [c for c in df.columns if c not in {"spot_id", "Total cells"}]
    aligned = (
        df.set_index("spot_id")[type_cols]
        .reindex(columns=full_type_list, fill_value=0.0)
        .reset_index()
    )
    if "Total cells" in df.columns:
        aligned["Total cells"] = df["Total cells"].values
    return aligned


def _mass_diff(before: pd.DataFrame, after: pd.DataFrame, type_cols: list[str]) -> dict[str, float]:
    diffs: dict[str, float] = {}
    for t in type_cols:
        if t in before.columns and t in after.columns:
            diffs[t] = float(after[t].sum() - before[t].sum())
    return diffs




def build_stage4_outputs(
    output_dir: Path,
    relabel_filtered: pd.DataFrame,
    missing_types: Sequence[str] | None,
    n_before: int,
    n_after: int,
    n_filtered: int,
    filter_mode: str,
    weak_th: float | None,
    strong_th: float | None,
    sampling_sub_spots: bool,
    n_subspots: int,
    solver_method: str,
    seed: int,
    missing_truth: str | None,
    project_root: Path,
    sample: str,
    cps_override: int | None,
    cell_type_column: str,
    filter_scope: FilterScope,
    mark_stats: dict,
    full_type_list: list[str] | None,
    active_type_list: list[str] | None,
    alias_map: dict[str, str],
) -> dict:
    assigned_path = output_dir / "assigned_locations.csv"
    if assigned_path.exists():
        assigned = pd.read_csv(assigned_path)
    else:
        assigned = pd.DataFrame()

    # Build simplified cell_assignment
    if not assigned.empty:
        assignment = assigned[["OriginalCID", "SpotID", "CellType"]].copy()
        assignment.rename(
            columns={"OriginalCID": "cell_id", "SpotID": "assigned_spot", "CellType": "cell_type"},
            inplace=True,
        )
        assignment.to_csv(output_dir / "cell_assignment.csv", index=False)
    else:
        assignment = pd.DataFrame(columns=["cell_id", "assigned_spot", "cell_type"])
        assignment.to_csv(output_dir / "cell_assignment.csv", index=False)

    # read composition to check missing type presence
    comp_path = output_dir / "cell_type_assignments_by_spot.csv"
    frac_path = output_dir / "fractional_abundances_by_spot.csv"
    missing_present = 0
    comp_total = None
    raw_output_types = []
    raw_fractional_types = []
    type_alignment_ok = None
    type_alignment_diff_max = None
    missing_keys: list[str] = []
    for mt in (missing_types or []):
        canon = canonicalize_type_name(mt, alias_map)
        if canon and canon not in missing_keys:
            missing_keys.append(canon)
    missing_present_by_type: dict[str, int] = {}

    if comp_path.exists():
        comp_raw, raw_output_types, _ = _load_fraction_table(
            comp_path,
            active_type_list or [],
            alias_map,
        )
        comp_aligned = comp_raw
        if full_type_list:
            comp_aligned = _align_to_full_types(comp_raw, full_type_list)
            diffs = _mass_diff(comp_raw, comp_aligned, raw_output_types)
            diff_max = max((abs(v) for v in diffs.values()), default=0.0)
            type_alignment_diff_max = diff_max
            if diff_max > 1e-6:
                raise ValueError("type alignment mismatch in cell_type_assignments_by_spot.csv")
        comp_aligned.to_csv(comp_path, index=False)
        for mk in missing_keys:
            if mk in comp_aligned.columns:
                missing_present_by_type[mk] = int(comp_aligned[mk].sum())
        missing_present = int(sum(missing_present_by_type.values()))
        if "Total cells" in comp_aligned.columns:
            comp_total = int(comp_aligned["Total cells"].sum())

    if frac_path.exists():
        frac_raw, raw_fractional_types, _ = _load_fraction_table(
            frac_path,
            active_type_list or [],
            alias_map,
        )
        frac_aligned = frac_raw
        if full_type_list:
            frac_aligned = _align_to_full_types(frac_raw, full_type_list)
            diffs = _mass_diff(frac_raw, frac_aligned, raw_fractional_types)
            diff_max = max((abs(v) for v in diffs.values()), default=0.0)
            if type_alignment_diff_max is None:
                type_alignment_diff_max = diff_max
            else:
                type_alignment_diff_max = max(type_alignment_diff_max, diff_max)
            if diff_max > 1e-6:
                raise ValueError("type alignment mismatch in fractional_abundances_by_spot.csv")
        frac_aligned.to_csv(frac_path, index=False)

    if active_type_list is not None:
        active_set = set(active_type_list)
        raw_set = set(raw_output_types) if raw_output_types else None
        if raw_set is not None:
            type_alignment_ok = (raw_set == active_set)

    ca_rows = len(assignment)

    summary = {
        "n_cells_before_prefilter": int(n_before),
        "n_cells_after_prefilter": int(n_after),
        "n_filtered": int(n_filtered),
        "assignment_rows": int(ca_rows),
        "assigned_locations_rows": int(len(assigned)) if not assigned.empty else 0,
        "by_spot_total_cells": comp_total,
        "missing_present_in_assignment": int(missing_present),
        "missing_present_in_assignment_by_type": missing_present_by_type,
        "filter_mode": filter_mode,
        "stage3_weak_th": weak_th,
        "stage3_strong_th": strong_th,
        "sampling_sub_spots": bool(sampling_sub_spots),
        "n_subspots": int(n_subspots),
        "solver_method": solver_method,
        "seed": int(seed),
        "missing_type_truth": missing_truth,
        "missing_type_list": missing_keys,
        "cells_per_spot_override": cps_override,
        "cell_type_column": cell_type_column,
        "filter_scope": filter_scope,
        "type_list": full_type_list,
        "full_type_list": full_type_list,
        "active_type_list": active_type_list,
        "raw_output_types_by_spot": raw_output_types,
        "raw_output_types_fractional": raw_fractional_types,
        "type_alignment_ok": type_alignment_ok,
        "type_alignment_diff_max": type_alignment_diff_max,
        # 这些统计描述的是：在Stage3视角下被标记为“低支持/Unknown”的规模
        # 注意：与实际被Stage4丢弃的 n_filtered 相互独立
        "marked_unknown_total": mark_stats.get("marked_unknown_total"),
        "marked_unknown_missing": mark_stats.get("marked_unknown_missing"),
        "marked_unknown_non_missing": mark_stats.get("marked_unknown_non_missing"),
        "missing_type_total_in_pool": mark_stats.get("missing_type_total"),
        "missing_type_total_in_pool_by_type": mark_stats.get("missing_type_totals_by_type"),
        "missing_type_stage3_detected": mark_stats.get("missing_type_stage3_detected"),
        "oracle_only_filtered_missing": mark_stats.get("oracle_only_filtered_missing"),
        "sha1": {
            "cell_assignment": sha1(output_dir / "cell_assignment.csv"),
            "sim_info": sha1(project_root / "data" / "processed" / sample / "stage1_preprocess" / "exported" / "sim_info.json"),
            "sc_metadata": sha1(project_root / "data" / "processed" / sample / "stage1_preprocess" / "exported" / "sc_metadata.csv"),
        },
        "cell_assignment_path": str(output_dir / "cell_assignment.csv"),
        "cytospace_output": str(output_dir),
    }
    (output_dir / "stage4_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return summary


def main():
    args = parse_args()

    # 强校验：依赖 missing_type 的 scope 必须提供 non-empty missing_type

    project_root = Path(args.project_root) if args.project_root else detect_root()

    dataset_cfg = load_dataset_cfg(project_root, args.sample)
    stage3_cfg = (dataset_cfg.get("stage3") or {}) if isinstance(dataset_cfg, dict) else {}
    stage4_cfg = (dataset_cfg.get("stage4") or {}) if isinstance(dataset_cfg, dict) else {}

    cfg_truth_filter = load_stage4_truth_filter(project_root, args.sample)
    filter_to_sim_truth = bool(args.filter_to_sim_truth or cfg_truth_filter)

    hotspot_rescue_post = bool(args.hotspot_rescue_post or stage4_cfg.get("hotspot_rescue_post", False))
    qvalue_max = (
        args.hotspot_rescue_qvalue_max
        if args.hotspot_rescue_qvalue_max is not None
        else float(stage4_cfg.get("hotspot_rescue_qvalue_max", 0.1))
    )
    qvalue_min_cfg = (
        args.hotspot_rescue_qvalue_min
        if args.hotspot_rescue_qvalue_min is not None
        else stage4_cfg.get("hotspot_rescue_qvalue_min", "final_fdr")
    )
    stage3_summary_path = project_root / "result" / args.sample / "stage3_typematch" / "stage3_summary.json"
    final_fdr = resolve_final_fdr(stage3_cfg, stage3_summary_path)
    qvalue_min = resolve_qvalue_min(qvalue_min_cfg, final_fdr)
    rescue_min_genes = (
        args.hotspot_rescue_min_genes
        if args.hotspot_rescue_min_genes is not None
        else int(stage4_cfg.get("hotspot_min_genes", 3))
    )
    rescue_expr_quantile = (
        args.hotspot_rescue_expr_quantile
        if args.hotspot_rescue_expr_quantile is not None
        else float(stage4_cfg.get("hotspot_expr_quantile", 0.9))
    )
    rescue_min_cluster_fraction = (
        args.hotspot_rescue_min_cluster_fraction
        if args.hotspot_rescue_min_cluster_fraction is not None
        else float(stage4_cfg.get("hotspot_min_cluster_fraction", 0.01))
    )
    rescue_dedupe = bool(args.hotspot_rescue_dedupe or stage4_cfg.get("hotspot_rescue_dedupe", False))
    allow_nonsignificant = bool(
        args.hotspot_rescue_allow_nonsignificant
        or stage4_cfg.get("hotspot_rescue_allow_nonsignificant", False)
    )
    max_replace_fraction = (
        args.hotspot_rescue_max_replace_fraction
        if args.hotspot_rescue_max_replace_fraction is not None
        else stage4_cfg.get("hotspot_rescue_max_replace_fraction")
    )
    replace_low_confidence = bool(
        args.hotspot_rescue_replace_low_confidence
        or stage4_cfg.get("hotspot_rescue_replace_low_confidence", False)
    )
    replace_support_max = (
        args.hotspot_rescue_replace_support_max
        if args.hotspot_rescue_replace_support_max is not None
        else stage4_cfg.get("hotspot_rescue_replace_support_max")
    )
    if replace_support_max is not None:
        replace_support_max = float(replace_support_max)

    rescue_ctrl_cfg = stage3_cfg.get("V5_rescue_control") or stage3_cfg.get("v5_rescue_control") or {}
    rescue_ctrl_enable = bool(rescue_ctrl_cfg.get("enable", False))
    try:
        rescue_weight = float(rescue_ctrl_cfg.get("prior_weight", 1.0))
    except (TypeError, ValueError):
        rescue_weight = 1.0
    if rescue_weight < 0:
        rescue_weight = 0.0
    if rescue_weight > 1:
        rescue_weight = 1.0
    rescued_types = []
    if rescue_ctrl_enable and stage3_summary_path.exists():
        try:
            summary = json.loads(stage3_summary_path.read_text(encoding="utf-8"))
            rescued_types = summary.get("v5_entropy_qc", {}).get("rescued_types_final", []) or []
        except Exception:
            rescued_types = []

    # 加载别名映射，并对 missing_type 做 canonicalize
    alias_map_path = project_root / "configs" / "type_aliases.yaml"
    alias_map = load_alias_map(alias_map_path)
    cli_missing_norm = _normalize_cli_missing_type(args.missing_type)
    missing_type_canonical = canonicalize_type_name(cli_missing_norm, alias_map) if cli_missing_norm else None
    effective_missing_types = resolve_effective_missing_types(
        project_root,
        args.sample,
        args.missing_type,
        alias_map,
    )
    if args.filter_scope == "missing_detected_only":
        for mt in collect_stage3_auto_missing_types(stage3_summary_path, alias_map):
            if mt not in effective_missing_types:
                effective_missing_types.append(mt)

    # 日志中显式打印 canonicalized missing_type 以及是否命中 alias
    if cli_missing_norm is None:
        print(f"[stage4] missing_type input disabled (raw='{args.missing_type}')")
    else:
        orig_norm = normalize_type_name(cli_missing_norm)
        if missing_type_canonical != orig_norm:
            print(
                f"[stage4] missing_type canonicalized: "
                f"raw='{args.missing_type}' -> normalized='{orig_norm}' -> canonical='{missing_type_canonical}' (via alias map)"
            )
        else:
            print(
                f"[stage4] missing_type canonicalized: "
                f"raw='{args.missing_type}' -> normalized='{orig_norm}' (no alias mapping applied)"
            )
    print(f"[stage4] effective_missing_types={effective_missing_types if effective_missing_types else []}")
    if args.filter_scope == "missing_only" and not effective_missing_types:
        raise SystemExit(
            f"filter_scope={args.filter_scope} requires at least one missing type "
            f"(from --missing_type or sim_info(.missing_types/.missing_type) chain)."
        )

    sc_expr, st_expr, st_coords, sc_meta = load_stage1(
        project_root,
        args.sample,
        sc_expr_source=args.sc_expr_source,
    )
    print(f"[stage4] sc_expr_source={args.sc_expr_source}")

    full_type_list = []
    for t in sc_meta["cell_type"].dropna().tolist():
        canon = canonicalize_type_name(t, alias_map)
        if canon and canon not in full_type_list:
            full_type_list.append(canon)
    if stage3_summary_path.exists():
        try:
            summary = json.loads(stage3_summary_path.read_text(encoding="utf-8"))
            plugin_types = summary.get("plugin_types", []) or []
            for t in plugin_types:
                canon = canonicalize_type_name(t, alias_map)
                if canon and canon not in full_type_list:
                    full_type_list.append(canon)
        except Exception:
            pass
    print(f"[debug] sc_expr shape={sc_expr.shape}, st_expr shape={st_expr.shape}")
    print(f"[debug] st_expr index sample={list(st_expr.index[:3])}")
    print(f"[debug] st_expr columns sample={list(st_expr.columns[:3])}")
    print(f"[debug] st_coords index sample={list(st_coords.index[:3])}")
    
    # Baseline (filter_mode=none) 不需要Stage3输出，直接使用所有SC细胞
    type_support = None
    if args.filter_mode == "none":
        # 对于baseline，使用所有SC细胞，不进行任何过滤
        n_before = len(sc_meta)
        n_after = len(sc_meta)
        n_filtered = 0
        mark_stats = {
            "marked_unknown_total": 0,
            "marked_unknown_missing": 0,
            "marked_unknown_non_missing": 0,
            "missing_type_list": list(effective_missing_types),
            "missing_type_total": 0,
            "missing_type_totals_by_type": {},
            "missing_type_stage3_detected": False,
            "oracle_only_filtered_missing": 0,
        }
        # 创建一个虚拟的relabel_filtered，用于prepare_cytospace_inputs
        # 但实际上cell_type_column=sc_meta时，不会使用这个
        relabel_filtered = sc_meta[["cell_id"]].copy()
    else:
        # Route2需要Stage3输出进行过滤
        relabel, type_support = load_stage3(project_root, args.sample)

        # 对 orig_type 做 canonicalize，便于后续与 missing_type 一致匹配
        if "orig_type" in relabel.columns:
            relabel["orig_type_canonical"] = relabel["orig_type"].apply(
                lambda x: canonicalize_type_name(x, alias_map)
            )
        if type_support is not None and "orig_type" in type_support.columns:
            type_support["orig_type_canonical"] = type_support["orig_type"].apply(
                lambda x: canonicalize_type_name(x, alias_map)
            )
        n_before = len(relabel)

        relabel_filtered, n_filtered, mark_stats = apply_filter(
            relabel,
            type_support,
            args.filter_mode,
            missing_types=effective_missing_types,
            filter_scope=args.filter_scope,
        )
        n_after = len(relabel_filtered)

    truth_filter_removed = 0
    if filter_to_sim_truth:
        truth_cell_ids = load_sim_truth_cell_ids(project_root, args.sample)
        if truth_cell_ids:
            before_truth = len(relabel_filtered)
            relabel_filtered = relabel_filtered[relabel_filtered["cell_id"].astype(str).isin(truth_cell_ids)].copy()
            truth_filter_removed = before_truth - len(relabel_filtered)
            n_after = len(relabel_filtered)
            n_filtered = int(n_filtered) + int(truth_filter_removed)

    relabel_for_mapping = relabel_filtered
    mapping_cell_type_column = args.cell_type_column
    route2_restore_stats = {
        "restored_unknown_non_missing": 0,
        "kept_unknown_missing": 0,
        "effective_column": args.cell_type_column,
    }
    if args.filter_mode != "none":
        relabel_for_mapping, mapping_cell_type_column, route2_restore_stats = apply_effective_route2_labels(
            relabel_filtered,
            args.cell_type_column,
            args.filter_scope,
            effective_missing_types,
        )

    active_type_list = []
    if args.cell_type_column == "sc_meta":
        source_series = sc_meta["cell_type"] if "cell_type" in sc_meta.columns else []
    else:
        if mapping_cell_type_column in relabel_for_mapping.columns:
            source_series = relabel_for_mapping[mapping_cell_type_column]
        else:
            source_series = sc_meta.get("cell_type", [])
    if hasattr(source_series, "dropna"):
        source_series = source_series.dropna().tolist()
    for t in source_series:
        canon = canonicalize_type_name(t, alias_map)
        if canon and canon not in active_type_list:
            active_type_list.append(canon)

    stage4_dir_name = "stage4_cytospace" + (args.stage4_suffix if args.stage4_suffix else "")
    out_root = project_root / "result" / args.sample / stage4_dir_name
    prep_dir = out_root / "cytospace_input"
    cyto_out_dir = out_root / "cytospace_output"

    cps_override = load_cells_per_spot_override(project_root, args.sample)
    mapping_cps = args.mapping_cells_per_spot if args.mapping_cells_per_spot is not None else cps_override
    if mapping_cps is None:
        mapping_cps = args.mean_cell_numbers
    inputs = prepare_cytospace_inputs(
        sc_expr,
        st_expr,
        st_coords,
        sc_meta,
        relabel_for_mapping,
        prep_dir,
        cells_per_spot_override=mapping_cps,
        cell_type_column=mapping_cell_type_column,
        rescue_types=rescued_types if rescue_ctrl_enable else None,
        rescue_weight=rescue_weight if rescue_ctrl_enable else None,
    )
    run_cytospace(inputs, cyto_out_dir, project_root, args)
    weak_th, strong_th = load_stage3_thresholds(project_root, args.sample)
    missing_truth = load_missing_truth(project_root, args.sample)
    summary = build_stage4_outputs(
        cyto_out_dir,
        relabel_filtered,
        effective_missing_types,
        n_before,
        n_after,
        n_filtered,
        args.filter_mode,
        weak_th,
        strong_th,
        args.sampling_sub_spots,
        args.n_subspots,
        args.solver_method,
        args.seed,
        missing_truth,
        project_root,
        args.sample,
        cps_override,
        args.cell_type_column,
        args.filter_scope,
        mark_stats,
        full_type_list,
        active_type_list,
        alias_map,
    )
    summary["mapping_cells_per_spot"] = int(mapping_cps) if mapping_cps is not None else None
    summary["truth_filter_enabled"] = bool(filter_to_sim_truth)
    summary["truth_filter_removed"] = int(truth_filter_removed)
    summary["mapping_cell_type_column"] = mapping_cell_type_column
    summary["route2_restored_unknown_non_missing"] = int(route2_restore_stats.get("restored_unknown_non_missing", 0))
    summary["route2_kept_unknown_missing"] = int(route2_restore_stats.get("kept_unknown_missing", 0))

    if hotspot_rescue_post and args.filter_mode != "none":
        stage1_export = project_root / "data" / "processed" / args.sample / "stage1_preprocess" / "exported"
        rescue_cells_per_spot = (
            args.hotspot_rescue_cells_per_spot
            if args.hotspot_rescue_cells_per_spot is not None
            else cps_override
        )
        if rescue_cells_per_spot is None:
            rescue_cells_per_spot = mapping_cps or args.mean_cell_numbers
        cells_per_spot = rescue_cells_per_spot
        assignment_post, rescue_stats = prepare_post_rescue_assignment(
            cyto_out_dir / "cell_assignment.csv",
            relabel,
            sc_expr,
            sc_meta,
            st_expr,
            stage1_export,
            stage3_cfg,
            qvalue_min,
            qvalue_max,
            rescue_expr_quantile,
            rescue_min_genes,
            rescue_min_cluster_fraction,
            int(cells_per_spot),
            args.seed,
            rescue_dedupe,
            allow_nonsignificant,
            max_replace_fraction,
            replace_low_confidence,
            replace_support_max,
        )
        post_assign_path = cyto_out_dir / "cell_assignment_post_rescue.csv"
        assignment_post.to_csv(post_assign_path, index=False)
        post_by_spot_path = cyto_out_dir / "cell_type_assignments_by_spot_post_rescue.csv"
        post_by_spot = build_by_spot_counts_from_assignment(assignment_post, st_expr.index)
        post_by_spot.to_csv(post_by_spot_path)

        rescue_stats["qvalue_window"] = [qvalue_min, qvalue_max]
        rescue_stats["expr_quantile"] = rescue_expr_quantile
        rescue_stats["min_genes"] = int(rescue_min_genes)
        rescue_stats["min_cluster_fraction"] = float(rescue_min_cluster_fraction)
        rescue_stats["cells_per_spot"] = int(cells_per_spot)
        rescue_stats["mapping_cells_per_spot"] = int(mapping_cps) if mapping_cps is not None else None
        rescue_stats["replace_support_max"] = replace_support_max

        summary["hotspot_rescue_post"] = rescue_stats
        summary["post_rescue_assignment_path"] = str(post_assign_path)
        summary["post_rescue_by_spot_path"] = str(post_by_spot_path)

    (cyto_out_dir / "stage4_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    print("[Stage4 cytospace] summary:")
    print(json.dumps(summary, ensure_ascii=False, indent=2))

    # 清理冗余文件，避免重复占用存储（assigned_expression、sc/st_expression_for_cytospace 未被 Stage5 使用）
    if not args.keep_redundant:
        _cleanup_redundant_stage4_outputs(cyto_out_dir, prep_dir)


def _cleanup_redundant_stage4_outputs(cyto_out_dir: Path, prep_dir: Path) -> None:
    """删除 Stage5 不依赖的冗余输出，节省约 6GB/run。"""
    removed = []
    ae = cyto_out_dir / "assigned_expression"
    if ae.exists() and ae.is_dir():
        shutil.rmtree(ae)
        removed.append("assigned_expression")
    for name in ("sc_expression_for_cytospace.csv", "st_expression_for_cytospace.csv"):
        for d in (cyto_out_dir, prep_dir):
            if d and (d / name).exists():
                (d / name).unlink()
                removed.append(name)
                break
    if removed:
        print(f"[Stage4] Cleaned redundant: {removed}")


if __name__ == "__main__":
    main()
