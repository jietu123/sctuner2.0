"""
Stage6 real-data audit (no ground truth).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable

import numpy as np
import pandas as pd

try:
    import yaml
except Exception:  # pragma: no cover
    yaml = None

# ensure project root in sys.path
_HERE = Path(__file__).resolve()
_ROOT = _HERE.parents[2]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src.utils.type_name import normalize_type_name

EPS = 1e-12


@dataclass
class FractionData:
    df: pd.DataFrame
    spot_ids: pd.Series
    type_cols: list[str]
    original_type_cols: list[str]
    filled_type_cols: list[str]
    counts: np.ndarray
    counts_total: np.ndarray | None
    fractions: np.ndarray
    row_sum_counts: np.ndarray
    row_sum_fractions: np.ndarray
    fraction_source: str


def stats_dict(values: Iterable[float]) -> Dict[str, float]:
    arr = np.asarray(list(values), dtype=float)
    if arr.size == 0:
        return {
            "min": float("nan"),
            "p05": float("nan"),
            "median": float("nan"),
            "p95": float("nan"),
            "max": float("nan"),
        }
    return {
        "min": float(np.nanmin(arr)),
        "p05": float(np.nanpercentile(arr, 5)),
        "median": float(np.nanmedian(arr)),
        "p95": float(np.nanpercentile(arr, 95)),
        "max": float(np.nanmax(arr)),
    }


def compute_entropy(fractions: np.ndarray, eps: float = EPS) -> np.ndarray:
    safe = np.clip(fractions, eps, 1.0)
    return -(safe * np.log(safe)).sum(axis=1)


def js_divergence(p: np.ndarray, q: np.ndarray, eps: float = EPS) -> float:
    p = np.clip(p, eps, 1.0)
    q = np.clip(q, eps, 1.0)
    p = p / p.sum()
    q = q / q.sum()
    m = 0.5 * (p + q)
    kl_pm = np.sum(p * np.log(p / m))
    kl_qm = np.sum(q * np.log(q / m))
    return float(0.5 * (kl_pm + kl_qm))


def _read_expression_wide(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, nrows=0)
    columns = list(df.columns)
    if "spot_id" in columns:
        full = pd.read_csv(path)
        full["spot_id"] = full["spot_id"].astype(str)
        return full

    def _looks_like_spot_ids(cols: list[str]) -> bool:
        sample = [str(c) for c in cols[: min(len(cols), 50)]]
        hits = sum(1 for c in sample if "__" in c)
        return hits >= max(1, len(sample) // 4)

    if columns and columns[0] == "Unnamed: 0":
        if _looks_like_spot_ids(columns[1:]):
            mat = pd.read_csv(path, index_col=0)
            mat.index = mat.index.astype(str)
            mat = mat.T
            mat.insert(0, "spot_id", mat.index.astype(str))
            return mat.reset_index(drop=True)
        full = pd.read_csv(path)
        full = full.rename(columns={columns[0]: "spot_id"})
        full["spot_id"] = full["spot_id"].astype(str)
        return full

    if columns and columns[0] in {"gene", "Gene", "gene_id"}:
        mat = pd.read_csv(path, index_col=0)
        mat.index = mat.index.astype(str)
        mat = mat.T
        mat.insert(0, "spot_id", mat.index.astype(str))
        return mat.reset_index(drop=True)

    full = pd.read_csv(path)
    if full.columns.size:
        full = full.rename(columns={full.columns[0]: "spot_id"})
        full["spot_id"] = full["spot_id"].astype(str)
    return full


def _spearman_corr_series(a: np.ndarray, b: np.ndarray) -> float:
    if a.size < 2 or b.size < 2:
        return float("nan")
    ra = pd.Series(a).rank(method="average").to_numpy()
    rb = pd.Series(b).rank(method="average").to_numpy()
    if np.std(ra) == 0 or np.std(rb) == 0:
        return float("nan")
    return float(np.corrcoef(ra, rb)[0, 1])


def _select_top_variance_genes(expr: pd.DataFrame, n_genes: int) -> list[str]:
    cols = [c for c in expr.columns if c != "spot_id"]
    if not cols:
        return []
    mat = expr[cols].to_numpy(dtype=float)
    var = np.var(mat, axis=0)
    order = np.argsort(var)[::-1]
    top = [cols[i] for i in order[: min(n_genes, len(cols))]]
    return top


def _smooth_values(values: np.ndarray, coords: np.ndarray, k: int, alpha: float) -> np.ndarray:
    if k <= 0 or alpha <= 0:
        return values
    n = coords.shape[0]
    out = values.copy()
    # pairwise distances
    diff_x = coords[:, 0][:, None] - coords[:, 0][None, :]
    diff_y = coords[:, 1][:, None] - coords[:, 1][None, :]
    dist = diff_x * diff_x + diff_y * diff_y
    for i in range(n):
        order = np.argsort(dist[i])
        neigh = order[1 : k + 1]
        if neigh.size == 0:
            continue
        out[i] = (1 - alpha) * values[i] + alpha * float(np.mean(values[neigh]))
    return out


def _load_sc_expression(path: Path, cell_id_col: str | None = None) -> pd.DataFrame:
    df = pd.read_csv(path)
    if cell_id_col is None:
        cell_id_col = "cell_id" if "cell_id" in df.columns else df.columns[0]
    df = df.set_index(cell_id_col)
    return df


def _build_type_profiles(sc_expr: pd.DataFrame, sc_meta: pd.DataFrame, type_col: str) -> pd.DataFrame:
    profiles = []
    for t, group in sc_meta.groupby(type_col):
        ids = group["cell_id"].astype(str).tolist()
        ids = [cid for cid in ids if cid in sc_expr.index]
        if not ids:
            continue
        mean_expr = sc_expr.loc[ids].mean(axis=0)
        profiles.append((t, mean_expr))
    if not profiles:
        return pd.DataFrame()
    types = [p[0] for p in profiles]
    mat = pd.DataFrame([p[1].values for p in profiles], columns=sc_expr.columns, index=types)
    return mat


def _align_types_for_profiles(types: list[str], type_cols: list[str]) -> dict[str, str]:
    norm_map = {normalize_type_name(c): c for c in type_cols}
    mapping = {}
    for t in types:
        norm = normalize_type_name(t)
        if norm in norm_map:
            mapping[t] = norm_map[norm]
    return mapping


def run_reconstruction_checks(
    fraction_data: FractionData,
    st_expr: pd.DataFrame,
    sc_profiles: pd.DataFrame,
    gene_set: str,
    n_genes: int,
    null_shuffle: int,
    corr_delta_min: float,
    rmse_delta_min: float,
) -> dict[str, Any]:
    reasons: list[str] = []
    status = "WARN"
    if sc_profiles.empty or st_expr.empty:
        return {"status": "SKIPPED", "metrics": {}, "reasons": ["missing profiles or st expression"]}

    # align genes
    st_genes = [c for c in st_expr.columns if c != "spot_id"]
    common_genes = [g for g in st_genes if g in sc_profiles.columns]
    if gene_set == "st_hvg2000":
        st_sub = st_expr[["spot_id"] + common_genes]
        top_genes = _select_top_variance_genes(st_sub, n_genes)
        common_genes = [g for g in top_genes if g in sc_profiles.columns]

    if not common_genes:
        return {"status": "SKIPPED", "metrics": {}, "reasons": ["no common genes"]}

    # align types
    type_cols = fraction_data.type_cols
    type_map = _align_types_for_profiles(sc_profiles.index.tolist(), type_cols)
    if not type_map:
        return {"status": "SKIPPED", "metrics": {}, "reasons": ["no common types for profiles"]}

    prof = sc_profiles.loc[list(type_map.keys()), common_genes]
    prof.index = [type_map[t] for t in prof.index]
    prof = prof.groupby(prof.index).mean()

    # build F and E
    frac_df = fraction_data.df.set_index("spot_id")[type_cols]
    st_expr = st_expr.set_index("spot_id")
    common_spots = frac_df.index.intersection(st_expr.index)
    if common_spots.empty:
        return {"status": "SKIPPED", "metrics": {}, "reasons": ["no common spots"]}
    F = frac_df.loc[common_spots, prof.index].to_numpy(dtype=float)
    E = st_expr.loc[common_spots, common_genes].to_numpy(dtype=float)
    P = prof.to_numpy(dtype=float)

    # reconstruct
    E_hat = F @ P

    def _corr_rows(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        a_center = a - a.mean(axis=1, keepdims=True)
        b_center = b - b.mean(axis=1, keepdims=True)
        denom = np.sqrt((a_center ** 2).sum(axis=1) * (b_center ** 2).sum(axis=1))
        denom = np.where(denom == 0, np.nan, denom)
        corr = (a_center * b_center).sum(axis=1) / denom
        return corr

    def _corr_cols(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        a_center = a - a.mean(axis=0, keepdims=True)
        b_center = b - b.mean(axis=0, keepdims=True)
        denom = np.sqrt((a_center ** 2).sum(axis=0) * (b_center ** 2).sum(axis=0))
        denom = np.where(denom == 0, np.nan, denom)
        corr = (a_center * b_center).sum(axis=0) / denom
        return corr

    spot_corr_real = _corr_rows(E, E_hat)
    gene_corr_real = _corr_cols(E, E_hat)
    rmse_real = float(np.mean((E - E_hat) ** 2) / (np.mean(E ** 2) + EPS))

    # null
    spot_corr_null = []
    rmse_null = []
    rng = np.random.default_rng(42)
    for _ in range(int(null_shuffle)):
        perm = rng.permutation(F.shape[0])
        E_hat_null = F[perm] @ P
        sc = _corr_rows(E, E_hat_null)
        spot_corr_null.append(np.nanmedian(sc) if sc.size else float("nan"))
        rmse_null.append(float(np.mean((E - E_hat_null) ** 2) / (np.mean(E ** 2) + EPS)))

    spot_corr_null_med = float(np.nanmedian(spot_corr_null)) if spot_corr_null else float("nan")
    rmse_null_med = float(np.nanmedian(rmse_null)) if rmse_null else float("nan")

    spot_corr_real_med = float(np.nanmedian(spot_corr_real)) if spot_corr_real.size else float("nan")
    gene_corr_real_med = float(np.nanmedian(gene_corr_real)) if gene_corr_real.size else float("nan")

    corr_delta = spot_corr_real_med - spot_corr_null_med
    rmse_delta = rmse_null_med - rmse_real

    if (corr_delta >= corr_delta_min) or (rmse_delta >= rmse_delta_min):
        status = "OK"
    else:
        status = "WARN"
        reasons.append("reconstruction not better than null")

    metrics = {
        "spot_corr_median": spot_corr_real_med,
        "gene_corr_median": gene_corr_real_med,
        "rmse_relative": rmse_real,
        "null_spot_corr_median": spot_corr_null_med,
        "null_rmse_relative": rmse_null_med,
        "corr_delta": corr_delta,
        "rmse_delta": rmse_delta,
        "n_genes_used": int(len(common_genes)),
        "n_spots_used": int(len(common_spots)),
        "corr_method": "pearson",
    }

    return {"status": status, "metrics": metrics, "reasons": reasons}


def run_compartment_checks(
    fraction_data: FractionData,
    sc_expr: pd.DataFrame,
    sc_meta: pd.DataFrame,
    st_expr: pd.DataFrame,
    st_coords: pd.DataFrame,
    config: dict[str, Any],
) -> dict[str, Any]:
    if st_expr.empty or sc_expr.empty or sc_meta.empty:
        return {"status": "SKIPPED", "metrics": {}, "reasons": ["missing inputs for compartments"]}

    comp_cfg = config.get("compartments", {}) if isinstance(config, dict) else {}
    mapping = comp_cfg.get("mapping", []) if isinstance(comp_cfg, dict) else []
    if not mapping:
        return {"status": "SKIPPED", "metrics": {}, "reasons": ["compartments mapping missing"]}

    n_genes = int(comp_cfg.get("n_genes_per_compartment", 200))
    detect_min = float(comp_cfg.get("detect_rate_min", 0.05))
    k = int(comp_cfg.get("smoothing_k", 6))
    alpha = float(comp_cfg.get("smoothing_alpha", 0.5))

    # coords
    coord_col = pick_column(st_coords, ["spot_id", "SpotID", "spot", "Spot"]) if st_coords is not None else None
    if coord_col is None:
        return {"status": "SKIPPED", "metrics": {}, "reasons": ["st_coords missing"]}
    coords_df = st_coords.copy()
    coords_df[coord_col] = coords_df[coord_col].astype(str)
    if "row" in coords_df.columns and "col" in coords_df.columns:
        coords_df = coords_df[[coord_col, "row", "col"]].copy()
        coords_df["row"] = pd.to_numeric(coords_df["row"], errors="coerce")
        coords_df["col"] = pd.to_numeric(coords_df["col"], errors="coerce")
    else:
        coord_cols = [c for c in coords_df.columns if c != coord_col]
        if len(coord_cols) < 2:
            return {"status": "SKIPPED", "metrics": {}, "reasons": ["st_coords missing row/col"]}
        coords_df = coords_df[[coord_col, coord_cols[0], coord_cols[1]]].copy()
        coords_df[coord_cols[0]] = pd.to_numeric(coords_df[coord_cols[0]], errors="coerce")
        coords_df[coord_cols[1]] = pd.to_numeric(coords_df[coord_cols[1]], errors="coerce")
        coords_df.columns = [coord_col, "row", "col"]

    st_expr = st_expr.set_index("spot_id")
    frac_df = fraction_data.df.set_index("spot_id")
    common_spots = frac_df.index.intersection(st_expr.index)
    if common_spots.empty:
        return {"status": "SKIPPED", "metrics": {}, "reasons": ["no common spots"]}

    # align coords
    coords_idx = coords_df.set_index(coord_col).loc[common_spots]
    coords = coords_idx.to_numpy(dtype=float)

    # compute detect_rate
    st_genes = [c for c in st_expr.columns if c != "spot_id"]
    st_mat = st_expr.loc[common_spots, st_genes].to_numpy(dtype=float)
    detect_rate = (st_mat > 0).mean(axis=0)
    detect_map = dict(zip(st_genes, detect_rate))

    sc_meta = sc_meta.copy()
    type_col = None
    for cand in ["cell_type", "CellType", "sc_type", "type"]:
        if cand in sc_meta.columns:
            type_col = cand
            break
    if type_col is None:
        return {"status": "SKIPPED", "metrics": {}, "reasons": ["sc_metadata missing cell_type"]}

    sc_meta["_type"] = sc_meta[type_col].astype(str)

    results = []
    strong_count = 0
    strong_ok = 0
    weak_count = 0

    for entry in mapping:
        if not isinstance(entry, dict):
            continue
        name = entry.get("name")
        types = entry.get("types", [])
        if not name or not types:
            continue

        type_norm = [normalize_type_name(t) for t in types]
        # fraction sum
        cols = [c for c in fraction_data.type_cols if normalize_type_name(c) in type_norm]
        if not cols:
            results.append({"compartment": name, "status": "SKIPPED", "reason": "types not found"})
            continue
        x = frac_df.loc[common_spots, cols].sum(axis=1).to_numpy(dtype=float)

        # build signature genes from scRNA
        mask = sc_meta["_type"].map(normalize_type_name).isin(type_norm)
        if mask.sum() == 0:
            results.append({"compartment": name, "status": "SKIPPED", "reason": "no sc cells"})
            continue
        sc_in = sc_expr.loc[sc_meta.loc[mask, "cell_id"].astype(str)]
        sc_out = sc_expr.loc[sc_meta.loc[~mask, "cell_id"].astype(str)]
        mean_in = sc_in.mean(axis=0)
        mean_out = sc_out.mean(axis=0)
        score = np.log2((mean_in + 1e-8) / (mean_out + 1e-8))

        candidates = []
        for gene, val in score.items():
            if gene not in detect_map:
                continue
            if detect_map[gene] < detect_min:
                continue
            if val <= 0:
                continue
            candidates.append((gene, val))
        candidates.sort(key=lambda x: x[1], reverse=True)
        genes = [g for g, _ in candidates[:n_genes]]
        if not genes:
            results.append({"compartment": name, "status": "SKIPPED", "reason": "no genes after filter"})
            continue

        y = st_expr.loc[common_spots, genes].mean(axis=1).to_numpy(dtype=float)

        # smooth
        x_s = _smooth_values(x, coords, k, alpha)
        y_s = _smooth_values(y, coords, k, alpha)

        rho = _spearman_corr_series(x_s, y_s)
        # effect
        x_series = pd.Series(x_s)
        top_thr = float(x_series.quantile(0.8))
        bottom_thr = float(x_series.quantile(0.2))
        top_mask = x_s >= top_thr
        bottom_mask = x_s <= bottom_thr
        if top_mask.sum() == 0 or bottom_mask.sum() == 0:
            effect = float("nan")
        else:
            effect = float(np.mean(y_s[top_mask]) - np.mean(y_s[bottom_mask]))

        strong = bool(np.isfinite(rho) and abs(rho) >= 0.10) or bool(np.isfinite(effect) and abs(effect) >= 0.15)
        if strong:
            strong_count += 1
            if rho > 0 and effect > 0:
                strong_ok += 1
        else:
            weak_count += 1

        results.append({
            "compartment": name,
            "rho_spearman": float(rho) if np.isfinite(rho) else None,
            "effect_size": float(effect) if np.isfinite(effect) else None,
            "n_genes": int(len(genes)),
            "strong": bool(strong),
        })

    status = "INCONCLUSIVE"
    reasons = []
    if strong_count >= 2:
        ok_fraction = strong_ok / strong_count if strong_count else 0.0
        if ok_fraction < 0.5:
            status = "FAIL"
            reasons.append("compartment direction negative")
        elif ok_fraction < 0.6:
            status = "WARN"
            reasons.append("compartment direction below threshold")
        else:
            status = "OK"
    else:
        if weak_count == len(results):
            status = "INCONCLUSIVE"
            reasons.append("weak_compartment_signal")
        else:
            status = "WARN"
            reasons.append("too_few_strong_compartments")

    metrics = {
        "results": results,
        "summary": {
            "n_compartments": int(len(results)),
            "n_strong": int(strong_count),
            "n_strong_ok": int(strong_ok),
        },
    }

    return {"status": status, "metrics": metrics, "reasons": reasons}





def pick_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for name in candidates:
        if name in df.columns:
            return name
    return None


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Stage6 real-data audit (no truth)")
    p.add_argument("--sample", default="real_brca", help="sample id")
    p.add_argument("--baseline_id", default="baseline", help="baseline run id")
    p.add_argument("--route2_id", default="route2_v5_3_1", help="route2 run id")

    p.add_argument(
        "--fraction_baseline",
        default="result/real_brca/stage4_cytospace/runs/baseline/cell_type_assignments_by_spot.csv",
        help="baseline cell_type_assignments_by_spot.csv",
    )
    p.add_argument(
        "--assignment_baseline",
        default="result/real_brca/stage4_cytospace/runs/baseline/cell_assignment.csv",
        help="baseline cell_assignment.csv",
    )
    p.add_argument(
        "--meta_baseline",
        default="result/real_brca/stage4_cytospace/runs/baseline/stage4_summary.json",
        help="baseline stage4_summary.json",
    )

    p.add_argument(
        "--fraction_route2",
        default="result/real_brca/stage4_cytospace/runs/route2_v5_3_1/cell_type_assignments_by_spot.csv",
        help="route2 cell_type_assignments_by_spot.csv",
    )
    p.add_argument(
        "--assignment_route2",
        default="result/real_brca/stage4_cytospace/runs/route2_v5_3_1/cell_assignment.csv",
        help="route2 cell_assignment.csv",
    )
    p.add_argument(
        "--meta_route2",
        default="result/real_brca/stage4_cytospace/runs/route2_v5_3_1/stage4_summary.json",
        help="route2 stage4_summary.json",
    )

    p.add_argument(
        "--stage3_summary",
        default="result/real_brca/stage3_typematch/stage3_summary.json",
        help="route2 stage3_summary.json",
    )
    p.add_argument(
        "--sc_meta",
        default="data/processed/real_brca/stage1_preprocess/exported/sc_metadata.csv",
        help="override sc_metadata.csv for rescue ledger",
    )

    p.add_argument(
        "--st_expr",
        default="data/processed/real_brca/stage1_preprocess/exported/st_expression_normalized.csv",
        help="spot expression matrix (csv)",
    )
    p.add_argument(
        "--st_coords",
        default="data/processed/real_brca/stage1_preprocess/exported/st_coordinates.csv",
        help="spot coordinates csv",
    )
    p.add_argument(
        "--marker_config",
        default="configs/datasets/real_brca.yaml",
        help="marker config yaml (can be dataset config with stage6 section)",
    )

    p.add_argument("--out_dir", default=None, help="output directory")

    p.add_argument(
        "--assignment_mode",
        choices=["cell", "expanded"],
        default="expanded",
        help="cell_assignment semantics: cell=unique real cells; expanded=pseudo cells per spot",
    )
    p.add_argument(
        "--marker_fraction_mode",
        choices=["soft", "hard"],
        default="soft",
        help="use soft fraction (fractional_abundances_by_spot) for marker check when available",
    )

    p.add_argument("--row_sum_tol", type=float, default=1e-6)
    p.add_argument("--nonneg_eps", type=float, default=1e-9)
    p.add_argument("--active_eps", type=float, default=1e-6)
    p.add_argument("--entropy_eps", type=float, default=1e-12)

    p.add_argument("--collapse_top1_median_fail", type=float, default=0.98)
    p.add_argument("--collapse_top1_median_warn", type=float, default=0.95)
    p.add_argument("--entropy_median_fail", type=float, default=0.01)
    p.add_argument("--entropy_median_warn", type=float, default=0.05)

    p.add_argument("--marker_min_spots", type=int, default=200)
    p.add_argument("--marker_min_var", type=float, default=1e-12)
    p.add_argument("--marker_top_quantile", type=float, default=0.2)
    p.add_argument("--marker_bottom_quantile", type=float, default=0.2)
    p.add_argument("--marker_inconclusive_rho", type=float, default=0.05)
    p.add_argument("--marker_inconclusive_delta", type=float, default=0.01)

    p.add_argument("--strict", action="store_true", help="treat missing optional inputs as failure")

    return p.parse_args()


def sha1_file(path: Path) -> str:
    h = hashlib.sha1()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def clean_spot_id(value: Any) -> str:
    return str(value).split()[0].split("\t")[0]


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def resolve_path(path_value: str | None, root: Path) -> Path | None:
    if not path_value:
        return None
    p = Path(path_value)
    return p if p.is_absolute() else root / p


def extract_type_cols(path: Path) -> list[str]:
    header = pd.read_csv(path, nrows=0)
    columns = list(header.columns)
    spot_col = None
    for cand in ["spot_id", "SpotID", "spot", "Spot"]:
        if cand in columns:
            spot_col = cand
            break
    if spot_col is None and columns:
        spot_col = columns[0]
    total_cols = [c for c in columns if c.lower().replace(" ", "_") in {"total_cells", "total_cell", "total"}]
    total_col = total_cols[0] if total_cols else None
    type_cols = [c for c in columns if c not in {spot_col, total_col}]
    return type_cols


def find_soft_fraction_path(fraction_path: Path) -> Path | None:
    candidate = fraction_path.parent / "fractional_abundances_by_spot.csv"
    if candidate.exists():
        return candidate
    return None


def load_fraction(path: Path, type_list: list[str] | None = None) -> FractionData:
    df = pd.read_csv(path)
    if "spot_id" not in df.columns:
        if "SpotID" in df.columns:
            df = df.rename(columns={"SpotID": "spot_id"})
        else:
            df = df.rename(columns={df.columns[0]: "spot_id"})

    df["spot_id"] = df["spot_id"].astype(str).map(clean_spot_id)

    total_cols = [c for c in df.columns if c.lower().replace(" ", "_") in {"total_cells", "total_cell", "total"}]
    total_col = total_cols[0] if total_cols else None

    original_type_cols = [c for c in df.columns if c != "spot_id" and c != total_col]
    if not original_type_cols:
        raise ValueError(f"no type columns found in {path}")

    filled_type_cols: list[str] = []
    if type_list:
        for col in type_list:
            if col not in df.columns:
                df[col] = 0.0
                filled_type_cols.append(col)
        type_cols = list(type_list)
    else:
        type_cols = list(original_type_cols)

    df[type_cols] = df[type_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    counts = df[type_cols].to_numpy(dtype=float)

    counts_total = None
    if total_col is not None:
        counts_total = pd.to_numeric(df[total_col], errors="coerce").fillna(0.0).to_numpy(dtype=float)

    row_sum_counts = np.nansum(counts, axis=1)
    max_val = float(np.nanmax(counts)) if counts.size else 0.0
    median_sum = float(np.nanmedian(row_sum_counts)) if row_sum_counts.size else 0.0
    is_counts = (max_val > 1.0 + 1e-6) or (median_sum > 1.0 + 1e-6)

    if is_counts:
        denom = np.where(row_sum_counts == 0, 1.0, row_sum_counts)
        fractions = counts / denom[:, None]
        fraction_source = "counts"
    else:
        fractions = counts
        fraction_source = "fractions"

    row_sum_fractions = np.nansum(fractions, axis=1)

    return FractionData(
        df=df[["spot_id"] + type_cols].copy(),
        spot_ids=df["spot_id"],
        type_cols=type_cols,
        original_type_cols=original_type_cols,
        filled_type_cols=filled_type_cols,
        counts=counts,
        counts_total=counts_total,
        fractions=fractions,
        row_sum_counts=row_sum_counts,
        row_sum_fractions=row_sum_fractions,
        fraction_source=fraction_source,
    )


def load_assignment(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "spot_id" not in df.columns:
        if "assigned_spot" in df.columns:
            df = df.rename(columns={"assigned_spot": "spot_id"})
        elif "SpotID" in df.columns:
            df = df.rename(columns={"SpotID": "spot_id"})

    if "cell_id" not in df.columns:
        if "OriginalCID" in df.columns:
            df = df.rename(columns={"OriginalCID": "cell_id"})

    if "spot_id" in df.columns:
        df["spot_id"] = df["spot_id"].astype(str).map(clean_spot_id)
    if "cell_id" in df.columns:
        df["cell_id"] = df["cell_id"].astype(str)

    return df


def load_coords(path: Path | None) -> pd.DataFrame | None:
    if path is None or not path.exists():
        return None
    df = pd.read_csv(path)
    if len(df.columns) == 1 and "\t" in df.columns[0]:
        df = pd.read_csv(path, sep="\t")
    spot_col = pick_column(df, ["spot_id", "SpotID", "spot", "Spot"])
    if spot_col is None:
        return df
    df[spot_col] = df[spot_col].astype(str).map(clean_spot_id)
    return df


def load_sc_meta(path: Path | None, sample: str, project_root: Path) -> pd.DataFrame | None:
    if path is None:
        candidate = project_root / "data" / "processed" / sample / "stage1_preprocess" / "exported" / "sc_metadata.csv"
    else:
        candidate = path
    if candidate is None or not candidate.exists():
        return None
    df = pd.read_csv(candidate)
    if "cell_id" in df.columns:
        df["cell_id"] = df["cell_id"].astype(str)
    return df

def load_marker_config(path: Path | None) -> tuple[dict[str, Any] | None, str | None]:
    if path is None:
        return None, "marker_config not provided"
    if yaml is None:
        return None, "pyyaml not available"
    try:
        with path.open("r", encoding="utf-8") as f:
            data = yaml.safe_load(f) or {}
    except Exception as exc:
        return None, f"failed to load marker_config: {exc}"
    if not isinstance(data, dict):
        return None, "marker_config invalid"
    if "stage6" in data and isinstance(data["stage6"], dict):
        data = data["stage6"]
    return data, None


def _collect_expr_candidates(
    args: argparse.Namespace,
    marker_config: dict[str, Any] | None,
    project_root: Path,
    primary_path: Path | None,
) -> list[tuple[str, Path]]:
    candidates: list[tuple[str, Path]] = []

    def _add(path: Path | None):
        if path is None:
            return
        if not path.exists():
            return
        for _, p in candidates:
            if p.resolve() == path.resolve():
                return
        label = path.stem
        candidates.append((label, path))

    if primary_path is not None:
        _add(primary_path)

    if marker_config and isinstance(marker_config, dict):
        expr_source = marker_config.get("expr_source", {})
        expr_raw = expr_source.get("path") if isinstance(expr_source, dict) else None
        if expr_raw:
            expr_candidate = Path(expr_raw)
            expr_candidate = expr_candidate if expr_candidate.is_absolute() else project_root / expr_candidate
            _add(expr_candidate)
        expr_candidates = marker_config.get("expr_candidates")
        if isinstance(expr_candidates, list):
            for raw in expr_candidates:
                if not raw:
                    continue
                cand = Path(raw)
                cand = cand if cand.is_absolute() else project_root / cand
                _add(cand)

    if primary_path is not None:
        folder = primary_path.parent
        for name in ["st_expression_lognorm.csv", "st_expression_log1p.csv", "st_expression_counts.csv", "st_expression.csv"]:
            cand = folder / name
            _add(cand)

    return candidates


def _spearman_corr(x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 2:
        return float("nan")
    rx = pd.Series(x).rank(method="average").to_numpy()
    ry = pd.Series(y).rank(method="average").to_numpy()
    if np.std(rx) == 0 or np.std(ry) == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def _read_expression(path: Path, markers: list[str]) -> tuple[pd.DataFrame | None, list[str]]:
    header = pd.read_csv(path, nrows=0)
    columns = list(header.columns)
    spot_col = None
    for cand in ["spot_id", "SpotID", "spot", "Spot"]:
        if cand in columns:
            spot_col = cand
            break
    if spot_col is None:
        if columns and columns[0] == "Unnamed: 0":
            available = [m for m in markers if m in columns]
            if available:
                spot_col = columns[0]
            else:
                spot_col = None
        if spot_col is None and columns and columns[0] in {"gene", "Gene", "gene_id", "Unnamed: 0"}:
            keep = []
            missing = []
            marker_set = set(markers)
            for chunk in pd.read_csv(path, chunksize=5000, index_col=0):
                hits = [m for m in chunk.index if m in marker_set]
                if hits:
                    keep.append(chunk.loc[hits])
            if not keep:
                return None, markers
            mat = pd.concat(keep, axis=0)
            mat.index = mat.index.astype(str)
            missing = [m for m in markers if m not in mat.index]
            mat = mat.loc[~mat.index.duplicated(keep="first")]
            out = mat.T
            out.insert(0, "spot_id", out.index.astype(str))
            return out.reset_index(drop=True), missing
        if spot_col is None:
            return None, markers

    available = [m for m in markers if m in columns]
    missing = [m for m in markers if m not in columns]
    usecols = [spot_col] + available
    df = pd.read_csv(path, usecols=usecols)
    if spot_col != "spot_id":
        df.rename(columns={spot_col: "spot_id"}, inplace=True)
    return df, missing


def run_marker_checks(
    fraction_spot_ids: pd.Series,
    fraction_df: pd.DataFrame,
    marker_config: dict[str, Any] | None,
    st_expr_path: Path | None,
    min_spots: int,
    min_var: float,
    top_quantile: float,
    bottom_quantile: float,
    inconclusive_rho: float,
    inconclusive_delta: float,
) -> dict[str, Any]:
    reasons: list[str] = []
    if marker_config is None:
        return {"status": "SKIPPED", "metrics": {}, "reasons": ["marker_config missing"]}

    expr_path = st_expr_path
    if expr_path is None:
        expr_source = marker_config.get("expr_source", {}) if isinstance(marker_config, dict) else {}
        expr_raw = expr_source.get("path") if isinstance(expr_source, dict) else None
        if expr_raw:
            expr_candidate = Path(expr_raw)
            expr_path = expr_candidate if expr_candidate.is_absolute() else _ROOT / expr_candidate

    if not expr_path:
        return {"status": "SKIPPED", "metrics": {}, "reasons": ["st_expression missing"]}
    if not expr_path.exists():
        return {"status": "SKIPPED", "metrics": {}, "reasons": ["st_expression not found"]}

    marker_sets = marker_config.get("marker_sets", []) if isinstance(marker_config, dict) else []
    composites = marker_config.get("composites", []) if isinstance(marker_config, dict) else []

    def _normalize_markers(value: Any) -> list[str]:
        if not isinstance(value, list):
            return []
        return [m for m in value if isinstance(m, str)]

    entries = []
    for entry in marker_sets:
        if not isinstance(entry, dict):
            continue
        type_name = entry.get("type")
        markers = _normalize_markers(entry.get("markers", []))
        if not type_name or not markers:
            continue
        entries.append({
            "type": type_name,
            "markers": markers,
            "sum_types": None,
            "source": "direct",
        })

    for entry in composites:
        if not isinstance(entry, dict):
            continue
        name = entry.get("name")
        sum_types = _normalize_markers(entry.get("sum_types", []))
        markers = _normalize_markers(entry.get("markers", []))
        if not name or not sum_types or not markers:
            continue
        entries.append({
            "type": name,
            "markers": markers,
            "sum_types": sum_types,
            "source": "composite",
        })

    if not entries:
        return {"status": "SKIPPED", "metrics": {}, "reasons": ["marker_sets empty"]}

    all_markers: list[str] = []
    for entry in entries:
        all_markers.extend(entry["markers"])
    all_markers = sorted(set(all_markers))
    if not all_markers:
        return {"status": "SKIPPED", "metrics": {}, "reasons": ["no markers defined"]}

    expr_df, missing_markers = _read_expression(expr_path, all_markers)
    if expr_df is None or expr_df.empty:
        return {"status": "SKIPPED", "metrics": {}, "reasons": ["failed to load st_expression"]}

    expr_df["spot_id"] = expr_df["spot_id"].astype(str)
    fraction_spot_ids = fraction_spot_ids.astype(str)

    common_spots = set(fraction_spot_ids).intersection(set(expr_df["spot_id"]))

    expr_df = expr_df[expr_df["spot_id"].isin(common_spots)].copy()
    fraction_df = fraction_df[fraction_df["spot_id"].isin(common_spots)].copy()

    fraction_df.set_index("spot_id", inplace=True)
    expr_df.set_index("spot_id", inplace=True)

    test_cfg = marker_config.get("test", {}) if isinstance(marker_config, dict) else {}
    filter_cfg = marker_config.get("filters", {}) if isinstance(marker_config, dict) else {}
    strength_cfg = marker_config.get("strength_gate", {}) if isinstance(marker_config, dict) else {}
    direction_cfg = marker_config.get("direction_gate", {}) if isinstance(marker_config, dict) else {}

    def _as_float(value: Any, default: float) -> float:
        try:
            return float(value)
        except Exception:
            return default

    def _as_int(value: Any, default: int) -> int:
        try:
            return int(value)
        except Exception:
            return default

    min_spots = _as_int(test_cfg.get("min_spots", min_spots), min_spots)
    top_quantile = _as_float(test_cfg.get("top_quantile", top_quantile), top_quantile)
    bottom_quantile = _as_float(test_cfg.get("bottom_quantile", bottom_quantile), bottom_quantile)

    if len(common_spots) < min_spots:
        return {"status": "SKIPPED", "metrics": {}, "reasons": ["not enough common spots"]}

    x_var_min = _as_float(filter_cfg.get("x_var_min", 1.0e-12), 1.0e-12)
    y_var_min = _as_float(filter_cfg.get("y_var_min", min_var), min_var)
    detect_rate_min = _as_float(filter_cfg.get("detect_rate_min", 0.05), 0.05)
    inconclusive_rho = _as_float(filter_cfg.get("inconclusive_rho_abs", inconclusive_rho), inconclusive_rho)
    inconclusive_delta = _as_float(filter_cfg.get("inconclusive_delta_abs", inconclusive_delta), inconclusive_delta)

    strength_rho_min = _as_float(strength_cfg.get("strength_rho_min", 0.10), 0.10)
    strength_effect_min = _as_float(strength_cfg.get("strength_effect_min", 0.15), 0.15)
    min_strong_tests = _as_int(strength_cfg.get("min_strong_tests", 6), 6)
    strong_rho_min = _as_float(strength_cfg.get("strong_test_rho_min", strength_rho_min), strength_rho_min)
    strong_effect_min = _as_float(strength_cfg.get("strong_test_effect_min", strength_effect_min), strength_effect_min)

    ok_fraction_pass = _as_float(direction_cfg.get("ok_fraction_pass", 0.60), 0.60)
    ok_fraction_fail = _as_float(direction_cfg.get("ok_fraction_fail", 0.50), 0.50)

    norm_to_col = {normalize_type_name(c): c for c in fraction_df.columns}

    def _finite_or_none(value: Any) -> float | None:
        try:
            val = float(value)
        except Exception:
            return None
        if not np.isfinite(val):
            return None
        return val

    def _cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
        if a.size == 0 or b.size == 0:
            return float("nan")
        gt = 0
        lt = 0
        for av in a:
            gt += int(np.sum(av > b))
            lt += int(np.sum(av < b))
        n = a.size * b.size
        if n == 0:
            return float("nan")
        return float((gt - lt) / n)

    def _classify_x(x: np.ndarray) -> tuple[bool, str, dict[str, Any]]:
        meta = {
            "type_mass": float(np.sum(x)) if x.size else 0.0,
            "x_var": float(np.var(x)) if x.size else 0.0,
            "x_unique": int(np.unique(x).size) if x.size else 0,
        }
        if x.size == 0:
            return True, "empty_x", meta
        if not np.isfinite(x).all():
            return True, "x_nonfinite", meta
        if meta["x_unique"] < 2 or meta["x_var"] < x_var_min:
            return True, "x_constant", meta
        return False, "", meta

    def _classify_y(y: np.ndarray) -> tuple[bool, str, dict[str, Any]]:
        meta = {
            "detect_rate": float(np.mean(y > 0)) if y.size else 0.0,
            "y_var": float(np.var(y)) if y.size else 0.0,
        }
        if y.size == 0:
            return True, "empty_y", meta
        if not np.isfinite(y).all():
            return True, "y_nonfinite", meta
        if meta["detect_rate"] < detect_rate_min:
            return True, f"low_detect_rate({meta['detect_rate']:.3f})", meta
        if meta["y_var"] < y_var_min:
            return True, "y_constant", meta
        return False, "", meta

    tests: list[dict[str, Any]] = []
    agg_corr_top5: list[dict[str, Any]] = []
    strength_summary: list[dict[str, Any]] = []

    valid_count = 0
    inconclusive_count = 0
    skipped_count = 0
    strong_count = 0
    strong_ok = 0
    weak_global_count = 0

    type_arrays = {col: fraction_df[col].to_numpy(dtype=float) for col in fraction_df.columns}

    def _append_test(base: dict[str, Any]) -> None:
        nonlocal valid_count, inconclusive_count, skipped_count, strong_count, strong_ok
        status = base.get("status")
        if status == "OK":
            valid_count += 1
            if base.get("strong_test"):
                strong_count += 1
                if base.get("direction_ok"):
                    strong_ok += 1
        elif status == "INCONCLUSIVE":
            inconclusive_count += 1
        elif status == "SKIPPED":
            skipped_count += 1
        tests.append(base)

    def _compute_agg_stats(y_agg: np.ndarray) -> tuple[list[dict[str, Any]], float, float]:
        corr_list: list[dict[str, Any]] = []
        best_abs_rho = 0.0
        best_abs_effect = 0.0
        for col, x in type_arrays.items():
            if x.size < min_spots:
                continue
            rho = _spearman_corr(x, y_agg)
            if np.isfinite(rho):
                corr_list.append({"type": col, "rho_spearman": float(rho)})
                best_abs_rho = max(best_abs_rho, abs(float(rho)))
            x_inconclusive, _, _ = _classify_x(x)
            if x_inconclusive:
                continue
            x_series = pd.Series(x)
            top_thr = float(x_series.quantile(1.0 - top_quantile))
            bottom_thr = float(x_series.quantile(bottom_quantile))
            top_mask = x >= top_thr
            bottom_mask = x <= bottom_thr
            if top_mask.sum() == 0 or bottom_mask.sum() == 0:
                continue
            eff = _cliffs_delta(y_agg[top_mask], y_agg[bottom_mask])
            if np.isfinite(eff):
                best_abs_effect = max(best_abs_effect, abs(float(eff)))
        corr_list.sort(key=lambda d: d["rho_spearman"], reverse=True)
        return corr_list, best_abs_rho, best_abs_effect

    for entry in entries:
        type_name = entry["type"]
        markers = entry["markers"]
        source = entry["source"]
        sum_types = entry["sum_types"]

        available = [m for m in markers if m in expr_df.columns]
        agg_status = "OK"
        agg_reason = None
        best_abs_rho = 0.0
        best_abs_effect = 0.0
        corr_list: list[dict[str, Any]] = []
        if not available:
            agg_status = "SKIPPED"
            agg_reason = "marker not found"
        else:
            expr_vals = expr_df[available].to_numpy(dtype=float)
            means = np.mean(expr_vals, axis=0)
            stds = np.std(expr_vals, axis=0)
            valid_mask = stds > 0
            if not np.any(valid_mask):
                agg_status = "INCONCLUSIVE"
                agg_reason = "marker_constant"
            else:
                z = (expr_vals[:, valid_mask] - means[valid_mask]) / stds[valid_mask]
                y_agg = np.mean(z, axis=1)
                corr_list, best_abs_rho, best_abs_effect = _compute_agg_stats(y_agg)

        weak_global_signal = False
        if agg_status == "OK":
            if best_abs_rho < strength_rho_min and best_abs_effect < strength_effect_min:
                weak_global_signal = True
                agg_reason = "weak_global_signal"
        if weak_global_signal:
            weak_global_count += 1

        agg_corr_top5.append({
            "marker_set": type_name,
            "type_source": source,
            "markers_used": available,
            "n_markers_used": int(len(available)),
            "status": agg_status,
            "reason": agg_reason,
            "top5": corr_list[:5],
        })

        strength_summary.append({
            "marker_set": type_name,
            "type_source": source,
            "best_abs_rho": float(best_abs_rho),
            "best_abs_effect": float(best_abs_effect),
            "weak_global_signal": bool(weak_global_signal or agg_status != "OK"),
            "reason": agg_reason,
        })

        # determine target x for direction tests
        x = None
        type_cols_used: list[str] = []
        if sum_types:
            for t in sum_types:
                col = norm_to_col.get(normalize_type_name(t))
                if col:
                    type_cols_used.append(col)
            if type_cols_used:
                x = fraction_df[type_cols_used].sum(axis=1).to_numpy(dtype=float)
        else:
            col = norm_to_col.get(normalize_type_name(type_name))
            if col:
                type_cols_used.append(col)
                x = fraction_df[col].to_numpy(dtype=float)

        def _base_fields() -> dict[str, Any]:
            return {
                "type": type_name,
                "type_source": source,
                "type_cols_used": type_cols_used,
                "n_spots_used": int(len(fraction_df)),
                "weak_global_signal": bool(weak_global_signal),
            }

        if x is None:
            _append_test({
                **_base_fields(),
                "marker": None,
                "status": "SKIPPED",
                "reason": "type not found",
            })
            continue

        x_inconclusive, x_reason, x_meta = _classify_x(x)
        x_series = pd.Series(x)
        top_thr = float(x_series.quantile(1.0 - top_quantile))
        bottom_thr = float(x_series.quantile(bottom_quantile))

        for marker in markers:
            if marker not in expr_df.columns:
                _append_test({
                    **_base_fields(),
                    "marker": marker,
                    "status": "SKIPPED",
                    "reason": "marker not found",
                })
                continue

            if weak_global_signal:
                _append_test({
                    **_base_fields(),
                    "marker": marker,
                    "status": "INCONCLUSIVE",
                    "reason": "weak_global_signal",
                })
                continue

            if x_inconclusive:
                _append_test({
                    **_base_fields(),
                    "marker": marker,
                    "status": "INCONCLUSIVE",
                    "reason": x_reason,
                    "type_mass": _finite_or_none(x_meta.get("type_mass")),
                    "x_var": _finite_or_none(x_meta.get("x_var")),
                    "x_unique": x_meta.get("x_unique"),
                })
                continue

            y = expr_df[marker].to_numpy(dtype=float)
            y_inconclusive, y_reason, y_meta = _classify_y(y)
            base = {
                **_base_fields(),
                "marker": marker,
                "type_mass": _finite_or_none(x_meta.get("type_mass")),
                "x_var": _finite_or_none(x_meta.get("x_var")),
                "x_unique": x_meta.get("x_unique"),
                "detect_rate": _finite_or_none(y_meta.get("detect_rate")),
                "y_var": _finite_or_none(y_meta.get("y_var")),
                "top_threshold": top_thr,
                "bottom_threshold": bottom_thr,
            }

            if y_inconclusive:
                _append_test({
                    **base,
                    "status": "INCONCLUSIVE",
                    "reason": y_reason,
                })
                continue

            rho = _spearman_corr(x, y)
            top_mask = x >= top_thr
            bottom_mask = x <= bottom_thr
            if top_mask.sum() == 0 or bottom_mask.sum() == 0:
                delta_mean = float("nan")
                effect = float("nan")
            else:
                delta_mean = float(np.mean(y[top_mask]) - np.mean(y[bottom_mask]))
                effect = _cliffs_delta(y[top_mask], y[bottom_mask])

            if not np.isfinite(rho) or not np.isfinite(delta_mean):
                _append_test({
                    **base,
                    "status": "INCONCLUSIVE",
                    "reason": "nonfinite_stats",
                    "rho_spearman": _finite_or_none(rho),
                    "delta_mean_top_bottom": _finite_or_none(delta_mean),
                    "effect_size": _finite_or_none(effect),
                })
                continue

            if abs(rho) < inconclusive_rho and abs(delta_mean) < inconclusive_delta:
                _append_test({
                    **base,
                    "status": "INCONCLUSIVE",
                    "reason": "weak_signal",
                    "rho_spearman": _finite_or_none(rho),
                    "delta_mean_top_bottom": _finite_or_none(delta_mean),
                    "effect_size": _finite_or_none(effect),
                })
                continue

            ok = bool((rho > 0) and (delta_mean > 0))
            strong_test = bool((abs(rho) >= strong_rho_min) or (abs(effect) >= strong_effect_min))
            _append_test({
                **base,
                "status": "OK",
                "rho_spearman": _finite_or_none(rho),
                "delta_mean_top_bottom": _finite_or_none(delta_mean),
                "effect_size": _finite_or_none(effect),
                "direction_ok": ok,
                "strong_test": strong_test,
            })

    status = "WARN"
    direction_ok_fraction = float(strong_ok / strong_count) if strong_count else 0.0

    if strong_count >= min_strong_tests:
        if direction_ok_fraction < ok_fraction_fail:
            status = "FAIL"
            reasons.append("marker direction majority negative")
        elif direction_ok_fraction < ok_fraction_pass:
            status = "WARN"
            reasons.append("marker direction fraction below threshold")
        else:
            status = "OK"
    else:
        if weak_global_count >= len(entries):
            status = "INCONCLUSIVE"
            reasons.append("weak_global_signal")
        else:
            status = "WARN"
            reasons.append("too_few_strong_tests")

    metrics = {
        "tests": tests,
        "summary": {
            "n_tests_total": int(len(tests)),
            "n_tests_valid": int(valid_count),
            "n_tests_inconclusive": int(inconclusive_count),
            "n_tests_skipped": int(skipped_count),
            "n_strong_tests": int(strong_count),
            "n_strong_ok": int(strong_ok),
            "direction_ok_fraction": direction_ok_fraction,
            "min_strong_tests": int(min_strong_tests),
            "n_weak_global_signal": int(weak_global_count),
        },
        "strength_summary": strength_summary,
        "agg_marker_type_corr_top5": agg_corr_top5,
        "missing_markers": missing_markers,
        "spot_overlap_stats": stats_dict([len(common_spots)]),
        "config_effective": {
            "min_spots": int(min_spots),
            "top_quantile": float(top_quantile),
            "bottom_quantile": float(bottom_quantile),
            "x_var_min": float(x_var_min),
            "y_var_min": float(y_var_min),
            "detect_rate_min": float(detect_rate_min),
            "strength_rho_min": float(strength_rho_min),
            "strength_effect_min": float(strength_effect_min),
            "ok_fraction_pass": float(ok_fraction_pass),
            "ok_fraction_fail": float(ok_fraction_fail),
        },
        "mode": "strength-gated",
    }

    return {
        "status": status,
        "metrics": metrics,
        "reasons": reasons,
    }


def _count_outside_tol(values: np.ndarray, target: float, tol: float) -> int:
    return int(np.sum(np.abs(values - target) > tol))


def run_io_checks(
    fraction_data: FractionData,
    assignment_df: pd.DataFrame | None,
    coords_df: pd.DataFrame | None,
    row_sum_tol: float,
    nonneg_eps: float,
    assignment_mode: str,
    cells_per_spot_expected: int | None,
    strict: bool,
) -> dict[str, Any]:
    reasons: list[str] = []
    status = "OK"

    spot_ids = fraction_data.spot_ids.astype(str)
    type_cols = fraction_data.type_cols
    fractions = fraction_data.fractions
    row_sum = fraction_data.row_sum_fractions

    n_spots = int(len(spot_ids))
    n_types = int(len(type_cols))

    if fraction_data.filled_type_cols:
        status = "FAIL"
        reasons.append("spot_type_fraction missing type columns (filled with 0)")

    fraction_nonneg_ok = bool(np.all(fractions >= -nonneg_eps))
    if not fraction_nonneg_ok:
        status = "FAIL"
        reasons.append("fraction contains negative values")

    has_nan = int(np.isnan(fractions).sum())
    has_inf = int(np.isinf(fractions).sum())
    if has_nan > 0 or has_inf > 0:
        status = "FAIL"
        reasons.append("fraction contains NaN/Inf")

    n_outside_tol = _count_outside_tol(row_sum, 1.0, row_sum_tol)
    row_sum_stats = stats_dict(row_sum)
    if n_outside_tol > 0:
        if n_outside_tol > max(1, int(0.005 * n_spots)):
            status = "FAIL"
            reasons.append("row_sum outside tolerance for many spots")
        else:
            if status != "FAIL":
                status = "WARN"
            reasons.append("row_sum outside tolerance for some spots")

    counts_mismatch = None
    if fraction_data.counts_total is not None:
        counts_mismatch = int(np.sum(fraction_data.row_sum_counts != fraction_data.counts_total))
        if counts_mismatch > 0:
            if status != "FAIL":
                status = "WARN"
            reasons.append("row_sum counts mismatch with Total cells")

    assignment_rows = 0
    cell_id_unique_ok = None
    assignment_spot_id_subset_ok = None
    cells_per_spot_stats = None
    assignment_rank_unique_ok = None
    spot_counts_match_fraction = None
    spot_counts_mismatch_count = None
    total_expected_cells = None
    total_cells_match = None

    if assignment_df is not None and not assignment_df.empty:
        assignment_rows = int(len(assignment_df))
        if "spot_id" in assignment_df.columns:
            assignment_spots = assignment_df["spot_id"].astype(str)
            assignment_spot_id_subset_ok = bool(set(assignment_spots).issubset(set(spot_ids)))
            if not assignment_spot_id_subset_ok:
                status = "FAIL"
                reasons.append("assignment spot_id not subset of fraction spot_id")
            counts = assignment_spots.value_counts().reindex(spot_ids, fill_value=0).to_numpy(dtype=float)
            cells_per_spot_stats = stats_dict(counts)

            if assignment_mode == "expanded":
                expected = cells_per_spot_expected
                assignment_df = assignment_df.copy()
                assignment_df["_rank"] = assignment_df.groupby("spot_id").cumcount() + 1
                assignment_rank_unique_ok = bool(~assignment_df.duplicated(subset=["spot_id", "_rank"]).any())
                if not assignment_rank_unique_ok:
                    status = "FAIL"
                    reasons.append("(spot_id, rank) not unique")

                if expected is not None:
                    mismatch = (counts != expected)
                    spot_counts_mismatch_count = int(np.sum(mismatch))
                    spot_counts_match_fraction = float(1.0 - spot_counts_mismatch_count / max(1, n_spots))
                    if spot_counts_mismatch_count > 0:
                        status = "FAIL"
                        reasons.append("cells_per_spot mismatch with expected")
                    total_expected_cells = int(expected * n_spots)
                    total_cells_match = bool(assignment_rows == total_expected_cells)
                    if not total_cells_match:
                        status = "FAIL"
                        reasons.append("total assignment rows not equal to expected")
            else:
                if "cell_id" in assignment_df.columns:
                    cell_id_unique_ok = bool(assignment_df["cell_id"].is_unique)
                    if not cell_id_unique_ok:
                        status = "FAIL"
                        reasons.append("cell_id not unique")
        else:
            status = "FAIL"
            reasons.append("assignment missing spot_id")

    spot_id_intersection_ok = None
    if coords_df is not None and not coords_df.empty:
        coord_spot_col = None
        for cand in ["spot_id", "SpotID", "spot", "Spot"]:
            if cand in coords_df.columns:
                coord_spot_col = cand
                break
        if coord_spot_col:
            coord_spots = coords_df[coord_spot_col].astype(str)
            spot_id_intersection_ok = bool(set(spot_ids).issubset(set(coord_spots)))
            if not spot_id_intersection_ok:
                if status != "FAIL":
                    status = "WARN"
                reasons.append("fraction spot_id not fully covered by coords")

    metrics: dict[str, Any] = {
        "n_spots_fraction": n_spots,
        "n_types_fraction": n_types,
        "type_cols_original": fraction_data.original_type_cols,
        "type_cols_filled": fraction_data.filled_type_cols,
        "type_list": fraction_data.type_cols,
        "n_cells_assignment": assignment_rows,
        "assignment_mode": assignment_mode,
        "cells_per_spot_expected": cells_per_spot_expected,
        "cell_id_unique_ok": cell_id_unique_ok,
        "assignment_rank_unique_ok": assignment_rank_unique_ok,
        "spot_counts_match_fraction": spot_counts_match_fraction,
        "spot_counts_mismatch_count": spot_counts_mismatch_count,
        "total_expected_cells": total_expected_cells,
        "total_cells_match": total_cells_match,
        "spot_id_intersection_ok": spot_id_intersection_ok,
        "assignment_spot_id_subset_fraction_ok": assignment_spot_id_subset_ok,
        "fraction_nonneg_ok": fraction_nonneg_ok,
        "fraction_row_sum_stats": row_sum_stats,
        "fraction_row_sum_outside_tol": n_outside_tol,
        "fraction_has_nan_or_inf": {"nan": has_nan, "inf": has_inf},
        "cells_per_spot_stats": cells_per_spot_stats,
        "fraction_source": fraction_data.fraction_source,
        "counts_total_mismatch": counts_mismatch,
    }

    return {
        "status": status,
        "metrics": metrics,
        "reasons": reasons,
    }


def run_distribution_checks(
    fraction_data: FractionData,
    active_eps: float,
    entropy_eps: float,
    collapse_top1_median_fail: float,
    entropy_median_fail: float,
    collapse_top1_median_warn: float,
    entropy_median_warn: float,
) -> dict[str, Any]:
    reasons: list[str] = []
    status = "OK"

    fractions = fraction_data.fractions
    type_cols = fraction_data.type_cols
    spot_ids = fraction_data.spot_ids.astype(str)
    n_spots, n_types = fractions.shape

    entropy = compute_entropy(fractions, eps=entropy_eps)
    if n_types > 1:
        entropy_norm = entropy / np.log(float(n_types))
    else:
        entropy_norm = np.zeros_like(entropy)

    top1_idx = np.argmax(fractions, axis=1)
    top1_fraction = fractions[np.arange(n_spots), top1_idx]
    top1_type = [type_cols[i] for i in top1_idx]
    n_active_types = (fractions > active_eps).sum(axis=1)

    spot_entropy_stats = stats_dict(entropy_norm)
    spot_top1_stats = stats_dict(top1_fraction)
    spot_active_stats = stats_dict(n_active_types)

    top1_series = pd.Series(top1_type)
    top1_counts = top1_series.value_counts()
    top1_total = float(top1_counts.sum()) if not top1_counts.empty else 0.0
    top1_top10 = []
    for t, c in top1_counts.head(10).items():
        top1_top10.append({"type": t, "count": int(c), "fraction": float(c / top1_total) if top1_total else 0.0})
    top1_dominant_fraction = float(top1_counts.max() / top1_total) if top1_total else 0.0

    type_mass_table = []
    presence_thresholds = [0.01, 0.05, 0.1]
    for idx, t in enumerate(type_cols):
        col = fractions[:, idx]
        entry = {
            "type": t,
            "type_total_mass": float(col.sum()),
            "type_mean_mass": float(col.mean()),
            "type_max_fraction": float(col.max()) if col.size else float("nan"),
        }
        for thr in presence_thresholds:
            entry[f"type_presence_spots_{thr}"] = int(np.sum(col >= thr))
        type_mass_table.append(entry)

    median_top1 = spot_top1_stats.get("median", float("nan"))
    median_entropy = spot_entropy_stats.get("median", float("nan"))

    if median_top1 >= collapse_top1_median_fail:
        status = "FAIL"
        reasons.append("median top1_fraction indicates near-hard assignment")
    if median_entropy <= entropy_median_fail:
        status = "FAIL"
        reasons.append("median entropy_norm indicates collapse")
    if top1_dominant_fraction >= 0.95:
        status = "FAIL"
        reasons.append("top1_type dominates most spots")

    if status != "FAIL":
        if median_top1 >= collapse_top1_median_warn:
            status = "WARN"
            reasons.append("median top1_fraction is high")
        if median_entropy <= entropy_median_warn:
            status = "WARN"
            reasons.append("median entropy_norm is low")

    metrics = {
        "spot_entropy_stats": spot_entropy_stats,
        "spot_top1_fraction_stats": spot_top1_stats,
        "spot_active_type_count_stats": spot_active_stats,
        "spot_top1_type_frequency_top10": top1_top10,
        "spot_top1_type_dominant_fraction": top1_dominant_fraction,
        "type_mass_table": type_mass_table,
    }

    plotdata = pd.DataFrame({
        "spot_id": spot_ids.values,
        "top1_type": top1_type,
        "top1_fraction": top1_fraction,
        "entropy": entropy,
        "entropy_norm": entropy_norm,
        "n_active_types": n_active_types,
    })

    return {
        "status": status,
        "metrics": metrics,
        "reasons": reasons,
        "plotdata": plotdata,
    }


def _safe_get(d: dict, keys: list[str]) -> Any:
    cur: Any = d
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return None
        cur = cur[k]
    return cur


def _pick_column(df: pd.DataFrame | None, candidates: list[str]) -> str | None:
    if df is None:
        return None
    for name in candidates:
        if name in df.columns:
            return name
    return None


def run_ledger_checks(
    stage4_summary: dict | None,
    stage3_summary: dict | None,
    assignment_df: pd.DataFrame | None,
    sc_meta_df: pd.DataFrame | None,
) -> dict[str, Any]:
    if stage3_summary is None:
        return {
            "status": "SKIPPED",
            "metrics": {},
            "reasons": ["stage3_summary missing"],
        }

    reasons: list[str] = []
    status = "OK"

    expected = _safe_get(stage3_summary, ["action_overview", "by_cell_count", "Dropped"])
    observed = None
    if stage4_summary:
        observed = stage4_summary.get("n_filtered")

    ledger_check_ok = None
    if expected is not None and observed is not None:
        ledger_check_ok = bool(int(expected) == int(observed))
        if not ledger_check_ok:
            status = "FAIL"
            reasons.append("filter ledger mismatch")

    rescued_types = _safe_get(stage3_summary, ["v5_entropy_qc", "rescued_types_final"])
    if rescued_types is None:
        rescued_types = _safe_get(stage3_summary, ["v5_rescue_control", "rescued_types"]) or []
    if not isinstance(rescued_types, list):
        rescued_types = []

    rescued_types_norm = {normalize_type_name(t) for t in rescued_types}

    rescue_metrics = {
        "rescued_types": sorted(rescued_types),
        "rescued_types_count": int(len(rescued_types)),
        "rescued_cells_total": None,
        "rescued_cells_in_assignment": None,
        "rescued_cells_filtered_or_missing": None,
        "rescued_cells_marked_unknown": None,
        "by_type": [],
    }

    if rescued_types_norm and sc_meta_df is not None and assignment_df is not None:
        cell_id_col = _pick_column(sc_meta_df, ["cell_id", "cell_id.1", "CellID", "OriginalCID"])
        cell_type_col = _pick_column(sc_meta_df, ["cell_type", "CellType", "sc_type", "type"])
        assign_cell_id_col = _pick_column(assignment_df, ["cell_id", "CellID", "OriginalCID"])

        if cell_id_col and cell_type_col and assign_cell_id_col:
            sc_meta_df = sc_meta_df.copy()
            assignment_df = assignment_df.copy()
            sc_meta_df["_type_norm"] = sc_meta_df[cell_type_col].map(normalize_type_name)
            assignment_df["_cell_id"] = assignment_df[assign_cell_id_col].astype(str)
            sc_meta_df["_cell_id"] = sc_meta_df[cell_id_col].astype(str)

            rescued_cells_total = 0
            rescued_cells_in_assignment = 0
            by_type = []
            for t in sorted(rescued_types_norm):
                sc_mask = sc_meta_df["_type_norm"] == t
                total = int(sc_mask.sum())
                if total == 0:
                    by_type.append({
                        "cell_type": t,
                        "n_sc_cells_total": 0,
                        "n_in_assignment": 0,
                    })
                    continue
                ids = set(sc_meta_df.loc[sc_mask, "_cell_id"].tolist())
                in_assign = int(assignment_df["_cell_id"].isin(ids).sum())
                by_type.append({
                    "cell_type": t,
                    "n_sc_cells_total": total,
                    "n_in_assignment": in_assign,
                })
                rescued_cells_total += total
                rescued_cells_in_assignment += in_assign

            rescue_metrics.update({
                "rescued_cells_total": rescued_cells_total,
                "rescued_cells_in_assignment": rescued_cells_in_assignment,
                "rescued_cells_filtered_or_missing": rescued_cells_total - rescued_cells_in_assignment,
                "by_type": by_type,
            })
        else:
            reasons.append("sc_metadata columns missing for rescue ledger")
            if status != "FAIL":
                status = "WARN"
    else:
        if rescued_types_norm:
            reasons.append("missing sc_metadata or assignment for rescue ledger")
            if status != "FAIL":
                status = "WARN"

    metrics = {
        "filter_ledger": {
            "expected_n_filtered": expected,
            "observed_n_filtered": observed,
            "ledger_check_ok": ledger_check_ok,
        },
        "rescue_ledger": rescue_metrics,
    }

    return {
        "status": status,
        "metrics": metrics,
        "reasons": reasons,
    }

def _type_mass_dict(type_mass_table: list[dict[str, Any]]) -> dict[str, float]:
    out: dict[str, float] = {}
    for entry in type_mass_table:
        if not isinstance(entry, dict):
            continue
        t = entry.get("type")
        v = entry.get("type_total_mass")
        if isinstance(t, str) and v is not None:
            out[t] = float(v)
    return out


def run_compare(
    baseline_fraction: FractionData,
    route2_fraction: FractionData,
    baseline_dist: dict[str, Any],
    route2_dist: dict[str, Any],
    eps: float = EPS,
) -> dict[str, Any]:
    delta = {}

    base_entropy = baseline_dist.get("spot_entropy_stats", {}).get("median")
    route_entropy = route2_dist.get("spot_entropy_stats", {}).get("median")
    if base_entropy is not None and route_entropy is not None:
        delta["delta_spot_entropy_median"] = float(route_entropy - base_entropy)

    base_top1 = baseline_dist.get("spot_top1_fraction_stats", {}).get("median")
    route_top1 = route2_dist.get("spot_top1_fraction_stats", {}).get("median")
    if base_top1 is not None and route_top1 is not None:
        delta["delta_spot_top1_fraction_median"] = float(route_top1 - base_top1)

    base_mass = _type_mass_dict(baseline_dist.get("type_mass_table", []))
    route_mass = _type_mass_dict(route2_dist.get("type_mass_table", []))
    all_types = sorted(set(base_mass) | set(route_mass))

    delta_type_mass = {
        t: float(route_mass.get(t, 0.0) - base_mass.get(t, 0.0)) for t in all_types
    }

    base_total = sum(base_mass.values()) if base_mass else 0.0
    route_total = sum(route_mass.values()) if route_mass else 0.0
    base_dist = np.array([base_mass.get(t, 0.0) for t in all_types], dtype=float)
    route_dist = np.array([route_mass.get(t, 0.0) for t in all_types], dtype=float)

    global_js = None
    if base_total > 0 and route_total > 0:
        global_js = js_divergence(base_dist, route_dist, eps=eps)

    base_top1 = baseline_dist.get("spot_top1_type_frequency_top10", [])
    route_top1 = route2_dist.get("spot_top1_type_frequency_top10", [])
    base_top1_map = {row["type"]: row.get("fraction", 0.0) for row in base_top1 if isinstance(row, dict)}
    route_top1_map = {row["type"]: row.get("fraction", 0.0) for row in route_top1 if isinstance(row, dict)}
    top1_shift = {}
    for t in sorted(set(base_top1_map) | set(route_top1_map)):
        top1_shift[t] = float(route_top1_map.get(t, 0.0) - base_top1_map.get(t, 0.0))

    spot_kl_stats = None
    common_spots = set(baseline_fraction.spot_ids).intersection(set(route2_fraction.spot_ids))
    if common_spots:
        base_df = baseline_fraction.df.set_index("spot_id")
        route_df = route2_fraction.df.set_index("spot_id")
        base_df = base_df.loc[base_df.index.intersection(common_spots)]
        route_df = route_df.loc[route_df.index.intersection(common_spots)]
        base_df = base_df[baseline_fraction.type_cols]
        route_df = route_df[route2_fraction.type_cols]

        all_types = sorted(set(baseline_fraction.type_cols) | set(route2_fraction.type_cols))
        base_arr = np.zeros((len(base_df), len(all_types)), dtype=float)
        route_arr = np.zeros((len(route_df), len(all_types)), dtype=float)
        base_idx = {t: i for i, t in enumerate(all_types)}

        for t in baseline_fraction.type_cols:
            base_arr[:, base_idx[t]] = base_df[t].to_numpy(dtype=float)
        for t in route2_fraction.type_cols:
            route_arr[:, base_idx[t]] = route_df[t].to_numpy(dtype=float)

        p = np.clip(route_arr, eps, 1.0)
        q = np.clip(base_arr, eps, 1.0)
        p = p / p.sum(axis=1, keepdims=True)
        q = q / q.sum(axis=1, keepdims=True)
        kl = np.sum(p * np.log(p / q), axis=1)
        spot_kl_stats = stats_dict(kl)

    return {
        "delta": delta,
        "global_type_mass_baseline": base_mass,
        "global_type_mass_route2": route_mass,
        "delta_type_total_mass": delta_type_mass,
        "global_distribution_js": global_js,
        "spotwise_kl_stats": spot_kl_stats,
        "top1_type_shift": top1_shift,
    }


def build_plotdata(
    base_plot: pd.DataFrame,
    coords_df: pd.DataFrame | None,
    st_expr_path: Path | None,
    marker_config: dict[str, Any] | None,
) -> pd.DataFrame:
    df = base_plot.copy()
    if coords_df is not None and not coords_df.empty:
        spot_col = pick_column(coords_df, ["spot_id", "SpotID", "spot", "Spot"])
        if spot_col:
            coords_df = coords_df.rename(columns={spot_col: "spot_id"})
        df = df.merge(coords_df, on="spot_id", how="left")

    if st_expr_path is None or marker_config is None:
        return df

    marker_sets = marker_config.get("marker_sets", []) if isinstance(marker_config, dict) else []
    composites = marker_config.get("composites", []) if isinstance(marker_config, dict) else []
    markers = []
    for entry in list(marker_sets) + list(composites):
        if not isinstance(entry, dict):
            continue
        markers.extend([m for m in entry.get("markers", []) if isinstance(m, str)])
    markers = sorted(set(markers))
    if not markers:
        return df

    expr_df, _missing = _read_expression(st_expr_path, markers)
    if expr_df is None or expr_df.empty:
        return df
    expr_df["spot_id"] = expr_df["spot_id"].astype(str).map(clean_spot_id)
    df = df.merge(expr_df, on="spot_id", how="left")
    return df


def aggregate_status(parts: list[dict[str, Any]], strict: bool) -> tuple[str, list[str]]:
    reasons: list[str] = []
    statuses = [p.get("status") for p in parts if isinstance(p, dict)]
    for p in parts:
        reasons.extend(p.get("reasons", []) or [])

    if any(s == "FAIL" for s in statuses):
        return "FAIL", reasons
    if any(s in ("WARN", "INCONCLUSIVE") for s in statuses):
        return "WARN", reasons
    if any(s == "SKIPPED" for s in statuses):
        return ("FAIL" if strict else "WARN"), reasons
    return "PASS", reasons


def run_single(
    run_id: str,
    fraction_path: Path,
    assignment_path: Path,
    meta_path: Path,
    stage3_summary: dict | None,
    coords_df: pd.DataFrame | None,
    sc_meta_df: pd.DataFrame | None,
    marker_config: dict | None,
    st_expr_path: Path | None,
    args: argparse.Namespace,
    out_dir: Path,
    enable_ledger: bool,
    type_list: list[str],
) -> tuple[dict[str, Any], FractionData, dict[str, Any], pd.DataFrame]:
    fraction_data = load_fraction(fraction_path, type_list=type_list)
    assignment_df = load_assignment(assignment_path)
    stage4_summary = read_json(meta_path)

    cells_per_spot_expected = stage4_summary.get("mapping_cells_per_spot")
    if cells_per_spot_expected is not None:
        try:
            cells_per_spot_expected = int(cells_per_spot_expected)
        except Exception:
            cells_per_spot_expected = None

    io_result = run_io_checks(
        fraction_data,
        assignment_df,
        coords_df,
        row_sum_tol=args.row_sum_tol,
        nonneg_eps=args.nonneg_eps,
        assignment_mode=args.assignment_mode,
        cells_per_spot_expected=cells_per_spot_expected,
        strict=args.strict,
    )

    dist_result = run_distribution_checks(
        fraction_data,
        active_eps=args.active_eps,
        entropy_eps=args.entropy_eps,
        collapse_top1_median_fail=args.collapse_top1_median_fail,
        entropy_median_fail=args.entropy_median_fail,
        collapse_top1_median_warn=args.collapse_top1_median_warn,
        entropy_median_warn=args.entropy_median_warn,
    )

    marker_fraction = fraction_data
    marker_fraction_path = fraction_path
    if args.marker_fraction_mode == "soft":
        soft_path = find_soft_fraction_path(fraction_path)
        if soft_path is not None:
            marker_fraction = load_fraction(soft_path, type_list=type_list)
            marker_fraction_path = soft_path

    expr_candidates = _collect_expr_candidates(args, marker_config, _ROOT, st_expr_path)
    marker_choice = None
    marker_runs = []
    for label, path in expr_candidates:
        result = run_marker_checks(
            marker_fraction.spot_ids,
            marker_fraction.df,
            marker_config,
            path,
            min_spots=args.marker_min_spots,
            min_var=args.marker_min_var,
            top_quantile=args.marker_top_quantile,
            bottom_quantile=args.marker_bottom_quantile,
            inconclusive_rho=args.marker_inconclusive_rho,
            inconclusive_delta=args.marker_inconclusive_delta,
        )
        summary = result.get("metrics", {}).get("strength_summary", [])
        scores = []
        for row in summary:
            if not isinstance(row, dict):
                continue
            rho = row.get("best_abs_rho")
            eff = row.get("best_abs_effect")
            rho_val = abs(float(rho)) if rho is not None else 0.0
            eff_val = abs(float(eff)) if eff is not None else 0.0
            if rho is not None or eff is not None:
                scores.append(max(rho_val, eff_val))
        strength_score = float(np.median(scores)) if scores else 0.0
        marker_runs.append({
            "label": label,
            "path": path,
            "score": strength_score,
            "result": result,
        })

    if marker_runs:
        marker_runs.sort(key=lambda d: d["score"], reverse=True)
        marker_choice = marker_runs[0]
        marker_result = marker_choice["result"]
    else:
        marker_result = run_marker_checks(
            marker_fraction.spot_ids,
            marker_fraction.df,
            marker_config,
            st_expr_path,
            min_spots=args.marker_min_spots,
            min_var=args.marker_min_var,
            top_quantile=args.marker_top_quantile,
            bottom_quantile=args.marker_bottom_quantile,
            inconclusive_rho=args.marker_inconclusive_rho,
            inconclusive_delta=args.marker_inconclusive_delta,
        )

    marker_result.setdefault("metrics", {})
    marker_result["metrics"]["fraction_source"] = marker_fraction.fraction_source
    marker_result["metrics"]["fraction_path"] = str(marker_fraction_path)
    marker_result["metrics"]["expr_layer_used"] = marker_choice["label"] if marker_choice else None
    marker_result["metrics"]["expr_layer_candidates"] = [
        {"label": m["label"], "path": str(m["path"]), "strength_score": m["score"]}
        for m in marker_runs
    ]

    reconstruction_result = {"status": "SKIPPED", "metrics": {}, "reasons": ["missing inputs"]}
    compartment_result = {"status": "SKIPPED", "metrics": {}, "reasons": ["missing inputs"]}

    # reconstruction uses preferred expr layer if available
    recon_cfg = marker_config.get("reconstruction", {}) if isinstance(marker_config, dict) else {}
    recon_gene_set = recon_cfg.get("gene_set", "st_hvg2000")
    recon_n_genes = int(recon_cfg.get("n_genes", 2000))
    recon_null = int(recon_cfg.get("null_shuffle", 50))
    recon_corr_delta = float(recon_cfg.get("corr_delta_min", 0.05))
    recon_rmse_delta = float(recon_cfg.get("rmse_delta_min", 0.05))

    # pick reconstruction layer: prefer lognorm
    recon_candidates = [
        _ROOT / "data" / "processed" / args.sample / "stage1_preprocess" / "exported" / "st_expression_lognorm.csv",
        _ROOT / "data" / "processed" / args.sample / "stage1_preprocess" / "exported" / "st_expression_normalized.csv",
        _ROOT / "data" / "processed" / args.sample / "stage1_preprocess" / "exported" / "st_expression_counts.csv",
    ]
    recon_expr = None
    for cand in recon_candidates:
        if cand.exists():
            recon_expr = _read_expression_wide(cand)
            recon_expr["spot_id"] = recon_expr["spot_id"].astype(str)
            recon_expr_layer = cand.name
            break

    sc_expr_path = _ROOT / "data" / "processed" / args.sample / "stage1_preprocess" / "exported" / "sc_expression_normalized.csv"
    sc_expr = _load_sc_expression(sc_expr_path)

    if recon_expr is not None and sc_expr is not None and sc_meta_df is not None:
        sc_profiles = _build_type_profiles(sc_expr, sc_meta_df, "cell_type")
        reconstruction_result = run_reconstruction_checks(
            fraction_data,
            recon_expr,
            sc_profiles,
            recon_gene_set,
            recon_n_genes,
            recon_null,
            recon_corr_delta,
            recon_rmse_delta,
        )
        reconstruction_result.setdefault("metrics", {})["expr_layer_used"] = recon_expr_layer

    # compartments use same expr layer as reconstruction
    if recon_expr is not None and sc_expr is not None and sc_meta_df is not None and coords_df is not None:
        compartment_result = run_compartment_checks(
            fraction_data,
            sc_expr,
            sc_meta_df,
            recon_expr,
            coords_df,
            marker_config,
        )
        compartment_result.setdefault("metrics", {})["expr_layer_used"] = recon_expr_layer

    if enable_ledger:
        ledger_result = run_ledger_checks(stage4_summary, stage3_summary, assignment_df, sc_meta_df)
    else:
        ledger_result = {"status": "SKIPPED", "metrics": {}, "reasons": ["baseline ledger skipped"]}

    parts = [io_result, dist_result, ledger_result]
    status, reasons = aggregate_status(parts, args.strict)

    # override overall status with reconstruction/compartment logic
    reconstruction_status = reconstruction_result.get("status") if isinstance(reconstruction_result, dict) else "SKIPPED"
    compartment_status = compartment_result.get("status") if isinstance(compartment_result, dict) else "SKIPPED"

    if io_result.get("status") == "FAIL":
        status = "FAIL"
    else:
        if reconstruction_status == "OK":
            if compartment_status == "FAIL":
                status = "FAIL"
            elif compartment_status in ("WARN", "INCONCLUSIVE"):
                status = "WARN"
            else:
                status = "PASS"
        elif reconstruction_status == "WARN":
            status = "WARN"
        else:
            status = "WARN"

    audit = {
        "stage": "stage6_real_audit",
        "sample": args.sample,
        "run_id": run_id,
        "created_at": pd.Timestamp.utcnow().isoformat(),
        "input_fingerprint": {
            "spot_type_fraction_sha1": sha1_file(fraction_path),
            "cell_assignment_sha1": sha1_file(assignment_path),
            "meta_sha1": sha1_file(meta_path),
            "st_expr_sha1": sha1_file(st_expr_path) if st_expr_path else None,
            "stage3_summary_sha1": sha1_file(Path(args.stage3_summary)) if args.stage3_summary else None,
        },
        "meta_trace": {
            "module_sha1": stage4_summary.get("sha1"),
            "runner_sha1": stage4_summary.get("runner_sha1"),
            "config_effective_subset": stage4_summary.get("config_effective_subset"),
            "knn_mode": stage4_summary.get("knn_mode"),
            "harden_method": stage4_summary.get("harden_method"),
            "cells_per_spot": stage4_summary.get("cells_per_spot_override"),
            "type_list": type_list,
        },
        "io_integrity": io_result,
        "ledger_integrity": ledger_result,
        "distribution_sanity": {
            "status": dist_result.get("status"),
            "metrics": dist_result.get("metrics"),
            "reasons": dist_result.get("reasons"),
        },
        "weak_supervision_marker_check": marker_result,
        "reconstruction_consistency": reconstruction_result,
        "compartment_marker_check": compartment_result,
        "overall_status": status,
        "reasons": reasons,
    }

    run_dir = out_dir / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    debug_tests = marker_result.get("metrics", {}).get("tests") if isinstance(marker_result, dict) else None
    if debug_tests:
        debug_df = pd.DataFrame(debug_tests)
        debug_path = run_dir / f"stage6_marker_debug_{run_id}.csv"
        debug_df.to_csv(debug_path, index=False)
        marker_result.setdefault("metrics", {})["debug_path"] = str(debug_path)

    agg_top5 = marker_result.get("metrics", {}).get("agg_marker_type_corr_top5")
    if agg_top5:
        rows = []
        for entry in agg_top5:
            base = {
                "marker_set": entry.get("marker_set"),
                "type_source": entry.get("type_source"),
                "status": entry.get("status"),
                "reason": entry.get("reason"),
                "n_markers_used": entry.get("n_markers_used"),
                "markers_used": "|".join(entry.get("markers_used", []) or []),
            }
            top5 = entry.get("top5") or []
            if top5:
                for rank, item in enumerate(top5, start=1):
                    rows.append({
                        **base,
                        "rank": rank,
                        "corr_type": item.get("type"),
                        "rho_spearman": item.get("rho_spearman"),
                    })
            else:
                rows.append({
                    **base,
                    "rank": None,
                    "corr_type": None,
                    "rho_spearman": None,
                })
        agg_df = pd.DataFrame(rows)
        agg_path = run_dir / f"stage6_marker_top5_{run_id}.csv"
        agg_df.to_csv(agg_path, index=False)
        marker_result.setdefault("metrics", {})["agg_top5_path"] = str(agg_path)

    if reconstruction_result.get("metrics"):
        recon_path = run_dir / f"stage6_reconstruction_summary_{run_id}.json"
        recon_path.write_text(json.dumps(reconstruction_result, ensure_ascii=False, indent=2), encoding="utf-8")
        marker_result.setdefault("metrics", {})["reconstruction_path"] = str(recon_path)

    if compartment_result.get("metrics"):
        comp_path = run_dir / f"stage6_compartment_marker_{run_id}.csv"
        comp_df = pd.DataFrame(compartment_result.get("metrics", {}).get("results", []))
        comp_df.to_csv(comp_path, index=False)
        marker_result.setdefault("metrics", {})["compartment_path"] = str(comp_path)

    strength_summary = marker_result.get("metrics", {}).get("strength_summary")
    if strength_summary:
        strength_df = pd.DataFrame(strength_summary)
        strength_path = run_dir / f"stage6_marker_strength_summary_{run_id}.csv"
        strength_df.to_csv(strength_path, index=False)
        marker_result.setdefault("metrics", {})["strength_summary_path"] = str(strength_path)

    status_explain = {
        "marker": {
            "status": marker_result.get("status"),
            "reasons": marker_result.get("reasons"),
            "summary": marker_result.get("metrics", {}).get("summary"),
            "expr_layer_used": marker_result.get("metrics", {}).get("expr_layer_used"),
            "expr_layer_candidates": marker_result.get("metrics", {}).get("expr_layer_candidates"),
            "config_effective": marker_result.get("metrics", {}).get("config_effective"),
            "mode": marker_result.get("metrics", {}).get("mode"),
        },
        "reconstruction": reconstruction_result,
        "compartment": compartment_result,
        "overall": {
            "status": status,
            "reasons": reasons,
        },
    }
    status_path = run_dir / f"stage6_marker_status_explain_{run_id}.json"
    status_path.write_text(json.dumps(status_explain, ensure_ascii=False, indent=2), encoding="utf-8")
    marker_result.setdefault("metrics", {})["status_explain_path"] = str(status_path)

    audit_path = run_dir / f"stage6_audit_{run_id}.json"
    audit_path.write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")

    plotdata = build_plotdata(dist_result["plotdata"], coords_df, st_expr_path, marker_config)
    plot_path = run_dir / f"stage6_plotdata_{run_id}.csv"
    plotdata.to_csv(plot_path, index=False)

    return audit, fraction_data, dist_result.get("metrics", {}), plotdata


def main() -> None:
    args = parse_args()
    project_root = _ROOT

    baseline_fraction = resolve_path(args.fraction_baseline, project_root)
    baseline_assignment = resolve_path(args.assignment_baseline, project_root)
    baseline_meta = resolve_path(args.meta_baseline, project_root)
    route2_fraction = resolve_path(args.fraction_route2, project_root)
    route2_assignment = resolve_path(args.assignment_route2, project_root)
    route2_meta = resolve_path(args.meta_route2, project_root)

    stage3_summary = None
    stage3_path = resolve_path(args.stage3_summary, project_root) if args.stage3_summary else None
    if stage3_path is not None and stage3_path.exists():
        stage3_summary = read_json(stage3_path)

    coords_path = resolve_path(args.st_coords, project_root) if args.st_coords else None
    coords_df = load_coords(coords_path)

    sc_meta_path = resolve_path(args.sc_meta, project_root) if args.sc_meta else None
    sc_meta_df = load_sc_meta(sc_meta_path, args.sample, project_root)

    marker_path = resolve_path(args.marker_config, project_root) if args.marker_config else None
    marker_config, marker_err = load_marker_config(marker_path)
    if marker_err and args.strict:
        raise RuntimeError(marker_err)

    st_expr_path = resolve_path(args.st_expr, project_root) if args.st_expr else None

    if baseline_fraction is None or baseline_assignment is None or baseline_meta is None:
        raise RuntimeError("baseline inputs missing")
    if route2_fraction is None or route2_assignment is None or route2_meta is None:
        raise RuntimeError("route2 inputs missing")

    type_list = []
    try:
        type_list.extend(extract_type_cols(baseline_fraction))
    except Exception:
        pass
    try:
        type_list.extend(extract_type_cols(route2_fraction))
    except Exception:
        pass
    if stage3_summary and isinstance(stage3_summary.get("plugin_types"), list):
        type_list.extend(stage3_summary["plugin_types"])
    type_list = sorted({t for t in type_list if t})

    out_dir = Path(args.out_dir) if args.out_dir else project_root / "result" / args.sample / "stage6_real_audit"
    out_dir.mkdir(parents=True, exist_ok=True)

    _, baseline_fraction_data, baseline_dist, _ = run_single(
        args.baseline_id,
        baseline_fraction,
        baseline_assignment,
        baseline_meta,
        stage3_summary=None,
        coords_df=coords_df,
        sc_meta_df=sc_meta_df,
        marker_config=marker_config,
        st_expr_path=st_expr_path,
        args=args,
        out_dir=out_dir,
        enable_ledger=False,
        type_list=type_list,
    )

    _, route2_fraction_data, route2_dist, _ = run_single(
        args.route2_id,
        route2_fraction,
        route2_assignment,
        route2_meta,
        stage3_summary=stage3_summary,
        coords_df=coords_df,
        sc_meta_df=sc_meta_df,
        marker_config=marker_config,
        st_expr_path=st_expr_path,
        args=args,
        out_dir=out_dir,
        enable_ledger=True,
        type_list=type_list,
    )

    compare = run_compare(
        baseline_fraction_data,
        route2_fraction_data,
        baseline_dist,
        route2_dist,
    )

    compare_path = out_dir / "stage6_compare_baseline_vs_route2.json"
    compare_path.write_text(json.dumps(compare, ensure_ascii=False, indent=2), encoding="utf-8")

    baseline_plot_path = out_dir / args.baseline_id / f"stage6_plotdata_{args.baseline_id}.csv"
    route2_plot_path = out_dir / args.route2_id / f"stage6_plotdata_{args.route2_id}.csv"
    if baseline_plot_path.exists() and route2_plot_path.exists():
        baseline_plot = pd.read_csv(baseline_plot_path)
        route2_plot = pd.read_csv(route2_plot_path)
        common = set(baseline_plot["spot_id"]).intersection(set(route2_plot["spot_id"]))
        if common:
            base_plot = baseline_plot[baseline_plot["spot_id"].isin(common)].copy().set_index("spot_id")
            route_plot = route2_plot[route2_plot["spot_id"].isin(common)].copy().set_index("spot_id")
            delta_df = pd.DataFrame({"spot_id": list(common)})
            delta_df["delta_entropy_norm"] = route_plot.loc[delta_df["spot_id"], "entropy_norm"].values - base_plot.loc[delta_df["spot_id"], "entropy_norm"].values
            delta_df["delta_top1_fraction"] = route_plot.loc[delta_df["spot_id"], "top1_fraction"].values - base_plot.loc[delta_df["spot_id"], "top1_fraction"].values
            delta_df["delta_n_active_types"] = route_plot.loc[delta_df["spot_id"], "n_active_types"].values - base_plot.loc[delta_df["spot_id"], "n_active_types"].values
            compare_plot_path = out_dir / "stage6_plotdata_compare.csv"
            delta_df.to_csv(compare_plot_path, index=False)

    print(f"[Stage6] audit outputs -> {out_dir}")


if __name__ == "__main__":
    main()
