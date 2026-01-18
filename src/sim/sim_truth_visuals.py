"""
Sim truth visualization generator (ground-truth only).
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd

try:  # pragma: no cover
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # type: ignore
    from matplotlib.colors import LinearSegmentedColormap  # type: ignore
except Exception:  # pragma: no cover
    plt = None

_HERE = Path(__file__).resolve()
_ROOT = _HERE.parents[2]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src.utils.type_name import canonicalize_type_name, normalize_type_name, load_alias_map

STYLE = {
    "title": 16,
    "subplot": 12,
    "label": 10,
    "tick": 9,
}


def _maybe_plot() -> bool:
    return plt is not None


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _clean_spot_id(value: Any) -> str:
    return str(value).split()[0].split("\t")[0]


def _safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", value)
    return cleaned.strip("_") or "unknown"


def _pick_spot_col(columns: Iterable[str]) -> str:
    for cand in ["spot_id", "SpotID", "spot", "Spot"]:
        if cand in columns:
            return cand
    return next(iter(columns))


def _load_coords(path: Path | None) -> Tuple[pd.DataFrame | None, str | None, str | None]:
    if path is None or not path.exists():
        return None, None, None
    df = pd.read_csv(path)
    if len(df.columns) == 1 and "\t" in df.columns[0]:
        df = pd.read_csv(path, sep="\t")
    spot_col = _pick_spot_col(df.columns)
    if spot_col != "spot_id":
        df.rename(columns={spot_col: "spot_id"}, inplace=True)
    df["spot_id"] = df["spot_id"].astype(str).map(_clean_spot_id)

    candidates = [("x", "y"), ("X", "Y"), ("row", "col"), ("Row", "Col")]
    for x_col, y_col in candidates:
        if x_col in df.columns and y_col in df.columns:
            return df, x_col, y_col

    numeric_cols = [c for c in df.columns if c != "spot_id" and pd.api.types.is_numeric_dtype(df[c])]
    if len(numeric_cols) >= 2:
        return df, numeric_cols[0], numeric_cols[1]
    return df, None, None


def _resolve_truth_fraction_path(stage1_export: Path, override: str | None) -> Path | None:
    if override:
        path = Path(override)
        return path if path.exists() else None
    for name in [
        "sim_truth_spot_type_fraction.csv",
        "spot_type_fraction_truth.csv",
        "truth_spot_type_fraction.csv",
    ]:
        cand = stage1_export / name
        if cand.exists():
            return cand
    return None


def _resolve_stage4_run_dir(stage4_root: Path, run_id: str) -> Path | None:
    run_dir = stage4_root / "runs" / run_id
    if run_dir.exists():
        return run_dir
    fallback = stage4_root / "cytospace_output"
    if fallback.exists():
        return fallback
    return None


def _load_fraction_table(path: Path) -> Tuple[pd.DataFrame, List[str], bool]:
    df = pd.read_csv(path)
    spot_col = _pick_spot_col(df.columns)
    if spot_col != "spot_id":
        df.rename(columns={spot_col: "spot_id"}, inplace=True)
    df["spot_id"] = df["spot_id"].astype(str).map(_clean_spot_id)

    total_cols = [c for c in df.columns if c.lower().replace(" ", "_") in {"total_cells", "total_cell", "total"}]
    total_col = total_cols[0] if total_cols else None
    type_cols = [c for c in df.columns if c not in {"spot_id", total_col}]
    if not type_cols:
        return df[["spot_id"]].copy(), [], False

    df[type_cols] = df[type_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    counts = df[type_cols].to_numpy(dtype=float)
    row_sum = np.nansum(counts, axis=1)
    max_val = float(np.nanmax(counts)) if counts.size else 0.0
    median_sum = float(np.nanmedian(row_sum)) if row_sum.size else 0.0
    is_counts = (max_val > 1.0 + 1e-6) or (median_sum > 1.0 + 1e-6)

    if is_counts:
        denom = np.where(row_sum == 0, 1.0, row_sum)
        fractions = counts / denom[:, None]
    else:
        fractions = counts

    out = pd.DataFrame(fractions, columns=type_cols)
    out.insert(0, "spot_id", df["spot_id"].values)
    return out, type_cols, is_counts


def _canonicalize_fraction_df(df: pd.DataFrame, alias_map: Dict[str, str]) -> Tuple[pd.DataFrame, List[str]]:
    if "spot_id" not in df.columns:
        return df, []
    type_cols = [c for c in df.columns if c != "spot_id"]
    groups: Dict[str, List[str]] = {}
    for col in type_cols:
        canon = canonicalize_type_name(col, alias_map) if alias_map else normalize_type_name(col)
        if not canon:
            continue
        groups.setdefault(canon, []).append(col)
    out = pd.DataFrame({"spot_id": df["spot_id"].values})
    for canon, cols in groups.items():
        out[canon] = df[cols].sum(axis=1)
    return out, list(groups.keys())


def _build_category_colors(types: List[str]) -> Dict[str, Any]:
    if not _maybe_plot() or not types:
        return {}
    cmap = plt.get_cmap("tab20", len(types))
    return {t: cmap(i) for i, t in enumerate(types)}


def _rgba_to_hex(color: Any) -> str:
    if isinstance(color, str):
        return color
    r, g, b = color[:3]
    return "#{:02x}{:02x}{:02x}".format(int(r * 255), int(g * 255), int(b * 255))


def _major_type_labels(df: pd.DataFrame, type_cols: List[str], unknown_label: str) -> List[str]:
    if not type_cols:
        return [unknown_label] * len(df)
    vals = df[type_cols].to_numpy(dtype=float)
    sums = np.sum(vals, axis=1)
    idx = np.argmax(vals, axis=1)
    labels = []
    for i, j in enumerate(idx):
        labels.append(type_cols[j] if sums[i] > 0 else unknown_label)
    return labels


def _select_missing_types(
    truth_df: pd.DataFrame,
    stage4_summary: dict | None,
    explicit: List[str],
    eps: float = 1e-6,
) -> List[str]:
    if explicit:
        return explicit
    summary = stage4_summary or {}
    missing = summary.get("missing_type_truth")
    if isinstance(missing, str) and missing:
        return [missing]
    if isinstance(missing, list):
        return [m for m in missing if isinstance(m, str)]
    sums = truth_df.sum(axis=0)
    return [t for t in sums.index if float(sums[t]) <= eps]


def _resolve_truth_cell_path(stage1_export: Path, override: str | None) -> Path | None:
    if override:
        path = Path(override)
        return path if path.exists() else None
    for name in [
        "sim_truth_query_cell_spot.csv",
        "sim_truth_cell_assignment.csv",
        "truth_cell_assignment.csv",
    ]:
        cand = stage1_export / name
        if cand.exists():
            return cand
    return None


def _load_truth_cells(path: Path, alias_map: Dict[str, str] | None) -> pd.DataFrame:
    df = pd.read_csv(path)
    if len(df.columns) == 1 and "\t" in df.columns[0]:
        df = pd.read_csv(path, sep="\t")
    if "true_spot_id" in df.columns and "spot_id" not in df.columns:
        df = df.rename(columns={"true_spot_id": "spot_id"})
    if "cell_type" not in df.columns:
        for cand in ["CellType", "celltype", "cell_type_label"]:
            if cand in df.columns:
                df = df.rename(columns={cand: "cell_type"})
                break
    if "spot_id" not in df.columns or "cell_type" not in df.columns:
        return pd.DataFrame()
    df["spot_id"] = df["spot_id"].astype(str).map(_clean_spot_id)
    df["cell_type"] = df["cell_type"].astype(str).apply(
        lambda t: canonicalize_type_name(t, alias_map) if alias_map else normalize_type_name(t)
    )
    return df[["spot_id", "cell_type"]].copy()


def _estimate_spot_spacing(coords_df: pd.DataFrame, x_col: str, y_col: str) -> float:
    xs = np.unique(coords_df[x_col].to_numpy(dtype=float))
    ys = np.unique(coords_df[y_col].to_numpy(dtype=float))
    dx = np.diff(np.sort(xs))
    dy = np.diff(np.sort(ys))
    dx = dx[dx > 0]
    dy = dy[dy > 0]
    vals = []
    if dx.size:
        vals.append(float(np.median(dx)))
    if dy.size:
        vals.append(float(np.median(dy)))
    return min(vals) if vals else 1.0


def _jitter_points(x: np.ndarray, y: np.ndarray, radius: float, seed: int = 0) -> Tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    angles = rng.uniform(0, 2 * np.pi, size=x.shape[0])
    r = np.sqrt(rng.uniform(0, 1.0, size=x.shape[0])) * radius
    return x + r * np.cos(angles), y + r * np.sin(angles)


def _convex_hull(points: np.ndarray) -> np.ndarray:
    if points.shape[0] < 3:
        return points
    pts = np.unique(points, axis=0)
    if pts.shape[0] < 3:
        return pts
    pts = pts[np.lexsort((pts[:, 1], pts[:, 0]))]

    def cross(o, a, b) -> float:
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower: List[np.ndarray] = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper: List[np.ndarray] = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    hull = np.array(lower[:-1] + upper[:-1])
    return hull


def _save_fig(fig, path: Path, dpi: int) -> None:
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _add_tissue_outline(ax, coords_df: pd.DataFrame, x_col: str, y_col: str) -> None:
    points = coords_df[[x_col, y_col]].to_numpy(dtype=float)
    hull = _convex_hull(points)
    if hull.shape[0] < 3:
        return
    closed = np.vstack([hull, hull[0]])
    ax.plot(closed[:, 0], closed[:, 1], color="#c7c7c7", linewidth=0.6, zorder=1)


def _scatter_spots(
    ax,
    coords_df: pd.DataFrame,
    x_col: str,
    y_col: str,
    values: List[Any],
    size: float,
    marker: str,
    edgecolor: str,
    linewidth: float,
) -> None:
    ax.scatter(
        coords_df[x_col],
        coords_df[y_col],
        c=values,
        s=size,
        marker=marker,
        edgecolor=edgecolor,
        linewidth=linewidth,
        zorder=2,
    )
    _add_tissue_outline(ax, coords_df, x_col, y_col)
    ax.set_aspect("equal")
    ax.set_axis_off()


def _plot_type_legend(
    out_path: Path,
    types: List[str],
    colors: Dict[str, Any],
    dpi: int,
    title: str | None = None,
) -> None:
    if not _maybe_plot() or not types:
        return
    ncol = min(4, max(1, math.ceil(len(types) / 6)))
    fig, ax = plt.subplots(figsize=(6.5, 2.5))
    ax.axis("off")
    handles = [
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=colors[t], markersize=6, label=t)
        for t in types
    ]
    legend = ax.legend(
        handles=handles,
        loc="center",
        ncol=ncol,
        frameon=False,
        fontsize=9,
        title=title,
        title_fontsize=10,
    )
    legend.set_zorder(5)
    _save_fig(fig, out_path, dpi)


def _write_meta(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def plot_dominant_type_map(
    out_path: Path,
    coords_df: pd.DataFrame,
    x_col: str,
    y_col: str,
    truth_df: pd.DataFrame,
    types: List[str],
    colors: Dict[str, Any],
    alias_map: Dict[str, str] | None,
    dpi: int,
    spot_size: float,
    marker: str,
) -> None:
    if not _maybe_plot():
        return
    unknown_label = "Unknown"
    labels = _major_type_labels(truth_df[types], types, unknown_label)
    color_vals = [colors.get(lbl, "#999999") for lbl in labels]
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    _scatter_spots(
        ax,
        coords_df,
        x_col,
        y_col,
        color_vals,
        size=spot_size,
        marker=marker,
        edgecolor="#ffffff",
        linewidth=0.3,
    )
    ax.set_title("Dominant type map", fontsize=STYLE["subplot"])
    _save_fig(fig, out_path, dpi)


def plot_fraction_map(
    out_path: Path,
    coords_df: pd.DataFrame,
    x_col: str,
    y_col: str,
    truth_df: pd.DataFrame,
    cell_type: str,
    base_color: Any,
    dpi: int,
    spot_size: float,
    marker: str,
) -> None:
    if not _maybe_plot() or cell_type not in truth_df.columns:
        return
    vals = truth_df[cell_type].to_numpy(dtype=float)
    cmap = LinearSegmentedColormap.from_list("fraction", ["#ffffff", base_color])
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    sc = ax.scatter(
        coords_df[x_col],
        coords_df[y_col],
        c=vals,
        cmap=cmap,
        vmin=0,
        vmax=1,
        s=spot_size,
        marker=marker,
        edgecolor="#ffffff",
        linewidth=0.3,
        zorder=2,
    )
    _add_tissue_outline(ax, coords_df, x_col, y_col)
    ax.set_title(f"Fraction heatmap: {cell_type}", fontsize=STYLE["subplot"])
    ax.set_aspect("equal")
    ax.set_axis_off()
    cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("fraction (0-1)")
    _save_fig(fig, out_path, dpi)


def plot_scalar_map(
    out_path: Path,
    coords_df: pd.DataFrame,
    x_col: str,
    y_col: str,
    values: np.ndarray,
    title: str,
    cmap: str,
    label: str,
    dpi: int,
    spot_size: float,
    marker: str,
) -> None:
    if not _maybe_plot():
        return
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    sc = ax.scatter(
        coords_df[x_col],
        coords_df[y_col],
        c=values,
        cmap=cmap,
        vmin=0,
        vmax=1,
        s=spot_size,
        marker=marker,
        edgecolor="#ffffff",
        linewidth=0.3,
        zorder=2,
    )
    _add_tissue_outline(ax, coords_df, x_col, y_col)
    ax.set_title(title, fontsize=STYLE["subplot"])
    ax.set_aspect("equal")
    ax.set_axis_off()
    cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(label)
    _save_fig(fig, out_path, dpi)


def plot_cells_type_map(
    out_path: Path,
    coords_df: pd.DataFrame,
    x_col: str,
    y_col: str,
    cell_df: pd.DataFrame,
    cell_types: List[str],
    colors: Dict[str, Any],
    dpi: int,
) -> None:
    if not _maybe_plot() or cell_df.empty:
        return
    merged = cell_df.merge(coords_df[["spot_id", x_col, y_col]], on="spot_id", how="inner")
    if merged.empty:
        return
    spacing = _estimate_spot_spacing(coords_df, x_col, y_col)
    jitter_radius = spacing * 0.35
    xs = merged[x_col].to_numpy(dtype=float)
    ys = merged[y_col].to_numpy(dtype=float)
    jx, jy = _jitter_points(xs, ys, jitter_radius, seed=0)
    labels = merged["cell_type"].astype(str).tolist()
    color_vals = [colors.get(lbl, "#999999") for lbl in labels]

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    ax.scatter(jx, jy, c=color_vals, s=4, edgecolor="none", alpha=0.9, zorder=3)
    _add_tissue_outline(ax, coords_df, x_col, y_col)
    ax.set_title("GT cell-level type map", fontsize=STYLE["subplot"])
    ax.set_aspect("equal")
    ax.set_axis_off()
    _save_fig(fig, out_path, dpi)


def _compute_mixing_metrics(truth_df: pd.DataFrame) -> Dict[str, np.ndarray]:
    vals = truth_df.to_numpy(dtype=float)
    row_sum = np.sum(vals, axis=1)
    denom = np.where(row_sum == 0, 1.0, row_sum)
    p = vals / denom[:, None]
    purity = np.max(p, axis=1)
    sorted_vals = np.sort(p, axis=1)
    top1 = sorted_vals[:, -1]
    top2 = sorted_vals[:, -2] if p.shape[1] > 1 else np.zeros_like(top1)
    gap = top1 - top2
    top2mix = 1.0 - gap
    p_safe = np.where(p > 0, p, 1.0)
    entropy = -np.sum(p * np.log(p_safe), axis=1)
    max_entropy = math.log(p.shape[1]) if p.shape[1] > 1 else 1.0
    entropy_norm = entropy / max_entropy if max_entropy > 0 else entropy
    return {
        "purity": np.clip(purity, 0.0, 1.0),
        "entropy": np.clip(entropy_norm, 0.0, 1.0),
        "top2gap": np.clip(top2mix, 0.0, 1.0),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--baseline_id", default="baseline")
    parser.add_argument("--truth_fraction_path", default=None)
    parser.add_argument("--truth_cell_path", default=None)
    parser.add_argument("--missing_types", default="")
    parser.add_argument("--out_dir", default=None)
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()

    if not _maybe_plot():
        print("[sim_truth_visuals] matplotlib not available.")
        return 1

    project_root = _ROOT
    alias_map = load_alias_map(project_root / "configs" / "aliases.yaml")
    stage1_export = project_root / "data" / "processed" / args.sample / "stage1_preprocess" / "exported"
    coords_df, x_col, y_col = _load_coords(stage1_export / "st_coordinates.csv")
    if coords_df is None or x_col is None or y_col is None:
        print("[sim_truth_visuals] missing coordinates.")
        return 1

    truth_path = _resolve_truth_fraction_path(stage1_export, args.truth_fraction_path)
    if truth_path is None or not truth_path.exists():
        print("[sim_truth_visuals] missing truth fraction table.")
        return 1

    truth_frac, _, _ = _load_fraction_table(truth_path)
    truth_frac, truth_types = _canonicalize_fraction_df(truth_frac, alias_map)
    truth_types = sorted(truth_types)
    if not truth_types:
        print("[sim_truth_visuals] no truth types found.")
        return 1

    merged = coords_df.merge(truth_frac, on="spot_id", how="inner")
    if merged.empty:
        print("[sim_truth_visuals] no shared spots between coords and truth.")
        return 1
    coords_df = merged[["spot_id", x_col, y_col]].copy()
    truth_frac = merged.drop(columns=[x_col, y_col])

    stage4_root = project_root / "result" / args.sample / "stage4_cytospace"
    base_run_dir = _resolve_stage4_run_dir(stage4_root, args.baseline_id)
    stage4_summary = {}
    if base_run_dir is not None:
        summary_path = base_run_dir / "stage4_summary.json"
        if summary_path.exists():
            stage4_summary = json.loads(summary_path.read_text(encoding="utf-8"))

    explicit_missing = [m.strip() for m in args.missing_types.split(",") if m.strip()]
    missing_types = _select_missing_types(truth_frac.drop(columns=["spot_id"]), stage4_summary, explicit_missing)
    missing_types = [canonicalize_type_name(t, alias_map) if alias_map else normalize_type_name(t) for t in missing_types]
    missing_not_found = [t for t in missing_types if t not in truth_frac.columns]
    if missing_not_found:
        print(f"[sim_truth_visuals] missing types not in truth fraction: {missing_not_found}")
        return 1
    missing_types = [t for t in missing_types if t in truth_frac.columns]

    out_dir = Path(args.out_dir) if args.out_dir else project_root / "result" / args.sample / "ground_truth_visuals"
    _ensure_dir(out_dir)

    truth_cell_path = _resolve_truth_cell_path(stage1_export, args.truth_cell_path)
    truth_cells = pd.DataFrame()
    cell_types: List[str] = []
    if truth_cell_path is not None and truth_cell_path.exists():
        truth_cells = _load_truth_cells(truth_cell_path, alias_map)
        if not truth_cells.empty:
            cell_types = sorted(truth_cells["cell_type"].dropna().unique().tolist())
        else:
            print("[sim_truth_visuals] truth cell assignment empty or missing columns; skip cell-level map.")
    else:
        print("[sim_truth_visuals] truth cell assignment not found; skip cell-level map.")

    all_types = list(truth_types)
    for t in sorted(set(cell_types) - set(all_types)):
        all_types.append(t)
    colors = _build_category_colors(all_types)

    spacing = _estimate_spot_spacing(coords_df, x_col, y_col)
    spot_size = max(18.0, min(40.0, spacing * 8.0))
    marker = "h"

    if not truth_cells.empty:
        plot_cells_type_map(
            out_dir / "gt_cells_type_map.png",
            coords_df,
            x_col,
            y_col,
            truth_cells,
            cell_types,
            colors,
            dpi=args.dpi,
        )
        _plot_type_legend(out_dir / "gt_cells_type_map_legend.png", cell_types, colors, args.dpi, title="Cell types")
        _write_meta(out_dir / "gt_cells_type_map.meta.json", {
            "type": "gt_cells_type_map",
            "sample": args.sample,
            "truth_cell_path": str(truth_cell_path),
            "coordinate_path": str(stage1_export / "st_coordinates.csv"),
            "colormap": {k: _rgba_to_hex(v) for k, v in colors.items() if k in cell_types},
            "cell_types": cell_types,
        })

    plot_dominant_type_map(
        out_dir / "dominant_type_map.png",
        coords_df,
        x_col,
        y_col,
        truth_frac,
        truth_types,
        colors,
        alias_map,
        dpi=args.dpi,
        spot_size=spot_size,
        marker=marker,
    )
    _plot_type_legend(out_dir / "dominant_type_map_legend.png", truth_types, colors, args.dpi, title="Cell types")
    _write_meta(out_dir / "dominant_type_map.meta.json", {
        "type": "dominant_type_map",
        "sample": args.sample,
        "truth_fraction_path": str(truth_path),
        "coordinate_path": str(stage1_export / "st_coordinates.csv"),
        "colormap": {k: _rgba_to_hex(v) for k, v in colors.items()},
        "cell_types": truth_types,
        "missing_types": missing_types,
    })

    if missing_types:
        frac_dir = out_dir / "fractions"
        _ensure_dir(frac_dir)
        for t in missing_types:
            slug = _safe_name(t)
            plot_fraction_map(
                frac_dir / f"fraction_{slug}.png",
                coords_df,
                x_col,
                y_col,
                truth_frac,
                t,
                colors.get(t, "#ff7f0e"),
                dpi=args.dpi,
                spot_size=spot_size,
                marker=marker,
            )
            _write_meta(frac_dir / f"fraction_{slug}.meta.json", {
                "type": "fraction_heatmap",
                "sample": args.sample,
                "truth_fraction_path": str(truth_path),
                "coordinate_path": str(stage1_export / "st_coordinates.csv"),
                "cell_type": t,
                "value_range": [0.0, 1.0],
                "base_color": _rgba_to_hex(colors.get(t, "#ff7f0e")),
                "missing_in_ref": True,
            })
    else:
        print("[sim_truth_visuals] no missing types detected; skip fraction heatmaps.")

    mixing_dir = out_dir / "mixing"
    _ensure_dir(mixing_dir)
    mixing_metrics = _compute_mixing_metrics(truth_frac.drop(columns=["spot_id"]))
    plot_scalar_map(
        mixing_dir / "mixing_entropy.png",
        coords_df,
        x_col,
        y_col,
        mixing_metrics["entropy"],
        title="Mixing map (entropy)",
        cmap="viridis",
        label="entropy (0-1)",
        dpi=args.dpi,
        spot_size=spot_size,
        marker=marker,
    )
    _write_meta(mixing_dir / "mixing_entropy.meta.json", {
        "type": "mixing_map",
        "sample": args.sample,
        "truth_fraction_path": str(truth_path),
        "coordinate_path": str(stage1_export / "st_coordinates.csv"),
        "metric": "entropy",
        "normalized": True,
        "value_range": [0.0, 1.0],
    })
    plot_scalar_map(
        mixing_dir / "mixing_purity.png",
        coords_df,
        x_col,
        y_col,
        mixing_metrics["purity"],
        title="Mixing map (purity)",
        cmap="viridis",
        label="purity (0-1)",
        dpi=args.dpi,
        spot_size=spot_size,
        marker=marker,
    )
    _write_meta(mixing_dir / "mixing_purity.meta.json", {
        "type": "mixing_map",
        "sample": args.sample,
        "truth_fraction_path": str(truth_path),
        "coordinate_path": str(stage1_export / "st_coordinates.csv"),
        "metric": "purity",
        "value_range": [0.0, 1.0],
    })
    plot_scalar_map(
        mixing_dir / "mixing_top2gap.png",
        coords_df,
        x_col,
        y_col,
        mixing_metrics["top2gap"],
        title="Mixing map (1 - top2 gap)",
        cmap="viridis",
        label="mixing (0-1)",
        dpi=args.dpi,
        spot_size=spot_size,
        marker=marker,
    )
    _write_meta(mixing_dir / "mixing_top2gap.meta.json", {
        "type": "mixing_map",
        "sample": args.sample,
        "truth_fraction_path": str(truth_path),
        "coordinate_path": str(stage1_export / "st_coordinates.csv"),
        "metric": "top2gap",
        "note": "1 - (top1 - top2)",
        "value_range": [0.0, 1.0],
    })

    print(f"[sim_truth_visuals] outputs -> {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
