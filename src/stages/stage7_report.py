"""
Stage7 report generator for real data (baseline vs route2).
Outputs figures, tables, JSON summary, and a lightweight HTML/Markdown report.
"""
from __future__ import annotations

import argparse
import base64
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd

try:
    import yaml  # type: ignore
except Exception:  # pragma: no cover
    yaml = None

try:  # pragma: no cover
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # type: ignore
    from matplotlib.patches import Polygon, Rectangle  # type: ignore
    from matplotlib.colors import TwoSlopeNorm  # type: ignore
except Exception:  # pragma: no cover
    plt = None

# ensure project root in sys.path
_HERE = Path(__file__).resolve()
_ROOT = _HERE.parents[2]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src.utils.type_name import normalize_type_name, canonicalize_type_name, load_alias_map

EPS = 1e-12

STYLE = {
    "title": 17,
    "subplot": 13,
    "label": 11,
    "tick": 9,
}

FIGSIZE_SINGLE = (6.5, 4.5)
FIGSIZE_WIDE = (12.0, 4.8)
FIGSIZE_GRID = (12.0, 8.0)
PAD_INCHES = 0.15
OUTPUT_FORMATS = ["png"]

TYPE_SHORT_MAP = {
    "Monocytes and Macrophages": "Mono/Mφ",
    "Endothelial cells": "Endo",
    "Epithelial cells": "Epi",
    "T cells CD4": "CD4 T",
    "T cells CD8": "CD8 T",
    "Unknown_sc_only": "Unknown(sc-only)",
}


def _apply_style() -> None:
    if not _maybe_plot():
        return
    plt.rcParams.update({
        "axes.titlesize": STYLE["subplot"],
        "axes.labelsize": STYLE["label"],
        "xtick.labelsize": STYLE["tick"],
        "ytick.labelsize": STYLE["tick"],
        "figure.titlesize": STYLE["title"],
        "legend.fontsize": STYLE["tick"],
    })


def _short_type(name: str, alias_map: Dict[str, str] | None) -> str:
    if not name:
        return name
    if name in TYPE_SHORT_MAP:
        return TYPE_SHORT_MAP[name]
    canon = canonicalize_type_name(name, alias_map) if alias_map else name
    if canon in TYPE_SHORT_MAP:
        return TYPE_SHORT_MAP[canon]
    norm = normalize_type_name(name)
    if norm in TYPE_SHORT_MAP:
        return TYPE_SHORT_MAP[norm]
    return name


def _build_type_abbrev(types: Iterable[str], alias_map: Dict[str, str] | None) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    for t in types:
        short = _short_type(t, alias_map)
        if short != t:
            mapping[t] = short
    return mapping


def _wrap_label(label: str) -> str:
    if len(label) <= 8:
        return label
    if "/" in label:
        return label.replace("/", "/\n")
    if " " in label:
        parts = label.split()
        if len(parts) >= 2:
            return parts[0] + "\n" + " ".join(parts[1:])
    return label


def _collect_must_include(
    base_summary: dict,
    route_summary: dict,
    base_mass: pd.Series,
    route_mass: pd.Series,
) -> List[str]:
    must: List[str] = []
    for summary in (base_summary, route_summary):
        missing = summary.get("missing_type_truth")
        if isinstance(missing, str) and missing:
            must.append(missing)
    for t in base_mass.index:
        if base_mass.get(t, 0.0) > 0 and route_mass.get(t, 0.0) == 0:
            must.append(t)
    for t in route_mass.index:
        if route_mass.get(t, 0.0) > 0 and base_mass.get(t, 0.0) == 0:
            must.append(t)
    out = []
    for t in must:
        if t not in out:
            out.append(t)
    return out


def detect_root() -> Path:
    here = Path(__file__).resolve()
    try:
        return here.parents[2]
    except IndexError:
        return here.parent.parent


def _read_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _load_yaml(path: Path) -> dict:
    if yaml is None or not path.exists():
        return {}
    try:
        with path.open("r", encoding="utf-8") as f:
            return yaml.safe_load(f) or {}
    except Exception:
        return {}


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _clean_spot_id(value: Any) -> str:
    return str(value).split()[0].split("\t")[0]


def _pick_spot_col(columns: Iterable[str]) -> str:
    for cand in ["spot_id", "SpotID", "spot", "Spot"]:
        if cand in columns:
            return cand
    return next(iter(columns))


def _resolve_stage4_run_dir(stage4_root: Path, run_id: str) -> Path | None:
    run_dir = stage4_root / "runs" / run_id
    if run_dir.exists():
        return run_dir
    fallback = stage4_root / "cytospace_output"
    if fallback.exists():
        return fallback
    return None


def _resolve_stage6_run_dir(stage6_root: Path, run_id: str) -> Path | None:
    run_dir = stage6_root / run_id
    return run_dir if run_dir.exists() else None


def _resolve_stage5_root(result_root: Path) -> Path | None:
    preferred = result_root / "stage5_route2_s0"
    if preferred.exists():
        return preferred
    candidates = sorted([p for p in result_root.iterdir() if p.is_dir() and p.name.startswith("stage5")])
    return candidates[0] if candidates else None


def _resolve_stage5_run_dir(stage5_root: Path | None, run_id: str) -> Path | None:
    if stage5_root is None:
        return None
    direct = stage5_root / run_id
    if direct.exists():
        return direct
    alt = run_id
    for prefix in ("route2_s0_", "route2_"):
        if alt.startswith(prefix):
            alt = alt[len(prefix):]
    if alt != run_id:
        alt_dir = stage5_root / alt
        if alt_dir.exists():
            return alt_dir
    matches = [p for p in stage5_root.iterdir() if p.is_dir() and alt in p.name]
    return sorted(matches)[0] if matches else None


def _find_stage5_json(run_dir: Path | None) -> Path | None:
    if run_dir is None or not run_dir.exists():
        return None
    files = sorted(run_dir.glob("stage5_route2_s0__*.json"))
    if not files:
        return None
    for f in files:
        if run_dir.name in f.name:
            return f
    return files[0]


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


def _align_fraction_tables(
    base_df: pd.DataFrame, route_df: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame, List[str]]:
    base = base_df.set_index("spot_id")
    route = route_df.set_index("spot_id")
    all_types = sorted(set(base.columns) | set(route.columns))
    base = base.reindex(columns=all_types, fill_value=0.0)
    route = route.reindex(columns=all_types, fill_value=0.0)
    return base, route, all_types


def _align_three_tables(
    base_df: pd.DataFrame, route_df: pd.DataFrame, truth_df: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, List[str]]:
    base = base_df.set_index("spot_id")
    route = route_df.set_index("spot_id")
    truth = truth_df.set_index("spot_id")
    common = base.index.intersection(route.index).intersection(truth.index)
    base = base.loc[common]
    route = route.loc[common]
    truth = truth.loc[common]
    all_types = sorted(set(base.columns) | set(route.columns) | set(truth.columns))
    base = base.reindex(columns=all_types, fill_value=0.0)
    route = route.reindex(columns=all_types, fill_value=0.0)
    truth = truth.reindex(columns=all_types, fill_value=0.0)
    return base, route, truth, all_types


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


def _collect_marker_candidates(cfg: dict) -> List[str]:
    preferred = ["PTPRC", "EPCAM", "KRT19", "COL1A1", "VWF", "PECAM1"]
    markers: List[str] = []
    stage6 = cfg.get("stage6", {}) if isinstance(cfg, dict) else {}
    for entry in stage6.get("marker_sets", []) or []:
        if isinstance(entry, dict):
            markers.extend([m for m in entry.get("markers", []) if isinstance(m, str)])
    for entry in stage6.get("composites", []) or []:
        if isinstance(entry, dict):
            markers.extend([m for m in entry.get("markers", []) if isinstance(m, str)])
    ordered = []
    for m in preferred:
        if m not in ordered:
            ordered.append(m)
    for m in markers:
        if m not in ordered:
            ordered.append(m)
    return ordered


def _collect_compartments(cfg: dict) -> List[dict]:
    stage6 = cfg.get("stage6", {}) if isinstance(cfg, dict) else {}
    compartments = stage6.get("compartments", {}) if isinstance(stage6, dict) else {}
    mapping = compartments.get("mapping") if isinstance(compartments, dict) else None
    if isinstance(mapping, list):
        return [m for m in mapping if isinstance(m, dict)]
    return []


def _extract_stage5_metrics(stage5_json: dict | None) -> Dict[str, float]:
    stage5_json = stage5_json or {}
    comp = stage5_json.get("composition", {})
    leak = stage5_json.get("leakage", {})
    cell_metrics = stage5_json.get("cell_spot_eval", {}).get("metrics", {})
    out: Dict[str, float] = {}
    for key in ("L1_mean", "JS_mean", "corr_mean"):
        val = comp.get(key)
        if val is not None:
            out[key] = float(val)
    miss = leak.get("missing_leak_rate")
    if miss is not None:
        out["missing_leak_rate"] = float(miss)
    extra = cell_metrics.get("extra_pred_fraction")
    if extra is not None:
        out["extra_pred_fraction"] = float(extra)
    acc = cell_metrics.get("acc_top1_cond_nonmissing")
    if acc is not None:
        out["acc_top1_cond_nonmissing"] = float(acc)
    return out


def _maybe_plot() -> bool:
    return plt is not None


def _save_fig(fig, path: Path, dpi: int, tight: bool = False) -> None:
    if tight:
        fig.tight_layout()
    for fmt in OUTPUT_FORMATS:
        out_path = path.with_suffix(f".{fmt}")
        fig.savefig(out_path, dpi=dpi if fmt == "png" else None, bbox_inches="tight", pad_inches=PAD_INCHES)
    plt.close(fig)


def plot_type_composition(
    out_path: Path,
    base_mass: pd.Series,
    route_mass: pd.Series,
    types: List[str],
    alias_map: Dict[str, str] | None,
    must_include: List[str],
    label_fmt: str,
    break_ratio: float,
    break_min: float,
    dpi: int,
) -> List[str]:
    if not _maybe_plot():
        return []
    all_types = sorted(set(base_mass.index) | set(route_mass.index))
    df = pd.DataFrame({
        "type": all_types,
        "baseline": [base_mass.get(t, 0.0) for t in all_types],
        "route2": [route_mass.get(t, 0.0) for t in all_types],
    })
    df["delta"] = df["route2"] - df["baseline"]
    df["abs_delta"] = df["delta"].abs()
    if types:
        df = df[df["type"].isin(types)]
    df = df.sort_values("delta", ascending=True)

    labels = [_short_type(t, alias_map) for t in df["type"].tolist()]
    y = np.arange(len(df))
    colors = ["#e15759" if v < 0 else "#59a14f" for v in df["delta"]]

    abs_vals = df["abs_delta"].to_numpy(dtype=float)
    broken = False
    if abs_vals.size >= 2:
        sorted_abs = np.sort(abs_vals)
        broken = (sorted_abs[-1] > break_min) and (sorted_abs[-1] > break_ratio * sorted_abs[-2])

    def _fmt(val: float) -> str:
        if label_fmt == "int":
            return f"{val:+.0f}"
        if label_fmt == "1f":
            return f"{val:+.1f}"
        return f"{val:+.2f}"

    xlabel = "delta mass (route2 - baseline)\n(positive = increase, negative = decrease)"
    if broken:
        fig, (ax1, ax2) = plt.subplots(
            ncols=2,
            figsize=(10.5, 4.5),
            gridspec_kw={"width_ratios": [4, 1], "wspace": 0.06},
            sharey=True,
        )
        second_max = sorted_abs[-2] if abs_vals.size >= 2 else sorted_abs[-1]
        limit = max(break_min, second_max * 1.2)
        for ax in (ax1, ax2):
            ax.axvline(0, color="#111111", linewidth=1.6)
            ax.grid(axis="x", linestyle="--", alpha=0.25)

        for i, (delta, color) in enumerate(zip(df["delta"], colors)):
            if abs(delta) <= limit:
                ax = ax1
            else:
                ax = ax2
            ax.hlines(i, 0, delta, color=color, linewidth=2, alpha=0.85)
            ax.scatter(delta, i, color=color, s=60, zorder=3, edgecolor="white", linewidth=0.6)
            offset = 1.5 if delta >= 0 else -1.5
            ax.text(delta + offset, i, _fmt(delta), va="center", ha="left" if delta >= 0 else "right", fontsize=9)

        ax1.set_xlim(-limit, limit)
        ax2_min = float(df["delta"].min())
        ax2_max = float(df["delta"].max())
        pad = max(5.0, 0.15 * (ax2_max - ax2_min))
        ax2.set_xlim(ax2_min - pad, ax2_max + pad)

        ax1.set_yticks(y)
        ax1.set_yticklabels(labels, fontsize=STYLE["tick"])
        ax1.set_xlabel(xlabel, labelpad=12)
        ax2.set_xlabel("outlier", labelpad=12)
        ax1.set_title("Top delta types\nΔmass = sum_spot_fractions(route2) - sum_spot_fractions(baseline)\nunit: spot-equivalent mass")
        ax2.set_title("Outlier")

        ax1.spines["right"].set_visible(False)
        ax2.spines["left"].set_visible(False)
        ax1.tick_params(labelright=False)
        ax2.tick_params(labelleft=False)
        d = 0.015
        ax1.plot((1 - d, 1 + d), (-d, +d), transform=ax1.transAxes, color="k", clip_on=False, linewidth=1)
        ax1.plot((1 - d, 1 + d), (1 - d, 1 + d), transform=ax1.transAxes, color="k", clip_on=False, linewidth=1)
        ax2.plot((-d, +d), (-d, +d), transform=ax2.transAxes, color="k", clip_on=False, linewidth=1)
        ax2.plot((-d, +d), (1 - d, 1 + d), transform=ax2.transAxes, color="k", clip_on=False, linewidth=1)
    else:
        fig, ax = plt.subplots(figsize=FIGSIZE_SINGLE)
        ax.axvline(0, color="#111111", linewidth=1.6)
        ax.grid(axis="x", linestyle="--", alpha=0.25)
        ax.hlines(y, 0, df["delta"], color=colors, linewidth=2, alpha=0.85)
        ax.scatter(df["delta"], y, color=colors, s=60, zorder=3, edgecolor="white", linewidth=0.6)
        for i, delta in enumerate(df["delta"]):
            offset = 1.5 if delta >= 0 else -1.5
            ax.text(delta + offset, i, _fmt(delta), va="center", ha="left" if delta >= 0 else "right", fontsize=9)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=STYLE["tick"])
        ax.set_xlabel(xlabel, labelpad=12)
        ax.set_title("Top delta types\nΔmass = sum_spot_fractions(route2) - sum_spot_fractions(baseline)\nunit: spot-equivalent mass")

    if must_include:
        note = "must-include: " + ", ".join([_short_type(t, alias_map) for t in must_include])
        fig.text(0.01, 0.02, note, fontsize=9, ha="left")
    fig.subplots_adjust(bottom=0.18)

    _save_fig(fig, out_path, dpi, tight=False)
    return df["type"].tolist()


def _extract_flow_counts(stage4_summary: dict | None, audit: dict | None) -> dict[str, int]:
    stage4_summary = stage4_summary or {}
    audit = audit or {}
    n_before = int(stage4_summary.get("n_cells_before_prefilter") or 0)
    n_filtered = int(stage4_summary.get("n_filtered") or 0)
    n_after = int(stage4_summary.get("n_cells_after_prefilter") or max(n_before - n_filtered, 0))

    rescue_ledger = (
        audit.get("ledger_integrity", {})
        .get("metrics", {})
        .get("rescue_ledger", {})
    )
    rescued = rescue_ledger.get("rescued_cells_total")
    if rescued is None:
        rescued = rescue_ledger.get("rescued_cells_count")
    if rescued is None:
        rescued = 0
    rescued = int(rescued)
    rescued = max(0, min(rescued, n_filtered))
    unknown = max(n_filtered - rescued, 0)

    filter_ledger = (
        audit.get("ledger_integrity", {})
        .get("metrics", {})
        .get("filter_ledger", {})
    )
    expected = filter_ledger.get("expected_n_filtered")
    observed = filter_ledger.get("observed_n_filtered")
    if expected is None:
        expected = n_filtered
    if observed is None:
        observed = n_filtered
    expected = int(expected)
    observed = int(observed)
    ledger_ok = filter_ledger.get("ledger_check_ok")
    if ledger_ok is None:
        ledger_ok = expected == observed

    return {
        "n_before": n_before,
        "n_filtered": n_filtered,
        "n_after": n_after,
        "n_rescued": rescued,
        "n_unknown": unknown,
        "expected_filtered": expected,
        "observed_filtered": observed,
        "ledger_ok": bool(ledger_ok),
    }


def _add_flow(ax, x0: float, x1: float, y0: float, y1: float, y2: float, y3: float, color: str) -> None:
    verts = [(x0, y0), (x0, y1), (x1, y3), (x1, y2)]
    ax.add_patch(Polygon(verts, closed=True, facecolor=color, edgecolor=color, linewidth=0.6, alpha=0.4))


def _add_zero_flow(ax, x0: float, x1: float, y: float, color: str) -> None:
    ax.plot([x0, x1], [y, y], color=color, linewidth=1.2, alpha=0.7)


def _draw_flow_panel(
    ax,
    title: str,
    counts: dict[str, int],
    scale: float,
    show_zero_flows: bool,
    filtered_label: str,
) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title(title, fontsize=10, pad=10)

    x0, x1, x2 = 0.05, 0.42, 0.79
    w = 0.15
    base_y = 0.10

    n_before = counts["n_before"]
    n_filtered = counts["n_filtered"]
    n_after = counts["n_after"]
    n_rescued = counts["n_rescued"]
    n_unknown = counts["n_unknown"]

    h_before = n_before * scale
    h_filtered = n_filtered * scale
    h_after = n_after * scale
    h_rescued = n_rescued * scale
    h_unknown = n_unknown * scale

    y0_start = (1.0 - h_before) / 2.0 if h_before > 0 else base_y
    y0_end = y0_start + h_before

    y_after_start = y0_start
    y_after_end = y_after_start + h_after
    gap = 0.008
    y_filtered_start = y_after_end + gap if h_filtered > 0 and h_after > 0 else y_after_end
    y_filtered_end = y_filtered_start + h_filtered

    y_unknown_start = y_filtered_start
    y_unknown_end = y_unknown_start + h_unknown
    y_rescued_start = y_unknown_end + gap if h_rescued > 0 and h_unknown > 0 else y_unknown_end
    y_rescued_end = y_rescued_start + h_rescued

    colors = {
        "input": "#9aa0a6",
        "filtered": "#f28e2b",
        "kept": "#4e79a7",
        "rescued": "#59a14f",
        "unknown": "#e15759",
    }

    if h_filtered > 0:
        _add_flow(ax, x0 + w, x1, y_filtered_start, y_filtered_end, y_filtered_start, y_filtered_end, colors["filtered"])
    elif show_zero_flows:
        _add_zero_flow(ax, x0 + w, x1, y_filtered_start, "#bbbbbb")
    if h_after > 0:
        _add_flow(ax, x0 + w, x1, y_after_start, y_after_end, y_after_start, y_after_end, colors["kept"])
    elif show_zero_flows:
        _add_zero_flow(ax, x0 + w, x1, y_after_start, "#bbbbbb")
    if h_rescued > 0:
        _add_flow(ax, x1 + w, x2, y_rescued_start, y_rescued_end, y_rescued_start, y_rescued_end, colors["rescued"])
    elif show_zero_flows:
        _add_zero_flow(ax, x1 + w, x2, y_rescued_start, "#bbbbbb")
    if h_unknown > 0:
        _add_flow(ax, x1 + w, x2, y_unknown_start, y_unknown_end, y_unknown_start, y_unknown_end, colors["unknown"])
    elif show_zero_flows:
        _add_zero_flow(ax, x1 + w, x2, y_unknown_start, "#bbbbbb")

    if h_before > 0:
        ax.add_patch(Rectangle((x0, y0_start), w, h_before, facecolor=colors["input"], edgecolor="#333333", linewidth=0.8))
    if h_filtered > 0 or show_zero_flows:
        if h_filtered > 0:
            ax.add_patch(Rectangle((x1, y_filtered_start), w, h_filtered, facecolor=colors["filtered"], edgecolor="#333333", linewidth=0.8))
    if h_after > 0 or show_zero_flows:
        if h_after > 0:
            ax.add_patch(Rectangle((x1, y_after_start), w, h_after, facecolor=colors["kept"], edgecolor="#333333", linewidth=0.8))
    if h_rescued > 0 or show_zero_flows:
        if h_rescued > 0:
            ax.add_patch(Rectangle((x2, y_rescued_start), w, h_rescued, facecolor=colors["rescued"], edgecolor="#333333", linewidth=0.8))
    if h_unknown > 0 or show_zero_flows:
        if h_unknown > 0:
            ax.add_patch(Rectangle((x2, y_unknown_start), w, h_unknown, facecolor=colors["unknown"], edgecolor="#333333", linewidth=0.8))

    def _label(x: float, y: float, text: str) -> None:
        ax.text(x, y, text, ha="center", va="center", fontsize=8)

    _label(x0 + w / 2, y0_start + max(h_before, 0.02) / 2, f"input\n{n_before}")
    if h_filtered > 0:
        _label(x1 + w / 2, (y_filtered_start + y_filtered_end) / 2, f"{filtered_label}\n{n_filtered}")
    if h_after > 0:
        _label(x1 + w / 2, (y_after_start + y_after_end) / 2, f"kept\n{n_after}")
    if h_rescued > 0:
        _label(x2 + w / 2, (y_rescued_start + y_rescued_end) / 2, f"rescued\n{n_rescued}")
    if h_unknown > 0:
        _label(x2 + w / 2, (y_unknown_start + y_unknown_end) / 2, f"unknown\n{n_unknown}")

    def _flow_label(xa: float, xb: float, y0: float, y1: float, text: str) -> None:
        y = (y0 + y1) / 2
        if y1 <= y0:
            y = y0 + 0.01
        ax.text((xa + xb) / 2, y, text, ha="center", va="center",
                fontsize=8, bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#666666", lw=0.6))

    if h_after > 0:
        _flow_label(x0 + w, x1, y_after_start, y_after_end, str(n_after))
    if h_filtered > 0:
        _flow_label(x0 + w, x1, y_filtered_start, y_filtered_end, str(n_filtered))
    if h_rescued > 0:
        _flow_label(x1 + w, x2, y_rescued_start, y_rescued_end, str(n_rescued))
    if h_unknown > 0:
        _flow_label(x1 + w, x2, y_unknown_start, y_unknown_end, str(n_unknown))


def plot_filter_rescue(
    out_path: Path,
    baseline_summary: dict | None,
    route_summary: dict | None,
    baseline_audit: dict | None,
    route_audit: dict | None,
    show_zero_flows: bool,
    dpi: int,
) -> None:
    if not _maybe_plot():
        return
    base_counts = _extract_flow_counts(baseline_summary, baseline_audit)
    route_counts = _extract_flow_counts(route_summary, route_audit)
    max_before = max(base_counts["n_before"], route_counts["n_before"], 1)
    scale = 0.80 / max_before

    fig, axes = plt.subplots(ncols=2, figsize=FIGSIZE_WIDE)
    if not isinstance(axes, np.ndarray):
        axes = np.array([axes])
    axes = axes.flatten()
    _draw_flow_panel(axes[0], "Baseline (no filter/rescue)", base_counts, scale, show_zero_flows, "filtered")
    _draw_flow_panel(axes[1], "Route2 (filter + rescue)", route_counts, scale, show_zero_flows, "filtered\n(marked unknown)")

    def _badge_text(counts: dict[str, int]) -> tuple[str, str]:
        status = "PASS" if counts.get("ledger_ok") else "FAIL"
        diff = counts.get("observed_filtered", 0) - counts.get("expected_filtered", 0)
        text = (
            f"ledger: {status}\n"
            f"expected vs observed: {counts.get('expected_filtered', 0)} vs {counts.get('observed_filtered', 0)}\n"
            f"Δ = {diff:+d}"
        )
        fc = "#e8f5e9" if status == "PASS" else "#ffebee"
        return text, fc

    fig.suptitle("Cell flow / ledger overview", fontsize=STYLE["title"], y=0.98)
    def _pct(val: int, total: int) -> str:
        if total <= 0:
            return "0.0%"
        return f"{(val / total) * 100:.1f}%"

    base_summary = (
        f"Input={base_counts['n_before']}; Filtered={base_counts['n_filtered']} ({_pct(base_counts['n_filtered'], base_counts['n_before'])}); "
        f"Kept={base_counts['n_after']} ({_pct(base_counts['n_after'], base_counts['n_before'])}); "
        f"Rescued={base_counts['n_rescued']} ({_pct(base_counts['n_rescued'], base_counts['n_before'])}); "
        f"Unknown={base_counts['n_unknown']} ({_pct(base_counts['n_unknown'], base_counts['n_before'])})"
    )
    route_summary = (
        f"Input={route_counts['n_before']}; Filtered={route_counts['n_filtered']} ({_pct(route_counts['n_filtered'], route_counts['n_before'])}); "
        f"Kept={route_counts['n_after']} ({_pct(route_counts['n_after'], route_counts['n_before'])}); "
        f"Rescued={route_counts['n_rescued']} ({_pct(route_counts['n_rescued'], route_counts['n_before'])}); "
        f"Unknown={route_counts['n_unknown']} ({_pct(route_counts['n_unknown'], route_counts['n_before'])})"
    )
    base_summary = base_summary.replace("; ", ";\n")
    route_summary = route_summary.replace("; ", ";\n")
    axes[0].text(0.5, -0.14, base_summary, transform=axes[0].transAxes, ha="center", va="top", fontsize=9)
    axes[1].text(0.5, -0.14, route_summary, transform=axes[1].transAxes, ha="center", va="top", fontsize=9)

    fig.subplots_adjust(bottom=0.32, top=0.86)
    base_text, base_fc = _badge_text(base_counts)
    route_text, route_fc = _badge_text(route_counts)
    axes[0].text(
        0.98,
        0.92,
        base_text,
        transform=axes[0].transAxes,
        ha="right",
        va="top",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.3", fc=base_fc, ec="#666666", lw=0.8),
    )
    axes[1].text(
        0.98,
        0.92,
        route_text,
        transform=axes[1].transAxes,
        ha="right",
        va="top",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.3", fc=route_fc, ec="#666666", lw=0.8),
    )
    _save_fig(fig, out_path, dpi, tight=False)


def _select_key_type(
    must_include: List[str],
    available_types: List[str],
    delta_mass: pd.Series | None,
) -> str | None:
    for t in must_include:
        if t in available_types:
            return t
    if delta_mass is not None:
        for t in delta_mass.abs().sort_values(ascending=False).index:
            if t in available_types:
                return t
    return available_types[0] if available_types else None


def _compute_spot_l1(
    pred_df: pd.DataFrame,
    truth_df: pd.DataFrame,
    type_cols: List[str],
) -> np.ndarray:
    pred = pred_df.reindex(columns=type_cols, fill_value=0.0).to_numpy(dtype=float)
    truth = truth_df.reindex(columns=type_cols, fill_value=0.0).to_numpy(dtype=float)
    return np.abs(pred - truth).sum(axis=1)


def _pick_missing_types(truth_df: pd.DataFrame, eps: float = 1e-9) -> List[str]:
    sums = truth_df.sum(axis=0)
    return [t for t in sums.index if float(sums[t]) <= eps]


def _ecdf(values: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    if values.size == 0:
        return np.array([]), np.array([])
    x = np.sort(values)
    y = np.arange(1, x.size + 1) / x.size
    return x, y


def _choose_panel_c_target(
    truth_df: pd.DataFrame,
    pred_base_df: pd.DataFrame,
    pred_route_df: pd.DataFrame,
    eps: float = 1e-6,
) -> Tuple[str, str | None, str | None]:
    type_cols = sorted(set(truth_df.columns) | set(pred_base_df.columns) | set(pred_route_df.columns))
    truth_sum = truth_df.reindex(columns=type_cols, fill_value=0.0).sum(axis=0)
    base_sum = pred_base_df.reindex(columns=type_cols, fill_value=0.0).sum(axis=0)
    candidates = [t for t in type_cols if (float(truth_sum[t]) <= eps) and (float(base_sum[t]) >= 1.0)]
    if candidates:
        t_star = max(candidates, key=lambda t: float(base_sum[t]))
        return "missing_leak", t_star, None

    best = None
    best_gain = -np.inf
    for t in type_cols:
        if t not in truth_df.columns:
            continue
        y_true = truth_df[t].to_numpy(dtype=float)
        y_b = pred_base_df.get(t, pd.Series(0.0, index=truth_df.index)).to_numpy(dtype=float)
        y_r = pred_route_df.get(t, pd.Series(0.0, index=truth_df.index)).to_numpy(dtype=float)
        rmse_b = float(np.sqrt(np.mean((y_b - y_true) ** 2))) if y_true.size else np.inf
        rmse_r = float(np.sqrt(np.mean((y_r - y_true) ** 2))) if y_true.size else np.inf
        if rmse_r >= rmse_b:
            continue
        gain = rmse_b - rmse_r
        if gain > best_gain:
            best_gain = gain
            best = t
    if best is not None:
        return "scatter", best, None
    best = type_cols[0] if type_cols else None
    return "scatter", best, "no rmse improvement found"


def _plot_lollipop_axis(
    ax,
    delta_mass: pd.Series,
    alias_map: Dict[str, str] | None,
    label_fmt: str,
    max_types: int = 6,
) -> None:
    if delta_mass.empty:
        ax.text(0.5, 0.5, "no delta data", ha="center", va="center", fontsize=10)
        ax.set_axis_off()
        return
    top = delta_mass.abs().sort_values(ascending=False).head(max_types).index.tolist()
    vals = delta_mass.loc[top].sort_values()
    y = np.arange(len(vals))
    colors = ["#e15759" if v < 0 else "#59a14f" for v in vals]
    ax.axvline(0, color="#111111", linewidth=1.2)
    ax.hlines(y, 0, vals.values, color=colors, linewidth=2.0, alpha=0.85)
    ax.scatter(vals.values, y, color=colors, s=45, zorder=3, edgecolor="white", linewidth=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels([_short_type(t, alias_map) for t in vals.index], fontsize=STYLE["tick"])
    ax.set_title("Δmass (route2 - baseline)", fontsize=STYLE["subplot"])
    ax.grid(axis="x", linestyle="--", alpha=0.25)
    if label_fmt in {"1f", "2f", "int"}:
        for i, v in enumerate(vals.values):
            if label_fmt == "int":
                lab = f"{v:+.0f}"
            elif label_fmt == "1f":
                lab = f"{v:+.1f}"
            else:
                lab = f"{v:+.2f}"
            ax.text(v, i, lab, va="center", ha="left" if v >= 0 else "right", fontsize=8)


def plot_perf_contrast(
    out_path: Path,
    base_df: pd.DataFrame,
    route_df: pd.DataFrame,
    truth_df: pd.DataFrame | None,
    base_metrics: Dict[str, float],
    route_metrics: Dict[str, float],
    delta_mass: pd.Series,
    alias_map: Dict[str, str] | None,
    must_include: List[str],
    label_fmt: str,
    dpi: int,
) -> None:
    if not _maybe_plot():
        return

    def _idx(df: pd.DataFrame) -> pd.DataFrame:
        if "spot_id" in df.columns:
            return df.set_index("spot_id")
        return df

    fig = plt.figure(figsize=(12.5, 9.0))
    grid = fig.add_gridspec(2, 2, wspace=0.28, hspace=0.35)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    sub_c = grid[1, 0].subgridspec(1, 2, wspace=0.15)
    ax_c1 = fig.add_subplot(sub_c[0, 0])
    ax_c2 = fig.add_subplot(sub_c[0, 1])
    ax_d = fig.add_subplot(grid[1, 1])

    ax_a.text(-0.12, 1.05, "A", transform=ax_a.transAxes, fontsize=13, fontweight="bold")
    ax_b.text(-0.12, 1.05, "B", transform=ax_b.transAxes, fontsize=13, fontweight="bold")
    ax_c1.text(-0.12, 1.05, "C", transform=ax_c1.transAxes, fontsize=13, fontweight="bold")
    ax_d.text(-0.12, 1.05, "D", transform=ax_d.transAxes, fontsize=13, fontweight="bold")

    metric_specs = [
        ("L1_mean", "L1_mean ↓"),
        ("JS_mean", "JS_mean ↓"),
        ("corr_mean", "corr_mean ↑"),
        ("missing_leak_rate", "missing_leak_rate ↓"),
        ("extra_pred_fraction", "extra_pred_fraction ↓"),
    ]
    labels = []
    base_vals = []
    route_vals = []
    for key, lab in metric_specs:
        if key in base_metrics and key in route_metrics:
            labels.append(lab)
            base_vals.append(base_metrics[key])
            route_vals.append(route_metrics[key])
    if labels:
        x = np.arange(len(labels))
        ax_a.vlines(x, base_vals, route_vals, color="#888888", linewidth=1.6, zorder=1)
        ax_a.scatter(x - 0.06, base_vals, s=50, color="#4e79a7", label="Baseline", zorder=2)
        ax_a.scatter(x + 0.06, route_vals, s=50, color="#f28e2b", label="Route2", zorder=2)
        for i, (b, r) in enumerate(zip(base_vals, route_vals)):
            delta = r - b
            ax_a.text(i + 0.10, max(b, r), f"{delta:+.3f}", fontsize=8, va="bottom")
        ax_a.set_xticks(x)
        ax_a.set_xticklabels(labels, rotation=25, ha="right", fontsize=9)
        ax_a.set_title("Global performance (truth-based)", fontsize=STYLE["subplot"])
        ax_a.grid(axis="y", linestyle="--", alpha=0.25)
        ax_a.legend(frameon=False, fontsize=9, loc="best")
    else:
        ax_a.text(0.5, 0.5, "metrics not available", ha="center", va="center", fontsize=10)
        ax_a.set_axis_off()

    truth_ok = truth_df is not None and not truth_df.empty
    if truth_ok:
        base_idx = _idx(base_df)
        route_idx = _idx(route_df)
        truth_idx = _idx(truth_df)
        common = base_idx.index.intersection(route_idx.index).intersection(truth_idx.index)
        base_idx = base_idx.loc[common]
        route_idx = route_idx.loc[common]
        truth_idx = truth_idx.loc[common]
        if common.size > 0:
            all_types = sorted(set(base_idx.columns) | set(route_idx.columns) | set(truth_idx.columns))
            base_idx = base_idx.reindex(columns=all_types, fill_value=0.0)
            route_idx = route_idx.reindex(columns=all_types, fill_value=0.0)
            truth_idx = truth_idx.reindex(columns=all_types, fill_value=0.0)
            err_base = _compute_spot_l1(base_idx, truth_idx, all_types)
            err_route = _compute_spot_l1(route_idx, truth_idx, all_types)
            improvement = err_base - err_route
            improvement = improvement[np.isfinite(improvement)]
            if improvement.size > 0:
                missing_types = _pick_missing_types(truth_idx, eps=1e-6)
                affected = None
                tau = 0.01
                if missing_types:
                    missing = missing_types[0]
                    if missing in base_idx.columns:
                        affected = base_idx[missing].to_numpy(dtype=float) > tau
                if affected is None:
                    thresh = float(np.nanquantile(err_base, 0.9)) if err_base.size else 0.0
                    affected = err_base >= thresh

                x_all, y_all = _ecdf(improvement)
                x_aff, y_aff = _ecdf(improvement[affected]) if affected is not None else (np.array([]), np.array([]))
                x_min = float(np.min(x_all)) if x_all.size else 0.0
                x_max = float(np.max(x_all)) if x_all.size else 0.0
                x_min = min(x_min, 0.0)
                x_max = max(x_max, 0.0)
                pad = max(0.02, 0.05 * (x_max - x_min))
                ax_b.axvspan(0.0, x_max, color="#e8f5e9", alpha=0.25, zorder=0)
                ax_b.plot(x_all, y_all, color="#4e79a7", linewidth=1.6, label=f"All spots (n={len(improvement)})")
                if x_aff.size:
                    ax_b.plot(x_aff, y_aff, color="#f28e2b", linewidth=2.6, label=f"Affected (n={int(np.sum(affected))})")
                ax_b.axvline(0, color="#333333", linewidth=1.2, linestyle="--")
                frac_imp_all = float(np.mean(improvement > 0)) if improvement.size else 0.0
                frac_imp_aff = float(np.mean(improvement[affected] > 0)) if np.any(affected) else np.nan
                frac_non_worse = float(np.mean(improvement >= 0)) if improvement.size else 0.0
                med = float(np.median(improvement)) if improvement.size else 0.0
                mean = float(np.mean(improvement)) if improvement.size else 0.0
                y0_all = float(np.mean(improvement <= 0)) if improvement.size else 0.0
                y0_aff = float(np.mean(improvement[affected] <= 0)) if np.any(affected) else np.nan
                ax_b.set_title("Spot-level improvement ECDF", fontsize=STYLE["subplot"])
                ax_b.set_xlabel("Improvement = err(baseline) - err(route2) (↑ better)")
                ax_b.set_ylabel("ECDF")
                ax_b.grid(alpha=0.25, linestyle="--")
                ax_b.set_xlim(x_min - pad, x_max + pad)
                if np.isfinite(y0_all):
                    ax_b.scatter([0.0], [y0_all], color="#4e79a7", s=24, zorder=3)
                if np.isfinite(y0_aff):
                    ax_b.scatter([0.0], [y0_aff], color="#f28e2b", s=24, zorder=3)
                ax_b.legend(frameon=False, fontsize=9, loc="upper left")
                ax_b.text(
                    0.02,
                    0.18,
                    "Affected: baseline missing-type pred > 0.01 (truth=0)",
                    transform=ax_b.transAxes,
                    ha="left",
                    va="bottom",
                    fontsize=8,
                )
                ax_b.text(
                    0.98,
                    0.98,
                    f"All improved = {frac_imp_all:.2f}\n"
                    f"Affected improved = {frac_imp_aff:.2f}\n"
                    f"non-worse = {frac_non_worse:.2f}\n"
                    f"median = {med:.3f} | mean = {mean:.3f}",
                    transform=ax_b.transAxes,
                    ha="right",
                    va="top",
                    fontsize=9,
                )
            else:
                ax_b.text(0.5, 0.5, "no spot-level errors", ha="center", va="center", fontsize=10)
                ax_b.set_axis_off()
        else:
            ax_b.text(0.5, 0.5, "no shared spots with truth", ha="center", va="center", fontsize=10)
            ax_b.set_axis_off()
    else:
        ax_b.text(0.5, 0.5, "truth not available", ha="center", va="center", fontsize=10)
        ax_b.set_axis_off()

    if truth_ok:
        truth_idx = _idx(truth_df)
        base_idx = _idx(base_df)
        route_idx = _idx(route_df)
        common = base_idx.index.intersection(route_idx.index).intersection(truth_idx.index)
        base_idx = base_idx.loc[common]
        route_idx = route_idx.loc[common]
        truth_idx = truth_idx.loc[common]
        if common.size > 0:
            type_cols = sorted(set(truth_idx.columns) | set(base_idx.columns) | set(route_idx.columns))
            base_idx = base_idx.reindex(columns=type_cols, fill_value=0.0)
            route_idx = route_idx.reindex(columns=type_cols, fill_value=0.0)
            truth_idx = truth_idx.reindex(columns=type_cols, fill_value=0.0)
            mode, key_type, note = _choose_panel_c_target(truth_idx, base_idx, route_idx, eps=1e-6)
            if mode == "missing_leak" and key_type:
                tau = 0.01
                base_vals = base_idx[key_type].to_numpy(dtype=float)
                route_vals = route_idx[key_type].to_numpy(dtype=float)
                base_plot = np.clip(base_vals, 0.0, 1.0)
                route_plot = np.clip(route_vals, 0.0, 1.0)
                base_plot[np.isclose(base_plot, 0.0, atol=1e-12)] = 0.0
                route_plot[np.isclose(route_plot, 0.0, atol=1e-12)] = 0.0
                p95_base = float(np.quantile(base_plot, 0.95)) if base_plot.size else 0.0
                x_max = max(0.05, tau * 5.0, p95_base)
                x_max = min(0.2, x_max)
                xb, yb = _ecdf(base_plot)
                xr, yr = _ecdf(route_plot)
                ax_c1.plot(xb, yb, color="#4e79a7", linewidth=2.0)
                if route_plot.size and float(np.max(route_plot)) <= 1e-12:
                    ax_c1.axvspan(0.0, tau, color="#f28e2b", alpha=0.08, zorder=1)
                    ax_c2.axvspan(0.0, tau, color="#f28e2b", alpha=0.15, zorder=1)
                    ax_c2.vlines(0.0, 0.0, 1.0, color="#f28e2b", linewidth=3.2, zorder=4)
                    ax_c2.plot([0.0, 0.0], [0.0, 1.0], transform=ax_c2.transAxes, color="#f28e2b", linewidth=3.2, zorder=5)
                    ax_c2.scatter([0.0], [1.0], s=28, color="#f28e2b", zorder=6)
                    ax_c2.spines["left"].set_visible(False)
                    ax_c2.tick_params(axis="y", left=False)
                else:
                    ax_c1.axvspan(0.0, tau, color="#f28e2b", alpha=0.08, zorder=1)
                    ax_c2.axvspan(0.0, tau, color="#f28e2b", alpha=0.15, zorder=1)
                    ax_c2.plot(xr, yr, color="#f28e2b", linewidth=2.0)
                ax_c1.axvline(0.0, color="#333333", linewidth=1.0, linestyle="--", alpha=0.18, zorder=1)
                ax_c2.axvline(0.0, color="#333333", linewidth=1.0, linestyle="--", alpha=0.18, zorder=1)
                ax_c1.axvline(tau, color="#333333", linewidth=1.0, linestyle="--", alpha=0.6)
                ax_c2.axvline(tau, color="#333333", linewidth=1.0, linestyle="--", alpha=0.6)
                ax_c1.set_title(f"Baseline (leak present)\n{_short_type(key_type, alias_map)}", fontsize=STYLE["subplot"])
                ax_c2.set_title(f"Route2 (leak removed)\n{_short_type(key_type, alias_map)}", fontsize=STYLE["subplot"])
                for ax in (ax_c1, ax_c2):
                    ax.set_xlabel(f"pred fraction ({_short_type(key_type, alias_map)})")
                    ax.set_ylabel("ECDF")
                    ax.grid(alpha=0.25, linestyle="--")

                def _stats(vals: np.ndarray) -> Dict[str, float]:
                    return {
                        "sum": float(np.sum(vals)),
                        "p99": float(np.quantile(vals, 0.99)) if vals.size else 0.0,
                        "frac": float(np.mean(vals > tau)) if vals.size else 0.0,
                    }

                sb = _stats(base_vals)
                sr = _stats(route_vals)
                ax_c1.set_xlim(0.0, x_max)
                ax_c2.set_xlim(0.0, x_max)
                ax_c1.text(
                    0.98,
                    0.02,
                    f"sum={sb['sum']:.2f}\n"
                    f"p99={sb['p99']:.3f}\n"
                    f"frac>{tau}={sb['frac']:.2f}",
                    transform=ax_c1.transAxes,
                    ha="right",
                    va="bottom",
                    fontsize=9,
                )
                ax_c2.text(
                    0.98,
                    0.02,
                    f"sum={sr['sum']:.2f}\n"
                    f"p99={sr['p99']:.3f}\n"
                    f"frac>{tau}={sr['frac']:.2f}",
                    transform=ax_c2.transAxes,
                    ha="right",
                    va="bottom",
                    fontsize=9,
                )
                ax_c1.set_ylabel("ECDF")
                ax_c2.set_ylabel("")
            elif key_type:
                all_types = sorted(set(base_idx.columns) | set(route_idx.columns) | set(truth_idx.columns))
                base_idx = base_idx.reindex(columns=all_types, fill_value=0.0)
                route_idx = route_idx.reindex(columns=all_types, fill_value=0.0)
                truth_idx = truth_idx.reindex(columns=all_types, fill_value=0.0)
                truth_vals = truth_idx[key_type].to_numpy(dtype=float)
                base_vals = base_idx[key_type].to_numpy(dtype=float)
                route_vals = route_idx[key_type].to_numpy(dtype=float)
                max_val = float(np.nanmax([truth_vals, base_vals, route_vals])) if truth_vals.size else 1.0
                max_val = max(max_val, 1e-6)
                for ax, pred, title_prefix in [
                    (ax_c1, base_vals, "Baseline"),
                    (ax_c2, route_vals, "Route2"),
                ]:
                    ax.scatter(truth_vals, pred, s=8, alpha=0.6, color="#4e79a7")
                    ax.plot([0, max_val], [0, max_val], color="#333333", linewidth=1.0)
                    ax.set_xlim(0, max_val)
                    ax.set_ylim(0, max_val)
                    ax.set_xlabel("Truth fraction")
                    ax.set_ylabel("Pred fraction")
                    if truth_vals.size and np.std(truth_vals) > 0 and np.std(pred) > 0:
                        corr = float(np.corrcoef(truth_vals, pred)[0, 1])
                    else:
                        corr = np.nan
                    rmse = float(np.sqrt(np.mean((pred - truth_vals) ** 2))) if truth_vals.size else np.nan
                    corr_str = f"{corr:.2f}" if np.isfinite(corr) else "NA"
                    rmse_str = f"{rmse:.3f}" if np.isfinite(rmse) else "NA"
                    ax.set_title(f"{title_prefix}: {_short_type(key_type, alias_map)}\ncorr={corr_str}, rmse={rmse_str}", fontsize=STYLE["subplot"])
                if note:
                    ax_c1.text(0.02, 0.95, note, transform=ax_c1.transAxes, ha="left", va="top", fontsize=8)
            else:
                for ax in (ax_c1, ax_c2):
                    ax.text(0.5, 0.5, "no key type", ha="center", va="center", fontsize=10)
                    ax.set_axis_off()
        else:
            for ax in (ax_c1, ax_c2):
                ax.text(0.5, 0.5, "no shared spots with truth", ha="center", va="center", fontsize=10)
                ax.set_axis_off()
    else:
        for ax in (ax_c1, ax_c2):
            ax.text(0.5, 0.5, "truth not available", ha="center", va="center", fontsize=10)
            ax.set_axis_off()

    _plot_lollipop_axis(ax_d, delta_mass, alias_map, label_fmt, max_types=6)

    fig.suptitle("Figure 1 | Performance contrast (truth-aware)", fontsize=STYLE["title"], y=0.98)
    fig.subplots_adjust(top=0.92, left=0.06, right=0.98, bottom=0.07)
    _save_fig(fig, out_path, dpi, tight=False)


def plot_relabel_heatmap(
    out_path: Path,
    relabel_path: Path,
    max_types: int,
    alias_map: Dict[str, str] | None,
    normalize: bool,
    drop_empty_types: bool,
    dpi: int,
) -> float | None:
    if not _maybe_plot() or not relabel_path.exists():
        return None
    df = pd.read_csv(relabel_path)
    if "orig_type" not in df.columns or "plugin_type" not in df.columns:
        return None
    ct = pd.crosstab(df["orig_type"], df["plugin_type"])
    if ct.empty:
        return None
    counts = ct.sum(axis=1).add(ct.sum(axis=0), fill_value=0.0)
    types = counts.sort_values(ascending=False).index.tolist()
    if max_types and len(types) > max_types:
        types = types[:max_types]
    for t in ["Unknown_sc_only", "Unknown"]:
        if t in counts.index and t not in types:
            types.append(t)
            if len(types) > max_types:
                types = types[:max_types]
            break
    ct = ct.reindex(index=types, columns=types, fill_value=0)
    empty_types: list[str] = []
    row_sum = ct.sum(axis=1)
    col_sum = ct.sum(axis=0)
    if drop_empty_types:
        keep = (row_sum > 0) & (col_sum > 0)
        empty_types = [t for t in ct.index if not keep.get(t, False)]
        ct = ct.loc[keep, keep]
        row_sum = ct.sum(axis=1)
    else:
        empty_types = [t for t in ct.index if row_sum.get(t, 0) == 0 and col_sum.get(t, 0) == 0]
    if normalize:
        row_sum = ct.sum(axis=1)
        valid = row_sum[row_sum > 0].index
        mat = ct.astype(float).copy()
        denom = row_sum.replace(0, np.nan)
        mat = (mat.div(denom, axis=0) * 100.0).fillna(0.0)
        vmin, vmax = 0.0, 100.0
        cmap = "Blues"
        title = "Type relabel consistency\nmax off-diagonal = {val:.1f}% (OK if < 5%)"
        show_diag = True
        if valid.size == 0:
            offdiag_max = 0.0
        else:
            frac = ct.loc[valid, valid].astype(float)
            frac = frac.div(frac.sum(axis=1), axis=0).fillna(0.0)
            mask = ~np.eye(frac.shape[0], dtype=bool)
            offdiag_vals = frac.to_numpy()[mask]
            offdiag_max = float(np.nanmax(offdiag_vals)) * 100.0 if offdiag_vals.size else 0.0
        if not np.isfinite(offdiag_max):
            offdiag_max = 0.0
        title = title.format(val=offdiag_max)
    else:
        mat = ct.astype(float)
        vmin, vmax = None, None
        cmap = "viridis"
        title = "Type relabel counts (supplement)"
        show_diag = False
        offdiag_max = None

    fig, ax = plt.subplots(figsize=FIGSIZE_SINGLE)
    im = ax.imshow(mat.values, cmap=cmap, vmin=vmin, vmax=vmax)
    labels = [_wrap_label(_short_type(t, alias_map)) for t in mat.index.tolist()]
    rotation = 90 if len(labels) > 10 else 45
    tick_size = 8 if rotation == 90 else STYLE["tick"]
    ax.set_xticks(np.arange(len(labels)))
    ax.set_yticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=rotation, ha="right", fontsize=tick_size)
    ax.set_yticklabels(labels, fontsize=STYLE["tick"])
    ax.set_title(title, fontsize=STYLE["subplot"])
    if empty_types and not drop_empty_types:
        ax.text(0.01, -0.18, "empty types: " + ", ".join([_short_type(t, alias_map) for t in empty_types]),
                transform=ax.transAxes, fontsize=8, ha="left")
    fig.subplots_adjust(bottom=0.32 if normalize else 0.45)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat.iat[i, j]
            if normalize:
                if i == j:
                    ax.text(j, i, f"{val:.0f}", ha="center", va="center", fontsize=8, color="#111111")
                elif val >= 1.0:
                    ax.text(j, i, f"{val:.0f}", ha="center", va="center", fontsize=7, color="#111111")
            else:
                if val > 0 and val >= mat.values.max() * 0.05:
                    ax.text(j, i, f"{int(val)}", ha="center", va="center", fontsize=7, color="#111111")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    _save_fig(fig, out_path, dpi, tight=False)
    return offdiag_max


def _knn_smooth(values: np.ndarray, coords: np.ndarray, k: int, alpha: float) -> np.ndarray:
    if k <= 0 or alpha <= 0:
        return values
    n = coords.shape[0]
    diff_x = coords[:, 0][:, None] - coords[:, 0][None, :]
    diff_y = coords[:, 1][:, None] - coords[:, 1][None, :]
    dist = diff_x * diff_x + diff_y * diff_y
    idx = np.argsort(dist, axis=1)[:, 1 : k + 1]
    neigh_mean = np.nanmean(values[idx], axis=1)
    return (1 - alpha) * values + alpha * neigh_mean


def plot_delta_maps(
    out_path: Path,
    delta_df: pd.DataFrame,
    coords_df: pd.DataFrame,
    x_col: str,
    y_col: str,
    types: List[str],
    delta_mass: pd.Series,
    global_vmax: float,
    clip_quantile: float,
    smooth_k: int,
    smooth_alpha: float,
    alias_map: Dict[str, str] | None,
    dpi: int,
) -> None:
    if not _maybe_plot() or not types:
        return
    merged = coords_df.merge(delta_df, on="spot_id", how="inner")
    if merged.empty:
        return
    coords = merged[[x_col, y_col]].to_numpy(dtype=float)
    n = len(types)
    ncols = 2 if n > 1 else 1
    nrows = math.ceil(n / ncols)
    fig = plt.figure(figsize=(6.5 * ncols + 0.8, 4.5 * nrows))
    grid = fig.add_gridspec(
        nrows=nrows,
        ncols=ncols + 1,
        width_ratios=[1] * ncols + [0.06],
        wspace=0.05,
        hspace=0.15,
    )
    axes = []
    for r in range(nrows):
        for c in range(ncols):
            axes.append(fig.add_subplot(grid[r, c]))
    cax = fig.add_subplot(grid[:, -1])
    sc_last = None
    for idx, t in enumerate(types):
        ax = axes[idx]
        vals = merged[t].to_numpy(dtype=float)
        vals = _knn_smooth(vals, coords, smooth_k, smooth_alpha)
        if not np.isfinite(vals).any():
            continue
        vmax = max(global_vmax, 1e-6)
        norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
        sc = ax.scatter(
            merged[x_col],
            merged[y_col],
            c=vals,
            cmap="coolwarm",
            s=6,
            norm=norm,
            edgecolor="none",
        )
        delta_val = float(delta_mass.get(t, 0.0))
        abs_vals = np.abs(vals)
        thresh = np.nanquantile(abs_vals, 0.9) if np.isfinite(abs_vals).any() else 0.0
        top_vals = vals[abs_vals >= thresh] if thresh > 0 else vals
        top_mean = float(np.nanmean(top_vals)) if top_vals.size else 0.0
        title = f"Δ {_short_type(t, alias_map)} (route2 - baseline)\nΔmass={delta_val:.1f} | Δtop10%={top_mean:.2f}"
        ax.set_title(title, fontsize=STYLE["subplot"])
        ax.set_xticks([])
        ax.set_yticks([])
        sc_last = sc
    for j in range(len(types), len(axes)):
        axes[j].axis("off")
    if sc_last is not None:
        cbar = fig.colorbar(sc_last, cax=cax)
        cbar.set_label(f"delta (clipped at p{int(clip_quantile * 100)})")
        cax.tick_params(labelsize=STYLE["tick"], pad=2)
    else:
        cax.axis("off")
    _save_fig(fig, out_path, dpi)


def plot_marker_maps(
    out_path: Path,
    plotdata_df: pd.DataFrame,
    markers: List[str],
    preferred: List[str],
    x_col: str,
    y_col: str,
    max_markers: int,
    clip_quantile: float,
    log1p: bool,
    dpi: int,
) -> Tuple[List[str], List[str]]:
    if not _maybe_plot():
        return [], []
    if plotdata_df.empty or x_col is None or y_col is None:
        return [], []
    if x_col not in plotdata_df.columns or y_col not in plotdata_df.columns:
        return [], []
    available = [m for m in markers if m in plotdata_df.columns]
    if not available:
        return [], preferred
    preferred = [m for m in preferred if m]
    used = [m for m in preferred if m in available]
    if len(used) < max_markers:
        for m in available:
            if m not in used:
                used.append(m)
            if len(used) >= max_markers:
                break
    missing = [m for m in preferred if m not in available]
    n = len(used)
    ncols = 2 if n > 1 else 1
    nrows = math.ceil(n / ncols)
    fig = plt.figure(figsize=(6.5 * ncols + 0.6, 4.5 * nrows))
    grid = fig.add_gridspec(nrows=nrows, ncols=ncols + 1, width_ratios=[1] * ncols + [0.05], wspace=0.12, hspace=0.25)
    axes = []
    for r in range(nrows):
        for c in range(ncols):
            axes.append(fig.add_subplot(grid[r, c]))
    cax = fig.add_subplot(grid[:, -1])
    sc_last = None
    for idx, m in enumerate(used):
        ax = axes[idx]
        vals = plotdata_df[m].to_numpy(dtype=float)
        raw = np.clip(vals, 0, None)
        if log1p:
            vals = np.log1p(raw)
        vmax = float(np.nanquantile(vals, clip_quantile)) if np.isfinite(vals).any() else 1.0
        vmax = max(vmax, 1e-6)
        vals = np.clip(vals, 0, vmax)
        scaled = vals / vmax if vmax > 0 else vals
        detect_rate = float(np.mean(raw > 0)) if raw.size else 0.0
        sc = ax.scatter(plotdata_df[x_col], plotdata_df[y_col], c=scaled, cmap="viridis", s=6, vmin=0, vmax=1, edgecolor="none")
        ax.set_title(f"{m} (detect_rate={detect_rate:.2f})", fontsize=STYLE["subplot"])
        ax.set_xticks([])
        ax.set_yticks([])
        sc_last = sc
    for j in range(len(used), len(axes)):
        axes[j].axis("off")
    if sc_last is not None:
        cbar = fig.colorbar(sc_last, cax=cax)
        cbar.set_label("log1p(counts), clipped@p99, min-max per gene")
    else:
        cax.axis("off")
    _save_fig(fig, out_path, dpi)
    return used, missing


def plot_compartment_maps(
    out_path: Path,
    comp_df: pd.DataFrame,
    coords_df: pd.DataFrame,
    x_col: str,
    y_col: str,
    compartments: List[str],
    comp_metrics: Dict[str, Dict[str, float]] | None,
    clip_quantile: float,
    smooth_k: int,
    smooth_alpha: float,
    dpi: int,
) -> None:
    if not _maybe_plot() or not compartments:
        return
    merged = coords_df.merge(comp_df, on="spot_id", how="inner")
    if merged.empty:
        return
    coords = merged[[x_col, y_col]].to_numpy(dtype=float)
    nrows = len(compartments)
    fig = plt.figure(figsize=(15.5, 4 * nrows))
    grid = fig.add_gridspec(
        nrows=nrows,
        ncols=3,
        width_ratios=[1, 1, 1],
        wspace=0.08,
        hspace=0.18,
        left=0.05,
        right=0.86,
        top=0.96,
        bottom=0.06,
    )
    row_axes = []
    comp_records = []
    for i, c in enumerate(compartments):
        base_col = f"{c}__baseline"
        route_col = f"{c}__route2"
        delta_col = f"{c}__delta"
        if base_col not in merged.columns or route_col not in merged.columns or delta_col not in merged.columns:
            continue
        base_vals = _knn_smooth(merged[base_col].to_numpy(dtype=float), coords, smooth_k, smooth_alpha)
        route_vals = _knn_smooth(merged[route_col].to_numpy(dtype=float), coords, smooth_k, smooth_alpha)
        delta_vals = _knn_smooth(merged[delta_col].to_numpy(dtype=float), coords, smooth_k, smooth_alpha)
        comp_records.append((c, base_vals, route_vals, delta_vals))

    if not comp_records:
        return

    all_frac = np.concatenate([np.concatenate([b, r]) for _, b, r, _ in comp_records])
    frac_vmax = float(np.nanquantile(all_frac, clip_quantile)) if np.isfinite(all_frac).any() else 1.0
    frac_vmax = max(frac_vmax, 1e-6)
    all_delta = np.concatenate([np.abs(d) for _, _, _, d in comp_records])
    delta_vmax = float(np.nanquantile(all_delta, clip_quantile)) if np.isfinite(all_delta).any() else 1.0
    delta_vmax = max(delta_vmax, float(np.nanmax(all_delta)), 1e-6)
    norm_delta = TwoSlopeNorm(vmin=-delta_vmax, vcenter=0.0, vmax=delta_vmax)

    for i, (c, base_vals, route_vals, delta_vals) in enumerate(comp_records):
        ax_base = fig.add_subplot(grid[i, 0])
        ax_route = fig.add_subplot(grid[i, 1])
        ax_delta = fig.add_subplot(grid[i, 2])

        for ax in (ax_base, ax_route, ax_delta):
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_facecolor("#f2f2f2")

        sc_base = ax_base.scatter(
            merged[x_col],
            merged[y_col],
            c=base_vals,
            cmap="viridis",
            s=6,
            vmin=0,
            vmax=frac_vmax,
            edgecolor="none",
        )
        ax_base.set_title(f"{c} baseline", fontsize=STYLE["subplot"])

        sc_route = ax_route.scatter(
            merged[x_col],
            merged[y_col],
            c=route_vals,
            cmap="viridis",
            s=6,
            vmin=0,
            vmax=frac_vmax,
            edgecolor="none",
        )
        ax_route.set_title(f"{c} route2", fontsize=STYLE["subplot"])

        sc_delta = ax_delta.scatter(
            merged[x_col],
            merged[y_col],
            c=delta_vals,
            cmap="coolwarm",
            s=6,
            norm=norm_delta,
            edgecolor="none",
        )
        title = f"{c} delta"
        if comp_metrics and c in comp_metrics:
            metrics = comp_metrics[c]
            corr = metrics.get("corr_delta")
            rmse = metrics.get("rmse_delta")
            label = []
            if corr is not None and np.isfinite(corr):
                label.append(f"corr(base,Δ)={corr:.2f}")
            if rmse is not None and np.isfinite(rmse):
                label.append(f"rmse(Δ)={rmse:.3f}")
            if label:
                title = f"{title}\n" + " | ".join(label)
        ax_delta.set_title(title, fontsize=STYLE["subplot"])
        row_axes.append(ax_base)

    m_frac = plt.cm.ScalarMappable(norm=plt.Normalize(vmin=0, vmax=frac_vmax), cmap="viridis")
    m_frac.set_array([])
    m_delta = plt.cm.ScalarMappable(norm=norm_delta, cmap="coolwarm")
    m_delta.set_array([])

    right_start = 0.88
    cbar_width = 0.012
    gap = 0.006
    for ax_ref in row_axes:
        bbox = ax_ref.get_position()
        cax_frac = fig.add_axes([right_start, bbox.y0, cbar_width, bbox.height])
        cax_delta = fig.add_axes([right_start + cbar_width + gap, bbox.y0, cbar_width, bbox.height])
        cb1 = fig.colorbar(m_frac, cax=cax_frac)
        cb1.set_label(f"fraction (p{int(clip_quantile * 100)})", fontsize=9, labelpad=6)
        cb2 = fig.colorbar(m_delta, cax=cax_delta)
        cb2.set_label(f"delta (p{int(clip_quantile * 100)})", fontsize=9, labelpad=6)
        cax_frac.tick_params(labelsize=8, pad=2)
        cax_delta.tick_params(labelsize=8, pad=2)

    _save_fig(fig, out_path, dpi, tight=False)


def plot_rescue_map(
    out_path: Path,
    rescue_df: pd.DataFrame,
    coords_df: pd.DataFrame,
    x_col: str,
    y_col: str,
    dpi: int,
) -> None:
    if not _maybe_plot():
        return
    merged = coords_df.merge(rescue_df, on="spot_id", how="inner")
    if merged.empty or "rescue_score" not in merged.columns:
        return
    fig, ax = plt.subplots(figsize=(6, 5))
    vals = merged["rescue_score"].to_numpy(dtype=float)
    sc = ax.scatter(merged[x_col], merged[y_col], c=vals, cmap="inferno", s=6)
    ax.set_title("Rescue score map")
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    _save_fig(fig, out_path, dpi)


def plot_pipeline_diagram(out_path: Path, dpi: int) -> None:
    if not _maybe_plot():
        return
    fig, ax = plt.subplots(figsize=FIGSIZE_WIDE)
    ax.axis("off")
    box = dict(boxstyle="round,pad=0.3", facecolor="#f3f3f3", edgecolor="#333333")
    ax.text(0.02, 0.55, "Stage1-3 inputs\nst_expression_*.csv\ncell_type_relabel.csv\n(+1 file)", bbox=box, fontsize=8)
    ax.text(0.30, 0.72, "Stage4 baseline\nstage4_summary.json\ncell_assignment.csv\nspot_type_fraction.csv (+2)", bbox=box, fontsize=8)
    ax.text(0.30, 0.32, "Stage4 route2\nstage4_summary.json\ncell_assignment.csv\nspot_type_fraction.csv (+2)", bbox=box, fontsize=8)
    ax.text(0.60, 0.52, "Stage6 audit\nstage6_audit_*.json\nstage6_plotdata_*.csv\nstage6_compare_*.json (+1)", bbox=box, fontsize=8)
    ax.text(0.82, 0.52, "Stage7 report\nfigures/*.png\noverall_summary.json\nStage7_report.html", bbox=box, fontsize=8)

    badge = dict(boxstyle="round,pad=0.2", facecolor="#ffffff", edgecolor="#666666")
    ax.text(0.30, 0.25, "type_alignment_ok", bbox=badge, fontsize=7)
    ax.text(0.30, 0.20, "filter_ledger_ok", bbox=badge, fontsize=7)

    ax.plot([0.20, 0.26], [0.60, 0.60], color="#333333", linewidth=1)
    ax.plot([0.26, 0.26], [0.60, 0.78], color="#333333", linewidth=1)
    ax.plot([0.26, 0.30], [0.78, 0.78], color="#333333", linewidth=1)
    ax.plot([0.26, 0.26], [0.60, 0.40], color="#333333", linewidth=1)
    ax.plot([0.26, 0.30], [0.40, 0.40], color="#333333", linewidth=1)

    ax.plot([0.48, 0.56], [0.78, 0.62], color="#333333", linewidth=1)
    ax.plot([0.48, 0.56], [0.40, 0.56], color="#333333", linewidth=1)
    ax.plot([0.72, 0.80], [0.56, 0.56], color="#333333", linewidth=1)

    ax.annotate("", xy=(0.30, 0.78), xytext=(0.26, 0.78), arrowprops=dict(arrowstyle="->"))
    ax.annotate("", xy=(0.30, 0.40), xytext=(0.26, 0.40), arrowprops=dict(arrowstyle="->"))
    ax.annotate("", xy=(0.60, 0.58), xytext=(0.56, 0.62), arrowprops=dict(arrowstyle="->"))
    ax.annotate("", xy=(0.60, 0.54), xytext=(0.56, 0.56), arrowprops=dict(arrowstyle="->"))
    ax.annotate("", xy=(0.82, 0.56), xytext=(0.80, 0.56), arrowprops=dict(arrowstyle="->"))

    ax.set_title("Stage7 workflow and evidence trace", fontsize=STYLE["title"])
    _save_fig(fig, out_path, dpi)


def plot_compartment_summary_table(
    out_path: Path,
    rows: List[Dict[str, Any]],
    dpi: int,
) -> None:
    if not _maybe_plot() or not rows:
        return
    df = pd.DataFrame(rows)
    df = df.sort_values("abs_delta", ascending=False)
    table_df = df.copy()
    table_df["corr"] = table_df["corr"].map(lambda v: f"{v:.2f}")
    table_df["rmse"] = table_df["rmse"].map(lambda v: f"{v:.3f}")
    table_df["delta_mass"] = table_df["delta_mass"].map(lambda v: f"{v:+.1f}")
    table_df = table_df.rename(columns={
        "corr": "corr(base,Δ)",
        "rmse": "rmse(Δ)",
        "delta_mass": "Δmass",
    })
    table_df = table_df[["compartment", "corr(base,Δ)", "rmse(Δ)", "Δmass"]]

    nrows = len(table_df) + 1
    fig, ax = plt.subplots(figsize=(6.5, 0.2 + 0.3 * nrows))
    ax.axis("off")
    tbl = ax.table(
        cellText=table_df.values,
        colLabels=table_df.columns,
        cellLoc="center",
        colLoc="center",
        loc="center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(10)
    tbl.scale(1.0, 1.1)
    for (row, col), cell in tbl.get_celld().items():
        if row == 0:
            cell.set_facecolor("#f0f0f0")
            cell.set_text_props(weight="bold")
    _save_fig(fig, out_path, dpi)


def build_report_html(
    out_path: Path,
    sample: str,
    baseline_id: str,
    route2_id: str,
    figures: Dict[str, str],
    highlights: List[str],
    abbrev_map: Dict[str, str],
    html_mode: str,
    fig_root: Path,
    comp_summary_rows: List[Dict[str, Any]] | None,
) -> None:
    def _mime_for(path: Path) -> str:
        ext = path.suffix.lower()
        if ext == ".pdf":
            return "application/pdf"
        if ext == ".svg":
            return "image/svg+xml"
        return "image/png"

    def _img_src(rel_path: str) -> str:
        if html_mode != "embed":
            return rel_path
        img_path = fig_root / rel_path
        if not img_path.exists():
            return rel_path
        data = base64.b64encode(img_path.read_bytes()).decode("ascii")
        mime = _mime_for(img_path)
        return f"data:{mime};base64,{data}"

    def _fig_tag(rel_path: str, width: int) -> str:
        ext = Path(rel_path).suffix.lower()
        src = _img_src(rel_path)
        if ext == ".pdf":
            return f"<embed src=\"{src}\" type=\"application/pdf\" width=\"{width}\" height=\"650\">"
        return f"<img src=\"{src}\" width=\"{width}\">"

    sections = []
    if figures.get("fig1_main"):
        sections.append(f"<h2>Figure 1 (main)</h2>{_fig_tag(figures['fig1_main'], 1000)}")
    if figures.get("fig1"):
        sections.append(f"<h2>Figure 1</h2>{_fig_tag(figures['fig1'], 900)}")
    if figures.get("fig1b"):
        sections.append(_fig_tag(figures["fig1b"], 900))
    if figures.get("fig1c"):
        sections.append(_fig_tag(figures["fig1c"], 900))
    if figures.get("fig1d"):
        sections.append(_fig_tag(figures["fig1d"], 900))
    if figures.get("fig2"):
        sections.append(f"<h2>Figure 2</h2>{_fig_tag(figures['fig2'], 900)}")
    if figures.get("fig2b"):
        sections.append(_fig_tag(figures["fig2b"], 900))
    if figures.get("fig2c"):
        sections.append(_fig_tag(figures["fig2c"], 900))
    if figures.get("fig2d"):
        sections.append(_fig_tag(figures["fig2d"], 900))
    if figures.get("fig2e"):
        sections.append(_fig_tag(figures["fig2e"], 900))
    if figures.get("fig3"):
        sections.append(f"<h2>Figure 3</h2>{_fig_tag(figures['fig3'], 900)}")

    hl_list = "".join([f"<li>{h}</li>" for h in highlights]) if highlights else "<li>no highlights</li>"
    abbrev_list = "".join([f"<li>{v} = {k}</li>" for k, v in abbrev_map.items()]) if abbrev_map else "<li>none</li>"
    note = ""
    if html_mode == "folder":
        note = "<p><b>Note:</b> keep this HTML in the same folder as the <code>figures/</code> directory.</p>"
    comp_table = ""
    if comp_summary_rows:
        rows = sorted(comp_summary_rows, key=lambda r: abs(r.get("delta_mass", 0.0)), reverse=True)
        header = "<tr><th>compartment</th><th>corr(base,Δ)</th><th>rmse(Δ)</th><th>Δmass</th></tr>"
        body = []
        for r in rows:
            corr = r.get("corr")
            rmse = r.get("rmse")
            dm = r.get("delta_mass")
            corr_str = f"{corr:.2f}" if corr is not None and np.isfinite(corr) else "NA"
            rmse_str = f"{rmse:.3f}" if rmse is not None and np.isfinite(rmse) else "NA"
            dm_str = f"{dm:+.1f}" if dm is not None else "NA"
            body.append(f"<tr><td>{r.get('compartment','')}</td><td>{corr_str}</td><td>{rmse_str}</td><td>{dm_str}</td></tr>")
        comp_table = "<h2>Compartment Summary</h2><table border=\"1\" cellspacing=\"0\" cellpadding=\"4\">" + header + "".join(body) + "</table>"

    unknown_note = ""
    if "Unknown_sc_only" in abbrev_map:
        unknown_note = "<p><b>Note:</b> Unknown(sc-only) is reserved for filter ledger; it is not a relabel target.</p>"

    html = f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8"/>
  <title>Stage7 Report - {sample}</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 24px; }}
    img {{ display: block; margin: 16px 0; }}
  </style>
</head>
<body>
  <h1>Stage7 Report</h1>
  <p>Sample: {sample}</p>
  <p>Baseline: {baseline_id} | Route2: {route2_id}</p>
  {note}
  <h2>Highlights</h2>
  <ul>{hl_list}</ul>
  <h2>Type Abbreviations</h2>
  <ul>{abbrev_list}</ul>
  {unknown_note}
  <p><b>Marker maps:</b> detect_rate = fraction of spots with expression &gt; 0; values are per-gene minmax of log1p(counts).</p>
  {comp_table}
  {''.join(sections)}
</body>
</html>
"""
    out_path.write_text(html, encoding="utf-8")


def build_report_markdown(
    out_path: Path,
    sample: str,
    baseline_id: str,
    route2_id: str,
    figures: Dict[str, str],
    highlights: List[str],
    abbrev_map: Dict[str, str],
) -> None:
    lines = [
        "# Stage7_summary_CN",
        "",
        f"- 样本：{sample}",
        f"- Baseline：{baseline_id}",
        f"- Route2：{route2_id}",
        "",
        "## 关键结论",
    ]
    if highlights:
        lines.extend([f"- {h}" for h in highlights])
    else:
        lines.append("- 暂无摘要（需要完整输入后生成）")
    lines.append("")
    lines.append("## 图表索引")
    for key in ["fig1", "fig1b", "fig1c", "fig1d", "fig2", "fig2b", "fig2c", "fig2d", "fig2e", "fig3"]:
        if figures.get(key):
            lines.append(f"- {key}: {figures[key]}")
    if abbrev_map:
        lines.append("")
        lines.append("## Type Abbreviations")
        for full, short in abbrev_map.items():
            lines.append(f"- {short} = {full}")
    out_path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Stage7 report generator (real data)")
    p.add_argument("--sample", default="real_brca", help="sample id")
    p.add_argument("--baseline_id", default="baseline", help="baseline run id")
    p.add_argument("--route2_id", default="route2_v5_3_1", help="route2 run id")
    p.add_argument("--stage4_root", default=None, help="override stage4 root dir")
    p.add_argument("--stage6_root", default=None, help="override stage6 root dir")
    p.add_argument("--stage1_export", default=None, help="override stage1 export dir")
    p.add_argument("--dataset_config", default=None, help="override dataset config")
    p.add_argument("--out_dir", default=None, help="output directory")
    p.add_argument("--top_types", type=int, default=6, help="top delta types for summaries")
    p.add_argument("--delta_map_max_types", type=int, default=4, help="max delta-map types to render")
    p.add_argument("--delta_map_control_type", default="Epithelial cells", help="control type to force in delta maps")
    p.add_argument("--delta_clip_quantile", type=float, default=0.99, help="quantile clip for delta color scale")
    p.add_argument("--delta_smooth_k", type=int, default=6, help="kNN smoothing neighbors for delta maps")
    p.add_argument("--delta_smooth_alpha", type=float, default=0.4, help="smoothing strength for delta maps")
    p.add_argument("--lollipop_label_fmt", choices=["int", "1f", "2f"], default="int", help="label format for lollipop")
    p.add_argument("--lollipop_break_ratio", type=float, default=2.8, help="ratio threshold for broken axis")
    p.add_argument("--lollipop_break_min", type=float, default=80.0, help="min delta to trigger broken axis")
    p.add_argument("--max_markers", type=int, default=4, help="max markers to plot")
    p.add_argument("--marker_clip_quantile", type=float, default=0.99, help="quantile clip for marker maps")
    p.add_argument("--marker_log1p", default=True, action=argparse.BooleanOptionalAction, help="log1p marker maps")
    p.add_argument("--show_zero_flows", default=False, action=argparse.BooleanOptionalAction, help="show zero flows in flow plot")
    p.add_argument("--truth_fraction_path", default=None, help="override truth spot fraction csv path")
    p.add_argument("--drop_empty_types", default=True, action=argparse.BooleanOptionalAction, help="drop empty types in relabel heatmaps")
    p.add_argument("--html_mode", choices=["embed", "folder"], default="embed", help="html output mode")
    p.add_argument("--output_formats", default="png", help="comma-separated output formats (png,pdf,svg)")
    p.add_argument("--dpi", type=int, default=240, help="figure dpi")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    _apply_style()
    global OUTPUT_FORMATS
    formats = [f.strip().lower() for f in (args.output_formats or "").split(",") if f.strip()]
    if not formats:
        formats = ["pdf"]
    OUTPUT_FORMATS = formats
    project_root = detect_root()
    stage4_root = Path(args.stage4_root) if args.stage4_root else project_root / "result" / args.sample / "stage4_cytospace"
    stage6_root = Path(args.stage6_root) if args.stage6_root else project_root / "result" / args.sample / "stage6_real_audit"
    stage1_export = Path(args.stage1_export) if args.stage1_export else project_root / "data" / "processed" / args.sample / "stage1_preprocess" / "exported"
    dataset_cfg_path = Path(args.dataset_config) if args.dataset_config else project_root / "configs" / "datasets" / f"{args.sample}.yaml"

    out_dir = Path(args.out_dir) if args.out_dir else project_root / "result" / args.sample / "stage7_report"
    fig_dir = out_dir / "figures"
    table_dir = out_dir / "tables"
    summary_dir = out_dir / "summary_json"
    report_html_dir = out_dir / "report_html"
    report_md_dir = out_dir / "report_markdown"
    for p in [fig_dir, table_dir, summary_dir, report_html_dir, report_md_dir]:
        _ensure_dir(p)

    alias_map = load_alias_map(project_root / "configs" / "aliases.yaml")
    dataset_cfg = _load_yaml(dataset_cfg_path)

    base_run_dir = _resolve_stage4_run_dir(stage4_root, args.baseline_id)
    route_run_dir = _resolve_stage4_run_dir(stage4_root, args.route2_id)
    if base_run_dir is None or route_run_dir is None:
        print("[Stage7] missing Stage4 outputs.")
        return 1

    base_summary = _read_json(base_run_dir / "stage4_summary.json")
    route_summary = _read_json(route_run_dir / "stage4_summary.json")

    base_frac_path = None
    route_frac_path = None
    for name in ["fractional_abundances_by_spot.csv", "cell_type_assignments_by_spot.csv"]:
        if (base_run_dir / name).exists() and base_frac_path is None:
            base_frac_path = base_run_dir / name
        if (route_run_dir / name).exists() and route_frac_path is None:
            route_frac_path = route_run_dir / name
    if base_frac_path is None or route_frac_path is None:
        print("[Stage7] missing spot fraction tables.")
        return 1

    base_frac, _, _ = _load_fraction_table(base_frac_path)
    route_frac, _, _ = _load_fraction_table(route_frac_path)
    base_frac, _ = _canonicalize_fraction_df(base_frac, alias_map)
    route_frac, _ = _canonicalize_fraction_df(route_frac, alias_map)
    base_aligned, route_aligned, all_types = _align_fraction_tables(base_frac, route_frac)

    base_mass = base_aligned.sum(axis=0)
    route_mass = route_aligned.sum(axis=0)
    delta_mass = route_mass - base_mass

    top_types = (
        delta_mass.abs()
        .sort_values(ascending=False)
        .head(args.top_types)
        .index
        .tolist()
    )
    must_include = _collect_must_include(base_summary, route_summary, base_mass, route_mass)
    lollipop_types = top_types[:]
    for t in must_include:
        if t in all_types and t not in lollipop_types:
            lollipop_types.append(t)

    type_summary = pd.DataFrame({
        "type": all_types,
        "baseline_mass": [base_mass.get(t, 0.0) for t in all_types],
        "route2_mass": [route_mass.get(t, 0.0) for t in all_types],
        "delta_mass": [delta_mass.get(t, 0.0) for t in all_types],
    }).sort_values("delta_mass", ascending=False)
    type_summary.to_csv(table_dir / "type_composition.csv", index=False)

    coords_path = stage1_export / "st_coordinates.csv"
    coords_df, x_col, y_col = _load_coords(coords_path)

    delta_df = route_aligned.loc[route_aligned.index.intersection(base_aligned.index)] - \
        base_aligned.loc[route_aligned.index.intersection(base_aligned.index)]
    delta_df = delta_df.reset_index().rename(columns={"index": "spot_id"})

    stage6_compare = _read_json(stage6_root / "stage6_compare_baseline_vs_route2.json")
    base_audit = _read_json(stage6_root / args.baseline_id / f"stage6_audit_{args.baseline_id}.json")
    route_audit = _read_json(stage6_root / args.route2_id / f"stage6_audit_{args.route2_id}.json")

    figures: Dict[str, str] = {}
    abbrev_map = _build_type_abbrev(all_types, alias_map)
    missing_markers: List[str] = []
    comp_summary_rows: List[Dict[str, Any]] = []
    missing_compartments: List[str] = []
    if _maybe_plot():
        stage5_root = _resolve_stage5_root(project_root / "result" / args.sample)
        base_stage5_json = _read_json(_find_stage5_json(_resolve_stage5_run_dir(stage5_root, "baseline")))
        route_stage5_json = _read_json(_find_stage5_json(_resolve_stage5_run_dir(stage5_root, args.route2_id)))
        base_metrics = _extract_stage5_metrics(base_stage5_json)
        route_metrics = _extract_stage5_metrics(route_stage5_json)
        truth_path = _resolve_truth_fraction_path(stage1_export, args.truth_fraction_path)
        truth_frac = None
        if truth_path and truth_path.exists():
            truth_frac, _, _ = _load_fraction_table(truth_path)
            truth_frac, _ = _canonicalize_fraction_df(truth_frac, alias_map)
        plot_perf_contrast(
            fig_dir / "fig1_perf_contrast.png",
            base_frac,
            route_frac,
            truth_frac,
            base_metrics,
            route_metrics,
            delta_mass,
            alias_map,
            must_include,
            args.lollipop_label_fmt,
            dpi=args.dpi,
        )
        if (fig_dir / "fig1_perf_contrast.png").exists():
            figures["fig1_main"] = "figures/fig1_perf_contrast.png"

        comp_types = plot_type_composition(
            fig_dir / "fig1_type_composition.png",
            base_mass,
            route_mass,
            lollipop_types,
            alias_map,
            must_include,
            args.lollipop_label_fmt,
            args.lollipop_break_ratio,
            args.lollipop_break_min,
            dpi=args.dpi,
        )
        figures["fig1"] = "figures/fig1_type_composition.png"

        plot_filter_rescue(
            fig_dir / "fig1_filter_rescue.png",
            base_summary,
            route_summary,
            base_audit,
            route_audit,
            args.show_zero_flows,
            dpi=args.dpi,
        )
        figures["fig1b"] = "figures/fig1_filter_rescue.png"

        relabel_path = project_root / "data" / "processed" / args.sample / "stage3_typematch" / "cell_type_relabel.csv"
        plot_relabel_heatmap(
            fig_dir / "fig1_relabel_heatmap.png",
            relabel_path,
            max_types=10,
            alias_map=alias_map,
            normalize=True,
            drop_empty_types=args.drop_empty_types,
            dpi=args.dpi,
        )
        if (fig_dir / "fig1_relabel_heatmap.png").exists():
            figures["fig1c"] = "figures/fig1_relabel_heatmap.png"
        plot_relabel_heatmap(
            fig_dir / "fig1_relabel_heatmap_counts.png",
            relabel_path,
            max_types=10,
            alias_map=alias_map,
            normalize=False,
            drop_empty_types=args.drop_empty_types,
            dpi=args.dpi,
        )
        if (fig_dir / "fig1_relabel_heatmap_counts.png").exists():
            figures["fig1d"] = "figures/fig1_relabel_heatmap_counts.png"

        delta_map_types = []
        for t in top_types:
            if t in all_types and t not in delta_map_types:
                delta_map_types.append(t)
        for t in must_include:
            if t in all_types and t not in delta_map_types:
                delta_map_types.append(t)
        control_type = canonicalize_type_name(args.delta_map_control_type, alias_map) if alias_map else normalize_type_name(args.delta_map_control_type)
        if control_type in all_types and control_type not in delta_map_types:
            delta_map_types.append(control_type)
        if len(delta_map_types) > args.delta_map_max_types:
            keep = []
            for t in must_include:
                if t in delta_map_types and t not in keep:
                    keep.append(t)
            for t in top_types:
                if t in delta_map_types and t not in keep:
                    keep.append(t)
            delta_map_types = keep[: args.delta_map_max_types]

        if coords_df is not None and x_col and y_col and delta_map_types:
            abs_vals = np.abs(delta_df[delta_map_types].to_numpy(dtype=float))
            global_vmax = float(np.nanquantile(abs_vals, args.delta_clip_quantile)) if abs_vals.size else 1.0
            global_vmax = max(global_vmax, 1e-6)
            plot_delta_maps(
                fig_dir / "fig2_delta_maps.png",
                delta_df,
                coords_df,
                x_col,
                y_col,
                delta_map_types,
                delta_mass,
                global_vmax,
                args.delta_clip_quantile,
                smooth_k=args.delta_smooth_k,
                smooth_alpha=args.delta_smooth_alpha,
                alias_map=alias_map,
                dpi=args.dpi,
            )
            figures["fig2"] = "figures/fig2_delta_maps.png"

            rescue_types = (
                route_audit.get("ledger_integrity", {})
                .get("metrics", {})
                .get("rescue_ledger", {})
                .get("rescued_types", [])
            )
            rescue_types = [
                canonicalize_type_name(t, alias_map) if alias_map else normalize_type_name(t)
                for t in (rescue_types or [])
            ]
            rescue_types = [t for t in rescue_types if t in delta_df.columns]
            if rescue_types:
                rescue_df = pd.DataFrame({"spot_id": delta_df["spot_id"]})
                vals = delta_df[rescue_types].to_numpy(dtype=float)
                rescue_df["rescue_score"] = np.clip(vals, 0, None).sum(axis=1)
                plot_rescue_map(
                    fig_dir / "fig2_rescue_map.png",
                    rescue_df,
                    coords_df,
                    x_col,
                    y_col,
                    dpi=args.dpi,
                )
                if (fig_dir / "fig2_rescue_map.png").exists():
                    figures["fig2b"] = "figures/fig2_rescue_map.png"

        plotdata_path = stage6_root / args.route2_id / f"stage6_plotdata_{args.route2_id}.csv"
        if plotdata_path.exists():
            plot_df = pd.read_csv(plotdata_path)
            if "spot_id" in plot_df.columns:
                plot_df["spot_id"] = plot_df["spot_id"].astype(str).map(_clean_spot_id)
            if coords_df is not None and x_col and y_col:
                if x_col not in plot_df.columns or y_col not in plot_df.columns:
                    plot_df = plot_df.merge(
                        coords_df[["spot_id", x_col, y_col]],
                        on="spot_id",
                        how="left",
                    )
                    if x_col not in plot_df.columns:
                        alt = f"{x_col}_y"
                        if alt in plot_df.columns:
                            plot_df.rename(columns={alt: x_col}, inplace=True)
                    if y_col not in plot_df.columns:
                        alt = f"{y_col}_y"
                        if alt in plot_df.columns:
                            plot_df.rename(columns={alt: y_col}, inplace=True)
                marker_candidates = _collect_marker_candidates(dataset_cfg)
                preferred_markers = marker_candidates[: max(args.max_markers, 1)]
                used_markers, missing = plot_marker_maps(
                    fig_dir / "fig2_marker_maps.png",
                    plot_df,
                    marker_candidates,
                    preferred_markers,
                    x_col,
                    y_col,
                    max_markers=args.max_markers,
                    clip_quantile=args.marker_clip_quantile,
                    log1p=args.marker_log1p,
                    dpi=args.dpi,
                )
                if missing:
                    missing_markers.extend([m for m in missing if m not in missing_markers])
                if used_markers:
                    figures["fig2c"] = "figures/fig2_marker_maps.png"

        compartments = _collect_compartments(dataset_cfg)
        if compartments and coords_df is not None and x_col and y_col:
            comp_values = {}
            comp_metrics: Dict[str, Dict[str, float]] = {}
            for comp in compartments:
                name = comp.get("name")
                types = comp.get("types") or []
                if not name or not isinstance(types, list):
                    continue
                cols = []
                for t in types:
                    canon = canonicalize_type_name(t, alias_map) if alias_map else normalize_type_name(t)
                    if canon in all_types:
                        cols.append(canon)
                if not cols:
                    missing_compartments.append(name)
                    continue
                comp_values[name] = {
                    "baseline": base_aligned[cols].sum(axis=1),
                    "route2": route_aligned[cols].sum(axis=1),
                }
                base_vals = comp_values[name]["baseline"].to_numpy(dtype=float)
                route_vals = comp_values[name]["route2"].to_numpy(dtype=float)
                delta_vals = route_vals - base_vals
                mask = np.isfinite(base_vals) & np.isfinite(delta_vals)
                corr = float(np.corrcoef(base_vals[mask], delta_vals[mask])[0, 1]) if mask.any() else float("nan")
                rmse = float(np.sqrt(np.mean((delta_vals[mask]) ** 2))) if mask.any() else float("nan")
                comp_metrics[name] = {"corr_delta": corr, "rmse_delta": rmse}
                delta_mass_val = float(np.nansum(route_vals - base_vals))
                comp_summary_rows.append({
                    "compartment": name,
                    "corr": corr,
                    "rmse": rmse,
                    "delta_mass": delta_mass_val,
                    "abs_delta": abs(delta_mass_val),
                })
            if comp_values:
                comp_delta = pd.DataFrame({"spot_id": delta_df["spot_id"]})
                spot_ids = comp_delta["spot_id"]
                for name, vals in comp_values.items():
                    base_vals = vals["baseline"].reindex(spot_ids)
                    route_vals = vals["route2"].reindex(spot_ids)
                    comp_delta[f"{name}__baseline"] = base_vals.values
                    comp_delta[f"{name}__route2"] = route_vals.values
                    comp_delta[f"{name}__delta"] = route_vals.values - base_vals.values
                plot_compartment_maps(
                    fig_dir / "fig2_compartment_maps.png",
                    comp_delta,
                    coords_df,
                    x_col,
                    y_col,
                    list(comp_values.keys()),
                    comp_metrics,
                    clip_quantile=args.delta_clip_quantile,
                    smooth_k=args.delta_smooth_k,
                    smooth_alpha=args.delta_smooth_alpha,
                    dpi=args.dpi,
                )
                if (fig_dir / "fig2_compartment_maps.png").exists():
                    figures["fig2d"] = "figures/fig2_compartment_maps.png"
                plot_compartment_summary_table(
                    fig_dir / "compartment_summary_table.png",
                    comp_summary_rows,
                    dpi=args.dpi,
                )
                if (fig_dir / "compartment_summary_table.png").exists():
                    figures["fig2e"] = "figures/compartment_summary_table.png"

        plot_pipeline_diagram(fig_dir / "fig3_pipeline.png", dpi=args.dpi)
        figures["fig3"] = "figures/fig3_pipeline.png"

    def _extract_filtered(audit: dict, summary: dict) -> int:
        val = (
            audit.get("ledger_integrity", {})
            .get("metrics", {})
            .get("filter_ledger", {})
            .get("observed_n_filtered")
        )
        if val is None:
            val = summary.get("n_filtered", 0)
        return int(val or 0)

    def _extract_rescued(audit: dict) -> int:
        val = (
            audit.get("ledger_integrity", {})
            .get("metrics", {})
            .get("rescue_ledger", {})
            .get("rescued_cells_count")
        )
        return int(val or 0)

    filter_summary = {
        "baseline": {
            "n_filtered": _extract_filtered(base_audit, base_summary),
            "n_rescued": _extract_rescued(base_audit),
        },
        "route2": {
            "n_filtered": _extract_filtered(route_audit, route_summary),
            "n_rescued": _extract_rescued(route_audit),
        },
    }
    pd.DataFrame([
        {"run_id": "baseline", **filter_summary["baseline"]},
        {"run_id": "route2", **filter_summary["route2"]},
    ]).to_csv(table_dir / "filter_rescue_summary.csv", index=False)

    delta_mass.sort_values(ascending=False).to_csv(table_dir / "delta_type_total_mass.csv", header=["delta_mass"])

    highlights = []
    if stage6_compare.get("global_distribution_js") is not None:
        highlights.append(f"global JS divergence: {stage6_compare['global_distribution_js']:.4f}")
    if filter_summary["route2"]["n_filtered"]:
        highlights.append(f"route2 filtered cells: {filter_summary['route2']['n_filtered']}")
    if top_types:
        highlights.append("top delta types: " + ", ".join(top_types[:3]))
    if missing_markers:
        highlights.append("missing markers: " + ", ".join(missing_markers))
    if missing_compartments:
        highlights.append("missing compartments: " + ", ".join(missing_compartments))

    build_report_html(
        report_html_dir / "Stage7_report.html",
        args.sample,
        args.baseline_id,
        args.route2_id,
        figures,
        highlights,
        abbrev_map,
        args.html_mode,
        out_dir,
        comp_summary_rows,
    )
    build_report_markdown(
        report_md_dir / "Stage7_summary_CN.md",
        args.sample,
        args.baseline_id,
        args.route2_id,
        figures,
        highlights,
        abbrev_map,
    )

    summary = {
        "sample": args.sample,
        "baseline_id": args.baseline_id,
        "route2_id": args.route2_id,
        "paths": {
            "stage4_baseline": str(base_run_dir),
            "stage4_route2": str(route_run_dir),
            "stage6_root": str(stage6_root),
            "stage1_export": str(stage1_export),
        },
        "top_delta_types": top_types,
        "filter_summary": filter_summary,
        "stage6_compare": stage6_compare,
    }
    (summary_dir / "overall_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    print(f"[Stage7] outputs -> {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
