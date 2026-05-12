#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

_HERE = Path(__file__).resolve()
_ROOT = _HERE.parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src.utils.sample_paths import infer_sim_group, resolve_sample_dir


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Three-panel spatial visualization for missing-type simulation: "
            "truth (with removed-type shown in gray), baseline mapping, route2 mapping."
        )
    )
    p.add_argument("--project_root", default=".", help="Project root path.")
    p.add_argument(
        "--sample",
        required=True,
        help="Simulation sample id, e.g. real_brca_clustered_sim_mt_b_cells_d100_fill_endo",
    )
    p.add_argument(
        "--missing_type",
        default=None,
        help="Missing type label override (single or comma-separated). Default from sim_info.json (missing_types/missing_type).",
    )
    p.add_argument(
        "--truth_show_removed_gray",
        default=None,
        choices=["true", "false"],
        help=(
            "Whether to paint source missing-type dominant spots as gray in truth panel. "
            "Default: auto (false when replacement_type exists in sim_info, else true)."
        ),
    )
    p.add_argument("--point_size", type=float, default=24.0, help="Spot point size.")
    p.add_argument("--alpha", type=float, default=0.93, help="Point alpha.")
    p.add_argument(
        "--out_png",
        default=None,
        help="Output png path (default: visualizations/simulations/<sim_group>/<sample>/missing_type_triptych.png).",
    )
    return p.parse_args()


def _load_sim_info(raw_dir: Path) -> dict:
    p = raw_dir / "sim_info.json"
    if not p.exists():
        raise FileNotFoundError(f"sim_info.json not found: {p}")
    return json.loads(p.read_text(encoding="utf-8-sig"))


def _parse_missing_types(raw: str | None) -> list[str]:
    if raw is None:
        return []
    txt = str(raw).strip()
    if not txt:
        return []
    txt = txt.replace(";", ",")
    out: list[str] = []
    for part in txt.split(","):
        t = str(part).strip()
        if t and t not in out:
            out.append(t)
    return out


def _normalize_missing_types(types: list[str]) -> list[str]:
    out: list[str] = []
    sentinel = {"__NO_MISSING__", "NO_MISSING", "NONE", "NULL", "NA", "N/A"}
    for t in types:
        x = str(t).strip()
        if not x:
            continue
        if x.upper() in sentinel:
            continue
        if x not in out:
            out.append(x)
    return out


def _read_coords(raw_dir: Path) -> pd.DataFrame:
    p = raw_dir / "brca_STdata_coordinates.txt"
    if not p.exists():
        raise FileNotFoundError(f"coordinate file not found: {p}")
    d = pd.read_csv(p, sep="\t")
    if d.shape[1] < 3:
        raise ValueError(f"invalid coordinate file: {p}")
    d = d.rename(columns={d.columns[0]: "spot_id", d.columns[1]: "row", d.columns[2]: "col"})
    d["spot_id"] = d["spot_id"].astype(str)
    return d[["spot_id", "row", "col"]]


def _dominant_from_fraction(csv_path: Path) -> tuple[pd.DataFrame, list[str]]:
    d = pd.read_csv(csv_path)
    if d.shape[1] < 2:
        raise ValueError(f"invalid fraction file: {csv_path}")
    d = d.rename(columns={d.columns[0]: "spot_id"})
    type_cols = [c for c in d.columns if c != "spot_id"]
    d[type_cols] = d[type_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    row_sum = d[type_cols].sum(axis=1)
    d["dominant"] = d[type_cols].idxmax(axis=1)
    d.loc[row_sum <= 0.0, "dominant"] = "__NoType__"
    return d[["spot_id", "dominant"]], type_cols


def _dominant_from_assignment(csv_path: Path) -> tuple[pd.DataFrame, list[str]]:
    d = pd.read_csv(csv_path)
    if d.shape[1] < 3:
        raise ValueError(f"invalid assignment file: {csv_path}")
    d = d.rename(columns={d.columns[0]: "spot_id"})
    ignore_cols = {"spot_id", "Total cells", "Unknown_sc_only"}
    type_cols = [c for c in d.columns if c not in ignore_cols]
    d[type_cols] = d[type_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    row_sum = d[type_cols].sum(axis=1)
    d["dominant"] = d[type_cols].idxmax(axis=1)
    d.loc[row_sum <= 0.0, "dominant"] = "__NoType__"
    return d[["spot_id", "dominant"]], type_cols


def _palette(type_list: list[str]) -> dict[str, str]:
    # Keep exactly the same base palette/order logic as visualize_sim_spatial_types.py
    # so colors are directly comparable across figures.
    base = [
        "#e41a1c",  # red
        "#377eb8",  # blue
        "#4daf4a",  # green
        "#984ea3",  # purple
        "#ff7f00",  # orange
        "#b8860b",  # dark-gold
        "#a65628",  # brown
        "#f781bf",  # pink
        "#17becf",  # cyan
        "#1b9e77",  # teal-green
        "#6a3d9a",  # deep violet
        "#b15928",
        "#00b4d8",
        "#fb5607",
        "#3a86ff",
    ]
    pal = {t: base[i % len(base)] for i, t in enumerate(type_list)}
    if "B cells" in pal and "Epithelial cells" in pal:
        pal["B cells"], pal["Epithelial cells"] = pal["Epithelial cells"], pal["B cells"]
    pal["__RemovedMissing__"] = "#bdbdbd"
    pal["__NoType__"] = "#e5e5e5"
    return pal


def _plot_panel(ax, df: pd.DataFrame, title: str, palette: dict[str, str], point_size: float, alpha: float) -> None:
    ax.set_facecolor("#e6e6e6")
    order = ["__NoType__"] + [x for x in palette.keys() if x in set(df["dominant"].tolist()) and x != "__NoType__"]
    for t in order:
        sub = df[df["dominant"] == t]
        if sub.empty:
            continue
        c = palette.get(t, "#333333")
        z = 1 if t == "__NoType__" else 2
        a = 0.55 if t == "__NoType__" else alpha
        ax.scatter(
            sub["col"],
            sub["row"],
            s=point_size,
            c=c,
            alpha=a,
            edgecolors="none",
            linewidths=0.0,
            zorder=z,
            rasterized=True,
        )
    ax.set_title(title, fontsize=12, weight="bold")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.invert_yaxis()
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)


def main() -> int:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    sample = args.sample

    raw_dir = resolve_sample_dir(project_root, sample, sim_group="real_brca", must_exist=True)
    res_dir = project_root / "result" / sample
    sim_group = infer_sim_group(project_root, sample, sim_group="real_brca") or "ungrouped"
    out_png = (
        Path(args.out_png).resolve()
        if args.out_png
        else (project_root / "visualizations" / "simulations" / sim_group / sample / "missing_type_triptych.png")
    )
    out_png.parent.mkdir(parents=True, exist_ok=True)

    sim_info = _load_sim_info(raw_dir)
    source_sample = str(sim_info.get("source_sample") or "").strip()
    cli_missing_types = _parse_missing_types(args.missing_type)
    sim_missing_types = []
    sim_missing_types.extend(_parse_missing_types(",".join([str(x) for x in (sim_info.get("missing_types") or [])])))
    sim_missing_types.extend(_parse_missing_types(str(sim_info.get("missing_type") or "")))
    missing_types = _normalize_missing_types(cli_missing_types or sim_missing_types)
    no_missing_mode = len(missing_types) == 0
    missing_type = missing_types[0] if missing_types else None
    missing_type_label = ", ".join(missing_types) if missing_types else "none"
    missing_type_set = set(missing_types)
    replacement_type = str(sim_info.get("replacement_type") or "").strip() or None
    if args.truth_show_removed_gray is None:
        truth_show_removed_gray = (not no_missing_mode) and (replacement_type is None)
    else:
        truth_show_removed_gray = args.truth_show_removed_gray.lower() == "true"

    coords = _read_coords(raw_dir)

    target_truth_csv = raw_dir / "sim_truth_spot_type_fraction.csv"
    source_truth_csv = None
    if truth_show_removed_gray:
        if not source_sample:
            raise ValueError("source_sample is missing in sim_info.json but truth_show_removed_gray=True.")
        source_truth_csv = (
            resolve_sample_dir(project_root, source_sample, sim_group="real_brca", must_exist=True)
            / "sim_truth_spot_type_fraction.csv"
        )
    base_assign_csv = res_dir / "stage4_cytospace_baseline" / "cytospace_output" / "cell_type_assignments_by_spot.csv"
    route2_assign_csv = res_dir / "stage4_cytospace_route2" / "cytospace_output" / "cell_type_assignments_by_spot.csv"

    if not target_truth_csv.exists():
        raise FileNotFoundError(f"missing target truth file: {target_truth_csv}")
    if source_truth_csv is not None and not source_truth_csv.exists():
        raise FileNotFoundError(f"missing source truth file: {source_truth_csv}")
    if not base_assign_csv.exists():
        raise FileNotFoundError(f"missing baseline assignment file: {base_assign_csv}")
    if not route2_assign_csv.exists():
        raise FileNotFoundError(f"missing route2 assignment file: {route2_assign_csv}")

    truth_target, truth_types_target = _dominant_from_fraction(target_truth_csv)
    truth_source, truth_types_source = (
        _dominant_from_fraction(source_truth_csv)
        if source_truth_csv is not None
        else (pd.DataFrame(columns=["spot_id", "dominant"]), [])
    )
    base_dom, base_types = _dominant_from_assignment(base_assign_csv)
    route2_dom, route2_types = _dominant_from_assignment(route2_assign_csv)

    # Reference order from source truth columns -> same color mapping as preview_truth figure.
    type_order = list(truth_types_source) if truth_types_source else list(truth_types_target)
    for t in list(truth_types_target) + list(base_types) + list(route2_types):
        if t not in type_order:
            type_order.append(t)
    pal = _palette(type_order)

    # Panel-1: use target truth, but paint source-missing dominant spots as gray.
    truth_panel = coords.merge(truth_target, on="spot_id", how="inner")
    removed_mask = pd.Series(False, index=truth_panel.index)
    if truth_show_removed_gray:
        truth_panel = truth_panel.merge(
            truth_source.rename(columns={"dominant": "source_dominant"}),
            on="spot_id",
            how="left",
        )
        removed_mask = truth_panel["source_dominant"].isin(missing_type_set)
        truth_panel.loc[removed_mask, "dominant"] = "__RemovedMissing__"

    baseline_panel = coords.merge(base_dom, on="spot_id", how="left")
    route2_panel = coords.merge(route2_dom, on="spot_id", how="left")
    baseline_panel["dominant"] = baseline_panel["dominant"].fillna("__NoType__")
    route2_panel["dominant"] = route2_panel["dominant"].fillna("__NoType__")

    fig, axes = plt.subplots(1, 3, figsize=(24, 8), dpi=180)
    fig.patch.set_facecolor("#e6e6e6")
    _plot_panel(
        axes[0],
        truth_panel,
        title=(
            "Truth (no missing type)"
            if no_missing_mode
            else (
                f"Truth (missing {missing_type_label}; removed spots in gray)"
                if truth_show_removed_gray
                else f"Truth (missing {missing_type_label}; filled with {replacement_type or 'replacement'})"
            )
        ),
        palette=pal,
        point_size=args.point_size,
        alpha=args.alpha,
    )
    _plot_panel(
        axes[1],
        baseline_panel,
        title="CytoSPACE Baseline Mapping",
        palette=pal,
        point_size=args.point_size,
        alpha=args.alpha,
    )
    _plot_panel(
        axes[2],
        route2_panel,
        title="SVTuner + CytoSPACE Mapping (Route2)",
        palette=pal,
        point_size=args.point_size,
        alpha=args.alpha,
    )

    # Compact legend on the right of the last panel.
    legend_types = [t for t in type_order if t != "__NoType__"]
    if truth_show_removed_gray and "__RemovedMissing__" not in legend_types:
        legend_types.append("__RemovedMissing__")
    if "__NoType__" not in legend_types:
        legend_types.append("__NoType__")
    legend_name = {t: t for t in legend_types}
    legend_name["__RemovedMissing__"] = f"Removed {missing_type_label}"
    legend_name["__NoType__"] = "Unassigned spots"
    handles = [
        plt.Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            markersize=7,
            markerfacecolor=pal[t],
            markeredgecolor="none",
            label=legend_name[t],
        )
        for t in legend_types
    ]
    axes[2].legend(handles=handles, title="Type", loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=8)

    fig.suptitle(
        f"{sample}: Truth vs Baseline vs Route2",
        fontsize=14,
        weight="bold",
        y=0.965,
    )
    fig.subplots_adjust(left=0.04, right=0.86, bottom=0.08, top=0.90, wspace=0.14)
    fig.savefig(out_png)
    plt.close(fig)

    print(f"[OK] wrote: {out_png}")
    print(f"[INFO] sample={sample} source_sample={source_sample} missing_types={missing_types}")
    print(f"[INFO] replacement_type={replacement_type}")
    print(f"[INFO] truth_show_removed_gray={truth_show_removed_gray}")
    print(f"[INFO] removed_spots(gray)={int(removed_mask.sum())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
