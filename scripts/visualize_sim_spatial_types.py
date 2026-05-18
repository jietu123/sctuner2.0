#!/usr/bin/env python
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgba

_HERE = Path(__file__).resolve()
_ROOT = _HERE.parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src.utils.sample_paths import infer_sim_group, resolve_sample_dir


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Visualize simulated cell-type spatial distribution from truth spot-type fractions."
    )
    p.add_argument("--project_root", default=".", help="Project root path.")
    p.add_argument("--sample", default="real_brca_clustered_sim", help="Sample id under data/sim or data/raw.")
    p.add_argument(
        "--truth_spot_csv",
        default=None,
        help="Path to sim_truth_spot_type_fraction.csv (default: <sample_dir>/sim_truth_spot_type_fraction.csv).",
    )
    p.add_argument(
        "--coord_csv",
        default=None,
        help="Path to ST coordinate file (default: <sample_dir>/brca_STdata_coordinates.txt).",
    )
    p.add_argument(
        "--out_dir",
        default=None,
        help="Output dir (default: visualizations/simulations/<sim_group>/<sample>).",
    )
    p.add_argument("--point_size", type=float, default=10.0, help="Spot point size.")
    p.add_argument("--bg_point_size", type=float, default=8.0, help="Background point size.")
    p.add_argument("--presence_threshold", type=float, default=0.02, help="Presence threshold for per-type panels.")
    p.add_argument(
        "--style",
        choices=["default", "triptych_like"],
        default="default",
        help="Visual style preset. 'triptych_like' aligns style with missing_type_triptych.",
    )
    return p.parse_args()


def _prepare_inputs(args: argparse.Namespace) -> tuple[Path, Path, Path]:
    project_root = Path(args.project_root).resolve()
    sample_dir = resolve_sample_dir(project_root, args.sample, sim_group="real_brca", must_exist=True)
    sim_group = infer_sim_group(project_root, args.sample, sim_group="real_brca") or "ungrouped"
    truth_spot_csv = Path(args.truth_spot_csv).resolve() if args.truth_spot_csv else sample_dir / "sim_truth_spot_type_fraction.csv"
    coord_csv = Path(args.coord_csv).resolve() if args.coord_csv else sample_dir / "brca_STdata_coordinates.txt"
    out_dir = (
        Path(args.out_dir).resolve()
        if args.out_dir
        else project_root / "visualizations" / "simulations" / sim_group / args.sample
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    return truth_spot_csv, coord_csv, out_dir


def _load_data(truth_spot_csv: Path, coord_csv: Path) -> tuple[pd.DataFrame, list[str]]:
    truth = pd.read_csv(truth_spot_csv)
    if truth.shape[1] < 2:
        raise ValueError(f"Invalid truth spot file: {truth_spot_csv}")
    truth = truth.rename(columns={truth.columns[0]: "spot_id"})
    type_cols = [c for c in truth.columns if c != "spot_id"]

    coords = pd.read_csv(coord_csv, sep="\t")
    if coords.shape[1] < 3:
        raise ValueError(f"Invalid coordinate file: {coord_csv}")
    coords = coords.rename(columns={coords.columns[0]: "spot_id", coords.columns[1]: "row", coords.columns[2]: "col"})

    merged = coords.merge(truth, on="spot_id", how="inner")
    if merged.empty:
        raise ValueError("No overlapping spot_id between coordinates and truth fractions.")

    merged["dominant_type"] = merged[type_cols].idxmax(axis=1)
    return merged, type_cols


def _make_palette(type_cols: list[str]) -> dict[str, str]:
    base = [
        "#e41a1c",  # red
        "#377eb8",  # blue
        "#4daf4a",  # green
        "#984ea3",  # purple
        "#ff7f00",  # orange
        "#b8860b",  # dark-gold (avoid over-bright yellow washout)
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
    palette = {t: base[i % len(base)] for i, t in enumerate(type_cols)}
    # Keep real_brca color semantics aligned with missing_type_triptych:
    # red should represent Epithelial cells and light purple B cells.
    if "B cells" in palette and "Epithelial cells" in palette:
        palette["B cells"], palette["Epithelial cells"] = palette["Epithelial cells"], palette["B cells"]
    return palette


def plot_dominant_map(
    df: pd.DataFrame,
    type_cols: list[str],
    palette: dict[str, str],
    out_png: Path,
    point_size: float,
    bg_point_size: float,
    sample: str,
    style: str = "default",
) -> None:
    if style == "triptych_like":
        fig, ax = plt.subplots(figsize=(14.8, 10.4), dpi=180)
        fig.patch.set_facecolor("#E6E6E6")
        ax.set_facecolor("#E6E6E6")
        bg_color = "#D7D7D7"
        bg_alpha = 0.22
        fg_alpha = 0.94
        edgecolors = "none"
        edge_width = 0.0
    else:
        fig, ax = plt.subplots(figsize=(14, 10), dpi=180)
        ax.set_facecolor("#F2F2F2")
        bg_color = "#D9D9D9"
        bg_alpha = 0.28
        fg_alpha = 0.92
        edgecolors = "#202020"
        edge_width = 0.12

    ax.scatter(df["col"], df["row"], s=bg_point_size, c=bg_color, alpha=bg_alpha, linewidths=0)

    counts = df["dominant_type"].value_counts().reindex(type_cols, fill_value=0)
    for t in type_cols:
        sub = df[df["dominant_type"] == t]
        if sub.empty:
            continue
        ax.scatter(
            sub["col"],
            sub["row"],
            s=point_size,
            c=palette[t],
            alpha=fg_alpha,
            edgecolors=edgecolors,
            linewidths=edge_width,
            label=f"{t} ({int(counts[t])})",
        )

    ax.set_title(f"{sample}: Dominant Cell-Type Spatial Map", fontsize=18, weight="bold")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.invert_yaxis()
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)
    ax.legend(title="Dominant Type", loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)


def plot_type_panels(
    df: pd.DataFrame,
    type_cols: list[str],
    palette: dict[str, str],
    out_png: Path,
    point_size: float,
    bg_point_size: float,
    presence_threshold: float,
) -> None:
    n = len(type_cols)
    ncol = 4
    nrow = math.ceil(n / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.8 * ncol, 4.6 * nrow), dpi=180)
    axes = np.array(axes).reshape(-1)

    x = df["col"].to_numpy()
    y = df["row"].to_numpy()

    for i, t in enumerate(type_cols):
        ax = axes[i]
        ax.set_facecolor("#F2F2F2")
        ax.scatter(x, y, s=bg_point_size, c="#DCDCDC", alpha=0.25, linewidths=0)

        vals = df[t].to_numpy(dtype=np.float32)
        mask = vals >= presence_threshold
        if mask.any():
            base_rgb = np.array(to_rgba(palette[t]))
            colors = np.tile(base_rgb, (mask.sum(), 1))
            alpha = np.clip(vals[mask], 0.25, 1.0)
            colors[:, 3] = alpha
            ax.scatter(
                x[mask],
                y[mask],
                s=point_size,
                c=colors,
                edgecolors="#202020",
                linewidths=0.10,
            )
            shown = int(mask.sum())
        else:
            shown = 0

        ax.set_title(f"{t}\nspots≥{presence_threshold:.2f}: {shown}", fontsize=11, weight="bold")
        ax.invert_yaxis()
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([])
        ax.set_yticks([])

    for j in range(n, len(axes)):
        axes[j].axis("off")

    fig.suptitle("Per-Type Spatial Distribution", fontsize=18, weight="bold", y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    truth_spot_csv, coord_csv, out_dir = _prepare_inputs(args)
    print(f"[INPUT] truth={truth_spot_csv}")
    print(f"[INPUT] coord={coord_csv}")
    print(f"[OUTDIR] {out_dir}")

    df, type_cols = _load_data(truth_spot_csv, coord_csv)
    palette = _make_palette(type_cols)

    out1 = out_dir / "cell_type_dominant_map_preview.png"
    out2 = out_dir / "cell_type_panels_preview.png"

    plot_dominant_map(
        df=df,
        type_cols=type_cols,
        palette=palette,
        out_png=out1,
        point_size=args.point_size,
        bg_point_size=args.bg_point_size,
        sample=args.sample,
        style=args.style,
    )
    plot_type_panels(
        df=df,
        type_cols=type_cols,
        palette=palette,
        out_png=out2,
        point_size=args.point_size,
        bg_point_size=args.bg_point_size,
        presence_threshold=args.presence_threshold,
    )

    print(f"[OK] wrote: {out1}")
    print(f"[OK] wrote: {out2}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
