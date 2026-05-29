#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colormaps
from matplotlib.colors import LinearSegmentedColormap, Normalize
from PIL import Image, ImageDraw, ImageFont


HIGHLIGHT_COLOR = "#00CFE8"
BASE_BG = "#E4E0D8"
N_PANEL_GENES = 30


def _make_deep_purple_magma() -> LinearSegmentedColormap:
    base = colormaps["magma"]
    colors = base(np.linspace(0.18, 1.0, 256))
    return LinearSegmentedColormap.from_list("magma_deep_purple", colors)


SIGNATURE_CMAP = _make_deep_purple_magma()


def _slugify(value: str) -> str:
    return re.sub(r"_+", "_", re.sub(r"[^a-z0-9]+", "_", str(value).lower())).strip("_") or "target"


def _load_font(size: int, bold: bool = False) -> ImageFont.ImageFont:
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


def _processed_dir(root: Path, sample: str) -> Path:
    candidates = [
        root / "data" / "processed" / sample,
        root / "data" / "processed" / "cytospace_fig2c_melanoma" / sample,
        root / "data" / "processed" / "cytospace_fig2d_tme" / sample,
        root / "data" / "processed" / "cytospace_fig2d_profile_mask" / sample,
    ]
    for path in candidates:
        if (path / "stage1_preprocess" / "exported").exists():
            return path
    raise FileNotFoundError(f"cannot resolve processed sample: {sample}")


def _read_gene_panel(masked_processed: Path, source_expr: Path, masked_expr: Path) -> list[str]:
    panel_path = masked_processed / "stage1_preprocess" / "fig2d_profile_mask_gene_panel.csv"
    if not panel_path.exists():
        raise FileNotFoundError(panel_path)
    source_cols = set(pd.read_csv(source_expr, nrows=0).columns)
    mask_cols = set(pd.read_csv(masked_expr, nrows=0).columns)
    allowed = source_cols.intersection(mask_cols)
    panel = pd.read_csv(panel_path)
    genes: list[str] = []
    for gene in panel["gene"].astype(str):
        if gene in allowed and gene not in genes:
            genes.append(gene)
        if len(genes) >= N_PANEL_GENES:
            break
    if not genes:
        raise ValueError(f"no usable marker genes in {panel_path}")
    return genes


def _load_st_signature(processed: Path, genes: list[str]) -> pd.DataFrame:
    export = processed / "stage1_preprocess" / "exported"
    expr_path = export / "st_expression_normalized.csv"
    coord_path = export / "st_coordinates.csv"
    expr = pd.read_csv(expr_path, usecols=["spot_id", *genes])
    coords = pd.read_csv(coord_path, usecols=["spot_id", "row", "col"])
    expr["spot_id"] = expr["spot_id"].astype(str)
    coords["spot_id"] = coords["spot_id"].astype(str)
    df = coords.merge(expr, on="spot_id", how="inner")
    df["target_signature"] = df[genes].mean(axis=1)
    return df[["spot_id", "row", "col", "target_signature"]]


def _resolve_assignment_column(columns: list[str], target_type: str) -> str | None:
    if target_type in columns:
        return target_type
    norm = {str(c).strip().casefold(): c for c in columns}
    return norm.get(str(target_type).strip().casefold())


def _load_overlay(root: Path, sample: str, target_type: str, stage4_dir: str) -> pd.DataFrame:
    path = root / "result" / sample / stage4_dir / "cytospace_output" / "cell_type_assignments_by_spot.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path)
    df = df.rename(columns={df.columns[0]: "spot_id"})
    df["spot_id"] = df["spot_id"].astype(str)
    col = _resolve_assignment_column(df.columns.tolist(), target_type)
    if col is None:
        df["target_cells"] = 0.0
    else:
        df["target_cells"] = pd.to_numeric(df[col], errors="coerce").fillna(0.0)
    df["target_present"] = df["target_cells"] > 0
    return df[["spot_id", "target_cells", "target_present"]]


def _is_sparse_integer_grid(df: pd.DataFrame) -> bool:
    span = max(
        float(df["col"].max() - df["col"].min()),
        float(df["row"].max() - df["row"].min()),
    )
    return len(df) < 600 and span < 100


def _spot_style(df: pd.DataFrame) -> dict[str, float]:
    """Use tile-like rendered spots for sparse integer-grid ST slides."""
    if _is_sparse_integer_grid(df):
        return {
            "bg_s": 135.0,
            "signal_s": 135.0,
            "outline_s": 112.0,
            "outline_lw": 1.05,
            "alpha": 0.98,
            "marker": "h",
        }
    return {
        "bg_s": 5.0,
        "signal_s": 9.0,
        "outline_s": 26.0,
        "outline_lw": 0.70,
        "alpha": 0.95,
        "marker": "o",
    }


def _expand_sparse_grid_for_display(df: pd.DataFrame) -> pd.DataFrame:
    # Render each low-resolution ST spot as a compact local patch. This is
    # visual-only and keeps all spot-level values unchanged.
    offsets = np.array([-0.32, -0.16, 0.0, 0.16, 0.32])
    dx, dy = np.meshgrid(offsets, offsets)
    # Round mask makes the patch look like a filled Visium spot instead of a square.
    keep = (dx.ravel() ** 2 + dy.ravel() ** 2) <= 0.35**2
    dx = dx.ravel()[keep]
    dy = dy.ravel()[keep]
    pieces = []
    base = df.reset_index(drop=True)
    for ox, oy in zip(dx, dy):
        tmp = base.copy()
        tmp["col"] = tmp["col"].astype(float) + float(ox)
        tmp["row"] = tmp["row"].astype(float) + float(oy)
        pieces.append(tmp)
    return pd.concat(pieces, ignore_index=True)


def _draw_panel(ax, df: pd.DataFrame, title: str, subtitle: str, norm: Normalize, outline_col: str | None = None):
    style = _spot_style(df)
    display_df = df
    ax.scatter(
        display_df["col"],
        -display_df["row"],
        s=style["bg_s"],
        marker=style["marker"],
        c=BASE_BG,
        linewidths=0,
        alpha=0.72,
        rasterized=True,
    )
    sc = ax.scatter(
        display_df["col"],
        -display_df["row"],
        s=style["signal_s"],
        marker=style["marker"],
        c=display_df["target_signature"],
        cmap=SIGNATURE_CMAP,
        norm=norm,
        linewidths=0,
        alpha=style["alpha"],
        rasterized=True,
    )
    if outline_col is not None and outline_col in df.columns:
        hi = df[df[outline_col].astype(bool)]
        if not hi.empty:
            ax.scatter(
                hi["col"],
                -hi["row"],
                s=style["outline_s"],
                marker="o",
                facecolors="none",
                edgecolors=HIGHLIGHT_COLOR,
                linewidths=style["outline_lw"],
                alpha=1.0,
                rasterized=False,
                zorder=10,
            )
    ax.set_title(f"{title}\n{subtitle}", loc="left", fontweight="bold", fontsize=9.6, pad=5)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("spatial col", fontsize=7.6)
    return sc


def _render_one(root: Path, row: pd.Series, out_root: Path) -> Path:
    pair_id = str(row["pair_id"])
    source_sample = str(row["source_sample"])
    sample = str(row["profile_mask_sample"])
    target_type = str(row["masked_target_type"])
    source_processed = _processed_dir(root, source_sample)
    masked_processed = _processed_dir(root, sample)
    source_expr = source_processed / "stage1_preprocess" / "exported" / "st_expression_normalized.csv"
    masked_expr = masked_processed / "stage1_preprocess" / "exported" / "st_expression_normalized.csv"
    genes = _read_gene_panel(masked_processed, source_expr, masked_expr)

    source = _load_st_signature(source_processed, genes)
    masked = _load_st_signature(masked_processed, genes)
    source["expr_high"] = source["target_signature"] >= source["target_signature"].quantile(0.90)
    masked["expr_high"] = masked["target_signature"] >= source["target_signature"].quantile(0.90)

    baseline = masked.merge(
        _load_overlay(root, sample, target_type, "stage4_cytospace_baseline_profile300"),
        on="spot_id",
        how="left",
    )
    route2 = masked.merge(
        _load_overlay(root, sample, target_type, "stage4_cytospace_route2_profile300"),
        on="spot_id",
        how="left",
    )
    for frame in [baseline, route2]:
        frame["target_present"] = frame["target_present"].fillna(False).astype(bool)
        frame["target_cells"] = frame["target_cells"].fillna(0.0)

    combined = pd.concat([source["target_signature"], masked["target_signature"]], ignore_index=True)
    norm = Normalize(float(np.nanpercentile(combined, 1)), float(np.nanpercentile(combined, 99)))

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, axes = plt.subplots(1, 4, figsize=(15.6, 4.8), dpi=220)
    sc = _draw_panel(axes[0], source, "A  Original ST", "target marker signature", norm, "expr_high")
    _draw_panel(axes[1], masked, "B  Profile-mask ST", "target signature after masking", norm, "expr_high")
    _draw_panel(axes[2], baseline, "C  CytoSPACE baseline", "mapped target outline", norm, "target_present")
    _draw_panel(axes[3], route2, "D  SVTuner + CytoSPACE", "mapped target outline", norm, "target_present")
    axes[0].set_ylabel("spatial row", fontsize=11)

    fig.suptitle(
        f"High-resolution profile-mask mapping: {pair_id} ({target_type})",
        x=0.03,
        ha="left",
        fontsize=13.5,
        fontweight="bold",
    )
    fig.subplots_adjust(left=0.024, right=0.928, top=0.79, bottom=0.15, wspace=-0.06)
    cax = fig.add_axes([0.941, 0.22, 0.013, 0.54])
    cb = fig.colorbar(sc, cax=cax)
    cb.set_label(f"{target_type} marker signature", fontsize=7.2, labelpad=5)
    cb.ax.tick_params(labelsize=7.2)

    out_dir = out_root / sample
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / f"{_slugify(target_type)}_mask_ABCD_mapping.png"
    fig.savefig(out_png, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[OK] wrote: {out_png}")
    return out_png


def _stack(pngs: list[Path], manifest: pd.DataFrame, out_root: Path) -> Path:
    images = [Image.open(path).convert("RGB") for path in pngs]
    target_w = 2400
    resized = []
    for image in images:
        scale = target_w / image.size[0]
        resized.append(image.resize((target_w, int(round(image.size[1] * scale))), Image.Resampling.LANCZOS))

    outer = 30
    title_h = 72
    row_gap = 18
    canvas_w = target_w + outer * 2
    canvas_h = outer * 2 + title_h + sum(img.size[1] for img in resized) + row_gap * (len(resized) - 1)
    canvas = Image.new("RGB", (canvas_w, canvas_h), (255, 255, 255))
    draw = ImageDraw.Draw(canvas)
    title_font = _load_font(36, bold=True)

    title = "High-resolution Profile-mask Mapping Overview"
    bb = draw.textbbox((0, 0), title, font=title_font)
    draw.text(((canvas_w - (bb[2] - bb[0])) // 2, outer + 2), title, fill=(24, 24, 24), font=title_font)

    y = outer + title_h
    for image in resized:
        canvas.paste(image, (outer, y))
        y += image.size[1] + row_gap

    out_png = out_root / "cytospace_fig2d_profile_mask_mapping_stack_6x4.png"
    canvas.save(out_png, optimize=True)
    print(f"[OK] wrote: {out_png}")
    return out_png


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project_root", default=".", type=Path)
    parser.add_argument(
        "--manifest",
        default="visualizations/cytospace_fig2d_profile_mask_benchmark/fig2d_profile_mask_benchmark_profile_mask_manifest.csv",
        type=Path,
    )
    args = parser.parse_args()
    root = args.project_root.resolve()
    manifest_path = args.manifest if args.manifest.is_absolute() else root / args.manifest
    manifest = pd.read_csv(manifest_path)
    out_root = root / "visualizations" / "cytospace_fig2d_profile_mask_mapping"
    out_root.mkdir(parents=True, exist_ok=True)
    pngs = [_render_one(root, row, out_root) for _, row in manifest.iterrows()]
    _stack(pngs, manifest, out_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
