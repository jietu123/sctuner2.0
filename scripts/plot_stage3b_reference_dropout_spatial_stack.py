#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colormaps
from matplotlib.colors import LinearSegmentedColormap, Normalize
from PIL import Image, ImageDraw, ImageFont

HIGHLIGHT_COLOR = "#00CFE8"
HIGHLIGHT_LW = 0.55
BASE_BG = "#E4E0D8"
PANEL_BORDER_COLOR = "#9A9A9A"
PANEL_BORDER_LW = 0.72


def make_deep_purple_magma() -> LinearSegmentedColormap:
    base = colormaps["magma"]
    colors = base(np.linspace(0.18, 1.0, 256))
    return LinearSegmentedColormap.from_list("magma_deep_purple", colors)


SIGNATURE_CMAP = make_deep_purple_magma()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Render final stack-style Stage3B reference-dropout spatial validation panels.")
    p.add_argument("--project_root", default=".")
    p.add_argument(
        "--metadata_csv",
        default="visualizations/stage3b_realdata_candidate_scan/spatial_9x2/stage3b_reference_dropout_spatial_stack_recommended_6x2_metadata.csv",
    )
    p.add_argument("--out_dir", default="visualizations/stage3b_realdata_candidate_scan/spatial_9x2")
    p.add_argument("--out_prefix", default="stage3b_reference_dropout_spatial_stack_recommended_6x2")
    return p.parse_args()


def load_font(size: int, bold: bool = False) -> ImageFont.ImageFont:
    candidates = []
    if bold:
        candidates.extend([r"C:\Windows\Fonts\arialbd.ttf", r"C:\Windows\Fonts\segoeuib.ttf", r"C:\Windows\Fonts\calibrib.ttf"])
    candidates.extend([r"C:\Windows\Fonts\arial.ttf", r"C:\Windows\Fonts\segoeui.ttf", r"C:\Windows\Fonts\calibri.ttf"])
    for path in candidates:
        try:
            return ImageFont.truetype(path, size=size)
        except Exception:
            pass
    return ImageFont.load_default()


def read_coords(path: Path) -> pd.DataFrame:
    coords = pd.read_csv(path)
    id_col = "spot_id" if "spot_id" in coords.columns else coords.columns[0]
    coords[id_col] = coords[id_col].astype(str)
    coords = coords.drop_duplicates(id_col).set_index(id_col)
    lower = {c.lower(): c for c in coords.columns}
    if "col" in lower and "row" in lower:
        x_col, y_col = lower["col"], lower["row"]
    elif "x" in lower and "y" in lower:
        x_col, y_col = lower["x"], lower["y"]
    else:
        numeric = [c for c in coords.columns if pd.api.types.is_numeric_dtype(coords[c])]
        if len(numeric) < 2:
            raise ValueError(f"Cannot infer coordinate columns: {path}")
        y_col, x_col = numeric[0], numeric[1]
    out = coords[[x_col, y_col]].apply(pd.to_numeric, errors="coerce").dropna()
    out.columns = ["x", "y"]
    out["y_plot"] = -out["y"]
    return out


def read_marker_percentile(expr_path: Path, genes: list[str]) -> pd.Series:
    header = pd.read_csv(expr_path, nrows=0).columns.tolist()
    id_col = header[0]
    present = [g for g in genes if g in header]
    if not present:
        raise ValueError(f"No marker genes present in {expr_path}")
    expr = pd.read_csv(expr_path, index_col=0, usecols=[id_col, *present])
    expr.index = expr.index.astype(str)
    expr = expr.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return expr.mean(axis=1).rank(method="average", pct=True)


def display_label(pair_id: str, target_type: str) -> str:
    mapping = {
        "brca_tnbc_fresh_frozen": "BRCA TNBC",
        "brca_her2_ffpe": "BRCA HER2 FFPE",
        "crc_fresh_frozen": "CRC fresh frozen",
        "mouse_embryo_real": "Mouse embryo",
        "human_intestine_cancer_real": "Human intestine cancer",
    }
    return f"{mapping.get(pair_id, pair_id.replace('_', ' '))} ({target_type})"


def point_size(n_spots: int, pair_id: str = "") -> float:
    if n_spots <= 1500:
        return 5.63
    if n_spots <= 3500:
        return 5.80 if pair_id == "brca_her2_ffpe" else 4.45
    return 2.37


def draw_panel(
    ax,
    df: pd.DataFrame,
    mask: pd.Series,
    title: str,
    subtitle: str,
    norm: Normalize,
    compact_svg: bool = False,
):
    size = point_size(len(df), str(getattr(draw_panel, "pair_id", "")))
    ax.scatter(df["x"], df["y_plot"], s=size, c=BASE_BG, marker="h", linewidths=0, alpha=0.72, rasterized=True)
    sc = ax.scatter(
        df["x"],
        df["y_plot"],
        s=size * 1.00,
        c=df["marker_percentile"],
        cmap=SIGNATURE_CMAP,
        norm=norm,
        linewidths=0,
        alpha=0.95,
        rasterized=True,
    )
    hi = df[mask.to_numpy()]
    if not hi.empty:
        ax.scatter(
            hi["x"],
            hi["y_plot"],
            s=max(size * 2.48, size + 4.0),
            facecolors="none",
            marker="h",
            edgecolors=HIGHLIGHT_COLOR,
            linewidths=HIGHLIGHT_LW,
            alpha=1.0,
            rasterized=compact_svg,
            zorder=10,
        )
    ax.set_title(f"{title}\n{subtitle}", loc="left", fontweight="bold", fontsize=9.3, pad=5)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("spatial col", fontsize=7.5)
    ax.set_ylabel("spatial row" if title.startswith("A") else "", fontsize=9.0)
    for spine in ax.spines.values():
        spine.set_visible(compact_svg)
        if compact_svg:
            spine.set_color(PANEL_BORDER_COLOR)
            spine.set_linewidth(PANEL_BORDER_LW)
    return sc


def load_scene(root: Path, row: pd.Series) -> pd.DataFrame:
    source_export = root / "data" / "processed" / str(row["storage_group"]) / str(row["source_sample"]) / "stage1_preprocess" / "exported"
    genes = [g for g in str(row["marker_genes_list"]).split(";") if g]
    coords = read_coords(source_export / "st_coordinates.csv")
    marker = read_marker_percentile(source_export / "st_expression_normalized.csv", genes)
    scores_path = root / "data" / "processed" / str(row["stage3b_group"]) / str(row["sample"]) / "stage3b_st_unsupported" / "spot_unsupported_scores.csv"
    scores = pd.read_csv(scores_path, index_col=0)
    scores.index = scores.index.astype(str)
    blank = scores["is_unsupported_region"]
    if blank.dtype != bool:
        blank = blank.astype(str).str.lower().isin(["true", "1", "yes"])
    common = coords.index.intersection(marker.index).intersection(blank.index)
    df = coords.loc[common].copy()
    df["marker_percentile"] = marker.loc[common]
    df["blank"] = blank.loc[common].to_numpy()
    df["target_top15"] = df["marker_percentile"] >= float(df["marker_percentile"].quantile(0.85))
    return df


def render_one(root: Path, out_dir: Path, row: pd.Series, idx: int) -> Path:
    pair_id = str(row["pair_id"])
    target_type = str(row["target_type"])
    df = load_scene(root, row)
    norm = Normalize(0.0, 1.0)
    plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.55), dpi=220)
    draw_panel.pair_id = pair_id
    sc = draw_panel(axes[0], df, df["target_top15"], "A  Real ST target signal", "cyan outline: marker top15 region", norm)
    draw_panel(axes[1], df, df["blank"], "B  Stage3B blank result", "cyan outline: blanked unsupported region", norm)
    fig.suptitle(
        f"Stage3B reference-dropout spatial validation: {display_label(pair_id, target_type)}",
        x=0.03,
        ha="left",
        fontsize=12.2,
        fontweight="bold",
    )
    fig.subplots_adjust(left=0.04, right=0.91, top=0.78, bottom=0.15, wspace=-0.05)
    cax = fig.add_axes([0.928, 0.24, 0.014, 0.50])
    cb = fig.colorbar(sc, cax=cax)
    cb.set_label("target marker\npercentile", fontsize=6.8, labelpad=4)
    cb.ax.tick_params(labelsize=6.8)
    safe_target = target_type.replace("/", "_").replace(" ", "_")
    out = out_dir / f"scene_{idx:02d}_{pair_id}_{safe_target}.png"
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out


def stack_images(paths: list[Path], out_png: Path, title: str) -> None:
    images = [Image.open(p).convert("RGB") for p in paths]
    target_w = 2280
    resized = []
    for image in images:
        scale = target_w / image.size[0]
        resized.append(image.resize((target_w, int(round(image.size[1] * scale))), Image.Resampling.LANCZOS))
    outer = 30
    title_h = 70
    row_gap = 16
    canvas_w = target_w + outer * 2
    canvas_h = outer * 2 + title_h + sum(img.size[1] for img in resized) + row_gap * (len(resized) - 1)
    canvas = Image.new("RGB", (canvas_w, canvas_h), (255, 255, 255))
    draw = ImageDraw.Draw(canvas)
    font = load_font(34, bold=True)
    bb = draw.textbbox((0, 0), title, font=font)
    draw.text(((canvas_w - (bb[2] - bb[0])) // 2, outer + 2), title, fill=(24, 24, 24), font=font)
    y = outer + title_h
    for image in resized:
        canvas.paste(image, (outer, y))
        y += image.size[1] + row_gap
    canvas.save(out_png, optimize=True)


def render_compact_editable_svg(
    root: Path,
    final_dir: Path,
    metadata: pd.DataFrame,
    out_prefix: str,
) -> Path:
    """Render a hybrid SVG: editable typography with compact raster spot layers."""
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "svg.fonttype": "none",
            "svg.image_inline": True,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    norm = Normalize(0.0, 1.0)
    n_rows = len(metadata)
    fig = plt.figure(figsize=(9.2, 3.45 * n_rows + 0.75), dpi=240, facecolor="white")
    outer = fig.add_gridspec(
        n_rows,
        1,
        left=0.055,
        right=0.955,
        top=0.955,
        bottom=0.025,
        hspace=0.30,
    )

    for row_index, (_, row) in enumerate(metadata.iterrows()):
        pair_id = str(row["pair_id"])
        target_type = str(row["target_type"])
        df = load_scene(root, row)
        inner = outer[row_index, 0].subgridspec(
            2,
            3,
            height_ratios=[0.30, 0.70],
            width_ratios=[1.0, 1.0, 0.035],
            wspace=0.045,
            hspace=0.12,
        )
        title_ax = fig.add_subplot(inner[0, :2])
        title_ax.set_axis_off()
        title_ax.text(
            0.0,
            0.72,
            display_label(pair_id, target_type),
            ha="left",
            va="center",
            fontsize=10.5,
            fontweight="bold",
        )
        ax_left = fig.add_subplot(inner[1, 0])
        ax_right = fig.add_subplot(inner[1, 1])
        cax = fig.add_subplot(inner[1, 2])
        draw_panel.pair_id = pair_id
        sc = draw_panel(
            ax_left,
            df,
            df["target_top15"],
            "A  Real ST target signal",
            "cyan outline: marker top15 region",
            norm,
            compact_svg=True,
        )
        draw_panel(
            ax_right,
            df,
            df["blank"],
            "B  Stage3B blank result",
            "cyan outline: blanked unsupported region",
            norm,
            compact_svg=True,
        )
        colorbar = fig.colorbar(sc, cax=cax)
        colorbar.set_label("target marker\npercentile", fontsize=6.8, labelpad=4)
        colorbar.ax.tick_params(labelsize=6.8, length=2)

    fig.suptitle(
        "Stage3B Reference-dropout Spatial Validation",
        fontsize=14,
        fontweight="bold",
        y=0.992,
    )
    out_svg = final_dir / f"{out_prefix}.svg"
    fig.savefig(out_svg, format="svg", facecolor="white")
    plt.close(fig)
    return out_svg


def main() -> int:
    args = parse_args()
    root = Path(args.project_root).resolve()
    final_dir = root / args.out_dir
    scene_dir = final_dir / f"{args.out_prefix}_scenes"
    scene_dir.mkdir(parents=True, exist_ok=True)
    metadata = pd.read_csv(root / args.metadata_csv)
    pngs = [render_one(root, scene_dir, row, i + 1) for i, (_, row) in enumerate(metadata.iterrows())]
    out_png = final_dir / f"{args.out_prefix}.png"
    stack_images(pngs, out_png, "Stage3B Reference-dropout Spatial Validation")
    out_svg = render_compact_editable_svg(root, final_dir, metadata, args.out_prefix)
    print(f"[done] {out_png}")
    print(f"[done] {out_svg}")
    print(f"[done] scenes={len(pngs)} dir={scene_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
