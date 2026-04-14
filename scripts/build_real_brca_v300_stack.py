#!/usr/bin/env python
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont


def _load_font(size: int, bold: bool = False) -> ImageFont.ImageFont:
    candidates = []
    if bold:
        candidates.extend(
            [
                r"C:\Windows\Fonts\arialbd.ttf",
                r"C:\Windows\Fonts\segoeuib.ttf",
                r"C:\Windows\Fonts\calibrib.ttf",
            ]
        )
    candidates.extend(
        [
            r"C:\Windows\Fonts\arial.ttf",
            r"C:\Windows\Fonts\segoeui.ttf",
            r"C:\Windows\Fonts\calibri.ttf",
        ]
    )
    for p in candidates:
        try:
            return ImageFont.truetype(p, size=size)
        except Exception:
            continue
    return ImageFont.load_default()


def _find_segments(signal: np.ndarray, threshold: float, min_len: int) -> list[tuple[int, int]]:
    mask = signal > threshold
    out: list[tuple[int, int]] = []
    s = None
    for i, v in enumerate(mask):
        if v and s is None:
            s = i
        elif (not v) and s is not None:
            if i - s >= min_len:
                out.append((s, i - 1))
            s = None
    if s is not None and (len(mask) - s) >= min_len:
        out.append((s, len(mask) - 1))
    return out


def _extract_triptych(row_png: Path) -> list[Image.Image]:
    arr = np.array(Image.open(row_png).convert("RGB"))
    # Triptych figures are drawn on light gray axes/canvas.
    bg = np.array([230, 230, 230], dtype=np.int16)
    diff = np.max(np.abs(arr.astype(np.int16) - bg[None, None, :]), axis=2)
    fg = diff > 8

    x_counts = fg.sum(axis=0)
    x_segs = _find_segments(x_counts, threshold=60, min_len=800)
    if len(x_segs) < 3:
        raise RuntimeError(f"Cannot locate three panels in: {row_png}")
    panel_x = x_segs[:3]

    panels: list[Image.Image] = []
    for x0, x1 in panel_x:
        sub = fg[:, x0 : x1 + 1]
        y_counts = sub.sum(axis=1)
        y_segs = _find_segments(y_counts, threshold=80, min_len=400)
        if y_segs:
            y0, y1 = y_segs[0]
        else:
            # Fallback for sparse hex-grid panels (e.g. kidney): smooth row counts first.
            sm = np.convolve(y_counts.astype(np.float32), np.ones(7, dtype=np.float32) / 7.0, mode="same")
            y_segs2 = _find_segments(sm, threshold=60, min_len=220)
            if not y_segs2:
                y_segs2 = _find_segments(sm, threshold=35, min_len=180)
            if not y_segs2:
                raise RuntimeError(f"Cannot locate panel y-range in: {row_png}")
            # Use the widest segment to avoid clipping panel body.
            y0, y1 = max(y_segs2, key=lambda t: t[1] - t[0])
        panel = Image.fromarray(arr[y0 : y1 + 1, x0 : x1 + 1, :], mode="RGB")
        panels.append(_whiten_gray_background(panel))

    return panels


def _whiten_gray_background(img: Image.Image) -> Image.Image:
    arr = np.array(img.convert("RGB"))
    # Matplotlib canvas/axes gray in source triptychs is around (230,230,230).
    d = np.max(np.abs(arr.astype(np.int16) - np.array([230, 230, 230], dtype=np.int16)[None, None, :]), axis=2)
    arr[d <= 16] = np.array([255, 255, 255], dtype=np.uint8)
    return Image.fromarray(arr, mode="RGB")


def _load_type_order(truth_csv: Path) -> list[str]:
    with truth_csv.open("r", encoding="utf-8", newline="") as f:
        reader = csv.reader(f)
        header = next(reader, [])
    return [x for x in header if x and x != "spot_id"]


def _base_palette(type_order: list[str]) -> dict[str, str]:
    # Keep exactly the same palette sequence as visualize_missing_type_triptych.py
    base = [
        "#e41a1c",
        "#377eb8",
        "#4daf4a",
        "#984ea3",
        "#ff7f00",
        "#b8860b",
        "#a65628",
        "#f781bf",
        "#17becf",
        "#1b9e77",
        "#6a3d9a",
        "#b15928",
        "#00b4d8",
        "#fb5607",
        "#3a86ff",
    ]
    pal = {t: base[i % len(base)] for i, t in enumerate(type_order)}
    pal["__NoType__"] = "#e5e5e5"
    return pal


def _build_legend_image(type_order: list[str]) -> Image.Image:
    palette = _base_palette(type_order)
    labels = type_order + ["Unassigned spots"]
    colors = [palette[t] for t in type_order] + [palette["__NoType__"]]

    title_font = _load_font(36, bold=True)
    item_font = _load_font(24, bold=False)
    pad = 16
    dot = 12
    row_h = 30
    gap = 2
    title = "Type"

    dummy = Image.new("RGB", (10, 10), "white")
    d0 = ImageDraw.Draw(dummy)
    max_tw = 0
    for txt in labels:
        b = d0.textbbox((0, 0), txt, font=item_font)
        max_tw = max(max_tw, b[2] - b[0])
    tb = d0.textbbox((0, 0), title, font=title_font)
    title_w = tb[2] - tb[0]

    w = pad * 2 + dot + 10 + max(max_tw, title_w)
    h = pad * 2 + (tb[3] - tb[1]) + 8 + len(labels) * (row_h + gap)
    legend = Image.new("RGB", (w, h), (255, 255, 255))
    d = ImageDraw.Draw(legend)
    d.text((pad, pad), title, fill=(32, 32, 32), font=title_font)

    y = pad + (tb[3] - tb[1]) + 8
    for txt, col in zip(labels, colors):
        cx = pad + dot // 2
        cy = y + row_h // 2
        d.ellipse((cx - dot // 2, cy - dot // 2, cx + dot // 2, cy + dot // 2), fill=col, outline=(120, 120, 120), width=1)
        d.text((pad + dot + 10, y + 2), txt, fill=(28, 28, 28), font=item_font)
        y += row_h + gap
    return legend


def _build():
    root = Path(__file__).resolve().parents[1]
    vis_root = root / "visualizations" / "simulations" / "real_brca"
    out = vis_root / "V300_real_brca_triptych_stack_vertical.png"

    rows = [
        vis_root / "real_brca_clustered_sim" / "mapping_triptych_no_missing.png",
        vis_root / "real_brca_clustered_sim_mt_b_cells_d100_fill_endo" / "missing_type_triptych.png",
        vis_root / "real_brca_clustered_sim_mt_b_cells_d100_fill_endo_mt_pcs_d100_fill_endo" / "missing_type_triptych.png",
        vis_root / "real_brca_clustered_sim_mt_b_cells_d100_fill_endo_mt_pcs_d100_fill_endo_mt_pvl_d100_fill_endo" / "missing_type_triptych.png",
        vis_root / "real_brca_clustered_sim_mt_b_cells_d100_fill_endo_mt_pcs_d100_fill_endo_mt_pvl_d100_fill_endo_mt_nk_cells_d100_fill_endo" / "missing_type_triptych.png",
        vis_root / "real_brca_chain_mt_b_pcs_pvl_nk_tcd8_fill_endo" / "missing_type_triptych.png",
        vis_root / "real_brca_chain_mt_b_pcs_pvl_nk_tcd8_fib_fill_endo" / "missing_type_triptych.png",
    ]
    labels = [
        "No missing",
        "- B cells",
        "- PCs",
        "- PVL",
        "- NK cells",
        "- T cells CD8",
        "- Fibroblasts",
    ]
    for p in rows:
        if not p.exists():
            raise FileNotFoundError(p)

    truth_csv = root / "data" / "sim" / "real_brca" / "real_brca_clustered_sim" / "sim_truth_spot_type_fraction.csv"
    if not truth_csv.exists():
        raise FileNotFoundError(truth_csv)
    type_order = _load_type_order(truth_csv)
    legend = _build_legend_image(type_order)

    all_panels: list[list[Image.Image]] = []
    for p in rows:
        panels = _extract_triptych(p)
        all_panels.append(panels)

    target_w = max(img.size[0] for row in all_panels for img in row)
    target_h = max(img.size[1] for row in all_panels for img in row)
    all_panels = [[img.resize((target_w, target_h), Image.Resampling.LANCZOS) for img in row] for row in all_panels]

    left_w = 350
    right_w = legend.size[0] + 30
    top_title_h = 74
    top_header_h = 50
    col_gap = 28
    row_gap = 6
    outer = 20

    canvas_w = left_w + (target_w * 3) + (col_gap * 2) + right_w + outer * 2
    canvas_h = outer + top_title_h + top_header_h + (target_h * len(all_panels)) + row_gap * (len(all_panels) - 1) + outer
    canvas = Image.new("RGB", (canvas_w, canvas_h), (255, 255, 255))
    draw = ImageDraw.Draw(canvas)

    title_font = _load_font(38, bold=True)
    header_font = _load_font(26, bold=True)
    label_font = _load_font(42, bold=True)

    title = "Real BRCA: Progressive Missing-Type Mapping"
    tw = draw.textbbox((0, 0), title, font=title_font)[2]
    draw.text(((canvas_w - tw) // 2, outer + 6), title, fill=(24, 24, 24), font=title_font)

    x_start = outer + left_w
    y_start = outer + top_title_h + top_header_h
    col_titles = ["Truth", "CytoSPACE Baseline", "SVTuner + CytoSPACE (Route2)"]
    for c in range(3):
        cx = x_start + c * (target_w + col_gap)
        t = col_titles[c]
        cw = draw.textbbox((0, 0), t, font=header_font)[2]
        draw.text((cx + (target_w - cw) // 2, outer + top_title_h + 6), t, fill=(30, 30, 30), font=header_font)

    for r, row_imgs in enumerate(all_panels):
        y = y_start + r * (target_h + row_gap)
        for c, img in enumerate(row_imgs):
            x = x_start + c * (target_w + col_gap)
            canvas.paste(img, (x, y))
        label = labels[r]
        lb = draw.textbbox((0, 0), label, font=label_font)
        ly = y + (target_h - (lb[3] - lb[1])) // 2
        draw.text((outer + 8, ly), label, fill=(30, 30, 30), font=label_font)

    lx = x_start + 3 * target_w + 2 * col_gap + 18
    ly = y_start + 10
    canvas.paste(legend, (lx, ly))

    out.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(out, optimize=True)
    print(f"[OK] wrote {out}")
    print(f"[INFO] size={canvas_w}x{canvas_h} panel={target_w}x{target_h}")


if __name__ == "__main__":
    _build()
