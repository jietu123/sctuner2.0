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
    start = None
    for i, v in enumerate(mask):
        if v and start is None:
            start = i
        elif (not v) and start is not None:
            if i - start >= min_len:
                out.append((start, i - 1))
            start = None
    if start is not None and (len(mask) - start) >= min_len:
        out.append((start, len(mask) - 1))
    return out


def _whiten_gray_background(img: Image.Image) -> Image.Image:
    arr = np.array(img.convert("RGB"))
    bg = np.array([230, 230, 230], dtype=np.int16)
    d = np.max(np.abs(arr.astype(np.int16) - bg[None, None, :]), axis=2)
    arr[d <= 16] = np.array([255, 255, 255], dtype=np.uint8)
    return Image.fromarray(arr, mode="RGB")


def _extract_triptych(row_png: Path) -> list[Image.Image]:
    arr = np.array(Image.open(row_png).convert("RGB"))
    bg = np.array([230, 230, 230], dtype=np.int16)
    diff = np.max(np.abs(arr.astype(np.int16) - bg[None, None, :]), axis=2)
    fg = diff > 8

    x_counts = fg.sum(axis=0)
    x_segs = _find_segments(x_counts, threshold=60, min_len=800)
    if len(x_segs) < 3:
        raise RuntimeError(f"Cannot locate three panels in: {row_png}")

    panels: list[Image.Image] = []
    for x0, x1 in x_segs[:3]:
        sub = fg[:, x0 : x1 + 1]
        y_counts = sub.sum(axis=1)
        y_segs = _find_segments(y_counts, threshold=80, min_len=400)
        if y_segs:
            y0, y1 = y_segs[0]
        else:
            sm = np.convolve(y_counts.astype(np.float32), np.ones(7, dtype=np.float32) / 7.0, mode="same")
            y_segs2 = _find_segments(sm, threshold=60, min_len=220)
            if not y_segs2:
                y_segs2 = _find_segments(sm, threshold=35, min_len=180)
            if not y_segs2:
                raise RuntimeError(f"Cannot locate panel y-range in: {row_png}")
            y0, y1 = max(y_segs2, key=lambda t: t[1] - t[0])
        panel = Image.fromarray(arr[y0 : y1 + 1, x0 : x1 + 1, :], mode="RGB")
        panels.append(_whiten_gray_background(panel))
    return panels


def _load_type_order(truth_csv: Path) -> list[str]:
    with truth_csv.open("r", encoding="utf-8", newline="") as f:
        reader = csv.reader(f)
        header = next(reader, [])
    return [x for x in header if x and x != "spot_id"]


def _base_palette(type_order: list[str]) -> dict[str, str]:
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
    palette = {t: base[i % len(base)] for i, t in enumerate(type_order)}
    if "B cells" in palette and "Epithelial cells" in palette:
        palette["B cells"], palette["Epithelial cells"] = palette["Epithelial cells"], palette["B cells"]
    palette["__NoType__"] = "#e5e5e5"
    return palette


def _build_legend(type_order: list[str]) -> Image.Image:
    palette = _base_palette(type_order)
    labels = type_order + ["Unassigned spots"]
    colors = [palette[t] for t in type_order] + [palette["__NoType__"]]

    title_font = _load_font(28, bold=True)
    item_font = _load_font(18)
    pad = 14
    dot = 10
    row_h = 24
    title = "Type"

    dummy = Image.new("RGB", (10, 10), "white")
    draw = ImageDraw.Draw(dummy)
    max_tw = 0
    for txt in labels:
        b = draw.textbbox((0, 0), txt, font=item_font)
        max_tw = max(max_tw, b[2] - b[0])
    tb = draw.textbbox((0, 0), title, font=title_font)

    w = pad * 2 + dot + 10 + max(max_tw, tb[2] - tb[0])
    h = pad * 2 + (tb[3] - tb[1]) + 8 + len(labels) * row_h
    legend = Image.new("RGB", (w, h), (255, 255, 255))
    draw = ImageDraw.Draw(legend)
    draw.text((pad, pad), title, fill=(32, 32, 32), font=title_font)

    y = pad + (tb[3] - tb[1]) + 10
    for txt, col in zip(labels, colors):
        cx = pad + dot // 2
        cy = y + row_h // 2
        draw.ellipse((cx - dot // 2, cy - dot // 2, cx + dot // 2, cy + dot // 2), fill=col, outline=(120, 120, 120), width=1)
        draw.text((pad + dot + 10, y + 2), txt, fill=(28, 28, 28), font=item_font)
        y += row_h
    return legend


def _build() -> None:
    root = Path(__file__).resolve().parents[1]
    vis_root = root / "visualizations" / "simulations" / "real_brca"
    out = vis_root / "real_brca_triptych_stack_3rows.png"

    rows = [
        vis_root / "real_brca_clustered_sim" / "mapping_triptych_no_missing.png",
        vis_root / "real_brca_clustered_sim_missing_epithelial_cells" / "missing_type_triptych.png",
        vis_root / "real_brca_clustered_sim_missing_epithelial_monocytes_macrophages" / "missing_type_triptych.png",
    ]
    row_labels = [
        "No missing",
        "- Epithelial cells",
        "- Epithelial cells\n- Monocytes & Macrophages",
    ]
    col_titles = [
        "Truth (no missing type)",
        "CytoSPACE Baseline Mapping",
        "SVTuner + CytoSPACE Mapping (Route2)",
    ]

    for p in rows:
        if not p.exists():
            raise FileNotFoundError(p)

    truth_csv = root / "data" / "sim" / "real_brca" / "real_brca_clustered_sim" / "sim_truth_spot_type_fraction.csv"
    if not truth_csv.exists():
        raise FileNotFoundError(truth_csv)

    type_order = _load_type_order(truth_csv)
    legend = _build_legend(type_order)
    all_panels = [_extract_triptych(p) for p in rows]

    target_w = max(img.size[0] for row in all_panels for img in row)
    target_h = max(img.size[1] for row in all_panels for img in row)
    all_panels = [[img.resize((target_w, target_h), Image.Resampling.LANCZOS) for img in row] for row in all_panels]

    left_w = 330
    right_w = legend.size[0] + 36
    top_title_h = 68
    top_header_h = 44
    col_gap = 24
    row_gap = 22
    outer = 24

    canvas_w = left_w + 3 * target_w + 2 * col_gap + right_w + outer * 2
    canvas_h = outer + top_title_h + top_header_h + len(all_panels) * target_h + (len(all_panels) - 1) * row_gap + outer
    canvas = Image.new("RGB", (canvas_w, canvas_h), (255, 255, 255))
    draw = ImageDraw.Draw(canvas)

    title_font = _load_font(36, bold=True)
    header_font = _load_font(22, bold=True)
    row_font = _load_font(28, bold=True)

    title = "Real BRCA Mapping Comparison"
    tb = draw.textbbox((0, 0), title, font=title_font)
    draw.text(((canvas_w - (tb[2] - tb[0])) // 2, outer + 4), title, fill=(24, 24, 24), font=title_font)

    x0 = outer + left_w
    y0 = outer + top_title_h + top_header_h

    for i, txt in enumerate(col_titles):
        cx = x0 + i * (target_w + col_gap)
        bb = draw.textbbox((0, 0), txt, font=header_font)
        draw.text((cx + (target_w - (bb[2] - bb[0])) // 2, outer + top_title_h), txt, fill=(30, 30, 30), font=header_font)

    for r, panels in enumerate(all_panels):
        y = y0 + r * (target_h + row_gap)
        label = row_labels[r]
        bb = draw.multiline_textbbox((0, 0), label, font=row_font, spacing=6)
        draw.multiline_text(
            (outer + 10, y + (target_h - (bb[3] - bb[1])) // 2),
            label,
            fill=(30, 30, 30),
            font=row_font,
            spacing=6,
        )
        for c, panel in enumerate(panels):
            x = x0 + c * (target_w + col_gap)
            canvas.paste(panel, (x, y))

    legend_x = x0 + 3 * target_w + 2 * col_gap + 18
    legend_y = y0 + 8
    canvas.paste(legend, (legend_x, legend_y))

    out.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(out, optimize=True)
    print(f"[OK] wrote {out}")
    print(f"[INFO] canvas={canvas_w}x{canvas_h} panel={target_w}x{target_h}")


if __name__ == "__main__":
    _build()
