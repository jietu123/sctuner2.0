#!/usr/bin/env python
from __future__ import annotations

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


def _crop_panel_grid(png_path: Path) -> list[Image.Image]:
    arr = np.array(Image.open(png_path).convert("RGB"))
    bg = np.array([255, 255, 255], dtype=np.int16)
    diff = np.max(np.abs(arr.astype(np.int16) - bg[None, None, :]), axis=2)
    fg = diff > 12

    x_counts = fg.sum(axis=0)
    x_segments = _find_segments(x_counts, threshold=40, min_len=80)
    if len(x_segments) < 4:
        raise RuntimeError(f"Cannot locate four panels in: {png_path}")

    panels: list[Image.Image] = []
    for x0, x1 in x_segments[:4]:
        sub = fg[:, x0 : x1 + 1]
        y_counts = sub.sum(axis=1)
        y_segments = _find_segments(y_counts, threshold=40, min_len=120)
        if not y_segments:
            raise RuntimeError(f"Cannot locate panel y-range in: {png_path}")
        y0, y1 = max(y_segments, key=lambda t: t[1] - t[0])
        panels.append(Image.fromarray(arr[y0 : y1 + 1, x0 : x1 + 1, :], mode="RGB"))
    return panels


def _label_for_dir(name: str) -> str:
    mapping = {
        "adult_mouse_kidney_real_profile_mask_endo": "Adult mouse kidney\nMask: Endo",
        "ffpe_mouse_brain_sagittal_real_profile_mask_microglia": "FFPE mouse brain\nMask: Microglia",
        "human_breast_cancer_real_profile_mask_basal_cell": "Human breast cancer\nMask: Basal cell",
        "human_breast_cancer_visium_ff_wta_real_profile_mask_macrophage": "Human breast cancer FFPE\nMask: Macrophage",
        "human_breast_cancer_wta_120_real_profile_mask_endothelial_cell": "Human breast cancer WTA120\nMask: Endothelial cell",
        "human_cervical_cancer_real_profile_mask_epithelial_cell": "Human cervical cancer\nMask: Epithelial cell",
        "human_heart_ff_real_profile_mask_endothelial_cell": "Human heart FF\nMask: Endothelial cell",
        "human_intestine_cancer_real_profile_mask_endothelial_cell": "Human intestine cancer\nMask: Endothelial cell",
        "human_lymph_node_real_profile_mask_b_cell": "Human lymph node\nMask: B cell",
        "mouse_embryo_real_profile_mask_erythroid": "Mouse embryo\nMask: Erythroid",
    }
    return mapping.get(name, name.replace("_", " "))


def _build() -> None:
    root = Path(__file__).resolve().parents[1]
    vis_root = root / "visualizations" / "masked_scenarios" / "real"
    out = root / "visualizations" / "masked_scenarios" / "masked_scenarios_real_stack_10x4.png"

    pngs = sorted(vis_root.glob("*/*_mask_ABCD_expression_identical_cyan_outline.png"))
    if not pngs:
        raise FileNotFoundError(f"No ABCD masked-scenario PNGs under {vis_root}")

    rows = []
    for png in pngs:
        rows.append((_label_for_dir(png.parent.name), _crop_panel_grid(png)))

    target_w = max(panel.size[0] for _, panels in rows for panel in panels)
    target_h = max(panel.size[1] for _, panels in rows for panel in panels)
    resized_rows = []
    for label, panels in rows:
        resized_rows.append(
            (
                label,
                [panel.resize((target_w, target_h), Image.Resampling.LANCZOS) for panel in panels],
            )
        )

    col_titles = [
        "Original real ST",
        "Profile-mask real ST",
        "CytoSPACE Baseline Mapping",
        "SVTuner + CytoSPACE Mapping (Route2)",
    ]

    outer = 28
    left_w = 360
    top_title_h = 72
    top_header_h = 48
    col_gap = 24
    row_gap = 18

    canvas_w = outer * 2 + left_w + 4 * target_w + 3 * col_gap
    canvas_h = outer * 2 + top_title_h + top_header_h + len(resized_rows) * target_h + (len(resized_rows) - 1) * row_gap
    canvas = Image.new("RGB", (canvas_w, canvas_h), (255, 255, 255))
    draw = ImageDraw.Draw(canvas)

    title_font = _load_font(38, bold=True)
    header_font = _load_font(22, bold=True)
    row_font = _load_font(24, bold=True)

    title = "Real Masked-scenario Visualization Overview"
    tb = draw.textbbox((0, 0), title, font=title_font)
    draw.text(((canvas_w - (tb[2] - tb[0])) // 2, outer + 4), title, fill=(24, 24, 24), font=title_font)

    x0 = outer + left_w
    y0 = outer + top_title_h + top_header_h

    for i, txt in enumerate(col_titles):
        cx = x0 + i * (target_w + col_gap)
        bb = draw.textbbox((0, 0), txt, font=header_font)
        draw.multiline_text(
            (cx + (target_w - (bb[2] - bb[0])) // 2, outer + top_title_h),
            txt,
            fill=(30, 30, 30),
            font=header_font,
            align="center",
        )

    for r, (label, panels) in enumerate(resized_rows):
        y = y0 + r * (target_h + row_gap)
        bb = draw.multiline_textbbox((0, 0), label, font=row_font, spacing=5)
        draw.multiline_text(
            (outer + 8, y + (target_h - (bb[3] - bb[1])) // 2),
            label,
            fill=(30, 30, 30),
            font=row_font,
            spacing=5,
        )
        for c, panel in enumerate(panels):
            x = x0 + c * (target_w + col_gap)
            canvas.paste(panel, (x, y))

    out.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(out, optimize=True)
    print(f"[OK] wrote {out}")
    print(f"[INFO] rows={len(resized_rows)} canvas={canvas_w}x{canvas_h} panel={target_w}x{target_h}")


if __name__ == "__main__":
    _build()
