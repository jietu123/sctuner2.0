#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

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


def _build() -> None:
    root = Path(__file__).resolve().parents[1]
    rows = [
        (
            "Real BRCA",
            root / "visualizations" / "simulations" / "real_brca" / "real_brca_triptych_stack_3rows.png",
            1.0,
        ),
        (
            "Mouse brain refined",
            root / "visualizations" / "simulations" / "mouse_brain_refined" / "mouse_brain_refined_triptych_stack_3rows.png",
            0.5,
        ),
        (
            "Human lung 5loc",
            root / "visualizations" / "simulations" / "human_lung_5loc" / "human_lung_5loc_triptych_stack_3rows.png",
            1.0,
        ),
    ]
    out = root / "visualizations" / "simulations" / "simulation_triptych_overview_stack_3datasets.png"

    images = []
    for label, path, y_scale in rows:
        if not path.exists():
            raise FileNotFoundError(path)
        img = Image.open(path).convert("RGB")
        if y_scale != 1.0:
            new_h = max(1, int(round(img.size[1] * y_scale)))
            img = img.resize((img.size[0], new_h), Image.Resampling.LANCZOS)
        images.append((label, img))

    max_w = max(img.size[0] for _, img in images)
    left_w = 330
    outer = 24
    row_gap = 28
    top_title_h = 72

    total_h = outer * 2 + top_title_h + sum(img.size[1] for _, img in images) + row_gap * (len(images) - 1)
    total_w = outer * 2 + left_w + max_w
    canvas = Image.new("RGB", (total_w, total_h), (255, 255, 255))
    draw = ImageDraw.Draw(canvas)

    title_font = _load_font(38, bold=True)
    row_font = _load_font(28, bold=True)

    title = "Simulation Mapping Visualization Overview"
    tb = draw.textbbox((0, 0), title, font=title_font)
    draw.text(((total_w - (tb[2] - tb[0])) // 2, outer + 4), title, fill=(24, 24, 24), font=title_font)

    x_img = outer + left_w
    y = outer + top_title_h
    for label, img in images:
        label_box = draw.multiline_textbbox((0, 0), label, font=row_font, spacing=6)
        draw.multiline_text(
            (outer + 8, y + (img.size[1] - (label_box[3] - label_box[1])) // 2),
            label,
            fill=(30, 30, 30),
            font=row_font,
            spacing=6,
        )
        x = x_img + (max_w - img.size[0]) // 2
        canvas.paste(img, (x, y))
        y += img.size[1] + row_gap

    out.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(out, optimize=True)
    print(f"[OK] wrote {out}")
    print(f"[INFO] canvas={total_w}x{total_h}")


if __name__ == "__main__":
    _build()
