from __future__ import annotations

import math
import random
from pathlib import Path

from rebuild_editable_svg_references import SVG, diverging_color, gaussian_clusters


ROOT = Path(__file__).resolve().parents[1]

BLUE = "#074bb5"
BLUE_LIGHT = "#edf5ff"
BLUE_LINE = "#2f69b8"
NAVY = "#07398e"
ORANGE = "#ff6200"
ORANGE_LIGHT = "#fff7ef"
PURPLE = "#6518aa"
PURPLE_LIGHT = "#f8f3fc"
GREEN = "#197b1e"
GREEN_LIGHT = "#f4fbf1"
CHARCOAL = "#30343b"
MUTED = "#68717c"
GRID = "#d5d9de"
TYPE_COLORS = ["#ff4f79", "#69bd28", "#7137c4", "#2d7de0", "#ff8514", "#929292"]


def begin_layer(svg: SVG, layer_id: str, label: str) -> None:
    svg.start(
        "g",
        id=layer_id,
        inkscape_groupmode="layer",
        inkscape_label=label,
    )


def end_layer(svg: SVG) -> None:
    svg.end("g")


def add_defs(svg: SVG) -> None:
    svg.raw(
        """<defs>
  <marker id="arrow-blue" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">
    <path d="M 0 0 L 10 5 L 0 10 z" fill="#074bb5"/>
  </marker>
  <marker id="arrow-purple" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">
    <path d="M 0 0 L 10 5 L 0 10 z" fill="#6518aa"/>
  </marker>
  <marker id="arrow-grey" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">
    <path d="M 0 0 L 10 5 L 0 10 z" fill="#4b5563"/>
  </marker>
</defs>"""
    )


def line_arrow(svg: SVG, x1, y1, x2, y2, color=BLUE, width=2.4, marker="arrow-blue", **attrs) -> None:
    args = dict(x1=x1, y1=y1, x2=x2, y2=y2, stroke=color, stroke_width=width, fill="none", marker_end=f"url(#{marker})")
    args.update(attrs)
    svg.el("line", **args)


def path_arrow(svg: SVG, d: str, color=BLUE, width=2.1, marker="arrow-blue", **attrs) -> None:
    args = dict(d=d, stroke=color, stroke_width=width, fill="none", marker_end=f"url(#{marker})")
    args.update(attrs)
    svg.el("path", **args)


def rounded_box(svg: SVG, x, y, w, h, stroke=BLUE_LINE, fill="#ffffff", radius=14, width=1.25, **attrs) -> None:
    args = dict(x=x, y=y, width=w, height=h, rx=radius, fill=fill, stroke=stroke, stroke_width=width)
    args.update(attrs)
    svg.el("rect", **args)


def footer_box(svg: SVG, x, y, w, h, lines, stroke=BLUE_LINE, fill=BLUE_LIGHT, color=BLUE, size=14) -> None:
    rounded_box(svg, x, y, w, h, stroke=stroke, fill=fill, radius=8, width=1)
    start_y = y + h / 2 - (len(lines) - 1) * size * .58 + size * .35
    svg.multiline(x + w / 2, start_y, lines, size=size, weight=600, anchor="middle", fill=color, line_height=1.25)


def badge(svg: SVG, x, y, label, fill=BLUE, radius=20, size=22) -> None:
    svg.el("circle", cx=x, cy=y, r=radius, fill=fill)
    svg.text(x, y + size * .34, label, size=size, weight=700, anchor="middle", fill="#ffffff")


def stage_title(svg: SVG, x, y, lines, color=BLUE, size=19) -> None:
    svg.multiline(x, y, lines, size=size, weight=700, anchor="middle", fill=color, line_height=1.2)


def axes(svg: SVG, x, y, scale=1.0) -> None:
    line_arrow(svg, x, y, x + 43 * scale, y, color="#111111", width=1.7, marker="arrow-grey")
    line_arrow(svg, x, y, x, y - 48 * scale, color="#111111", width=1.7, marker="arrow-grey")
    svg.text(x + 48 * scale, y + 7, "x", size=12 * scale, weight=600)
    svg.text(x - 10, y - 52 * scale, "y", size=12 * scale, weight=600)


def heatmap(svg: SVG, x, y, rows=7, cols=7, cell=14, seed=1, mode="diagonal", border="#59616d") -> None:
    rng = random.Random(seed)
    for r in range(rows):
        for c in range(cols):
            if mode == "diagonal":
                val = max(-1, min(1, .85 - abs(r - c) * .37 + rng.uniform(-.55, .35)))
            elif mode == "residual":
                val = max(-1, min(1, rng.gauss(0, .73) + (.5 if (r + 2 * c) % 7 == 0 else 0)))
            else:
                val = rng.uniform(-1, 1)
            svg.el(
                "rect",
                x=x + c * cell,
                y=y + r * cell,
                width=cell,
                height=cell,
                fill=diverging_color(val),
                stroke="#ffffff",
                stroke_width=.4,
            )
    svg.el("rect", x=x, y=y, width=cols * cell, height=rows * cell, fill="none", stroke=border, stroke_width=1)


def cluster_reference(svg: SVG, centers, scale=1.0, seed=1, points=42) -> None:
    gaussian_clusters(
        svg,
        centers,
        TYPE_COLORS[: len(centers)],
        points_per=points,
        sx=15 * scale,
        sy=16 * scale,
        radius=3.0 * scale,
        seed=seed,
    )


def tissue_points(svg: SVG, x, y, w, h, seed=1, faded=False, mode="purple", dot=3.2, spacing=7.0) -> None:
    rng = random.Random(seed)
    points = []
    row = 0
    yy = y
    while yy <= y + h:
        xx = x + (spacing / 2 if row % 2 else 0)
        while xx <= x + w:
            nx = (xx - (x + .51 * w)) / (.50 * w)
            ny = (yy - (y + .49 * h)) / (.47 * h)
            edge = nx * nx + ny * ny
            top_wave = .11 * math.sin((nx + 1) * 4.4) - .07 * math.cos(nx * 8)
            notch = xx > x + .70 * w and yy > y + .61 * h and yy < y + .82 * h
            lower_cut = xx < x + .18 * w and yy > y + .80 * h
            if edge < 1 + top_wave and not notch and not lower_cut:
                points.append((xx, yy, nx, ny))
            xx += spacing
        yy += spacing * .866
        row += 1

    svg.start("g", opacity=.44 if faded else 1)
    for px, py, nx, ny in points:
        if mode == "spatial":
            score = [
                -1.2 * nx - .4 * ny,
                nx - .8 * ny,
                -nx + 1.25 * ny,
                nx + .8 * ny,
                math.sin(nx * 4) + math.cos(ny * 3),
                -abs(nx) + .65 * ny,
            ]
            color = TYPE_COLORS[max(range(6), key=lambda i: score[i])]
        elif mode == "mixed":
            field = math.sin(px * .08) + math.cos(py * .07) + rng.gauss(0, .9)
            color = diverging_color(max(-1, min(1, field / 2.3)))
        else:
            t = max(0, min(1, .35 + .26 * math.sin(nx * 3) + .25 * rng.random()))
            color = diverging_color(-.9 + 1.35 * t)
        svg.el("circle", cx=f"{px:.2f}", cy=f"{py:.2f}", r=dot, fill=color, stroke="#ffffff", stroke_width=.35)
    svg.end("g")


def spot_grid(svg: SVG, cx, cy, radius, spacing=14, dot=5.5, filled=False, seed=1) -> None:
    rng = random.Random(seed)
    row = 0
    yy = cy - radius
    while yy <= cy + radius:
        xx = cx - radius + (spacing / 2 if row % 2 else 0)
        while xx <= cx + radius:
            if (xx - cx) ** 2 + (yy - cy) ** 2 <= radius ** 2:
                fill = TYPE_COLORS[rng.randrange(5)] if filled else "#ffffff"
                svg.el("circle", cx=xx, cy=yy, r=dot, fill=fill, stroke="#6d5bc5", stroke_width=1.25)
            xx += spacing
        yy += spacing * .866
        row += 1


def colorful_assignment(svg: SVG, x, y, w, rows=3, cols=15, seed=1, r=5.5) -> None:
    rng = random.Random(seed)
    dx = w / max(1, cols - 1)
    for row in range(rows):
        for col in range(cols):
            px = x + col * dx + (row % 2) * dx * .22
            py = y + row * r * 1.7
            svg.el("circle", cx=px, cy=py, r=r, fill=TYPE_COLORS[rng.randrange(5)], stroke="#ffffff", stroke_width=.6)


def cell_spot_icons(svg: SVG, x, y, cols=5, scale=1.0) -> None:
    for col in range(cols):
        cx = x + col * 42 * scale
        svg.el("path", d=f"M {cx-3} {y-18*scale} q {-8*scale} {7*scale} 0 {15*scale} q {7*scale} {6*scale} {-2*scale} {14*scale}", fill="none", stroke="#111111", stroke_width=1.4)
        svg.el("circle", cx=cx, cy=y + 22 * scale, r=6 * scale, fill="none", stroke=PURPLE, stroke_width=1.7)
        svg.el("line", x1=cx - 4 * scale, y1=y + 26 * scale, x2=cx + 4 * scale, y2=y + 18 * scale, stroke=PURPLE, stroke_width=1.2)


def mini_stacked_bars(svg: SVG, x, y, width=145, height=145, bars=6) -> None:
    rng = random.Random(44)
    bw = width / (bars * 1.45)
    gap = (width - bw * bars) / (bars - 1)
    for i in range(bars):
        vals = [rng.uniform(.05, .4) for _ in range(5)]
        total = sum(vals)
        yy = y + height
        for value, color in zip(vals, TYPE_COLORS[:5]):
            h = height * value / total
            yy -= h
            svg.el("rect", x=x + i * (bw + gap), y=yy, width=bw, height=h, fill=color, stroke="#ffffff", stroke_width=.4)
    svg.el("line", x1=x, y1=y, x2=x, y2=y + height, stroke="#111111", stroke_width=1)
    svg.el("line", x1=x, y1=y + height, x2=x + width, y2=y + height, stroke="#111111", stroke_width=1)
    for frac, label in [(0, "1.0"), (.5, "0.5"), (1, "0")]:
        yy = y + height * frac
        svg.text(x - 8, yy + 5, label, size=13, anchor="end")


def draw_input_panel(svg: SVG, x, y, w, h, compact=False) -> None:
    rounded_box(svg, x, y, w, h, stroke=BLUE_LINE, radius=18, width=1.4)
    badge(svg, x + 34, y + 41, "1", radius=20 if not compact else 19)
    stage_title(svg, x + w * .61, y + 48, ["Input data"], size=18 if compact else 20)
    svg.multiline(x + w / 2, y + 115, ["Spatial transcriptomics", "expression + coordinates"], size=15 if compact else 16, weight=600, anchor="middle", line_height=1.35)
    tissue_points(svg, x + 40, y + 166, w - 75, 160 if not compact else 145, seed=11, mode="purple", dot=3.0, spacing=7.1)
    axes(svg, x + 30, y + 335 if not compact else y + 315, .78)
    svg.text(x + w / 2, y + 390 if not compact else y + 365, "+", size=34, weight=500, anchor="middle", fill="#111111")
    svg.multiline(x + w / 2, y + 432 if not compact else y + 405, ["scRNA–seq reference", "expression + cell labels"], size=15 if compact else 16, weight=600, anchor="middle", line_height=1.35)
    centers = [
        (x + 68, y + 530), (x + 145, y + 515), (x + 112, y + 565),
        (x + 65, y + 625), (x + 113, y + 625), (x + 178, y + 595),
    ]
    if compact:
        centers = [(cx, cy - 45) for cx, cy in centers]
    cluster_reference(svg, centers, scale=.85, seed=22, points=30)
    footer_box(svg, x + 15, y + h - 102, w - 30, 70, ["Matched sample universe", "(spot set & coordinate system)"], size=13)


def draw_stage1_panel(svg: SVG, x, y, w, h, compact=False) -> None:
    rounded_box(svg, x, y, w, h, stroke=BLUE_LINE, radius=18, width=1.4)
    badge(svg, x + 34, y + 41, "2", radius=20 if not compact else 19)
    stage_title(svg, x + w * .60, y + 45, ["Stage 1", "harmonization"], size=17 if compact else 19)
    labels = ["QC and annotation", "Gene–ID matching", "Shared–gene matrices", "Coordinate alignment"]
    by = y + 107
    bh = 47 if not compact else 40
    gap = 12 if not compact else 10
    for label in labels:
        rounded_box(svg, x + 19, by, w - 38, bh, stroke="#99bce3", fill=BLUE_LIGHT, radius=7, width=.9)
        svg.text(x + w / 2, by + bh * .62, label, size=14 if compact else 15, weight=500, anchor="middle", fill="#111111")
        by += bh + gap
    svg.el("line", x1=x + 19, y1=by + 5, x2=x + w - 19, y2=by + 5, stroke="#75a8dd", stroke_width=1)
    svg.text(x + w / 2, by + 50, "Unified exported inputs", size=15 if compact else 17, weight=700, anchor="middle", fill=BLUE)
    heatmap(svg, x + w / 2 - 61, by + 80, rows=6, cols=6, cell=20, seed=30, mode="residual")
    files = ["sc_expression.csv", "st_expression.csv", "sc_metadata.csv", "st_coordinates.csv", "(ann_truth_*.csv)"]
    fy = by + 225
    for name in files:
        svg.text(x + w / 2, fy, name, size=12 if compact else 13, weight=500, anchor="middle")
        fy += 25 if not compact else 22
    footer_box(svg, x + 15, y + h - 102, w - 30, 70, ["Ready for diagnosis", "and matched mapping"], size=13)


def draw_simple_stage3(svg: SVG, x, y, w, h) -> None:
    rounded_box(svg, x, y, w, h, stroke="#f0a15c", radius=18, width=1.4)
    badge(svg, x + 66, y + 41, "3", fill=ORANGE, radius=20)
    stage_title(svg, x + w * .61, y + 44, ["Stage 3A / 3B", "mismatch diagnosis"], color=ORANGE, size=19)
    mid = x + w / 2
    svg.el("line", x1=mid, y1=y + 130, x2=mid, y2=y + 388, stroke="#b6b6b6", stroke_width=1.3, stroke_dasharray="8 7")
    svg.multiline(x + w * .26, y + 145, ["3A: SC–only", "mismatch", "(ref has, ST lacks)"], size=16, weight=700, anchor="middle", fill=ORANGE, line_height=1.35)
    svg.multiline(x + w * .75, y + 145, ["3B: ST–only", "mismatch", "(ST has, ref lacks)"], size=16, weight=700, anchor="middle", fill=PURPLE, line_height=1.35)
    tissue_points(svg, x + 28, y + 224, w * .40, 125, seed=41, mode="purple", dot=2.7, spacing=6.1)
    tissue_points(svg, mid + 26, y + 224, w * .40, 125, seed=42, mode="mixed", faded=True, dot=2.7, spacing=6.1)
    line_arrow(svg, x + w * .25, y + 364, x + w * .25, y + 408, width=2.1)
    line_arrow(svg, x + w * .75, y + 364, x + w * .75, y + 408, width=2.1)
    svg.text(mid, y + 449, "Type–level support evaluation", size=16, weight=700, anchor="middle")
    heatmap(svg, x + 30, y + 482, rows=5, cols=5, cell=20, seed=51, mode="diagonal")
    cluster_reference(svg, [(x + 225, y + 497), (x + 281, y + 496), (x + 230, y + 554), (x + 286, y + 554)], scale=.70, seed=53, points=25)
    svg.el("ellipse", cx=x + 225, cy=y + 497, rx=28, ry=31, fill="none", stroke=PURPLE, stroke_width=1.2, stroke_dasharray="6 5")
    svg.el("ellipse", cx=x + 281, cy=y + 496, rx=28, ry=31, fill="none", stroke="#53a739", stroke_width=1.2, stroke_dasharray="6 5")
    footer_box(svg, x + 15, y + h - 151, w - 30, 119, ["Outputs:", "•  plugin_type / cleaned reference", "•  unsupported score / masks"], stroke="#f0a15c", fill=ORANGE_LIGHT, color="#a51f11", size=14)


def draw_stage4_panel(svg: SVG, x, y, w, h) -> None:
    rounded_box(svg, x, y, w, h, stroke="#5c6168", radius=18, width=1.3)
    badge(svg, x + 33, y + 41, "4", fill=CHARCOAL, radius=20)
    stage_title(svg, x + w * .60, y + 45, ["Stage 4", "CytoSPACE mapping"], color="#111111", size=18)
    svg.multiline(x + w / 2, y + 143, ["Same mapping backend for", "baseline and SVTuner"], size=15, weight=600, anchor="middle", line_height=1.35)
    rounded_box(svg, x + 16, y + 239, w - 32, 215, stroke="#9ea2a7", fill="#f8f8f8", radius=9, width=1)
    svg.text(x + w / 2, y + 282, "Matched inputs", size=15, weight=700, anchor="middle")
    svg.multiline(x + w / 2, y + 331, ["ST expression + coordinates", "(common genes)", "+", "Reference cells", "(baseline / cleaned)"], size=14, weight=500, anchor="middle", line_height=1.55)
    svg.multiline(x + w / 2, y + 520, ["CytoSPACE optimization", "(cell–spot assignment)"], size=15, weight=600, anchor="middle", line_height=1.3)
    cell_spot_icons(svg, x + 38, y + 619, cols=5, scale=.83)


def draw_simple_stage5(svg: SVG, x, y, w, h) -> None:
    rounded_box(svg, x, y, w, h, stroke=BLUE_LINE, radius=18, width=1.4)
    badge(svg, x + 34, y + 41, "5", radius=20)
    stage_title(svg, x + w * .60, y + 45, ["Mapping outputs and", "abstention–aware results"], size=18)
    svg.text(x + 28, y + 112, "Baseline route", size=15, weight=700, fill=BLUE)
    svg.text(x + w / 2, y + 147, "Original reference", size=14, weight=600, anchor="middle")
    line_arrow(svg, x + w / 2, y + 160, x + w / 2, y + 196, width=1.8)
    svg.multiline(x + w / 2, y + 223, ["CytoSPACE assignment", "(forced at every spot)"], size=13, weight=500, anchor="middle", line_height=1.25)
    colorful_assignment(svg, x + 52, y + 267, w - 105, rows=3, cols=15, seed=61, r=5.5)
    svg.el("line", x1=x + 18, y1=y + 335, x2=x + w - 18, y2=y + 335, stroke="#73a8dd", stroke_width=1.3, stroke_dasharray="7 6")
    svg.text(x + 28, y + 374, "SVTuner route", size=15, weight=700, fill=BLUE)
    svg.text(x + w / 2, y + 408, "Stage3A cleaned pool", size=13, weight=600, anchor="middle")
    line_arrow(svg, x + w / 2, y + 420, x + w / 2, y + 455, width=1.8)
    svg.multiline(x + w / 2, y + 481, ["CytoSPACE assignment", "(respect same backend)"], size=13, weight=500, anchor="middle", line_height=1.25)
    svg.text(x + w / 2, y + 548, "Stage3B withheld mask", size=13, weight=600, anchor="middle")
    for i in range(9):
        cx = x + 49 + i * (w - 98) / 8
        svg.el("circle", cx=cx, cy=y + 580, r=5.5, fill="none", stroke=PURPLE, stroke_width=1.4)
        svg.el("line", x1=cx - 4, y1=y + 584, x2=cx + 4, y2=y + 576, stroke=PURPLE, stroke_width=1.2)
    line_arrow(svg, x + w / 2, y + 596, x + w / 2, y + 630, width=1.8)
    svg.text(x + w / 2, y + 662, "Withheld–aware output", size=13, weight=600, anchor="middle")
    colorful_assignment(svg, x + 52, y + 687, w - 105, rows=3, cols=15, seed=62, r=5.5)
    footer_box(svg, x + 15, y + h - 102, w - 30, 70, ["Downstream analysis", "(with abstention awareness)"], size=13)


def build_a_simple() -> None:
    svg = SVG(1536, 1024, "SVTuner workflow – concise editable reconstruction")
    svg.el("rect", x=0, y=0, width=1536, height=1024, fill="#ffffff")
    add_defs(svg)

    begin_layer(svg, "top-brackets", "Top workflow brackets")
    svg.text(598, 63, "Reference–spatial compatibility diagnosis (pre–mapping)", size=20, weight=700, anchor="middle", fill=BLUE)
    svg.text(1190, 63, "Matched mapping and abstention–aware output", size=20, weight=700, anchor="middle", fill=BLUE)
    svg.el("line", x1=304, y1=88, x2=884, y2=88, stroke=BLUE, stroke_width=2)
    for xx in (304, 688, 884):
        line_arrow(svg, xx, 88, xx, 107, width=2)
    svg.el("line", x1=958, y1=88, x2=1404, y2=88, stroke=BLUE, stroke_width=2)
    for xx in (958, 1404):
        line_arrow(svg, xx, 88, xx, 107, width=2)
    end_layer(svg)

    y, h = 116, 808
    begin_layer(svg, "stage-1-input", "1 Input data")
    draw_input_panel(svg, 23, y, 236, h)
    end_layer(svg)
    begin_layer(svg, "stage-2-harmonization", "2 Stage 1 harmonization")
    draw_stage1_panel(svg, 282, y, 215, h)
    end_layer(svg)
    begin_layer(svg, "stage-3-diagnosis", "3 Stage 3A and 3B mismatch diagnosis")
    draw_simple_stage3(svg, 521, y, 338, h)
    end_layer(svg)
    begin_layer(svg, "stage-4-mapping", "4 CytoSPACE mapping")
    draw_stage4_panel(svg, 887, y, 254, h)
    end_layer(svg)
    begin_layer(svg, "stage-5-results", "5 Mapping outputs")
    draw_simple_stage5(svg, 1166, y, 285, h)
    end_layer(svg)

    begin_layer(svg, "workflow-arrows", "Workflow arrows")
    for x1, x2 in [(259, 280), (497, 519), (859, 885), (1141, 1164)]:
        line_arrow(svg, x1, 480, x2, 480, width=3.4)
    end_layer(svg)
    svg.finish(ROOT / "A简洁版_editable.svg")


def build_b_simulation() -> None:
    svg = SVG(1536, 1024, "Simulated data generation – editable reconstruction")
    svg.el("rect", x=0, y=0, width=1536, height=1024, fill="#ffffff")
    add_defs(svg)

    begin_layer(svg, "title", "Title")
    rounded_box(svg, 459, 20, 524, 65, stroke=BLUE_LINE, fill=BLUE_LIGHT, radius=10, width=1.2)
    svg.text(721, 64, "Simulated data generation", size=32, weight=700, anchor="middle", fill=NAVY)
    end_layer(svg)

    begin_layer(svg, "inputs-panel", "Inputs")
    rounded_box(svg, 81, 105, 1276, 337, stroke=BLUE_LINE, radius=10, width=1.3)
    svg.el("rect", x=82, y=106, width=1274, height=39, rx=9, fill=BLUE_LIGHT)
    svg.el("line", x1=81, y1=145, x2=1357, y2=145, stroke=BLUE_LINE, stroke_width=1.2)
    svg.text(719, 134, "Inputs", size=22, weight=700, anchor="middle", fill=BLUE)

    titles = [
        (216, ["scRNA–seq", "expression"]),
        (526, ["Cell labels", "(cell types)"]),
        (835, ["Cell coordinates", "(2D/3D)"]),
        (1178, ["Spot grid", "(Visium–like)"]),
    ]
    for tx, lines in titles:
        svg.multiline(tx, 181, lines, size=19, weight=700, anchor="middle", fill=BLUE, line_height=1.16)

    cluster_reference(svg, [(167, 274), (260, 268), (226, 327), (153, 375), (213, 387), (291, 366)], scale=1.0, seed=71, points=44)
    svg.text(383, 324, "+", size=42, weight=500, anchor="middle")
    for i, (label, color) in enumerate(zip(["Type 1", "Type 2", "Type 3", "Type 4", "Type 5", "…"], TYPE_COLORS)):
        yy = 246 + i * 33
        svg.el("circle", cx=491, cy=yy, r=8.5, fill=color)
        svg.text(513, yy + 6, label, size=17, weight=500)
    svg.text(671, 324, "+", size=42, weight=500, anchor="middle")
    cluster_reference(svg, [(789, 272), (879, 270), (836, 329), (778, 375), (844, 383), (904, 365)], scale=1.0, seed=72, points=44)
    axes(svg, 736, 414, 1.2)
    svg.text(990, 324, "+", size=42, weight=500, anchor="middle")
    spot_grid(svg, 1180, 328, 93, spacing=17, dot=6.4, filled=False, seed=73)
    end_layer(svg)

    begin_layer(svg, "input-output-arrow", "Input to output arrow")
    line_arrow(svg, 704, 442, 704, 486, width=5)
    end_layer(svg)

    begin_layer(svg, "outputs-panel", "Outputs")
    rounded_box(svg, 81, 491, 1276, 330, stroke="#4a9949", radius=10, width=1.3)
    svg.el("rect", x=82, y=492, width=1274, height=43, rx=9, fill=GREEN_LIGHT)
    svg.el("line", x1=81, y1=535, x2=1357, y2=535, stroke="#4a9949", stroke_width=1.2)
    svg.text(719, 523, "Outputs", size=22, weight=700, anchor="middle", fill=GREEN)
    out_titles = [
        (195, ["Simulated sc", "expression"]),
        (440, ["Simulated ST", "expression"]),
        (683, ["Spot coordinates"]),
        (914, ["Cell–spot", "truth mapping"]),
        (1168, ["Spot × type", "truth fractions"]),
    ]
    for tx, lines in out_titles:
        svg.multiline(tx, 565, lines, size=18, weight=700, anchor="middle", fill=GREEN, line_height=1.18)

    heatmap(svg, 117, 609, rows=7, cols=7, cell=22, seed=81, mode="residual")
    svg.text(194, 783, "Genes", size=17, weight=600, anchor="middle")
    svg.text(105, 692, "Cells", size=17, weight=600, anchor="middle", transform="rotate(-90 105 692)")
    heatmap(svg, 350, 609, rows=7, cols=8, cell=22, seed=82, mode="mixed")
    svg.text(438, 783, "Genes", size=17, weight=600, anchor="middle")
    svg.text(338, 692, "Spots", size=17, weight=600, anchor="middle", transform="rotate(-90 338 692)")
    spot_grid(svg, 682, 692, 78, spacing=17, dot=4.4, filled=False, seed=83)
    heatmap(svg, 842, 613, rows=7, cols=8, cell=22, seed=84, mode="diagonal")
    svg.text(930, 790, "Spots", size=17, weight=600, anchor="middle")
    svg.text(830, 696, "Cells", size=17, weight=600, anchor="middle", transform="rotate(-90 830 696)")
    mini_stacked_bars(svg, 1101, 617, width=140, height=145, bars=6)
    svg.text(1171, 790, "Spots", size=17, weight=600, anchor="middle")
    svg.text(1081, 694, "Fraction", size=16, weight=600, anchor="middle", transform="rotate(-90 1081 694)")
    for i, (label, color) in enumerate(zip(["Type 1", "Type 2", "Type 3", "Type 4", "Type 5", "…"], TYPE_COLORS)):
        yy = 625 + i * 29
        svg.el("rect", x=1264, y=yy - 11, width=17, height=17, fill=color)
        svg.text(1289, yy + 3, label, size=14)
    end_layer(svg)

    begin_layer(svg, "generation-notes", "Generation controls")
    rounded_box(svg, 220, 856, 969, 101, stroke="#353a40", fill="#ffffff", radius=14, width=1.5, stroke_dasharray="10 7")
    svg.text(269, 897, "•  Cell sampling", size=17, weight=500)
    svg.text(438, 897, "•  Biologically informed spatial placement", size=17, weight=500)
    svg.text(845, 897, "•  Expression simulation with noise", size=17, weight=500)
    svg.text(269, 934, "•  Mapping generation", size=17, weight=500)
    svg.text(497, 934, "•  Coverage control", size=17, weight=500)
    end_layer(svg)
    svg.finish(ROOT / "B_editable.svg")


def draw_detailed_stage3a(svg: SVG, x, y, w, h) -> None:
    rounded_box(svg, x, y, w, h, stroke="#f19a4d", radius=16, width=1.4)
    badge(svg, x + 28, y + 34, "3A", fill=ORANGE, radius=19, size=18)
    stage_title(svg, x + w * .59, y + 35, ["Stage 3A: SC–only mismatch", "(reference has, ST lacks)"], color=ORANGE, size=16)
    footer_box(svg, x + 24, y + 87, w - 48, 57, ["Reference contains a type that", "ST does not support"], stroke="#e9ad72", fill=ORANGE_LIGHT, color="#222222", size=13)
    svg.text(x + w / 2, y + 184, "Type–level support evaluation", size=14, weight=700, anchor="middle")
    svg.text(x + 70, y + 231, "Marker evidence", size=12, weight=600, anchor="middle")
    svg.text(x + 205, y + 231, "Similarity evidence", size=12, weight=600, anchor="middle")
    heatmap(svg, x + 28, y + 254, rows=5, cols=5, cell=18, seed=101, mode="diagonal")
    svg.text(x + 13, y + 305, "Genes", size=10, weight=600, anchor="middle", transform=f"rotate(-90 {x+13} {y+305})")
    svg.text(x + 73, y + 363, "Types", size=11, weight=600, anchor="middle")
    cluster_reference(svg, [(x + 174, y + 272), (x + 224, y + 271), (x + 176, y + 324), (x + 227, y + 323)], scale=.6, seed=102, points=20)
    for cx, cy, color in [(x + 174, y + 272, PURPLE), (x + 224, y + 271, "#4f9e3b"), (x + 176, y + 324, "#3177d0"), (x + 227, y + 323, ORANGE)]:
        svg.el("ellipse", cx=cx, cy=cy, rx=24, ry=28, fill="none", stroke=color, stroke_width=1, stroke_dasharray="5 4")
    svg.text(x + 151, y + 306, "UMAP", size=9, weight=600, anchor="middle", transform=f"rotate(-90 {x+151} {y+306})")

    svg.text(x + 18, y + 433, "Type support scoring", size=13, weight=700)
    scores = [
        ("Relative support", "#5cac4b"),
        ("Marker identity", "#9ccf62"),
        ("Neighbour similarity", "#f6a21a"),
        ("Replacement pressure", "#ec3268"),
        ("Robust support / FDR", "#d38b00"),
    ]
    for i, (label, color) in enumerate(scores):
        yy = y + 466 + i * 27
        svg.el("circle", cx=x + 18, cy=yy - 4, r=5.5, fill=color)
        svg.text(x + 31, yy, label, size=10.5)
    svg.el("line", x1=x + 139, y1=y + 452, x2=x + 139, y2=y + 594, stroke="#c4c4c4", stroke_width=1)
    decisions = [("Supported", "keep"), ("Grey zone", "protect"), ("Missing–like", "unknown"), ("Unsupported", "filter")]
    for i, (left, right) in enumerate(decisions):
        yy = y + 478 + i * 34
        svg.text(x + 148, yy, left, size=9.2)
        svg.text(x + 218, yy, "→", size=10, anchor="middle")
        svg.text(x + 266, yy, right, size=9.2, anchor="end")
    footer_box(svg, x + 24, y + h - 112, w - 48, 80, ["Output: plugin_type", "+ cleaned reference pool"], stroke="#e9ad72", fill=ORANGE_LIGHT, color="#a52218", size=14)


def draw_detailed_stage3b(svg: SVG, x, y, w, h) -> None:
    rounded_box(svg, x, y, w, h, stroke=PURPLE, radius=16, width=1.35)
    badge(svg, x + 28, y + 34, "3B", fill=PURPLE, radius=19, size=18)
    stage_title(svg, x + w * .59, y + 35, ["Stage 3B: ST–only mismatch", "(ST has, reference lacks)"], color=PURPLE, size=16)
    footer_box(svg, x + 20, y + 87, w - 40, 57, ["ST contains structure that", "the reference cannot explain"], stroke="#a98bbd", fill=PURPLE_LIGHT, color=PURPLE, size=13)
    svg.text(x + 68, y + 185, "Observed ST", size=12, weight=600, anchor="middle")
    svg.multiline(x + 214, y + 183, ["NNLS reference", "reconstruction"], size=12, weight=600, anchor="middle", line_height=1.2)
    tissue_points(svg, x + 20, y + 211, 112, 90, seed=111, mode="mixed", dot=2.0, spacing=4.9)
    tissue_points(svg, x + 165, y + 211, 112, 90, seed=112, mode="purple", faded=True, dot=2.0, spacing=4.9)
    line_arrow(svg, x + 138, y + 257, x + 158, y + 257, color="#4b5563", width=1.6, marker="arrow-grey")
    svg.text(x + w / 2, y + 338, "Observed − reconstructed residuals", size=12.5, weight=700, anchor="middle")
    heatmap(svg, x + 56, y + 360, rows=4, cols=9, cell=18, seed=113, mode="residual")

    svg.text(x + 13, y + 461, "Multi–feature spot anomaly", size=11.5, weight=700)
    left_items = ["Relative reconstruction error", "Cosine deficit", "Positive residual fraction", "Residual concentration"]
    for i, item in enumerate(left_items):
        svg.text(x + 14, y + 492 + i * 27, "•  " + item, size=9.5)
    svg.el("line", x1=x + 153, y1=y + 450, x2=x + 153, y2=y + 661, stroke="#bababa", stroke_width=1)
    svg.text(x + 165, y + 461, "Statistical calibration", size=11.5, weight=700)
    for i, item in enumerate(["Pseudo–ST (supported)", "One–sided Stouffer test", "BH–FDR control"]):
        svg.text(x + 165, y + 492 + i * 27, "•  " + item, size=9.5)
    svg.text(x + 165, y + 583, "Spatial region detection", size=11.5, weight=700)
    for i, item in enumerate(["Contiguity graph", "Permutation test (region)", "Residual–group branch"]):
        svg.text(x + 165, y + 611 + i * 23, "•  " + item, size=9.2)
    footer_box(svg, x + 20, y + h - 112, w - 40, 80, ["Output: unsupported score", "+ withheld / blank mask"], stroke="#a98bbd", fill=PURPLE_LIGHT, color=PURPLE, size=14)


def draw_detailed_stage4(svg: SVG, x, y, w, h) -> None:
    rounded_box(svg, x, y, w, h, stroke="#60646a", radius=16, width=1.3)
    badge(svg, x + 28, y + 34, "4", fill=CHARCOAL, radius=19)
    stage_title(svg, x + w * .60, y + 34, ["Stage 4: CytoSPACE", "mapping backend"], color="#111111", size=15)
    svg.multiline(x + w / 2, y + 116, ["Same mapping backend", "for baseline and SVTuner"], size=13, weight=600, anchor="middle", line_height=1.3)
    rounded_box(svg, x + 16, y + 168, w - 32, 342, stroke="#a7aaad", fill="#f8f8f8", radius=8, width=1)
    svg.text(x + w / 2, y + 203, "Matched inputs", size=13, weight=700, anchor="middle")
    svg.multiline(x + w / 2, y + 248, ["ST expression + coordinates", "(common genes)", "+", "Reference cells"], size=12, weight=500, anchor="middle", line_height=1.55)
    svg.text(x + w * .31, y + 354, "Original", size=11, weight=600, anchor="middle")
    svg.text(x + w * .31, y + 371, "reference", size=11, weight=600, anchor="middle")
    svg.text(x + w * .72, y + 354, "Stage3A", size=11, weight=600, anchor="middle")
    svg.text(x + w * .72, y + 371, "cleaned pool", size=11, weight=600, anchor="middle")
    cluster_reference(svg, [(x + 64, y + 425), (x + 84, y + 446), (x + 52, y + 457), (x + 91, y + 475)], scale=.42, seed=121, points=12)
    cluster_reference(svg, [(x + 154, y + 425), (x + 174, y + 446), (x + 142, y + 457), (x + 181, y + 475)], scale=.42, seed=122, points=12)
    svg.multiline(x + w / 2, y + 565, ["CytoSPACE optimization", "(cell–spot assignment)"], size=13, weight=600, anchor="middle", line_height=1.3)
    cell_spot_icons(svg, x + 32, y + 646, cols=5, scale=.76)


def draw_detailed_stage5(svg: SVG, x, y, w, h) -> None:
    rounded_box(svg, x, y, w, h, stroke=BLUE_LINE, radius=16, width=1.35)
    badge(svg, x + 28, y + 34, "5", radius=19)
    stage_title(svg, x + w * .60, y + 34, ["Mapping outputs and", "abstention–aware results"], size=15)
    svg.text(x + 16, y + 93, "Baseline route", size=13, weight=700, fill=BLUE)
    svg.text(x + 17, y + 127, "•", size=16, fill=BLUE)
    svg.text(x + 37, y + 127, "Original reference", size=11.5)
    line_arrow(svg, x + w / 2, y + 142, x + w / 2, y + 170, width=1.8)
    svg.multiline(x + w / 2, y + 197, ["CytoSPACE assignment", "(forced at every spot)"], size=11, anchor="middle", line_height=1.25)
    tissue_points(svg, x + 22, y + 224, w - 43, 73, seed=131, mode="spatial", dot=2.5, spacing=5.6)
    svg.el("line", x1=x + 12, y1=y + 310, x2=x + w - 12, y2=y + 310, stroke="#76a7d8", stroke_width=1.2, stroke_dasharray="6 5")
    svg.text(x + 16, y + 345, "SVTuner route", size=13, weight=700, fill=BLUE)
    svg.text(x + 17, y + 378, "•", size=16, fill=BLUE)
    svg.text(x + 37, y + 378, "Stage3A cleaned pool", size=11.5)
    line_arrow(svg, x + w / 2, y + 391, x + w / 2, y + 421, width=1.8)
    svg.multiline(x + w / 2, y + 449, ["CytoSPACE assignment", "(respect same backend)"], size=11, anchor="middle", line_height=1.25)
    svg.text(x + w / 2, y + 502, "Stage3B withheld mask", size=11.5, weight=600, anchor="middle")
    for i in range(8):
        cx = x + 24 + i * (w - 48) / 7
        svg.el("circle", cx=cx, cy=y + 533, r=5, fill="none", stroke=PURPLE, stroke_width=1.3)
        svg.el("line", x1=cx - 3.5, y1=y + 536.5, x2=cx + 3.5, y2=y + 529.5, stroke=PURPLE, stroke_width=1.1)
    line_arrow(svg, x + w / 2, y + 546, x + w / 2, y + 579, width=1.8)
    svg.text(x + w / 2, y + 606, "Withheld–aware output", size=11.5, weight=600, anchor="middle")
    tissue_points(svg, x + 22, y + 625, w - 43, 73, seed=132, mode="spatial", dot=2.5, spacing=5.6)
    footer_box(svg, x + 14, y + h - 99, w - 28, 67, ["Downstream analysis", "(with abstention awareness)"], size=12)


def build_realistic_workflow() -> None:
    svg = SVG(1692, 929, "SVTuner 2.0 realistic workflow – editable reconstruction")
    svg.el("rect", x=0, y=0, width=1692, height=929, fill="#ffffff")
    add_defs(svg)
    begin_layer(svg, "top-brackets", "Top workflow brackets")
    svg.el("line", x1=400, y1=33, x2=1180, y2=33, stroke=BLUE, stroke_width=2)
    for xx in (400, 1180):
        line_arrow(svg, xx, 33, xx, 72, width=2)
    svg.el("line", x1=1182, y1=33, x2=1585, y2=33, stroke=BLUE, stroke_width=2)
    for xx in (1182, 1585):
        line_arrow(svg, xx, 33, xx, 72, width=2)
    svg.el("rect", x=490, y=20, width=497, height=25, fill="#ffffff")
    svg.el("rect", x=1194, y=20, width=390, height=25, fill="#ffffff")
    svg.text(738, 39, "Reference–spatial compatibility diagnosis (pre–mapping)", size=17, weight=700, anchor="middle", fill=BLUE)
    svg.text(1389, 39, "Matched mapping and abstention–aware output", size=17, weight=700, anchor="middle", fill=BLUE)
    end_layer(svg)

    y, h = 84, 782
    begin_layer(svg, "stage-1-input", "1 Input data")
    draw_input_panel(svg, 25, y, 235, h, compact=True)
    end_layer(svg)
    begin_layer(svg, "stage-2-harmonization", "2 Stage 1 harmonization")
    draw_stage1_panel(svg, 285, y, 209, h, compact=True)
    end_layer(svg)
    begin_layer(svg, "stage-3a", "3A SC-only mismatch")
    draw_detailed_stage3a(svg, 518, y, 284, h)
    end_layer(svg)
    begin_layer(svg, "stage-3b", "3B ST-only mismatch")
    draw_detailed_stage3b(svg, 820, y, 295, h)
    end_layer(svg)
    begin_layer(svg, "stage-4-mapping", "4 CytoSPACE mapping backend")
    draw_detailed_stage4(svg, 1138, y, 234, h)
    end_layer(svg)
    begin_layer(svg, "stage-5-results", "5 Mapping outputs")
    draw_detailed_stage5(svg, 1394, y, 274, h)
    end_layer(svg)

    begin_layer(svg, "workflow-arrows", "Workflow arrows")
    for x1, x2, color, marker in [
        (260, 282, BLUE, "arrow-blue"),
        (494, 516, BLUE, "arrow-blue"),
        (802, 818, PURPLE, "arrow-purple"),
        (1115, 1136, BLUE, "arrow-blue"),
        (1372, 1392, BLUE, "arrow-blue"),
    ]:
        line_arrow(svg, x1, 430, x2, 430, color=color, width=2.9, marker=marker)
    end_layer(svg)
    svg.finish(ROOT / "SVTuner_2_0_workflow_realistic_style_editable.svg")


if __name__ == "__main__":
    build_a_simple()
    build_b_simulation()
    build_realistic_workflow()
    for name in [
        "A简洁版_editable.svg",
        "B_editable.svg",
        "SVTuner_2_0_workflow_realistic_style_editable.svg",
    ]:
        print(ROOT / name)
