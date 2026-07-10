#!/usr/bin/env python
"""Export the legacy simulation overview as a fully editable, geometry-matched SVG.

The PNG assembly chain tightly cropped every source panel, resized crops to a
dataset-specific common size, and vertically compressed the mouse-brain stack.
This exporter reproduces those transforms from the underlying spot data instead
of embedding any raster image.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from xml.sax.saxutils import escape, quoteattr

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image


CANVAS_WIDTH = 4044.0
CANVAS_HEIGHT = 5354.0
BASE_RADIUS = 6.12  # 24 pt^2 scatter marker rendered at the legacy 180 dpi.
BASE_COLORS = (
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
)
COLUMN_TITLES = (
    "Truth (no missing type)",
    "CytoSPACE Baseline Mapping",
    "SVTuner + CytoSPACE Mapping (Route2)",
)


@dataclass(frozen=True)
class DatasetSpec:
    group: str
    outer_label: str
    section_title: str
    samples: tuple[str, str, str]
    row_labels: tuple[str, str, str]
    panel_width: int
    panel_height: int
    source_crop_sizes: tuple[
        tuple[tuple[int, int], tuple[int, int], tuple[int, int]],
        tuple[tuple[int, int], tuple[int, int], tuple[int, int]],
        tuple[tuple[int, int], tuple[int, int], tuple[int, int]],
    ]
    stack_width: int
    stack_height: int
    stack_left_width: int
    legend_width: int
    legend_height: int
    row_font_size: int
    global_x: float
    global_y: float
    global_y_scale: float


DATASETS = (
    DatasetSpec(
        group="real_brca",
        outer_label="Real BRCA",
        section_title="Real BRCA Mapping Comparison",
        samples=(
            "real_brca7_candidate_stable_control",
            "real_brca7_candidate_stable_control_missing_epithelial_cells",
            "real_brca7_candidate_stable_control_missing_epithelial_cells_pcs",
        ),
        row_labels=("No missing", "- Epithelial cells", "- Epithelial cells\n- PCs"),
        panel_width=967,
        panel_height=562,
        source_crop_sizes=(
            ((967, 562), (967, 562), (967, 562)),
            ((967, 562), (967, 562), (967, 562)),
            ((967, 562), (967, 562), (967, 562)),
        ),
        stack_width=3654,
        stack_height=1890,
        stack_left_width=430,
        legend_width=191,
        legend_height=254,
        row_font_size=28,
        global_x=360,
        global_y=96,
        global_y_scale=1.0,
    ),
    DatasetSpec(
        group="mouse_brain_refined",
        outer_label="Mouse brain refined",
        section_title="Mouse Brain Refined Mapping Comparison",
        samples=(
            "mouse_brain_refined7_balanced_clustered_sim",
            "mouse_brain_refined7_balanced_clustered_sim_missing_micro_fill_ext_l56",
            "mouse_brain_refined7_balanced_clustered_sim_missing_micro_oligo_2_fill_ext_l56",
        ),
        row_labels=("No missing", "- Microglia", "- Microglia\n- Oligo_2"),
        panel_width=989,
        panel_height=1050,
        source_crop_sizes=(
            ((988, 1050), (988, 1050), (988, 1050)),
            ((988, 1050), (988, 1050), (988, 1050)),
            ((988, 1050), (989, 1050), (988, 1050)),
        ),
        stack_width=3640,
        stack_height=3354,
        stack_left_width=350,
        legend_width=191,
        legend_height=254,
        row_font_size=26,
        global_x=367,
        global_y=2014,
        global_y_scale=0.5,
    ),
    DatasetSpec(
        group="human_lung_5loc",
        outer_label="Human lung 5loc",
        section_title="Human Lung 5-Location Mapping Comparison",
        samples=(
            "human_lung_5loc_fine9_clustered_sim",
            "human_lung_5loc_fine9_clustered_sim_missing_at2",
            "human_lung_5loc_fine9_clustered_sim_missing_at2_fibroblast",
        ),
        row_labels=("No missing", "- AT2", "- AT2\n- Fibroblast"),
        panel_width=994,
        panel_height=469,
        source_crop_sizes=(
            ((993, 442), (993, 469), (994, 469)),
            ((993, 442), (993, 469), (994, 469)),
            ((993, 442), (993, 469), (994, 442)),
        ),
        stack_width=3666,
        stack_height=1611,
        stack_left_width=330,
        legend_width=222,
        legend_height=302,
        row_font_size=28,
        global_x=354,
        global_y=3719,
        global_y_scale=1.0,
    ),
)


def _safe_id(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def _read_coordinates(raw_dir: Path) -> pd.DataFrame:
    frame = pd.read_csv(raw_dir / "brca_STdata_coordinates.txt", sep="\t", usecols=[0, 1, 2])
    frame.columns = ["spot_id", "row", "col"]
    frame["spot_id"] = frame["spot_id"].astype(str)
    return frame


def _find_segments(signal: np.ndarray, threshold: float, min_len: int) -> list[tuple[int, int]]:
    mask = signal > threshold
    segments: list[tuple[int, int]] = []
    start: int | None = None
    for index, value in enumerate(mask):
        if value and start is None:
            start = index
        elif not value and start is not None:
            if index - start >= min_len:
                segments.append((start, index - 1))
            start = None
    if start is not None and len(mask) - start >= min_len:
        segments.append((start, len(mask) - 1))
    return segments


def _source_crop_bounds(path: Path) -> list[tuple[int, int, int, int]]:
    pixels = np.array(Image.open(path).convert("RGB"))
    background = np.array([230, 230, 230], dtype=np.int16)
    difference = np.max(
        np.abs(pixels.astype(np.int16) - background[None, None, :]), axis=2
    )
    foreground = difference > 8
    x_segments = _find_segments(foreground.sum(axis=0), threshold=60, min_len=800)
    if len(x_segments) < 3:
        raise RuntimeError(f"Cannot locate three source panels in {path}")

    bounds: list[tuple[int, int, int, int]] = []
    for x0, x1 in x_segments[:3]:
        y_counts = foreground[:, x0 : x1 + 1].sum(axis=1)
        y_segments = _find_segments(y_counts, threshold=80, min_len=400)
        if y_segments:
            y0, y1 = y_segments[0]
        else:
            smoothed = np.convolve(
                y_counts.astype(np.float32),
                np.ones(7, dtype=np.float32) / 7.0,
                mode="same",
            )
            fallback = _find_segments(smoothed, threshold=60, min_len=220)
            if not fallback:
                fallback = _find_segments(smoothed, threshold=35, min_len=180)
            if not fallback:
                raise RuntimeError(f"Cannot locate source panel y-range in {path}")
            y0, y1 = max(fallback, key=lambda item: item[1] - item[0])
        bounds.append((x0, x1, y0, y1))
    return bounds


def _source_pixel_coordinates(coordinates: pd.DataFrame) -> list[pd.DataFrame]:
    figure, axes = plt.subplots(1, 3, figsize=(24, 8), dpi=180)
    for axis in axes:
        axis.scatter(
            coordinates["col"],
            coordinates["row"],
            s=24.0,
            edgecolors="none",
            linewidths=0.0,
        )
        axis.invert_yaxis()
        axis.set_aspect("equal", adjustable="box")
    figure.subplots_adjust(left=0.04, right=0.86, bottom=0.08, top=0.90, wspace=0.14)
    figure.canvas.draw()
    figure_height = float(figure.canvas.get_width_height()[1])

    projected: list[pd.DataFrame] = []
    values = coordinates[["col", "row"]].to_numpy(dtype=float)
    for axis in axes:
        display = axis.transData.transform(values)
        frame = coordinates[["spot_id"]].copy()
        frame["source_x"] = display[:, 0]
        frame["source_y"] = figure_height - display[:, 1]
        projected.append(frame)
    plt.close(figure)
    return projected


def _dominant_truth(path: Path) -> tuple[pd.DataFrame, list[str]]:
    frame = pd.read_csv(path)
    frame = frame.rename(columns={frame.columns[0]: "spot_id"})
    type_columns = [column for column in frame.columns if column != "spot_id"]
    values = frame[type_columns].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    frame["dominant"] = values.idxmax(axis=1)
    frame.loc[values.sum(axis=1) <= 0.0, "dominant"] = "__NoType__"
    frame["spot_id"] = frame["spot_id"].astype(str)
    return frame[["spot_id", "dominant"]], type_columns


def _dominant_mapping(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    frame = frame.rename(columns={frame.columns[0]: "spot_id"})
    ignored = {"spot_id", "Total cells", "Unknown_sc_only"}
    type_columns = [column for column in frame.columns if column not in ignored]
    values = frame[type_columns].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    frame["dominant"] = values.idxmax(axis=1)
    frame.loc[values.sum(axis=1) <= 0.0, "dominant"] = "__NoType__"
    frame["spot_id"] = frame["spot_id"].astype(str)
    return frame[["spot_id", "dominant"]]


def _load_panels(root: Path, spec: DatasetSpec, sample: str) -> tuple[list[pd.DataFrame], list[str]]:
    raw_dir = root / "data" / "sim" / spec.group / sample
    result_dir = root / "result" / sample
    coordinates = _read_coordinates(raw_dir)
    truth, type_order = _dominant_truth(raw_dir / "sim_truth_spot_type_fraction.csv")
    baseline = _dominant_mapping(
        result_dir
        / "stage4_cytospace_baseline"
        / "cytospace_output"
        / "cell_type_assignments_by_spot.csv"
    )
    route2 = _dominant_mapping(
        result_dir
        / "stage4_cytospace_route2"
        / "cytospace_output"
        / "cell_type_assignments_by_spot.csv"
    )
    projected = _source_pixel_coordinates(coordinates)
    panels = []
    for labels, positions in zip((truth, baseline, route2), projected):
        panel = positions.merge(labels, on="spot_id", how="left")
        panel["dominant"] = panel["dominant"].fillna("__NoType__")
        panels.append(panel)
    return panels, type_order


def _palette(type_order: list[str]) -> dict[str, str]:
    palette = {
        cell_type: BASE_COLORS[index % len(BASE_COLORS)]
        for index, cell_type in enumerate(type_order)
    }
    if "B cells" in palette and "Epithelial cells" in palette:
        palette["B cells"], palette["Epithelial cells"] = (
            palette["Epithelial cells"],
            palette["B cells"],
        )
    palette["__NoType__"] = "#e5e5e5"
    return palette


def _text(
    x: float,
    y: float,
    value: str,
    size: float,
    *,
    anchor: str = "start",
    weight: str = "normal",
    line_spacing: float = 1.16,
) -> str:
    lines = value.split("\n")
    attrs = (
        f'x="{x:.3f}" y="{y:.3f}" text-anchor="{anchor}" '
        f'font-size="{size:.3f}" font-weight="{weight}" '
        'font-family="Arial, sans-serif" fill="#1e1e1e"'
    )
    if len(lines) == 1:
        return f"<text {attrs} dominant-baseline=\"middle\">{escape(value)}</text>"
    half = (len(lines) - 1) * size * line_spacing / 2.0
    tspans = []
    for index, line in enumerate(lines):
        dy = 0.0 if index == 0 else size * line_spacing
        first_y = y - half if index == 0 else None
        y_attr = f' y="{first_y:.3f}"' if first_y is not None else ""
        tspans.append(
            f'<tspan x="{x:.3f}"{y_attr} dy="{dy:.3f}">{escape(line)}</tspan>'
        )
    return f"<text {attrs}>{''.join(tspans)}</text>"


def _panel_svg(
    frame: pd.DataFrame,
    type_order: list[str],
    palette: dict[str, str],
    x: float,
    y: float,
    target_width: float,
    target_height: float,
    crop_bounds: tuple[int, int, int, int],
    panel_id: str,
) -> list[str]:
    assigned = frame.loc[frame["dominant"] != "__NoType__"].copy()
    if assigned.empty:
        return []

    x0, x1, y0, y1 = crop_bounds
    source_width = float(x1 - x0 + 1)
    source_height = float(y1 - y0 + 1)
    scale_x = target_width / source_width
    scale_y = target_height / source_height
    radius_x = BASE_RADIUS * scale_x
    radius_y = BASE_RADIUS * scale_y
    clip_id = panel_id + "__clip"

    lines = [
        f'<defs><clipPath id={quoteattr(clip_id)}><rect x="{x:.3f}" y="{y:.3f}" '
        f'width="{target_width:.3f}" height="{target_height:.3f}"/></clipPath></defs>',
        f'<g id={quoteattr(panel_id)} clip-path="url(#{clip_id})">',
    ]
    for cell_type in type_order:
        subset = assigned.loc[assigned["dominant"] == cell_type]
        if subset.empty:
            continue
        color = palette[cell_type]
        path_commands: list[str] = []
        for row in subset.itertuples(index=False):
            cx = x + (float(row.source_x) - x0) * scale_x
            cy = y + (float(row.source_y) - y0) * scale_y
            if (
                cx + radius_x < x
                or cx - radius_x > x + target_width
                or cy + radius_y < y
                or cy - radius_y > y + target_height
            ):
                continue
            path_commands.append(
                f"M {cx - radius_x:.2f} {cy:.2f} "
                f"a {radius_x:.2f} {radius_y:.2f} 0 1 0 {2.0 * radius_x:.2f} 0 "
                f"a {radius_x:.2f} {radius_y:.2f} 0 1 0 {-2.0 * radius_x:.2f} 0 z"
            )
        if path_commands:
            lines.append(
                f'<path id={quoteattr(panel_id + "__" + _safe_id(cell_type))} '
                f'd={quoteattr(" ".join(path_commands))} fill="{color}" '
                'fill-opacity="0.93" stroke="none"/>'
            )
    lines.append("</g>")
    return lines


def _legend_svg(
    spec: DatasetSpec,
    type_order: list[str],
    palette: dict[str, str],
    x: float,
    y: float,
) -> list[str]:
    labels = type_order + ["Unassigned spots"]
    colors = [palette[cell_type] for cell_type in type_order] + [palette["__NoType__"]]
    lines = [f'<g id={quoteattr(spec.group + "__legend")}>']
    lines.append(_text(x + 14, y + 27, "Type", 28, weight="bold"))
    item_y = y + 54
    for label, color in zip(labels, colors):
        cy = item_y + 12
        lines.append(
            f'<circle cx="{x + 19:.3f}" cy="{cy:.3f}" r="5" fill="{color}" '
            'stroke="#787878" stroke-width="1"/>'
        )
        lines.append(_text(x + 34, cy, label, 18))
        item_y += 24
    lines.append("</g>")
    return lines


def _section_svg(root: Path, spec: DatasetSpec) -> list[str]:
    lines = [
        f'<g id={quoteattr(spec.group)} transform="translate({spec.global_x:.3f} '
        f'{spec.global_y:.3f}) scale(1 {spec.global_y_scale:.6f})">'
    ]
    x0 = 24 + spec.stack_left_width
    y0 = 24 + 68 + 44
    column_gap = 24
    row_gap = 22

    lines.append(
        _text(
            spec.stack_width / 2.0,
            47,
            spec.section_title,
            36,
            anchor="middle",
            weight="bold",
        )
    )
    for column_index, title in enumerate(COLUMN_TITLES):
        cx = x0 + column_index * (spec.panel_width + column_gap) + spec.panel_width / 2.0
        lines.append(_text(cx, 103, title, 22, anchor="middle", weight="bold"))

    control_truth = (
        root
        / "data"
        / "sim"
        / spec.group
        / spec.samples[0]
        / "sim_truth_spot_type_fraction.csv"
    )
    _, type_order = _dominant_truth(control_truth)
    palette = _palette(type_order)

    for row_index, (sample, row_label) in enumerate(zip(spec.samples, spec.row_labels)):
        panels, sample_types = _load_panels(root, spec, sample)
        if sample_types != type_order:
            raise ValueError(f"Type order changed within {spec.group}: {sample}")
        panel_y = y0 + row_index * (spec.panel_height + row_gap)
        lines.append(
            _text(
                34,
                panel_y + spec.panel_height / 2.0,
                row_label,
                spec.row_font_size,
                weight="bold",
            )
        )
        source_png_name = "mapping_triptych_no_missing.png" if row_index == 0 else "missing_type_triptych.png"
        source_png = (
            root
            / "visualizations"
            / "simulations"
            / spec.group
            / sample
            / source_png_name
        )
        crop_bounds = _source_crop_bounds(source_png)
        for column_index, panel in enumerate(panels):
            panel_x = x0 + column_index * (spec.panel_width + column_gap)
            source_width, source_height = spec.source_crop_sizes[row_index][column_index]
            crop_x0, crop_x1, crop_y0, crop_y1 = crop_bounds[column_index]
            actual_size = (crop_x1 - crop_x0 + 1, crop_y1 - crop_y0 + 1)
            if actual_size != (source_width, source_height):
                raise ValueError(
                    f"Unexpected source crop for {sample} panel {column_index + 1}: "
                    f"{actual_size} != {(source_width, source_height)}"
                )
            panel_id = f"{spec.group}__row_{row_index + 1}__panel_{column_index + 1}"
            lines.extend(
                _panel_svg(
                    panel,
                    type_order,
                    palette,
                    panel_x,
                    panel_y,
                    spec.panel_width,
                    spec.panel_height,
                    crop_bounds[column_index],
                    panel_id,
                )
            )

    legend_x = x0 + 3 * spec.panel_width + 2 * column_gap + 18
    legend_y = y0 + 8
    lines.extend(_legend_svg(spec, type_order, palette, legend_x, legend_y))
    lines.append("</g>")
    return lines


def _build(root: Path, output: Path) -> None:
    lines = [
        '<?xml version="1.0" encoding="UTF-8" standalone="no"?>',
        (
            '<svg xmlns="http://www.w3.org/2000/svg" '
            'xmlns:xlink="http://www.w3.org/1999/xlink" '
            'width="2022" height="2677" viewBox="0 0 4044 5354">'
        ),
        '<title>Simulation Mapping Visualization Overview</title>',
        '<desc>Fully editable vector reconstruction from spot coordinates and mapping assignments.</desc>',
        '<rect id="background" x="0" y="0" width="4044" height="5354" fill="#ffffff"/>',
        _text(
            CANVAS_WIDTH / 2.0,
            51,
            "Simulation Mapping Visualization Overview",
            38,
            anchor="middle",
            weight="bold",
        ),
    ]

    for spec in DATASETS:
        displayed_height = spec.stack_height * spec.global_y_scale
        lines.append(
            _text(
                32,
                spec.global_y + displayed_height / 2.0,
                spec.outer_label,
                28,
                weight="bold",
            )
        )
        lines.extend(_section_svg(root, spec))

    lines.append("</svg>")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    output = (
        root
        / "visualizations"
        / "simulations"
        / "simulation_triptych_overview_stack_3datasets.svg"
    )
    _build(root, output)
    print(f"[OK] wrote geometry-matched editable SVG: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
