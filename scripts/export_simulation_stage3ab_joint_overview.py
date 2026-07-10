#!/usr/bin/env python
"""Build the ordered Stage3A+Stage3B simulation overview as SVG and PNG.

The SVG keeps text editable and combines spots of the same type within each
panel into compound paths so Adobe applications do not need to manage more
than one hundred thousand independent objects.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import re
import subprocess
import tempfile
from xml.sax.saxutils import escape, quoteattr

import pandas as pd


CANVAS_WIDTH = 4044
CANVAS_HEIGHT = 5960
PANEL_X = 420
PANEL_WIDTH = 930
PANEL_GAP = 24
SECTION_HEADER_HEIGHT = 105
ROW_GAP = 24
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
BLANK_COLOR = "#bdbdbd"
BLANK_EDGE = "#969696"
COLUMN_TITLES = (
    "Simulation Truth",
    "CytoSPACE Baseline (SC target absent)",
    "SVTuner + CytoSPACE (Stage3B blank-preserving)",
)


@dataclass(frozen=True)
class Scenario:
    sample: str
    row_label: str


@dataclass(frozen=True)
class DatasetSpec:
    group: str
    outer_label: str
    section_title: str
    target_type: str
    scenarios: tuple[Scenario, Scenario, Scenario]
    section_y: int
    panel_height: int
    spot_radius: float
    visible_row_min: float | None = None
    visible_row_max: float | None = None


DATASETS = (
    DatasetSpec(
        group="real_brca",
        outer_label="Real BRCA",
        section_title="Real BRCA Joint Stage3A + Stage3B Comparison",
        target_type="Endothelial cells",
        scenarios=(
            Scenario(
                "real_brca7_endothelial_marker_control_sc_missing_endothelial_cells",
                "Control ST\n- Endothelial cells (SC)",
            ),
            Scenario(
                "real_brca7_endothelial_marker_missing_epithelial_cells_sc_missing_endothelial_cells",
                "- Epithelial cells (ST)\n- Endothelial cells (SC)",
            ),
            Scenario(
                "real_brca7_endothelial_marker_missing_epithelial_cells_pcs_sc_missing_endothelial_cells",
                "- Epithelial cells, PCs (ST)\n- Endothelial cells (SC)",
            ),
        ),
        section_y=90,
        panel_height=520,
        spot_radius=5.90,
    ),
    DatasetSpec(
        group="mouse_brain_refined",
        outer_label="Mouse brain refined",
        section_title="Mouse Brain Refined Joint Stage3A + Stage3B Comparison",
        target_type="Ext_L56",
        scenarios=(
            Scenario(
                "mouse_brain_refined7_balanced_clustered_sim_sc_missing_ext_l56",
                "Control ST\n- Ext_L56 (SC)",
            ),
            Scenario(
                "mouse_brain_refined7_balanced_clustered_sim_missing_micro_fill_inh_pvalb_sc_missing_ext_l56",
                "- Micro (ST)\n- Ext_L56 (SC)",
            ),
            Scenario(
                "mouse_brain_refined7_balanced_clustered_sim_missing_micro_oligo_2_fill_inh_pvalb_sc_missing_ext_l56",
                "- Micro, Oligo_2 (ST)\n- Ext_L56 (SC)",
            ),
        ),
        section_y=1855,
        panel_height=760,
        spot_radius=3.45,
    ),
    DatasetSpec(
        group="human_lung_5loc",
        outer_label="Human lung 5loc",
        section_title="Human Lung 5-Location Joint Stage3A + Stage3B Comparison",
        target_type="B_cell",
        scenarios=(
            Scenario(
                "human_lung_5loc_fine9_clustered_sim_sc_missing_b_cell",
                "Control ST\n- B_cell (SC)",
            ),
            Scenario(
                "human_lung_5loc_fine9_clustered_sim_missing_at2_sc_missing_b_cell",
                "- AT2 (ST)\n- B_cell (SC)",
            ),
            Scenario(
                "human_lung_5loc_fine9_clustered_sim_missing_at2_fibroblast_sc_missing_b_cell",
                "- AT2, Fibroblast (ST)\n- B_cell (SC)",
            ),
        ),
        section_y=4340,
        panel_height=470,
        spot_radius=5.70,
        visible_row_min=10874.0,
        visible_row_max=17396.0,
    ),
)


def _safe_id(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def _text(
    x: float,
    y: float,
    value: str,
    size: float,
    *,
    anchor: str = "start",
    weight: str = "normal",
    fill: str = "#1e1e1e",
    line_spacing: float = 1.16,
) -> str:
    lines = value.split("\n")
    attributes = (
        f'x="{x:.2f}" y="{y:.2f}" text-anchor="{anchor}" '
        f'font-size="{size:.2f}" font-weight="{weight}" '
        f'font-family="Arial, sans-serif" fill="{fill}"'
    )
    if len(lines) == 1:
        return f'<text {attributes} dominant-baseline="middle">{escape(value)}</text>'
    half = (len(lines) - 1) * size * line_spacing / 2.0
    tspans = []
    for index, line in enumerate(lines):
        if index == 0:
            tspans.append(
                f'<tspan x="{x:.2f}" y="{y - half:.2f}">{escape(line)}</tspan>'
            )
        else:
            tspans.append(
                f'<tspan x="{x:.2f}" dy="{size * line_spacing:.2f}">{escape(line)}</tspan>'
            )
    return f"<text {attributes}>{''.join(tspans)}</text>"


def _read_coordinates(raw_dir: Path) -> pd.DataFrame:
    frame = pd.read_csv(raw_dir / "brca_STdata_coordinates.txt", sep="\t", usecols=[0, 1, 2])
    frame.columns = ["spot_id", "row", "col"]
    frame["spot_id"] = frame["spot_id"].astype(str)
    return frame


def _dominant_truth(path: Path) -> tuple[pd.DataFrame, list[str]]:
    frame = pd.read_csv(path)
    frame = frame.rename(columns={frame.columns[0]: "spot_id"})
    type_columns = [column for column in frame.columns if column != "spot_id"]
    values = frame[type_columns].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    frame["dominant"] = values.idxmax(axis=1)
    frame.loc[values.sum(axis=1) <= 0.0, "dominant"] = "__NoType__"
    frame["spot_id"] = frame["spot_id"].astype(str)
    return frame[["spot_id", "dominant"]], type_columns


def _dominant_mapping(path: Path) -> tuple[pd.DataFrame, list[str]]:
    frame = pd.read_csv(path)
    frame = frame.rename(columns={frame.columns[0]: "spot_id"})
    ignored = {"spot_id", "Total cells", "Unknown_sc_only"}
    type_columns = [column for column in frame.columns if column not in ignored]
    values = frame[type_columns].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    frame["dominant"] = values.idxmax(axis=1)
    frame.loc[values.sum(axis=1) <= 0.0, "dominant"] = "__NoType__"
    frame["spot_id"] = frame["spot_id"].astype(str)
    return frame[["spot_id", "dominant"]], type_columns


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
    palette["__BlankUnsupported__"] = BLANK_COLOR
    return palette


def _load_scenario(
    root: Path, spec: DatasetSpec, scenario: Scenario
) -> tuple[list[pd.DataFrame], list[str], int]:
    raw_dir = root / "data" / "sim" / spec.group / scenario.sample
    result_dir = root / "result" / scenario.sample
    coordinates = _read_coordinates(raw_dir)
    truth, truth_types = _dominant_truth(raw_dir / "sim_truth_spot_type_fraction.csv")
    baseline, baseline_types = _dominant_mapping(
        result_dir
        / "stage4_cytospace_baseline"
        / "cytospace_output"
        / "cell_type_assignments_by_spot.csv"
    )
    stage3b_dir = result_dir / "stage4_cytospace_stage3b_blank" / "cytospace_output"
    stage3b_path = stage3b_dir / "cell_type_assignments_by_spot.csv"
    stage3b, stage3b_types = _dominant_mapping(stage3b_path)
    blank_manifest = pd.read_csv(stage3b_dir / "stage3b_blank_spots.csv", index_col=0)
    blank_ids = set(blank_manifest.index.astype(str))

    raw_stage3b = pd.read_csv(stage3b_path, index_col=0)
    raw_stage3b.index = raw_stage3b.index.astype(str)
    numeric_stage3b = raw_stage3b.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    missing_blank_rows = blank_ids.difference(numeric_stage3b.index)
    if missing_blank_rows:
        raise ValueError(
            f"{scenario.sample}: {len(missing_blank_rows)} blank spots missing from Stage3B output"
        )
    nonzero_blank_rows = numeric_stage3b.loc[list(blank_ids)].abs().sum(axis=1) > 0
    if nonzero_blank_rows.any():
        raise ValueError(
            f"{scenario.sample}: {int(nonzero_blank_rows.sum())} blank spots have nonzero mapping"
        )

    type_order = list(truth_types)
    for cell_type in baseline_types + stage3b_types:
        if cell_type not in type_order:
            type_order.append(cell_type)

    panels = []
    for labels in (truth, baseline, stage3b):
        panel = coordinates.merge(labels, on="spot_id", how="left")
        panel["dominant"] = panel["dominant"].fillna("__NoType__")
        panels.append(panel)
    panels[2].loc[
        panels[2]["spot_id"].isin(blank_ids), "dominant"
    ] = "__BlankUnsupported__"
    return panels, type_order, len(blank_ids)


def _circle_compound_path(
    points: pd.DataFrame,
    radius: float,
    map_x,
    map_y,
) -> str:
    commands = []
    for point in points.itertuples(index=False):
        cx = map_x(float(point.col))
        cy = map_y(float(point.row))
        commands.append(
            f"M {cx - radius:.2f} {cy:.2f} "
            f"a {radius:.2f} {radius:.2f} 0 1 0 {2.0 * radius:.2f} 0 "
            f"a {radius:.2f} {radius:.2f} 0 1 0 {-2.0 * radius:.2f} 0 z"
        )
    return " ".join(commands)


def _panel_svg(
    frame: pd.DataFrame,
    type_order: list[str],
    palette: dict[str, str],
    panel_x: float,
    panel_y: float,
    panel_height: float,
    radius: float,
    panel_id: str,
    visible_row_min: float | None = None,
    visible_row_max: float | None = None,
) -> list[str]:
    display_frame = frame
    if visible_row_min is not None:
        display_frame = display_frame.loc[display_frame["row"] >= visible_row_min]
    if visible_row_max is not None:
        display_frame = display_frame.loc[display_frame["row"] <= visible_row_max]
    if display_frame.empty:
        raise ValueError(f"No spots remain in display window for {panel_id}")

    min_col = float(display_frame["col"].min())
    max_col = float(display_frame["col"].max())
    min_row = float(display_frame["row"].min())
    max_row = float(display_frame["row"].max())
    col_range = max(1.0, max_col - min_col)
    row_range = max(1.0, max_row - min_row)
    margin = max(10.0, radius * 2.0)
    scale = min(
        (PANEL_WIDTH - 2.0 * margin) / col_range,
        (panel_height - 2.0 * margin) / row_range,
    )
    tissue_width = col_range * scale
    tissue_height = row_range * scale
    offset_x = panel_x + (PANEL_WIDTH - tissue_width) / 2.0
    offset_y = panel_y + (panel_height - tissue_height) / 2.0
    map_x = lambda value: offset_x + (value - min_col) * scale
    map_y = lambda value: offset_y + (value - min_row) * scale
    clip_id = panel_id + "__clip"

    lines = [
        f'<defs><clipPath id={quoteattr(clip_id)}><rect x="{panel_x:.2f}" '
        f'y="{panel_y:.2f}" width="{PANEL_WIDTH:.2f}" height="{panel_height:.2f}"/>'
        "</clipPath></defs>",
        f'<rect x="{panel_x:.2f}" y="{panel_y:.2f}" width="{PANEL_WIDTH:.2f}" '
        f'height="{panel_height:.2f}" fill="#ffffff" stroke="#d1d1d1" stroke-width="1.2"/>',
        f'<g id={quoteattr(panel_id)} clip-path="url(#{clip_id})">',
    ]
    present = set(display_frame["dominant"].astype(str))
    draw_order = [cell_type for cell_type in type_order if cell_type in present]
    if "__BlankUnsupported__" in present:
        draw_order.append("__BlankUnsupported__")
    for cell_type in draw_order:
        if cell_type == "__NoType__":
            continue
        subset = display_frame.loc[
            display_frame["dominant"] == cell_type, ["col", "row"]
        ]
        if subset.empty:
            continue
        path_data = _circle_compound_path(subset, radius, map_x, map_y)
        stroke = BLANK_EDGE if cell_type == "__BlankUnsupported__" else "none"
        stroke_width = "0.65" if cell_type == "__BlankUnsupported__" else "0"
        opacity = "1" if cell_type == "__BlankUnsupported__" else "0.93"
        lines.append(
            f'<path id={quoteattr(panel_id + "__" + _safe_id(cell_type))} '
            f'd={quoteattr(path_data)} fill="{palette[cell_type]}" fill-opacity="{opacity}" '
            f'stroke="{stroke}" stroke-width="{stroke_width}"/>'
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
    labels = type_order + ["__BlankUnsupported__"]
    display = {cell_type: cell_type for cell_type in labels}
    display["__BlankUnsupported__"] = "Stage3B unsupported / blank"
    lines = [f'<g id={quoteattr(spec.group + "__legend")}>']
    lines.append(_text(x, y, "Type", 25, weight="bold"))
    item_y = y + 38
    for cell_type in labels:
        lines.append(
            f'<circle cx="{x + 6:.2f}" cy="{item_y:.2f}" r="5.2" '
            f'fill="{palette[cell_type]}" '
            f'stroke="{BLANK_EDGE if cell_type == "__BlankUnsupported__" else "none"}" '
            'stroke-width="0.7"/>'
        )
        lines.append(_text(x + 22, item_y, display[cell_type], 16.5))
        item_y += 27
    lines.append("</g>")
    return lines


def _build_svg(root: Path, output: Path) -> dict:
    lines = [
        '<?xml version="1.0" encoding="UTF-8" standalone="no"?>',
        (
            '<svg xmlns="http://www.w3.org/2000/svg" '
            f'width="{CANVAS_WIDTH}" height="{CANVAS_HEIGHT}" '
            f'viewBox="0 0 {CANVAS_WIDTH} {CANVAS_HEIGHT}">'
        ),
        '<title>Stage3A + Stage3B Joint Simulation Mapping Overview</title>',
        '<desc>Ordered joint Stage3A and Stage3B simulation comparison across three datasets.</desc>',
        f'<rect x="0" y="0" width="{CANVAS_WIDTH}" height="{CANVAS_HEIGHT}" fill="#ffffff"/>',
        _text(
            CANVAS_WIDTH / 2,
            43,
            "Stage3A + Stage3B Joint Simulation Mapping Overview",
            38,
            anchor="middle",
            weight="bold",
        ),
    ]
    summary = {"datasets": [], "order_preserved": True}

    for spec in DATASETS:
        section_height = SECTION_HEADER_HEIGHT + 3 * spec.panel_height + 2 * ROW_GAP
        lines.append(
            _text(
                32,
                spec.section_y + section_height / 2.0,
                spec.outer_label,
                27,
                weight="bold",
            )
        )
        lines.append(
            _text(
                PANEL_X + (3 * PANEL_WIDTH + 2 * PANEL_GAP) / 2.0,
                spec.section_y + 27,
                spec.section_title,
                31,
                anchor="middle",
                weight="bold",
            )
        )
        for column_index, title in enumerate(COLUMN_TITLES):
            panel_x = PANEL_X + column_index * (PANEL_WIDTH + PANEL_GAP)
            lines.append(
                _text(
                    panel_x + PANEL_WIDTH / 2.0,
                    spec.section_y + 77,
                    title,
                    20.5,
                    anchor="middle",
                    weight="bold",
                )
            )

        dataset_summary = {
            "group": spec.group,
            "target_type": spec.target_type,
            "scenarios": [],
        }
        control_type_order: list[str] | None = None
        control_palette: dict[str, str] | None = None
        for row_index, scenario in enumerate(spec.scenarios):
            panels, type_order, blank_count = _load_scenario(root, spec, scenario)
            if control_type_order is None:
                control_type_order = type_order
                control_palette = _palette(type_order)
            elif type_order != control_type_order:
                raise ValueError(f"Type order changed within {spec.group}: {scenario.sample}")
            assert control_palette is not None
            panel_y = spec.section_y + SECTION_HEADER_HEIGHT + row_index * (
                spec.panel_height + ROW_GAP
            )
            lines.append(
                _text(
                    138,
                    panel_y + spec.panel_height / 2.0,
                    scenario.row_label,
                    21.5,
                    weight="bold",
                )
            )
            for column_index, panel in enumerate(panels):
                panel_x = PANEL_X + column_index * (PANEL_WIDTH + PANEL_GAP)
                lines.extend(
                    _panel_svg(
                        panel,
                        control_type_order,
                        control_palette,
                        panel_x,
                        panel_y,
                        spec.panel_height,
                        spec.spot_radius,
                        f"{spec.group}__row_{row_index + 1}__panel_{column_index + 1}",
                        spec.visible_row_min,
                        spec.visible_row_max,
                    )
                )
            dataset_summary["scenarios"].append(
                {
                    "sample": scenario.sample,
                    "row_label": scenario.row_label,
                    "stage3b_blank_spots": blank_count,
                }
            )

        assert control_type_order is not None and control_palette is not None
        legend_x = PANEL_X + 3 * PANEL_WIDTH + 2 * PANEL_GAP + 24
        legend_y = spec.section_y + SECTION_HEADER_HEIGHT + 28
        lines.extend(
            _legend_svg(spec, control_type_order, control_palette, legend_x, legend_y)
        )
        summary["datasets"].append(dataset_summary)

    lines.append("</svg>")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return summary


def _edge_executable() -> Path:
    candidates = (
        Path(r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe"),
        Path(r"C:\Program Files\Microsoft\Edge\Application\msedge.exe"),
        Path(r"C:\Program Files\Google\Chrome\Application\chrome.exe"),
        Path(r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe"),
    )
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError("No supported headless Edge/Chrome executable found")


def _render_png(svg: Path, png: Path) -> None:
    browser = _edge_executable()
    with tempfile.TemporaryDirectory(prefix="svtuner_stage3ab_overview_") as profile:
        command = [
            str(browser),
            "--headless=new",
            "--disable-gpu",
            "--hide-scrollbars",
            "--allow-file-access-from-files",
            "--force-device-scale-factor=1",
            f"--user-data-dir={profile}",
            f"--window-size={CANVAS_WIDTH},{CANVAS_HEIGHT}",
            f"--screenshot={png.resolve()}",
            svg.resolve().as_uri(),
        ]
        subprocess.run(
            command,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    if not png.exists():
        raise RuntimeError(f"PNG render did not produce {png}")


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    output_base = (
        root
        / "visualizations"
        / "simulations"
        / "simulation_stage3ab_joint_triptych_overview_stack_3datasets"
    )
    svg = output_base.with_suffix(".svg")
    png = output_base.with_suffix(".png")
    summary_path = output_base.with_suffix(".json")
    summary = _build_svg(root, svg)
    _render_png(svg, png)
    summary.update({"svg": str(svg), "png": str(png)})
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"[OK] wrote {svg}")
    print(f"[OK] wrote {png}")
    print(f"[OK] wrote {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
