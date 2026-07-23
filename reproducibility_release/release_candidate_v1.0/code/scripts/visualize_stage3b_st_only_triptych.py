#!/usr/bin/env python
"""Three-panel visualization for an ST-only unsupported-type simulation."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

_HERE = Path(__file__).resolve()
_ROOT = _HERE.parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src.utils.sample_paths import infer_sim_group, resolve_sample_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot truth, CytoSPACE, and Stage3B blank-preserving mapping."
    )
    parser.add_argument("--project_root", default=".")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--target_type", required=True)
    parser.add_argument("--point_size", type=float, default=24.0)
    parser.add_argument("--alpha", type=float, default=0.93)
    parser.add_argument("--out_png", default=None)
    return parser.parse_args()


def _read_coordinates(raw_dir: Path) -> pd.DataFrame:
    path = raw_dir / "brca_STdata_coordinates.txt"
    coords = pd.read_csv(path, sep="\t")
    if coords.shape[1] < 3:
        raise ValueError(f"Invalid coordinate file: {path}")
    coords = coords.iloc[:, :3].copy()
    coords.columns = ["spot_id", "row", "col"]
    coords["spot_id"] = coords["spot_id"].astype(str)
    return coords


def _dominant_truth(path: Path) -> tuple[pd.DataFrame, list[str]]:
    truth = pd.read_csv(path)
    truth = truth.rename(columns={truth.columns[0]: "spot_id"})
    types = [column for column in truth.columns if column != "spot_id"]
    truth[types] = truth[types].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    truth["dominant"] = truth[types].idxmax(axis=1)
    truth.loc[truth[types].sum(axis=1) <= 0, "dominant"] = "__NoType__"
    truth["spot_id"] = truth["spot_id"].astype(str)
    return truth[["spot_id", "dominant"]], types


def _dominant_mapping(path: Path) -> tuple[pd.DataFrame, list[str]]:
    mapping = pd.read_csv(path)
    mapping = mapping.rename(columns={mapping.columns[0]: "spot_id"})
    ignore = {"spot_id", "Total cells", "Unknown_sc_only"}
    types = [column for column in mapping.columns if column not in ignore]
    mapping[types] = mapping[types].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    mapping["dominant"] = mapping[types].idxmax(axis=1)
    mapping.loc[mapping[types].sum(axis=1) <= 0, "dominant"] = "__NoType__"
    mapping["spot_id"] = mapping["spot_id"].astype(str)
    return mapping[["spot_id", "dominant"]], types


def _palette(type_order: list[str]) -> dict[str, str]:
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
    palette = {cell_type: base[i % len(base)] for i, cell_type in enumerate(type_order)}
    if "B cells" in palette and "Epithelial cells" in palette:
        palette["B cells"], palette["Epithelial cells"] = (
            palette["Epithelial cells"],
            palette["B cells"],
        )
    palette["__NoType__"] = "#e5e5e5"
    palette["__BlankUnsupported__"] = "#bdbdbd"
    return palette


def _dataset_label(group: str) -> str:
    labels = {
        "real_brca": "BRCA",
        "human_lung_5loc": "Human Lung",
        "mouse_brain_refined": "Mouse Brain Refined",
    }
    return labels.get(group, group.replace("_", " ").title())


def _plot_panel(
    axis,
    frame: pd.DataFrame,
    title: str,
    palette: dict[str, str],
    point_size: float,
    alpha: float,
) -> None:
    axis.set_facecolor("#e6e6e6")
    present = set(frame["dominant"].astype(str))
    order = ["__NoType__"] + [
        cell_type
        for cell_type in palette
        if cell_type in present and cell_type != "__NoType__"
    ]
    for cell_type in order:
        subset = frame.loc[frame["dominant"] == cell_type]
        if subset.empty:
            continue
        if cell_type == "__BlankUnsupported__":
            axis.scatter(
                subset["col"],
                subset["row"],
                s=point_size,
                facecolors=palette["__BlankUnsupported__"],
                edgecolors="#969696",
                linewidths=0.25,
                alpha=1.0,
                zorder=3,
                rasterized=True,
            )
        else:
            axis.scatter(
                subset["col"],
                subset["row"],
                s=point_size,
                c=palette.get(cell_type, "#333333"),
                edgecolors="none",
                linewidths=0,
                alpha=0.55 if cell_type == "__NoType__" else alpha,
                zorder=1 if cell_type == "__NoType__" else 2,
                rasterized=True,
            )
    axis.set_title(title, fontsize=12, weight="bold")
    axis.set_xlabel("x")
    axis.set_ylabel("y")
    axis.invert_yaxis()
    axis.set_aspect("equal", adjustable="box")
    axis.grid(False)


def main() -> int:
    args = parse_args()
    root = Path(args.project_root).resolve()
    sample = args.sample
    raw_dir = resolve_sample_dir(root, sample)
    group = infer_sim_group(root, sample) or "ungrouped"
    output = (
        Path(args.out_png).resolve()
        if args.out_png
        else root
        / "visualizations"
        / "simulations"
        / group
        / sample
        / "st_only_stage3b_triptych.png"
    )
    output.parent.mkdir(parents=True, exist_ok=True)

    coords = _read_coordinates(raw_dir)
    truth, truth_types = _dominant_truth(
        raw_dir / "sim_truth_spot_type_fraction.csv"
    )
    baseline_path = (
        root
        / "result"
        / sample
        / "stage4_cytospace_baseline"
        / "cytospace_output"
        / "cell_type_assignments_by_spot.csv"
    )
    baseline, baseline_types = _dominant_mapping(baseline_path)
    stage3b_output_dir = (
        root
        / "result"
        / sample
        / "stage4_cytospace_stage3b_blank"
        / "cytospace_output"
    )
    stage3b_mapping_path = (
        stage3b_output_dir / "cell_type_assignments_by_spot.csv"
    )
    stage3b_mapping, stage3b_types = _dominant_mapping(stage3b_mapping_path)
    blank_manifest_path = stage3b_output_dir / "stage3b_blank_spots.csv"
    blank_manifest = pd.read_csv(blank_manifest_path, index_col=0)
    blank_manifest.index = blank_manifest.index.astype(str)
    blank_ids = set(blank_manifest.index)

    stage3b_raw = pd.read_csv(stage3b_mapping_path, index_col=0)
    stage3b_raw.index = stage3b_raw.index.astype(str)
    stage3b_numeric = stage3b_raw.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    missing_blank_rows = blank_ids.difference(stage3b_numeric.index)
    if missing_blank_rows:
        raise ValueError(
            f"{len(missing_blank_rows)} Stage3B blank spots are absent from mapping output"
        )
    blank_nonzero = (
        stage3b_numeric.loc[list(blank_ids)].abs().sum(axis=1) > 0
    )
    if blank_nonzero.any():
        raise ValueError(
            f"{int(blank_nonzero.sum())} Stage3B blank spots have nonzero mapping output"
        )

    type_order = list(truth_types)
    for cell_type in baseline_types + stage3b_types:
        if cell_type not in type_order:
            type_order.append(cell_type)
    palette = _palette(type_order)

    truth_panel = coords.merge(truth, on="spot_id", how="left")
    baseline_panel = coords.merge(baseline, on="spot_id", how="left")
    baseline_panel["dominant"] = baseline_panel["dominant"].fillna("__NoType__")
    svtuner_panel = coords.merge(stage3b_mapping, on="spot_id", how="left")
    svtuner_panel["dominant"] = svtuner_panel["dominant"].fillna("__NoType__")
    svtuner_panel.loc[
        svtuner_panel["spot_id"].isin(blank_ids),
        "dominant",
    ] = "__BlankUnsupported__"

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    figure, axes = plt.subplots(1, 3, figsize=(24, 8), dpi=180)
    figure.patch.set_facecolor("#e6e6e6")
    _plot_panel(
        axes[0],
        truth_panel,
        f"Simulation Truth ({args.target_type} present in ST)",
        palette,
        args.point_size,
        args.alpha,
    )
    _plot_panel(
        axes[1],
        baseline_panel,
        f"CytoSPACE Mapping ({args.target_type} absent from SC)",
        palette,
        args.point_size,
        args.alpha,
    )
    _plot_panel(
        axes[2],
        svtuner_panel,
        "SVTuner + CytoSPACE (unsupported region left blank)",
        palette,
        args.point_size,
        args.alpha,
    )

    legend_types = [cell_type for cell_type in type_order if cell_type != "__NoType__"]
    legend_types.extend(["__BlankUnsupported__", "__NoType__"])
    labels = {cell_type: cell_type for cell_type in legend_types}
    labels["__BlankUnsupported__"] = "Stage3B unsupported / blank"
    labels["__NoType__"] = "Unassigned spot"
    handles = []
    for cell_type in legend_types:
        if cell_type == "__BlankUnsupported__":
            handles.append(
                plt.Line2D(
                    [],
                    [],
                    marker="o",
                    linestyle="",
                    markersize=7,
                    markerfacecolor=palette["__BlankUnsupported__"],
                    markeredgecolor="#888888",
                    markeredgewidth=0.7,
                    label=labels[cell_type],
                )
            )
        else:
            handles.append(
                plt.Line2D(
                    [],
                    [],
                    marker="o",
                    linestyle="",
                    markersize=7,
                    markerfacecolor=palette[cell_type],
                    markeredgecolor="none",
                    label=labels[cell_type],
                )
            )
    figure.legend(
        handles=handles,
        loc="center right",
        bbox_to_anchor=(0.995, 0.5),
        frameon=False,
        fontsize=9,
    )
    figure.suptitle(
        f"{_dataset_label(group)} ST-only unsupported-type simulation",
        fontsize=15,
        weight="bold",
        y=0.98,
    )
    figure.subplots_adjust(left=0.055, right=0.84, bottom=0.1, top=0.86, wspace=0.22)
    figure.savefig(output, dpi=300, bbox_inches="tight", facecolor=figure.get_facecolor())
    figure.savefig(output.with_suffix(".pdf"), bbox_inches="tight", facecolor=figure.get_facecolor())
    plt.close(figure)

    summary = {
        "sample": sample,
        "target_type": args.target_type,
        "truth_spots": int(len(truth_panel)),
        "baseline_spots": int(len(baseline_panel)),
        "stage3b_blank_spots": int(len(blank_ids)),
        "stage3b_mapping_source": str(stage3b_mapping_path),
        "stage3b_blank_manifest": str(blank_manifest_path),
        "stage3b_blank_nonzero_rows": int(blank_nonzero.sum()),
        "png": str(output),
        "pdf": str(output.with_suffix(".pdf")),
    }
    output.with_suffix(".json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
