#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


SCENARIOS: list[tuple[str, str]] = [
    ("real_brca", "real_brca_clustered_sim"),
    ("real_brca", "real_brca_clustered_sim_missing_epithelial_cells"),
    ("real_brca", "real_brca_clustered_sim_missing_epithelial_monocytes_macrophages"),
    ("real_brca", "real_brca_clustered_sim_missing_epithelial_monocytes_endothelial"),
    ("real_brca", "real_brca_clustered_sim_missing_epithelial_monocytes_endothelial_fibroblasts"),
    ("human_lung_5loc", "human_lung_5loc_fine9_clustered_sim"),
    ("human_lung_5loc", "human_lung_5loc_fine9_clustered_sim_missing_ciliated"),
    ("human_lung_5loc", "human_lung_5loc_fine9_clustered_sim_missing_ciliated_endothelia_vascular"),
    ("mouse_brain_refined", "mouse_brain_refined8_balanced_clustered_sim"),
    ("mouse_brain_refined", "mouse_brain_refined8_balanced_clustered_sim_missing_micro_fill_ext_l56"),
    ("mouse_brain_refined", "mouse_brain_refined8_balanced_clustered_sim_missing_micro_astro_ctx_fill_ext_l56"),
    ("mouse_brain_refined", "mouse_brain_refined8_balanced_clustered_sim_missing_micro_astro_ctx_oligo_2_fill_ext_l56"),
]


METHODS = [
    ("CytoSPACE", "cytospace_baseline"),
    ("SVTuner + CytoSPACE", "cytospace_route2"),
    ("Tangram\n(all genes)", "tangram_all"),
    ("Tangram\n(marker genes)", "tangram_marker"),
    ("novoSpaRc", "novosparc"),
    ("SpaOTsc", "spaotsc"),
    ("CellTrek", "celltrek"),
]


PALETTE = {
    "CytoSPACE": "#E83E78",
    "SVTuner + CytoSPACE": "#2A9D8F",
    "Tangram\n(all genes)": "#3AA1B8",
    "Tangram\n(marker genes)": "#6EA35C",
    "novoSpaRc": "#7A4DA0",
    "SpaOTsc": "#FF7F0E",
    "CellTrek": "#B87333",
}


def _truth_path(project_root: Path, group: str, sample: str) -> Path:
    path = project_root / "data" / "sim" / group / sample / "sim_truth_spot_type_fraction.csv"
    if not path.exists():
        raise FileNotFoundError(f"truth fraction not found: {path}")
    return path


def _prediction_path(project_root: Path, sample: str, method_dir: str) -> Path:
    if method_dir == "cytospace_baseline":
        return project_root / "result" / sample / "stage4_cytospace_baseline" / "cytospace_output" / "fractional_abundances_by_spot.csv"
    if method_dir == "cytospace_route2":
        return project_root / "result" / sample / "stage4_cytospace_route2" / "cytospace_output" / "fractional_abundances_by_spot.csv"
    return project_root / "result" / sample / "stage4_mapping" / method_dir / "spot_type_fraction.csv"


def _load_fraction(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"fraction file not found: {path}")
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str)
    return df.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def _composition_recovery(pred: pd.DataFrame, truth: pd.DataFrame) -> tuple[float, pd.DataFrame]:
    common_spots = pred.index.intersection(truth.index)
    if common_spots.empty:
        raise ValueError("no common spots between prediction and truth.")
    type_cols = sorted(set(pred.columns).union(set(truth.columns)))
    p = pred.reindex(index=common_spots, columns=type_cols, fill_value=0.0).astype(float)
    t = truth.reindex(index=common_spots, columns=type_cols, fill_value=0.0).astype(float)
    p = p.div(p.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    t = t.div(t.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)

    rows: list[dict[str, Any]] = []
    weights = []
    values = []
    for cell_type in type_cols:
        truth_mass = float(t[cell_type].sum())
        if truth_mass <= 0:
            continue
        recovered = float(np.minimum(p[cell_type].to_numpy(), t[cell_type].to_numpy()).sum())
        value = recovered / truth_mass
        rows.append(
            {
                "cell_type": cell_type,
                "truth_mass": truth_mass,
                "recovered_mass": recovered,
                "composition_recovery": value,
            }
        )
        weights.append(truth_mass)
        values.append(value)
    if not values:
        return float("nan"), pd.DataFrame(rows)
    scenario_score = float(np.average(np.asarray(values), weights=np.asarray(weights)))
    return scenario_score, pd.DataFrame(rows)


def build_tables(project_root: Path, sample_suffix: str = "") -> tuple[pd.DataFrame, pd.DataFrame]:
    scenario_rows: list[dict[str, Any]] = []
    cell_type_rows: list[pd.DataFrame] = []
    for group, sample in SCENARIOS:
        eval_sample = f"{sample}{sample_suffix}"
        truth = _load_fraction(_truth_path(project_root, group, eval_sample))
        for method_label, method_dir in METHODS:
            pred = _load_fraction(_prediction_path(project_root, eval_sample, method_dir))
            score, per_type = _composition_recovery(pred, truth)
            scenario_rows.append(
                {
                    "group": group,
                    "sample": eval_sample,
                    "method": method_label,
                    "composition_recovery": score,
                }
            )
            per_type["group"] = group
            per_type["sample"] = eval_sample
            per_type["method"] = method_label
            cell_type_rows.append(per_type)
    scenario_df = pd.DataFrame(scenario_rows)
    per_type_df = pd.concat(cell_type_rows, ignore_index=True)
    order = [m[0] for m in METHODS]
    scenario_df["method"] = pd.Categorical(scenario_df["method"], order, ordered=True)
    per_type_df["method"] = pd.Categorical(per_type_df["method"], order, ordered=True)
    return scenario_df, per_type_df


def plot(df: pd.DataFrame, out_png: Path, out_pdf: Path | None = None, title: str | None = None) -> None:
    sns.set_theme(style="whitegrid")
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.edgecolor": "#333333",
            "axes.linewidth": 1.0,
            "axes.titleweight": "bold",
            "axes.labelweight": "bold",
        }
    )
    order = [m[0] for m in METHODS]
    fig, ax = plt.subplots(figsize=(10.8, 5.8), dpi=220, constrained_layout=True)
    sns.boxplot(
        data=df,
        x="method",
        hue="method",
        y="composition_recovery",
        order=order,
        hue_order=order,
        palette=PALETTE,
        width=0.62,
        linewidth=1.6,
        fliersize=0,
        legend=False,
        ax=ax,
        boxprops={"alpha": 0.95, "edgecolor": "#333333"},
        medianprops={"color": "#333333", "linewidth": 1.8},
        whiskerprops={"color": "#333333", "linewidth": 1.3, "linestyle": "--"},
        capprops={"color": "#333333", "linewidth": 1.3},
    )
    sns.stripplot(
        data=df,
        x="method",
        y="composition_recovery",
        order=order,
        color="#333333",
        size=3.2,
        jitter=0.16,
        alpha=0.55,
        ax=ax,
        zorder=5,
    )
    ax.set_title(title or "Spatial cell-type composition recovery, no sc noise", fontsize=15, pad=12)
    ax.set_xlabel("")
    ax.set_ylabel("Composition recovery score", fontsize=12)
    ax.set_ylim(0, 1.02)
    ax.grid(axis="y", color="#E2E2E2", linestyle="--", linewidth=0.9)
    ax.grid(axis="x", visible=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="x", labelsize=10)
    for tick in ax.get_xticklabels():
        tick.set_rotation(0)
        tick.set_ha("center")

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, bbox_inches="tight")
    if out_pdf is not None:
        fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot scenario-level spatial composition recovery across 7 methods.")
    p.add_argument("--project_root", default=".", help="Project root.")
    p.add_argument(
        "--out_dir",
        default="visualizations/method_comparison/no_noise",
        help="Output directory.",
    )
    p.add_argument(
        "--output_prefix",
        default="composition_recovery_7mapping_methods_no_noise",
        help="Output filename prefix.",
    )
    p.add_argument(
        "--sample_suffix",
        default="",
        help="Suffix appended to each base simulation sample before loading truth/results.",
    )
    p.add_argument(
        "--title",
        default="Spatial cell-type composition recovery, no sc noise",
        help="Figure title.",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    out_dir = project_root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    scenario_df, per_type_df = build_tables(project_root, args.sample_suffix)
    prefix = args.output_prefix
    scenario_df.to_csv(out_dir / f"{prefix}_scenario.csv", index=False, encoding="utf-8")
    per_type_df.to_csv(out_dir / f"{prefix}_cell_type.csv", index=False, encoding="utf-8")
    summary = (
        scenario_df.groupby("method", observed=True)["composition_recovery"]
        .agg(["count", "mean", "median", "std"])
        .reset_index()
    )
    summary.to_csv(out_dir / f"{prefix}_summary.csv", index=False, encoding="utf-8")
    plot(
        scenario_df,
        out_dir / f"{prefix}_boxplot.png",
        out_dir / f"{prefix}_boxplot.pdf",
        args.title,
    )
    print(f"[OK] wrote: {out_dir / f'{prefix}_boxplot.png'}")
    print(f"[OK] wrote: {out_dir / f'{prefix}_scenario.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
