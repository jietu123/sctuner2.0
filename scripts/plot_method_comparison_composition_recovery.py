#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


SCENARIOS: list[tuple[str, str]] = [
    ("real_brca", "real_brca7_candidate_stable_control"),
    ("real_brca", "real_brca7_candidate_stable_control_missing_epithelial_cells"),
    (
        "real_brca",
        "real_brca7_candidate_stable_control_missing_epithelial_cells_pcs",
    ),
    ("human_lung_5loc", "human_lung_5loc_fine9_clustered_sim"),
    ("human_lung_5loc", "human_lung_5loc_fine9_clustered_sim_missing_at2"),
    ("human_lung_5loc", "human_lung_5loc_fine9_clustered_sim_missing_at2_fibroblast"),
    ("mouse_brain_refined", "mouse_brain_refined7_balanced_clustered_sim"),
    ("mouse_brain_refined", "mouse_brain_refined7_balanced_clustered_sim_missing_micro_fill_ext_l56"),
    ("mouse_brain_refined", "mouse_brain_refined7_balanced_clustered_sim_missing_micro_oligo_2_fill_ext_l56"),
]

COMPOSITE_SCENARIOS: list[tuple[str, str]] = [
    ("real_brca", "real_brca7_endothelial_marker_control_sc_missing_endothelial_cells"),
    (
        "real_brca",
        "real_brca7_endothelial_marker_missing_epithelial_cells_sc_missing_endothelial_cells",
    ),
    (
        "real_brca",
        "real_brca7_endothelial_marker_missing_epithelial_cells_pcs_sc_missing_endothelial_cells",
    ),
    ("human_lung_5loc", "human_lung_5loc_fine9_clustered_sim_sc_missing_b_cell"),
    (
        "human_lung_5loc",
        "human_lung_5loc_fine9_clustered_sim_missing_at2_sc_missing_b_cell",
    ),
    (
        "human_lung_5loc",
        "human_lung_5loc_fine9_clustered_sim_missing_at2_fibroblast_sc_missing_b_cell",
    ),
    ("mouse_brain_refined", "mouse_brain_refined7_balanced_clustered_sim_sc_missing_ext_l56"),
    (
        "mouse_brain_refined",
        "mouse_brain_refined7_balanced_clustered_sim_missing_micro_fill_inh_pvalb_sc_missing_ext_l56",
    ),
    (
        "mouse_brain_refined",
        "mouse_brain_refined7_balanced_clustered_sim_missing_micro_oligo_2_fill_inh_pvalb_sc_missing_ext_l56",
    ),
]


METHODS = [
    ("CytoSPACE", "cytospace_baseline"),
    ("SVTuner", "cytospace_route2"),
    ("Tangram\n(all genes)", "tangram_all"),
    ("Tangram\n(marker genes)", "tangram_marker"),
    ("novoSpaRc", "novosparc"),
    ("SpaOTsc", "spaotsc"),
    ("CellTrek", "celltrek"),
]


PALETTE = {
    "CytoSPACE": "#E83E78",
    "SVTuner": "#2A9D8F",
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


def _prediction_path(
    project_root: Path,
    sample: str,
    method_dir: str,
    route2_stage4_dir: str = "stage4_cytospace_route2",
) -> Path:
    if method_dir == "cytospace_baseline":
        return project_root / "result" / sample / "stage4_cytospace_baseline" / "cytospace_output" / "fractional_abundances_by_spot.csv"
    if method_dir == "cytospace_route2":
        return project_root / "result" / sample / route2_stage4_dir / "cytospace_output" / "fractional_abundances_by_spot.csv"
    return project_root / "result" / sample / "stage4_mapping" / method_dir / "spot_type_fraction.csv"


def _load_fraction(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"fraction file not found: {path}")
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str)
    return df.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def _stage3b_supported_spots(
    project_root: Path,
    sample: str,
    route2_stage4_dir: str,
) -> set[str]:
    blank_path = (
        project_root
        / "result"
        / sample
        / route2_stage4_dir
        / "cytospace_output"
        / "stage3b_blank_spots.csv"
    )
    if not blank_path.exists():
        raise FileNotFoundError(f"Stage3B blank manifest not found: {blank_path}")
    blank = pd.read_csv(blank_path, index_col=0)
    blank_ids = set(blank.index.astype(str))
    truth = _load_fraction(
        project_root
        / "data"
        / "sim"
        / _group_for_sample(sample)
        / sample
        / "sim_truth_spot_type_fraction.csv"
    )
    return set(truth.index.astype(str)).difference(blank_ids)


def _group_for_sample(sample: str) -> str:
    for group, scenario_sample in SCENARIOS + COMPOSITE_SCENARIOS:
        if sample == scenario_sample or sample.startswith(f"{scenario_sample}_"):
            return group
    raise KeyError(f"Could not infer simulation group for sample: {sample}")


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


def _abstention_aware_composition_recovery(
    pred: pd.DataFrame,
    truth: pd.DataFrame,
    blank_spots: set[str],
    unsupported_types: list[str],
    truth_rule: str = "target_dominant",
) -> tuple[float, dict[str, Any]]:
    """Score mapping overlap on every spot and reward only truth-supported abstention."""
    all_spots = truth.index.astype(str)
    type_cols = sorted(set(pred.columns).union(set(truth.columns)))
    p = pred.reindex(index=all_spots, columns=type_cols, fill_value=0.0).astype(float)
    t = truth.reindex(index=all_spots, columns=type_cols, fill_value=0.0).astype(float)
    p = p.div(p.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    t = t.div(t.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)

    unsupported_types = [cell_type for cell_type in unsupported_types if cell_type in t.columns]
    if not unsupported_types:
        raise ValueError(
            "None of the reference-dropped cell types are present in simulation truth: "
            f"{unsupported_types}"
        )
    if truth_rule == "target_positive":
        truth_unsupported = t[unsupported_types].sum(axis=1).gt(0.0)
    elif truth_rule == "target_dominant":
        truth_unsupported = t.idxmax(axis=1).isin(unsupported_types)
    else:
        raise ValueError(f"Unknown abstention truth rule: {truth_rule}")

    predicted_blank = pd.Series(all_spots.isin(blank_spots), index=all_spots)
    correct_blank = predicted_blank & truth_unsupported
    incorrect_blank = predicted_blank & ~truth_unsupported
    spot_overlap = pd.Series(
        np.minimum(p.to_numpy(), t.to_numpy()).sum(axis=1),
        index=all_spots,
        dtype=float,
    )
    spot_score = spot_overlap.copy()
    spot_score.loc[correct_blank] = 1.0
    spot_score.loc[incorrect_blank] = 0.0
    return float(spot_score.mean()), {
        "all_spots": int(len(all_spots)),
        "predicted_blank_spots": int(predicted_blank.sum()),
        "correct_abstention_spots": int(correct_blank.sum()),
        "incorrect_abstention_spots": int(incorrect_blank.sum()),
        "truth_unsupported_spots": int(truth_unsupported.sum()),
        "raw_composition_overlap": float(spot_overlap.mean()),
        "abstention_credit": float(correct_blank.sum() / max(len(all_spots), 1)),
        "abstention_truth_rule": truth_rule,
    }


def build_tables(
    project_root: Path,
    sample_suffix: str = "",
    methods: list[tuple[str, str]] | None = None,
    groups: set[str] | None = None,
    scenarios: list[tuple[str, str]] | None = None,
    route2_stage4_dir: str = "stage4_cytospace_route2",
    exclude_stage3b_blank_spots: bool = False,
    reward_correct_stage3b_abstention: bool = False,
    abstention_truth_rule: str = "target_dominant",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if exclude_stage3b_blank_spots and reward_correct_stage3b_abstention:
        raise ValueError(
            "Supported-region filtering and whole-space abstention credit are mutually exclusive."
        )
    selected_methods = methods or METHODS
    selected_scenarios = scenarios or SCENARIOS
    scenario_rows: list[dict[str, Any]] = []
    cell_type_rows: list[pd.DataFrame] = []
    for group, sample in selected_scenarios:
        if groups is not None and group not in groups:
            continue
        eval_sample = f"{sample}{sample_suffix}"
        truth = _load_fraction(_truth_path(project_root, group, eval_sample))
        unsupported_types: list[str] = []
        if reward_correct_stage3b_abstention:
            sim_info_path = _truth_path(project_root, group, eval_sample).parent / "sim_info.json"
            if not sim_info_path.exists():
                raise FileNotFoundError(f"simulation metadata not found: {sim_info_path}")
            sim_info = json.loads(sim_info_path.read_text(encoding="utf-8"))
            unsupported_types = [
                str(cell_type)
                for cell_type in sim_info.get("sc_reference_drop_types", [])
                if str(cell_type)
            ]
            if not unsupported_types:
                raise ValueError(
                    f"sc_reference_drop_types absent from simulation metadata: {sim_info_path}"
                )
        if exclude_stage3b_blank_spots:
            supported_spots = _stage3b_supported_spots(
                project_root,
                eval_sample,
                route2_stage4_dir,
            )
            truth = truth.loc[truth.index.intersection(supported_spots)].copy()
        for method_label, method_dir in selected_methods:
            pred = _load_fraction(
                _prediction_path(
                    project_root,
                    eval_sample,
                    method_dir,
                    route2_stage4_dir=route2_stage4_dir,
                )
            )
            if exclude_stage3b_blank_spots:
                pred = pred.loc[pred.index.intersection(truth.index)].copy()
            score, per_type = _composition_recovery(pred, truth)
            audit: dict[str, Any] = {
                "all_spots": int(len(truth)),
                "predicted_blank_spots": 0,
                "correct_abstention_spots": 0,
                "incorrect_abstention_spots": 0,
                "truth_unsupported_spots": 0,
                "raw_composition_overlap": score,
                "abstention_credit": 0.0,
                "abstention_truth_rule": "not_applicable",
            }
            if reward_correct_stage3b_abstention and method_dir == "cytospace_route2":
                blank_spots = set(
                    pd.read_csv(
                        project_root
                        / "result"
                        / eval_sample
                        / route2_stage4_dir
                        / "cytospace_output"
                        / "stage3b_blank_spots.csv",
                        index_col=0,
                    ).index.astype(str)
                )
                score, audit = _abstention_aware_composition_recovery(
                    pred,
                    truth,
                    blank_spots,
                    unsupported_types,
                    truth_rule=abstention_truth_rule,
                )
            scenario_rows.append(
                {
                    "group": group,
                    "sample": eval_sample,
                    "method": method_label,
                    "composition_recovery": score,
                    **audit,
                }
            )
            per_type["group"] = group
            per_type["sample"] = eval_sample
            per_type["method"] = method_label
            cell_type_rows.append(per_type)
    scenario_df = pd.DataFrame(scenario_rows)
    per_type_df = pd.concat(cell_type_rows, ignore_index=True)
    order = [m[0] for m in selected_methods]
    scenario_df["method"] = pd.Categorical(scenario_df["method"], order, ordered=True)
    per_type_df["method"] = pd.Categorical(per_type_df["method"], order, ordered=True)
    return scenario_df, per_type_df


def plot(
    df: pd.DataFrame,
    out_png: Path,
    out_pdf: Path | None = None,
    title: str | None = None,
    ylabel: str = "Composition recovery score",
) -> None:
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
    order = [str(x) for x in df["method"].dropna().drop_duplicates().tolist()]
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
    ax.set_ylabel(ylabel, fontsize=12)
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
        default="visualizations/method_comparison/composition_recovery",
        help="Output directory.",
    )
    p.add_argument(
        "--output_prefix",
        default="composition_recovery_7mapping_methods",
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
    p.add_argument(
        "--groups",
        default="",
        help="Comma-separated simulation groups to include. Empty uses all configured groups.",
    )
    p.add_argument(
        "--scenario_preset",
        choices=["standard", "composite"],
        default="standard",
        help="Scenario set to plot.",
    )
    p.add_argument(
        "--route2_stage4_dir",
        default="stage4_cytospace_route2",
        help="Stage4 directory used for the SVTuner + CytoSPACE method.",
    )
    p.add_argument(
        "--exclude_stage3b_blank_spots",
        action="store_true",
        help="Evaluate composition recovery only on spots not flagged as Stage3B unsupported.",
    )
    p.add_argument(
        "--reward_correct_stage3b_abstention",
        action="store_true",
        help=(
            "Evaluate all spots, replacing SVTuner blank-spot overlap with 1 only when "
            "simulation truth confirms a correct abstention."
        ),
    )
    p.add_argument(
        "--abstention_truth_rule",
        choices=["target_dominant", "target_positive"],
        default="target_dominant",
        help="Simulation-truth rule used to validate a Stage3B abstention.",
    )
    p.add_argument(
        "--ylabel",
        default="Composition recovery score",
        help="Y-axis label.",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    project_root = Path(args.project_root).absolute()
    out_dir = project_root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    selected_groups = {x.strip() for x in args.groups.split(",") if x.strip()} or None
    scenarios = COMPOSITE_SCENARIOS if args.scenario_preset == "composite" else SCENARIOS
    scenario_df, per_type_df = build_tables(
        project_root,
        args.sample_suffix,
        methods=METHODS,
        groups=selected_groups,
        scenarios=scenarios,
        route2_stage4_dir=args.route2_stage4_dir,
        exclude_stage3b_blank_spots=args.exclude_stage3b_blank_spots,
        reward_correct_stage3b_abstention=args.reward_correct_stage3b_abstention,
        abstention_truth_rule=args.abstention_truth_rule,
    )
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
        args.ylabel,
    )
    print(f"[OK] wrote: {out_dir / f'{prefix}_boxplot.png'}")
    print(f"[OK] wrote: {out_dir / f'{prefix}_scenario.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
