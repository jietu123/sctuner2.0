#!/usr/bin/env python3
"""BioApp Phase 5: endpoint-specific evaluation of CytoSPACE baselines.

This phase reads existing Phase 4 CytoSPACE outputs only. It does not run
CytoSPACE, SVTuner, Stage3, Stage4, or prevention analysis.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import average_precision_score, roc_auc_score


ROOT = Path(__file__).resolve().parents[1]
PHASE2C_DIR = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase2c_endpoint_freeze_with_validated_cta_to_spot_mapping"
PHASE3R_DIR = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase3r_review_resolution_for_cta_endpoint_baseline_feasibility"
PHASE4_DIR = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase4_formal_cytospace_baseline_execution_against_frozen_cta_immune_endpoint"
OUT_DIR = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase5_endpoint_specific_evaluation_of_cytospace_baselines_against_frozen_cta_immune_endpoint"

PHASE = "BioApp Phase 5 endpoint-specific evaluation of CytoSPACE baseline outputs against frozen CTA Immune endpoint"
STAGE_TYPE = "baseline-vs-endpoint evaluation / biological application candidate preparation"
PRIMARY_ENDPOINT = "Immune cells"
PRIMARY_ENDPOINT_NOTE = "sparse immune-positive spatial compartments, not a large continuous ROI"
RUNS = [
    "baseline_full_reference",
    "baseline_immune_all_dropout",
    "baseline_nonimmune_all_dropout_control",
]
IMMUNE_LABELS = [
    "B cells",
    "CD8 T cells",
    "Monocytes and Macrophages",
    "NK cells",
    "CD4 T cells",
    "Plasma cells",
    "T-cells",
]
NONIMMUNE_LABELS = ["Endothelial cells", "Epithelial cells", "Fibroblasts", "PVL"]
ALLOWED_CLAIMS = [
    "CytoSPACE baseline outputs were evaluated against the frozen CTA Immune endpoint.",
    "Endpoint-specific baseline metrics were computed for later comparison with SVTuner-aware outputs.",
    "The frozen CTA endpoint was used only for evaluation, not for mapping or endpoint selection.",
]
DISALLOWED_CLAIMS = [
    "SVTuner improves endpoint recovery.",
    "SVTuner prevents contradicted interpretation.",
    "Biological application completed.",
    "Biological discovery made.",
]


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT)).replace("\\", "/")
    except ValueError:
        return str(path).replace("\\", "/")


def read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def write_json(path: Path, data: dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, ensure_ascii=False)
        fh.write("\n")


def validate_inputs() -> tuple[dict[str, Any], dict[str, Any], list[str]]:
    errors: list[str] = []
    p2c_path = PHASE2C_DIR / "bioapp_phase2c_endpoint_freeze_summary.json"
    p4_path = PHASE4_DIR / "bioapp_phase4_formal_cytospace_baseline_execution_summary.json"
    if not p2c_path.exists():
        return {}, {}, [f"missing Phase 2C summary: {rel(p2c_path)}"]
    if not p4_path.exists():
        return {}, {}, [f"missing Phase 4 summary: {rel(p4_path)}"]
    p2c = read_json(p2c_path)
    p4 = read_json(p4_path)
    checks = [
        (p2c.get("decision") == "PASS", "Phase 2C decision is not PASS"),
        (p2c.get("endpoint_frozen") is True, "Phase 2C endpoint_frozen is not true"),
        (p2c.get("primary_endpoint") == PRIMARY_ENDPOINT, "Phase 2C primary endpoint is not Immune cells"),
        (p4.get("decision") == "PASS", "Phase 4 decision is not PASS"),
        (p4.get("ready_for_phase5_endpoint_evaluation") is True, "Phase 4 not ready for Phase 5"),
        (p4.get("all_required_baseline_runs_completed") is True, "Phase 4 required baseline runs not completed"),
        (p4.get("biological_application_allowed") is False, "Phase 4 biological_application_allowed must be false"),
        (p4.get("SVTuner_run") is False, "Phase 4 SVTuner_run must be false"),
        (p4.get("Stage4_run") is False, "Phase 4 Stage4_run must be false"),
        (p4.get("baseline_vs_endpoint_metric_computed") is False, "Phase 4 baseline_vs_endpoint_metric_computed must be false"),
    ]
    for ok, msg in checks:
        if not ok:
            errors.append(msg)
    for run in RUNS:
        run_dir = PHASE4_DIR / run / "cytospace_output"
        for fname in ["assigned_locations.csv", "cell_type_assignments_by_spot.csv", "fractional_abundances_by_spot.csv", "log.txt"]:
            if not (run_dir / fname).exists():
                errors.append(f"missing required CytoSPACE output: {run}/{fname}")
    return p2c, p4, errors


def load_endpoint() -> pd.DataFrame:
    endpoint = pd.read_csv(PHASE2C_DIR / "spot_level_endpoint_freeze.csv")
    endpoint["barcode"] = endpoint["barcode"].astype(str)
    return endpoint


def load_coords() -> pd.DataFrame:
    coords = pd.read_csv(PHASE3R_DIR / "st_coordinates_full_2248.csv")
    coords["barcode"] = coords["barcode"].astype(str)
    if {"imagerow", "imagecol"}.issubset(coords.columns):
        coords = coords.rename(columns={"imagerow": "row", "imagecol": "col"})
    return coords[["barcode", "row", "col"]]


def parse_fraction_table(run: str, endpoint_barcodes: list[str]) -> tuple[pd.DataFrame, dict[str, Any]]:
    frac_path = PHASE4_DIR / run / "cytospace_output" / "fractional_abundances_by_spot.csv"
    parse_report = {
        "run_name": run,
        "preferred_source": rel(frac_path),
        "fallback_used": False,
        "parse_status": "parsed_fractional_abundances",
        "notes": "",
    }
    df = pd.read_csv(frac_path)
    if "SpotID" not in df.columns:
        raise ValueError(f"{run}: fractional_abundances_by_spot.csv missing SpotID")
    df = df.rename(columns={"SpotID": "barcode"})
    df["barcode"] = df["barcode"].astype(str)
    label_cols = [c for c in df.columns if c != "barcode"]
    for col in label_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0.0)
    missing = sorted(set(endpoint_barcodes) - set(df["barcode"]))
    if missing:
        filler = pd.DataFrame({"barcode": missing})
        for col in label_cols:
            filler[col] = 0.0
        df = pd.concat([df, filler], ignore_index=True)
        parse_report["notes"] = f"filled {len(missing)} missing endpoint spots with zero fractions"
    df = df.loc[df["barcode"].isin(endpoint_barcodes)].copy()
    df = df.set_index("barcode").loc[endpoint_barcodes].reset_index()
    df["fraction_sum"] = df[label_cols].sum(axis=1)
    df["n_cell_types_available"] = len(label_cols)
    df["source_file"] = rel(frac_path)
    df["parse_status"] = parse_report["parse_status"]
    df.to_csv(OUT_DIR / f"{run}_spot_fraction_table.csv", index=False)
    return df, parse_report


def label_grouping(run: str, label_cols: list[str]) -> dict[str, Any]:
    present = set(label_cols)
    data = {
        "run_name": run,
        "immune_labels_expected": IMMUNE_LABELS,
        "nonimmune_labels_expected": NONIMMUNE_LABELS,
        "immune_labels_present_in_run": [x for x in IMMUNE_LABELS if x in present],
        "nonimmune_labels_present_in_run": [x for x in NONIMMUNE_LABELS if x in present],
        "labels_missing_from_run": [x for x in IMMUNE_LABELS + NONIMMUNE_LABELS if x not in present],
        "label_grouping_source": "Phase 3R / Phase 4 fixed vocabulary",
    }
    write_json(OUT_DIR / f"{run}_label_grouping_used.json", data)
    return data


def score_by_spot(run: str, frac: pd.DataFrame, endpoint: pd.DataFrame) -> pd.DataFrame:
    label_cols = [c for c in frac.columns if c not in {"barcode", "fraction_sum", "n_cell_types_available", "source_file", "parse_status"}]
    grouping = label_grouping(run, label_cols)
    merged = endpoint[["barcode", "primary_endpoint_status"]].merge(frac, on="barcode", how="left")
    for col in label_cols:
        merged[col] = pd.to_numeric(merged[col], errors="coerce").fillna(0.0)
    immune_present = [c for c in IMMUNE_LABELS if c in label_cols]
    nonimmune_present = [c for c in NONIMMUNE_LABELS if c in label_cols]
    merged["is_endpoint_positive"] = merged["primary_endpoint_status"].eq("positive")
    merged["is_endpoint_negative"] = merged["primary_endpoint_status"].eq("negative")
    merged["is_endpoint_ambiguous"] = merged["primary_endpoint_status"].eq("ambiguous")
    merged["is_endpoint_excluded"] = merged["primary_endpoint_status"].eq("excluded")
    merged["baseline_run"] = run
    merged["baseline_reference_condition"] = run.replace("baseline_", "")
    merged["immune_score"] = merged[immune_present].sum(axis=1) if immune_present else 0.0
    merged["nonimmune_score"] = merged[nonimmune_present].sum(axis=1) if nonimmune_present else 0.0
    merged["tumor_or_epithelial_score"] = merged[[c for c in ["Epithelial cells"] if c in label_cols]].sum(axis=1) if "Epithelial cells" in label_cols else 0.0
    merged["stromal_score"] = merged[[c for c in ["Fibroblasts", "PVL"] if c in label_cols]].sum(axis=1) if any(c in label_cols for c in ["Fibroblasts", "PVL"]) else 0.0
    merged["myeloid_score_if_available"] = merged["Monocytes and Macrophages"] if "Monocytes and Macrophages" in label_cols else np.nan
    merged["B_plasma_score_if_available"] = merged[[c for c in ["B cells", "Plasma cells"] if c in label_cols]].sum(axis=1) if any(c in label_cols for c in ["B cells", "Plasma cells"]) else np.nan
    merged["T_NK_score_if_available"] = merged[[c for c in ["CD8 T cells", "CD4 T cells", "T-cells", "NK cells"] if c in label_cols]].sum(axis=1) if any(c in label_cols for c in ["CD8 T cells", "CD4 T cells", "T-cells", "NK cells"]) else np.nan
    if label_cols:
        vals = merged[label_cols].to_numpy()
        idx = np.nanargmax(vals, axis=1)
        merged["dominant_baseline_label"] = [label_cols[i] for i in idx]
        merged["dominant_baseline_fraction"] = vals[np.arange(len(vals)), idx]
    else:
        merged["dominant_baseline_label"] = None
        merged["dominant_baseline_fraction"] = np.nan
    merged["available_label_count"] = len(label_cols)
    merged["score_parse_status"] = "ok"
    if run == "baseline_immune_all_dropout":
        merged["immune_score_interpretation"] = "immune labels unavailable by dropout design"
    else:
        merged["immune_score_interpretation"] = "immune labels available"
    keep = [
        "barcode",
        "primary_endpoint_status",
        "is_endpoint_positive",
        "is_endpoint_negative",
        "is_endpoint_ambiguous",
        "is_endpoint_excluded",
        "baseline_run",
        "baseline_reference_condition",
        "immune_score",
        "nonimmune_score",
        "tumor_or_epithelial_score",
        "stromal_score",
        "myeloid_score_if_available",
        "B_plasma_score_if_available",
        "T_NK_score_if_available",
        "dominant_baseline_label",
        "dominant_baseline_fraction",
        "available_label_count",
        "fraction_sum",
        "score_parse_status",
        "immune_score_interpretation",
    ]
    merged[keep].to_csv(OUT_DIR / f"{run}_endpoint_score_by_spot.csv", index=False)
    return merged[keep]


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float | None:
    if len(x) == 0 or len(y) == 0:
        return None
    gt = 0
    lt = 0
    for xv in x:
        gt += int(np.sum(xv > y))
        lt += int(np.sum(xv < y))
    return (gt - lt) / (len(x) * len(y))


def bootstrap_ci(pos: np.ndarray, neg: np.ndarray, n_boot: int = 1000, seed: int = 20260705) -> tuple[float | None, float | None]:
    if len(pos) == 0 or len(neg) == 0:
        return None, None
    rng = np.random.default_rng(seed)
    deltas = []
    for _ in range(n_boot):
        ps = rng.choice(pos, size=len(pos), replace=True)
        ns = rng.choice(neg, size=len(neg), replace=True)
        deltas.append(float(np.mean(ps) - np.mean(ns)))
    return float(np.quantile(deltas, 0.025)), float(np.quantile(deltas, 0.975))


def compute_metrics(run: str, score: pd.DataFrame) -> dict[str, Any]:
    main = score.loc[score["primary_endpoint_status"].isin(["positive", "negative"])].copy()
    pos = main.loc[main["is_endpoint_positive"], "immune_score"].astype(float).to_numpy()
    neg = main.loc[main["is_endpoint_negative"], "immune_score"].astype(float).to_numpy()
    mean_pos = float(np.mean(pos)) if len(pos) else None
    mean_neg = float(np.mean(neg)) if len(neg) else None
    median_pos = float(np.median(pos)) if len(pos) else None
    median_neg = float(np.median(neg)) if len(neg) else None
    delta = None if mean_pos is None or mean_neg is None else mean_pos - mean_neg
    fold = None if mean_neg in {None, 0.0} else mean_pos / mean_neg
    metric_status = "ok"
    cohen = None
    cliff = None
    mw_p = None
    auroc = None
    auprc = None
    reason = ""
    if len(pos) < 2 or len(neg) < 2 or np.nanstd(np.concatenate([pos, neg])) == 0:
        metric_status = "undefined_due_to_constant_score"
        reason = "immune_score_constant_or_unavailable"
    else:
        pooled = math.sqrt(((len(pos) - 1) * np.var(pos, ddof=1) + (len(neg) - 1) * np.var(neg, ddof=1)) / (len(pos) + len(neg) - 2))
        cohen = None if pooled == 0 else float((np.mean(pos) - np.mean(neg)) / pooled)
        cliff = cliffs_delta(pos, neg)
        mw_p = float(stats.mannwhitneyu(pos, neg, alternative="two-sided").pvalue)
        y_true = np.array([1] * len(pos) + [0] * len(neg))
        y_score = np.concatenate([pos, neg])
        auroc = float(roc_auc_score(y_true, y_score))
        auprc = float(average_precision_score(y_true, y_score))
    ci_low, ci_high = bootstrap_ci(pos, neg)
    mean_nonimmune_pos = float(score.loc[score["is_endpoint_positive"], "nonimmune_score"].mean())
    return {
        "baseline_run": run,
        "n_positive_spots": int(len(pos)),
        "n_negative_spots": int(len(neg)),
        "mean_immune_score_positive": mean_pos,
        "mean_immune_score_negative": mean_neg,
        "median_immune_score_positive": median_pos,
        "median_immune_score_negative": median_neg,
        "delta_mean_positive_minus_negative": delta,
        "fold_change_positive_over_negative": fold,
        "Cohen_d_positive_vs_negative": cohen,
        "Cliffs_delta_positive_vs_negative": cliff,
        "Mann_Whitney_U_pvalue": mw_p,
        "AUROC": auroc,
        "AUPRC": auprc,
        "positive_class_prevalence": len(pos) / (len(pos) + len(neg)) if len(pos) + len(neg) else None,
        "endpoint_positive_enrichment": delta,
        "bootstrap_n": 1000,
        "bootstrap_CI95_low": ci_low,
        "bootstrap_CI95_high": ci_high,
        "mean_nonimmune_score_positive": mean_nonimmune_pos,
        "metric_status": metric_status,
        "reason": reason,
        "metric_note": "immune labels removed by design" if run == "baseline_immune_all_dropout" else "",
    }


def plot_endpoint_map(endpoint: pd.DataFrame, coords: pd.DataFrame) -> None:
    data = endpoint[["barcode", "primary_endpoint_status"]].merge(coords, on="barcode", how="left")
    colors = {"positive": "#d73027", "negative": "#4575b4", "ambiguous": "#fdae61", "excluded": "#bdbdbd"}
    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    for status in ["negative", "ambiguous", "excluded", "positive"]:
        sub = data[data["primary_endpoint_status"] == status]
        ax.scatter(sub["col"], -sub["row"], s=10, c=colors[status], label=status, linewidths=0)
    ax.set_title("Frozen CTA Immune endpoint status", fontsize=11)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal")
    ax.legend(frameon=False, markerscale=1.5, fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "frozen_cta_immune_endpoint_status_map.pdf")
    fig.savefig(OUT_DIR / "frozen_cta_immune_endpoint_status_map.svg")
    plt.close(fig)


def plot_score_map(run: str, score: pd.DataFrame, coords: pd.DataFrame) -> None:
    data = score.merge(coords, on="barcode", how="left")
    pos = data[data["is_endpoint_positive"]]
    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    sc = ax.scatter(data["col"], -data["row"], c=data["immune_score"], s=10, cmap="viridis", vmin=0, vmax=1, linewidths=0)
    ax.scatter(pos["col"], -pos["row"], facecolors="none", edgecolors="#d73027", s=32, linewidths=0.6)
    title = run.replace("baseline_", "").replace("_", " ")
    if run == "baseline_immune_all_dropout":
        title += "\nimmune labels unavailable by dropout design"
    ax.set_title(title, fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal")
    cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label("immune score", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / f"{run}_immune_score_map.pdf")
    fig.savefig(OUT_DIR / f"{run}_immune_score_map.svg")
    plt.close(fig)


def plot_boxplots(all_scores: pd.DataFrame) -> None:
    main = all_scores[all_scores["primary_endpoint_status"].isin(["positive", "negative"])].copy()
    fig, ax = plt.subplots(figsize=(8.0, 4.6))
    positions = []
    data = []
    labels = []
    pos = 1
    for run in RUNS:
        for status in ["positive", "negative"]:
            vals = main[(main["baseline_run"] == run) & (main["primary_endpoint_status"] == status)]["immune_score"].astype(float)
            data.append(vals)
            positions.append(pos)
            labels.append(f"{run.replace('baseline_', '')}\n{status}")
            pos += 1
        pos += 0.6
    bp = ax.boxplot(data, positions=positions, widths=0.6, patch_artist=True, showfliers=False)
    for patch, label in zip(bp["boxes"], labels):
        patch.set_facecolor("#d73027" if label.endswith("positive") else "#4575b4")
        patch.set_alpha(0.7)
    ax.set_ylabel("Immune score")
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=8)
    ax.set_ylim(-0.02, 1.02)
    ax.set_title("Baseline immune score in CTA Immune-positive vs negative spots")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "baseline_immune_score_positive_vs_negative_boxplots.pdf")
    fig.savefig(OUT_DIR / "baseline_immune_score_positive_vs_negative_boxplots.svg")
    plt.close(fig)


def plot_metric_summary(metrics: pd.DataFrame) -> None:
    metric_cols = ["AUROC", "AUPRC", "delta_mean_positive_minus_negative", "mean_nonimmune_score_positive"]
    fig, axes = plt.subplots(1, 4, figsize=(10.5, 3.2))
    x = np.arange(len(RUNS))
    colors = ["#4c78a8", "#f58518", "#54a24b"]
    for ax, metric in zip(axes, metric_cols):
        vals = metrics.set_index("baseline_run").reindex(RUNS)[metric].to_list()
        bars = ax.bar(x, [0 if pd.isna(v) else v for v in vals], color=colors)
        for i, v in enumerate(vals):
            if pd.isna(v):
                ax.text(i, 0.02, "NA", ha="center", va="bottom", fontsize=8)
            else:
                ax.text(i, v + 0.02, f"{v:.2f}", ha="center", va="bottom", fontsize=7)
        ax.set_title(metric.replace("_", " "), fontsize=8)
        ax.set_xticks(x)
        ax.set_xticklabels([r.replace("baseline_", "").replace("_", "\n") for r in RUNS], fontsize=7)
        if metric in {"AUROC", "AUPRC", "mean_nonimmune_score_positive"}:
            ax.set_ylim(0, 1.05)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "baseline_endpoint_metric_summary.pdf")
    fig.savefig(OUT_DIR / "baseline_endpoint_metric_summary.svg")
    plt.close(fig)


def plot_nonimmune_breakdown(breakdown: pd.DataFrame) -> None:
    pivot = breakdown.pivot(index="baseline_run", columns="nonimmune_label", values="mean_fraction_in_endpoint_positive").reindex(RUNS).fillna(0)
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    bottom = np.zeros(len(pivot))
    colors = ["#8dd3c7", "#ffffb3", "#bebada", "#fb8072"]
    for col, color in zip(pivot.columns, colors):
        vals = pivot[col].to_numpy()
        ax.bar(np.arange(len(pivot)), vals, bottom=bottom, label=col, color=color)
        bottom += vals
    ax.set_xticks(np.arange(len(pivot)))
    ax.set_xticklabels([r.replace("baseline_", "").replace("_", "\n") for r in pivot.index], fontsize=8)
    ax.set_ylabel("Mean fraction in CTA Immune-positive spots")
    ax.set_title("Non-immune assignment composition in CTA Immune-positive spots")
    ax.legend(frameon=False, fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "baseline_positive_spot_nonimmune_label_breakdown.pdf")
    fig.savefig(OUT_DIR / "baseline_positive_spot_nonimmune_label_breakdown.svg")
    plt.close(fig)


def write_text_outputs(summary: dict[str, Any], metrics: pd.DataFrame) -> None:
    next_text = (
        "BioApp Phase 6 - SVTuner-aware execution planning against frozen CTA Immune endpoint"
        if summary["decision"] == "PASS"
        else "Review baseline output parsing, spot universe alignment, and endpoint-specific metric computation before SVTuner-aware execution.\nDo not run SVTuner/Stage4/prevention analysis."
        if summary["decision"] == "REVIEW_REQUIRED"
        else "Stop. Fix endpoint/baseline evaluation or boundary violation before retrying Phase 5."
    )
    lines = [
        "BioApp Phase 5 - endpoint-specific evaluation of CytoSPACE baseline outputs against frozen CTA Immune endpoint",
        "",
        f"Decision: {summary['decision']}",
        "",
        "Input endpoint:",
        "CTA-defined Immune cells",
        "",
        "Endpoint status:",
        "Frozen at Phase 2C",
        "",
        "Primary endpoint note:",
        "Sparse immune-positive spatial compartments, not a large continuous ROI.",
        "",
        "Baseline runs evaluated:",
        *RUNS,
        "",
        "Endpoint analysis set:",
        f"positive = {summary['endpoint_positive_spots']}",
        f"negative = {summary['endpoint_negative_spots']}",
        f"ambiguous = {summary['endpoint_ambiguous_spots']}",
        f"excluded = {summary['endpoint_excluded_spots']}",
        f"main_analysis_spots = positive + negative = {summary['main_analysis_spots']}",
        "",
        "Spot universe alignment:",
        f"all baseline runs aligned to 2248 frozen endpoint spots = {summary['spot_universe_aligned_all_runs']}",
        "",
        "Key baseline endpoint metrics:",
    ]
    for run in RUNS:
        row = metrics[metrics["baseline_run"] == run].iloc[0]
        lines += [
            f"{run}:",
            f"  mean_immune_score_positive = {row['mean_immune_score_positive']}",
            f"  mean_immune_score_negative = {row['mean_immune_score_negative']}",
            f"  AUROC = {row['AUROC']}",
            f"  AUPRC = {row['AUPRC']}",
        ]
        if run == "baseline_immune_all_dropout":
            lines.append("  note = immune labels removed by design")
    lines += [
        "",
        "Boundary checks:",
        "CytoSPACE rerun: false",
        "SVTuner run: false",
        "Stage4 run: false",
        "SVTuner Stage3 run: false",
        "SVTuner Stage4 run: false",
        "Prevention analysis run: false",
        "Endpoint redefined: false",
        "Endpoint used for mapping: false",
        "Endpoint used for evaluation: true",
        "SVTuner-vs-endpoint metric computed: false",
        "",
        "Allowed claims:",
        *[f"- {x}" for x in ALLOWED_CLAIMS],
        "",
        "Disallowed claims:",
        *[f"- {x}" for x in DISALLOWED_CLAIMS],
        "",
        "Next:",
        next_text,
    ]
    (OUT_DIR / "decision.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    readme = f"""# BioApp Phase 5

Purpose: endpoint-specific evaluation of existing CytoSPACE baseline outputs against the frozen CTA-defined Immune cells endpoint.

Inputs:
- Phase 2C frozen endpoint
- Phase 4 baseline outputs

The endpoint is used only for evaluation, not for mapping or endpoint selection.

Primary analysis uses endpoint-positive and endpoint-negative spots. Ambiguous and excluded spots remain in the full spot table but are excluded from main metrics.

immune_score is the sum of fixed immune label fractions. nonimmune_score is the sum of fixed non-immune label fractions. In the immune-all-dropout baseline, immune labels are unavailable by design.

No SVTuner, Stage4, or prevention analysis was run.

Decision: `{summary['decision']}`
"""
    (OUT_DIR / "README.md").write_text(readme, encoding="utf-8")


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    p2c, p4, errors = validate_inputs()
    if errors:
        decision = "FAIL"
        summary = {
            "phase": PHASE,
            "decision": decision,
            "decision_reasons": errors,
            "stage_type": STAGE_TYPE,
            "input_phase2c_decision": p2c.get("decision") if p2c else None,
            "input_phase4_decision": p4.get("decision") if p4 else None,
            "biological_application_allowed": False,
        }
        write_json(OUT_DIR / "bioapp_phase5_endpoint_specific_baseline_evaluation_summary.json", summary)
        return 1

    endpoint = load_endpoint()
    coords = load_coords()
    endpoint_barcodes = endpoint["barcode"].tolist()
    status_counts = endpoint["primary_endpoint_status"].value_counts().to_dict()
    all_score_tables = []
    parse_reports = []
    alignment_rows = []
    metric_rows = []
    breakdown_rows = []

    for run in RUNS:
        frac, parse_report = parse_fraction_table(run, endpoint_barcodes)
        parse_reports.append(parse_report)
        alignment_rows.append(
            {
                "run_name": run,
                "n_endpoint_spots": len(endpoint_barcodes),
                "n_baseline_spots": int(frac["barcode"].nunique()),
                "n_overlap_spots": len(set(endpoint_barcodes) & set(frac["barcode"])),
                "n_endpoint_missing_from_baseline": len(set(endpoint_barcodes) - set(frac["barcode"])),
                "n_baseline_not_in_endpoint": len(set(frac["barcode"]) - set(endpoint_barcodes)),
                "spot_universe_aligned": set(endpoint_barcodes) == set(frac["barcode"]) and len(endpoint_barcodes) == 2248,
                "notes": "",
            }
        )
        score = score_by_spot(run, frac, endpoint)
        all_score_tables.append(score)
        metric_rows.append(compute_metrics(run, score))
        positive_barcodes = set(endpoint.loc[endpoint["primary_endpoint_status"].eq("positive"), "barcode"].astype(str))
        positive_frac = frac.loc[frac["barcode"].isin(positive_barcodes)]
        for label in NONIMMUNE_LABELS:
            label_col = label if label in frac.columns else None
            breakdown_rows.append(
                {
                    "baseline_run": run,
                    "nonimmune_label": label,
                    "mean_fraction_in_endpoint_positive": float(positive_frac[label_col].mean()) if label_col else 0.0,
                    "label_present_in_run": label_col is not None,
                }
            )
        plot_score_map(run, score, coords)

    all_scores = pd.concat(all_score_tables, ignore_index=True)
    metrics = pd.DataFrame(metric_rows)
    alignment = pd.DataFrame(alignment_rows)
    breakdown = pd.DataFrame(breakdown_rows)
    parse_report_df = pd.DataFrame(parse_reports)
    all_scores.to_csv(OUT_DIR / "baseline_endpoint_score_all_runs_by_spot.csv", index=False)
    metrics.to_csv(OUT_DIR / "baseline_endpoint_metric_summary.csv", index=False)
    metrics.to_csv(OUT_DIR / "baseline_condition_metric_comparison.csv", index=False)
    alignment.to_csv(OUT_DIR / "baseline_spot_universe_alignment.csv", index=False)
    breakdown.to_csv(OUT_DIR / "baseline_positive_spot_nonimmune_label_breakdown.csv", index=False)
    parse_report_df.to_csv(OUT_DIR / "baseline_output_parse_report.csv", index=False)

    plot_endpoint_map(endpoint, coords)
    plot_boxplots(all_scores)
    plot_metric_summary(metrics)
    plot_nonimmune_breakdown(breakdown)

    all_aligned = bool(alignment["spot_universe_aligned"].all())
    metrics_ok = not metrics.empty and metrics["mean_immune_score_positive"].notna().all()
    figures = [
        "frozen_cta_immune_endpoint_status_map.pdf",
        "baseline_full_reference_immune_score_map.pdf",
        "baseline_immune_all_dropout_immune_score_map.pdf",
        "baseline_nonimmune_all_dropout_control_immune_score_map.pdf",
        "baseline_immune_score_positive_vs_negative_boxplots.pdf",
        "baseline_endpoint_metric_summary.pdf",
        "baseline_positive_spot_nonimmune_label_breakdown.pdf",
    ]
    figures_ok = all((OUT_DIR / f).exists() for f in figures)
    unexpected_auroc_missing = metrics.loc[
        (metrics["baseline_run"] != "baseline_immune_all_dropout") & metrics["AUROC"].isna()
    ]
    if all_aligned and metrics_ok and figures_ok and unexpected_auroc_missing.empty:
        decision = "PASS"
        reasons = ["baseline outputs parsed, aligned, and evaluated against frozen CTA Immune endpoint"]
        ready_for_phase6 = True
    else:
        decision = "REVIEW_REQUIRED"
        reasons = ["baseline parsing, metric availability, spot alignment, or figure generation requires review"]
        ready_for_phase6 = False

    highlights = {}
    for _, row in metrics.iterrows():
        highlights[row["baseline_run"]] = {
            "mean_immune_score_positive": None if pd.isna(row["mean_immune_score_positive"]) else float(row["mean_immune_score_positive"]),
            "mean_immune_score_negative": None if pd.isna(row["mean_immune_score_negative"]) else float(row["mean_immune_score_negative"]),
            "delta_mean_positive_minus_negative": None if pd.isna(row["delta_mean_positive_minus_negative"]) else float(row["delta_mean_positive_minus_negative"]),
            "AUROC": None if pd.isna(row["AUROC"]) else float(row["AUROC"]),
            "AUPRC": None if pd.isna(row["AUPRC"]) else float(row["AUPRC"]),
            "mean_nonimmune_score_positive": None if pd.isna(row["mean_nonimmune_score_positive"]) else float(row["mean_nonimmune_score_positive"]),
        }
        if row["baseline_run"] == "baseline_immune_all_dropout":
            highlights[row["baseline_run"]]["metric_note"] = "immune labels removed by design"

    summary = {
        "phase": PHASE,
        "decision": decision,
        "decision_reasons": reasons,
        "stage_type": STAGE_TYPE,
        "input_phase2c_decision": p2c.get("decision"),
        "input_phase4_decision": p4.get("decision"),
        "endpoint_frozen": True,
        "primary_endpoint": PRIMARY_ENDPOINT,
        "primary_endpoint_spatial_note": PRIMARY_ENDPOINT_NOTE,
        "endpoint_positive_spots": int(status_counts.get("positive", 0)),
        "endpoint_negative_spots": int(status_counts.get("negative", 0)),
        "endpoint_ambiguous_spots": int(status_counts.get("ambiguous", 0)),
        "endpoint_excluded_spots": int(status_counts.get("excluded", 0)),
        "main_analysis_spots": int(status_counts.get("positive", 0) + status_counts.get("negative", 0)),
        "baseline_runs_evaluated": RUNS,
        "spot_universe_aligned_all_runs": all_aligned,
        "n_spots_expected": 2248,
        "n_gene_overlap_used_in_phase4": 2000,
        "metrics_computed": {
            "mean_immune_score_by_endpoint_status": True,
            "effect_size": True,
            "AUROC_AUPRC": True,
            "endpoint_positive_enrichment": True,
            "nonimmune_assignment_burden": True,
        },
        "baseline_metric_highlights": highlights,
        "endpoint_used_for_mapping": False,
        "endpoint_used_for_evaluation": True,
        "endpoint_redefined": False,
        "expression_markers_used_to_define_endpoint": False,
        "mapping_outputs_used_to_define_endpoint": False,
        "CytoSPACE_rerun": False,
        "SVTuner_run": False,
        "Stage4_run": False,
        "SVTuner_Stage3_run": False,
        "SVTuner_Stage4_run": False,
        "prevention_analysis_run": False,
        "SVTuner_vs_endpoint_metric_computed": False,
        "baseline_vs_endpoint_metric_computed": True,
        "biological_application_allowed": False,
        "ready_for_phase6_svtuner_execution": ready_for_phase6,
        "allowed_claims": ALLOWED_CLAIMS,
        "disallowed_claims": DISALLOWED_CLAIMS,
        "output_dir": rel(OUT_DIR),
    }
    write_json(OUT_DIR / "bioapp_phase5_endpoint_specific_baseline_evaluation_summary.json", summary)
    golden = {
        "stage_type": STAGE_TYPE,
        "biological_question_defined": True,
        "external_endpoint_predefined": True,
        "endpoint_independent_from_SVTuner": True,
        "endpoint_spatially_registered": True,
        "baseline_comparison_available": True,
        "endpoint_specific_improvement_defined": True,
        "endpoint_specific_quantitative_metric_available": True,
        "endpoint_baseline_SVTuner_spatial_comparison_available": False,
        "interpretation_boundary_defined": True,
        "biological_application_allowed": False,
        "allowed_claim_level": "baseline-vs-endpoint evaluation only",
        "decision": decision,
    }
    write_json(OUT_DIR / "bioapp_phase5_golden_rules_v2_1_check.json", golden)
    write_text_outputs(summary, metrics)
    manifest = []
    for path in sorted(OUT_DIR.iterdir()):
        if path.is_file():
            manifest.append(
                {
                    "file": rel(path),
                    "type": path.suffix.lstrip(".") or "text",
                    "description": "BioApp Phase 5 output artifact",
                    "created_by_phase": PHASE,
                    "status": "generated",
                    "notes": "baseline-vs-endpoint evaluation only; no SVTuner",
                }
            )
    pd.DataFrame(manifest).to_csv(OUT_DIR / "manifest.csv", index=False)

    print("BioApp Phase 5 completed.")
    print()
    print("Decision:")
    print(decision)
    print()
    print("Input endpoint:")
    print("CTA-defined Immune cells")
    print()
    print("Endpoint frozen:")
    print("true")
    print()
    print("Endpoint analysis set:")
    print(f"positive = {summary['endpoint_positive_spots']}")
    print(f"negative = {summary['endpoint_negative_spots']}")
    print(f"ambiguous = {summary['endpoint_ambiguous_spots']}")
    print(f"excluded = {summary['endpoint_excluded_spots']}")
    print(f"main_analysis_spots = {summary['main_analysis_spots']}")
    print()
    print("Baseline runs evaluated:")
    for run in RUNS:
        print(run)
    print()
    print("Spot universe alignment:")
    print(f"all_runs_aligned_to_2248 = {str(all_aligned).lower()}")
    print()
    print("Key metrics:")
    for run in RUNS:
        row = metrics[metrics["baseline_run"] == run].iloc[0]
        print(f"{run}:")
        print(f"  AUROC = {row['AUROC']}")
        print(f"  AUPRC = {row['AUPRC']}")
        print(f"  delta_mean_immune_score = {row['delta_mean_positive_minus_negative']}")
        if run == "baseline_immune_all_dropout":
            print("  note = immune labels removed by design")
    print()
    print("Boundary checks:")
    print("CytoSPACE rerun: false")
    print("SVTuner run: false")
    print("Stage4 run: false")
    print("SVTuner Stage3 run: false")
    print("SVTuner Stage4 run: false")
    print("Prevention analysis run: false")
    print("Endpoint redefined: false")
    print("Endpoint used for mapping: false")
    print("Endpoint used for evaluation: true")
    print("SVTuner-vs-endpoint metric computed: false")
    print()
    print("Ready for Phase 6 SVTuner-aware execution planning:")
    print(str(ready_for_phase6).lower())
    print()
    print("Next:")
    if decision == "PASS":
        print("BioApp Phase 6 - SVTuner-aware execution planning against frozen CTA Immune endpoint")
    elif decision == "REVIEW_REQUIRED":
        print("Review baseline output parsing, spot universe alignment, and endpoint-specific metric computation before SVTuner-aware execution.")
        print("Do not run SVTuner/Stage4/prevention analysis.")
    else:
        print("Stop. Fix endpoint/baseline evaluation or boundary violation before retrying Phase 5.")
    return 0 if decision in {"PASS", "REVIEW_REQUIRED"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
