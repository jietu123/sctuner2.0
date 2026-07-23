from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = PROJECT_ROOT / "visualizations" / "cytospace_fig2e_stage3_profile_mask"
RESULT_DIR = PROJECT_ROOT / "result" / "cytospace_fig2e_stage3_profile_mask"
GENE_SET_DIR = PROJECT_ROOT / "data" / "raw" / "cytospace_fig2c_melanoma" / "prepared" / "gene_sets" / "generated"
EPS = 1.0e-12


@dataclass(frozen=True)
class Scenario:
    sample: str
    dataset_label: str
    scrna_label: str
    masked_label: str
    tumor_labels: tuple[str, ...]
    platform_color: str


SCENARIOS = [
    Scenario(
        sample="cytospace_fig2c_melanoma_mel1_rep2_screen_mask_macrophages",
        dataset_label="Melanoma 1",
        scrna_label="Tirosh et al.",
        masked_label="Macrophages",
        tumor_labels=("Melanoma", "Melanoma cells"),
        platform_color="#4c78a8",
    ),
    Scenario(
        sample="cytospace_fig2c_melanoma_mel2_rep1_profile_mask_nk_cells",
        dataset_label="Melanoma 2",
        scrna_label="Tirosh et al.",
        masked_label="NK cells",
        tumor_labels=("Melanoma", "Melanoma cells"),
        platform_color="#4c78a8",
    ),
    Scenario(
        sample="cytospace_fig2d_tme_brca_er_her2_fresh_frozen_profile_mask_t_cells",
        dataset_label="HER2+ BRCA\nFF",
        scrna_label="Wu et al.",
        masked_label="T-cells",
        tumor_labels=("Epithelial cells", "Cancer epithelial cells", "Malignant cells"),
        platform_color="#54a24b",
    ),
    Scenario(
        sample="cytospace_fig2d_tme_brca_her2_ffpe_profile_mask_t_cells",
        dataset_label="HER2+ BRCA\nFFPE",
        scrna_label="Wu et al.",
        masked_label="T-cells",
        tumor_labels=("Epithelial cells", "Cancer epithelial cells", "Malignant cells"),
        platform_color="#e45756",
    ),
    Scenario(
        sample="cytospace_fig2d_tme_brca_tnbc_fresh_frozen_profile_mask_fibroblasts",
        dataset_label="TNBC BRCA",
        scrna_label="Wu et al.",
        masked_label="Fibroblasts",
        tumor_labels=("Epithelial cells", "Cancer epithelial cells", "Malignant cells"),
        platform_color="#54a24b",
    ),
    Scenario(
        sample="cytospace_fig2d_tme_crc_fresh_frozen_profile_mask_monocytes_and_macrophages",
        dataset_label="CRC",
        scrna_label="Lee et al.",
        masked_label="Mono/Mac",
        tumor_labels=("Epithelial cells", "Cancer epithelial cells", "Malignant cells"),
        platform_color="#54a24b",
    ),
]

METHODS = [
    ("CytoSPACE", "CytoSPACE"),
    ("SVTuner + CytoSPACE", "SVTuner +\nCytoSPACE"),
]

CELL_TYPES = [
    ("CD4 T cells", ("CD4 T cells",), "CD4 T cells", "ecotyper_cd4_t_cells", "#6f3fb1"),
    ("CD8 T cells", ("CD8 T cells",), "CD8 T cells", "ecotyper_cd8_t_cells", "#4aa6d8"),
    (
        "Mono/Mac",
        ("Monocytes and Macrophages", "Macrophages", "Mono/Mac"),
        "Mono/Mac",
        "ecotyper_monocytes_and_macrophages",
        "#7a6a28",
    ),
    ("B cells", ("B cells",), "B cells", "ecotyper_b_cells", "#f2c339"),
    ("Plasma cells", ("Plasma cells", "PCs"), "Plasma cells", "ecotyper_pcs", "#2aa6a1"),
]
FEATURES = ["CE9", "CE10"]
GENE_SET_FILE_OVERRIDES = {("ecotyper_cd8_t_cells", "CE10"): "ecotyper_cd8_cells_ce10.txt"}


def _read_gene_set(prefix: str, feature: str) -> set[str]:
    filename = GENE_SET_FILE_OVERRIDES.get((prefix, feature), f"{prefix}_{feature.lower()}.txt")
    path = GENE_SET_DIR / filename
    return {x.strip() for x in path.read_text(encoding="utf-8-sig").splitlines() if x.strip()}


def _processed_dir(sample: str) -> Path:
    direct = PROJECT_ROOT / "data" / "processed" / sample
    if (direct / "stage1_preprocess" / "exported" / "sc_expression_normalized.csv").exists():
        return direct
    matches = list((PROJECT_ROOT / "data" / "processed").glob(f"*/{sample}"))
    for p in matches:
        if (p / "stage1_preprocess" / "exported" / "sc_expression_normalized.csv").exists():
            return p
    raise FileNotFoundError(f"processed directory not found for {sample}")


def _cytospace_dir(sample: str, method: str) -> Path:
    base = PROJECT_ROOT / "result" / sample
    if method == "CytoSPACE":
        patterns = ["stage4_cytospace_baseline_strictfig2c", "stage4_cytospace_baseline_profile300", "stage4_cytospace_baseline*"]
    else:
        patterns = ["stage4_cytospace_route2_strictfig2c", "stage4_cytospace_route2_profile300", "stage4_cytospace_route2*"]
    for pattern in patterns:
        for p in sorted(base.glob(pattern)):
            out = p / "cytospace_output"
            if (out / "assigned_locations.csv").exists():
                return out
    raise FileNotFoundError(f"assigned_locations not found for {sample} {method}")


def _load_stage3_evidence(scenario: Scenario) -> dict:
    summary_path = PROJECT_ROOT / "result" / scenario.sample / "stage3_typematch" / "stage3_summary.json"
    data = json.loads(summary_path.read_text(encoding="utf-8"))
    auto = data.get("auto_missing_detection", {})
    masked = data.get("masked_missing_detection", {})
    return {
        "sample": scenario.sample,
        "dataset_label": scenario.dataset_label,
        "masked_label": scenario.masked_label,
        "summary_path": str(summary_path.relative_to(PROJECT_ROOT)),
        "auto_missing_types": auto.get("auto_missing_types", []),
        "masked_missing_types": masked.get("masked_missing_types", []),
        "white_list_used": False,
    }


def _mean_nearest5(target_xy: np.ndarray, tumor_xy: np.ndarray) -> np.ndarray:
    diff = target_xy[:, None, :] - tumor_xy[None, :, :]
    dist = np.sqrt(np.sum(diff * diff, axis=2))
    k = min(5, dist.shape[1])
    return np.partition(dist, kth=k - 1, axis=1)[:, :k].mean(axis=1)


def _running_es(stats: pd.Series, gene_set: set[str]) -> tuple[float, int]:
    stats = stats.replace([np.inf, -np.inf], np.nan).dropna().sort_values(ascending=False)
    hits = stats.index.to_series().isin(gene_set).to_numpy()
    n_hits = int(hits.sum())
    n_miss = int(len(stats) - n_hits)
    if n_hits == 0 or n_miss == 0:
        return float("nan"), 0
    weights = np.abs(stats.to_numpy(dtype=np.float64))
    hit_norm = float(weights[hits].sum())
    if hit_norm <= EPS:
        increments = np.where(hits, 1.0 / n_hits, -1.0 / n_miss)
    else:
        increments = np.where(hits, weights / hit_norm, -1.0 / n_miss)
    running = np.cumsum(increments)
    peak_idx = int(np.nanargmax(running))
    return float(running[peak_idx]), int(peak_idx + 1)


def _null_nes(stats: pd.Series, n_genes: int, observed_es: float, nperm: int, seed: int) -> tuple[float, float]:
    if not np.isfinite(observed_es):
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    genes = stats.replace([np.inf, -np.inf], np.nan).dropna().index.to_numpy()
    if n_genes >= len(genes) or n_genes == 0:
        return float("nan"), float("nan")
    null = np.empty(nperm, dtype=np.float64)
    for i in range(nperm):
        selected = set(rng.choice(genes, size=n_genes, replace=False).tolist())
        null[i], _ = _running_es(stats, selected)
    denom = float(np.mean(np.abs(null))) + EPS
    nes = observed_es / denom
    pval = (float(np.sum(null >= observed_es)) + 1.0) / (float(nperm) + 1.0)
    return float(nes), float(pval)


def _select_cell_type(loc: pd.DataFrame, aliases: tuple[str, ...]) -> tuple[str, pd.DataFrame]:
    for alias in aliases:
        target = loc[loc["CellType"].eq(alias)].copy()
        if len(target) > 0:
            return alias, target
    return aliases[0], loc.iloc[0:0].copy()


def compute_metrics(nperm: int = 1000) -> pd.DataFrame:
    rows = []
    gene_sets = {
        (prefix, feature): _read_gene_set(prefix, feature)
        for _, _, _, prefix, _ in CELL_TYPES
        for feature in FEATURES
    }
    for scenario_idx, scenario in enumerate(SCENARIOS):
        processed = _processed_dir(scenario.sample)
        sc_expr = pd.read_csv(processed / "stage1_preprocess" / "exported" / "sc_expression_normalized.csv", index_col=0)
        genes = pd.Index(sc_expr.columns.astype(str))
        for method_idx, (method, display_method) in enumerate(METHODS):
            loc = pd.read_csv(_cytospace_dir(scenario.sample, method) / "assigned_locations.csv")
            loc["OriginalCID"] = loc["OriginalCID"].astype(str)
            loc["CellType"] = loc["CellType"].astype(str)
            coord_cols = ["row", "col"] if {"row", "col"}.issubset(loc.columns) else ["X", "Y"]
            tumor = loc[loc["CellType"].isin(scenario.tumor_labels)]
            if len(tumor) < 5:
                raise RuntimeError(f"{scenario.sample} {method}: too few tumor mapped cells ({len(tumor)})")
            tumor_xy = tumor[coord_cols].to_numpy(dtype=np.float64)
            for cell_idx, (canonical, aliases, cell_label, gs_prefix, _) in enumerate(CELL_TYPES):
                actual_cell_type, target = _select_cell_type(loc, aliases)
                if len(target) < 20:
                    for feature in FEATURES:
                        rows.append(
                            {
                                "scenario": scenario.sample,
                                "dataset_label": scenario.dataset_label,
                                "scrna_label": scenario.scrna_label,
                                "masked_label": scenario.masked_label,
                                "method": method,
                                "display_method": display_method,
                                "cell_type": canonical,
                                "actual_cell_type": actual_cell_type,
                                "cell_label": cell_label,
                                "feature": feature,
                                "NES": np.nan,
                                "P-value": np.nan,
                                "ES": np.nan,
                                "peak_rank": 0,
                                "n_cells": int(len(target)),
                                "n_close": 0,
                                "n_far": 0,
                                "n_genes_used": 0,
                            }
                        )
                    continue
                target_xy = target[coord_cols].to_numpy(dtype=np.float64)
                dist = _mean_nearest5(target_xy, tumor_xy)
                close = dist <= float(np.median(dist))
                close_ids = target.loc[close, "OriginalCID"].to_numpy()
                far_ids = target.loc[~close, "OriginalCID"].to_numpy()
                close_expr = sc_expr.reindex(close_ids).fillna(0.0)
                far_expr = sc_expr.reindex(far_ids).fillna(0.0)
                stats = pd.Series(
                    close_expr.mean(axis=0).to_numpy(dtype=np.float64)
                    - far_expr.mean(axis=0).to_numpy(dtype=np.float64),
                    index=genes,
                )
                for feature_idx, feature in enumerate(FEATURES):
                    gene_set = pd.Index(list(gene_sets[(gs_prefix, feature)])).intersection(genes)
                    es, peak_rank = _running_es(stats, set(gene_set.tolist()))
                    nes, pval = _null_nes(
                        stats,
                        len(gene_set),
                        es,
                        nperm=nperm,
                        seed=20260526 + scenario_idx * 1000 + method_idx * 100 + cell_idx * 10 + feature_idx,
                    )
                    rows.append(
                        {
                            "scenario": scenario.sample,
                            "dataset_label": scenario.dataset_label,
                            "scrna_label": scenario.scrna_label,
                            "masked_label": scenario.masked_label,
                            "method": method,
                            "display_method": display_method,
                            "cell_type": canonical,
                            "actual_cell_type": actual_cell_type,
                            "cell_label": cell_label,
                            "feature": feature,
                            "NES": nes,
                            "P-value": pval,
                            "ES": es,
                            "peak_rank": peak_rank,
                            "n_cells": int(len(target)),
                            "n_close": int(close.sum()),
                            "n_far": int((~close).sum()),
                            "n_genes_used": int(len(gene_set)),
                        }
                    )
    return pd.DataFrame(rows)


def _p_size(p: float) -> float:
    if not np.isfinite(p):
        return 0
    return 28 if p < 0.01 else 14


def _nes_color(nes: float) -> str:
    if not np.isfinite(nes):
        return "#e8e8e8"
    return "#f28e2b" if nes >= 0 else "#2f6fab"


def plot(metrics: pd.DataFrame) -> Path:
    mpl.rcParams["svg.fonttype"] = "none"
    mpl.rcParams["font.family"] = "Arial"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    RESULT_DIR.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(4.9, 2.55), dpi=300)
    fig.patch.set_facecolor("white")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(
        0.50,
        0.965,
        "CE9 spatial enrichment after Stage3-detected profile masking",
        fontsize=6.8,
        weight="bold",
        ha="center",
        va="top",
    )
    ax.text(0.075, 0.800, "ST dataset", fontsize=4.3, ha="left", va="bottom")
    ax.text(0.215, 0.800, "scRNA-seq", fontsize=4.3, ha="left", va="bottom")
    ax.text(0.452, 0.835, "CytoSPACE", fontsize=4.5, ha="center", va="bottom")
    ax.text(0.628, 0.835, "SVTuner + CytoSPACE", fontsize=4.5, ha="center", va="bottom")

    x_dataset = 0.075
    x_scrna = 0.215
    x_platform_bar = 0.330
    method_x0 = {
        "CytoSPACE": 0.380,
        "SVTuner + CytoSPACE": 0.556,
    }
    dx = 0.029
    y_top = 0.728
    dy_scenario = 0.052
    cell_w = dx
    cell_h = dy_scenario

    for method, x0 in method_x0.items():
        for cidx, (_, _, _, _, color) in enumerate(CELL_TYPES):
            ax.add_patch(
                plt.Rectangle(
                    (x0 + cidx * dx - 0.013, 0.792),
                    0.026,
                    0.022,
                    facecolor=color,
                    edgecolor="white",
                    linewidth=0.25,
                )
        )

    for sidx, scenario in enumerate(SCENARIOS):
        y_mid = y_top - sidx * dy_scenario
        ax.text(x_dataset, y_mid, scenario.dataset_label, fontsize=4.1, ha="left", va="center", linespacing=0.72)
        ax.add_patch(
            plt.Rectangle(
                (x_platform_bar - 0.006, y_mid - cell_h * 0.34),
                0.012,
                cell_h * 0.68,
                facecolor=scenario.platform_color,
                edgecolor="none",
            )
        )
        ax.text(x_scrna, y_mid, scenario.scrna_label, fontsize=4.1, ha="left", va="center")
        for method, _ in METHODS:
            y = y_mid
            feature = "CE9"
            for cidx, (canonical, _, _, _, _) in enumerate(CELL_TYPES):
                row = metrics[
                    (metrics["scenario"] == scenario.sample)
                    & (metrics["method"] == method)
                    & (metrics["cell_type"] == canonical)
                    & (metrics["feature"] == feature)
                ].iloc[0]
                x = method_x0[method] + cidx * dx
                nes = float(row["NES"]) if pd.notna(row["NES"]) else np.nan
                pval = float(row["P-value"]) if pd.notna(row["P-value"]) else np.nan
                if not np.isfinite(nes):
                    ax.add_patch(
                        plt.Rectangle(
                            (x - cell_w / 2 + 0.002, y - cell_h / 2 + 0.002),
                            cell_w - 0.004,
                            cell_h - 0.004,
                            facecolor="#d9d9d9",
                            edgecolor="none",
                            zorder=1,
                        )
                    )
                else:
                    ax.scatter([x], [y], s=_p_size(pval), facecolor=_nes_color(nes), edgecolor="none", zorder=4)

    row_bottom = y_top - (len(SCENARIOS) - 1) * dy_scenario
    for x0 in method_x0.values():
        left = x0 - cell_w / 2
        right = x0 + (len(CELL_TYPES) - 0.5) * cell_w
        bottom = row_bottom - cell_h / 2
        top = y_top + cell_h / 2
        ax.add_patch(plt.Rectangle((left, bottom), right - left, top - bottom, fill=False, edgecolor="#222222", linewidth=0.65, zorder=5))
        for k in range(1, len(CELL_TYPES)):
            xline = x0 + (k - 0.5) * cell_w
            ax.plot([xline, xline], [bottom, top], color="#d0d0d0", linewidth=0.35, zorder=2)
        for k in range(1, len(SCENARIOS)):
            yline = y_top - (k - 0.5) * cell_h
            ax.plot([left, right], [yline, yline], color="#d0d0d0", linewidth=0.35, zorder=2)

    legend_x = 0.755
    legend_y = 0.690
    ax.text(legend_x, legend_y + 0.055, "Cell types", fontsize=4.7, ha="left", va="bottom")
    for i, (_, _, label, _, color) in enumerate(CELL_TYPES):
        y = legend_y + 0.018 - i * 0.037
        ax.add_patch(plt.Rectangle((legend_x, y - 0.010), 0.020, 0.020, facecolor=color, edgecolor="none"))
        ax.text(legend_x + 0.026, y, label, fontsize=4.2, ha="left", va="center")

    ax.text(0.160, 0.285, "ST platform", fontsize=4.8, ha="left", va="bottom")
    platform_items = [("Legacy ST", "#4c78a8"), ("Visium FFPE", "#e45756"), ("Visium fresh-frozen", "#54a24b")]
    for i, (label, color) in enumerate(platform_items):
        y = 0.260 - i * 0.030
        ax.add_patch(plt.Rectangle((0.165, y - 0.009), 0.018, 0.018, facecolor=color, edgecolor="none"))
        ax.text(0.190, y, label, fontsize=4.4, ha="left", va="center")

    ax.text(0.430, 0.285, "Enrichment score", fontsize=4.8, ha="center", va="bottom")
    ax.text(0.382, 0.247, "Close to\n tumor", fontsize=4.1, ha="right", va="center", linespacing=0.85)
    ax.scatter([0.405], [0.247], s=31, color="#f28e2b", edgecolor="none")
    ax.plot([0.428, 0.465], [0.247, 0.247], color="#555555", linewidth=0.75)
    ax.scatter([0.488], [0.247], s=31, color="#2f6fab", edgecolor="none")
    ax.text(0.512, 0.247, "Far from\n tumor", fontsize=4.1, ha="left", va="center", linespacing=0.85)

    ax.text(0.630, 0.285, "P value scale", fontsize=4.8, ha="left", va="bottom")
    for i, (label, size) in enumerate([("P < 0.01", 42), ("P > 0.01", 24), ("N/A", 13)]):
        y = 0.260 - i * 0.030
        face = "white" if label != "N/A" else "#e8e8e8"
        edge = "#111111" if label != "N/A" else "#8f8f8f"
        ax.scatter([0.645], [y], s=size, facecolor=face, edgecolor=edge, linewidth=0.55)
        ax.text(0.668, y, label, fontsize=4.4, ha="left", va="center")

    out_png = OUT_DIR / "fig2e_stage3_profile_mask_baseline_vs_route2.png"
    out_pdf = OUT_DIR / "fig2e_stage3_profile_mask_baseline_vs_route2.pdf"
    fig.savefig(out_png, dpi=600, bbox_inches="tight", pad_inches=0.005)
    fig.savefig(out_pdf, bbox_inches="tight", pad_inches=0.005)
    plt.close(fig)
    return out_png


def main() -> int:
    evidence = [_load_stage3_evidence(s) for s in SCENARIOS]
    metrics = compute_metrics(nperm=1000)
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    (RESULT_DIR / "stage3_detection_evidence.json").write_text(
        json.dumps(evidence, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    metrics.to_csv(RESULT_DIR / "fig2e_stage3_profile_mask_metrics.csv", index=False)
    out = plot(metrics)
    print(f"[OK] wrote: {out}")
    print(metrics[["dataset_label", "masked_label", "method", "cell_type", "feature", "NES", "P-value", "n_cells"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


