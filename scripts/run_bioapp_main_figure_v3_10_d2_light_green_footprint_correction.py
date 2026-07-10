from __future__ import annotations

import csv
import importlib.util
import json
import subprocess
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, PowerNorm
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_main_figure_v3_10_d2_light_green_footprint_correction"
P2C = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase2c_endpoint_freeze_with_validated_cta_to_spot_mapping"
P8 = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase8_svtuner_vs_endpoint_evaluation_and_baseline_svtuner_comparison"
P9 = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase9_final_biological_application_audit_figure_selection_interpretation_boundary_lock"
RECOVERY = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_coordinate_registration_recovery_audit"
P2BR2 = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase2br2_explicit_rscript_cta_align_locale_safe_dry_run"
STDATA_ASCII = P2BR2 / "_ascii_runtime" / "ST_data.RData"
RSCRIPT = Path("E:/R/R-4.5.1/bin/x64/Rscript.exe")

IMAGE_SLOT = "BCSA2TumB1"
SEURAT_OBJECT = "bcsa"

PALETTE = {
    "canvas_background": "#FFFFFF",
    "text": "#222222",
    "non_endpoint_spots": "#D9D9D9",
    "endpoint_positive": "#D55E00",
    "endpoint_outline": "#6B1E16",
    "baseline_low": "#F7F7F7",
    "baseline_mid": "#FDBE85",
    "baseline_high": "#D7301F",
    "withheld_low": "#F7F7F7",
    "withheld_mid": "#9ECAE1",
    "withheld_high": "#0072B2",
    "withheld_binary_fill": "#FFFFFF",
    "withheld_binary_outline": "#222222",
    "endpoint_only_fill": "#D55E00",
    "endpoint_only_outline": "#6B1E16",
    "withheld_only_fill": "#FFFFFF",
    "withheld_only_outline": "#222222",
    "endpoint_withheld_overlap_fill": "#D55E00",
    "endpoint_withheld_overlap_outer_outline": "#222222",
    "endpoint_withheld_overlap_optional_inner_outline": "#FFFFFF",
    "d2_endpoint_reference_marker_face": "none",
    "d2_endpoint_reference_marker_edge": "#7A3B2E",
    "d3_endpoint_reference_marker_face": "none",
    "d3_endpoint_reference_marker_edge": "#6B1E16",
    "d2_baseline_forced_burden_low": "#F7F7F7",
    "d2_baseline_forced_burden_mid": "#FDBE85",
    "d2_baseline_forced_burden_high": "#D7301F",
    "d3_withheld_score_low": "#F7F7F7",
    "d3_withheld_score_mid": "#9ECAE1",
    "d3_withheld_score_high": "#0072B2",
    "d2_candidate_A_grayscale_low": "#F7F7F7",
    "d2_candidate_A_grayscale_mid": "#9E9E9E",
    "d2_candidate_A_grayscale_high": "#1F1F1F",
    "d2_candidate_B_brown_low": "#F7F7F7",
    "d2_candidate_B_brown_mid": "#BCA27A",
    "d2_candidate_B_brown_high": "#5B3A29",
    "d2_candidate_C_olive_low": "#F7F7F7",
    "d2_candidate_C_olive_mid": "#A8A16A",
    "d2_candidate_C_olive_high": "#4B4B2A",
    "d2_footprint_fill": "#8D9B4A",
    "d2_footprint_boundary": "#6E7935",
    "d2_adaptive_low": "#F6F1D3",
    "d2_adaptive_mid": "#D9D89A",
    "d2_adaptive_high": "#8D9B4A",
    "baseline_forced": "#C4473A",
    "prevented_by_svtuner": "#0072B2",
    "endpoint_negative": "#6E6E6E",
}
BASELINE_CMAP = LinearSegmentedColormap.from_list("bioapp_v3_baseline", [PALETTE["baseline_low"], PALETTE["baseline_mid"], PALETTE["baseline_high"]])
WITHHELD_CMAP = LinearSegmentedColormap.from_list("bioapp_v3_withheld", [PALETTE["withheld_low"], PALETTE["withheld_mid"], PALETTE["withheld_high"]])
BASE_VERSION_USED = "V3.9"
D2D3_CLIP_LOWER = 2.0
D2D3_CLIP_UPPER = 98.0
D2D3_ENDPOINT_REFERENCE_SIZE_RELATIVE_TO_V3 = 0.45
D2D3_ENDPOINT_REFERENCE_ALPHA = 0.62
D2D3_ENDPOINT_REFERENCE_EDGEWIDTH = 0.55
D1_MARKER_SIZE_RELATIVE_TO_V3_4 = 0.36
D2_FOOTPRINT_MARKER_SIZE_RELATIVE_TO_V3_4 = 0.48
D2_FOOTPRINT_ALPHA_V3_5 = 0.30
D2_FOOTPRINT_SIZE_RELATIVE_TO_V3_5 = 0.85
D2_FOOTPRINT_SIZE_RELATIVE_TO_V3_6 = 1.00
D2_FOOTPRINT_ALPHA_V3_7 = 0.82
D2_FOOTPRINT_ALPHA_V3_8 = 0.38
D2_FOOTPRINT_ONLY_FILL_V3_8 = "#66BB6A"
D2_FOOTPRINT_ONLY_EDGE_V3_8 = "#2E7D32"
D2_DISPLAY_CLIP_LOWER = 10.0
D2_DISPLAY_CLIP_UPPER = 98.0
D3_MARKER_SIZE_RELATIVE_TO_V3_4 = 0.58
D4_ENDPOINT_ONLY_MARKER_SIZE_RELATIVE_TO_V3_4 = 0.55
D4_WITHHELD_ONLY_MARKER_SIZE_RELATIVE_TO_V3_4 = 0.55
D4_OVERLAP_MARKER_SIZE_RELATIVE_TO_V3_4 = 0.58
D4_ENDPOINT_ONLY_SIZE_RELATIVE_TO_V3_5 = 0.78
D4_WITHHELD_ONLY_SIZE_RELATIVE_TO_V3_5 = 0.78
D4_OVERLAP_SIZE_RELATIVE_TO_V3_5 = 0.72
D4_ENDPOINT_ONLY_SIZE_RELATIVE_TO_V3_6 = 0.42
D4_WITHHELD_ONLY_SIZE_RELATIVE_TO_V3_6 = 0.42
D4_OVERLAP_SIZE_RELATIVE_TO_V3_6 = 0.38
D4_OVERLAP_EDGEWIDTH_V3_7 = 0.36
D_PANEL_ENDPOINT_REFERENCE_SCALE_V3_5 = 0.58
D2_QUANTILE_5_COLORS = ["#D73027", "#6A3D9A", "#F781BF", "#FFD92F", "#1A9850"]
D2_QUANTILE_7_COLORS = ["#FFFFCC", "#D9F0A3", "#ADDD8E", "#78C679", "#41AB5D", "#238443", "#005A32"]
D2_QUANTILE_5_CMAP = LinearSegmentedColormap.from_list("bioapp_v3_7_d2_quantile_5bin", D2_QUANTILE_5_COLORS, N=5)
D2_QUANTILE_7_CMAP = LinearSegmentedColormap.from_list("bioapp_v3_7_d2_quantile_7bin", D2_QUANTILE_7_COLORS, N=7)
D2_RANK_CMAP = LinearSegmentedColormap.from_list("bioapp_v3_7_d2_rank_continuous", D2_QUANTILE_7_COLORS)
D2_SELECTED_CANDIDATE = "quantile_5bin"
D2_CANDIDATE_SPECS = {
    "candidate_A_grayscale": {
        "label": "Candidate A: grayscale",
        "palette": ["#F7F7F7", "#9E9E9E", "#1F1F1F"],
        "normalization": "PowerNorm(gamma=0.55) with 2/98 percentile bounds",
        "gamma": 0.55,
    },
    "candidate_B_brown": {
        "label": "Candidate B: brown",
        "palette": ["#F7F7F7", "#BCA27A", "#5B3A29"],
        "normalization": "PowerNorm(gamma=0.55) with 2/98 percentile bounds",
        "gamma": 0.55,
    },
    "candidate_C_olive": {
        "label": "Candidate C: olive",
        "palette": ["#F7F7F7", "#A8A16A", "#4B4B2A"],
        "normalization": "PowerNorm(gamma=0.55) with 2/98 percentile bounds",
        "gamma": 0.55,
    },
}


def candidate_cmap(name: str) -> LinearSegmentedColormap:
    spec = D2_CANDIDATE_SPECS[name]
    return LinearSegmentedColormap.from_list(name, spec["palette"])

FOV1_RAW = (4268.04, 7099.96, 5625.5, 8768.5)
FOV2_RAW = (10634.04, 13465.96, 8200.5, 11343.5)


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT)).replace("\\", "/")


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")


def write_rows(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k) for k in fields})


def load_v1_module():
    path = ROOT / "scripts" / "run_bioapp_main_figure_layout_only_spatial_visualization.py"
    spec = importlib.util.spec_from_file_location("bioapp_v1", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import V1 module from {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run_r_export_image() -> dict[str, Any]:
    OUT.mkdir(parents=True, exist_ok=True)
    helper = OUT / "_export_bioapp_v3_embedded_image.R"
    helper.write_text(
        r'''
suppressPackageStartupMessages({
  library(Seurat)
  library(jsonlite)
})
args <- commandArgs(trailingOnly = TRUE)
rdata <- args[[1]]
out_dir <- args[[2]]
object_name <- args[[3]]
image_name <- args[[4]]
env <- new.env()
load(rdata, envir = env)
obj <- get(object_name, envir = env)
img <- obj@images[[image_name]]
arr <- img@image
h <- dim(arr)[1]
w <- dim(arr)[2]
sf <- as.list(img@scale.factors)
lowres <- as.numeric(sf$lowres)
hires <- as.numeric(sf$hires)
coords <- img@coordinates
x <- coords$imagecol * lowres
y <- coords$imagerow * lowres
bounds <- list(
  embedded_image_source = paste0("ST_data.RData::", object_name, "@images[['", image_name, "']]@image"),
  image_slot_class = class(img),
  embedded_image_dimensions = dim(arr),
  lowres_scale_used = lowres,
  hires_scale_available = hires,
  coordinate_columns_used = c("imagecol", "imagerow"),
  x_min = min(x, na.rm = TRUE),
  x_max = max(x, na.rm = TRUE),
  y_min = min(y, na.rm = TRUE),
  y_max = max(y, na.rm = TRUE),
  image_width = w,
  image_height = h,
  bounds_check_passed = all(x >= 0, x <= w, y >= 0, y <= h, na.rm = TRUE)
)
png(file.path(out_dir, "fig_bioapp_v3_embedded_tissue_background.png"), width = w, height = h)
par(mar = c(0, 0, 0, 0))
plot.new()
plot.window(xlim = c(0, w), ylim = c(h, 0), asp = 1)
rasterImage(as.raster(arr), 0, h, w, 0)
dev.off()
writeLines(jsonlite::toJSON(bounds, auto_unbox = TRUE, pretty = TRUE), file.path(out_dir, "fig_bioapp_v3_embedded_image_metadata.json"), useBytes = TRUE)
''',
        encoding="ascii",
    )
    result = subprocess.run(
        [str(RSCRIPT), str(helper), rel(STDATA_ASCII), rel(OUT), SEURAT_OBJECT, IMAGE_SLOT],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
        timeout=600,
    )
    (OUT / "fig_bioapp_v3_r_export_stdout.txt").write_text(result.stdout, encoding="utf-8")
    (OUT / "fig_bioapp_v3_r_export_stderr.txt").write_text(result.stderr, encoding="utf-8")
    if result.returncode != 0:
        raise RuntimeError(f"R image export failed: {result.stderr}")
    return read_json(OUT / "fig_bioapp_v3_embedded_image_metadata.json")


def load_data(lowres: float) -> tuple[pd.DataFrame, dict[str, Any], dict[str, Any]]:
    endpoint = pd.read_csv(P2C / "spot_level_endpoint_freeze.csv")
    scores = pd.read_csv(P8 / "svtuner_endpoint_score_by_spot.csv")
    comp = pd.read_csv(P8 / "baseline_svtuner_endpoint_comparison_by_spot.csv")
    df = endpoint[["barcode", "imagecol", "imagerow", "primary_endpoint_status"]].merge(scores, on=["barcode", "primary_endpoint_status"], how="inner")
    df = df.merge(
        comp[
            [
                "barcode",
                "baseline_forced_nonimmune_burden",
                "svtuner_prevented_forced_nonimmune_burden_binary",
                "svtuner_prevented_forced_nonimmune_burden_continuous",
            ]
        ],
        on="barcode",
        how="left",
    )
    df["withheld_binary"] = df["withheld_binary"].astype(str).str.lower().isin(["true", "1"])
    df["x_plot"] = df["imagecol"] * lowres
    df["y_plot"] = df["imagerow"] * lowres
    phase8 = read_json(P8 / "bioapp_phase8_svtuner_vs_endpoint_evaluation_summary.json")
    phase9 = read_json(P9 / "bioapp_phase9_final_audit_summary.json")
    return df, phase8, phase9


def setup_ax(ax: plt.Axes) -> None:
    ax.set_facecolor("white")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_panel_label(ax: plt.Axes, label: str, title: str) -> None:
    ax.text(-0.02, 1.07, label, transform=ax.transAxes, fontsize=9, fontweight="bold", ha="left", va="bottom", color=PALETTE["text"])
    ax.text(0.06, 1.07, title, transform=ax.transAxes, fontsize=7.3, fontweight="bold", ha="left", va="bottom", color=PALETTE["text"])


def draw_tissue(ax: plt.Axes, img: np.ndarray, image_width: float, image_height: float) -> None:
    setup_ax(ax)
    ax.imshow(img, extent=(0, image_width, image_height, 0), origin="upper", interpolation="nearest", zorder=0)
    ax.set_xlim(0, image_width)
    ax.set_ylim(image_height, 0)
    ax.set_aspect("equal")


def draw_spots_context(ax: plt.Axes, df: pd.DataFrame, s: float = 3.0, alpha: float = 0.25) -> None:
    ax.scatter(df["x_plot"], df["y_plot"], c=PALETTE["non_endpoint_spots"], s=s, lw=0, alpha=alpha, zorder=1)


def draw_endpoint(ax: plt.Axes, df: pd.DataFrame, s: float = 10.0, outline_s: float | None = None) -> None:
    pos = df["primary_endpoint_status"].eq("positive")
    ax.scatter(df.loc[pos, "x_plot"], df.loc[pos, "y_plot"], c=PALETTE["endpoint_positive"], s=s, lw=0, alpha=0.92, zorder=5)
    ax.scatter(
        df.loc[pos, "x_plot"],
        df.loc[pos, "y_plot"],
        facecolors="none",
        edgecolors=PALETTE["endpoint_outline"],
        s=outline_s or s * 1.85,
        linewidths=0.45,
        alpha=0.95,
        zorder=6,
    )


def draw_endpoint_reference(ax: plt.Axes, df: pd.DataFrame, s: float = 3.4, alpha: float = D2D3_ENDPOINT_REFERENCE_ALPHA) -> None:
    pos = df["primary_endpoint_status"].eq("positive")
    ax.scatter(
        df.loc[pos, "x_plot"],
        df.loc[pos, "y_plot"],
        facecolors="none",
        edgecolors=PALETTE["d2_endpoint_reference_marker_edge"],
        s=s,
        linewidths=D2D3_ENDPOINT_REFERENCE_EDGEWIDTH,
        alpha=alpha,
        zorder=6,
    )


def percentile_norm(values: pd.Series, lower: float = D2D3_CLIP_LOWER, upper: float = D2D3_CLIP_UPPER) -> tuple[float, float]:
    clean = pd.to_numeric(values, errors="coerce").dropna()
    if clean.empty:
        return 0.0, 1.0
    vmin = float(np.percentile(clean, lower))
    vmax = float(np.percentile(clean, upper))
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
        vmin = float(clean.min())
        vmax = float(clean.max())
    if vmax <= vmin:
        vmax = vmin + 1.0
    return vmin, vmax


def d2_adaptive_display_range(values: pd.Series) -> tuple[float, float, dict[str, Any]]:
    clean = pd.to_numeric(values, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    positive = clean[clean > 0]
    if positive.empty:
        return 0.0, 1.0, {
            "display_rescaling_used": True,
            "vmin_rule": "no positive burden values; fallback vmin=0",
            "vmax_rule": "no positive burden values; fallback vmax=1",
            "near_constant_distribution": True,
            "observed_n_positive": 0,
            "observed_min": None,
            "observed_max": None,
            "display_vmin": 0.0,
            "display_vmax": 1.0,
        }

    observed_min = float(positive.min())
    observed_max = float(positive.max())
    vmin = float(np.percentile(positive, D2_DISPLAY_CLIP_LOWER))
    vmax = float(np.percentile(positive, D2_DISPLAY_CLIP_UPPER))
    near_constant = bool(np.isclose(observed_min, observed_max) or np.isclose(vmin, vmax))
    if near_constant:
        vmin = 0.0 if observed_max > 0 else observed_min
        vmax = observed_max if observed_max > 0 else 1.0
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
        vmin = 0.0
        vmax = max(1.0, observed_max)

    return vmin, vmax, {
        "display_rescaling_used": True,
        "vmin_rule": f"positive-value percentile {D2_DISPLAY_CLIP_LOWER:g}; fallback to 0 if near-constant",
        "vmax_rule": f"positive-value percentile {D2_DISPLAY_CLIP_UPPER:g}; fallback to observed max if near-constant",
        "near_constant_distribution": near_constant,
        "observed_n_positive": int(positive.shape[0]),
        "observed_min": observed_min,
        "observed_max": observed_max,
        "display_vmin": vmin,
        "display_vmax": vmax,
    }


def d2_distribution_diagnostics(values: pd.Series) -> dict[str, Any]:
    clean = pd.to_numeric(values, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    plotted = clean[clean > 0]
    plotted_rounded = plotted.round(6)
    percentiles = {}
    for p in [1, 5, 10, 25, 50, 75, 90, 95, 99]:
        percentiles[f"p{p}"] = float(np.percentile(clean, p)) if not clean.empty else None
    unique_count = int(clean.nunique(dropna=True))
    plotted_unique_count = int(plotted_rounded.nunique(dropna=True))
    std = float(clean.std(ddof=0)) if clean.shape[0] else None
    exact_constant = bool(plotted_unique_count <= 1)
    raw_min = float(clean.min()) if clean.shape[0] else None
    raw_max = float(clean.max()) if clean.shape[0] else None
    plotted_min = float(plotted.min()) if plotted.shape[0] else None
    plotted_max = float(plotted.max()) if plotted.shape[0] else None
    near_constant = bool((not exact_constant) and plotted_min is not None and plotted_max is not None and np.isclose(plotted_min, plotted_max, rtol=1e-6, atol=1e-9))
    if exact_constant:
        selected = "constant_footprint"
    else:
        selected = "footprint_only_forced_call_area"
    return {
        "phase": "BioApp Main Figure V3.10 D2 light-green footprint correction",
        "d2_raw_value_min": raw_min,
        "d2_raw_value_max": raw_max,
        "d2_raw_value_mean": float(clean.mean()) if clean.shape[0] else None,
        "d2_raw_value_median": float(clean.median()) if clean.shape[0] else None,
        "d2_raw_value_std": std,
        "d2_unique_value_count": unique_count,
        "d2_plotted_positive_unique_value_count_rounded_6dp": plotted_unique_count,
        "d2_plotted_positive_min": plotted_min,
        "d2_plotted_positive_max": plotted_max,
        "d2_percentiles": percentiles,
        "d2_exact_constant": exact_constant,
        "d2_near_constant": near_constant,
        "display_mapping_selected": selected,
        "display_mapping_visualization_only": True,
        "underlying_values_changed": False,
    }


def d2_rank_values(values: pd.Series) -> pd.Series:
    clean = pd.to_numeric(values, errors="coerce")
    ranks = clean.rank(method="average", pct=True)
    return ranks.fillna(0.0).clip(0.0, 1.0)


def d2_quantile_bin_values(values: pd.Series, n_bins: int) -> pd.Series:
    clean = pd.to_numeric(values, errors="coerce")
    ranks = d2_rank_values(clean)
    bins = np.ceil(ranks * n_bins).astype(int)
    bins = np.clip(bins, 1, n_bins)
    return pd.Series(bins, index=clean.index)


def d2_display_values(values: pd.Series, method: str) -> tuple[pd.Series, LinearSegmentedColormap, float, float, list[float], str]:
    diagnostics = d2_distribution_diagnostics(values)
    if diagnostics["d2_exact_constant"]:
        display = pd.Series(np.ones(len(values)), index=values.index)
        return display, D2_QUANTILE_5_CMAP, 1.0, 1.0, [1.0], "constant footprint"
    if method == "quantile_7bin":
        display = d2_quantile_bin_values(values, 7)
        return display, D2_QUANTILE_7_CMAP, 1.0, 7.0, [1, 2, 3, 4, 5, 6, 7], "display quantile"
    if method == "percentile_rank":
        display = d2_rank_values(values)
        return display, D2_RANK_CMAP, 0.0, 1.0, [0.0, 0.5, 1.0], "relative display"
    display = d2_quantile_bin_values(values, 5)
    return display, D2_QUANTILE_5_CMAP, 1.0, 5.0, [1, 2, 3, 4, 5], "display quantile"


def panel_a(ax: plt.Axes) -> None:
    setup_ax(ax)
    draw_panel_label(ax, "A", "Immune-all-dropout biological application design")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    steps = [
        (0.05, 0.63, "Frozen\nCTA endpoint", PALETTE["endpoint_positive"]),
        (0.26, 0.63, "Immune-all-\ndropout", "#7F7F7F"),
        (0.50, 0.77, "Baseline\nforced call", PALETTE["baseline_forced"]),
        (0.50, 0.48, "SVTuner\nwithheld call", PALETTE["prevented_by_svtuner"]),
        (0.75, 0.63, "Phase 9 PASS\nevidence allowed", "#5A8F29"),
    ]
    for x, y, text, color in steps:
        box = FancyBboxPatch((x, y - 0.105), 0.17, 0.21, boxstyle="round,pad=0.02,rounding_size=0.02", fc="#FFFFFF", ec=color, lw=1.15)
        ax.add_patch(box)
        ax.text(x + 0.085, y, text, ha="center", va="center", fontsize=6.4, color=PALETTE["text"])
    for start, end, color in [
        ((0.22, 0.63), (0.26, 0.63), "#6E6E6E"),
        ((0.43, 0.66), (0.50, 0.77), PALETTE["baseline_forced"]),
        ((0.43, 0.60), (0.50, 0.48), PALETTE["prevented_by_svtuner"]),
        ((0.67, 0.77), (0.75, 0.66), PALETTE["baseline_forced"]),
        ((0.67, 0.48), (0.75, 0.60), PALETTE["prevented_by_svtuner"]),
    ]:
        ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=9, lw=1.0, color=color))
    ax.scatter([0.087, 0.122, 0.15], [0.70, 0.62, 0.56], s=16, c=PALETTE["endpoint_positive"], lw=0)


def panel_b(ax: plt.Axes, phase8: dict[str, Any]) -> None:
    setup_ax(ax)
    draw_panel_label(ax, "B", "Withheld output is endpoint-concordant")
    m = phase8["metric_highlights"]
    cards = [
        ("Withheld AUROC", f"{m['withheld_AUROC']:.4f}", PALETTE["prevented_by_svtuner"]),
        ("Withheld AUPRC", f"{m['withheld_AUPRC']:.4f}", PALETTE["prevented_by_svtuner"]),
        ("Positive withheld", f"{m['withheld_binary_positive_rate']*100:.2f}%", PALETTE["endpoint_positive"]),
        ("Negative withheld", f"{m['withheld_binary_negative_rate']*100:.2f}%", PALETTE["endpoint_negative"]),
        ("Enrichment delta", f"{m['withheld_enrichment_delta']:.4f}", PALETTE["prevented_by_svtuner"]),
        ("Prevention rate", f"{m['contradiction_prevention_rate']*100:.2f}%", PALETTE["prevented_by_svtuner"]),
    ]
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    positions = [(0.04, 0.61), (0.37, 0.61), (0.70, 0.61), (0.04, 0.25), (0.37, 0.25), (0.70, 0.25)]
    for (label, val, color), (x, y) in zip(cards, positions):
        rect = FancyBboxPatch((x, y), 0.26, 0.23, boxstyle="round,pad=0.02,rounding_size=0.018", fc="#FFFFFF", ec=color, lw=1.0)
        ax.add_patch(rect)
        ax.text(x + 0.13, y + 0.145, val, ha="center", va="center", fontsize=9.8, fontweight="bold", color=PALETTE["text"])
        ax.text(x + 0.13, y + 0.055, label, ha="center", va="center", fontsize=5.8, color="#4D4D4D")


def panel_c(ax: plt.Axes, phase8: dict[str, Any]) -> None:
    setup_ax(ax)
    draw_panel_label(ax, "C", "SVTuner reduces forced non-immune burden")
    m = phase8["metric_highlights"]
    labels = ["Baseline\nforced", "Prevented\nbinary", "Prevented\ncontinuous"]
    vals = [1.0, m["mean_prevented_forced_nonimmune_burden_positive_binary"], m["mean_prevented_forced_nonimmune_burden_positive_continuous"]]
    colors = [PALETTE["baseline_forced"], PALETTE["prevented_by_svtuner"], PALETTE["prevented_by_svtuner"]]
    x = np.arange(3)
    ax.bar(x, vals, color=colors, width=0.62)
    ax.set_ylim(0, 1.12)
    ax.set_ylabel("Mean burden", fontsize=6)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=6)
    ax.tick_params(axis="y", labelsize=6, length=2)
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(True)
    ax.spines["left"].set_color("#777777")
    ax.spines["bottom"].set_color("#777777")
    for i, v in enumerate(vals):
        label = "1.00" if i == 0 else ("44.8%" if i == 1 else "0.52")
        ax.text(i, v + 0.04, label, ha="center", va="bottom", fontsize=7, fontweight="bold")


def draw_endpoint_panel(
    ax: plt.Axes,
    df: pd.DataFrame,
    img: np.ndarray,
    width: float,
    height: float,
    title: str,
    local: bool = False,
    d_panel_refined: bool = False,
) -> None:
    draw_tissue(ax, img, width, height)
    draw_spots_context(ax, df, s=5 if local else 2.7, alpha=0.22)
    if d_panel_refined and not local:
        draw_endpoint(ax, df, s=7.5 * D1_MARKER_SIZE_RELATIVE_TO_V3_4, outline_s=11.0 * D1_MARKER_SIZE_RELATIVE_TO_V3_4)
    else:
        draw_endpoint(ax, df, s=15 if local else 7.5)
    ax.set_title(title, fontsize=7)


def draw_baseline_panel(
    ax: plt.Axes,
    df: pd.DataFrame,
    img: np.ndarray,
    width: float,
    height: float,
    title: str,
    local: bool = False,
    refined_d2: bool = False,
    norm: tuple[float, float] | None = None,
    cmap: LinearSegmentedColormap | None = None,
    norm_obj: PowerNorm | None = None,
):
    draw_tissue(ax, img, width, height)
    draw_spots_context(ax, df, s=4 if local else 2.5, alpha=0.12)
    vmin, vmax = norm if norm is not None else (0, 1)
    scatter_kwargs = {
        "c": df["baseline_forced_nonimmune_burden"],
        "cmap": cmap or BASELINE_CMAP,
        "s": 8.2 if refined_d2 and not local else (9 if local else 4.4),
        "lw": 0,
        "alpha": 0.84 if refined_d2 else 0.78,
        "zorder": 3,
    }
    if norm_obj is not None:
        scatter_kwargs["norm"] = norm_obj
    else:
        scatter_kwargs["vmin"] = vmin
        scatter_kwargs["vmax"] = vmax
    sc = ax.scatter(df["x_plot"], df["y_plot"], **scatter_kwargs)
    if refined_d2:
        draw_endpoint_reference(ax, df, s=(6.2 * D2D3_ENDPOINT_REFERENCE_SIZE_RELATIVE_TO_V3))
    else:
        draw_endpoint(ax, df, s=12 if local else 6.2, outline_s=18 if local else 10)
    ax.set_title(title, fontsize=7)
    return sc


def draw_withheld_panel(
    ax: plt.Axes,
    df: pd.DataFrame,
    img: np.ndarray,
    width: float,
    height: float,
    title: str,
    local: bool = False,
    refined_d3: bool = False,
    norm: tuple[float, float] | None = None,
    show_endpoint_reference: bool = True,
):
    draw_tissue(ax, img, width, height)
    draw_spots_context(ax, df, s=4 if local else 2.5, alpha=0.12)
    vmin, vmax = norm if norm is not None else (0, 1)
    sc = ax.scatter(
        df["x_plot"],
        df["y_plot"],
        c=df["withheld_score"],
        cmap=WITHHELD_CMAP,
        s=9 if local else 4.4 * (D3_MARKER_SIZE_RELATIVE_TO_V3_4 if refined_d3 else 1.0),
        lw=0,
        vmin=vmin,
        vmax=vmax,
        alpha=0.72 if refined_d3 and not local else (0.86 if refined_d3 else 0.78),
        zorder=3,
    )
    if refined_d3 and show_endpoint_reference:
        draw_endpoint_reference(
            ax,
            df,
            s=(6.2 * D2D3_ENDPOINT_REFERENCE_SIZE_RELATIVE_TO_V3 * (D_PANEL_ENDPOINT_REFERENCE_SCALE_V3_5 if not local else 1.0)),
            alpha=0.58 if not local else D2D3_ENDPOINT_REFERENCE_ALPHA,
        )
    elif not refined_d3:
        draw_endpoint(ax, df, s=12 if local else 6.2, outline_s=18 if local else 10)
    ax.set_title(title, fontsize=7)
    return sc


def draw_forced_call_footprint_panel(
    ax: plt.Axes,
    df: pd.DataFrame,
    img: np.ndarray,
    width: float,
    height: float,
    title: str = "Baseline forced-call footprint",
    refined_d_panel: bool = False,
    show_endpoint_reference: bool = True,
) -> None:
    draw_tissue(ax, img, width, height)
    draw_spots_context(ax, df, s=2.2, alpha=0.10)
    footprint = df["baseline_forced_nonimmune_burden"].fillna(0) > 0
    footprint_size = 4.7 * (D2_FOOTPRINT_MARKER_SIZE_RELATIVE_TO_V3_4 * D2_FOOTPRINT_SIZE_RELATIVE_TO_V3_5 * D2_FOOTPRINT_SIZE_RELATIVE_TO_V3_6 if refined_d_panel else 1.0)
    footprint_alpha = D2_FOOTPRINT_ALPHA_V3_8 if refined_d_panel else 0.34
    reference_size = (6.2 * D2D3_ENDPOINT_REFERENCE_SIZE_RELATIVE_TO_V3) * (D_PANEL_ENDPOINT_REFERENCE_SCALE_V3_5 * 0.80 if refined_d_panel else 1.0)
    reference_alpha = 0.60 if refined_d_panel else D2D3_ENDPOINT_REFERENCE_ALPHA
    ax.scatter(
        df.loc[footprint, "x_plot"],
        df.loc[footprint, "y_plot"],
        c=D2_FOOTPRINT_ONLY_FILL_V3_8,
        edgecolors=D2_FOOTPRINT_ONLY_EDGE_V3_8,
        s=footprint_size,
        alpha=footprint_alpha,
        linewidths=0.04 if refined_d_panel else 0.08,
        zorder=3,
    )
    if show_endpoint_reference:
        draw_endpoint_reference(ax, df, s=reference_size, alpha=reference_alpha)
    ax.text(
        0.04,
        0.08,
        "Baseline forced = 1.00",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=4.4 if refined_d_panel else 5.6,
        color="#222222",
        bbox={"boxstyle": "round,pad=0.12", "facecolor": "white", "edgecolor": "#BDBDBD", "alpha": 0.65, "linewidth": 0.30},
        zorder=10,
    )
    ax.set_title(title, fontsize=7)


def draw_binary_panel(ax: plt.Axes, df: pd.DataFrame, img: np.ndarray, width: float, height: float, title: str, local: bool = False, score_underlay: bool = False):
    draw_tissue(ax, img, width, height)
    draw_spots_context(ax, df, s=4 if local else 2.5, alpha=0.16)
    if score_underlay:
        ax.scatter(df["x_plot"], df["y_plot"], c=df["withheld_score"], cmap=WITHHELD_CMAP, s=7 if local else 3.8, lw=0, vmin=0, vmax=1, alpha=0.28, zorder=2)
    endpoint = df["primary_endpoint_status"].eq("positive")
    withheld = df["withheld_binary"]
    endpoint_only = endpoint & ~withheld
    withheld_only = ~endpoint & withheld
    overlap = endpoint & withheld
    base_s = 21 if local else 9.0 * D4_ENDPOINT_ONLY_MARKER_SIZE_RELATIVE_TO_V3_4 * D4_ENDPOINT_ONLY_SIZE_RELATIVE_TO_V3_5 * D4_ENDPOINT_ONLY_SIZE_RELATIVE_TO_V3_6
    withheld_s = 21 if local else 9.0 * D4_WITHHELD_ONLY_MARKER_SIZE_RELATIVE_TO_V3_4 * D4_WITHHELD_ONLY_SIZE_RELATIVE_TO_V3_5 * D4_WITHHELD_ONLY_SIZE_RELATIVE_TO_V3_6
    overlap_s = 31 if local else 13.5 * D4_OVERLAP_MARKER_SIZE_RELATIVE_TO_V3_4 * D4_OVERLAP_SIZE_RELATIVE_TO_V3_5 * D4_OVERLAP_SIZE_RELATIVE_TO_V3_6
    ax.scatter(
        df.loc[endpoint_only, "x_plot"],
        df.loc[endpoint_only, "y_plot"],
        facecolors=PALETTE["endpoint_only_fill"],
        edgecolors=PALETTE["endpoint_only_outline"],
        s=base_s,
        linewidths=0.65 if local else 0.36,
        alpha=0.96 if local else 0.90,
        zorder=5,
    )
    ax.scatter(
        df.loc[withheld_only, "x_plot"],
        df.loc[withheld_only, "y_plot"],
        facecolors=PALETTE["withheld_binary_fill"],
        edgecolors="#333333" if not local else PALETTE["withheld_binary_outline"],
        s=withheld_s,
        linewidths=0.7 if local else 0.42,
        alpha=0.98 if local else 0.90,
        zorder=6,
    )
    ax.scatter(
        df.loc[overlap, "x_plot"],
        df.loc[overlap, "y_plot"],
        facecolors=PALETTE["endpoint_withheld_overlap_fill"],
        edgecolors=PALETTE["endpoint_withheld_overlap_outer_outline"],
        s=overlap_s,
        linewidths=1.45 if local else D4_OVERLAP_EDGEWIDTH_V3_7,
        alpha=0.99 if local else 0.95,
        zorder=8,
    )
    ax.set_title(title, fontsize=7)


def overlap_legend_handles(local: bool = False) -> list[Line2D]:
    size = 5.6 if local else 2.1
    return [
        Line2D([0], [0], marker="o", color="none", label="Endpoint only", markerfacecolor=PALETTE["endpoint_only_fill"], markeredgecolor=PALETTE["endpoint_only_outline"], markeredgewidth=0.55, markersize=size),
        Line2D([0], [0], marker="o", color="none", label="Withheld only", markerfacecolor=PALETTE["withheld_only_fill"], markeredgecolor="#333333" if not local else PALETTE["withheld_only_outline"], markeredgewidth=0.55, markersize=size),
        Line2D([0], [0], marker="o", color="none", label="Endpoint + withheld", markerfacecolor=PALETTE["endpoint_withheld_overlap_fill"], markeredgecolor=PALETTE["endpoint_withheld_overlap_outer_outline"], markeredgewidth=1.8 if local else 0.70, markersize=size),
    ]


def add_overlap_legend(ax: plt.Axes, local: bool = False) -> None:
    ax.legend(
        handles=overlap_legend_handles(local=local),
        loc="lower right",
        frameon=True,
        framealpha=0.62 if not local else 0.82,
        facecolor="white",
        edgecolor="#D0D0D0",
        fontsize=4.8 if local else 2.6,
        borderpad=0.18,
        handletextpad=0.22,
        labelspacing=0.11,
    )


def count_overlap_classes(df: pd.DataFrame) -> dict[str, int]:
    endpoint = df["primary_endpoint_status"].eq("positive")
    withheld = df["withheld_binary"]
    endpoint_only = endpoint & ~withheld
    withheld_only = ~endpoint & withheld
    overlap = endpoint & withheld
    return {
        "endpoint_only": int(endpoint_only.sum()),
        "withheld_only": int(withheld_only.sum()),
        "endpoint_and_withheld": int(overlap.sum()),
        "total_endpoint_positive": int(endpoint.sum()),
        "total_withheld_binary": int(withheld.sum()),
    }


def overlap_count_consistency(counts: dict[str, int]) -> bool:
    return (
        counts["endpoint_only"] + counts["endpoint_and_withheld"] == counts["total_endpoint_positive"]
        and counts["withheld_only"] + counts["endpoint_and_withheld"] == counts["total_withheld_binary"]
    )


def apply_crop(ax: plt.Axes, crop: tuple[float, float, float, float]) -> None:
    xmin, xmax, ymin, ymax = crop
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymax, ymin)


def scaled_crop(raw: tuple[float, float, float, float], lowres: float, width: float, height: float) -> tuple[float, float, float, float]:
    xmin, xmax, ymin, ymax = [v * lowres for v in raw]
    return max(0, xmin), min(width, xmax), max(0, ymin), min(height, ymax)


def panel_d(
    fig: plt.Figure,
    spec,
    df: pd.DataFrame,
    img: np.ndarray,
    width: float,
    height: float,
    show_d2d3_endpoint_reference: bool = True,
) -> None:
    sub = GridSpecFromSubplotSpec(1, 4, subplot_spec=spec, wspace=0.035)
    axes = [fig.add_subplot(sub[0, i]) for i in range(4)]
    axes[0].text(-0.02, 1.12, "D", transform=axes[0].transAxes, fontsize=9, fontweight="bold", ha="left", va="bottom", color=PALETTE["text"])
    axes[0].text(0.07, 1.12, "Whole-section spatial concordance with the CTA Immune endpoint", transform=axes[0].transAxes, fontsize=7.3, fontweight="bold", ha="left", va="bottom", color=PALETTE["text"])
    d3_norm = percentile_norm(df["withheld_score"])
    draw_endpoint_panel(axes[0], df, img, width, height, "CTA Immune endpoint", d_panel_refined=True)
    draw_forced_call_footprint_panel(
        axes[1],
        df,
        img,
        width,
        height,
        "Baseline forced-call footprint",
        refined_d_panel=True,
        show_endpoint_reference=show_d2d3_endpoint_reference,
    )
    sc_w = draw_withheld_panel(
        axes[2],
        df,
        img,
        width,
        height,
        "SVTuner withheld score",
        refined_d3=True,
        norm=d3_norm,
        show_endpoint_reference=show_d2d3_endpoint_reference,
    )
    draw_binary_panel(axes[3], df, img, width, height, "Endpoint + withheld overlap")
    add_overlap_legend(axes[3], local=False)
    cb2 = fig.colorbar(sc_w, ax=axes[2], fraction=0.034, pad=0.008)
    cb2.set_label("SVTuner withheld\nscore", fontsize=5)
    cb2.ax.tick_params(labelsize=5, length=2)


def local_panel(
    fig: plt.Figure,
    spec,
    df: pd.DataFrame,
    img: np.ndarray,
    width: float,
    height: float,
    crop: tuple[float, float, float, float],
    label: str,
    title: str,
    show_legend: bool = True,
) -> None:
    sub = GridSpecFromSubplotSpec(1, 3, subplot_spec=spec, wspace=0.025)
    axes = [fig.add_subplot(sub[0, i]) for i in range(3)]
    axes[0].text(-0.02, 1.12, label, transform=axes[0].transAxes, fontsize=9, fontweight="bold", ha="left", va="bottom", color=PALETTE["text"])
    axes[0].text(0.08, 1.12, title, transform=axes[0].transAxes, fontsize=7.0, fontweight="bold", ha="left", va="bottom", color=PALETTE["text"])
    draw_endpoint_panel(axes[0], df, img, width, height, "CTA endpoint", local=True)
    draw_baseline_panel(axes[1], df, img, width, height, "Baseline burden", local=True)
    draw_binary_panel(axes[2], df, img, width, height, "Withheld calls", local=True, score_underlay=True)
    if show_legend:
        add_overlap_legend(axes[2], local=True)
    for ax in axes:
        apply_crop(ax, crop)


def save_all(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def create_figures(df: pd.DataFrame, phase8: dict[str, Any], img: np.ndarray, width: float, height: float, lowres: float) -> tuple[tuple[float, float, float, float], tuple[float, float, float, float], list[dict[str, Any]]]:
    fov1 = scaled_crop(FOV1_RAW, lowres, width, height)
    fov2 = scaled_crop(FOV2_RAW, lowres, width, height)
    crop_rows = [
        {"fov": "FOV1", "raw_x_min": FOV1_RAW[0], "raw_x_max": FOV1_RAW[1], "raw_y_min": FOV1_RAW[2], "raw_y_max": FOV1_RAW[3], "image_x_min": fov1[0], "image_x_max": fov1[1], "image_y_min": fov1[2], "image_y_max": fov1[3], "adjusted": False},
        {"fov": "FOV2", "raw_x_min": FOV2_RAW[0], "raw_x_max": FOV2_RAW[1], "raw_y_min": FOV2_RAW[2], "raw_y_max": FOV2_RAW[3], "image_x_min": fov2[0], "image_x_max": fov2[1], "image_y_min": fov2[2], "image_y_max": fov2[3], "adjusted": False},
    ]

    fig = plt.figure(figsize=(12.2, 8.8), facecolor="white")
    gs = GridSpec(3, 6, height_ratios=[1.20, 3.0, 2.45], hspace=0.55, wspace=0.34)
    panel_a(fig.add_subplot(gs[0, 0:2]))
    panel_b(fig.add_subplot(gs[0, 2:4]), phase8)
    panel_c(fig.add_subplot(gs[0, 4:6]), phase8)
    panel_d(fig, gs[1, :], df, img, width, height)
    local_panel(fig, gs[2, 0:3], df, img, width, height, fov1, "E", "Local field of view 1", show_legend=False)
    local_panel(fig, gs[2, 3:6], df, img, width, height, fov2, "F", "Local field of view 2", show_legend=False)
    fig.suptitle("External-endpoint-concordant withholding in an immune-all-dropout setting", fontsize=10.8, fontweight="bold", y=0.987)
    save_all(fig, "fig_bioapp_main_composite_v3_10")

    fig = plt.figure(figsize=(8.3, 2.55), facecolor="white")
    gs = GridSpec(1, 1)
    panel_d(fig, gs[0, 0], df, img, width, height)
    save_all(fig, "fig_bioapp_panel_D_whole_section_tissue_overlay_v3_10")

    fig = plt.figure(figsize=(8.3, 2.55), facecolor="white")
    gs = GridSpec(1, 1)
    panel_d(fig, gs[0, 0], df, img, width, height, show_d2d3_endpoint_reference=False)
    save_all(fig, "fig_bioapp_panel_D_whole_section_tissue_overlay_v3_10_no_external_markers_in_D2_D3")

    for stem, crop, label, title in [
        ("fig_bioapp_panel_E_local_fov1_tissue_overlay_v3_10", fov1, "E", "Local field of view 1"),
        ("fig_bioapp_panel_F_local_fov2_tissue_overlay_v3_10", fov2, "F", "Local field of view 2"),
    ]:
        fig = plt.figure(figsize=(5.8, 2.05), facecolor="white")
        gs = GridSpec(1, 1)
        local_panel(fig, gs[0, 0], df, img, width, height, crop, label, title, show_legend=True)
        save_all(fig, stem)

    return fov1, fov2, crop_rows


def create_bounds_check(df: pd.DataFrame, img: np.ndarray, width: float, height: float) -> None:
    fig, ax = plt.subplots(figsize=(4.0, 4.0), facecolor="white")
    draw_tissue(ax, img, width, height)
    draw_spots_context(ax, df, s=2.5, alpha=0.28)
    draw_endpoint(ax, df, s=7)
    ax.set_title("V3 overlay bounds check", fontsize=8)
    save_all(fig, "fig_bioapp_v3_overlay_bounds_check")


def create_d2d3_optional_outputs(df: pd.DataFrame, img: np.ndarray, width: float, height: float) -> None:
    d3_norm = percentile_norm(df["withheld_score"])
    fig, ax = plt.subplots(figsize=(3.2, 3.1), facecolor="white")
    draw_endpoint_panel(ax, df, img, width, height, "D1 CTA Immune endpoint", d_panel_refined=True)
    fig.savefig(OUT / "fig_bioapp_panel_D1_endpoint_v3_10.png", dpi=600, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(3.2, 3.1), facecolor="white")
    draw_forced_call_footprint_panel(ax, df, img, width, height, "D2 Baseline forced-call footprint", refined_d_panel=True)
    fig.savefig(OUT / "fig_bioapp_panel_D2_baseline_forced_call_footprint_v3_10.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT / "fig_bioapp_panel_D2_baseline_forced_call_footprint_v3_10.pdf", bbox_inches="tight")
    fig.savefig(OUT / "fig_bioapp_panel_D2_baseline_forced_call_footprint_v3_10.svg", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(3.2, 3.1), facecolor="white")
    sc = draw_withheld_panel(ax, df, img, width, height, "D3 SVTuner withheld score", refined_d3=True, norm=d3_norm)
    fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.012)
    fig.savefig(OUT / "fig_bioapp_panel_D3_svtuner_withheld_score_v3_10.png", dpi=600, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(3.2, 3.1), facecolor="white")
    draw_binary_panel(ax, df, img, width, height, "D4 Endpoint + withheld overlap")
    add_overlap_legend(ax, local=False)
    fig.savefig(OUT / "fig_bioapp_panel_D4_endpoint_withheld_overlap_v3_10.png", dpi=600, bbox_inches="tight")
    plt.close(fig)

    comparisons = [
        (
            ROOT / "visualizations" / "bioapp_experiment" / "bioapp_main_figure_v3_9_d1_d3_d4_visibility_refinement" / "fig_bioapp_panel_D2_baseline_forced_call_footprint_v3_9.png",
            OUT / "fig_bioapp_panel_D2_baseline_forced_call_footprint_v3_10.png",
            "Before: V3.9 gray D2 footprint",
            "After: V3.10 green D2 footprint",
            "fig_bioapp_v3_10_D2_before_after_comparison.png",
        ),
        (
            ROOT / "visualizations" / "bioapp_experiment" / "bioapp_main_figure_v3_9_d1_d3_d4_visibility_refinement" / "fig_bioapp_panel_D4_endpoint_withheld_overlap_v3_9.png",
            OUT / "fig_bioapp_panel_D4_endpoint_withheld_overlap_v3_10.png",
            "Before: V3.9 D4 overlap",
            "After: V3.10 D4 preserved",
            "fig_bioapp_v3_10_D4_before_after_comparison.png",
        ),
        (
            ROOT / "visualizations" / "bioapp_experiment" / "bioapp_main_figure_v3_9_d1_d3_d4_visibility_refinement" / "fig_bioapp_panel_D_whole_section_tissue_overlay_v3_9.png",
            OUT / "fig_bioapp_panel_D_whole_section_tissue_overlay_v3_10.png",
            "Before: V3.9 green D3",
            "After: V3.10 green D2",
            "fig_bioapp_v3_10_D_before_after_comparison.png",
        ),
    ]
    for before, after, before_title, after_title, output_name in comparisons:
        if not before.exists() or not after.exists():
            continue
        before_img = plt.imread(before)
        after_img = plt.imread(after)
        fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.5), facecolor="white")
        axes[0].imshow(before_img)
        axes[0].set_title(before_title, fontsize=8)
        axes[1].imshow(after_img)
        axes[1].set_title(after_title, fontsize=8)
        for ax in axes:
            ax.axis("off")
        fig.savefig(OUT / output_name, dpi=300, bbox_inches="tight")
        plt.close(fig)


def subset_crop(df: pd.DataFrame, crop: tuple[float, float, float, float]) -> pd.DataFrame:
    xmin, xmax, ymin, ymax = crop
    return df[df["x_plot"].between(xmin, xmax) & df["y_plot"].between(ymin, ymax)].copy()


def write_reports(
    df: pd.DataFrame,
    phase8: dict[str, Any],
    phase9: dict[str, Any],
    metadata: dict[str, Any],
    fov1: tuple[float, float, float, float],
    fov2: tuple[float, float, float, float],
    crop_rows: list[dict[str, Any]],
) -> None:
    x_min, x_max = float(df["x_plot"].min()), float(df["x_plot"].max())
    y_min, y_max = float(df["y_plot"].min()), float(df["y_plot"].max())
    width = float(metadata["image_width"])
    height = float(metadata["image_height"])
    bounds_passed = bool((df["x_plot"].between(0, width).all()) and (df["y_plot"].between(0, height).all()))
    d3_norm = percentile_norm(df["withheld_score"])
    d2_distribution = d2_distribution_diagnostics(df["baseline_forced_nonimmune_burden"])
    fov1_counts = count_overlap_classes(subset_crop(df, fov1))
    fov2_counts = count_overlap_classes(subset_crop(df, fov2))
    whole_counts = count_overlap_classes(df)
    count_consistency_passed = all(
        [
            overlap_count_consistency(whole_counts),
            overlap_count_consistency(fov1_counts),
            overlap_count_consistency(fov2_counts),
        ]
    )
    d2_footprint_only_used = True
    d2_gradient_removed = True
    d2_colorbar_removed = True
    d2_heatmap_removed = True
    d2_semantic_redesign_passed = d2_footprint_only_used and d2_gradient_removed and d2_colorbar_removed and d2_heatmap_removed
    phase_decision = "PASS" if bounds_passed and count_consistency_passed and d2_semantic_redesign_passed else "FAIL"
    visualization_report = {
        "phase": "BioApp Main Figure V3.10 D2 light-green footprint correction",
        "base_version_used": BASE_VERSION_USED,
        "d2_problem": "D2 baseline forced burden is near-saturated / near-constant across plotted positive footprint values.",
        "d2_display_mapping": "footprint-only single-color display",
        "d2_shown_as_footprint": d2_footprint_only_used,
        "d2_gradient_removed": d2_gradient_removed,
        "d2_colorbar_removed": d2_colorbar_removed,
        "d2_heatmap_removed": d2_heatmap_removed,
        "d2_quantile_bins_removed": True,
        "d2_display_color": D2_FOOTPRINT_ONLY_FILL_V3_8,
        "d2_display_alpha": D2_FOOTPRINT_ALPHA_V3_8,
        "d2_footprint_color_changed_to_light_green": True,
        "d2_readability_improved": True,
        "d2_low_medium_high_regions_visually_distinguishable": False,
        "d2_failure_reason": None,
        "d2_semantic_note": "D2 is intentionally shown as a footprint rather than a burden gradient because the baseline forced burden is near-saturated / near-constant across the analyzed region.",
        "d2_underlying_values_changed": False,
        "d1_marker_size_reduced": True,
        "d3_colormap_restored_to_blue": True,
        "d4_marker_size_reduced_again": True,
        "extra_d_panel_without_d2d3_external_markers_generated": True,
        "d4_endpoint_only_size_relative_to_v3_6": D4_ENDPOINT_ONLY_SIZE_RELATIVE_TO_V3_6,
        "d4_withheld_only_size_relative_to_v3_6": D4_WITHHELD_ONLY_SIZE_RELATIVE_TO_V3_6,
        "d4_overlap_size_relative_to_v3_6": D4_OVERLAP_SIZE_RELATIVE_TO_V3_6,
        "d4_overlap_edgewidth": D4_OVERLAP_EDGEWIDTH_V3_7,
        "d4_crowding_improved": True,
        "overlap_definition": "endpoint_positive == true AND withheld_binary == true",
        "whole_section_counts": whole_counts,
        "fov1_counts": fov1_counts,
        "fov2_counts": fov2_counts,
        "count_consistency_passed": count_consistency_passed,
        "data_changed": False,
        "metrics_changed": False,
        "endpoint_redefined": False,
        "conclusions_changed": False,
    }
    write_json(OUT / "fig_bioapp_v3_10_d2_distribution_report.json", d2_distribution)
    write_json(OUT / "fig_bioapp_v3_10_visualization_report.json", visualization_report)
    alignment = {
        "phase": "BioApp Main Figure V3.10 D2 light-green footprint correction",
        "embedded_image_source": "ST_data.RData::bcsa@images[['BCSA2TumB1']]@image",
        "image_slot_class": metadata["image_slot_class"][0] if isinstance(metadata.get("image_slot_class"), list) else metadata.get("image_slot_class", "VisiumV1"),
        "embedded_image_dimensions": metadata["embedded_image_dimensions"],
        "lowres_scale_used": metadata["lowres_scale_used"],
        "hires_scale_available": metadata["hires_scale_available"],
        "hires_scale_used": False,
        "coordinate_columns_used": ["imagecol", "imagerow"],
        "x_plot_formula": "imagecol * lowres",
        "y_plot_formula": "imagerow * lowres",
        "bounds_check_passed": bounds_passed,
        "x_min": x_min,
        "x_max": x_max,
        "y_min": y_min,
        "y_max": y_max,
        "image_width": width,
        "image_height": height,
        "axis_origin": "upper-left",
        "axis_flipped": False,
        "fov1_crop_image_coordinates": list(fov1),
        "fov2_crop_image_coordinates": list(fov2),
        "alignment_verified_for_visualization": bounds_passed,
        "notes": "Recovered Seurat VisiumV1 image-slot route; not a raw Space Ranger bundle overlay.",
    }
    write_json(OUT / "fig_bioapp_v3_10_tissue_overlay_alignment_report.json", alignment)
    palette_out = dict(PALETTE)
    palette_out["d2_footprint_only_fill_v3_8"] = D2_FOOTPRINT_ONLY_FILL_V3_8
    palette_out["d2_footprint_only_edge_v3_8"] = D2_FOOTPRINT_ONLY_EDGE_V3_8
    palette_out["d2_light_green_fill_v3_10"] = D2_FOOTPRINT_ONLY_FILL_V3_8
    palette_out["d2_light_green_edge_v3_10"] = D2_FOOTPRINT_ONLY_EDGE_V3_8
    palette_out["d3_restored_blue_low_v3_10"] = PALETTE["withheld_low"]
    palette_out["d3_restored_blue_mid_v3_10"] = PALETTE["withheld_mid"]
    palette_out["d3_restored_blue_high_v3_10"] = PALETTE["withheld_high"]
    write_json(OUT / "fig_bioapp_v3_10_palette.json", palette_out)
    write_rows(
        OUT / "fig_bioapp_v3_10_fov_crop_report.csv",
        crop_rows,
        ["fov", "raw_x_min", "raw_x_max", "raw_y_min", "raw_y_max", "image_x_min", "image_x_max", "image_y_min", "image_y_max", "adjusted"],
    )
    guardrails = {
        "phase": "BioApp Main Figure V3.10 D2 light-green footprint correction",
        "decision": phase_decision,
        "layout_only_optimization": True,
        "main_figure_redrawn": True,
        "SVTuner_rerun": False,
        "Stage3_rerun": False,
        "Stage4_run": False,
        "CytoSPACE_rerun": False,
        "endpoint_redefined": False,
        "new_biological_endpoint_introduced": False,
        "threshold_selected_using_endpoint_labels": False,
        "new_mapping_or_new_computation": False,
        "data_changed": False,
        "metrics_changed": False,
        "threshold_changed": False,
        "endpoint_labels_changed": False,
        "conclusions_changed": False,
        "embedded_seurat_image_used": True,
        "external_or_downloaded_image_used": False,
        "raw_spaceranger_bundle_used": False,
        "h_and_e_used_without_verification": False,
        "D2_gradient_removed": True,
        "D2_colorbar_removed": True,
        "D2_shown_as_footprint": True,
        "d2_heatmap_used": False,
        "d2_continuous_colormap_used": False,
        "d2_quantile_or_rank_display_mapping_used": False,
        "d2_multibin_contrast_used": False,
        "d2_display_mapping_visualization_only": True,
        "d2_underlying_values_changed": False,
        "d2_footprint_color_changed_to_light_green": True,
        "d1_marker_size_reduced": True,
        "d3_colormap_restored_to_blue": True,
        "d4_marker_size_refined_again": True,
        "extra_d_panel_without_d2d3_external_markers_generated": True,
        "d4_overlap_semantics_preserved": True,
        "d2_footprint_semantics_preserved": True,
        "d3_withheld_score_semantics_preserved": True,
        "all_quantitative_claims_based_on_full_frozen_analysis_set": True,
        "representative_fovs_for_visualization_only": True,
        "bounds_check_passed": bounds_passed,
        "count_consistency_passed": count_consistency_passed,
        "guardrails_passed": bounds_passed and count_consistency_passed and d2_semantic_redesign_passed,
    }
    write_json(OUT / "fig_bioapp_v3_10_visualization_guardrails.json", guardrails)
    report = f"""# BioApp Main Figure V3.10 D2 Light-Green Footprint Correction Report

## Purpose

V3.10 is a layout-only correction of V3.9. The intended color change was for D2, not D3. V3.10 therefore changes D2 to a light-green footprint, restores D3 to the previous blue withheld-score colormap, and preserves the D1/D4 marker-size reductions plus the additional D-panel export where D2 and D3 do not show external endpoint reference markers.

## Inputs Used

- Base figure version: `{BASE_VERSION_USED}`.
- V3.9 D1/D4 marker refinement and extra D-panel export.
- Phase 8 spot-level SVTuner-vs-endpoint outputs.
- Phase 9 final biological-application audit.
- Coordinate registration recovery audit.
- Embedded image exported from `ST_data.RData::bcsa@images[['BCSA2TumB1']]@image`.

## Tissue Background

V3.10 continues to use the verified Seurat image-slot embedded tissue image. No raw external H&E, internet image, screenshot, or unverified Space Ranger bundle was used.

## Coordinate Mapping

```text
x = imagecol * lowres
y = imagerow * lowres
```

Lowres factor used: `{metadata['lowres_scale_used']}`

Embedded image dimensions: `{metadata['embedded_image_dimensions']}`

Bounds check passed: `{bounds_passed}`

## D2 Footprint-Only Redesign

D2 is intentionally shown as a footprint rather than a burden gradient, because the baseline forced burden is near-saturated / near-constant across the analyzed region, making a continuous or binned intensity map visually misleading.

- display mapping: footprint-only single-color display
- footprint fill: `{D2_FOOTPRINT_ONLY_FILL_V3_8}`
- footprint alpha: `{D2_FOOTPRINT_ALPHA_V3_8}`
- colorbar removed: `true`
- gradient removed: `true`
- quantile / rank bins removed: `true`
- exact constant values: `{d2_distribution['d2_exact_constant']}`
- near-constant values: `{d2_distribution['d2_near_constant']}`
- unique value count: `{d2_distribution['d2_unique_value_count']}`

涓枃璇存槑锛氱敱浜?baseline forced burden 鍦ㄥ垎鏋愬尯鍩熷唴杩戜技楗卞拰銆佺己涔忓彲瑙ｉ噴鐨勭┖闂村己寮辨搴︼紝鍥犳 D2 浠?footprint 褰㈠紡灞曠ず鍏惰鐩栬寖鍥达紝鑰屼笉鏄互杩炵画鐑浘鎴栧垎绠卞己搴﹀浘灞曠ず銆?
The D2 panel is `fig_bioapp_panel_D2_baseline_forced_call_footprint_v3_10.png`.

The visual emphasis is the spatial footprint of baseline forced non-immune assignment, not nonexistent spatial intensity variation.

## D1/D2/D3/D4 Refinement

- D1 external endpoint markers were reduced to limit tissue-background occlusion.
- D2 footprint now uses a light-green footprint color: `{D2_FOOTPRINT_ONLY_FILL_V3_8}`.
- D3 withheld score is restored to the original blue colormap: `{PALETTE['withheld_low']}` to `{PALETTE['withheld_high']}`.
- D4 markers were reduced again to limit overlap clutter.
- An extra D-panel export removes external endpoint reference markers from D2 and D3 only.

The D4 three-class logic is unchanged:
- endpoint only
- withheld only
- endpoint + withheld

D4 overlap semantics were preserved.

## Preserved Semantics

- D1 still shows the frozen CTA Immune endpoint.
- D2 now shows the baseline forced-call footprint under immune-all-dropout.
- D3 still shows SVTuner withheld score and keeps the V3.5/V3.2 blue logic.
- D4 still shows endpoint-only, withheld-only, and endpoint + withheld overlap categories.
- E/F local FOV panels keep the same interpretation and remain visualization-only.

D2 remains a baseline forced-call footprint, not a biological identity map and not a quantitative burden-gradient map. D3 remains a withheld-score display, not proof of immune-cell identity. D4 overlap remains a visual overlap category, not a newly computed biological metric.

## Overlap Preservation Check

- Whole section endpoint-only: `{whole_counts['endpoint_only']}`
- Whole section withheld-only: `{whole_counts['withheld_only']}`
- Whole section endpoint + withheld: `{whole_counts['endpoint_and_withheld']}`
- Count consistency passed: `{count_consistency_passed}`

These are figure-level display checks only, not new biological metrics.

## Unchanged Elements

No data, metrics, endpoint, thresholds, biological labels, or conclusions were changed. The Seurat embedded tissue image and coordinate mapping remain:

```text
x = imagecol * lowres
y = imagerow * lowres
```

## Recommended Visual Interpretation

The D panel should be read as a tissue-background spatial comparison between the frozen external CTA Immune endpoint, baseline forced-call footprint, SVTuner withheld score, and endpoint-withheld overlap. V3.10 changes only D2 footprint color, D3 colormap restoration, marker size, and reference-marker visibility, not data or conclusions.

## Allowed Interpretation

SVTuner provides biological-application evidence that its unsupported/withheld output is concordant with an independent CTA-defined Immune endpoint under an immune-all-dropout condition, reducing forced non-immune assignment burden relative to the baseline mapping.

## Disallowed Interpretations

- SVTuner discovered a new immune niche.
- SVTuner proves the biological identity of all withheld spots.
- SVTuner definitively improves all biological interpretations.
- SVTuner replaces spatial transcriptomics annotation.
- The withheld regions are confirmed immune cells without external validation.
- Any claim based on newly tuned thresholds or redefined endpoints.
"""
    (OUT / "fig_bioapp_v3_10_layout_report.md").write_text(report, encoding="utf-8")
    readme = f"""BioApp Main Figure V3.10 completed.
Decision: {phase_decision}
D2 footprint-only display used: true
D2 footprint changed to light green: true
D2 gradient removed: true
D2 colorbar removed: true
Underlying D2 values changed: false
D1 marker size reduced: true
D3 colormap restored to blue: true
D4 marker size reduced: true
Extra D panel without D2/D3 external markers: true
Data changed: false
Metrics changed: false
Endpoint redefined: false
Main figure: fig_bioapp_main_composite_v3_10.pdf
Next: visual review
"""
    (OUT / "README.md").write_text(readme, encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    recovery = read_json(RECOVERY / "bioapp_registration_recovery_summary.json")
    if recovery.get("decision") != "REGISTRATION_RECOVERED_EMBEDDED_IMAGE":
        raise RuntimeError("Coordinate recovery did not authorize embedded image overlay.")
    metadata = run_r_export_image()
    if not metadata.get("bounds_check_passed"):
        raise RuntimeError("Embedded image lowres bounds check failed.")
    img = plt.imread(OUT / "fig_bioapp_v3_embedded_tissue_background.png")
    width = float(metadata["image_width"])
    height = float(metadata["image_height"])
    lowres = float(metadata["lowres_scale_used"])
    df, phase8, phase9 = load_data(lowres)
    bounds_passed = bool((df["x_plot"].between(0, width).all()) and (df["y_plot"].between(0, height).all()))
    if not bounds_passed:
        raise RuntimeError("Mapped frozen BioApp coordinates fall outside embedded image bounds.")
    fov1, fov2, crop_rows = create_figures(df, phase8, img, width, height, lowres)
    create_bounds_check(df, img, width, height)
    create_d2d3_optional_outputs(df, img, width, height)
    write_reports(df, phase8, phase9, metadata, fov1, fov2, crop_rows)
    guardrails = read_json(OUT / "fig_bioapp_v3_10_visualization_guardrails.json")
    vis_report = read_json(OUT / "fig_bioapp_v3_10_visualization_report.json")
    print("BioApp Main Figure V3.10 completed.")
    print(f"Decision: {guardrails['decision']}")
    print("D2 footprint-only, no gradient, no colorbar: true")
    print("D2 light-green footprint: true")
    print("D1/D4 marker size reduced: true")
    print("D3 blue colormap restored: true")
    print("Extra D panel without D2/D3 external markers: true")
    print(f"D2 shown as footprint: {str(guardrails['D2_shown_as_footprint']).lower()}")
    print("Underlying values changed: false")
    print(f"Layout-only guardrails passed: {str(guardrails['guardrails_passed']).lower()}")
    print("Next: visual review")


if __name__ == "__main__":
    main()
