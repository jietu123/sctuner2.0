from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.build_fig2e_stage3_profile_mask_benchmark import (
    CELL_TYPES,
    OUT_DIR,
    RESULT_DIR,
    SCENARIOS,
    _nes_color,
    _p_size,
)


METHOD = "SVTuner + CytoSPACE"
FEATURES = ["CE9", "CE10"]


def plot_route2_only(metrics: pd.DataFrame) -> Path:
    mpl.rcParams["svg.fonttype"] = "none"
    mpl.rcParams["font.family"] = "Arial"
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(4.85, 2.55), dpi=300)
    fig.patch.set_facecolor("white")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(
        0.50,
        0.965,
        "CE9/CE10 spatial enrichment with SVTuner + CytoSPACE",
        fontsize=6.8,
        weight="bold",
        ha="center",
        va="top",
    )
    ax.text(0.075, 0.800, "ST dataset", fontsize=4.3, ha="left", va="bottom")
    ax.text(0.215, 0.800, "scRNA-seq", fontsize=4.3, ha="left", va="bottom")

    x_dataset = 0.075
    x_scrna = 0.215
    x_platform_bar = 0.330
    block_x0 = {"CE9": 0.395, "CE10": 0.565}
    dx = 0.029
    y_top = 0.728
    dy_scenario = 0.052
    cell_w = dx
    cell_h = dy_scenario

    for feature, x0 in block_x0.items():
        ax.text(x0 + 2 * dx, 0.835, feature, fontsize=5.0, ha="center", va="bottom")
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
        ax.text(x_scrna, y_mid, scenario.scrna_label, fontsize=4.1, ha="left", va="center")
        ax.add_patch(
            plt.Rectangle(
                (x_platform_bar - 0.006, y_mid - cell_h * 0.34),
                0.012,
                cell_h * 0.68,
                facecolor=scenario.platform_color,
                edgecolor="none",
            )
        )

        for feature, x0 in block_x0.items():
            for cidx, (canonical, _, _, _, _) in enumerate(CELL_TYPES):
                row = metrics[
                    (metrics["scenario"] == scenario.sample)
                    & (metrics["method"] == METHOD)
                    & (metrics["cell_type"] == canonical)
                    & (metrics["feature"] == feature)
                ].iloc[0]
                x = x0 + cidx * dx
                nes = float(row["NES"]) if pd.notna(row["NES"]) else float("nan")
                pval = float(row["P-value"]) if pd.notna(row["P-value"]) else float("nan")
                if not pd.notna(nes):
                    ax.add_patch(
                        plt.Rectangle(
                            (x - cell_w / 2 + 0.002, y_mid - cell_h / 2 + 0.002),
                            cell_w - 0.004,
                            cell_h - 0.004,
                            facecolor="#d9d9d9",
                            edgecolor="none",
                            zorder=1,
                        )
                    )
                else:
                    ax.scatter([x], [y_mid], s=_p_size(pval), facecolor=_nes_color(nes), edgecolor="none", zorder=4)

    for x0 in block_x0.values():
        left = x0 - cell_w / 2
        right = x0 + (len(CELL_TYPES) - 0.5) * cell_w
        bottom = y_top - (len(SCENARIOS) - 1) * dy_scenario - cell_h / 2
        top = y_top + cell_h / 2
        ax.add_patch(
            plt.Rectangle((left, bottom), right - left, top - bottom, fill=False, edgecolor="#222222", linewidth=0.65, zorder=5)
        )
        for k in range(1, len(CELL_TYPES)):
            xline = x0 + (k - 0.5) * cell_w
            ax.plot([xline, xline], [bottom, top], color="#d0d0d0", linewidth=0.35, zorder=2)
        for k in range(1, len(SCENARIOS)):
            yline = y_top - (k - 0.5) * cell_h
            ax.plot([left, right], [yline, yline], color="#d0d0d0", linewidth=0.35, zorder=2)

    legend_x = 0.750
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
    for i, (label, size) in enumerate([("P < 0.01", 28), ("P > 0.01", 14), ("N/A", 13)]):
        y = 0.260 - i * 0.030
        if label == "N/A":
            ax.add_patch(plt.Rectangle((0.640, y - 0.010), 0.020, 0.020, facecolor="#d9d9d9", edgecolor="none"))
        else:
            ax.scatter([0.650], [y], s=size, facecolor="white", edgecolor="#111111", linewidth=0.55)
        ax.text(0.672, y, label, fontsize=4.4, ha="left", va="center")

    out_png = OUT_DIR / "fig2e_stage3_profile_mask_route2_ce9_ce10.png"
    out_pdf = OUT_DIR / "fig2e_stage3_profile_mask_route2_ce9_ce10.pdf"
    fig.savefig(out_png, dpi=600, bbox_inches="tight", pad_inches=0.005)
    fig.savefig(out_pdf, bbox_inches="tight", pad_inches=0.005)
    plt.close(fig)
    return out_png


def main() -> int:
    metrics = pd.read_csv(RESULT_DIR / "fig2e_stage3_profile_mask_metrics.csv")
    out = plot_route2_only(metrics)
    print(f"[OK] wrote: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
