from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


PROJECT_ROOT = Path(__file__).resolve().parents[1]
FOUNDATION_DIR = PROJECT_ROOT / "result" / "real_profile_mask_foundation"
VIS_DIR = PROJECT_ROOT / "visualizations" / "simulations" / "real_profile_mask_fig2c_only"

# Use the previously identified representative scenario where route2 shows cleaner low-support suppression.
REPRESENTATIVE_SAMPLE = "ffpe_mouse_brain_sagittal_real_profile_mask_astrocyte"
SCENARIO_LABEL = "Mouse brain  Microglia"

# Warm orange/yellow palette closer to paper Fig.2c.
CURVE_COLOR = "#f0b323"
TICK_COLOR = "#f28e2b"
TRACK_COLOR = "#d8d8d8"
TRACK_EDGE = "#bcbcbc"


def _running_enrichment_curve(scores: pd.Series, frac: float = 0.10) -> tuple[np.ndarray, np.ndarray]:
    hits = scores <= scores.quantile(frac)
    hits_np = hits.to_numpy(dtype=bool)
    n = len(hits_np)
    nh = int(hits_np.sum())
    if nh == 0 or nh == n:
        return np.zeros(n, dtype=float), hits_np
    hit_w = 1.0 / nh
    miss_w = 1.0 / (n - nh)
    running = np.cumsum(np.where(hits_np, hit_w, -miss_w))
    return running, hits_np


def _draw_row(fig: plt.Figure, cell, df: pd.DataFrame, score_col: str, title: str) -> None:
    sub_gs = cell.subgridspec(3, 1, height_ratios=[2.35, 0.42, 0.34], hspace=0.025)
    ax_curve = fig.add_subplot(sub_gs[0, 0])
    ax_hits = fig.add_subplot(sub_gs[1, 0], sharex=ax_curve)
    ax_track = fig.add_subplot(sub_gs[2, 0], sharex=ax_curve)

    running, hits = _running_enrichment_curve(df[score_col], frac=0.10)
    x = df["rank"].to_numpy(dtype=float)

    ax_curve.plot(x, running, color=CURVE_COLOR, linewidth=1.75)
    ax_curve.axhline(0, color="#9a9a9a", linestyle="--", linewidth=0.7)
    ax_curve.set_xlim(1, len(df))
    ax_curve.set_yticks([0.0])
    ax_curve.set_yticklabels(["0"], fontsize=7.8)
    ax_curve.set_ylabel("Enrichment", fontsize=7.4)
    ax_curve.text(0.965, 0.90, f"NES = {running[np.argmax(np.abs(running))]:.3f}", transform=ax_curve.transAxes, ha="right", va="top", fontsize=7.0)
    sns.despine(ax=ax_curve, top=True, right=True)
    ax_curve.tick_params(axis="x", bottom=False, labelbottom=False)

    ax_hits.vlines(df.loc[hits, "rank"], 0.05, 0.95, color=TICK_COLOR, linewidth=0.65)
    ax_hits.set_ylim(0, 1)
    ax_hits.set_yticks([])
    sns.despine(ax=ax_hits, left=True, bottom=True, right=True, top=True)
    ax_hits.tick_params(axis="x", bottom=False, labelbottom=False)

    support = df["masked_support_norm"].to_numpy(dtype=float)
    support = (support - support.min()) / (support.max() - support.min() + 1e-12)
    ax_track.fill_between(x, 0, support, color=TRACK_COLOR, linewidth=0)
    ax_track.plot(x, support, color=TRACK_EDGE, linewidth=0.75)
    ax_track.set_ylim(0, 1.0)
    ax_track.set_yticks([])
    sns.despine(ax=ax_track, left=True, right=True, top=True)
    ax_track.tick_params(axis="x", bottom=False, labelbottom=False)


def main() -> int:
    sns.set_theme(style="white", font="DejaVu Sans")
    VIS_DIR.mkdir(parents=True, exist_ok=True)

    spot_df = pd.read_csv(FOUNDATION_DIR / REPRESENTATIVE_SAMPLE / "spot_foundation.csv")
    spot_df = spot_df.sort_values("masked_support_norm", ascending=True).reset_index(drop=True).copy()
    spot_df["rank"] = np.arange(1, len(spot_df) + 1)

    fig = plt.figure(figsize=(8.45, 5.0), dpi=220)
    outer = fig.add_gridspec(2, 1, hspace=0.22)

    _draw_row(fig, outer[0, 0], spot_df, "baseline_reconstructed_score", "CytoSPACE")
    _draw_row(fig, outer[1, 0], spot_df, "route2_reconstructed_score", "SVTuner + CytoSPACE")

    fig.text(0.015, 0.975, "c", ha="left", va="top", fontsize=12.0, weight="bold")
    fig.text(0.048, 0.976, "Spatial enrichment of mapped target-like signal", ha="left", va="top", fontsize=8.7, weight="bold")
    fig.text(0.048, 0.952, SCENARIO_LABEL, ha="left", va="top", fontsize=7.2, color="#444444")
    fig.text(0.075, 0.900, "CytoSPACE", ha="left", va="bottom", fontsize=7.9, color="#555555")
    fig.text(0.075, 0.443, "SVTuner + CytoSPACE", ha="left", va="bottom", fontsize=7.9, color="#555555")
    fig.text(0.08, 0.041, "Lower masked support", ha="left", va="bottom", fontsize=7.6)
    fig.text(0.50, 0.041, "Relative support rank", ha="center", va="bottom", fontsize=7.8)
    fig.text(0.92, 0.041, "Higher masked support", ha="right", va="bottom", fontsize=7.6)
    fig.subplots_adjust(left=0.075, right=0.99, top=0.910, bottom=0.085)

    out_png = VIS_DIR / "fig2_panel_c_real_profile_mask.png"
    out_pdf = VIS_DIR / "fig2_panel_c_real_profile_mask.pdf"
    fig.savefig(out_png, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

    print(f"[OK] wrote: {out_png}")
    print(f"[OK] wrote: {out_pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
