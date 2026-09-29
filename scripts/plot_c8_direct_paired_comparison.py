from __future__ import annotations

from pathlib import Path
import textwrap

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "result" / "c8_direct_paired_comparison"
OUTPUT = ROOT / "visualizations" / "c8_direct_paired_comparison"

PALETTE = {
    "cytospace": "#C47A2C",       # warm amber
    "svtuner": "#087F8C",         # deep teal
    "profile": "#7563A8",         # profile-mask violet
    "profile_light": "#DCD5EA",
    "merscope": "#4F8A72",        # muted green
    "merscope_light": "#D7E7DF",
    "unfavorable": "#9AA3AA",
    "neutral": "#D8DDE1",
    "grid": "#E7EAEC",
    "text": "#26323A",
    "secondary": "#68737B",
    "negative": "#B86E62",
    "white": "#FFFFFF",
}

METRICS = [
    "A1_reciprocal_suppression",
    "A2_cosine_similarity",
    "A3_ecotyper_experiment_mean",
    "B1_merscope_suppression",
    "B2_merscope_peak_es",
]

SHORT_TITLES = {
    "A1_reciprocal_suppression": "Reciprocal suppression",
    "A2_cosine_similarity": "Expression cosine similarity",
    "A3_ecotyper_experiment_mean": "EcoTyper enrichment",
    "B1_merscope_suppression": "Low-support suppression",
    "B2_merscope_peak_es": "Peak ES",
}

ROW_LABELS = {
    "A1_reciprocal_suppression": "Reciprocal\nsuppression",
    "A2_cosine_similarity": "Expression cosine\nsimilarity",
    "A3_ecotyper_experiment_mean": "EcoTyper\nenrichment",
    "B1_merscope_suppression": "MERSCOPE\nlow-support\nsuppression",
    "B2_merscope_peak_es": "MERSCOPE\nPeak ES",
}

DATASET_LABELS = {
    "highres_humanbreastcancerpatient1_profile_mask_monocytes_and_macrophages": "Breast cancer P1\nMono./macro.",
    "highres_humancoloncancerpatient1_profile_mask_fibroblasts": "Colon cancer P1\nFibroblasts",
    "highres_humanlungcancerpatient1_profile_mask_plasma_cells": "Lung cancer P1\nPlasma cells",
    "highres_humanmelanomapatient1_profile_mask_fibroblasts": "Melanoma P1\nFibroblasts",
    "highres_humanmelanomapatient2_profile_mask_b_cells": "Melanoma P2\nB cells",
}


def _style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7.4,
        "axes.titlesize": 8.2,
        "axes.labelsize": 7.5,
        "xtick.labelsize": 6.8,
        "ytick.labelsize": 6.8,
        "legend.fontsize": 6.8,
        "axes.edgecolor": PALETTE["secondary"],
        "axes.linewidth": 0.65,
        "xtick.color": PALETTE["text"],
        "ytick.color": PALETTE["text"],
        "axes.labelcolor": PALETTE["text"],
        "text.color": PALETTE["text"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def _clean_axis(ax: plt.Axes, grid: bool = True) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid:
        ax.grid(axis="y", color=PALETTE["grid"], linewidth=0.55, zorder=0)
    ax.tick_params(length=2.4, width=0.55)


def _fmt_num(x: float, digits: int = 4, signed: bool = False) -> str:
    return f"{x:+.{digits}f}" if signed else f"{x:.{digits}f}"


def _pair_columns(df: pd.DataFrame) -> tuple[str, str]:
    candidates = [
        ("cytospace", "svtuner"),
        ("cytospace_value", "svtuner_value"),
        ("cytospace_score", "svtuner_score"),
    ]
    for pair in candidates:
        if all(c in df.columns for c in pair):
            return pair
    raise KeyError(f"Cannot identify method columns in {list(df.columns)}")


def _unit_column(df: pd.DataFrame) -> str:
    for c in ("independent_unit_id", "pair_id", "dataset"):
        if c in df.columns:
            return c
    raise KeyError("No independent-unit identifier found")


def _delta_column(df: pd.DataFrame) -> str:
    for c in ("delta_raw", "delta", "delta_favorable"):
        if c in df.columns:
            return c
    cy, sv = _pair_columns(df)
    df["delta_raw"] = df[sv] - df[cy]
    return "delta_raw"


def _metric_subset(canonical: pd.DataFrame, metric_id: str) -> pd.DataFrame:
    out = canonical.loc[canonical["metric_id"] == metric_id].copy()
    if out.empty:
        raise ValueError(f"No rows found for {metric_id}")
    return out


def _draw_panel_a(fig: plt.Figure, spec, canonical: pd.DataFrame, stats: pd.DataFrame) -> None:
    gs = spec.subgridspec(1, 3, wspace=0.40)
    for idx, metric_id in enumerate(METRICS[:3]):
        ax = fig.add_subplot(gs[0, idx])
        sub = _metric_subset(canonical, metric_id)
        cy_col, sv_col = _pair_columns(sub)
        delta_col = _delta_column(sub)
        y_all = np.r_[sub[cy_col].to_numpy(float), sub[sv_col].to_numpy(float)]
        span = max(np.ptp(y_all), max(abs(np.mean(y_all)), 1.0) * 0.025)
        pad = span * 0.19
        for _, row in sub.iterrows():
            favorable = float(row[delta_col]) > 0
            line_color = PALETTE["profile"] if favorable else PALETTE["unfavorable"]
            ax.plot([0, 1], [row[cy_col], row[sv_col]], color=line_color,
                    alpha=0.58 if favorable else 0.50, linewidth=0.80, zorder=1)
        ax.scatter(np.zeros(len(sub)), sub[cy_col], s=18, color=PALETTE["cytospace"],
                   edgecolor="white", linewidth=0.35, zorder=3)
        ax.scatter(np.ones(len(sub)), sub[sv_col], s=18, color=PALETTE["svtuner"],
                   edgecolor="white", linewidth=0.35, zorder=3)
        ax.set_xlim(-0.32, 1.32)
        ax.set_ylim(float(y_all.min() - pad), float(y_all.max() + pad * 1.90))
        ax.set_xticks([0, 1], ["CytoSPACE", "SVTuner"])
        ax.set_title(SHORT_TITLES[metric_id], pad=5.5, fontweight="semibold")
        st = stats.loc[stats["metric_id"] == metric_id].iloc[0]
        annotation = (
            f"Δmean {_fmt_num(st.mean_delta_raw, 4, True)} "
            f"[{st.mean_delta_ci_low:.4f}, {st.mean_delta_ci_high:.4f}]\n"
            f"P={st.wilcoxon_p_two_sided:.4g} · {int(st.wins)}/{int(st.n_pairs)} favorable"
        )
        ax.text(0.5, 0.985, annotation, transform=ax.transAxes, ha="center", va="top",
                color=PALETTE["secondary"], fontsize=6.3, linespacing=1.25)
        _clean_axis(ax, grid=True)
        ax.spines["left"].set_color(PALETTE["neutral"])
        ax.spines["bottom"].set_color(PALETTE["neutral"])
        if idx == 0:
            ax.text(-0.27, 1.18, "A", transform=ax.transAxes, fontsize=12,
                    fontweight="bold", va="top")
            ax.text(-0.12, 1.18, "Low-resolution paired performance",
                    transform=ax.transAxes, fontsize=8.5, fontweight="semibold", va="top")


def _draw_panel_b(fig: plt.Figure, spec, stats: pd.DataFrame) -> None:
    gs = spec.subgridspec(1, 2, width_ratios=[0.96, 1.55], wspace=0.06)
    ax = fig.add_subplot(gs[0, 0])
    tx = fig.add_subplot(gs[0, 1], sharey=ax)
    ordered = stats.set_index("metric_id").loc[METRICS].reset_index()
    y = np.arange(len(ordered))[::-1]
    colors = [PALETTE["profile"]] * 3 + [PALETTE["merscope"]] * 2
    ax.axvline(0, color=PALETTE["unfavorable"], linestyle=(0, (3, 2)), linewidth=0.75, zorder=0)
    ax.scatter(ordered["cohen_dz"], y, s=30, c=colors, edgecolor="white",
               linewidth=0.5, zorder=3)
    ax.set_yticks(y, [ROW_LABELS[m] for m in ordered["metric_id"]])
    ax.set_xlabel("Cohen $d_z$")
    xmin = min(-0.12, float(ordered["cohen_dz"].min()) - 0.10)
    xmax = float(ordered["cohen_dz"].max()) + 0.18
    ax.set_xlim(xmin, xmax)
    ax.grid(axis="x", color=PALETTE["grid"], linewidth=0.55)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.text(-0.40, 1.16, "B", transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")
    ax.text(-0.19, 1.16, "Direct effect summary", transform=ax.transAxes,
            fontsize=8.5, fontweight="semibold", va="top")

    tx.set_xlim(0, 1)
    tx.set_ylim(ax.get_ylim())
    tx.axis("off")
    headers = [(0.00, "Mean Δ [95% CI]"), (0.69, "Fav."), (0.87, "P")]
    for x, label in headers:
        tx.text(x, y.max() + 0.65, label, ha="left" if x == 0 else "center",
                va="center", fontsize=6.2, fontweight="semibold", color=PALETTE["secondary"])
    for yi, row in zip(y, ordered.itertuples(index=False)):
        ci_text = f"{row.mean_delta_raw:+.4f} [{row.mean_delta_ci_low:.4f}, {row.mean_delta_ci_high:.4f}]"
        tx.text(0.00, yi, ci_text, ha="left", va="center", fontsize=5.85)
        tx.text(0.69, yi, f"{int(row.wins)}/{int(row.n_pairs)}", ha="center", va="center", fontsize=6.0)
        tx.text(0.87, yi, f"{row.wilcoxon_p_two_sided:.4g}", ha="center", va="center", fontsize=6.0)
    tx.plot([0, 0.99], [y.max() + 0.38] * 2, color=PALETTE["neutral"], linewidth=0.6)


def _short_dataset(row: pd.Series) -> str:
    if row.get("dataset") in DATASET_LABELS:
        return DATASET_LABELS[row["dataset"]]
    unit = str(row.get("independent_unit_id", row.get("pair_id", row.get("dataset", ""))))
    for key, label in DATASET_LABELS.items():
        if key in unit:
            return label
    return unit.replace("highres_", "").replace("_profile_mask_", "\n").replace("_", " ")


def _draw_dumbbell(ax: plt.Axes, sub: pd.DataFrame, title: str, win_text: str) -> None:
    cy_col, sv_col = _pair_columns(sub)
    delta_col = "delta_favorable" if "delta_favorable" in sub.columns else _delta_column(sub)
    sub = sub.sort_values(delta_col, ascending=True).copy()
    y = np.arange(len(sub))
    labels = [_short_dataset(r) for _, r in sub.iterrows()]
    for yi, (_, row) in enumerate(sub.iterrows()):
        ax.plot([row[cy_col], row[sv_col]], [yi, yi], color=PALETTE["neutral"], linewidth=1.6, zorder=1)
    ax.scatter(sub[cy_col], y, s=24, color=PALETTE["cytospace"], edgecolor="white", linewidth=0.4,
               label="CytoSPACE", zorder=3)
    ax.scatter(sub[sv_col], y, s=24, color=PALETTE["svtuner"], edgecolor="white", linewidth=0.4,
               label="SVTuner", zorder=3)
    ax.set_yticks(y, labels)
    ax.set_title(title, loc="left", pad=7, fontweight="semibold")
    ax.text(1.0, 1.03, win_text, transform=ax.transAxes, ha="right", va="bottom",
            color=PALETTE["merscope"], fontsize=6.6, fontweight="semibold")
    _clean_axis(ax, grid=True)
    ax.grid(axis="x", color=PALETTE["grid"], linewidth=0.55)
    ax.grid(axis="y", visible=False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)


def _draw_panel_c(fig: plt.Figure, spec, canonical: pd.DataFrame) -> None:
    gs = spec.subgridspec(1, 2, wspace=0.52)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    _draw_dumbbell(ax1, _metric_subset(canonical, METRICS[3]), "Low-support suppression", "5/5 favorable")
    _draw_dumbbell(ax2, _metric_subset(canonical, METRICS[4]), "Peak ES", "2/5 favorable")
    ax1.set_xlabel("Score")
    ax2.set_xlabel("Peak ES")
    handles, labels = ax1.get_legend_handles_labels()
    ax2.legend(handles, labels, frameon=False, ncol=2, loc="lower right",
               bbox_to_anchor=(1.0, -0.29), handletextpad=0.35, columnspacing=0.9)
    ax1.text(-0.28, 1.18, "C", transform=ax1.transAxes, fontsize=12, fontweight="bold", va="top")
    ax1.text(-0.11, 1.18, "MERSCOPE paired outcomes", transform=ax1.transAxes,
             fontsize=8.5, fontweight="semibold", va="top")


def _readout_short_label(pair_id: str) -> str:
    p = pair_id.lower()
    mappings = [
        ("melanoma_slide1", "Melanoma S1"),
        ("melanoma_slide2", "Melanoma S2"),
        ("brca_er_her2_fresh_frozen", "BRCA ER/HER2 FF"),
        ("brca_her2_ffpe", "BRCA HER2 FFPE"),
        ("brca_tnbc_fresh_frozen", "BRCA TNBC FF"),
        ("crc_fresh_frozen", "CRC FF"),
    ]
    for key, label in mappings:
        if key in p:
            return label
    return pair_id.replace("_", " ")


def _draw_panel_d(fig: plt.Figure, spec, readouts: pd.DataFrame) -> None:
    ax = fig.add_subplot(spec)
    unit_col = "pair_id" if "pair_id" in readouts.columns else _unit_column(readouts)
    readout_col = "readout" if "readout" in readouts.columns else "cell_type"
    delta_col = "delta_favorable" if "delta_favorable" in readouts.columns else _delta_column(readouts)
    pivot = readouts.pivot(index=unit_col, columns=readout_col, values=delta_col)
    desired_cols = [c for c in ("CD4 T cells", "CD8 T cells") if c in pivot.columns]
    if len(desired_cols) == 2:
        pivot = pivot[desired_cols]
    pivot = pivot.sort_index()
    vals = pivot.to_numpy(float)
    vmax = max(abs(np.nanmin(vals)), abs(np.nanmax(vals)), 1e-6)
    cmap = LinearSegmentedColormap.from_list("ecotyper_delta", [PALETTE["negative"], "#F4F2EF", PALETTE["profile"]])
    im = ax.imshow(vals, aspect="auto", cmap=cmap, norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax))
    ax.set_xticks(np.arange(pivot.shape[1]), pivot.columns, rotation=0)
    ax.set_yticks(np.arange(pivot.shape[0]), [_readout_short_label(str(x)) for x in pivot.index])
    ax.tick_params(length=0)
    for i in range(vals.shape[0]):
        for j in range(vals.shape[1]):
            v = vals[i, j]
            rgba = cmap(TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)(v))
            luminance = 0.2126 * rgba[0] + 0.7152 * rgba[1] + 0.0722 * rgba[2]
            ax.text(j, i, f"{v:+.3f}", ha="center", va="center",
                    fontsize=6.0, color="white" if luminance < 0.58 else PALETTE["text"])
    for edge in np.arange(-0.5, pivot.shape[0], 1):
        ax.axhline(edge, color="white", linewidth=1.0)
    for edge in np.arange(-0.5, pivot.shape[1], 1):
        ax.axvline(edge, color="white", linewidth=1.0)
    for s in ax.spines.values():
        s.set_visible(False)
    cbar = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.075, pad=0.35, aspect=25)
    cbar.set_label("SVTuner − CytoSPACE (Δ NES)", fontsize=6.3)
    cbar.ax.tick_params(labelsize=5.9, length=2)
    cbar.outline.set_linewidth(0.45)
    ax.text(-0.25, 1.16, "D", transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")
    ax.text(-0.04, 1.16, "EcoTyper readout differences", transform=ax.transAxes,
            fontsize=8.5, fontweight="semibold", va="top")
    ax.text(0.5, -0.15, "12 readouts shown descriptively;\npaired inference uses 6 experiment means.",
            transform=ax.transAxes, ha="center", va="top", fontsize=5.7,
            linespacing=1.15, color=PALETTE["secondary"])


def _build_source_table(canonical: pd.DataFrame, readouts: pd.DataFrame, stats: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    cy_col, sv_col = _pair_columns(canonical)
    unit_col = _unit_column(canonical)
    delta_col = _delta_column(canonical)
    for _, row in canonical.iterrows():
        metric_id = row["metric_id"]
        panel = "A" if metric_id in METRICS[:3] else "C"
        stat = stats.loc[stats["metric_id"] == metric_id].iloc[0]
        for method, col in (("CytoSPACE", cy_col), ("SVTuner", sv_col)):
            rows.append({
                "panel": panel,
                "metric": stat["metric_label"],
                "independent_unit": row[unit_col],
                "method": method,
                "value": row[col],
                "delta": row[delta_col],
                "inferential_unit": "Yes",
                "ci_low": np.nan,
                "ci_high": np.nan,
                "cohen_dz": np.nan,
                "favorable_pairs": np.nan,
                "n_pairs": np.nan,
                "wilcoxon_p_two_sided": np.nan,
            })
    for _, stat in stats.iterrows():
        rows.append({
            "panel": "B",
            "metric": stat["metric_label"],
            "independent_unit": "summary",
            "method": "Cohen dz",
            "value": stat["cohen_dz"],
            "delta": stat["mean_delta_raw"],
            "inferential_unit": "Summary of independent pairs",
            "ci_low": stat["mean_delta_ci_low"],
            "ci_high": stat["mean_delta_ci_high"],
            "cohen_dz": stat["cohen_dz"],
            "favorable_pairs": stat["wins"],
            "n_pairs": stat["n_pairs"],
            "wilcoxon_p_two_sided": stat["wilcoxon_p_two_sided"],
        })
    rcy, rsv = _pair_columns(readouts)
    runit = "pair_id" if "pair_id" in readouts.columns else _unit_column(readouts)
    rdelta = _delta_column(readouts)
    for _, row in readouts.iterrows():
        metric = f"EcoTyper {row.get('readout', 'readout')} normalized enrichment"
        for method, col in (("CytoSPACE", rcy), ("SVTuner", rsv)):
            rows.append({
                "panel": "D",
                "metric": metric,
                "independent_unit": row[runit],
                "method": method,
                "value": row[col],
                "delta": row[rdelta],
                "inferential_unit": "No; descriptive readout",
                "ci_low": np.nan,
                "ci_high": np.nan,
                "cohen_dz": np.nan,
                "favorable_pairs": np.nan,
                "n_pairs": np.nan,
                "wilcoxon_p_two_sided": np.nan,
            })
    return pd.DataFrame(rows)


def _build_stats_table(stats: pd.DataFrame) -> pd.DataFrame:
    ordered = stats.set_index("metric_id").loc[METRICS].reset_index()
    return pd.DataFrame({
        "Metric": ordered["metric_label"],
        "Evidence": ordered["evidence_layer"],
        "n": ordered["n_pairs"].astype(int),
        "CytoSPACE mean": ordered["cytospace_mean"].map(lambda x: f"{x:.4f}"),
        "SVTuner mean": ordered["svtuner_mean"].map(lambda x: f"{x:.4f}"),
        "Mean paired difference": ordered["mean_delta_raw"].map(lambda x: f"{x:+.5f}"),
        "95% CI": ordered.apply(lambda r: f"[{r.mean_delta_ci_low:.5f}, {r.mean_delta_ci_high:.5f}]", axis=1),
        "Cohen dz": ordered["cohen_dz"].map(lambda x: f"{x:.4f}"),
        "Favorable pairs": ordered.apply(lambda r: f"{int(r.wins)}/{int(r.n_pairs)}", axis=1),
        "Two-sided paired Wilcoxon P": ordered["wilcoxon_p_two_sided"].map(lambda x: f"{x:.6g}"),
    })


def main() -> None:
    _style()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    canonical = pd.read_csv(INPUT / "c8_canonical_pairs.csv")
    readouts = pd.read_csv(INPUT / "c8_ecotyper_readout_pairs.csv")
    stats = pd.read_csv(INPUT / "c8_paired_statistics.csv")

    expected_n = dict(zip(METRICS, [10, 10, 6, 5, 5]))
    for metric_id, n in expected_n.items():
        if len(_metric_subset(canonical, metric_id)) != n:
            raise ValueError(f"Unexpected pair count for {metric_id}")
    if len(readouts) != 12 or len(stats) != 5:
        raise ValueError("Frozen C8 input dimensions do not match the expected design")

    fig = plt.figure(figsize=(7.25, 8.15), facecolor="white")
    outer = fig.add_gridspec(3, 1, height_ratios=[1.72, 1.64, 1.75],
                             left=0.105, right=0.975, top=0.955, bottom=0.075, hspace=0.51)
    _draw_panel_a(fig, outer[0], canonical, stats)
    middle = outer[1].subgridspec(1, 2, width_ratios=[1.50, 1.02], wspace=0.34)
    _draw_panel_b(fig, middle[0], stats)
    _draw_panel_d(fig, middle[1], readouts)
    _draw_panel_c(fig, outer[2], canonical)

    png = OUTPUT / "C8_direct_paired_comparison_main.png"
    pdf = OUTPUT / "C8_direct_paired_comparison_main.pdf"
    svg = OUTPUT / "C8_direct_paired_comparison_main.svg"
    fig.savefig(png, dpi=600, facecolor="white")
    fig.savefig(pdf, facecolor="white")
    fig.savefig(svg, facecolor="white")
    plt.close(fig)

    source = _build_source_table(canonical, readouts, stats)
    source.to_csv(OUTPUT / "C8_direct_paired_comparison_source_values.csv", index=False)
    table = _build_stats_table(stats)
    table.to_csv(OUTPUT / "C8_direct_paired_comparison_stats_table.csv", index=False)

    caption = (
        "Direct paired comparison of CytoSPACE and SVTuner across frozen profile-masking benchmarks. "
        "(A) Paired experiment-level values for reciprocal suppression and reconstructed-expression cosine "
        "similarity (n=10 independent sample–target experiments each) and EcoTyper normalized enrichment "
        "(n=6 independent experiments). Twelve CD4/CD8 EcoTyper readouts are shown descriptively in panel D, "
        "whereas inference used the six experiment-level means. (B) Standardized paired effects (Cohen’s dz; "
        "points) with raw mean paired differences, bootstrap 95% confidence intervals, favorable-pair counts, "
        "and two-sided Wilcoxon signed-rank P values; confidence intervals were not estimated for dz. "
        "(C) Dataset–target-level MERSCOPE contrasts for low-support suppression and Peak ES (n=5 independent "
        "dataset–target pairs each). (D) Descriptive EcoTyper readout-level differences. Mean-difference "
        "confidence intervals used 10,000 paired bootstrap resamples (seed 20260927). All hypothesis tests were "
        "two-sided Wilcoxon signed-rank tests. Spots and cells were not treated as independent replicates. "
        "MERSCOPE Peak ES showed no detectable between-method difference (P=1.0)."
    )
    (OUTPUT / "C8_direct_paired_comparison_caption.txt").write_text(caption + "\n", encoding="utf-8")

    print(f"figure_size_inches=7.25x8.15")
    print(f"source_rows={len(source)}")
    print(f"stats_rows={len(table)}")
    print(str(png))
    print(str(pdf))
    print(str(svg))


if __name__ == "__main__":
    main()
