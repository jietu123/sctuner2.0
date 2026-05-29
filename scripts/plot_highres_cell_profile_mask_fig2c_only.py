from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MANIFEST = (
    PROJECT_ROOT
    / "visualizations"
    / "highres_profile_mask_mapping"
    / "highres_profile_mask_mapping_manifest.csv"
)
VIS_DIR = PROJECT_ROOT / "visualizations" / "highres_profile_mask_fig2c_only"
RESULT_DIR = PROJECT_ROOT / "result" / "highres_profile_mask_fig2c_only"
N_PANEL_GENES = 30

CURVE_COLOR = "#f0b323"
TICK_COLOR = "#f28e2b"
TRACK_COLOR = "#d8d8d8"
TRACK_EDGE = "#bcbcbc"


def _processed_dir(sample: str) -> Path:
    candidates = [
        PROJECT_ROOT / "data" / "processed" / sample,
        PROJECT_ROOT / "data" / "processed" / "highres_profile_mask" / sample,
    ]
    for path in candidates:
        if (path / "stage1_preprocess" / "exported").exists():
            return path
    raise FileNotFoundError(f"Cannot resolve processed directory for {sample}")


def _read_genes(processed: Path) -> list[str]:
    panel_path = processed / "stage1_preprocess" / "fig2d_profile_mask_gene_panel.csv"
    panel = pd.read_csv(panel_path)
    expr_cols = set(
        pd.read_csv(
            processed / "stage1_preprocess" / "exported" / "st_expression_normalized.csv",
            nrows=0,
        ).columns
    )
    genes: list[str] = []
    for gene in panel["gene"].astype(str):
        if gene in expr_cols and gene not in genes:
            genes.append(gene)
        if len(genes) >= N_PANEL_GENES:
            break
    if not genes:
        raise ValueError(f"No usable profile-mask marker genes in {processed}")
    return genes


def _masked_support(processed: Path, genes: list[str]) -> pd.DataFrame:
    expr = pd.read_csv(
        processed / "stage1_preprocess" / "exported" / "st_expression_normalized.csv",
        usecols=["spot_id", *genes],
    )
    expr["spot_id"] = expr["spot_id"].astype(str)
    expr["masked_support"] = expr[genes].mean(axis=1)
    support = expr[["spot_id", "masked_support"]].copy()
    lo = float(support["masked_support"].min())
    hi = float(support["masked_support"].max())
    support["masked_support_norm"] = (support["masked_support"] - lo) / (hi - lo + 1e-12)
    return support


def _cell_type_marker_scores(processed: Path, genes: list[str]) -> dict[str, float]:
    export = processed / "stage1_preprocess" / "exported"
    meta = pd.read_csv(export / "sc_metadata.csv", usecols=["cell_id", "cell_type"])
    expr = pd.read_csv(export / "sc_expression_normalized.csv", usecols=["cell_id", *genes])
    expr["target_marker_score"] = expr[genes].mean(axis=1)
    joined = expr[["cell_id", "target_marker_score"]].merge(meta, on="cell_id", how="inner")
    return joined.groupby("cell_type")["target_marker_score"].mean().to_dict()


def _assignment_path(sample: str, suffix: str) -> Path:
    return (
        PROJECT_ROOT
        / "result"
        / sample
        / f"stage4_cytospace_{suffix}_highres"
        / "cytospace_output"
        / "cell_type_assignments_by_spot.csv"
    )


def _predicted_target_like_score(sample: str, suffix: str, type_scores: dict[str, float]) -> pd.DataFrame:
    path = _assignment_path(sample, suffix)
    df = pd.read_csv(path)
    df = df.rename(columns={df.columns[0]: "spot_id"})
    df["spot_id"] = df["spot_id"].astype(str)

    score = np.zeros(len(df), dtype=float)
    total = np.zeros(len(df), dtype=float)
    for cell_type, marker_score in type_scores.items():
        if cell_type not in df.columns:
            continue
        counts = pd.to_numeric(df[cell_type], errors="coerce").fillna(0.0).to_numpy(dtype=float)
        score += counts * float(marker_score)
        total += counts

    if "Total cells" in df.columns:
        reported_total = pd.to_numeric(df["Total cells"], errors="coerce").to_numpy(dtype=float)
        total = np.where(np.isfinite(reported_total), reported_total, total)
    score = score / np.maximum(total, 1.0)
    return df[["spot_id"]].assign(predicted_target_like_score=score)


def _running_enrichment_curve(scores: pd.Series, frac: float = 0.10) -> tuple[np.ndarray, np.ndarray, float]:
    hits = scores <= scores.quantile(frac)
    hits_np = hits.to_numpy(dtype=bool)
    n = len(hits_np)
    nh = int(hits_np.sum())
    if nh == 0 or nh == n:
        running = np.zeros(n, dtype=float)
        return running, hits_np, 0.0
    hit_w = 1.0 / nh
    miss_w = 1.0 / (n - nh)
    running = np.cumsum(np.where(hits_np, hit_w, -miss_w))
    return running, hits_np, float(running.max())


def _build_spot_table(row: pd.Series) -> tuple[pd.DataFrame, dict[str, float]]:
    sample = str(row["profile_mask_sample"])
    processed = _processed_dir(sample)
    genes = _read_genes(processed)
    support = _masked_support(processed, genes)
    type_scores = _cell_type_marker_scores(processed, genes)
    baseline = _predicted_target_like_score(sample, "baseline", type_scores).rename(
        columns={"predicted_target_like_score": "baseline_reconstructed_score"}
    )
    route2 = _predicted_target_like_score(sample, "route2", type_scores).rename(
        columns={"predicted_target_like_score": "route2_reconstructed_score"}
    )
    spot_df = support.merge(baseline, on="spot_id", how="inner").merge(route2, on="spot_id", how="inner")
    spot_df = spot_df.sort_values("masked_support_norm", ascending=True).reset_index(drop=True)
    spot_df["rank"] = np.arange(1, len(spot_df) + 1)

    _, _, baseline_peak = _running_enrichment_curve(spot_df["baseline_reconstructed_score"])
    _, _, route2_peak = _running_enrichment_curve(spot_df["route2_reconstructed_score"])
    metrics = {
        "baseline_peak_es": baseline_peak,
        "route2_peak_es": route2_peak,
        "delta_peak_es": route2_peak - baseline_peak,
        "n_spots": float(len(spot_df)),
        "n_genes": float(len(genes)),
        "support_sd": float(spot_df["masked_support_norm"].std(ddof=0)),
    }
    return spot_df, metrics


def _draw_row(fig: plt.Figure, cell, df: pd.DataFrame, score_col: str, method_label: str) -> float:
    sub_gs = cell.subgridspec(3, 1, height_ratios=[2.35, 0.42, 0.34], hspace=0.025)
    ax_curve = fig.add_subplot(sub_gs[0, 0])
    ax_hits = fig.add_subplot(sub_gs[1, 0], sharex=ax_curve)
    ax_track = fig.add_subplot(sub_gs[2, 0], sharex=ax_curve)

    running, hits, peak_es = _running_enrichment_curve(df[score_col], frac=0.10)
    x = df["rank"].to_numpy(dtype=float)

    if len(running) >= 50:
        win = max(5, int(round(len(running) * 0.006)))
        smooth = pd.Series(running).rolling(win, center=True, min_periods=1).mean().to_numpy()
    else:
        smooth = running

    ax_curve.plot(x, smooth, color=CURVE_COLOR, linewidth=1.95)
    ax_curve.axhline(0, color="#9a9a9a", linestyle="--", linewidth=0.7)
    ax_curve.set_xlim(1, len(df))
    ax_curve.set_yticks([0.0])
    ax_curve.set_yticklabels(["0"], fontsize=7.8)
    ax_curve.set_ylabel("Enrichment", fontsize=7.4)
    ax_curve.text(0.965, 0.90, f"NES = {peak_es:.3f}", transform=ax_curve.transAxes, ha="right", va="top", fontsize=7.0)
    sns.despine(ax=ax_curve, top=True, right=True)
    ax_curve.tick_params(axis="x", bottom=False, labelbottom=False)

    ax_hits.vlines(df.loc[hits, "rank"], 0.05, 0.95, color=TICK_COLOR, linewidth=0.65)
    ax_hits.set_ylim(0, 1)
    ax_hits.set_yticks([])
    sns.despine(ax=ax_hits, left=True, bottom=True, right=True, top=True)
    ax_hits.tick_params(axis="x", bottom=False, labelbottom=False)

    support = df["masked_support_norm"].to_numpy(dtype=float)
    ax_track.fill_between(x, 0, support, color=TRACK_COLOR, linewidth=0)
    ax_track.plot(x, support, color=TRACK_EDGE, linewidth=0.75)
    ax_track.set_ylim(0, 1.0)
    ax_track.set_yticks([])
    sns.despine(ax=ax_track, left=True, right=True, top=True)
    ax_track.tick_params(axis="x", bottom=False, labelbottom=False)

    ax_curve.text(0.0, 1.02, method_label, transform=ax_curve.transAxes, ha="left", va="bottom", fontsize=7.9, color="#555555")
    return peak_es


def main() -> int:
    plt.rcParams["svg.fonttype"] = "none"
    sns.set_theme(style="white", font="DejaVu Sans")
    VIS_DIR.mkdir(parents=True, exist_ok=True)
    RESULT_DIR.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(MANIFEST)
    tables: dict[str, pd.DataFrame] = {}
    metric_rows = []
    for _, row in manifest.iterrows():
        sample = str(row["profile_mask_sample"])
        spot_df, metrics = _build_spot_table(row)
        tables[sample] = spot_df
        metric_rows.append(
            {
                "raw_sample": row["raw_sample"],
                "profile_mask_sample": sample,
                "masked_target_type": row["masked_target_type"],
                **metrics,
            }
        )

    metrics_df = pd.DataFrame(metric_rows).sort_values(
        ["delta_peak_es", "route2_peak_es", "support_sd"], ascending=[False, False, False]
    )
    metrics_df.to_csv(RESULT_DIR / "highres_cell_profile_mask_fig2c_candidate_metrics.csv", index=False)
    selected = metrics_df.iloc[0]
    selected_sample = str(selected["profile_mask_sample"])
    spot_df = tables[selected_sample]
    spot_df.to_csv(RESULT_DIR / f"{selected_sample}_spot_fig2c_foundation.csv", index=False)

    scenario_label = f"{selected['raw_sample']}  {selected['masked_target_type']}"
    fig = plt.figure(figsize=(8.45, 5.0), dpi=220)
    outer = fig.add_gridspec(2, 1, hspace=0.22)
    _draw_row(fig, outer[0, 0], spot_df, "baseline_reconstructed_score", "CytoSPACE")
    _draw_row(fig, outer[1, 0], spot_df, "route2_reconstructed_score", "SVTuner + CytoSPACE")

    fig.text(0.015, 0.975, "c", ha="left", va="top", fontsize=12.0, weight="bold")
    fig.text(0.048, 0.976, "Spatial enrichment of mapped target-like signal", ha="left", va="top", fontsize=8.7, weight="bold")
    fig.text(0.048, 0.952, scenario_label, ha="left", va="top", fontsize=7.2, color="#444444")
    fig.text(0.08, 0.041, "Lower masked support", ha="left", va="bottom", fontsize=7.6)
    fig.text(0.50, 0.041, "Relative support rank", ha="center", va="bottom", fontsize=7.8)
    fig.text(0.92, 0.041, "Higher masked support", ha="right", va="bottom", fontsize=7.6)
    fig.subplots_adjust(left=0.075, right=0.99, top=0.910, bottom=0.085)

    for ext in ("png", "pdf", "svg"):
        out = VIS_DIR / f"fig2_panel_c_highres_cell_profile_mask.{ext}"
        fig.savefig(out, bbox_inches="tight")
        print(f"[OK] wrote: {out}")
    plt.close(fig)
    print(
        "[SELECTED] "
        f"{selected_sample}: target={selected['masked_target_type']}, "
        f"baseline_ES={selected['baseline_peak_es']:.4f}, "
        f"route2_ES={selected['route2_peak_es']:.4f}, "
        f"delta={selected['delta_peak_es']:.4f}"
    )
    print(f"[OK] wrote: {RESULT_DIR / 'highres_cell_profile_mask_fig2c_candidate_metrics.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
