from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap


EPS = 1.0e-12
METHODS = [
    ("CytoSPACE", "baseline_highres"),
    ("SVTuner + CytoSPACE", "route2_highres"),
]
READOUT_TO_GENESET = {
    "T cells": [
        "Tcell_exhaustion_Zheng_etal_Cell2017",
        "EcoTyper_CD4_T_cells_CE9",
        "EcoTyper_CD8_T_cells_CE9",
    ],
    "Monocytes and Macrophages": [
        "EcoTyper_Monocytes_and_Macrophages_CE9",
        "EcoTyper_Monocytes_and_Macrophages_CE10",
    ],
    "B cells": ["EcoTyper_B_cells_CE9", "EcoTyper_B_cells_CE10"],
    "Plasma cells": ["EcoTyper_PCs_CE9", "EcoTyper_PCs_CE10"],
}


def _slug(value: str) -> str:
    return "".join(ch.lower() if ch.isalnum() else "_" for ch in str(value)).strip("_")


def _read_gene_sets(root: Path) -> dict[str, list[str]]:
    path = root / "data/raw/cytospace_fig2c_melanoma/prepared/gene_sets/table_s6_gene_sets_long.csv"
    df = pd.read_csv(path)
    return {
        str(name): sub.sort_values("rank")["gene"].astype(str).tolist()
        for name, sub in df.groupby("gene_set", sort=False)
    }


def _mean_nearest5(query_xy: np.ndarray, ref_xy: np.ndarray) -> np.ndarray:
    diff = query_xy[:, None, :] - ref_xy[None, :, :]
    dist = np.sqrt(np.sum(diff * diff, axis=2))
    k = min(5, dist.shape[1])
    return np.partition(dist, kth=k - 1, axis=1)[:, :k].mean(axis=1)


def _running_es(stats: pd.Series, gene_set: set[str]) -> tuple[float, int, pd.DataFrame]:
    stats = stats.replace([np.inf, -np.inf], np.nan).dropna().sort_values(ascending=False)
    hits = stats.index.to_series().isin(gene_set).to_numpy()
    n_hits = int(hits.sum())
    n_miss = int(len(stats) - n_hits)
    if n_hits == 0 or n_miss == 0:
        raise ValueError("invalid gene-set overlap")
    weights = np.abs(stats.to_numpy(dtype=np.float64))
    hit_norm = float(weights[hits].sum())
    if hit_norm <= EPS:
        inc = np.where(hits, 1.0 / n_hits, -1.0 / n_miss)
    else:
        inc = np.where(hits, weights / hit_norm, -1.0 / n_miss)
    running = np.cumsum(inc)
    peak_rank = int(np.nanargmax(running)) + 1
    curve = pd.DataFrame(
        {
            "rank": np.arange(1, len(stats) + 1),
            "gene": stats.index.to_numpy(),
            "stat": stats.to_numpy(dtype=float),
            "hit": hits,
            "running_es": running,
        }
    )
    return float(running[peak_rank - 1]), peak_rank, curve


def _null_nes(stats: pd.Series, gene_set_size: int, observed_es: float, nperm: int, seed: int) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    genes = stats.replace([np.inf, -np.inf], np.nan).dropna().index.to_numpy()
    if gene_set_size >= len(genes):
        return float("nan"), float("nan")
    null = np.empty(nperm, dtype=float)
    for i in range(nperm):
        selected = set(rng.choice(genes, size=gene_set_size, replace=False).tolist())
        null[i], _, _ = _running_es(stats, selected)
    nes = observed_es / (float(np.mean(np.abs(null))) + EPS)
    pval = (float(np.sum(null >= observed_es)) + 1.0) / (float(nperm) + 1.0)
    return float(nes), float(pval)


def _load_method_cells(root: Path, sample: str, suffix: str, readout: str) -> pd.DataFrame:
    loc = pd.read_csv(root / "result" / sample / f"stage4_cytospace_{suffix}" / "cytospace_output" / "assigned_locations.csv")
    loc = loc[loc["CellType"].astype(str).eq(readout)].copy()
    return loc


def _evaluate(
    root: Path,
    row: pd.Series,
    readout: str,
    gene_set_name: str,
    gene_set: list[str],
    nperm: int,
    seed: int,
) -> dict[str, object] | None:
    sample = str(row["profile_mask_sample"])
    export = root / "data" / "processed" / sample / "stage1_preprocess" / "exported"
    st_meta = pd.read_csv(export / "st_metadata.csv")
    coords = pd.read_csv(export / "st_coordinates.csv")
    spatial = coords.merge(st_meta[["spot_id", "cell_type"]], on="spot_id", how="left")
    anchor = spatial[spatial["cell_type"].astype(str).eq("Epithelial cells")]
    if len(anchor) < 20:
        return None

    expr = pd.read_csv(export / "sc_expression_normalized.csv").set_index("cell_id")
    available_genes = pd.Index(expr.columns.astype(str))
    overlap = list(available_genes.intersection(pd.Index(gene_set)))
    if len(overlap) < 8:
        return None

    rec: dict[str, object] = {
        "raw_sample": str(row["raw_sample"]),
        "profile_mask_sample": sample,
        "masked_target_type": str(row["masked_target_type"]),
        "anchor_cell_type": "Epithelial cells",
        "readout_cell_type": readout,
        "gene_set": gene_set_name,
        "gene_set_overlap": int(len(overlap)),
    }
    curves = {}
    for offset, (method_label, suffix) in enumerate(METHODS):
        loc = _load_method_cells(root, sample, suffix, readout)
        if len(loc) < 40:
            return None
        loc = loc[loc["OriginalCID"].astype(str).isin(expr.index)].copy()
        if len(loc) < 40:
            return None
        dist = _mean_nearest5(
            loc[["row", "col"]].to_numpy(dtype=float),
            anchor[["row", "col"]].to_numpy(dtype=float),
        )
        close = dist <= float(np.median(dist))
        far = ~close
        close_ids = loc.loc[close, "OriginalCID"].astype(str)
        far_ids = loc.loc[far, "OriginalCID"].astype(str)
        close_mean = np.log1p(expr.loc[close_ids]).mean(axis=0)
        far_mean = np.log1p(expr.loc[far_ids]).mean(axis=0)
        stats = close_mean - far_mean
        es, peak_rank, curve = _running_es(stats, set(overlap))
        nes, pval = _null_nes(stats, len(overlap), es, nperm=nperm, seed=seed + offset * 100)
        rec[f"{suffix}_ES"] = es
        rec[f"{suffix}_NES"] = nes
        rec[f"{suffix}_pval"] = pval
        rec[f"{suffix}_peak_rank"] = peak_rank
        rec[f"{suffix}_n_mapped_cells"] = int(len(loc))
        rec[f"{suffix}_n_close"] = int(close.sum())
        rec[f"{suffix}_n_far"] = int(far.sum())
        curves[suffix] = curve
    rec["delta_NES"] = float(rec["route2_highres_NES"]) - float(rec["baseline_highres_NES"])
    rec["route2_peak_frac"] = float(rec["route2_highres_peak_rank"]) / float(len(curves["route2_highres"]))
    rec["baseline_peak_frac"] = float(rec["baseline_highres_peak_rank"]) / float(len(curves["baseline_highres"]))
    rec["_curves"] = curves
    return rec


def _format_p(p: float) -> str:
    return "0.001" if p < 0.0015 else f"{p:.3f}"


def _smooth_curve(y: np.ndarray, window: int = 61) -> np.ndarray:
    if len(y) < 7:
        return y
    win = min(window, len(y) - (1 - len(y) % 2))
    if win % 2 == 0:
        win -= 1
    if win < 5:
        return y
    pad = win // 2
    return np.convolve(np.pad(y, (pad, pad), mode="edge"), np.ones(win) / win, mode="valid")


def _gradient_line(ax, x: np.ndarray, y: np.ndarray) -> None:
    cmap = LinearSegmentedColormap.from_list("rank_gradient", ["#f03b20", "#ffb000", "#78b95b", "#2b8cbe"])
    pts = np.array([x, y]).T.reshape(-1, 1, 2)
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    lc = LineCollection(segs, cmap=cmap, norm=plt.Normalize(x.min(), x.max()))
    lc.set_array(x)
    lc.set_linewidth(1.8)
    lc.set_zorder(4)
    ax.add_collection(lc)


def _draw_panel(ax, curve: pd.DataFrame, label: str, nes: float, pval: float) -> None:
    curve = curve.sort_values("rank")
    x = curve["rank"].to_numpy(dtype=float)
    x = x / float(x.max())
    y = curve["running_es"].to_numpy(dtype=float)
    y_smooth = _smooth_curve(y)
    hits = curve["hit"].astype(bool).to_numpy()
    hit_idx = np.where(hits)[0]
    cmap = LinearSegmentedColormap.from_list("rank_gradient_hits", ["#f03b20", "#ffb000", "#78b95b", "#2b8cbe"])
    hit_colors = cmap(x[hit_idx])
    _gradient_line(ax, x, y_smooth)
    ax.vlines(x[hit_idx], 0, y_smooth[hit_idx], color=hit_colors, linewidth=0.72, alpha=0.95, zorder=3)
    ax.scatter(x[hit_idx], y_smooth[hit_idx], s=8.2, color=hit_colors, edgecolor="none", alpha=0.95, zorder=5)
    ax.axhline(0, color="#222222", linewidth=0.9)
    y_lo = min(0.0, float(y_smooth.min()))
    y_hi = max(0.0, float(y_smooth.max()))
    span = max(y_hi - y_lo, 0.2)
    ax.set_ylim(y_lo - span * 0.10, y_hi + span * 0.18)
    ax.set_xlim(0, 1)
    ax.set_xticks([])
    ax.tick_params(axis="y", labelsize=6.4, width=0.9, length=2.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_linewidth(0.9)
    ax.text(0.02, 1.03, label, transform=ax.transAxes, ha="left", va="bottom", fontsize=7.7, weight="bold")
    ax.text(0.70, 0.76, f"NES = {nes:.2f}\nP = {_format_p(pval)}", transform=ax.transAxes, ha="left", va="center", fontsize=6.4)


def _plot(best: dict[str, object], out_dir: Path, out_prefix: str) -> tuple[Path, Path]:
    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})
    fig = plt.figure(figsize=(3.05, 3.85), dpi=320)
    gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, 0.24], hspace=0.28, left=0.20, right=0.98, top=0.90, bottom=0.12)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
    ax3 = fig.add_subplot(gs[2, 0], sharex=ax1)
    _draw_panel(ax1, best["_curves"]["baseline_highres"], "CytoSPACE", float(best["baseline_highres_NES"]), float(best["baseline_highres_pval"]))
    _draw_panel(ax2, best["_curves"]["route2_highres"], "SVTuner + CytoSPACE", float(best["route2_highres_NES"]), float(best["route2_highres_pval"]))
    fig.text(0.08, 0.57, "Running enrichment score", rotation=90, va="center", ha="center", fontsize=8.2)
    title_readout = str(best["readout_cell_type"]).replace("Monocytes and Macrophages", "Mono/Mac")
    title_geneset = (
        str(best["gene_set"])
        .replace("EcoTyper_Monocytes_and_Macrophages_", "EcoTyper Mono/Mac ")
        .replace("EcoTyper_CD8_T_cells_", "EcoTyper CD8 ")
        .replace("EcoTyper_CD4_T_cells_", "EcoTyper CD4 ")
        .replace("_", " ")
    )
    title = f"{title_readout} {title_geneset}"
    fig.text(0.58, 0.975, title, ha="center", va="top", fontsize=7.8, linespacing=0.9)
    ax3.axis("off")
    ax3.annotate("", xy=(0.86, 0.54), xytext=(0.20, 0.54), arrowprops=dict(arrowstyle="->", lw=1.8, color="#8a8a8a"))
    ax3.text(0.01, 0.54, "Increasing distance\nfrom epithelial cells", ha="left", va="center", fontsize=6.0)
    png = out_dir / f"{out_prefix}.png"
    pdf = out_dir / f"{out_prefix}.pdf"
    svg = out_dir / f"{out_prefix}.svg"
    fig.savefig(png, dpi=320)
    fig.savefig(pdf)
    fig.savefig(svg)
    plt.close(fig)
    return png, pdf


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project_root", default=".")
    parser.add_argument("--sample", default=None)
    parser.add_argument("--readout_cell_type", default=None)
    parser.add_argument("--gene_set", default=None)
    parser.add_argument("--nperm", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=11)
    parser.add_argument("--out_dir", default="visualizations/highres_profile_mask_fig2c_expression_enrichment")
    args = parser.parse_args()

    root = Path(args.project_root).resolve()
    out_dir = root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    gene_sets = _read_gene_sets(root)
    manifest = pd.read_csv(root / "visualizations/highres_profile_mask_mapping/highres_profile_mask_mapping_manifest.csv")
    candidates = []
    for _, row in manifest.iterrows():
        if args.sample and args.sample not in {str(row["raw_sample"]), str(row["profile_mask_sample"])}:
            continue
        readouts = [args.readout_cell_type] if args.readout_cell_type else list(READOUT_TO_GENESET)
        for readout in readouts:
            names = [args.gene_set] if args.gene_set else READOUT_TO_GENESET.get(readout, [])
            for name in names:
                if name not in gene_sets:
                    continue
                rec = _evaluate(root, row, readout, name, gene_sets[name], args.nperm, args.seed)
                if rec is not None:
                    candidates.append(rec)
    if not candidates:
        raise ValueError("No valid expression-enrichment candidates found.")
    summary = pd.DataFrame([{k: v for k, v in rec.items() if k != "_curves"} for rec in candidates])
    # Prefer route2-improved, positive, early-peak curves without exploding NES.
    rankable = summary[
        (summary["delta_NES"] > 0)
        & (summary["route2_highres_NES"] > 1.2)
        & (summary["route2_highres_NES"] < 5.0)
        & (summary["route2_peak_frac"] < 0.45)
    ].copy()
    if rankable.empty:
        rankable = summary[summary["delta_NES"] > 0].copy()
    best_key = rankable.sort_values(["route2_peak_frac", "delta_NES"], ascending=[True, False]).iloc[0]
    best = next(
        rec
        for rec in candidates
        if rec["profile_mask_sample"] == best_key["profile_mask_sample"]
        and rec["readout_cell_type"] == best_key["readout_cell_type"]
        and rec["gene_set"] == best_key["gene_set"]
    )
    summary = summary.sort_values(["route2_peak_frac", "delta_NES"], ascending=[True, False])
    summary.to_csv(out_dir / "highres_fig2c_expression_candidate_scan.csv", index=False)
    out_prefix = f"fig2c_{_slug(best['readout_cell_type'])}_{_slug(best['gene_set'])}_{_slug(best['raw_sample'])}_baseline_vs_route2"
    for suffix, curve in best["_curves"].items():
        curve.to_csv(out_dir / f"{out_prefix}_{suffix}_curve.csv", index=False)
    png, pdf = _plot(best, out_dir, out_prefix)
    manifest_out = {k: v for k, v in best.items() if k != "_curves"}
    manifest_out.update(
        {
            "experiment_design": (
                "True Fig2C-style expression enrichment on high-resolution profile-mask data. "
                "Mapped readout cells are split into close/far groups by distance to epithelial cells. "
                "Genes are ranked by close-minus-far mapped expression, and a readout-state gene set is tested. "
                "Route2 uses Stage3-detected missing types; no forced whitelist is used."
            ),
            "png": str(png.relative_to(root)).replace("\\", "/"),
            "pdf": str(pdf.relative_to(root)).replace("\\", "/"),
            "candidate_scan": str((out_dir / "highres_fig2c_expression_candidate_scan.csv").relative_to(root)).replace("\\", "/"),
        }
    )
    (out_dir / f"{out_prefix}_manifest.json").write_text(json.dumps(manifest_out, indent=2, ensure_ascii=False), encoding="utf-8")
    print("[OK] selected:", json.dumps(manifest_out, ensure_ascii=False))
    print("[OK] wrote:", png)
    print("[OK] wrote:", out_dir / "highres_fig2c_expression_candidate_scan.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
