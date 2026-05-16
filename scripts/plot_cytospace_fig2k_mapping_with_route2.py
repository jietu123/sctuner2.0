from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr


STATE_MARKERS = {
    "Temra": ["GZMB", "GZMH", "NKG7", "PRF1", "KLRG1", "EOMES"],
    "CREM+ Tm": ["CREM", "CCR7", "IL7R"],
    "TNF+": ["TNF"],
    "AREG+ Tm": ["AREG", "IL7R"],
    "CCL5+ Tm": ["CCL5", "IL7R"],
    "Tn": ["CCR7", "SELL", "IL7R", "TCF7"],
    "CCR6+ Th17": ["CCR6", "RORC", "IL17A"],
    "ADSL+ Tn": ["ADSL", "CCR7", "IL7R"],
    "CXCR5+ Tfh": ["CXCR5", "ICOS", "IL21"],
    "CAPG+CREM- Tm": ["CAPG", "IL7R"],
    "IL26+ Th17": ["IL26", "CCR6", "RORC"],
    "TIMP1+ Tm": ["TIMP1", "IL7R"],
    "GZMK+ Tem": ["GZMK", "CCL5"],
    "IL21+ Tfh": ["IL21", "CXCR5", "ICOS"],
    "NME1+CCR4+": ["NME1", "CCR4"],
    "CAPG+ Tm": ["CAPG", "IL7R"],
    "TNFRSF9- Treg": ["FOXP3", "IL2RA", "CTLA4"],
    "S1PR1+ Treg": ["S1PR1", "FOXP3", "IL2RA", "CTLA4"],
    "NME1+CCR4-": ["NME1"],
    "TNFRSF9+ Treg": ["TNFRSF9", "FOXP3", "IL2RA", "CTLA4"],
    "IFNG+ Th": ["IFNG"],
    "ISG+ Th": ["STAT1", "IFNAR1", "IFNGR1", "IFNGR2", "NLRC5", "TAP1"],
    "ISG+ Treg": ["STAT1", "IFNAR1", "IFNGR1", "IFNGR2", "NLRC5", "TAP1", "FOXP3", "IL2RA"],
}

METHODS = [
    ("stage4_cytospace_baseline", "CytoSPACE", "#ef6a5b"),
    ("stage4_cytospace_route2", "SVTuner + CytoSPACE", "#2f9a8f"),
]


def _load_table_s9(path: Path) -> pd.DataFrame:
    raw = pd.read_excel(path, sheet_name="Table S9", header=None)
    df = raw.iloc[8:, 0:4].copy()
    df.columns = ["known_rank", "state", "official_cytospace_nes", "official_merscope_nes"]
    df = df.dropna(subset=["known_rank", "state"]).copy()
    df["known_rank"] = pd.to_numeric(df["known_rank"], errors="coerce")
    df["official_cytospace_nes"] = pd.to_numeric(df["official_cytospace_nes"], errors="coerce")
    df["official_merscope_nes"] = pd.to_numeric(df["official_merscope_nes"], errors="coerce")
    df = df.dropna(subset=["known_rank"]).copy()
    df["known_rank"] = df["known_rank"].astype(int)
    df["state"] = df["state"].astype(str).str.strip()
    return df.sort_values("known_rank").reset_index(drop=True)


def _format_p(pval: float) -> str:
    if pval >= 0.001:
        return rf"$P$ = {pval:.3f}".rstrip("0").rstrip(".")
    mantissa, exponent = f"{pval:.1e}".split("e")
    return rf"$P$ = {mantissa} x 10$^{{{int(exponent)}}}$"


def _fit_line(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    slope, intercept = np.polyfit(x, y, deg=1)
    x_grid = np.linspace(1, 23, 100)
    return x_grid, slope * x_grid + intercept


def _compute_values(root: Path, sample: str, table_s9: pd.DataFrame) -> pd.DataFrame:
    export = root / "data" / "processed" / sample / "stage1_preprocess" / "exported"
    sc_expr = pd.read_csv(export / "sc_expression_normalized.csv", index_col=0)
    spot_meta = pd.read_csv(export / "merscope_spot_metadata.csv")
    spot_meta["spot_id"] = spot_meta["spot_id"].astype(str)
    spot_zone = spot_meta.set_index("spot_id")["tumor_zone"]
    rows: list[dict[str, object]] = []
    for stage4_dir, method, _ in METHODS:
        assigned_path = root / "result" / sample / stage4_dir / "cytospace_output" / "assigned_locations.csv"
        assigned = pd.read_csv(assigned_path)
        assigned["OriginalCID"] = assigned["OriginalCID"].astype(str)
        assigned["SpotID"] = assigned["SpotID"].astype(str)
        assigned = assigned[assigned["OriginalCID"].isin(sc_expr.index)].copy()
        assigned["tumor_zone"] = assigned["SpotID"].map(spot_zone)
        expr = sc_expr.loc[assigned["OriginalCID"].to_list()].copy()
        expr.index = assigned.index
        for _, row in table_s9.iterrows():
            genes = [g for g in STATE_MARKERS.get(row["state"], []) if g in expr.columns]
            if genes:
                score = expr[genes].astype("float32").mean(axis=1)
                high = score[assigned["tumor_zone"].eq("tumor_high")].mean()
                low = score[assigned["tumor_zone"].eq("tumor_low")].mean()
                enrichment = float(high - low)
            else:
                enrichment = np.nan
            rows.append(
                {
                    "method": method,
                    "stage4_dir": stage4_dir,
                    "known_rank": int(row["known_rank"]),
                    "state": row["state"],
                    "state_marker_genes_used": ";".join(genes),
                    "n_state_marker_genes_used": len(genes),
                    "tumor_high_minus_low_score": enrichment,
                    "n_assigned_cells": int(len(assigned)),
                    "n_tumor_high_assignments": int(assigned["tumor_zone"].eq("tumor_high").sum()),
                    "n_tumor_low_assignments": int(assigned["tumor_zone"].eq("tumor_low").sum()),
                }
            )
    out = pd.DataFrame(rows)
    out["predicted_rank"] = out.groupby("method")["tumor_high_minus_low_score"].rank(method="average", ascending=True)
    return out


def _panel(ax: plt.Axes, df: pd.DataFrame, method: str, color: str, state_colors: dict[str, object]) -> tuple[float, float]:
    sub = df[df["method"].eq(method)].dropna(subset=["predicted_rank"]).copy()
    x = sub["known_rank"].to_numpy(float)
    y = sub["predicted_rank"].to_numpy(float)
    r, pval = pearsonr(x, y)
    grid, fit = _fit_line(x, y)
    colors = [state_colors.get(state, "#111111") for state in sub["state"].astype(str)]
    ax.scatter(x, y, s=8.5, color=colors, edgecolor="none", zorder=3)
    ax.plot(grid, fit, color=color, linewidth=1.1, zorder=2)
    ax.set_title(method, fontsize=6.8, pad=4)
    ax.text(0.07, 0.92, rf"$r$ = {r:.2f}" + "\n" + _format_p(pval), transform=ax.transAxes, va="top", ha="left", fontsize=5.8)
    ax.set_xlim(0.5, 23.5)
    ax.set_ylim(0.5, 23.5)
    ax.set_aspect("equal", adjustable="box")
    order_x = df.drop_duplicates("known_rank").sort_values("known_rank")
    ax.set_xticks(order_x["known_rank"])
    ax.set_xticklabels(order_x["state"], rotation=90, fontsize=2.1)
    order_y = sub.sort_values("predicted_rank")
    ax.set_yticks(order_y["predicted_rank"])
    ax.set_yticklabels(order_y["state"], fontsize=2.1)
    ax.tick_params(axis="both", width=0.55, length=1.1, pad=0.5)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_linewidth(0.75)
    ax.spines["bottom"].set_linewidth(0.75)
    return float(r), float(pval)


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot Fig.2k-style CytoSPACE vs Route2 mapped result.")
    parser.add_argument("--project_root", default=".")
    parser.add_argument("--sample", default="cytospace_fig2k_breast_forced_unsupported")
    parser.add_argument("--source_xlsx", default="data/raw/cytospace_fig2c_melanoma/41587_2023_1697_MOESM3_ESM.xlsx")
    parser.add_argument("--out_dir", default="visualizations/cytospace_fig2k_tcell_states_route2")
    parser.add_argument("--out_prefix", default="fig2k_mapping_baseline_vs_route2")
    args = parser.parse_args()

    root = Path(args.project_root).resolve()
    out_dir = root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    table_s9 = _load_table_s9(root / args.source_xlsx)
    values = _compute_values(root, args.sample, table_s9)
    values.to_csv(out_dir / f"{args.out_prefix}_source_values.csv", index=False)

    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(4.15, 2.45), dpi=300, sharey=True)
    stats = {}
    state_order = table_s9.sort_values("known_rank")["state"].astype(str).tolist()
    cmap = plt.get_cmap("viridis")
    state_colors = {state: cmap(i) for state, i in zip(state_order, np.linspace(0.05, 0.95, len(state_order)))}
    for ax, (_, method, color) in zip(axes, METHODS):
        stats[method] = _panel(ax, values, method, color, state_colors)
    axes[0].set_ylabel("Predicted tumor enrichment rank", fontsize=6.8)
    fig.text(0.56, 0.16, "Known tumor enrichment rank\n(Zheng et al.)", ha="center", va="center", fontsize=6.6)
    fig.suptitle("Tumor/normal enrichment of CD4 T cell states", fontsize=7.8, y=0.965)

    overlay = fig.add_axes([0, 0, 1, 1], frameon=False)
    overlay.set_axis_off()
    overlay.annotate("", xy=(0.087, 0.66), xytext=(0.087, 0.43), xycoords="figure fraction",
                     arrowprops=dict(arrowstyle="<->", color="#8a8a8a", lw=0.9), annotation_clip=False)
    fig.text(0.052, 0.68, "Higher in\ntumor", ha="center", va="bottom", fontsize=5.4)
    fig.text(0.052, 0.39, "Higher in\nnormal-like\nregions", ha="center", va="top", fontsize=5.4)
    overlay.annotate("", xy=(0.78, 0.105), xytext=(0.34, 0.105), xycoords="figure fraction",
                     arrowprops=dict(arrowstyle="<->", color="#8a8a8a", lw=0.9), annotation_clip=False)
    fig.text(0.28, 0.075, "Normal-like\nregions", ha="center", va="top", fontsize=5.4)
    fig.text(0.82, 0.096, "Tumor", ha="center", va="center", fontsize=5.4)

    fig.subplots_adjust(left=0.22, right=0.98, top=0.78, bottom=0.36, wspace=0.42)
    png = out_dir / f"{args.out_prefix}.png"
    pdf = out_dir / f"{args.out_prefix}.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)

    manifest = {
        "sample": args.sample,
        "figure": str(png),
        "source_values": str(out_dir / f"{args.out_prefix}_source_values.csv"),
        "stats": {method: {"r": r, "p": p} for method, (r, p) in stats.items()},
    }
    (out_dir / f"{args.out_prefix}_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[OK] wrote: {png}")
    print(f"[OK] wrote: {pdf}")
    for method, (r, p) in stats.items():
        n = int(values[values["method"].eq(method)]["predicted_rank"].notna().sum())
        print(f"[INFO] {method}: n={n} r={r:.4f} p={p:.4g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
