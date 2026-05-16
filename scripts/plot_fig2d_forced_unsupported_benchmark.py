from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon


SELECTED = {
    "melanoma_slide1": "macrophages",
    "melanoma_slide2": "non_malignant_unresolved",
    "brca_er_her2_fresh_frozen": "t_cells",
    "brca_her2_ffpe": "plasma_cells",
    "brca_tnbc_fresh_frozen": "t_cells",
    "crc_fresh_frozen": "fibroblasts",
}
METHODS = ["CytoSPACE", "SVTuner + CytoSPACE", "Tangram", "CellTrek"]
COLORS = {
    "CytoSPACE": "#ef6a5b",
    "SVTuner + CytoSPACE": "#2f9a8f",
    "Tangram": "#d9d9d9",
    "CellTrek": "#d9d9d9",
}


def _sig_label(p: float) -> str:
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "n.s."


def _load_summary(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path)
    return df[df["cell_type"].isin(["CD4 T cells", "CD8 T cells"])].copy()


def collect(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    result_root = root / "result" / "cytospace_fig2d_forced_unsupported"
    rows = []
    selected_rows = []
    for pair, forced_type in SELECTED.items():
        method_sources = [
            ("CytoSPACE", result_root / pair / "baseline_fast300" / "fig2c_official_enrichment_summary.csv", "none"),
            (
                "SVTuner + CytoSPACE",
                result_root / pair / f"route2_fast300_{forced_type}" / "fig2c_official_enrichment_summary.csv",
                forced_type,
            ),
            ("Tangram", result_root / pair / "enrichment_tangram_marker" / "fig2c_official_enrichment_summary.csv", "none"),
            ("CellTrek", result_root / pair / "enrichment_celltrek" / "fig2c_official_enrichment_summary.csv", "none"),
        ]
        selected_rows.append({"pair_id": pair, "forced_unsupported_type": forced_type})
        for method, path, forced in method_sources:
            df = _load_summary(path)
            for _, row in df.iterrows():
                rows.append(
                    {
                        "pair_id": pair,
                        "cell_type": row["cell_type"],
                        "method": method,
                        "nes": float(row["NES"]),
                        "pval": float(row["pval"]),
                        "forced_unsupported_type": forced,
                    }
                )
    return pd.DataFrame(rows), pd.DataFrame(selected_rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot Fig.2d-style forced unsupported perturbation benchmark.")
    parser.add_argument("--project_root", default=".")
    parser.add_argument("--out_dir", default="visualizations/cytospace_fig2d_forced_unsupported")
    parser.add_argument("--out_prefix", default="fig2d_forced_unsupported_benchmark")
    args = parser.parse_args()

    root = Path(args.project_root).resolve()
    out_dir = root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    df, selected = collect(root)
    df.to_csv(out_dir / f"{args.out_prefix}_source_values.csv", index=False)
    selected.to_csv(out_dir / f"{args.out_prefix}_selected_forced_types.csv", index=False)

    wide = df.pivot_table(index=["pair_id", "cell_type"], columns="method", values="nes", aggfunc="first")
    paired = wide[METHODS].dropna()
    pvalues = {}
    for method in METHODS[1:]:
        pvalues[method] = float(wilcoxon(paired[method], paired["CytoSPACE"], alternative="greater").pvalue)

    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, ax = plt.subplots(figsize=(3.45, 3.05), dpi=300)
    positions = np.arange(1, len(METHODS) + 1)
    values = [df.loc[df["method"].eq(m), "nes"].to_numpy(float) for m in METHODS]
    bp = ax.boxplot(
        values,
        positions=positions,
        widths=0.52,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#222222", "linewidth": 0.95},
        boxprops={"color": "#8a8a8a", "linewidth": 0.85},
        whiskerprops={"color": "#555555", "linewidth": 0.85},
        capprops={"color": "#555555", "linewidth": 0.85},
    )
    for patch, method in zip(bp["boxes"], METHODS):
        patch.set_facecolor(COLORS[method])
        patch.set_alpha(0.95)
        patch.set_zorder(1)
    for key in ["medians", "whiskers", "caps"]:
        for artist in bp[key]:
            artist.set_zorder(2)

    rng = np.random.default_rng(19)
    for i, method in enumerate(METHODS, start=1):
        y = df.loc[df["method"].eq(method), "nes"].to_numpy(float)
        ax.scatter(
            np.full(len(y), i) + rng.normal(0, 0.04, len(y)),
            y,
            s=8,
            color="#222222",
            alpha=0.82,
            linewidth=0,
            zorder=5,
        )

    ax.axhline(0, color="#555555", linewidth=0.75, linestyle=(0, (2.2, 1.8)))
    ymin = min(-2.2, float(np.nanmin(np.concatenate(values))) - 0.15)
    ymax = max(3.05, float(np.nanmax(np.concatenate(values))) + 0.55)
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(0.45, len(METHODS) + 0.75)
    ax.set_ylabel("Normalized\nenrichment score", fontsize=7.2)
    ax.set_xticks(positions)
    ax.set_xticklabels(["CytoSPACE", "SVTuner +\nCytoSPACE", "Tangram", "CellTrek"], rotation=45, ha="right", rotation_mode="anchor", fontsize=6.2)
    ax.tick_params(axis="y", labelsize=6.2, width=0.85, length=3)
    ax.tick_params(axis="x", width=0.85, length=2.5, pad=0)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_linewidth(0.9)
    ax.spines["bottom"].set_linewidth(0.9)

    ax.annotate(
        "",
        xy=(4.50, 1.8),
        xytext=(4.50, -1.8),
        arrowprops=dict(arrowstyle="<->", color="#8a8a8a", lw=1.1),
        annotation_clip=False,
    )
    ax.text(4.62, 1.45, "Enriched\nclose to\ntumor", ha="left", va="center", fontsize=5.7)
    ax.text(4.62, -1.50, "Enriched\nfar from\ntumor", ha="left", va="center", fontsize=5.7)
    fig.text(
        0.190,
        0.965,
        "T cell exhaustion enrichment under\nunsupported-type perturbation",
        fontsize=6.4,
        fontweight="bold",
        va="top",
    )

    fig.subplots_adjust(left=0.18, right=0.78, top=0.84, bottom=0.30)
    png = out_dir / f"{args.out_prefix}.png"
    pdf = out_dir / f"{args.out_prefix}.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)
    print(f"[OK] wrote: {png}")
    print(f"[OK] wrote: {pdf}")
    print(f"[INFO] n per method: {df.groupby('method').size().to_dict()}")
    print(f"[INFO] mean NES: {df.groupby('method')['nes'].mean().to_dict()}")
    print(f"[INFO] route2 > baseline: {int((paired['SVTuner + CytoSPACE'] > paired['CytoSPACE']).sum())}/{len(paired)}")
    print(f"[INFO] Wilcoxon greater pvals vs CytoSPACE={pvalues}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
