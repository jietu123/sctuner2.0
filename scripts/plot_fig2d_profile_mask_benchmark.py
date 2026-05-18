from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


METHODS = ["CytoSPACE", "SVTuner + CytoSPACE", "Tangram", "CellTrek"]
METHOD_SOURCES = {
    "CytoSPACE": "baseline_profile300",
    "SVTuner + CytoSPACE": "route2_profile300",
    "Tangram": "enrichment_tangram_marker",
    "CellTrek": "enrichment_celltrek",
}
COLORS = {
    "CytoSPACE": "#ef6a5b",
    "SVTuner + CytoSPACE": "#2f9a8f",
    "Tangram": "#d9d9d9",
    "CellTrek": "#d9d9d9",
}


def _read_manifest(path: Path) -> list[dict]:
    if not path.exists():
        raise FileNotFoundError(path)
    return json.loads(path.read_text(encoding="utf-8-sig"))


def collect(root: Path, require_stage3_detected: bool = True) -> tuple[pd.DataFrame, pd.DataFrame]:
    result_root = root / "result" / "cytospace_fig2d_profile_mask_benchmark"
    manifest = _read_manifest(result_root / "profile_mask_manifest.json")
    meta = pd.DataFrame(manifest)
    fallback_values = result_root / "fig2d_profile_mask_benchmark_source_values.csv"
    if fallback_values.exists() and not any((result_root / str(row["pair_id"])).exists() for row in manifest):
        return pd.read_csv(fallback_values), meta
    rows = []
    for entry in manifest:
        if require_stage3_detected and not bool(entry.get("target_detected_by_stage3")):
            continue
        pair = str(entry["pair_id"])
        for method, subdir in METHOD_SOURCES.items():
            path = result_root / pair / subdir / "fig2c_official_enrichment_summary.csv"
            if not path.exists():
                continue
            df = pd.read_csv(path)
            df = df[df["cell_type"].isin(["CD4 T cells", "CD8 T cells"])].copy()
            for _, row in df.iterrows():
                rows.append(
                    {
                        "pair_id": pair,
                        "profile_mask_sample": entry["profile_mask_sample"],
                        "masked_target_type": entry["masked_target_type"],
                        "stage3_missing_types": "|".join(entry.get("stage3_missing_types") or []),
                        "target_detected_by_stage3": bool(entry.get("target_detected_by_stage3")),
                        "cell_type": row["cell_type"],
                        "method": method,
                        "nes": float(row["NES"]),
                        "pval": float(row["pval"]),
                    }
                )
    if not rows:
        raise RuntimeError("No enrichment rows found for profile-mask benchmark.")
    return pd.DataFrame(rows), meta


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot non-forced profile-mask Fig.2d benchmark.")
    parser.add_argument("--project_root", default=".")
    parser.add_argument("--out_dir", default="visualizations/cytospace_fig2d_profile_mask_benchmark")
    parser.add_argument("--out_prefix", default="fig2d_profile_mask_benchmark")
    parser.add_argument("--include_stage3_failures", action="store_true")
    args = parser.parse_args()

    root = Path(args.project_root).resolve()
    out_dir = root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    df, meta = collect(root, require_stage3_detected=not args.include_stage3_failures)
    df.to_csv(out_dir / f"{args.out_prefix}_source_values.csv", index=False)
    meta.to_csv(out_dir / f"{args.out_prefix}_profile_mask_manifest.csv", index=False)

    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, ax = plt.subplots(figsize=(3.55, 3.25), dpi=300)
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

    rng = np.random.default_rng(23)
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

    flat = np.concatenate([v for v in values if len(v)])
    ymin = min(-2.2, float(np.nanmin(flat)) - 0.20)
    ymax = max(3.05, float(np.nanmax(flat)) + 0.55)
    ax.axhline(0, color="#555555", linewidth=0.75, linestyle=(0, (2.2, 1.8)))
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(0.45, len(METHODS) + 0.75)
    ax.set_ylabel("Normalized\nenrichment score", fontsize=7.2)
    ax.set_xticks(positions)
    ax.set_xticklabels(
        ["CytoSPACE", "SVTuner +\nCytoSPACE", "Tangram", "CellTrek"],
        rotation=45,
        ha="right",
        rotation_mode="anchor",
        fontsize=6.2,
    )
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
        "T cell exhaustion enrichment under\nStage3-detected profile-mask perturbation",
        fontsize=6.2,
        fontweight="bold",
        va="top",
    )
    n_pairs = df[["pair_id", "cell_type"]].drop_duplicates().shape[0]
    fig.text(0.190, 0.905, f"n = {n_pairs} CD4/CD8 profiles from auto-detected masked targets", fontsize=4.9, va="top")

    fig.subplots_adjust(left=0.18, right=0.78, top=0.82, bottom=0.30)
    png = out_dir / f"{args.out_prefix}.png"
    pdf = out_dir / f"{args.out_prefix}.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)
    print(f"[OK] wrote: {png}")
    print(f"[OK] wrote: {pdf}")
    print(f"[INFO] rows per method: {df.groupby('method').size().to_dict()}")
    print(f"[INFO] mean NES: {df.groupby('method')['nes'].mean().to_dict()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
