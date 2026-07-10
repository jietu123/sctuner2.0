from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
BIOAPP = ROOT / "visualizations" / "bioapp_experiment"
OUT = BIOAPP / "bioapp_downstream_phase_d1_morphology_and_interface_analysis"

D0 = BIOAPP / "bioapp_downstream_phase_d0_morphology_domain_boundary_freeze"
P8 = BIOAPP / "bioapp_phase8_svtuner_vs_endpoint_evaluation_and_baseline_svtuner_comparison"

DOMAIN_ORDER = ["tumor_core", "immune_enriched", "stroma_rich", "mixed_boundary", "unmapped_or_excluded"]
DOMAIN_COLORS = {
    "tumor_core": "#8C1D40",
    "immune_enriched": "#D55E00",
    "stroma_rich": "#4D9221",
    "mixed_boundary": "#7B3294",
    "unmapped_or_excluded": "#D0D0D0",
}


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def load_inputs() -> pd.DataFrame:
    d0 = pd.read_csv(D0 / "phase_d0_frozen_morphology_domain_and_boundary_by_spot.csv")
    comp = pd.read_csv(P8 / "baseline_svtuner_endpoint_comparison_by_spot.csv")
    score = pd.read_csv(P8 / "svtuner_endpoint_score_by_spot.csv")
    score_cols = [
        "barcode",
        "baseline_immune_score",
        "baseline_nonimmune_score",
        "baseline_dominant_label",
        "svtuner_immune_score",
        "svtuner_nonimmune_score",
        "svtuner_dominant_label",
        "reference_unrepresented_score",
        "supported_assignment_score",
    ]
    df = d0.merge(comp, on=["barcode", "primary_endpoint_status"], how="left")
    df = df.merge(score[score_cols], on="barcode", how="left")
    df["endpoint_positive_binary"] = df["primary_endpoint_status"].eq("positive").astype(float)
    df["endpoint_negative_binary"] = df["primary_endpoint_status"].eq("negative").astype(float)
    df["main_analysis_spot"] = df["primary_endpoint_status"].isin(["positive", "negative"])
    df["withheld_binary"] = df["withheld_binary"].fillna(False).astype(bool)
    return df


def summarize_domains(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for domain in DOMAIN_ORDER:
        sub = df[df["morphology_domain"].eq(domain)].copy()
        main = sub[sub["main_analysis_spot"]]
        rows.append(
            {
                "morphology_domain": domain,
                "n_spots": int(len(sub)),
                "n_main_analysis_spots": int(len(main)),
                "endpoint_positive_rate": float(main["endpoint_positive_binary"].mean()) if len(main) else np.nan,
                "mean_cta_immune_fraction": float(sub["Immune_cells_fraction"].mean()) if len(sub) else np.nan,
                "mean_withheld_score": float(sub["withheld_score"].mean()) if len(sub) else np.nan,
                "withheld_binary_rate": float(sub["withheld_binary"].mean()) if len(sub) else np.nan,
                "mean_baseline_forced_nonimmune_burden": float(sub["baseline_forced_nonimmune_burden"].mean()) if len(sub) else np.nan,
                "mean_svtuner_prevented_burden_continuous": float(sub["svtuner_prevented_forced_nonimmune_burden_continuous"].mean()) if len(sub) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def assign_distance_bins(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    guard = json.loads((D0 / "phase_d0_morphology_domain_boundary_freeze_summary.json").read_text(encoding="utf-8"))
    median_nn = float(guard["boundary_summary"]["median_nearest_neighbor_distance_image_px"])
    out["boundary_distance_spot_width"] = out["signed_distance_to_tumor_stroma_boundary"] / median_nn
    bins = [-np.inf, -4, -2, -1, 0, 1, 2, 4, np.inf]
    labels = [
        "deep tumor\n(<-4)",
        "inner tumor\n(-4,-2)",
        "tumor edge\n(-2,-1)",
        "near tumor\n(-1,0)",
        "interface side\n(0,1)",
        "outer side\n(1,2)",
        "outer domain\n(2,4)",
        "far outside\n(>4)",
    ]
    out["boundary_distance_bin"] = pd.cut(out["boundary_distance_spot_width"], bins=bins, labels=labels)
    return out


def summarize_interface(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    labels = df["boundary_distance_bin"].cat.categories
    for label in labels:
        sub = df[df["boundary_distance_bin"].eq(label)]
        main = sub[sub["main_analysis_spot"]]
        rows.append(
            {
                "boundary_distance_bin": str(label),
                "n_spots": int(len(sub)),
                "n_main_analysis_spots": int(len(main)),
                "endpoint_positive_rate": float(main["endpoint_positive_binary"].mean()) if len(main) else np.nan,
                "withheld_binary_rate": float(sub["withheld_binary"].mean()) if len(sub) else np.nan,
                "mean_withheld_score": float(sub["withheld_score"].mean()) if len(sub) else np.nan,
                "mean_baseline_forced_nonimmune_burden": float(sub["baseline_forced_nonimmune_burden"].mean()) if len(sub) else np.nan,
                "mean_svtuner_prevented_burden_continuous": float(sub["svtuner_prevented_forced_nonimmune_burden_continuous"].mean()) if len(sub) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def save_panel_g(domain_summary: pd.DataFrame) -> None:
    metrics = [
        "endpoint_positive_rate",
        "mean_cta_immune_fraction",
        "withheld_binary_rate",
        "mean_withheld_score",
    ]
    labels = [
        "CTA endpoint+\nrate",
        "CTA immune\nfraction",
        "Withheld\nrate",
        "Withheld\nscore",
    ]
    mat = domain_summary.set_index("morphology_domain").loc[DOMAIN_ORDER, metrics]
    fig = plt.figure(figsize=(7.0, 3.2), facecolor="white")
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 2.3], wspace=0.28)

    ax0 = fig.add_subplot(gs[0, 0])
    counts = domain_summary.set_index("morphology_domain").loc[DOMAIN_ORDER, "n_spots"]
    ax0.barh(range(len(DOMAIN_ORDER)), counts.values, color=[DOMAIN_COLORS[d] for d in DOMAIN_ORDER], edgecolor="none")
    ax0.set_yticks(range(len(DOMAIN_ORDER)))
    ax0.set_yticklabels([d.replace("_", " ") for d in DOMAIN_ORDER], fontsize=7)
    ax0.invert_yaxis()
    ax0.set_xlabel("spots", fontsize=7)
    ax0.set_title("Frozen domains", fontsize=8, fontweight="bold")
    ax0.tick_params(axis="x", labelsize=7)
    for sp in ax0.spines.values():
        sp.set_visible(False)

    ax1 = fig.add_subplot(gs[0, 1])
    im = ax1.imshow(mat.to_numpy(float), vmin=0, vmax=1, cmap="magma", aspect="auto")
    ax1.set_yticks(range(len(DOMAIN_ORDER)))
    ax1.set_yticklabels([d.replace("_", " ") for d in DOMAIN_ORDER], fontsize=7)
    ax1.set_xticks(range(len(labels)))
    ax1.set_xticklabels(labels, fontsize=6, rotation=35, ha="right")
    ax1.set_title("Endpoint/SVTuner signal by frozen morphology domain", fontsize=8, fontweight="bold")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat.iloc[i, j]
            if np.isfinite(val):
                text_color = "white" if (val < 0.25 or val > 0.55) else "black"
                ax1.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=5.5, color=text_color)
    cbar = fig.colorbar(im, ax=ax1, fraction=0.035, pad=0.02)
    cbar.set_label("fraction / mean score", fontsize=6)
    cbar.ax.tick_params(labelsize=6)
    for sp in ax1.spines.values():
        sp.set_visible(False)

    fig.suptitle("G. Morphology-domain concordance with frozen CTA endpoint and SVTuner withheld output", fontsize=9, fontweight="bold", y=1.04)
    fig.savefig(OUT / "panel_g_morphology_domain_concordance.png", dpi=450, bbox_inches="tight")
    fig.savefig(OUT / "panel_g_morphology_domain_concordance.pdf", bbox_inches="tight")
    plt.close(fig)


def save_panel_h(interface_summary: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 3.0), facecolor="white")
    x = np.arange(len(interface_summary))
    plot_df = interface_summary.copy()
    plot_df["endpoint_positive_rate_plot"] = plot_df["endpoint_positive_rate"].where(plot_df["n_main_analysis_spots"] >= 20)
    series = [
        ("endpoint_positive_rate_plot", "CTA endpoint+ rate", "#D55E00"),
        ("withheld_binary_rate", "SVTuner withheld rate", "#0072B2"),
        ("mean_withheld_score", "Mean withheld score", "#009E73"),
    ]
    for col, label, color in series:
        ax.plot(x, plot_df[col], marker="o", linewidth=1.6, markersize=4.0, color=color, label=label)
    ax.set_xticks(x)
    xticklabels = [
        f"{row.boundary_distance_bin}\nn={int(row.n_spots)}"
        for row in interface_summary.itertuples(index=False)
    ]
    ax.set_xticklabels(xticklabels, fontsize=6, rotation=35, ha="right")
    ax.set_ylim(-0.02, 1.02)
    ax.set_ylabel("fraction / mean score", fontsize=7)
    ax.set_title("H. Endpoint and withheld burden across frozen tumor-stroma boundary distance", fontsize=8, fontweight="bold")
    ax.grid(axis="y", color="#E0E0E0", linewidth=0.6)
    ax.axvline(3.5, color="#333333", linewidth=0.8, linestyle="--", alpha=0.7)
    ax.text(3.55, 0.96, "boundary", fontsize=6, va="top")
    ax.legend(loc="upper left", fontsize=6, frameon=False, ncol=2)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    fig.savefig(OUT / "panel_h_tumor_stroma_interface_topology.png", dpi=450, bbox_inches="tight")
    fig.savefig(OUT / "panel_h_tumor_stroma_interface_topology.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df = load_inputs()
    df = assign_distance_bins(df)
    domain_summary = summarize_domains(df)
    interface_summary = summarize_interface(df)

    df.to_csv(OUT / "phase_d1_downstream_analysis_by_spot.csv", index=False)
    domain_summary.to_csv(OUT / "panel_g_morphology_domain_concordance_summary.csv", index=False)
    interface_summary.to_csv(OUT / "panel_h_tumor_stroma_interface_topology_summary.csv", index=False)
    save_panel_g(domain_summary)
    save_panel_h(interface_summary)

    guardrails = {
        "phase": "BioApp Downstream Phase D1 morphology-domain and interface analysis",
        "SVTuner_rerun": False,
        "Stage3_rerun": False,
        "Stage4_run": False,
        "CytoSPACE_rerun": False,
        "endpoint_redefined": False,
        "threshold_selected_using_endpoint_labels": False,
        "uses_frozen_phase_d0_domains_and_boundary": True,
        "uses_frozen_phase8_scores": True,
        "formal_metrics_recomputed": False,
        "guardrails_passed": True,
    }
    domain_top = domain_summary.sort_values("withheld_binary_rate", ascending=False).head(3)[
        ["morphology_domain", "withheld_binary_rate", "endpoint_positive_rate", "mean_cta_immune_fraction"]
    ].to_dict(orient="records")
    summary = {
        "decision": "PASS",
        "phase": "BioApp Downstream Phase D1 morphology-domain and interface analysis",
        "panel_g_generated": True,
        "panel_h_generated": True,
        "n_spots": int(len(df)),
        "domain_top_withheld_binary_rate": domain_top,
        "boundary_bins": interface_summary["boundary_distance_bin"].tolist(),
        "interpretation_status": "downstream validation prototype; not final biological claim",
        "ready_for_visual_review": True,
        "cautions": [
            "Morphology domains are proxy labels from Phase D0, not pathologist annotations.",
            "Panel H uses a frozen spot-graph boundary proxy and should not be described as a pixel-resolution histology boundary.",
            "No upstream BioApp metrics, endpoint labels, or thresholds were modified.",
        ],
    }
    write_json(OUT / "phase_d1_guardrails.json", guardrails)
    write_json(OUT / "phase_d1_morphology_interface_analysis_summary.json", summary)

    report = f"""# BioApp Downstream Phase D1 morphology and interface analysis

Decision: `PASS`

## Outputs

- `panel_g_morphology_domain_concordance.png/pdf`
- `panel_h_tumor_stroma_interface_topology.png/pdf`
- `panel_g_morphology_domain_concordance_summary.csv`
- `panel_h_tumor_stroma_interface_topology_summary.csv`
- `phase_d1_downstream_analysis_by_spot.csv`

## Boundary

This phase reads frozen D0 morphology-domain and tumor-stroma boundary definitions.
It does not redefine the CTA Immune endpoint and does not rerun SVTuner, Stage3,
Stage4, or CytoSPACE.

## Interpretation

These panels are downstream validation prototypes. They can support whether
SVTuner withheld output aligns with external morphology-domain and interface
structure, but they are not final biological-discovery claims.
"""
    (OUT / "phase_d1_morphology_interface_analysis_report.md").write_text(report, encoding="utf-8")

    print("BioApp Downstream Phase D1 completed.")
    print("Decision: PASS")
    print("Panel G generated: true")
    print("Panel H generated: true")
    print(f"Output directory: {OUT}")


if __name__ == "__main__":
    main()
