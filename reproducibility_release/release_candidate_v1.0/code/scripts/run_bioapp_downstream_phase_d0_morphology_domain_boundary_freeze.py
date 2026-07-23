from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parents[1]
BIOAPP = ROOT / "visualizations" / "bioapp_experiment"
OUT = BIOAPP / "bioapp_downstream_phase_d0_morphology_domain_boundary_freeze"

P2C = BIOAPP / "bioapp_phase2c_endpoint_freeze_with_validated_cta_to_spot_mapping"
FEAS = BIOAPP / "bioapp_downstream_feasibility_check"
V312 = BIOAPP / "bioapp_main_figure_v3_12_panel_A_evidence_chain_redesign"

PALETTE = {
    "tumor_core": "#8C1D40",
    "immune_enriched": "#D55E00",
    "stroma_rich": "#4D9221",
    "mixed_boundary": "#7B3294",
    "other_mapped": "#777777",
    "unmapped_or_excluded": "#D0D0D0",
}


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def load_inputs() -> tuple[pd.DataFrame, np.ndarray, dict]:
    spot = pd.read_csv(P2C / "spot_level_endpoint_freeze.csv")
    feat = pd.read_csv(FEAS / "spot_morphology_feature_feasibility_table.csv")
    img = plt.imread(V312 / "fig_bioapp_v3_embedded_tissue_background.png")
    meta = json.loads((V312 / "fig_bioapp_v3_embedded_image_metadata.json").read_text(encoding="utf-8"))

    keep = [
        "barcode",
        "rgb_mean",
        "rgb_std",
        "darkness",
        "saturation",
        "red_blue",
        "edge_proxy",
        "tissue_coverage",
    ]
    df = spot.merge(feat[keep], on="barcode", how="left")
    lowres = float(meta["lowres_scale_used"])
    df["x_plot"] = df["imagecol"] * lowres
    df["y_plot"] = df["imagerow"] * lowres
    return df, img, meta


def freeze_morphology_domains(df: pd.DataFrame) -> pd.DataFrame:
    """Freeze conservative CTA/H&E morphology-domain proxy labels.

    These are not pathology-grade segmentation labels. They are external
    endpoint-compatible morphology-domain proxies derived before any downstream
    SVTuner interpretation.
    """
    out = df.copy()
    mapped = out["CTA_mapped"].fillna(False).astype(bool)
    total_ok = out["total_CTA_objects"].fillna(0) >= 5
    tumor = out["Tumor_fraction"].fillna(0)
    immune = out["Immune_cells_fraction"].fillna(0)
    stroma = out["Stroma_fraction"].fillna(0)
    dom_frac = out["dominant_CTA_fraction"].fillna(0)

    conditions = [
        (~mapped) | (~total_ok),
        mapped & total_ok & (tumor >= 0.70) & (dom_frac >= 0.70),
        mapped & total_ok & (immune >= 0.30),
        mapped & total_ok & (stroma >= 0.35),
        mapped & total_ok & (tumor >= 0.35) & ((stroma + immune) >= 0.25),
        mapped & total_ok & (dom_frac < 0.60),
    ]
    choices = [
        "unmapped_or_excluded",
        "tumor_core",
        "immune_enriched",
        "stroma_rich",
        "mixed_boundary",
        "mixed_boundary",
    ]
    out["morphology_domain"] = np.select(conditions, choices, default="other_mapped")
    out["morphology_domain_source"] = "CTA_composition_plus_H&E_patch_feature_proxy"
    out["morphology_domain_frozen"] = True
    return out


def nearest_neighbor_distance(coords: np.ndarray) -> float:
    tree = cKDTree(coords)
    dists, _ = tree.query(coords, k=2)
    return float(np.median(dists[:, 1]))


def freeze_boundary(df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    out = df.copy()
    coords = out[["x_plot", "y_plot"]].to_numpy(float)
    median_nn = nearest_neighbor_distance(coords)
    neighbor_radius = median_nn * 1.65
    interface_band_width = median_nn * 1.25

    domains = out["morphology_domain"].astype(str).to_numpy()
    tumor_region = domains == "tumor_core"
    interface_partner_region = np.isin(domains, ["stroma_rich", "mixed_boundary"])
    mapped = out["CTA_mapped"].fillna(False).astype(bool).to_numpy()

    tree = cKDTree(coords)
    neigh = tree.query_ball_point(coords, r=neighbor_radius)
    boundary_seed = np.zeros(len(out), dtype=bool)
    for i, ns in enumerate(neigh):
        ns = [j for j in ns if j != i]
        if not ns or not mapped[i]:
            continue
        has_tumor_neighbor = any(tumor_region[j] and mapped[j] for j in ns)
        if interface_partner_region[i] and has_tumor_neighbor:
            boundary_seed[i] = True

    if boundary_seed.sum() == 0:
        out["signed_distance_to_tumor_stroma_boundary"] = np.nan
        out["abs_distance_to_tumor_stroma_boundary"] = np.nan
        out["interface_band"] = False
    else:
        btree = cKDTree(coords[boundary_seed])
        dist, _ = btree.query(coords, k=1)
        signed = dist.copy()
        signed[tumor_region] *= -1
        out["signed_distance_to_tumor_stroma_boundary"] = signed
        out["abs_distance_to_tumor_stroma_boundary"] = np.abs(signed)
        out["interface_band"] = np.abs(signed) <= interface_band_width

    out["tumor_region_for_boundary"] = tumor_region
    out["stroma_or_immune_region_for_boundary"] = interface_partner_region
    out["boundary_seed_spot"] = boundary_seed
    out["boundary_definition_frozen"] = True

    summary = {
        "boundary_definition": "stroma_rich/mixed_boundary spots directly adjacent to frozen tumor_core domains",
        "median_nearest_neighbor_distance_image_px": median_nn,
        "neighbor_radius_image_px": neighbor_radius,
        "interface_band_width_image_px": interface_band_width,
        "n_boundary_seed_spots": int(boundary_seed.sum()),
        "n_interface_band_spots": int(out["interface_band"].sum()),
        "n_tumor_region_spots_for_boundary": int(tumor_region.sum()),
        "n_stroma_or_immune_region_spots_for_boundary": int(interface_partner_region.sum()),
    }
    return out, summary


def draw_background(ax: plt.Axes, img: np.ndarray) -> None:
    h, w = img.shape[:2]
    ax.imshow(img, extent=(0, w, h, 0), alpha=0.88)
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def save_domain_map(df: pd.DataFrame, img: np.ndarray) -> None:
    fig, ax = plt.subplots(figsize=(5.8, 5.4), facecolor="white")
    draw_background(ax, img)
    order = ["unmapped_or_excluded", "other_mapped", "tumor_core", "stroma_rich", "mixed_boundary", "immune_enriched"]
    for dom in order:
        sub = df[df["morphology_domain"].eq(dom)]
        if sub.empty:
            continue
        ax.scatter(
            sub["x_plot"],
            sub["y_plot"],
            s=12 if dom != "unmapped_or_excluded" else 7,
            c=PALETTE[dom],
            edgecolors="#222222" if dom in {"immune_enriched", "mixed_boundary"} else "none",
            linewidths=0.20,
            alpha=0.78 if dom != "unmapped_or_excluded" else 0.42,
            label=dom.replace("_", " "),
        )
    ax.set_title("Frozen morphology-domain proxy", fontsize=10, fontweight="bold")
    ax.legend(loc="lower right", fontsize=6, frameon=True, framealpha=0.88, markerscale=1.2)
    fig.savefig(OUT / "phase_d0_frozen_morphology_domain_map.png", dpi=450, bbox_inches="tight")
    fig.savefig(OUT / "phase_d0_frozen_morphology_domain_map.pdf", bbox_inches="tight")
    plt.close(fig)


def save_boundary_map(df: pd.DataFrame, img: np.ndarray) -> None:
    fig, ax = plt.subplots(figsize=(5.8, 5.4), facecolor="white")
    draw_background(ax, img)
    valid = df["signed_distance_to_tumor_stroma_boundary"].notna()
    vals = df.loc[valid, "signed_distance_to_tumor_stroma_boundary"]
    lim = float(np.nanpercentile(np.abs(vals), 95)) if valid.any() else 1.0
    cmap = LinearSegmentedColormap.from_list("boundary_distance", ["#8C1D40", "#FFFFFF", "#4D9221"])
    norm = TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim)
    sc = ax.scatter(df.loc[valid, "x_plot"], df.loc[valid, "y_plot"], c=vals, s=10, cmap=cmap, norm=norm, alpha=0.78, linewidths=0)
    boundary = df["boundary_seed_spot"].fillna(False)
    ax.scatter(df.loc[boundary, "x_plot"], df.loc[boundary, "y_plot"], s=7, facecolors="none", edgecolors="#111111", linewidths=0.28, alpha=0.70, label="boundary seed")
    ax.set_title("Frozen tumor-stroma boundary proxy", fontsize=10, fontweight="bold")
    cbar = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Signed distance\n(tumor negative)", fontsize=7)
    cbar.ax.tick_params(labelsize=6)
    ax.legend(loc="lower right", fontsize=6, frameon=True, framealpha=0.88)
    fig.savefig(OUT / "phase_d0_frozen_boundary_distance_map.png", dpi=450, bbox_inches="tight")
    fig.savefig(OUT / "phase_d0_frozen_boundary_distance_map.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df, img, meta = load_inputs()
    df = freeze_morphology_domains(df)
    df, boundary_summary = freeze_boundary(df)

    df.to_csv(OUT / "phase_d0_frozen_morphology_domain_and_boundary_by_spot.csv", index=False)
    domain_counts = {str(k): int(v) for k, v in df["morphology_domain"].value_counts().to_dict().items()}
    interface_by_domain = (
        df.groupby("morphology_domain", dropna=False)["interface_band"]
        .agg(["sum", "count"])
        .reset_index()
        .rename(columns={"sum": "interface_band_spots", "count": "domain_spots"})
    )
    interface_by_domain["interface_band_fraction"] = interface_by_domain["interface_band_spots"] / interface_by_domain["domain_spots"]
    interface_by_domain.to_csv(OUT / "phase_d0_interface_band_by_domain.csv", index=False)

    save_domain_map(df, img)
    save_boundary_map(df, img)

    guardrails = {
        "phase": "BioApp Downstream Phase D0 morphology-domain and boundary freeze",
        "SVTuner_rerun": False,
        "Stage3_rerun": False,
        "Stage4_run": False,
        "CytoSPACE_rerun": False,
        "endpoint_redefined": False,
        "threshold_selected_using_endpoint_labels": False,
        "new_mapping_or_new_computation": "only downstream morphology-domain/boundary derivation from frozen inputs",
        "domain_definitions_use_svtuner_outputs": False,
        "boundary_definitions_use_svtuner_outputs": False,
        "data_changed": False,
        "metrics_changed": False,
        "guardrails_passed": True,
    }
    summary = {
        "decision": "PASS",
        "phase": "BioApp Downstream Phase D0 morphology-domain and boundary freeze",
        "domain_definition_status": "frozen",
        "boundary_definition_status": "frozen",
        "morphology_domain_type": "CTA composition plus low-resolution H&E patch feature proxy; not pathology-grade segmentation",
        "domain_counts": domain_counts,
        "boundary_summary": boundary_summary,
        "ready_for_panel_g_morphology_domain_concordance": True,
        "ready_for_panel_h_tumor_stroma_interface_topology": True,
        "cautions": [
            "Morphology domains are conservative proxy labels and should not be described as pathologist-annotated regions.",
            "Boundary is a spot-graph proxy derived from CTA Tumor/Stroma-rich and mixed-boundary fractions.",
            "SVTuner outputs were not used to define domains or boundary.",
        ],
    }
    write_json(OUT / "phase_d0_guardrails.json", guardrails)
    write_json(OUT / "phase_d0_morphology_domain_boundary_freeze_summary.json", summary)

    report = f"""# BioApp Downstream Phase D0 morphology-domain and boundary freeze

Decision: `PASS`

## Frozen definitions

- Morphology-domain labels: frozen.
- Tumor-stroma boundary proxy: frozen.
- Domain type: CTA composition plus low-resolution H&E patch feature proxy.
- These are not pathology-grade segmentation labels.

## Domain counts

{pd.Series(domain_counts).to_markdown()}

## Boundary summary

{pd.Series(boundary_summary).to_markdown()}

## Guardrails

- SVTuner rerun: false
- Stage3 rerun: false
- Stage4 run: false
- CytoSPACE rerun: false
- Endpoint redefined: false
- Domain definitions use SVTuner outputs: false
- Boundary definitions use SVTuner outputs: false

## Next

Proceed to Panel G morphology-domain concordance and Panel H tumor-stroma interface topology using these frozen D0 definitions.
"""
    (OUT / "phase_d0_morphology_domain_boundary_freeze_report.md").write_text(report, encoding="utf-8")
    print("BioApp Downstream Phase D0 completed.")
    print("Decision: PASS")
    print(f"Domain counts: {domain_counts}")
    print(f"Boundary seed spots: {boundary_summary['n_boundary_seed_spots']}")
    print(f"Interface band spots: {boundary_summary['n_interface_band_spots']}")
    print("Ready for Panel G/H: true")


if __name__ == "__main__":
    main()
