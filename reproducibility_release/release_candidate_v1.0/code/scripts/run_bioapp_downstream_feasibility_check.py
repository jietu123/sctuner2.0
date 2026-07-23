from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
BIOAPP = ROOT / "visualizations" / "bioapp_experiment"
OUT = BIOAPP / "bioapp_downstream_feasibility_check"

P2C = BIOAPP / "bioapp_phase2c_endpoint_freeze_with_validated_cta_to_spot_mapping"
P8 = BIOAPP / "bioapp_phase8_svtuner_vs_endpoint_evaluation_and_baseline_svtuner_comparison"
V312 = BIOAPP / "bioapp_main_figure_v3_12_panel_A_evidence_chain_redesign"


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def cohens_d(a: pd.Series, b: pd.Series) -> float | None:
    a = pd.Series(a).dropna().astype(float)
    b = pd.Series(b).dropna().astype(float)
    if len(a) < 2 or len(b) < 2:
        return None
    pooled = math.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return float((a.mean() - b.mean()) / pooled) if pooled else None


def spearman(a: pd.Series | np.ndarray, b: pd.Series | np.ndarray) -> float | None:
    s = pd.DataFrame({"a": a, "b": b}).dropna()
    if len(s) < 10:
        return None
    return float(s["a"].corr(s["b"], method="spearman"))


def extract_patch_features(df: pd.DataFrame, img: np.ndarray, radius: int = 9) -> pd.DataFrame:
    h, w = img.shape[:2]
    rows: list[dict] = []
    for _, r in df.iterrows():
        x = int(round(float(r["x_plot"])))
        y = int(round(float(r["y_plot"])))
        x0, x1 = max(0, x - radius), min(w, x + radius + 1)
        y0, y1 = max(0, y - radius), min(h, y + radius + 1)
        patch = img[y0:y1, x0:x1, :]
        vals: dict[str, float | str | None]
        if patch.size == 0:
            vals = {
                "rgb_mean": None,
                "rgb_std": None,
                "darkness": None,
                "saturation": None,
                "red_blue": None,
                "edge_proxy": None,
                "tissue_coverage": None,
            }
        else:
            gray = patch.mean(axis=2)
            mx = patch.max(axis=2)
            mn = patch.min(axis=2)
            gy = np.abs(np.diff(gray, axis=0)).mean() if gray.shape[0] > 1 else 0.0
            gx = np.abs(np.diff(gray, axis=1)).mean() if gray.shape[1] > 1 else 0.0
            vals = {
                "rgb_mean": float(patch.mean()),
                "rgb_std": float(patch.std()),
                "darkness": float((1.0 - gray).mean()),
                "saturation": float((mx - mn).mean()),
                "red_blue": float((patch[:, :, 0] - patch[:, :, 2]).mean()),
                "edge_proxy": float(gx + gy),
                "tissue_coverage": float((gray < 0.92).mean()),
            }
        vals["barcode"] = r["barcode"]
        rows.append(vals)
    return pd.DataFrame(rows)


def interface_proxy(df: pd.DataFrame) -> dict:
    coords = df[["x_plot", "y_plot"]].to_numpy(float)
    tumor_pos = df["Tumor_endpoint_status"].eq("positive").to_numpy()
    stroma_signal = (df["Stroma_fraction"].fillna(0).to_numpy() >= 0.20) | df["Stroma_endpoint_status"].eq("positive").to_numpy()
    withheld = df["withheld_binary"].eq(True).to_numpy()
    not_withheld = df["withheld_binary"].eq(False).to_numpy()
    endpoint_pos = df["primary_endpoint_status"].eq("positive").to_numpy()
    try:
        from scipy.spatial import cKDTree

        tree_all = cKDTree(coords)
        neigh = tree_all.query_ball_point(coords, r=18.0)
        boundary = []
        for i, ns in enumerate(neigh):
            if not tumor_pos[i]:
                boundary.append(False)
                continue
            ns = [j for j in ns if j != i]
            has_stroma = any(stroma_signal[j] for j in ns)
            has_non_tumor = any(not tumor_pos[j] for j in ns)
            boundary.append(has_stroma or has_non_tumor)
        boundary_arr = np.array(boundary, dtype=bool)
        if boundary_arr.sum() == 0:
            return {"boundary_proxy_available": False, "n_boundary_proxy_spots": 0}
        btree = cKDTree(coords[boundary_arr])
        dist, _ = btree.query(coords, k=1)
        signed = dist.copy()
        signed[tumor_pos] *= -1
        df["signed_distance_to_tumor_boundary_proxy"] = signed
        return {
            "boundary_proxy_available": True,
            "n_boundary_proxy_spots": int(boundary_arr.sum()),
            "median_abs_distance_endpoint_positive": float(np.median(np.abs(signed[endpoint_pos]))) if endpoint_pos.any() else None,
            "median_abs_distance_withheld": float(np.median(np.abs(signed[withheld]))) if withheld.any() else None,
            "median_abs_distance_not_withheld": float(np.median(np.abs(signed[not_withheld]))) if not_withheld.any() else None,
            "spearman_abs_distance_withheld_score": spearman(np.abs(signed), df["withheld_score"]),
        }
    except Exception as exc:  # noqa: BLE001
        return {"boundary_proxy_available": False, "error": repr(exc)}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    spot = pd.read_csv(P2C / "spot_level_endpoint_freeze.csv")
    sv = pd.read_csv(P8 / "svtuner_endpoint_score_by_spot.csv")
    img_path = V312 / "fig_bioapp_v3_embedded_tissue_background.png"
    meta = read_json(V312 / "fig_bioapp_v3_embedded_image_metadata.json")
    img = np.asarray(Image.open(img_path).convert("RGB")).astype(float) / 255.0
    h, w = img.shape[:2]
    lowres = float(meta["lowres_scale_used"])

    spot["x_plot"] = spot["imagecol"] * lowres
    spot["y_plot"] = spot["imagerow"] * lowres
    df = spot.merge(
        sv[
            [
                "barcode",
                "withheld_score",
                "reference_unrepresented_score",
                "withheld_binary",
                "baseline_nonimmune_score",
                "baseline_dominant_label",
                "svtuner_dominant_label",
            ]
        ],
        on="barcode",
        how="left",
    )

    feat = extract_patch_features(df, img)
    merged = df.merge(feat, on="barcode", how="left")
    merged["withheld_binary"] = merged["withheld_binary"].astype(bool)

    counts = {
        "n_spots": int(len(merged)),
        "image_width": int(w),
        "image_height": int(h),
        "coords_within_image": bool(merged["x_plot"].between(0, w).all() and merged["y_plot"].between(0, h).all()),
        "primary_endpoint_counts": {str(k): int(v) for k, v in merged["primary_endpoint_status"].value_counts(dropna=False).to_dict().items()},
        "tumor_endpoint_counts": {str(k): int(v) for k, v in merged["Tumor_endpoint_status"].value_counts(dropna=False).to_dict().items()},
        "stroma_endpoint_counts": {str(k): int(v) for k, v in merged["Stroma_endpoint_status"].value_counts(dropna=False).to_dict().items()},
        "dominant_CTA_counts": {str(k): int(v) for k, v in merged["dominant_CTA_class"].value_counts(dropna=False).to_dict().items()},
        "withheld_binary_counts": {str(k): int(v) for k, v in merged["withheld_binary"].value_counts(dropna=False).to_dict().items()},
    }

    pos = merged["primary_endpoint_status"].eq("positive")
    neg = merged["primary_endpoint_status"].eq("negative")
    withheld = merged["withheld_binary"].eq(True)
    not_withheld = merged["withheld_binary"].eq(False)
    features = ["rgb_mean", "rgb_std", "darkness", "saturation", "red_blue", "edge_proxy", "tissue_coverage"]
    feature_tests = []
    for f in features:
        feature_tests.append(
            {
                "feature": f,
                "endpoint_pos_vs_neg_d": cohens_d(merged.loc[pos, f], merged.loc[neg, f]),
                "withheld_vs_not_d": cohens_d(merged.loc[withheld, f], merged.loc[not_withheld, f]),
                "spearman_withheld_score": spearman(merged[f], merged["withheld_score"]),
                "spearman_immune_fraction": spearman(merged[f], merged["Immune_cells_fraction"]),
                "spearman_tumor_fraction": spearman(merged[f], merged["Tumor_fraction"]),
                "spearman_stroma_fraction": spearman(merged[f], merged["Stroma_fraction"]),
            }
        )
    feature_tests_df = pd.DataFrame(feature_tests)
    interface_summary = interface_proxy(merged)

    merged.to_csv(OUT / "spot_morphology_feature_feasibility_table.csv", index=False)
    feature_tests_df.to_csv(OUT / "he_patch_feature_association_screen.csv", index=False)

    summary = {
        "decision": "FEASIBLE_WITH_CAUTION",
        "experiment_1_he_morphology_domain_concordance": "feasible as low-cost patch-feature and domain-proxy screen",
        "experiment_2_tumor_stroma_interface_topology": "feasible as CTA-fraction boundary proxy; stronger if H&E morphology domains are frozen",
        "inputs_present": {
            "embedded_tissue_image": img_path.exists(),
            "spot_coordinates": True,
            "cta_composition": True,
            "endpoint_status": True,
            "baseline_svtuner_scores": True,
        },
        "counts": counts,
        "feature_screen_top_abs_associations": feature_tests_df.assign(max_abs=feature_tests_df.drop(columns=["feature"]).abs().max(axis=1))
        .sort_values("max_abs", ascending=False)
        .head(5)
        .drop(columns=["max_abs"])
        .to_dict(orient="records"),
        "interface_summary": interface_summary,
        "cautions": [
            "Current embedded H&E is low-resolution 600x594; patch features are feasible but limited.",
            "CTA Tumor/Stroma/Immune fractions provide strong domain proxies, but pure H&E morphology-domain labels are not yet frozen.",
            "Tumor-stroma boundary can be approximated from CTA fractions and the spot graph; a pathology-grade H&E segmentation would be stronger.",
            "Downstream panels should be treated as validation/preparation until morphology-domain labels are explicitly frozen.",
        ],
    }
    (OUT / "bioapp_downstream_feasibility_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    report = [
        "# BioApp downstream feasibility check",
        "",
        f"Decision: `{summary['decision']}`",
        "",
        "## Available inputs",
        "",
    ]
    for k, v in summary["inputs_present"].items():
        report.append(f"- {k}: `{v}`")
    report += ["", "## Spot/domain counts", ""]
    for k, v in counts.items():
        report.append(f"- {k}: `{v}`")
    report += ["", "## H&E patch feature screen", "", feature_tests_df.to_markdown(index=False), "", "## Tumor-stroma interface proxy", ""]
    for k, v in interface_summary.items():
        report.append(f"- {k}: `{v}`")
    report += [
        "",
        "## Recommendation",
        "",
        "Proceed first with a lightweight downstream Phase D0/D1 design: freeze H&E/CTA morphology-domain labels, then implement Panel G morphology-domain concordance and Panel H tumor-stroma interface topology. Do not yet run all four downstream experiments.",
    ]
    (OUT / "bioapp_downstream_feasibility_report.md").write_text("\n".join(report), encoding="utf-8")

    print("BioApp downstream feasibility check completed.")
    print(f"Decision: {summary['decision']}")
    print(f"Experiment 1 feasible: {summary['experiment_1_he_morphology_domain_concordance']}")
    print(f"Experiment 2 feasible: {summary['experiment_2_tumor_stroma_interface_topology']}")
    print(f"Output directory: {OUT}")


if __name__ == "__main__":
    main()
