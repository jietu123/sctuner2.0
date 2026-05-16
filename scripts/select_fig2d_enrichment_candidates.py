from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.metrics import silhouette_score


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SIM_ROOT = PROJECT_ROOT / "data" / "sim"
OUT_DIR = PROJECT_ROOT / "result" / "fig2d_enrichment_candidates"

K_NEIGHBORS = 8
MIN_POS_SPOTS = 50
MAX_POS_FRAC = 0.60
TARGET_POS_FRAC = 0.20


def _simulation_dirs() -> list[Path]:
    dirs: list[Path] = []
    for group_dir in sorted(SIM_ROOT.glob("*")):
        if not group_dir.is_dir():
            continue
        for sample_dir in sorted(group_dir.glob("*")):
            if not sample_dir.is_dir():
                continue
            if sample_dir.name.endswith("_scnoise10"):
                continue
            if (sample_dir / "sim_info.json").exists():
                dirs.append(sample_dir)
    return dirs


def _load_coords(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=None, engine="python")
    cols = {c.lower(): c for c in df.columns}
    spot_col = cols.get("spot_id") or cols.get("spotid")
    row_col = cols.get("row")
    col_col = cols.get("col")
    if not (spot_col and row_col and col_col):
        raise ValueError(f"Unrecognized coordinate columns in {path}")
    out = df[[spot_col, row_col, col_col]].copy()
    out.columns = ["spot_id", "row", "col"]
    return out


def _load_truth_fraction(sample_dir: Path) -> pd.DataFrame:
    preferred = sample_dir / "sim_truth_spot_type_fraction_from_cells.csv"
    fallback = sample_dir / "sim_truth_spot_type_fraction.csv"
    path = preferred if preferred.exists() else fallback
    df = pd.read_csv(path)
    if "spot_id" not in df.columns:
        raise ValueError(f"spot_id missing in {path}")
    return df


def _load_dominant(sample_dir: Path) -> pd.DataFrame:
    path = sample_dir / "sim_truth_spot_dominant_type.csv"
    df = pd.read_csv(path)
    if "spot_id" not in df.columns or "dominant_type" not in df.columns:
        raise ValueError(f"dominant_type table malformed: {path}")
    return df


def _neighbor_purity(coords: np.ndarray, mask: np.ndarray, k: int = K_NEIGHBORS) -> float:
    pos_idx = np.flatnonzero(mask)
    if len(pos_idx) <= 1:
        return float("nan")
    tree = cKDTree(coords)
    query_k = min(k + 1, len(coords))
    _, nn = tree.query(coords[pos_idx], k=query_k)
    if nn.ndim == 1:
        nn = nn[:, None]
    nn = nn[:, 1:]
    return float(mask[nn].mean())


def _balance_score(pos_frac: float, target: float = TARGET_POS_FRAC) -> float:
    if pos_frac <= 0 or pos_frac >= 1:
        return 0.0
    return float(max(0.0, 1.0 - abs(pos_frac - target) / target))


def _scenario_penalty(missing_type: object) -> float:
    return 0.15 if pd.notna(missing_type) and str(missing_type).strip() not in {"", "None", "null"} else 0.0


def _evaluate_sample(sample_dir: Path) -> list[dict[str, object]]:
    info = json.loads((sample_dir / "sim_info.json").read_text(encoding="utf-8"))
    coords = _load_coords(sample_dir / "brca_STdata_coordinates.txt")
    truth = _load_truth_fraction(sample_dir)
    dominant = _load_dominant(sample_dir)

    merged = coords.merge(truth, on="spot_id", how="inner").merge(dominant[["spot_id", "dominant_type", "assigned_cells"]], on="spot_id", how="inner")
    type_cols = [c for c in truth.columns if c != "spot_id"]
    xy = merged[["row", "col"]].to_numpy(dtype=float)

    rows: list[dict[str, object]] = []
    for cell_type in type_cols:
        frac = merged[cell_type].to_numpy(dtype=float)
        pos_mask = (merged["dominant_type"].to_numpy(dtype=object) == cell_type) & (frac > 0)
        n_pos = int(pos_mask.sum())
        pos_frac = n_pos / max(len(merged), 1)
        if n_pos < MIN_POS_SPOTS or pos_frac >= MAX_POS_FRAC:
            continue

        dominant_conf = float(frac[pos_mask].mean())
        contrast = float(frac[pos_mask].mean() - frac[~pos_mask].mean())
        purity = _neighbor_purity(xy, pos_mask, k=K_NEIGHBORS)
        sil = float(silhouette_score(xy, pos_mask.astype(int))) if 1 < n_pos < (len(merged) - 1) else float("nan")
        sil01 = 0.5 if not np.isfinite(sil) else float((sil + 1.0) / 2.0)
        balance = _balance_score(pos_frac)
        penalty = _scenario_penalty(info.get("missing_type"))

        composite = (
            0.35 * purity
            + 0.25 * sil01
            + 0.20 * contrast
            + 0.10 * dominant_conf
            + 0.10 * balance
            - penalty
        )

        rows.append(
            {
                "group": sample_dir.parent.name,
                "sample": sample_dir.name,
                "source_sample": info.get("source_sample"),
                "simulation_type": info.get("simulation_type"),
                "missing_type": info.get("missing_type"),
                "cell_type": cell_type,
                "n_spots": int(len(merged)),
                "n_pos_spots": n_pos,
                "pos_spot_fraction": pos_frac,
                "dominant_confidence": dominant_conf,
                "fraction_contrast": contrast,
                "neighbor_purity_k8": purity,
                "silhouette": sil,
                "silhouette_01": sil01,
                "balance_score": balance,
                "scenario_penalty": penalty,
                "composite_score": float(composite),
            }
        )
    return rows


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    all_rows: list[dict[str, object]] = []
    for sample_dir in _simulation_dirs():
        try:
            all_rows.extend(_evaluate_sample(sample_dir))
        except Exception as exc:
            all_rows.append(
                {
                    "group": sample_dir.parent.name,
                    "sample": sample_dir.name,
                    "cell_type": "__ERROR__",
                    "composite_score": float("nan"),
                    "error": str(exc),
                }
            )

    long_df = pd.DataFrame(all_rows)
    candidate_df = long_df[long_df["cell_type"] != "__ERROR__"].copy()
    candidate_df = candidate_df.sort_values(["composite_score", "neighbor_purity_k8", "fraction_contrast"], ascending=[False, False, False])

    per_sample = (
        candidate_df.groupby(["group", "sample"], as_index=False)
        .first()[[
            "group",
            "sample",
            "source_sample",
            "simulation_type",
            "missing_type",
            "cell_type",
            "n_pos_spots",
            "pos_spot_fraction",
            "dominant_confidence",
            "fraction_contrast",
            "neighbor_purity_k8",
            "silhouette",
            "balance_score",
            "scenario_penalty",
            "composite_score",
        ]]
        .rename(columns={"cell_type": "recommended_cell_type"})
        .sort_values("composite_score", ascending=False)
    )

    recommended = per_sample.head(5).to_dict(orient="records")
    summary = {
        "selection_goal": "Identify simulation scenarios with the cleanest known spatial enrichment structure for a paper-style Fig.2d NES benchmark.",
        "n_candidate_rows": int(len(candidate_df)),
        "n_samples_ranked": int(len(per_sample)),
        "top_recommendations": recommended,
    }

    candidate_df.to_csv(OUT_DIR / "scenario_celltype_candidates.csv", index=False)
    per_sample.to_csv(OUT_DIR / "scenario_recommendations.csv", index=False)
    (OUT_DIR / "selection_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(f"[OK] wrote: {OUT_DIR / 'scenario_celltype_candidates.csv'}")
    print(f"[OK] wrote: {OUT_DIR / 'scenario_recommendations.csv'}")
    print(f"[OK] wrote: {OUT_DIR / 'selection_summary.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
