from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

EPS = 1e-12


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Stage5 real-data metrics summary (baseline vs route2)")
    p.add_argument("--sample", default="real_brca")
    p.add_argument("--project_root", default=".")
    p.add_argument("--baseline_dir", required=True, help=".../stage4_cytospace_baseline/cytospace_output")
    p.add_argument("--route2_dir", required=True, help=".../stage4_cytospace_route2/cytospace_output")
    p.add_argument("--out_dir", default=None, help="default: result/<sample>/stage5_real_summary")
    p.add_argument("--top_k", type=int, default=10)
    return p.parse_args()


def _clean_spot_id(x: object) -> str:
    return str(x).split()[0].split("\t")[0]


def _read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _load_by_spot(stage4_dir: Path) -> pd.DataFrame:
    candidates = [
        stage4_dir / "cell_type_assignments_by_spot.csv",
        stage4_dir / "cell_type_assignment_by_spot.csv",
    ]
    by_spot_path = None
    for c in candidates:
        if c.exists():
            by_spot_path = c
            break
    if by_spot_path is None:
        raise FileNotFoundError(f"missing by-spot csv under {stage4_dir}")

    df = pd.read_csv(by_spot_path, low_memory=False)
    if "spot_id" not in df.columns:
        first = df.columns[0]
        df = df.rename(columns={first: "spot_id"})
    df["spot_id"] = df["spot_id"].map(_clean_spot_id)
    non_type = {"spot_id", "Total cells"}
    type_cols = [c for c in df.columns if c not in non_type]
    for c in type_cols + (["Total cells"] if "Total cells" in df.columns else []):
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)
    if "Total cells" not in df.columns:
        df["Total cells"] = df[type_cols].sum(axis=1)
    df = df[["spot_id"] + type_cols + ["Total cells"]].copy()
    return df


def _align_types(a: pd.DataFrame, b: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    type_cols = sorted(set(a.columns).union(set(b.columns)) - {"spot_id", "Total cells"})
    aa = a.copy()
    bb = b.copy()
    for c in type_cols:
        if c not in aa.columns:
            aa[c] = 0.0
        if c not in bb.columns:
            bb[c] = 0.0
    aa = aa[["spot_id"] + type_cols + ["Total cells"]]
    bb = bb[["spot_id"] + type_cols + ["Total cells"]]
    return aa, bb, type_cols


def _row_normalize(x: np.ndarray) -> np.ndarray:
    s = x.sum(axis=1, keepdims=True)
    s[s <= 0] = 1.0
    return x / s


def _js_divergence_rows(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    p = np.clip(p, EPS, None)
    q = np.clip(q, EPS, None)
    p = p / p.sum(axis=1, keepdims=True)
    q = q / q.sum(axis=1, keepdims=True)
    m = 0.5 * (p + q)
    kl_pm = np.sum(p * np.log(p / m), axis=1)
    kl_qm = np.sum(q * np.log(q / m), axis=1)
    return 0.5 * (kl_pm + kl_qm)


def _entropy_rows(p: np.ndarray) -> np.ndarray:
    p = np.clip(p, EPS, None)
    p = p / p.sum(axis=1, keepdims=True)
    ent = -np.sum(p * np.log(p), axis=1)
    k = max(p.shape[1], 2)
    return ent / np.log(k)


def _row_corr(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    out = np.full(p.shape[0], np.nan, dtype=float)
    for i in range(p.shape[0]):
        x = p[i]
        y = q[i]
        if np.std(x) <= EPS or np.std(y) <= EPS:
            continue
        out[i] = float(np.corrcoef(x, y)[0, 1])
    return out


def main() -> int:
    args = parse_args()
    root = Path(args.project_root).resolve()
    baseline_dir = Path(args.baseline_dir).resolve()
    route2_dir = Path(args.route2_dir).resolve()
    out_dir = Path(args.out_dir).resolve() if args.out_dir else (root / "result" / args.sample / "stage5_real_summary")
    out_dir.mkdir(parents=True, exist_ok=True)

    base_summary = _read_json(baseline_dir / "stage4_summary.json")
    route_summary = _read_json(route2_dir / "stage4_summary.json")

    base = _load_by_spot(baseline_dir)
    route = _load_by_spot(route2_dir)
    base, route, type_cols = _align_types(base, route)

    b = base.set_index("spot_id")
    r = route.set_index("spot_id")
    common_spots = b.index.intersection(r.index)
    if len(common_spots) == 0:
        raise ValueError("no overlapping spots between baseline and route2")

    b = b.loc[common_spots]
    r = r.loc[common_spots]

    b_counts = b[type_cols].to_numpy(dtype=float)
    r_counts = r[type_cols].to_numpy(dtype=float)
    b_frac = _row_normalize(b_counts)
    r_frac = _row_normalize(r_counts)

    spot_l1 = np.sum(np.abs(r_frac - b_frac), axis=1)
    spot_js = _js_divergence_rows(r_frac, b_frac)
    spot_corr = _row_corr(r_frac, b_frac)

    b_entropy = _entropy_rows(b_frac)
    r_entropy = _entropy_rows(r_frac)
    b_top1 = np.max(b_frac, axis=1)
    r_top1 = np.max(r_frac, axis=1)
    b_active = np.sum(b_frac > EPS, axis=1)
    r_active = np.sum(r_frac > EPS, axis=1)

    base_mass = pd.Series(b_counts.sum(axis=0), index=type_cols, dtype=float)
    route_mass = pd.Series(r_counts.sum(axis=0), index=type_cols, dtype=float)
    delta_mass = route_mass - base_mass

    base_global = (base_mass / max(base_mass.sum(), EPS)).to_numpy(dtype=float)[None, :]
    route_global = (route_mass / max(route_mass.sum(), EPS)).to_numpy(dtype=float)[None, :]
    global_js = float(_js_divergence_rows(route_global, base_global)[0])

    unknown_name = "Unknown_sc_only" if "Unknown_sc_only" in type_cols else None
    unknown_base = float(base_mass.get(unknown_name, 0.0)) if unknown_name else 0.0
    unknown_route = float(route_mass.get(unknown_name, 0.0)) if unknown_name else 0.0

    top_delta = delta_mass.abs().sort_values(ascending=False).head(max(args.top_k, 1))
    top_delta_types = top_delta.index.tolist()

    type_mass_df = pd.DataFrame(
        {
            "type": type_cols,
            "baseline_mass": [float(base_mass.get(t, 0.0)) for t in type_cols],
            "route2_mass": [float(route_mass.get(t, 0.0)) for t in type_cols],
            "delta_mass": [float(delta_mass.get(t, 0.0)) for t in type_cols],
        }
    ).sort_values("delta_mass", ascending=False)
    type_mass_df.to_csv(out_dir / "stage5_real_type_mass.csv", index=False)

    spot_metrics_df = pd.DataFrame(
        {
            "spot_id": common_spots.astype(str),
            "l1": spot_l1,
            "js": spot_js,
            "corr": spot_corr,
            "entropy_baseline": b_entropy,
            "entropy_route2": r_entropy,
            "top1_baseline": b_top1,
            "top1_route2": r_top1,
            "n_active_baseline": b_active,
            "n_active_route2": r_active,
        }
    )
    spot_metrics_df.to_csv(out_dir / "stage5_real_spot_metrics.csv", index=False)

    summary = {
        "sample": args.sample,
        "stage": "stage5_real_summary",
        "paths": {
            "baseline_dir": str(baseline_dir),
            "route2_dir": str(route2_dir),
            "out_dir": str(out_dir),
        },
        "stage4": {
            "baseline_n_filtered": int(base_summary.get("n_filtered", 0)),
            "route2_n_filtered": int(route_summary.get("n_filtered", 0)),
            "route2_marked_unknown_total": int(route_summary.get("marked_unknown_total", 0)),
            "route2_restored_unknown_non_missing": int(route_summary.get("route2_restored_unknown_non_missing", 0)),
        },
        "distribution": {
            "n_common_spots": int(len(common_spots)),
            "n_types": int(len(type_cols)),
            "global_distribution_js": global_js,
            "spot_l1_mean": float(np.nanmean(spot_l1)),
            "spot_js_mean": float(np.nanmean(spot_js)),
            "spot_corr_mean": float(np.nanmean(spot_corr)),
            "spot_entropy_median_baseline": float(np.nanmedian(b_entropy)),
            "spot_entropy_median_route2": float(np.nanmedian(r_entropy)),
            "spot_entropy_delta_median": float(np.nanmedian(r_entropy - b_entropy)),
            "spot_top1_fraction_median_baseline": float(np.nanmedian(b_top1)),
            "spot_top1_fraction_median_route2": float(np.nanmedian(r_top1)),
            "spot_top1_fraction_delta_median": float(np.nanmedian(r_top1 - b_top1)),
            "spot_n_active_types_median_baseline": float(np.nanmedian(b_active)),
            "spot_n_active_types_median_route2": float(np.nanmedian(r_active)),
            "spot_n_active_types_delta_median": float(np.nanmedian(r_active - b_active)),
        },
        "mass": {
            "baseline_total_mass": float(base_mass.sum()),
            "route2_total_mass": float(route_mass.sum()),
            "unknown_mass_baseline": unknown_base,
            "unknown_mass_route2": unknown_route,
            "unknown_mass_delta": float(unknown_route - unknown_base),
        },
        "top_delta_types": top_delta_types,
    }

    out_json = out_dir / "stage5_real_summary.json"
    out_json.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[OK] wrote {out_json}")
    print(f"[OK] wrote {out_dir / 'stage5_real_type_mass.csv'}")
    print(f"[OK] wrote {out_dir / 'stage5_real_spot_metrics.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

