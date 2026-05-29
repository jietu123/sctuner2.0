from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon


PROJECT_ROOT = Path(__file__).resolve().parents[1]
FOUNDATION_DIR = PROJECT_ROOT / "result" / "real_profile_mask_foundation"
PROCESSED_ROOT = PROJECT_ROOT / "data" / "processed" / "low_resolution_experiments"
RESULT_ROOT = PROJECT_ROOT / "result"
OUT_DIR = PROJECT_ROOT / "result" / "real_profile_mask_expression_recovery"

BASELINE = "CytoSPACE"
ROUTE2 = "SVTuner + CytoSPACE"
METHOD_ASSIGNMENTS = {
    BASELINE: Path("stage4_cytospace_baseline") / "cytospace_output" / "cell_assignment.csv",
    ROUTE2: Path("stage4_cytospace_route2") / "cytospace_output" / "cell_assignment.csv",
}


def _read_header(path: Path) -> list[str]:
    return pd.read_csv(path, nrows=0).columns.tolist()


def _select_hvgs(st_expr: pd.DataFrame, common_genes: list[str], max_genes: int) -> list[str]:
    if not common_genes:
        raise ValueError("no shared genes between sc and ST expression matrices")
    variances = st_expr[common_genes].var(axis=0).sort_values(ascending=False)
    return variances.head(min(max_genes, len(variances))).index.tolist()


def _recover_spot_expression(
    assignment_csv: Path,
    sc_expr: pd.DataFrame,
    genes: list[str],
    target_spots: pd.Index,
) -> pd.DataFrame:
    assign = pd.read_csv(assignment_csv, usecols=["cell_id", "assigned_spot"])
    assign["cell_id"] = assign["cell_id"].astype(str)
    assign["assigned_spot"] = assign["assigned_spot"].astype(str)

    assigned_cells = set(assign["cell_id"])
    sc_expr = sc_expr[sc_expr.index.isin(assigned_cells)]

    assign = assign[assign["cell_id"].isin(sc_expr.index)].set_index("cell_id")
    expr = sc_expr.loc[assign.index, genes].copy()
    expr["assigned_spot"] = assign["assigned_spot"].to_numpy()
    recovered = expr.groupby("assigned_spot", sort=False)[genes].sum()
    recovered = recovered.reindex(target_spots.astype(str)).fillna(0.0)
    return recovered


def _row_cosine(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    numerator = np.sum(a * b, axis=1)
    denom = np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1)
    out = np.full(a.shape[0], np.nan, dtype=float)
    valid = denom > 0
    out[valid] = numerator[valid] / denom[valid]
    return out


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(FOUNDATION_DIR / "foundation_manifest.csv")

    spot_rows: list[dict[str, object]] = []
    scenario_rows: list[dict[str, object]] = []

    for row in manifest.itertuples(index=False):
        sample = str(row.sample)
        processed = PROCESSED_ROOT / sample / "stage1_preprocess" / "exported"
        st_expr_csv = processed / "st_expression_normalized.csv"
        sc_expr_csv = processed / "sc_expression_normalized.csv"

        st_expr = pd.read_csv(st_expr_csv).set_index("spot_id")
        st_expr.index = st_expr.index.astype(str)
        st_genes = [c for c in st_expr.columns]
        sc_header = _read_header(sc_expr_csv)
        sc_genes = [c for c in sc_header if c != "cell_id"]
        common_genes = [g for g in st_genes if g in set(sc_genes)]
        genes = _select_hvgs(st_expr, common_genes, max_genes=2000)
        st_selected = st_expr[genes].astype(float)
        sc_selected = pd.read_csv(sc_expr_csv, usecols=["cell_id", *genes])
        sc_selected["cell_id"] = sc_selected["cell_id"].astype(str)
        sc_selected = sc_selected.set_index("cell_id")

        method_values: dict[str, float] = {}
        for method, rel_assignment in METHOD_ASSIGNMENTS.items():
            assignment_csv = RESULT_ROOT / sample / rel_assignment
            if not assignment_csv.exists():
                raise FileNotFoundError(f"missing assignment for {sample} {method}: {assignment_csv}")
            recovered = _recover_spot_expression(assignment_csv, sc_selected, genes, st_selected.index)
            common_spots = st_selected.index.intersection(recovered.index)
            st_mat = st_selected.loc[common_spots, genes].to_numpy(dtype=float)
            rec_mat = recovered.loc[common_spots, genes].to_numpy(dtype=float)
            cosine = _row_cosine(rec_mat, st_mat)
            mean_cosine = float(np.nanmean(cosine))
            median_cosine = float(np.nanmedian(cosine))
            method_values[method] = mean_cosine

            scenario_rows.append(
                {
                    "sample": sample,
                    "source_sample": row.source_sample,
                    "target_type": row.target_type,
                    "scenario_label": row.scenario_label,
                    "method": method,
                    "metric": "expression_recovery_cosine",
                    "mean_cosine": mean_cosine,
                    "median_cosine": median_cosine,
                    "n_spots": int(len(common_spots)),
                    "n_genes": int(len(genes)),
                }
            )
            for spot_id, value in zip(common_spots, cosine, strict=True):
                spot_rows.append(
                    {
                        "sample": sample,
                        "scenario_label": row.scenario_label,
                        "method": method,
                        "spot_id": spot_id,
                        "cosine_similarity": float(value) if np.isfinite(value) else np.nan,
                    }
                )

        scenario_wide_row = {
            "sample": sample,
            "source_sample": row.source_sample,
            "target_type": row.target_type,
            "scenario_label": row.scenario_label,
        }
        for method in METHOD_ASSIGNMENTS:
            scenario_wide_row[method] = method_values.get(method, np.nan)
        scenario_wide_row["delta_route2_minus_baseline"] = scenario_wide_row[ROUTE2] - scenario_wide_row[BASELINE]
        scenario_wide_row["route2_better"] = scenario_wide_row["delta_route2_minus_baseline"] > 0
        scenario_rows.append({**scenario_wide_row, "method": "__paired__", "metric": "expression_recovery_cosine"})

    scenario_long = pd.DataFrame([r for r in scenario_rows if r.get("method") != "__paired__"])
    scenario_wide = pd.DataFrame([r for r in scenario_rows if r.get("method") == "__paired__"]).drop(columns=["method", "metric"])
    spot_df = pd.DataFrame(spot_rows)

    pair = scenario_wide[[BASELINE, ROUTE2]].dropna()
    if len(pair) >= 2 and np.any(np.abs(pair[ROUTE2] - pair[BASELINE]) > 1e-12):
        stat = wilcoxon(pair[ROUTE2], pair[BASELINE], alternative="greater", zero_method="wilcox")
        pvalue = float(stat.pvalue)
    else:
        pvalue = float("nan")

    cfg = {
        "metric": "expression_recovery_cosine",
        "metric_full_name": "Masked-ST expression recovery cosine similarity",
        "definition": "Cosine similarity between masked ST expression and spot expression recovered by summing mapped single-cell expression.",
        "reference_expression": "masked ST expression used by each real profile-mask scenario",
        "gene_selection": "top 2000 shared genes by variance in masked ST expression",
        "n_scenarios": int(scenario_wide["sample"].nunique()),
        "n_route2_better": int(scenario_wide["route2_better"].sum()),
        "method_order": list(METHOD_ASSIGNMENTS.keys()),
        "mean_baseline": float(scenario_wide[BASELINE].mean()),
        "mean_route2": float(scenario_wide[ROUTE2].mean()),
        "method_means": {method: float(scenario_wide[method].mean()) for method in METHOD_ASSIGNMENTS},
        "mean_delta_route2_minus_baseline": float(scenario_wide["delta_route2_minus_baseline"].mean()),
        "wilcoxon_pvalue_route2_gt_baseline": pvalue,
    }

    scenario_long.to_csv(OUT_DIR / "expression_recovery_long.csv", index=False)
    scenario_wide.to_csv(OUT_DIR / "expression_recovery_by_scenario.csv", index=False)
    spot_df.to_csv(OUT_DIR / "expression_recovery_by_spot.csv", index=False)
    (OUT_DIR / "expression_recovery_config.json").write_text(json.dumps(cfg, indent=2, ensure_ascii=False), encoding="utf-8")

    print(f"[OK] wrote: {OUT_DIR / 'expression_recovery_by_scenario.csv'}")
    print(json.dumps(cfg, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
