from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "result" / "c5_decomposition"
ROUTES = OUT / "c5_merscope_four_route_summary.csv"
FROZEN_SOURCE = ROOT / "visualizations" / "highres_profile_mask_fig3b_enrichment" / "fig3b_highres_profile_mask_enrichment_source_values.csv"
FROZEN_METRICS = ROOT / "visualizations" / "highres_profile_mask_fig3b_enrichment" / "fig3b_highres_profile_mask_enrichment_metrics.csv"
LOW_SUPPORT_QUANTILE = 0.20


def processed_dir(sample: str) -> Path:
    candidates = [
        ROOT / "data" / "processed" / sample,
        ROOT / "data" / "processed" / "highres_profile_mask" / sample,
    ]
    for path in candidates:
        if (path / "stage1_preprocess" / "exported").exists():
            return path
    raise FileNotFoundError(f"Cannot resolve processed directory for {sample}")


def reconstructed_spot_score(assignment_path: Path, sc_expr: pd.DataFrame,
                             genes: list[str], ordered_spots: list[str]) -> pd.Series:
    assignment = pd.read_csv(assignment_path, usecols=["cell_id", "assigned_spot"])
    assignment["cell_id"] = assignment["cell_id"].astype(str)
    assignment["assigned_spot"] = assignment["assigned_spot"].astype(str)
    assignment = assignment[assignment["cell_id"].isin(sc_expr.index)].copy()
    merged = assignment.join(sc_expr, on="cell_id", how="left")
    spot_score = merged.groupby("assigned_spot")[genes].mean().mean(axis=1)
    return spot_score.reindex(ordered_spots).fillna(0.0).astype(float)


def read_mask(path: Path, ordered_spots: list[str]) -> pd.Series:
    mask = pd.read_csv(path, usecols=["spot_id", "is_unsupported_region"])
    mask["spot_id"] = mask["spot_id"].astype(str)
    mask["withheld"] = mask["is_unsupported_region"].astype(str).str.lower().eq("true")
    result = mask.set_index("spot_id")["withheld"].reindex(ordered_spots)
    if result.isna().any() or len(mask) != len(ordered_spots):
        raise RuntimeError(f"Mask does not preserve the full spatial-unit index: {path}")
    return result.astype(bool)


route_manifest = pd.read_csv(ROUTES)
frozen_source = pd.read_csv(FROZEN_SOURCE)
frozen_metrics = pd.read_csv(FROZEN_METRICS)
long_rows = []

route_specs = [
    ("Baseline", "baseline_path", False, False, None),
    ("Stage3A-only", "stage3a_only_path", True, False, None),
    ("Stage3B-only", "stage3b_only_path", False, True, "standalone_mask_path"),
    ("Full", "full_path", True, True, "sequential_mask_path"),
]

for record in route_manifest.itertuples(index=False):
    sample = str(record.dataset)
    fixed = frozen_source[frozen_source["profile_mask_sample"].astype(str).eq(sample)].copy()
    if fixed.empty:
        raise RuntimeError(f"Frozen residual-support table missing {sample}")
    fixed = fixed.sort_values("rank", kind="mergesort")
    ordered_spots = fixed["spot_id"].astype(str).tolist()
    low_cut = float(fixed["masked_support_norm"].quantile(LOW_SUPPORT_QUANTILE))
    low_mask = fixed["masked_support_norm"].le(low_cut).to_numpy()
    n_low = int(low_mask.sum())

    metric_row = frozen_metrics[frozen_metrics["profile_mask_sample"].astype(str).eq(sample)]
    if len(metric_row) != 1:
        raise RuntimeError(f"Frozen marker panel metadata missing or duplicated for {sample}")
    genes = [g for g in str(metric_row.iloc[0]["marker_genes"]).split(";") if g]
    processed = processed_dir(sample)
    sc_path = processed / "stage1_preprocess" / "exported" / "sc_expression_normalized.csv"
    sc_cols = set(pd.read_csv(sc_path, nrows=0).columns.astype(str))
    genes = [gene for gene in genes if gene in sc_cols]
    if len(genes) != int(metric_row.iloc[0]["n_marker_genes"]):
        raise RuntimeError(f"Frozen marker panel cannot be reproduced for {sample}")
    sc_expr = pd.read_csv(sc_path, usecols=["cell_id", *genes]).set_index("cell_id")
    sc_expr.index = sc_expr.index.astype(str)
    sc_expr = sc_expr.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    standalone_mask_path = Path(str(record.standalone_mask_path))
    sequential_mask_path = Path(str(record.sequential_mask_path))
    if standalone_mask_path.resolve() == sequential_mask_path.resolve():
        raise RuntimeError(f"Full does not have an independent sequential mask for {sample}")

    for route, path_field, stage3a_on, stage3b_on, mask_field in route_specs:
        output_dir = Path(str(getattr(record, path_field)))
        assignment_path = output_dir / "cell_assignment.csv"
        if not assignment_path.exists():
            raise FileNotFoundError(f"Missing existing mapping assignment: {assignment_path}")
        scores = reconstructed_spot_score(assignment_path, sc_expr, genes, ordered_spots)
        if stage3b_on:
            mask_path = Path(str(getattr(record, mask_field)))
            withheld = read_mask(mask_path, ordered_spots)
            n_withheld = int(withheld.sum())
            n_withheld_low = int(withheld.to_numpy()[low_mask].sum())
        else:
            mask_path = None
            n_withheld = 0
            n_withheld_low = 0
        low_mean = float(scores.to_numpy()[low_mask].mean())
        suppression = 1.0 / (1.0 + low_mean)
        long_rows.append({
            "dataset": sample,
            "target": record.target,
            "capacity": int(record.capacity),
            "route": route,
            "stage3a_on": stage3a_on,
            "stage3b_on": stage3b_on,
            "n_total_units": len(ordered_spots),
            "n_low_support_units": n_low,
            "n_withheld_total": n_withheld,
            "n_withheld_in_low_support": n_withheld_low,
            "mean_mapped_target_signal_low_support": low_mean,
            "low_support_suppression_score": suppression,
            "low_support_quantile": LOW_SUPPORT_QUANTILE,
            "low_support_cutoff": low_cut,
            "n_marker_genes": len(genes),
            "marker_genes": ";".join(genes),
            "assignment_path": str(assignment_path),
            "stage3b_mask_path": "" if mask_path is None else str(mask_path),
            "fixed_support_source": str(FROZEN_SOURCE),
        })

long_df = pd.DataFrame(long_rows)
if len(long_df) != 20 or long_df["dataset"].nunique() != 5:
    raise RuntimeError("Expected exactly 5 datasets x 4 routes")
if long_df.isna().any().any():
    raise RuntimeError("Undefined values found in four-route endpoint table")
for _, group in long_df.groupby("dataset"):
    if len(group) != 4 or group["n_total_units"].nunique() != 1 or group["n_low_support_units"].nunique() != 1:
        raise RuntimeError("Evaluation denominator differs across routes")

long_path = OUT / "c5_merscope_four_route_suppression_long.csv"
long_df.to_csv(long_path, index=False)

wide = long_df.pivot(index=["dataset", "target", "capacity"], columns="route",
                     values="low_support_suppression_score").reset_index()
wide = wide.rename(columns={"Baseline": "baseline", "Stage3A-only": "stage3a_only",
                            "Stage3B-only": "stage3b_only", "Full": "full"})
wide["effect_stage3a"] = wide["stage3a_only"] - wide["baseline"]
wide["effect_stage3b"] = wide["stage3b_only"] - wide["baseline"]
wide["effect_full"] = wide["full"] - wide["baseline"]
wide["increment_stage3b_after_stage3a"] = wide["full"] - wide["stage3a_only"]
wide["increment_stage3a_after_stage3b"] = wide["full"] - wide["stage3b_only"]
wide_path = OUT / "c5_merscope_four_route_suppression_wide.csv"
wide.to_csv(wide_path, index=False)
print("SUCCESS")
