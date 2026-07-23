from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon


PROJECT_ROOT = Path(__file__).resolve().parents[1]
FOUNDATION_DIR = PROJECT_ROOT / "result" / "real_profile_mask_foundation"
OUT_DIR = PROJECT_ROOT / "result" / "real_profile_mask_fig2d_nes"

BASELINE = "CytoSPACE"
ROUTE2 = "SVTuner + CytoSPACE"
METHODS = [
    ("baseline_reconstructed_score", BASELINE),
    ("route2_reconstructed_score", ROUTE2),
]

LOW_SUPPORT_QUANTILE = 0.15
N_PERMUTATIONS = 1000
RNG_SEED = 42


def _running_es(rank_scores: np.ndarray, hits: np.ndarray) -> float:
    scores = np.asarray(rank_scores, dtype=float)
    hits = np.asarray(hits, dtype=bool)
    n = len(scores)
    nh = int(hits.sum())
    if nh == 0 or nh == n:
        return float("nan")
    weights = np.abs(scores)
    hit_denom = float(weights[hits].sum())
    if hit_denom <= 0:
        weights = np.ones(n, dtype=float)
        hit_denom = float(nh)
    miss_w = 1.0 / float(n - nh)
    running = np.cumsum(np.where(hits, weights / hit_denom, -miss_w))
    pos = float(running.max())
    neg = float(running.min())
    return pos if abs(pos) >= abs(neg) else neg


def _signal_from_values(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    return -(arr - arr.mean()) / (arr.std(ddof=0) + 1e-12)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(FOUNDATION_DIR / "foundation_manifest.csv")
    rng = np.random.default_rng(RNG_SEED)

    rows: list[dict[str, object]] = []

    for item in manifest.itertuples(index=False):
        spot_df = pd.read_csv(FOUNDATION_DIR / item.sample / "spot_foundation.csv")
        low_thr = float(spot_df["masked_support_norm"].quantile(LOW_SUPPORT_QUANTILE))
        hit_mask = (spot_df["masked_support_norm"] <= low_thr).to_numpy(dtype=bool)
        nh = int(hit_mask.sum())

        for score_col, method in METHODS:
            signal = _signal_from_values(spot_df[score_col].to_numpy(dtype=float))
            order = np.argsort(-signal, kind="mergesort")
            ordered_signal = signal[order]
            ordered_hits = hit_mask[order]

            es = _running_es(ordered_signal, ordered_hits)
            null_es = np.empty(N_PERMUTATIONS, dtype=float)
            for i in range(N_PERMUTATIONS):
                perm_hits = np.zeros(len(ordered_hits), dtype=bool)
                perm_hits[rng.choice(len(ordered_hits), size=nh, replace=False)] = True
                null_es[i] = _running_es(ordered_signal, perm_hits)

            # Use absolute-null normalization to stabilize scenario-level direction
            # across masked-real samples while retaining a proper permutation-based scale.
            denom = float(np.mean(np.abs(null_es)))
            nes = float(es / denom) if denom > 0 else float("nan")

            rows.append(
                {
                    "sample": item.sample,
                    "source_sample": item.source_sample,
                    "target_type": item.target_type,
                    "scenario_label": item.scenario_label,
                    "method": method,
                    "score_column": score_col,
                    "n_spots": int(len(spot_df)),
                    "n_hits": nh,
                    "low_support_quantile": LOW_SUPPORT_QUANTILE,
                    "signal_definition": "z_inverse_reconstructed_target_like_score",
                    "es": es,
                    "nes": nes,
                }
            )

    long_df = pd.DataFrame(rows)
    by_scenario = long_df.pivot_table(index=["sample", "source_sample", "target_type", "scenario_label"], columns="method", values="nes", aggfunc="first").reset_index()
    by_scenario["delta_route2_minus_baseline"] = by_scenario[ROUTE2] - by_scenario[BASELINE]
    by_scenario["route2_better"] = by_scenario["delta_route2_minus_baseline"] > 0

    summary = (
        long_df.groupby("method", as_index=False)
        .agg(
            mean_nes=("nes", "mean"),
            median_nes=("nes", "median"),
            min_nes=("nes", "min"),
            max_nes=("nes", "max"),
            mean_es=("es", "mean"),
        )
    )

    pair = by_scenario[[BASELINE, ROUTE2]].dropna()
    stat = wilcoxon(pair[ROUTE2], pair[BASELINE], alternative="greater", zero_method="wilcox")
    config = {
        "metric": "low_support_suppression_nes",
        "metric_label": "Suppression NES",
        "meaning": "NES for enrichment of low-support masked-ST spots among locations ranked by per-method suppression signal, normalized by mean absolute null ES. Higher is better.",
        "low_support_quantile": LOW_SUPPORT_QUANTILE,
        "signal_definition": "z_inverse_reconstructed_target_like_score",
        "n_permutations": N_PERMUTATIONS,
        "rng_seed": RNG_SEED,
        "n_scenarios": int(len(pair)),
        "n_route2_better": int(by_scenario["route2_better"].sum()),
        "mean_nes_baseline": float(summary.loc[summary["method"] == BASELINE, "mean_nes"].iloc[0]),
        "mean_nes_route2": float(summary.loc[summary["method"] == ROUTE2, "mean_nes"].iloc[0]),
        "wilcoxon_pvalue_route2_gt_baseline": float(stat.pvalue),
    }

    long_df.to_csv(OUT_DIR / "fig2d_nes_long.csv", index=False)
    by_scenario.to_csv(OUT_DIR / "fig2d_nes_by_scenario.csv", index=False)
    summary.to_csv(OUT_DIR / "fig2d_nes_summary.csv", index=False)
    (OUT_DIR / "fig2d_nes_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")

    print(f"[OK] wrote: {OUT_DIR / 'fig2d_nes_long.csv'}")
    print(f"[OK] wrote: {OUT_DIR / 'fig2d_nes_by_scenario.csv'}")
    print(f"[OK] wrote: {OUT_DIR / 'fig2d_nes_summary.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
