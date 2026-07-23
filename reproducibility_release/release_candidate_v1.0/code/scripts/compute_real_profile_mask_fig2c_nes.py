from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon


PROJECT_ROOT = Path(__file__).resolve().parents[1]
FOUNDATION_DIR = PROJECT_ROOT / "result" / "real_profile_mask_foundation"
OUT_DIR = PROJECT_ROOT / "result" / "real_profile_mask_fig2c_nes"

METHOD_SPECS = [
    ("baseline_reconstructed_score", "CytoSPACE"),
    ("route2_reconstructed_score", "SVTuner + CytoSPACE"),
]

RANK_DIRECTION = "low_support_to_high_support"
HIT_FRACTION = 0.10
N_PERMUTATIONS = 1000
PERM_BATCH_SIZE = 200
RNG_SEED = 42


def enrichment_score_from_hits(hits: np.ndarray) -> float:
    hits = hits.astype(bool, copy=False)
    n = hits.size
    nh = int(hits.sum())
    if nh == 0 or nh == n:
        return 0.0
    hit_w = 1.0 / nh
    miss_w = 1.0 / (n - nh)
    running = np.cumsum(np.where(hits, hit_w, -miss_w))
    max_v = float(running.max())
    min_v = float(running.min())
    return max_v if abs(max_v) >= abs(min_v) else min_v


def permuted_es_distribution(n: int, nh: int, n_perm: int, batch_size: int, rng: np.random.Generator) -> np.ndarray:
    if nh <= 0 or nh >= n:
        return np.zeros(n_perm, dtype=float)
    out = np.empty(n_perm, dtype=float)
    done = 0
    while done < n_perm:
        cur = min(batch_size, n_perm - done)
        rand = rng.random((cur, n))
        hit_pos = np.argpartition(rand, nh - 1, axis=1)[:, :nh]
        hits = np.zeros((cur, n), dtype=np.int8)
        row_idx = np.arange(cur)[:, None]
        hits[row_idx, hit_pos] = 1
        hit_w = 1.0 / nh
        miss_w = 1.0 / (n - nh)
        running = np.cumsum(np.where(hits.astype(bool), hit_w, -miss_w), axis=1)
        max_v = running.max(axis=1)
        min_v = running.min(axis=1)
        out[done : done + cur] = np.where(np.abs(max_v) >= np.abs(min_v), max_v, min_v)
        done += cur
    return out


def compute_observed_hits(scores: pd.Series, frac: float) -> np.ndarray:
    threshold = float(scores.quantile(frac))
    hits = (scores <= threshold).to_numpy(dtype=bool)
    if hits.sum() == 0:
        # Fallback in degenerate case.
        n = len(scores)
        keep = max(1, int(round(n * frac)))
        order = np.argsort(scores.to_numpy(dtype=float), kind="mergesort")
        hits = np.zeros(n, dtype=bool)
        hits[order[:keep]] = True
    return hits


def normalized_enrichment_score(observed_es: float, null_es: np.ndarray) -> tuple[float, float]:
    if observed_es >= 0:
        pos = null_es[null_es >= 0]
        scale = float(np.mean(pos)) if pos.size else float(np.mean(np.abs(null_es)))
        tail = np.sum(null_es >= observed_es)
    else:
        neg = null_es[null_es < 0]
        scale = float(np.mean(np.abs(neg))) if neg.size else float(np.mean(np.abs(null_es)))
        tail = np.sum(null_es <= observed_es)
    if scale <= 0:
        nes = 0.0
    else:
        nes = float(observed_es / scale)
    pvalue = float((tail + 1) / (len(null_es) + 1))
    return nes, pvalue


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(FOUNDATION_DIR / "foundation_manifest.csv")
    rng = np.random.default_rng(RNG_SEED)

    long_rows: list[dict[str, object]] = []
    null_rows: list[dict[str, object]] = []

    for row in manifest.itertuples(index=False):
        sample = row.sample
        spot_df = pd.read_csv(FOUNDATION_DIR / sample / "spot_foundation.csv")
        spot_df = spot_df.sort_values("masked_support_norm", ascending=True).reset_index(drop=True)
        n_spots = int(len(spot_df))

        for score_col, method in METHOD_SPECS:
            hits = compute_observed_hits(spot_df[score_col], HIT_FRACTION)
            nh = int(hits.sum())
            observed_es = enrichment_score_from_hits(hits)
            null_es = permuted_es_distribution(
                n=n_spots,
                nh=nh,
                n_perm=N_PERMUTATIONS,
                batch_size=PERM_BATCH_SIZE,
                rng=rng,
            )
            nes, pvalue = normalized_enrichment_score(observed_es, null_es)

            long_rows.append(
                {
                    "sample": sample,
                    "source_sample": row.source_sample,
                    "target_type": row.target_type,
                    "scenario_label": row.scenario_label,
                    "method": method,
                    "score_column": score_col,
                    "rank_direction": RANK_DIRECTION,
                    "hit_definition": f"bottom_{int(HIT_FRACTION*100)}pct_reconstructed_score",
                    "n_spots": n_spots,
                    "n_hits": nh,
                    "observed_es": observed_es,
                    "nes": nes,
                    "permutation_pvalue": pvalue,
                }
            )

            null_rows.append(
                {
                    "sample": sample,
                    "method": method,
                    "null_es_mean": float(np.mean(null_es)),
                    "null_es_abs_mean": float(np.mean(np.abs(null_es))),
                    "null_es_std": float(np.std(null_es, ddof=1)),
                    "null_es_pos_mean": float(np.mean(null_es[null_es >= 0])) if np.any(null_es >= 0) else np.nan,
                    "null_es_neg_abs_mean": float(np.mean(np.abs(null_es[null_es < 0]))) if np.any(null_es < 0) else np.nan,
                }
            )

    long_df = pd.DataFrame(long_rows)
    null_df = pd.DataFrame(null_rows)
    by_scenario = (
        long_df.pivot(
            index=["sample", "source_sample", "target_type", "scenario_label"],
            columns="method",
            values="nes",
        )
        .reset_index()
        .rename_axis(columns=None)
    )
    by_scenario["delta_route2_minus_baseline"] = by_scenario["SVTuner + CytoSPACE"] - by_scenario["CytoSPACE"]
    by_scenario["route2_better"] = by_scenario["delta_route2_minus_baseline"] > 0

    summary = (
        long_df.groupby("method", as_index=False)
        .agg(
            mean_nes=("nes", "mean"),
            median_nes=("nes", "median"),
            min_nes=("nes", "min"),
            max_nes=("nes", "max"),
            mean_es=("observed_es", "mean"),
            mean_perm_pvalue=("permutation_pvalue", "mean"),
        )
    )

    pair = by_scenario[["CytoSPACE", "SVTuner + CytoSPACE"]].dropna()
    stat = wilcoxon(
        pair["SVTuner + CytoSPACE"],
        pair["CytoSPACE"],
        alternative="greater",
        zero_method="wilcox",
    )

    config = {
        "rank_direction": RANK_DIRECTION,
        "hit_definition": f"bottom_{int(HIT_FRACTION*100)}pct_reconstructed_score",
        "hit_fraction": HIT_FRACTION,
        "n_permutations": N_PERMUTATIONS,
        "perm_batch_size": PERM_BATCH_SIZE,
        "rng_seed": RNG_SEED,
        "n_scenarios": int(by_scenario["sample"].nunique()),
        "n_route2_better": int(by_scenario["route2_better"].sum()),
        "mean_nes_baseline": float(summary.loc[summary["method"] == "CytoSPACE", "mean_nes"].iloc[0]),
        "mean_nes_route2": float(summary.loc[summary["method"] == "SVTuner + CytoSPACE", "mean_nes"].iloc[0]),
        "wilcoxon_pvalue_route2_gt_baseline": float(stat.pvalue),
    }

    long_df.to_csv(OUT_DIR / "fig2c_nes_long.csv", index=False)
    by_scenario.to_csv(OUT_DIR / "fig2c_nes_by_scenario.csv", index=False)
    summary.to_csv(OUT_DIR / "fig2c_nes_summary.csv", index=False)
    null_df.to_csv(OUT_DIR / "fig2c_nes_null_summary.csv", index=False)
    (OUT_DIR / "fig2c_nes_config.json").write_text(
        json.dumps(config, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    print(f"[OK] wrote: {OUT_DIR / 'fig2c_nes_long.csv'}")
    print(f"[OK] wrote: {OUT_DIR / 'fig2c_nes_by_scenario.csv'}")
    print(f"[OK] wrote: {OUT_DIR / 'fig2c_nes_summary.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
