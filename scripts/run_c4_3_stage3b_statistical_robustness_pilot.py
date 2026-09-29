from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SAMPLE = "cytospace_fig2d_tme_brca_tnbc_fresh_frozen_sc_missing_plasma_cells"
GROUP = "stage3b_realdata_reference_dropout"
FDR = 0.05
SEED = 42


def adjust_pvalues(pvalues: np.ndarray, dependence_factor: float) -> np.ndarray:
    pvalues = np.asarray(pvalues, dtype=np.float64)
    n = len(pvalues)
    order = np.argsort(pvalues, kind="stable")
    ranked = pvalues[order]
    adjusted = ranked * n * dependence_factor / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty(n, dtype=np.float64)
    result[order] = np.clip(adjusted, 0.0, 1.0)
    return result


def region_signatures(scores: pd.DataFrame) -> dict[tuple[str, ...], int]:
    signatures: dict[tuple[str, ...], int] = {}
    for region_id, sub in scores.loc[scores["region_id"] > 0].groupby("region_id"):
        signatures[tuple(sorted(sub["spot_id"].astype(str)))] = int(region_id)
    return signatures


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", default=SAMPLE)
    args = parser.parse_args()
    sample = args.sample
    original_data = (
        ROOT / "data" / "processed" / GROUP / sample / "stage3b_st_unsupported"
    )
    original_result = ROOT / "result" / GROUP / sample / "stage3b_st_unsupported"
    scores_path = original_data / "spot_unsupported_scores.csv"
    regions_path = original_data / "unsupported_regions.csv"
    summary_path = original_result / "stage3b_summary.json"
    dataset_config = ROOT / "configs" / "datasets" / f"{sample}.yaml"
    for path in (scores_path, regions_path, summary_path, dataset_config):
        if not path.is_file():
            raise FileNotFoundError(path)

    original_summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if int(original_summary["config"]["n_spatial_permutations"]) != 200:
        raise ValueError("Original Stage3B result did not use 200 permutations")
    if int(original_summary["config"]["random_seed"]) != SEED:
        raise ValueError("Original Stage3B result did not use seed 42")

    scores = pd.read_csv(scores_path)
    pvalues = scores["whole_profile_pvalue"].to_numpy(dtype=np.float64)
    bh_q = adjust_pvalues(pvalues, 1.0)
    harmonic = float(np.sum(1.0 / np.arange(1, len(pvalues) + 1)))
    by_q = adjust_pvalues(pvalues, harmonic)
    saved_bh_q = scores["whole_profile_qvalue"].to_numpy(dtype=np.float64)
    max_bh_difference = float(np.max(np.abs(bh_q - saved_bh_q)))
    if max_bh_difference > 1e-12:
        raise ValueError(
            f"Recomputed BH q-values do not match saved values: {max_bh_difference}"
        )
    bh_mask = bh_q <= FDR
    saved_candidate_mask = scores["is_whole_profile_candidate"].astype(bool).to_numpy()
    if not np.array_equal(bh_mask, saved_candidate_mask):
        raise ValueError("Saved whole-profile candidate mask does not match BH q-values")
    by_mask = by_q <= FDR
    intersection = int(np.sum(bh_mask & by_mask))
    union = int(np.sum(bh_mask | by_mask))
    jaccard = float(intersection / union) if union else 1.0

    output_dir = ROOT / "result" / "c4_stage3b_statistical_robustness" / sample
    output_dir.mkdir(parents=True, exist_ok=True)
    spot_comparison = pd.DataFrame(
        {
            "spot_id": scores["spot_id"].astype(str),
            "whole_profile_pvalue": pvalues,
            "bh_qvalue": bh_q,
            "by_qvalue": by_q,
            "bh_significant": bh_mask,
            "by_significant": by_mask,
        }
    )
    spot_comparison.to_csv(output_dir / "c4_3_spot_bh_by_comparison.csv", index=False)

    command = [
        sys.executable,
        "-m",
        "src.stages.stage3b_st_unsupported",
        "--project_root",
        str(ROOT),
        "--sample",
        sample,
        "--dataset_config",
        str(dataset_config),
        "--output_suffix",
        "_c4_p1000",
        "--n_spatial_permutations",
        "1000",
        "--random_seed",
        str(SEED),
    ]
    print("[C4.3] Running Stage3B with 1000 spatial permutations", flush=True)
    completed = subprocess.run(command, cwd=ROOT, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"1000-permutation Stage3B failed: exit {completed.returncode}")

    new_data = (
        ROOT
        / "data"
        / "processed"
        / GROUP
        / sample
        / "stage3b_st_unsupported_c4_p1000"
    )
    new_result = (
        ROOT
        / "result"
        / GROUP
        / sample
        / "stage3b_st_unsupported_c4_p1000"
    )
    scores_1000 = pd.read_csv(new_data / "spot_unsupported_scores.csv")
    regions_200 = pd.read_csv(regions_path)
    regions_1000 = pd.read_csv(new_data / "unsupported_regions.csv")
    summary_1000 = json.loads(
        (new_result / "stage3b_summary.json").read_text(encoding="utf-8")
    )
    if int(summary_1000["config"]["n_spatial_permutations"]) != 1000:
        raise ValueError("New Stage3B result did not use 1000 permutations")
    if int(summary_1000["config"]["random_seed"]) != SEED:
        raise ValueError("New Stage3B result did not use seed 42")

    signatures_200 = region_signatures(scores)
    signatures_1000 = region_signatures(scores_1000)
    if set(signatures_200) != set(signatures_1000):
        raise ValueError("The 200- and 1000-permutation region memberships differ")
    regions_200_by_id = regions_200.set_index("region_id")
    regions_1000_by_id = regions_1000.set_index("region_id")
    region_rows: list[dict[str, object]] = []
    for signature, region_200_id in signatures_200.items():
        region_1000_id = signatures_1000[signature]
        old = regions_200_by_id.loc[region_200_id]
        new = regions_1000_by_id.loc[region_1000_id]
        region_rows.append(
            {
                "region_id_200": region_200_id,
                "region_id_1000": region_1000_id,
                "n_spots": int(old["n_spots"]),
                "region_mass": float(old["region_mass"]),
                "pvalue_200": float(old["region_pvalue"]),
                "pvalue_1000": float(new["region_pvalue"]),
                "absolute_pvalue_difference": abs(
                    float(old["region_pvalue"]) - float(new["region_pvalue"])
                ),
                "significant_200": bool(old["is_unsupported_region"]),
                "significant_1000": bool(new["is_unsupported_region"]),
            }
        )
    region_comparison = pd.DataFrame(region_rows).sort_values("region_id_200")
    region_comparison.to_csv(
        output_dir / "c4_3_region_permutation_comparison.csv", index=False
    )

    old_final_mask = scores["is_unsupported_region"].astype(bool).to_numpy()
    new_by_spot = scores_1000.set_index("spot_id").reindex(scores["spot_id"].astype(str))
    if new_by_spot.isna().all(axis=1).any():
        raise ValueError("Could not align 1000-permutation spot output")
    new_final_mask = new_by_spot["is_unsupported_region"].astype(bool).to_numpy()
    final_intersection = int(np.sum(old_final_mask & new_final_mask))
    final_union = int(np.sum(old_final_mask | new_final_mask))
    final_jaccard = float(final_intersection / final_union) if final_union else 1.0

    method_summary = pd.DataFrame(
        [
            {
                "method": "BH",
                "significant_spots": int(bh_mask.sum()),
                "significant_fraction": float(bh_mask.mean()),
                "overlap_spots": intersection,
                "jaccard_vs_other": jaccard,
            },
            {
                "method": "BY",
                "significant_spots": int(by_mask.sum()),
                "significant_fraction": float(by_mask.mean()),
                "overlap_spots": intersection,
                "jaccard_vs_other": jaccard,
            },
        ]
    )
    permutation_summary = pd.DataFrame(
        [
            {
                "permutations": 200,
                "regions": len(regions_200),
                "significant_regions": int(regions_200["is_unsupported_region"].sum()),
                "final_withheld_spots": int(old_final_mask.sum()),
                "region_pvalue_min": float(regions_200["region_pvalue"].min()),
                "region_pvalue_max": float(regions_200["region_pvalue"].max()),
            },
            {
                "permutations": 1000,
                "regions": len(regions_1000),
                "significant_regions": int(regions_1000["is_unsupported_region"].sum()),
                "final_withheld_spots": int(new_final_mask.sum()),
                "region_pvalue_min": float(regions_1000["region_pvalue"].min()),
                "region_pvalue_max": float(regions_1000["region_pvalue"].max()),
            },
        ]
    )
    method_summary.to_csv(output_dir / "c4_3_bh_by_summary.csv", index=False)
    permutation_summary.to_csv(
        output_dir / "c4_3_permutation_summary.csv", index=False
    )
    run_summary = {
        "sample": sample,
        "fdr": FDR,
        "seed": SEED,
        "bh_by_intersection_spots": intersection,
        "bh_by_union_spots": union,
        "bh_by_jaccard": jaccard,
        "max_recomputed_bh_qvalue_difference": max_bh_difference,
        "max_region_pvalue_absolute_difference": float(
            region_comparison["absolute_pvalue_difference"].max()
        ),
        "changed_region_significance_count": int(
            (region_comparison["significant_200"] != region_comparison["significant_1000"]).sum()
        ),
        "final_mask_jaccard_200_vs_1000": final_jaccard,
    }
    (output_dir / "c4_3_pilot_summary.json").write_text(
        json.dumps(run_summary, indent=2), encoding="utf-8"
    )
    print("\nBH vs BY")
    print(method_summary.to_string(index=False))
    print("\n200 vs 1000 spatial permutations")
    print(permutation_summary.to_string(index=False))
    print(json.dumps(run_summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
