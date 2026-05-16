from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from scipy.stats import wilcoxon


PROJECT_ROOT = Path(__file__).resolve().parents[1]
FOUNDATION_DIR = PROJECT_ROOT / "result" / "real_profile_mask_foundation"
OUT_DIR = PROJECT_ROOT / "result" / "real_profile_mask_fig2d_score"

METRIC = "low_support_suppression_score"
METRIC_LABEL = "Low-support suppression score"
METRIC_MEANING = (
    "Reciprocal transform of the mean reconstructed target-like score "
    "inside the lowest-support 20% of masked-ST spots. Higher is better."
)
LOW_QUANTILE = 0.20
SCORE_COLUMNS = {
    "baseline_reconstructed_score": "CytoSPACE",
    "route2_reconstructed_score": "SVTuner + CytoSPACE",
    "tangram_reconstructed_score": "Tangram",
    "celltrek_reconstructed_score": "CellTrek",
}


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(FOUNDATION_DIR / "foundation_manifest.csv")

    long_rows: list[dict[str, object]] = []

    for row in manifest.itertuples(index=False):
        sample = row.sample
        df = pd.read_csv(FOUNDATION_DIR / sample / "spot_foundation.csv")
        low_mask = df["masked_support_norm"] <= df["masked_support_norm"].quantile(LOW_QUANTILE)

        available = [(col, method) for col, method in SCORE_COLUMNS.items() if col in df.columns]
        missing = [col for col in SCORE_COLUMNS if col not in df.columns]
        if missing:
            raise ValueError(f"{sample} foundation missing score columns: {missing}")
        for score_col, method in available:
            low_mean = float(df.loc[low_mask, score_col].mean())
            suppression_score = 1.0 / (1.0 + low_mean)
            long_rows.append(
                {
                    "sample": sample,
                    "source_sample": row.source_sample,
                    "target_type": row.target_type,
                    "scenario_label": row.scenario_label,
                    "method": method,
                    "metric": METRIC,
                    "metric_label": METRIC_LABEL,
                    "direction": "higher",
                    "low_support_quantile": LOW_QUANTILE,
                    "raw_low_support_mean": low_mean,
                    "value": suppression_score,
                    "meaning": METRIC_MEANING,
                }
            )

    long_df = pd.DataFrame(long_rows)
    by_scenario = (
        long_df.pivot(
            index=["sample", "source_sample", "target_type", "scenario_label"],
            columns="method",
            values="value",
        )
        .reset_index()
        .rename_axis(columns=None)
    )
    by_scenario["delta_route2_minus_baseline"] = by_scenario["SVTuner + CytoSPACE"] - by_scenario["CytoSPACE"]
    by_scenario["route2_better"] = by_scenario["delta_route2_minus_baseline"] > 0

    summary = (
        long_df.groupby("method", as_index=False)
        .agg(
            mean_score=("value", "mean"),
            median_score=("value", "median"),
            min_score=("value", "min"),
            max_score=("value", "max"),
            mean_raw_low_support_mean=("raw_low_support_mean", "mean"),
        )
    )

    pair = by_scenario[["CytoSPACE", "SVTuner + CytoSPACE"]].dropna()
    stat = wilcoxon(pair["SVTuner + CytoSPACE"], pair["CytoSPACE"], alternative="greater", zero_method="wilcox")
    pairwise_vs_cytospace = {}
    for method in ["SVTuner + CytoSPACE", "Tangram", "CellTrek"]:
        paired = by_scenario[["CytoSPACE", method]].dropna()
        if len(paired) < 2:
            pairwise_vs_cytospace[method] = None
            continue
        pairwise_vs_cytospace[method] = {
            "n": int(len(paired)),
            "wilcoxon_pvalue_method_gt_cytospace": float(
                wilcoxon(paired[method], paired["CytoSPACE"], alternative="greater", zero_method="wilcox").pvalue
            ),
            "mean_delta_method_minus_cytospace": float((paired[method] - paired["CytoSPACE"]).mean()),
        }

    config = {
        "metric": METRIC,
        "metric_label": METRIC_LABEL,
        "meaning": METRIC_MEANING,
        "low_support_quantile": LOW_QUANTILE,
        "n_scenarios": int(by_scenario["sample"].nunique()),
        "n_route2_better": int(by_scenario["route2_better"].sum()),
        "mean_score_baseline": float(summary.loc[summary["method"] == "CytoSPACE", "mean_score"].iloc[0]),
        "mean_score_route2": float(summary.loc[summary["method"] == "SVTuner + CytoSPACE", "mean_score"].iloc[0]),
        "wilcoxon_pvalue_route2_gt_baseline": float(stat.pvalue),
        "pairwise_vs_cytospace": pairwise_vs_cytospace,
    }

    long_df.to_csv(OUT_DIR / "fig2d_score_long.csv", index=False)
    by_scenario.to_csv(OUT_DIR / "fig2d_score_by_scenario.csv", index=False)
    summary.to_csv(OUT_DIR / "fig2d_score_summary.csv", index=False)
    (OUT_DIR / "fig2d_score_config.json").write_text(
        json.dumps(config, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    print(f"[OK] wrote: {OUT_DIR / 'fig2d_score_long.csv'}")
    print(f"[OK] wrote: {OUT_DIR / 'fig2d_score_by_scenario.csv'}")
    print(f"[OK] wrote: {OUT_DIR / 'fig2d_score_summary.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
