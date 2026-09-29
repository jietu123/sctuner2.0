"""Build the frozen C8 canonical pair tables and direct paired statistics."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "result" / "c8_direct_paired_comparison"

BOOTSTRAP_ITERATIONS = 10_000
BOOTSTRAP_SEED = 20_260_927

METRIC_ORDER = [
    "A1_reciprocal_suppression",
    "A2_cosine_similarity",
    "A3_ecotyper_experiment_mean",
    "B1_merscope_suppression",
    "B2_merscope_peak_es",
]

METADATA = {
    "A1_reciprocal_suppression": {
        "metric_label": "Reciprocal suppression",
        "evidence_layer": "Low-resolution profile masking",
        "metric_direction": "higher_is_better",
    },
    "A2_cosine_similarity": {
        "metric_label": "Reconstructed-expression cosine similarity",
        "evidence_layer": "Low-resolution profile masking",
        "metric_direction": "higher_is_better",
    },
    "A3_ecotyper_experiment_mean": {
        "metric_label": "EcoTyper normalized enrichment (experiment mean)",
        "evidence_layer": "EcoTyper profile-masking enrichment",
        "metric_direction": "higher_is_better",
    },
    "B1_merscope_suppression": {
        "metric_label": "MERSCOPE low-support suppression score",
        "evidence_layer": "MERSCOPE cell-resolution",
        "metric_direction": "higher_is_better",
    },
    "B2_merscope_peak_es": {
        "metric_label": "MERSCOPE peak ES / low-support enrichment",
        "evidence_layer": "MERSCOPE cell-resolution",
        "metric_direction": "higher_is_better",
    },
}


def _canonical_frame(
    metric_id: str,
    frame: pd.DataFrame,
    *,
    independent_unit_id: pd.Series,
    dataset: pd.Series,
    target: pd.Series,
    readout: pd.Series | str,
    cytospace: pd.Series,
    svtuner: pd.Series,
) -> pd.DataFrame:
    meta = METADATA[metric_id]
    n = len(frame)
    readout_values = pd.Series([readout] * n, index=frame.index) if isinstance(readout, str) else readout
    result = pd.DataFrame(
        {
            "evidence_layer": meta["evidence_layer"],
            "metric_id": metric_id,
            "metric_label": meta["metric_label"],
            "independent_unit_id": independent_unit_id.astype(str),
            "dataset": dataset.astype(str),
            "target": target.astype(str),
            "readout": readout_values.astype(str),
            "cytospace": pd.to_numeric(cytospace),
            "svtuner": pd.to_numeric(svtuner),
            "metric_direction": meta["metric_direction"],
        }
    ).reset_index(drop=True)
    result["delta_raw"] = result["svtuner"] - result["cytospace"]
    if meta["metric_direction"] == "higher_is_better":
        result["delta_favorable"] = result["delta_raw"]
    elif meta["metric_direction"] == "lower_is_better":
        result["delta_favorable"] = -result["delta_raw"]
    else:
        raise ValueError(f"Unknown metric direction for {metric_id}")
    return result[
        [
            "evidence_layer",
            "metric_id",
            "metric_label",
            "independent_unit_id",
            "dataset",
            "target",
            "readout",
            "cytospace",
            "svtuner",
            "delta_raw",
            "metric_direction",
            "delta_favorable",
        ]
    ]


def load_a1() -> pd.DataFrame:
    path = (
        ROOT
        / "visualizations"
        / "simulations"
        / "real_profile_mask_fig2d_only"
        / "fig2_panel_d_real_profile_mask_source_values.csv"
    )
    data = pd.read_csv(path)
    data = data[
        (data["metric"] == "low_support_suppression_score")
        & data["method"].isin(["CytoSPACE", "SVTuner + CytoSPACE"])
    ]
    keys = ["sample", "target_type"]
    if data.duplicated(keys + ["method"]).any():
        raise ValueError("Duplicate A1 pairing key")
    wide = data.pivot(index=keys, columns="method", values="value").reset_index()
    return _canonical_frame(
        "A1_reciprocal_suppression",
        wide,
        independent_unit_id=wide["sample"] + " | " + wide["target_type"],
        dataset=wide["sample"],
        target=wide["target_type"],
        readout="reciprocal suppression",
        cytospace=wide["CytoSPACE"],
        svtuner=wide["SVTuner + CytoSPACE"],
    )


def load_a2() -> pd.DataFrame:
    path = (
        ROOT
        / "reproducibility_release"
        / "release_candidate_v1.0"
        / "source_values"
        / "SV-0042_expression_recovery_long.csv"
    )
    data = pd.read_csv(path)
    data = data[data["method"].isin(["CytoSPACE", "SVTuner + CytoSPACE"])]
    keys = ["sample", "target_type"]
    if data.duplicated(keys + ["method"]).any():
        raise ValueError("Duplicate A2 pairing key")
    wide = data.pivot(index=keys, columns="method", values="mean_cosine").reset_index()
    return _canonical_frame(
        "A2_cosine_similarity",
        wide,
        independent_unit_id=wide["sample"] + " | " + wide["target_type"],
        dataset=wide["sample"],
        target=wide["target_type"],
        readout="mean cosine similarity",
        cytospace=wide["CytoSPACE"],
        svtuner=wide["SVTuner + CytoSPACE"],
    )


def load_a3() -> tuple[pd.DataFrame, pd.DataFrame]:
    path = (
        ROOT
        / "visualizations"
        / "cytospace_fig2d_profile_mask_benchmark"
        / "fig2d_profile_mask_benchmark_source_values.csv"
    )
    data = pd.read_csv(path)
    data = data[data["method"].isin(["CytoSPACE", "SVTuner + CytoSPACE"])]
    readout_keys = ["pair_id", "profile_mask_sample", "masked_target_type", "cell_type"]
    if data.duplicated(readout_keys + ["method"]).any():
        raise ValueError("Duplicate A3 readout pairing key")
    readout_wide = data.pivot(index=readout_keys, columns="method", values="nes").reset_index()
    readout_pairs = _canonical_frame(
        "A3_ecotyper_experiment_mean",
        readout_wide,
        independent_unit_id=readout_wide["pair_id"] + " | " + readout_wide["cell_type"],
        dataset=readout_wide["profile_mask_sample"],
        target=readout_wide["masked_target_type"],
        readout=readout_wide["cell_type"],
        cytospace=readout_wide["CytoSPACE"],
        svtuner=readout_wide["SVTuner + CytoSPACE"],
    )
    readout_pairs["pair_id"] = readout_wide["pair_id"].astype(str).to_numpy()
    readout_pairs["inferential_unit"] = False

    experiment = (
        readout_wide.groupby("pair_id", sort=False, as_index=False)
        .agg(
            dataset=("profile_mask_sample", "first"),
            target=("masked_target_type", "first"),
            n_readouts=("cell_type", "size"),
            cytospace=("CytoSPACE", "mean"),
            svtuner=("SVTuner + CytoSPACE", "mean"),
        )
    )
    canonical = _canonical_frame(
        "A3_ecotyper_experiment_mean",
        experiment,
        independent_unit_id=experiment["pair_id"],
        dataset=experiment["dataset"],
        target=experiment["target"],
        readout="mean of CD4/CD8 readouts",
        cytospace=experiment["cytospace"],
        svtuner=experiment["svtuner"],
    )
    return canonical, readout_pairs


def load_b1() -> pd.DataFrame:
    path = ROOT / "result" / "c5_decomposition" / "c5_merscope_four_route_suppression_wide.csv"
    data = pd.read_csv(path)
    if data.duplicated(["dataset", "target"]).any():
        raise ValueError("Duplicate B1 pairing key")
    return _canonical_frame(
        "B1_merscope_suppression",
        data,
        independent_unit_id=data["dataset"] + " | " + data["target"],
        dataset=data["dataset"],
        target=data["target"],
        readout="low-support suppression score",
        cytospace=data["baseline"],
        svtuner=data["full"],
    )


def load_b2() -> pd.DataFrame:
    path = ROOT / "result" / "c5_decomposition" / "c5_final_quantitative_source.csv"
    data = pd.read_csv(path)
    data = data[data["resolution"] == "MERSCOPE"].copy()
    expected_metric = "peak_es_low_masked_support_enrichment"
    if set(data["mapping_quality_metric"].dropna()) != {expected_metric}:
        raise ValueError("Unexpected B2 metric definition")
    if data.duplicated(["dataset", "target"]).any():
        raise ValueError("Duplicate B2 pairing key")
    return _canonical_frame(
        "B2_merscope_peak_es",
        data,
        independent_unit_id=data["dataset"] + " | " + data["target"],
        dataset=data["dataset"],
        target=data["target"],
        readout=expected_metric,
        cytospace=data["baseline_peak_es"],
        svtuner=data["full_peak_es"],
    )


def bootstrap_mean_ci(differences: np.ndarray) -> tuple[float, float]:
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    indices = rng.integers(0, len(differences), size=(BOOTSTRAP_ITERATIONS, len(differences)))
    means = differences[indices].mean(axis=1)
    low, high = np.percentile(means, [2.5, 97.5])
    return float(low), float(high)


def calculate_statistics(canonical: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for metric_id in METRIC_ORDER:
        group = canonical[canonical["metric_id"] == metric_id].copy()
        raw = group["delta_raw"].to_numpy(dtype=float)
        favorable = group["delta_favorable"].to_numpy(dtype=float)
        low, high = bootstrap_mean_ci(raw)
        sample_sd = float(np.std(raw, ddof=1))
        dz = np.nan if sample_sd == 0 else float(np.mean(raw) / sample_sd)
        test = wilcoxon(raw, alternative="two-sided", method="auto")
        meta = METADATA[metric_id]
        rows.append(
            {
                "metric_id": metric_id,
                "metric_label": meta["metric_label"],
                "evidence_layer": meta["evidence_layer"],
                "n_pairs": len(group),
                "metric_direction": meta["metric_direction"],
                "cytospace_mean": group["cytospace"].mean(),
                "cytospace_median": group["cytospace"].median(),
                "svtuner_mean": group["svtuner"].mean(),
                "svtuner_median": group["svtuner"].median(),
                "mean_delta_raw": raw.mean(),
                "median_delta_raw": np.median(raw),
                "mean_delta_ci_low": low,
                "mean_delta_ci_high": high,
                "mean_delta_favorable": favorable.mean(),
                "median_delta_favorable": np.median(favorable),
                "cohen_dz": dz,
                "wins": int(np.sum(favorable > 0)),
                "ties": int(np.sum(favorable == 0)),
                "losses": int(np.sum(favorable < 0)),
                "wilcoxon_statistic": float(test.statistic),
                "wilcoxon_p_two_sided": float(test.pvalue),
                "bootstrap_iterations": BOOTSTRAP_ITERATIONS,
                "bootstrap_seed": BOOTSTRAP_SEED,
            }
        )
    return pd.DataFrame(rows)


def run_qc(canonical: pd.DataFrame, readouts: pd.DataFrame, stats: pd.DataFrame) -> None:
    expected_n = {
        "A1_reciprocal_suppression": 10,
        "A2_cosine_similarity": 10,
        "A3_ecotyper_experiment_mean": 6,
        "B1_merscope_suppression": 5,
        "B2_merscope_peak_es": 5,
    }
    observed = canonical.groupby("metric_id").size().to_dict()
    if observed != expected_n:
        raise ValueError(f"Unexpected inferential pair counts: {observed}")
    if len(readouts) != 12 or readouts["inferential_unit"].any():
        raise ValueError("A3 readout table must contain 12 non-inferential pairs")
    if canonical.duplicated(["metric_id", "independent_unit_id"]).any():
        raise ValueError("Duplicate canonical pairing key")
    numeric = ["cytospace", "svtuner", "delta_raw", "delta_favorable"]
    if canonical[numeric].isna().any().any() or not np.isfinite(canonical[numeric].to_numpy()).all():
        raise ValueError("NA/Inf in canonical metric values")
    readout_numeric = ["cytospace", "svtuner", "delta_raw", "delta_favorable"]
    if readouts[readout_numeric].isna().any().any() or not np.isfinite(readouts[readout_numeric].to_numpy()).all():
        raise ValueError("NA/Inf in A3 readout values")
    if len(stats) != 5 or stats["n_pairs"].to_list() != [10, 10, 6, 5, 5]:
        raise ValueError("Statistics table failed row-count QC")


def fmt(value: float) -> str:
    return f"{value:.6f}"


def write_summary(stats: pd.DataFrame) -> str:
    lines = ["C8-1 DIRECT PAIRED STATISTICS", ""]
    headings = {
        "A1_reciprocal_suppression": "A1 Reciprocal suppression",
        "A2_cosine_similarity": "A2 Cosine similarity",
        "A3_ecotyper_experiment_mean": "A3 EcoTyper enrichment",
        "B1_merscope_suppression": "B1 MERSCOPE suppression",
        "B2_merscope_peak_es": "B2 MERSCOPE peak ES",
    }
    by_id = stats.set_index("metric_id")
    for metric_id in METRIC_ORDER:
        row = by_id.loc[metric_id]
        lines.append(headings[metric_id])
        if metric_id == "A3_ecotyper_experiment_mean":
            lines.append("12 readouts shown descriptively; inference based on 6 independent experiments.")
        lines.extend(
            [
                f"n = {int(row['n_pairs'])}",
                f"CytoSPACE mean = {fmt(row['cytospace_mean'])}",
                f"SVTuner mean = {fmt(row['svtuner_mean'])}",
                f"mean delta = {fmt(row['mean_delta_raw'])}",
                f"95% CI = [{fmt(row['mean_delta_ci_low'])}, {fmt(row['mean_delta_ci_high'])}]",
                f"Cohen dz = {fmt(row['cohen_dz'])}",
                f"wins/ties/losses = {int(row['wins'])}/{int(row['ties'])}/{int(row['losses'])}",
                f"two-sided Wilcoxon P = {fmt(row['wilcoxon_p_two_sided'])}",
                "",
            ]
        )
    lines.extend(["QC:", "PASS", "", "STATUS:", "C8-1 = PASS", ""])
    text = "\n".join(lines)
    (OUT / "C8_1_summary.txt").write_text(text, encoding="utf-8")
    return text


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    a1 = load_a1()
    a2 = load_a2()
    a3, readouts = load_a3()
    b1 = load_b1()
    b2 = load_b2()
    canonical = pd.concat([a1, a2, a3, b1, b2], ignore_index=True)
    stats = calculate_statistics(canonical)
    run_qc(canonical, readouts, stats)
    canonical.to_csv(OUT / "c8_canonical_pairs.csv", index=False)
    readouts.to_csv(OUT / "c8_ecotyper_readout_pairs.csv", index=False)
    stats.to_csv(OUT / "c8_paired_statistics.csv", index=False)
    write_summary(stats)
    print("C8-1 = PASS")


if __name__ == "__main__":
    main()
