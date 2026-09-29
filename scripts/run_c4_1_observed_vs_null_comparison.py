#!/usr/bin/env python
"""C4.1-C observed-versus-pseudo-null comparison for three frozen real ST cases."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.config import load_project_config_yaml
from src.stages.stage3b_st_unsupported import (
    FEATURE_NAMES,
    _resolve_dataset_config,
    upper_tail_pvalues,
)
from src.stages.storage import processed_dir


CASES = (
    (
        "cytospace_fig2d_tme_brca_tnbc_fresh_frozen_sc_missing_plasma_cells",
        "BRCA TNBC – Plasma cells",
    ),
    (
        "cytospace_fig2d_tme_crc_fresh_frozen_sc_missing_b_cells",
        "CRC – B cells",
    ),
    (
        "cytospace_fig2d_tme_brca_her2_ffpe_sc_missing_plasma_cells",
        "BRCA HER2 FFPE – Plasma cells",
    ),
)

FEATURE_LABELS = {
    "relative_reconstruction_error": "Relative reconstruction error",
    "cosine_deficit": "Cosine deficit",
    "positive_residual_fraction": "Positive residual fraction",
    "residual_concentration": "Residual concentration",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project_root", default=str(ROOT))
    parser.add_argument(
        "--result_dir",
        default="result/c4_stage3b_calibration",
    )
    parser.add_argument(
        "--figure_dir",
        default="visualizations/c4_stage3b_calibration",
    )
    return parser.parse_args()


def resolve_under_root(root: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else root / path


def export_null(root: Path, result_root: Path, sample: str) -> Path:
    sample_dir = result_root / sample
    command = [
        sys.executable,
        str(root / "scripts" / "export_c4_1_stage3b_raw_calibration_null.py"),
        "--project_root",
        str(root),
        "--sample",
        sample,
        "--output_dir",
        str(sample_dir),
    ]
    completed = subprocess.run(
        command,
        cwd=root,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"Null export failed for {sample}:\n{completed.stderr.strip()}"
        )
    output = sample_dir / "c4_1_raw_null_statistics.csv"
    if not output.exists():
        raise FileNotFoundError(f"Null export missing for {sample}: {output}")
    return output


def load_case(
    root: Path,
    project_cfg: dict,
    result_root: Path,
    sample: str,
) -> tuple[pd.DataFrame, pd.DataFrame, float, Path]:
    null_path = export_null(root, result_root, sample)
    _, dataset_cfg = _resolve_dataset_config(root, sample, project_cfg, None)
    scores_path = (
        processed_dir(root, sample, dataset_cfg)
        / "stage3b_st_unsupported"
        / "spot_unsupported_scores.csv"
    )
    observed = pd.read_csv(scores_path, index_col=0)
    null = pd.read_csv(null_path, index_col=0)
    observed_values = observed.loc[:, FEATURE_NAMES].apply(pd.to_numeric, errors="raise")
    null_values = null.loc[:, FEATURE_NAMES].apply(pd.to_numeric, errors="raise")
    if observed_values.isna().any().any() or null_values.isna().any().any():
        raise ValueError(f"NaN found in observed/null statistics for {sample}")

    max_pvalue_difference = 0.0
    for feature in FEATURE_NAMES:
        reproduced = upper_tail_pvalues(
            observed_values[feature].to_numpy(dtype=float),
            null_values[feature].to_numpy(dtype=float),
        )
        expected = pd.to_numeric(
            observed[f"{feature}_pvalue"],
            errors="raise",
        ).to_numpy(dtype=float)
        max_pvalue_difference = max(
            max_pvalue_difference,
            float(np.max(np.abs(reproduced - expected))),
        )
    if max_pvalue_difference > 1e-12:
        raise ValueError(
            f"Stage3B p-value reproduction failed for {sample}: "
            f"max_abs_diff={max_pvalue_difference}"
        )
    return observed_values, null_values, max_pvalue_difference, null_path


def plot_qq(case_data: list[dict[str, object]], figure_dir: Path) -> list[Path]:
    fig, axes = plt.subplots(3, 4, figsize=(11.2, 7.2), dpi=300)
    probabilities = np.linspace(0.01, 0.99, 99)
    for row_index, case in enumerate(case_data):
        observed = case["observed"]
        null = case["null"]
        for column_index, feature in enumerate(FEATURE_NAMES):
            ax = axes[row_index, column_index]
            observed_quantiles = np.quantile(observed[feature], probabilities)
            null_quantiles = np.quantile(null[feature], probabilities)
            lower = float(min(observed_quantiles.min(), null_quantiles.min()))
            upper = float(max(observed_quantiles.max(), null_quantiles.max()))
            margin = max((upper - lower) * 0.04, np.finfo(float).eps)
            lower -= margin
            upper += margin
            ax.plot([lower, upper], [lower, upper], color="#B8BEC1", linewidth=0.8, zorder=1)
            ax.scatter(
                null_quantiles,
                observed_quantiles,
                s=9,
                color="#176F75",
                alpha=0.78,
                edgecolors="none",
                zorder=2,
            )
            ax.set_xlim(lower, upper)
            ax.set_ylim(lower, upper)
            ax.set_aspect("equal", adjustable="box")
            if row_index == 0:
                ax.set_title(FEATURE_LABELS[feature], fontsize=9, fontweight="normal")
            if row_index == len(case_data) - 1:
                ax.set_xlabel("Pseudo-null quantiles", fontsize=8)
            if column_index == 0:
                ax.set_ylabel(f"{case['label']}\nObserved quantiles", fontsize=8)
            ax.grid(color="#E6E8E9", linewidth=0.5)
            ax.tick_params(labelsize=7)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
    fig.tight_layout(pad=0.8, w_pad=1.0, h_pad=1.1)
    outputs = [
        figure_dir / "c4_1_observed_vs_null_qq.png",
        figure_dir / "c4_1_observed_vs_null_qq.pdf",
    ]
    fig.savefig(outputs[0], dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(outputs[1], bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return outputs


def ecdf(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.sort(np.asarray(values, dtype=float))
    y = np.arange(1, len(x) + 1, dtype=float) / len(x)
    return x, y


def plot_ecdf(case_data: list[dict[str, object]], figure_dir: Path) -> list[Path]:
    fig, axes = plt.subplots(3, 4, figsize=(11.2, 7.2), dpi=300)
    for row_index, case in enumerate(case_data):
        observed = case["observed"]
        null = case["null"]
        for column_index, feature in enumerate(FEATURE_NAMES):
            ax = axes[row_index, column_index]
            null_x, null_y = ecdf(null[feature].to_numpy(dtype=float))
            observed_x, observed_y = ecdf(observed[feature].to_numpy(dtype=float))
            ax.plot(null_x, null_y, color="#AEB5B8", linewidth=1.0, label="Pseudo-null")
            ax.plot(observed_x, observed_y, color="#176F75", linewidth=1.2, label="Observed")
            if row_index == 0:
                ax.set_title(FEATURE_LABELS[feature], fontsize=9, fontweight="normal")
            if row_index == len(case_data) - 1:
                ax.set_xlabel("Statistic value", fontsize=8)
            if column_index == 0:
                ax.set_ylabel(f"{case['label']}\nECDF", fontsize=8)
            ax.set_ylim(0, 1.01)
            ax.grid(axis="y", color="#E6E8E9", linewidth=0.5)
            ax.tick_params(labelsize=7)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            if row_index == 0 and column_index == 0:
                ax.legend(frameon=False, fontsize=7, loc="lower right")
    fig.tight_layout(pad=0.8, w_pad=1.0, h_pad=1.1)
    outputs = [
        figure_dir / "c4_1_observed_vs_null_ecdf.png",
        figure_dir / "c4_1_observed_vs_null_ecdf.pdf",
    ]
    fig.savefig(outputs[0], dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(outputs[1], bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return outputs


def main() -> int:
    args = parse_args()
    root = Path(args.project_root).resolve()
    result_root = resolve_under_root(root, args.result_dir)
    figure_dir = resolve_under_root(root, args.figure_dir)
    result_root.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    project_cfg = load_project_config_yaml(root, "configs/project_config.yaml")

    summary_rows: list[dict[str, object]] = []
    case_data: list[dict[str, object]] = []
    null_paths: dict[str, str] = {}
    reproduction: dict[str, float] = {}
    for sample, label in CASES:
        observed, null, max_diff, null_path = load_case(
            root,
            project_cfg,
            result_root,
            sample,
        )
        case_data.append(
            {"sample": sample, "label": label, "observed": observed, "null": null}
        )
        null_paths[sample] = str(null_path.resolve())
        reproduction[sample] = max_diff
        for feature in FEATURE_NAMES:
            observed_values = observed[feature].to_numpy(dtype=float)
            null_values = null[feature].to_numpy(dtype=float)
            observed_median = float(np.median(observed_values))
            null_median = float(np.median(null_values))
            summary_rows.append(
                {
                    "sample": sample,
                    "dataset_label": label,
                    "statistic": feature,
                    "n_observed": len(observed_values),
                    "n_null": len(null_values),
                    "observed_median": observed_median,
                    "null_median": null_median,
                    "observed_95th_percentile": float(np.quantile(observed_values, 0.95)),
                    "null_95th_percentile": float(np.quantile(null_values, 0.95)),
                    "observed_null_median_ratio": (
                        observed_median / null_median if null_median != 0 else np.nan
                    ),
                }
            )

    summary_path = result_root / "c4_1_observed_vs_null_summary.csv"
    pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
    qq_paths = plot_qq(case_data, figure_dir)
    ecdf_paths = plot_ecdf(case_data, figure_dir)
    print(
        json.dumps(
            {
                "summary": str(summary_path.resolve()),
                "qq": [str(path.resolve()) for path in qq_paths],
                "ecdf": [str(path.resolve()) for path in ecdf_paths],
                "raw_null": null_paths,
                "pvalue_reproduction_max_abs_difference": reproduction,
            },
            indent=2,
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
