#!/usr/bin/env python
"""Export the frozen Stage3B pseudo-null statistics for C4.1-B."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.config import load_project_config_yaml
from src.stages.stage3b_st_unsupported import (
    FEATURE_NAMES,
    Stage3BConfig,
    _compositional_normalize_rows,
    _expression_path,
    _linearize_expression,
    _read_expression,
    _resolve_dataset_config,
    _select_common_genes,
    anomaly_features,
    build_type_profiles,
    fit_nonnegative_mixtures,
    generate_supported_calibration,
    upper_tail_pvalues,
)
from src.stages.storage import processed_dir, result_dir


SAMPLE = "cytospace_fig2d_tme_brca_tnbc_fresh_frozen_sc_missing_plasma_cells"
EXPECTED_SEED = 42
EXPECTED_CALIBRATION = 1211


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Reproduce and export the Stage3B raw calibration null for C4.1-B."
    )
    parser.add_argument("--project_root", default=str(ROOT))
    parser.add_argument("--sample", default=SAMPLE)
    parser.add_argument(
        "--output_dir",
        default=None,
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    sample = str(args.sample)
    project_cfg = load_project_config_yaml(project_root, "configs/project_config.yaml")
    _, dataset_cfg = _resolve_dataset_config(
        project_root,
        sample,
        project_cfg,
        None,
    )

    stage3b_data = processed_dir(project_root, sample, dataset_cfg) / "stage3b_st_unsupported"
    summary_path = result_dir(project_root, sample, dataset_cfg) / "stage3b_st_unsupported" / "stage3b_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    config = Stage3BConfig(**summary["config"])
    n_calibration = int(summary["dimensions"]["calibration_spots"])
    if config.random_seed != EXPECTED_SEED:
        raise ValueError(f"Expected random_seed={EXPECTED_SEED}, found {config.random_seed}")
    if sample == SAMPLE and n_calibration != EXPECTED_CALIBRATION:
        raise ValueError(
            f"Expected n_calibration={EXPECTED_CALIBRATION}, found {n_calibration}"
        )

    export_dir = processed_dir(project_root, sample, dataset_cfg) / "stage1_preprocess" / "exported"
    profile_source = config.sc_profile_source or config.sc_expr_source
    profile_scale = config.sc_profile_scale or config.expression_scale
    sc_path = _expression_path(export_dir, "sc", profile_source)
    st_path = _expression_path(export_dir, "st", config.sc_expr_source)
    metadata_path = export_dir / "sc_metadata.csv"

    sc = _read_expression(sc_path)
    st = _read_expression(st_path)
    metadata = pd.read_csv(metadata_path, index_col=0)
    metadata.index = metadata.index.astype(str)
    shared_cells = sc.index.intersection(metadata.index, sort=False)
    sc = sc.loc[shared_cells]
    labels = metadata.loc[shared_cells, "cell_type"].astype(str)
    genes = _select_common_genes(sc, st, config.max_genes)
    sc_values = _compositional_normalize_rows(
        _linearize_expression(
            sc.loc[:, genes].to_numpy(dtype=np.float64),
            profile_scale,
        )
    )
    st_values = _compositional_normalize_rows(
        _linearize_expression(
            st.loc[:, genes].to_numpy(dtype=np.float64),
            config.expression_scale,
        )
    )

    profiles, type_names, cells_by_type = build_type_profiles(sc_values, labels)
    fitted_weights, reconstruction = fit_nonnegative_mixtures(st_values, profiles)
    observed_features = anomaly_features(st_values, reconstruction)

    rng = np.random.default_rng(config.random_seed)
    pseudo = generate_supported_calibration(
        fitted_weights,
        cells_by_type,
        type_names,
        n_calibration,
        rng,
    )
    _, pseudo_reconstruction = fit_nonnegative_mixtures(pseudo, profiles)
    calibration_features = anomaly_features(pseudo, pseudo_reconstruction)

    output_dir = Path(args.output_dir or f"result/c4_stage3b_calibration/{sample}")
    if not output_dir.is_absolute():
        output_dir = project_root / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "c4_1_raw_null_statistics.csv"
    null_frame = pd.DataFrame(calibration_features, columns=FEATURE_NAMES)
    null_frame.index = np.arange(1, len(null_frame) + 1)
    null_frame.index.name = "pseudo_null_id"
    null_frame.to_csv(output_path)

    existing_scores = pd.read_csv(
        stage3b_data / "spot_unsupported_scores.csv",
        index_col=0,
    )
    existing_scores.index = existing_scores.index.astype(str)
    existing_scores = existing_scores.reindex(st.index)
    pvalue_differences: dict[str, float] = {}
    observed_statistic_differences: dict[str, float] = {}
    for feature_index, feature_name in enumerate(FEATURE_NAMES):
        reproduced_p = upper_tail_pvalues(
            observed_features[:, feature_index],
            calibration_features[:, feature_index],
        )
        expected_p = existing_scores[f"{feature_name}_pvalue"].to_numpy(dtype=float)
        pvalue_differences[feature_name] = float(
            np.max(np.abs(reproduced_p - expected_p))
        )
        expected_observed = existing_scores[feature_name].to_numpy(dtype=float)
        observed_statistic_differences[feature_name] = float(
            np.max(np.abs(observed_features[:, feature_index] - expected_observed))
        )

    print(
        json.dumps(
            {
                "sample": sample,
                "random_seed": config.random_seed,
                "n_calibration": n_calibration,
                "output": str(output_path.resolve()),
                "pvalue_max_absolute_difference": pvalue_differences,
                "observed_statistic_max_absolute_difference": observed_statistic_differences,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
