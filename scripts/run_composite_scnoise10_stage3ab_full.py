#!/usr/bin/env python
"""Run and audit the isolated full Stage3A + Stage3B 10% noise benchmark."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.plot_method_comparison_composition_recovery import (  # noqa: E402
    COMPOSITE_SCENARIOS,
    METHODS,
    _prediction_path,
    build_tables,
    plot,
)
from src.stages.storage import (  # noqa: E402
    processed_dir,
    read_dataset_config,
    result_dir,
    stage1_export_dir,
)


SUFFIX = "_scnoise10_stage3ab_full"
STAGE3_DIR = f"stage3_typematch{SUFFIX}"
STAGE3B_DIR = f"stage3b_st_unsupported{SUFFIX}"
STAGE4_DIR = f"stage4_cytospace{SUFFIX}"
NOISE_SUFFIX = "_scnoise10"
OUT_REL = Path("visualizations/method_comparison/composite_scnoise10_stage3ab_full")


@dataclass(frozen=True)
class Scenario:
    group: str
    base_sample: str
    condition: str
    stage3b_target: str

    @property
    def sample(self) -> str:
        return f"{self.base_sample}{NOISE_SUFFIX}"


SCENARIOS = [
    Scenario(*COMPOSITE_SCENARIOS[0], "control", "Endothelial cells"),
    Scenario(*COMPOSITE_SCENARIOS[1], "single_missing", "Endothelial cells"),
    Scenario(*COMPOSITE_SCENARIOS[2], "double_missing", "Endothelial cells"),
    Scenario(*COMPOSITE_SCENARIOS[3], "control", "B_cell"),
    Scenario(*COMPOSITE_SCENARIOS[4], "single_missing", "B_cell"),
    Scenario(*COMPOSITE_SCENARIOS[5], "double_missing", "B_cell"),
    Scenario(*COMPOSITE_SCENARIOS[6], "control", "Ext_L56"),
    Scenario(*COMPOSITE_SCENARIOS[7], "single_missing", "Ext_L56"),
    Scenario(*COMPOSITE_SCENARIOS[8], "double_missing", "Ext_L56"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    # Resolve from cwd on Windows; some isolated Conda builds misdecode non-ASCII __file__ paths.
    parser.add_argument("--project_root", default=".")
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--execute", action="store_true", help="run Stage3A/B/4 after audit")
    parser.add_argument("--resume", action="store_true", help="reuse complete new outputs")
    parser.add_argument("--n_processors", type=int, default=1)
    parser.add_argument("--n_subspots", type=int, default=800)
    parser.add_argument("--mapping_cells_per_spot", type=int, default=5)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def csv_index(path: Path) -> set[str]:
    return set(pd.read_csv(path, usecols=[0]).iloc[:, 0].astype(str))


def truth_path(root: Path, scenario: Scenario) -> Path:
    return root / "data" / "sim" / scenario.group / scenario.sample / "sim_truth_spot_type_fraction.csv"


def method_output(root: Path, sample: str, method_dir: str) -> Path:
    return _prediction_path(root, sample, method_dir, route2_stage4_dir=STAGE4_DIR)


def run_command(command: list[str], log_path: Path, root: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["PYTHONNOUSERSITE"] = "1"
    env["CYTOSPACE_SKIP_ASSIGNED_EXPRESSION"] = "1"
    with log_path.open("w", encoding="utf-8") as log:
        log.write("COMMAND\n" + subprocess.list2cmdline(command) + "\n\nOUTPUT\n")
        log.flush()
        completed = subprocess.run(
            command,
            cwd=root,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    if completed.returncode:
        raise RuntimeError(f"Command failed ({completed.returncode}); see {log_path}")


def input_and_reuse_audit(root: Path, out_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    input_rows: list[dict[str, Any]] = []
    reuse_rows: list[dict[str, Any]] = []
    reusable_methods = [item for item in METHODS if item[1] != "cytospace_route2"]
    for scenario in SCENARIOS:
        cfg = read_dataset_config(root, scenario.sample)
        export = stage1_export_dir(root, scenario.sample, cfg)
        info_path = truth_path(root, scenario).parent / "sim_info.json"
        noise_summary_path = root / "result" / scenario.sample / "sc_noise_generation_summary.json"
        info = read_json(info_path)
        source_sample = str(info["source_sample"])
        source_cfg = read_dataset_config(root, source_sample)
        source_export = stage1_export_dir(root, source_sample, source_cfg)
        current_truth = truth_path(root, scenario)
        source_truth = root / "data" / "sim" / scenario.group / source_sample / current_truth.name
        paths = {
            "sc_expression_normalized": export / "sc_expression_normalized.csv",
            "st_expression_normalized": export / "st_expression_normalized.csv",
            "st_coordinates": export / "st_coordinates.csv",
            "simulation_truth": current_truth,
        }
        source_paths = {
            "sc_expression_normalized": source_export / "sc_expression_normalized.csv",
            "st_expression_normalized": source_export / "st_expression_normalized.csv",
            "st_coordinates": source_export / "st_coordinates.csv",
            "simulation_truth": source_truth,
        }
        for kind, path in paths.items():
            source = source_paths[kind]
            input_rows.append(
                {
                    "group": scenario.group,
                    "sample": scenario.sample,
                    "condition": scenario.condition,
                    "input_kind": kind,
                    "path": str(path),
                    "sha256": sha256(path),
                    "source_sample": source_sample,
                    "source_path": str(source),
                    "source_sha256": sha256(source),
                    "same_as_no_noise_source": sha256(path) == sha256(source),
                    "expected_same_as_source": kind != "sc_expression_normalized",
                    "noise_fraction": info.get("noise_fraction"),
                    "noise_seed": info.get("noise_seed"),
                }
            )
        truth_spots = csv_index(current_truth)
        stage1_st_spots = csv_index(paths["st_expression_normalized"])
        if not stage1_st_spots.issubset(truth_spots):
            raise RuntimeError(f"Stage1 ST spots are not a subset of truth: {scenario.sample}")
        generated_at = max(
            paths["sc_expression_normalized"].stat().st_mtime,
            noise_summary_path.stat().st_mtime if noise_summary_path.exists() else 0,
        )
        noise_summary = read_json(noise_summary_path) if noise_summary_path.exists() else {}
        normalized_noise = noise_summary.get("normalized", {})
        deterministic_recreation = bool(
            noise_summary.get("source_sample") == source_sample
            and noise_summary.get("target_sample") == scenario.sample
            and float(noise_summary.get("noise_fraction", -1)) == 0.1
            and int(noise_summary.get("seed", -1)) == 42
            and int(normalized_noise.get("perturbed_entries", 0)) > 0
            and int(normalized_noise.get("total_entries", 0)) > 0
        )
        for method_label, method_dir in reusable_methods:
            output = method_output(root, scenario.sample, method_dir)
            exists = output.exists() and output.stat().st_size > 0
            output_spots = csv_index(output) if exists else set()
            spot_match = output_spots == stage1_st_spots
            postdates_noise = exists and output.stat().st_mtime >= generated_at
            verified = bool(exists and spot_match and (postdates_noise or deterministic_recreation))
            reuse_rows.append(
                {
                    "group": scenario.group,
                    "sample": scenario.sample,
                    "method": method_label.replace("\n", " "),
                    "method_key": method_dir,
                    "output_path": str(output),
                    "output_sha256": sha256(output) if exists else "",
                    "input_sc_sha256": sha256(paths["sc_expression_normalized"]),
                    "input_st_sha256": sha256(paths["st_expression_normalized"]),
                    "coordinates_sha256": sha256(paths["st_coordinates"]),
                    "truth_sha256": sha256(current_truth),
                    "output_exists": exists,
                    "spot_universe_matches_stage1_st": spot_match,
                    "output_spots": len(output_spots),
                    "stage1_st_spots": len(stage1_st_spots),
                    "truth_spots": len(truth_spots),
                    "stage1_st_is_truth_subset": stage1_st_spots.issubset(truth_spots),
                    "output_postdates_noisy_input": postdates_noise,
                    "deterministic_noise_recreation_verified": deterministic_recreation,
                    "provenance_verified": verified,
                    "rerun_required": not verified,
                    "evidence": (
                        "exact noisy-sample output path; runner consumes that sample's Stage1 normalized "
                        "matrices; output spot IDs match the QC-filtered Stage1 ST input; "
                        + ("output postdates noisy input" if postdates_noise else
                           "input was deterministically recreated later from the same source with p=0.1 and seed=42")
                    ),
                }
            )
    input_df = pd.DataFrame(input_rows)
    reuse_df = pd.DataFrame(reuse_rows)
    input_df.to_csv(out_dir / "composite_scnoise10_stage3ab_full_input_audit.csv", index=False)
    reuse_df.to_csv(out_dir / "composite_scnoise10_stage3ab_full_reused_method_audit.csv", index=False)
    bad_inputs = input_df["same_as_no_noise_source"] != input_df["expected_same_as_source"]
    if bad_inputs.any():
        raise RuntimeError("Noisy-input invariants failed; see input audit CSV")
    if not reuse_df["provenance_verified"].all():
        raise RuntimeError("At least one reused method failed provenance verification")
    return input_df, reuse_df


def scenario_manifest(root: Path, out_dir: Path) -> pd.DataFrame:
    rows = []
    for scenario in SCENARIOS:
        info = read_json(truth_path(root, scenario).parent / "sim_info.json")
        cfg = read_dataset_config(root, scenario.sample)
        export = stage1_export_dir(root, scenario.sample, cfg)
        missing = [str(x) for x in info.get("missing_types", [])]
        rows.append(
            {
                "group": scenario.group,
                "sample": scenario.sample,
                "condition": scenario.condition,
                "expected_stage3a_missing_types_audit_only": ";".join(missing),
                "stage3b_target": scenario.stage3b_target,
                "stage3a_filter_mode": "none" if scenario.condition == "control" else "plugin_unknown",
                "filter_scope": "unsupported_all",
                "cell_type_column": "sc_meta" if scenario.condition == "control" else "plugin_type",
                "noise_type": info.get("noise_type"),
                "noise_fraction": info.get("noise_fraction"),
                "noise_seed": info.get("noise_seed"),
                "sc_input_path": str(export / "sc_expression_normalized.csv"),
                "st_input_path": str(export / "st_expression_normalized.csv"),
                "coordinates_path": str(export / "st_coordinates.csv"),
                "truth_path": str(truth_path(root, scenario)),
                "output_path": str(root / "result" / scenario.sample / STAGE4_DIR),
            }
        )
    frame = pd.DataFrame(rows)
    frame.to_csv(out_dir / "composite_scnoise10_stage3ab_full_scenario_manifest.csv", index=False)
    return frame


def execute_pipeline(
    root: Path,
    python: str,
    out_dir: Path,
    resume: bool,
    n_processors: int,
    n_subspots: int,
    mapping_cells_per_spot: int,
) -> None:
    logs = out_dir / "execution_logs"
    for scenario in SCENARIOS:
        cfg = read_dataset_config(root, scenario.sample)
        stage3_summary = result_dir(root, scenario.sample, cfg) / STAGE3_DIR / "stage3_summary.json"
        if not (resume and stage3_summary.exists()):
            run_command(
                [python, "-m", "src.stages.stage3_type_plugin", "--sample", scenario.sample,
                 "--project_root", str(root), "--output_suffix", SUFFIX],
                logs / f"{scenario.sample}_stage3a.log",
                root,
            )
        stage3b_summary = result_dir(root, scenario.sample, cfg) / STAGE3B_DIR / "stage3b_summary.json"
        if not (resume and stage3b_summary.exists()):
            run_command(
                [python, "-m", "src.stages.stage3b_st_unsupported", "--sample", scenario.sample,
                 "--project_root", str(root), "--output_suffix", SUFFIX],
                logs / f"{scenario.sample}_stage3b.log",
                root,
            )
        new_summary = root / "result" / scenario.sample / STAGE4_DIR / "cytospace_output" / "stage4_summary.json"
        if resume and new_summary.exists():
            continue
        command = [
            python, "-m", "src.stages.stage4_cytospace", "--sample", scenario.sample,
            "--project_root", str(root), "--stage3_suffix", SUFFIX, "--stage4_suffix", SUFFIX,
            "--filter_scope", "unsupported_all", "--stage3b_blank_regions",
            "--stage3b_scores_path", str(processed_dir(root, scenario.sample, cfg) / STAGE3B_DIR / "spot_unsupported_scores.csv"),
            "--sc_expr_source", "normalized", "--n_processors", str(n_processors),
            "--n_subspots", str(n_subspots),
            "--mapping_cells_per_spot", str(mapping_cells_per_spot),
        ]
        if scenario.condition == "control":
            command.extend(["--filter_mode", "none", "--cell_type_column", "sc_meta"])
        else:
            command.extend(["--filter_mode", "plugin_unknown", "--cell_type_column", "plugin_type"])
        run_command(command, logs / f"{scenario.sample}_stage4.log", root)


def stage_audits(root: Path, out_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    stage3a_rows = []
    stage3b_rows = []
    resolved_dir = out_dir / "resolved_configs"
    resolved_dir.mkdir(parents=True, exist_ok=True)
    for scenario in SCENARIOS:
        cfg = read_dataset_config(root, scenario.sample)
        relabel_path = processed_dir(root, scenario.sample, cfg) / STAGE3_DIR / "cell_type_relabel.csv"
        relabel = pd.read_csv(relabel_path)
        unknown = relabel["plugin_type"].astype(str).str.startswith("Unknown")
        detected_types = sorted(relabel.loc[unknown, "orig_type"].astype(str).unique())
        summary_path = result_dir(root, scenario.sample, cfg) / STAGE3_DIR / "stage3_summary.json"
        stage3_summary = read_json(summary_path)
        stage4_path = root / "result" / scenario.sample / STAGE4_DIR / "cytospace_output" / "stage4_summary.json"
        stage4 = read_json(stage4_path)
        expected_types = read_json(truth_path(root, scenario).parent / "sim_info.json").get("missing_types", [])
        stage3a_rows.append(
            {
                "group": scenario.group,
                "sample": scenario.sample,
                "condition": scenario.condition,
                "stage3_summary_path": str(summary_path),
                "detected_unsupported_types": ";".join(detected_types),
                "n_cells_marked_unknown": int(unknown.sum()),
                "expected_types_audit_only": ";".join(str(x) for x in expected_types),
                "detection_status": (
                    "not_applicable_control" if scenario.condition == "control" else
                    "detected" if set(map(str, expected_types)).issubset(set(detected_types)) else
                    "partial_detection" if detected_types else "detection_failure"
                ),
                "n_cells_before": stage4.get("n_cells_before_prefilter"),
                "n_cells_after": stage4.get("n_cells_after_prefilter"),
                "stage4_filter_mode": stage4.get("filter_mode"),
                "stage4_filter_scope": stage4.get("filter_scope"),
                "stage4_cell_type_column": stage4.get("cell_type_column"),
                "n_filtered": stage4.get("n_filtered"),
                "truth_filter_enabled": stage4.get("truth_filter_enabled"),
                "oracle_only_filtered_missing": stage4.get("oracle_only_filtered_missing"),
                "forced_stage3_outputs_written": "not_available",
                "stage3_relabel_path": str(relabel_path),
                "stage3_relabel_sha256": sha256(relabel_path),
                "stage4_summary_path": str(stage4_path),
                "stage4_summary_sha256": sha256(stage4_path),
                "execution_log": str(out_dir / "execution_logs" / f"{scenario.sample}_stage3a.log"),
            }
        )
        score_path = processed_dir(root, scenario.sample, cfg) / STAGE3B_DIR / "spot_unsupported_scores.csv"
        scores = pd.read_csv(score_path, index_col=0)
        truth = pd.read_csv(truth_path(root, scenario), index_col=0).apply(pd.to_numeric, errors="coerce").fillna(0)
        scores.index = scores.index.astype(str)
        truth.index = truth.index.astype(str)
        truth = truth.reindex(scores.index)
        fraction = truth[scenario.stage3b_target]
        detected = scores["is_unsupported_region"].astype(str).str.lower().eq("true")
        dominant = truth.idxmax(axis=1).eq(scenario.stage3b_target)
        positive = fraction.gt(0)
        zero = fraction.eq(0)
        stage3b_rows.append(
            {
                "group": scenario.group,
                "sample": scenario.sample,
                "target_type": scenario.stage3b_target,
                "n_spots": len(scores),
                "blank_spots": int(detected.sum()),
                "dominant_precision": float((detected & dominant).sum() / max(int(detected.sum()), 1)),
                "dominant_recall": float((detected & dominant).sum() / max(int(dominant.sum()), 1)),
                "target_positive_recall": float((detected & positive).sum() / max(int(positive.sum()), 1)),
                "target_mass_recall": float(fraction.loc[detected].sum() / max(float(fraction.sum()), np.finfo(float).eps)),
                "zero_target_false_positive_rate": float((detected & zero).sum() / max(int(zero.sum()), 1)),
                "score_path": str(score_path),
                "score_sha256": sha256(score_path),
                "stage3b_enabled_in_stage4": bool(stage4.get("stage3b_blank_regions", {}).get("enabled")),
                "blank_mask_integrated": bool(stage4.get("stage3b_blank_regions", {}).get("applied_before_mapping")),
                "execution_log": str(out_dir / "execution_logs" / f"{scenario.sample}_stage3b.log"),
            }
        )
        stage3b_summary_path = result_dir(root, scenario.sample, cfg) / STAGE3B_DIR / "stage3b_summary.json"
        resolved = {
            "sample": scenario.sample,
            "condition": scenario.condition,
            "dataset_config_path": str(root / "configs" / "datasets" / f"{scenario.sample}.yaml"),
            "dataset_config_sha256": sha256(root / "configs" / "datasets" / f"{scenario.sample}.yaml"),
            "stage3a_output_suffix": SUFFIX,
            "stage3a_params": stage3_summary.get("params", {}),
            "stage3b_output_suffix": SUFFIX,
            "stage3b_config": read_json(stage3b_summary_path).get("config", {}),
            "stage4_output_suffix": SUFFIX,
            "stage4_runtime": {
                key: stage4.get(key) for key in [
                    "filter_mode", "filter_scope", "cell_type_column", "n_filtered",
                    "sampling_sub_spots", "n_subspots", "solver_method", "seed",
                    "mapping_cells_per_spot", "truth_filter_enabled",
                ]
            },
            "stage3b_blank_regions": stage4.get("stage3b_blank_regions", {}),
        }
        (resolved_dir / f"{scenario.sample}.json").write_text(
            json.dumps(resolved, indent=2), encoding="utf-8"
        )
    stage3a_df = pd.DataFrame(stage3a_rows)
    stage3b_df = pd.DataFrame(stage3b_rows)
    stage3a_df.to_csv(out_dir / "composite_scnoise10_stage3ab_full_stage3a_audit.csv", index=False)
    stage3b_df.to_csv(out_dir / "composite_scnoise10_stage3ab_full_stage3b_audit.csv", index=False)
    return stage3a_df, stage3b_df


def evaluate(root: Path, out_dir: Path, reuse_df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    scenario_df, cell_df = build_tables(
        root,
        sample_suffix=NOISE_SUFFIX,
        methods=METHODS,
        scenarios=COMPOSITE_SCENARIOS,
        route2_stage4_dir=STAGE4_DIR,
        reward_correct_stage3b_abstention=True,
        abstention_truth_rule="target_dominant",
    )
    scenario_df["method"] = scenario_df["method"].astype(str)
    cell_df["method"] = cell_df["method"].astype(str)
    reuse_lookup = {(row["sample"], row["method_key"]): row for _, row in reuse_df.iterrows()}
    method_key = {label: key for label, key in METHODS}
    scenario_df["method_key"] = scenario_df["method"].map(method_key)
    scenario_df["prediction_path"] = scenario_df.apply(
        lambda row: str(method_output(root, row["sample"], row["method_key"])), axis=1
    )
    scenario_df["reused_mapping"] = scenario_df["method_key"].ne("cytospace_route2")
    scenario_df["provenance_verified"] = scenario_df.apply(
        lambda row: True if row["method_key"] == "cytospace_route2" else bool(
            reuse_lookup[(row["sample"], row["method_key"])]["provenance_verified"]
        ), axis=1
    )
    prefix = "composite_scnoise10_stage3ab_full"
    scenario_df.to_csv(out_dir / f"{prefix}_source_values.csv", index=False)
    cell_df.to_csv(out_dir / f"{prefix}_cell_type_values.csv", index=False)
    summary = scenario_df.groupby("method", sort=False)["composition_recovery"].agg(["count", "mean", "median", "std"]).reset_index()
    summary.to_csv(out_dir / f"{prefix}_method_summary.csv", index=False)
    wide = scenario_df.pivot(index=["group", "sample"], columns="method", values="composition_recovery").reset_index()
    wide["svtuner_minus_cytospace"] = wide["SVTuner"] - wide["CytoSPACE"]
    wide["svtuner_win"] = wide["svtuner_minus_cytospace"] > 0
    wide.to_csv(out_dir / f"{prefix}_paired_comparison.csv", index=False)
    new_blank = scenario_df[scenario_df["method"] == "SVTuner"].copy()
    new_blank["blank_precision"] = new_blank["correct_abstention_spots"] / new_blank["predicted_blank_spots"].replace(0, np.nan)
    new_blank["blank_recall"] = new_blank["correct_abstention_spots"] / new_blank["truth_unsupported_spots"].replace(0, np.nan)
    blank_columns = [
        "group", "sample", "predicted_blank_spots", "correct_abstention_spots",
        "incorrect_abstention_spots", "truth_unsupported_spots", "blank_precision", "blank_recall",
    ]
    new_blank[blank_columns].to_csv(out_dir / f"{prefix}_blank_summary.csv", index=False)
    plot(
        scenario_df,
        out_dir / "fig1d_scnoise10_stage3ab_full.png",
        out_dir / "fig1d_scnoise10_stage3ab_full.pdf",
        "Spatial cell-type composition recovery, 10% sc-reference noise",
        "Abstention-aware recovery score",
    )
    cy_mean = float(summary.loc[summary["method"] == "CytoSPACE", "mean"].iloc[0])
    sv_mean = float(summary.loc[summary["method"] == "SVTuner", "mean"].iloc[0])
    result = {
        "svtuner_mean": sv_mean,
        "cytospace_mean": cy_mean,
        "absolute_improvement": sv_mean - cy_mean,
        "relative_improvement_percent": 100 * (sv_mean - cy_mean) / cy_mean,
        "paired_wins": int(wide["svtuner_win"].sum()),
        "paired_total": int(len(wide)),
        "scoring_rule": "whole-space target-dominant abstention-aware composition overlap",
        "blank_precision": float(new_blank["correct_abstention_spots"].sum() / max(int(new_blank["predicted_blank_spots"].sum()), 1)),
        "blank_recall": float(new_blank["correct_abstention_spots"].sum() / max(int(new_blank["truth_unsupported_spots"].sum()), 1)),
    }
    method_comparisons = []
    for method in [label for label, _ in METHODS if label != "SVTuner"]:
        delta = wide["SVTuner"] - wide[method]
        method_comparisons.append({
            "comparison_method": method.replace("\n", " "),
            "mean_difference": float(delta.mean()),
            "wins": int((delta > 0).sum()),
            "losses": int((delta < 0).sum()),
            "ties": int((delta == 0).sum()),
        })
    pd.DataFrame(method_comparisons).to_csv(out_dir / f"{prefix}_svtuner_method_comparisons.csv", index=False)
    result["comparisons_to_other_methods"] = method_comparisons
    return scenario_df, result


def main() -> int:
    args = parse_args()
    # Do not resolve a subst drive back to its non-ASCII physical path.
    root = Path(args.project_root).absolute()
    out_dir = root / OUT_REL
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = scenario_manifest(root, out_dir)
    input_df, reuse_df = input_and_reuse_audit(root, out_dir)
    if not args.execute:
        print(f"[AUDIT PASS] {len(reuse_df)} reused outputs verified. Use --execute to run.")
        return 0
    execute_pipeline(
        root,
        args.python,
        out_dir,
        args.resume,
        args.n_processors,
        args.n_subspots,
        args.mapping_cells_per_spot,
    )
    stage3a_df, stage3b_df = stage_audits(root, out_dir)
    scenario_df, metrics = evaluate(root, out_dir, reuse_df)
    controls = stage3a_df["condition"].eq("control")
    missing = ~controls
    guardrails = {
        "standalone_from_retired_stage3b_only_outputs": True,
        "all_reused_method_provenance_verified": bool(reuse_df["provenance_verified"].all()),
        "controls_filter_none": bool((stage3a_df.loc[controls, "stage4_filter_mode"] == "none").all()),
        "controls_n_filtered_zero": bool((stage3a_df.loc[controls, "n_filtered"] == 0).all()),
        "missing_filter_plugin_unknown": bool((stage3a_df.loc[missing, "stage4_filter_mode"] == "plugin_unknown").all()),
        "missing_stage3a_outputs_and_logs_present": bool(
            stage3a_df.loc[missing, "stage3_summary_path"].map(lambda x: Path(x).exists()).all()
            and stage3a_df.loc[missing, "execution_log"].map(lambda x: Path(x).exists()).all()
        ),
        "missing_n_filtered_positive_quality_check": bool((stage3a_df.loc[missing, "n_filtered"] > 0).all()),
        "filter_scope_unsupported_all": bool((stage3a_df["stage4_filter_scope"] == "unsupported_all").all()),
        "truth_filter_disabled": bool((stage3a_df["truth_filter_enabled"] == False).all()),  # noqa: E712
        "oracle_filtering_disabled": bool(
            (stage3a_df["truth_filter_enabled"] == False).all()  # noqa: E712
            and (stage3a_df["stage4_filter_scope"] != "missing_only").all()
        ),
        "all_expected_stage3a_types_detected_quality_check": bool(
            stage3a_df.loc[missing, "detection_status"].eq("detected").all()
        ),
        "nine_scenarios_seven_methods": len(scenario_df) == 63,
        "stage3b_all_scenarios_present": len(stage3b_df) == 9,
        "stage3b_enabled_and_integrated": bool(
            stage3b_df["stage3b_enabled_in_stage4"].all()
            and stage3b_df["blank_mask_integrated"].all()
        ),
    }
    critical_keys = [
        "standalone_from_retired_stage3b_only_outputs", "all_reused_method_provenance_verified", "controls_filter_none",
        "controls_n_filtered_zero", "missing_filter_plugin_unknown", "missing_stage3a_outputs_and_logs_present",
        "filter_scope_unsupported_all", "truth_filter_disabled", "oracle_filtering_disabled",
        "nine_scenarios_seven_methods", "stage3b_all_scenarios_present", "stage3b_enabled_and_integrated",
    ]
    quality_keys = ["missing_n_filtered_positive_quality_check", "all_expected_stage3a_types_detected_quality_check"]
    if not all(guardrails[key] for key in critical_keys):
        guardrails["decision"] = "FAIL"
    elif not all(guardrails[key] for key in quality_keys):
        guardrails["decision"] = "PARTIAL"
    else:
        guardrails["decision"] = "PASS"
    (out_dir / "composite_scnoise10_stage3ab_full_guardrails.json").write_text(json.dumps(guardrails, indent=2), encoding="utf-8")
    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "decision": guardrails["decision"],
        "route": "Stage3A unsupported-reference filtering + Stage3B blanking + CytoSPACE",
        "noise": "10% independent entry-wise within-gene replacement, seed 42",
        "scenario_count": len(manifest),
        "input_audit_rows": len(input_df),
        "n_subspots": args.n_subspots,
        "mapping_cells_per_spot": args.mapping_cells_per_spot,
        **metrics,
    }
    (out_dir / "composite_scnoise10_stage3ab_full_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0 if guardrails["decision"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
