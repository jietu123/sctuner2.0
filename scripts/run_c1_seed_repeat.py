#!/usr/bin/env python
"""Run one existing composite scenario with one new simulation seed."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.plot_method_comparison_composition_recovery import (  # noqa: E402
    COMPOSITE_SCENARIOS,
    _abstention_aware_composition_recovery,
    _composition_recovery,
    _load_fraction,
)

SCENARIO_IDS = ["BC-1", "BC-2", "BC-3", "Lung-1", "Lung-2", "Lung-3", "Brain-1", "Brain-2", "Brain-3"]
ALGORITHM_SEED = 42


def command_for_step(python: str, root: Path, group: str, source: str, target: str, info: dict) -> list[str]:
    common = ["--project_root", str(root), "--source_sample", source, "--target_sample", target, "--sim_group", group]
    kind = str(info["simulation_type"])
    if kind.startswith("spatial_"):
        params = info["params"]
        cmd = [python, str(root / "scripts/generate_real_brca_clustered_sim.py"), *common,
               "--simulation_type_label", kind, "--centers_per_type", str(params["centers_per_type"]),
               "--distance_scale", str(params["distance_scale"]), "--temperature", str(params["temperature"]),
               "--noise_sigma", str(params["noise_sigma"]), "--mix_alpha", str(params["mix_alpha"]),
               "--depth_scale", str(params["depth_scale"])]
        for old, value in (params.get("replace_cell_type") or {}).items():
            cmd += ["--replace_cell_type", f"{old}={value['replacement_type']}"]
        if params.get("exclude_cell_types"):
            cmd += ["--exclude_cell_types", *params["exclude_cell_types"]]
        return cmd
    if kind.startswith("missing_type_from_existing_sim"):
        cmd = [python, str(root / "scripts/generate_missing_type_from_sim.py"), *common,
               "--drop_cell_type", str(info["missing_type"]), "--drop_fraction", str(info["drop_fraction"])]
        if info.get("replacement_type"):
            cmd += ["--replacement_cell_type", str(info["replacement_type"])]
        return cmd
    if kind == "purified_target_region_from_existing_sim":
        return [python, str(root / "scripts/purify_sim_target_region.py"), *common,
                "--target_type", str(info["target_type"]), "--threshold", str(info["threshold"]),
                "--replacement_type", str(info["replacement_type"]),
                "--marker_boost_top_n", str(info["marker_boost_top_n"]),
                "--marker_boost_factor", str(info["marker_boost_factor"]), "--hvg_nfeatures", "0"]
    if kind == "st_only_type_missing_from_sc_reference":
        cmd = [python, str(root / "scripts/generate_st_only_reference_dropout_from_sim.py"), *common]
        for cell_type in info["sc_reference_drop_types"]:
            cmd += ["--drop_cell_type", str(cell_type)]
        return cmd
    raise ValueError(f"Unsupported historical simulation step: {kind}")


def lineage(root: Path, group: str, final_sample: str) -> list[tuple[str, dict]]:
    steps = []
    sample = final_sample
    while True:
        path = root / "data/sim" / group / sample / "sim_info.json"
        if not path.exists():
            break
        info = json.loads(path.read_text(encoding="utf-8"))
        steps.append((sample, info))
        sample = str(info["source_sample"])
    steps.reverse()
    return steps


def run(cmd: list[str], root: Path) -> None:
    print("[C1]", subprocess.list2cmdline(cmd), flush=True)
    subprocess.run(cmd, cwd=root, check=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scenario", choices=SCENARIO_IDS, required=True)
    parser.add_argument("--simulation-seed", type=int, required=True)
    parser.add_argument("--project_root", default=".")
    parser.add_argument("--python", default=sys.executable)
    args = parser.parse_args()
    root = Path(args.project_root).absolute()
    index = SCENARIO_IDS.index(args.scenario)
    group, historical_final = COMPOSITE_SCENARIOS[index]
    historical_steps = lineage(root, group, historical_final)
    if not historical_steps:
        raise FileNotFoundError(f"No historical generation lineage for {historical_final}")

    tag = f"c1_{args.scenario.lower().replace('-', '')}_s{args.simulation_seed}"
    renamed = {sample: f"{tag}_{i + 1}" for i, (sample, _) in enumerate(historical_steps)}
    source = str(historical_steps[0][1]["source_sample"])
    generated = []
    for old_sample, info in historical_steps:
        target = renamed[old_sample]
        mapped_source = renamed.get(str(info["source_sample"]), source)
        cmd = command_for_step(args.python, root, group, mapped_source, target, info)
        if not str(info["simulation_type"]).startswith("st_only_"):
            cmd += ["--seed", str(args.simulation_seed)]
        if (root / "data/sim" / group / target).exists():
            raise FileExistsError(f"Refusing to overwrite existing sample: {target}")
        run(cmd, root)
        generated.append(target)

    base_sample, sample = generated[0], generated[-1]
    run([args.python, str(root / "scripts/prepare_stage1_from_sim_source.py"), "--sample", base_sample,
         "--project_root", str(root), "--rebuild_sc_from_raw"], root)
    run([args.python, str(root / "scripts/prepare_stage1_from_sim_source.py"), "--sample", sample,
         "--project_root", str(root)], root)
    run([args.python, "-m", "src.stages.stage3_type_plugin", "--sample", sample,
         "--project_root", str(root), "--sc_expr_source", "normalized"], root)
    run([args.python, "-m", "src.stages.stage3b_st_unsupported", "--sample", sample,
         "--project_root", str(root), "--random_seed", str(ALGORITHM_SEED)], root)

    stage4 = [args.python, "-m", "src.stages.stage4_cytospace", "--sample", sample,
              "--project_root", str(root), "--seed", str(ALGORITHM_SEED), "--n_processors", "1",
              "--n_subspots", "800", "--sc_expr_source", "normalized"]
    run([*stage4, "--stage4_suffix", "_baseline", "--filter_mode", "none",
         "--cell_type_column", "sc_meta"], root)
    control = index % 3 == 0
    run([*stage4, "--stage4_suffix", "_stage3b_blank", "--filter_mode", "none" if control else "plugin_unknown",
         "--cell_type_column", "sc_meta" if control else "plugin_type", "--filter_scope", "unsupported_all",
         "--stage3b_blank_regions"], root)

    truth = _load_fraction(root / "data/sim" / group / sample / "sim_truth_spot_type_fraction.csv")
    result = root / "result" / sample
    baseline = _load_fraction(result / "stage4_cytospace_baseline/cytospace_output/fractional_abundances_by_spot.csv")
    svtuner_dir = result / "stage4_cytospace_stage3b_blank/cytospace_output"
    svtuner = _load_fraction(svtuner_dir / "fractional_abundances_by_spot.csv")
    baseline_score, _ = _composition_recovery(baseline, truth)
    blank_ids = set(pd.read_csv(svtuner_dir / "stage3b_blank_spots.csv", index_col=0).index.astype(str))
    drop_types = json.loads((root / "data/sim" / group / sample / "sim_info.json").read_text(encoding="utf-8"))["sc_reference_drop_types"]
    svtuner_score, _ = _abstention_aware_composition_recovery(
        svtuner, truth, blank_ids, drop_types, truth_rule="target_dominant"
    )
    row = pd.DataFrame([{"scenario_id": args.scenario, "simulation_seed": args.simulation_seed,
                         "sample": sample, "cytospace_score": baseline_score, "svtuner_score": svtuner_score}])
    score_path = root / "result/c1_seed_scores.csv"
    if score_path.exists():
        old = pd.read_csv(score_path)
        if ((old["scenario_id"] == args.scenario) & (old["simulation_seed"] == args.simulation_seed)).any():
            raise FileExistsError(f"Score already exists for {args.scenario}, seed {args.simulation_seed}")
        row = pd.concat([old, row], ignore_index=True)
    row.to_csv(score_path, index=False)
    print(row.tail(1).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
