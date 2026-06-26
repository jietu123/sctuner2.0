#!/usr/bin/env python
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path


SCENARIOS = [
    ("real_brca", "real_brca7_candidate_stable_control"),
    ("real_brca", "real_brca7_candidate_stable_control_missing_epithelial_cells"),
    ("real_brca", "real_brca7_candidate_stable_control_missing_epithelial_cells_pcs"),
    ("human_lung_5loc", "human_lung_5loc_fine9_clustered_sim"),
    ("human_lung_5loc", "human_lung_5loc_fine9_clustered_sim_missing_at2"),
    ("human_lung_5loc", "human_lung_5loc_fine9_clustered_sim_missing_at2_fibroblast"),
    ("mouse_brain_refined", "mouse_brain_refined7_balanced_clustered_sim"),
    (
        "mouse_brain_refined",
        "mouse_brain_refined7_balanced_clustered_sim_missing_micro_fill_ext_l56",
    ),
    (
        "mouse_brain_refined",
        "mouse_brain_refined7_balanced_clustered_sim_missing_micro_oligo_2_fill_ext_l56",
    ),
]

COMPOSITE_SCENARIOS = [
    ("real_brca", "real_brca7_endothelial_marker_control_sc_missing_endothelial_cells"),
    (
        "real_brca",
        "real_brca7_endothelial_marker_missing_epithelial_cells_sc_missing_endothelial_cells",
    ),
    (
        "real_brca",
        "real_brca7_endothelial_marker_missing_epithelial_cells_pcs_sc_missing_endothelial_cells",
    ),
    ("human_lung_5loc", "human_lung_5loc_fine9_clustered_sim_sc_missing_b_cell"),
    (
        "human_lung_5loc",
        "human_lung_5loc_fine9_clustered_sim_missing_at2_sc_missing_b_cell",
    ),
    (
        "human_lung_5loc",
        "human_lung_5loc_fine9_clustered_sim_missing_at2_fibroblast_sc_missing_b_cell",
    ),
    ("mouse_brain_refined", "mouse_brain_refined7_balanced_clustered_sim_sc_missing_ext_l56"),
    (
        "mouse_brain_refined",
        "mouse_brain_refined7_balanced_clustered_sim_missing_micro_fill_inh_pvalb_sc_missing_ext_l56",
    ),
    (
        "mouse_brain_refined",
        "mouse_brain_refined7_balanced_clustered_sim_missing_micro_oligo_2_fill_inh_pvalb_sc_missing_ext_l56",
    ),
]

METHODS = [
    "tangram_all",
    "tangram_marker",
    "novosparc",
    "spaotsc",
    "celltrek",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Run retained external mapping methods on the nine 0% noise scenarios."
    )
    p.add_argument("--project_root", default=".")
    p.add_argument(
        "--python",
        default=sys.executable,
        help="Python executable containing Tangram, novoSpaRc, and SpaOTsc.",
    )
    p.add_argument(
        "--methods",
        default=",".join(METHODS),
        help="Comma-separated methods to run.",
    )
    p.add_argument(
        "--sample_suffix",
        default="",
        help="Suffix appended to all scenario ids, for example _scnoise10.",
    )
    p.add_argument(
        "--scenario_preset",
        choices=["standard", "composite"],
        default="standard",
        help="Scenario set to run.",
    )
    p.add_argument("--force", action="store_true")
    return p.parse_args()


def _output_path(root: Path, sample: str, method: str) -> Path:
    return (
        root
        / "result"
        / sample
        / "stage4_mapping"
        / method
        / "spot_type_fraction.csv"
    )


def _command(
    python: str, root: Path, group: str, sample: str, method: str
) -> list[str]:
    common = ["--project_root", str(root), "--group", group, "--sample", sample]
    if method.startswith("tangram_"):
        mode = "all" if method == "tangram_all" else "marker"
        return [
            python,
            str(root / "scripts" / "run_tangram_marker_mapping.py"),
            *common,
            "--gene_mode",
            mode,
            "--top_n_marker",
            "50",
            "--num_epochs",
            "200",
            "--device",
            "cpu",
        ]
    if method == "celltrek":
        return [
            python,
            str(root / "scripts" / "run_celltrek_mapping.py"),
            *common,
            "--max_genes",
            "2000",
            "--n_pcs",
            "30",
            "--ntree",
            "500",
        ]
    return [
        python,
        str(root / "scripts" / "run_ot_mapping.py"),
        "--method",
        method,
        *common,
        "--max_genes",
        "500",
        "--n_pcs",
        "30",
    ]


def main() -> int:
    args = parse_args()
    root = Path(args.project_root).resolve()
    selected = [x.strip() for x in args.methods.split(",") if x.strip()]
    unknown = sorted(set(selected) - set(METHODS))
    if unknown:
        raise ValueError(f"Unknown methods: {unknown}")

    env = os.environ.copy()
    env.pop("PYTHONNOUSERSITE", None)
    env.setdefault("OMP_NUM_THREADS", "4")
    env.setdefault("MKL_NUM_THREADS", "4")

    failures: list[tuple[str, str, int]] = []
    completed = 0
    skipped = 0
    started = time.perf_counter()
    scenarios = COMPOSITE_SCENARIOS if args.scenario_preset == "composite" else SCENARIOS
    for method in selected:
        for group, base_sample in scenarios:
            sample = f"{base_sample}{args.sample_suffix}"
            output = _output_path(root, sample, method)
            if output.exists() and not args.force:
                print(f"[SKIP] {method}: {sample}", flush=True)
                skipped += 1
                continue
            print(f"[RUN] {method}: {sample}", flush=True)
            result = subprocess.run(
                _command(args.python, root, group, sample, method),
                cwd=root,
                env=env,
                check=False,
            )
            if result.returncode:
                failures.append((method, sample, result.returncode))
                print(
                    f"[FAIL] {method}: {sample} exit={result.returncode}",
                    flush=True,
                )
            else:
                completed += 1
    elapsed = time.perf_counter() - started
    print(
        f"[SUMMARY] completed={completed} skipped={skipped} "
        f"failed={len(failures)} elapsed_seconds={elapsed:.1f}",
        flush=True,
    )
    for method, sample, code in failures:
        print(f"[FAILED_ITEM] {method}\t{sample}\t{code}", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
