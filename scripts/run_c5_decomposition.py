from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PYTHON = Path(r"E:\ANACONDA\python.exe")
CYTOSPACE_PYTHON = Path(r"E:\ANACONDA\envs\cytospace_v1.1.0_py310\python.exe")
OUT = ROOT / "result" / "c5_decomposition"
LOGS = OUT / "logs"
STANDALONE_CFG = OUT / "c5_standalone_stage3b_dataset_config.yaml"
SEQUENTIAL_CFG = (
    ROOT / "result" / "c5_stage3a_decomposition" / "c5_sequential_stage3b_dataset_config.yaml"
)
C5_STAGE3A = (
    ROOT / "result" / "c5_stage3a_decomposition" / "c5_stage3a_exclusion_by_dataset.csv"
)
PILOTS = {
    "adult_mouse_kidney_real_profile_mask_endo",
    "highres_humanbreastcancerpatient1_profile_mask_monocytes_and_macrophages",
}

STAGE4_CORE_REQUIRED = (
    "assigned_locations.csv",
    "cell_assignment.csv",
    "cell_type_assignments_by_spot.csv",
    "fractional_abundances_by_spot.csv",
    "stage4_summary.json",
    "log.txt",
)


def source_root(dataset: str, resolution: str) -> Path:
    if resolution == "spot-resolution":
        return ROOT / "data" / "processed" / "low_resolution_experiments" / dataset
    return ROOT / "data" / "processed" / dataset


def stage1_export(base: Path) -> Path:
    return base / "stage1_preprocess" / "exported"


def hardlink(src: Path, dst: Path) -> None:
    if not src.is_file():
        raise FileNotFoundError(src)
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        if os.path.samefile(src, dst):
            return
        raise FileExistsError(f"Existing non-matching destination: {dst}")
    os.link(src, dst)


def complete_stage3b(mask: Path, summary: Path) -> bool:
    if not (mask.is_file() and summary.is_file()):
        return False
    try:
        scores = pd.read_csv(mask)
        json.loads(summary.read_text(encoding="utf-8"))
        return "is_unsupported_region" in scores.columns and len(scores) > 0
    except Exception:
        return False


def complete_stage4(output_dir: Path, requires_stage3b_mask: bool = True) -> bool:
    required = STAGE4_CORE_REQUIRED + (("stage3b_blank_spots.csv",) if requires_stage3b_mask else ())
    return all((output_dir / name).is_file() for name in required)


def run_logged(command: list[str], log_path: Path, env: dict[str, str] | None = None) -> tuple[bool, str]:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as handle:
        handle.write("COMMAND: " + subprocess.list2cmdline(command) + "\n\n")
        handle.flush()
        proc = subprocess.run(
            command,
            cwd=ROOT,
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
            env=env,
            check=False,
        )
    if proc.returncode == 0:
        return True, ""
    tail = log_path.read_text(encoding="utf-8", errors="replace").splitlines()[-15:]
    return False, " | ".join(tail)


def prepare_reference(row: pd.Series, mode: str) -> tuple[Path, int, int, dict[str, int]]:
    dataset = str(row["dataset"])
    original = source_root(dataset, str(row["resolution"]))
    original_export = stage1_export(original)
    group = "c5_stage3b_standalone" if mode == "standalone" else "c5_stage3a_stage3b_sequential"
    target_root = ROOT / "data" / "processed" / group / dataset
    target_export = stage1_export(target_root)

    for name in ("sc_expression_normalized.csv", "st_expression_normalized.csv", "st_coordinates.csv"):
        hardlink(original_export / name, target_export / name)

    original_meta = original_export / "sc_metadata.csv"
    if mode == "standalone":
        hardlink(original_meta, target_export / "sc_metadata.csv")
        metadata = pd.read_csv(original_meta)
        return target_root, len(metadata), 0, {}

    relabel_path = original / "stage3_typematch" / "cell_type_relabel.csv"
    if not relabel_path.is_file():
        raise FileNotFoundError(relabel_path)
    relabel = pd.read_csv(relabel_path)
    metadata = pd.read_csv(original_meta)
    if not {"cell_id", "action", "orig_type"}.issubset(relabel.columns):
        raise ValueError(f"Unexpected relabel schema: {relabel_path}")
    dropped = relabel.loc[relabel["action"].astype(str).eq("Dropped")].copy()
    keep_ids = set(relabel.loc[~relabel["action"].astype(str).eq("Dropped"), "cell_id"].astype(str))
    retained = metadata.loc[metadata["cell_id"].astype(str).isin(keep_ids)].copy()
    if len(retained) != len(keep_ids):
        raise ValueError(
            f"Retained metadata mismatch for {dataset}: rows={len(retained)}, ids={len(keep_ids)}"
        )
    target_meta = target_export / "sc_metadata.csv"
    if target_meta.exists():
        prior = pd.read_csv(target_meta)
        if set(prior["cell_id"].astype(str)) != set(retained["cell_id"].astype(str)):
            raise ValueError(f"Existing retained metadata differs: {target_meta}")
    else:
        target_meta.parent.mkdir(parents=True, exist_ok=True)
        retained.to_csv(target_meta, index=False)

    dropped_by_type = {
        str(k): int(v) for k, v in dropped["orig_type"].astype(str).value_counts().items()
    }
    info = {
        "dataset": dataset,
        "basis": "Stage3A action != Dropped",
        "source_stage1": str(original_export),
        "source_relabel": str(relabel_path),
        "full_reference_cells": int(len(metadata)),
        "retained_reference_cells": int(len(retained)),
        "dropped_reference_cells": int(len(dropped)),
        "dropped_populations": dropped_by_type,
        "large_expression_files": "hardlinked from original Stage1",
    }
    info_path = target_root / "stage1_preprocess" / "retained_reference_info.json"
    info_path.parent.mkdir(parents=True, exist_ok=True)
    info_path.write_text(json.dumps(info, indent=2, ensure_ascii=False), encoding="utf-8")
    return target_root, len(retained), len(dropped), dropped_by_type


def stage3b_paths(dataset: str, mode: str, pilot: bool) -> tuple[Path, Path]:
    if mode == "standalone" and pilot:
        if dataset == "adult_mouse_kidney_real_profile_mask_endo":
            data_base = ROOT / "data" / "processed" / "low_resolution_experiments" / dataset
            result_base = ROOT / "result" / "low_resolution_experiments" / dataset
        else:
            data_base = ROOT / "data" / "processed" / dataset
            result_base = ROOT / "result" / dataset
        suffix = "stage3b_st_unsupported_c5_pilot"
    elif mode == "sequential" and pilot:
        data_base = ROOT / "data" / "processed" / "c5_stage3a_stage3b_sequential" / dataset
        result_base = ROOT / "result" / "c5_stage3a_stage3b_sequential" / dataset
        suffix = "stage3b_st_unsupported_sequential_full"
    elif mode == "standalone":
        data_base = ROOT / "data" / "processed" / "c5_stage3b_standalone" / dataset
        result_base = ROOT / "result" / "c5_stage3b_standalone" / dataset
        suffix = "stage3b_st_unsupported_c5_standalone"
    else:
        data_base = ROOT / "data" / "processed" / "c5_stage3a_stage3b_sequential" / dataset
        result_base = ROOT / "result" / "c5_stage3a_stage3b_sequential" / dataset
        suffix = "stage3b_st_unsupported_sequential_full"
    return data_base / suffix / "spot_unsupported_scores.csv", result_base / suffix / "stage3b_summary.json"


def run_stage3b(dataset: str, mode: str, pilot: bool) -> tuple[bool, Path, Path, str]:
    mask, summary = stage3b_paths(dataset, mode, pilot)
    if complete_stage3b(mask, summary):
        return True, mask, summary, "reused"
    if pilot:
        return False, mask, summary, "pilot output incomplete"
    cfg = STANDALONE_CFG if mode == "standalone" else SEQUENTIAL_CFG
    suffix = "_c5_standalone" if mode == "standalone" else "_sequential_full"
    command = [
        str(PYTHON), "-m", "src.stages.stage3b_st_unsupported",
        "--sample", dataset,
        "--project_root", str(ROOT),
        "--dataset_config", str(cfg),
        "--output_suffix", suffix,
    ]
    ok, error = run_logged(command, LOGS / f"{dataset}.{mode}.stage3b.log")
    ok = ok and complete_stage3b(mask, summary)
    return ok, mask, summary, error if not ok else "completed"


def score_mask(path: Path) -> tuple[set[str], int]:
    frame = pd.read_csv(path)
    spot_col = next(c for c in ("spot_id", "spot", "SpotID", "spot_name") if c in frame.columns)
    flag = frame["is_unsupported_region"].astype(str).str.lower().isin(("true", "1"))
    return set(frame.loc[flag, spot_col].astype(str)), len(frame)


def prepare_mapping_project(dataset: str, mode: str, row: pd.Series) -> Path:
    project = OUT / "mapping_projects" / mode
    cfg_dir = project / "configs" / "datasets"
    cfg_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT / "configs" / "project_config.yaml", project / "configs" / "project_config.yaml")
    alias_config = ROOT / "configs" / "type_aliases.yaml"
    if alias_config.is_file():
        shutil.copy2(alias_config, project / "configs" / "type_aliases.yaml")
    shutil.copy2(ROOT / "configs" / "datasets" / f"{dataset}.yaml", cfg_dir / f"{dataset}.yaml")

    source = source_root(dataset, str(row["resolution"]))
    src_export = stage1_export(source)
    dst_export = project / "data" / "processed" / dataset / "stage1_preprocess" / "exported"
    for name in ("sc_expression_normalized.csv", "st_expression_normalized.csv", "st_coordinates.csv"):
        hardlink(src_export / name, dst_export / name)
    if mode == "standalone":
        metadata = src_export / "sc_metadata.csv"
    else:
        metadata = (
            ROOT / "data" / "processed" / "c5_stage3a_stage3b_sequential" / dataset
            / "stage1_preprocess" / "exported" / "sc_metadata.csv"
        )
    hardlink(metadata, dst_export / "sc_metadata.csv")
    return project


def run_mapping(dataset: str, mode: str, mask: Path, row: pd.Series) -> tuple[bool, Path, str]:
    if dataset == "highres_humanbreastcancerpatient1_profile_mask_monocytes_and_macrophages":
        suffix = "stage4_cytospace_c5_stage3b_only" if mode == "standalone" else "stage4_cytospace_c5_sequential_full"
        output = ROOT / "result" / dataset / suffix / "cytospace_output"
        return complete_stage4(output), output, "reused"

    project = prepare_mapping_project(dataset, mode, row)
    suffix_arg = "_c5_stage3b_only" if mode == "standalone" else "_c5_sequential_full"
    output = project / "result" / dataset / f"stage4_cytospace{suffix_arg}" / "cytospace_output"
    if complete_stage4(output):
        return True, output, "reused"
    command = [
        str(CYTOSPACE_PYTHON), "-m", "src.stages.stage4_cytospace",
        "--sample", dataset,
        "--project_root", str(project),
        "--stage4_suffix", suffix_arg,
        "--filter_mode", "none",
        "--cell_type_column", "sc_meta",
        "--stage3b_blank_regions",
        "--stage3b_scores_path", str(mask.resolve()),
        "--mapping_cells_per_spot", "1",
        "--no_sampling_sub_spots",
        "--seed", "42",
    ]
    env = os.environ.copy()
    additions = [str(ROOT), str(ROOT / "external" / "cytospace")]
    env["PYTHONPATH"] = os.pathsep.join(additions + ([env["PYTHONPATH"]] if env.get("PYTHONPATH") else []))
    ok, error = run_logged(command, LOGS / f"{dataset}.{mode}.stage4.log", env=env)
    ok = ok and complete_stage4(output)
    return ok, output, error if not ok else "completed"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    source_table = pd.read_csv(C5_STAGE3A)
    decomposition: list[dict] = []
    runtime: dict[str, dict] = {}

    for index, row in source_table.iterrows():
        dataset = str(row["dataset"])
        pilot = dataset in PILOTS
        print(f"[{index + 1:02d}/15] {dataset}: preparing references", flush=True)
        state: dict = {"errors": []}
        try:
            if pilot:
                # Validate the already prepared sequential reference; do not rewrite pilot assets.
                seq_info = (
                    ROOT / "data" / "processed" / "c5_stage3a_stage3b_sequential" / dataset
                    / "stage1_preprocess" / "retained_reference_info.json"
                )
                info = json.loads(seq_info.read_text(encoding="utf-8"))
                state["retained_cells"] = int(info["retained_reference_cells"])
                state["dropped_cells"] = int(
                    info.get("dropped_reference_cells", info.get("dropped_cells", 0))
                )
                state["dropped_populations"] = info["dropped_populations"]
            else:
                prepare_reference(row, "standalone")
                _, retained, dropped, dropped_types = prepare_reference(row, "sequential")
                state["retained_cells"] = int(retained)
                state["dropped_cells"] = int(dropped)
                state["dropped_populations"] = dropped_types
        except Exception as exc:
            state["errors"].append(f"reference preparation: {exc}")

        masks: dict[str, Path] = {}
        for mode in ("standalone", "sequential"):
            try:
                ok, mask, summary, note = run_stage3b(dataset, mode, pilot)
                state[f"{mode}_ok"] = bool(ok)
                state[f"{mode}_mask"] = str(mask)
                state[f"{mode}_summary"] = str(summary)
                state[f"{mode}_note"] = note
                if ok:
                    masks[mode] = mask
            except Exception as exc:
                state[f"{mode}_ok"] = False
                state["errors"].append(f"{mode} Stage3B: {exc}")

        record = {
            "dataset": dataset,
            "target": row["target"],
            "resolution": row["resolution"],
            "capacity": int(row["capacity"]),
            "total_units": pd.NA,
            "target_exclusion_fraction": float(row["target_exclusion_fraction"]),
            "standalone_withheld": pd.NA,
            "standalone_fraction": pd.NA,
            "sequential_withheld": pd.NA,
            "sequential_fraction": pd.NA,
            "delta_withheld": pd.NA,
            "delta_fraction": pd.NA,
            "mask_intersection": pd.NA,
            "mask_union": pd.NA,
            "jaccard": pd.NA,
            "standalone_mask_path": state.get("standalone_mask", ""),
            "sequential_mask_path": state.get("sequential_mask", ""),
            "stage3b_success": False,
            "error": "",
        }
        if set(masks) == {"standalone", "sequential"}:
            stand, n1 = score_mask(masks["standalone"])
            seq, n2 = score_mask(masks["sequential"])
            if n1 != n2:
                state["errors"].append(f"spot count mismatch: {n1} vs {n2}")
            else:
                inter, union = stand & seq, stand | seq
                record.update(
                    total_units=n1,
                    standalone_withheld=len(stand),
                    standalone_fraction=len(stand) / n1,
                    sequential_withheld=len(seq),
                    sequential_fraction=len(seq) / n1,
                    delta_withheld=len(seq) - len(stand),
                    delta_fraction=(len(seq) - len(stand)) / n1,
                    mask_intersection=len(inter),
                    mask_union=len(union),
                    jaccard=(len(inter) / len(union)) if union else 1.0,
                    stage3b_success=True,
                )
        record["error"] = " | ".join(state["errors"])
        decomposition.append(record)
        runtime[dataset] = state
        print(
            f"[{index + 1:02d}/15] {dataset}: "
            f"standalone={state.get('standalone_ok', False)}, sequential={state.get('sequential_ok', False)}",
            flush=True,
        )

    decomposition_path = OUT / "c5_stage3b_decomposition_by_dataset.csv"
    pd.DataFrame(decomposition).to_csv(decomposition_path, index=False)

    merscope_rows: list[dict] = []
    for _, row in source_table.loc[source_table["resolution"].eq("MERSCOPE")].iterrows():
        dataset = str(row["dataset"])
        state = runtime[dataset]
        baseline = ROOT / "result" / dataset / "stage4_cytospace_baseline_highres" / "cytospace_output"
        stage3a = ROOT / "result" / dataset / "stage4_cytospace_route2_highres" / "cytospace_output"
        stand_mask = Path(state["standalone_mask"]) if state.get("standalone_ok") else None
        seq_mask = Path(state["sequential_mask"]) if state.get("sequential_ok") else None
        if stand_mask is not None:
            stand_ok, stand_route, stand_note = run_mapping(dataset, "standalone", stand_mask, row)
        else:
            stand_ok, stand_route, stand_note = False, Path(), "standalone Stage3B missing"
        if seq_mask is not None:
            full_ok, full_route, full_note = run_mapping(dataset, "sequential", seq_mask, row)
        else:
            full_ok, full_route, full_note = False, Path(), "sequential Stage3B missing"
        merscope_rows.append(
            {
                "dataset": dataset,
                "target": row["target"],
                "resolution": row["resolution"],
                "capacity": int(row["capacity"]),
                "baseline_complete": complete_stage4(baseline, requires_stage3b_mask=False),
                "baseline_path": str(baseline),
                "stage3a_only_complete": complete_stage4(stage3a, requires_stage3b_mask=False),
                "stage3a_only_path": str(stage3a),
                "stage3b_only_complete": bool(stand_ok),
                "stage3b_only_path": str(stand_route),
                "full_complete": bool(full_ok),
                "full_path": str(full_route),
                "standalone_mask_path": str(stand_mask) if stand_mask else "",
                "sequential_mask_path": str(seq_mask) if seq_mask else "",
                "notes": f"Stage3B-only: {stand_note}; Full: {full_note}",
            }
        )
        print(f"[mapping] {dataset}: Stage3B-only={stand_ok}, Full={full_ok}", flush=True)

    routes_path = OUT / "c5_merscope_four_route_summary.csv"
    pd.DataFrame(merscope_rows).to_csv(routes_path, index=False)
    (OUT / "c5_run_runtime.json").write_text(
        json.dumps(runtime, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    n_stage3b = int(pd.DataFrame(decomposition)["stage3b_success"].sum())
    n_routes = int(
        pd.DataFrame(merscope_rows)[
            ["baseline_complete", "stage3a_only_complete", "stage3b_only_complete", "full_complete"]
        ].all(axis=1).sum()
    )
    print(f"DONE stage3b={n_stage3b}/15 merscope_routes={n_routes}/5", flush=True)
    return 0 if n_stage3b == 15 and n_routes == 5 else 1


if __name__ == "__main__":
    raise SystemExit(main())
