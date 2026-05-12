#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Create a simulation sample with scRNA transcriptome noise by copying an existing "
            "raw simulation directory and Stage1 processed exports."
        )
    )
    p.add_argument("--project_root", default=".", help="Project root.")
    p.add_argument("--sim_group", required=True, help="Simulation group under data/sim and processed simulation_experiments.")
    p.add_argument("--source_sample", required=True, help="Existing source sample.")
    p.add_argument("--target_sample", required=True, help="New noisy target sample.")
    p.add_argument("--noise_fraction", type=float, default=0.10, help="Fraction of sc expression entries to perturb.")
    p.add_argument("--seed", type=int, default=42, help="Random seed.")
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing target outputs.")
    return p.parse_args()


def _copytree(src: Path, dst: Path, overwrite: bool) -> None:
    if not src.exists():
        raise FileNotFoundError(f"source not found: {src}")
    if dst.exists():
        if not overwrite:
            raise FileExistsError(f"target exists: {dst} (use --overwrite)")
        shutil.rmtree(dst)
    shutil.copytree(src, dst)


def _replace_config_paths(obj, source_sample: str, target_sample: str):
    if isinstance(obj, dict):
        return {k: _replace_config_paths(v, source_sample, target_sample) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_replace_config_paths(v, source_sample, target_sample) for v in obj]
    if isinstance(obj, str):
        return obj.replace(source_sample, target_sample)
    return obj


def _write_dataset_config(root: Path, source_sample: str, target_sample: str) -> None:
    src_cfg = root / "configs" / "datasets" / f"{source_sample}.yaml"
    dst_cfg = root / "configs" / "datasets" / f"{target_sample}.yaml"
    if not src_cfg.exists():
        raise FileNotFoundError(f"source dataset config not found: {src_cfg}")
    cfg = yaml.safe_load(src_cfg.read_text(encoding="utf-8")) or {}
    cfg = _replace_config_paths(cfg, source_sample, target_sample)
    dst_cfg.parent.mkdir(parents=True, exist_ok=True)
    dst_cfg.write_text(yaml.safe_dump(cfg, sort_keys=False, allow_unicode=True), encoding="utf-8")


def _perturb_matrix(path: Path, rng: np.random.Generator, noise_fraction: float) -> dict:
    df = pd.read_csv(path, index_col=0)
    arr = df.to_numpy(dtype=np.float32, copy=True)
    mask = rng.random(arr.shape) < float(noise_fraction)

    n_rows, n_cols = arr.shape
    for j in range(n_cols):
        rows = np.flatnonzero(mask[:, j])
        if rows.size == 0:
            continue
        donor_rows = rng.integers(0, n_rows, size=rows.size)
        arr[rows, j] = arr[donor_rows, j]

    out = pd.DataFrame(arr, index=df.index, columns=df.columns)
    out.index.name = df.index.name
    out.to_csv(path)
    info = {
        "path": str(path),
        "shape": [int(n_rows), int(n_cols)],
        "perturbed_entries": int(mask.sum()),
        "total_entries": int(mask.size),
        "observed_noise_fraction": float(mask.mean()),
    }
    return info


def main() -> int:
    args = parse_args()
    if not (0.0 <= args.noise_fraction <= 1.0):
        raise ValueError("--noise_fraction must be in [0,1]")
    root = Path(args.project_root).resolve()
    rng = np.random.default_rng(args.seed)

    raw_src = root / "data" / "sim" / args.sim_group / args.source_sample
    raw_dst = root / "data" / "sim" / args.sim_group / args.target_sample
    proc_src = root / "data" / "processed" / "simulation_experiments" / args.sim_group / args.source_sample
    proc_dst = root / "data" / "processed" / "simulation_experiments" / args.sim_group / args.target_sample
    result_dst = root / "result" / args.target_sample

    _copytree(raw_src, raw_dst, args.overwrite)
    _copytree(proc_src, proc_dst, args.overwrite)
    if result_dst.exists() and args.overwrite:
        shutil.rmtree(result_dst)
    result_dst.mkdir(parents=True, exist_ok=True)
    _write_dataset_config(root, args.source_sample, args.target_sample)

    export = proc_dst / "stage1_preprocess" / "exported"
    if not export.exists():
        raise FileNotFoundError(f"processed Stage1 export not found: {export}")

    norm_info = _perturb_matrix(export / "sc_expression_normalized.csv", rng, args.noise_fraction)
    data_info = _perturb_matrix(export / "sc_expression_data.csv", rng, args.noise_fraction)
    counts_info = None
    counts_path = export / "sc_expression_counts.csv"
    if counts_path.exists():
        counts_info = _perturb_matrix(counts_path, rng, args.noise_fraction)

    for sim_info_path in [raw_dst / "sim_info.json", export / "sim_info.json"]:
        if sim_info_path.exists():
            info = json.loads(sim_info_path.read_text(encoding="utf-8-sig"))
        else:
            info = {}
        info.update(
            {
                "sample": args.target_sample,
                "source_sample": args.source_sample,
                "noise_type": "sc_reference_entry_permutation",
                "noise_fraction": float(args.noise_fraction),
                "noise_percent": float(args.noise_fraction * 100.0),
                "noise_seed": int(args.seed),
                "noise_note": (
                    "10% scRNA transcriptome noise: selected sc expression entries are replaced "
                    "by values sampled from other cells for the same gene; ST truth/expression unchanged."
                ),
            }
        )
        sim_info_path.write_text(json.dumps(info, ensure_ascii=False, indent=2), encoding="utf-8")

    summary = {
        "source_sample": args.source_sample,
        "target_sample": args.target_sample,
        "sim_group": args.sim_group,
        "noise_fraction": float(args.noise_fraction),
        "seed": int(args.seed),
        "normalized": norm_info,
        "data": data_info,
        "counts": counts_info,
        "raw_dir": str(raw_dst),
        "processed_dir": str(proc_dst),
    }
    (result_dst / "sc_noise_generation_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(f"[OK] generated noisy sample: {args.target_sample}")
    print(f"[OK] raw: {raw_dst}")
    print(f"[OK] processed: {proc_dst}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
