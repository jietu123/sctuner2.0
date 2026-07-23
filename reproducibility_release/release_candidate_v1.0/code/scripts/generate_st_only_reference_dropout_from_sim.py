#!/usr/bin/env python
"""Construct an ST-only scenario by removing a type from the SC reference.

ST expression, coordinates, and truth are copied unchanged. The target type is
removed only from SC metadata; Stage1 uses the recorded reference-drop metadata
to filter every SC expression export consistently.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
from pathlib import Path

import pandas as pd
import yaml

_HERE = Path(__file__).resolve()
_ROOT = _HERE.parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src.utils.sample_paths import resolve_sample_dir


UNCHANGED_FILES = (
    "brca_scRNA_GEP.txt",
    "brca_STdata_GEP.txt",
    "brca_STdata_coordinates.txt",
    "sim_truth_query_cell_spot.csv",
    "sim_truth_spot_type_fraction.csv",
    "sim_truth_spot_type_fraction_from_cells.csv",
    "sim_truth_spot_dominant_type.csv",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Remove one or more cell types from SC while leaving ST unchanged."
    )
    parser.add_argument("--project_root", default=".")
    parser.add_argument("--source_sample", required=True)
    parser.add_argument("--target_sample", required=True)
    parser.add_argument("--sim_group", default="real_brca")
    parser.add_argument(
        "--drop_cell_type",
        action="append",
        required=True,
        help="SC reference type to remove; repeat for multiple types",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def link_or_copy(source: Path, target: Path) -> str:
    try:
        os.link(source, target)
        return "hardlink"
    except OSError:
        shutil.copy2(source, target)
        return "copy"


def main() -> int:
    args = parse_args()
    root = Path(args.project_root).resolve()
    source_dir = resolve_sample_dir(
        root,
        args.source_sample,
        sim_group=args.sim_group,
        must_exist=True,
    )
    target_dir = root / "data" / "sim" / args.sim_group / args.target_sample
    drop_types = list(dict.fromkeys(str(x).strip() for x in args.drop_cell_type if str(x).strip()))
    if not drop_types:
        raise ValueError("At least one non-empty --drop_cell_type is required")

    if target_dir.exists():
        if not args.overwrite:
            raise FileExistsError(f"Target already exists: {target_dir}")
        shutil.rmtree(target_dir)
    target_dir.mkdir(parents=True)

    metadata_path = source_dir / "brca_scRNA_celllabels.txt"
    metadata = pd.read_csv(metadata_path, sep="\t")
    if metadata.shape[1] < 2:
        raise ValueError(f"Invalid SC metadata: {metadata_path}")
    metadata = metadata.iloc[:, :2].copy()
    metadata.columns = ["cell_id", "cell_type"]
    metadata["cell_id"] = metadata["cell_id"].astype(str)
    metadata["cell_type"] = metadata["cell_type"].astype(str)

    counts_before = metadata["cell_type"].value_counts().sort_index()
    missing = sorted(set(drop_types).difference(counts_before.index))
    if missing:
        raise ValueError(f"Types not present in source SC metadata: {missing}")
    keep = ~metadata["cell_type"].isin(drop_types)
    filtered = metadata.loc[keep].copy()
    if filtered.empty or filtered["cell_type"].nunique() < 2:
        raise ValueError("Reference dropout leaves fewer than two SC types")
    filtered.to_csv(
        target_dir / "brca_scRNA_celllabels.txt",
        sep="\t",
        index=False,
    )

    transfer_modes: dict[str, str] = {}
    for name in UNCHANGED_FILES:
        source = source_dir / name
        if not source.exists():
            if name == "sim_truth_spot_type_fraction_from_cells.csv":
                continue
            raise FileNotFoundError(f"Required source file missing: {source}")
        transfer_modes[name] = link_or_copy(source, target_dir / name)

    source_info_path = source_dir / "sim_info.json"
    source_info = (
        json.loads(source_info_path.read_text(encoding="utf-8"))
        if source_info_path.exists()
        else {}
    )
    source_missing_types: list[str] = []
    if isinstance(source_info.get("missing_types"), list):
        for value in source_info.get("missing_types") or []:
            cell_type = str(value).strip()
            if cell_type and cell_type not in source_missing_types:
                source_missing_types.append(cell_type)
    source_missing_type = str(source_info.get("missing_type") or "").strip()
    if source_missing_type and source_missing_type not in source_missing_types:
        source_missing_types.append(source_missing_type)
    removed_counts = {
        cell_type: int(counts_before.get(cell_type, 0)) for cell_type in drop_types
    }
    info = {
        "sample": args.target_sample,
        "source_sample": args.source_sample,
        "simulation_type": "st_only_type_missing_from_sc_reference",
        "sc_reference_drop_types": drop_types,
        "removed_sc_cells_by_type": removed_counts,
        "n_sc_cells_before": int(len(metadata)),
        "n_sc_cells_after": int(len(filtered)),
        "n_sc_types_before": int(metadata["cell_type"].nunique()),
        "n_sc_types_after": int(filtered["cell_type"].nunique()),
        "st_unchanged_from_source": True,
        "truth_unchanged_from_source": True,
        "source_simulation_type": source_info.get("simulation_type"),
        "file_transfer_modes": transfer_modes,
    }
    if source_missing_type:
        info["missing_type"] = source_missing_type
    if source_missing_types:
        info["missing_types"] = source_missing_types
    (target_dir / "sim_info.json").write_text(
        json.dumps(info, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    source_cfg = root / "configs" / "datasets" / f"{args.source_sample}.yaml"
    target_cfg = root / "configs" / "datasets" / f"{args.target_sample}.yaml"
    if not source_cfg.exists():
        raise FileNotFoundError(f"Source dataset config missing: {source_cfg}")
    config = yaml.safe_load(source_cfg.read_text(encoding="utf-8")) or {}
    plugin_path = ((config.get("stage3") or {}).get("plugin_genes_path"))
    if isinstance(plugin_path, str):
        config.setdefault("stage3", {})
        config["stage3"]["plugin_genes_path"] = plugin_path.replace(
            args.source_sample,
            args.target_sample,
        )
    config["stage3b"] = {
        "fdr": 0.05,
        "n_calibration": 0,
        "n_spatial_permutations": 200,
        "random_seed": 42,
        "sc_expr_source": "normalized",
        "expression_scale": "log1p",
        "max_genes": 0,
    }
    target_cfg.write_text(
        yaml.safe_dump(config, sort_keys=False, allow_unicode=True),
        encoding="utf-8",
    )

    print(f"[DONE] ST-only scenario: {target_dir}")
    print(f"[SC] removed: {removed_counts}")
    print(f"[SC] cells: {len(metadata)} -> {len(filtered)}")
    print("[ST] expression, coordinates, and truth unchanged")
    print(f"[CFG] {target_cfg}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
