"""
SimGen M0: weak-missing + spatial bias for a target cell type (e.g., CD8).

Goal:
- Keep a small fraction of the missing_type in ST (keep_fraction)
- Add spatial bias so missing_type concentrates in a hotspot region
"""

from __future__ import annotations

import argparse
import json
import hashlib
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd
import yaml


def load_yaml(path: Path) -> dict:
    if not path.exists() or path.stat().st_size == 0:
        return {}
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="SimGen M0 generator (weak missing + spatial bias)")
    p.add_argument("--sample", default="real_brca", help="sample id")
    p.add_argument("--project_root", default=None, help="override project root")
    p.add_argument("--config", default="configs/project_config.yaml", help="project config path")
    p.add_argument("--dataset_config", default=None, help="dataset config path (overrides auto-detection)")
    p.add_argument("--sim_config", default=None, help="simgen config path for M0 (optional)")
    p.add_argument("--missing_type", default=None, help="cell type to downsample in ST (M0)")
    p.add_argument("--cells_per_spot", type=int, default=None, help="cells per spot for synthetic ST")
    p.add_argument("--seed", type=int, default=None, help="random seed")
    p.add_argument("--keep_fraction", type=float, default=None, help="overall missing_type fraction in ST")
    p.add_argument("--hotspot_fraction", type=float, default=None, help="fraction of spots in hotspot region")
    p.add_argument("--hotspot_cd8_fraction", type=float, default=None, help="missing_type fraction inside hotspot spots")
    p.add_argument("--scenario_tag", default=None, help="output tag under data/sim/<sample>/M0/")
    return p.parse_args()


def detect_project_root() -> Path:
    here = Path(__file__).resolve()
    try:
        return here.parents[2]
    except IndexError:
        return here.parent.parent


def load_paths(project_root: Path, sample: str, args: argparse.Namespace) -> Dict[str, str]:
    project_cfg = load_yaml(project_root / args.config)
    dataset_cfg_map = (project_cfg or {}).get("dataset_config_map", {}) or {}
    if args.dataset_config:
        dataset_cfg_path = Path(args.dataset_config)
    else:
        mapped_name = dataset_cfg_map.get(sample)
        dataset_cfg_path = project_root / "configs" / "datasets" / (mapped_name or f"{sample}.yaml")
    dataset_cfg = load_yaml(dataset_cfg_path)
    paths = (dataset_cfg or {}).get("paths", {}) or {}
    return {
        "sc_expr": paths.get("sc_expr", f"{sample}_scRNA_GEP.txt"),
        "sc_meta": paths.get("sc_meta", f"{sample}_scRNA_celllabels.txt"),
        "st_expr": paths.get("st_expr", f"{sample}_STdata_GEP.txt"),
        "st_meta": paths.get("st_meta", f"{sample}_STdata_coordinates.txt"),
    }


def read_expression(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=None, engine="python")
    genes = df.iloc[:, 0].astype(str)
    mat = df.iloc[:, 1:]
    mat.index = genes
    mat = mat.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return mat


def choose_missing_type(sc_meta: pd.DataFrame, user_type: str | None) -> str:
    type_col = "cell_type" if "cell_type" in sc_meta.columns else sc_meta.columns[1]
    if user_type:
        if user_type not in set(sc_meta[type_col].astype(str)):
            raise ValueError(f"missing_type={user_type} not found in sc_meta")
        return user_type
    counts = sc_meta[type_col].astype(str).value_counts()
    return counts.index[0]


def pick_hotspot_spots(coords: pd.DataFrame, rng: np.random.Generator, hotspot_fraction: float) -> Tuple[list[str], str]:
    spot_ids = coords.index.astype(str).tolist()
    if hotspot_fraction <= 0:
        return [], ""
    n_hot = max(1, int(len(spot_ids) * hotspot_fraction))
    coord_vals = coords.iloc[:, :2].to_numpy(dtype=float)
    center_idx = int(rng.integers(0, len(spot_ids)))
    center = coord_vals[center_idx]
    dists = np.sum((coord_vals - center) ** 2, axis=1)
    hot_idx = np.argsort(dists)[:n_hot]
    hotspot_spots = [spot_ids[i] for i in hot_idx]
    return hotspot_spots, spot_ids[center_idx]


def build_st_expression_m0(
    sc_expr: pd.DataFrame,
    sc_meta: pd.DataFrame,
    coords: pd.DataFrame,
    missing_type: str,
    cells_per_spot: int,
    seed: int,
    keep_fraction: float,
    hotspot_fraction: float,
    hotspot_cd8_fraction: float,
) -> Tuple[pd.DataFrame, Dict[str, float], pd.DataFrame, pd.DataFrame, list[str], str]:
    rng = np.random.default_rng(seed)
    type_col = "cell_type" if "cell_type" in sc_meta.columns else sc_meta.columns[1]
    sc_meta = sc_meta.copy()
    sc_meta[type_col] = sc_meta[type_col].astype(str)

    target_cells = sc_meta.loc[sc_meta[type_col] == missing_type, "cell_id"].tolist()
    other_cells = sc_meta.loc[sc_meta[type_col] != missing_type, "cell_id"].tolist()
    if len(target_cells) == 0 or len(other_cells) == 0:
        raise ValueError("missing_type pool or other pool empty, cannot build ST")

    hotspot_spots, hotspot_center = pick_hotspot_spots(coords, rng, hotspot_fraction)
    hotspot_set = set(hotspot_spots)

    if hotspot_fraction >= 1.0:
        base_fraction = keep_fraction
        hotspot_fraction_eff = 1.0
    else:
        base_fraction = (keep_fraction - hotspot_fraction * hotspot_cd8_fraction) / (1.0 - hotspot_fraction)
        hotspot_fraction_eff = hotspot_fraction
    if base_fraction < 0:
        raise ValueError("Computed base_fraction < 0. Reduce hotspot_cd8_fraction or hotspot_fraction.")
    if hotspot_cd8_fraction > 1.0 or base_fraction > 1.0:
        raise ValueError("CD8 fraction cannot exceed 1.0")

    cell2type = sc_meta.set_index("cell_id")[type_col].to_dict()

    gene_names = sc_expr.index
    spot_ids = coords.index.astype(str).tolist()
    st_mat = pd.DataFrame(0.0, index=gene_names, columns=spot_ids)
    truth_rows = []
    spot_type_counts = {spot: {} for spot in spot_ids}
    qid = 0
    cd8_count = 0
    total_count = 0

    for spot in spot_ids:
        p_cd8 = hotspot_cd8_fraction if spot in hotspot_set else base_fraction
        pick = []
        for _ in range(cells_per_spot):
            if rng.random() < p_cd8:
                cid = rng.choice(target_cells)
                cd8_count += 1
            else:
                cid = rng.choice(other_cells)
            pick.append(cid)
            total_count += 1
        sub = sc_expr[pick]
        st_mat[spot] = sub.sum(axis=1)
        for cid in pick:
            qid += 1
            ctype = cell2type.get(cid, "unknown")
            truth_rows.append({"query_id": qid, "cell_id": cid, "true_spot_id": spot, "cell_type": ctype})
            spot_type_counts[spot][ctype] = spot_type_counts[spot].get(ctype, 0) + 1

    truth_query = pd.DataFrame(truth_rows)
    spot_df = pd.DataFrame(spot_type_counts).T.fillna(0)
    spot_df.index.name = "spot_id"
    counts_sum = spot_df.sum(axis=1)
    frac_df = spot_df.div(counts_sum, axis=0).fillna(0)

    stats = {
        "keep_fraction_target": float(keep_fraction),
        "keep_fraction_actual": float(cd8_count / total_count) if total_count else 0.0,
        "hotspot_fraction": float(hotspot_fraction_eff),
        "hotspot_cd8_fraction": float(hotspot_cd8_fraction),
        "base_cd8_fraction": float(base_fraction),
    }
    return st_mat, stats, truth_query, frac_df, hotspot_spots, hotspot_center


def sha1_of_file(path: Path) -> str:
    h = hashlib.sha1()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    args = parse_args()
    project_root = Path(args.project_root) if args.project_root else detect_project_root()

    sim_cfg = load_yaml(project_root / (args.sim_config or "configs/simgen/m0.yaml"))
    cells_per_spot = args.cells_per_spot or sim_cfg.get("cells_per_spot", 5)
    seed = args.seed if args.seed is not None else sim_cfg.get("seed", 42)
    missing_type_cli = args.missing_type or sim_cfg.get("missing_type")
    keep_fraction = args.keep_fraction if args.keep_fraction is not None else sim_cfg.get("keep_fraction", 0.05)
    hotspot_fraction = args.hotspot_fraction if args.hotspot_fraction is not None else sim_cfg.get("hotspot_fraction", 0.2)
    hotspot_cd8_fraction = (
        args.hotspot_cd8_fraction if args.hotspot_cd8_fraction is not None else sim_cfg.get("hotspot_cd8_fraction", 0.2)
    )
    scenario_tag = args.scenario_tag or sim_cfg.get("scenario_tag", "cd8_weak_hotspot")

    input_dir = project_root / "data" / "raw" / args.sample
    paths = load_paths(project_root, args.sample, args)

    sc_expr = read_expression(input_dir / paths["sc_expr"])
    sc_meta = pd.read_csv(input_dir / paths["sc_meta"], sep=None, engine="python")
    if "cell_id" not in sc_meta.columns:
        cols = sc_meta.columns.tolist()
        if len(cols) < 2:
            raise ValueError("sc_meta requires at least 2 columns (cell_id, cell_type)")
        sc_meta = sc_meta.iloc[:, :2]
        sc_meta.columns = ["cell_id", "cell_type"]
    st_coords = pd.read_csv(input_dir / paths["st_meta"], index_col=0)

    missing_type = choose_missing_type(sc_meta, missing_type_cli)
    print(f"[SimGen M0] missing_type: {missing_type}")

    st_expr_m0, stats, truth_query, truth_frac, hotspot_spots, hotspot_center = build_st_expression_m0(
        sc_expr=sc_expr,
        sc_meta=sc_meta,
        coords=st_coords,
        missing_type=missing_type,
        cells_per_spot=cells_per_spot,
        seed=seed,
        keep_fraction=keep_fraction,
        hotspot_fraction=hotspot_fraction,
        hotspot_cd8_fraction=hotspot_cd8_fraction,
    )

    out_base = project_root / "data" / "sim" / args.sample / "M0" / scenario_tag
    out_dir = out_base / f"seed_{seed}"
    out_dir.mkdir(parents=True, exist_ok=True)

    sc_expr_path = out_dir / "sc_expression.csv"
    sc_meta_path = out_dir / "sc_metadata.csv"
    st_expr_path = out_dir / "st_expression.csv"
    st_coords_path = out_dir / "st_coordinates.csv"
    truth_query_path = out_dir / "sim_truth_query_cell_spot.csv"
    truth_frac_path = out_dir / "sim_truth_spot_type_fraction.csv"

    sc_expr.to_csv(sc_expr_path)
    sc_meta.to_csv(sc_meta_path, index=False)
    st_expr_m0.to_csv(st_expr_path)
    st_coords.to_csv(st_coords_path)
    truth_query.to_csv(truth_query_path, index=False)
    truth_frac.to_csv(truth_frac_path)

    info = {
        "sample": args.sample,
        "scenario": "M0",
        "scenario_tag": scenario_tag,
        "missing_type": missing_type,
        "cells_per_spot": cells_per_spot,
        "seed": seed,
        "spots": int(len(st_coords)),
        "cells_pool": int(len(sc_meta)),
        "keep_fraction": stats,
        "hotspot": {
            "fraction": float(hotspot_fraction),
            "center_spot": hotspot_center,
            "spot_count": int(len(hotspot_spots)),
        },
        "input_paths": paths,
        "sha1": {
            "sc_expression.csv": sha1_of_file(sc_expr_path),
            "sc_metadata.csv": sha1_of_file(sc_meta_path),
            "st_expression.csv": sha1_of_file(st_expr_path),
            "st_coordinates.csv": sha1_of_file(st_coords_path),
            "sim_truth_query_cell_spot.csv": sha1_of_file(truth_query_path),
            "sim_truth_spot_type_fraction.csv": sha1_of_file(truth_frac_path),
        },
        "output_dir": str(out_base),
    }
    out_info = out_dir / "sim_info.json"
    out_info.write_text(json.dumps(info, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[SimGen M0] done. outputs in {out_dir}")


if __name__ == "__main__":
    main()
