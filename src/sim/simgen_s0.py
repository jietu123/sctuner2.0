"""
SimGen S0: 构造“ST 缺失某类型、SC 保留”的模拟数据。

输出目录: data/sim/<sample>/S0/
- sc_expression.csv   (genes x cells)
- sc_metadata.csv     (cell_id, cell_type)
- st_expression.csv   (genes x spots)  — 生成的 ST，移除 missing_type
- st_coordinates.csv  (复制原始坐标)
- sim_info.json       (记录 missing_type, cells_per_spot, seed 等)
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd
import yaml
import hashlib


def load_yaml(path: Path) -> dict:
    if not path.exists() or path.stat().st_size == 0:
        return {}
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="SimGen S0 generator (remove one type from ST)")
    p.add_argument("--sample", default="real_brca", help="sample id")
    p.add_argument("--project_root", default=None, help="override project root")
    p.add_argument("--config", default="configs/project_config.yaml", help="project config path")
    p.add_argument("--dataset_config", default=None, help="dataset config path (overrides auto-detection)")
    p.add_argument("--sim_config", default=None, help="simgen config path for S0 (optional)")
    p.add_argument("--missing_type", default=None, help="cell type to remove from ST (S0). If not set, pick the most frequent type")
    p.add_argument("--cells_per_spot", type=int, default=None, help="cells per spot for synthetic ST")
    p.add_argument("--seed", type=int, default=None, help="random seed")
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
    # 输入格式: 第一列基因名，后续列为细胞/spot
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
            raise ValueError(f"指定的 missing_type={user_type} 不存在于 sc_meta")
        return user_type
    counts = sc_meta[type_col].astype(str).value_counts()
    return counts.index[0]


def build_st_expression(
    sc_expr: pd.DataFrame,
    sc_meta: pd.DataFrame,
    coords: pd.DataFrame,
    missing_type: str,
    cells_per_spot: int,
    seed: int,
) -> Tuple[pd.DataFrame, Dict[str, int], pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    type_col = "cell_type" if "cell_type" in sc_meta.columns else sc_meta.columns[1]
    sc_meta = sc_meta.copy()
    sc_meta[type_col] = sc_meta[type_col].astype(str)
    pool_cells = sc_meta.loc[sc_meta[type_col] != missing_type, "cell_id"].tolist()
    if len(pool_cells) == 0:
        raise ValueError(f"所有细胞都是 missing_type={missing_type}，无法生成 ST")

    # map cell_id -> cell_type for truth export
    cell2type = sc_meta.set_index("cell_id")[type_col].to_dict()

    gene_names = sc_expr.index
    spot_ids = coords.index.astype(str).tolist()
    st_mat = pd.DataFrame(0.0, index=gene_names, columns=spot_ids)
    used_counts = {"spots": len(spot_ids), "cells_pool": len(pool_cells)}
    truth_rows = []
    spot_type_counts = {spot: {} for spot in spot_ids}
    qid = 0

    for spot in spot_ids:
        pick = rng.choice(pool_cells, size=cells_per_spot, replace=True)
        sub = sc_expr[pick]  # genes x picked_cells
        st_mat[spot] = sub.sum(axis=1)
        for cid in pick:
            qid += 1
            ctype = cell2type.get(cid, "unknown")
            truth_rows.append({"query_id": qid, "cell_id": cid, "true_spot_id": spot, "cell_type": ctype})
            spot_type_counts[spot][ctype] = spot_type_counts[spot].get(ctype, 0) + 1

    truth_query = pd.DataFrame(truth_rows)

    # spot-level fractions
    spot_df = pd.DataFrame(spot_type_counts).T.fillna(0)
    spot_df.index.name = "spot_id"
    counts_sum = spot_df.sum(axis=1)
    frac_df = spot_df.div(counts_sum, axis=0).fillna(0)

    return st_mat, used_counts, truth_query, frac_df


def sha1_of_file(path: Path) -> str:
    h = hashlib.sha1()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    args = parse_args()
    project_root = Path(args.project_root) if args.project_root else detect_project_root()

    # simgen config (optional override)
    sim_cfg = load_yaml(project_root / (args.sim_config or "configs/simgen/s0.yaml"))
    cells_per_spot = args.cells_per_spot or sim_cfg.get("cells_per_spot", 5)
    seed = args.seed if args.seed is not None else sim_cfg.get("seed", 42)
    missing_type_cli = args.missing_type or sim_cfg.get("missing_type")

    input_dir = project_root / "data" / "raw" / args.sample
    paths = load_paths(project_root, args.sample, args)

    sc_expr = read_expression(input_dir / paths["sc_expr"])
    sc_meta = pd.read_csv(input_dir / paths["sc_meta"], sep=None, engine="python")
    if "cell_id" not in sc_meta.columns:
        cols = sc_meta.columns.tolist()
        if len(cols) < 2:
            raise ValueError("sc_meta 至少需要两列 (cell_id, cell_type)")
        sc_meta = sc_meta.iloc[:, :2]
        sc_meta.columns = ["cell_id", "cell_type"]
    st_coords = pd.read_csv(input_dir / paths["st_meta"], index_col=0)

    missing_type = choose_missing_type(sc_meta, missing_type_cli)
    print(f"[SimGen S0] missing_type: {missing_type}")

    st_expr_s0, used_counts, truth_query, truth_frac = build_st_expression(
        sc_expr=sc_expr,
        sc_meta=sc_meta,
        coords=st_coords,
        missing_type=missing_type,
        cells_per_spot=cells_per_spot,
        seed=seed,
    )

    out_dir = project_root / "data" / "sim" / args.sample / "S0"
    out_dir.mkdir(parents=True, exist_ok=True)

    # 输出
    sc_expr_path = out_dir / "sc_expression.csv"
    sc_meta_path = out_dir / "sc_metadata.csv"
    st_expr_path = out_dir / "st_expression.csv"
    st_coords_path = out_dir / "st_coordinates.csv"
    truth_query_path = out_dir / "sim_truth_query_cell_spot.csv"
    truth_frac_path = out_dir / "sim_truth_spot_type_fraction.csv"

    sc_expr.to_csv(sc_expr_path)
    sc_meta.to_csv(sc_meta_path, index=False)
    st_expr_s0.to_csv(st_expr_path)
    st_coords.to_csv(st_coords_path)
    truth_query.to_csv(truth_query_path, index=False)
    truth_frac.to_csv(truth_frac_path)

    info = {
        "sample": args.sample,
        "missing_type": missing_type,
        "cells_per_spot": cells_per_spot,
        "seed": seed,
        "spots": used_counts["spots"],
        "cells_pool": used_counts["cells_pool"],
        "input_paths": paths,
        "sha1": {
            "sc_expression.csv": sha1_of_file(sc_expr_path),
            "sc_metadata.csv": sha1_of_file(sc_meta_path),
            "st_expression.csv": sha1_of_file(st_expr_path),
            "st_coordinates.csv": sha1_of_file(st_coords_path),
            "sim_truth_query_cell_spot.csv": sha1_of_file(truth_query_path),
            "sim_truth_spot_type_fraction.csv": sha1_of_file(truth_frac_path),
        },
        "output_dir": str(out_dir),
    }
    out_info = out_dir / "sim_info.json"
    out_info.write_text(json.dumps(info, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[SimGen S0] done. outputs in {out_dir}")


if __name__ == "__main__":
    main()

