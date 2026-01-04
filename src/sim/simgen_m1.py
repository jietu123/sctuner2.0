"""
SimGen M1: remove one cell type from sc reference while keeping ST intact.

Outputs:
data/sim/<sample>/M1/<scenario_tag>/seed_<seed>/
  - sc_expression.csv
  - sc_metadata.csv
  - st_expression.csv
  - st_coordinates.csv
  - sim_truth_query_cell_spot.csv
  - sim_truth_spot_type_fraction.csv
  - sim_info.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
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
    p = argparse.ArgumentParser(description="SimGen M1 generator (sc missing type, ST intact)")
    p.add_argument("--sample", default="real_brca", help="sample id")
    p.add_argument("--project_root", default=None, help="override project root")
    p.add_argument("--config", default="configs/project_config.yaml", help="project config path")
    p.add_argument("--dataset_config", default=None, help="dataset config path (overrides auto-detection)")
    p.add_argument("--sim_config", default=None, help="simgen config path for M1 (optional)")
    p.add_argument("--missing_type", default=None, help="cell type to remove from sc reference (M1)")
    p.add_argument("--cells_per_spot", type=int, default=None, help="cells per spot for synthetic ST")
    p.add_argument("--seed", type=int, default=None, help="random seed")
    p.add_argument("--scenario_tag", default=None, help="output tag under data/sim/<sample>/M1/")
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


def build_st_expression_full(
    sc_expr: pd.DataFrame,
    sc_meta: pd.DataFrame,
    coords: pd.DataFrame,
    cells_per_spot: int,
    seed: int,
) -> Tuple[pd.DataFrame, Dict[str, int], pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    type_col = "cell_type" if "cell_type" in sc_meta.columns else sc_meta.columns[1]
    sc_meta = sc_meta.copy()
    sc_meta[type_col] = sc_meta[type_col].astype(str)
    pool_cells = sc_meta["cell_id"].tolist()
    if len(pool_cells) == 0:
        raise ValueError("sc_meta empty, cannot build ST")

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

    return st_mat, used_counts, truth_query, frac_df


def sha1_of_file(path: Path) -> str:
    h = hashlib.sha1()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root) if args.project_root else detect_project_root()

    sim_cfg = load_yaml(project_root / (args.sim_config or "configs/simgen/m1.yaml"))
    cells_per_spot = args.cells_per_spot or sim_cfg.get("cells_per_spot", 5)
    seed = args.seed if args.seed is not None else sim_cfg.get("seed", 42)
    missing_type_cli = args.missing_type or sim_cfg.get("missing_type")
    scenario_tag = args.scenario_tag or sim_cfg.get("scenario_tag", "sc_missing_type")

    input_dir = project_root / "data" / "raw" / args.sample
    paths = load_paths(project_root, args.sample, args)

    sc_expr_full = read_expression(input_dir / paths["sc_expr"])
    sc_meta_full = pd.read_csv(input_dir / paths["sc_meta"], sep=None, engine="python")
    if "cell_id" not in sc_meta_full.columns:
        cols = sc_meta_full.columns.tolist()
        if len(cols) < 2:
            raise ValueError("sc_meta needs at least two columns (cell_id, cell_type)")
        sc_meta_full = sc_meta_full.iloc[:, :2]
        sc_meta_full.columns = ["cell_id", "cell_type"]
    st_coords = pd.read_csv(input_dir / paths["st_meta"], sep=None, engine="python", index_col=0)
    st_coords.index = st_coords.index.astype(str).str.split(r"\s|\t").str[0]

    missing_type = choose_missing_type(sc_meta_full, missing_type_cli)
    print(f"[SimGen M1] missing_type: {missing_type}")

    st_expr, used_counts, truth_query, truth_frac = build_st_expression_full(
        sc_expr=sc_expr_full,
        sc_meta=sc_meta_full,
        coords=st_coords,
        cells_per_spot=cells_per_spot,
        seed=seed,
    )

    type_col = "cell_type" if "cell_type" in sc_meta_full.columns else sc_meta_full.columns[1]
    sc_meta_full[type_col] = sc_meta_full[type_col].astype(str)
    keep_mask = sc_meta_full[type_col] != missing_type
    sc_meta = sc_meta_full.loc[keep_mask].copy()
    sc_expr = sc_expr_full[sc_meta["cell_id"].tolist()]

    sc_total = int(sc_meta_full.shape[0])
    sc_missing_count = int((~keep_mask).sum())
    sc_missing_fraction = float(sc_missing_count / sc_total) if sc_total else 0.0

    st_missing_count = int((truth_query["cell_type"] == missing_type).sum())
    st_total = int(truth_query.shape[0])
    st_missing_fraction = float(st_missing_count / st_total) if st_total else 0.0

    out_base = project_root / "data" / "sim" / args.sample / "M1" / scenario_tag / f"seed_{seed}"
    out_base.mkdir(parents=True, exist_ok=True)

    sc_expr_path = out_base / "sc_expression.csv"
    sc_meta_path = out_base / "sc_metadata.csv"
    st_expr_path = out_base / "st_expression.csv"
    st_coords_path = out_base / "st_coordinates.csv"
    truth_query_path = out_base / "sim_truth_query_cell_spot.csv"
    truth_frac_path = out_base / "sim_truth_spot_type_fraction.csv"

    sc_expr.to_csv(sc_expr_path)
    sc_meta.to_csv(sc_meta_path, index=False)
    st_expr.to_csv(st_expr_path)
    st_coords.to_csv(st_coords_path)
    truth_query.to_csv(truth_query_path, index=False)
    truth_frac.to_csv(truth_frac_path)

    info = {
        "sample": args.sample,
        "scenario": "M1",
        "scenario_tag": scenario_tag,
        "missing_type": missing_type,
        "cells_per_spot": cells_per_spot,
        "seed": seed,
        "spots": used_counts["spots"],
        "cells_pool": used_counts["cells_pool"],
        "stats": {
            "sc_missing_count": sc_missing_count,
            "sc_missing_fraction": sc_missing_fraction,
            "st_missing_count": st_missing_count,
            "st_missing_fraction": st_missing_fraction,
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
    out_info = out_base / "sim_info.json"
    out_info.write_text(json.dumps(info, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[SimGen M1] done. outputs in {out_base}")


if __name__ == "__main__":
    main()
