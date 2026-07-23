#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

_HERE = Path(__file__).resolve()
_ROOT = _HERE.parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src.utils.sample_paths import resolve_sample_dir


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Create a simulation variant with a purified spatial region for one target type."
    )
    p.add_argument("--project_root", default=".")
    p.add_argument("--source_sample", required=True)
    p.add_argument("--target_sample", required=True)
    p.add_argument("--sim_group", default="real_brca")
    p.add_argument("--target_type", required=True)
    p.add_argument("--threshold", type=float, default=0.25)
    p.add_argument("--replacement_type", default="PVL")
    p.add_argument(
        "--marker_boost_top_n",
        type=int,
        default=0,
        help="Boost top target-specific genes in purified spots; <=0 disables.",
    )
    p.add_argument(
        "--marker_boost_factor",
        type=float,
        default=1.0,
        help="Multiplicative boost for target-specific genes in purified spots.",
    )
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--depth_scale", type=float, default=1.0)
    p.add_argument(
        "--hvg_nfeatures",
        type=int,
        default=None,
        help="Override qc.hvg_nfeatures in the generated dataset config.",
    )
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def read_expr_tsv(path: Path) -> tuple[np.ndarray, list[str], np.ndarray]:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    genes = df.iloc[:, 0].astype(str).to_numpy()
    ids = [str(c) for c in df.columns[1:]]
    mat = df.iloc[:, 1:].to_numpy(dtype=np.float32, copy=False)
    return genes, ids, mat


def write_expr_tsv(path: Path, genes: np.ndarray, spot_ids: list[str], mat: np.ndarray) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Gene", *spot_ids])
        for i, gene in enumerate(genes):
            writer.writerow([gene, *mat[i, :].tolist()])
            if (i + 1) % 500 == 0 or (i + 1) == len(genes):
                print(f"[WRITE] rows: {i + 1}/{len(genes)}")


def main() -> int:
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    root = Path(args.project_root).resolve()
    src = resolve_sample_dir(root, args.source_sample, sim_group=args.sim_group, must_exist=True)
    dst = root / "data" / "sim" / args.sim_group / args.target_sample
    if dst.exists():
        if not args.overwrite:
            raise FileExistsError(f"target exists: {dst}")
        shutil.rmtree(dst)
    dst.mkdir(parents=True)

    truth_path = src / "sim_truth_spot_type_fraction.csv"
    sc_meta_path = src / "brca_scRNA_celllabels.txt"
    sc_expr_path = src / "brca_scRNA_GEP.txt"
    st_expr_path = src / "brca_STdata_GEP.txt"
    st_coord_path = src / "brca_STdata_coordinates.txt"
    for path in [truth_path, sc_meta_path, sc_expr_path, st_expr_path, st_coord_path]:
        if not path.exists():
            raise FileNotFoundError(path)

    truth = pd.read_csv(truth_path)
    truth = truth.rename(columns={truth.columns[0]: "spot_id"})
    truth["spot_id"] = truth["spot_id"].astype(str)
    type_cols = [c for c in truth.columns if c != "spot_id"]
    if args.target_type not in type_cols:
        raise ValueError(f"target type not in truth table: {args.target_type}")
    if args.replacement_type not in type_cols:
        raise ValueError(f"replacement type not in truth table: {args.replacement_type}")

    weights = truth.set_index("spot_id")[type_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    target = weights[args.target_type]
    dominant = weights.idxmax(axis=1) == args.target_type
    purified = dominant | (target >= float(args.threshold))
    if int(purified.sum()) < 20:
        raise ValueError(f"too few purified spots: {int(purified.sum())}")

    # Make the target region explicit and remove diffuse target leakage elsewhere.
    weights.loc[purified, :] = 0.0
    weights.loc[purified, args.target_type] = 1.0
    leaked = weights.loc[~purified, args.target_type].copy()
    weights.loc[~purified, args.target_type] = 0.0
    weights.loc[~purified, args.replacement_type] = (
        weights.loc[~purified, args.replacement_type] + leaked
    )
    weights = weights.div(np.clip(weights.sum(axis=1), 1e-8, None), axis=0).astype(np.float32)

    print("[STEP] Build ST expression from purified truth")
    sc_meta = pd.read_csv(sc_meta_path, sep="\t")
    sc_meta = sc_meta.iloc[:, :2].copy()
    sc_meta.columns = ["cell_id", "cell_type"]
    sc_meta["cell_id"] = sc_meta["cell_id"].astype(str)
    sc_meta["cell_type"] = sc_meta["cell_type"].astype(str)
    type_names = list(type_cols)
    type_to_idx = {t: i for i, t in enumerate(type_names)}
    missing_meta_types = sorted(set(type_names).difference(sc_meta["cell_type"].unique()))
    if missing_meta_types:
        raise ValueError(f"truth types absent from SC metadata: {missing_meta_types}")

    sc_genes, sc_cells, sc_mat = read_expr_tsv(sc_expr_path)
    st_genes, st_spots, st_mat = read_expr_tsv(st_expr_path)
    meta_map = dict(zip(sc_meta["cell_id"], sc_meta["cell_type"]))
    keep_cells = np.array([cell in meta_map for cell in sc_cells], dtype=bool)
    sc_cells = [cell for cell, keep in zip(sc_cells, keep_cells) if keep]
    sc_mat = sc_mat[:, keep_cells]
    cell_type_idx = np.array([type_to_idx[meta_map[cell]] for cell in sc_cells], dtype=np.int32)

    sc_g2i = {g: i for i, g in enumerate(sc_genes)}
    st_g2i = {g: i for i, g in enumerate(st_genes)}
    common = [g for g in st_genes if g in sc_g2i]
    if len(common) < 1000:
        raise ValueError(f"too few common genes: {len(common)}")
    sc_idx = np.array([sc_g2i[g] for g in common], dtype=np.int32)
    st_idx = np.array([st_g2i[g] for g in common], dtype=np.int32)
    sc_common = sc_mat[sc_idx, :]
    st_common = st_mat[st_idx, :]
    one_hot = np.zeros((len(sc_cells), len(type_names)), dtype=np.float32)
    one_hot[np.arange(len(sc_cells)), cell_type_idx] = 1.0
    type_sum = sc_common @ one_hot
    type_profile = type_sum + 1e-3
    type_profile = type_profile / np.clip(type_profile.sum(axis=0, keepdims=True), 1e-8, None)

    w = weights.reindex(index=st_spots, columns=type_names, fill_value=0.0).to_numpy(dtype=np.float32)
    lib_sizes = np.clip(st_common.sum(axis=0, dtype=np.float64), 1.0, None) * float(args.depth_scale)
    expected = (type_profile @ w.T) * lib_sizes[np.newaxis, :]

    boosted_genes: list[str] = []
    boost_top_n = max(0, int(args.marker_boost_top_n))
    boost_factor = float(args.marker_boost_factor)
    if boost_top_n > 0 and boost_factor > 1.0:
        target_idx = type_to_idx[args.target_type]
        other_idx = [idx for idx in range(len(type_names)) if idx != target_idx]
        target_profile = type_profile[:, target_idx]
        other_profile = type_profile[:, other_idx].max(axis=1)
        specificity = np.log2((target_profile + 1e-8) / (other_profile + 1e-8))
        score = (target_profile / (other_profile + 0.02)) * target_profile
        eligible = (
            (target_profile > np.quantile(target_profile, 0.5))
            & (specificity > 1.0)
        )
        order = np.argsort(np.where(eligible, score, -np.inf))[::-1]
        keep = [idx for idx in order if np.isfinite(score[idx]) and score[idx] > 0][:boost_top_n]
        if keep:
            target_spot_mask = purified.reindex(index=st_spots).to_numpy(dtype=bool)
            boost = np.ones_like(expected, dtype=np.float32)
            target_spot_idx = np.flatnonzero(target_spot_mask)
            boost[np.ix_(np.array(keep, dtype=np.int32), target_spot_idx)] = boost_factor
            expected = expected * boost
            boosted_genes = [str(common[idx]) for idx in keep]
            print(
                f"[BOOST] target={args.target_type} genes={len(boosted_genes)} "
                f"factor={boost_factor:g} spots={int(target_spot_mask.sum())}"
            )

    sim_st = rng.poisson(np.clip(expected, 0.0, None)).astype(np.int32)

    shutil.copy2(sc_expr_path, dst / "brca_scRNA_GEP.txt")
    shutil.copy2(sc_meta_path, dst / "brca_scRNA_celllabels.txt")
    shutil.copy2(st_coord_path, dst / "brca_STdata_coordinates.txt")
    for name in ["sim_truth_query_cell_spot.csv", "sim_truth_spot_type_fraction_from_cells.csv"]:
        source = src / name
        if source.exists():
            shutil.copy2(source, dst / name)
    write_expr_tsv(dst / "brca_STdata_GEP.txt", np.array(common, dtype=object), st_spots, sim_st)

    out_truth = weights.copy()
    out_truth.insert(0, "spot_id", out_truth.index)
    out_truth.to_csv(dst / "sim_truth_spot_type_fraction.csv", index=False, encoding="utf-8")
    dominant_out = pd.DataFrame(
        {
            "spot_id": weights.index,
            "dominant_type": weights.idxmax(axis=1).to_numpy(),
            "assigned_cells": 0,
        }
    )
    dominant_out.to_csv(dst / "sim_truth_spot_dominant_type.csv", index=False, encoding="utf-8")

    source_info_path = src / "sim_info.json"
    source_info = json.loads(source_info_path.read_text(encoding="utf-8")) if source_info_path.exists() else {}
    source_missing_types: list[str] = []
    if isinstance(source_info.get("missing_types"), list):
        for value in source_info.get("missing_types") or []:
            cell_type = str(value).strip()
            if cell_type and cell_type not in source_missing_types:
                source_missing_types.append(cell_type)
    source_missing_type = str(source_info.get("missing_type") or "").strip()
    if source_missing_type and source_missing_type not in source_missing_types:
        source_missing_types.append(source_missing_type)
    info = {
        "sample": args.target_sample,
        "source_sample": args.source_sample,
        "simulation_type": "purified_target_region_from_existing_sim",
        "target_type": args.target_type,
        "purified_spots": int(purified.sum()),
        "threshold": float(args.threshold),
        "replacement_type": args.replacement_type,
        "marker_boost_top_n": int(args.marker_boost_top_n),
        "marker_boost_factor": float(args.marker_boost_factor),
        "boosted_genes": boosted_genes,
        "seed": int(args.seed),
        "source_simulation_type": source_info.get("simulation_type"),
        "missing_type": source_missing_types[0] if len(source_missing_types) == 1 else None,
        "missing_types": source_missing_types,
    }
    (dst / "sim_info.json").write_text(json.dumps(info, ensure_ascii=False, indent=2), encoding="utf-8")

    cfg_src = root / "configs" / "datasets" / f"{args.source_sample}.yaml"
    cfg_dst = root / "configs" / "datasets" / f"{args.target_sample}.yaml"
    if cfg_src.exists():
        cfg = yaml.safe_load(cfg_src.read_text(encoding="utf-8")) or {}
        plugin = ((cfg.get("stage3") or {}).get("plugin_genes_path"))
        if isinstance(plugin, str):
            cfg.setdefault("stage3", {})
            cfg["stage3"]["plugin_genes_path"] = plugin.replace(args.source_sample, args.target_sample)
        if args.hvg_nfeatures is not None:
            cfg.setdefault("qc", {})
            cfg["qc"]["hvg_nfeatures"] = int(args.hvg_nfeatures)
        cfg_dst.write_text(yaml.safe_dump(cfg, sort_keys=False, allow_unicode=True), encoding="utf-8")
        print(f"[CFG] wrote {cfg_dst}")

    print(
        f"[DONE] {dst} target={args.target_type} purified_spots={int(purified.sum())}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
