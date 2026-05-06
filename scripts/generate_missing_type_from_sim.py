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

from src.utils.sample_paths import resolve_sample_dir, sample_dir_candidates


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate a new simulation sample by dropping one cell type from ST truth."
    )
    p.add_argument("--project_root", default=".", help="Project root path.")
    p.add_argument("--source_sample", required=True, help="Source sample id.")
    p.add_argument("--target_sample", required=True, help="Target sample id.")
    p.add_argument(
        "--sim_group",
        default="real_brca",
        help="simulation group under data/sim/<sim_group>/ (default: real_brca)",
    )
    p.add_argument("--drop_cell_type", required=True, help="Cell type to remove from ST truth.")
    p.add_argument("--drop_fraction", type=float, default=1.0, help="Fraction to remove from ST truth.")
    p.add_argument(
        "--replacement_cell_type",
        default=None,
        help=(
            "Optional replacement type for dropped cells. "
            "If set, dropped cells are replaced at the same spots (no spatial holes)."
        ),
    )
    p.add_argument("--seed", type=int, default=42, help="Random seed.")
    p.add_argument("--depth_scale", type=float, default=1.0, help="Scale for simulated ST library size.")
    p.add_argument("--overwrite", action="store_true", help="Overwrite target dir if exists.")
    return p.parse_args()


def read_expr_tsv(path: Path) -> tuple[np.ndarray, list[str], np.ndarray]:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    genes = df.iloc[:, 0].astype(str).to_numpy()
    sample_ids = [str(c) for c in df.columns[1:]]
    mat = df.iloc[:, 1:].to_numpy(dtype=np.float32, copy=False)
    return genes, sample_ids, mat


def write_expr_tsv(path: Path, genes: np.ndarray, spot_ids: list[str], mat: np.ndarray) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Gene", *spot_ids])
        n = mat.shape[0]
        for i, g in enumerate(genes):
            w.writerow([g, *mat[i, :].tolist()])
            if (i + 1) % 500 == 0 or (i + 1) == n:
                print(f"[WRITE] rows: {i + 1}/{n}")


def main() -> int:
    args = parse_args()
    if not (0.0 <= args.drop_fraction <= 1.0):
        raise ValueError("--drop_fraction must be in [0,1]")
    rng = np.random.default_rng(args.seed)

    root = Path(args.project_root).resolve()
    src_dir = resolve_sample_dir(root, args.source_sample, sim_group=args.sim_group, must_exist=True)
    dst_dir = root / "data" / "sim" / args.sim_group / args.target_sample

    if not src_dir.exists():
        raise FileNotFoundError(f"source sample dir not found: {src_dir}")
    if dst_dir.exists():
        if not args.overwrite:
            raise FileExistsError(f"target sample exists: {dst_dir} (use --overwrite)")
        shutil.rmtree(dst_dir)
    dst_dir.mkdir(parents=True, exist_ok=True)

    sc_expr_src = src_dir / "brca_scRNA_GEP.txt"
    sc_meta_src = src_dir / "brca_scRNA_celllabels.txt"
    st_expr_src = src_dir / "brca_STdata_GEP.txt"
    st_meta_src = src_dir / "brca_STdata_coordinates.txt"
    truth_query_src = src_dir / "sim_truth_query_cell_spot.csv"
    truth_spot_src = src_dir / "sim_truth_spot_type_fraction.csv"
    for p in [sc_expr_src, sc_meta_src, st_expr_src, st_meta_src, truth_query_src, truth_spot_src]:
        if not p.exists():
            raise FileNotFoundError(f"missing file: {p}")

    def collect_missing_types_from_chain(start_sample: str, start_dir: Path | None = None, max_depth: int = 12) -> list[str]:
        out: list[str] = []
        cur = str(start_sample).strip()
        seen: set[str] = set()
        for _ in range(max_depth):
            if not cur or cur in seen:
                break
            seen.add(cur)
            info = None
            bases = []
            if start_dir is not None:
                bases.append(start_dir)
                start_dir = None
            bases.extend(sample_dir_candidates(root, cur, sim_group=args.sim_group))
            for base in bases:
                p = base / "sim_info.json"
                if not p.exists():
                    continue
                try:
                    info = json.loads(p.read_text(encoding="utf-8"))
                    break
                except Exception:
                    info = None
            if info is None:
                break
            mts = info.get("missing_types")
            if isinstance(mts, list):
                for x in mts:
                    xx = str(x).strip()
                    if xx and xx not in out:
                        out.append(xx)
            mt = info.get("missing_type")
            mt = str(mt).strip() if mt is not None else ""
            if mt and mt not in out:
                out.append(mt)
            parent = info.get("source_sample")
            cur = str(parent).strip() if parent is not None else ""
        return out

    inherited_missing_types: list[str] = collect_missing_types_from_chain(args.source_sample, src_dir)

    print("[STEP] load truth query")
    truth_q = pd.read_csv(truth_query_src)
    required = {"query_id", "cell_id", "true_spot_id", "cell_type"}
    miss = required.difference(truth_q.columns)
    if miss:
        raise ValueError(f"truth query missing columns: {sorted(miss)}")

    total_drop_type = int((truth_q["cell_type"] == args.drop_cell_type).sum())
    if total_drop_type == 0:
        raise ValueError(f"drop type not found in source truth query: {args.drop_cell_type}")

    print("[STEP] load SC metadata")
    sc_meta = pd.read_csv(sc_meta_src, sep="\t")
    sc_meta = sc_meta.rename(columns={sc_meta.columns[0]: "cell_id", sc_meta.columns[1]: "cell_type"})
    sc_meta = sc_meta[["cell_id", "cell_type"]].drop_duplicates("cell_id")
    sc_meta["cell_id"] = sc_meta["cell_id"].astype(str)
    sc_meta["cell_type"] = sc_meta["cell_type"].astype(str)

    replacement_type = None
    donor_ids = None
    if args.replacement_cell_type is not None and str(args.replacement_cell_type).strip():
        replacement_type = str(args.replacement_cell_type).strip()
        if replacement_type == args.drop_cell_type:
            raise ValueError("--replacement_cell_type cannot be the same as --drop_cell_type.")
        donor_ids = sc_meta.loc[sc_meta["cell_type"] == replacement_type, "cell_id"].to_numpy(dtype=object)
        if donor_ids.size == 0:
            raise ValueError(f"replacement type not found in SC metadata: {replacement_type}")

    drop_mask = (truth_q["cell_type"] == args.drop_cell_type) & (rng.random(len(truth_q)) < args.drop_fraction)
    truth_q_new = truth_q.loc[~drop_mask].copy()
    dropped_n = int(drop_mask.sum())
    print(
        f"[DROP] type={args.drop_cell_type} total={total_drop_type} dropped={dropped_n} "
        f"fraction={dropped_n / max(total_drop_type,1):.4f}"
    )

    replaced_n = 0
    n_after_drop_before_replacement = int(len(truth_q_new))
    if replacement_type is not None and dropped_n > 0:
        repl = truth_q.loc[drop_mask].copy()
        repl["cell_type"] = replacement_type
        # Keep spatial location (true_spot_id) unchanged; swap to donor cell IDs from replacement type.
        repl["cell_id"] = rng.choice(donor_ids, size=len(repl), replace=True).astype(str)
        truth_q_new = pd.concat([truth_q_new, repl], ignore_index=True)
        replaced_n = int(len(repl))
        print(f"[REPLACE] dropped cells replaced with {replacement_type}: {replaced_n}")

    if truth_q_new.empty:
        raise ValueError("all cells dropped; invalid scenario")

    all_types = sorted(sc_meta["cell_type"].unique().tolist())
    type_to_idx = {t: i for i, t in enumerate(all_types)}

    print("[STEP] recompute spot fractions from filtered truth query")
    coords = pd.read_csv(st_meta_src, sep="\t")
    coords = coords.rename(columns={coords.columns[0]: "spot_id", coords.columns[1]: "row", coords.columns[2]: "col"})
    spot_ids = coords["spot_id"].astype(str).to_list()

    ctab = pd.crosstab(
        pd.Categorical(truth_q_new["true_spot_id"].astype(str), categories=spot_ids),
        pd.Categorical(truth_q_new["cell_type"].astype(str), categories=all_types),
    )
    ctab.index = ctab.index.astype(str)
    ctab.columns = ctab.columns.astype(str)
    ctab = ctab.reindex(index=spot_ids, fill_value=0).reindex(columns=all_types, fill_value=0)
    frac_from_query = ctab.div(np.clip(ctab.sum(axis=1), 1, None), axis=0).astype(np.float32)

    # When replacement type is provided, keep spot coverage from source truth and move dropped-type mass
    # to replacement type on the same spots (more biologically plausible than leaving holes).
    if replacement_type is not None:
        truth_spot_src_df = pd.read_csv(truth_spot_src)
        if truth_spot_src_df.shape[1] < 2:
            raise ValueError(f"invalid source truth spot fraction file: {truth_spot_src}")
        truth_spot_src_df = truth_spot_src_df.rename(columns={truth_spot_src_df.columns[0]: "spot_id"})
        truth_spot_src_df["spot_id"] = truth_spot_src_df["spot_id"].astype(str)
        for t in all_types:
            if t not in truth_spot_src_df.columns:
                truth_spot_src_df[t] = 0.0
        frac = truth_spot_src_df.set_index("spot_id").reindex(index=spot_ids, columns=all_types, fill_value=0.0)
        frac = frac.apply(pd.to_numeric, errors="coerce").fillna(0.0)
        moved = frac[args.drop_cell_type] * float(args.drop_fraction)
        frac[args.drop_cell_type] = frac[args.drop_cell_type] * (1.0 - float(args.drop_fraction))
        frac[replacement_type] = frac[replacement_type] + moved
        frac = frac.astype(np.float32)
        print("[SPOT] mode=source_truth_transfer (keep coverage, move missing mass to replacement type)")
    else:
        frac = frac_from_query
        print("[SPOT] mode=query_recount")

    print("[STEP] build ST expression from SC profiles + new spot fractions")
    sc_genes, sc_cells, sc_mat = read_expr_tsv(sc_expr_src)
    st_genes, st_spots, st_mat_src = read_expr_tsv(st_expr_src)
    sc_cells_arr = np.array(sc_cells, dtype=object)

    # Align cells present in sc metadata
    meta_map = dict(zip(sc_meta["cell_id"], sc_meta["cell_type"]))
    keep_cell = np.array([c in meta_map for c in sc_cells_arr], dtype=bool)
    sc_cells_arr = sc_cells_arr[keep_cell]
    sc_mat = sc_mat[:, keep_cell]
    sc_type_idx = np.array([type_to_idx[meta_map[c]] for c in sc_cells_arr], dtype=np.int32)

    # Align genes between SC and source ST
    sc_g2i = {g: i for i, g in enumerate(sc_genes)}
    common_genes = [g for g in st_genes if g in sc_g2i]
    if len(common_genes) < 1000:
        raise ValueError(f"too few common genes: {len(common_genes)}")
    st_g2i = {g: i for i, g in enumerate(st_genes)}
    sc_idx = np.array([sc_g2i[g] for g in common_genes], dtype=np.int32)
    st_idx = np.array([st_g2i[g] for g in common_genes], dtype=np.int32)
    sc_common = sc_mat[sc_idx, :]
    st_common = st_mat_src[st_idx, :]
    del sc_mat
    del st_mat_src

    n_cells = len(sc_cells_arr)
    n_types = len(all_types)
    onehot = np.zeros((n_cells, n_types), dtype=np.float32)
    onehot[np.arange(n_cells), sc_type_idx] = 1.0

    # Type profile: gene x type
    type_sum = sc_common @ onehot
    type_profile = type_sum + 1e-3
    type_profile = type_profile / np.clip(type_profile.sum(axis=0, keepdims=True), 1e-8, None)

    weights = frac.reindex(index=st_spots, columns=all_types, fill_value=0.0).to_numpy(dtype=np.float32)
    lib_sizes = st_common.sum(axis=0, dtype=np.float64)
    lib_sizes = np.clip(lib_sizes, 1.0, None) * float(args.depth_scale)
    expected = type_profile @ weights.T
    expected = expected * lib_sizes[np.newaxis, :]
    sim_st = rng.poisson(np.clip(expected, 0.0, None)).astype(np.int32)

    print("[STEP] write new sample files")
    shutil.copy2(sc_expr_src, dst_dir / "brca_scRNA_GEP.txt")
    shutil.copy2(sc_meta_src, dst_dir / "brca_scRNA_celllabels.txt")
    shutil.copy2(st_meta_src, dst_dir / "brca_STdata_coordinates.txt")
    write_expr_tsv(dst_dir / "brca_STdata_GEP.txt", np.array(common_genes, dtype=object), st_spots, sim_st)

    # truth outputs
    truth_q_new.to_csv(dst_dir / "sim_truth_query_cell_spot.csv", index=False, encoding="utf-8")

    truth_spot_out = frac.copy()
    truth_spot_out.insert(0, "spot_id", truth_spot_out.index)
    truth_spot_out.to_csv(dst_dir / "sim_truth_spot_type_fraction.csv", index=False, encoding="utf-8")

    dominant = pd.DataFrame(
        {
            "spot_id": frac.index.to_numpy(),
            "dominant_type": np.array(all_types, dtype=object)[np.argmax(frac.to_numpy(), axis=1)],
            "assigned_cells": ctab.sum(axis=1).to_numpy(dtype=np.int32),
        }
    )
    dominant.to_csv(dst_dir / "sim_truth_spot_dominant_type.csv", index=False, encoding="utf-8")

    query_id_new = f"{args.target_sample}_q0"
    truth_q_new["query_id"] = query_id_new
    truth_q_new.to_csv(dst_dir / "sim_truth_query_cell_spot.csv", index=False, encoding="utf-8")

    sim_info = {
        "sample": args.target_sample,
        "source_sample": args.source_sample,
        "simulation_type": (
            "missing_type_from_existing_sim_with_replacement"
            if replacement_type is not None
            else "missing_type_from_existing_sim"
        ),
        "seed": int(args.seed),
        "query_id": query_id_new,
        "missing_type": args.drop_cell_type,
        "missing_types": list(
            dict.fromkeys(inherited_missing_types + ([args.drop_cell_type] if dropped_n > 0 else []))
        ),
        "drop_fraction": float(args.drop_fraction),
        "n_cells_before_drop": int(len(truth_q)),
        "n_cells_after_drop_before_replacement": int(n_after_drop_before_replacement),
        "n_cells_after_drop": int(len(truth_q_new)),
        "replacement_type": replacement_type,
        "replacement_cells_added": int(replaced_n),
        "dropped_cells": int(dropped_n),
        "cells_per_spot": float(ctab.sum(axis=1).mean()),
        "n_spots": int(len(st_spots)),
        "n_types": int(len(all_types)),
        "n_genes": int(len(common_genes)),
    }
    (dst_dir / "sim_info.json").write_text(json.dumps(sim_info, indent=2, ensure_ascii=False), encoding="utf-8")

    # Optional: create dataset config for the target sample by copying source config.
    cfg_src = root / "configs" / "datasets" / f"{args.source_sample}.yaml"
    cfg_dst = root / "configs" / "datasets" / f"{args.target_sample}.yaml"
    if cfg_src.exists():
        cfg = yaml.safe_load(cfg_src.read_text(encoding="utf-8")) or {}
        try:
            plugin_genes = (
                (cfg.get("stage3") or {}).get("plugin_genes_path")
                if isinstance(cfg, dict)
                else None
            )
            if isinstance(plugin_genes, str) and args.source_sample in plugin_genes:
                cfg.setdefault("stage3", {})
                cfg["stage3"]["plugin_genes_path"] = plugin_genes.replace(args.source_sample, args.target_sample)
        except Exception:
            pass
        cfg_dst.parent.mkdir(parents=True, exist_ok=True)
        cfg_dst.write_text(yaml.safe_dump(cfg, sort_keys=False, allow_unicode=True), encoding="utf-8")
        print(f"[CFG] wrote dataset config: {cfg_dst}")
    else:
        print(f"[CFG] source dataset config not found, skip: {cfg_src}")

    print(f"[DONE] generated: {dst_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
