#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import yaml
from scipy import sparse

_HERE = Path(__file__).resolve()
_ROOT = _HERE.parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src.stages.storage import stage1_export_dir


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Prepare a cell2location mouse brain reference-drop case for Stage3B."
    )
    p.add_argument("--project_root", default=".")
    p.add_argument("--section", default="ST8059048")
    p.add_argument("--target_type", default="Oligo_2")
    p.add_argument(
        "--drop_types",
        default=None,
        help="comma-separated annotation_1 types to remove; defaults to target_type",
    )
    p.add_argument("--target_label", default=None, help="display label for a multi-type target")
    p.add_argument("--sample", default=None)
    p.add_argument("--max_cells_per_type", type=int, default=80)
    p.add_argument("--n_genes", type=int, default=1500)
    p.add_argument("--target_marker_genes", type=int, default=300)
    p.add_argument("--random_seed", type=int, default=42)
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def read_10x_h5(path: Path) -> tuple[sparse.csr_matrix, np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as handle:
        matrix = handle["matrix"]
        data = matrix["data"][:]
        indices = matrix["indices"][:]
        indptr = matrix["indptr"][:]
        shape = tuple(matrix["shape"][:])
        counts = sparse.csc_matrix((data, indices, indptr), shape=shape).T.tocsr()
        barcodes = np.array(
            [x.decode() if isinstance(x, bytes) else str(x) for x in matrix["barcodes"][:]]
        )
        genes = np.array(
            [
                x.decode() if isinstance(x, bytes) else str(x)
                for x in matrix["features"]["name"][:]
            ]
        )
    return counts, barcodes, genes


def log_normalize_rows(values: sparse.spmatrix | np.ndarray, scale_factor: float = 10000.0) -> np.ndarray:
    if sparse.issparse(values):
        arr = values.astype(np.float32).toarray()
    else:
        arr = np.asarray(values, dtype=np.float32)
    totals = arr.sum(axis=1, keepdims=True)
    totals[totals <= 0] = 1.0
    return np.log1p(arr / totals * scale_factor).astype(np.float32)


def select_cells(labels: pd.Series, max_cells_per_type: int, rng: np.random.Generator) -> list[str]:
    selected: list[str] = []
    for cell_type, ids in labels.groupby(labels).groups.items():
        cell_type = str(cell_type)
        if cell_type.startswith(("LowQ", "Unk")):
            continue
        ids = np.array(list(ids), dtype=object)
        if len(ids) > max_cells_per_type:
            ids = rng.choice(ids, size=max_cells_per_type, replace=False)
        selected.extend(str(x) for x in ids)
    return selected


def unique_gene_positions(genes: np.ndarray) -> dict[str, int]:
    out: dict[str, int] = {}
    for i, gene in enumerate(genes.astype(str)):
        if gene not in out:
            out[gene] = i
    return out


def target_markers_from_profiles(raw_var: pd.DataFrame, target_type: str, st_gene_pos: dict[str, int]) -> list[str]:
    column = f"means_cov_effect_annotation_1_{target_type}"
    if column not in raw_var.columns:
        return []
    genes = (
        raw_var["SYMBOL"].astype(str).to_numpy()
        if "SYMBOL" in raw_var.columns
        else raw_var.index.astype(str).to_numpy()
    )
    profile = raw_var[column].to_numpy(dtype=float)
    order = np.argsort(profile)[::-1]
    markers: list[str] = []
    for idx in order:
        gene = str(genes[idx])
        if gene in st_gene_pos and gene not in markers:
            markers.append(gene)
    return markers


def target_markers_from_counts(
    raw_counts: sparse.spmatrix,
    raw_var: pd.DataFrame,
    labels: pd.Series,
    target_types: list[str],
    st_gene_pos: dict[str, int],
) -> list[str]:
    genes = (
        raw_var["SYMBOL"].astype(str).to_numpy()
        if "SYMBOL" in raw_var.columns
        else raw_var.index.astype(str).to_numpy()
    )
    target_mask = labels.astype(str).isin(target_types).to_numpy()
    if int(target_mask.sum()) == 0:
        raise ValueError(f"Target types are absent from sc.h5ad: {target_types}")
    common_idx = np.array(
        [i for i, gene in enumerate(genes) if str(gene) in st_gene_pos],
        dtype=int,
    )
    common_genes = np.array([str(genes[i]) for i in common_idx])
    x = raw_counts[:, common_idx].tocsr()
    xt = x[target_mask, :]
    xo = x[~target_mask, :]
    target_mean = np.asarray(xt.mean(axis=0)).ravel()
    other_mean = np.asarray(xo.mean(axis=0)).ravel()
    target_detect = np.asarray((xt > 0).mean(axis=0)).ravel()
    other_detect = np.asarray((xo > 0).mean(axis=0)).ravel()
    score = np.log2((target_mean + 0.05) / (other_mean + 0.05))
    score *= np.maximum(target_detect - 0.5 * other_detect, 0.0)
    order = np.argsort(score)[::-1]
    markers: list[str] = []
    for idx in order:
        if score[idx] <= 0:
            break
        gene = str(common_genes[idx])
        if gene not in markers:
            markers.append(gene)
    return markers


def build_gene_set(
    st_counts: sparse.csr_matrix,
    st_genes: np.ndarray,
    sc_genes: np.ndarray,
    target_markers: list[str],
    n_genes: int,
    n_target_markers: int,
) -> list[str]:
    st_pos = unique_gene_positions(st_genes)
    sc_pos = unique_gene_positions(sc_genes)
    common = [gene for gene in st_pos if gene in sc_pos]
    if not common:
        raise ValueError("No common genes between SC and ST")
    common_idx = np.array([st_pos[g] for g in common], dtype=int)
    x = st_counts[:, common_idx].astype(np.float32).tocsr()
    lib = np.asarray(x.sum(axis=1)).ravel()
    lib[lib <= 0] = 1.0
    norm = x.multiply(10000.0 / lib[:, None]).log1p().tocsr()
    mean = np.asarray(norm.mean(axis=0)).ravel()
    mean_sq = np.asarray(norm.power(2).mean(axis=0)).ravel()
    variance = mean_sq - mean * mean
    variable = [common[i] for i in np.argsort(variance)[::-1]]

    selected: list[str] = []
    for gene in target_markers[:n_target_markers]:
        if gene in sc_pos and gene in st_pos and gene not in selected:
            selected.append(gene)
    for gene in variable:
        if gene not in selected:
            selected.append(gene)
        if len(selected) >= n_genes:
            break
    return selected


def write_config(project_root: Path, sample: str) -> None:
    cfg = {
        "paths": {
            "sc_expr": "unused_stage1_direct.csv",
            "sc_meta": "unused_stage1_direct.csv",
            "st_expr": "unused_stage1_direct.csv",
            "st_meta": "unused_stage1_direct.csv",
            "svg_marker_whitelist": None,
        },
        "qc": {
            "sc_min_genes": 5,
            "sc_max_genes": 10000,
            "sc_max_mt": 100,
            "st_min_genes": 5,
            "st_max_genes": float("inf"),
            "st_max_mt": 100,
            "hvg_nfeatures": 0,
            "mt_pattern": "^(mt-|MT-)",
        },
        "gene_filter": {"min_cells_sc": 0, "min_cells_st": 0},
        "stage3b": {
            "fdr": 0.05,
            "n_calibration": 0,
            "n_spatial_permutations": 200,
            "random_seed": 42,
            "sc_expr_source": "counts",
            "sc_profile_source": "counts",
            "expression_scale": "linear",
            "sc_profile_scale": "linear",
            "max_genes": 0,
        },
        "storage": {"group": "cell2location_mouse_brain"},
    }
    path = project_root / "configs" / "datasets" / f"{sample}.yaml"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(cfg, sort_keys=False, allow_unicode=True), encoding="utf-8")


def main() -> int:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    rng = np.random.default_rng(args.random_seed)
    drop_types = [
        item.strip()
        for item in (args.drop_types if args.drop_types is not None else args.target_type).split(",")
        if item.strip()
    ]
    target_label = args.target_label or ("+".join(drop_types) if len(drop_types) > 1 else args.target_type)
    sample = args.sample or f"cell2loc_mouse_brain_{args.section.lower()}_sc_missing_{target_label.lower()}"

    source = project_root / "data" / "raw" / "cell2location_mouse_brain"
    sc_path = source / "sc.h5ad"
    section_dir = source / "mouse_brain_visium_wo_cloupe_data" / "rawdata" / args.section
    st_h5 = section_dir / "filtered_feature_bc_matrix.h5"
    pos_csv = section_dir / "spatial" / "tissue_positions_list.csv"

    sc_data = ad.read_h5ad(sc_path)
    if sc_data.raw is None:
        raise ValueError("sc.h5ad lacks raw counts")
    labels = sc_data.obs["annotation_1"].astype(str)
    selected_cells = select_cells(labels, args.max_cells_per_type, rng)
    selected_cells = [cell for cell in selected_cells if labels.loc[cell] not in set(drop_types)]
    if not selected_cells:
        raise ValueError("SC reference dropout removed every selected cell")

    st_counts, st_barcodes, st_genes = read_10x_h5(st_h5)
    positions = pd.read_csv(
        pos_csv,
        header=None,
        names=["spot_id", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"],
    )
    positions["spot_id"] = positions["spot_id"].astype(str)
    barcode_pos = {barcode: i for i, barcode in enumerate(st_barcodes.astype(str))}
    tissue_spots = [
        spot for spot in positions.loc[positions["in_tissue"] == 1, "spot_id"].tolist() if spot in barcode_pos
    ]
    st_rows = np.array([barcode_pos[spot] for spot in tissue_spots], dtype=int)
    st_counts = st_counts[st_rows, :].tocsr()

    sc_genes = (
        sc_data.raw.var["SYMBOL"].astype(str).to_numpy()
        if "SYMBOL" in sc_data.raw.var.columns
        else sc_data.raw.var.index.astype(str).to_numpy()
    )
    st_gene_pos = unique_gene_positions(st_genes)
    sc_gene_pos = unique_gene_positions(sc_genes)
    raw_counts = sc_data.raw.X
    if not sparse.issparse(raw_counts):
        raw_counts = sparse.csr_matrix(raw_counts)
    target_markers = target_markers_from_counts(
        raw_counts,
        sc_data.raw.var,
        labels,
        drop_types,
        st_gene_pos,
    )
    fallback_markers: list[str] = []
    for drop_type in drop_types:
        fallback_markers.extend(target_markers_from_profiles(sc_data.raw.var, drop_type, st_gene_pos))
    for gene in fallback_markers:
        if gene not in target_markers:
            target_markers.append(gene)
    genes = build_gene_set(
        st_counts,
        st_genes,
        sc_genes,
        target_markers,
        args.n_genes,
        args.target_marker_genes,
    )
    sc_gene_idx = np.array([sc_gene_pos[g] for g in genes], dtype=int)
    st_gene_idx = np.array([st_gene_pos[g] for g in genes], dtype=int)

    sc_cell_idx = np.array([sc_data.obs_names.get_loc(cell) for cell in selected_cells], dtype=int)
    sc_counts = raw_counts[sc_cell_idx, :][:, sc_gene_idx]
    if sparse.issparse(sc_counts):
        sc_counts = sc_counts.toarray()
    sc_counts = np.asarray(sc_counts, dtype=np.float32)
    st_counts_sel = st_counts[:, st_gene_idx].astype(np.float32)

    sc_norm = log_normalize_rows(sc_counts)
    st_norm = log_normalize_rows(st_counts_sel)
    st_counts_dense = st_counts_sel.toarray().astype(np.float32)

    export_dir = stage1_export_dir(project_root, sample, {"storage": {"group": "cell2location_mouse_brain"}})
    if export_dir.exists() and not args.overwrite:
        raise FileExistsError(f"Stage1 export already exists: {export_dir}")
    export_dir.mkdir(parents=True, exist_ok=True)

    sc_meta = pd.DataFrame(
        {
            "cell_id": selected_cells,
            "cell_type": labels.loc[selected_cells].astype(str).to_numpy(),
        }
    )
    st_meta = positions.set_index("spot_id").loc[tissue_spots, ["array_row", "array_col", "pxl_row", "pxl_col"]]
    st_meta = st_meta.rename(columns={"array_row": "row", "array_col": "col"})

    pd.DataFrame(sc_counts, index=selected_cells, columns=genes).to_csv(
        export_dir / "sc_expression_counts.csv", index_label="cell_id"
    )
    pd.DataFrame(sc_norm, index=selected_cells, columns=genes).to_csv(
        export_dir / "sc_expression_normalized.csv", index_label="cell_id"
    )
    pd.DataFrame(st_counts_dense, index=tissue_spots, columns=genes).to_csv(
        export_dir / "st_expression_counts.csv", index_label="spot_id"
    )
    pd.DataFrame(st_norm, index=tissue_spots, columns=genes).to_csv(
        export_dir / "st_expression_normalized.csv", index_label="spot_id"
    )
    sc_meta.to_csv(export_dir / "sc_metadata.csv", index=False)
    st_meta.to_csv(export_dir / "st_coordinates.csv", index_label="spot_id")

    stage1_dir = export_dir.parent
    (stage1_dir / "common_genes.txt").write_text("\n".join(genes) + "\n", encoding="utf-8")
    (stage1_dir / "hvg_genes.txt").write_text("\n".join(genes) + "\n", encoding="utf-8")
    marker_path = stage1_dir / "target_marker_genes.txt"
    marker_path.write_text("\n".join([g for g in target_markers if g in genes]) + "\n", encoding="utf-8")
    info = {
        "source": "cell2location mouse brain tutorial data",
        "section": args.section,
        "target_type": target_label,
        "drop_types": drop_types,
        "sample": sample,
        "sc_cells": len(selected_cells),
        "st_spots": len(tissue_spots),
        "genes": len(genes),
        "max_cells_per_type": args.max_cells_per_type,
        "target_reference_cells_removed": int(labels.isin(drop_types).sum()),
    }
    (stage1_dir / "stage3b_case_info.json").write_text(json.dumps(info, indent=2), encoding="utf-8")
    write_config(project_root, sample)
    print(json.dumps(info, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
