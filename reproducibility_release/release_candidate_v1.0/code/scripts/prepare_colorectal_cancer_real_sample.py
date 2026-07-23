#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import yaml
from scipy import sparse


CLUSTER_TO_CELL_TYPE = {
    1: "Epithelial cell",
    2: "Epithelial cell",
    3: "Epithelial cell",
    4: "Epithelial cell",
    5: "Epithelial cell",
    6: "Fibroblast",
    7: "Plasma B cell",
    8: "Myeloid cell",
    9: "Myeloid cell",
    10: "Epithelial cell",
    11: "T cell",
    12: "Endothelial cell",
}

CLUSTER_EVIDENCE = {
    1: "EPCAM/AGR2/KRT8/MUC2 epithelial marker score dominates",
    2: "EPCAM/AGR2/KRT8/MUC2/REG4 epithelial marker score dominates",
    3: "AGR2/EPCAM/KRT8/MUC2/TFF3/KRT19 epithelial marker score dominates",
    4: "EPCAM/AGR2/KRT8/MUC2 epithelial marker score dominates",
    5: "AGR2/EPCAM/KRT8/MUC2/TFF3/KRT19 epithelial marker score dominates",
    6: "COL1A1/COL3A1/COL1A2/COL6A2/LUM/DCN fibroblast marker score dominates",
    7: "IGKC/IGHG1/XBP1/IGHG3/MZB1/JCHAIN plasma-cell marker score dominates",
    8: "CD68/LYZ/TYROBP/C1QA/C1QC/C1QB myeloid marker score dominates",
    9: "LYZ/CD68/C1QA/TYROBP/C1QB/C1QC myeloid marker score dominates",
    10: "EPCAM/AGR2/KRT8/MUC2/TFF3/REG4 epithelial marker score dominates",
    11: "TRBC2/TRAC/CCL5/CD3D/CD2/TRBC1 T-cell marker score dominates",
    12: "PECAM1/PLVAP/FLT1/VWF/KDR/RAMP2 endothelial marker score dominates",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Prepare 10x Human Colorectal Cancer low-resolution data as a standard Stage1 real sample."
    )
    p.add_argument("--project_root", default=".")
    p.add_argument(
        "--source_dir",
        default="data/raw/low_resolution/Human Colorectal Cancer",
        help="Directory containing Cell Ranger sc/Flex and Visium outputs.",
    )
    p.add_argument("--target_sample", default="human_colorectal_cancer_real")
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def _decode_arr(arr: np.ndarray) -> list[str]:
    out: list[str] = []
    for x in arr.tolist():
        out.append(x.decode("utf-8") if isinstance(x, (bytes, bytearray)) else str(x))
    return out


def _make_unique_gene_names(feature_ids: list[str], gene_names: list[str]) -> tuple[list[str], dict[str, int]]:
    seen: dict[str, int] = {}
    duplicate_counts: dict[str, int] = {}
    unique: list[str] = []
    for feature_id, gene_name in zip(feature_ids, gene_names, strict=True):
        n = seen.get(gene_name, 0) + 1
        seen[gene_name] = n
        if n == 1:
            unique.append(gene_name)
        else:
            duplicate_counts[gene_name] = n
            unique.append(f"{gene_name}__dup{n}_{feature_id}")
    return unique, duplicate_counts


def _read_10x_h5(path: Path) -> tuple[list[str], list[str], list[str], sparse.csc_matrix]:
    with h5py.File(path, "r") as f:
        g = f["matrix"]
        shape = tuple(int(x) for x in g["shape"][()])
        data = g["data"][()]
        indices = g["indices"][()]
        indptr = g["indptr"][()]
        barcodes = _decode_arr(g["barcodes"][()])
        feature_ids = _decode_arr(g["features"]["id"][()])
        genes_raw = _decode_arr(g["features"]["name"][()])
        genes, duplicate_counts = _make_unique_gene_names(feature_ids, genes_raw)
    if duplicate_counts:
        print(f"[INFO] renamed duplicate gene symbols in {path.name}: {duplicate_counts}")
    mat = sparse.csc_matrix((data, indices, indptr), shape=shape)
    return feature_ids, genes, barcodes, mat


def _write_expr_tsv(path: Path, genes: list[str], sample_ids: list[str], mat_g_by_s: sparse.spmatrix) -> None:
    mat_csr = mat_g_by_s.tocsr()
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["Gene", *sample_ids])
        for i, gene in enumerate(genes):
            row = np.asarray(mat_csr.getrow(i).toarray()).ravel()
            writer.writerow([gene, *row.tolist()])
            if (i + 1) % 1000 == 0 or i + 1 == len(genes):
                print(f"[WRITE] {path.name}: {i + 1}/{len(genes)}")


def _write_dataset_yaml(project_root: Path, sample: str) -> None:
    cfg = {
        "paths": {
            "sc_expr": "brca_scRNA_GEP.txt",
            "sc_meta": "brca_scRNA_celllabels.txt",
            "st_expr": "brca_STdata_GEP.txt",
            "st_meta": "brca_STdata_coordinates.txt",
            "svg_marker_whitelist": None,
        },
        "qc": {
            "sc_min_genes": 100,
            "sc_max_genes": 10000,
            "sc_max_mt": 15,
            "st_min_genes": 50,
            "st_max_genes": "Inf",
            "st_max_mt": 30,
            "hvg_nfeatures": 2000,
            "mt_pattern": "^(MT-|mt-)",
        },
        "gene_filter": {"min_cells_sc": 0, "min_cells_st": 0},
        "stage3": {
            "strong_th": 0.7,
            "weak_th": 0.4,
            "st_cluster_k": 30,
            "unknown_floor": 0.3,
            "min_cells_rare_type": 20,
            "eps": 1.0e-8,
            "plugin_genes_path": f"data/processed/{sample}/stage1_preprocess/hvg_genes.txt",
            "gene_weights_path": None,
            "auto_missing_detection": {
                "enable": True,
                "method": "adaptive_low_support",
                "min_cells": 50,
                "robust_z_th": -2.0,
                "soft_z_th": -0.9,
                "require_masked_for_soft": False,
                "require_masked_for_hard": True,
                "require_confirmation": True,
                "confirmation_max_support_score": None,
                "confirmation_use_masked_missing": True,
                "confirmation_use_marker_identity": True,
                "confirmation_marker_identity_z_th": -1.5,
                "confirmation_marker_support_score_th": 0.5,
                "max_fraction_types": 0.3,
                "max_types": 2,
                "action": "mark_unknown",
            },
            "masked_missing_detection": {
                "enable": True,
                "apply_to_auto_missing": True,
                "neighbor_cosine_th": 0.9,
                "neighbor_cell_ratio_min": 1.0,
                "marker_top_n": 20,
                "min_identity_markers": 3,
                "min_marker_specificity": 1.2,
                "min_marker_type_mean": 0.001,
                "min_marker_st_detect_frac": 0.005,
                "st_presence_quantile": 0.9,
                "identity_z_th": -0.8,
                "pressure_z_th": 1.0,
                "min_support_score_for_apply": 0.8,
                "max_types": 2,
            },
            "marker_identity_diagnostics": {
                "enable": True,
                "marker_top_n": 80,
                "min_identity_markers": 5,
                "min_all_specificity": 1.3,
                "min_neighbor_specificity": 1.1,
                "min_marker_type_mean": 0.001,
                "min_marker_st_detect_frac": 0.005,
                "st_presence_quantile": 0.9,
                "depleted_z_th": -0.9,
            },
        },
    }
    out = project_root / "configs" / "datasets" / f"{sample}.yaml"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(yaml.safe_dump(cfg, sort_keys=False, allow_unicode=True), encoding="utf-8")
    print(f"[CFG] wrote: {out}")


def main() -> int:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    source_dir = (project_root / args.source_dir).resolve()
    target_dir = project_root / "data" / "raw" / args.target_sample
    if target_dir.exists():
        if not args.overwrite:
            raise FileExistsError(f"Target exists: {target_dir} (use --overwrite)")
        shutil.rmtree(target_dir)
    target_dir.mkdir(parents=True, exist_ok=True)

    sc_h5 = source_dir / "4plex_human_colorectal_kidney_scFFPE_multiplex_Colorectal_Manual_BC1_count_sample_filtered_feature_bc_matrix.h5"
    st_h5 = source_dir / "CytAssist_11mm_FFPE_Human_Colorectal_Cancer_filtered_feature_bc_matrix.h5"
    clusters_csv = source_dir / "analysis" / "clustering" / "gene_expression_graphclust" / "clusters.csv"
    positions_csv = source_dir / "spatial" / "tissue_positions.csv"
    for path in (sc_h5, st_h5, clusters_csv, positions_csv):
        if not path.exists():
            raise FileNotFoundError(path)

    print("[STEP] load 10x h5 matrices")
    sc_feature_ids, sc_genes, sc_cells, sc_mat = _read_10x_h5(sc_h5)
    st_feature_ids, st_genes, st_spots, st_mat = _read_10x_h5(st_h5)

    clusters = pd.read_csv(clusters_csv)
    clusters["Barcode"] = clusters["Barcode"].astype(str)
    clusters["Cluster"] = clusters["Cluster"].astype(int)
    cell_meta = pd.DataFrame({"cell_id": sc_cells})
    cell_meta = cell_meta.merge(clusters, left_on="cell_id", right_on="Barcode", how="left")
    if cell_meta["Cluster"].isna().any():
        missing = int(cell_meta["Cluster"].isna().sum())
        raise ValueError(f"{missing} sc cells are missing graphclust cluster labels.")
    cell_meta["Cluster"] = cell_meta["Cluster"].astype(int)
    cell_meta["cell_type"] = cell_meta["Cluster"].map(CLUSTER_TO_CELL_TYPE)
    if cell_meta["cell_type"].isna().any():
        bad = sorted(cell_meta.loc[cell_meta["cell_type"].isna(), "Cluster"].unique().tolist())
        raise ValueError(f"Unmapped clusters: {bad}")

    print("[STEP] align genes")
    st_feature_to_idx = {feature_id: i for i, feature_id in enumerate(st_feature_ids)}
    common_feature_ids = [feature_id for feature_id in sc_feature_ids if feature_id in st_feature_to_idx]
    sc_feature_to_idx = {feature_id: i for i, feature_id in enumerate(sc_feature_ids)}
    common_genes = [sc_genes[sc_feature_to_idx[feature_id]] for feature_id in common_feature_ids]
    if len(common_genes) != len(set(common_genes)):
        raise ValueError("Gene names must be unique after duplicate-symbol renaming.")
    sc_idx = np.array([sc_feature_to_idx[feature_id] for feature_id in common_feature_ids], dtype=np.int32)
    st_idx = np.array([st_feature_to_idx[feature_id] for feature_id in common_feature_ids], dtype=np.int32)
    if len(common_genes) < 500:
        raise ValueError(f"Too few common genes: {len(common_genes)}")
    sc_common = sc_mat[sc_idx, :]
    st_common = st_mat[st_idx, :]
    print(f"[INFO] common_genes={len(common_genes)} cells={sc_common.shape[1]} spots={st_common.shape[1]}")

    print("[STEP] write expression and metadata")
    _write_expr_tsv(target_dir / "brca_scRNA_GEP.txt", common_genes, sc_cells, sc_common)
    cell_meta[["cell_id", "cell_type"]].to_csv(target_dir / "brca_scRNA_celllabels.txt", sep="\t", index=False)

    positions = pd.read_csv(positions_csv)
    positions["barcode"] = positions["barcode"].astype(str)
    coord_df = positions.set_index("barcode").reindex(st_spots)
    if coord_df["pxl_row_in_fullres"].isna().any():
        missing = int(coord_df["pxl_row_in_fullres"].isna().sum())
        raise ValueError(f"{missing} ST spots are missing tissue positions.")
    pd.DataFrame(
        {
            "spot_id": st_spots,
            "row": coord_df["pxl_row_in_fullres"].to_numpy(dtype=int),
            "col": coord_df["pxl_col_in_fullres"].to_numpy(dtype=int),
        }
    ).to_csv(target_dir / "brca_STdata_coordinates.txt", sep="\t", index=False)
    _write_expr_tsv(target_dir / "brca_STdata_GEP.txt", common_genes, st_spots, st_common)

    audit = (
        cell_meta.groupby(["Cluster", "cell_type"], as_index=False)
        .size()
        .rename(columns={"size": "n_cells", "Cluster": "cluster"})
        .sort_values("cluster")
    )
    audit["annotation_evidence"] = audit["cluster"].map(CLUSTER_EVIDENCE)
    audit.to_csv(target_dir / "cell_type_annotation_audit.csv", index=False)

    info = {
        "sample": args.target_sample,
        "source_dir": str(source_dir),
        "n_cells": int(sc_common.shape[1]),
        "n_spots": int(st_common.shape[1]),
        "n_genes": int(len(common_genes)),
        "n_renamed_duplicate_gene_symbols": int(len(common_genes) - len({g.split('__dup', 1)[0] for g in common_genes})),
        "annotation_source": "Cell Ranger gene_expression_graphclust clusters mapped to broad cell types by canonical marker scores.",
        "cell_type_counts": cell_meta["cell_type"].value_counts().to_dict(),
        "params": vars(args),
    }
    (target_dir / "real_input_info.json").write_text(json.dumps(info, indent=2, ensure_ascii=False), encoding="utf-8")
    _write_dataset_yaml(project_root, args.target_sample)
    print(f"[DONE] prepared: {target_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
