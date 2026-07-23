from __future__ import annotations

import argparse
import gzip
import io as pyio
import json
import tarfile
import tempfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import yaml
from scipy import io, sparse


GROUP = "cytospace_fig2d_tme"

PAIR_SPECS = {
    "brca_er_her2_fresh_frozen": {
        "kind": "brca",
        "sc_subtype": "HER2+",
        "st_h5": "brca_ST_10x_fresh_frozen_block_a_section1/V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5",
        "st_spatial_tar": "brca_ST_10x_fresh_frozen_block_a_section1/V1_Breast_Cancer_Block_A_Section_1_spatial.tar.gz",
        "tumor_label": "Epithelial cells",
    },
    "brca_her2_ffpe": {
        "kind": "brca",
        "sc_subtype": "HER2+",
        "st_h5": "brca_ST_10x_ffpe/Visium_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix (1).h5",
        "st_spatial_tar": "brca_ST_10x_ffpe/Visium_FFPE_Human_Breast_Cancer_spatial.tar.gz",
        "tumor_label": "Epithelial cells",
    },
    "brca_tnbc_fresh_frozen": {
        "kind": "brca",
        "sc_subtype": "TNBC",
        "st_mex_tar": "brca_ST_Wu_TNBC_Zenodo/filtered_count_matrices.tar.gz",
        "st_mex_prefix": "filtered_count_matrices/CID4465_filtered_count_matrix",
        "st_spatial_tar": "brca_ST_Wu_TNBC_Zenodo/spatial.tar.gz",
        "st_spatial_prefix": "spatial/CID4465_spatial",
        "tumor_label": "Epithelial cells",
    },
    "crc_fresh_frozen": {
        "kind": "crc",
        "st_h5": "crc_ST_10x_parent_visium/Parent_Visium_Human_ColorectalCancer_filtered_feature_bc_matrix.h5",
        "st_spatial_tar": "crc_ST_10x_parent_visium/Parent_Visium_Human_ColorectalCancer_spatial.tar.gz",
        "tumor_label": "Epithelial cells",
    },
}


BRCA_LABEL_MAP = {
    "B-cells": "B cells",
    "Plasmablasts": "Plasma cells",
    "Myeloid": "Monocytes and Macrophages",
    "Endothelial": "Endothelial cells",
    "CAFs": "Fibroblasts",
    "PVL": "PVL",
    "Normal Epithelial": "Epithelial cells",
    "Cancer Epithelial": "Epithelial cells",
}

CRC_LABEL_MAP = {
    "B cells": "B cells",
    "Epithelial cells": "Epithelial cells",
    "Mast cells": "Mast cells",
    "Myeloids": "Monocytes and Macrophages",
    "Stromal cells": "Fibroblasts",
}


def _read_10x_h5(path: Path) -> tuple[sparse.csc_matrix, list[str], list[str]]:
    with h5py.File(path, "r") as f:
        g = f["matrix"]
        mat = sparse.csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=tuple(g["shape"][:]))
        genes = [x.decode() if isinstance(x, bytes) else str(x) for x in g["features"]["name"][:]]
        barcodes = [x.decode() if isinstance(x, bytes) else str(x) for x in g["barcodes"][:]]
    return mat, genes, barcodes


def _read_10x_mex_from_tar(tar_path: Path, prefix: str) -> tuple[sparse.csc_matrix, list[str], list[str]]:
    with tarfile.open(tar_path, "r:gz") as tf:
        def read_member(name: str) -> bytes:
            handle = tf.extractfile(name)
            if handle is None:
                raise FileNotFoundError(name)
            data = handle.read()
            return gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data

        mat = io.mmread(pyio.BytesIO(read_member(f"{prefix}/matrix.mtx.gz"))).tocsc()
        features = pd.read_csv(pyio.BytesIO(read_member(f"{prefix}/features.tsv.gz")), sep="\t", header=None)
        barcodes = pd.read_csv(pyio.BytesIO(read_member(f"{prefix}/barcodes.tsv.gz")), sep="\t", header=None)[0].astype(str).tolist()
    genes = features.iloc[:, 1 if features.shape[1] > 1 else 0].astype(str).tolist()
    return mat, genes, barcodes


def _read_spatial_from_tar(tar_path: Path, suffix: str | None = None) -> pd.DataFrame:
    with tarfile.open(tar_path, "r:gz") as tf:
        members = [m.name for m in tf.getmembers()]
        candidates = [m for m in members if m.endswith("tissue_positions_list.csv") or m.endswith("tissue_positions.csv")]
        if suffix:
            candidates = [m for m in candidates if m.startswith(suffix)]
        if not candidates:
            raise FileNotFoundError(f"No tissue_positions file in {tar_path}")
        f = tf.extractfile(candidates[0])
        coords = pd.read_csv(f, header=None)
    coords = coords.iloc[:, :6]
    coords.columns = ["spot_id", "in_tissue", "array_row", "array_col", "row", "col"]
    return coords[["spot_id", "row", "col"]]


def _collapse_duplicate_genes(mat: sparse.spmatrix, genes: list[str]) -> tuple[sparse.csr_matrix, list[str]]:
    gene_index: dict[str, list[int]] = {}
    for i, gene in enumerate(genes):
        gene_index.setdefault(str(gene), []).append(i)
    unique = sorted(gene_index)
    rows = [sparse.csr_matrix(mat[idxs, :].sum(axis=0)) for idxs in (gene_index[g] for g in unique)]
    return sparse.vstack(rows, format="csr"), unique


def _log_norm_cells(mat: sparse.spmatrix) -> sparse.csr_matrix:
    mat = mat.astype(np.float32).tocsc(copy=True)
    totals = np.asarray(mat.sum(axis=0)).ravel()
    totals[totals <= 0] = 1.0
    mat = mat @ sparse.diags((10000.0 / totals).astype(np.float32), format="csc")
    mat.data = np.log1p(mat.data)
    return mat.tocsr()


def _select_hvgs(sc_norm: sparse.spmatrix, genes: list[str], n: int) -> list[str]:
    if n <= 0 or n >= len(genes):
        return list(genes)
    mean = np.asarray(sc_norm.mean(axis=1)).ravel()
    mean_sq = np.asarray(sc_norm.multiply(sc_norm).mean(axis=1)).ravel()
    var = mean_sq - mean * mean
    idx = np.argsort(-var)[:n]
    return sorted([genes[i] for i in idx])


def _downsample_cells(meta: pd.DataFrame, max_per_type: int, seed: int) -> pd.DataFrame:
    if max_per_type <= 0:
        return meta
    rng = np.random.default_rng(seed)
    keep = []
    for _, sub in meta.groupby("cell_type", sort=False):
        if len(sub) <= max_per_type:
            keep.extend(sub.index.tolist())
        else:
            keep.extend(rng.choice(sub.index.to_numpy(), size=max_per_type, replace=False).tolist())
    return meta.loc[keep].sort_index()


def _load_brca_sc(root: Path, subtype: str, max_per_type: int, seed: int):
    extracted = root / "data/raw/cytospace_fig2d_tme/brca_scRNA_GSE176078/Wu_etal_2021_BRCA_scRNASeq"
    if extracted.exists():
        print(f"[STEP] reading extracted BRCA scRNA: {extracted}", flush=True)
        meta = pd.read_csv(extracted / "metadata.csv", index_col=0)
        genes = (extracted / "count_matrix_genes.tsv").read_text(encoding="utf-8").splitlines()
        barcodes = (extracted / "count_matrix_barcodes.tsv").read_text(encoding="utf-8").splitlines()
        mat = io.mmread(extracted / "count_matrix_sparse.mtx").tocsc()
    else:
        print("[STEP] reading compressed BRCA scRNA tar", flush=True)
        tar_path = root / "data/raw/cytospace_fig2d_tme/brca_scRNA_GSE176078/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz"
        with tarfile.open(tar_path, "r:gz") as tf:
            meta = pd.read_csv(tf.extractfile("Wu_etal_2021_BRCA_scRNASeq/metadata.csv"), index_col=0)
            genes = [line.decode().strip() for line in tf.extractfile("Wu_etal_2021_BRCA_scRNASeq/count_matrix_genes.tsv")]
            barcodes = [line.decode().strip() for line in tf.extractfile("Wu_etal_2021_BRCA_scRNASeq/count_matrix_barcodes.tsv")]
            mat = io.mmread(tf.extractfile("Wu_etal_2021_BRCA_scRNASeq/count_matrix_sparse.mtx")).tocsc()
    meta = meta.loc[barcodes].copy()
    meta = meta[meta["subtype"].astype(str).eq(subtype)].copy()
    labels = []
    for _, row in meta.iterrows():
        major = str(row["celltype_major"])
        minor = str(row["celltype_minor"])
        if minor == "T cells CD4+":
            labels.append("CD4 T cells")
        elif minor == "T cells CD8+":
            labels.append("CD8 T cells")
        elif minor in {"NK cells", "NKT cells"}:
            labels.append("NK cells")
        else:
            labels.append(BRCA_LABEL_MAP.get(major, major))
    meta["cell_type"] = labels
    meta = _downsample_cells(meta, max_per_type, seed)
    idx = [barcodes.index(x) for x in meta.index.astype(str)]
    return mat[:, idx], genes, meta.index.astype(str).tolist(), meta[["cell_type"]].reset_index(names="cell_id")


def _load_crc_sc(root: Path, max_per_type: int, seed: int):
    base = root / "data/raw/cytospace_fig2d_tme/crc_scRNA_GSE132465"
    anno = pd.read_csv(base / "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz", sep="\t", index_col=0)
    labels = []
    for _, row in anno.iterrows():
        cell_type = str(row["Cell_type"])
        subtype = str(row["Cell_subtype"])
        if subtype == "CD4+ T cells":
            labels.append("CD4 T cells")
        elif subtype == "CD8+ T cells":
            labels.append("CD8 T cells")
        elif subtype == "NK cells":
            labels.append("NK cells")
        else:
            labels.append(CRC_LABEL_MAP.get(cell_type, cell_type))
    anno["cell_type"] = labels
    anno = _downsample_cells(anno, max_per_type, seed)
    use_cells = set(anno.index.astype(str))
    header = pd.read_csv(base / "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz", sep="\t", nrows=0).columns.tolist()
    cells = header[1:]
    keep_cols = [i for i, c in enumerate(cells) if c in use_cells]
    usecols = [0] + [i + 1 for i in keep_cols]
    expr = pd.read_csv(
        base / "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz",
        sep="\t",
        index_col=0,
        usecols=usecols,
    )
    genes = expr.index.astype(str).tolist()
    kept_cells = expr.columns.astype(str).tolist()
    meta = anno.loc[kept_cells, ["cell_type"]].reset_index(names="cell_id")
    return sparse.csr_matrix(expr.to_numpy(dtype=np.float32)), genes, kept_cells, meta


def _load_st(root: Path, spec: dict):
    base = root / "data/raw/cytospace_fig2d_tme"
    if "st_h5" in spec:
        mat, genes, spots = _read_10x_h5(base / spec["st_h5"])
        coords = _read_spatial_from_tar(base / spec["st_spatial_tar"])
    else:
        mat, genes, spots = _read_10x_mex_from_tar(base / spec["st_mex_tar"], spec["st_mex_prefix"])
        coords = _read_spatial_from_tar(base / spec["st_spatial_tar"], spec["st_spatial_prefix"])
    return mat, genes, spots, coords


def _write_config(root: Path, sample: str, hvg_path: Path) -> None:
    cfg = {
        "paths": {"sc_expr": "unused", "sc_meta": None, "st_expr": None, "st_meta": None, "svg_marker_whitelist": None},
        "qc": {
            "sc_min_genes": 0,
            "sc_max_genes": 100000,
            "sc_max_mt": 100,
            "st_min_genes": 0,
            "st_max_genes": "Inf",
            "st_max_mt": 100,
            "hvg_nfeatures": 2000,
            "mt_pattern": "^(MT-|mt-)",
        },
        "gene_filter": {"min_cells_sc": 0, "min_cells_st": 0},
        "stage3": {
            "strong_th": 0.7,
            "weak_th": 0.4,
            "st_cluster_k": 0,
            "unknown_floor": 0.3,
            "min_cells_rare_type": 20,
            "eps": 1.0e-8,
            "plugin_genes_path": str(hvg_path.relative_to(root)).replace("\\", "/"),
            "gene_weights_path": None,
            "auto_missing_detection": {"enable": False, "method": "adaptive_low_support", "action": "mark_unknown"},
            "masked_missing_detection": {"enable": False},
            "marker_identity_diagnostics": {"enable": True, "marker_top_n": 80},
        },
        "storage": {"group": GROUP},
    }
    path = root / "configs/datasets" / f"{sample}.yaml"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")


def prepare_pair(root: Path, pair_id: str, max_genes: int, max_per_type: int, seed: int, overwrite: bool) -> None:
    spec = PAIR_SPECS[pair_id]
    print(f"[STEP] preparing pair: {pair_id}", flush=True)
    sample = f"cytospace_fig2d_tme_{pair_id}"
    stage1 = root / "data/processed" / GROUP / sample / "stage1_preprocess"
    export = stage1 / "exported"
    if export.exists() and not overwrite:
        raise FileExistsError(f"{export} exists; pass --overwrite")
    export.mkdir(parents=True, exist_ok=True)

    if spec["kind"] == "brca":
        sc_mat, sc_genes, sc_cells, sc_meta = _load_brca_sc(root, spec["sc_subtype"], max_per_type, seed)
    elif spec["kind"] == "crc":
        sc_mat, sc_genes, sc_cells, sc_meta = _load_crc_sc(root, max_per_type, seed)
    else:
        raise ValueError(spec["kind"])
    st_mat, st_genes, spots, coords = _load_st(root, spec)

    sc_mat, sc_genes = _collapse_duplicate_genes(sc_mat, sc_genes)
    st_mat, st_genes = _collapse_duplicate_genes(st_mat, st_genes)
    sc_lookup = {g: i for i, g in enumerate(sc_genes)}
    st_lookup = {g: i for i, g in enumerate(st_genes)}
    common = sorted(set(sc_lookup) & set(st_lookup))
    if len(common) < 100:
        raise ValueError(f"{pair_id}: too few common genes: {len(common)}")
    sc_idx = [sc_lookup[g] for g in common]
    st_idx = [st_lookup[g] for g in common]
    sc_norm = _log_norm_cells(sc_mat[sc_idx, :])
    st_norm = _log_norm_cells(st_mat[st_idx, :])
    hvgs = _select_hvgs(sc_norm, common, max_genes)
    hvg_set = set(hvgs)
    keep_common_idx = [i for i, g in enumerate(common) if g in hvg_set]
    genes = [common[i] for i in keep_common_idx]
    sc_norm = sc_norm[keep_common_idx, :]
    st_norm = st_norm[keep_common_idx, :]

    pd.DataFrame(sc_norm.T.toarray(), index=sc_cells, columns=genes).to_csv(export / "sc_expression_normalized.csv", index_label="cell_id")
    pd.DataFrame(st_norm.T.toarray(), index=spots, columns=genes).to_csv(export / "st_expression_normalized.csv", index_label="spot_id")
    sc_meta.to_csv(export / "sc_metadata.csv", index=False)
    coords = coords[coords["spot_id"].astype(str).isin(set(spots))].set_index("spot_id").reindex(spots)
    coords.to_csv(export / "st_coordinates.csv", index_label="spot_id")
    (stage1 / "common_genes.txt").write_text("\n".join(genes) + "\n", encoding="utf-8")
    (stage1 / "hvg_genes.txt").write_text("\n".join(genes) + "\n", encoding="utf-8")
    _write_config(root, sample, stage1 / "hvg_genes.txt")
    summary = {
        "pair_id": pair_id,
        "sample": sample,
        "n_cells": len(sc_cells),
        "n_spots": len(spots),
        "n_genes": len(genes),
        "cell_type_counts": sc_meta["cell_type"].value_counts().to_dict(),
        "tumor_label": spec["tumor_label"],
    }
    (stage1 / "fig2d_stage1_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"[OK] prepared {sample}: cells={len(sc_cells)} spots={len(spots)} genes={len(genes)}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare Stage1 exports for CytoSPACE Fig.2d TME route2 pairs.")
    parser.add_argument("--project_root", default=".")
    parser.add_argument("--pair", action="append", choices=sorted(PAIR_SPECS), default=[])
    parser.add_argument("--max_genes", type=int, default=2000)
    parser.add_argument("--max_sc_cells_per_type", type=int, default=1500)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    root = Path(args.project_root).resolve()
    pairs = args.pair or list(PAIR_SPECS)
    for pair_id in pairs:
        prepare_pair(root, pair_id, args.max_genes, args.max_sc_cells_per_type, args.seed, args.overwrite)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
