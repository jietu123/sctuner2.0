#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import re
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

from src.stages.storage import read_dataset_config, result_dir, stage1_dir
from src.utils.sample_paths import resolve_sample_dir


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Prepare Stage1 exported files for a simulation-derived sample "
            "without running R (fallback for Windows Rscript runtime issues)."
        )
    )
    p.add_argument("--sample", required=True, help="target sample id")
    p.add_argument("--project_root", default=".", help="project root")
    p.add_argument("--source_sample", default=None, help="optional source sample id override")
    p.add_argument("--scale_factor", type=float, default=10000.0, help="log-normalize scale factor")
    return p.parse_args()


def _load_source_sample(raw_dir: Path, override: str | None) -> str | None:
    if override:
        return override
    sim_info = raw_dir / "sim_info.json"
    if not sim_info.exists():
        required = [
            raw_dir / "brca_scRNA_GEP.txt",
            raw_dir / "brca_scRNA_celllabels.txt",
            raw_dir / "brca_STdata_GEP.txt",
            raw_dir / "brca_STdata_coordinates.txt",
        ]
        if all(p.exists() for p in required):
            return None
        raise FileNotFoundError(f"sim_info.json not found and raw text inputs are incomplete: {sim_info}")
    info = json.loads(sim_info.read_text(encoding="utf-8"))
    source_sample = str(info.get("source_sample") or "").strip()
    if not source_sample:
        raise ValueError("sim_info.json missing non-empty source_sample")
    return source_sample


def _copy_sc_exports(src_export: Path, dst_export: Path) -> None:
    copied = []
    for name in ("sc_metadata.csv", "sc_expression_normalized.csv", "sc_expression_data.csv", "sc_expression_counts.csv"):
        s = src_export / name
        if s.exists():
            shutil.copy2(s, dst_export / name)
            copied.append(name)
    if not (dst_export / "sc_metadata.csv").exists():
        raise FileNotFoundError(f"required source file missing: {src_export / 'sc_metadata.csv'}")
    if not (dst_export / "sc_expression_normalized.csv").exists():
        raise FileNotFoundError(f"required source file missing: {src_export / 'sc_expression_normalized.csv'}")
    print(f"[Stage1-fallback] copied SC exports: {copied}")


def _load_qc_config(project_root: Path, sample: str) -> dict:
    defaults = {
        "sc_min_genes": 200,
        "sc_max_genes": 6000,
        "sc_max_mt": 10,
        "st_min_genes": 100,
        "st_max_genes": float("inf"),
        "st_max_mt": 20,
        "hvg_nfeatures": 2000,
        "mt_pattern": "^MT-",
    }
    cfg_path = project_root / "configs" / "datasets" / f"{sample}.yaml"
    if not cfg_path.exists():
        return defaults
    with cfg_path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}
    out = defaults.copy()
    out.update(cfg.get("qc") or {})
    return out


def _mt_mask(genes: list[str], pattern: str) -> np.ndarray:
    rx = re.compile(pattern or "^MT-")
    return np.array([bool(rx.search(str(g))) for g in genes], dtype=bool)


def _build_sc_exports_from_raw(
    raw_sc_expr: Path,
    raw_sc_meta: Path,
    qc: dict,
    scale_factor: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str], list[str]]:
    print(f"[Stage1-fallback] building SC exports from raw: {raw_sc_expr}")
    sc_raw = pd.read_csv(raw_sc_expr, sep="\t", low_memory=False)
    if sc_raw.shape[1] < 2:
        raise ValueError(f"invalid SC expression file: {raw_sc_expr}")
    gene_col = sc_raw.columns[0]
    sc_raw = sc_raw.set_index(gene_col)
    sc_raw.index = sc_raw.index.astype(str)
    sc_raw.columns = [str(c) for c in sc_raw.columns]
    sc_raw = sc_raw.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    meta = pd.read_csv(raw_sc_meta, sep="\t")
    if meta.shape[1] < 2:
        raise ValueError(f"invalid SC metadata file: {raw_sc_meta}")
    meta = meta.iloc[:, :2].copy()
    meta.columns = ["cell_id", "cell_type"]
    meta["cell_id"] = meta["cell_id"].astype(str)
    meta["cell_type"] = meta["cell_type"].astype(str)
    meta = meta.drop_duplicates("cell_id")

    common_cells = [c for c in sc_raw.columns if c in set(meta["cell_id"])]
    if len(common_cells) == 0:
        raise ValueError("No overlapping cells between SC expression and SC metadata.")
    sc_raw = sc_raw[common_cells]
    genes = sc_raw.index.astype(str).tolist()

    mat = sc_raw.to_numpy(dtype=np.float32, copy=False)  # genes x cells
    lib = mat.sum(axis=0, dtype=np.float64)
    n_count = lib.astype(np.int64)
    n_feature = (mat > 0).sum(axis=0).astype(np.int64)
    mt_mask = _mt_mask(genes, str(qc.get("mt_pattern") or "^MT-"))
    if mt_mask.any():
        mt_sum = mat[mt_mask, :].sum(axis=0, dtype=np.float64)
        pct_mt = np.where(lib > 0, mt_sum / np.where(lib > 0, lib, 1.0) * 100.0, 0.0).astype(np.float32)
    else:
        pct_mt = np.zeros(shape=(mat.shape[1],), dtype=np.float32)

    keep_mask = (
        (n_feature >= int(qc.get("sc_min_genes", 200)))
        & (n_feature <= float(qc.get("sc_max_genes", 6000)))
        & (pct_mt <= float(qc.get("sc_max_mt", 10)))
    )
    kept_cells = [c for c, keep in zip(common_cells, keep_mask) if bool(keep)]
    if not kept_cells:
        raise ValueError("SC QC removed all cells in Stage1 fallback.")
    mat = mat[:, keep_mask]
    lib = lib[keep_mask]
    n_count = n_count[keep_mask]
    n_feature = n_feature[keep_mask]
    pct_mt = pct_mt[keep_mask]
    denom = np.where(lib > 0.0, lib, 1.0)
    norm = np.log1p((mat / denom[np.newaxis, :]) * float(scale_factor)).astype(np.float32)

    hvg_n = int(qc.get("hvg_nfeatures", 2000))
    if hvg_n > 0 and hvg_n < len(genes):
        vars_ = np.var(norm, axis=1)
        hvg_idx = np.argsort(vars_)[::-1][:hvg_n]
        hvg_genes = [genes[i] for i in hvg_idx]
    else:
        hvg_genes = list(genes)

    export_genes = hvg_genes if hvg_n > 0 and hvg_n < len(genes) else genes
    export_idx = np.array([genes.index(g) for g in export_genes], dtype=np.int32)
    norm_export = norm[export_idx, :]
    mat_export = mat[export_idx, :]

    meta = meta.set_index("cell_id").reindex(kept_cells)
    sc_meta = pd.DataFrame(
        {
            "orig.ident": "SeuratProject",
            "nCount_RNA": n_count,
            "nFeature_RNA": n_feature,
            "cell_type": meta["cell_type"].fillna("Unknown_sc_only").astype(str).to_numpy(),
            "percent.mt": pct_mt,
        },
        index=kept_cells,
    )
    sc_meta.index.name = "cell_id"
    sc_expr_norm = pd.DataFrame(norm_export.T, index=kept_cells, columns=export_genes)
    sc_expr_norm.index.name = "cell_id"

    sc_expr_counts = pd.DataFrame(mat_export.T, index=kept_cells, columns=export_genes)
    sc_expr_counts.index.name = "cell_id"

    print(
        f"[Stage1-fallback] SC QC kept {len(kept_cells)}/{len(common_cells)} cells; "
        f"export_genes={len(export_genes)}/{len(genes)}"
    )
    return sc_expr_norm, sc_expr_counts, sc_meta, export_genes, hvg_genes


def _read_gene_order_from_sc(sc_expr_csv: Path) -> list[str]:
    cols = pd.read_csv(sc_expr_csv, nrows=0).columns.tolist()
    if not cols or cols[0] != "cell_id":
        raise ValueError(f"unexpected sc_expression header: first column is not 'cell_id' in {sc_expr_csv}")
    genes = [str(c) for c in cols[1:]]
    if len(genes) == 0:
        raise ValueError(f"no gene columns found in {sc_expr_csv}")
    return genes


def _build_st_normalized(
    raw_st_expr: Path,
    gene_order: list[str],
    qc: dict,
    scale_factor: float,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray]:
    print(f"[Stage1-fallback] loading ST counts: {raw_st_expr}")
    st_raw = pd.read_csv(raw_st_expr, sep="\t", low_memory=False)
    if st_raw.shape[1] < 2:
        raise ValueError(f"invalid ST expression file: {raw_st_expr}")
    gene_col = st_raw.columns[0]
    st_raw = st_raw.set_index(gene_col)
    st_raw.index = st_raw.index.astype(str)
    st_raw.columns = [str(c) for c in st_raw.columns]
    st_raw = st_raw.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    # align to SC gene order; fill missing genes with zero
    st_raw = st_raw.reindex(gene_order).fillna(0.0)
    spot_ids = st_raw.columns.to_numpy(dtype=str)
    mat = st_raw.to_numpy(dtype=np.float32, copy=False)  # genes x spots
    del st_raw

    lib = mat.sum(axis=0, dtype=np.float64)
    n_count = lib.astype(np.int64)
    n_feature = (mat > 0).sum(axis=0).astype(np.int64)
    mt_mask = _mt_mask(gene_order, str(qc.get("mt_pattern") or "^MT-"))
    if mt_mask.any():
        mt_sum = mat[mt_mask, :].sum(axis=0, dtype=np.float64)
        pct_mt = np.where(lib > 0, mt_sum / np.where(lib > 0, lib, 1.0) * 100.0, 0.0).astype(np.float32)
    else:
        pct_mt = np.zeros(shape=(mat.shape[1],), dtype=np.float32)

    keep_mask = (
        (n_feature >= int(qc.get("st_min_genes", 100)))
        & (n_feature <= float(qc.get("st_max_genes", float("inf"))))
        & (pct_mt <= float(qc.get("st_max_mt", 20)))
    )
    if not bool(keep_mask.any()):
        raise ValueError("ST QC removed all spots in Stage1 fallback.")
    spot_ids = spot_ids[keep_mask]
    mat = mat[:, keep_mask]
    lib = lib[keep_mask]
    n_count = n_count[keep_mask]
    n_feature = n_feature[keep_mask]
    pct_mt = pct_mt[keep_mask]
    denom = np.where(lib > 0.0, lib, 1.0)
    norm = np.log1p((mat / denom[np.newaxis, :]) * float(scale_factor)).astype(np.float32)

    st_expr_norm = pd.DataFrame(norm.T, index=spot_ids, columns=gene_order)
    st_expr_norm.index.name = "spot_id"
    print(f"[Stage1-fallback] ST QC kept {len(spot_ids)}/{len(keep_mask)} spots")
    return st_expr_norm, n_count, n_feature, pct_mt


def _build_st_coords(raw_st_meta: Path, spot_ids: list[str], n_count: np.ndarray, n_feature: np.ndarray, pct_mt: np.ndarray) -> pd.DataFrame:
    coords = pd.read_csv(raw_st_meta, sep="\t")
    if coords.shape[1] < 3:
        raise ValueError(f"invalid ST coordinate file: {raw_st_meta}")
    coords = coords.iloc[:, :3].copy()
    coords.columns = ["spot_id", "row", "col"]
    coords["spot_id"] = coords["spot_id"].astype(str)
    coords = coords.set_index("spot_id").reindex(spot_ids)
    if coords[["row", "col"]].isna().any().any():
        missing = int(coords[["row", "col"]].isna().any(axis=1).sum())
        raise ValueError(f"spot coordinates missing for {missing} spots in {raw_st_meta}")

    coords["orig.ident"] = "SeuratProject"
    coords["nCount_RNA"] = n_count
    coords["nFeature_RNA"] = n_feature
    coords["spot_id"] = coords.index
    coords["percent.mt"] = pct_mt
    coords["UMI_total"] = n_count
    cols = ["orig.ident", "nCount_RNA", "nFeature_RNA", "spot_id", "row", "col", "percent.mt", "UMI_total"]
    return coords[cols]


def main() -> int:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    sample = args.sample

    raw_dir = resolve_sample_dir(project_root, sample, sim_group="real_brca", must_exist=True)
    source_sample = _load_source_sample(raw_dir, args.source_sample)

    cfg = read_dataset_config(project_root, sample)
    dst_stage1 = stage1_dir(project_root, sample, cfg)
    dst_export = dst_stage1 / "exported"
    dst_export.mkdir(parents=True, exist_ok=True)
    src_cfg = read_dataset_config(project_root, source_sample) if source_sample else {}
    src_stage1 = stage1_dir(project_root, source_sample, src_cfg) if source_sample else None
    src_export = src_stage1 / "exported" if src_stage1 else None
    if src_export is not None and src_export.exists() and src_export.resolve() != dst_export.resolve():
        _copy_sc_exports(src_export, dst_export)
        # keep Stage3 plugin gene path consistent
        for name in ("hvg_genes.txt", "common_genes.txt"):
            s = src_stage1 / name
            if s.exists():
                shutil.copy2(s, dst_stage1 / name)
    else:
        print(f"[Stage1-fallback] source stage1 export not found, rebuild SC from raw: {src_export}")
        raw_sc_expr = raw_dir / "brca_scRNA_GEP.txt"
        raw_sc_meta = raw_dir / "brca_scRNA_celllabels.txt"
        if not raw_sc_expr.exists() or not raw_sc_meta.exists():
            raise FileNotFoundError(
                "source stage1 missing and target raw SC files missing: "
                "brca_scRNA_GEP.txt / brca_scRNA_celllabels.txt"
            )
        qc = _load_qc_config(project_root, sample)
        sc_expr_norm, sc_expr_counts, sc_meta, genes, hvg_genes = _build_sc_exports_from_raw(
            raw_sc_expr,
            raw_sc_meta,
            qc,
            args.scale_factor,
        )
        sc_expr_norm_out = sc_expr_norm.copy()
        sc_expr_norm_out.insert(0, "cell_id", sc_expr_norm_out.index)
        sc_expr_norm_out.to_csv(dst_export / "sc_expression_normalized.csv", index=False)

        sc_expr_data_out = sc_expr_norm.copy()
        sc_expr_data_out.insert(0, "cell_id", sc_expr_data_out.index)
        sc_expr_data_out.to_csv(dst_export / "sc_expression_data.csv", index=False)

        sc_expr_counts_out = sc_expr_counts.copy()
        sc_expr_counts_out.insert(0, "cell_id", sc_expr_counts_out.index)
        sc_expr_counts_out.to_csv(dst_export / "sc_expression_counts.csv", index=False)

        sc_meta_out = sc_meta.copy()
        sc_meta_out.insert(0, "cell_id", sc_meta_out.index)
        sc_meta_out.to_csv(dst_export / "sc_metadata.csv", index=False)

        (dst_stage1 / "hvg_genes.txt").write_text("\n".join(hvg_genes) + "\n", encoding="utf-8")
        (dst_stage1 / "common_genes.txt").write_text("\n".join(genes) + "\n", encoding="utf-8")

    # keep simulation info available for stage4 summary hash and audit
    sim_info = raw_dir / "sim_info.json"
    if sim_info.exists():
        shutil.copy2(sim_info, dst_export / "sim_info.json")

    raw_st_expr = raw_dir / "brca_STdata_GEP.txt"
    raw_st_meta = raw_dir / "brca_STdata_coordinates.txt"
    if not raw_st_expr.exists() or not raw_st_meta.exists():
        raise FileNotFoundError("target raw ST files missing: brca_STdata_GEP.txt / brca_STdata_coordinates.txt")

    gene_order = _read_gene_order_from_sc(dst_export / "sc_expression_normalized.csv")
    qc = _load_qc_config(project_root, sample)
    st_expr_norm, n_count, n_feature, pct_mt = _build_st_normalized(raw_st_expr, gene_order, qc, args.scale_factor)
    st_expr_norm.to_csv(dst_export / "st_expression_normalized.csv")

    spot_ids = st_expr_norm.index.astype(str).tolist()
    st_coords = _build_st_coords(raw_st_meta, spot_ids, n_count, n_feature, pct_mt)
    st_coords.to_csv(dst_export / "st_coordinates.csv", index=True, index_label="spot_id")

    summary_dir = result_dir(project_root, sample, cfg) / "stage1_preprocess"
    summary_dir.mkdir(parents=True, exist_ok=True)
    summary = {
        "sample": sample,
        "mode": "python_fallback_from_sim_source",
        "source_sample": source_sample,
        "scale_factor": float(args.scale_factor),
        "n_spots": int(len(spot_ids)),
        "n_genes": int(len(gene_order)),
        "output_dir": str(dst_export),
    }
    (summary_dir / "stage1_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    print(f"[Stage1-fallback] done -> {dst_export}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
