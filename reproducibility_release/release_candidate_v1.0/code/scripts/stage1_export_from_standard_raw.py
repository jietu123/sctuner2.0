#!/usr/bin/env python
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from src.stages.storage import raw_dir, stage1_dir


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Export Stage1-compatible normalized CSVs from standard raw TSV inputs without R/Seurat."
    )
    p.add_argument("--project_root", default=".")
    p.add_argument("--sample", required=True)
    p.add_argument("--hvg_nfeatures", type=int, default=2000)
    p.add_argument(
        "--export_gene_mode",
        choices=("hvg", "all_common"),
        default="hvg",
        help="Default hvg keeps exported CSVs small enough for Stage3 on large h5ad-derived datasets.",
    )
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def _read_yaml(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(path)
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def _resolve_raw_path(project_root: Path, sample: str, cfg: dict, rel_path: str) -> Path:
    p = Path(rel_path)
    if p.is_absolute():
        return p
    return raw_dir(project_root, sample, cfg) / p


def _read_expr(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    if df.index.has_duplicates:
        raise ValueError(f"Duplicate gene names are not allowed: {path}")
    return df.astype(np.float32, copy=False)


def _log_normalize_cells(mat: pd.DataFrame, scale_factor: float = 10000.0) -> pd.DataFrame:
    arr = mat.to_numpy(dtype=np.float32, copy=True)
    if np.nanmin(arr) < 0:
        raise ValueError("Expression matrix contains negative values; use raw counts/nonnegative expression.")
    totals = arr.sum(axis=0, keepdims=True)
    totals[totals <= 0] = 1.0
    arr = arr / totals * scale_factor
    np.log1p(arr, out=arr)
    return pd.DataFrame(arr, index=mat.index, columns=mat.columns)


def _select_hvg(sc_norm: pd.DataFrame, n: int) -> list[str]:
    if n <= 0 or n >= sc_norm.shape[0]:
        return sc_norm.index.astype(str).tolist()
    variances = sc_norm.var(axis=1)
    return variances.sort_values(ascending=False).head(n).index.astype(str).tolist()


def _write_transposed_csv(path: Path, mat: pd.DataFrame, index_label: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    mat.T.to_csv(path, index_label=index_label)


def main() -> int:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    cfg_path = project_root / "configs" / "datasets" / f"{args.sample}.yaml"
    cfg = _read_yaml(cfg_path)
    paths = cfg.get("paths", {})

    sc_expr_path = _resolve_raw_path(project_root, args.sample, cfg, paths.get("sc_expr", "brca_scRNA_GEP.txt"))
    sc_meta_path = _resolve_raw_path(project_root, args.sample, cfg, paths.get("sc_meta", "brca_scRNA_celllabels.txt"))
    st_expr_path = _resolve_raw_path(project_root, args.sample, cfg, paths.get("st_expr", "brca_STdata_GEP.txt"))
    st_meta_path = _resolve_raw_path(project_root, args.sample, cfg, paths.get("st_meta", "brca_STdata_coordinates.txt"))

    out_dir = stage1_dir(project_root, args.sample, cfg)
    if out_dir.exists() and args.overwrite:
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    export_dir = out_dir / "exported"
    export_dir.mkdir(parents=True, exist_ok=True)

    print(f"[STEP] read sc expression: {sc_expr_path}")
    sc = _read_expr(sc_expr_path)
    print(f"[STEP] read st expression: {st_expr_path}")
    st = _read_expr(st_expr_path)

    st_gene_set = set(st.index.astype(str).tolist())
    common_genes = [g for g in sc.index.astype(str).tolist() if g in st_gene_set]
    if len(common_genes) < 100:
        raise ValueError(f"Too few common genes: {len(common_genes)}")
    sc = sc.loc[common_genes]
    st = st.loc[common_genes]

    print("[STEP] log-normalize sc/st")
    sc_norm = _log_normalize_cells(sc)
    st_norm = _log_normalize_cells(st)

    hvg_genes = _select_hvg(sc_norm, min(args.hvg_nfeatures, len(common_genes)))
    export_genes = hvg_genes if args.export_gene_mode == "hvg" else common_genes
    sc_norm = sc_norm.loc[export_genes]
    st_norm = st_norm.loc[export_genes]

    print(f"[STEP] export normalized CSVs: genes={len(export_genes)} mode={args.export_gene_mode}")
    _write_transposed_csv(export_dir / "sc_expression_normalized.csv", sc_norm, "cell_id")
    _write_transposed_csv(export_dir / "st_expression_normalized.csv", st_norm, "spot_id")

    sc_meta = pd.read_csv(sc_meta_path, sep="\t")
    if list(sc_meta.columns[:2]) != ["cell_id", "cell_type"]:
        sc_meta = sc_meta.iloc[:, :2]
        sc_meta.columns = ["cell_id", "cell_type"]
    sc_meta.to_csv(export_dir / "sc_metadata.csv", index=False)

    st_coords = pd.read_csv(st_meta_path, sep="\t")
    if list(st_coords.columns[:3]) != ["spot_id", "row", "col"]:
        st_coords = st_coords.iloc[:, :3]
        st_coords.columns = ["spot_id", "row", "col"]
    st_coords.to_csv(export_dir / "st_coordinates.csv", index=False)

    (out_dir / "common_genes.txt").write_text("\n".join(export_genes) + "\n", encoding="utf-8")
    (out_dir / "hvg_genes.txt").write_text("\n".join(hvg_genes) + "\n", encoding="utf-8")

    print(f"[DONE] exported Python Stage1-compatible files: {out_dir}")
    print(
        f"[INFO] common_genes={len(common_genes)} export_genes={len(export_genes)} "
        f"hvg_genes={len(hvg_genes)} sc_cells={sc_norm.shape[1]} st_spots={st_norm.shape[1]}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
