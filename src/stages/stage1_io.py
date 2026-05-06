"""
Shared loader for Stage1 exports. Normalizes to canonical format:
  sc_expr: cells x genes
  st_expr: spots x genes
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from src.stages.storage import read_dataset_config, stage1_export_dir


def _clean_spot_index(idx: pd.Index) -> pd.Index:
    return pd.Index([str(x).split("\t")[0] for x in idx])


def _resolve_sc_expr_path(base: Path, sc_expr_source: str) -> Path:
    source = str(sc_expr_source or "normalized").strip().lower()
    if source in ("counts", "count"):
        candidates = ["sc_expression_counts.csv", "sc_expression_normalized.csv"]
    elif source in ("data", "normalized_data"):
        candidates = ["sc_expression_data.csv", "sc_expression_normalized.csv"]
    elif source in ("auto", "best"):
        candidates = [
            "sc_expression_counts.csv",
            "sc_expression_data.csv",
            "sc_expression_normalized.csv",
        ]
    else:
        candidates = ["sc_expression_normalized.csv"]

    for name in candidates:
        path = base / name
        if path.exists():
            return path
    raise FileNotFoundError(
        f"Missing SC expression file under {base}. Tried: {', '.join(candidates)}"
    )


def load_stage1(
    root: Path,
    sample: str,
    sc_expr_source: str = "normalized",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load Stage1 exports. Normalizes to canonical format: sc_expr=cells x genes, st_expr=spots x genes."""
    cfg = read_dataset_config(root, sample)
    base = stage1_export_dir(root, sample, cfg)
    sc_expr_path = _resolve_sc_expr_path(base, sc_expr_source)
    sc_expr = pd.read_csv(sc_expr_path, index_col=0, sep=None, engine="python")
    st_expr = pd.read_csv(base / "st_expression_normalized.csv", index_col=0, sep=None, engine="python")
    st_coords = pd.read_csv(base / "st_coordinates.csv", index_col=0, sep=None, engine="python")
    sc_meta = pd.read_csv(base / "sc_metadata.csv", sep=None, engine="python")
    if "cell_id" in sc_expr.columns:
        sc_expr = sc_expr.drop(columns=["cell_id"])
    if "cell_id" in st_expr.columns:
        st_expr = st_expr.drop(columns=["cell_id"])
    sc_expr = sc_expr.apply(pd.to_numeric, errors="coerce").dropna(axis=1, how="all").astype("float32")
    st_expr = st_expr.apply(pd.to_numeric, errors="coerce").dropna(axis=1, how="all").astype("float32")
    st_expr.index = _clean_spot_index(st_expr.index)
    st_coords.index = _clean_spot_index(st_coords.index)

    # Normalize orientation to cells x genes and spots x genes for Stage3/Stage4 compatibility.
    cell_ids = set(sc_meta["cell_id"].astype(str))
    n_cells_in_idx = sum(1 for c in cell_ids if c in sc_expr.index)
    n_cells_in_col = sum(1 for c in cell_ids if c in sc_expr.columns)
    if n_cells_in_col > n_cells_in_idx:
        sc_expr = sc_expr.T
    spot_ids = set(st_coords.index.astype(str))
    n_spots_in_idx = sum(1 for s in spot_ids if s in st_expr.index)
    n_spots_in_col = sum(1 for s in spot_ids if s in st_expr.columns)
    if n_spots_in_col > n_spots_in_idx:
        st_expr = st_expr.T
        st_expr.index = _clean_spot_index(st_expr.index)
    return sc_expr, st_expr, st_coords, sc_meta
