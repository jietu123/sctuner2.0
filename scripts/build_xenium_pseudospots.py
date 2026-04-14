from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build minimal Xenium pseudo-spots for SVTuner Stage1 inputs")
    p.add_argument("--sample", default="xenium_minimal")
    p.add_argument("--project_root", default=".")
    p.add_argument(
        "--xenium_dir",
        default="data/raw/xenium_outs",
    )
    p.add_argument("--output_dir", default=None)
    p.add_argument("--gene_list", default=None, help="Optional gene list file. Defaults to Xenium panel genes.")
    p.add_argument("--grid_size_um", type=float, default=250.0)
    p.add_argument("--min_cells_per_spot", type=int, default=25)
    p.add_argument("--max_spots", type=int, default=1800)
    p.add_argument(
        "--spot_keep_strategy",
        choices=["top_cells", "all"],
        default="top_cells",
        help="If max_spots is exceeded, keep the densest spots or keep all.",
    )
    p.add_argument("--write_cell_map", action="store_true")
    p.add_argument(
        "--cell_annotation_csv",
        default=None,
        help=(
            "Optional cell-level annotation CSV/CSV.GZ containing `cell_id` and type columns, "
            "e.g. data/sim/<sample>/preview_truth/xenium_cell_major_type_truth.csv.gz"
        ),
    )
    p.add_argument(
        "--annotation_type_column",
        default="pred_major_type",
        help="Type column in --cell_annotation_csv used for optional dropout filtering.",
    )
    p.add_argument(
        "--drop_cell_type",
        default=None,
        help="If set with --cell_annotation_csv, randomly drop this type before pseudo-spot aggregation.",
    )
    p.add_argument(
        "--drop_fraction",
        type=float,
        default=1.0,
        help="Fraction of --drop_cell_type cells to drop (0-1, default 1.0 means full dropout).",
    )
    p.add_argument("--seed", type=int, default=42, help="Random seed for dropout sampling.")
    return p.parse_args()


def resolve_path(project_root: Path, value: str | None) -> Path | None:
    if value is None:
        return None
    p = Path(value)
    return p if p.is_absolute() else (project_root / p).resolve()


def load_panel_genes(xenium_dir: Path) -> list[str]:
    panel = json.loads((xenium_dir / "gene_panel.json").read_text(encoding="utf-8"))
    genes: list[str] = []
    for target in panel.get("payload", {}).get("targets", []):
        t = target.get("type", {})
        if t.get("descriptor") == "gene":
            name = t.get("data", {}).get("name")
            if name:
                genes.append(str(name))
    seen = set()
    ordered: list[str] = []
    for gene in genes:
        if gene not in seen:
            ordered.append(gene)
            seen.add(gene)
    return ordered


def load_cells_table(xenium_dir: Path) -> pd.DataFrame:
    cols = ["cell_id", "x_centroid", "y_centroid", "transcript_counts", "cell_area", "nucleus_area"]
    parquet_path = xenium_dir / "cells.parquet"
    csv_path = xenium_dir / "cells.csv.gz"
    if parquet_path.exists():
        try:
            return pd.read_parquet(parquet_path, columns=cols)
        except Exception:
            pass
    if csv_path.exists():
        return pd.read_csv(csv_path, usecols=cols, compression="infer")
    raise FileNotFoundError(f"Neither cells.parquet nor cells.csv.gz found under {xenium_dir}")


def load_h5_matrix(h5_path: Path) -> tuple[sparse.csc_matrix, np.ndarray, np.ndarray]:
    with h5py.File(h5_path, "r") as f:
        grp = f["matrix"]
        data = grp["data"][:]
        indices = grp["indices"][:]
        indptr = grp["indptr"][:]
        shape = tuple(int(x) for x in grp["shape"][:])
        barcodes = grp["barcodes"][:]
        features = grp["features"]["name"][:]
    mat = sparse.csc_matrix((data, indices, indptr), shape=shape, dtype=np.float32)
    barcode_ids = np.array([x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in barcodes], dtype=object)
    feature_names = np.array([x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in features], dtype=object)
    return mat, barcode_ids, feature_names


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    xenium_dir = resolve_path(project_root, args.xenium_dir)
    if xenium_dir is None or not xenium_dir.exists():
        raise FileNotFoundError(f"Xenium dir not found: {xenium_dir}")

    output_dir = resolve_path(
        project_root,
        args.output_dir or f"data/processed/{args.sample}/stage1_preprocess/exported",
    )
    assert output_dir is not None
    output_dir.mkdir(parents=True, exist_ok=True)
    stage1_dir = output_dir.parent
    stage1_dir.mkdir(parents=True, exist_ok=True)

    if args.gene_list:
        gene_list_path = resolve_path(project_root, args.gene_list)
        if gene_list_path is None or not gene_list_path.exists():
            raise FileNotFoundError(f"gene_list not found: {gene_list_path}")
        gene_list = [line.strip() for line in gene_list_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    else:
        gene_list = load_panel_genes(xenium_dir)

    cells = load_cells_table(xenium_dir)
    if cells.empty:
        raise ValueError("cells.parquet is empty")
    n_cells_total = int(len(cells))
    dropout_stats = {
        "enabled": False,
        "annotation_path": None,
        "type_column": None,
        "drop_cell_type": None,
        "drop_fraction": None,
        "target_cells_before": 0,
        "dropped_cells": 0,
        "target_cells_after": 0,
    }

    if args.cell_annotation_csv and args.drop_cell_type:
        ann_path = resolve_path(project_root, args.cell_annotation_csv)
        if ann_path is None or not ann_path.exists():
            raise FileNotFoundError(f"cell_annotation_csv not found: {ann_path}")
        ann = pd.read_csv(ann_path, compression="infer")
        if "cell_id" not in ann.columns:
            raise ValueError("cell_annotation_csv must include column `cell_id`")
        if args.annotation_type_column not in ann.columns:
            raise ValueError(
                f"annotation_type_column '{args.annotation_type_column}' not found in {ann_path}"
            )
        ann = ann[["cell_id", args.annotation_type_column]].copy()
        ann["cell_id"] = ann["cell_id"].astype(str)
        ann = ann.drop_duplicates(subset=["cell_id"], keep="first")

        cells = cells.copy()
        cells["cell_id"] = cells["cell_id"].astype(str)
        cells = cells.merge(ann, on="cell_id", how="left")

        drop_fraction = float(args.drop_fraction)
        drop_fraction = min(max(drop_fraction, 0.0), 1.0)
        target_mask = cells[args.annotation_type_column].astype(str) == str(args.drop_cell_type)
        target_idx = cells.index[target_mask].to_numpy()
        n_target_before = int(target_mask.sum())
        n_drop = int(round(n_target_before * drop_fraction))

        if n_drop > 0:
            rng = np.random.default_rng(int(args.seed))
            drop_idx = rng.choice(target_idx, size=n_drop, replace=False)
            cells = cells.drop(index=drop_idx).copy()

        n_target_after = int((cells[args.annotation_type_column].astype(str) == str(args.drop_cell_type)).sum())
        cells = cells.drop(columns=[args.annotation_type_column])

        dropout_stats = {
            "enabled": True,
            "annotation_path": str(ann_path),
            "type_column": args.annotation_type_column,
            "drop_cell_type": str(args.drop_cell_type),
            "drop_fraction": float(drop_fraction),
            "target_cells_before": int(n_target_before),
            "dropped_cells": int(n_drop),
            "target_cells_after": int(n_target_after),
        }
        if cells.empty:
            raise ValueError("No cells left after dropout filtering; lower drop_fraction or change drop_cell_type.")

    xmin = float(cells["x_centroid"].min())
    ymin = float(cells["y_centroid"].min())
    grid = float(args.grid_size_um)

    gx = np.floor((cells["x_centroid"].to_numpy(dtype=np.float64) - xmin) / grid).astype(np.int32)
    gy = np.floor((cells["y_centroid"].to_numpy(dtype=np.float64) - ymin) / grid).astype(np.int32)
    cells["grid_x"] = gx
    cells["grid_y"] = gy
    cells["pseudo_spot_id"] = pd.Series(gx).astype(str).radd("gx") + "_gy" + pd.Series(gy).astype(str)

    spot_stats = (
        cells.groupby("pseudo_spot_id", sort=False)
        .agg(
            cell_count=("cell_id", "size"),
            row=("grid_y", "mean"),
            col=("grid_x", "mean"),
            row_um=("y_centroid", "mean"),
            col_um=("x_centroid", "mean"),
            transcript_sum=("transcript_counts", "sum"),
        )
        .reset_index()
    )
    spot_stats = spot_stats.loc[spot_stats["cell_count"] >= int(args.min_cells_per_spot)].copy()
    if args.spot_keep_strategy == "top_cells" and args.max_spots and len(spot_stats) > int(args.max_spots):
        spot_stats = spot_stats.nlargest(int(args.max_spots), "cell_count").copy()
    selected_spots = set(spot_stats["pseudo_spot_id"].astype(str))
    cells = cells.loc[cells["pseudo_spot_id"].isin(selected_spots)].copy()
    spot_stats = spot_stats.loc[spot_stats["pseudo_spot_id"].isin(selected_spots)].copy()
    spot_stats = spot_stats.sort_values(["row_um", "col_um"], kind="mergesort").reset_index(drop=True)

    matrix, barcode_ids, feature_names = load_h5_matrix(xenium_dir / "cell_feature_matrix.h5")
    feature_index = {name: i for i, name in enumerate(feature_names.tolist())}
    keep_genes = [g for g in gene_list if g in feature_index]
    if not keep_genes:
        raise ValueError("No requested genes found in Xenium matrix")

    cell_to_spot = pd.Series(cells["pseudo_spot_id"].to_numpy(), index=cells["cell_id"].astype(str))
    spot_order = spot_stats["pseudo_spot_id"].astype(str).tolist()
    spot_to_col = {spot: idx for idx, spot in enumerate(spot_order)}

    assigned_spots = pd.Series(barcode_ids, copy=False).map(cell_to_spot)
    keep_mask = assigned_spots.notna().to_numpy()
    sel_cols = np.flatnonzero(keep_mask)
    if len(sel_cols) == 0:
        raise ValueError("No Xenium cells mapped to selected pseudo-spots")

    selected_spot_ids = assigned_spots.iloc[sel_cols].astype(str).to_numpy()
    selected_spot_cols = np.fromiter((spot_to_col[s] for s in selected_spot_ids), dtype=np.int32, count=len(selected_spot_ids))

    indicator = sparse.csr_matrix(
        (
            np.ones(len(sel_cols), dtype=np.float32),
            (np.arange(len(sel_cols), dtype=np.int32), selected_spot_cols),
        ),
        shape=(len(sel_cols), len(spot_order)),
    )

    gene_idx = np.array([feature_index[g] for g in keep_genes], dtype=np.int32)
    sub = matrix[gene_idx][:, sel_cols]
    aggregated = (sub @ indicator).astype(np.float32).toarray().T

    st_expr = pd.DataFrame(aggregated, index=spot_order, columns=keep_genes)
    st_coords = spot_stats.set_index("pseudo_spot_id")[["row_um", "col_um", "cell_count", "transcript_sum"]].copy()
    st_coords = st_coords.rename(columns={"row_um": "row", "col_um": "col"})
    st_coords = st_coords.reindex(spot_order)

    st_expr.to_csv(output_dir / "st_expression_normalized.csv")
    st_coords.to_csv(output_dir / "st_coordinates.csv", index_label="spot_id")
    (stage1_dir / "hvg_genes.txt").write_text("\n".join(keep_genes) + "\n", encoding="utf-8")
    (stage1_dir / "common_genes.txt").write_text("\n".join(keep_genes) + "\n", encoding="utf-8")

    spot_stats.to_csv(output_dir / "pseudo_spot_stats.csv", index=False)
    if args.write_cell_map:
        cells[["cell_id", "x_centroid", "y_centroid", "pseudo_spot_id", "transcript_counts"]].to_csv(
            output_dir / "xenium_cell_to_pseudospot.csv",
            index=False,
        )

    summary = {
        "sample": args.sample,
        "xenium_dir": str(xenium_dir),
        "output_dir": str(output_dir),
        "grid_size_um": float(args.grid_size_um),
        "min_cells_per_spot": int(args.min_cells_per_spot),
        "max_spots": int(args.max_spots),
        "spot_keep_strategy": args.spot_keep_strategy,
        "n_cells_total": n_cells_total,
        "n_cells_selected": int(len(cells)),
        "n_pseudospots": int(len(spot_order)),
        "n_genes": int(len(keep_genes)),
        "selected_genes": keep_genes,
        "dropout": dropout_stats,
    }
    (stage1_dir / "xenium_pseudospot_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(f"[Xenium] spots={len(spot_order)} cells={len(cells)} genes={len(keep_genes)}")
    print(f"[Xenium] output_dir={output_dir}")


if __name__ == "__main__":
    main()
