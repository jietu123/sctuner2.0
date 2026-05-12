#!/usr/bin/env python
from __future__ import annotations

import argparse
import contextlib
import json
import sys
import threading
import time
import traceback
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

try:
    import psutil
except Exception:  # pragma: no cover
    psutil = None

_HERE = Path(__file__).resolve()
_ROOT = _HERE.parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from scripts.run_similarity_mapping import (  # noqa: E402
    _common_genes,
    _eval_against_truth,
    _json_write,
    _load_inputs,
    _processed_export_dir,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Run a CellTrek-style spatial charting baseline on one simulation sample. "
            "This Python runner is used when the R CellTrek/Seurat stack is unavailable."
        )
    )
    p.add_argument("--project_root", default=".", help="Project root.")
    p.add_argument("--group", required=True, help="Simulation group.")
    p.add_argument("--sample", required=True, help="Simulation sample id.")
    p.add_argument("--cell_type_column", default="cell_type", help="Cell type column in sc_metadata.csv.")
    p.add_argument("--max_genes", type=int, default=2000, help="Maximum variable common genes used for PCA.")
    p.add_argument("--n_pcs", type=int, default=30, help="Number of PCA dimensions.")
    p.add_argument("--ntree", type=int, default=500, help="Number of random forest trees.")
    p.add_argument(
        "--cells_per_spot",
        type=int,
        default=0,
        help="Number of charted cells summarized per spot. 0 means round(n_cells / n_spots).",
    )
    p.add_argument("--out_dir", default=None, help="Override output directory.")
    return p.parse_args()


class ResourceMonitor:
    def __init__(self, interval_seconds: float = 0.5) -> None:
        self.interval_seconds = interval_seconds
        self.process = psutil.Process() if psutil is not None else None
        self.peak_rss_mb = 0.0
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)

    def start(self) -> None:
        if self.process is not None:
            self._thread.start()

    def stop(self) -> float:
        if self.process is None:
            return self.peak_rss_mb
        self._stop.set()
        self._thread.join(timeout=2.0)
        try:
            self.peak_rss_mb = max(self.peak_rss_mb, self.process.memory_info().rss / (1024**2))
        except Exception:
            pass
        return self.peak_rss_mb

    def _run(self) -> None:
        while not self._stop.is_set():
            try:
                self.peak_rss_mb = max(self.peak_rss_mb, self.process.memory_info().rss / (1024**2))
            except Exception:
                pass
            self._stop.wait(self.interval_seconds)


def _load_coordinates(export_dir: Path, st_index: pd.Index) -> pd.DataFrame:
    coord_path = export_dir / "st_coordinates.csv"
    if not coord_path.exists():
        raise FileNotFoundError(f"required input missing: {coord_path}")
    coords = pd.read_csv(coord_path)
    spot_cols = [c for c in coords.columns if c == "spot_id" or c.startswith("spot_id.")]
    if not spot_cols:
        raise ValueError("st_coordinates.csv must contain a spot_id column.")
    row_col = "row" if "row" in coords.columns else "spatial_row"
    col_col = "col" if "col" in coords.columns else "spatial_col"
    if row_col not in coords.columns or col_col not in coords.columns:
        raise ValueError("st_coordinates.csv must contain row/col or spatial_row/spatial_col columns.")
    coords = coords[[spot_cols[0], row_col, col_col]].copy()
    coords.columns = ["spot_id", "coord_x", "coord_y"]
    coords["spot_id"] = coords["spot_id"].astype(str)
    coords = coords.set_index("spot_id")
    common = st_index.astype(str).intersection(coords.index)
    if common.empty:
        raise ValueError("no common spots between ST expression and coordinates.")
    return coords.loc[common].astype(float)


def _select_variable_genes(sc_expr: pd.DataFrame, st_expr: pd.DataFrame, genes: list[str], max_genes: int) -> list[str]:
    if len(genes) <= max_genes:
        return genes
    sc_var = sc_expr[genes].astype("float32").var(axis=0)
    st_var = st_expr[genes].astype("float32").var(axis=0)
    score = (sc_var.rank(pct=True) + st_var.rank(pct=True)).sort_values(ascending=False)
    return score.index[:max_genes].tolist()


def _spot_type_fraction_from_cells(
    spot_coords: pd.DataFrame,
    cell_coords: pd.DataFrame,
    cell_types: pd.Series,
    cells_per_spot: int,
) -> pd.DataFrame:
    cell_types = cell_types.astype(str)
    type_order = sorted(cell_types.unique().tolist())
    k = max(1, min(int(cells_per_spot), len(cell_coords)))
    nn = NearestNeighbors(n_neighbors=k, algorithm="auto")
    nn.fit(cell_coords[["coord_x", "coord_y"]].to_numpy(dtype=float))
    _, idx = nn.kneighbors(spot_coords[["coord_x", "coord_y"]].to_numpy(dtype=float))

    out = pd.DataFrame(0.0, index=spot_coords.index.astype(str), columns=type_order)
    labels = cell_types.to_numpy()
    for i, spot_id in enumerate(out.index):
        vals, counts = np.unique(labels[idx[i]], return_counts=True)
        out.loc[spot_id, vals] = counts.astype(float) / float(k)
    out.index.name = "spot_id"
    return out


def _cell_assignment_from_coordinates(spot_coords: pd.DataFrame, cell_coords: pd.DataFrame) -> pd.DataFrame:
    nn = NearestNeighbors(n_neighbors=1, algorithm="auto")
    nn.fit(spot_coords[["coord_x", "coord_y"]].to_numpy(dtype=float))
    dist, idx = nn.kneighbors(cell_coords[["coord_x", "coord_y"]].to_numpy(dtype=float))
    spot_ids = spot_coords.index.astype(str).to_numpy()
    out = cell_coords[["cell_id", "cell_type"]].copy()
    out["assigned_spot"] = spot_ids[idx[:, 0]]
    out["assignment_distance"] = dist[:, 0]
    return out[["cell_id", "assigned_spot", "cell_type", "assignment_distance"]]


def _write_mapping_outputs(
    frac: pd.DataFrame,
    cell_coords: pd.DataFrame,
    cell_assignment: pd.DataFrame,
    out_dir: Path,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    frac.to_csv(out_dir / "spot_type_fraction.csv")
    dominant = pd.DataFrame(
        {
            "spot_id": frac.index.astype(str),
            "dominant_type": frac.idxmax(axis=1).to_numpy(),
            "dominant_score": frac.max(axis=1).to_numpy(),
        }
    )
    dominant.to_csv(out_dir / "spot_type_dominant.csv", index=False)
    cell_coords.to_csv(out_dir / "celltrek_cell_coordinates.csv", index=False)
    cell_assignment.to_csv(out_dir / "cell_assignment.csv", index=False)


def main() -> int:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    out_dir = (
        Path(args.out_dir).resolve()
        if args.out_dir
        else project_root / "result" / args.sample / "stage4_mapping" / "celltrek"
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "run.log"
    start = time.perf_counter()
    monitor = ResourceMonitor()
    monitor.start()
    status = "failed"

    try:
        with log_path.open("w", encoding="utf-8") as log_fh, contextlib.redirect_stdout(log_fh), contextlib.redirect_stderr(log_fh):
            export_dir = _processed_export_dir(project_root, args.group, args.sample)
            sc_expr, st_expr, sc_meta, sim_info = _load_inputs(export_dir, args.cell_type_column)
            cell_type_column = str(sc_meta.attrs.get("cell_type_column", args.cell_type_column))
            st_coords = _load_coordinates(export_dir, st_expr.index)
            st_expr = st_expr.loc[st_coords.index]

            genes = _select_variable_genes(sc_expr, st_expr, _common_genes(sc_expr, st_expr), args.max_genes)
            x_sc = sc_expr[genes].apply(pd.to_numeric, errors="coerce").fillna(0.0).astype("float32")
            x_st = st_expr[genes].apply(pd.to_numeric, errors="coerce").fillna(0.0).astype("float32")
            x_all = pd.concat([x_st, x_sc], axis=0)

            scaler = StandardScaler(with_mean=True, with_std=True)
            x_scaled = scaler.fit_transform(x_all.to_numpy(dtype=np.float32))
            n_pcs = max(2, min(args.n_pcs, x_scaled.shape[0] - 1, x_scaled.shape[1]))
            pca = PCA(n_components=n_pcs, svd_solver="randomized", random_state=42)
            emb = pca.fit_transform(x_scaled)
            st_emb = emb[: len(x_st), :]
            sc_emb = emb[len(x_st) :, :]

            rf = RandomForestRegressor(
                n_estimators=args.ntree,
                random_state=42,
                n_jobs=1,
                min_samples_leaf=2,
            )
            rf.fit(st_emb, st_coords[["coord_x", "coord_y"]].to_numpy(dtype=float))
            pred_cell_coords = rf.predict(sc_emb)
            cell_coords = pd.DataFrame(
                {
                    "cell_id": x_sc.index.astype(str),
                    "cell_type": sc_meta.loc[x_sc.index, cell_type_column].astype(str).to_numpy(),
                    "coord_x": pred_cell_coords[:, 0],
                    "coord_y": pred_cell_coords[:, 1],
                }
            )

            cps = args.cells_per_spot if args.cells_per_spot > 0 else int(round(len(x_sc) / max(len(x_st), 1)))
            frac = _spot_type_fraction_from_cells(
                st_coords,
                cell_coords.set_index("cell_id")[["coord_x", "coord_y"]],
                sc_meta.loc[x_sc.index, cell_type_column],
                cells_per_spot=max(1, cps),
            )
            cell_assignment = _cell_assignment_from_coordinates(st_coords, cell_coords)
            _write_mapping_outputs(frac, cell_coords, cell_assignment, out_dir)
            (out_dir / "genes_used.txt").write_text("\n".join(genes) + "\n", encoding="utf-8")
            metrics = _eval_against_truth(project_root, args.group, args.sample, frac, sim_info)
            _json_write(out_dir / "metrics_simulation.json", metrics)
            _json_write(
                out_dir / "celltrek_summary.json",
                {
                    "sample": args.sample,
                    "group": args.group,
                    "method": "celltrek",
                    "implementation": "python_fallback_celltrek_core",
                    "status": "ok",
                    "note": (
                        "Official R CellTrek source was downloaded to external/CellTrek, but the current "
                        "Windows conda R/Seurat stack fails to load with a Mingw-w64 runtime error. "
                        "This runner keeps the CellTrek spatial-charting idea: expression PCA -> RF coordinate "
                        "regression on ST -> chart sc cells back to ST space -> spot-level type fractions."
                    ),
                    "processed_export_dir": str(export_dir),
                    "cell_type_column": cell_type_column,
                    "n_cells": int(sc_expr.shape[0]),
                    "n_spots": int(st_expr.shape[0]),
                    "n_genes_used": int(len(genes)),
                    "n_pcs": int(n_pcs),
                    "ntree": int(args.ntree),
                    "cells_per_spot": int(max(1, cps)),
                    "output_files": {
                        "spot_type_fraction": str(out_dir / "spot_type_fraction.csv"),
                        "spot_type_dominant": str(out_dir / "spot_type_dominant.csv"),
                        "cell_assignment": str(out_dir / "cell_assignment.csv"),
                        "celltrek_cell_coordinates": str(out_dir / "celltrek_cell_coordinates.csv"),
                        "metrics_simulation": str(out_dir / "metrics_simulation.json"),
                    },
                },
            )
            status = "ok"
    except Exception as exc:
        _json_write(
            out_dir / "celltrek_summary.json",
            {
                "sample": args.sample,
                "group": args.group,
                "method": "celltrek",
                "status": "failed",
                "error_type": type(exc).__name__,
                "error": str(exc),
                "traceback": traceback.format_exc(),
            },
        )
        raise
    finally:
        peak = monitor.stop()
        elapsed = time.perf_counter() - start
        _json_write(
            out_dir / "resource_metrics.json",
            {
                "sample": args.sample,
                "group": args.group,
                "method": "celltrek",
                "status": status,
                "wall_time_seconds": float(elapsed),
                "peak_memory_mb": float(peak),
                "log_path": str(log_path),
            },
        )

    print(f"[OK] celltrek {status}: {args.sample} -> {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
