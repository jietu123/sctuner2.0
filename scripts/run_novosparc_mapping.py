#!/usr/bin/env python
from __future__ import annotations

import argparse
import contextlib
import json
import os
import sys
import threading
import time
import traceback
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

try:
    import psutil
except Exception:  # pragma: no cover
    psutil = None

_HERE = Path(__file__).resolve()
_ROOT = _HERE.parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from scripts.run_celltrek_mapping import _load_coordinates  # noqa: E402
from scripts.run_similarity_mapping import (  # noqa: E402
    _common_genes,
    _eval_against_truth,
    _json_write,
    _load_inputs,
    _processed_export_dir,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Run novoSpaRc on one simulation sample and export method-agnostic mapping outputs."
    )
    p.add_argument("--project_root", default=".", help="Project root.")
    p.add_argument("--group", required=True, help="Simulation group.")
    p.add_argument("--sample", required=True, help="Simulation sample id.")
    p.add_argument("--cell_type_column", default="cell_type", help="Cell type column in sc_metadata.csv.")
    p.add_argument("--max_genes", type=int, default=2000, help="Maximum variable common genes used.")
    p.add_argument("--n_pcs", type=int, default=30, help="PCA dimensions for novoSpaRc smooth expression cost.")
    p.add_argument("--num_neighbors_s", type=int, default=5, help="Source-cell kNN size.")
    p.add_argument("--num_neighbors_t", type=int, default=5, help="Target-location kNN size.")
    p.add_argument("--alpha_linear", type=float, default=0.5, help="Weight of ST atlas expression cost.")
    p.add_argument("--epsilon", type=float, default=5e-3, help="Entropy regularization.")
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


def _ensure_novosparc_available() -> None:
    if os.environ.get("PYTHONNOUSERSITE") == "1":
        raise RuntimeError(
            "PYTHONNOUSERSITE=1 is set, but novosparc was installed in the user site-packages. "
            "Unset PYTHONNOUSERSITE before running novoSpaRc."
        )
    try:
        import novosparc  # noqa: F401
    except Exception as exc:
        raise RuntimeError("Package 'novosparc' is not available in the active Python environment.") from exc


def _select_variable_genes(sc_expr: pd.DataFrame, st_expr: pd.DataFrame, genes: list[str], max_genes: int) -> list[str]:
    if len(genes) <= max_genes:
        return genes
    sc_var = sc_expr[genes].apply(pd.to_numeric, errors="coerce").fillna(0.0).astype("float32").var(axis=0)
    st_var = st_expr[genes].apply(pd.to_numeric, errors="coerce").fillna(0.0).astype("float32").var(axis=0)
    score = (sc_var.rank(pct=True) + st_var.rank(pct=True)).sort_values(ascending=False)
    return score.index[:max_genes].tolist()


def _pca_representation(x: pd.DataFrame, n_pcs: int) -> np.ndarray:
    arr = x.apply(pd.to_numeric, errors="coerce").fillna(0.0).astype("float32").to_numpy()
    n_components = max(2, min(int(n_pcs), arr.shape[0] - 1, arr.shape[1]))
    arr = StandardScaler(with_mean=True, with_std=True).fit_transform(arr)
    return PCA(n_components=n_components, random_state=42).fit_transform(arr).astype("float32")


def _fraction_from_coupling(gw: np.ndarray, cell_types: pd.Series, spot_ids: pd.Index) -> pd.DataFrame:
    labels = cell_types.astype(str).to_numpy()
    type_order = sorted(pd.unique(labels).tolist())
    onehot = np.zeros((len(labels), len(type_order)), dtype="float64")
    type_to_idx = {t: i for i, t in enumerate(type_order)}
    for i, label in enumerate(labels):
        onehot[i, type_to_idx[label]] = 1.0
    frac = gw.T @ onehot
    frac = frac / np.maximum(frac.sum(axis=1, keepdims=True), 1e-12)
    out = pd.DataFrame(frac, index=spot_ids.astype(str), columns=type_order)
    out.index.name = "spot_id"
    return out


def _cell_assignment_from_coupling(
    gw: np.ndarray,
    cell_ids: pd.Index,
    cell_types: pd.Series,
    spot_ids: pd.Index,
) -> pd.DataFrame:
    best_idx = np.argmax(gw, axis=1)
    row_sum = np.maximum(gw.sum(axis=1), 1e-12)
    score = gw[np.arange(gw.shape[0]), best_idx] / row_sum
    return pd.DataFrame(
        {
            "cell_id": cell_ids.astype(str),
            "assigned_spot": spot_ids.astype(str).to_numpy()[best_idx],
            "cell_type": cell_types.astype(str).to_numpy(),
            "assignment_score": score,
        }
    )


def _write_mapping_outputs(
    frac: pd.DataFrame,
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
    cell_assignment.to_csv(out_dir / "cell_assignment.csv", index=False)


def main() -> int:
    args = parse_args()
    _ensure_novosparc_available()
    import novosparc

    project_root = Path(args.project_root).resolve()
    out_dir = (
        Path(args.out_dir).resolve()
        if args.out_dir
        else project_root / "result" / args.sample / "stage4_mapping" / "novosparc"
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
            coords = _load_coordinates(export_dir, st_expr.index)
            common_spots = st_expr.index.astype(str).intersection(coords.index.astype(str))
            if common_spots.empty:
                raise ValueError("no common spots between ST expression and coordinates.")
            st_expr = st_expr.loc[common_spots]
            coords = coords.loc[common_spots]

            genes = _select_variable_genes(sc_expr, st_expr, _common_genes(sc_expr, st_expr), args.max_genes)
            sc_x = sc_expr[genes].apply(pd.to_numeric, errors="coerce").fillna(0.0).astype("float32")
            st_x = st_expr[genes].apply(pd.to_numeric, errors="coerce").fillna(0.0).astype("float32")
            adata = AnnData(X=sc_x.to_numpy(dtype="float32"))
            adata.obs_names = sc_x.index.astype(str)
            adata.var_names = genes
            locations = coords[["coord_x", "coord_y"]].to_numpy(dtype="float64")
            atlas_matrix = st_x.to_numpy(dtype="float64")
            markers_to_use = np.arange(len(genes), dtype=int)

            tissue = novosparc.cm.Tissue(
                dataset=adata,
                locations=locations,
                atlas_matrix=atlas_matrix,
                markers_to_use=markers_to_use,
                output_folder=str(out_dir),
            )
            dge_rep = _pca_representation(sc_x, args.n_pcs)
            tissue.setup_linear_cost(markers_to_use=markers_to_use, atlas_matrix=atlas_matrix)
            tissue.setup_smooth_costs(
                dge_rep=dge_rep,
                num_neighbors_s=args.num_neighbors_s,
                num_neighbors_t=args.num_neighbors_t,
                verbose=True,
            )
            tissue.reconstruct(
                alpha_linear=float(args.alpha_linear),
                epsilon=float(args.epsilon),
                search_epsilon=True,
                random_ini=False,
                verbose=True,
            )

            gw = np.asarray(tissue.gw, dtype="float64")
            frac = _fraction_from_coupling(gw, sc_meta.loc[sc_x.index, cell_type_column], st_x.index)
            cell_assignment = _cell_assignment_from_coupling(
                gw,
                sc_x.index,
                sc_meta.loc[sc_x.index, cell_type_column],
                st_x.index,
            )
            _write_mapping_outputs(frac, cell_assignment, out_dir)
            (out_dir / "genes_used.txt").write_text("\n".join(genes) + "\n", encoding="utf-8")
            metrics = _eval_against_truth(project_root, args.group, args.sample, frac, sim_info)
            _json_write(out_dir / "metrics_simulation.json", metrics)
            _json_write(
                out_dir / "novosparc_summary.json",
                {
                    "sample": args.sample,
                    "group": args.group,
                    "method": "novosparc",
                    "status": "ok",
                    "processed_export_dir": str(export_dir),
                    "cell_type_column": cell_type_column,
                    "n_cells": int(sc_x.shape[0]),
                    "n_spots": int(st_x.shape[0]),
                    "n_genes_used": int(len(genes)),
                    "n_cell_types": int(sc_meta.loc[sc_x.index, cell_type_column].nunique()),
                    "alpha_linear": float(args.alpha_linear),
                    "epsilon": float(tissue.epsilon),
                    "num_neighbors_s": int(args.num_neighbors_s),
                    "num_neighbors_t": int(args.num_neighbors_t),
                    "output_files": {
                        "spot_type_fraction": str(out_dir / "spot_type_fraction.csv"),
                        "spot_type_dominant": str(out_dir / "spot_type_dominant.csv"),
                        "cell_assignment": str(out_dir / "cell_assignment.csv"),
                        "metrics_simulation": str(out_dir / "metrics_simulation.json"),
                    },
                },
            )
            status = "ok"
    except Exception as exc:
        _json_write(
            out_dir / "novosparc_summary.json",
            {
                "sample": args.sample,
                "group": args.group,
                "method": "novosparc",
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
                "method": "novosparc",
                "status": status,
                "wall_time_seconds": float(elapsed),
                "peak_memory_mb": float(peak),
                "log_path": str(log_path),
            },
        )

    print(f"[OK] novosparc {status}: {args.sample} -> {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
