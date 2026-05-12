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
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances
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
from scripts.run_novosparc_mapping import (  # noqa: E402
    _cell_assignment_from_coupling,
    _fraction_from_coupling,
    _select_variable_genes,
    _write_mapping_outputs,
)
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
            "Run a SpaOTsc-style structured optimal transport mapping on one simulation sample. "
            "The runner uses SpaOTsc's official usot solver while avoiding the incompatible top-level import."
        )
    )
    p.add_argument("--project_root", default=".", help="Project root.")
    p.add_argument("--group", required=True, help="Simulation group.")
    p.add_argument("--sample", required=True, help="Simulation sample id.")
    p.add_argument("--cell_type_column", default="cell_type", help="Cell type column in sc_metadata.csv.")
    p.add_argument("--max_genes", type=int, default=500, help="Maximum variable common genes used.")
    p.add_argument("--n_pcs", type=int, default=30, help="PCA dimensions for expression costs.")
    p.add_argument("--alpha", type=float, default=0.1, help="Weight for SpaOTsc structured GW term.")
    p.add_argument("--epsilon", type=float, default=0.1, help="Entropy regularization.")
    p.add_argument("--rho", default="inf", help="Unbalanced OT KL weight. Use 'inf' for balanced OT.")
    p.add_argument("--niter", type=int, default=10, help="SpaOTsc structured OT outer iterations.")
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


def _ensure_spaotsc_solver_available() -> Any:
    if os.environ.get("PYTHONNOUSERSITE") == "1":
        raise RuntimeError(
            "PYTHONNOUSERSITE=1 is set, but spaotsc was installed in the user site-packages. "
            "Unset PYTHONNOUSERSITE before running SpaOTsc."
        )
    try:
        from spaotsc.utils import usot
    except Exception as exc:
        raise RuntimeError("SpaOTsc solver module 'spaotsc.utils.usot' is not available.") from exc
    return usot


def _rho_value(raw: str) -> float:
    return float("inf") if str(raw).lower() in {"inf", "infinity", "np.inf"} else float(raw)


def _pca_joint(sc_x: pd.DataFrame, st_x: pd.DataFrame, n_pcs: int) -> tuple[np.ndarray, np.ndarray]:
    combined = pd.concat([sc_x, st_x], axis=0)
    arr = combined.apply(pd.to_numeric, errors="coerce").fillna(0.0).astype("float32").to_numpy()
    arr = StandardScaler(with_mean=True, with_std=True).fit_transform(arr)
    n_components = max(2, min(int(n_pcs), arr.shape[0] - 1, arr.shape[1]))
    z = PCA(n_components=n_components, random_state=42).fit_transform(arr).astype("float64")
    return z[: sc_x.shape[0]], z[sc_x.shape[0] :]


def _normalize_cost(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype="float64")
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
    x -= float(np.min(x))
    mx = float(np.max(x))
    if mx > 0:
        x /= mx
    return x


def _transport_plan(
    sc_z: np.ndarray,
    st_z: np.ndarray,
    coords: pd.DataFrame,
    alpha: float,
    epsilon: float,
    rho: float,
    niter: int,
) -> np.ndarray:
    usot = _ensure_spaotsc_solver_available()
    cost = _normalize_cost(pairwise_distances(sc_z, st_z, metric="euclidean"))
    weight = np.exp(1.0 - cost)
    mu = weight.sum(axis=1)
    nu = weight.sum(axis=0)
    mu = mu / np.sum(mu)
    nu = nu / np.sum(nu)

    if alpha <= 0:
        gamma = usot.uot(mu, nu, cost, epsilon, rho=rho)
    else:
        g_sc = _normalize_cost(pairwise_distances(sc_z, metric="euclidean"))
        g_is = _normalize_cost(pairwise_distances(coords[["coord_x", "coord_y"]].to_numpy(dtype="float64"), metric="euclidean"))
        gamma = usot.usot(mu, nu, cost, g_sc, g_is, float(alpha), epsilon=float(epsilon), rho=rho, niter=int(niter))
    return np.asarray(gamma, dtype="float64")


def main() -> int:
    args = parse_args()
    _ensure_spaotsc_solver_available()
    project_root = Path(args.project_root).resolve()
    out_dir = (
        Path(args.out_dir).resolve()
        if args.out_dir
        else project_root / "result" / args.sample / "stage4_mapping" / "spaotsc"
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
            sc_z, st_z = _pca_joint(sc_x, st_x, args.n_pcs)
            rho = _rho_value(args.rho)
            gamma = _transport_plan(sc_z, st_z, coords, args.alpha, args.epsilon, rho, args.niter)

            frac = _fraction_from_coupling(gamma, sc_meta.loc[sc_x.index, cell_type_column], st_x.index)
            cell_assignment = _cell_assignment_from_coupling(
                gamma,
                sc_x.index,
                sc_meta.loc[sc_x.index, cell_type_column],
                st_x.index,
            )
            _write_mapping_outputs(frac, cell_assignment, out_dir)
            (out_dir / "genes_used.txt").write_text("\n".join(genes) + "\n", encoding="utf-8")
            metrics = _eval_against_truth(project_root, args.group, args.sample, frac, sim_info)
            _json_write(out_dir / "metrics_simulation.json", metrics)
            _json_write(
                out_dir / "spaotsc_summary.json",
                {
                    "sample": args.sample,
                    "group": args.group,
                    "method": "spaotsc",
                    "status": "ok",
                    "processed_export_dir": str(export_dir),
                    "cell_type_column": cell_type_column,
                    "n_cells": int(sc_x.shape[0]),
                    "n_spots": int(st_x.shape[0]),
                    "n_genes_used": int(len(genes)),
                    "n_cell_types": int(sc_meta.loc[sc_x.index, cell_type_column].nunique()),
                    "alpha": float(args.alpha),
                    "epsilon": float(args.epsilon),
                    "rho": "inf" if np.isinf(rho) else float(rho),
                    "niter": int(args.niter),
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
            out_dir / "spaotsc_summary.json",
            {
                "sample": args.sample,
                "group": args.group,
                "method": "spaotsc",
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
                "method": "spaotsc",
                "status": status,
                "wall_time_seconds": float(elapsed),
                "peak_memory_mb": float(peak),
                "log_path": str(log_path),
            },
        )

    print(f"[OK] spaotsc {status}: {args.sample} -> {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
