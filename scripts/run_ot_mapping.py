#!/usr/bin/env python
from __future__ import annotations

import argparse
import contextlib
import importlib.util
import json
import os
import sys
import threading
import time
import traceback
import types
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import cdist

try:
    import psutil
except Exception:  # pragma: no cover
    psutil = None

_HERE = Path(__file__).resolve()
_ROOT = _HERE.parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from scripts.run_celltrek_mapping import _common_genes, _load_coordinates  # noqa: E402
from scripts.run_tangram_marker_mapping import (  # noqa: E402
    _eval_against_truth,
    _json_write,
    _load_inputs,
    _processed_export_dir,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Run novoSpaRc or SpaOTsc on one simulation sample."
    )
    p.add_argument("--method", required=True, choices=["novosparc", "spaotsc"])
    p.add_argument("--project_root", default=".")
    p.add_argument("--group", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--cell_type_column", default="cell_type")
    p.add_argument("--max_genes", type=int, default=500)
    p.add_argument("--n_pcs", type=int, default=30)
    p.add_argument("--num_neighbors_s", type=int, default=5)
    p.add_argument("--num_neighbors_t", type=int, default=5)
    p.add_argument("--alpha_linear", type=float, default=0.5)
    p.add_argument("--novosparc_epsilon", type=float, default=5e-3)
    p.add_argument("--spaotsc_alpha", type=float, default=0.1)
    p.add_argument("--spaotsc_epsilon", type=float, default=0.1)
    p.add_argument("--rho", default="inf")
    p.add_argument("--niter", type=int, default=10)
    p.add_argument("--out_dir", default=None)
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
            return 0.0
        self._stop.set()
        self._thread.join(timeout=2.0)
        try:
            self.peak_rss_mb = max(
                self.peak_rss_mb, self.process.memory_info().rss / (1024**2)
            )
        except Exception:
            pass
        return self.peak_rss_mb

    def _run(self) -> None:
        while not self._stop.is_set():
            try:
                self.peak_rss_mb = max(
                    self.peak_rss_mb, self.process.memory_info().rss / (1024**2)
                )
            except Exception:
                pass
            self._stop.wait(self.interval_seconds)


def _select_variable_genes(
    sc_expr: pd.DataFrame,
    st_expr: pd.DataFrame,
    genes: list[str],
    max_genes: int,
) -> list[str]:
    if len(genes) <= max_genes:
        return genes
    sc_var = sc_expr[genes].astype("float32").var(axis=0)
    st_var = st_expr[genes].astype("float32").var(axis=0)
    score = (sc_var.rank(pct=True) + st_var.rank(pct=True)).sort_values(
        ascending=False
    )
    return score.index[:max_genes].tolist()


def _pca(x: pd.DataFrame, n_pcs: int) -> np.ndarray:
    arr = x.astype("float32").to_numpy()
    arr = StandardScaler(with_mean=True, with_std=True).fit_transform(arr)
    n_components = max(2, min(n_pcs, arr.shape[0] - 1, arr.shape[1]))
    return PCA(n_components=n_components, random_state=42).fit_transform(arr)


def _joint_pca(
    sc_x: pd.DataFrame, st_x: pd.DataFrame, n_pcs: int
) -> tuple[np.ndarray, np.ndarray]:
    z = _pca(pd.concat([sc_x, st_x], axis=0), n_pcs).astype("float64")
    return z[: len(sc_x)], z[len(sc_x) :]


def _normalize_cost(x: np.ndarray) -> np.ndarray:
    out = np.nan_to_num(np.asarray(x, dtype="float64"))
    out -= float(out.min())
    maximum = float(out.max())
    if maximum > 0:
        out /= maximum
    return out


def _fraction_from_coupling(
    coupling: np.ndarray, cell_types: pd.Series, spot_ids: pd.Index
) -> pd.DataFrame:
    labels = cell_types.astype(str).to_numpy()
    type_order = sorted(pd.unique(labels).tolist())
    onehot = np.zeros((len(labels), len(type_order)), dtype="float64")
    lookup = {label: idx for idx, label in enumerate(type_order)}
    for row, label in enumerate(labels):
        onehot[row, lookup[label]] = 1.0
    frac = coupling.T @ onehot
    frac /= np.maximum(frac.sum(axis=1, keepdims=True), 1e-12)
    result = pd.DataFrame(
        frac, index=spot_ids.astype(str), columns=type_order
    )
    result.index.name = "spot_id"
    return result


def _cell_assignment(
    coupling: np.ndarray,
    cell_ids: pd.Index,
    cell_types: pd.Series,
    spot_ids: pd.Index,
) -> pd.DataFrame:
    best = np.argmax(coupling, axis=1)
    row_sum = np.maximum(coupling.sum(axis=1), 1e-12)
    return pd.DataFrame(
        {
            "cell_id": cell_ids.astype(str),
            "assigned_spot": spot_ids.astype(str).to_numpy()[best],
            "cell_type": cell_types.astype(str).to_numpy(),
            "assignment_score": coupling[np.arange(len(best)), best] / row_sum,
        }
    )


def _write_outputs(
    frac: pd.DataFrame, assignment: pd.DataFrame, out_dir: Path
) -> None:
    frac.to_csv(out_dir / "spot_type_fraction.csv")
    pd.DataFrame(
        {
            "spot_id": frac.index,
            "dominant_type": frac.idxmax(axis=1),
            "dominant_score": frac.max(axis=1),
        }
    ).to_csv(out_dir / "spot_type_dominant.csv", index=False)
    assignment.to_csv(out_dir / "cell_assignment.csv", index=False)


def _run_novosparc(
    sc_x: pd.DataFrame,
    st_x: pd.DataFrame,
    coords: pd.DataFrame,
    args: argparse.Namespace,
    out_dir: Path,
) -> tuple[np.ndarray, dict[str, Any]]:
    package_spec = importlib.util.find_spec("novosparc")
    if package_spec is None or package_spec.submodule_search_locations is None:
        raise RuntimeError("Package 'novosparc' is not available.")
    package_root = Path(next(iter(package_spec.submodule_search_locations)))

    reconstruction_spec = importlib.util.spec_from_file_location(
        "_svtuner_novosparc_reconstruction",
        package_root / "reconstruction" / "_reconstruction.py",
    )
    if reconstruction_spec is None or reconstruction_spec.loader is None:
        raise RuntimeError("Could not load novoSpaRc reconstruction core.")
    reconstruction = importlib.util.module_from_spec(reconstruction_spec)
    reconstruction_spec.loader.exec_module(reconstruction)

    fake_package = types.ModuleType("novosparc")
    fake_package.analysis = types.SimpleNamespace()
    previous_package = sys.modules.get("novosparc")
    sys.modules["novosparc"] = fake_package
    try:
        gw_spec = importlib.util.spec_from_file_location(
            "_svtuner_novosparc_gw",
            package_root / "reconstruction" / "_GWadjusted.py",
        )
        if gw_spec is None or gw_spec.loader is None:
            raise RuntimeError("Could not load novoSpaRc GW core.")
        gw_module = importlib.util.module_from_spec(gw_spec)
        gw_spec.loader.exec_module(gw_module)
    finally:
        if previous_package is None:
            sys.modules.pop("novosparc", None)
        else:
            sys.modules["novosparc"] = previous_package

    locations = coords[["coord_x", "coord_y"]].to_numpy(dtype="float64")
    expression_cost, location_cost = reconstruction.setup_for_OT_reconstruction(
        _pca(sc_x, args.n_pcs).astype("float32"),
        locations,
        num_neighbors_source=args.num_neighbors_s,
        num_neighbors_target=args.num_neighbors_t,
        verbose=True,
    )
    sc_values = sc_x.to_numpy(dtype="float64")
    st_values = st_x.to_numpy(dtype="float64")
    sc_values /= max(float(np.max(sc_values)), 1e-12)
    st_values /= max(float(np.max(st_values)), 1e-12)
    marker_cost = cdist(sc_values, st_values, metric="minkowski", p=2)
    p_expression = np.full(len(sc_x), 1.0 / len(sc_x), dtype="float64")
    p_locations = np.full(len(st_x), 1.0 / len(st_x), dtype="float64")
    coupling = gw_module.gromov_wasserstein_adjusted_norm(
        marker_cost,
        expression_cost,
        location_cost,
        args.alpha_linear,
        p_expression,
        p_locations,
        "square_loss",
        epsilon=args.novosparc_epsilon,
        verbose=True,
        random_ini=False,
    )
    return np.asarray(coupling, dtype="float64"), {
        "alpha_linear": args.alpha_linear,
        "epsilon": args.novosparc_epsilon,
        "num_neighbors_s": args.num_neighbors_s,
        "num_neighbors_t": args.num_neighbors_t,
        "import_mode": "official_reconstruction_core_without_top_level_init",
    }


def _run_spaotsc(
    sc_x: pd.DataFrame,
    st_x: pd.DataFrame,
    coords: pd.DataFrame,
    args: argparse.Namespace,
) -> tuple[np.ndarray, dict[str, Any]]:
    from spaotsc.utils import usot

    sc_z, st_z = _joint_pca(sc_x, st_x, args.n_pcs)
    cost = _normalize_cost(pairwise_distances(sc_z, st_z))
    weight = np.exp(1.0 - cost)
    mu = weight.sum(axis=1)
    nu = weight.sum(axis=0)
    mu /= mu.sum()
    nu /= nu.sum()
    rho = (
        float("inf")
        if str(args.rho).lower() in {"inf", "infinity", "np.inf"}
        else float(args.rho)
    )
    g_sc = _normalize_cost(pairwise_distances(sc_z))
    g_st = _normalize_cost(
        pairwise_distances(coords[["coord_x", "coord_y"]].to_numpy())
    )
    coupling = usot.usot(
        mu,
        nu,
        cost,
        g_sc,
        g_st,
        args.spaotsc_alpha,
        epsilon=args.spaotsc_epsilon,
        rho=rho,
        niter=args.niter,
    )
    return np.asarray(coupling, dtype="float64"), {
        "alpha": args.spaotsc_alpha,
        "epsilon": args.spaotsc_epsilon,
        "rho": "inf" if np.isinf(rho) else rho,
        "niter": args.niter,
    }


def main() -> int:
    args = parse_args()
    if os.environ.get("PYTHONNOUSERSITE") == "1":
        raise RuntimeError("Unset PYTHONNOUSERSITE to access mapping packages.")
    project_root = Path(args.project_root).absolute()
    out_dir = (
        Path(args.out_dir).resolve()
        if args.out_dir
        else project_root
        / "result"
        / args.sample
        / "stage4_mapping"
        / args.method
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "run.log"
    monitor = ResourceMonitor()
    monitor.start()
    started = time.perf_counter()
    status = "failed"
    try:
        with log_path.open("w", encoding="utf-8") as log_fh, contextlib.redirect_stdout(
            log_fh
        ), contextlib.redirect_stderr(log_fh):
            export_dir = _processed_export_dir(
                project_root, args.group, args.sample
            )
            sc_expr, st_expr, sc_meta, sim_info = _load_inputs(
                export_dir, args.cell_type_column
            )
            type_col = str(
                sc_meta.attrs.get("cell_type_column", args.cell_type_column)
            )
            coords = _load_coordinates(export_dir, st_expr.index)
            common_spots = st_expr.index.astype(str).intersection(coords.index)
            st_expr = st_expr.loc[common_spots]
            coords = coords.loc[common_spots]
            genes = _select_variable_genes(
                sc_expr,
                st_expr,
                _common_genes(sc_expr, st_expr),
                args.max_genes,
            )
            sc_x = sc_expr[genes].astype("float32")
            st_x = st_expr[genes].astype("float32")
            if args.method == "novosparc":
                coupling, method_params = _run_novosparc(
                    sc_x, st_x, coords, args, out_dir
                )
            else:
                coupling, method_params = _run_spaotsc(
                    sc_x, st_x, coords, args
                )
            types = sc_meta.loc[sc_x.index, type_col]
            frac = _fraction_from_coupling(coupling, types, st_x.index)
            assignment = _cell_assignment(
                coupling, sc_x.index, types, st_x.index
            )
            _write_outputs(frac, assignment, out_dir)
            (out_dir / "genes_used.txt").write_text(
                "\n".join(genes) + "\n", encoding="utf-8"
            )
            _json_write(
                out_dir / "metrics_simulation.json",
                _eval_against_truth(
                    project_root, args.group, args.sample, frac, sim_info
                ),
            )
            summary = {
                "sample": args.sample,
                "group": args.group,
                "method": args.method,
                "status": "ok",
                "processed_export_dir": str(export_dir),
                "cell_type_column": type_col,
                "n_cells": len(sc_x),
                "n_spots": len(st_x),
                "n_genes_used": len(genes),
                "n_cell_types": int(types.nunique()),
                **method_params,
            }
            _json_write(out_dir / f"{args.method}_summary.json", summary)
            status = "ok"
    except Exception as exc:
        _json_write(
            out_dir / f"{args.method}_summary.json",
            {
                "sample": args.sample,
                "group": args.group,
                "method": args.method,
                "status": "failed",
                "error_type": type(exc).__name__,
                "error": str(exc),
                "traceback": traceback.format_exc(),
            },
        )
        raise
    finally:
        peak_memory = monitor.stop()
        _json_write(
            out_dir / "resource_metrics.json",
            {
                "sample": args.sample,
                "group": args.group,
                "method": args.method,
                "status": status,
                "wall_time_seconds": time.perf_counter() - started,
                "peak_memory_mb": peak_memory,
                "log_path": str(log_path),
            },
        )
    print(f"[OK] {args.method}: {args.sample} -> {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
