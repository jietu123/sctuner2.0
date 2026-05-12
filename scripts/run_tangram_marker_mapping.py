#!/usr/bin/env python
from __future__ import annotations

import argparse
import contextlib
import json
import math
import os
import sys
import threading
import time
import traceback
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import psutil

_HERE = Path(__file__).resolve()
_ROOT = _HERE.parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src.utils.sample_paths import resolve_sample_dir


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Run Tangram(marker/all genes) on one simulation sample and export method-agnostic mapping outputs."
    )
    p.add_argument("--project_root", default=".", help="Project root.")
    p.add_argument("--group", default="real_brca", help="Simulation group under data/sim and processed simulation_experiments.")
    p.add_argument("--sample", required=True, help="Simulation sample id.")
    p.add_argument("--cell_type_column", default="cell_type", help="Cell type column in sc_metadata.csv.")
    p.add_argument("--gene_mode", choices=["marker", "all"], default="marker", help="Gene set used by Tangram.")
    p.add_argument("--top_n_marker", type=int, default=50, help="Top marker genes per cell type.")
    p.add_argument("--min_marker_genes", type=int, default=100, help="Fallback to HVGs until at least this many genes.")
    p.add_argument("--max_marker_genes", type=int, default=1000, help="Maximum marker genes used for Tangram.")
    p.add_argument("--num_epochs", type=int, default=500, help="Tangram optimization epochs.")
    p.add_argument("--device", default="cpu", help="Tangram device, e.g. cpu or cuda:0.")
    p.add_argument("--out_dir", default=None, help="Override output directory.")
    p.add_argument(
        "--prepare_only",
        action="store_true",
        help="Prepare marker genes and input summary without importing/running Tangram.",
    )
    return p.parse_args()


class ResourceMonitor:
    def __init__(self, interval_seconds: float = 0.5) -> None:
        self.interval_seconds = interval_seconds
        self.process = psutil.Process()
        self.peak_rss_mb = 0.0
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)

    def start(self) -> None:
        self._thread.start()

    def stop(self) -> float:
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
                rss = self.process.memory_info().rss / (1024**2)
                self.peak_rss_mb = max(self.peak_rss_mb, rss)
            except Exception:
                pass
            self._stop.wait(self.interval_seconds)


def _json_write(path: Path, obj: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _processed_export_dir(project_root: Path, group: str, sample: str) -> Path:
    candidates = [
        project_root / "data" / "processed" / "simulation_experiments" / group / sample / "stage1_preprocess" / "exported",
        project_root / "data" / "processed" / sample / "stage1_preprocess" / "exported",
    ]
    for cand in candidates:
        if cand.exists():
            return cand
    raise FileNotFoundError("Stage1 exported directory not found. Tried: " + ", ".join(str(x) for x in candidates))


def _ensure_tangram_available() -> None:
    if os.environ.get("PYTHONNOUSERSITE") == "1":
        raise RuntimeError(
            "PYTHONNOUSERSITE=1 is set, but tangram-sc was installed in the user site-packages. "
            "Unset PYTHONNOUSERSITE before running Tangram."
        )
    try:
        import tangram  # noqa: F401
    except Exception as exc:
        raise RuntimeError(
            "Tangram is not installed in the active environment. Install package 'tangram-sc' "
            "before running this script without --prepare_only."
        ) from exc


def _load_inputs(export_dir: Path, cell_type_column: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    sc_expr_path = export_dir / "sc_expression_normalized.csv"
    st_expr_path = export_dir / "st_expression_normalized.csv"
    sc_meta_path = export_dir / "sc_metadata.csv"
    sim_info_path = export_dir / "sim_info.json"
    for p in [sc_expr_path, st_expr_path, sc_meta_path, sim_info_path]:
        if not p.exists():
            raise FileNotFoundError(f"required input missing: {p}")

    sc_expr = pd.read_csv(sc_expr_path, index_col=0)
    st_expr = pd.read_csv(st_expr_path, index_col=0)
    sc_meta = pd.read_csv(sc_meta_path)
    sim_info = _read_json(sim_info_path)

    if cell_type_column not in sc_meta.columns:
        candidates = [c for c in ["cell_type", "sc_meta", "orig_type", "type", "label"] if c in sc_meta.columns]
        if not candidates:
            raise ValueError(f"cell type column '{cell_type_column}' not found and no fallback column is available.")
        cell_type_column = candidates[0]

    cell_id_cols = [c for c in sc_meta.columns if c == "cell_id" or c.startswith("cell_id.")]
    if not cell_id_cols:
        raise ValueError("sc_metadata.csv must contain a cell_id column.")
    sc_meta = sc_meta.copy()
    sc_meta["__cell_id__"] = sc_meta[cell_id_cols[0]].astype(str)
    sc_meta[cell_type_column] = sc_meta[cell_type_column].astype(str)
    sc_meta = sc_meta.set_index("__cell_id__", drop=False)
    sc_meta = sc_meta.loc[sc_meta.index.intersection(sc_expr.index)]
    sc_expr = sc_expr.loc[sc_meta.index]
    if sc_expr.empty or st_expr.empty:
        raise ValueError("empty sc or ST expression matrix after alignment.")
    sc_meta.attrs["cell_type_column"] = cell_type_column
    return sc_expr, st_expr, sc_meta, sim_info


def _select_marker_genes(
    sc_expr: pd.DataFrame,
    st_expr: pd.DataFrame,
    sc_meta: pd.DataFrame,
    export_dir: Path,
    cell_type_column: str,
    top_n: int,
    min_genes: int,
    max_genes: int,
) -> tuple[list[str], dict[str, Any]]:
    common_genes = [g for g in sc_expr.columns if g in set(st_expr.columns)]
    if not common_genes:
        raise ValueError("no common genes between sc and ST expression matrices.")

    x = sc_expr[common_genes].apply(pd.to_numeric, errors="coerce").fillna(0.0).astype("float32")
    labels = sc_meta.loc[x.index, cell_type_column].astype(str)
    type_counts = labels.value_counts().to_dict()
    total_sum = x.sum(axis=0)
    n_total = float(x.shape[0])

    selected: dict[str, float] = {}
    per_type_counts: dict[str, int] = {}
    for cell_type, n_cells in sorted(type_counts.items()):
        mask = labels == cell_type
        if int(n_cells) <= 0:
            continue
        type_mean = x.loc[mask].mean(axis=0)
        if n_total - float(n_cells) > 0:
            rest_mean = (total_sum - x.loc[mask].sum(axis=0)) / (n_total - float(n_cells))
        else:
            rest_mean = pd.Series(0.0, index=x.columns)
        score = (type_mean - rest_mean).sort_values(ascending=False)
        score = score[score > 0]
        top = score.head(max(1, top_n))
        per_type_counts[cell_type] = int(top.shape[0])
        for gene, val in top.items():
            selected[str(gene)] = max(float(val), selected.get(str(gene), -math.inf))

    hvg_used = 0
    hvg_path = export_dir.parent / "hvg_genes.txt"
    if len(selected) < min_genes and hvg_path.exists():
        hvg_genes = [g.strip() for g in hvg_path.read_text(encoding="utf-8").splitlines() if g.strip()]
        for gene in hvg_genes:
            if gene in common_genes and gene not in selected:
                selected[gene] = 0.0
                hvg_used += 1
                if len(selected) >= min_genes:
                    break

    if not selected:
        raise ValueError("failed to select marker genes.")

    ranked = sorted(selected.items(), key=lambda kv: kv[1], reverse=True)
    marker_genes = [g for g, _ in ranked[:max_genes]]
    stats = {
        "n_common_genes": int(len(common_genes)),
        "n_marker_genes": int(len(marker_genes)),
        "top_n_marker": int(top_n),
        "min_marker_genes": int(min_genes),
        "max_marker_genes": int(max_genes),
        "hvg_fallback_genes_used": int(hvg_used),
        "cell_type_counts": {str(k): int(v) for k, v in type_counts.items()},
        "marker_genes_per_type_positive": per_type_counts,
    }
    return marker_genes, stats


def _select_all_genes(sc_expr: pd.DataFrame, st_expr: pd.DataFrame) -> tuple[list[str], dict[str, Any]]:
    st_gene_set = set(st_expr.columns)
    genes = [g for g in sc_expr.columns if g in st_gene_set]
    if not genes:
        raise ValueError("no common genes between sc and ST expression matrices.")
    return genes, {
        "n_common_genes": int(len(genes)),
        "n_genes_used": int(len(genes)),
        "gene_mode": "all",
    }


def _as_dense_mapping(ad_map_x: Any, spot_ids: list[str]) -> np.ndarray:
    if hasattr(ad_map_x, "toarray"):
        arr = ad_map_x.toarray()
    else:
        arr = np.asarray(ad_map_x)
    if arr.shape[1] != len(spot_ids):
        raise ValueError(f"unexpected Tangram map shape {arr.shape}; expected second dimension {len(spot_ids)} spots.")
    return arr


def _matrix_to_type_fraction(arr: np.ndarray, cell_types: pd.Series, spot_ids: list[str]) -> pd.DataFrame:
    types = list(pd.Series(cell_types.astype(str)).drop_duplicates())
    out = {}
    for t in types:
        idx = np.where(cell_types.to_numpy(dtype=str) == t)[0]
        out[t] = arr[idx, :].sum(axis=0) if idx.size else np.zeros(arr.shape[1], dtype="float32")
    frac = pd.DataFrame(out, index=spot_ids)
    row_sum = frac.sum(axis=1).replace(0.0, np.nan)
    frac = frac.div(row_sum, axis=0).fillna(0.0)
    frac.index.name = "spot_id"
    return frac


def _matrix_to_cell_assignment(
    arr: np.ndarray,
    cell_ids: list[str],
    cell_types: pd.Series,
    spot_ids: list[str],
) -> pd.DataFrame:
    best_idx = np.argmax(arr, axis=1)
    best_score = arr[np.arange(arr.shape[0]), best_idx]
    return pd.DataFrame(
        {
            "cell_id": cell_ids,
            "assigned_spot": [spot_ids[i] for i in best_idx],
            "cell_type": cell_types.astype(str).to_numpy(),
            "assignment_score": best_score,
        }
    )


def _write_mapping_outputs(frac: pd.DataFrame, out_dir: Path, cell_assignment: pd.DataFrame | None = None) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    frac.to_csv(out_dir / "spot_type_fraction.csv")
    dominant = pd.DataFrame(
        {
            "spot_id": frac.index.astype(str),
            "dominant_type": frac.idxmax(axis=1).to_numpy(),
            "dominant_score": frac.max(axis=1).to_numpy(),
        }
    )
    no_signal = frac.sum(axis=1) <= 0
    dominant.loc[no_signal.to_numpy(), "dominant_type"] = "__NoType__"
    dominant.to_csv(out_dir / "spot_type_dominant.csv", index=False)
    if cell_assignment is not None:
        cell_assignment.to_csv(out_dir / "cell_assignment.csv", index=False)


def _eval_against_truth(project_root: Path, group: str, sample: str, pred_frac: pd.DataFrame, sim_info: dict[str, Any]) -> dict[str, Any]:
    try:
        raw_dir = resolve_sample_dir(project_root, sample, sim_group=group, must_exist=True)
    except Exception:
        return {"status": "skipped", "reason": "raw simulation directory not found"}
    truth_path = raw_dir / "sim_truth_spot_type_fraction.csv"
    if not truth_path.exists():
        return {"status": "skipped", "reason": f"truth file not found: {truth_path}"}
    truth = pd.read_csv(truth_path, index_col=0)
    common_spots = pred_frac.index.intersection(truth.index)
    if common_spots.empty:
        return {"status": "skipped", "reason": "no common spots with truth"}
    type_cols = sorted(set(pred_frac.columns).union(set(truth.columns)))
    pred = pred_frac.reindex(index=common_spots, columns=type_cols, fill_value=0.0).astype(float)
    obs = truth.reindex(index=common_spots, columns=type_cols, fill_value=0.0).astype(float)
    pred = pred.div(pred.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    obs = obs.div(obs.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)

    eps = 1e-12
    p = pred.to_numpy() + eps
    q = obs.to_numpy() + eps
    p = p / p.sum(axis=1, keepdims=True)
    q = q / q.sum(axis=1, keepdims=True)
    m = 0.5 * (p + q)
    js = 0.5 * np.sum(p * np.log(p / m), axis=1) + 0.5 * np.sum(q * np.log(q / m), axis=1)
    l1 = np.abs(p - q).sum(axis=1)

    pred_dom = pred.idxmax(axis=1)
    truth_dom = obs.idxmax(axis=1)
    missing_types = []
    mts = sim_info.get("missing_types")
    if isinstance(mts, list):
        missing_types.extend([str(x) for x in mts if str(x).strip()])
    mt = sim_info.get("missing_type")
    if mt is not None and str(mt).strip() and str(mt) not in missing_types:
        missing_types.append(str(mt))
    missing_cols = [t for t in missing_types if t in pred.columns]
    leakage_by_type = {t: float(pred[t].sum()) for t in missing_cols}

    per_type_corr = {}
    for t in type_cols:
        a = pred[t].to_numpy()
        b = obs[t].to_numpy()
        if np.std(a) <= 0 or np.std(b) <= 0:
            per_type_corr[t] = None
        else:
            per_type_corr[t] = float(np.corrcoef(a, b)[0, 1])

    return {
        "status": "ok",
        "n_common_spots": int(len(common_spots)),
        "mean_js_divergence": float(np.mean(js)),
        "mean_l1_distance": float(np.mean(l1)),
        "dominant_type_accuracy": float((pred_dom == truth_dom).mean()),
        "missing_types": missing_types,
        "missing_type_total_leakage": float(sum(leakage_by_type.values())),
        "missing_type_leakage_by_type": leakage_by_type,
        "per_type_pearson": per_type_corr,
    }


def _run_tangram(
    sc_expr: pd.DataFrame,
    st_expr: pd.DataFrame,
    sc_meta: pd.DataFrame,
    marker_genes: list[str],
    cell_type_column: str,
    num_epochs: int,
    device: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    try:
        import anndata as ad
        import tangram as tg
    except Exception as exc:
        raise RuntimeError(
            "Tangram is not installed in the active environment. Install package 'tangram-sc' "
            "before running this script without --prepare_only."
        ) from exc

    genes = [g for g in marker_genes if g in sc_expr.columns and g in st_expr.columns]
    if not genes:
        raise ValueError("no marker genes remain after intersecting sc/ST matrices.")
    adata_sc = ad.AnnData(sc_expr[genes].astype("float32"))
    adata_sp = ad.AnnData(st_expr[genes].astype("float32"))
    adata_sc.obs[cell_type_column] = sc_meta.loc[adata_sc.obs_names, cell_type_column].astype(str).to_numpy()

    tg.pp_adatas(adata_sc, adata_sp, genes=genes)
    ad_map = tg.map_cells_to_space(
        adata_sc,
        adata_sp,
        mode="cells",
        device=device,
        num_epochs=num_epochs,
    )
    spot_ids = list(adata_sp.obs_names.astype(str))
    arr = _as_dense_mapping(ad_map.X, spot_ids)
    frac = _matrix_to_type_fraction(arr, adata_sc.obs[cell_type_column], spot_ids)
    cell_assignment = _matrix_to_cell_assignment(
        arr,
        list(adata_sc.obs_names.astype(str)),
        adata_sc.obs[cell_type_column],
        spot_ids,
    )
    return frac, cell_assignment


def main() -> int:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    os.environ.setdefault("NUMBA_CACHE_DIR", str(project_root / ".numba_cache"))
    Path(os.environ["NUMBA_CACHE_DIR"]).mkdir(parents=True, exist_ok=True)
    method_name = "tangram_marker" if args.gene_mode == "marker" else "tangram_all"
    out_dir = (
        Path(args.out_dir).resolve()
        if args.out_dir
        else project_root / "result" / args.sample / "stage4_mapping" / method_name
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "run.log"
    resource_path = out_dir / "resource_metrics.json"
    start = time.perf_counter()
    monitor = ResourceMonitor()
    monitor.start()
    status = "failed"

    try:
        with log_path.open("w", encoding="utf-8") as log_fh, contextlib.redirect_stdout(log_fh), contextlib.redirect_stderr(log_fh):
            if not args.prepare_only:
                _ensure_tangram_available()
            export_dir = _processed_export_dir(project_root, args.group, args.sample)
            sc_expr, st_expr, sc_meta, sim_info = _load_inputs(export_dir, args.cell_type_column)
            cell_type_column = str(sc_meta.attrs.get("cell_type_column", args.cell_type_column))
            if args.gene_mode == "marker":
                genes_used, gene_stats = _select_marker_genes(
                    sc_expr,
                    st_expr,
                    sc_meta,
                    export_dir,
                    cell_type_column,
                    args.top_n_marker,
                    args.min_marker_genes,
                    args.max_marker_genes,
                )
                (out_dir / "marker_genes.txt").write_text("\n".join(genes_used) + "\n", encoding="utf-8")
            else:
                genes_used, gene_stats = _select_all_genes(sc_expr, st_expr)
            (out_dir / "genes_used.txt").write_text("\n".join(genes_used) + "\n", encoding="utf-8")

            summary = {
                "sample": args.sample,
                "group": args.group,
                "method": method_name,
                "gene_mode": args.gene_mode,
                "status": "prepared_only" if args.prepare_only else "running",
                "processed_export_dir": str(export_dir),
                "cell_type_column": cell_type_column,
                "n_cells": int(sc_expr.shape[0]),
                "n_spots": int(st_expr.shape[0]),
                "n_sc_genes": int(sc_expr.shape[1]),
                "n_st_genes": int(st_expr.shape[1]),
                "num_epochs": int(args.num_epochs),
                "device": args.device,
                "sim_info_missing_type": sim_info.get("missing_type"),
                "sim_info_missing_types": sim_info.get("missing_types"),
                "gene_stats": gene_stats,
            }
            _json_write(out_dir / "tangram_summary.json", summary)

            if args.prepare_only:
                status = "prepared_only"
            else:
                pred_frac, cell_assignment = _run_tangram(
                    sc_expr,
                    st_expr,
                    sc_meta,
                    genes_used,
                    cell_type_column,
                    args.num_epochs,
                    args.device,
                )
                _write_mapping_outputs(pred_frac, out_dir, cell_assignment)
                metrics = _eval_against_truth(project_root, args.group, args.sample, pred_frac, sim_info)
                _json_write(out_dir / "metrics_simulation.json", metrics)
                summary["status"] = "ok"
                summary["output_files"] = {
                    "spot_type_fraction": str(out_dir / "spot_type_fraction.csv"),
                    "spot_type_dominant": str(out_dir / "spot_type_dominant.csv"),
                    "cell_assignment": str(out_dir / "cell_assignment.csv"),
                    "metrics_simulation": str(out_dir / "metrics_simulation.json"),
                }
                _json_write(out_dir / "tangram_summary.json", summary)
                status = "ok"
    except Exception as exc:
        error = {
            "status": "failed",
            "sample": args.sample,
            "group": args.group,
            "method": method_name,
            "gene_mode": args.gene_mode,
            "error_type": type(exc).__name__,
            "error": str(exc),
            "traceback": traceback.format_exc(),
        }
        _json_write(out_dir / "tangram_summary.json", error)
        raise
    finally:
        peak = monitor.stop()
        elapsed = time.perf_counter() - start
        _json_write(
            resource_path,
            {
                "sample": args.sample,
                "group": args.group,
                "method": method_name,
                "gene_mode": args.gene_mode,
                "status": status,
                "wall_time_seconds": float(elapsed),
                "peak_memory_mb": float(peak),
                "log_path": str(log_path),
            },
        )

    print(f"[OK] {method_name} {status}: {args.sample} -> {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
