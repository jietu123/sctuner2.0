"""
Stage 5: Simulated-scenario evaluation (with ground truth from SimGen).

This runner is intentionally lightweight and focuses on "contract + guardrail" metrics:
- Use `unique_cid` (cell instance id) as the unique evaluation unit when available.
- Exclude Query cells with `true_spot_id = NA` from cell-level accuracy denominator.
- Report duplicate-cell-id rate and NA-truth fractions for traceability.

Outputs (per SVTuner.md conventions):
- result/<scenario_id>/stage5_eval/generic/metrics_<backend>.json
- result/<scenario_id>/stage5_eval/generic/metrics_<backend>_runs.csv
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "src"
for p in (PROJECT_ROOT, SRC_ROOT):
    if str(p) not in sys.path:
        sys.path.append(str(p))


def detect_project_root() -> Path:
    here = Path(__file__).resolve()
    try:
        return here.parents[2]
    except IndexError:
        return here.parent.parent


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _read_json(path: Path) -> Dict[str, Any]:
    if not path.exists() or path.stat().st_size == 0:
        return {}
    with path.open("r", encoding="utf-8") as f:
        return json.load(f) or {}


def _read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(str(path))
    return pd.read_csv(path)


def _as_str_series(s: pd.Series) -> pd.Series:
    # Use pandas "string" dtype to preserve missing values as <NA> (avoid "nan" strings).
    return s.astype("string")


def load_query_truth(sim_dir: Path) -> Tuple[pd.DataFrame, str]:
    candidates = [
        ("sim_truth_query_cell_spot.csv", "canonical"),
        ("cell_true_spot.csv", "alias"),
    ]
    for fname, _kind in candidates:
        p = sim_dir / fname
        if not p.exists():
            continue
        df = _read_csv(p)
        required = {"cell_id", "true_spot_id"}
        missing = required - set(df.columns)
        if missing:
            raise KeyError(f"Query truth missing columns {sorted(list(missing))}: {p}")
        df = df.copy()
        df["cell_id"] = _as_str_series(df["cell_id"])
        if "truth_status" in df.columns:
            df["truth_status"] = _as_str_series(df["truth_status"])
        if "celltype" in df.columns:
            df["celltype"] = _as_str_series(df["celltype"])
        return df, fname
    raise FileNotFoundError(f"Query truth not found under: {sim_dir} (expected sim_truth_query_cell_spot.csv or cell_true_spot.csv)")


@dataclass(frozen=True)
class Stage4Run:
    backend: str
    run_dir: Path
    run_id: str
    mode: str


def iter_stage4_runs(stage4_root: Path, backend: str) -> Iterable[Stage4Run]:
    bdir = stage4_root / backend
    if not bdir.exists():
        return []
    out: List[Stage4Run] = []
    for d in sorted([p for p in bdir.iterdir() if p.is_dir()]):
        meta = _read_json(d / "meta.json")
        run_id = str(meta.get("run_id") or d.name)
        mode = str(meta.get("mode") or ("baseline" if "baseline" in d.name else "plus"))
        out.append(Stage4Run(backend=backend, run_dir=d, run_id=run_id, mode=mode))
    return out


def _resolve_cell_assignment_path(run: Stage4Run) -> Path:
    expected = run.run_dir / f"cell_assignment_{run.mode}.csv"
    if expected.exists():
        return expected
    cand = sorted(list(run.run_dir.glob("cell_assignment_*.csv")))
    if cand:
        return cand[0]
    raise FileNotFoundError(f"No cell_assignment_*.csv found under {run.run_dir}")


def evaluate_cell_level(
    *,
    stage4_assignment: pd.DataFrame,
    truth_query: pd.DataFrame,
    strict_cell_id_alignment: bool,
) -> Tuple[Dict[str, Any], pd.DataFrame]:
    required = {"cell_id", "spot_id"}
    missing = required - set(stage4_assignment.columns)
    if missing:
        raise KeyError(f"Stage4 cell_assignment missing columns {sorted(list(missing))}")

    df = stage4_assignment.copy()
    df["cell_id"] = _as_str_series(df["cell_id"])
    df["spot_id"] = _as_str_series(df["spot_id"])

    eval_id_col = "unique_cid"
    eval_id_source = "stage4"
    if eval_id_col in df.columns:
        df[eval_id_col] = _as_str_series(df[eval_id_col])
        if df[eval_id_col].isna().any():
            raise ValueError("Stage4 unique_cid contains NA")
        dup = df[eval_id_col].duplicated().sum()
        if dup:
            raise ValueError(f"Stage4 unique_cid is not unique (duplicated={int(dup)})")
    else:
        eval_id_source = "generated_row_id"
        df = df.reset_index(drop=True)
        df[eval_id_col] = [f"row_{i:06d}" for i in range(len(df))]

    truth = truth_query.copy()
    truth_ids = set(truth["cell_id"].astype(str).tolist())
    stage4_ids = set(df["cell_id"].astype(str).tolist())
    not_in_truth = sorted(list(stage4_ids - truth_ids))
    if strict_cell_id_alignment and not_in_truth:
        raise ValueError(f"Stage4 cell_id not in Query truth (n={len(not_in_truth)}), examples={not_in_truth[:5]}")

    truth_missing_in_stage4 = sorted(list(truth_ids - stage4_ids))

    n_rows = int(len(df))
    n_unique_cell_id = int(df["cell_id"].nunique())
    duplicate_cell_id_rate = float(1.0 - (n_unique_cell_id / n_rows)) if n_rows > 0 else 0.0

    truth_key_cols = ["cell_id", "true_spot_id"]
    for opt in ["truth_status", "celltype"]:
        if opt in truth.columns:
            truth_key_cols.append(opt)
    truth_small = truth[truth_key_cols].copy()
    truth_small["cell_id"] = _as_str_series(truth_small["cell_id"])

    merged = df.merge(truth_small, on="cell_id", how="left", validate="many_to_one", indicator=True)
    merged["true_spot_id"] = _as_str_series(merged["true_spot_id"])

    n_missing_truth_row = int((merged["_merge"] == "left_only").sum())
    if strict_cell_id_alignment and n_missing_truth_row:
        raise ValueError(f"Stage4->truth join has left_only rows: n={n_missing_truth_row}")
    merged = merged.drop(columns=["_merge"])

    evaluable_mask = ~merged["true_spot_id"].isna()
    n_evaluable = int(evaluable_mask.sum())
    n_na_truth = int((~evaluable_mask).sum())
    na_truth_fraction = float(n_na_truth / n_rows) if n_rows > 0 else 0.0
    evaluable_query_cell_fraction = float(n_evaluable / n_rows) if n_rows > 0 else 0.0

    if n_evaluable > 0:
        correct = (merged.loc[evaluable_mask, "spot_id"].astype(str) == merged.loc[evaluable_mask, "true_spot_id"].astype(str))
        cell_level_accuracy = float(correct.mean())
        n_correct = int(correct.sum())
    else:
        cell_level_accuracy = None
        n_correct = 0

    metrics = {
        "eval_unit": eval_id_col,
        "eval_unit_source": eval_id_source,
        "n_rows": n_rows,
        "n_unique_cell_id": n_unique_cell_id,
        "duplicate_cell_id_rate": duplicate_cell_id_rate,
        "cell_id_not_in_truth_count": int(len(not_in_truth)),
        "truth_cells_not_in_stage4_count": int(len(truth_missing_in_stage4)),
        "truth_cells_not_in_stage4_fraction": (float(len(truth_missing_in_stage4) / len(truth_ids)) if len(truth_ids) else 0.0),
        "n_evaluable": n_evaluable,
        "n_na_truth": n_na_truth,
        "na_truth_fraction": na_truth_fraction,
        "evaluable_query_cell_fraction": evaluable_query_cell_fraction,
        "n_correct": n_correct,
        "cell_level_accuracy": cell_level_accuracy,
    }
    return metrics, merged


def run_stage5_for_scenario(
    *,
    scenario_id: str,
    project_root: Path,
    backends: List[str],
    strict_cell_id_alignment: bool,
) -> None:
    sim_dir = project_root / "data" / "sim" / scenario_id
    stage4_root = project_root / "result" / scenario_id / "stage4_mapping"
    out_root = project_root / "result" / scenario_id / "stage5_eval"
    out_generic = out_root / "generic"
    _ensure_dir(out_generic)

    scenario_meta = _read_json(sim_dir / "scenario_meta.json")
    truth_query, truth_file_used = load_query_truth(sim_dir)

    truth_total = int(len(truth_query))
    truth_na = int(truth_query["true_spot_id"].isna().sum())
    truth_evaluable = truth_total - truth_na
    truth_summary = {
        "truth_file_used": truth_file_used,
        "truth_total": truth_total,
        "truth_na": truth_na,
        "truth_na_fraction": (float(truth_na / truth_total) if truth_total else 0.0),
        "truth_evaluable": truth_evaluable,
        "truth_evaluable_fraction": (float(truth_evaluable / truth_total) if truth_total else 0.0),
    }

    for backend in backends:
        runs = list(iter_stage4_runs(stage4_root, backend))
        if not runs:
            continue

        run_rows = []
        metrics_by_run: Dict[str, Any] = {}
        for run in runs:
            ca_path = _resolve_cell_assignment_path(run)
            ca = _read_csv(ca_path)
            m_cell, _merged = evaluate_cell_level(
                stage4_assignment=ca,
                truth_query=truth_query,
                strict_cell_id_alignment=strict_cell_id_alignment,
            )
            meta = _read_json(run.run_dir / "meta.json")
            row = {
                "scenario_id": scenario_id,
                "backend": backend,
                "run_id": run.run_id,
                "mode": run.mode,
                "cell_assignment_path": str(ca_path.as_posix()),
                "truth_file_used": truth_file_used,
                **m_cell,
            }
            if meta:
                row["stage4_meta_status"] = meta.get("status")
                row["cell_id_space"] = meta.get("cell_id_space")
                row["cell_instance_id_column"] = meta.get("cell_instance_id_column")
            run_rows.append(row)
            metrics_by_run[run.run_id] = row

        out_json = out_generic / f"metrics_{backend}.json"
        payload = {
            "scenario_id": scenario_id,
            "backend": backend,
            "sim_dir": str(sim_dir.as_posix()),
            "stage4_root": str(stage4_root.as_posix()),
            "out_dir": str(out_root.as_posix()),
            "scenario_meta": scenario_meta,
            "truth_summary": truth_summary,
            "runs": metrics_by_run,
            "strict_cell_id_alignment": bool(strict_cell_id_alignment),
        }
        out_json.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

        out_csv = out_generic / f"metrics_{backend}_runs.csv"
        pd.DataFrame(run_rows).sort_values(["backend", "run_id"]).to_csv(out_csv, index=False)

        print(f"[Stage5] Wrote: {out_json}")
        print(f"[Stage5] Wrote: {out_csv}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Stage5: evaluate simulated scenarios using SimGen truth")
    p.add_argument("--sample", required=True, help="scenario_id, e.g. S0_matched / M1_sc_missing_Bcell")
    p.add_argument("--project_root", default=None)
    p.add_argument("--backends", default="cytospace", help="comma-separated backends (default: cytospace)")
    p.add_argument("--non_strict_cell_id_alignment", action="store_true", default=False)
    return p.parse_args()


def main():
    args = parse_args()
    project_root = Path(args.project_root) if args.project_root else detect_project_root()
    backends = [b.strip() for b in str(args.backends).split(",") if b.strip()]
    strict = not bool(args.non_strict_cell_id_alignment)
    run_stage5_for_scenario(
        scenario_id=str(args.sample),
        project_root=project_root,
        backends=backends or ["cytospace"],
        strict_cell_id_alignment=strict,
    )


if __name__ == "__main__":
    main()
