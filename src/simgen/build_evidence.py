"""
SimGen evidence printer (stdout-only).

This script reads SimGen + Stage1-exported + (optional) Stage4 outputs and prints
an evidence report to stdout. It does NOT write any evidence files under `result/`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import pandas as pd


def _now_iso() -> str:
    return datetime.now().isoformat(timespec="seconds")


def _detect_project_root() -> Path:
    here = Path(__file__).resolve()
    try:
        return here.parents[2]
    except IndexError:
        return here.parent.parent


def _read_text_head(path: Path, *, n_lines: int = 20) -> str:
    lines: List[str] = []
    try:
        with path.open("r", encoding="utf-8", errors="replace") as f:
            for _ in range(int(n_lines)):
                ln = f.readline()
                if not ln:
                    break
                lines.append(ln.rstrip("\n\r"))
    except Exception as e:
        return f"[read_head_error] {e}"
    return "\n".join(lines) + ("\n" if lines else "")


def _read_lines(path: Path) -> List[str]:
    return [ln.strip() for ln in path.read_text(encoding="utf-8", errors="replace").splitlines() if ln.strip()]


def _sha1_file(path: Path) -> Optional[str]:
    try:
        h = hashlib.sha1()
        with path.open("rb") as f:
            for chunk in iter(lambda: f.read(1024 * 1024), b""):
                h.update(chunk)
        return h.hexdigest().lower()
    except Exception:
        return None


def _mtime_iso(path: Path) -> Optional[str]:
    try:
        ts = path.stat().st_mtime
        return datetime.fromtimestamp(ts).isoformat(timespec="seconds")
    except Exception:
        return None


def _relpath_str(path: Path, project_root: Path) -> str:
    try:
        return path.resolve().relative_to(project_root.resolve()).as_posix()
    except Exception:
        return str(path.resolve())


def fingerprint_file(path: Path, *, project_root: Path) -> Dict[str, Any]:
    out: Dict[str, Any] = {"path": _relpath_str(path, project_root), "exists": bool(path.exists())}
    if not path.exists() or not path.is_file():
        out.update({"sha1": None, "bytes": None, "mtime_iso": None})
        return out
    try:
        out["bytes"] = int(path.stat().st_size)
    except Exception:
        out["bytes"] = None
    out["mtime_iso"] = _mtime_iso(path)
    out["sha1"] = _sha1_file(path)
    return out


def _required_files_check(dir_path: Path, required: Sequence[str]) -> Tuple[bool, List[str]]:
    missing = [name for name in required if not (dir_path / name).exists()]
    return (len(missing) == 0), missing


def _csv_stats_series(values: pd.Series) -> Dict[str, float]:
    v = pd.to_numeric(values, errors="coerce").dropna().astype(float)
    if v.empty:
        return {"min": float("nan"), "p10": float("nan"), "median": float("nan"), "mean": float("nan"), "p90": float("nan"), "max": float("nan")}
    return {
        "min": float(v.min()),
        "p10": float(v.quantile(0.10)),
        "median": float(v.quantile(0.50)),
        "mean": float(v.mean()),
        "p90": float(v.quantile(0.90)),
        "max": float(v.max()),
    }


def _run_py_compile(py_files: Sequence[Path]) -> Dict[str, Any]:
    results: List[Dict[str, Any]] = []
    ok = True
    for f in py_files:
        cmd = [sys.executable, "-m", "py_compile", str(f)]
        p = subprocess.run(cmd, capture_output=True, text=True)
        rec = {
            "cmd": " ".join(cmd),
            "file": str(f),
            "returncode": int(p.returncode),
            "stdout": (p.stdout or "").strip(),
            "stderr": (p.stderr or "").strip(),
        }
        results.append(rec)
        if p.returncode != 0:
            ok = False
    return {"ok": ok, "results": results}


def _find_stage4_cell_assignment_files(stage4_dir: Path) -> List[Path]:
    if not stage4_dir.exists():
        return []
    out: List[Path] = []
    for p in stage4_dir.rglob("cell_assignment_*.csv"):
        if p.is_file():
            out.append(p)
    return sorted(out, key=lambda x: str(x).lower())


def _stage4_truth_alignment(
    *,
    cell_assignment_files: Sequence[Path],
    truth_query_path: Path,
) -> Dict[str, Any]:
    if not cell_assignment_files:
        return {"checked": False, "files": []}
    truth = pd.read_csv(truth_query_path)
    truth_ids = set(truth["cell_id"].astype(str))
    out_rows: List[Dict[str, Any]] = []
    max_not_in_truth = 0
    for f in cell_assignment_files:
        df = pd.read_csv(f)
        if "cell_id" not in df.columns:
            out_rows.append({"file": str(f), "error": "missing cell_id column"})
            max_not_in_truth = max(max_not_in_truth, 1)
            continue
        cids = df["cell_id"].astype(str)
        n_total = int(len(cids))
        n_unique = int(cids.nunique(dropna=False))
        duplicate_rate = float(1.0 - (n_unique / n_total)) if n_total > 0 else 0.0
        not_in_truth = sorted(set(cids) - truth_ids)
        n_not_in_truth = int(len(not_in_truth))
        max_not_in_truth = max(max_not_in_truth, n_not_in_truth)
        out_rows.append(
            {
                "file": str(f),
                "n_rows": n_total,
                "n_unique_cell_id": n_unique,
                "duplicate_cell_id_rate": duplicate_rate,
                "cell_id_not_in_truth": n_not_in_truth,
            }
        )
    return {"checked": True, "files": out_rows, "max_cell_id_not_in_truth": int(max_not_in_truth)}


def build_evidence_for_scenario(*, project_root: Path, scenario_id: str, sim_root: Path) -> Dict[str, Any]:
    sim_out_dir = project_root / sim_root / scenario_id
    stage1_export_dir = project_root / "data" / "processed" / scenario_id / "stage1_preprocess" / "exported"
    stage1_dir = stage1_export_dir.parent
    stage4_dir = project_root / "result" / scenario_id / "stage4_mapping"

    hard_errors: List[str] = []
    warnings: List[str] = []

    required_sim = [
        "sim_sc_expression.csv",
        "sim_sc_metadata.csv",
        "sim_st_expression.csv",
        "sim_st_coordinates.csv",
        "sim_truth_spot_type_fraction.csv",
        "spot_true_type_fraction.csv",
        "sim_truth_query_cell_spot.csv",
        "cell_true_spot.csv",
        "sim_truth_cell_spot.csv",
        "scenario_meta.json",
    ]
    required_stage1 = [
        "sc_expression_normalized.csv",
        "st_expression_normalized.csv",
        "sc_metadata.csv",
        "st_coordinates.csv",
    ]

    ok_sim, missing_sim = _required_files_check(sim_out_dir, required_sim)
    ok_stage1, missing_stage1 = _required_files_check(stage1_export_dir, required_stage1)
    if not sim_out_dir.exists():
        hard_errors.append(f"sim_out_dir not found: {sim_out_dir}")
    if missing_sim:
        hard_errors.append(f"missing sim outputs: {missing_sim}")
    if missing_stage1:
        hard_errors.append(f"missing stage1 exported outputs: {missing_stage1}")

    scenario_meta_path = sim_out_dir / "scenario_meta.json"
    scenario_meta: Dict[str, Any] = {}
    if scenario_meta_path.exists():
        try:
            scenario_meta = json.loads(scenario_meta_path.read_text(encoding="utf-8"))
        except Exception as e:
            hard_errors.append(f"failed to parse scenario_meta.json: {e}")

    # ---- py_compile
    py_files = [
        project_root / "src" / "simgen" / "make_scenario.py",
        project_root / "src" / "simgen" / "build_evidence.py",
        project_root / "src" / "stages" / "backends" / "cytospace_backend.py",
        project_root / "src" / "stages" / "stage2_svg_plugin.py",
        project_root / "src" / "stages" / "stage3_type_plugin.py",
        project_root / "src" / "stages" / "stage4_mapping.py",
    ]
    py_compile = _run_py_compile(py_files)
    if not py_compile.get("ok", False):
        hard_errors.append("py_compile failed")

    # ---- schema / alignment checks
    schema_checks: Dict[str, Any] = {
        "scenario_id": scenario_id,
        "generated_at": _now_iso(),
        "paths": {
            "sim_out_dir": str(sim_out_dir),
            "stage1_export_dir": str(stage1_export_dir),
        },
        "required_files": {
            "sim": {"ok": ok_sim, "missing": missing_sim, "required": required_sim},
            "stage1_exported": {"ok": ok_stage1, "missing": missing_stage1, "required": required_stage1},
        },
        "schema": {},
        "alignment": {},
        "uniqueness": {},
        "scenario_semantics": {},
        "stage4_alignment": {},
        "hard_errors": [],
        "warnings": [],
    }

    # Core reads (guarded)
    sim_sc_meta = None
    sim_truth_query = None
    sim_truth_spot = None
    sim_truth_world = None
    sim_st_expr = None
    sim_st_coords = None
    stage1_sc_meta = None
    stage1_st_expr = None
    stage1_st_coords = None

    try:
        sim_sc_meta = pd.read_csv(sim_out_dir / "sim_sc_metadata.csv")
        schema_checks["schema"]["sim_sc_metadata.csv"] = {"columns": list(sim_sc_meta.columns)}
        for col in ("cell_id", "celltype"):
            if col not in sim_sc_meta.columns:
                hard_errors.append(f"sim_sc_metadata.csv missing required column: {col}")
    except Exception as e:
        hard_errors.append(f"failed to read sim_sc_metadata.csv: {e}")

    try:
        sim_truth_query = pd.read_csv(sim_out_dir / "sim_truth_query_cell_spot.csv")
        schema_checks["schema"]["sim_truth_query_cell_spot.csv"] = {"columns": list(sim_truth_query.columns)}
        for col in ("cell_id", "true_spot_id"):
            if col not in sim_truth_query.columns:
                hard_errors.append(f"sim_truth_query_cell_spot.csv missing required column: {col}")
    except Exception as e:
        hard_errors.append(f"failed to read sim_truth_query_cell_spot.csv: {e}")

    try:
        sim_truth_spot = pd.read_csv(sim_out_dir / "sim_truth_spot_type_fraction.csv", index_col=0)
        schema_checks["schema"]["sim_truth_spot_type_fraction.csv"] = {
            "n_spots": int(sim_truth_spot.shape[0]),
            "n_types": int(sim_truth_spot.shape[1]),
        }
    except Exception as e:
        hard_errors.append(f"failed to read sim_truth_spot_type_fraction.csv: {e}")

    try:
        sim_truth_world = pd.read_csv(sim_out_dir / "sim_truth_cell_spot.csv")
        schema_checks["schema"]["sim_truth_cell_spot.csv"] = {"columns": list(sim_truth_world.columns)}
    except Exception as e:
        warnings.append(f"failed to read sim_truth_cell_spot.csv (audit-only): {e}")

    try:
        sim_st_expr = pd.read_csv(sim_out_dir / "sim_st_expression.csv", index_col=0)
        sim_st_coords = pd.read_csv(sim_out_dir / "sim_st_coordinates.csv", index_col=0)
        schema_checks["schema"]["sim_st_coordinates.csv"] = {"columns": list(sim_st_coords.columns)}
        if "in_tissue" not in sim_st_coords.columns:
            warnings.append("sim_st_coordinates.csv missing in_tissue column")
        if set(sim_st_expr.index.astype(str)) != set(sim_st_coords.index.astype(str)):
            hard_errors.append("sim_st_expression vs sim_st_coordinates spot_id set mismatch")
        if not sim_st_expr.index.astype(str).equals(sim_st_coords.index.astype(str)):
            hard_errors.append("sim_st_expression vs sim_st_coordinates spot_id order mismatch")
    except Exception as e:
        hard_errors.append(f"failed to read/align sim_st_expression/sim_st_coordinates: {e}")

    try:
        stage1_sc_meta = pd.read_csv(stage1_export_dir / "sc_metadata.csv")
        stage1_st_expr = pd.read_csv(stage1_export_dir / "st_expression_normalized.csv", index_col=0)
        stage1_st_coords = pd.read_csv(stage1_export_dir / "st_coordinates.csv", index_col=0)
        schema_checks["schema"]["stage1_exported/st_coordinates.csv"] = {"columns": list(stage1_st_coords.columns)}
        if "in_tissue" not in stage1_st_coords.columns:
            warnings.append("stage1 exported st_coordinates.csv missing in_tissue column")
        if set(stage1_st_expr.index.astype(str)) != set(stage1_st_coords.index.astype(str)):
            hard_errors.append("stage1 st_expression_normalized vs st_coordinates spot_id set mismatch")
        if not stage1_st_expr.index.astype(str).equals(stage1_st_coords.index.astype(str)):
            hard_errors.append("stage1 st_expression_normalized vs st_coordinates spot_id order mismatch")
        if "cell_id" not in stage1_sc_meta.columns or "celltype" not in stage1_sc_meta.columns:
            hard_errors.append("stage1 sc_metadata.csv must include cell_id and celltype")
    except Exception as e:
        hard_errors.append(f"failed to read/align stage1 exported: {e}")

    # Uniqueness / duplicate rates
    if sim_sc_meta is not None and "cell_id" in sim_sc_meta.columns:
        total = int(len(sim_sc_meta))
        uniq = int(sim_sc_meta["cell_id"].astype(str).nunique(dropna=False))
        dup_rate = float(1.0 - (uniq / total)) if total > 0 else 0.0
        schema_checks["uniqueness"]["sim_sc_metadata.cell_id"] = {
            "n_total": total,
            "n_unique": uniq,
            "duplicate_rate": dup_rate,
        }
    if stage1_sc_meta is not None and "cell_id" in stage1_sc_meta.columns:
        total = int(len(stage1_sc_meta))
        uniq = int(stage1_sc_meta["cell_id"].astype(str).nunique(dropna=False))
        dup_rate = float(1.0 - (uniq / total)) if total > 0 else 0.0
        schema_checks["uniqueness"]["stage1_sc_metadata.cell_id"] = {
            "n_total": total,
            "n_unique": uniq,
            "duplicate_rate": dup_rate,
        }

    # Query truth NA fraction
    if sim_truth_query is not None and "true_spot_id" in sim_truth_query.columns:
        na_mask = sim_truth_query["true_spot_id"].isna() | (sim_truth_query["true_spot_id"].astype(str).str.strip().str.lower().isin({"na", "nan", ""}))
        schema_checks["alignment"]["query_truth_na_fraction"] = float(na_mask.mean()) if len(na_mask) else 0.0

    # Coverage (from scenario_meta preferred; fallback to exported st_coordinates)
    cov = (scenario_meta.get("coverage_audit") or {}) if isinstance(scenario_meta, dict) else {}
    cov_summary: Dict[str, Any] = {}
    try:
        steps = list(cov.get("fallback_steps") or [])
        actions = [str(s.get("action")) for s in steps if isinstance(s, dict) and "action" in s]
        action_counts: Dict[str, int] = {}
        for a in actions:
            action_counts[a] = int(action_counts.get(a, 0) + 1)
        cov_summary = {
            "min_cells_per_spot": cov.get("min_cells_per_spot"),
            "max_empty_spot_fraction": cov.get("max_empty_spot_fraction"),
            "max_resample_rounds": cov.get("max_resample_rounds"),
            "max_radius_multiplier": cov.get("max_radius_multiplier"),
            "radius_final": cov.get("radius_final"),
            "empty_spot_fraction": cov.get("empty_spot_fraction"),
            "min_cells_per_spot_satisfied_fraction": cov.get("min_cells_per_spot_satisfied_fraction"),
            "min_cells_per_spot_satisfied_fraction_warn": cov.get("min_cells_per_spot_satisfied_fraction_warn"),
            "fallback_steps_n": int(len(steps)),
            "fallback_action_counts": action_counts,
            "fallback_actions_tail": actions[-10:],
        }
        if (cov_summary.get("empty_spot_fraction") is None) or (cov_summary.get("min_cells_per_spot_satisfied_fraction") is None):
            for s in reversed(list(cov.get("fallback_steps") or [])):
                if str(s.get("action")) == "coverage_check":
                    cov_summary["empty_spot_fraction"] = cov_summary.get("empty_spot_fraction") or s.get("empty_spot_fraction")
                    cov_summary["min_cells_per_spot_satisfied_fraction"] = cov_summary.get("min_cells_per_spot_satisfied_fraction") or s.get("min_cells_per_spot_satisfied_fraction")
                    break
        if stage1_st_coords is not None and "world_cells_in_spot" in stage1_st_coords.columns:
            world_counts = pd.to_numeric(stage1_st_coords["world_cells_in_spot"], errors="coerce").fillna(0.0)
            cov_summary.setdefault("empty_spot_fraction_from_exported", float((world_counts <= 0).mean()))
            try:
                m = cov_summary.get("min_cells_per_spot")
                m = int(m) if m is not None else None
            except Exception:
                m = None
            if m is not None:
                cov_summary.setdefault("min_cells_per_spot_satisfied_fraction_from_exported", float((world_counts >= float(m)).mean()))
        if stage1_st_coords is not None and "in_tissue" in stage1_st_coords.columns:
            it = pd.to_numeric(stage1_st_coords["in_tissue"], errors="coerce").fillna(0.0)
            cov_summary["in_tissue_fraction"] = float((it > 0.5).mean())
            cov_summary["n_in_tissue_spots"] = int((it > 0.5).sum())
    except Exception as e:
        warnings.append(f"failed to compute coverage summary: {e}")
    schema_checks["alignment"]["coverage_summary"] = cov_summary

    # Stage1 non-exported gene lists (minimal compat)
    try:
        common_path = stage1_dir / "common_genes.txt"
        hvg_path = stage1_dir / "hvg_genes.txt"
        manifest_path = stage1_dir / "genes_manifest.json"
        common_genes = _read_lines(common_path) if common_path.exists() else []
        hvg_genes = _read_lines(hvg_path) if hvg_path.exists() else []
        hvg_subset = bool(set(hvg_genes).issubset(set(common_genes))) if (common_genes or hvg_genes) else True
        manifest_summary: Dict[str, Any] = {}
        if manifest_path.exists():
            try:
                m = json.loads(manifest_path.read_text(encoding="utf-8"))
                manifest_summary = {
                    "hvg_source": m.get("hvg_source"),
                    "hvg_cap_n_top": m.get("hvg_cap_n_top"),
                    "hvg_cap_applied": m.get("hvg_cap_applied"),
                    "min_hvg_after_filter": m.get("min_hvg_after_filter"),
                }
            except Exception:
                manifest_summary = {"parse_error": True}
        schema_checks["alignment"]["stage1_gene_lists"] = {
            "common_genes_exists": bool(common_path.exists()),
            "hvg_genes_exists": bool(hvg_path.exists()),
            "genes_manifest_exists": bool(manifest_path.exists()),
            "n_common_genes": int(len(common_genes)),
            "n_hvg_genes": int(len(hvg_genes)),
            "hvg_subset_of_common": bool(hvg_subset),
            "manifest": manifest_summary,
        }
        if (common_path.exists() and hvg_path.exists()) and not hvg_subset:
            hard_errors.append("stage1 gene lists: hvg_genes.txt is not a subset of common_genes.txt")
    except Exception as e:
        warnings.append(f"failed to check stage1 gene lists: {e}")

    # Truth alias integrity (canonical vs aliases must be identical)
    try:
        alias_checks: Dict[str, Any] = {}
        q_can = sim_out_dir / "sim_truth_query_cell_spot.csv"
        q_alias = sim_out_dir / "cell_true_spot.csv"
        if q_can.exists() and q_alias.exists():
            alias_checks["cell_true_spot_matches_canonical"] = bool(_sha1_file(q_can) == _sha1_file(q_alias))
        s_can = sim_out_dir / "sim_truth_spot_type_fraction.csv"
        s_alias = sim_out_dir / "spot_true_type_fraction.csv"
        if s_can.exists() and s_alias.exists():
            alias_checks["spot_true_type_fraction_matches_canonical"] = bool(_sha1_file(s_can) == _sha1_file(s_alias))
        schema_checks["alignment"]["truth_alias_checks"] = alias_checks
        for k, v in alias_checks.items():
            if v is False:
                hard_errors.append(f"truth alias mismatch: {k}")
    except Exception as e:
        warnings.append(f"failed to check truth aliases: {e}")

    # Stage1-exported sc MUST be Query set (contract)
    try:
        if sim_sc_meta is not None and stage1_sc_meta is not None and "cell_id" in sim_sc_meta.columns and "cell_id" in stage1_sc_meta.columns:
            sim_set = set(sim_sc_meta["cell_id"].astype(str))
            stg_set = set(stage1_sc_meta["cell_id"].astype(str))
            schema_checks["alignment"]["stage1_sc_vs_sim_sc"] = {
                "n_sim": int(len(sim_set)),
                "n_stage1": int(len(stg_set)),
                "n_intersection": int(len(sim_set & stg_set)),
                "n_sim_only": int(len(sim_set - stg_set)),
                "n_stage1_only": int(len(stg_set - sim_set)),
            }
            if sim_set != stg_set:
                hard_errors.append("stage1 sc_metadata cell_id set != sim_sc_metadata cell_id set (stage1_sc must be Query)")
    except Exception as e:
        warnings.append(f"failed to check stage1_sc vs sim_sc: {e}")

    # Scenario semantics audit (S0/M1/M2)
    sc_missing = []
    st_missing = []
    try:
        cfg_eff = scenario_meta.get("config_effective_subset") or {}
        sc_missing = list(cfg_eff.get("sc_missing_types") or [])
        st_missing = list(cfg_eff.get("st_missing_types") or [])
    except Exception:
        sc_missing = []
        st_missing = []

    sem_rows: List[Dict[str, Any]] = []
    if sim_sc_meta is not None and sim_truth_query is not None and sim_truth_spot is not None:
        q_type_col = "celltype" if "celltype" in sim_sc_meta.columns else ("type" if "type" in sim_sc_meta.columns else None)
        if q_type_col is None:
            hard_errors.append("sim_sc_metadata.csv missing type column (celltype/type)")
        else:
            q_counts = sim_sc_meta[q_type_col].astype(str).value_counts()
            tq = sim_truth_query[["cell_id", "true_spot_id"]].copy()
            tq["cell_id"] = tq["cell_id"].astype(str)
            meta2 = sim_sc_meta[["cell_id", q_type_col]].copy()
            meta2["cell_id"] = meta2["cell_id"].astype(str)
            joined = meta2.merge(tq, on="cell_id", how="left", validate="m:1")
            joined["is_na_truth"] = joined["true_spot_id"].isna() | (
                joined["true_spot_id"].astype(str).str.strip().str.lower().isin({"na", "nan", ""})
            )
            na_by_type = joined.groupby(q_type_col)["is_na_truth"].mean()

            for t in sorted(set(list(q_counts.index.astype(str)) + list(sim_truth_spot.columns.astype(str)))):
                spot_sum = float(pd.to_numeric(sim_truth_spot.get(t, pd.Series(dtype=float)), errors="coerce").fillna(0.0).sum()) if t in sim_truth_spot.columns else 0.0
                spot_nonzero = int((pd.to_numeric(sim_truth_spot.get(t, 0.0), errors="coerce").fillna(0.0) > 0).sum()) if t in sim_truth_spot.columns else 0
                sem_rows.append(
                    {
                        "type": str(t),
                        "query_cell_count": int(q_counts.get(t, 0)),
                        "spot_truth_col_sum": spot_sum,
                        "spot_truth_nonzero_spots": spot_nonzero,
                        "query_truth_na_fraction": float(na_by_type.get(t, float("nan"))),
                        "is_sc_missing_type": bool(t in sc_missing),
                        "is_st_missing_type": bool(t in st_missing),
                    }
                )
    schema_checks["scenario_semantics"]["sc_missing_types"] = sc_missing
    schema_checks["scenario_semantics"]["st_missing_types"] = st_missing

    # Enforce semantics contracts as hard errors (minimal but useful)
    for t in sc_missing:
        row = next((r for r in sem_rows if r["type"] == str(t)), None)
        if row is None:
            warnings.append(f"sc_missing_type not found in type tables: {t}")
            continue
        if int(row.get("query_cell_count", 0)) != 0:
            hard_errors.append(f"sc_missing_type {t}: query_cell_count expected 0, got {row.get('query_cell_count')}")
        if float(row.get("spot_truth_col_sum", 0.0)) <= 0.0:
            hard_errors.append(f"sc_missing_type {t}: spot_truth_col_sum expected >0, got {row.get('spot_truth_col_sum')}")
    for t in st_missing:
        row = next((r for r in sem_rows if r["type"] == str(t)), None)
        if row is None:
            warnings.append(f"st_missing_type not found in type tables: {t}")
            continue
        if float(row.get("spot_truth_col_sum", 0.0)) != 0.0:
            hard_errors.append(f"st_missing_type {t}: spot_truth_col_sum expected 0, got {row.get('spot_truth_col_sum')}")
        # NA truth for missing-in-ST should be high; record as warning if not.
        na_frac = row.get("query_truth_na_fraction")
        if isinstance(na_frac, (int, float)) and pd.notna(na_frac) and float(na_frac) < 0.95:
            warnings.append(f"st_missing_type {t}: query_truth_na_fraction expected ~1.0, got {na_frac}")

    # Stage4 ↔ truth alignment
    stage4_alignment: Dict[str, Any] = {"checked": False, "files": []}
    truth_query_path = sim_out_dir / "sim_truth_query_cell_spot.csv"
    if truth_query_path.exists():
        cell_assignment_files = _find_stage4_cell_assignment_files(stage4_dir)
        stage4_alignment = _stage4_truth_alignment(cell_assignment_files=cell_assignment_files, truth_query_path=truth_query_path)
        if stage4_alignment.get("checked") and int(stage4_alignment.get("max_cell_id_not_in_truth", 0)) != 0:
            hard_errors.append(f"stage4 cell_id not in truth: max={stage4_alignment.get('max_cell_id_not_in_truth')}")
    schema_checks["stage4_alignment"] = stage4_alignment

    schema_checks["hard_errors"] = list(hard_errors)
    schema_checks["warnings"] = list(warnings)

    # ---- stats tables
    counts_summary_rows: List[Dict[str, Any]] = []
    per_type_rows: List[Dict[str, Any]] = list(sem_rows)
    cells_per_spot_rows: List[Dict[str, Any]] = []

    n_query = int(len(sim_sc_meta)) if sim_sc_meta is not None else None
    n_world = int(len(sim_truth_world)) if sim_truth_world is not None else None
    n_spots = int(sim_truth_spot.shape[0]) if sim_truth_spot is not None else None
    n_types = int(sim_truth_spot.shape[1]) if sim_truth_spot is not None else None
    query_na = schema_checks.get("alignment", {}).get("query_truth_na_fraction")
    cov_sum = schema_checks.get("alignment", {}).get("coverage_summary") or {}

    counts_summary_rows.append(
        {
            "scenario_id": scenario_id,
            "n_query_cells": n_query,
            "n_world_cells": n_world,
            "n_spots": n_spots,
            "n_types": n_types,
            "query_truth_na_fraction": query_na,
            "coverage_min_cells_per_spot": cov_sum.get("min_cells_per_spot"),
            "coverage_max_empty_spot_fraction": cov_sum.get("max_empty_spot_fraction"),
            "coverage_empty_spot_fraction": cov_sum.get("empty_spot_fraction"),
            "coverage_min_cells_satisfied_fraction": cov_sum.get("min_cells_per_spot_satisfied_fraction"),
            "coverage_min_cells_satisfied_warn": cov_sum.get("min_cells_per_spot_satisfied_fraction_warn"),
            "coverage_in_tissue_fraction": cov_sum.get("in_tissue_fraction"),
            "n_sc_missing_types": int(len(sc_missing)),
            "n_st_missing_types": int(len(st_missing)),
            "stage4_checked": bool(stage4_alignment.get("checked", False)),
            "stage4_max_cell_id_not_in_truth": stage4_alignment.get("max_cell_id_not_in_truth"),
        }
    )

    # cells per spot stats
    try:
        if sim_truth_world is not None and "true_spot_id" in sim_truth_world.columns and sim_truth_spot is not None:
            spot_ids = sim_truth_spot.index.astype(str)
            w_counts = sim_truth_world["true_spot_id"].astype(str).value_counts()
            w_per_spot = pd.Series([int(w_counts.get(sid, 0)) for sid in spot_ids], index=spot_ids)
            w_stats = _csv_stats_series(w_per_spot)
            w_empty_frac = float((w_per_spot == 0).mean()) if len(w_per_spot) else 0.0

            if sim_truth_query is not None and "true_spot_id" in sim_truth_query.columns:
                q_valid = sim_truth_query[~sim_truth_query["true_spot_id"].isna()].copy()
                q_valid = q_valid[q_valid["true_spot_id"].astype(str).str.strip() != ""]
                q_counts = q_valid["true_spot_id"].astype(str).value_counts()
                q_per_spot = pd.Series([int(q_counts.get(sid, 0)) for sid in spot_ids], index=spot_ids)
                q_stats = _csv_stats_series(q_per_spot)
                q_empty_frac = float((q_per_spot == 0).mean()) if len(q_per_spot) else 0.0
            else:
                q_stats = {"min": float("nan"), "p10": float("nan"), "median": float("nan"), "mean": float("nan"), "p90": float("nan"), "max": float("nan")}
                q_empty_frac = float("nan")

            cells_per_spot_rows.append({"group": "world_cells_per_spot", **w_stats, "empty_fraction": w_empty_frac})
            cells_per_spot_rows.append({"group": "query_cells_per_spot(non_na_truth)", **q_stats, "empty_fraction": q_empty_frac})
    except Exception as e:
        warnings.append(f"failed to compute cells_per_spot stats: {e}")

    # ---- file fingerprints
    paths_to_fingerprint: List[Path] = []
    paths_to_fingerprint.extend(py_files)
    paths_to_fingerprint.extend(
        [
            project_root / "dosc" / "SVTuner.md",
            project_root / "dosc" / "插件统一接口规范api.md",
            project_root / "configs" / "simgen" / f"{scenario_id}.yaml",
        ]
    )
    for name in required_sim:
        paths_to_fingerprint.append(sim_out_dir / name)
    for name in required_stage1:
        paths_to_fingerprint.append(stage1_export_dir / name)

    for f in _find_stage4_cell_assignment_files(stage4_dir):
        paths_to_fingerprint.append(f)
        meta_json = f.parent / "meta.json"
        if meta_json.exists():
            paths_to_fingerprint.append(meta_json)

    # Dedup while preserving order
    seen = set()
    uniq_paths: List[Path] = []
    for p in paths_to_fingerprint:
        key = str(p.resolve()).lower()
        if key in seen:
            continue
        seen.add(key)
        uniq_paths.append(p)

    fingerprints = [fingerprint_file(p, project_root=project_root) for p in uniq_paths]

    # ---- evidence summary
    status = "success" if not hard_errors else "failed"
    evidence = {
        "status": status,
        "scenario_id": scenario_id,
        "generated_at": _now_iso(),
        "paths": {
            "sim_out_dir": str(sim_out_dir),
            "stage1_export_dir": str(stage1_export_dir),
            "stage4_dir": str(stage4_dir),
        },
        "py_compile": {"ok": bool(py_compile.get("ok", False))},
        "hard_errors": list(hard_errors),
        "warnings": list(warnings),
    }
    snippets: Dict[str, str] = {
        "head_sim_sc_metadata": _read_text_head(sim_out_dir / "sim_sc_metadata.csv", n_lines=8),
        "head_sim_st_coordinates": _read_text_head(sim_out_dir / "sim_st_coordinates.csv", n_lines=8),
        "head_truth_query_cell_spot": _read_text_head(sim_out_dir / "sim_truth_query_cell_spot.csv", n_lines=8),
        "head_truth_spot_type_fraction": _read_text_head(sim_out_dir / "sim_truth_spot_type_fraction.csv", n_lines=3),
    }

    # Also fingerprint input_files recorded in scenario_meta (if any)
    input_paths: List[Path] = []
    try:
        for rec in (scenario_meta.get("input_files") or []):
            p = rec.get("path")
            if p:
                input_paths.append(Path(str(p)))
    except Exception:
        input_paths = []
    input_paths = [p for p in input_paths if p.exists() and p.is_file()]
    if input_paths:
        extra = [fingerprint_file(p, project_root=project_root) for p in input_paths]
    else:
        extra = []

    return {
        "evidence": evidence,
        "schema_checks": schema_checks,
        "py_compile": py_compile,
        "stats_tables": {
            "counts_summary": counts_summary_rows,
            "per_type_counts": per_type_rows,
            "cells_per_spot_stats": cells_per_spot_rows,
        },
        "snippets": snippets,
        "fingerprints": fingerprints,
        "input_fingerprints": extra,
    }


def _render_markdown(result: Dict[str, Any], *, project_root: Path) -> str:
    ev = result.get("evidence", {})
    checks = result.get("schema_checks", {})
    tables = result.get("stats_tables", {})
    counts_summary = (tables.get("counts_summary") or [{}])[0]
    stage4 = (checks.get("stage4_alignment") or {})
    snippets = result.get("snippets", {})

    lines: List[str] = []
    lines.append(f"# SimGen Evidence (stdout): {ev.get('scenario_id')}")
    lines.append("")
    lines.append(f"- status: `{ev.get('status')}`")
    lines.append(f"- generated_at: `{ev.get('generated_at')}`")
    lines.append(f"- sim_out_dir: `{_relpath_str(Path(ev.get('paths', {}).get('sim_out_dir', '')), project_root)}`")
    lines.append(f"- stage1_export_dir: `{_relpath_str(Path(ev.get('paths', {}).get('stage1_export_dir', '')), project_root)}`")
    lines.append("")

    lines.append("## Summary")
    lines.append(
        f"- n_query_cells={counts_summary.get('n_query_cells')}, n_world_cells={counts_summary.get('n_world_cells')}, "
        f"n_spots={counts_summary.get('n_spots')}, n_types={counts_summary.get('n_types')}, "
        f"query_truth_na_fraction={counts_summary.get('query_truth_na_fraction')}"
    )
    lines.append(
        f"- coverage: min_cells_per_spot={counts_summary.get('coverage_min_cells_per_spot')}, "
        f"empty_spot_fraction={counts_summary.get('coverage_empty_spot_fraction')} (max={counts_summary.get('coverage_max_empty_spot_fraction')}), "
        f"min_cells_satisfied_fraction={counts_summary.get('coverage_min_cells_satisfied_fraction')} (warn<{counts_summary.get('coverage_min_cells_satisfied_warn')}), "
        f"in_tissue_fraction={counts_summary.get('coverage_in_tissue_fraction')}"
    )
    cov_sum = checks.get("alignment", {}).get("coverage_summary") or {}
    if cov_sum.get("fallback_action_counts") is not None:
        lines.append(
            f"- coverage_fallback_steps: n={cov_sum.get('fallback_steps_n')}, action_counts={cov_sum.get('fallback_action_counts')}, "
            f"actions_tail={cov_sum.get('fallback_actions_tail')}"
        )
    lines.append(f"- sc_missing_types={checks.get('scenario_semantics', {}).get('sc_missing_types')}")
    lines.append(f"- st_missing_types={checks.get('scenario_semantics', {}).get('st_missing_types')}")

    req = checks.get("required_files", {})
    lines.append("")
    lines.append("## Contracts")
    lines.append(f"- required sim outputs ok: `{req.get('sim', {}).get('ok')}` missing={req.get('sim', {}).get('missing')}")
    lines.append(f"- required stage1 exported ok: `{req.get('stage1_exported', {}).get('ok')}` missing={req.get('stage1_exported', {}).get('missing')}")
    lines.append(f"- py_compile ok: `{ev.get('py_compile', {}).get('ok')}`")
    if checks.get("alignment", {}).get("truth_alias_checks") is not None:
        lines.append(f"- truth alias checks: `{checks.get('alignment', {}).get('truth_alias_checks')}`")
    if checks.get("alignment", {}).get("stage1_sc_vs_sim_sc") is not None:
        lines.append(f"- stage1_sc_vs_sim_sc: `{checks.get('alignment', {}).get('stage1_sc_vs_sim_sc')}`")
    if checks.get("uniqueness"):
        lines.append(f"- cell_id duplicates: `{checks.get('uniqueness')}`")
    if checks.get("alignment", {}).get("stage1_gene_lists") is not None:
        lines.append(f"- stage1 gene lists: `{checks.get('alignment', {}).get('stage1_gene_lists')}`")

    # Missing-type focused rows
    per_type = tables.get("per_type_counts") or []
    sc_missing = list(checks.get("scenario_semantics", {}).get("sc_missing_types") or [])
    st_missing = list(checks.get("scenario_semantics", {}).get("st_missing_types") or [])
    if sc_missing or st_missing:
        lines.append("")
        lines.append("## Missing-Type Audit (types only)")
        want = [str(t) for t in (sc_missing + st_missing)]
        for t in want:
            row = next((r for r in per_type if str(r.get("type")) == t), None)
            if not row:
                lines.append(f"- {t}: (not found)")
                continue
            lines.append(
                f"- {t}: query_cell_count={row.get('query_cell_count')}, spot_truth_col_sum={row.get('spot_truth_col_sum')}, "
                f"query_truth_na_fraction={row.get('query_truth_na_fraction')}"
            )

    # Stage4 alignment
    lines.append("")
    lines.append("## Stage4 <-> Query truth alignment")
    if stage4.get("checked"):
        lines.append(f"- files_checked: `{len(stage4.get('files') or [])}`")
        lines.append(f"- max_cell_id_not_in_truth: `{stage4.get('max_cell_id_not_in_truth')}`")
        for r in stage4.get("files") or []:
            if "error" in r:
                lines.append(f"- {r.get('file')}: error={r.get('error')}")
            else:
                lines.append(
                    f"- {Path(str(r.get('file'))).name}: n_rows={r.get('n_rows')}, dup_rate={r.get('duplicate_cell_id_rate')}, "
                    f"not_in_truth={r.get('cell_id_not_in_truth')}"
                )
    else:
        lines.append("- stage4 outputs not found; skipped")

    # Snippets
    lines.append("")
    lines.append("## Snippets (head)")
    for k in ["head_sim_sc_metadata", "head_sim_st_coordinates", "head_truth_query_cell_spot", "head_truth_spot_type_fraction"]:
        lines.append(f"### {k}")
        lines.append("```")
        lines.append((snippets.get(k) or "").rstrip("\n"))
        lines.append("```")
        lines.append("")

    # Fingerprints
    lines.append("## File fingerprints (sha1/bytes/mtime)")
    for rec in result.get("fingerprints") or []:
        if not rec.get("exists"):
            lines.append(f"- {rec.get('path')}: (missing)")
        else:
            lines.append(f"- {rec.get('path')}: sha1={rec.get('sha1')} bytes={rec.get('bytes')} mtime={rec.get('mtime_iso')}")
    if result.get("input_fingerprints"):
        lines.append("")
        lines.append("## Input files (from scenario_meta.json)")
        for rec in result.get("input_fingerprints") or []:
            if not rec.get("exists"):
                lines.append(f"- {rec.get('path')}: (missing)")
            else:
                lines.append(f"- {rec.get('path')}: sha1={rec.get('sha1')} bytes={rec.get('bytes')} mtime={rec.get('mtime_iso')}")

    # Errors / warnings
    if ev.get("hard_errors"):
        lines.append("")
        lines.append("## Hard errors")
        for e in ev.get("hard_errors") or []:
            lines.append(f"- {e}")
    if ev.get("warnings"):
        lines.append("")
        lines.append("## Warnings")
        for w in ev.get("warnings") or []:
            lines.append(f"- {w}")

    # py_compile details if failed
    if not ev.get("py_compile", {}).get("ok"):
        lines.append("")
        lines.append("## py_compile details")
        for rec in (result.get("py_compile", {}).get("results") or []):
            if int(rec.get("returncode", 0)) == 0:
                continue
            lines.append(f"- {rec.get('file')}: returncode={rec.get('returncode')}")
            if rec.get("stderr"):
                lines.append(f"  stderr: {rec.get('stderr')}")
            if rec.get("stdout"):
                lines.append(f"  stdout: {rec.get('stdout')}")

    return "\n".join(lines).rstrip() + "\n"


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Print SimGen evidence report to stdout (no files written).")
    p.add_argument("--scenario_id", nargs="+", required=True, help="Scenario id(s), e.g. S0_matched")
    p.add_argument("--sim_root", default="data/sim", help="SimGen root dir (default: data/sim)")
    p.add_argument("--project_root", default=None, help="Project root (auto-detect by default)")
    p.add_argument("--format", default="md", choices=["md", "json"], help="Output format to stdout")
    return p.parse_args(list(argv) if argv is not None else None)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    project_root = Path(args.project_root).resolve() if args.project_root else _detect_project_root()
    sim_root = Path(args.sim_root)
    scenario_ids = list(args.scenario_id)

    script_path = Path(__file__).resolve()
    print(f"[evidence] script={script_path}")
    print(f"[evidence] script_sha1={_sha1_file(script_path)}")
    print(f"[evidence] project_root={project_root}")
    print(f"[evidence] sim_root={sim_root.as_posix()}")
    failed = 0
    for sid in scenario_ids:
        print(f"\n[evidence] scenario_id={sid}")
        result = build_evidence_for_scenario(project_root=project_root, scenario_id=sid, sim_root=sim_root)
        status = (result.get("evidence") or {}).get("status")
        if status != "success":
            failed += 1
        if str(args.format).lower() == "json":
            print(json.dumps(result, ensure_ascii=False, indent=2))
        else:
            print(_render_markdown(result, project_root=project_root))
    return 0 if failed == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
