"""
Stage4 run manifest utilities.

Stage5/6 must enumerate runs ONLY from this manifest (no directory guessing),
per SVTuner.md planning docs.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


MANIFEST_VERSION_STAGE4 = "stage4_run_manifest_v1"


def _now_iso() -> str:
    return datetime.now().isoformat(timespec="seconds")


def _read_json(path: Path) -> Dict[str, Any]:
    if not path.exists() or path.stat().st_size == 0:
        return {}
    return json.loads(path.read_text(encoding="utf-8")) or {}


def _relpath_posix(path: Path, project_root: Path) -> str:
    try:
        return path.resolve().relative_to(project_root.resolve()).as_posix()
    except Exception:
        return path.resolve().as_posix()


def stage4_manifest_paths(*, project_root: Path, scenario_id: str) -> Tuple[Path, Path]:
    scenario_root = project_root / "result" / scenario_id
    return (scenario_root / "stage4_run_manifest.json", scenario_root / "stage4_run_manifest.csv")


@dataclass(frozen=True)
class Stage4RunRecord:
    scenario_id: str
    backend: str
    mode: str
    run_id: str
    variant: Optional[str]
    run_dir: str
    meta_path: Optional[str]
    cell_assignment_path: str
    cell_spot_matrix_path: str
    spot_type_fraction_path: str
    config_id: Optional[str]
    seed: Optional[int]
    svg_refine_lambda: Optional[float]
    lambda_svg: Optional[float]
    refine_enabled: bool

    cell_assignment_exists: bool
    cell_spot_matrix_exists: bool
    spot_type_fraction_exists: bool
    meta_exists: bool

    def to_row(self) -> Dict[str, Any]:
        return {
            "scenario_id": self.scenario_id,
            "backend": self.backend,
            "mode": self.mode,
            "run_id": self.run_id,
            "variant": self.variant,
            "run_dir": self.run_dir,
            "meta_path": self.meta_path,
            "cell_assignment_path": self.cell_assignment_path,
            "cell_spot_matrix_path": self.cell_spot_matrix_path,
            "spot_type_fraction_path": self.spot_type_fraction_path,
            "config_id": self.config_id,
            "seed": self.seed,
            "svg_refine_lambda": self.svg_refine_lambda,
            "lambda_svg": self.lambda_svg,
            "refine_enabled": self.refine_enabled,
            "cell_assignment_exists": self.cell_assignment_exists,
            "cell_spot_matrix_exists": self.cell_spot_matrix_exists,
            "spot_type_fraction_exists": self.spot_type_fraction_exists,
            "meta_exists": self.meta_exists,
        }


def _coerce_int(v: Any) -> Optional[int]:
    if v in (None, ""):
        return None
    try:
        return int(v)
    except Exception:
        return None


def _coerce_float(v: Any) -> Optional[float]:
    if v in (None, ""):
        return None
    try:
        return float(v)
    except Exception:
        return None


def _coerce_str(v: Any) -> Optional[str]:
    if v is None:
        return None
    s = str(v)
    return s


def _resolve_run_outputs(run_dir: Path, mode: str) -> Tuple[Path, Path, Path]:
    return (
        run_dir / f"cell_assignment_{mode}.csv",
        run_dir / f"cell_spot_matrix_{mode}.npz",
        run_dir / f"spot_type_fraction_{mode}.csv",
    )


def build_stage4_run_manifest(
    *,
    project_root: Path,
    scenario_id: str,
    stage4_mapping_root: Path,
) -> Dict[str, Any]:
    """
    Build a Stage4 run manifest by scanning stage4_mapping outputs.

    Stage5/6 MUST consume this manifest only; Stage4 is the only writer.
    """
    runs: List[Stage4RunRecord] = []
    if stage4_mapping_root.exists():
        for backend_dir in sorted([p for p in stage4_mapping_root.iterdir() if p.is_dir()], key=lambda p: p.name.lower()):
            backend_name = backend_dir.name
            for run_dir in sorted([p for p in backend_dir.iterdir() if p.is_dir()], key=lambda p: p.name.lower()):
                meta_path = run_dir / "meta.json"
                meta = _read_json(meta_path)

                mode = str(meta.get("mode") or ("baseline" if "baseline" in run_dir.name else "plus"))
                run_id = str(meta.get("run_id") or run_dir.name)
                variant = meta.get("variant")
                config_id = _coerce_str(meta.get("config_id"))
                seed = _coerce_int(meta.get("seed"))
                svg_refine_lambda = _coerce_float(meta.get("svg_refine_lambda"))
                if svg_refine_lambda is None:
                    svg_refine_lambda = _coerce_float((meta.get("config_effective_subset") or {}).get("svg_refine_lambda"))
                if mode == "baseline":
                    svg_refine_lambda = None
                    lambda_svg = None
                    refine_enabled = False
                else:
                    lambda_svg = svg_refine_lambda
                    refine_enabled = bool(lambda_svg is not None and float(lambda_svg) > 0.0)

                ca_path, mat_path, stf_path = _resolve_run_outputs(run_dir, mode)

                runs.append(
                    Stage4RunRecord(
                        scenario_id=str(meta.get("sample") or meta.get("scenario_id") or scenario_id),
                        backend=str(meta.get("backend") or backend_name),
                        mode=mode,
                        run_id=run_id,
                        variant=_coerce_str(variant),
                        run_dir=_relpath_posix(run_dir, project_root),
                        meta_path=_relpath_posix(meta_path, project_root) if meta_path.exists() else None,
                        cell_assignment_path=_relpath_posix(ca_path, project_root),
                        cell_spot_matrix_path=_relpath_posix(mat_path, project_root),
                        spot_type_fraction_path=_relpath_posix(stf_path, project_root),
                        config_id=config_id,
                        seed=seed,
                        svg_refine_lambda=svg_refine_lambda,
                        lambda_svg=lambda_svg,
                        refine_enabled=refine_enabled,
                        cell_assignment_exists=bool(ca_path.exists()),
                        cell_spot_matrix_exists=bool(mat_path.exists()),
                        spot_type_fraction_exists=bool(stf_path.exists()),
                        meta_exists=bool(meta_path.exists()),
                    )
                )

    payload = {
        "manifest_version": MANIFEST_VERSION_STAGE4,
        "generated_at": _now_iso(),
        "scenario_id": scenario_id,
        "stage4_mapping_root": _relpath_posix(stage4_mapping_root, project_root),
        "runs": [r.to_row() for r in runs],
    }
    return payload


def write_stage4_run_manifest(
    *,
    project_root: Path,
    scenario_id: str,
    stage4_mapping_root: Path,
) -> Tuple[Path, Path]:
    """
    Write stage4_run_manifest.{json,csv} under result/<scenario_id>/.
    """
    out_json, out_csv = stage4_manifest_paths(project_root=project_root, scenario_id=scenario_id)
    out_json.parent.mkdir(parents=True, exist_ok=True)

    payload = build_stage4_run_manifest(project_root=project_root, scenario_id=scenario_id, stage4_mapping_root=stage4_mapping_root)
    out_json.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    runs = payload.get("runs") or []
    fieldnames = [
        "scenario_id",
        "backend",
        "mode",
        "run_id",
        "variant",
        "run_dir",
        "meta_path",
        "cell_assignment_path",
        "cell_spot_matrix_path",
        "spot_type_fraction_path",
        "config_id",
        "seed",
        "svg_refine_lambda",
        "lambda_svg",
        "refine_enabled",
        "cell_assignment_exists",
        "cell_spot_matrix_exists",
        "spot_type_fraction_exists",
        "meta_exists",
    ]
    with out_csv.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in runs:
            row = {k: r.get(k) for k in fieldnames}
            w.writerow(row)

    return out_json, out_csv


def load_stage4_run_manifest(
    *,
    project_root: Path,
    scenario_id: str,
) -> Dict[str, Any]:
    """
    Load stage4_run_manifest for Stage5/6 consumption. Prefers JSON, falls back to CSV.
    """
    manifest_json, manifest_csv = stage4_manifest_paths(project_root=project_root, scenario_id=scenario_id)
    if manifest_json.exists():
        payload = _read_json(manifest_json)
        if payload.get("manifest_version") != MANIFEST_VERSION_STAGE4:
            raise ValueError(
                f"Unsupported Stage4 manifest_version={payload.get('manifest_version')!r}; expected {MANIFEST_VERSION_STAGE4!r}"
            )
        return payload

    if manifest_csv.exists():
        # Minimal CSV fallback: convert to JSON-like structure
        rows: List[Dict[str, Any]] = []
        with manifest_csv.open("r", encoding="utf-8", newline="") as f:
            r = csv.DictReader(f)
            for row in r:
                rows.append(dict(row))
        return {
            "manifest_version": MANIFEST_VERSION_STAGE4,
            "generated_at": None,
            "scenario_id": scenario_id,
            "stage4_mapping_root": _relpath_posix(project_root / "result" / scenario_id / "stage4_mapping", project_root),
            "runs": rows,
        }

    raise FileNotFoundError(f"Stage4 run manifest not found under result/{scenario_id}/ (expected stage4_run_manifest.json or .csv)")
