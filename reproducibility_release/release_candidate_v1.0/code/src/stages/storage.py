from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


DEFAULT_LOW_RES_GROUP = "low_resolution_experiments"


def read_dataset_config(project_root: Path, sample: str) -> dict[str, Any]:
    cfg_path = project_root / "configs" / "datasets" / f"{sample}.yaml"
    if not cfg_path.exists():
        return {}
    return yaml.safe_load(cfg_path.read_text(encoding="utf-8")) or {}


def sample_group(cfg: dict[str, Any] | None) -> str | None:
    storage = (cfg or {}).get("storage", {}) or {}
    group = storage.get("group")
    if group is None:
        return None
    group = str(group).strip().strip("/\\")
    return group or None


def sample_dir(project_root: Path, kind: str, sample: str, cfg: dict[str, Any] | None = None) -> Path:
    group = sample_group(cfg)
    if kind == "raw":
        base = project_root / "data" / "raw"
    elif kind == "processed":
        base = project_root / "data" / "processed"
    elif kind == "result":
        base = project_root / "result"
    else:
        raise ValueError(f"Unknown storage kind: {kind}")
    return base / group / sample if group else base / sample


def raw_dir(project_root: Path, sample: str, cfg: dict[str, Any] | None = None) -> Path:
    return sample_dir(project_root, "raw", sample, cfg)


def processed_dir(project_root: Path, sample: str, cfg: dict[str, Any] | None = None) -> Path:
    return sample_dir(project_root, "processed", sample, cfg)


def result_dir(project_root: Path, sample: str, cfg: dict[str, Any] | None = None) -> Path:
    return sample_dir(project_root, "result", sample, cfg)


def stage1_dir(project_root: Path, sample: str, cfg: dict[str, Any] | None = None) -> Path:
    return processed_dir(project_root, sample, cfg) / "stage1_preprocess"


def stage1_export_dir(project_root: Path, sample: str, cfg: dict[str, Any] | None = None) -> Path:
    return stage1_dir(project_root, sample, cfg) / "exported"

