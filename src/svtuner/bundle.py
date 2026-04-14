from __future__ import annotations

import json
import shutil
import zipfile
from datetime import datetime
from pathlib import Path
from typing import Iterable


DEFAULT_INCLUDE_FILES = [
    "README.md",
    "pyproject.toml",
    ".gitignore",
]

DEFAULT_INCLUDE_DIRS = [
    "configs",
    "docs",
    "src",
    "scripts",
    "r_scripts",
    "external/cytospace",
]

SKIP_DIR_NAMES = {
    ".git",
    ".idea",
    ".vscode",
    "__pycache__",
}

SKIP_FILE_SUFFIXES = {
    ".pyc",
    ".pyo",
    ".pyd",
}


def _iter_files(base: Path) -> Iterable[Path]:
    for p in base.rglob("*"):
        if not p.is_file():
            continue
        parts = set(p.parts)
        if parts & SKIP_DIR_NAMES:
            continue
        if p.suffix.lower() in SKIP_FILE_SUFFIXES:
            continue
        yield p


def _copy_file(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def _should_skip_external_large(rel_from_project: Path) -> bool:
    s = str(rel_from_project).replace("\\", "/")
    if s.startswith("external/cytospace/data/"):
        return True
    if s.startswith("external/cytospace/images/"):
        return True
    if "cytospace_results" in s:
        return True
    return False


def create_bundle(
    project_root: Path,
    out_dir: Path,
    bundle_name: str,
    *,
    include_raw_data: bool = False,
    include_results: bool = False,
    include_external: bool = True,
) -> dict:
    project_root = project_root.resolve()
    out_dir = out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    stage_dir = out_dir / bundle_name
    if stage_dir.exists():
        shutil.rmtree(stage_dir)
    stage_dir.mkdir(parents=True, exist_ok=True)

    include_files = list(DEFAULT_INCLUDE_FILES)
    include_dirs = [d for d in DEFAULT_INCLUDE_DIRS if include_external or d != "external/cytospace"]
    if include_raw_data:
        include_dirs.append("data/raw")
    if include_results:
        include_dirs.append("result")

    copied_files = 0
    copied_bytes = 0

    for rel in include_files:
        src = project_root / rel
        if not src.exists() or not src.is_file():
            continue
        dst = stage_dir / rel
        _copy_file(src, dst)
        copied_files += 1
        copied_bytes += src.stat().st_size

    for rel in include_dirs:
        src_dir = project_root / rel
        if not src_dir.exists() or not src_dir.is_dir():
            continue
        for src_file in _iter_files(src_dir):
            rel_file = src_file.relative_to(project_root)
            if _should_skip_external_large(rel_file):
                continue
            dst_file = stage_dir / rel_file
            _copy_file(src_file, dst_file)
            copied_files += 1
            copied_bytes += src_file.stat().st_size

    manifest = {
        "bundle_name": bundle_name,
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "project_root": str(project_root),
        "options": {
            "include_raw_data": include_raw_data,
            "include_results": include_results,
            "include_external": include_external,
        },
        "copied_files": copied_files,
        "copied_bytes": copied_bytes,
        "included_files": include_files,
        "included_dirs": include_dirs,
    }
    manifest_path = stage_dir / "BUNDLE_MANIFEST.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    zip_path = out_dir / f"{bundle_name}.zip"
    if zip_path.exists():
        zip_path.unlink()

    with zipfile.ZipFile(zip_path, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:
        for fp in _iter_files(stage_dir):
            arcname = fp.relative_to(stage_dir)
            zf.write(fp, arcname=str(arcname))

    return {
        "stage_dir": str(stage_dir),
        "zip_path": str(zip_path),
        "manifest": str(manifest_path),
        "copied_files": copied_files,
        "copied_bytes": copied_bytes,
    }

