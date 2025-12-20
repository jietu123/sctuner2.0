"""
Stage 4: Mapping orchestrator
- 跑 baseline / plus 两条分支
- 支持多 backend（首个 cytospace），统一输出目录与文件命名
- baseline 不用插件，plus 必须用 Stage2/3 插件
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path
from typing import Dict, Any, List, Optional

import sys

# 确保可以从项目根导入 src 包
PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "src"
for p in (PROJECT_ROOT, SRC_ROOT):
    if str(p) not in sys.path:
        sys.path.append(str(p))

from stages.backends.base_backend import MappingBackend
from stages.backends.cytospace_backend import CytoSPACEBackend
from config import load_project_config, ProjectConfig
from utils.run_manifest import write_stage4_run_manifest


BACKEND_REGISTRY = {
    "cytospace": CytoSPACEBackend,
    # "tangram": TangramBackend,  # 预留
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Stage 4: mapping orchestrator")
    p.add_argument("--sample", default="real_brca")
    p.add_argument("--project_root", default=None)
    p.add_argument("--backends", default=None, help="comma-separated backend names; default uses project_config")
    p.add_argument("--enable_plus", action="store_true", default=False)
    p.add_argument("--config_id", default=None)
    return p.parse_args()


def detect_project_root() -> Path:
    here = Path(__file__).resolve()
    try:
        return here.parents[2]
    except IndexError:
        return here.parent.parent


def runner_fingerprint() -> Dict[str, Optional[str]]:
    p = Path(__file__).resolve()
    try:
        sha1 = hashlib.sha1(p.read_bytes()).hexdigest().lower()
    except Exception:
        sha1 = None
    return {"runner_file": str(p), "runner_sha1": sha1}


def _format_run_token(value: Any) -> str:
    if value is None or value == "":
        return "NA"
    s = str(value).strip()
    if not s:
        return "NA"
    for ch in ["\\", "/", ":", "|"]:
        s = s.replace(ch, "-")
    return s.replace(" ", "_")


def _format_float_token(value: Any) -> str:
    if value is None or value == "":
        return "NA"
    try:
        f = float(value)
    except Exception:
        return _format_run_token(value)
    if not math.isfinite(f):
        return "NA"
    return _format_run_token(f"{f:.6g}")


def _resolve_refine_enabled(lambda_svg: Any) -> bool:
    try:
        return lambda_svg is not None and float(lambda_svg) > 0.0
    except Exception:
        return False


def _build_run_id(
    *,
    mode: str,
    variant: Optional[str],
    seed: Any,
    svg_refine_lambda: Any,
    config_id: Any,
    refine_enabled: bool,
) -> str:
    prefix = f"{mode}_{variant}" if variant else mode
    return (
        f"{prefix}"
        f"__seed_{_format_run_token(seed)}"
        f"__lam_{_format_float_token(svg_refine_lambda)}"
        f"__ref_{'1' if refine_enabled else '0'}"
        f"__cfg_{_format_run_token(config_id)}"
    )


def ensure_inputs(stage1_dir: Path, stage2_dir: Path, stage3_dir: Path, enable_plus: bool):
    # 基本存在性检查
    exp_dir = stage1_dir / "exported"
    for fname in ["sc_expression_normalized.csv", "st_expression_normalized.csv", "sc_metadata.csv", "st_coordinates.csv"]:
        f = exp_dir / fname
        if not f.exists():
            raise FileNotFoundError(f"缺少 Stage1 文件: {f}")
    if enable_plus:
        if not stage2_dir.exists():
            raise FileNotFoundError(f"缺少 Stage2 目录: {stage2_dir}")
        if not stage3_dir.exists():
            raise FileNotFoundError(f"缺少 Stage3 目录: {stage3_dir}")
        for f in [stage2_dir / "plugin_genes.txt", stage2_dir / "gene_weights.csv"]:
            if not f.exists():
                raise FileNotFoundError(f"缺少 Stage2 文件: {f}")
        for f in [stage3_dir / "cell_type_relabel.csv", stage3_dir / "type_prior_matrix.csv"]:
            if not f.exists():
                raise FileNotFoundError(f"缺少 Stage3 文件: {f}")


def run_backend(backend: MappingBackend, name: str, stage1_dir: Path, stage2_dir: Path, stage3_dir: Path, out_root: Path, cfg: Dict[str, Any], enable_plus: bool):
    # baseline
    base_cfg = cfg | {"mode": "baseline", "variant": None, "svg_refine_lambda": None}
    base_run_id = _build_run_id(
        mode="baseline",
        variant=None,
        seed=base_cfg.get("seed"),
        svg_refine_lambda=None,
        config_id=base_cfg.get("config_id"),
        refine_enabled=False,
    )
    base_cfg = base_cfg | {"run_id": base_run_id}
    out_base = out_root / name / base_run_id
    out_base.mkdir(parents=True, exist_ok=True)
    print(f"[Stage4] Running {name} baseline -> {out_base}")
    log_base = out_base / "log.txt"
    with log_base.open("w", encoding="utf-8") as f:
        f.write(f"baseline backend={name}\n")
        f.write(f"mode={base_cfg['mode']} run_id={base_cfg['run_id']} variant={base_cfg['variant']}\n")
        f.write(json.dumps(base_cfg, ensure_ascii=False, indent=2) + "\n")
        try:
            backend.run_baseline(stage1_dir=stage1_dir, out_dir=out_base, config=base_cfg)
            f.write("status=success\n")
        except Exception as e:
            f.write(f"status=failed\nerror={e}\n")
            raise

    # Ensure meta.json exists even for backends that don't emit it.
    meta_path = out_base / "meta.json"
    if not meta_path.exists():
        meta = {
            "backend": name,
            "mode": base_cfg.get("mode"),
            "run_id": base_cfg.get("run_id"),
            "variant": base_cfg.get("variant"),
            "sample": base_cfg.get("sample"),
            "config_id": base_cfg.get("config_id"),
            "seed": base_cfg.get("seed"),
            "svg_refine_lambda": base_cfg.get("svg_refine_lambda"),
            "status": "success",
            "resolved_mapping_config": base_cfg,
        }
        meta_path.write_text(json.dumps(meta, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    if enable_plus:
        plus_cfg = cfg | {"mode": "plus", "variant": "svg_type"}
        lambda_svg = plus_cfg.get("svg_refine_lambda")
        plus_run_id = _build_run_id(
            mode="plus",
            variant="svg_type",
            seed=plus_cfg.get("seed"),
            svg_refine_lambda=lambda_svg,
            config_id=plus_cfg.get("config_id"),
            refine_enabled=_resolve_refine_enabled(lambda_svg),
        )
        plus_cfg = plus_cfg | {"run_id": plus_run_id}
        out_plus = out_root / name / plus_run_id
        out_plus.mkdir(parents=True, exist_ok=True)
        print(f"[Stage4] Running {name} plus_svg_type -> {out_plus}")
        log_plus = out_plus / "log.txt"
        with log_plus.open("w", encoding="utf-8") as f:
            f.write(f"plus backend={name}\n")
            f.write(f"mode={plus_cfg['mode']} run_id={plus_cfg['run_id']} variant={plus_cfg['variant']}\n")
            f.write(json.dumps(plus_cfg, ensure_ascii=False, indent=2) + "\n")
            try:
                backend.run_plus(
                    stage1_dir=stage1_dir,
                    stage2_dir=stage2_dir,
                    stage3_dir=stage3_dir,
                    out_dir=out_plus,
                    config=plus_cfg,
                )
                try:
                    meta = json.loads((out_plus / "meta.json").read_text(encoding="utf-8"))
                    norm_applied = meta.get("norm_applied")
                    norm_fallback = meta.get("norm_fallback")
                    norm_reason = meta.get("norm_reason")
                    f.write(f"norm_applied={norm_applied} norm_fallback={norm_fallback} norm_reason={norm_reason}\n")
                except Exception as meta_err:
                    f.write(f"norm_meta_read_error={meta_err}\n")
                f.write("status=success\n")
            except Exception as e:
                f.write(f"status=failed\nerror={e}\n")
                raise

        meta_path = out_plus / "meta.json"
        if not meta_path.exists():
            meta = {
                "backend": name,
                "mode": plus_cfg.get("mode"),
                "run_id": plus_cfg.get("run_id"),
                "variant": plus_cfg.get("variant"),
                "sample": plus_cfg.get("sample"),
                "config_id": plus_cfg.get("config_id"),
                "seed": plus_cfg.get("seed"),
                "svg_refine_lambda": plus_cfg.get("svg_refine_lambda"),
                "status": "success",
                "resolved_mapping_config": plus_cfg,
            }
            meta_path.write_text(json.dumps(meta, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def main():
    args = parse_args()
    project_root = Path(args.project_root) if args.project_root else detect_project_root()
    cfg = load_project_config(project_root)

    if args.backends:
        backend_names = [b.strip() for b in args.backends.split(",") if b.strip()]
    else:
        backend_names = None

    run_stage4_for_sample(
        cfg=cfg,
        sample_id=args.sample,
        backends=backend_names,
        enable_plus=args.enable_plus,
        config_id=args.config_id,
    )

    print("[Stage4] All done.")


def run_stage4_for_sample(
    cfg: ProjectConfig,
    sample_id: str,
    backends: Optional[List[str]] = None,
    enable_plus: bool = False,
    config_id: Optional[str] = None,
) -> None:
    ds_cfg = cfg.load_dataset_cfg(sample_id) or {}
    stage1_dir = cfg.get_stage_dir(sample_id, "stage1_preprocess")
    stage2_dir = cfg.get_stage_dir(sample_id, "stage2_svg_plugin")
    stage3_dir = cfg.get_stage_dir(sample_id, "stage3_typematch")
    out_root = cfg.get_stage_dir(sample_id, "stage4_mapping")

    ensure_inputs(stage1_dir, stage2_dir, stage3_dir, enable_plus=enable_plus)

    backend_list = backends if backends else cfg.enabled_backends(sample_id)
    run_fp = runner_fingerprint()
    for bname in backend_list:
        if bname not in BACKEND_REGISTRY:
            raise ValueError(f"未知 backend: {bname}")
        ds_map = ds_cfg.get("mapping", {}) if isinstance(ds_cfg, dict) else {}
        ds_backend = ds_map.get(bname, {}) if isinstance(ds_map, dict) else {}
        default_cfg_id = ds_backend.get("default_config_id")
        effective_config_id = config_id if config_id is not None else default_cfg_id
        backend_cfg = cfg.get_mapping_config(sample_id, bname, config_id=effective_config_id)
        backend_cfg = backend_cfg | {
            "sample": sample_id,
            "backend": bname,
            "config_id": effective_config_id,
            "project_root": str(cfg.project_root),
            "project_config_path": str(cfg.project_root / "configs" / "project_config.yaml"),
            "dataset_config_path": str(cfg.dataset_cfg_path(sample_id)),
        } | run_fp
        backend = BACKEND_REGISTRY[bname]()
        run_backend(
            backend=backend,
            name=bname,
            stage1_dir=stage1_dir,
            stage2_dir=stage2_dir,
            stage3_dir=stage3_dir,
            out_root=out_root,
            cfg=backend_cfg,
            enable_plus=enable_plus,
        )

    # Stage4 manifest is the single source of truth for Stage5/6 run enumeration.
    out_json, out_csv = write_stage4_run_manifest(
        project_root=cfg.project_root,
        scenario_id=sample_id,
        stage4_mapping_root=out_root,
    )
    print(f"[Stage4] Wrote manifest: {out_json}")
    print(f"[Stage4] Wrote manifest: {out_csv}")


if __name__ == "__main__":
    main()
