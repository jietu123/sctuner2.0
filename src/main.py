"""
SVTuner pipeline orchestrator.

Design goals (per SVTuner.md):
- CLI 控制样本 / 阶段 / 后端 / R 调用方式。
- 调度 Stage0（环境检测，Python）和 Stage1（预处理，R/Seurat）。
- 其余阶段预留接口，待后续接入（2/3/4/5/6/7）。
"""
from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys
from pathlib import Path

import config


# 基础工具
def parse_stage_list(stage_str: str) -> list[int]:
    if stage_str.lower() == "all":
        return list(range(0, 8))
    parts = [s.strip() for s in stage_str.split(",") if s.strip()]
    return [int(p) for p in parts]


def run_cmd(cmd, env=None, cwd=None) -> int:
    print(f"[RUN] {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=cwd, env=env)
    if result.returncode != 0:
        print(f"[ERR] Command failed with code {result.returncode}")
    return result.returncode


# Stage 调度
def run_stage0(project_root: Path) -> int:
    """Stage0: 极简占位检查，不做实际检测。"""
    print("[Stage0] 环境配置成功（占位检查，目前不执行实际环境检测）。")
    return 0


def run_stage1(project_root: Path, sample: str, r_cmd: str) -> int:
    """调用 R 预处理脚本（Stage1：Seurat 预处理）。"""
    script = str(project_root / "r_scripts" / "stage1_preprocess.R")
    cmd = shlex.split(r_cmd) + [script, "--sample", sample]
    env = dict(os.environ)
    env.setdefault("CONDA_DLL_SEARCH_MODIFICATION_ENABLE", "1")
    return run_cmd(cmd, env=env, cwd=str(project_root))


def main(argv=None):
    parser = argparse.ArgumentParser(description="SVTuner pipeline orchestrator")
    parser.add_argument("--sample", default="real_brca", help="Sample ID (e.g., real_brca, S0_matched)")
    parser.add_argument(
        "--stages",
        default=config.DEFAULT_STAGES,
        help=(
            "Comma-separated stage numbers (0-7) or 'all'. "
            "Implemented: 1 (preprocess). "
            "Stage 0 is an optional dummy env-check placeholder."
        ),
    )
    parser.add_argument(
        "--r-cmd",
        default=config.DEFAULT_R_CMD,
        help='Command to invoke Rscript (e.g., "Rscript" or "conda run -n <env> Rscript")',
    )
    parser.add_argument(
        "--backends",
        default="cytospace",
        help="Comma-separated backends for later stages (placeholder for future).",
    )
    parser.add_argument(
        "--config",
        help="Optional extra YAML config (placeholder; Stage1 自行读取 configs/*.yaml).",
    )
    args = parser.parse_args(argv)

    project_root = config.detect_project_root()
    print(f"[Main] Project root: {project_root}")
    stages = parse_stage_list(args.stages)

    for st in stages:
        if st == 0:
            print("[Main] Running Stage 0: Env check")
            code = run_stage0(project_root)
        elif st == 1:
            print("[Main] Running Stage 1: Preprocess (R/Seurat)")
            code = run_stage1(project_root, args.sample, args.r_cmd)
        else:
            print(f"[Main] Stage {st} not implemented in main.py. Skipping.")
            code = 0

        if code != 0:
            print(f"[Main] Stage {st} failed with code {code}. Aborting pipeline.")
            sys.exit(code)

    print("[Main] Pipeline completed.")


if __name__ == "__main__":
    main()
