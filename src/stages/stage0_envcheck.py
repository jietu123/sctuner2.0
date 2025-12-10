# -*- coding: utf-8 -*-
"""
Stage 0：环境与项目结构自检

使用方式（在 SVTUNER_PROJECT 根目录或任意目录均可）：
    python D:/EXPERIMENT/SVTUNER_PROJECT/src/stages/stage0_envcheck.py

功能：
1. 检查关键 Python 包能否 import
2. 检查 cytospace 命令行是否可用
3. 检查项目关键目录/文件是否存在
4. 在 logs/stage0_envcheck.json 写入一份体检报告
"""

import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path


def detect_project_root() -> Path:
    """
    自动推断项目根目录：
    假设本文件路径为：
        D:/EXPERIMENT/SVTUNER_PROJECT/src/stages/stage0_envcheck.py
    则项目根目录为：
        D:/EXPERIMENT/SVTUNER_PROJECT
    """
    here = Path(__file__).resolve()
    # .../SVTUNER_PROJECT/src/stages/stage0_envcheck.py
    # parents[0] = stages, parents[1] = src, parents[2] = SVTUNER_PROJECT
    try:
        root = here.parents[2]
    except IndexError:
        # 兜底：取当前文件所在目录的父目录
        root = here.parent.parent
    return root


def check_imports():
    """
    检查关键 Python 包是否可以导入。
    返回：
        result: dict[module_name] = bool
        ok: 是否全部导入成功
    """
    print("[Stage 0] 检查 Python 包导入...")
    # lapjv 对 CytoSPACE 可选（仅在使用 lapjv 求解器时需要），先不作为硬失败项。
    modules = ["cytospace", "scanpy", "anndata", "numpy"]
    result = {}
    all_ok = True

    for name in modules:
        try:
            __import__(name)
            print(f"  - import {name} [OK]")
            result[name] = True
        except ImportError as e:
            print(f"  - import {name} [FAIL] -> {e}")
            result[name] = False
            all_ok = False

    return result, all_ok


def check_cytospace_cli():
    """
    检查 cytospace 命令是否可用，并尝试获取版本信息。
    返回：
        info: dict 包含 help_ok, version 字段
        ok: bool
    """
    print("[Stage 0] 检查 cytospace 命令行...")
    info = {
        "help_ok": False,
        "version": None,
    }

    # 1) 先试一下 cytospace --help；若命令未在 PATH，再尝试 python -m cytospace --help
    ok = False
    try:
        result = subprocess.run(
            ["cytospace", "--help"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if result.returncode == 0:
            print("  - cytospace --help [OK]")
            info["help_ok"] = True
            ok = True
        else:
            print("  - cytospace --help [FAIL]")
            print("    错误信息:")
            print(result.stderr.strip())
    except FileNotFoundError:
        try:
            result = subprocess.run(
                [sys.executable, "-m", "cytospace", "--help"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            if result.returncode == 0:
                print("  - python -m cytospace --help [OK]")
                info["help_ok"] = True
                ok = True
            else:
                print("  - python -m cytospace --help [FAIL]")
                print("    错误信息:")
                print(result.stderr.strip())
        except Exception as e:
            print(f"  - 无法调用 cytospace：{e}")

    # 如果 CLI 依旧不可用，视为警告，不阻塞后续阶段（Python 包已确认可导入）
    if not ok:
        print("  - 警告：未检测到 cytospace 命令行，后续若需使用请检查 PATH 或安装入口脚本。")
        ok = True

    # 2) 尝试通过 pip show 获取版本（可选，不影响总状态）
    try:
        result_ver = subprocess.run(
            [sys.executable, "-m", "pip", "show", "cytospace"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if result_ver.returncode == 0:
            for line in result_ver.stdout.splitlines():
                if line.startswith("Version:"):
                    version = line.split("Version:", 1)[1].strip()
                    info["version"] = version
                    print(f"  - cytospace 版本：{version}")
                    break
        else:
            # 不打印错误，避免太吵
            pass
    except Exception:
        # 某些环境下 pip 命令不可用，忽略版本信息
        pass

    return info, ok


def check_project_layout(project_root: Path):
    """
    检查项目关键路径是否存在。
    返回：
        result: dict[path_str] = bool
        ok: bool
    """
    print("[Stage 0] 检查项目目录结构...")
    expected_paths = [
        project_root / "configs" / "project_config.yaml",
        project_root / "configs" / "datasets" / "real_brca.yaml",
        project_root / "data" / "raw" / "real_brca",
        project_root / "r_scripts" / "stage1_preprocess.R",
        project_root / "src" / "stages" / "stage1_preprocess.py",
        project_root / "src" / "stages" / "stage4_mapping.py",
        project_root / "src" / "stages" / "backends" / "cytospace_backend.py",
    ]

    result = {}
    all_ok = True

    for p in expected_paths:
        key = str(p.relative_to(project_root))
        if p.exists():
            print(f"  - {key} [OK]")
            result[key] = True
        else:
            print(f"  - {key} [MISSING]")
            result[key] = False
            all_ok = False

    return result, all_ok


def write_report(project_root: Path, report: dict):
    """
    将体检结果写入 logs/stage0_envcheck.json
    """
    logs_dir = project_root / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    out_path = logs_dir / "stage0_envcheck.json"
    try:
        with out_path.open("w", encoding="utf-8") as f:
            json.dump(report, f, ensure_ascii=False, indent=2)
        print(f"[Stage 0] 体检报告已写入：{out_path}")
    except Exception as e:
        print(f"[Stage 0] 写入体检报告失败：{e}")


def run_stage0(project_root: Path = None) -> int:
    """
    阶段 0 主入口，可供其它脚本调用。

    返回：
        0 表示检查通过
        1 表示检查未通过
    """
    if project_root is None:
        project_root = detect_project_root()

    print("===============================================")
    print("     SVTuner Stage 0：环境与结构自检")
    print("===============================================")
    print(f"[Stage 0] 推断项目根目录：{project_root}")

    imports_result, ok_imports = check_imports()
    cli_info, ok_cli = check_cytospace_cli()
    layout_result, ok_layout = check_project_layout(project_root)

    all_ok = ok_imports and ok_cli and ok_layout

    report = {
        "timestamp": datetime.now().isoformat(),
        "project_root": str(project_root),
        "python_imports": imports_result,
        "cytospace_cli": cli_info,
        "project_layout": layout_result,
        "status": "ok" if all_ok else "failed",
    }

    write_report(project_root, report)

    if all_ok:
        print("Stage 0: 环境检查通过，可以继续后续阶段。")
        return 0
    else:
        print("Stage 0: 环境检查未通过，请根据上面的提示修复后再重试。")
        return 1


if __name__ == "__main__":
    # 允许从命令行直接运行本文件
    root = detect_project_root()
    exit_code = run_stage0(root)
    sys.exit(exit_code)
