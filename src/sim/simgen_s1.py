"""
SimGen S1 占位：暂不实现，预留接口。

目标（未来）：构造 SC 缺失类型B、ST 缺失类型A 的双向缺失场景。
当前：仅打印提示并退出。
"""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="SimGen S1 placeholder (not implemented)")
    p.add_argument("--sample", default="real_brca", help="sample id")
    p.add_argument("--project_root", default=None, help="override project root")
    p.add_argument("--config", default="configs/project_config.yaml", help="project config path")
    p.add_argument("--dataset_config", default=None, help="dataset config path")
    p.add_argument("--sim_config", default=None, help="simgen config path for S1 (optional)")
    return p.parse_args()


def main():
    args = parse_args()
    here = Path(__file__).resolve()
    project_root = Path(args.project_root) if args.project_root else here.parents[2]
    sim_cfg = project_root / (args.sim_config or "configs/simgen/s1.yaml")
    print(f"[SimGen S1] 占位，未实现。项目根: {project_root}")
    print(f"[SimGen S1] 可在此扩展，读取配置: {sim_cfg}")


if __name__ == "__main__":
    main()

