"""
Route 2 Stage4 (MVP):
- 读取 Stage3 输出 (cell_type_relabel.csv, type_support.csv)
- 在构建 cell pool 前过滤掉 plugin_type == Unknown_sc_only（或 orig_type 对应 unsupported）
- 输出简单的 assignment 占位（未做线性指派），仅用于验收过滤效果

验收指标写入 stage4_summary.json:
- n_cells_before_prefilter
- n_cells_after_prefilter
- n_filtered_unknown
- missing_present_in_assignment（指定 missing_type 的细胞在分配表中的计数）
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Optional

import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Route2 Stage4 MVP (filter Unknown/unsupported before assignment)")
    p.add_argument("--sample", default="real_brca_simS0_seed42", help="sample id")
    p.add_argument("--project_root", default=None, help="override project root")
    p.add_argument("--missing_type", default="T cells CD8", help="target missing type for reporting")
    p.add_argument("--filter_mode", choices=["plugin_unknown", "unsupported"], default="plugin_unknown",
                   help="过滤规则: plugin_unknown=移除 plugin_type==Unknown_sc_only; unsupported=移除 support_category==unsupported")
    return p.parse_args()


def detect_root() -> Path:
    here = Path(__file__).resolve()
    try:
        return here.parents[2]
    except IndexError:
        return here.parent.parent


def load_stage3_outputs(root: Path, sample: str) -> tuple[pd.DataFrame, Optional[pd.DataFrame]]:
    base = root / "data" / "processed" / sample / "stage3_typematch"
    relabel = pd.read_csv(base / "cell_type_relabel.csv")
    ts_path = base / "type_support.csv"
    ts = pd.read_csv(ts_path) if ts_path.exists() else None
    return relabel, ts


def main():
    args = parse_args()
    project_root = Path(args.project_root) if args.project_root else detect_root()

    relabel, type_support = load_stage3_outputs(project_root, args.sample)
    n_before = len(relabel)

    if args.filter_mode == "plugin_unknown":
        mask_drop = relabel["plugin_type"] == "Unknown_sc_only"
    else:
        unsupported_types = set(type_support.loc[type_support["support_category"] == "unsupported", "orig_type"]) if type_support is not None else set()
        mask_drop = relabel["orig_type"].isin(unsupported_types)

    filtered = relabel.loc[~mask_drop].copy()
    n_after = len(filtered)
    n_filtered = mask_drop.sum()

    # assignment 占位：未做指派，仅保留过滤后的细胞
    assignment = filtered[["cell_id", "plugin_type", "orig_type"]].copy()
    assignment["assigned_spot"] = None

    # 统计 missing_type 是否仍在 assignment
    missing_in_assign = (assignment["orig_type"] == args.missing_type).sum()

    out_dir = project_root / "result" / args.sample / "stage4_route2"
    out_dir.mkdir(parents=True, exist_ok=True)
    assignment.to_csv(out_dir / "cell_assignment.csv", index=False)

    summary = {
        "sample": args.sample,
        "filter_mode": args.filter_mode,
        "missing_type": args.missing_type,
        "n_cells_before_prefilter": int(n_before),
        "n_cells_after_prefilter": int(n_after),
        "n_filtered_unknown": int(n_filtered),
        "missing_present_in_assignment": int(missing_in_assign),
    }
    (out_dir / "stage4_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print("[Stage4 route2] summary:")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

