#!/usr/bin/env python3
"""
清理冗余存储：删除 pipeline 未使用的 Stage4 大文件

可安全删除（Stage5 及后续不依赖）：
- assigned_expression/ (matrix.mtx 等) ~146 GB
- sc_expression_for_cytospace.csv ~22 GB
- st_expression_for_cytospace.csv ~14 GB

Usage:
    python scripts/cleanup_redundant_storage.py --dry-run   # 仅列出，不删除
    python scripts/cleanup_redundant_storage.py            # 执行删除
"""
from __future__ import annotations

import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "result"


def find_redundant() -> list[tuple[Path, str]]:
    """返回 (路径, 类型) 列表。"""
    items: list[tuple[Path, str]] = []
    if not RESULT.exists():
        return items

    for sample_dir in RESULT.iterdir():
        if not sample_dir.is_dir():
            continue
        for stage4 in sample_dir.glob("stage4_cytospace_*"):
            cyto_out = stage4 / "cytospace_output"
            cyto_in = stage4 / "cytospace_input"
            if not cyto_out.exists():
                continue
            # assigned_expression 目录 (cytospace_output)
            ae = cyto_out / "assigned_expression"
            if ae.exists() and ae.is_dir():
                items.append((ae, "assigned_expression"))
            # 两个大 CSV：可能在 cytospace_output 或 cytospace_input
            for name in ("sc_expression_for_cytospace.csv", "st_expression_for_cytospace.csv"):
                for base in (cyto_out, cyto_in):
                    if base and base.exists():
                        p = base / name
                        if p.exists():
                            items.append((p, name))
                            break  # 每个 name 只加一次
    return items


def dir_size(p: Path) -> int:
    if not p.exists():
        return 0
    if p.is_file():
        return p.stat().st_size
    total = 0
    for f in p.rglob("*"):
        if f.is_file():
            total += f.stat().st_size
    return total


def main():
    ap = argparse.ArgumentParser(description="Clean redundant Stage4 files")
    ap.add_argument("--dry-run", action="store_true", help="Only list, do not delete")
    ap.add_argument("-y", "--yes", action="store_true", help="Skip confirmation prompt")
    args = ap.parse_args()

    items = find_redundant()
    if not items:
        print("[INFO] No redundant files found.")
        return

    total_bytes = 0
    by_type: dict[str, list[Path]] = {}
    for p, typ in items:
        total_bytes += dir_size(p)
        by_type.setdefault(typ, []).append(p)

    print("=" * 60)
    print("Redundant files to remove:")
    print("=" * 60)
    for typ, paths in sorted(by_type.items()):
        n = len(paths)
        sz = sum(dir_size(x) for x in paths)
        print(f"  {typ}: {n} items, {sz / 1e9:.2f} GB")
    print("-" * 60)
    print(f"  TOTAL: {total_bytes / 1e9:.2f} GB")
    print("=" * 60)

    if args.dry_run:
        print("\n[DRY-RUN] No files deleted. Run without --dry-run to delete.")
        return

    if not args.yes:
        confirm = input("\nProceed with deletion? [y/N]: ").strip().lower()
        if confirm != "y":
            print("Aborted.")
            return

    deleted = 0
    for p, typ in items:
        try:
            if p.is_dir():
                import shutil
                shutil.rmtree(p)
            else:
                p.unlink()
            deleted += 1
            print(f"  Deleted: {p.relative_to(ROOT)}")
        except Exception as e:
            print(f"  ERROR: {p}: {e}")

    print(f"\n[DONE] Deleted {deleted} items.")


if __name__ == "__main__":
    main()
