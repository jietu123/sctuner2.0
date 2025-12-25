"""小任务3：自动调整weak_th直到Route2清零泄漏

固定流程：
1. weak_th = 0.50：跑Stage4+Stage5(route2)
   - 若Route2 missing_present_in_assignment == 0：直接定为0.50（收工）
2. 若不为0：按步长上调（0.02）
   - 0.52 → 0.54 → 0.56 → 0.58 → ...
   每次只看2件事：
   - Route2 missing_present_in_assignment是否清零
   - filter_audit.missing_type_contribution.count是否已>0（说明missing_type已开始被过滤/生效）
3. 一旦Route2清零立刻停，记录该阈值为weak_th*，同时记录代价三件套

Usage:
  python scripts/auto_tune_weak_th.py \
    --sample real_brca_simS0_mt_b_cells_seed42 \
    --missing_type "B cells" \
    --seed 42 \
    --max_weak_th 0.70
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

# 设置输出编码
if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')


def detect_root() -> Path:
    here = Path(__file__).resolve()
    try:
        return here.parents[1]
    except IndexError:
        return here.parent


def load_yaml(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        import yaml
        with path.open("r", encoding="utf-8") as f:
            return yaml.safe_load(f) or {}
    except Exception:
        return {}


def save_yaml(path: Path, data: dict):
    try:
        import yaml
        with path.open("w", encoding="utf-8") as f:
            yaml.dump(data, f, default_flow_style=False, allow_unicode=True)
    except Exception:
        pass


def run_stage3(sample: str, weak_th: float, strong_th: float = 0.80):
    """运行Stage3"""
    root = detect_root()
    cfg_path = root / "configs" / "datasets" / f"{sample}.yaml"
    
    # 更新配置
    cfg = load_yaml(cfg_path)
    if "stage3" not in cfg:
        cfg["stage3"] = {}
    cfg["stage3"]["weak_th"] = weak_th
    cfg["stage3"]["strong_th"] = strong_th
    cfg["stage3"]["abstain_unknown_sc_only"] = True
    save_yaml(cfg_path, cfg)
    
    # 运行Stage3
    cmd = [
        "conda", "run", "-n", "cytospace_v1.1.0_py310", "python",
        "src/stages/stage3_type_plugin.py",
        "--sample", sample,
        "--dataset_config", str(cfg_path)
    ]
    result = subprocess.run(cmd, cwd=root, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[ERROR] Stage3 failed with weak_th={weak_th}")
        print(result.stderr)
        return False
    return True


def run_stage4_baseline(sample: str, missing_type: str, seed: int = 42):
    """运行Stage4 Baseline（官方CytoSPACE，无过滤）"""
    root = detect_root()
    cmd = [
        "conda", "run", "-n", "cytospace_v1.1.0_py310", "python",
        "src/stages/stage4_cytospace.py",
        "--sample", sample,
        "--filter_mode", "none",
        "--cell_type_column", "sc_meta",
        "--missing_type", missing_type,
        "--mean_cell_numbers", "5",
        "--solver_method", "lap_CSPR",
        "--seed", str(seed),
        "--sampling_sub_spots",
        "--n_subspots", "800",
    ]
    result = subprocess.run(cmd, cwd=root, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[ERROR] Stage4 Baseline failed")
        print(result.stderr)
        return False
    return True


def run_stage4_route2(sample: str, missing_type: str, seed: int = 42, filter_scope: str = "unsupported_all"):
    """运行Stage4 Route2"""
    root = detect_root()
    cmd = [
        "conda", "run", "-n", "cytospace_v1.1.0_py310", "python",
        "src/stages/stage4_cytospace.py",
        "--sample", sample,
        "--filter_mode", "plugin_unknown",
        "--cell_type_column", "plugin_type",
        "--missing_type", missing_type,
        "--mean_cell_numbers", "5",
        "--solver_method", "lap_CSPR",
        "--seed", str(seed),
        "--sampling_sub_spots",
        "--n_subspots", "800",
        "--filter_scope", filter_scope,
    ]
    result = subprocess.run(cmd, cwd=root, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[ERROR] Stage4 Route2 failed")
        print(result.stderr)
        return False
    return True


def run_stage5_baseline(sample: str, sim_dir: Path, out_dir: Path):
    """运行Stage5 Baseline"""
    root = detect_root()
    stage4_dir = root / "result" / sample / "stage4_cytospace" / "cytospace_output"
    cmd = [
        "conda", "run", "-n", "cytospace_v1.1.0_py310", "python",
        "src/stages/stage5_route2_s0.py",
        "--sample", sample,
        "--run_tag", "baseline",
        "--stage4_dir", str(stage4_dir),
        "--sim_dir", str(sim_dir),
        "--out_dir", str(out_dir),
        "--strict_config",
    ]
    result = subprocess.run(cmd, cwd=root, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[ERROR] Stage5 Baseline failed")
        print(result.stderr)
        return None
    return json.loads((out_dir / "stage5_route2_s0__baseline.json").read_text(encoding="utf-8"))


def run_stage5_route2(sample: str, sim_dir: Path, out_dir: Path):
    """运行Stage5 Route2"""
    root = detect_root()
    stage4_dir = root / "result" / sample / "stage4_cytospace" / "cytospace_output"
    cmd = [
        "conda", "run", "-n", "cytospace_v1.1.0_py310", "python",
        "src/stages/stage5_route2_s0.py",
        "--sample", sample,
        "--run_tag", "route2",
        "--stage4_dir", str(stage4_dir),
        "--sim_dir", str(sim_dir),
        "--out_dir", str(out_dir),
        "--strict_config",
    ]
    result = subprocess.run(cmd, cwd=root, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[ERROR] Stage5 Route2 failed")
        print(result.stderr)
        return None
    return json.loads((out_dir / "stage5_route2_s0__route2.json").read_text(encoding="utf-8"))


def check_route2_cleared(stage5_json: dict) -> tuple[bool, int]:
    """检查Route2是否清零泄漏"""
    leakage = stage5_json.get("leakage", {})
    missing_in_assignment = leakage.get("missing_in_assignment", -1)
    return missing_in_assignment == 0, missing_in_assignment


def get_filter_audit_summary(stage5_json: dict) -> dict:
    """获取filter_audit摘要"""
    audit = stage5_json.get("filter_audit", {})
    if not audit:
        return {}
    
    return {
        "missing_type_contrib_count": audit.get("missing_type_contribution", {}).get("count", 0),
        "non_missing_filtered_fraction": audit.get("non_missing_filtered", {}).get("fraction_of_total", 0.0),
        "top3_filtered_types": audit.get("top10_filtered_types", [])[:3],
        "total_filtered_fraction": audit.get("total_filtered", 0) / audit.get("total_sc_cells", 1) if audit.get("total_sc_cells", 0) > 0 else 0.0
    }


def main():
    parser = argparse.ArgumentParser(description="自动调整weak_th直到Route2清零泄漏")
    parser.add_argument("--sample", required=True, help="样本名（如real_brca_simS0_mt_b_cells_seed42）")
    parser.add_argument("--missing_type", required=True, help="缺失类型（如'B cells'）")
    parser.add_argument("--seed", type=int, default=42, help="随机种子")
    parser.add_argument("--max_weak_th", type=float, default=0.70, help="最大weak_th（防止无限上调）")
    parser.add_argument("--step", type=float, default=0.02, help="上调步长")
    parser.add_argument("--sim_dir", help="SimGen输出目录（如data/sim/real_brca/S0/b_cells/seed_42）")
    parser.add_argument(
        "--filter_scope",
        choices=["unsupported_all", "missing_only"],
        default="missing_only",
        help="Route2过滤范围：missing_only=精准止漏（推荐，默认）；unsupported_all=全局过滤（上界参考）",
    )
    parser.add_argument(
        "--max_non_missing_filtered_fraction",
        type=float,
        default=0.30,
        help="允许的最大非missing误伤比例（相对于SC总量），默认0.30；超出则视为代价不可接受",
    )
    parser.add_argument(
        "--max_total_filtered_fraction",
        type=float,
        default=0.50,
        help="允许的最大总过滤比例（相对于SC总量），默认0.50；超出则视为代价不可接受",
    )
    args = parser.parse_args()
    
    root = detect_root()
    
    # 确定sim_dir
    if args.sim_dir:
        sim_dir = Path(args.sim_dir)
    else:
        # 从sample推断
        mt_slug = args.missing_type.lower().replace(" ", "_")
        sim_dir = root / "data" / "sim" / "real_brca" / "S0" / mt_slug / f"seed_{args.seed}"
    
    out_dir = root / "result" / args.sample / "stage5_route2_s0"
    out_dir.mkdir(parents=True, exist_ok=True)
    
    print("="*80)
    print(f"自动调整weak_th直到Route2清零泄漏")
    print("="*80)
    print(f"Sample: {args.sample}")
    print(f"Missing type: {args.missing_type}")
    print(f"Seed: {args.seed}")
    print(f"Sim dir: {sim_dir}")
    print(f"Max weak_th: {args.max_weak_th}")
    print(f"Step: {args.step}")
    print(f"Filter_scope: {args.filter_scope} (推荐: missing_only)")
    print(f"Max non-missing filtered fraction: {args.max_non_missing_filtered_fraction:.3f}")
    print(f"Max total filtered fraction: {args.max_total_filtered_fraction:.3f}")
    print()
    
    # 首先运行 baseline（只运行一次，在开始 tuning 之前）
    print("[Baseline] Running baseline (official CytoSPACE, no filtering)...")
    print("  Running Stage4 Baseline...")
    if not run_stage4_baseline(args.sample, args.missing_type, args.seed):
        print("[ERROR] Baseline Stage4 failed, but continuing with Route2 tuning...")
    else:
        print("  Running Stage5 Baseline...")
        baseline_json = run_stage5_baseline(args.sample, sim_dir, out_dir)
        if baseline_json is None:
            print("[WARN] Baseline Stage5 failed, but continuing with Route2 tuning...")
        else:
            baseline_leak = baseline_json.get("leakage", {}).get("missing_leak_rate", -1)
            print(f"  [OK] Baseline leak_rate: {baseline_leak:.4f}")
    print()
    
    weak_th = 0.50
    weak_th_star = None
    results = []
    
    while weak_th <= args.max_weak_th:
        print(f"[Step] Testing weak_th = {weak_th:.2f}")
        
        # 运行Stage3
        print(f"  Running Stage3...")
        if not run_stage3(args.sample, weak_th):
            print(f"  [FAIL] Stage3 failed, skipping")
            weak_th += args.step
            continue
        
        # 运行Stage4 Route2
        print(f"  Running Stage4 Route2...")
        if not run_stage4_route2(args.sample, args.missing_type, args.seed, filter_scope=args.filter_scope):
            print(f"  [FAIL] Stage4 Route2 failed, skipping")
            weak_th += args.step
            continue
        
        # 运行Stage5 Route2
        print(f"  Running Stage5 Route2...")
        stage5_json = run_stage5_route2(args.sample, sim_dir, out_dir)
        if stage5_json is None:
            print(f"  [FAIL] Stage5 Route2 failed, skipping")
            weak_th += args.step
            continue
        
        # 检查是否清零
        cleared, missing_count = check_route2_cleared(stage5_json)
        audit_summary = get_filter_audit_summary(stage5_json)
        missing_contrib_count = audit_summary.get("missing_type_contrib_count", 0)

        total_filtered_fraction = audit_summary.get("total_filtered_fraction", 0.0)
        non_missing_filtered_fraction = audit_summary.get("non_missing_filtered_fraction", 0.0)

        # 约束判断：是否在可接受误伤范围内
        over_non_missing = (
            args.max_non_missing_filtered_fraction is not None
            and non_missing_filtered_fraction > args.max_non_missing_filtered_fraction
        )
        over_total = (
            args.max_total_filtered_fraction is not None
            and total_filtered_fraction > args.max_total_filtered_fraction
        )
        acceptable = cleared and not (over_non_missing or over_total)

        print(f"  Result: missing_present_in_assignment   = {missing_count}")
        print(f"          missing_type_contrib_count     = {missing_contrib_count}")
        print(f"          total_filtered_fraction        = {total_filtered_fraction:.4f}")
        print(f"          non_missing_filtered_fraction = {non_missing_filtered_fraction:.4f}")
        if cleared and not acceptable:
            print("  [WARN] 已清零，但误伤超过可接受上限（建议调整filter_scope或双阈值策略）")
        elif cleared and acceptable:
            print("  [OK] 已清零且误伤在可接受范围内")

        stop_reason = ""
        # unsupported_all 下，清零但超限时提前停止，建议切换策略
        if args.filter_scope == "unsupported_all" and cleared and (over_non_missing or over_total):
            stop_reason = "cleared_but_exceeds_limits_early_stop_unsupported_all"
        # missing_only 下，若清零但 total_filtered_fraction 已超过上限，也可以提前停止
        elif args.filter_scope == "missing_only" and cleared and over_total:
            stop_reason = "cleared_but_exceeds_total_limit_early_stop_missing_only"
        elif acceptable:
            stop_reason = "cleared_within_limits"

        results.append({
            "weak_th": weak_th,
            "missing_present_in_assignment": missing_count,
            "missing_type_contrib_count": missing_contrib_count,
            "cleared": cleared,
            "acceptable": acceptable,
            "over_total_limit": over_total,
            "over_non_missing_limit": over_non_missing,
            "stop_reason": stop_reason,
            **audit_summary
        })

        if stop_reason == "cleared_but_exceeds_limits_early_stop_unsupported_all":
            print("  [EARLY-STOP] 在 unsupported_all 下已清零但误伤爆炸，继续上调 weak_th 几乎只会更糟。")
            print("              建议改用 --filter_scope missing_only 或引入 weak_th_missing / weak_th_other 双阈值策略。")
            break
        if stop_reason == "cleared_but_exceeds_total_limit_early_stop_missing_only":
            print("  [EARLY-STOP] 在 missing_only 下已清零但总体过滤比例超过上限，继续上调 weak_th 只会过滤更多细胞。")
            print("              建议放宽 max_total_filtered_fraction 或改用更温和的配置。")
            break

        if acceptable and weak_th_star is None:
            weak_th_star = weak_th
            print(f"  [SUCCESS] Route2 cleared at weak_th = {weak_th:.2f} (within injury budget)")
            break
        
        # 检查missing_type是否已开始被过滤
        if missing_contrib_count == 0:
            print(f"  [WARNING] missing_type not yet filtered, continuing...")
        
        weak_th += args.step
    
    # 输出最终总结
    print()
    print("="*80)
    print("最终总结")
    print("="*80)
    
    if weak_th_star is None:
        # 尝试区分两种失败模式：
        any_cleared = any(r.get("cleared") for r in results)
        if any_cleared:
            print(f"[FAIL] 存在可以清零的weak_th，但都超出了误伤上限（已尝试到{weak_th-args.step:.2f}）")
            print("       建议：使用 --filter_scope missing_only 或引入 weak_th_missing / weak_th_other 双阈值策略。")
        else:
            print(f"[FAIL] 未能找到使Route2清零的weak_th（已尝试到{weak_th-args.step:.2f}）")
        # 仍然照常写出summary CSV
    else:
        # 选取第一个 acceptable 结果作为 final_result
        final_result = next(r for r in results if r.get("acceptable"))
        print(f"Missing type: {args.missing_type}")
        print(f"weak_th*: {weak_th_star:.2f}")
        print(f"Route2_missing_present: {final_result['missing_present_in_assignment']} (should be 0)")
        print(f"total_filtered_fraction: {final_result['total_filtered_fraction']*100:.2f}%")
        print(f"non_missing_filtered_fraction: {final_result['non_missing_filtered_fraction']*100:.2f}%")
        print(f"top3_filtered_types:")
        for item in final_result.get("top3_filtered_types", []):
            print(f"  - {item['type']}: {item['count']} ({item['fraction']*100:.2f}%)")
    
    # 保存总结到CSV
    summary_path = out_dir / "weak_th_tuning_summary.csv"
    import pandas as pd
    df = pd.DataFrame(results)
    df.to_csv(summary_path, index=False, encoding="utf-8")
    print(f"\n[OK] Summary saved to: {summary_path}")


if __name__ == "__main__":
    main()

