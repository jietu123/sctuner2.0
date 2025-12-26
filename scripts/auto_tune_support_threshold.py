"""自动调整 V4.1 support_threshold 以优化检测精度

目标：
1. 确保 missing_type 被正确检测（FlagMismatch=Yes）
2. 最小化非 missing 类型的误标记（false positives）
3. 优化组成指标（L1, JS, corr）

策略：
- 从较高的阈值开始（如 2.5），逐步下调（步长 0.1）
- 对于每个阈值，评估：
  - missing_type 是否被检测到
  - 非 missing 类型的误标记数量
  - 组成指标（L1, JS, corr）
- 选择在检测到 missing_type 的同时，误标记最少的阈值

Usage:
  python scripts/auto_tune_support_threshold.py \
    --sample real_brca_simS0_seed42 \
    --missing_type "T cells CD8" \
    --seed 42 \
    --min_support_threshold 1.0 \
    --max_support_threshold 2.5
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


def run_stage3(sample: str, support_threshold: float, marker_gene_count: int = 30):
    """运行Stage3 V4.1"""
    root = detect_root()
    cfg_path = root / "configs" / "datasets" / f"{sample}.yaml"
    
    # 更新配置
    cfg = load_yaml(cfg_path)
    if "stage3" not in cfg:
        cfg["stage3"] = {}
    cfg["stage3"]["enable_mismatch_detection"] = True
    cfg["stage3"]["marker_gene_count"] = marker_gene_count
    cfg["stage3"]["score_method"] = "fisher_z"
    cfg["stage3"]["support_threshold"] = support_threshold
    # 保持其他配置不变
    if "strong_th" not in cfg["stage3"]:
        cfg["stage3"]["strong_th"] = 0.8
    if "weak_th" not in cfg["stage3"]:
        cfg["stage3"]["weak_th"] = 0.5
    save_yaml(cfg_path, cfg)
    
    # 运行Stage3
    cmd = [
        "conda", "run", "-n", "cytospace_v1.1.0_py310", "python",
        "src/stages/stage3_type_plugin.py",
        "--sample", sample,
        "--dataset_config", str(cfg_path)
    ]
    result = subprocess.run(
        cmd, 
        cwd=root, 
        capture_output=True, 
        text=True,
        encoding='utf-8',
        errors='replace'
    )
    if result.returncode != 0:
        print(f"[ERROR] Stage3 failed with support_threshold={support_threshold}")
        if result.stderr:
            stderr_lines = result.stderr.split('\n')[:5]
            for line in stderr_lines:
                if line.strip():
                    print(f"  {line}")
        return False
    return True


def run_stage4_route2(sample: str, missing_type: str, seed: int = 42, filter_scope: str = "missing_only"):
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
    result = subprocess.run(
        cmd, 
        cwd=root, 
        capture_output=True, 
        text=True,
        encoding='utf-8',
        errors='replace'
    )
    if result.returncode != 0:
        print(f"[ERROR] Stage4 Route2 failed")
        if result.stderr:
            stderr_lines = result.stderr.split('\n')[:5]
            for line in stderr_lines:
                if line.strip():
                    print(f"  {line}")
        return False
    return True


def run_stage5_route2(sample: str, sim_dir: Path, out_dir: Path):
    """运行Stage5 Route2"""
    root = detect_root()
    stage4_dir = root / "result" / sample / "stage4_cytospace" / "cytospace_output"
    json_path = out_dir / "stage5_route2_s0__route2.json"
    
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
    # 使用 errors='replace' 来避免编码问题
    result = subprocess.run(
        cmd, 
        cwd=root, 
        capture_output=True, 
        text=True,
        encoding='utf-8',
        errors='replace'
    )
    
    # 即使返回码非 0，也尝试读取 JSON（可能已经生成了）
    if json_path.exists():
        try:
            return json.loads(json_path.read_text(encoding="utf-8"))
        except Exception as e:
            if result.returncode != 0:
                print(f"[WARN] Stage5 returned non-zero, but trying to read JSON anyway...")
                print(f"[ERROR] Failed to read JSON: {e}")
            return None
    
    if result.returncode != 0:
        print(f"[ERROR] Stage5 Route2 failed (returncode={result.returncode})")
        # 只打印错误的关键部分，避免编码问题
        if result.stderr:
            stderr_lines = result.stderr.split('\n')[:5]  # 只显示前5行
            for line in stderr_lines:
                if line.strip():
                    print(f"  {line}")
        return None
    
    return json.loads(json_path.read_text(encoding="utf-8"))


def check_missing_detected(sample: str, missing_type: str) -> tuple[bool, float]:
    """检查 missing_type 是否被检测到（FlagMismatch=Yes）"""
    root = detect_root()
    type_support_path = root / "data" / "processed" / sample / "stage3_typematch" / "type_support.csv"
    if not type_support_path.exists():
        return False, 0.0
    
    import pandas as pd
    df = pd.read_csv(type_support_path)
    missing_row = df[df["CellType"] == missing_type]
    if len(missing_row) == 0:
        return False, 0.0
    
    flag_mismatch = missing_row.iloc[0]["FlagMismatch"]
    support_score = missing_row.iloc[0]["SupportScore"]
    return flag_mismatch == "Yes", support_score


def count_false_positives(sample: str, missing_type: str) -> int:
    """统计非 missing 类型被误标记为 FlagMismatch=Yes 的数量"""
    root = detect_root()
    type_support_path = root / "data" / "processed" / sample / "stage3_typematch" / "type_support.csv"
    if not type_support_path.exists():
        return 0
    
    import pandas as pd
    df = pd.read_csv(type_support_path)
    # 排除 missing_type 本身
    false_positives = df[(df["CellType"] != missing_type) & (df["FlagMismatch"] == "Yes")]
    return len(false_positives)


def get_composition_metrics(stage5_json: dict) -> dict:
    """获取组成指标"""
    comp = stage5_json.get("composition", {})
    return {
        "L1_mean": comp.get("L1_mean", float('inf')),
        "JS_mean": comp.get("JS_mean", float('inf')),
        "corr_mean": comp.get("corr_mean", -float('inf')),
    }


def get_filter_audit_summary(stage5_json: dict) -> dict:
    """获取filter_audit摘要"""
    audit = stage5_json.get("filter_audit", {})
    if not audit:
        return {}
    
    return {
        "missing_type_contrib_count": audit.get("missing_type_contribution", {}).get("count", 0),
        "non_missing_filtered_fraction": audit.get("non_missing_filtered", {}).get("fraction_of_total", 0.0),
        "marked_unknown_non_missing": audit.get("marked_unknown", {}).get("non_missing", 0),
        "total_filtered_fraction": audit.get("total_filtered", 0) / audit.get("total_sc_cells", 1) if audit.get("total_sc_cells", 0) > 0 else 0.0
    }


def main():
    parser = argparse.ArgumentParser(description="自动调整 V4.1 support_threshold 以优化检测精度")
    parser.add_argument("--sample", required=True, help="样本名（如real_brca_simS0_seed42）")
    parser.add_argument("--missing_type", required=True, help="缺失类型（如'T cells CD8'）")
    parser.add_argument("--seed", type=int, default=42, help="随机种子")
    parser.add_argument("--min_support_threshold", type=float, default=1.0, help="最小 support_threshold")
    parser.add_argument("--max_support_threshold", type=float, default=2.5, help="最大 support_threshold")
    parser.add_argument("--step", type=float, default=0.1, help="下调步长")
    parser.add_argument("--sim_dir", help="SimGen输出目录")
    parser.add_argument("--marker_gene_count", type=int, default=30, help="标记基因数量")
    parser.add_argument(
        "--max_non_missing_filtered_fraction",
        type=float,
        default=0.20,
        help="允许的最大非missing误标记比例（相对于SC总量），默认0.20",
    )
    args = parser.parse_args()
    
    root = detect_root()
    
    # 确定sim_dir
    if args.sim_dir:
        sim_dir = Path(args.sim_dir)
    else:
        # 从sample推断
        sim_dir = root / "data" / "sim" / "real_brca" / "S0" / f"seed_{args.seed}"
    
    out_dir = root / "result" / args.sample / "stage5_route2_s0"
    out_dir.mkdir(parents=True, exist_ok=True)
    
    print("="*80)
    print(f"自动调整 V4.1 support_threshold 以优化检测精度")
    print("="*80)
    print(f"Sample: {args.sample}")
    print(f"Missing type: {args.missing_type}")
    print(f"Seed: {args.seed}")
    print(f"Sim dir: {sim_dir}")
    print(f"Threshold range: [{args.min_support_threshold}, {args.max_support_threshold}]")
    print(f"Step: {args.step}")
    print(f"Max non-missing filtered fraction: {args.max_non_missing_filtered_fraction:.3f}")
    print()
    
    support_threshold = args.max_support_threshold
    best_threshold = None
    best_score = float('inf')
    results = []
    
    total_steps = int((args.max_support_threshold - args.min_support_threshold) / args.step) + 1
    current_step = 0
    
    while support_threshold >= args.min_support_threshold:
        current_step += 1
        print(f"[Step {current_step}/{total_steps}] Testing support_threshold = {support_threshold:.2f}")
        
        # 运行Stage3
        print(f"  Running Stage3...")
        if not run_stage3(args.sample, support_threshold, args.marker_gene_count):
            print(f"  [FAIL] Stage3 failed, skipping")
            support_threshold -= args.step
            continue
        
        # 检查 missing_type 是否被检测到
        missing_detected, missing_support_score = check_missing_detected(args.sample, args.missing_type)
        n_false_positives = count_false_positives(args.sample, args.missing_type)
        
        print(f"  Missing type detection: detected={missing_detected}, support_score={missing_support_score:.4f}, false_positives={n_false_positives}")
        
        if not missing_detected:
            print(f"  [SKIP] Missing type not detected (support_score={missing_support_score:.4f} >= threshold={support_threshold:.2f}), skipping pipeline")
            results.append({
                "support_threshold": support_threshold,
                "missing_detected": False,
                "missing_support_score": missing_support_score,
                "n_false_positives": n_false_positives,
                "skipped": True,
            })
            support_threshold -= args.step
            continue
        
        # 运行Stage4 Route2
        print(f"  Running Stage4 Route2...")
        if not run_stage4_route2(args.sample, args.missing_type, args.seed):
            print(f"  [FAIL] Stage4 Route2 failed, skipping")
            support_threshold -= args.step
            continue
        
        # 运行Stage5 Route2
        print(f"  Running Stage5 Route2...")
        stage5_json = run_stage5_route2(args.sample, sim_dir, out_dir)
        if stage5_json is None:
            print(f"  [FAIL] Stage5 Route2 failed, skipping")
            support_threshold -= args.step
            continue
        
        # 获取指标
        metrics = get_composition_metrics(stage5_json)
        audit_summary = get_filter_audit_summary(stage5_json)
        non_missing_filtered_fraction = audit_summary.get("non_missing_filtered_fraction", 0.0)
        
        # 检查是否在可接受范围内
        acceptable = non_missing_filtered_fraction <= args.max_non_missing_filtered_fraction
        
        # 计算综合评分（越小越好）
        # 评分 = L1_mean + JS_mean - corr_mean + 10 * non_missing_filtered_fraction + false_positives 惩罚
        # 权重可以根据需要调整
        # 对 false_positives 进行惩罚：每个 false positive 惩罚 0.5 分（增加权重以减少误标记）
        false_positive_penalty = n_false_positives * 0.5
        score = (
            metrics["L1_mean"] +
            metrics["JS_mean"] * 2 +  # JS 权重稍高
            (1.0 - metrics["corr_mean"]) * 3 +  # corr 越高越好，增加权重以更重视 corr_mean
            non_missing_filtered_fraction * 10 +  # 误标记惩罚
            false_positive_penalty  # false positives 惩罚（增加权重）
        )
        
        print(f"  Result: missing_detected          = {missing_detected}")
        print(f"          missing_support_score     = {missing_support_score:.4f}")
        print(f"          n_false_positives        = {n_false_positives} (penalty: {false_positive_penalty:.2f})")
        print(f"          L1_mean                  = {metrics['L1_mean']:.4f}")
        print(f"          JS_mean                  = {metrics['JS_mean']:.4f}")
        print(f"          corr_mean                = {metrics['corr_mean']:.4f}")
        print(f"          non_missing_filtered_frac = {non_missing_filtered_fraction:.4f}")
        print(f"          score                    = {score:.4f}")
        if acceptable:
            print("  [OK] Within acceptable limits")
        else:
            print(f"  [WARN] Non-missing fraction exceeds limit ({args.max_non_missing_filtered_fraction:.3f})")
        
        results.append({
            "support_threshold": support_threshold,
            "missing_detected": missing_detected,
            "missing_support_score": missing_support_score,
            "n_false_positives": n_false_positives,
            "acceptable": acceptable,
            "score": score,
            "skipped": False,
            **metrics,
            **audit_summary
        })
        
        # 更新最佳阈值
        if acceptable and score < best_score:
            best_score = score
            best_threshold = support_threshold
            print(f"  [BEST] New best threshold: {support_threshold:.2f} (score={score:.4f})")
        
        support_threshold -= args.step
    
    # 输出最终总结
    print()
    print("="*80)
    print("最终总结")
    print("="*80)
    
    if best_threshold is None:
        # 检查是否有任何可接受的结果
        acceptable_results = [r for r in results if r.get("acceptable") and not r.get("skipped", False)]
        if acceptable_results:
            # 选择评分最低的
            best_result = min(acceptable_results, key=lambda x: x.get("score", float('inf')))
            best_threshold = best_result["support_threshold"]
            print(f"[WARN] 未找到完全可接受的结果，但选择了评分最低的阈值")
        else:
            print(f"[FAIL] 未能找到可接受的 support_threshold")
            # 仍然选择评分最低的
            valid_results = [r for r in results if not r.get("skipped", False)]
            if valid_results:
                best_result = min(valid_results, key=lambda x: x.get("score", float('inf')))
                best_threshold = best_result["support_threshold"]
                print(f"[INFO] 选择评分最低的阈值: {best_threshold:.2f}")
    
    if best_threshold is not None:
        best_result = next(r for r in results if r.get("support_threshold") == best_threshold and not r.get("skipped", False))
        print(f"Missing type: {args.missing_type}")
        print(f"Best support_threshold*: {best_threshold:.2f}")
        print(f"Missing detected: {best_result['missing_detected']}")
        print(f"Missing support_score: {best_result['missing_support_score']:.4f}")
        print(f"False positives: {best_result['n_false_positives']}")
        print(f"L1_mean: {best_result['L1_mean']:.4f}")
        print(f"JS_mean: {best_result['JS_mean']:.4f}")
        print(f"corr_mean: {best_result['corr_mean']:.4f}")
        print(f"non_missing_filtered_fraction: {best_result['non_missing_filtered_fraction']:.4f}")
        print(f"score: {best_result['score']:.4f}")
    
    # 保存总结到CSV
    summary_path = out_dir / "support_threshold_tuning_summary.csv"
    import pandas as pd
    df = pd.DataFrame(results)
    df.to_csv(summary_path, index=False, encoding="utf-8")
    print(f"\n[OK] Summary saved to: {summary_path}")
    
    # 更新配置文件
    if best_threshold is not None:
        cfg_path = root / "configs" / "datasets" / f"{args.sample}.yaml"
        cfg = load_yaml(cfg_path)
        if "stage3" not in cfg:
            cfg["stage3"] = {}
        cfg["stage3"]["support_threshold"] = best_threshold
        save_yaml(cfg_path, cfg)
        print(f"[OK] Updated config: {cfg_path} (support_threshold = {best_threshold:.2f})")


if __name__ == "__main__":
    main()
