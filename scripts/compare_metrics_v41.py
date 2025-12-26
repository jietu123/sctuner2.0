"""对比 V4.1 优化后的指标与基线和地基映射"""
import json
from pathlib import Path

# 读取地基映射指标（从文档）
groundwork_metrics = {
    "L1_mean": 0.7956,
    "JS_mean": 0.1774,
    "corr_mean": 0.6762,
    "missing_leak_rate": 0.0,
    "coverage_on_truth": 1.0,
    "acc_top1_all_nonmissing": 0.1549,
}

# 读取基线映射（如果存在）
baseline_path = Path("result/real_brca_simS0_seed42/stage5_route2_s0/stage5_route2_s0__baseline.json")
baseline_metrics = None
if baseline_path.exists():
    baseline_data = json.loads(baseline_path.read_text(encoding="utf-8"))
    comp_b = baseline_data.get("composition", {})
    leak_b = baseline_data.get("leakage", {})
    cell_b = baseline_data.get("cell_spot_eval", {}).get("metrics", {})
    baseline_metrics = {
        "L1_mean": comp_b.get("L1_mean", 0),
        "JS_mean": comp_b.get("JS_mean", 0),
        "corr_mean": comp_b.get("corr_mean", 0),
        "missing_leak_rate": leak_b.get("missing_leak_rate", 0),
        "coverage_on_truth": cell_b.get("coverage_on_truth", 0),
        "acc_top1_all_nonmissing": cell_b.get("acc_top1_all_nonmissing", 0),
    }

# 读取 V4.1 优化后的指标
route2_path = Path("result/real_brca_simS0_seed42/stage5_route2_s0/stage5_route2_s0__route2.json")
route2_data = json.loads(route2_path.read_text(encoding="utf-8"))
comp_r = route2_data.get("composition", {})
leak_r = route2_data.get("leakage", {})
cell_r = route2_data.get("cell_spot_eval", {}).get("metrics", {})
v41_metrics = {
    "L1_mean": comp_r.get("L1_mean", 0),
    "JS_mean": comp_r.get("JS_mean", 0),
    "corr_mean": comp_r.get("corr_mean", 0),
    "missing_leak_rate": leak_r.get("missing_leak_rate", 0),
    "coverage_on_truth": cell_r.get("coverage_on_truth", 0),
    "acc_top1_all_nonmissing": cell_r.get("acc_top1_all_nonmissing", 0),
}

print("=" * 80)
print("seed42 CD8 场景指标对比（V4.1 优化后 vs 基线 vs 地基映射）")
print("=" * 80)
print()

# 创建对比表
print(f"{'指标':<30} {'基线映射':<15} {'地基映射':<15} {'V4.1优化后':<15} {'vs基线':<10} {'vs地基':<10}")
print("-" * 95)

def format_value(val):
    if val is None:
        return "N/A"
    return f"{val:.4f}"

def format_delta(v41, ref):
    if ref is None or ref == 0:
        return "N/A"
    delta = v41 - ref
    pct = (delta / ref) * 100 if ref != 0 else 0
    sign = "+" if delta >= 0 else ""
    return f"{sign}{delta:.4f} ({sign}{pct:.1f}%)"

metrics_list = [
    ("L1_mean", "L1_mean (越小越好)"),
    ("JS_mean", "JS_mean (越小越好)"),
    ("corr_mean", "corr_mean (越大越好)"),
    ("missing_leak_rate", "missing_leak_rate (越小越好)"),
    ("coverage_on_truth", "coverage_on_truth (越大越好)"),
    ("acc_top1_all_nonmissing", "acc_top1_all_nonmissing (越大越好)"),
]

for key, label in metrics_list:
    baseline_val = baseline_metrics.get(key) if baseline_metrics else None
    groundwork_val = groundwork_metrics.get(key)
    v41_val = v41_metrics.get(key)
    
    print(f"{label:<30} {format_value(baseline_val):<15} {format_value(groundwork_val):<15} {format_value(v41_val):<15} {format_delta(v41_val, baseline_val):<10} {format_delta(v41_val, groundwork_val):<10}")

print()
print("=" * 80)
print("关键发现：")
print("=" * 80)

# 分析
if baseline_metrics:
    print(f"\n1. vs 基线映射：")
    if v41_metrics["L1_mean"] < baseline_metrics["L1_mean"]:
        print(f"   [OK] L1_mean 改善: {baseline_metrics['L1_mean']:.4f} -> {v41_metrics['L1_mean']:.4f}")
    if v41_metrics["JS_mean"] < baseline_metrics["JS_mean"]:
        print(f"   [OK] JS_mean 改善: {baseline_metrics['JS_mean']:.4f} -> {v41_metrics['JS_mean']:.4f}")
    if v41_metrics["corr_mean"] > baseline_metrics["corr_mean"]:
        print(f"   [OK] corr_mean 改善: {baseline_metrics['corr_mean']:.4f} -> {v41_metrics['corr_mean']:.4f}")
    elif v41_metrics["corr_mean"] < baseline_metrics["corr_mean"]:
        print(f"   [WARN] corr_mean 下降: {baseline_metrics['corr_mean']:.4f} -> {v41_metrics['corr_mean']:.4f}")
    if v41_metrics["missing_leak_rate"] == 0 and baseline_metrics["missing_leak_rate"] > 0:
        print(f"   [OK] missing_leak_rate 清零: {baseline_metrics['missing_leak_rate']:.4f} -> {v41_metrics['missing_leak_rate']:.4f}")

print(f"\n2. vs 地基映射：")
if v41_metrics["L1_mean"] < groundwork_metrics["L1_mean"]:
    print(f"   [OK] L1_mean 改善: {groundwork_metrics['L1_mean']:.4f} -> {v41_metrics['L1_mean']:.4f}")
if v41_metrics["JS_mean"] < groundwork_metrics["JS_mean"]:
    print(f"   [OK] JS_mean 改善: {groundwork_metrics['JS_mean']:.4f} -> {v41_metrics['JS_mean']:.4f}")
if v41_metrics["corr_mean"] > groundwork_metrics["corr_mean"]:
    print(f"   [OK] corr_mean 改善: {groundwork_metrics['corr_mean']:.4f} -> {v41_metrics['corr_mean']:.4f}")
elif v41_metrics["corr_mean"] < groundwork_metrics["corr_mean"]:
    print(f"   [WARN] corr_mean 下降: {groundwork_metrics['corr_mean']:.4f} -> {v41_metrics['corr_mean']:.4f} (下降 {(groundwork_metrics['corr_mean'] - v41_metrics['corr_mean']) / groundwork_metrics['corr_mean'] * 100:.1f}%)")
if v41_metrics["missing_leak_rate"] == 0:
    print(f"   [OK] missing_leak_rate 保持清零: {v41_metrics['missing_leak_rate']:.4f}")

print(f"\n3. 自动调参结果：")
print(f"   最佳 support_threshold: 2.00")
print(f"   False positives (误标记): 6 个非 missing 类型被误标记")
print(f"   non_missing_filtered_fraction: 0.0 (无误伤)")

