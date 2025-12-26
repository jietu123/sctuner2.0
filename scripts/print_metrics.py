"""快速输出 Stage5 指标"""
import json
from pathlib import Path

json_path = Path("result/real_brca_simS0_seed42/stage5_route2_s0/stage5_route2_s0__route2.json")
if not json_path.exists():
    print(f"[ERROR] File not found: {json_path}")
    exit(1)

data = json.loads(json_path.read_text(encoding="utf-8"))
comp = data.get("composition", {})
leak = data.get("leakage", {})
cell = data.get("cell_spot_eval", {}).get("metrics", {})

print("=" * 80)
print("V4.1 优化后 (support_threshold=2.00) - Route2 指标")
print("=" * 80)
print("\n=== 组成指标 ===")
print(f"L1_mean: {comp.get('L1_mean', 0):.4f}")
print(f"JS_mean: {comp.get('JS_mean', 0):.4f}")
print(f"corr_mean: {comp.get('corr_mean', 0):.4f}")

print("\n=== 泄漏指标 ===")
print(f"missing_leak_rate: {leak.get('missing_leak_rate', 0):.4f}")
print(f"missing_in_assignment: {leak.get('missing_in_assignment', 0)}")

print("\n=== 细胞级指标 ===")
print(f"coverage_on_truth: {cell.get('coverage_on_truth', 0):.4f}")
print(f"extra_pred_fraction: {cell.get('extra_pred_fraction', 0):.4f}")
print(f"acc_top1_cond_nonmissing: {cell.get('acc_top1_cond_nonmissing', 0):.4f}")
print(f"acc_top1_all_nonmissing: {cell.get('acc_top1_all_nonmissing', 0):.4f}")

