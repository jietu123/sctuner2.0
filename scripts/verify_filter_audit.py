"""验证filter_audit的三个强制检查指标"""
import json
import sys

json_path = sys.argv[1] if len(sys.argv) > 1 else "result/real_brca_simS0_mt_t_cells_cd4_seed42/stage5_route2_s0/stage5_route2_s0__route2.json"

with open(json_path, "r", encoding="utf-8") as f:
    data = json.load(f)

audit = data["filter_audit"]["validation_checks"]

print("="*60)
print("Filter Audit Validation Checks")
print("="*60)
print()
print(f"Check1 (sum_top10 == total):")
print(f"  Passed: {audit['check1_sum_top10_equals_total']['passed']}")
print(f"  sum_top10: {audit['check1_sum_top10_equals_total']['sum_top10']}")
print(f"  n_filtered_total: {audit['check1_sum_top10_equals_total']['n_filtered_total']}")
print()
print(f"Check2 (missing_type_contrib):")
print(f"  Value: {audit['check2_missing_type_contrib']['value']*100:.2f}%")
print(f"  Description: {audit['check2_missing_type_contrib']['description']}")
print()
print(f"Check3 (non_missing_filtered_fraction):")
print(f"  Value: {audit['check3_non_missing_filtered_fraction']['value']*100:.2f}%")
print(f"  Description: {audit['check3_non_missing_filtered_fraction']['description']}")
print()
print("="*60)

