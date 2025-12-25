"""生成类型过滤审计报告"""
import pandas as pd
import json
from pathlib import Path

sample = "real_brca_simS0_mt_t_cells_cd4_seed42"
missing_type = "T cells CD4"

# 读取真值
truth_path = Path(f"data/sim/real_brca/S0/t_cells_cd4/seed_42/sim_truth_spot_type_fraction.csv")
truth = pd.read_csv(truth_path, index_col=0)

# 读取Stage3支持度
ts_path = Path(f"data/processed/{sample}/stage3_typematch/type_support.csv")
ts = pd.read_csv(ts_path)

# 读取filter_audit
audit_path = Path(f"result/{sample}/stage4_cytospace/filter_audit.json")
with audit_path.open("r", encoding="utf-8") as f:
    audit = json.load(f)

types_to_check = ["T cells CD8", "NK cells", "Endothelial cells"]

print("="*80)
print("类型过滤审计报告")
print("="*80)
print(f"\n样本: {sample}")
print(f"Missing type (应该被过滤): {missing_type}")
print(f"当前 weak_th: 0.58\n")

print("-"*80)
print("检查结果：")
print("-"*80)

for t in types_to_check:
    print(f"\n【{t}】")
    
    # 检查真值
    if t in truth.columns:
        col = truth[t]
        n_spots_with = (col > 0).sum()
        n_spots_total = len(col)
        mean_fraction = col[col > 0].mean() if (col > 0).sum() > 0 else 0
        print(f"  [EXISTS] 真值中存在: {n_spots_with}/{n_spots_total} spots ({n_spots_with/n_spots_total*100:.2f}%)")
        print(f"    平均fraction (当存在时): {mean_fraction:.4f}")
    else:
        print(f"  [NOT_EXISTS] 真值中不存在")
    
    # 检查Stage3支持度
    ts_row = ts[ts.orig_type == t]
    if len(ts_row) > 0:
        support_score = ts_row.iloc[0]['support_score']
        support_category = ts_row.iloc[0]['support_category']
        print(f"  Stage3 support_score: {support_score:.4f}")
        print(f"  Stage3 support_category: {support_category}")
        
        if support_category == "unsupported":
            print(f"  [WARNING] 被判定为 unsupported (会被过滤)")
        else:
            print(f"  [OK] 未被判定为 unsupported")
    
    # 检查是否被过滤
    filtered_count = audit["all_filtered_types"].get(t, 0)
    if filtered_count > 0:
        print(f"  [FILTERED] 被过滤: {filtered_count} cells")
        if t in truth.columns and (col > 0).sum() > 0:
            print(f"  [ERROR] 错误过滤：该类型在ST真值中存在，但被Stage3判为unsupported并过滤")
        else:
            print(f"  [OK] 正确过滤：该类型在ST真值中不存在或极少")
    else:
        print(f"  [OK] 未被过滤")

print("\n" + "="*80)
print("总结：")
print("="*80)

# 统计
correctly_filtered = []
incorrectly_filtered = []

for t in types_to_check:
    if t in truth.columns:
        col = truth[t]
        n_spots_with = (col > 0).sum()
        filtered_count = audit["all_filtered_types"].get(t, 0)
        
        if n_spots_with > 0 and filtered_count > 0:
            incorrectly_filtered.append((t, filtered_count, n_spots_with))
        elif n_spots_with == 0 and filtered_count > 0:
            correctly_filtered.append((t, filtered_count))

print(f"\n错误过滤的类型（在ST真值中存在但被过滤）:")
if incorrectly_filtered:
    for t, count, n_spots in incorrectly_filtered:
        print(f"  - {t}: {count} cells (在{n_spots}个spots中存在)")
else:
    print("  无")

print(f"\n正确过滤的类型（在ST真值中不存在）:")
if correctly_filtered:
    for t, count in correctly_filtered:
        print(f"  - {t}: {count} cells")
else:
    print("  无")

print(f"\nMissing type ({missing_type}):")
if missing_type in truth.columns:
    col = truth[missing_type]
    n_spots_with = (col > 0).sum()
    print(f"  [WARNING] missing_type在真值中存在 {n_spots_with} spots")
else:
    print(f"  [OK] 正确：missing_type不在真值中（符合S0设计）")

