"""分析Stage4过滤的细胞类型构成"""
import pandas as pd
import json
import sys
from pathlib import Path

sample = sys.argv[1] if len(sys.argv) > 1 else "real_brca_simS0_mt_t_cells_cd4_seed42"
missing_type = sys.argv[2] if len(sys.argv) > 2 else "T cells CD4"

relabel_path = Path(f"data/processed/{sample}/stage3_typematch/cell_type_relabel.csv")
sim_info_path = Path(f"data/processed/{sample}/stage1_preprocess/exported/sim_info.json")

relabel = pd.read_csv(relabel_path)

# 找出所有被过滤的细胞（plugin_type == Unknown_sc_only）
filtered = relabel[relabel["plugin_type"] == "Unknown_sc_only"].copy()
total_filtered = len(filtered)
total_cells = len(relabel)

# 统计被过滤细胞的orig_type分布
type_counts = filtered["orig_type"].value_counts().to_dict()

# 读取sim_info获取missing_type
if sim_info_path.exists():
    with sim_info_path.open("r", encoding="utf-8") as f:
        sim_info = json.load(f)
    missing_type = sim_info.get("missing_type", missing_type)

# 计算missing_type的贡献
missing_count = type_counts.get(missing_type, 0)
missing_fraction = missing_count / total_filtered if total_filtered > 0 else 0.0

# 计算非missing类型的误伤
non_missing_filtered = total_filtered - missing_count
non_missing_fraction = non_missing_filtered / total_filtered if total_filtered > 0 else 0.0

# Top10被过滤类型
top10_types = list(type_counts.items())[:10]
top10_list = [{"type": t, "count": c, "fraction": c/total_filtered} for t, c in top10_types]

audit = {
    "sample": sample,
    "missing_type": missing_type,
    "total_cells": int(total_cells),
    "total_filtered": int(total_filtered),
    "filtered_fraction": float(total_filtered / total_cells) if total_cells > 0 else 0.0,
    "missing_type_contribution": {
        "count": int(missing_count),
        "fraction_of_filtered": float(missing_fraction),
        "fraction_of_total": float(missing_count / total_cells) if total_cells > 0 else 0.0
    },
    "non_missing_filtered": {
        "count": int(non_missing_filtered),
        "fraction_of_filtered": float(non_missing_fraction),
        "fraction_of_total": float(non_missing_filtered / total_cells) if total_cells > 0 else 0.0
    },
    "top10_filtered_types": top10_list,
    "all_filtered_types": {k: int(v) for k, v in type_counts.items()}
}

# 输出JSON
out_path = Path(f"result/{sample}/stage4_cytospace/filter_audit.json")
out_path.parent.mkdir(parents=True, exist_ok=True)
out_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8")

print(f"[OK] filter_audit written to: {out_path}")
print(f"\nSummary:")
print(f"  Total filtered: {total_filtered} / {total_cells} ({total_filtered/total_cells*100:.2f}%)")
print(f"  Missing type ({missing_type}): {missing_count} ({missing_fraction*100:.2f}% of filtered)")
print(f"  Non-missing filtered: {non_missing_filtered} ({non_missing_fraction*100:.2f}% of filtered)")
print(f"\nTop 10 filtered types:")
for item in top10_list:
    print(f"  {item['type']}: {item['count']} ({item['fraction']*100:.2f}%)")

