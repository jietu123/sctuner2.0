"""检查真值文件中CD8/NK/Endo在ST中的实际分布"""
import pandas as pd
import sys

truth_path = sys.argv[1] if len(sys.argv) > 1 else "data/sim/real_brca/S0/t_cells_cd4/seed_42/sim_truth_spot_type_fraction.csv"

truth = pd.read_csv(truth_path, index_col=0)

# 检查三个类型
types_to_check = ["T cells CD8", "NK cells", "Endothelial cells"]
missing_type = "T cells CD4"

print(f"Checking truth file: {truth_path}")
print(f"Missing type (intentionally removed): {missing_type}\n")

for t in types_to_check:
    if t in truth.columns:
        col = truth[t]
        total_fraction = col.sum()
        n_spots_with_type = (col > 0).sum()
        n_spots_total = len(col)
        mean_fraction = col.mean()
        max_fraction = col.max()
        
        print(f"{t}:")
        print(f"  Total fraction across all spots: {total_fraction:.4f}")
        print(f"  Spots with this type: {n_spots_with_type} / {n_spots_total} ({n_spots_with_type/n_spots_total*100:.2f}%)")
        print(f"  Mean fraction (when present): {mean_fraction:.4f}")
        print(f"  Max fraction in a spot: {max_fraction:.4f}")
        print()
    else:
        print(f"{t}: NOT FOUND in truth columns")
        print()

# 检查missing_type是否真的不在真值中
if missing_type in truth.columns:
    col = truth[missing_type]
    total_fraction = col.sum()
    print(f"{missing_type} (should be missing):")
    print(f"  Total fraction: {total_fraction:.4f} (should be 0)")
    print(f"  Spots with this type: {(col > 0).sum()} (should be 0)")
else:
    print(f"{missing_type} (should be missing): NOT FOUND in truth columns (correct)")

