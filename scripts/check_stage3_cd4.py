"""快速检查Stage3中CD4的support_score和support_category"""
import pandas as pd
import sys

sample = sys.argv[1] if len(sys.argv) > 1 else "real_brca_simS0_mt_t_cells_cd4_seed42"
ts_path = f"data/processed/{sample}/stage3_typematch/type_support.csv"
relabel_path = f"data/processed/{sample}/stage3_typematch/cell_type_relabel.csv"

ts = pd.read_csv(ts_path)
cd4_row = ts[ts.orig_type == "T cells CD4"]
if len(cd4_row) > 0:
    print(f"CD4 support_score: {cd4_row.iloc[0]['support_score']:.4f}")
    print(f"CD4 support_category: {cd4_row.iloc[0]['support_category']}")
else:
    print("CD4 not found in type_support.csv")

relabel = pd.read_csv(relabel_path)
unknown_count = (relabel["plugin_type"] == "Unknown_sc_only").sum()
total = len(relabel)
print(f"Unknown_sc_only cells: {unknown_count} / {total} ({unknown_count/total*100:.2f}%)")

