"""检查ST数据中CD8/NK/Endo的表达情况"""
import pandas as pd
import numpy as np
import sys

st_expr_path = sys.argv[1] if len(sys.argv) > 1 else "data/processed/real_brca_simS0_mt_t_cells_cd4_seed42/stage1_preprocess/exported/st_expression_normalized.csv"
sc_expr_path = sys.argv[2] if len(sys.argv) > 2 else "data/processed/real_brca_simS0_mt_t_cells_cd4_seed42/stage1_preprocess/exported/sc_expression_normalized.csv"
sc_meta_path = sys.argv[3] if len(sys.argv) > 3 else "data/processed/real_brca_simS0_mt_t_cells_cd4_seed42/stage1_preprocess/exported/sc_metadata.csv"

print("Loading data...")
st_expr = pd.read_csv(st_expr_path, index_col=0, sep=None, engine="python")
sc_expr = pd.read_csv(sc_expr_path, index_col=0, sep=None, engine="python")
sc_meta = pd.read_csv(sc_meta_path, sep=None, engine="python")

# 获取cell_type列
cell_type_col = sc_meta.columns[1] if len(sc_meta.columns) > 1 else "cell_type"
if "cell_type" in sc_meta.columns:
    cell_type_col = "cell_type"

types_to_check = ["T cells CD8", "NK cells", "Endothelial cells"]

print(f"\nST expression shape: {st_expr.shape}")
print(f"SC expression shape: {sc_expr.shape}")
print(f"SC metadata shape: {sc_meta.shape}")
print(f"Cell type column: {cell_type_col}\n")

# 检查每个类型
for t in types_to_check:
    print(f"{'='*60}")
    print(f"{t}:")
    
    # SC中该类型的细胞
    sc_cells_of_type = sc_meta[sc_meta[cell_type_col] == t]
    n_sc_cells = len(sc_cells_of_type)
    print(f"  SC cells of this type: {n_sc_cells}")
    
    if n_sc_cells > 0:
        # 计算SC中该类型的平均表达
        sc_cell_ids = sc_cells_of_type.iloc[:, 0].values  # 第一列是cell_id
        sc_cell_ids_in_expr = [cid for cid in sc_cell_ids if cid in sc_expr.index]
        if len(sc_cell_ids_in_expr) > 0:
            sc_type_expr = sc_expr.loc[sc_cell_ids_in_expr]
            sc_mean_expr = sc_type_expr.mean(axis=0)
            sc_top_genes = sc_mean_expr.nlargest(5).index.tolist()
            print(f"  SC top 5 genes (mean): {sc_top_genes}")
        
        # 计算ST中所有spots的平均表达
        st_mean_expr = st_expr.mean(axis=0)
        st_top_genes = st_mean_expr.nlargest(5).index.tolist()
        print(f"  ST top 5 genes (mean): {st_top_genes}")
        
        # 检查共同高表达基因
        common_genes = set(sc_top_genes) & set(st_top_genes)
        print(f"  Common top genes: {len(common_genes)} / 5")
    
    print()

