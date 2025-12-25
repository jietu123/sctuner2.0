"""为自动调参准备样本数据（Stage1转置和配置文件）"""
import argparse
import shutil
from pathlib import Path
import json

def prepare_sample(sample: str, missing_type: str, seed: int = 42):
    root = Path(".")
    sim_dir = root / "data" / "sim" / "real_brca" / "S0" / missing_type.lower().replace(" ", "_") / f"seed_{seed}"
    processed_dir = root / "data" / "processed" / sample
    exported_dir = processed_dir / "stage1_preprocess" / "exported"
    exported_dir.mkdir(parents=True, exist_ok=True)
    
    # 复制SimGen输出到processed（先复制原始文件）
    for f in ["sc_expression.csv", "st_expression.csv", "sc_metadata.csv", "st_coordinates.csv", "sim_info.json", 
              "sim_truth_query_cell_spot.csv", "sim_truth_spot_type_fraction.csv"]:
        src = sim_dir / f
        if src.exists():
            shutil.copy2(src, exported_dir / f)
    
    # 转置表达矩阵（SimGen输出是cell/spot x gene，需要转成gene x cell/spot）
    import pandas as pd
    sc = pd.read_csv(exported_dir / "sc_expression.csv", sep=None, engine="python")
    sc_t = sc.set_index(sc.columns[0]).T
    sc_t.to_csv(exported_dir / "sc_expression_normalized.csv")
    
    st = pd.read_csv(exported_dir / "st_expression.csv", sep=None, engine="python")
    st_t = st.set_index(st.columns[0]).T
    st_t.to_csv(exported_dir / "st_expression_normalized.csv")
    
    # 复制hvg_genes.txt（从已有的样本）
    hvg_src = root / "data" / "processed" / "real_brca_simS0_seed42" / "stage1_preprocess" / "hvg_genes.txt"
    if hvg_src.exists():
        shutil.copy2(hvg_src, processed_dir / "stage1_preprocess" / "hvg_genes.txt")
    
    # 创建配置文件
    cfg_path = root / "configs" / "datasets" / f"{sample}.yaml"
    cfg_content = f"""stage3:
  weak_th: 0.50
  strong_th: 0.80
  abstain_unknown_sc_only: true
"""
    cfg_path.write_text(cfg_content, encoding="utf-8")
    
    print(f"[OK] Prepared sample: {sample}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--missing_type", required=True)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    prepare_sample(args.sample, args.missing_type, args.seed)

