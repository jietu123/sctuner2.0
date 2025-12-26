#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""检查被误标记为 FlagMismatch=Yes 的非 missing 类型"""

import pandas as pd
import sys
from pathlib import Path

root = Path(__file__).parent.parent
sys.path.insert(0, str(root))

def main():
    sample = "real_brca_simS0_seed42"
    missing_type = "T cells CD8"
    
    # 读取 type_support.csv
    support_path = root / "data" / "processed" / sample / "stage3_typematch" / "type_support.csv"
    if not support_path.exists():
        print(f"[ERROR] File not found: {support_path}")
        return
    
    df = pd.read_csv(support_path)
    
    # 找出所有被标记为 FlagMismatch=Yes 的类型
    flagged = df[df['FlagMismatch'] == 'Yes'].copy()
    
    # 排除 missing type
    false_positives = flagged[flagged['CellType'] != missing_type].copy()
    
    print("="*80)
    print("误标记分析 (False Positives)")
    print("="*80)
    print(f"\nMissing type: {missing_type}")
    missing_row = df[df['CellType'] == missing_type].iloc[0]
    print(f"Missing type SupportScore: {missing_row['SupportScore']:.4f}")
    print(f"Current threshold: 1.85")
    print(f"\n被误标记的非 missing 类型数量: {len(false_positives)}")
    
    if len(false_positives) > 0:
        print("\n误标记的类型详情:")
        print(false_positives[['CellType', 'SupportScore', 'FlagMismatch']].to_string(index=False))
        
        print("\n误标记类型的 SupportScore 分布:")
        print(f"  最小值: {false_positives['SupportScore'].min():.4f}")
        print(f"  最大值: {false_positives['SupportScore'].max():.4f}")
        print(f"  平均值: {false_positives['SupportScore'].mean():.4f}")
        print(f"  中位数: {false_positives['SupportScore'].median():.4f}")
        
        print(f"\nMissing type SupportScore: {missing_row['SupportScore']:.4f}")
        print(f"\n问题分析:")
        print(f"  - Missing type 的 SupportScore ({missing_row['SupportScore']:.4f}) 与误标记类型的 SupportScore 很接近")
        print(f"  - 如果提高阈值以减少误标记，可能会漏掉 missing type")
        print(f"  - 如果降低阈值以确保检测 missing type，会增加误标记")
        
        # 检查是否有阈值可以同时检测到 missing 且减少误标记
        print(f"\n阈值分析:")
        all_scores = sorted(df['SupportScore'].unique(), reverse=True)
        missing_score = missing_row['SupportScore']
        
        # 找出所有非 missing 类型的 SupportScore
        non_missing_scores = df[df['CellType'] != missing_type]['SupportScore'].values
        non_missing_below_threshold = [s for s in non_missing_scores if s < 1.85]
        
        print(f"  当前阈值 1.85 下，非 missing 类型中 SupportScore < 1.85 的数量: {len(non_missing_below_threshold)}")
        if len(non_missing_below_threshold) > 0:
            print(f"  这些类型的 SupportScore 范围: [{min(non_missing_below_threshold):.4f}, {max(non_missing_below_threshold):.4f}]")
        
        # 尝试找到一个更好的阈值
        print(f"\n可能的优化方向:")
        print(f"  1. 如果这些误标记类型的 SupportScore 都 < missing type 的 SupportScore，可以设置阈值在它们之间")
        fp_scores = false_positives['SupportScore'].values
        if len(fp_scores) > 0 and all(s < missing_score for s in fp_scores):
            optimal_threshold = (max(fp_scores) + missing_score) / 2
            print(f"    建议阈值: {optimal_threshold:.4f} (在误标记最高分 {max(fp_scores):.4f} 和 missing {missing_score:.4f} 之间)")
        elif len(fp_scores) > 0 and any(s > missing_score for s in fp_scores):
            print(f"    问题: 部分误标记类型的 SupportScore 高于 missing type")
            print(f"    这表明 Fisher z-transform 方法可能无法完全区分这些类型")
            print(f"    建议: 考虑使用更复杂的统计方法（如 V4.2 的显著性检验）")
    else:
        print("\n[OK] 没有误标记！")

if __name__ == "__main__":
    main()
