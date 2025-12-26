#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""测试 V4.2 p 值计算"""

import pandas as pd
import numpy as np
from scipy.stats import norm

# 读取结果
df = pd.read_csv('data/processed/real_brca_simS0_seed42/stage3_typematch/type_support.csv')

print("="*80)
print("V4.2 p 值计算结果分析")
print("="*80)

print("\n所有类型的统计信息:")
print(df[['CellType', 'SupportScore', 'PValue', 'QValue', 'Significant', 'FlagMismatch']].to_string())

print("\n" + "="*80)
print("重点检查:")
print("="*80)

# 检查 T cells CD8 (missing type)
cd8 = df[df['CellType'] == 'T cells CD8'].iloc[0]
print(f"\nT cells CD8 (missing type):")
print(f"  SupportScore: {cd8['SupportScore']:.4f}")
print(f"  PValue: {cd8['PValue']}")
print(f"  QValue: {cd8['QValue']}")
print(f"  Significant: {cd8['Significant']}")
print(f"  FlagMismatch: {cd8['FlagMismatch']}")

# 检查其他类型
print(f"\n其他类型（应该存在）:")
for _, row in df[df['CellType'] != 'T cells CD8'].iterrows():
    if row['PValue'] > 0:
        print(f"  {row['CellType']}: SupportScore={row['SupportScore']:.4f}, PValue={row['PValue']:.2e}, QValue={row['QValue']:.2e}")

print(f"\n问题分析:")
print(f"  - 大部分类型的 PValue 为 0.0，这可能表示 z 值太大，导致 norm.cdf 接近 1")
print(f"  - 对于 missing type (T cells CD8)，PValue 也应该是较大的值（接近 1），但现在也是 0.0")
print(f"  - 这可能表示 p 值计算逻辑有问题")

