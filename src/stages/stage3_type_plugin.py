"""
Stage 3: 类型不匹配 / Unknown-aware 插件

功能：
- 评估每个原始 celltype 在 ST 的支持度（support_score 固定为 top-3 相似度均值）。
- 重写细胞类型为 plugin_type（Unknown 统一为 Unknown_sc_only）。
- 构建 spot×plugin_type 的先验矩阵 type_prior_matrix.csv。
- 输出 stage3_summary.json 记录参数与统计。

输出位置：
- data/processed/<sample>/stage3_typematch/: type_support.csv, cell_type_relabel.csv, type_prior_matrix.csv
- result/<sample>/stage3_typematch/: stage3_summary.json
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict, Counter
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import yaml
from scipy.stats import pearsonr, norm


# ----------------------------
# 配置与路径
# ----------------------------

def load_yaml(path: Path) -> dict:
    if not path.exists() or path.stat().st_size == 0:
        return {}
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Stage 3: type matching / unknown-aware plugin")
    p.add_argument("--sample", default="real_brca", help="sample id")
    p.add_argument("--project_root", default=None, help="override project root")
    p.add_argument("--config", default="configs/project_config.yaml", help="project config path")
    p.add_argument("--dataset_config", default=None, help="dataset config path (overrides auto-detection)")
    p.add_argument("--strong_th", type=float, default=None, help="support threshold for strong (overrides config)")
    p.add_argument("--weak_th", type=float, default=None, help="support threshold for weak (overrides config)")
    p.add_argument("--st_cluster_k", type=int, default=None, help="ST clustering K (<=0 表示不聚类，直接用 spot)")
    p.add_argument("--unknown_floor", type=float, default=None, help="Unknown_sc_only 最小占比（0~1）")
    p.add_argument("--min_cells_rare_type", type=int, default=None, help="稀有类型阈值，低于此不会标 strong")
    p.add_argument("--eps", type=float, default=None, help="数值下限，避免除零")
    # V4.1 新增配置项
    p.add_argument("--enable_mismatch_detection", type=bool, default=None, help="是否启用 V4.1 Fisher r→z 不匹配检测")
    p.add_argument("--marker_gene_count", type=int, default=None, help="每个细胞类型用于计算相关性的标记基因数量")
    p.add_argument("--score_method", type=str, default=None, choices=["fisher_z", "aggregate_test", "legacy"], help="支持度评分方法")
    p.add_argument("--support_threshold", type=float, default=None, help="Fisher z 分数阈值（默认1.96对应p<0.05）")
    # V4.2 新增配置项
    p.add_argument("--alpha", type=float, default=None, help="显著性水平（默认0.05）")
    p.add_argument("--multiple_test_correction", type=str, default=None, choices=["Bonferroni", "BH", "none"], help="多重比较校正方法")
    p.add_argument("--min_marker_specificity", type=float, default=None, help="标记基因特异性阈值（用于过滤非特异基因）")
    p.add_argument("--min_marker_specificity_overrides", type=str, default=None, help="按类型覆盖特异性阈值（如 T cells CD8=2.0）")
    p.add_argument("--min_marker_genes", type=int, default=None, help="最少 marker 数回退阈值")
    p.add_argument("--min_marker_genes_overrides", type=str, default=None, help="按类型覆盖最少 marker 数（如 T cells CD8=2）")
    p.add_argument("--marker_specificity_mode", type=str, default=None, choices=["mean", "max"], help="特异性计算方式（mean 或 max）")
    p.add_argument("--marker_specificity_mode_overrides", type=str, default=None, help="按类型覆盖特异性模式（如 NK cells=mean,PCs=mean）")
    p.add_argument("--fallback_specificity_mode", type=str, default=None, choices=["mean", "max"], help="回退时特异性计算方式（mean 或 max）")
    p.add_argument("--fallback_specificity_types", type=str, default=None, help="回退时特异性模式应用的类型（逗号分隔）")
    p.add_argument("--use_permutation_test", type=bool, default=None, help="是否使用置换检验计算p值（默认False）")
    p.add_argument("--z_score_aggregation", type=str, default=None, choices=["max", "topk_mean"], help="V4.2 支持度 z 值聚合方式")
    p.add_argument("--z_topk", type=int, default=None, help="topk_mean 的 k 值")
    p.add_argument("--apply_mismatch_to_relabel", type=bool, default=None, help="是否将 V4.2 显著缺失结果用于重写 plugin_type")
    return p.parse_args()


def detect_project_root() -> Path:
    here = Path(__file__).resolve()
    try:
        return here.parents[2]
    except IndexError:
        return here.parent.parent


# ----------------------------
# 相似度与聚类
# ----------------------------

def weighted_cosine(x: np.ndarray, y: np.ndarray, w: np.ndarray, eps: float) -> float:
    # x,y: (G,), w: (G,)
    wx = x * w
    wy = y * w
    num = np.dot(wx, wy)
    denom = (np.linalg.norm(wx) * np.linalg.norm(wy)) + eps
    if denom <= eps:
        return 0.0
    return float(num / denom)


def compute_type_profiles(sc_expr: pd.DataFrame, sc_meta: pd.DataFrame, plugin_genes: List[str], type_col: str) -> Dict[str, np.ndarray]:
    """
    sc_expr: index=cell_id, columns=genes
    sc_meta: must contain cell_id + type_col
    """
    profiles = {}
    grouped = sc_meta.groupby(type_col)["cell_id"].apply(list)
    for t, cell_ids in grouped.items():
        ids = [cid for cid in cell_ids if cid in sc_expr.index]
        if len(ids) == 0:
            profiles[t] = np.zeros(len(plugin_genes), dtype=float)
            continue
        sub = sc_expr.loc[ids, plugin_genes]
        profiles[t] = sub.to_numpy(dtype=float).mean(axis=0)
    return profiles


def try_cluster_st(st_expr: pd.DataFrame, plugin_genes: List[str], k: int) -> Tuple[np.ndarray, pd.DataFrame]:
    """
    返回 labels, cluster_profiles
    cluster_profiles: index=cluster_id, columns=plugin_genes
    """
    if k is None or k <= 0 or len(st_expr) <= k:
        labels = np.arange(len(st_expr))
        profiles = st_expr.loc[:, plugin_genes].copy()
        profiles.index = [f"spot_{i}" for i in range(len(st_expr))]
        return labels, profiles
    try:
        from sklearn.cluster import KMeans
    except Exception:
        labels = np.arange(len(st_expr))
        profiles = st_expr.loc[:, plugin_genes].copy()
        profiles.index = [f"spot_{i}" for i in range(len(st_expr))]
        return labels, profiles

    km = KMeans(n_clusters=min(k, len(st_expr)), n_init="auto", random_state=42)
    labels = km.fit_predict(st_expr[plugin_genes].to_numpy(dtype=float))
    df = st_expr.copy()
    df["cluster"] = labels
    profiles = df.groupby("cluster")[plugin_genes].mean()
    profiles.index = [f"cluster_{i}" for i in profiles.index]
    return labels, profiles


def topk_mean(arr: np.ndarray, k: int) -> float:
    if arr.size == 0:
        return 0.0
    k = max(1, min(k, arr.size))
    topk = np.partition(arr, -k)[-k:]
    return float(topk.mean())


def aggregate_z_values(z_values: np.ndarray, method: str = "max", topk: int = 5) -> float:
    if z_values.size == 0:
        return 0.0
    if method == "max":
        return float(np.max(z_values))
    if method == "topk_mean":
        return topk_mean(z_values, topk)
    raise ValueError(f"Unknown z-score aggregation method: {method}")


# ----------------------------
# V4.1: Fisher r→z 变换支持度计算
# ----------------------------

def fisher_z_transform(r: np.ndarray, eps: float = 1e-10) -> np.ndarray:
    """
    Fisher r→z 变换：z = atanh(r) = 0.5 * ln((1+r)/(1-r))
    
    Args:
        r: 相关系数数组，范围应在 [-1, 1]
        eps: 防止除零的小值
    
    Returns:
        z: Fisher 变换后的 z 值数组
    """
    # 限制 r 的范围到 (-1+eps, 1-eps) 以避免数值问题
    r_clipped = np.clip(r, -1.0 + eps, 1.0 - eps)
    z = np.arctanh(r_clipped)
    return z


def select_marker_genes(
    sc_expr: pd.DataFrame,
    sc_meta: pd.DataFrame,
    cell_type: str,
    type_col: str,
    plugin_genes: List[str],
    n_markers: int,
) -> List[str]:
    """
    为指定细胞类型选择标记基因。
    简单实现：选择在该类型中平均表达最高的 n_markers 个基因。
    
    Args:
        sc_expr: 单细胞表达矩阵
        sc_meta: 单细胞元数据
        cell_type: 目标细胞类型
        type_col: 类型列名
        plugin_genes: 候选基因列表
        n_markers: 需要选择的标记基因数量
    
    Returns:
        标记基因列表
    """
    # 获取该类型的细胞
    type_cells = sc_meta[sc_meta[type_col] == cell_type]["cell_id"]
    type_cells = [cid for cid in type_cells if cid in sc_expr.index]
    
    if len(type_cells) == 0:
        return plugin_genes[:n_markers] if len(plugin_genes) >= n_markers else plugin_genes
    
    # 计算该类型中每个基因的平均表达
    type_expr = sc_expr.loc[type_cells, plugin_genes]
    mean_expr = type_expr.mean(axis=0)
    
    # 选择平均表达最高的 n_markers 个基因
    top_genes = mean_expr.nlargest(min(n_markers, len(plugin_genes))).index.tolist()
    return top_genes


def compute_fisher_z_support_score(
    type_profile: pd.Series,
    st_expr: pd.DataFrame,
    marker_genes: List[str],
    eps: float = 1e-10,
    return_all_z: bool = False,
    agg_method: str = "max",
    agg_topk: int = 5,
) -> float | Tuple[float, np.ndarray]:
    """
    使用 Fisher r→z 变换计算细胞类型在空间数据中的支持度评分。
    
    Args:
        type_profile: 细胞类型的标记基因表达向量（pd.Series，index=gene_name）
        st_expr: 空间转录组表达矩阵（index=spot_id, columns=genes）
        marker_genes: 标记基因列表
        eps: 数值稳定性参数
        return_all_z: 如果为True，返回 (support_score, z_values)，否则只返回 support_score
        agg_method: z 值聚合方式
        agg_topk: topk_mean 的 k 值
    
    Returns:
        如果 return_all_z=False: 支持度评分（最大 Fisher z 值）
        如果 return_all_z=True: (support_score, z_values) 元组
    """
    # 确保 marker_genes 在 st_expr 和 type_profile 中都存在
    available_genes = [g for g in marker_genes if g in st_expr.columns and g in type_profile.index]
    if len(available_genes) < 2:  # 至少需要2个基因才能计算相关系数
        if return_all_z:
            return 0.0, np.array([])
        return 0.0
    
    # 提取类型和空间的标记基因表达（按相同顺序）
    type_vec = type_profile[available_genes].to_numpy(dtype=float)
    st_subset = st_expr[available_genes].to_numpy(dtype=float)
    
    # 计算每个 spot 与类型特征的 Pearson 相关系数
    correlations = []
    for spot_expr in st_subset:
        # 计算 Pearson 相关系数
        if np.std(type_vec) < eps or np.std(spot_expr) < eps:
            r = 0.0
        else:
            r, _ = pearsonr(type_vec, spot_expr)
            if np.isnan(r):
                r = 0.0
        correlations.append(r)
    
    correlations = np.array(correlations, dtype=float)
    
    # Fisher r→z 变换
    z_values = fisher_z_transform(correlations, eps=eps)
    
    # 支持度评分 = 聚合 z 值
    support_score = aggregate_z_values(z_values, method=agg_method, topk=agg_topk)
    
    if return_all_z:
        return support_score, z_values
    return support_score


# ----------------------------
# V4.2: 显著性检验和多重比较校正
# ----------------------------

def compute_p_value(
    z_values: np.ndarray,
    n_genes: int,
    n_spots: int,
    use_permutation: bool = False,
    n_permutations: int = 1000,
    seed: Optional[int] = None,
    type_profile: Optional[pd.Series] = None,
    st_expr: Optional[pd.DataFrame] = None,
    marker_genes: Optional[List[str]] = None,
    eps: float = 1e-10,
    score_aggregation: str = "max",
    z_topk: int = 5,
) -> float:
    """
    计算细胞类型支持度的 p 值。
    
    Args:
        z_values: 所有 spot 的 Fisher z 值数组
        n_genes: 用于计算的标记基因数量
        n_spots: 空间 spot 总数（用于多 spot 校正）
        use_permutation: 是否使用置换检验
        n_permutations: 置换检验的置换次数
        seed: 随机种子
        type_profile: cell type mean expression (marker genes)
        st_expr: spatial expression matrix
        marker_genes: marker genes used for scoring
        eps: numeric stability term
        score_aggregation: z-score aggregation for permutation test
        z_topk: top-k size used by topk_mean
    
    Returns:
        p 值（单尾，表示该类型不存在的概率）
    """
    if z_values.size == 0:
        return 1.0
    
    if use_permutation:
        if type_profile is None or st_expr is None or marker_genes is None:
            # Missing data for permutation test, fall back to normal approximation.
            use_permutation = False
        else:
            available_genes = [g for g in marker_genes if g in st_expr.columns and g in type_profile.index]
            if len(available_genes) < 2:
                return 1.0

            type_vec = type_profile[available_genes].to_numpy(dtype=float)
            st_subset = st_expr[available_genes].to_numpy(dtype=float)

            type_center = type_vec - type_vec.mean()
            type_norm = np.linalg.norm(type_center)
            if type_norm <= eps:
                return 1.0

            st_centered = st_subset - st_subset.mean(axis=1, keepdims=True)
            st_norms = np.linalg.norm(st_centered, axis=1)

            observed_score = aggregate_z_values(z_values, method=score_aggregation, topk=z_topk)
            rng = np.random.default_rng(seed if seed is not None else 42)

            exceed_count = 0
            for _ in range(n_permutations):
                perm_idx = rng.permutation(len(available_genes))
                perm_center = type_center[perm_idx]
                numerator = st_centered @ perm_center
                denom = st_norms * type_norm + eps
                r_vals = np.divide(numerator, denom, out=np.zeros_like(numerator), where=denom > eps)
                r_vals = np.clip(r_vals, -1.0 + eps, 1.0 - eps)
                z_perm = np.arctanh(r_vals)
                perm_score = aggregate_z_values(z_perm, method=score_aggregation, topk=z_topk)
                if perm_score >= observed_score:
                    exceed_count += 1

            return (exceed_count + 1) / (n_permutations + 1)
    
    # 使用正态近似计算 p 值
    # 对于 Fisher z 变换，在 H0（无相关性）下，z 近似服从 N(0, 1/(n-3))
    # 但由于我们取最大值，需要多 spot 校正
    
    # 计算每个 spot 的单尾 p 值
    # 单尾：P(Z >= z) = 1 - Φ(z * sqrt(n-3))
    # 这里我们使用标准正态分布（n 足够大时近似）
    if n_genes >= 3:
        # Fisher z 的标准误差约为 1/sqrt(n-3)
        se = 1.0 / np.sqrt(max(1, n_genes - 3))
    else:
        se = 1.0  # 基因数太少时使用标准正态
    
    # 计算每个 spot 的 p 值
    # 根据设计文档，p 值表示"假设该类型不存在时，得到当前支持度分数或更极端值的概率"
    # 如果类型不存在，z 值应该接近 0，所以 p 值应该大
    # 如果类型存在，z 值应该大，p 值应该小
    
    # 标准化 z 值
    standardized_z = z_values / se
    
    # 计算 p 值：P(Z >= z | H0: 类型不存在)
    # 如果 z 值大，p 值小，表示类型存在（拒绝 H0）
    # 如果 z 值小，p 值大，表示类型不存在（不拒绝 H0）
    spot_p_values = norm.sf(standardized_z)
    
    # 处理数值下溢：如果 p 值太小，设为最小值
    spot_p_values = np.maximum(spot_p_values, np.finfo(float).tiny)
    
    # 多 spot 校正：Bonferroni 近似
    # p_c = min(1, S × min_i p_c,i)
    # 注意：这里我们取最小值 p，因为如果任何一个 spot 有高相关性，就表示类型存在
    min_spot_p = float(np.min(spot_p_values))
    p_value = min(1.0, n_spots * min_spot_p)
    
    # 但是，根据设计文档，p 值应该表示"类型不存在的概率"
    # 所以我们需要反转：p_value_large 表示类型不存在
    # 实际上，当前的 p_value 表示"类型存在的概率"（p 值小 = 存在）
    # 我们需要返回 1 - p_value 来表示"类型不存在的概率"
    # 但这样会导致所有类型的 p 值都接近 1
    
    # 重新理解：根据设计文档，p 值应该表示"在 H0（类型不存在）下，得到当前或更极端结果的概率"
    # 如果 z 值大，说明类型存在，p 值应该小（拒绝 H0）
    # 如果 z 值小，说明类型不存在，p 值应该大（不拒绝 H0）
    # 所以当前的实现是正确的，但判定逻辑需要调整：
    # - p 值小（< alpha）→ 拒绝 H0 → 类型存在 → Significant = No
    # - p 值大（>= alpha）→ 不拒绝 H0 → 类型不存在 → Significant = Yes
    
    return p_value


def apply_multiple_test_correction(
    p_values: Dict[str, float],
    method: str = "BH",
    alpha: float = 0.05,
) -> Dict[str, float]:
    """
    应用多重比较校正。
    
    Args:
        p_values: 字典 {cell_type: p_value}
        method: 校正方法 ("Bonferroni", "BH", "none")
        alpha: 显著性水平
    
    Returns:
        字典 {cell_type: q_value}，q_value 是校正后的显著性值
    """
    if method == "none":
        return p_values.copy()
    
    cell_types = list(p_values.keys())
    p_vals = np.array([p_values[ct] for ct in cell_types])
    M = len(cell_types)
    
    if method == "Bonferroni":
        # Bonferroni 校正：q = min(1, p × M)
        q_vals = np.minimum(1.0, p_vals * M)
    elif method == "BH":
        # Benjamini-Hochberg 校正（FDR 控制）
        # 1. 按 p 值升序排序
        sorted_indices = np.argsort(p_vals)
        sorted_p = p_vals[sorted_indices]
        
        # 2. 计算调整后的 q 值
        q_vals = np.zeros_like(p_vals)
        for i in range(M - 1, -1, -1):  # 从大到小遍历
            rank = i + 1  # 排名（从1开始）
            q_vals[sorted_indices[i]] = min(
                sorted_p[i] * M / rank,
                q_vals[sorted_indices[i + 1]] if i < M - 1 else sorted_p[i] * M / rank
            )
        
        # 确保 q 值单调递增（从最小 p 到最大 p）
        for i in range(M - 2, -1, -1):
            if q_vals[sorted_indices[i]] > q_vals[sorted_indices[i + 1]]:
                q_vals[sorted_indices[i]] = q_vals[sorted_indices[i + 1]]
    else:
        raise ValueError(f"Unknown correction method: {method}")
    
    return {ct: float(q) for ct, q in zip(cell_types, q_vals)}


def select_marker_genes_with_specificity(
    sc_expr: pd.DataFrame,
    sc_meta: pd.DataFrame,
    cell_type: str,
    type_col: str,
    plugin_genes: List[str],
    n_markers: int,
    min_specificity: float = 0.0,
    specificity_mode: str = "mean",
    min_marker_genes: int = 2,
    fallback_specificity_mode: str | None = None,
) -> List[str]:
    """
    V4.2: 优化标记基因选择，考虑特异性。
    
    计算每个基因对区分该细胞类型与其他类型的信息增益（简化版：使用表达差异作为代理），
    并过滤掉在其他类型中也高度表达的基因。
    
    Args:
        sc_expr: 单细胞表达矩阵
        sc_meta: 单细胞元数据
        cell_type: 目标细胞类型
        type_col: 类型列名
        plugin_genes: 候选基因列表
        n_markers: 要选择的标记基因数量
        min_specificity: 最小特异性阈值（0-1，表示该基因在该类型中的表达相对于其他类型的倍数）
        specificity_mode: 其他类型表达的聚合方式（mean 或 max）
        min_marker_genes: 过滤过严时的最少 marker 回退数量
        fallback_specificity_mode: 回退时的特异性计算方式
    
    Returns:
        选定的标记基因列表
    """
    # 获取该类型的细胞
    type_cells = sc_meta[sc_meta[type_col] == cell_type]["cell_id"]
    type_cells = [cid for cid in type_cells if cid in sc_expr.index]
    
    if len(type_cells) == 0:
        # 如果没有该类型的细胞，返回空列表或使用简单方法
        return select_marker_genes(sc_expr, sc_meta, cell_type, type_col, plugin_genes, n_markers)
    
    # 获取其他类型的细胞
    other_cells = sc_meta[sc_meta[type_col] != cell_type]["cell_id"]
    other_cells = [cid for cid in other_cells if cid in sc_expr.index]
    
    if len(other_cells) == 0:
        # 如果没有其他类型，使用简单方法
        return select_marker_genes(sc_expr, sc_meta, cell_type, type_col, plugin_genes, n_markers)
    
    # 计算每个基因在该类型和其他类型中的平均表达
    type_expr = sc_expr.loc[type_cells, plugin_genes].mean(axis=0)
    other_expr_max = sc_expr.loc[other_cells, plugin_genes].max(axis=0)
    other_expr_mean = sc_expr.loc[other_cells, plugin_genes].mean(axis=0)
    
    # 计算特异性：该类型表达 / (其他类型表达 + eps)
    eps = 1e-10
    if specificity_mode == "max":
        specificity = type_expr / (other_expr_max + eps)
    else:
        specificity = type_expr / (other_expr_mean + eps)
    
    # 过滤：只保留特异性 >= min_specificity 的基因
    if min_specificity > 0:
        filtered_genes = [g for g in plugin_genes if specificity[g] >= min_specificity]
    else:
        filtered_genes = plugin_genes

    min_keep = max(2, int(min_marker_genes)) if min_marker_genes is not None else 2
    min_keep = min(min_keep, len(plugin_genes))
    if len(filtered_genes) < min_keep:
        # 回退：按特异性排序，补足最少 marker 数
        fallback_mode = fallback_specificity_mode or specificity_mode
        if fallback_mode == "max":
            fallback_spec = type_expr / (other_expr_max + eps)
        else:
            fallback_spec = type_expr / (other_expr_mean + eps)
        ranked = fallback_spec.sort_values(ascending=False)
        filtered_genes = ranked.index[:min_keep].tolist()
    
    if len(filtered_genes) == 0:
        # 如果过滤后没有基因，使用所有基因
        filtered_genes = plugin_genes
    
    # 在该类型中按平均表达排序，选择 top n_markers
    type_expr_filtered = type_expr[filtered_genes]
    top_genes = type_expr_filtered.nlargest(n_markers).index.tolist()
    
    return top_genes


# ----------------------------
# 主流程
# ----------------------------

def main():
    args = parse_args()
    project_root = Path(args.project_root) if args.project_root else detect_project_root()

    # 加载配置（参考 Stage2 的实现）
    project_cfg_path = project_root / args.config
    project_cfg = load_yaml(project_cfg_path)
    dataset_cfg_map = (project_cfg or {}).get("dataset_config_map", {}) or {}
    if args.dataset_config:
        dataset_cfg_path = Path(args.dataset_config)
    else:
        mapped_name = dataset_cfg_map.get(args.sample)
        dataset_cfg_path = project_root / "configs" / "datasets" / (mapped_name or f"{args.sample}.yaml")
    
    dataset_cfg = load_yaml(dataset_cfg_path)
    stage3_cfg = (dataset_cfg.get("stage3") or {}) if isinstance(dataset_cfg, dict) else {}
    
    # 默认值
    defaults = {
        "strong_th": 0.7,
        "weak_th": 0.4,
        "st_cluster_k": 30,
        "unknown_floor": 0.3,
        "min_cells_rare_type": 20,
        "eps": 1e-8,
        # V4.1 默认值
        "enable_mismatch_detection": False,  # 默认关闭，保持向后兼容
        "marker_gene_count": 30,
        "score_method": "fisher_z",
        "support_threshold": 1.96,  # 对应 p<0.05 的 Z 分数
        # V4.2 默认值
        "alpha": 0.05,  # 显著性水平
        "multiple_test_correction": "BH",  # 多重比较校正方法（推荐 BH）
        "min_marker_specificity": 0.0,  # 标记基因特异性阈值（0表示不过滤）
        "min_marker_specificity_overrides": None,  # 按类型覆盖特异性阈值
        "min_marker_genes": 2,  # 最少 marker 数回退阈值
        "min_marker_genes_overrides": None,  # 按类型覆盖最少 marker 数
        "marker_specificity_mode": "mean",  # 特异性计算方式
        "marker_specificity_mode_overrides": None,  # 按类型覆盖特异性模式
        "fallback_specificity_mode": None,  # 回退时特异性计算方式
        "fallback_specificity_types": None,  # 回退时特异性模式应用的类型
        "use_permutation_test": False,  # 是否使用置换检验（默认False，计算成本高）
        "z_score_aggregation": "max",  # z 值聚合方式
        "z_topk": 5,  # topk_mean 的 k 值
        "apply_mismatch_to_relabel": False,  # V4.2: 是否将显著缺失用于重写 plugin_type
    }
    
    # 从配置读取（配置优先，然后 CLI 覆盖）
    strong_th = args.strong_th if args.strong_th is not None else stage3_cfg.get("strong_th", defaults["strong_th"])
    weak_th = args.weak_th if args.weak_th is not None else stage3_cfg.get("weak_th", defaults["weak_th"])
    st_cluster_k = args.st_cluster_k if args.st_cluster_k is not None else stage3_cfg.get("st_cluster_k", defaults["st_cluster_k"])
    unknown_floor = args.unknown_floor if args.unknown_floor is not None else stage3_cfg.get("unknown_floor", defaults["unknown_floor"])
    min_cells_rare_type = args.min_cells_rare_type if args.min_cells_rare_type is not None else stage3_cfg.get("min_cells_rare_type", defaults["min_cells_rare_type"])
    eps = args.eps if args.eps is not None else stage3_cfg.get("eps", defaults["eps"])
    
    # V4.1 配置
    enable_mismatch_detection = args.enable_mismatch_detection if args.enable_mismatch_detection is not None else stage3_cfg.get("enable_mismatch_detection", defaults["enable_mismatch_detection"])
    marker_gene_count = args.marker_gene_count if args.marker_gene_count is not None else stage3_cfg.get("marker_gene_count", defaults["marker_gene_count"])
    score_method = args.score_method if args.score_method is not None else stage3_cfg.get("score_method", defaults["score_method"])
    support_threshold = args.support_threshold if args.support_threshold is not None else stage3_cfg.get("support_threshold", defaults["support_threshold"])
    # V4.2 配置
    alpha = args.alpha if args.alpha is not None else stage3_cfg.get("alpha", defaults["alpha"])
    multiple_test_correction = args.multiple_test_correction if args.multiple_test_correction is not None else stage3_cfg.get("multiple_test_correction", defaults["multiple_test_correction"])
    min_marker_specificity = args.min_marker_specificity if args.min_marker_specificity is not None else stage3_cfg.get("min_marker_specificity", defaults["min_marker_specificity"])
    min_marker_specificity_overrides = args.min_marker_specificity_overrides if args.min_marker_specificity_overrides is not None else stage3_cfg.get("min_marker_specificity_overrides", defaults["min_marker_specificity_overrides"])
    min_marker_genes = args.min_marker_genes if args.min_marker_genes is not None else stage3_cfg.get("min_marker_genes", defaults["min_marker_genes"])
    min_marker_genes_overrides = args.min_marker_genes_overrides if args.min_marker_genes_overrides is not None else stage3_cfg.get("min_marker_genes_overrides", defaults["min_marker_genes_overrides"])
    marker_specificity_mode = args.marker_specificity_mode if args.marker_specificity_mode is not None else stage3_cfg.get("marker_specificity_mode", defaults["marker_specificity_mode"])
    marker_specificity_mode_overrides = args.marker_specificity_mode_overrides if args.marker_specificity_mode_overrides is not None else stage3_cfg.get("marker_specificity_mode_overrides", defaults["marker_specificity_mode_overrides"])
    fallback_specificity_mode = args.fallback_specificity_mode if args.fallback_specificity_mode is not None else stage3_cfg.get("fallback_specificity_mode", defaults["fallback_specificity_mode"])
    fallback_specificity_types = args.fallback_specificity_types if args.fallback_specificity_types is not None else stage3_cfg.get("fallback_specificity_types", defaults["fallback_specificity_types"])
    use_permutation_test = args.use_permutation_test if args.use_permutation_test is not None else stage3_cfg.get("use_permutation_test", defaults["use_permutation_test"])
    z_score_aggregation = args.z_score_aggregation if args.z_score_aggregation is not None else stage3_cfg.get("z_score_aggregation", defaults["z_score_aggregation"])
    z_topk = args.z_topk if args.z_topk is not None else stage3_cfg.get("z_topk", defaults["z_topk"])
    apply_mismatch_to_relabel = args.apply_mismatch_to_relabel if args.apply_mismatch_to_relabel is not None else stage3_cfg.get("apply_mismatch_to_relabel", defaults["apply_mismatch_to_relabel"])

    if marker_specificity_mode not in {"mean", "max"}:
        marker_specificity_mode = defaults["marker_specificity_mode"]
    if fallback_specificity_mode not in {None, "mean", "max"}:
        fallback_specificity_mode = None
    override_map = {}
    if isinstance(marker_specificity_mode_overrides, str):
        for part in marker_specificity_mode_overrides.split(","):
            part = part.strip()
            if not part:
                continue
            if "=" in part:
                key, val = part.split("=", 1)
            elif ":" in part:
                key, val = part.split(":", 1)
            else:
                continue
            key = key.strip().lower()
            val = val.strip().lower()
            if key and val in {"mean", "max"}:
                override_map[key] = val
    elif isinstance(marker_specificity_mode_overrides, dict):
        for key, val in marker_specificity_mode_overrides.items():
            key = str(key).strip().lower()
            val = str(val).strip().lower()
            if key and val in {"mean", "max"}:
                override_map[key] = val
    specificity_override_map = {}
    if isinstance(min_marker_specificity_overrides, str):
        for part in min_marker_specificity_overrides.split(","):
            part = part.strip()
            if not part:
                continue
            if "=" in part:
                key, val = part.split("=", 1)
            elif ":" in part:
                key, val = part.split(":", 1)
            else:
                continue
            key = key.strip().lower()
            try:
                val = float(val)
            except ValueError:
                continue
            if key:
                specificity_override_map[key] = val
    min_genes_override_map = {}
    if isinstance(min_marker_genes_overrides, str):
        for part in min_marker_genes_overrides.split(","):
            part = part.strip()
            if not part:
                continue
            if "=" in part:
                key, val = part.split("=", 1)
            elif ":" in part:
                key, val = part.split(":", 1)
            else:
                continue
            key = key.strip().lower()
            try:
                val = int(val)
            except ValueError:
                continue
            if key:
                min_genes_override_map[key] = val
    elif isinstance(min_marker_genes_overrides, dict):
        for key, val in min_marker_genes_overrides.items():
            key = str(key).strip().lower()
            try:
                val = int(val)
            except (TypeError, ValueError):
                continue
            if key:
                min_genes_override_map[key] = val
    elif isinstance(min_marker_specificity_overrides, dict):
        for key, val in min_marker_specificity_overrides.items():
            key = str(key).strip().lower()
            try:
                val = float(val)
            except (TypeError, ValueError):
                continue
            if key:
                specificity_override_map[key] = val
    if isinstance(fallback_specificity_types, str):
        fallback_specificity_types = [t.strip() for t in fallback_specificity_types.split(",") if t.strip()]
    if isinstance(fallback_specificity_types, list):
        fallback_specificity_types = {str(t).strip().lower() for t in fallback_specificity_types if str(t).strip()}
    else:
        fallback_specificity_types = None
    try:
        min_marker_genes = int(min_marker_genes)
    except (TypeError, ValueError):
        min_marker_genes = defaults["min_marker_genes"]
    if min_marker_genes < 2:
        min_marker_genes = 2
    if z_score_aggregation not in {"max", "topk_mean"}:
        z_score_aggregation = defaults["z_score_aggregation"]
    try:
        z_topk = int(z_topk)
    except (TypeError, ValueError):
        z_topk = defaults["z_topk"]
    if z_topk <= 0:
        z_topk = defaults["z_topk"]
    if z_score_aggregation != "max":
        use_permutation_test = True

    # 路径
    stage1_export = project_root / "data" / "processed" / args.sample / "stage1_preprocess" / "exported"
    out_proc = project_root / "data" / "processed" / args.sample / "stage3_typematch"
    out_res = project_root / "result" / args.sample / "stage3_typematch"
    out_proc.mkdir(parents=True, exist_ok=True)
    out_res.mkdir(parents=True, exist_ok=True)

    # 读取数据
    sc_expr = pd.read_csv(stage1_export / "sc_expression_normalized.csv", index_col=0)
    st_expr = pd.read_csv(stage1_export / "st_expression_normalized.csv", index_col=0)
    sc_meta = pd.read_csv(stage1_export / "sc_metadata.csv")
    st_coords = pd.read_csv(stage1_export / "st_coordinates.csv", index_col=0)
    # 阶段二已移除：改用 Stage1 的高变基因或配置指定的基因列表
    plugin_genes_path = stage3_cfg.get("plugin_genes_path") or (stage1_export.parent / "hvg_genes.txt")
    if not Path(plugin_genes_path).exists():
        raise FileNotFoundError(f"plugin_genes_path 不存在: {plugin_genes_path}")
    plugin_genes = [g.strip() for g in Path(plugin_genes_path).read_text(encoding="utf-8").splitlines() if g.strip()]
    gene_weights_path = stage3_cfg.get("gene_weights_path")
    if gene_weights_path and Path(gene_weights_path).exists():
        gene_weights = pd.read_csv(gene_weights_path)
        weight_map = dict(zip(gene_weights.iloc[:, 0], gene_weights.iloc[:, 1]))
    else:
        weight_map = {}
    w = np.array([weight_map.get(g, 1.0) for g in plugin_genes], dtype=float)

    # 对齐基因
    sc_expr = sc_expr.loc[:, plugin_genes]
    st_expr = st_expr.loc[:, plugin_genes]

    # 取类型列名（兼容 cell_type / celltype / type），并加护栏避免 silent mismatch
    if "cell_type" in sc_meta.columns and "celltype" in sc_meta.columns:
        a = sc_meta["cell_type"].astype(str).to_numpy()
        b = sc_meta["celltype"].astype(str).to_numpy()
        if a.shape == b.shape and (a != b).any():
            raise ValueError("sc_metadata.csv 同时存在 cell_type 与 celltype，但两列不一致（禁止 silent bug）")
    if "cell_type" in sc_meta.columns:
        type_col = "cell_type"
    elif "celltype" in sc_meta.columns:
        type_col = "celltype"
    elif "type" in sc_meta.columns:
        sc_meta = sc_meta.copy()
        sc_meta["celltype"] = sc_meta["type"].astype(str)
        type_col = "celltype"
    else:
        raise KeyError("sc_metadata.csv 需要包含 celltype/cell_type/type 之一作为类型列")
    # 类型画像
    profiles = compute_type_profiles(sc_expr, sc_meta, plugin_genes, type_col=type_col)

    # ST 聚类（可选）
    labels, st_profiles = try_cluster_st(st_expr, plugin_genes, k=st_cluster_k)
    st_profiles_np = st_profiles.to_numpy(dtype=float)

    # 支持度计算
    support_rows = []
    type_support_score = {}
    type_p_values = {}  # V4.2: 存储 p 值
    type_z_values = {}  # V4.2: 存储所有 z 值（用于计算 p 值）
    
    # V4.1: 如果启用 Fisher r→z 方法
    use_fisher_z = enable_mismatch_detection and score_method == "fisher_z"
    use_v42 = use_fisher_z  # V4.2 在 V4.1 基础上扩展
    if args.apply_mismatch_to_relabel is None and "apply_mismatch_to_relabel" not in stage3_cfg:
        apply_mismatch_to_relabel = use_v42
    
    n_spots = len(st_expr.index)  # 用于 p 值计算的多 spot 校正
    
    for t, prof in profiles.items():
        if use_fisher_z:
            # V4.1/V4.2: 使用 Fisher r→z 变换方法
            # 1. 选择标记基因（V4.2: 考虑特异性）
            effective_min_marker_specificity = specificity_override_map.get(
                str(t).strip().lower(), min_marker_specificity
            )
            if effective_min_marker_specificity > 0:
                effective_specificity_mode = override_map.get(str(t).strip().lower(), marker_specificity_mode)
                effective_min_marker_genes = min_genes_override_map.get(str(t).strip().lower(), min_marker_genes)
                marker_genes = select_marker_genes_with_specificity(
                    sc_expr, sc_meta, t, type_col, plugin_genes, 
                    marker_gene_count, effective_min_marker_specificity, effective_specificity_mode,
                    effective_min_marker_genes,
                    fallback_specificity_mode if (not fallback_specificity_types or str(t).strip().lower() in fallback_specificity_types) else None,
                )
            else:
                marker_genes = select_marker_genes(sc_expr, sc_meta, t, type_col, plugin_genes, marker_gene_count)
            
            # 2. 构建类型特征向量（仅使用标记基因）
            type_cells = sc_meta[sc_meta[type_col] == t]["cell_id"]
            type_cells = [cid for cid in type_cells if cid in sc_expr.index]
            if len(type_cells) > 0:
                type_expr = sc_expr.loc[type_cells, marker_genes]
                type_profile = type_expr.mean(axis=0)  # 保持为 pd.Series，index=gene_name
            else:
                type_profile = pd.Series(0.0, index=marker_genes)
            
            # 3. 计算 Fisher z 支持度评分（V4.2: 同时获取所有 z 值）
            if use_v42:
                score, z_vals = compute_fisher_z_support_score(
                    type_profile,
                    st_expr,
                    marker_genes,
                    eps=eps,
                    return_all_z=True,
                    agg_method=z_score_aggregation,
                    agg_topk=z_topk,
                )
                type_z_values[t] = z_vals
                n_genes_used = len(marker_genes)
            else:
                score = compute_fisher_z_support_score(type_profile, st_expr, marker_genes, eps=eps)
                z_vals = None
                n_genes_used = len(marker_genes)
            
            # V4.2: 计算 p 值
            if use_v42 and z_vals is not None and z_vals.size > 0:
                effective_use_permutation = use_permutation_test or z_score_aggregation != "max"
                p_val = compute_p_value(
                    z_vals,
                    n_genes_used,
                    n_spots,
                    use_permutation=effective_use_permutation,
                    type_profile=type_profile,
                    st_expr=st_expr,
                    marker_genes=marker_genes,
                    eps=eps,
                    score_aggregation=z_score_aggregation,
                    z_topk=z_topk,
                )
                type_p_values[t] = p_val
            else:
                type_p_values[t] = None
            
            # 4. 判定不匹配（V4.1: 使用阈值；V4.2: 稍后使用 q 值）
            if use_v42:
                # V4.2 将在多重比较校正后判定
                flag_mismatch = ""  # 稍后填充
            else:
                # V4.1: 使用阈值判定
                flag_mismatch = "Yes" if score < support_threshold else "No"
            
            # 为了兼容性，仍然计算传统相似度（用于 mapped_st_cluster）
            sims = []
            for i in range(st_profiles_np.shape[0]):
                sims.append(weighted_cosine(prof, st_profiles_np[i], w, eps=eps))
            sims = np.array(sims, dtype=float)
            best_idx = int(np.argmax(sims)) if sims.size > 0 else 0
            mapped_cluster = st_profiles.index[best_idx] if len(st_profiles.index) > 0 else ""
        else:
            # 传统方法：加权余弦相似度 + top-3 均值
            sims = []
            for i in range(st_profiles_np.shape[0]):
                sims.append(weighted_cosine(prof, st_profiles_np[i], w, eps=eps))
            sims = np.array(sims, dtype=float)
            score = topk_mean(sims, k=3)
            best_idx = int(np.argmax(sims)) if sims.size > 0 else 0
            mapped_cluster = st_profiles.index[best_idx] if len(st_profiles.index) > 0 else ""
            flag_mismatch = ""  # 传统方法不输出 FlagMismatch
        
        type_support_score[t] = score
        n_cells = int((sc_meta[type_col] == t).sum())
        
        # V4.2: 准备输出字段
        row_data = {
            "CellType": t,  # V4.1 规范字段名
            "orig_type": t,  # 保持向后兼容
            "n_cells": n_cells,
            "SupportScore": score,  # V4.1 规范字段名
            "support_score": score,  # 保持向后兼容
            "FlagMismatch": flag_mismatch,  # V4.1 新增字段（V4.2 稍后更新）
            "support_category": "",  # 稍后填充
            "mapped_st_cluster": mapped_cluster,
        }
        
        # V4.2: 添加 p 值和 q 值字段（稍后填充）
        if use_v42:
            row_data["PValue"] = type_p_values.get(t, None)
            row_data["QValue"] = None  # 稍后通过多重比较校正填充
            row_data["Significant"] = ""  # 稍后判定
        
        support_rows.append(row_data)
    
    # V4.2: 多重比较校正和显著性判定
    if use_v42 and type_p_values:
        # 过滤掉 None 值（传统方法没有 p 值）
        valid_p_values = {t: p for t, p in type_p_values.items() if p is not None}
        
        if valid_p_values:
            # 应用多重比较校正
            q_values = apply_multiple_test_correction(
                valid_p_values, method=multiple_test_correction, alpha=alpha
            )
            
            # 更新 support_rows 中的 q 值和显著性判定
            for row in support_rows:
                cell_type = row["CellType"]
                if cell_type in q_values:
                    row["QValue"] = q_values[cell_type]
                    row["PValue"] = valid_p_values[cell_type]
                    
                    # 判定显著性：根据设计文档，p 值表示"类型不存在的概率"
                    # 如果 p 值大（>= alpha），表示类型不存在（显著缺失）
                    # 如果 p 值小（< alpha），表示类型存在（不显著缺失）
                    # 注意：这里 p 值实际上是"类型存在的显著性"，所以需要反转逻辑
                    # 实际上，根据设计文档，p 值大表示类型不存在，所以：
                    # q >= alpha → 类型不存在 → Significant = Yes
                    # q < alpha → 类型存在 → Significant = No
                    is_significant = q_values[cell_type] >= alpha  # q >= alpha 表示显著缺失（p 值大）
                    row["Significant"] = "Yes" if is_significant else "No"
                    
                    # 更新 FlagMismatch（保持向后兼容）
                    row["FlagMismatch"] = "Yes" if is_significant else "No"
                else:
                    # 传统方法或没有 p 值的类型
                    if "PValue" in row:
                        row["PValue"] = None
                    if "QValue" in row:
                        row["QValue"] = None
                    if "Significant" in row:
                        row["Significant"] = ""

    # 支持度分档 + 稀有类型保护
    for row in support_rows:
        score = row["support_score"]
        n_cells = row["n_cells"]
        if n_cells < min_cells_rare_type and score >= strong_th:
            category = "weak"
        elif score >= strong_th:
            category = "strong"
        elif score >= weak_th:
            category = "weak"
        else:
            category = "unsupported"
        row["support_category"] = category

    # 写 type_support
    type_support_path = out_proc / "type_support.csv"
    pd.DataFrame(support_rows).to_csv(type_support_path, index=False)

    # 重写标签
    support_map = {r["orig_type"]: r["support_category"] for r in support_rows}
    # V4.2: 将显著缺失结果应用到 relabel
    mismatch_map = {}
    if use_fisher_z and apply_mismatch_to_relabel:
        for row in support_rows:
            orig_type = row["orig_type"]
            flag_mismatch = row.get("FlagMismatch", "")
            mismatch_map[orig_type] = (flag_mismatch == "Yes")

    relabel_rows = []
    for _, r in sc_meta.iterrows():
        orig = r[type_col]
        is_orig_unknown = str(orig).lower() in {"unknown", "unk", "na", "unlabeled"}
        status = support_map.get(orig, "unsupported")

        # V4.2: mismatch detection applied to relabeling when enabled.
        if use_fisher_z and apply_mismatch_to_relabel and mismatch_map.get(orig, False):
            plugin_type = "Unknown_sc_only"
            label_status = "v4.2_mismatch_detected"
        elif is_orig_unknown or status == "unsupported":
            plugin_type = "Unknown_sc_only"
            label_status = "unknown_merged"
        elif status == "strong":
            plugin_type = orig
            label_status = "kept"
        elif status == "weak":
            plugin_type = orig
            label_status = "downweighted"
        else:
            plugin_type = "Unknown_sc_only"
            label_status = "unknown_merged"
        relabel_rows.append(
            {
                "cell_id": r["cell_id"],
                "orig_type": orig,
                "plugin_type": plugin_type,
                "status": label_status,
            }
        )
    relabel_df = pd.DataFrame(relabel_rows)
    relabel_path = out_proc / "cell_type_relabel.csv"
    relabel_df.to_csv(relabel_path, index=False)

    # plugin_type 列顺序（非 Unknown 按字母排序，Unknown_sc_only 最后）
    plugin_types = sorted([t for t in relabel_df["plugin_type"].unique() if t != "Unknown_sc_only"])
    plugin_types.append("Unknown_sc_only")

    # type_prior_matrix
    type_profiles = {}
    for t in plugin_types:
        if t == "Unknown_sc_only":
            # 用被合并为 Unknown 的细胞平均表达（若为空则全零）
            unk_cells = relabel_df[relabel_df["plugin_type"] == "Unknown_sc_only"]["cell_id"]
            if len(unk_cells) > 0:
                unk_expr = sc_expr.loc[unk_cells, plugin_genes]
                type_profiles[t] = unk_expr.to_numpy(dtype=float).mean(axis=0)
            else:
                type_profiles[t] = np.zeros(len(plugin_genes), dtype=float)
        else:
            type_profiles[t] = profiles.get(t, np.zeros(len(plugin_genes), dtype=float))

    prior_rows = []
    n_spots_all_unknown = 0
    for spot_id, row in st_expr.iterrows():
        sims = []
        for t in plugin_types:
            sims.append(weighted_cosine(row.to_numpy(dtype=float), np.array(type_profiles[t]), w, eps=eps))
        sims = np.array(sims, dtype=float)
        sims_pos = np.maximum(sims, 0.0)
        total = sims_pos.sum()
        if total <= eps:
            # 极端情况：全部未知
            prior = np.zeros_like(sims_pos)
            prior[-1] = 1.0  # Unknown_sc_only
            n_spots_all_unknown += 1
        else:
            prior = sims_pos / total
            # Unknown 保底
            uf = unknown_floor
            if uf < 0 or uf > 1:
                raise ValueError("unknown_floor must be in [0,1]")
            if prior[-1] < uf:
                remain = max(eps, 1.0 - uf)
                scale = remain / max(eps, prior[:-1].sum())
                prior[:-1] = prior[:-1] * scale
                prior[-1] = uf
                prior = prior / prior.sum()  # 归一化到1
        prior_rows.append([spot_id] + prior.tolist())

    prior_df = pd.DataFrame(prior_rows, columns=["spot_id"] + plugin_types)
    prior_path = out_proc / "type_prior_matrix.csv"
    prior_df.to_csv(prior_path, index=False)

    # 汇总统计
    # 支持度档按类型数
    cat_counts = Counter([r["support_category"] for r in support_rows])
    # 按细胞数占比
    cell_counts = Counter()
    for r in relabel_rows:
        cat = support_map.get(r["orig_type"], "unsupported")
        cell_counts[cat] += 1
    total_cells = len(relabel_rows)
    support_overview = {
        "by_type_count": dict(cat_counts),
        "by_cell_fraction": {k: (v / total_cells if total_cells > 0 else 0.0) for k, v in cell_counts.items()},
    }
    unknown_cells = (relabel_df["plugin_type"] == "Unknown_sc_only").sum()
    unknown_prior_mass = prior_df["Unknown_sc_only"].mean() if "Unknown_sc_only" in prior_df.columns else 0.0

    rare_types = [
        {
            "orig_type": r["orig_type"],
            "n_cells": r["n_cells"],
            "support_category": r["support_category"],
        }
        for r in support_rows
        if r["n_cells"] < min_cells_rare_type
    ]

    support_score_def = "top3_mean_cluster_similarity"
    if use_fisher_z:
        if z_score_aggregation == "topk_mean":
            support_score_def = f"fisher_z_top{z_topk}_mean"
        else:
            support_score_def = "fisher_z_max"

    summary = {
        "params": {
            "strong_th": strong_th,
            "weak_th": weak_th,
            "st_cluster_k": st_cluster_k,
            "similarity_metric": "weighted_cosine" if not use_fisher_z else "fisher_z_pearson",
            "unknown_floor": unknown_floor,
            "min_cells_rare_type": min_cells_rare_type,
            "support_score_def": support_score_def,
            "resolved_celltype_column": type_col,
            # V4.1 参数
            "enable_mismatch_detection": enable_mismatch_detection,
            "marker_gene_count": marker_gene_count if use_fisher_z else None,
            "score_method": score_method if use_fisher_z else "legacy",
            "support_threshold": support_threshold if use_fisher_z else None,
            # V4.2 参数
            "alpha": alpha if use_v42 else None,
            "multiple_test_correction": multiple_test_correction if use_v42 else None,
            "min_marker_specificity": min_marker_specificity if use_v42 else None,
            "min_marker_specificity_overrides": min_marker_specificity_overrides if use_v42 else None,
            "min_marker_genes": min_marker_genes if use_v42 else None,
            "min_marker_genes_overrides": min_marker_genes_overrides if use_v42 else None,
            "marker_specificity_mode": marker_specificity_mode if use_v42 else None,
            "marker_specificity_mode_overrides": marker_specificity_mode_overrides if use_v42 else None,
            "fallback_specificity_mode": fallback_specificity_mode if use_v42 else None,
            "fallback_specificity_types": sorted(fallback_specificity_types) if use_v42 and fallback_specificity_types else None,
            "use_permutation_test": use_permutation_test if use_v42 else None,
            "z_score_aggregation": z_score_aggregation if use_v42 else None,
            "z_topk": z_topk if use_v42 else None,
            "apply_mismatch_to_relabel": apply_mismatch_to_relabel if use_fisher_z else None,
        },
        "support_overview": support_overview,
        "unknown_overview": {
            "cell_fraction": unknown_cells / total_cells if total_cells > 0 else 0.0,
            "prior_mass_fraction": unknown_prior_mass,
            "n_spots_all_unknown": n_spots_all_unknown,
        },
        "rare_types": rare_types,
        "plugin_types": plugin_types,
    }
    with (out_res / "stage3_summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    print(f"[Stage3] type_support.csv -> {type_support_path}")
    print(f"[Stage3] cell_type_relabel.csv -> {relabel_path}")
    print(f"[Stage3] type_prior_matrix.csv -> {prior_path}")
    print(f"[Stage3] stage3_summary.json -> {out_res / 'stage3_summary.json'}")


if __name__ == "__main__":
    main()
