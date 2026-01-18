"""
Stage 3: 类型不匹配 / Unknown-aware 插件

功能：
- 评估每个原始 celltype 在 ST 的支持度（support_score 固定为 top-3 相似度均值）。
- 重写细胞类型为 plugin_type（Unknown 统一为 Unknown_sc_only）。
- 构建 spot×plugin_type 的先验矩阵 type_prior_matrix.csv。
- 输出 stage3_summary.json 记录参数与统计。

输出位置：
- data/processed/<sample>/stage3_typematch/: type_support.csv, cell_type_relabel.csv, type_prior_matrix.csv, stage3_adjusted_annotations.csv
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
from scipy.stats import pearsonr, norm, spearmanr


# 模块流程概览：
# 1) 读取项目与数据集配置，统一 Stage3 参数来源（defaults + dataset + CLI）。
# 2) 计算类型支持度：V4.1 Fisher r->z 与 V4.2 显著性检验/多重校正。
# 3) 生成 Unknown / Dropped / Relabel 与 type_prior_matrix。
# 4) V5.x 进行去噪、生态位拯救、熵质控与拯救控制。
# 5) 输出 stage3_summary.json 与 data/processed 下的中间文件。
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
    p.add_argument("--min_effect_size", type=float, default=None, help="灰区效应量阈值（支持度过低也进入灰区）")
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
    p.add_argument(
        "--marker_selection_strategy",
        type=str,
        default=None,
        choices=["expression", "specificity"],
        help="marker selection strategy: expression (default) or specificity",
    )
    p.add_argument("--use_permutation_test", type=bool, default=None, help="是否使用置换检验计算p值（默认False）")
    p.add_argument("--z_score_aggregation", type=str, default=None, choices=["max", "topk_mean"], help="V4.2 支持度 z 值聚合方式")
    p.add_argument("--z_topk", type=int, default=None, help="topk_mean 的 k 值")
    p.add_argument("--apply_mismatch_to_relabel", type=bool, default=None, help="是否将 V4.2 显著缺失结果用于重写 plugin_type")
    # V4.3 新增配置项
    p.add_argument("--mismatch_action", type=str, default=None, choices=["relabel", "mark_unknown", "ignore"], help="缺失类型处理策略")
    p.add_argument("--relabel_similarity_threshold", type=float, default=None, help="重标注相似度阈值")
    p.add_argument("--merge_target_map", type=str, default=None, help="手工指定合并映射（JSON或A=B,B=C）")
    p.add_argument("--unknown_label_prefix", type=str, default=None, help="unknown类型标签前缀")
    p.add_argument("--drop_unknown", type=bool, default=None, help="unknown类型是否从映射中剔除")
    p.add_argument("--support_margin_min", type=float, default=None, help="strong ?????????top1-top2 z?")
    p.add_argument("--support_margin_enable", type=bool, default=None, help="????????????")
    p.add_argument("--conflict_demotion_enable", type=bool, default=None, help="????? strong ???")
    p.add_argument(
        "--protect_strong_from_missing",
        type=bool,
        default=None,
        help="?????????? strong????? Keep???? missing_types",
    )
    return p.parse_args()


def detect_project_root() -> Path:
    here = Path(__file__).resolve()
    try:
        return here.parents[2]
    except IndexError:
        return here.parent.parent


# ----------------------------
# 相似度与聚类
# - 计算 ST/SC 在 marker 基因空间中的相似度。
# - 可选 KMeans 将 ST spot 聚类，降低噪声与计算量；失败时回退到逐 spot。
# ----------------------------

def weighted_cosine(x: np.ndarray, y: np.ndarray, w: np.ndarray, eps: float) -> float:
    # x,y: (G,), w: (G,)
    # Weighted cosine with stability guard for near-zero norms.
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


def normalize_type_key(name: str | None) -> str:
    return " ".join(str(name or "").strip().lower().split())


def parse_merge_target_map(raw: object) -> Dict[str, str]:
    if raw is None:
        return {}
    if isinstance(raw, dict):
        return {str(k).strip(): str(v).strip() for k, v in raw.items() if str(k).strip() and str(v).strip()}
    if isinstance(raw, str):
        raw = raw.strip()
        if not raw:
            return {}
        try:
            data = json.loads(raw)
            if isinstance(data, dict):
                return {str(k).strip(): str(v).strip() for k, v in data.items() if str(k).strip() and str(v).strip()}
        except Exception:
            pass
        mapping = {}
        for part in raw.split(","):
            part = part.strip()
            if not part:
                continue
            if "=" in part:
                key, val = part.split("=", 1)
            elif ":" in part:
                key, val = part.split(":", 1)
            else:
                continue
            key = key.strip()
            val = val.strip()
            if key and val:
                mapping[key] = val
        return mapping
    return {}


def safe_pearson(a: np.ndarray, b: np.ndarray, eps: float = 1e-12) -> float:
    if a.size == 0 or b.size == 0:
        return 0.0
    if np.std(a) < eps or np.std(b) < eps:
        return 0.0
    r = pearsonr(a, b)[0]
    if not np.isfinite(r):
        return 0.0
    return float(r)


def safe_spearman(a: np.ndarray, b: np.ndarray, eps: float = 1e-12) -> float:
    if a.size == 0 or b.size == 0:
        return 0.0
    if np.std(a) < eps or np.std(b) < eps:
        return 0.0
    r = spearmanr(a, b)[0]
    if not np.isfinite(r):
        return 0.0
    return float(r)


def try_cluster_st(st_expr: pd.DataFrame, plugin_genes: List[str], k: int) -> Tuple[np.ndarray, pd.DataFrame]:
    """
    返回 labels, cluster_profiles
    cluster_profiles: index=cluster_id, columns=plugin_genes
    """
    # Fall back to no clustering when k is invalid or too large.
    if k is None or k <= 0 or len(st_expr) <= k:
        labels = np.arange(len(st_expr))
        profiles = st_expr.loc[:, plugin_genes].copy()
        profiles.index = [f"spot_{i}" for i in range(len(st_expr))]
        return labels, profiles
    try:
        from sklearn.cluster import KMeans
    except Exception:
        # If sklearn is missing, treat each spot as its own cluster.
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


def compute_normalized_entropy(z_values: np.ndarray, temperature: float = 1.0, eps: float = 1e-12) -> float:
    if z_values.size <= 1:
        return 0.0
    # Softmax with temperature, then normalize entropy by log(K).
    temp = max(eps, float(temperature))
    logits = z_values / temp
    logits = logits - np.max(logits)
    exp_vals = np.exp(logits)
    probs = exp_vals / (np.sum(exp_vals) + eps)
    entropy = -np.sum(probs * np.log(probs + eps))
    entropy_max = np.log(float(z_values.size))
    if entropy_max <= eps:
        return 0.0
    return float(entropy / entropy_max)


# ----------------------------
# V4.1: Fisher r→z 变换支持度计算
# - 将相关系数转为 z 值，得到更稳定的支持度统计量。
# - 后续用于显著性判定与灰区/缺失类型识别。
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
    # Simple marker selection: top mean expression within the target type.
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
    # Compute per-spot correlations, then aggregate Fisher z scores.
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
# V5.1: Signal denoising (marker pruning)
# - 通过 CV 评估 marker 噪声，按比例剪枝，提升弱信号的可辨识性。
# ----------------------------

def compute_marker_noise_scores_cv(
    st_expr: pd.DataFrame,
    marker_genes: List[str],
    eps: float = 1e-10,
) -> pd.Series:
    # Use inverse CV as a stability score (higher = less noisy).
    available = [g for g in marker_genes if g in st_expr.columns]
    if not available:
        return pd.Series(dtype=float)
    sub = st_expr[available]
    mean = sub.mean(axis=0)
    std = sub.std(axis=0)
    cv = std / (mean + eps)
    noise = 1.0 / (cv + eps)
    return noise.sort_values(ascending=False)


def prune_noisy_markers(
    marker_genes: List[str],
    st_expr: pd.DataFrame,
    max_pruning_ratio: float,
    min_markers_left: int,
    eps: float = 1e-10,
) -> Tuple[List[str], List[str]]:
    # Drop the noisiest markers up to the configured pruning ratio.
    n_total = len(marker_genes)
    if n_total <= min_markers_left:
        return marker_genes, []
    max_pruning_ratio = max(0.0, min(1.0, max_pruning_ratio))
    drop_n = int(n_total * max_pruning_ratio)
    drop_n = min(drop_n, n_total - min_markers_left)
    if drop_n <= 0:
        return marker_genes, []
    noise_scores = compute_marker_noise_scores_cv(st_expr, marker_genes, eps=eps)
    if noise_scores.empty:
        return marker_genes, []
    drop_genes = noise_scores.head(drop_n).index.tolist()
    drop_set = set(drop_genes)
    kept = [g for g in marker_genes if g not in drop_set]
    return kept, drop_genes


# ----------------------------
# V4.2: 显著性检验和多重比较校正
# - 计算 p 值并进行 BH/Bonferroni 校正，为“灰区准入”提供统计依据。
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
    
    # Prefer permutation test when requested and inputs are available.
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
    # Using min spot p-value yields a conservative multi-spot correction.
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
        
        # Enforce monotonicity of BH-adjusted q-values.
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
    selection_strategy: str = "expression",
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
    # Filter by specificity, then fall back to a minimum marker count if needed.
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
    if selection_strategy not in {"expression", "specificity"}:
        selection_strategy = "expression"

    # ????????? >= min_specificity ???
    if min_specificity > 0:
        filtered_genes = [g for g in plugin_genes if specificity[g] >= min_specificity]
    else:
        filtered_genes = plugin_genes

    min_keep = max(2, int(min_marker_genes)) if min_marker_genes is not None else 2
    min_keep = min(min_keep, len(plugin_genes))
    specificity_rank = specificity
    if len(filtered_genes) < min_keep:
        # ?????????????? marker ?
        fallback_mode = fallback_specificity_mode or specificity_mode
        if fallback_mode == "max":
            fallback_spec = type_expr / (other_expr_max + eps)
        else:
            fallback_spec = type_expr / (other_expr_mean + eps)
        ranked = fallback_spec.sort_values(ascending=False)
        filtered_genes = ranked.index[:min_keep].tolist()
        specificity_rank = fallback_spec

    if len(filtered_genes) == 0:
        # ????????????????
        filtered_genes = plugin_genes

    if selection_strategy == "specificity":
        ranked = specificity_rank.loc[filtered_genes].sort_values(ascending=False)
        top_genes = ranked.index[:n_markers].tolist()
    else:
        # ??????????????? top n_markers
        type_expr_filtered = type_expr[filtered_genes]
        top_genes = type_expr_filtered.nlargest(n_markers).index.tolist()

    return top_genes


# ----------------------------
# 主流程
# - 加载配置与数据。
# - 计算支持度与显著性，并生成 Unknown/Drop/Relabel。
# - 应用 V5 拯救链，输出汇总与中间文件。
# ----------------------------

def main():
    args = parse_args()
    project_root = Path(args.project_root) if args.project_root else detect_project_root()

    # Config precedence: dataset config -> defaults, with CLI overrides.
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
    # Stage3 子配置（若缺省则为空字典）
    
    # 默认值
    # 统一的默认参数，用于缺省字段和向后兼容
    defaults = {
        "strong_th": 0.7,
        "weak_th": 0.4,
        "min_effect_size": 0.15,
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
        "marker_selection_strategy": "expression",
        "use_permutation_test": False,  # 是否使用置换检验（默认False，计算成本高）
        "z_score_aggregation": "max",  # z 值聚合方式
        "z_topk": 5,  # topk_mean 的 k 值
        "apply_mismatch_to_relabel": False,  # V4.2: 是否将显著缺失用于重写 plugin_type
        # V4.3 默认值
        "mismatch_action": "mark_unknown",
        "relabel_similarity_threshold": 0.8,
        "merge_target_map": None,
        "unknown_label_prefix": "Unknown_",
        "drop_unknown": True,
        "support_margin_min": None,
        "support_margin_enable": False,
        "conflict_demotion_enable": False,
        "protect_strong_from_missing": False,
    }

    # V5 系列子模块配置分组：去噪 / 生态位 / 熵质控 / 拯救控制
    v5_cfg = {}
    if isinstance(stage3_cfg, dict):
        v5_cfg = stage3_cfg.get("V5_denoising") or stage3_cfg.get("v5_denoising") or {}

    v5_niche_cfg = {}
    if isinstance(stage3_cfg, dict):
        v5_niche_cfg = stage3_cfg.get("V5_niche_rescue") or stage3_cfg.get("v5_niche_rescue") or {}

    v5_entropy_cfg = {}
    if isinstance(stage3_cfg, dict):
        v5_entropy_cfg = stage3_cfg.get("V5_entropy_qc") or stage3_cfg.get("v5_entropy_qc") or {}

    v5_rescue_ctrl_cfg = {}
    if isinstance(stage3_cfg, dict):
        v5_rescue_ctrl_cfg = stage3_cfg.get("V5_rescue_control") or stage3_cfg.get("v5_rescue_control") or {}

    # 将 V5 参数做类型转换与范围裁剪，避免非法值影响流程
    v5_enable = bool(v5_cfg.get("enable", False))
    try:
        v5_p_upper = float(v5_cfg.get("p_value_upper_limit", 0.2))
    except (TypeError, ValueError):
        v5_p_upper = 0.2
    if v5_p_upper > 1.0:
        v5_p_upper = 1.0
    try:
        v5_prune_ratio = float(v5_cfg.get("max_pruning_ratio", 0.3))
    except (TypeError, ValueError):
        v5_prune_ratio = 0.3
    try:
        v5_min_markers_left = int(v5_cfg.get("min_markers_left", 10))
    except (TypeError, ValueError):
        v5_min_markers_left = 10
    if v5_min_markers_left < 2:
        v5_min_markers_left = 2

    v5_niche_enable = bool(v5_niche_cfg.get("enable", False))
    try:
        v5_niche_anchor_p = float(v5_niche_cfg.get("anchor_p_threshold", 0.01))
    except (TypeError, ValueError):
        v5_niche_anchor_p = 0.01
    try:
        v5_niche_corr_th = float(v5_niche_cfg.get("correlation_threshold", 0.35))
    except (TypeError, ValueError):
        v5_niche_corr_th = 0.35
    v5_niche_metric = str(v5_niche_cfg.get("correlation_metric", "spearman")).strip().lower()
    if v5_niche_metric not in {"spearman", "pearson"}:
        v5_niche_metric = "spearman"
    try:
        v5_niche_p_upper = float(v5_niche_cfg.get("p_value_upper_limit", v5_p_upper))
    except (TypeError, ValueError):
        v5_niche_p_upper = v5_p_upper
    if v5_niche_p_upper > 1.0:
        v5_niche_p_upper = 1.0
    if v5_niche_p_upper < 0.0:
        v5_niche_p_upper = 0.0

    v5_entropy_enable = bool(v5_entropy_cfg.get("enable", False))
    try:
        v5_entropy_threshold = float(v5_entropy_cfg.get("entropy_threshold", 0.9))
    except (TypeError, ValueError):
        v5_entropy_threshold = 0.9
    try:
        v5_entropy_temperature = float(v5_entropy_cfg.get("temperature", 1.0))
    except (TypeError, ValueError):
        v5_entropy_temperature = 1.0
    if v5_entropy_temperature <= 0:
        v5_entropy_temperature = 1.0

    v5_rescue_ctrl_enable = bool(v5_rescue_ctrl_cfg.get("enable", False))
    try:
        v5_rescue_prior_weight = float(v5_rescue_ctrl_cfg.get("prior_weight", 1.0))
    except (TypeError, ValueError):
        v5_rescue_prior_weight = 1.0
    if v5_rescue_prior_weight < 0:
        v5_rescue_prior_weight = 0.0
    if v5_rescue_prior_weight > 1:
        v5_rescue_prior_weight = 1.0

    # 从配置读取（配置优先，然后 CLI 覆盖）
    # 读取 V4/V5 参数：dataset 配置优先，CLI 覆盖
    strong_th = args.strong_th if args.strong_th is not None else stage3_cfg.get("strong_th", defaults["strong_th"])
    weak_th = args.weak_th if args.weak_th is not None else stage3_cfg.get("weak_th", defaults["weak_th"])
    min_effect_size = (
        args.min_effect_size
        if args.min_effect_size is not None
        else stage3_cfg.get("min_effect_size", defaults["min_effect_size"])
    )
    st_cluster_k = args.st_cluster_k if args.st_cluster_k is not None else stage3_cfg.get("st_cluster_k", defaults["st_cluster_k"])
    unknown_floor = args.unknown_floor if args.unknown_floor is not None else stage3_cfg.get("unknown_floor", defaults["unknown_floor"])
    min_cells_rare_type = args.min_cells_rare_type if args.min_cells_rare_type is not None else stage3_cfg.get("min_cells_rare_type", defaults["min_cells_rare_type"])
    eps = args.eps if args.eps is not None else stage3_cfg.get("eps", defaults["eps"])
    try:
        min_effect_size = float(min_effect_size)
    except (TypeError, ValueError):
        min_effect_size = defaults["min_effect_size"]
    if min_effect_size < 0:
        min_effect_size = 0.0
    
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
    marker_selection_strategy = (
        args.marker_selection_strategy
        if args.marker_selection_strategy is not None
        else stage3_cfg.get("marker_selection_strategy", defaults["marker_selection_strategy"])
    )
    use_permutation_test = args.use_permutation_test if args.use_permutation_test is not None else stage3_cfg.get("use_permutation_test", defaults["use_permutation_test"])
    z_score_aggregation = args.z_score_aggregation if args.z_score_aggregation is not None else stage3_cfg.get("z_score_aggregation", defaults["z_score_aggregation"])
    z_topk = args.z_topk if args.z_topk is not None else stage3_cfg.get("z_topk", defaults["z_topk"])
    apply_mismatch_to_relabel = args.apply_mismatch_to_relabel if args.apply_mismatch_to_relabel is not None else stage3_cfg.get("apply_mismatch_to_relabel", defaults["apply_mismatch_to_relabel"])
    # V4.3 配置
    mismatch_action = args.mismatch_action if args.mismatch_action is not None else stage3_cfg.get("mismatch_action", defaults["mismatch_action"])
    relabel_similarity_threshold = (
        args.relabel_similarity_threshold
        if args.relabel_similarity_threshold is not None
        else stage3_cfg.get("relabel_similarity_threshold", defaults["relabel_similarity_threshold"])
    )
    merge_target_map = parse_merge_target_map(
        args.merge_target_map if args.merge_target_map is not None else stage3_cfg.get("merge_target_map", defaults["merge_target_map"])
    )
    unknown_label_prefix = args.unknown_label_prefix if args.unknown_label_prefix is not None else stage3_cfg.get("unknown_label_prefix", defaults["unknown_label_prefix"])
    drop_unknown = args.drop_unknown if args.drop_unknown is not None else stage3_cfg.get("drop_unknown", defaults["drop_unknown"])
    support_margin_min = (
        args.support_margin_min
        if args.support_margin_min is not None
        else stage3_cfg.get("support_margin_min", defaults["support_margin_min"])
    )
    support_margin_enable = (
        args.support_margin_enable
        if args.support_margin_enable is not None
        else stage3_cfg.get("support_margin_enable", defaults["support_margin_enable"])
    )
    conflict_demotion_enable = (
        args.conflict_demotion_enable
        if args.conflict_demotion_enable is not None
        else stage3_cfg.get("conflict_demotion_enable", defaults["conflict_demotion_enable"])
    )
    protect_strong_from_missing = (
        args.protect_strong_from_missing
        if args.protect_strong_from_missing is not None
        else stage3_cfg.get("protect_strong_from_missing", defaults["protect_strong_from_missing"])
    )

    if mismatch_action not in {"relabel", "mark_unknown", "ignore"}:
        mismatch_action = defaults["mismatch_action"]
    if relabel_similarity_threshold is None:
        relabel_similarity_threshold = defaults["relabel_similarity_threshold"]
    if unknown_label_prefix is None:
        unknown_label_prefix = defaults["unknown_label_prefix"]
    if drop_unknown is None:
        drop_unknown = defaults["drop_unknown"]

    if marker_specificity_mode not in {"mean", "max"}:
        marker_specificity_mode = defaults["marker_specificity_mode"]
    if fallback_specificity_mode not in {None, "mean", "max"}:
        fallback_specificity_mode = None
    if marker_selection_strategy not in {"expression", "specificity"}:
        marker_selection_strategy = defaults["marker_selection_strategy"]
    # 解析按类型覆盖的特异性/marker 下限配置
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
    # Stage1 导出与 Stage3 输出目录
    stage1_export = project_root / "data" / "processed" / args.sample / "stage1_preprocess" / "exported"
    out_proc = project_root / "data" / "processed" / args.sample / "stage3_typematch"
    out_res = project_root / "result" / args.sample / "stage3_typematch"
    out_proc.mkdir(parents=True, exist_ok=True)
    out_res.mkdir(parents=True, exist_ok=True)

    # 读取数据
    # 读取 Stage1 导出的表达矩阵与元数据
    sc_expr = pd.read_csv(stage1_export / "sc_expression_normalized.csv", index_col=0)
    st_expr = pd.read_csv(stage1_export / "st_expression_normalized.csv", index_col=0)
    sc_meta = pd.read_csv(stage1_export / "sc_metadata.csv")
    st_coords = pd.read_csv(stage1_export / "st_coordinates.csv", index_col=0)
    # 阶段二已移除：改用 Stage1 的高变基因或配置指定的基因列表
    # 选择用于匹配的基因集合（HVG 或显式路径）
    plugin_genes_path = stage3_cfg.get("plugin_genes_path") or (stage1_export.parent / "hvg_genes.txt")
    if not Path(plugin_genes_path).exists():
        raise FileNotFoundError(f"plugin_genes_path 不存在: {plugin_genes_path}")
    plugin_genes = [g.strip() for g in Path(plugin_genes_path).read_text(encoding="utf-8").splitlines() if g.strip()]
    # 可选基因权重，用于加权相似度
    gene_weights_path = stage3_cfg.get("gene_weights_path")
    if gene_weights_path and Path(gene_weights_path).exists():
        gene_weights = pd.read_csv(gene_weights_path)
        weight_map = dict(zip(gene_weights.iloc[:, 0], gene_weights.iloc[:, 1]))
    else:
        weight_map = {}
    w = np.array([weight_map.get(g, 1.0) for g in plugin_genes], dtype=float)

    # 对齐基因
    # 对齐 gene 列顺序，确保 SC/ST 使用同一集合
    sc_expr = sc_expr.loc[:, plugin_genes]
    st_expr = st_expr.loc[:, plugin_genes]

    # 取类型列名（兼容 cell_type / celltype / type），并加护栏避免 silent mismatch
    # 解析细胞类型列并检查不一致
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
    # 计算每个 cell type 的平均表达 profile，用于相似度与支持度评估
    profiles = compute_type_profiles(sc_expr, sc_meta, plugin_genes, type_col=type_col)

    # ST 聚类（可选）
    # 将 spot 聚类为若干代表性簇，提高稳定性与计算效率
    labels, st_profiles = try_cluster_st(st_expr, plugin_genes, k=st_cluster_k)
    st_profiles_np = st_profiles.to_numpy(dtype=float)

    # 支持度计算
    # 准备容器记录支持度、p 值、marker 等中间结果
    support_rows = []
    type_support_score = {}
    type_p_values = {}  # V4.2: 存储 p 值
    type_z_values = {}  # V4.2: 存储所有 z 值（用于计算 p 值）
    marker_genes_map: dict[str, List[str]] = {}
    type_cells_map: dict[str, List[str]] = {}
    
    # V4.1: 如果启用 Fisher r→z 方法
    # use_v42 代表启用 V4.2 的显著性检验与多重校正
    use_fisher_z = enable_mismatch_detection and score_method == "fisher_z"
    use_v42 = use_fisher_z  # V4.2 在 V4.1 基础上扩展
    if args.apply_mismatch_to_relabel is None and "apply_mismatch_to_relabel" not in stage3_cfg:
        apply_mismatch_to_relabel = use_v42
    if mismatch_action != "ignore" and not use_v42:
        print("[Stage3] V4.3 mismatch_action requires V4.2 results; fallback to ignore.")
        mismatch_action = "ignore"
    
    n_spots = len(st_expr.index)  # 用于 p 值计算的多 spot 校正
    
    # 逐类型计算支持度与显著性
    for t, prof in profiles.items():
        if use_fisher_z:
            # V4.1/V4.2: 使用 Fisher r→z 变换方法
            # 1. 选择标记基因（V4.2: 考虑特异性）
            effective_min_marker_specificity = specificity_override_map.get(
                str(t).strip().lower(), min_marker_specificity
            )
            use_specificity_ranking = marker_selection_strategy == "specificity"
            if effective_min_marker_specificity > 0 or use_specificity_ranking:
                effective_specificity_mode = override_map.get(str(t).strip().lower(), marker_specificity_mode)
                effective_min_marker_genes = min_genes_override_map.get(str(t).strip().lower(), min_marker_genes)
                marker_genes = select_marker_genes_with_specificity(
                    sc_expr, sc_meta, t, type_col, plugin_genes, 
                    marker_gene_count, effective_min_marker_specificity, effective_specificity_mode,
                    effective_min_marker_genes,
                    fallback_specificity_mode if (not fallback_specificity_types or str(t).strip().lower() in fallback_specificity_types) else None,
                    selection_strategy=marker_selection_strategy,
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
            marker_genes_map[t] = list(marker_genes)
            type_cells_map[t] = list(type_cells)
            
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
            # 传统方法：加权余弦相似度 + top-3 均值（兼容旧版本）
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
    # 基于 q 值决定显著缺失（Significant=Yes 表示缺失）
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
                    # 这里用 q 值（校正后）来判定缺失，避免多重比较带来的假阳性
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

    # V5 灰区定义：
    # 1) p 在 [alpha, upper_limit) 之间：统计不确定，进入拯救候选
    # 2) p < alpha 但支持度 <= min_effect_size：统计显著但效应太弱
    # upper_limit 由 V5 配置控制，用于限制“过于不显著”的类型
    def _in_v5_grey_zone(p_val: float, support_score: float, upper_limit: float) -> bool:
        if p_val is None:
            return False
        if p_val >= alpha and p_val < upper_limit:
            return True
        if p_val < alpha and support_score <= min_effect_size:
            return True
        return False

    # V5.1: 信号提纯拯救（只针对灰区 P 值）
    # 思路：对灰区类型做 marker 剪枝，重新计算支持度/显著性；若转为显著则拯救
    v5_denoising_rescued: set[str] = set()
    if v5_enable and use_v42:
        if v5_p_upper < alpha:
            v5_p_upper = alpha
        # 遍历所有类型，筛出灰区候选
        for row in support_rows:
            p_val = row.get("PValue", None)
            if p_val is None:
                continue
            try:
                p_val = float(p_val)
            except (TypeError, ValueError):
                continue
            if not _in_v5_grey_zone(p_val, row.get("support_score", 0.0), v5_p_upper):
                continue
            t = row["orig_type"]
            marker_genes = marker_genes_map.get(t, [])
            if not marker_genes:
                continue
            # 计算噪声并剪枝 marker
            cleaned_genes, dropped_genes = prune_noisy_markers(
                marker_genes,
                st_expr,
                max_pruning_ratio=v5_prune_ratio,
                min_markers_left=v5_min_markers_left,
                eps=eps,
            )
            if len(cleaned_genes) < 2:
                continue
            type_cells = type_cells_map.get(t, [])
            type_cells = [cid for cid in type_cells if cid in sc_expr.index]
            if len(type_cells) == 0:
                continue
            # 用剪枝后的 marker 重新评估支持度与 p 值
            type_profile = sc_expr.loc[type_cells, cleaned_genes].mean(axis=0)
            score_new, z_vals = compute_fisher_z_support_score(
                type_profile,
                st_expr,
                cleaned_genes,
                eps=eps,
                return_all_z=True,
                agg_method=z_score_aggregation,
                agg_topk=z_topk,
            )
            effective_use_permutation = use_permutation_test or z_score_aggregation != "max"
            p_new = compute_p_value(
                z_vals,
                len(cleaned_genes),
                n_spots,
                use_permutation=effective_use_permutation,
                type_profile=type_profile,
                st_expr=st_expr,
                marker_genes=cleaned_genes,
                eps=eps,
                score_aggregation=z_score_aggregation,
                z_topk=z_topk,
            )
            row["V5_denoising"] = True
            row["V5_p_value"] = p_new
            row["V5_support_score"] = score_new
            row["V5_marker_pruned"] = len(dropped_genes)
            row["V5_marker_left"] = len(cleaned_genes)
            row["V5_rescued"] = "Yes" if p_new < alpha else "No"
            # 记录是否被成功拯救
            if p_new < alpha:
                v5_denoising_rescued.add(t)

    # V5.2: niche co-occurrence rescue
    # 以高置信 anchor 类型为参照，基于 z 向量相关性判定拯救
    v5_niche_rescued: set[str] = set()
    v5_niche_anchor_map: dict[str, str] = {}
    v5_niche_corr_map: dict[str, float] = {}
    if v5_niche_enable and use_v42:
        if v5_niche_p_upper < alpha:
            v5_niche_p_upper = alpha
        anchor_types = []
        # 先筛出 anchor（极显著存在的类型）
        for row in support_rows:
            p_val = row.get("PValue", None)
            if p_val is None:
                continue
            try:
                p_val = float(p_val)
            except (TypeError, ValueError):
                continue
            if p_val < v5_niche_anchor_p:
                t = row["orig_type"]
                z_vals = type_z_values.get(t)
                if z_vals is not None and z_vals.size > 0:
                    anchor_types.append(t)
        # 再对灰区类型计算与 anchor 的相关性
        for row in support_rows:
            p_val = row.get("PValue", None)
            if p_val is None:
                continue
            try:
                p_val = float(p_val)
            except (TypeError, ValueError):
                continue
            if not _in_v5_grey_zone(p_val, row.get("support_score", 0.0), v5_niche_p_upper):
                continue
            t = row["orig_type"]
            if t in v5_denoising_rescued:
                continue
            z_c = type_z_values.get(t)
            if z_c is None or z_c.size == 0:
                continue
            best_anchor = None
            best_corr = -2.0
            for anchor in anchor_types:
                z_a = type_z_values.get(anchor)
                if z_a is None or z_a.size == 0:
                    continue
                if v5_niche_metric == "pearson":
                    corr = safe_pearson(z_c, z_a, eps=eps)
                else:
                    corr = safe_spearman(z_c, z_a, eps=eps)
                if corr > best_corr:
                    best_corr = corr
                    best_anchor = anchor
            if best_anchor is None:
                continue
            # 记录 anchor 与相关系数，用于审计与 summary
            row["V5_niche"] = True
            row["V5_niche_anchor"] = best_anchor
            row["V5_niche_corr"] = best_corr
            row["V5_niche_rescued"] = "Yes" if best_corr > v5_niche_corr_th else "No"
            if best_corr > v5_niche_corr_th:
                v5_niche_rescued.add(t)
                v5_niche_anchor_map[t] = best_anchor
                v5_niche_corr_map[t] = best_corr

    # V5.3: 熵质控（过滤分布过于弥散的拯救类型）
    # 熵越低表示信号集中（更可能真实），熵过高则撤销拯救
    v5_rescued_types = v5_denoising_rescued | v5_niche_rescued
    v5_final_rescued = set(v5_rescued_types)
    v5_entropy_values: dict[str, float] = {}
    v5_entropy_pass: dict[str, bool] = {}
    if v5_entropy_enable and use_v42 and v5_rescued_types:
        # 对每个已拯救类型计算归一化熵
        for t in sorted(v5_rescued_types):
            z_vals = type_z_values.get(t)
            if z_vals is None or z_vals.size == 0:
                v5_entropy_pass[t] = False
                continue
            entropy = compute_normalized_entropy(
                z_vals, temperature=v5_entropy_temperature, eps=eps
            )
            v5_entropy_values[t] = entropy
            v5_entropy_pass[t] = entropy < v5_entropy_threshold
        v5_final_rescued = {t for t in v5_rescued_types if v5_entropy_pass.get(t, False)}
    # 回填熵质控结果，便于审计与可视化
    for row in support_rows:
        t = row["orig_type"]
        if t in v5_entropy_values:
            row["V5_entropy"] = v5_entropy_values[t]
            row["V5_entropy_pass"] = "Yes" if v5_entropy_pass.get(t, False) else "No"
            row["V5_final_rescued"] = "Yes" if t in v5_final_rescued else "No"

    # 支持度分档 + 稀有类型保护
    # 根据 strong/weak/unsupported 将类型分层，稀有类型避免被判 strong
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
    # ???????????
    support_margin_min_value = None
    if support_margin_min is not None:
        try:
            support_margin_min_value = float(support_margin_min)
        except (TypeError, ValueError):
            support_margin_min_value = None
    for row in support_rows:
        t = row["orig_type"]
        margin = None
        if use_fisher_z:
            z_vals = type_z_values.get(t)
            if z_vals is not None and len(z_vals) >= 2:
                top2 = np.partition(z_vals, -2)[-2:]
                margin = float(np.max(top2) - np.min(top2))
        row["support_margin"] = margin
        if support_margin_min_value is not None and margin is not None:
            row["support_margin_pass"] = "Yes" if margin >= support_margin_min_value else "No"
        else:
            row["support_margin_pass"] = ""
        row["support_margin_downgraded"] = "No"
        row["conflict_downgraded"] = "No"
        if support_margin_enable and support_margin_min_value is not None and margin is not None:
            if row.get("support_category") == "strong" and margin < support_margin_min_value:
                row["support_category"] = "weak"
                row["support_margin_downgraded"] = "Yes"
        if conflict_demotion_enable and use_v42:
            if str(row.get("Significant", "")).strip().lower() == "yes" and row.get("support_category") == "strong":
                row["support_category"] = "weak"
                row["conflict_downgraded"] = "Yes"
    unsupported_types = {row["orig_type"] for row in support_rows if row["support_category"] == "unsupported"}
    # 拯救成功的类型不再视为 unsupported
    if v5_final_rescued:
        unsupported_types -= v5_final_rescued

    # V4.3: 根据显著缺失结果进行处理决策
    # missing_types 用于 relabel 或标记 Unknown
    missing_types = set()
    missing_conflicts = []
    if use_v42:
        for row in support_rows:
            if str(row.get("Significant", "")).strip().lower() == "yes":
                if protect_strong_from_missing and row.get("support_category") == "strong":
                    missing_conflicts.append(
                        {
                            "orig_type": row.get("orig_type"),
                            "support_category": row.get("support_category"),
                            "support_score": row.get("support_score"),
                            "PValue": row.get("PValue"),
                            "QValue": row.get("QValue"),
                        }
                    )
                    continue
                missing_types.add(row["orig_type"])
    # 拯救成功的类型不再视为缺失
    if v5_final_rescued:
        missing_types -= v5_final_rescued

    # 生成动作映射：Keep / Dropped / Unknown / Relabel
    # 初始全部 Keep，再应用 unsupported 与 missing_types 规则
    action_map = {row["orig_type"]: "Keep" for row in support_rows}
    for t in unsupported_types:
        action_map[t] = "Dropped" if drop_unknown else "Unknown"
    for t in v5_final_rescued:
        action_map[t] = "Keep"
    merge_sources: dict[str, List[str]] = defaultdict(list)
    similarity_target_map: dict[str, str] = {}

    # 缺失类型的处理策略：重标注或直接 Unknown
    # relabel 优先使用 merge_target_map 的人工映射，其次按表达相似度自动匹配
    if mismatch_action in {"relabel", "mark_unknown"} and missing_types:
        type_by_norm = {normalize_type_key(t): t for t in profiles.keys()}
        manual_map = {normalize_type_key(k): v for k, v in merge_target_map.items()}
        for miss in sorted(missing_types):
            if miss in unsupported_types:
                continue
            if mismatch_action == "mark_unknown":
                action_map[miss] = "Dropped" if drop_unknown else "Unknown"
                continue
            target = None
            manual_target = manual_map.get(normalize_type_key(miss))
            if manual_target:
                resolved = type_by_norm.get(normalize_type_key(manual_target))
                target = resolved or manual_target
                if target == miss or target in missing_types:
                    target = None
            if target is None:
                best_target = None
                best_sim = -1.0
                for cand in profiles.keys():
                    if cand == miss or cand in missing_types:
                        continue
                    sim = safe_pearson(profiles[miss], profiles[cand])
                    if sim > best_sim:
                        best_sim = sim
                        best_target = cand
                if best_target is not None and best_sim >= relabel_similarity_threshold:
                    target = best_target
                    similarity_target_map[miss] = best_target
            if target is not None:
                action_map[miss] = f"Relabel->{target}"
                merge_sources[target].append(miss)
            else:
                action_map[miss] = "Dropped" if drop_unknown else "Unknown"

    for target, sources in merge_sources.items():
        if not sources:
            continue
        merged = ",".join(sorted(sources))
        prev = action_map.get(target, "Keep")
        if prev.startswith("Merged<-"):
            prev_sources = prev.replace("Merged<-", "").strip()
            merged = ",".join(sorted({*prev_sources.split(","), *sources}))
        action_map[target] = f"Merged<-{merged}"

    # 将 action_map 写回 support_rows，供输出与后续分析
    for row in support_rows:
        row["Action"] = action_map.get(row["orig_type"], "Keep")

    # 写 type_support
    # 记录支持度/显著性/处理动作的明细表
    type_support_path = out_proc / "type_support.csv"
    pd.DataFrame(support_rows).to_csv(type_support_path, index=False)

    # 重写标签
    # 生成 cell_type_relabel.csv，提供给 Stage4 过滤与映射
    # 同时保留原始类型与 action，便于回溯
    support_map = {r["orig_type"]: r["support_category"] for r in support_rows}
    # V4.2: 将显著缺失结果应用到 relabel
    mismatch_map = {}
    if use_fisher_z and apply_mismatch_to_relabel:
        for row in support_rows:
            orig_type = row["orig_type"]
            flag_mismatch = row.get("FlagMismatch", "")
            mismatch_map[orig_type] = (flag_mismatch == "Yes")

    # 生成每个细胞的 plugin_type 与状态标签
    # 依赖 action_map / mismatch 结果决定是否 Unknown 或 Relabel
    relabel_rows = []
    for _, r in sc_meta.iterrows():
        orig = r[type_col]
        is_orig_unknown = str(orig).lower() in {"unknown", "unk", "na", "unlabeled"}
        status = support_map.get(orig, "unsupported")
        action = action_map.get(orig, "Keep")

        if mismatch_action in {"relabel", "mark_unknown"} and use_v42:
            if is_orig_unknown:
                plugin_type = "Unknown_sc_only"
                label_status = "unknown_merged"
            elif action.startswith("Relabel->"):
                plugin_type = action.split("->", 1)[1].strip()
                label_status = "v4.3_relabel"
            elif action in {"Unknown", "Dropped"}:
                plugin_type = "Unknown_sc_only"
                label_status = "v4.3_unknown"
            else:
                plugin_type = orig
                label_status = "kept"
        else:
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
                "action": action,
            }
        )
    relabel_df = pd.DataFrame(relabel_rows)
    relabel_path = out_proc / "cell_type_relabel.csv"
    relabel_df.to_csv(relabel_path, index=False)

    # 生成最终注释（可选择去除 Unknown）
    # 这里将 Relabel/Unknown 的最终标签写入 stage3_adjusted_annotations.csv
    adjusted_rows = []
    for row in relabel_rows:
        orig = row["orig_type"]
        action = action_map.get(orig, "Keep")
        if action.startswith("Relabel->"):
            final_type = action.split("->", 1)[1].strip()
        elif action in {"Unknown", "Dropped"}:
            final_type = f"{unknown_label_prefix}{orig}"
        else:
            final_type = orig
        # drop_unknown=True 时直接丢弃 Unknown/Dropped 细胞
        if drop_unknown and action in {"Unknown", "Dropped"}:
            continue
        adjusted_rows.append(
            {
                "cell_id": row["cell_id"],
                "cell_type_final": final_type,
                "orig_type": orig,
                "action": action,
            }
        )
    adjusted_df = pd.DataFrame(adjusted_rows)
    adjusted_path = out_proc / "stage3_adjusted_annotations.csv"
    adjusted_df.to_csv(adjusted_path, index=False)

    # plugin_type 列顺序（非 Unknown 按字母排序，Unknown_sc_only 最后）
    # rescued_plugin_types 用于后续 prior 的权重控制
    plugin_types = sorted([t for t in relabel_df["plugin_type"].unique() if t != "Unknown_sc_only"])
    plugin_types.append("Unknown_sc_only")
    rescued_plugin_types = [t for t in plugin_types if t in v5_final_rescued]

    # type_prior_matrix
    # 构建 spot × plugin_type 的先验分布，供 Stage4 使用
    # 注意：Unknown_sc_only 作为兜底类型参与归一化
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
            # V5 拯救控制：对被拯救类型下调 prior 权重，避免过拟合
            if v5_rescue_ctrl_enable and rescued_plugin_types and v5_rescue_prior_weight < 1.0:
                for idx, t in enumerate(plugin_types):
                    if t in rescued_plugin_types:
                        prior[idx] *= v5_rescue_prior_weight
                prior_sum = prior.sum()
                if prior_sum > eps:
                    prior = prior / prior_sum
            # Unknown 保底：保证 Unknown_sc_only 有最低占比
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
    # 汇总支持度分布、动作分布、Unknown 占比等
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

    def _action_bucket(action: str) -> str:
        if action.startswith("Relabel->"):
            return "Relabel"
        if action.startswith("Merged<-"):
            return "Merged"
        if action in {"Unknown", "Dropped"}:
            return action
        return "Keep"

    action_by_type = Counter(_action_bucket(v) for v in action_map.values())
    action_by_cell = Counter()
    for r in relabel_rows:
        action_by_cell[_action_bucket(action_map.get(r["orig_type"], "Keep"))] += 1

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

    # 汇总输出（stage3_summary.json）
    # 结构包含：params / support_overview / action_overview / v5_* / unknown_overview / rare_types / plugin_types
    unknown_label_effective = "Unknown_sc_only"
    output_paths = {
        "stage3_summary": str((out_res / "stage3_summary.json").relative_to(project_root)),
        "cell_type_relabel": str(relabel_path.relative_to(project_root)),
        "type_support": str(type_support_path.relative_to(project_root)),
        "stage3_adjusted_annotations": str(adjusted_path.relative_to(project_root)),
        "type_prior_matrix": str(prior_path.relative_to(project_root)),
    }

    summary = {
        # 参数快照：用于复现实验与排查差异
        "params": {
            "strong_th": strong_th,
            "weak_th": weak_th,
            "min_effect_size": min_effect_size,
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
            "marker_selection_strategy": marker_selection_strategy if use_v42 else None,
            "use_permutation_test": use_permutation_test if use_v42 else None,
            "z_score_aggregation": z_score_aggregation if use_v42 else None,
            "z_topk": z_topk if use_v42 else None,
            "apply_mismatch_to_relabel": apply_mismatch_to_relabel if use_fisher_z else None,
            # V4.3 参数
            "mismatch_action": mismatch_action if use_v42 else None,
            "relabel_similarity_threshold": relabel_similarity_threshold if use_v42 else None,
            "merge_target_map": merge_target_map if use_v42 else None,
            "unknown_label_prefix": unknown_label_prefix if use_v42 else None,
            "unknown_label_effective": unknown_label_effective if use_v42 else None,
            "drop_unknown": drop_unknown if use_v42 else None,
            "support_margin_min": support_margin_min if use_fisher_z else None,
            "support_margin_enable": support_margin_enable if use_fisher_z else None,
            "conflict_demotion_enable": conflict_demotion_enable if use_v42 else None,
            "protect_strong_from_missing": protect_strong_from_missing if use_v42 else None,
            # V5.1 参数
            "v5_denoising_enable": v5_enable if use_v42 else None,
            "v5_p_value_upper_limit": v5_p_upper if use_v42 else None,
            "v5_max_pruning_ratio": v5_prune_ratio if use_v42 else None,
            "v5_min_markers_left": v5_min_markers_left if use_v42 else None,
            "v5_niche_enable": v5_niche_enable if use_v42 else None,
            "v5_niche_anchor_p_threshold": v5_niche_anchor_p if use_v42 else None,
            "v5_niche_correlation_metric": v5_niche_metric if use_v42 else None,
            "v5_niche_correlation_threshold": v5_niche_corr_th if use_v42 else None,
            "v5_niche_p_value_upper_limit": v5_niche_p_upper if use_v42 else None,
            "v5_entropy_enable": v5_entropy_enable if use_v42 else None,
            "v5_entropy_threshold": v5_entropy_threshold if use_v42 else None,
            "v5_entropy_temperature": v5_entropy_temperature if use_v42 else None,
            "v5_rescue_control_enable": v5_rescue_ctrl_enable if use_v42 else None,
            "v5_rescue_prior_weight": v5_rescue_prior_weight if use_v42 else None,
        },
        # 支持度分布与细胞占比
        "support_overview": support_overview,
        # 动作汇总（Keep/Dropped/Unknown/Relabel）与缺失类型列表
        "action_overview": {
            "by_type_count": dict(action_by_type),
            "by_cell_count": dict(action_by_cell),
            "missing_types": sorted(missing_types),
            "missing_types_conflicts": missing_conflicts,
        },
        # V5.1 去噪拯救统计
        "v5_denoising": {
            "rescued_types": sorted(v5_denoising_rescued),
            "rescued_count": len(v5_denoising_rescued),
        },
        # V5.2 生态位拯救统计
        "v5_niche_rescue": {
            "rescued_types": sorted(v5_niche_rescued),
            "rescued_count": len(v5_niche_rescued),
            "anchor_map": v5_niche_anchor_map,
            "correlations": {k: v5_niche_corr_map[k] for k in sorted(v5_niche_corr_map)},
        },
        # V5.3 熵质控统计
        "v5_entropy_qc": {
            "enabled": v5_entropy_enable,
            "threshold": v5_entropy_threshold,
            "temperature": v5_entropy_temperature,
            "rescued_types_before": sorted(v5_rescued_types),
            "rescued_types_final": sorted(v5_final_rescued),
            "rejected_types": sorted(set(v5_rescued_types) - set(v5_final_rescued)),
            "entropy": {k: v5_entropy_values[k] for k in sorted(v5_entropy_values)},
        },
        # V5 拯救控制（prior 权重调整）
        "v5_rescue_control": {
            "enabled": v5_rescue_ctrl_enable,
            "prior_weight": v5_rescue_prior_weight,
            "rescued_types": rescued_plugin_types,
        },
        # Unknown 占比与 prior 统计
        "unknown_overview": {
            "cell_fraction": unknown_cells / total_cells if total_cells > 0 else 0.0,
            "prior_mass_fraction": unknown_prior_mass,
            "n_spots_all_unknown": n_spots_all_unknown,
            "unknown_label_effective": unknown_label_effective,
        },
        # 稀有类型清单与 plugin_type 列表（输出顺序）
        "rare_types": rare_types,
        "plugin_types": plugin_types,
        "output_paths": output_paths,
    }
    # 写入 summary 结果，作为后续指标与审计的权威来源
    with (out_res / "stage3_summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    print(f"[Stage3] type_support.csv -> {type_support_path}")
    print(f"[Stage3] cell_type_relabel.csv -> {relabel_path}")
    print(f"[Stage3] stage3_adjusted_annotations.csv -> {adjusted_path}")
    print(f"[Stage3] type_prior_matrix.csv -> {prior_path}")
    print(f"[Stage3] stage3_summary.json -> {out_res / 'stage3_summary.json'}")


if __name__ == "__main__":
    main()
