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
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import yaml


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
    }
    
    # 从配置读取（配置优先，然后 CLI 覆盖）
    strong_th = args.strong_th if args.strong_th is not None else stage3_cfg.get("strong_th", defaults["strong_th"])
    weak_th = args.weak_th if args.weak_th is not None else stage3_cfg.get("weak_th", defaults["weak_th"])
    st_cluster_k = args.st_cluster_k if args.st_cluster_k is not None else stage3_cfg.get("st_cluster_k", defaults["st_cluster_k"])
    unknown_floor = args.unknown_floor if args.unknown_floor is not None else stage3_cfg.get("unknown_floor", defaults["unknown_floor"])
    min_cells_rare_type = args.min_cells_rare_type if args.min_cells_rare_type is not None else stage3_cfg.get("min_cells_rare_type", defaults["min_cells_rare_type"])
    eps = args.eps if args.eps is not None else stage3_cfg.get("eps", defaults["eps"])

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
    for t, prof in profiles.items():
        sims = []
        for i in range(st_profiles_np.shape[0]):
            sims.append(weighted_cosine(prof, st_profiles_np[i], w, eps=eps))
        sims = np.array(sims, dtype=float)
        score = topk_mean(sims, k=3)
        type_support_score[t] = score
        best_idx = int(np.argmax(sims)) if sims.size > 0 else 0
        mapped_cluster = st_profiles.index[best_idx] if len(st_profiles.index) > 0 else ""
        n_cells = int((sc_meta[type_col] == t).sum())
        support_rows.append(
            {
                "orig_type": t,
                "n_cells": n_cells,
                "support_score": score,
                "support_category": "",  # 稍后填充
                "mapped_st_cluster": mapped_cluster,
            }
        )

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
    relabel_rows = []
    for _, r in sc_meta.iterrows():
        orig = r[type_col]
        is_orig_unknown = str(orig).lower() in {"unknown", "unk", "na", "unlabeled"}
        status = support_map.get(orig, "unsupported")
        if is_orig_unknown or status == "unsupported":
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

    summary = {
        "params": {
            "strong_th": strong_th,
            "weak_th": weak_th,
            "st_cluster_k": st_cluster_k,
            "similarity_metric": "weighted_cosine",
            "unknown_floor": unknown_floor,
            "min_cells_rare_type": min_cells_rare_type,
            "support_score_def": "top3_mean_cluster_similarity",
            "resolved_celltype_column": type_col,
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
