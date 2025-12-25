"""
Stage5 route2 S0 minimal evaluator:
- Computes leakage of missing type
- Computes spot-level composition errors (L1 / JS / corr)
- Writes a single JSON per run

Usage (single run):
  python src/stages/stage5_route2_s0.py \
    --sample real_brca_simS0_seed42 \
    --run_tag route2 \
    --stage4_dir result/real_brca_simS0_seed42/stage4_cytospace_filtered/cytospace_output \
    --sim_dir data/sim/real_brca/S0 \
    --out_dir result/real_brca_simS0_seed42/stage5_route2_s0/default

For baseline, pass run_tag=baseline and stage4_dir pointing to the baseline Stage4 output.
"""

from __future__ import annotations

import argparse
import json
import hashlib
import sys
from pathlib import Path
from typing import Tuple, List

import numpy as np
import pandas as pd

try:
    import yaml
except ImportError:
    yaml = None

# 添加项目根目录到路径，以便导入 src.utils
def detect_root() -> Path:
    here = Path(__file__).resolve()
    try:
        return here.parents[2]
    except IndexError:
        return here.parent.parent

_root = detect_root()
if str(_root) not in sys.path:
    sys.path.insert(0, str(_root))

from src.utils.type_name import normalize_type_name, load_alias_map, canonicalize_type_name

EPS = 1e-12


def sha1_file(p: Path) -> str:
    h = hashlib.sha1()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def normalize_rows(mat: np.ndarray) -> np.ndarray:
    row_sum = mat.sum(axis=1, keepdims=True)
    row_sum[row_sum == 0] = 1.0
    return mat / row_sum


def js_div(p: np.ndarray, q: np.ndarray) -> float:
    p = np.clip(p, EPS, 1.0)
    q = np.clip(q, EPS, 1.0)
    p = p / p.sum()
    q = q / q.sum()
    m = 0.5 * (p + q)
    return 0.5 * (np.sum(p * np.log(p / m)) + np.sum(q * np.log(q / m)))


def row_pearson(p: np.ndarray, q: np.ndarray) -> float:
    p0 = p - p.mean()
    q0 = q - q.mean()
    denom = np.sqrt((p0 * p0).sum()) * np.sqrt((q0 * q0).sum())
    if denom < EPS:
        return 0.0
    return float((p0 * q0).sum() / denom)


def load_by_spot_counts(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=None, engine="python")
    # rename first column to spot_id
    first_col = df.columns[0]
    df = df.rename(columns={first_col: "spot_id"})
    if "Total cells" in df.columns:
        df = df.drop(columns=["Total cells"])
    # clean spot_id (remove tab/space artifacts)
    df["spot_id"] = df["spot_id"].astype(str).str.split(r"\s|\t").str[0]
    df = df.set_index("spot_id")
    return df


def compute_composition_metrics(pred: pd.DataFrame, truth: pd.DataFrame) -> Tuple[float, float, float, int]:
    common_spots = truth.index.intersection(pred.index)
    truth = truth.loc[common_spots]
    pred = pred.loc[common_spots]

    common_types = [c for c in truth.columns if c in pred.columns]
    truth = truth[common_types].fillna(0.0)
    pred = pred[common_types].fillna(0.0)

    truth_np = normalize_rows(truth.to_numpy(dtype=float))
    pred_np = normalize_rows(pred.to_numpy(dtype=float))

    l1_list, js_list, corr_list = [], [], []
    for i in range(truth_np.shape[0]):
        p = pred_np[i]
        q = truth_np[i]
        l1_list.append(float(np.sum(np.abs(p - q))))
        js_list.append(js_div(p, q))
        corr_list.append(row_pearson(p, q))

    return float(np.mean(l1_list)), float(np.mean(js_list)), float(np.mean(corr_list)), len(common_spots)


def detect_root() -> Path:
    here = Path(__file__).resolve()
    try:
        return here.parents[2]
    except IndexError:
        return here.parent.parent


def load_stage3_weak_th(root: Path, sample: str) -> float | None:
    """从stage3_summary.json或配置文件读取weak_th"""
    # 优先从stage3_summary.json读取
    stage3_summary_path = root / "result" / sample / "stage3_typematch" / "stage3_summary.json"
    if stage3_summary_path.exists():
        try:
            with stage3_summary_path.open("r", encoding="utf-8") as f:
                summary = json.load(f)
            return summary.get("params", {}).get("weak_th")
        except Exception:
            pass
    
    # 从配置文件读取
    cfg_path = root / "configs" / "datasets" / f"{sample}.yaml"
    if cfg_path.exists() and yaml is not None:
        try:
            with cfg_path.open("r", encoding="utf-8") as f:
                data = yaml.safe_load(f)
            stage3 = (data or {}).get("stage3") or {}
            return stage3.get("weak_th")
        except Exception:
            pass
    
    return None


def eval_cell_spot_accuracy(
    truth_path: Path,
    pred_path: Path,
    missing_types: List[str],
    out_dir: Path,
    strict_config: bool = True,
    query_id: str | None = None,
    enable_spot_f1: bool = True,
    project_root: Path | None = None,
) -> tuple[dict, pd.DataFrame | None, pd.DataFrame | None]:
    """
    单细胞 spot 位置准确性评估（cell-level spot accuracy）
    
    Args:
        truth_path: truth_query_cell_spot.csv 路径
        pred_path: cell_assignment.csv 路径
        missing_types: missing/noise 类型列表（如 ["T cells CD8"]）
        out_dir: 输出目录
        strict_config: 是否启用严格检查
        query_id: 如果 truth 中有多个 query_id，指定要评估的 query_id
        enable_spot_f1: 是否计算 spot-level F1
        project_root: 项目根目录（用于加载 type_aliases.yaml）
    
    Returns:
        (summary_dict, by_type_df, spot_f1_df)
    """
    if not truth_path.exists():
        return None, None, None
    if not pred_path.exists():
        return None, None, None
    
    # 加载别名映射（用于类型规范化）
    if project_root is None:
        project_root = detect_root()
    alias_map_path = project_root / "configs" / "type_aliases.yaml"
    alias_map = load_alias_map(alias_map_path)
    
    # 读取 truth
    truth = pd.read_csv(truth_path)
    required_truth_cols = ["query_id", "cell_id", "true_spot_id", "cell_type"]
    missing_cols = [c for c in required_truth_cols if c not in truth.columns]
    if missing_cols:
        raise ValueError(f"truth_query_cell_spot.csv missing columns: {missing_cols}")
    
    # 读取 pred
    pred = pd.read_csv(pred_path)
    required_pred_cols = ["cell_id", "assigned_spot", "cell_type"]
    missing_cols = [c for c in required_pred_cols if c not in pred.columns]
    if missing_cols:
        raise ValueError(f"cell_assignment.csv missing columns: {missing_cols}")
    
    # 严格检查1: truth (query_id, cell_id) 唯一性
    truth_dup = truth.groupby(["query_id", "cell_id"]).size()
    if (truth_dup > 1).any():
        dup_count = int((truth_dup > 1).sum())
        msg = f"truth_query_cell_spot.csv: {dup_count} duplicate (query_id, cell_id) pairs"
        if strict_config:
            raise ValueError(msg)
        else:
            print(f"[WARN] {msg}")
    
    # 严格检查2: 处理多个 query_id 的情况
    unique_query_ids = truth["query_id"].unique()
    if len(unique_query_ids) > 1:
        # 检查是否每个 cell_id 都有唯一的 query_id
        cell_id_counts = truth.groupby("cell_id")["query_id"].nunique()
        has_overlap = (cell_id_counts > 1).any()
        
        if has_overlap:
            # 有 cell_id 对应多个 query_id，合并所有数据但去重 cell_id（保留第一次出现）
            # 这样我们可以使用所有 cell 的真值，而不仅仅是第一个 query_id
            print(f"[INFO] truth_query_cell_spot.csv contains {len(unique_query_ids)} query_ids with overlapping cell_id. Merging all query_ids and deduplicating by cell_id.")
            query_id = f"merged_{len(unique_query_ids)}_queries"
            # 去重 cell_id，保留第一次出现（这样可以保留所有唯一的 cell）
            truth = truth.drop_duplicates(subset=["cell_id"], keep="first").copy()
        else:
            # 每个 cell_id 只有一个 query_id，可能是索引，合并所有 query_id
            print(f"[INFO] truth_query_cell_spot.csv contains {len(unique_query_ids)} query_ids, but each cell_id has unique query_id. Merging all query_ids.")
            query_id = f"merged_{len(unique_query_ids)}_queries"
            # 不需要过滤，直接使用所有数据
    elif len(unique_query_ids) == 1:
        query_id = unique_query_ids[0]
    
    # 严格检查3: pred cell_id 唯一性
    pred_dup = pred["cell_id"].duplicated()
    if pred_dup.any():
        dup_count = int(pred_dup.sum())
        # 对于 cell_assignment，一个 cell 可能被分配到多个 spot（在 CytoSPACE 中可能发生）
        # 我们选择保留第一个分配（或可以改为保留所有，但需要修改后续逻辑）
        print(f"[INFO] cell_assignment.csv: {dup_count} duplicate cell_id (one cell assigned to multiple spots). Keeping first assignment.")
        pred = pred.drop_duplicates(subset=["cell_id"], keep="first")
    
    # 规范化类型名
    truth["cell_type_norm"] = truth["cell_type"].apply(lambda x: normalize_type_name(x))
    pred["cell_type_norm"] = pred["cell_type"].apply(lambda x: normalize_type_name(x))
    
    # 规范化 missing_types
    missing_types_norm = [normalize_type_name(mt) for mt in missing_types]
    
    # 核心 join（outer join 以统计缺失）
    merged = truth.merge(
        pred,
        on="cell_id",
        how="outer",
        suffixes=("_truth", "_pred")
    )
    
    # 标记列
    merged["pred_present"] = merged["assigned_spot"].notna()
    # truth_present: true_spot_id 非空（对于 non-missing cells，这应该总是 true）
    merged["truth_present"] = merged["true_spot_id"].notna()
    merged["type_mismatch"] = (
        merged["cell_type_norm_truth"] != merged["cell_type_norm_pred"]
    ) & merged["pred_present"] & merged["truth_present"]
    
    # 标记 missing/noise 类型（以 truth 为准）
    merged["is_missing"] = merged["cell_type_norm_truth"].isin(missing_types_norm)
    
    # 对于 missing cells，true_spot_id 可能为空（这是正常的）
    # 但我们需要确保在统计时正确处理
    
    # 严格检查4: 大量 cell_id 不在 truth（可能用错文件）
    pred_only = merged[merged["truth_present"].isna()]
    if len(pred_only) > 0:
        pred_only_frac = len(pred_only) / len(pred)
        if pred_only_frac > 0.1:  # 超过 10%
            msg = f"{pred_only_frac*100:.1f}% of prediction cell_id not in truth (possible wrong file)"
            if strict_config:
                raise ValueError(msg)
            else:
                print(f"[WARN] {msg}")
    
    # 严格检查5: type_mismatch 比例
    type_mismatch_frac = merged["type_mismatch"].sum() / len(merged[merged["pred_present"] & merged["truth_present"]]) if (merged["pred_present"] & merged["truth_present"]).any() else 0.0
    if type_mismatch_frac > 0.01:  # 超过 1%
        msg = f"type_mismatch_fraction = {type_mismatch_frac:.4f} > 0.01 (possible type column confusion)"
        if strict_config:
            raise ValueError(msg)
        else:
            print(f"[WARN] {msg}")
    
    # 计算覆盖率
    n_truth_total = len(truth)
    n_pred_present = int(merged["pred_present"].sum())
    # n_pred_present_in_truth: pred 中在 truth 中存在的 cell_id 数量
    n_pred_present_in_truth = int((merged["pred_present"] & merged["truth_present"]).sum())
    pred_cell_id_unique = len(pred) == pred["cell_id"].nunique()
    truth_cell_id_unique = len(truth) == truth["cell_id"].nunique()
    
    # 打印关键统计信息
    print(f"[cell_spot_eval_stats] n_truth_total = {n_truth_total}")
    print(f"[cell_spot_eval_stats] n_pred_present = {n_pred_present}")
    print(f"[cell_spot_eval_stats] n_pred_present_in_truth = {n_pred_present_in_truth}")
    print(f"[cell_spot_eval_stats] pred.cell_id unique = {pred_cell_id_unique} (n_pred={len(pred)}, n_unique={pred['cell_id'].nunique()})")
    print(f"[cell_spot_eval_stats] truth.cell_id unique = {truth_cell_id_unique} (n_truth={len(truth)}, n_unique={truth['cell_id'].nunique()})")
    
    # 自检1: 检查 extra_pred (pred中有但truth中没有的cell) 的 cell_type 分布
    # truth_present是布尔值（notna()的结果），所以应该用~merged["truth_present"]而不是isna()
    extra_pred_mask = merged["pred_present"] & (~merged["truth_present"])
    n_extra_pred = int(extra_pred_mask.sum())
    if n_extra_pred > 0:
        extra_pred_df = merged[extra_pred_mask].copy()
        extra_pred_type_counts = extra_pred_df["cell_type_norm_pred"].value_counts()
        print(f"[self_check_1] extra_pred count = {n_extra_pred} (pred中有但truth中没有的cell)")
        print(f"[self_check_1] extra_pred cell_type 分布:")
        for cell_type, count in extra_pred_type_counts.items():
            print(f"  {cell_type}: {count}")
        # 确认不包含CD8
        missing_types_in_extra = [mt for mt in missing_types_norm if mt in extra_pred_type_counts.index]
        if missing_types_in_extra:
            print(f"[self_check_1] ⚠️ WARNING: extra_pred中包含missing_types: {missing_types_in_extra}")
        else:
            print(f"[self_check_1] ✅ OK: extra_pred中不包含missing_types (CD8等)")
    else:
        print(f"[self_check_1] extra_pred count = 0 (所有pred都在truth中)")
    
    # 改1: 删除 coverage_all，新增 coverage_on_truth 和 extra_pred_fraction
    # coverage_all = n_pred_present / n_truth_total 在 truth 是 query 子集时会 >100%，容易误解
    # 改为：coverage_on_truth = n_pred_present_in_truth / n_truth_total（交集覆盖率）
    coverage_on_truth = n_pred_present_in_truth / n_truth_total if n_truth_total > 0 else 0.0
    extra_pred_fraction = (n_pred_present - n_pred_present_in_truth) / n_pred_present if n_pred_present > 0 else 0.0
    
    # non-missing cells
    non_missing_mask = ~merged["is_missing"]
    n_truth_nonmissing = int(non_missing_mask.sum())
    n_pred_nonmissing_present = int((non_missing_mask & merged["pred_present"]).sum())
    coverage_nonmissing = n_pred_nonmissing_present / n_truth_nonmissing if n_truth_nonmissing > 0 else 0.0
    
    # 改2: cell-level 模块里不再输出 leak_rate（truth_query 不含 missing-type）
    # missing cells 的 leak_rate 只在 Stage5 主输出中计算（基于 by_spot 统计）
    # 这里保留统计量计算，但 leak_rate/reject_rate 设为 None
    missing_mask = merged["is_missing"]
    n_truth_missing = int(missing_mask.sum())
    n_pred_missing_present = int((missing_mask & merged["pred_present"]).sum())
    leak_rate = None  # truth_query_cell_spot 不含 missing-type，泄漏评估见 Stage5 主输出
    reject_rate = None  # 同上
    
    # 单细胞 spot Top-1 命中率（只对 non-missing 算）
    # 自检2: 确保只使用交集（truth_query cells），即 truth_present=True 的cell
    non_missing_merged = merged[non_missing_mask & merged["truth_present"]].copy()
    
    # 验证：Acc@1计算只用了truth中存在的cell
    n_non_missing_in_truth = len(non_missing_merged)
    print(f"[self_check_2] Acc@1计算使用的cell数 = {n_non_missing_in_truth} (应该等于n_truth_nonmissing，且只包含truth中存在的cell)")
    print(f"[self_check_2] ✅ 确认: Acc@1只基于truth_query cells (truth_present=True)，不包含extra_pred")
    
    # 清理 spot_id（去除可能的 tab/space artifacts）
    non_missing_merged["true_spot_id_clean"] = non_missing_merged["true_spot_id"].astype(str).str.split(r"\s|\t").str[0]
    non_missing_merged["assigned_spot_clean"] = non_missing_merged["assigned_spot"].astype(str).str.split(r"\s|\t").str[0]
    
    # 计算 hit
    non_missing_merged["hit"] = (
        (non_missing_merged["true_spot_id_clean"] == non_missing_merged["assigned_spot_clean"])
        & non_missing_merged["pred_present"]
    )
    
    # 口径1: 只在"有预测"的 non-missing 上算
    non_missing_with_pred = non_missing_merged[non_missing_merged["pred_present"]]
    n_non_missing_with_pred = len(non_missing_with_pred)
    print(f"[self_check_2] acc_top1_cond计算使用的cell数 = {n_non_missing_with_pred} (有预测的non-missing cells)")
    acc_top1_cond = float(non_missing_with_pred["hit"].mean()) if len(non_missing_with_pred) > 0 else 0.0
    
    # 口径2: 在全部 non-missing 上算（未分配算错）
    non_missing_merged["hit_all"] = non_missing_merged["hit"].fillna(False)
    acc_top1_all = float(non_missing_merged["hit_all"].mean()) if len(non_missing_merged) > 0 else 0.0
    
    # 分层准确性（按 cell_type_truth）
    by_type_rows = []
    for cell_type in non_missing_merged["cell_type_norm_truth"].dropna().unique():
        type_mask = non_missing_merged["cell_type_norm_truth"] == cell_type
        type_merged = non_missing_merged[type_mask]
        n_truth_type = len(type_merged)
        n_pred_present_type = int(type_merged["pred_present"].sum())
        coverage_type = n_pred_present_type / n_truth_type if n_truth_type > 0 else 0.0
        
        type_with_pred = type_merged[type_merged["pred_present"]]
        acc_top1_cond_type = float(type_with_pred["hit"].mean()) if len(type_with_pred) > 0 else 0.0
        
        acc_top1_all_type = float(type_merged["hit_all"].mean()) if len(type_merged) > 0 else 0.0
        
        by_type_rows.append({
            "cell_type": cell_type,
            "n_truth": n_truth_type,  # truth_query 中该类型的 cell 数
            "n_pred_present": n_pred_present_type,  # pred 中该类型且有预测的 cell 数（交集内）
            "coverage": coverage_type,  # 交集覆盖率：n_pred_present / n_truth
            "acc_top1_cond": acc_top1_cond_type,  # 仅在 truth_query 交集上算
            "acc_top1_all": acc_top1_all_type,  # 同上
        })
    
    by_type_df = pd.DataFrame(by_type_rows)
    
    # Spot-level F1（可选）
    spot_f1_df = None
    macro_f1 = None
    micro_f1 = None
    
    if enable_spot_f1:
        spot_rows = []
        all_spots = set(non_missing_merged["true_spot_id_clean"].dropna().unique()) | set(non_missing_merged["assigned_spot_clean"].dropna().unique())
        
        tp_total = 0
        fp_total = 0
        fn_total = 0
        
        for spot_id in all_spots:
            truth_cells = set(non_missing_merged[non_missing_merged["true_spot_id_clean"] == spot_id]["cell_id"].dropna())
            pred_cells = set(non_missing_merged[non_missing_merged["assigned_spot_clean"] == spot_id]["cell_id"].dropna())
            
            n_truth_cells = len(truth_cells)
            n_pred_cells = len(pred_cells)
            
            intersection = truth_cells & pred_cells
            tp = len(intersection)
            fp = len(pred_cells - truth_cells)
            fn = len(truth_cells - pred_cells)
            
            tp_total += tp
            fp_total += fp
            fn_total += fn
            
            precision = tp / n_pred_cells if n_pred_cells > 0 else 0.0
            recall = tp / n_truth_cells if n_truth_cells > 0 else 0.0
            f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
            
            spot_rows.append({
                "spot_id": spot_id,
                "n_truth_cells": n_truth_cells,
                "n_pred_cells": n_pred_cells,
                "precision": precision,
                "recall": recall,
                "f1": f1,
            })
        
        spot_f1_df = pd.DataFrame(spot_rows)
        
        # Macro F1
        macro_f1 = float(spot_f1_df["f1"].mean()) if len(spot_f1_df) > 0 else 0.0
        
        # Micro F1
        micro_f1 = (2 * tp_total) / (2 * tp_total + fp_total + fn_total) if (2 * tp_total + fp_total + fn_total) > 0 else 0.0
    
    # 保存 CSV
    by_type_csv = out_dir / "cell_spot_acc_by_type.csv"
    by_type_df.to_csv(by_type_csv, index=False, encoding="utf-8")
    
    spot_f1_csv = None
    if spot_f1_df is not None:
        spot_f1_csv = out_dir / "spot_f1.csv"
        spot_f1_df.to_csv(spot_f1_csv, index=False, encoding="utf-8")
    
    # 构建 summary dict（确保所有值都是 JSON 可序列化的）
    summary = {
        "inputs": {
            "truth_query_cell_spot_path": str(truth_path),
            "cell_assignment_path": str(pred_path),
            "query_id": str(query_id) if query_id is not None else None,
            "missing_types": [str(mt) for mt in missing_types_norm],
            "eval_mode": "hard_assignment",
        },
        "counts": {
            "n_truth_total": int(n_truth_total),
            "n_pred_present": int(n_pred_present),
            "n_pred_present_in_truth": int(n_pred_present_in_truth),  # 交集：pred 中在 truth 中存在的 cell 数
            "n_extra_pred": int(n_extra_pred),  # pred 中有但 truth 中没有的 cell 数（不参与 Acc@1）
            "n_truth_missing": int(n_truth_missing),
            "n_pred_missing_present": int(n_pred_missing_present),
            "n_truth_nonmissing": int(n_truth_nonmissing),
            "n_pred_nonmissing_present": int(n_pred_nonmissing_present),
        },
        "metrics": {
            # 改1: 删除 coverage_all，新增 coverage_on_truth 和 extra_pred_fraction
            "coverage_on_truth": float(coverage_on_truth),
            "extra_pred_fraction": float(extra_pred_fraction),
            "coverage_nonmissing": float(coverage_nonmissing),
            # 改2: cell-level 模块里不再输出 leak_rate（truth_query 不含 missing-type）
            "leak_rate": None,  # truth_query_cell_spot 不含 missing-type，泄漏评估见 Stage5 主输出
            "reject_rate": None,  # 同上
            # 改3: 明确 Acc@1 仅在 truth_query 交集上算
            "acc_top1_cond_nonmissing": float(acc_top1_cond),  # 使用 n_pred_present_in_truth 个 truth_query cells（交集）
            "acc_top1_all_nonmissing": float(acc_top1_all),  # 同上
        },
        "validation_checks": {
            "truth_unique_cell_id": len(truth) == truth["cell_id"].nunique(),
            "pred_unique_cell_id": len(pred) == pred["cell_id"].nunique(),
            "single_query_id": len(unique_query_ids) == 1,
            "type_mismatch_fraction": float(type_mismatch_frac),
        },
        "note": {
            "acc_top1_calculation": "Acc@1 仅在 truth_query cells 交集上计算（n_pred_present_in_truth 个 cell），不包含 extra_pred",
            "leak_rate_note": "truth_query_cell_spot 不含 missing-type，泄漏评估见 Stage5 主输出（leakage.missing_leak_rate）",
        },
        "artifacts": {
            "by_cell_type_csv": str(by_type_csv.relative_to(out_dir)) if by_type_csv.exists() else None,
            "spot_f1_csv": str(spot_f1_csv.relative_to(out_dir)) if spot_f1_csv and spot_f1_csv.exists() else None,
        },
    }
    
    if macro_f1 is not None:
        summary["metrics"]["macro_F1_spot"] = float(macro_f1)
    if micro_f1 is not None:
        summary["metrics"]["micro_F1_spot"] = float(micro_f1)
    
    return summary, by_type_df, spot_f1_df


def generate_filter_audit(
    sample: str,
    missing_type: str,
    run_tag: str,
    n_filtered_total: int,
    n_sc_total: int,
    out_dir: Path,
    filter_scope: str | None = None,
) -> dict | None:
    """生成filter_audit：统计被过滤类型的构成

    filter_scope:
        - None / "unsupported_all": 使用所有 plugin_type == Unknown_sc_only 的细胞
        - "missing_only": 仅统计 orig_type == missing_type 且 plugin_type == Unknown_sc_only 的细胞

    Returns:
        dict: filter_audit摘要（JSON中保留）
        CSV文件会保存到out_dir/filter_audit_all_types.csv
    """
    root = detect_root()
    relabel_path = root / "data" / "processed" / sample / "stage3_typematch" / "cell_type_relabel.csv"
    
    if not relabel_path.exists():
        return None
    
    relabel = pd.read_csv(relabel_path)

    # 找出所有被过滤的细胞（与Stage4 filter_scope保持一致）
    mask_unknown = relabel["plugin_type"] == "Unknown_sc_only"
    if filter_scope == "missing_only":
        # missing_only 模式：直接过滤所有 orig_type == missing_type 的细胞，不管 plugin_type 是什么
        # 这与 Stage4 的逻辑保持一致：missing_only 模式下直接过滤 missing_type，不依赖 base_mask
        mask_type = relabel["orig_type"] == missing_type
        filtered = relabel[mask_type].copy()
    else:
        # unsupported_all 模式：只过滤 plugin_type == Unknown_sc_only 的细胞
        filtered = relabel[mask_unknown].copy()

    # 实际filtered数量，用于和Stage4 n_filtered对比
    actual_filtered = int(len(filtered))

    # 统计被过滤细胞的orig_type分布
    type_counts = filtered["orig_type"].value_counts().to_dict()
    
    # 计算missing_type的贡献
    missing_count = type_counts.get(missing_type, 0)
    missing_fraction_of_filtered = missing_count / n_filtered_total if n_filtered_total > 0 else 0.0
    
    # 计算非missing类型的误伤
    non_missing_filtered = n_filtered_total - missing_count
    non_missing_fraction_of_total = non_missing_filtered / n_sc_total if n_sc_total > 0 else 0.0
    
    # Top10被过滤类型（基于实际统计到的filtered集合）
    top10_types = list(type_counts.items())[:10]
    denom = float(n_filtered_total) if n_filtered_total > 0 else 1.0
    top10_list = [{"type": t, "count": int(c), "fraction": float(c/denom)} for t, c in top10_types]

    # 强制检查1：sum(top_types.n_filtered) == n_filtered_total
    sum_top10 = sum(c for _, c in top10_types)
    check1_passed = abs(sum_top10 - n_filtered_total) < 1e-6
    
    # 强制检查2：missing_type_contrib
    check2_value = missing_fraction_of_filtered
    
    # 强制检查3：non_missing_filtered_fraction
    check3_value = non_missing_fraction_of_total
    
    # 读取weak_th
    weak_th_effective = load_stage3_weak_th(root, sample)

    # 额外统计：Stage3 视角下所有被标记为 Unknown_sc_only 的规模（不论是否被Stage4丢弃）
    marked_total = int(mask_unknown.sum())
    marked_missing = int((mask_unknown & (relabel["orig_type"] == missing_type)).sum())
    marked_non_missing = marked_total - marked_missing
    
    # 保存全量all_filtered_types到CSV
    denom = float(n_filtered_total) if n_filtered_total > 0 else 1.0
    all_types_df = pd.DataFrame([
        {"type": t, "count": int(c), "fraction": float(c/denom)}
        for t, c in sorted(type_counts.items(), key=lambda x: x[1], reverse=True)
    ])
    csv_path = out_dir / "filter_audit_all_types.csv"
    all_types_df.to_csv(csv_path, index=False, encoding="utf-8")
    
    # 计算相对路径
    try:
        csv_path_rel = str(csv_path.relative_to(root))
    except ValueError:
        # 如果不在同一路径下，使用绝对路径
        csv_path_rel = str(csv_path)
    
    return {
        "weak_th_effective": weak_th_effective,
        "missing_type": missing_type,
        "run_tag": run_tag,
        "source_counts_basis": "sc_total_before_prefilter",
        "total_filtered": int(n_filtered_total),
        "total_sc_cells": int(n_sc_total),
        "missing_type_contribution": {
            "count": int(missing_count),
            "fraction_of_filtered": float(missing_fraction_of_filtered),
            "fraction_of_total": float(missing_count / n_sc_total) if n_sc_total > 0 else 0.0
        },
        "marked_unknown": {
            "total": int(marked_total),
            "missing": int(marked_missing),
            "non_missing": int(marked_non_missing),
            "non_missing_fraction_of_total": float(marked_non_missing / n_sc_total) if n_sc_total > 0 else 0.0,
        },
        "non_missing_filtered": {
            "count": int(non_missing_filtered),
            "fraction_of_filtered": float(non_missing_filtered / n_filtered_total) if n_filtered_total > 0 else 0.0,
            "fraction_of_total": float(non_missing_fraction_of_total)
        },
        "top10_filtered_types": top10_list,
        "filter_audit_path": csv_path_rel if csv_path.exists() else None,
        "validation_checks": {
            "check1_sum_top10_equals_total": {
                "passed": check1_passed,
                "sum_top10": int(sum_top10),
                "n_filtered_total": int(n_filtered_total)
            },
            "check2_missing_type_contrib": {
                "value": float(check2_value),
                "description": "missing_type_contrib = n_filtered(missing_type)/n_filtered_total"
            },
            "check3_non_missing_filtered_fraction": {
                "value": float(check3_value),
                "description": "non_missing_filtered_fraction = (n_filtered_total - n_filtered(missing_type))/n_sc_total"
            },
            "check4_stage4_stage5_n_filtered_match": {
                "value": bool(actual_filtered == n_filtered_total),
                "actual_filtered": int(actual_filtered),
                "stage4_n_filtered": int(n_filtered_total),
                "description": "Ensure Stage4 n_filtered equals count of filtered rows reconstructed in Stage5"
            },
        }
    }


def run_eval(args: argparse.Namespace):
    stage4_dir = Path(args.stage4_dir)
    sim_dir = Path(args.sim_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_p = stage4_dir / "stage4_summary.json"
    assign_p = stage4_dir / "cell_assignment.csv"
    by_spot_p = stage4_dir / "cell_type_assignments_by_spot.csv"

    sim_info_p = sim_dir / "sim_info.json"
    truth_spot_p = sim_dir / "sim_truth_spot_type_fraction.csv"
    truth_query_p = sim_dir / "sim_truth_query_cell_spot.csv"

    with summary_p.open("r", encoding="utf-8") as f:
        s4 = json.load(f)

    with sim_info_p.open("r", encoding="utf-8") as f:
        sim_info = json.load(f)

    missing_type = sim_info.get("missing_type")

    ca = pd.read_csv(assign_p)
    leak_assign = int((ca["cell_type"] == missing_type).sum()) if "cell_type" in ca.columns else None

    by_spot = load_by_spot_counts(by_spot_p)
    leak_by_spot = int(by_spot[missing_type].sum()) if missing_type in by_spot.columns else 0

    assignment_rows = int(s4.get("assignment_rows", len(ca)))
    leak_rate = (leak_by_spot / assignment_rows) if assignment_rows else 0.0

    truth = load_by_spot_counts(truth_spot_p)
    l1_mean, js_mean, corr_mean, n_spots = compute_composition_metrics(by_spot, truth)

    # 生成filter_audit（仅对route2运行）
    filter_audit = None
    if args.run_tag == "route2":
        n_filtered = s4.get("n_filtered", 0)
        n_sc_total = s4.get("n_cells_before_prefilter", 0)
        filter_audit = generate_filter_audit(
            args.sample,
            missing_type,
            args.run_tag,
            n_filtered,
            n_sc_total,
            out_dir,
            filter_scope=s4.get("filter_scope"),
        )
        # 严格模式下，若Stage4与Stage5重建的filtered数量不一致则直接报错
        if args.strict_config and filter_audit is not None:
            vc = (filter_audit.get("validation_checks") or {}).get("check4_stage4_stage5_n_filtered_match")
            if isinstance(vc, dict) and not vc.get("value", True):
                msg = (
                    "[strict_config] Stage4 vs Stage5 n_filtered mismatch: "
                    f"Stage4={vc.get('stage4_n_filtered')} vs Stage5_reconstructed={vc.get('actual_filtered')}"
                )
                raise ValueError(msg)

    # 单细胞 spot 位置准确性评估
    cell_spot_eval = None
    if truth_query_p.exists():
        try:
            # 对于 route2，type_mismatch 是正常的（因为使用了 plugin_type），所以放宽检查
            eval_strict = args.strict_config and args.run_tag == "baseline"
            cell_spot_eval, _, _ = eval_cell_spot_accuracy(
                truth_path=truth_query_p,
                pred_path=assign_p,
                missing_types=[missing_type] if missing_type else [],
                out_dir=out_dir,
                strict_config=eval_strict,
                query_id=None,  # 自动检测
                enable_spot_f1=True,
                project_root=detect_root(),
            )
        except Exception as e:
            if args.strict_config and args.run_tag == "baseline":
                raise
            else:
                print(f"[WARN] cell_spot_accuracy evaluation failed: {e}")
                cell_spot_eval = None

    out = {
        "scenario": "S0",
        "sample": args.sample,
        "run_tag": args.run_tag,
        "missing_type_truth": missing_type,
        "leakage": {
            "missing_in_assignment": leak_assign,
            "missing_in_by_spot": leak_by_spot,
            "missing_leak_rate": leak_rate,
        },
        "composition": {
            "L1_mean": l1_mean,
            "JS_mean": js_mean,
            "corr_mean": corr_mean,
            "n_spots_eval": n_spots,
        },
        "coverage": {
            "n_cells_before_prefilter": s4.get("n_cells_before_prefilter"),
            "n_cells_after_prefilter": s4.get("n_cells_after_prefilter"),
            "n_filtered": s4.get("n_filtered"),
        },
        "filter_audit": filter_audit,
        "scale": {
            "spots": int(n_spots),
            "cells_per_spot": s4.get("cells_per_spot_override", sim_info.get("cells_per_spot")),
            "assignment_rows": assignment_rows,
        },
        "engine": {
            "solver": s4.get("solver_method"),
            "sampling_sub_spots": s4.get("sampling_sub_spots"),
            "n_subspots": s4.get("n_subspots"),
            "seed": s4.get("seed"),
        },
        "acc_slot": {
            "acc_top1": None,
            "reason": "Not computed: Stage4 output does not carry query_id mapping.",
        },
        "cell_spot_eval": cell_spot_eval,
        "sha1": {
            "stage4_summary": sha1_file(summary_p),
            "cell_assignment": sha1_file(assign_p),
            "by_spot": sha1_file(by_spot_p),
            "sim_info": sha1_file(sim_info_p),
            "truth_spot_type_fraction": sha1_file(truth_spot_p),
            "truth_query_cell_spot": sha1_file(truth_query_p) if truth_query_p.exists() else None,
        },
    }

    out_p = out_dir / f"stage5_route2_s0__{args.run_tag}.json"
    out_p.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"[OK] wrote: {out_p}")
    return out


def compare_runs(baseline_json: Path, route2_json: Path, out_path: Path):
    def load(p: Path):
        return json.loads(p.read_text(encoding="utf-8"))

    b = load(baseline_json)
    r = load(route2_json)

    def diff(a, b):
        if a is None or b is None:
            return None
        return b - a

    comp = {
        "scenario": "S0",
        "baseline": baseline_json.name,
        "route2": route2_json.name,
        "delta": {
            "leakage_missing_in_by_spot": diff(b["leakage"]["missing_in_by_spot"], r["leakage"]["missing_in_by_spot"]),
            "L1_mean": diff(b["composition"]["L1_mean"], r["composition"]["L1_mean"]),
            "JS_mean": diff(b["composition"]["JS_mean"], r["composition"]["JS_mean"]),
            "corr_mean": diff(b["composition"]["corr_mean"], r["composition"]["corr_mean"]),
        },
    }
    out_path.write_text(json.dumps(comp, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"[OK] wrote compare: {out_path}")


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=False, default="real_brca_simS0_seed42")
    ap.add_argument("--run_tag", help="baseline / route2", required=False)
    ap.add_argument("--stage4_dir", help="stage4 output dir for this run", required=False)
    ap.add_argument("--sim_dir", help="sim truth dir (data/sim/<sample_base>/S0)", required=False)
    ap.add_argument("--out_dir", help="output dir for stage5 jsons", required=False)
    ap.add_argument("--seeds", nargs="+", type=int, help="run multiple seeds (will append seed to sample name)", required=False)
    ap.add_argument("--stage4_dir_template", help="template with {seed} for stage4_dir when batch", required=False)
    ap.add_argument("--out_dir_template", help="template with {seed} for out_dir when batch", required=False)
    ap.add_argument("--compare_baseline", help="baseline json path", required=False)
    ap.add_argument("--compare_route2", help="route2 json path", required=False)
    ap.add_argument("--compare_out", help="output compare json path", required=False)
    ap.add_argument(
        "--strict_config",
        action="store_true",
        help="若为true，则在Stage4与Stage5 n_filtered不一致时立即报错；否则仅在JSON中标记validation_checks并打印警告",
    )
    return ap.parse_args()


def main():
    args = parse_args()

    if args.seeds:
        # batch mode: run baseline/route2 assumed already produced per seed
        for seed in args.seeds:
            stage4_dir = Path(args.stage4_dir_template.format(seed=seed))
            sim_dir = Path(args.sim_dir.format(seed=seed)) if "{seed}" in (args.sim_dir or "") else Path(args.sim_dir)
            out_dir = Path(args.out_dir_template.format(seed=seed))
            for run_tag in ["baseline", "route2"]:
                run_args = argparse.Namespace(
                    sample=f"{args.sample}_seed{seed}" if "seed" not in args.sample else args.sample,
                    run_tag=run_tag,
                    stage4_dir=stage4_dir / run_tag if stage4_dir.is_dir() and (stage4_dir / run_tag).exists() else stage4_dir,
                    sim_dir=sim_dir,
                    out_dir=out_dir,
                    strict_config=args.strict_config,
                )
                run_eval(run_args)
        return

    if args.compare_baseline and args.compare_route2 and args.compare_out:
        compare_runs(Path(args.compare_baseline), Path(args.compare_route2), Path(args.compare_out))
        return

    if not all([args.run_tag, args.stage4_dir, args.sim_dir, args.out_dir]):
        raise SystemExit("Missing required arguments for evaluation run.")

    run_eval(args)


if __name__ == "__main__":
    main()

