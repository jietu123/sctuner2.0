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
import re

import numpy as np
import pandas as pd

try:
    import yaml
except ImportError:
    yaml = None
try:
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.svm import SVC
    from sklearn.model_selection import StratifiedKFold
    from sklearn.metrics import roc_auc_score, accuracy_score
except ImportError:
    make_pipeline = None
    StandardScaler = None
    SVC = None
    StratifiedKFold = None
    roc_auc_score = None
    accuracy_score = None

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


def resolve_optional_path(value: str | None, base_dir: Path) -> Path | None:
    if not value:
        return None
    p = Path(value)
    return p if p.is_absolute() else base_dir / p


def resolve_eval_paths(
    stage4_dir: Path,
    s4_summary: dict,
    args: argparse.Namespace,
) -> tuple[Path, Path, str]:
    assign_p = stage4_dir / "cell_assignment.csv"
    by_spot_p = stage4_dir / "cell_type_assignments_by_spot.csv"
    source = "default"

    assign_override = resolve_optional_path(args.assignment_csv, stage4_dir)
    by_spot_override = resolve_optional_path(args.by_spot_csv, stage4_dir)

    if assign_override is not None:
        assign_p = assign_override
        source = "override"
    elif args.use_post_rescue:
        source = "post_rescue"
        summary_path = s4_summary.get("post_rescue_assignment_path")
        if summary_path:
            assign_p = Path(summary_path)
            if not assign_p.is_absolute():
                assign_p = stage4_dir / assign_p
        else:
            assign_p = stage4_dir / "cell_assignment_post_rescue.csv"

    post_hint = assign_override is not None and "post_rescue" in assign_override.name
    if by_spot_override is not None:
        by_spot_p = by_spot_override
    elif args.use_post_rescue or post_hint:
        summary_path = s4_summary.get("post_rescue_by_spot_path")
        if summary_path:
            by_spot_p = Path(summary_path)
            if not by_spot_p.is_absolute():
                by_spot_p = stage4_dir / by_spot_p
        else:
            by_spot_p = stage4_dir / "cell_type_assignments_by_spot_post_rescue.csv"

    return assign_p, by_spot_p, source


def resolve_confidence_inputs(
    args: argparse.Namespace,
    stage4_dir: Path,
    sim_dir: Path,
    project_root: Path,
) -> dict:
    stage1_export = project_root / "data" / "processed" / args.sample / "stage1_preprocess" / "exported"

    def resolve_override(value: str | None) -> Path | None:
        if not value:
            return None
        p = Path(value)
        return p if p.is_absolute() else project_root / p

    def pick_path(candidates: list[Path | None]) -> tuple[Path | None, str | None]:
        for p in candidates:
            if p is not None and p.exists():
                return p, "override" if p in overrides else sources.get(p, None)
        return None, None

    overrides = set()
    sources = {}

    assigned_override = resolve_override(args.confidence_assigned_locations)
    if assigned_override is not None:
        overrides.add(assigned_override)
    assigned_candidate = assigned_override or (stage4_dir / "assigned_locations.csv")
    assigned_p = assigned_candidate if assigned_candidate.exists() else None
    assigned_source = "override" if assigned_override is not None else ("stage4_dir" if assigned_p else None)

    sc_expr_override = resolve_override(args.confidence_sc_expr)
    st_expr_override = resolve_override(args.confidence_st_expr)
    sc_meta_override = resolve_override(args.confidence_sc_meta)
    if sc_expr_override is not None:
        overrides.add(sc_expr_override)
    if st_expr_override is not None:
        overrides.add(st_expr_override)
    if sc_meta_override is not None:
        overrides.add(sc_meta_override)

    sim_sc_expr = sim_dir / "sc_expression.csv"
    sim_st_expr = sim_dir / "st_expression.csv"
    sim_sc_meta = sim_dir / "sc_metadata.csv"

    stage1_sc_expr = stage1_export / "sc_expression_normalized.csv"
    stage1_st_expr = stage1_export / "st_expression_normalized.csv"
    stage1_sc_meta = stage1_export / "sc_metadata.csv"

    sources.update({
        sim_sc_expr: "sim_dir",
        sim_st_expr: "sim_dir",
        sim_sc_meta: "sim_dir",
        stage1_sc_expr: "stage1_export",
        stage1_st_expr: "stage1_export",
        stage1_sc_meta: "stage1_export",
    })

    sc_expr_p, sc_expr_source = pick_path([sc_expr_override, sim_sc_expr, stage1_sc_expr])
    st_expr_p, st_expr_source = pick_path([st_expr_override, sim_st_expr, stage1_st_expr])
    sc_meta_p, sc_meta_source = pick_path([sc_meta_override, sim_sc_meta, stage1_sc_meta])

    ready = all([assigned_p, sc_expr_p, st_expr_p, sc_meta_p])

    return {
        "assigned_locations": str(assigned_p) if assigned_p else None,
        "sc_expression": str(sc_expr_p) if sc_expr_p else None,
        "st_expression": str(st_expr_p) if st_expr_p else None,
        "sc_metadata": str(sc_meta_p) if sc_meta_p else None,
        "sources": {
            "assigned_locations": assigned_source,
            "sc_expression": sc_expr_source,
            "st_expression": st_expr_source,
            "sc_metadata": sc_meta_source,
        },
        "ready": bool(ready),
    }


def clean_spot_id(value: str) -> str:
    return str(value).split()[0].split("\t")[0]


def load_expression_matrix(path: Path, clean_cols: bool = False) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    if "Gene" not in df.columns:
        raise ValueError(f"Missing Gene column in expression matrix: {path}")
    expr = df.drop(columns=["Gene"])
    expr = expr.apply(pd.to_numeric, errors="coerce").fillna(0.0).astype(np.float32)
    expr.index = df["Gene"].astype(str)
    if clean_cols:
        expr.columns = pd.Index(expr.columns).astype(str).map(clean_spot_id)
    else:
        expr.columns = expr.columns.astype(str)
    return expr


def load_truth_spot_presence(path: Path) -> pd.DataFrame:
    truth = pd.read_csv(path, sep=None, engine="python")
    if "spot_id" not in truth.columns:
        first_col = truth.columns[0]
        truth = truth.rename(columns={first_col: "spot_id"})
    truth["spot_id"] = truth["spot_id"].astype(str).map(clean_spot_id)
    truth = truth.set_index("spot_id")
    return truth


def load_truth_query(path: Path) -> pd.DataFrame:
    truth = pd.read_csv(path)
    required = {"query_id", "cell_id", "true_spot_id", "cell_type"}
    missing = required - set(truth.columns)
    if missing:
        raise ValueError(f"truth_query_cell_spot.csv missing columns: {sorted(missing)}")
    truth = truth.drop_duplicates(subset=["cell_id"], keep="first").copy()
    truth["cell_id"] = truth["cell_id"].astype(str)
    truth["true_spot_id"] = truth["true_spot_id"].astype(str).map(clean_spot_id)
    return truth


def compute_marker_scores(
    sc_expr: pd.DataFrame,
    pos_cells: list[str],
    neg_cells: list[str],
    strategy: str,
) -> pd.Series:
    if not pos_cells or not neg_cells:
        return pd.Series(dtype=float)
    if strategy == "t_stat":
        pos_log = np.log1p(sc_expr[pos_cells])
        neg_log = np.log1p(sc_expr[neg_cells])
        pos_mean = pos_log.mean(axis=1)
        neg_mean = neg_log.mean(axis=1)
        pos_var = pos_log.var(axis=1, ddof=1)
        neg_var = neg_log.var(axis=1, ddof=1)
        denom = np.sqrt(pos_var / max(len(pos_cells), 1) + neg_var / max(len(neg_cells), 1))
        denom = denom.replace(0.0, np.nan)
        score = (pos_mean - neg_mean) / denom
        return score.fillna(0.0)
    pos_mean = sc_expr[pos_cells].mean(axis=1)
    neg_mean = sc_expr[neg_cells].mean(axis=1)
    return np.log1p(pos_mean) - np.log1p(neg_mean)


def select_marker_genes(
    sc_expr: pd.DataFrame,
    pos_cells: list[str],
    neg_cells: list[str],
    st_genes: set[str],
    n_markers: int,
    min_markers: int,
    strategy: str,
) -> tuple[list[str], pd.Series]:
    score = compute_marker_scores(sc_expr, pos_cells, neg_cells, strategy)
    if score.empty:
        return [], score
    score = score.loc[score.index.intersection(st_genes)]
    score = score.sort_values(ascending=False, kind="mergesort")
    markers = score.head(n_markers).index.tolist()
    return (markers if len(markers) >= min_markers else []), score


def build_pseudobulk(
    sc_expr: pd.DataFrame,
    markers: list[str],
    spot_cells: dict[str, list[str]],
    spots: list[str],
) -> np.ndarray:
    rows = []
    for spot in spots:
        cell_ids = spot_cells.get(spot, [])
        if not cell_ids:
            continue
        mat = sc_expr.loc[markers, cell_ids]
        rows.append(mat.mean(axis=1).to_numpy(dtype=float))
    return np.vstack(rows) if rows else np.empty((0, len(markers)))


def compute_spot_support(
    st_expr: pd.DataFrame,
    markers: list[str],
    weights: pd.Series,
) -> pd.Series:
    mat = np.log1p(st_expr.loc[markers].to_numpy(dtype=float)).T
    score = mat.dot(weights.to_numpy(dtype=float))
    if np.allclose(score.max(), score.min()):
        return pd.Series(0.5, index=st_expr.columns)
    score = (score - score.min()) / (score.max() - score.min())
    return pd.Series(score, index=st_expr.columns)


def auc_score(y_true: np.ndarray, y_score: np.ndarray) -> float | None:
    y_true = np.asarray(y_true)
    y_score = np.asarray(y_score)
    pos = y_true == 1
    neg = y_true == 0
    n_pos = int(pos.sum())
    n_neg = int(neg.sum())
    if n_pos == 0 or n_neg == 0:
        return None
    ranks = pd.Series(y_score).rank(method="average").to_numpy()
    sum_pos = float(ranks[pos].sum())
    return (sum_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def stratified_folds(labels: np.ndarray, n_splits: int, rng: np.random.Generator):
    pos_idx = np.where(labels == 1)[0]
    neg_idx = np.where(labels == 0)[0]
    rng.shuffle(pos_idx)
    rng.shuffle(neg_idx)
    folds = [[] for _ in range(n_splits)]
    for i, idx in enumerate(pos_idx):
        folds[i % n_splits].append(idx)
    for i, idx in enumerate(neg_idx):
        folds[i % n_splits].append(idx)
    for fold in folds:
        test_idx = np.array(fold, dtype=int)
        mask = np.zeros(len(labels), dtype=bool)
        mask[test_idx] = True
        train_idx = np.where(~mask)[0]
        yield train_idx, test_idx


def compute_cv_metrics(
    train: np.ndarray,
    labels: np.ndarray,
    method: str,
    rng: np.random.Generator,
    n_folds: int,
    svm_c: float | None = None,
) -> dict | None:
    labels = np.asarray(labels)
    n_pos = int((labels == 1).sum())
    n_neg = int((labels == 0).sum())
    max_folds = min(n_folds, n_pos, n_neg)
    if max_folds < 2:
        return None
    auc_list = []
    acc_list = []
    if StratifiedKFold is not None:
        splitter = StratifiedKFold(n_splits=max_folds, shuffle=True, random_state=int(rng.integers(0, 1_000_000)))
        splits = splitter.split(train, labels)
    else:
        splits = stratified_folds(labels, max_folds, rng)
    for train_idx, test_idx in splits:
        x_train = train[train_idx]
        y_train = labels[train_idx]
        x_test = train[test_idx]
        y_test = labels[test_idx]
        if method == "svm" and SVC is not None:
            c_value = svm_c if svm_c is not None else 1.0
            model = make_pipeline(StandardScaler(), SVC(kernel="linear", C=c_value, probability=True, random_state=1))
            model.fit(x_train, y_train)
            scores = model.predict_proba(x_test)[:, 1]
        else:
            pos_mean = np.log1p(x_train[y_train == 1].mean(axis=0))
            neg_mean = np.log1p(x_train[y_train == 0].mean(axis=0))
            weights = pos_mean - neg_mean
            scores = np.log1p(x_test).dot(weights)
            if np.allclose(scores.max(), scores.min()):
                scores = np.full_like(scores, 0.5, dtype=float)
            else:
                scores = (scores - scores.min()) / (scores.max() - scores.min())
        if roc_auc_score is not None:
            auc = float(roc_auc_score(y_test, scores))
        else:
            auc = auc_score(y_test, scores)
        acc = float(accuracy_score(y_test, scores >= 0.5)) if accuracy_score is not None else float((scores >= 0.5).mean())
        if auc is not None:
            auc_list.append(auc)
        acc_list.append(acc)
    return {
        "folds": int(max_folds),
        "auc_mean": float(np.mean(auc_list)) if auc_list else None,
        "acc_mean": float(np.mean(acc_list)) if acc_list else None,
    }


def parse_float_list(value: str | None, default: list[float]) -> list[float]:
    if not value:
        return default
    parts = [p for p in re.split(r"[,\s]+", value.strip()) if p]
    values = []
    for part in parts:
        try:
            values.append(float(part))
        except ValueError:
            continue
    values = sorted(set(values))
    return values if values else default


def choose_svm_c(
    train: np.ndarray,
    labels: np.ndarray,
    c_grid: list[float],
    rng: np.random.Generator,
    n_folds: int,
) -> tuple[float, dict | None]:
    best_c = c_grid[0]
    best_cv = None
    best_score = -1.0
    for c_value in c_grid:
        cv = compute_cv_metrics(train, labels, "svm", rng, n_folds, svm_c=c_value)
        if cv and cv.get("auc_mean") is not None and cv["auc_mean"] > best_score:
            best_score = cv["auc_mean"]
            best_c = c_value
            best_cv = cv
    return best_c, best_cv


def parse_thresholds(value: str | None, default: list[float]) -> list[float]:
    if not value:
        return default
    parts = [p for p in re.split(r"[,\s]+", value.strip()) if p]
    thresholds = []
    for part in parts:
        try:
            thresholds.append(float(part))
        except ValueError:
            continue
    thresholds = sorted(set(thresholds))
    return thresholds if thresholds else default


def build_segment_stats(values: pd.Series, low: float, high: float) -> dict:
    vals = values.dropna()
    total = len(vals)
    if total == 0:
        return {"low": None, "mid": None, "high": None}
    low_mask = vals < low
    high_mask = vals >= high
    mid_mask = (~low_mask) & (~high_mask)
    return {
        "low": float(low_mask.sum() / total),
        "mid": float(mid_mask.sum() / total),
        "high": float(high_mask.sum() / total),
    }


def compute_type_support_rate(assigned_conf: pd.DataFrame, threshold: float) -> dict:
    rates = {}
    for cell_type, group in assigned_conf.groupby("CellType"):
        vals = group["confidence"].dropna()
        if len(vals) == 0:
            continue
        rates[cell_type] = float((vals >= threshold).sum() / len(vals))
    return rates


def compute_confidence_metrics(
    confidence_inputs: dict,
    truth_spot_path: Path,
    truth_query_path: Path | None,
    out_dir: Path,
    args: argparse.Namespace,
) -> dict:
    if args.disable_confidence:
        return {"enabled": False, "ready": False, "reason": "disabled_by_flag"}
    if not confidence_inputs.get("ready"):
        return {"enabled": True, "ready": False, "reason": "missing_inputs"}
    if not truth_spot_path.exists():
        return {"enabled": True, "ready": False, "reason": "missing_truth_spot"}

    assigned_p = Path(confidence_inputs["assigned_locations"])
    sc_expr_p = Path(confidence_inputs["sc_expression"])
    st_expr_p = Path(confidence_inputs["st_expression"])
    sc_meta_p = Path(confidence_inputs["sc_metadata"])

    assigned = pd.read_csv(assigned_p)
    required = {"OriginalCID", "CellType", "SpotID"}
    if not required.issubset(set(assigned.columns)):
        raise ValueError(f"assigned_locations.csv missing columns: {sorted(required - set(assigned.columns))}")
    assigned["SpotID"] = assigned["SpotID"].astype(str).map(clean_spot_id)
    assigned["OriginalCID"] = assigned["OriginalCID"].astype(str)
    assigned["CellType"] = assigned["CellType"].astype(str)

    sc_meta = pd.read_csv(sc_meta_p)
    if "cell_id" not in sc_meta.columns or "cell_type" not in sc_meta.columns:
        raise ValueError("sc_metadata.csv missing cell_id or cell_type")
    sc_meta["cell_id"] = sc_meta["cell_id"].astype(str)
    sc_meta["cell_type"] = sc_meta["cell_type"].astype(str)

    sc_expr = load_expression_matrix(sc_expr_p, clean_cols=False)
    st_expr = load_expression_matrix(st_expr_p, clean_cols=True)
    truth = load_truth_spot_presence(truth_spot_path)

    sc_cells = [c for c in sc_meta["cell_id"].unique().tolist() if c in sc_expr.columns]
    if not sc_cells:
        return {"enabled": True, "ready": False, "reason": "no_sc_cells_matched"}
    sc_expr = sc_expr[sc_cells]
    sc_meta = sc_meta[sc_meta["cell_id"].isin(sc_cells)]

    assigned = assigned[assigned["OriginalCID"].isin(sc_cells)]
    if assigned.empty:
        return {"enabled": True, "ready": False, "reason": "no_assigned_cells_matched"}

    truth_types = [c for c in truth.columns if c in sc_meta["cell_type"].unique().tolist()]
    assigned_types = assigned["CellType"].unique().tolist()
    types = sorted(set(truth_types).intersection(assigned_types))
    if not types:
        return {"enabled": True, "ready": False, "reason": "no_common_cell_types"}

    spot_cells = (
        assigned.groupby("SpotID")["OriginalCID"]
        .apply(list)
        .to_dict()
    )
    st_genes = set(st_expr.index)
    rng = np.random.default_rng(args.confidence_seed)
    c_grid = parse_float_list(args.confidence_svm_c_grid, [1.0])

    support = pd.DataFrame(index=st_expr.columns)
    marker_summary = []
    marker_rows = []
    used_types = []
    method_used = None
    cv_auc_by_type = {}
    cv_acc_by_type = {}
    svm_c_by_type = {}

    for cell_type in types:
        pos_cells = sc_meta.loc[sc_meta["cell_type"] == cell_type, "cell_id"].tolist()
        neg_cells = sc_meta.loc[sc_meta["cell_type"] != cell_type, "cell_id"].tolist()
        markers, score = select_marker_genes(
            sc_expr,
            pos_cells,
            neg_cells,
            st_genes,
            args.confidence_n_markers,
            args.confidence_min_markers,
            args.confidence_marker_strategy,
        )
        if not markers:
            marker_summary.append({"cell_type": cell_type, "markers": 0, "status": "skip"})
            continue
        for rank, gene in enumerate(markers, start=1):
            marker_rows.append({
                "cell_type": cell_type,
                "gene": gene,
                "rank": rank,
                "score": float(score.get(gene, 0.0)),
            })

        pos_spots = truth.index[truth[cell_type] > 0].tolist()
        neg_spots = truth.index[truth[cell_type] <= 0].tolist()
        pos_spots = [s for s in pos_spots if s in spot_cells]
        neg_spots = [s for s in neg_spots if s in spot_cells]
        if not pos_spots or not neg_spots:
            marker_summary.append({"cell_type": cell_type, "markers": len(markers), "status": "skip_spots"})
            continue

        pos_sample = rng.choice(pos_spots, size=args.confidence_groupsize, replace=len(pos_spots) < args.confidence_groupsize).tolist()
        neg_sample = rng.choice(neg_spots, size=args.confidence_groupsize, replace=len(neg_spots) < args.confidence_groupsize).tolist()

        train_pos = build_pseudobulk(sc_expr, markers, spot_cells, pos_sample)
        train_neg = build_pseudobulk(sc_expr, markers, spot_cells, neg_sample)
        if train_pos.shape[0] == 0 or train_neg.shape[0] == 0:
            marker_summary.append({"cell_type": cell_type, "markers": len(markers), "status": "skip_train"})
            continue

        train = np.vstack([train_pos, train_neg])
        labels = np.array([1] * train_pos.shape[0] + [0] * train_neg.shape[0])
        used_types.append(cell_type)

        method = args.confidence_method
        if method == "auto":
            method = "svm" if SVC is not None else "weight_diff"
        if method == "svm" and SVC is None:
            method = "weight_diff"

        if method == "svm":
            if len(c_grid) > 1:
                best_c, best_cv = choose_svm_c(train, labels, c_grid, rng, args.confidence_cv_folds)
            else:
                best_c, best_cv = c_grid[0], None
            svm_c_by_type[cell_type] = float(best_c)
            model = make_pipeline(StandardScaler(), SVC(kernel="linear", C=best_c, probability=True, random_state=args.confidence_seed))
            model.fit(train, labels)
            test = np.log1p(st_expr.loc[markers].to_numpy(dtype=float)).T
            pred = model.predict_proba(test)[:, 1]
            support[cell_type] = pred
        else:
            pos_mean = np.log1p(sc_expr.loc[markers, pos_cells].mean(axis=1))
            neg_mean = np.log1p(sc_expr.loc[markers, neg_cells].mean(axis=1))
            weights = pos_mean - neg_mean
            support[cell_type] = compute_spot_support(st_expr, markers, weights)

        method_used = method
        cv = best_cv if (method == "svm" and best_cv is not None) else compute_cv_metrics(
            train,
            labels,
            method,
            rng,
            args.confidence_cv_folds,
            svm_c=svm_c_by_type.get(cell_type),
        )
        cv_auc = cv.get("auc_mean") if cv else None
        cv_acc = cv.get("acc_mean") if cv else None
        if cv_auc is not None:
            cv_auc_by_type[cell_type] = cv_auc
        if cv_acc is not None:
            cv_acc_by_type[cell_type] = cv_acc
        marker_summary.append({
            "cell_type": cell_type,
            "markers": len(markers),
            "status": "ok",
            "method": method,
            "cv_auc_mean": cv_auc,
            "cv_acc_mean": cv_acc,
            "cv_folds": cv.get("folds") if cv else None,
            "svm_c": svm_c_by_type.get(cell_type),
        })

    if support.empty:
        return {"enabled": True, "ready": False, "reason": "no_support_matrix"}

    marker_path = None
    if marker_rows:
        marker_path = out_dir / "confidence_markers.csv"
        pd.DataFrame(marker_rows).to_csv(marker_path, index=False, encoding="utf-8")

    support_path = out_dir / "spot_type_support.csv"
    support.to_csv(support_path, index_label="spot_id")

    support_long = support.stack().reset_index()
    support_long.columns = ["SpotID", "CellType", "confidence"]
    assigned_conf = assigned.merge(support_long, on=["SpotID", "CellType"], how="left")

    conf_path = out_dir / "cell_confidence.csv"
    assigned_conf.to_csv(conf_path, index=False, encoding="utf-8")

    seg_low = float(args.confidence_segment_low)
    seg_high = float(args.confidence_segment_high)
    if seg_low >= seg_high:
        seg_low, seg_high = min(seg_low, seg_high), max(seg_low, seg_high)

    validation = None
    validation_path = None
    if truth_query_path is not None and truth_query_path.exists():
        try:
            truth_query = load_truth_query(truth_query_path)
            truth_map = truth_query.set_index("cell_id")["true_spot_id"].to_dict()
            assigned_conf["true_spot_id"] = assigned_conf["OriginalCID"].map(truth_map)
            assigned_conf["has_truth"] = assigned_conf["true_spot_id"].notna()
            assigned_conf["is_correct"] = (
                assigned_conf["has_truth"] & (assigned_conf["true_spot_id"] == assigned_conf["SpotID"])
            )
            validation_path = out_dir / "confidence_validation.csv"
            assigned_conf.to_csv(validation_path, index=False, encoding="utf-8")

            truth_mask = assigned_conf["has_truth"]
            correct_vals = assigned_conf.loc[truth_mask & assigned_conf["is_correct"], "confidence"]
            incorrect_vals = assigned_conf.loc[truth_mask & (~assigned_conf["is_correct"]), "confidence"]
            validation = {
                "n_total_with_truth": int(truth_mask.sum()),
                "n_correct": int((truth_mask & assigned_conf["is_correct"]).sum()),
                "n_incorrect": int((truth_mask & (~assigned_conf["is_correct"])).sum()),
                "mean_confidence_correct": float(correct_vals.mean()) if len(correct_vals) else None,
                "mean_confidence_incorrect": float(incorrect_vals.mean()) if len(incorrect_vals) else None,
                "median_confidence_correct": float(correct_vals.median()) if len(correct_vals) else None,
                "median_confidence_incorrect": float(incorrect_vals.median()) if len(incorrect_vals) else None,
                "high_confidence_fraction_correct": float((correct_vals >= args.confidence_threshold).sum() / len(correct_vals)) if len(correct_vals) else None,
                "high_confidence_fraction_incorrect": float((incorrect_vals >= args.confidence_threshold).sum() / len(incorrect_vals)) if len(incorrect_vals) else None,
                "segment_correct": build_segment_stats(correct_vals, seg_low, seg_high),
                "segment_incorrect": build_segment_stats(incorrect_vals, seg_low, seg_high),
            }
        except Exception as e:
            print(f"[WARN] confidence validation failed: {e}")

    conf_values = assigned_conf["confidence"]
    n_total = int(len(assigned_conf))
    n_with = int(conf_values.notna().sum())
    if n_with > 0:
        mean_conf = float(conf_values.mean(skipna=True))
        high_conf = float((conf_values >= args.confidence_threshold).sum() / n_with)
    else:
        mean_conf = None
        high_conf = None

    type_support = compute_type_support_rate(assigned_conf, args.confidence_threshold)
    type_support_mean = float(np.mean(list(type_support.values()))) if type_support else None
    cv_auc_mean = float(np.mean(list(cv_auc_by_type.values()))) if cv_auc_by_type else None
    cv_acc_mean = float(np.mean(list(cv_acc_by_type.values()))) if cv_acc_by_type else None

    thresholds = parse_thresholds(args.confidence_thresholds, [0.1, 0.3, 0.5, 0.7, 0.9])
    threshold_scan = []
    for th in thresholds:
        th_rates = compute_type_support_rate(assigned_conf, th)
        th_mean = float(np.mean(list(th_rates.values()))) if th_rates else None
        th_high = float((conf_values.dropna() >= th).sum() / n_with) if n_with else None
        threshold_scan.append({
            "threshold": float(th),
            "high_confidence_fraction": th_high,
            "type_support_rate_mean": th_mean,
            "type_support_rate": th_rates,
        })

    segment_global = build_segment_stats(conf_values, seg_low, seg_high)
    segment_by_type = {}
    for cell_type, group in assigned_conf.groupby("CellType"):
        segment_by_type[cell_type] = build_segment_stats(group["confidence"], seg_low, seg_high)

    report = {
        "thresholds": thresholds,
        "threshold_scan": threshold_scan,
        "segments": {
            "low": seg_low,
            "high": seg_high,
            "global": segment_global,
            "by_type": segment_by_type,
        },
        "cv": {
            "folds": int(args.confidence_cv_folds),
            "auc_mean": cv_auc_mean,
            "acc_mean": cv_acc_mean,
        },
    }
    report_path = out_dir / "confidence_report.json"
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    md_path = out_dir / "confidence_report.md"
    md_lines = [
        "# Confidence Report",
        "",
        "## Threshold Scan",
        "",
        "| threshold | high_confidence_fraction | type_support_rate_mean |",
        "| --- | --- | --- |",
    ]
    for row in threshold_scan:
        md_lines.append(
            f"| {row['threshold']:.2f} | {row['high_confidence_fraction'] if row['high_confidence_fraction'] is not None else 'NA'} | {row['type_support_rate_mean'] if row['type_support_rate_mean'] is not None else 'NA'} |"
        )
    md_lines.extend([
        "",
        "## Confidence Segments",
        "",
        "| segment | fraction |",
        "| --- | --- |",
        f"| low (<{seg_low}) | {segment_global.get('low')} |",
        f"| mid ({seg_low}-{seg_high}) | {segment_global.get('mid')} |",
        f"| high (>= {seg_high}) | {segment_global.get('high')} |",
        "",
    ])
    md_path.write_text("\n".join(md_lines), encoding="utf-8")

    return {
        "enabled": True,
        "ready": True,
        "method": method_used,
        "marker_strategy": args.confidence_marker_strategy,
        "svm_c_grid": c_grid,
        "svm_c_by_type": svm_c_by_type,
        "marker_path": str(marker_path.relative_to(out_dir)) if marker_path else None,
        "threshold": float(args.confidence_threshold),
        "thresholds": thresholds,
        "threshold_scan": threshold_scan,
        "segments": {
            "low": seg_low,
            "high": seg_high,
            "global": segment_global,
            "by_type": segment_by_type,
        },
        "groupsize": int(args.confidence_groupsize),
        "n_markers": int(args.confidence_n_markers),
        "min_markers": int(args.confidence_min_markers),
        "seed": int(args.confidence_seed),
        "cv_folds": int(args.confidence_cv_folds),
        "cv_auc_mean": cv_auc_mean,
        "cv_acc_mean": cv_acc_mean,
        "cv_auc_by_type": cv_auc_by_type,
        "cv_acc_by_type": cv_acc_by_type,
        "n_cell_types": int(len(used_types)),
        "n_cells_total": n_total,
        "n_cells_with_confidence": n_with,
        "mean_confidence": mean_conf,
        "high_confidence_fraction": high_conf,
        "type_support_rate": type_support,
        "type_support_rate_mean": type_support_mean,
        "spot_type_support_path": str(support_path.relative_to(out_dir)),
        "cell_confidence_path": str(conf_path.relative_to(out_dir)),
        "confidence_validation_path": str(validation_path.relative_to(out_dir)) if validation_path else None,
        "confidence_validation": validation,
        "confidence_report_path": str(report_path.relative_to(out_dir)),
        "confidence_report_md_path": str(md_path.relative_to(out_dir)),
        "marker_summary": marker_summary,
    }


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
            print(f"[self_check_1] WARNING: extra_pred中包含missing_types: {missing_types_in_extra}")
        else:
            print(f"[self_check_1] OK: extra_pred中不包含missing_types (CD8等)")
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
    print("[self_check_2] OK: Acc@1只基于truth_query cells (truth_present=True)，不包含extra_pred")
    
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
    truth_filter_enabled: bool | None = None,
    truth_filter_removed: int | None = None,
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

    # 强制检查4：Stage4 n_filtered vs Stage5重建数量（含 truth-filter）
    truth_removed = int(truth_filter_removed or 0)
    expected_stage4_filtered = actual_filtered + truth_removed
    
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
        "truth_filter": {
            "enabled": bool(truth_filter_enabled),
            "removed": int(truth_removed),
        },
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
                "value": bool(expected_stage4_filtered == n_filtered_total),
                "actual_filtered": int(actual_filtered),
                "truth_filter_removed": int(truth_removed),
                "expected_stage4_n_filtered": int(expected_stage4_filtered),
                "stage4_n_filtered": int(n_filtered_total),
                "description": "Ensure Stage4 n_filtered equals Stage5 filtered + truth_filter_removed"
            },
        }
    }


def run_eval(args: argparse.Namespace):
    stage4_dir = Path(args.stage4_dir)
    sim_dir = Path(args.sim_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    project_root = detect_root()
    summary_p = stage4_dir / "stage4_summary.json"

    sim_info_p = sim_dir / "sim_info.json"
    truth_spot_p = sim_dir / "sim_truth_spot_type_fraction.csv"
    truth_query_p = sim_dir / "sim_truth_query_cell_spot.csv"

    with summary_p.open("r", encoding="utf-8") as f:
        s4 = json.load(f)

    assign_p, by_spot_p, assignment_source = resolve_eval_paths(stage4_dir, s4, args)
    if not assign_p.exists():
        raise FileNotFoundError(f"assignment_csv not found: {assign_p}")
    if not by_spot_p.exists():
        raise FileNotFoundError(f"by_spot_csv not found: {by_spot_p}")

    with sim_info_p.open("r", encoding="utf-8") as f:
        sim_info = json.load(f)

    confidence_inputs = resolve_confidence_inputs(args, stage4_dir, sim_dir, project_root)
    try:
        confidence_metrics = compute_confidence_metrics(
            confidence_inputs,
            truth_spot_p,
            truth_query_p if truth_query_p.exists() else None,
            out_dir,
            args,
        )
    except Exception as e:
        print(f"[WARN] confidence evaluation failed: {e}")
        confidence_metrics = {"enabled": True, "ready": False, "reason": f"error:{e}"}
    missing_type = sim_info.get("missing_type")

    ca = pd.read_csv(assign_p)
    leak_assign = int((ca["cell_type"] == missing_type).sum()) if "cell_type" in ca.columns else None

    by_spot = load_by_spot_counts(by_spot_p)
    leak_by_spot = int(by_spot[missing_type].sum()) if missing_type in by_spot.columns else 0

    assignment_rows = int(len(ca))
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
            truth_filter_enabled=s4.get("truth_filter_enabled"),
            truth_filter_removed=s4.get("truth_filter_removed"),
        )
        # 严格模式下，若Stage4与Stage5重建的filtered数量不一致则直接报错
        if args.strict_config and filter_audit is not None:
            vc = (filter_audit.get("validation_checks") or {}).get("check4_stage4_stage5_n_filtered_match")
            if isinstance(vc, dict) and not vc.get("value", True):
                msg = (
                    "[strict_config] Stage4 vs Stage5 n_filtered mismatch: "
                    f"Stage4={vc.get('stage4_n_filtered')} vs "
                    f"Stage5_filtered={vc.get('actual_filtered')} + truth_removed={vc.get('truth_filter_removed')}"
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
        "assignment": {
            "source": assignment_source,
            "cell_assignment_path": str(assign_p),
            "by_spot_path": str(by_spot_p),
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
        "confidence_inputs": confidence_inputs,
        "confidence_metrics": confidence_metrics,
        "sha1": {
            "stage4_summary": sha1_file(summary_p),
            "cell_assignment": sha1_file(assign_p) if assign_p.exists() else None,
            "by_spot": sha1_file(by_spot_p) if by_spot_p.exists() else None,
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
    ap.add_argument("--assignment_csv", help="override assignment csv (absolute or relative to stage4_dir)", required=False)
    ap.add_argument("--by_spot_csv", help="override by-spot csv (absolute or relative to stage4_dir)", required=False)
    ap.add_argument(
        "--use_post_rescue",
        action="store_true",
        help="use post-rescue outputs if available (cell_assignment_post_rescue.csv)",
    )
    ap.add_argument("--seeds", nargs="+", type=int, help="run multiple seeds (will append seed to sample name)", required=False)
    ap.add_argument("--stage4_dir_template", help="template with {seed} for stage4_dir when batch", required=False)
    ap.add_argument("--out_dir_template", help="template with {seed} for out_dir when batch", required=False)
    ap.add_argument("--compare_baseline", help="baseline json path", required=False)
    ap.add_argument("--compare_route2", help="route2 json path", required=False)
    ap.add_argument("--compare_out", help="output compare json path", required=False)
    ap.add_argument("--confidence_sc_expr", help="override sc expression path for confidence inputs", required=False)
    ap.add_argument("--confidence_st_expr", help="override st expression path for confidence inputs", required=False)
    ap.add_argument("--confidence_sc_meta", help="override sc metadata path for confidence inputs", required=False)
    ap.add_argument("--confidence_assigned_locations", help="override assigned_locations.csv path for confidence inputs", required=False)
    ap.add_argument("--disable_confidence", action="store_true", help="skip confidence evaluation")
    ap.add_argument(
        "--confidence_threshold",
        type=float,
        default=0.35,
        help="threshold for high-confidence fraction (recommended range: 0.35-0.40)",
    )
    ap.add_argument("--confidence_n_markers", type=int, default=50, help="number of marker genes per cell type")
    ap.add_argument("--confidence_min_markers", type=int, default=5, help="minimum markers required to score a type")
    ap.add_argument("--confidence_groupsize", type=int, default=50, help="pseudobulk group size per class")
    ap.add_argument("--confidence_seed", type=int, default=1, help="seed for confidence evaluation")
    ap.add_argument(
        "--confidence_marker_strategy",
        choices=["mean_logfc", "t_stat"],
        default="mean_logfc",
        help="marker score strategy for confidence evaluation",
    )
    ap.add_argument(
        "--confidence_cv_folds",
        type=int,
        default=5,
        help="number of CV folds for confidence model stability",
    )
    ap.add_argument(
        "--confidence_thresholds",
        default="0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90",
        help="threshold scan list, e.g. '0.35,0.4,0.5'",
    )
    ap.add_argument(
        "--confidence_segment_low",
        type=float,
        default=0.3,
        help="low confidence upper bound for segment stats",
    )
    ap.add_argument(
        "--confidence_segment_high",
        type=float,
        default=0.7,
        help="high confidence lower bound for segment stats",
    )
    ap.add_argument(
        "--confidence_method",
        choices=["auto", "svm", "weight_diff"],
        default="auto",
        help="confidence model: svm if available, otherwise weight_diff",
    )
    ap.add_argument(
        "--confidence_svm_c_grid",
        default="1",
        help="comma-separated C values for SVM grid search (e.g. '0.1,1,10')",
    )
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

