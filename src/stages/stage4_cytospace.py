"""
Stage4 (Route2, full): wrap external CytoSPACE mapping with Stage3 filtering.

Pipeline:
- Load Stage1 exports (sc_expression_normalized.csv, st_expression_normalized.csv, st_coordinates.csv)
- Load Stage3 outputs (cell_type_relabel.csv, type_support.csv)
- Filter cells: either drop plugin_type == Unknown_sc_only (default) or drop support_category == unsupported
- Prepare CytoSPACE input files (gene x cell, gene x spot, cell type table, type fractions)
- Run CytoSPACE main_cytospace
- Convert outputs to Stage4 format (cell_assignment.csv, stage4_summary.json)

Outputs are written to result/<sample>/stage4_cytospace/
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Literal

import pandas as pd
import hashlib
import sys as _sys

# 确保项目根目录在 sys.path 中，便于在多种运行方式下都能导入 src.*
_HERE = Path(__file__).resolve()
_ROOT = _HERE.parents[2]
if str(_ROOT) not in _sys.path:
    _sys.path.insert(0, str(_ROOT))

from src.utils.type_name import load_alias_map, canonicalize_type_name, normalize_type_name


FilterMode = Literal["plugin_unknown", "unsupported", "none"]
FilterScope = Literal["unsupported_all", "missing_only"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Route2 Stage4 full (CytoSPACE wrapper)")
    p.add_argument("--sample", default="real_brca_simS0_seed42", help="sample id")
    p.add_argument("--project_root", default=None, help="override project root")
    p.add_argument("--filter_mode", choices=["plugin_unknown", "unsupported", "none"], default="none",
                   help="Filter mode: 'none' for baseline (official CytoSPACE), 'plugin_unknown' or 'unsupported' for Route2")
    p.add_argument("--cell_type_column", choices=["plugin_type", "orig_type", "sc_meta"], default="sc_meta",
                   help="which column to use as cell_type: 'sc_meta' for baseline (original types), 'plugin_type' for Route2")
    p.add_argument("--missing_type", default="T cells CD8", help="target missing type for reporting")
    p.add_argument("--mean_cell_numbers", type=int, default=5, help="CytoSPACE mean_cell_numbers")
    p.add_argument("--solver_method", default="lap_CSPR", help="CytoSPACE solver_method (lap_CSPR avoids lapjv dependency)")
    p.add_argument("--distance_metric", default="Pearson_correlation", help="CytoSPACE distance_metric")
    p.add_argument("--seed", type=int, default=42, help="random seed")
    p.add_argument("--sampling_sub_spots", action="store_true", default=True, help="enable sub-spot sampling to reduce memory")
    p.add_argument("--n_subspots", type=int, default=800, help="number of selected sub-spots when sampling_sub_spots is on")
    p.add_argument(
        "--n_processors",
        type=int,
        default=1,
        help="number of processors (reduce to 1 to avoid memory issues)",
    )
    p.add_argument(
        "--filter_scope",
        choices=["unsupported_all", "missing_only"],
        default="unsupported_all",
        help=(
            "scope of Route2 filtering when filter_mode uses Stage3 plugin/unsupported info; "
            "'unsupported_all' (default) drops all unsupported/Unknown_sc_only cells; "
            "'missing_only' only drops cells whose orig_type matches the specified missing_type"
        ),
    )
    return p.parse_args()


def detect_root() -> Path:
    here = Path(__file__).resolve()
    try:
        return here.parents[2]
    except IndexError:
        return here.parent.parent


def _clean_spot_index(idx: pd.Index) -> pd.Index:
    return pd.Index([str(x).split("\t")[0] for x in idx])


def load_stage1(root: Path, sample: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    base = root / "data" / "processed" / sample / "stage1_preprocess" / "exported"
    sc_expr = pd.read_csv(base / "sc_expression_normalized.csv", index_col=0, sep=None, engine="python")
    st_expr = pd.read_csv(base / "st_expression_normalized.csv", index_col=0, sep=None, engine="python")
    st_coords = pd.read_csv(base / "st_coordinates.csv", index_col=0, sep=None, engine="python")
    sc_meta = pd.read_csv(base / "sc_metadata.csv", sep=None, engine="python")
    if "cell_id" in sc_expr.columns:
        sc_expr = sc_expr.drop(columns=["cell_id"])
    if "cell_id" in st_expr.columns:
        st_expr = st_expr.drop(columns=["cell_id"])
    sc_expr = sc_expr.apply(pd.to_numeric, errors="coerce").dropna(axis=1, how="all").astype("float32")
    st_expr = st_expr.apply(pd.to_numeric, errors="coerce").dropna(axis=1, how="all").astype("float32")
    st_expr.index = _clean_spot_index(st_expr.index)
    st_coords.index = _clean_spot_index(st_coords.index)
    return sc_expr, st_expr, st_coords, sc_meta


def load_stage3(root: Path, sample: str) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    base = root / "data" / "processed" / sample / "stage3_typematch"
    relabel = pd.read_csv(base / "cell_type_relabel.csv")
    ts_path = base / "type_support.csv"
    ts = pd.read_csv(ts_path) if ts_path.exists() else None
    return relabel, ts


def load_stage3_thresholds(root: Path, sample: str) -> tuple[float | None, float | None]:
    cfg_path = root / "configs" / "datasets" / f"{sample}.yaml"
    if not cfg_path.exists():
        return None, None
    try:
        import yaml  # type: ignore
    except Exception:
        return None, None
    with cfg_path.open("r", encoding="utf-8") as f:
        data = yaml.safe_load(f)
    stage3 = (data or {}).get("stage3") or {}
    return stage3.get("weak_th"), stage3.get("strong_th")


def load_missing_truth(root: Path, sample: str) -> str | None:
    sim_info_path = root / "data" / "processed" / sample / "stage1_preprocess" / "exported" / "sim_info.json"
    if not sim_info_path.exists():
        return None
    try:
        data = json.loads(sim_info_path.read_text(encoding="utf-8"))
        return data.get("missing_type")
    except Exception:
        return None


def load_cells_per_spot_override(root: Path, sample: str) -> int | None:
    sim_info_path = root / "data" / "processed" / sample / "stage1_preprocess" / "exported" / "sim_info.json"
    if not sim_info_path.exists():
        return None
    try:
        data = json.loads(sim_info_path.read_text(encoding="utf-8"))
        return data.get("cells_per_spot")
    except Exception:
        return None


def sha1(path: Path) -> str:
    h = hashlib.sha1()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def apply_filter(
    relabel: pd.DataFrame,
    type_support: pd.DataFrame | None,
    mode: FilterMode,
    missing_type: str | None = None,
    filter_scope: FilterScope = "unsupported_all",
) -> tuple[pd.DataFrame, int, dict]:
    """
    根据Stage3输出决定在Stage4中要丢弃哪些细胞。

    - mode == "plugin_unknown": 使用 plugin_type == Unknown_sc_only
    - mode == "unsupported": 使用 type_support.support_category == unsupported
    - mode == "none": 不做过滤

    filter_scope 进一步约束过滤范围：
    - "unsupported_all": 按上述规则丢弃所有命中的细胞（当前默认行为）
    - "missing_only": 只丢弃 orig_type 属于 missing_type 的细胞，用于“只针对缺失类型止漏”
    """
    if filter_scope == "missing_only":
        if missing_type is None or not str(missing_type).strip():
            raise ValueError("filter_scope=missing_only requires non-empty missing_type.")
    if mode == "plugin_unknown":
        base_mask = relabel["plugin_type"] == "Unknown_sc_only"
    elif mode == "unsupported":
        col_name = "orig_type_canonical" if "orig_type_canonical" in relabel.columns else "orig_type"
        ts_col = "orig_type_canonical" if type_support is not None and "orig_type_canonical" in type_support.columns else "orig_type"
        unsupported_types = (
            set(type_support.loc[type_support["support_category"] == "unsupported", ts_col])
            if type_support is not None
            else set()
        )
        base_mask = relabel[col_name].isin(unsupported_types)
    else:
        base_mask = pd.Series(False, index=relabel.index)

    # 统计所有被Stage3标记为"低支持/Unknown"的规模（与filter_scope无关）
    col_canonical = "orig_type_canonical" if "orig_type_canonical" in relabel.columns else "orig_type"
    missing_mask = relabel[col_canonical] == missing_type if missing_type is not None else pd.Series(False, index=relabel.index)
    marked_total = int(base_mask.sum())
    marked_missing = int((base_mask & missing_mask).sum())
    marked_non_missing = marked_total - marked_missing

    # 进一步按scope约束：missing_only 只动指定missing_type，其它类型保留
    if filter_scope == "missing_only" and missing_type is not None:
        # missing_only 模式：直接过滤指定的 missing_type，不管它是否被标记为 unknown
        # 这样可以处理 missing_type 被标记为 weak 但仍需过滤的情况
        mask_drop = missing_mask
    else:
        # 默认行为：所有unsupported/Unknown_sc_only都被丢弃
        mask_drop = base_mask

    filtered = relabel.loc[~mask_drop].copy()
    n_filtered = int(mask_drop.sum())

    mark_stats = {
        "marked_unknown_total": marked_total,
        "marked_unknown_missing": marked_missing,
        "marked_unknown_non_missing": marked_non_missing,
    }
    return filtered, n_filtered, mark_stats


def prepare_cytospace_inputs(
    sc_expr: pd.DataFrame,
    st_expr: pd.DataFrame,
    st_coords: pd.DataFrame,
    sc_meta_stage1: pd.DataFrame,
    relabel_filtered: pd.DataFrame,
    work_dir: Path,
    cells_per_spot_override: int | None = None,
    cell_type_column: str = "plugin_type",
) -> dict[str, Path]:
    work_dir.mkdir(parents=True, exist_ok=True)

    # align SC expression to filtered cells
    # 如果使用sc_meta（baseline），直接使用所有SC细胞（不过滤）
    if cell_type_column == "sc_meta":
        keep_ids = [cid for cid in sc_meta_stage1["cell_id"] if cid in sc_expr.index]
    else:
        keep_ids = [cid for cid in relabel_filtered["cell_id"] if cid in sc_expr.index]
    sc_expr_filtered = sc_expr.loc[keep_ids]

    # CytoSPACE expects genes as rows, samples as columns
    sc_expr_gene_by_cell = sc_expr_filtered.T
    st_expr_gene_by_spot = st_expr.T

    sc_expr_path = work_dir / "sc_expression_for_cytospace.csv"
    st_expr_path = work_dir / "st_expression_for_cytospace.csv"
    coords_path = work_dir / "st_coordinates_for_cytospace.csv"
    cell_type_path = work_dir / "cell_types_for_cytospace.csv"
    fractions_path = work_dir / "cell_type_fractions.csv"
    n_cells_path = work_dir / "n_cells_per_spot.csv"

    sc_expr_gene_by_cell.to_csv(sc_expr_path, sep=",")
    st_expr_gene_by_spot.to_csv(st_expr_path, sep=",")
    coords_for_cyto = st_coords.copy()
    coords_for_cyto.columns = ["row", "col"]
    coords_for_cyto.to_csv(coords_path, sep=",", index_label="SpotID")

    if cell_type_column == "sc_meta":
        if "cell_type" not in sc_meta_stage1.columns:
            raise ValueError("sc_metadata.csv missing cell_type column")
        cell_type_df = sc_meta_stage1[["cell_id", "cell_type"]].copy()
    else:
        if cell_type_column not in relabel_filtered.columns:
            raise ValueError(f"cell_type_column={cell_type_column} not found in relabel file")
        cell_type_df = relabel_filtered[["cell_id", cell_type_column]].copy()
        cell_type_df.rename(columns={cell_type_column: "cell_type"}, inplace=True)
    cell_type_df.set_index("cell_id", inplace=True)
    cell_type_df.to_csv(cell_type_path, header=["cell_type"], sep=",")

    counts = cell_type_df["cell_type"].value_counts()
    frac_series = counts / counts.sum()
    frac_df = pd.DataFrame([frac_series.values], columns=frac_series.index)
    frac_df.index = ["Fraction"]
    frac_df.to_csv(fractions_path, sep=",")

    inputs = {
        "sc_expr": sc_expr_path,
        "st_expr": st_expr_path,
        "coords": coords_path,
        "cell_types": cell_type_path,
        "fractions": fractions_path,
    }
    if cells_per_spot_override is not None:
        df_n = pd.DataFrame({"n_cells": [cells_per_spot_override] * len(st_coords)}, index=st_coords.index)
        df_n.to_csv(n_cells_path)
        inputs["n_cells_per_spot"] = n_cells_path
    return inputs


def run_cytospace(
    inputs: dict[str, Path],
    output_dir: Path,
    project_root: Path,
    args: argparse.Namespace,
):
    sys.path.insert(0, str(project_root / "external" / "cytospace"))
    from cytospace.cytospace import main_cytospace  # type: ignore

    output_dir.mkdir(parents=True, exist_ok=True)
    main_cytospace(
        scRNA_path=str(inputs["sc_expr"]),
        cell_type_path=str(inputs["cell_types"]),
        n_cells_per_spot_path=str(inputs.get("n_cells_per_spot")) if "n_cells_per_spot" in inputs else None,
        st_cell_type_path=None,
        cell_type_fraction_estimation_path=str(inputs["fractions"]),
        spaceranger_path=None,
        st_path=str(inputs["st_expr"]),
        coordinates_path=str(inputs["coords"]),
        output_folder=str(output_dir),
        output_prefix="",
        mean_cell_numbers=args.mean_cell_numbers,
        downsample_off=False,
        scRNA_max_transcripts_per_cell=1500,
        solver_method=args.solver_method,
        distance_metric=args.distance_metric,
        sampling_method="duplicates",
        single_cell=False,
        number_of_selected_spots=10000,
        sampling_sub_spots=args.sampling_sub_spots,
        number_of_selected_sub_spots=args.n_subspots,
        number_of_processors=args.n_processors,
        seed=args.seed,
        plot_off=True,
        geometry="honeycomb",
        max_num_cells_plot=50000,
        num_column=3,
    )


def build_stage4_outputs(
    output_dir: Path,
    relabel_filtered: pd.DataFrame,
    missing_type: str,
    n_before: int,
    n_after: int,
    n_filtered: int,
    filter_mode: str,
    weak_th: float | None,
    strong_th: float | None,
    sampling_sub_spots: bool,
    n_subspots: int,
    solver_method: str,
    seed: int,
    missing_truth: str | None,
    project_root: Path,
    sample: str,
    cps_override: int | None,
    cell_type_column: str,
    filter_scope: FilterScope,
    mark_stats: dict,
) -> dict:
    assigned_path = output_dir / "assigned_locations.csv"
    if assigned_path.exists():
        assigned = pd.read_csv(assigned_path)
    else:
        assigned = pd.DataFrame()

    # Build simplified cell_assignment
    if not assigned.empty:
        assignment = assigned[["OriginalCID", "SpotID", "CellType"]].copy()
        assignment.rename(
            columns={"OriginalCID": "cell_id", "SpotID": "assigned_spot", "CellType": "cell_type"},
            inplace=True,
        )
        assignment.to_csv(output_dir / "cell_assignment.csv", index=False)
    else:
        assignment = pd.DataFrame(columns=["cell_id", "assigned_spot", "cell_type"])
        assignment.to_csv(output_dir / "cell_assignment.csv", index=False)

    # read composition to check missing type presence
    comp_path = output_dir / "cell_type_assignments_by_spot.csv"
    missing_present = 0
    if comp_path.exists():
        comp = pd.read_csv(comp_path, index_col=0)
        if missing_type in comp.columns:
            missing_present = int(comp[missing_type].sum())
    ca_rows = len(assignment)
    comp_total = int(comp["Total cells"].sum()) if comp_path.exists() else None

    summary = {
        "n_cells_before_prefilter": int(n_before),
        "n_cells_after_prefilter": int(n_after),
        "n_filtered": int(n_filtered),
        "assignment_rows": int(ca_rows),
        "assigned_locations_rows": int(len(assigned)) if not assigned.empty else 0,
        "by_spot_total_cells": comp_total,
        "missing_present_in_assignment": int(missing_present),
        "filter_mode": filter_mode,
        "stage3_weak_th": weak_th,
        "stage3_strong_th": strong_th,
        "sampling_sub_spots": bool(sampling_sub_spots),
        "n_subspots": int(n_subspots),
        "solver_method": solver_method,
        "seed": int(seed),
        "missing_type_truth": missing_truth,
        "cells_per_spot_override": cps_override,
        "cell_type_column": cell_type_column,
        "filter_scope": filter_scope,
        # 这些统计描述的是：在Stage3视角下被标记为“低支持/Unknown”的规模
        # 注意：与实际被Stage4丢弃的 n_filtered 相互独立
        "marked_unknown_total": mark_stats.get("marked_unknown_total"),
        "marked_unknown_missing": mark_stats.get("marked_unknown_missing"),
        "marked_unknown_non_missing": mark_stats.get("marked_unknown_non_missing"),
        "sha1": {
            "cell_assignment": sha1(output_dir / "cell_assignment.csv"),
            "sim_info": sha1(project_root / "data" / "processed" / sample / "stage1_preprocess" / "exported" / "sim_info.json"),
            "sc_metadata": sha1(project_root / "data" / "processed" / sample / "stage1_preprocess" / "exported" / "sc_metadata.csv"),
        },
        "cell_assignment_path": str(output_dir / "cell_assignment.csv"),
        "cytospace_output": str(output_dir),
    }
    (output_dir / "stage4_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return summary


def main():
    args = parse_args()

    # 强校验：missing_only 必须提供 non-empty missing_type
    if args.filter_scope == "missing_only":
        if args.missing_type is None or not str(args.missing_type).strip():
            raise SystemExit("filter_scope=missing_only requires non-empty --missing_type.")

    project_root = Path(args.project_root) if args.project_root else detect_root()

    # 加载别名映射，并对 missing_type 做 canonicalize
    alias_map_path = project_root / "configs" / "type_aliases.yaml"
    alias_map = load_alias_map(alias_map_path)
    missing_type_canonical = canonicalize_type_name(args.missing_type, alias_map)

    # 日志中显式打印 canonicalized missing_type 以及是否命中 alias
    orig_norm = normalize_type_name(args.missing_type)
    if missing_type_canonical != orig_norm:
        print(
            f"[stage4] missing_type canonicalized: "
            f"raw='{args.missing_type}' -> normalized='{orig_norm}' -> canonical='{missing_type_canonical}' (via alias map)"
        )
    else:
        print(
            f"[stage4] missing_type canonicalized: "
            f"raw='{args.missing_type}' -> normalized='{orig_norm}' (no alias mapping applied)"
        )

    sc_expr, st_expr, st_coords, sc_meta = load_stage1(project_root, args.sample)
    print(f"[debug] sc_expr shape={sc_expr.shape}, st_expr shape={st_expr.shape}")
    print(f"[debug] st_expr index sample={list(st_expr.index[:3])}")
    print(f"[debug] st_expr columns sample={list(st_expr.columns[:3])}")
    print(f"[debug] st_coords index sample={list(st_coords.index[:3])}")
    
    # Baseline (filter_mode=none) 不需要Stage3输出，直接使用所有SC细胞
    if args.filter_mode == "none":
        # 对于baseline，使用所有SC细胞，不进行任何过滤
        n_before = len(sc_meta)
        n_after = len(sc_meta)
        n_filtered = 0
        mark_stats = {
            "marked_unknown_total": 0,
            "marked_unknown_missing": 0,
            "marked_unknown_non_missing": 0,
        }
        # 创建一个虚拟的relabel_filtered，用于prepare_cytospace_inputs
        # 但实际上cell_type_column=sc_meta时，不会使用这个
        relabel_filtered = sc_meta[["cell_id"]].copy()
    else:
        # Route2需要Stage3输出进行过滤
        relabel, type_support = load_stage3(project_root, args.sample)

        # 对 orig_type 做 canonicalize，便于后续与 missing_type 一致匹配
        if "orig_type" in relabel.columns:
            relabel["orig_type_canonical"] = relabel["orig_type"].apply(
                lambda x: canonicalize_type_name(x, alias_map)
            )
        if type_support is not None and "orig_type" in type_support.columns:
            type_support["orig_type_canonical"] = type_support["orig_type"].apply(
                lambda x: canonicalize_type_name(x, alias_map)
            )
        n_before = len(relabel)

        relabel_filtered, n_filtered, mark_stats = apply_filter(
            relabel,
            type_support,
            args.filter_mode,
            missing_type=missing_type_canonical,
            filter_scope=args.filter_scope,
        )
        n_after = len(relabel_filtered)

    out_root = project_root / "result" / args.sample / "stage4_cytospace"
    prep_dir = out_root / "cytospace_input"
    cyto_out_dir = out_root / "cytospace_output"

    cps_override = load_cells_per_spot_override(project_root, args.sample)
    inputs = prepare_cytospace_inputs(
        sc_expr,
        st_expr,
        st_coords,
        sc_meta,
        relabel_filtered,
        prep_dir,
        cells_per_spot_override=cps_override,
        cell_type_column=args.cell_type_column,
    )
    run_cytospace(inputs, cyto_out_dir, project_root, args)
    weak_th, strong_th = load_stage3_thresholds(project_root, args.sample)
    missing_truth = load_missing_truth(project_root, args.sample)
    summary = build_stage4_outputs(
        cyto_out_dir,
        relabel_filtered,
        args.missing_type,
        n_before,
        n_after,
        n_filtered,
        args.filter_mode,
        weak_th,
        strong_th,
        args.sampling_sub_spots,
        args.n_subspots,
        args.solver_method,
        args.seed,
        missing_truth,
        project_root,
        args.sample,
        cps_override,
        args.cell_type_column,
        args.filter_scope,
        mark_stats,
    )

    print("[Stage4 cytospace] summary:")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

