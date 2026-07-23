#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


EPS = 1.0e-8


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Build a melanoma Fig.2c sample with spatial support removed for one target type "
            "while keeping the scRNA reference unchanged."
        )
    )
    p.add_argument("--project_root", default=".")
    p.add_argument("--source_sample", default="cytospace_fig2c_melanoma_mel1_rep2")
    p.add_argument("--target_sample", default="cytospace_fig2c_melanoma_mel1_rep2_profile_mask_bcells")
    p.add_argument("--target_type", default="B cells")
    p.add_argument("--cell_type_col", default="cell_type")
    p.add_argument("--overwrite", action="store_true")

    p.add_argument("--max_genes", type=int, default=240)
    p.add_argument("--min_genes", type=int, default=25)
    p.add_argument("--min_target_mean", type=float, default=0.001)
    p.add_argument("--min_st_detect_frac", type=float, default=0.01)
    p.add_argument("--min_specificity_all", type=float, default=1.2)
    p.add_argument("--min_specificity_neighbor", type=float, default=1.15)
    p.add_argument("--min_specificity_max", type=float, default=1.05)
    p.add_argument("--specificity_cap", type=float, default=6.0)

    p.add_argument("--target_profile_impact", type=float, default=0.40)
    p.add_argument("--mask_strength", type=float, default=0.95)
    p.add_argument("--min_scale", type=float, default=0.01)
    p.add_argument("--hotspot_quantile", type=float, default=0.80)
    p.add_argument("--hotspot_mask_strength", type=float, default=0.98)
    p.add_argument("--hotspot_min_scale", type=float, default=0.01)
    p.add_argument("--disable_hotspot", action="store_true")
    p.add_argument("--component_profile_weight", type=float, default=0.0)
    p.add_argument("--component_hotspot_boost", type=float, default=1.5)
    p.add_argument("--component_score_quantile", type=float, default=0.95)

    p.add_argument("--min_impact_ratio", type=float, default=1.5)
    p.add_argument("--max_non_target_total_impact", type=float, default=0.18)
    p.add_argument("--allow_unsafe_panel", action="store_true")
    p.add_argument("--enable_auto_missing", action="store_true")
    return p.parse_args()


def _read_yaml(path: Path) -> dict:
    if not path.exists():
        return {}
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def _write_yaml(path: Path, cfg: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(cfg, sort_keys=False, allow_unicode=True), encoding="utf-8")


def _cosine_neighbor(type_means: pd.DataFrame, target_type: str) -> str:
    target = type_means.loc[target_type].to_numpy(dtype=np.float64)
    target_norm = float(np.linalg.norm(target)) + EPS
    best_type = ""
    best_score = -2.0
    for cell_type, row in type_means.iterrows():
        if str(cell_type) == target_type:
            continue
        vec = row.to_numpy(dtype=np.float64)
        score = float(np.dot(target, vec) / (target_norm * (np.linalg.norm(vec) + EPS)))
        if score > best_score:
            best_score = score
            best_type = str(cell_type)
    return best_type


def _build_gene_panel(
    sc: pd.DataFrame,
    cell_types: pd.Series,
    st: pd.DataFrame,
    target_type: str,
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    available_types = sorted(cell_types.astype(str).unique().tolist())
    if target_type not in available_types:
        raise ValueError(f"target_type not found: {target_type}. Available types: {available_types}")

    common_genes = [g for g in sc.columns if g in st.columns]
    if not common_genes:
        raise ValueError("No shared genes between source scRNA and ST exports.")
    sc = sc.loc[:, common_genes]
    st_common = st.loc[:, common_genes]

    type_means = sc.groupby(cell_types).mean(numeric_only=True).astype(np.float64)
    target_mean = type_means.loc[target_type]
    other_type_means = type_means.drop(index=target_type)
    all_other_mean = sc.loc[cell_types != target_type].mean(axis=0).astype(np.float64)
    max_other_mean = other_type_means.max(axis=0).astype(np.float64)
    max_other_type = other_type_means.idxmax(axis=0).astype(str)
    neighbor_type = _cosine_neighbor(type_means.loc[:, common_genes], target_type)
    neighbor_mean = type_means.loc[neighbor_type].astype(np.float64)
    st_detect_frac = (st_common > 0).mean(axis=0).astype(np.float64)

    specificity_all = (target_mean + EPS) / (all_other_mean + EPS)
    specificity_max = (target_mean + EPS) / (max_other_mean + EPS)
    specificity_neighbor = (target_mean + EPS) / (neighbor_mean + EPS)
    specificity_floor = pd.concat(
        [specificity_all, specificity_max, specificity_neighbor],
        axis=1,
    ).min(axis=1)

    panel = pd.DataFrame(
        {
            "gene": common_genes,
            "target_type": target_type,
            "neighbor_type": neighbor_type,
            "target_mean": target_mean.loc[common_genes].to_numpy(),
            "all_other_mean": all_other_mean.loc[common_genes].to_numpy(),
            "max_other_mean": max_other_mean.loc[common_genes].to_numpy(),
            "max_other_type": max_other_type.loc[common_genes].to_numpy(),
            "neighbor_mean": neighbor_mean.loc[common_genes].to_numpy(),
            "specificity_all": specificity_all.loc[common_genes].to_numpy(),
            "specificity_max": specificity_max.loc[common_genes].to_numpy(),
            "specificity_neighbor": specificity_neighbor.loc[common_genes].to_numpy(),
            "specificity_floor": specificity_floor.loc[common_genes].to_numpy(),
            "st_detect_frac": st_detect_frac.loc[common_genes].to_numpy(),
        }
    ).set_index("gene", drop=False)

    keep = (
        (panel["target_mean"] >= args.min_target_mean)
        & (panel["st_detect_frac"] >= args.min_st_detect_frac)
        & (panel["specificity_all"] >= args.min_specificity_all)
        & (panel["specificity_neighbor"] >= args.min_specificity_neighbor)
        & (panel["specificity_max"] >= args.min_specificity_max)
    )
    panel = panel.loc[keep].copy()
    if len(panel) < args.min_genes:
        raise ValueError(
            f"Only {len(panel)} target genes passed filters; need at least {args.min_genes}."
        )

    panel["score"] = panel["target_mean"] * np.log1p(panel["specificity_floor"])
    panel = panel.sort_values(["score", "specificity_floor", "target_mean"], ascending=False).head(args.max_genes)
    panel.attrs["target_profile_full"] = type_means.loc[target_type].copy()

    _assign_mask_scales(type_means.loc[:, common_genes], panel, target_type, args)

    impact = _build_impact_audit(type_means.loc[:, common_genes], panel, target_type)
    target_impact = float(impact.loc[target_type, "plugin_profile_impact_fraction"])
    non_target = impact.drop(index=target_type)
    max_non_target_impact = (
        float(non_target["plugin_profile_impact_fraction"].max()) if not non_target.empty else 0.0
    )
    impact_ratio = target_impact / (max_non_target_impact + EPS)
    summary = {
        "target_type": target_type,
        "neighbor_type": neighbor_type,
        "n_panel_genes": int(len(panel)),
        "target_total_impact": target_impact,
        "max_non_target_total_impact": max_non_target_impact,
        "impact_ratio": impact_ratio,
        "target_profile_impact_requested": float(args.target_profile_impact),
        "target_profile_impact_achieved": target_impact,
    }

    if not args.allow_unsafe_panel:
        if impact_ratio < args.min_impact_ratio:
            raise ValueError(
                f"Unsafe panel: target/non-target impact ratio {impact_ratio:.3f} < {args.min_impact_ratio:.3f}"
            )
        if max_non_target_impact > args.max_non_target_total_impact:
            raise ValueError(
                f"Unsafe panel: max non-target total impact {max_non_target_impact:.3f} > "
                f"{args.max_non_target_total_impact:.3f}"
            )

    return panel, impact, summary


def _assign_mask_scales(
    type_means: pd.DataFrame,
    panel: pd.DataFrame,
    target_type: str,
    args: argparse.Namespace,
) -> None:
    denom = max(args.specificity_cap - 1.0, EPS)
    weights = ((panel["specificity_floor"].astype(float) - 1.0) / denom).clip(0.0, 1.0)
    if float(weights.sum()) <= EPS:
        weights = pd.Series(1.0, index=panel.index)

    max_reduction = max(0.0, min(1.0, 1.0 - args.min_scale))
    if args.target_profile_impact and args.target_profile_impact > 0:
        reductions = _solve_reductions_for_target_impact(
            type_means=type_means,
            genes=panel.index.tolist(),
            target_type=target_type,
            weights=weights.to_numpy(dtype=np.float64),
            target_impact=float(args.target_profile_impact),
            max_reduction=max_reduction,
        )
    else:
        reductions = np.minimum(max_reduction, float(args.mask_strength) * weights.to_numpy(dtype=np.float64))

    hotspot_max_reduction = max(0.0, min(1.0, 1.0 - args.hotspot_min_scale))
    hotspot_reductions = np.maximum(
        reductions,
        np.minimum(
            hotspot_max_reduction,
            float(args.hotspot_mask_strength) * np.maximum(weights.to_numpy(dtype=np.float64), reductions),
        ),
    )

    panel["specificity_weight"] = weights.to_numpy(dtype=np.float64)
    panel["global_reduction"] = reductions
    panel["global_scale"] = 1.0 - reductions
    panel["hotspot_reduction"] = hotspot_reductions
    panel["hotspot_scale"] = 1.0 - hotspot_reductions


def _solve_reductions_for_target_impact(
    type_means: pd.DataFrame,
    genes: list[str],
    target_type: str,
    weights: np.ndarray,
    target_impact: float,
    max_reduction: float,
) -> np.ndarray:
    target_profile = type_means.loc[target_type, genes].to_numpy(dtype=np.float64)
    target_total = float(type_means.loc[target_type].sum()) + EPS

    def impact(multiplier: float) -> float:
        reductions = np.minimum(max_reduction, multiplier * weights)
        return float(np.dot(target_profile, reductions) / target_total)

    max_possible = impact(1.0e9)
    if max_possible <= target_impact:
        return np.minimum(max_reduction, 1.0e9 * weights)

    lo = 0.0
    hi = 1.0
    while impact(hi) < target_impact:
        hi *= 2.0
        if hi > 1.0e9:
            break
    for _ in range(80):
        mid = (lo + hi) / 2.0
        if impact(mid) < target_impact:
            lo = mid
        else:
            hi = mid
    return np.minimum(max_reduction, hi * weights)


def _build_impact_audit(type_means: pd.DataFrame, panel: pd.DataFrame, target_type: str) -> pd.DataFrame:
    genes = panel.index.tolist()
    reductions = 1.0 - panel["global_scale"].astype(float)
    rows = []
    for cell_type, profile in type_means.iterrows():
        profile = profile.astype(float)
        panel_expr = profile.loc[genes]
        panel_total = float(panel_expr.sum())
        plugin_total = float(profile.sum())
        reduction_mass = float((panel_expr * reductions).sum())
        rows.append(
            {
                "cell_type": cell_type,
                "is_target": str(cell_type) == target_type,
                "panel_expression_sum": panel_total,
                "plugin_expression_sum": plugin_total,
                "panel_impact_fraction": reduction_mass / (panel_total + EPS),
                "plugin_profile_impact_fraction": reduction_mass / (plugin_total + EPS),
            }
        )
    out = pd.DataFrame(rows).set_index("cell_type")
    target_total = float(out.loc[target_type, "plugin_profile_impact_fraction"])
    out["relative_to_target_total_impact"] = out["plugin_profile_impact_fraction"] / (target_total + EPS)
    return out.sort_values("plugin_profile_impact_fraction", ascending=False)


def _apply_mask(st: pd.DataFrame, panel: pd.DataFrame, args: argparse.Namespace) -> tuple[pd.DataFrame, dict]:
    masked = st.copy()
    genes = panel.index.tolist()
    expr = masked.loc[:, genes].astype(np.float64)

    if args.disable_hotspot:
        hotspot_spots: list[str] = []
        scales = panel["global_scale"].astype(float)
        masked.loc[:, genes] = expr.mul(scales, axis=1)
    else:
        weights = panel["specificity_weight"].to_numpy(dtype=np.float64)
        values = expr.to_numpy(dtype=np.float64)
        norm_values = values / (values.mean(axis=1, keepdims=True) + EPS)
        spot_scores = (norm_values @ weights) / (weights.sum() + EPS)
        threshold = float(np.quantile(spot_scores, args.hotspot_quantile))
        hotspot_mask = spot_scores >= threshold
        hotspot_spots = expr.index[hotspot_mask].astype(str).tolist()

        global_scales = panel["global_scale"].to_numpy(dtype=np.float64)[None, :]
        hotspot_scales = panel["hotspot_scale"].to_numpy(dtype=np.float64)[None, :]
        scale_matrix = np.repeat(global_scales, expr.shape[0], axis=0)
        scale_matrix[hotspot_mask, :] = np.repeat(hotspot_scales, int(hotspot_mask.sum()), axis=0)
        masked.loc[:, genes] = values * scale_matrix

    component_subtraction_total = 0.0
    if float(args.component_profile_weight) > 0:
        target_profile = panel.attrs.get("target_profile_full")
        if isinstance(target_profile, pd.Series) and not target_profile.empty:
            full_expr = masked.astype(np.float64)
            target_profile = target_profile.reindex(full_expr.columns).fillna(0.0).astype(np.float64)
            target_profile_sum = float(target_profile.sum())
            if target_profile_sum > EPS:
                target_profile = target_profile / target_profile_sum
                spot_score_scale = float(np.quantile(spot_scores, args.component_score_quantile)) + EPS
                spot_score_norm = np.clip(spot_scores / spot_score_scale, 0.0, 1.0)
                hotspot_boost = np.ones_like(spot_score_norm)
                hotspot_boost[np.asarray(spot_scores >= float(np.quantile(spot_scores, args.hotspot_quantile)))] = float(
                    args.component_hotspot_boost
                )
                component_strength = np.clip(
                    spot_score_norm * hotspot_boost * float(args.component_profile_weight),
                    0.0,
                    1.0,
                )
                row_totals = full_expr.sum(axis=1).to_numpy(dtype=np.float64)
                subtraction = np.outer(component_strength * row_totals, target_profile.to_numpy(dtype=np.float64))
                new_values = np.maximum(full_expr.to_numpy(dtype=np.float64) - subtraction, 0.0)
                component_subtraction_total = float(np.maximum(full_expr.to_numpy(dtype=np.float64) - new_values, 0.0).sum())
                masked = pd.DataFrame(new_values, index=full_expr.index, columns=full_expr.columns)

    info = {
        "hotspot_enabled": not args.disable_hotspot,
        "hotspot_quantile": None if args.disable_hotspot else args.hotspot_quantile,
        "n_hotspot_spots": int(len(hotspot_spots)),
        "component_profile_weight": float(args.component_profile_weight),
        "component_hotspot_boost": float(args.component_hotspot_boost),
        "component_score_quantile": float(args.component_score_quantile),
        "component_subtraction_total": component_subtraction_total,
    }
    return masked, info


def _write_target_config(project_root: Path, source_sample: str, target_sample: str, args: argparse.Namespace) -> None:
    source_cfg_path = project_root / "configs" / "datasets" / f"{source_sample}.yaml"
    cfg = _read_yaml(source_cfg_path)
    if not cfg:
        raise FileNotFoundError(source_cfg_path)

    stage3 = cfg.setdefault("stage3", {})
    stage3.pop("force_unsupported_types", None)
    stage3["plugin_genes_path"] = (
        f"data/processed/cytospace_fig2c_melanoma/{target_sample}/stage1_preprocess/hvg_genes.txt"
    )
    if args.enable_auto_missing:
        amd = stage3.setdefault("auto_missing_detection", {})
        amd["enable"] = True
        amd["method"] = "adaptive_low_support"
        amd["min_cells"] = 50
        amd["robust_z_th"] = -2.0
        amd["soft_z_th"] = -0.9
        amd["require_masked_for_soft"] = False
        amd["require_masked_for_hard"] = True
        amd["max_fraction_types"] = 0.3
        amd["max_types"] = 2
        amd["action"] = "mark_unknown"
        amd["require_confirmation"] = True
        amd["confirmation_max_support_score"] = 0.55
        amd["confirmation_use_masked_missing"] = True
        amd["confirmation_use_marker_identity"] = True
        amd["confirmation_marker_identity_z_th"] = -1.5
        amd["confirmation_marker_support_score_th"] = 0.5
        stage3["masked_missing_detection"] = {
            "enable": True,
            "apply_to_auto_missing": True,
            "neighbor_cosine_th": 0.89,
            "neighbor_cell_ratio_min": 0.50,
            "marker_top_n": 20,
            "min_identity_markers": 3,
            "min_marker_specificity": 1.2,
            "min_marker_type_mean": 0.001,
            "min_marker_st_detect_frac": 0.005,
            "st_presence_quantile": 0.90,
            "identity_z_th": -0.70,
            "pressure_z_th": 0.0,
            "min_support_score_for_apply": 0.52,
            "max_types": 2,
        }
        stage3["marker_identity_diagnostics"] = {
            "enable": True,
            "marker_top_n": 80,
            "min_identity_markers": 5,
            "min_all_specificity": 1.3,
            "min_neighbor_specificity": 1.1,
            "min_marker_type_mean": 0.001,
            "min_marker_st_detect_frac": 0.005,
            "st_presence_quantile": 0.9,
            "depleted_z_th": -0.9,
        }
    _write_yaml(project_root / "configs" / "datasets" / f"{target_sample}.yaml", cfg)


def main() -> int:
    args = parse_args()
    project_root = Path(args.project_root).resolve()

    source_dir = (
        project_root
        / "data"
        / "processed"
        / "cytospace_fig2c_melanoma"
        / args.source_sample
        / "stage1_preprocess"
    )
    target_dir = (
        project_root
        / "data"
        / "processed"
        / "cytospace_fig2c_melanoma"
        / args.target_sample
        / "stage1_preprocess"
    )
    source_export = source_dir / "exported"
    target_export = target_dir / "exported"
    if not source_export.exists():
        raise FileNotFoundError(source_export)
    if target_dir.exists():
        if not args.overwrite:
            raise FileExistsError(f"Target exists: {target_dir}")
        shutil.rmtree(target_dir.parent)
    shutil.copytree(source_dir.parent, target_dir.parent)

    sc = pd.read_csv(target_export / "sc_expression_normalized.csv", index_col=0)
    st = pd.read_csv(target_export / "st_expression_normalized.csv", index_col=0)
    sc_meta = pd.read_csv(target_export / "sc_metadata.csv")
    cell_types = sc_meta.set_index("cell_id")[args.cell_type_col].astype(str)
    sc.index = sc.index.astype(str)
    st.index = st.index.astype(str)
    sc.columns = sc.columns.astype(str)
    st.columns = st.columns.astype(str)
    cell_types.index = cell_types.index.astype(str)
    cell_types = cell_types.reindex(sc.index)

    panel, impact, impact_summary = _build_gene_panel(sc, cell_types, st, args.target_type, args)
    masked_st, hotspot_info = _apply_mask(st, panel, args)
    masked_st.to_csv(target_export / "st_expression_normalized.csv", index_label="spot_id")

    panel_path = target_dir / "st_profile_mask_gene_panel.csv"
    impact_path = target_dir / "st_profile_mask_impact_audit.csv"
    info_path = target_dir / "fig2c_profile_mask_info.json"
    panel.to_csv(panel_path, index=False)
    impact.to_csv(impact_path)

    stage1_summary_path = target_dir / "fig2c_stage1_summary.json"
    if stage1_summary_path.exists():
        stage1_summary = json.loads(stage1_summary_path.read_text(encoding="utf-8"))
    else:
        stage1_summary = {}
    stage1_summary["sample"] = args.target_sample
    stage1_summary["source_sample"] = args.source_sample
    stage1_summary["st_profile_mask"] = {
        "target_type": args.target_type,
        "params": vars(args),
        "impact_summary": impact_summary,
        **hotspot_info,
        "gene_panel_csv": str(panel_path.relative_to(project_root)).replace("\\", "/"),
        "impact_audit_csv": str(impact_path.relative_to(project_root)).replace("\\", "/"),
    }
    stage1_summary_path.write_text(json.dumps(stage1_summary, indent=2, ensure_ascii=False), encoding="utf-8")
    info_path.write_text(json.dumps(stage1_summary["st_profile_mask"], indent=2, ensure_ascii=False), encoding="utf-8")

    _write_target_config(project_root, args.source_sample, args.target_sample, args)
    print(f"[DONE] built profile-mask sample: {args.target_sample}")
    print(json.dumps(stage1_summary["st_profile_mask"]["impact_summary"], indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
