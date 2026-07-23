from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.build_cytospace_fig2c_profile_mask_sample import (
    _apply_mask,
    _build_gene_panel,
    _read_yaml,
    _write_yaml,
)


def _safe_name(label: str) -> str:
    return "".join(ch.lower() if ch.isalnum() else "_" for ch in str(label)).strip("_")


def _candidate_processed_dirs(root: Path, sample: str) -> list[Path]:
    return [
        root / "data" / "processed" / sample,
        root / "data" / "processed" / "cytospace_fig2c_melanoma" / sample,
        root / "data" / "processed" / "cytospace_fig2d_tme" / sample,
        root / "data" / "processed" / "cytospace_fig2d_profile_mask" / sample,
    ]


def _resolve_processed_dir(root: Path, sample: str) -> Path:
    for path in _candidate_processed_dirs(root, sample):
        if (path / "stage1_preprocess" / "exported").exists():
            return path.resolve()
    raise FileNotFoundError(f"Cannot find Stage1 export for sample: {sample}")


def _write_target_config(
    project_root: Path,
    source_sample: str,
    target_sample: str,
    target_group: str,
    hvg_path: Path,
) -> None:
    source_cfg_path = project_root / "configs" / "datasets" / f"{source_sample}.yaml"
    cfg = _read_yaml(source_cfg_path)
    if not cfg:
        raise FileNotFoundError(source_cfg_path)

    cfg["storage"] = {"group": target_group}
    stage3 = cfg.setdefault("stage3", {})
    stage3.pop("force_unsupported_types", None)
    stage3["plugin_genes_path"] = str(hvg_path.relative_to(project_root)).replace("\\", "/")

    amd = stage3.setdefault("auto_missing_detection", {})
    amd.update(
        {
            "enable": True,
            "method": "adaptive_low_support",
            "min_cells": 50,
            "robust_z_th": -2.0,
            "soft_z_th": -0.9,
            "require_masked_for_soft": False,
            "require_masked_for_hard": True,
            "max_fraction_types": 0.3,
            "max_types": 2,
            "action": "mark_unknown",
            "require_confirmation": True,
            "confirmation_max_support_score": 0.55,
            "confirmation_use_masked_missing": True,
            "confirmation_use_marker_identity": True,
            "confirmation_marker_identity_z_th": -1.5,
            "confirmation_marker_support_score_th": 0.5,
        }
    )
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


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Build a Fig.2d-style profile-mask sample by weakening ST support for one target type. "
            "No forced unsupported types are written; Stage3 must detect the missing/unsupported type automatically."
        )
    )
    p.add_argument("--project_root", default=".")
    p.add_argument("--source_sample", required=True)
    p.add_argument("--target_sample", default=None)
    p.add_argument("--target_group", default="cytospace_fig2d_profile_mask")
    p.add_argument("--target_type", required=True)
    p.add_argument("--cell_type_col", default="cell_type")
    p.add_argument("--overwrite", action="store_true")

    p.add_argument("--max_genes", type=int, default=240)
    p.add_argument("--min_genes", type=int, default=20)
    p.add_argument("--min_target_mean", type=float, default=0.001)
    p.add_argument("--min_st_detect_frac", type=float, default=0.005)
    p.add_argument("--min_specificity_all", type=float, default=1.15)
    p.add_argument("--min_specificity_neighbor", type=float, default=1.05)
    p.add_argument("--min_specificity_max", type=float, default=1.02)
    p.add_argument("--specificity_cap", type=float, default=6.0)
    p.add_argument("--target_profile_impact", type=float, default=0.45)
    p.add_argument("--mask_strength", type=float, default=0.95)
    p.add_argument("--min_scale", type=float, default=0.01)
    p.add_argument("--hotspot_quantile", type=float, default=0.80)
    p.add_argument("--hotspot_mask_strength", type=float, default=0.98)
    p.add_argument("--hotspot_min_scale", type=float, default=0.01)
    p.add_argument("--disable_hotspot", action="store_true")
    p.add_argument("--component_profile_weight", type=float, default=0.0)
    p.add_argument("--component_hotspot_boost", type=float, default=1.5)
    p.add_argument("--component_score_quantile", type=float, default=0.95)
    p.add_argument("--min_impact_ratio", type=float, default=1.25)
    p.add_argument("--max_non_target_total_impact", type=float, default=0.25)
    p.add_argument("--allow_unsafe_panel", action="store_true")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    root = Path(args.project_root).resolve()
    target_sample = args.target_sample or f"{args.source_sample}_profile_mask_{_safe_name(args.target_type)}"

    source_dir = _resolve_processed_dir(root, args.source_sample)
    target_root = root / "data" / "processed" / args.target_group / target_sample
    if target_root.exists():
        if not args.overwrite:
            raise FileExistsError(f"Target exists: {target_root}")
        shutil.rmtree(target_root)
    shutil.copytree(source_dir, target_root)

    stage1_dir = target_root / "stage1_preprocess"
    export = stage1_dir / "exported"
    sc = pd.read_csv(export / "sc_expression_normalized.csv", index_col=0)
    st = pd.read_csv(export / "st_expression_normalized.csv", index_col=0)
    sc_meta = pd.read_csv(export / "sc_metadata.csv")
    if args.cell_type_col not in sc_meta.columns:
        raise ValueError(f"cell type column not found: {args.cell_type_col}")

    cell_types = sc_meta.set_index("cell_id")[args.cell_type_col].astype(str)
    sc.index = sc.index.astype(str)
    st.index = st.index.astype(str)
    sc.columns = sc.columns.astype(str)
    st.columns = st.columns.astype(str)
    cell_types.index = cell_types.index.astype(str)
    cell_types = cell_types.reindex(sc.index)

    panel, impact, impact_summary = _build_gene_panel(sc, cell_types, st, args.target_type, args)
    masked_st, hotspot_info = _apply_mask(st, panel, args)
    masked_st.to_csv(export / "st_expression_normalized.csv", index_label="spot_id")

    panel_path = stage1_dir / "fig2d_profile_mask_gene_panel.csv"
    impact_path = stage1_dir / "fig2d_profile_mask_impact_audit.csv"
    info_path = stage1_dir / "fig2d_profile_mask_info.json"
    panel.to_csv(panel_path, index=False)
    impact.to_csv(impact_path)

    summary = {
        "sample": target_sample,
        "source_sample": args.source_sample,
        "profile_mask": {
            "target_type": args.target_type,
            "params": vars(args),
            "impact_summary": impact_summary,
            **hotspot_info,
            "gene_panel_csv": str(panel_path.relative_to(root)).replace("\\", "/"),
            "impact_audit_csv": str(impact_path.relative_to(root)).replace("\\", "/"),
            "force_unsupported_types": None,
        },
    }
    info_path.write_text(json.dumps(summary["profile_mask"], indent=2, ensure_ascii=False), encoding="utf-8")
    # Compatibility with existing profile-mask screening utilities.
    (stage1_dir / "fig2c_profile_mask_info.json").write_text(
        json.dumps(summary["profile_mask"], indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    (stage1_dir / "fig2d_profile_mask_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    _write_target_config(root, args.source_sample, target_sample, args.target_group, stage1_dir / "hvg_genes.txt")

    print(f"[DONE] built profile-mask sample: {target_sample}")
    print(json.dumps(summary["profile_mask"]["impact_summary"], indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
