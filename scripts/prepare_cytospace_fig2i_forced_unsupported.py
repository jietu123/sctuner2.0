from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


def _copy_if_exists(src: Path, dst: Path) -> None:
    if src.exists():
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create a Fig.2i kidney state-level forced-unsupported decoy sample."
    )
    parser.add_argument("--project_root", default=".")
    parser.add_argument("--source_sample", default="cytospace_fig2i_mouse_kidney_strict")
    parser.add_argument("--target_sample", default="cytospace_fig2i_mouse_kidney_forced_unsupported")
    parser.add_argument("--decoy_source_state", default="State_32")
    parser.add_argument("--decoy_type", default="Unsupported_decoy_State32like")
    parser.add_argument("--n_decoys", type=int, default=1200)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    root = Path(args.project_root).resolve()
    source_export = root / "data" / "processed" / args.source_sample / "stage1_preprocess" / "exported"
    target_stage1 = root / "data" / "processed" / args.target_sample / "stage1_preprocess"
    target_export = target_stage1 / "exported"
    target_stage3 = root / "data" / "processed" / args.target_sample / "stage3_typematch"
    target_result_stage3 = root / "result" / args.target_sample / "stage3_typematch"

    if not source_export.exists():
        raise FileNotFoundError(source_export)
    if target_stage1.exists() and not args.overwrite:
        raise FileExistsError(f"{target_stage1} exists; pass --overwrite")
    if target_stage1.exists():
        shutil.rmtree(target_stage1)
    if target_stage3.exists():
        shutil.rmtree(target_stage3)
    if target_result_stage3.exists():
        shutil.rmtree(target_result_stage3)
    target_export.mkdir(parents=True, exist_ok=True)
    target_stage3.mkdir(parents=True, exist_ok=True)
    target_result_stage3.mkdir(parents=True, exist_ok=True)

    for name in [
        "st_expression_normalized.csv",
        "st_coordinates.csv",
        "selected_genes.txt",
        "kidneycellexplorer_genes.txt",
        "kidneycellexplorer_gene_variance.csv",
        "r_export_stdout.log",
        "r_export_stderr.log",
    ]:
        _copy_if_exists(source_export / name, target_export / name)

    for name in ["common_genes.txt", "hvg_genes.txt"]:
        _copy_if_exists(source_export.parent / name, target_stage1 / name)

    meta = pd.read_csv(source_export / "sc_metadata.csv")
    expr = pd.read_csv(source_export / "sc_expression_normalized.csv", index_col=0)
    if "cell_type" not in meta.columns:
        raise ValueError("source sc_metadata.csv missing cell_type column")
    source_ids = meta.loc[meta["cell_type"].astype(str).eq(args.decoy_source_state), "cell_id"].astype(str).tolist()
    if not source_ids:
        raise ValueError(f"No cells found for decoy_source_state={args.decoy_source_state}")
    missing_expr_ids = sorted(set(source_ids).difference(expr.index.astype(str)))
    if missing_expr_ids:
        raise ValueError(f"{len(missing_expr_ids)} source state cells missing from expression matrix")

    rng = np.random.default_rng(args.seed)
    sampled = rng.choice(source_ids, size=args.n_decoys, replace=True)
    decoy_ids = [f"{args.decoy_type}_{i:04d}_{cid}" for i, cid in enumerate(sampled)]
    decoy_expr = expr.loc[sampled].copy()
    decoy_expr.index = decoy_ids
    expr_out = pd.concat([expr, decoy_expr], axis=0)

    template_meta = meta.set_index("cell_id").loc[sampled].reset_index(drop=True).copy()
    template_meta["cell_id"] = decoy_ids
    template_meta["cell_type"] = args.decoy_type
    if "state_id" in template_meta.columns:
        template_meta["state_id"] = args.decoy_type
    meta_out = pd.concat([meta, template_meta], axis=0, ignore_index=True)

    expr_out.to_csv(target_export / "sc_expression_normalized.csv", index_label="cell_id")
    meta_out.to_csv(target_export / "sc_metadata.csv", index=False)

    relabel = pd.DataFrame(
        {
            "cell_id": meta_out["cell_id"].astype(str),
            "orig_type": meta_out["cell_type"].astype(str),
            "plugin_type": meta_out["cell_type"].astype(str),
            "label_status": "kept",
        }
    )
    is_decoy = relabel["orig_type"].eq(args.decoy_type)
    relabel.loc[is_decoy, "plugin_type"] = "Unknown_sc_only"
    relabel.loc[is_decoy, "label_status"] = "forced_unsupported_decoy"
    relabel.to_csv(target_stage3 / "cell_type_relabel.csv", index=False)

    support_rows = []
    for cell_type, n in meta_out["cell_type"].value_counts().sort_index().items():
        support_rows.append(
            {
                "orig_type": cell_type,
                "n_cells": int(n),
                "support_category": "unsupported" if cell_type == args.decoy_type else "supported",
                "forced_decoy": bool(cell_type == args.decoy_type),
            }
        )
    type_support = pd.DataFrame(support_rows)
    type_support.to_csv(target_stage3 / "type_support.csv", index=False)

    summary = {
        "sample": args.target_sample,
        "source_sample": args.source_sample,
        "decoy_source_state": args.decoy_source_state,
        "decoy_type": args.decoy_type,
        "n_decoys": int(args.n_decoys),
        "n_cells_original": int(len(meta)),
        "n_cells_total": int(len(meta_out)),
        "n_forced_unknown": int(is_decoy.sum()),
        "plugin_types": sorted([x for x in relabel["plugin_type"].unique().tolist() if x != "Unknown_sc_only"])
        + ["Unknown_sc_only"],
        "unknown_overview": {
            "cell_fraction": float(is_decoy.mean()),
            "unknown_label_effective": "Unknown_sc_only",
        },
    }
    (target_result_stage3 / "stage3_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (target_stage3 / "stage3_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    cfg = {
        "paths": {"sc_expr": "unused", "sc_meta": "unused", "st_expr": "unused", "st_meta": "unused"},
        "qc": {
            "sc_min_genes": 0,
            "sc_max_genes": 100000,
            "sc_max_mt": 100,
            "st_min_genes": 0,
            "st_max_genes": "Inf",
            "st_max_mt": 100,
            "hvg_nfeatures": 2000,
            "mt_pattern": "^(MT-|mt-)",
        },
        "gene_filter": {"min_cells_sc": 0, "min_cells_st": 0},
        "stage3": {
            "strong_th": 0.7,
            "weak_th": 0.4,
            "unknown_floor": 0.3,
            "plugin_genes_path": str((target_stage1 / "hvg_genes.txt").relative_to(root)).replace("\\", "/"),
            "auto_missing_detection": {"enable": False, "action": "mark_unknown"},
            "masked_missing_detection": {"enable": False},
        },
    }
    cfg_path = root / "configs" / "datasets" / f"{args.target_sample}.yaml"
    cfg_path.parent.mkdir(parents=True, exist_ok=True)
    cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")

    print(f"[OK] wrote forced sample: {target_export}")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
