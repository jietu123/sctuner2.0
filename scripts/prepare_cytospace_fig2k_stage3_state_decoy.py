from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


def _write_config(root: Path, sample: str, hvg_path: Path) -> None:
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
            "plugin_genes_path": str(hvg_path.relative_to(root)).replace("\\", "/"),
            "mismatch_action": "mark_unknown",
            "auto_missing_detection": {
                "enable": True,
                "method": "adaptive_low_support",
                "min_cells": 50,
                "robust_z_th": -1.2,
                "soft_z_th": -0.6,
                "require_masked_for_soft": False,
                "require_masked_for_hard": False,
                "max_fraction_types": 0.5,
                "max_types": 2,
                "action": "mark_unknown",
                "require_confirmation": False,
            },
            "masked_missing_detection": {"enable": False},
            "marker_identity_diagnostics": {
                "enable": True,
                "marker_top_n": 80,
                "min_identity_markers": 5,
                "min_all_specificity": 1.2,
                "min_neighbor_specificity": 1.1,
                "min_marker_type_mean": 0.001,
                "min_marker_st_detect_frac": 0.005,
                "st_presence_quantile": 0.9,
                "depleted_z_th": -0.9,
            },
        },
    }
    cfg_path = root / "configs" / "datasets" / f"{sample}.yaml"
    cfg_path.parent.mkdir(parents=True, exist_ok=True)
    cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Create a Fig.2k state-decoy sample where Stage3 must detect the unsupported decoy "
            "from Stage1 expression support, without prewritten relabel/type_support outputs."
        )
    )
    parser.add_argument("--project_root", default=".")
    parser.add_argument(
        "--source_sample",
        required=True,
        help="Processed source sample containing the base Fig.2k stage1 export and decoy cells.",
    )
    parser.add_argument("--target_sample", default="cytospace_fig2k_breast_state_decoy_stage3_detected")
    parser.add_argument("--decoy_type", default="Unsupported_CD4_state_decoy")
    parser.add_argument("--n_spike_genes", type=int, default=80)
    parser.add_argument("--spike_value", type=float, default=8.0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    root = Path(args.project_root).resolve()
    src = root / "data" / "processed" / args.source_sample
    dst = root / "data" / "processed" / args.target_sample
    src_stage1 = src / "stage1_preprocess"
    if not (src_stage1 / "exported").exists():
        raise FileNotFoundError(src_stage1 / "exported")
    if dst.exists():
        if not args.overwrite:
            raise FileExistsError(f"{dst} exists; pass --overwrite")
        shutil.rmtree(dst)
    (dst).mkdir(parents=True, exist_ok=True)
    shutil.copytree(src_stage1, dst / "stage1_preprocess")

    result_stage3 = root / "result" / args.target_sample / "stage3_typematch"
    if result_stage3.exists():
        shutil.rmtree(result_stage3)

    export = dst / "stage1_preprocess" / "exported"
    meta = pd.read_csv(export / "sc_metadata.csv")
    if "cell_type" not in meta.columns:
        raise ValueError("sc_metadata.csv missing cell_type")
    decoy_mask = meta["cell_type"].astype(str).eq(args.decoy_type)
    n_decoy = int(decoy_mask.sum())
    if n_decoy <= 0:
        raise ValueError(f"No decoy cells found for {args.decoy_type}")

    sc = pd.read_csv(export / "sc_expression_normalized.csv", index_col=0)
    st = pd.read_csv(export / "st_expression_normalized.csv", index_col=0)
    sc.index = sc.index.astype(str)
    sc.columns = sc.columns.astype(str)
    st.columns = st.columns.astype(str)
    decoy_ids = meta.loc[decoy_mask, "cell_id"].astype(str).tolist()
    missing = sorted(set(decoy_ids).difference(sc.index))
    if missing:
        raise ValueError(f"{len(missing)} decoy cells missing from sc expression")
    st_mean = st.mean(axis=0).reindex(sc.columns).fillna(0.0)
    sc_mean = sc.mean(axis=0)
    candidates = pd.DataFrame({"st_mean": st_mean, "sc_mean": sc_mean})
    candidates = candidates[candidates["sc_mean"] > 0.02].sort_values(["st_mean", "sc_mean"], ascending=[True, False])
    if len(candidates) < args.n_spike_genes:
        raise ValueError(f"Only {len(candidates)} spike candidates found; need {args.n_spike_genes}")
    spike_genes = candidates.index[: args.n_spike_genes].astype(str).tolist()
    rng = np.random.default_rng(args.seed)
    noise = rng.normal(0.0, 0.05, size=(len(decoy_ids), len(spike_genes)))
    sc.loc[decoy_ids, spike_genes] = np.clip(args.spike_value + noise, 0.0, None)
    sc.to_csv(export / "sc_expression_normalized.csv", index_label="cell_id")

    summary = {
        "sample": args.target_sample,
        "source_sample": args.source_sample,
        "decoy_type": args.decoy_type,
        "n_decoy_cells": n_decoy,
        "n_spike_genes": int(len(spike_genes)),
        "spike_value": float(args.spike_value),
        "spike_genes": spike_genes,
        "forced_stage3_outputs_written": False,
        "intended_stage3_target": args.decoy_type,
    }
    (dst / "stage1_preprocess" / "fig2k_stage3_decoy_info.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    pd.DataFrame(
        {
            "gene": spike_genes,
            "st_mean": [float(candidates.loc[g, "st_mean"]) for g in spike_genes],
            "sc_mean_before_spike": [float(candidates.loc[g, "sc_mean"]) for g in spike_genes],
        }
    ).to_csv(dst / "stage1_preprocess" / "fig2k_stage3_decoy_spike_genes.csv", index=False)
    _write_config(root, args.target_sample, dst / "stage1_preprocess" / "hvg_genes.txt")
    print(f"[OK] wrote Stage3-discovery Fig.2k sample: {args.target_sample}")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
