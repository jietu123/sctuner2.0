from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from scripts.plot_cytospace_fig2k_mapping_with_route2 import STATE_MARKERS, _load_table_s9


DEFAULT_SAMPLE = "cytospace_fig2k_breast_state_decoy"
MERSCOPE_DIR = "data/raw/cytospace_fig2k_breast_merscope/HumanBreastCancerPatient1"
MERSCOPE_PREP = "data/processed/cytospace_fig2k_breast_merscope/HumanBreastCancerPatient1/prepared"
SUPP_XLSX = "data/raw/cytospace_fig2c_melanoma/41587_2023_1697_MOESM3_ESM.xlsx"


def _load_merscope_subset(expr_path: Path, cell_ids: list[str], genes: list[str], chunk_size: int) -> pd.DataFrame:
    wanted = set(map(str, cell_ids))
    chunks = []
    for chunk in pd.read_csv(expr_path, usecols=["cell", *genes], chunksize=chunk_size):
        chunk["cell_id"] = chunk["cell"].astype(str)
        sub = chunk[chunk["cell_id"].isin(wanted)].copy()
        if len(sub):
            chunks.append(sub.drop(columns=["cell"]).set_index("cell_id"))
    if not chunks:
        raise ValueError("No selected MERSCOPE cells found in expression matrix.")
    out = pd.concat(chunks, axis=0)
    return out.loc[[cid for cid in cell_ids if cid in out.index], genes].astype("float32")


def _write_config(root: Path, sample: str) -> None:
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
            "auto_missing_detection": {"enable": False, "action": "mark_unknown"},
            "masked_missing_detection": {"enable": False},
        },
    }
    path = root / "configs" / "datasets" / f"{sample}.yaml"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")


def _apply_state_signal(
    base: pd.DataFrame,
    zones: pd.Series,
    table_s9: pd.DataFrame,
    amplitude: float,
    invert: bool,
) -> pd.DataFrame:
    out = base.copy()
    ranks = table_s9.set_index("state")["known_rank"].astype(float)
    rank_z = (ranks - ranks.mean()) / ranks.std(ddof=0)
    for state, marker_genes in STATE_MARKERS.items():
        genes = [g for g in marker_genes if g in out.columns]
        if not genes or state not in rank_z.index:
            continue
        state_signal = float(rank_z.loc[state])
        if invert:
            state_signal *= -1.0
        zone_sign = zones.map({"tumor_high": 1.0, "tumor_low": -1.0}).fillna(0.0).to_numpy(dtype=float)
        delta = amplitude * state_signal * zone_sign
        out.loc[:, genes] = np.maximum(0.0, out[genes].to_numpy(dtype=float) + delta[:, None])
    return out.astype("float32")


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare a Fig.2k CD4 state-level unsupported decoy benchmark.")
    parser.add_argument("--project_root", default=".")
    parser.add_argument("--sample", default=DEFAULT_SAMPLE)
    parser.add_argument("--n_spots", type=int, default=1200)
    parser.add_argument("--state_signal_amplitude", type=float, default=1.8)
    parser.add_argument("--noise_sd", type=float, default=0.02)
    parser.add_argument("--decoy_copies_per_spot", type=int, default=1)
    parser.add_argument("--chunk_size", type=int, default=50000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    root = Path(args.project_root).resolve()
    rng = np.random.default_rng(args.seed)
    sample = args.sample
    out_root = root / "data" / "processed" / sample
    export = out_root / "stage1_preprocess" / "exported"
    stage3 = out_root / "stage3_typematch"
    result_stage3 = root / "result" / sample / "stage3_typematch"
    if out_root.exists() and args.overwrite:
        shutil.rmtree(out_root)
    if result_stage3.exists() and args.overwrite:
        shutil.rmtree(result_stage3)
    export.mkdir(parents=True, exist_ok=True)
    stage3.mkdir(parents=True, exist_ok=True)
    result_stage3.mkdir(parents=True, exist_ok=True)

    merscope = root / MERSCOPE_DIR
    merscope_prep = root / MERSCOPE_PREP
    table_s9 = _load_table_s9(root / SUPP_XLSX)
    header = pd.read_csv(merscope / "cell_by_gene.csv", nrows=0).columns.astype(str).tolist()
    genes = [g for g in header if g != "cell" and not g.startswith("Blank-")]

    candidates = pd.read_csv(merscope_prep / "merscope_tcell_candidate_scores.csv")
    candidates["cell_id"] = candidates["cell_id"].astype(str)
    cd4 = candidates[candidates["is_cd4_t_candidate"].astype(bool)].copy()
    high_thr = cd4["tumor_proximity_score"].quantile(0.75)
    low_thr = cd4["tumor_proximity_score"].quantile(0.25)
    high_pool = cd4[cd4["tumor_proximity_score"].ge(high_thr)]["cell_id"].astype(str).tolist()
    low_pool = cd4[cd4["tumor_proximity_score"].le(low_thr)]["cell_id"].astype(str).tolist()
    n_high = args.n_spots // 2
    selected = (
        rng.choice(high_pool, size=min(n_high, len(high_pool)), replace=False).astype(str).tolist()
        + rng.choice(low_pool, size=min(args.n_spots - n_high, len(low_pool)), replace=False).astype(str).tolist()
    )
    rng.shuffle(selected)
    spot_meta = candidates.set_index("cell_id").loc[selected].copy()
    zones = pd.Series(
        np.where(
            spot_meta["tumor_proximity_score"].ge(high_thr),
            "tumor_high",
            np.where(spot_meta["tumor_proximity_score"].le(low_thr), "tumor_low", "middle"),
        ),
        index=selected,
        name="tumor_zone",
    )

    base = np.log1p(_load_merscope_subset(merscope / "cell_by_gene.csv", selected, genes, args.chunk_size))
    st_expr = _apply_state_signal(base, zones, table_s9, args.state_signal_amplitude, invert=False)
    st_expr.index.name = "spot_id"
    st_expr.to_csv(export / "st_expression_normalized.csv")

    true_expr = st_expr.copy()
    true_expr.index = [f"true_{idx}" for idx in st_expr.index]
    decoy_base = _apply_state_signal(base, zones, table_s9, args.state_signal_amplitude, invert=True)
    decoy_exprs = []
    for copy_idx in range(max(1, args.decoy_copies_per_spot)):
        copied = decoy_base.copy()
        copied.index = [f"state_decoy{copy_idx}_{idx}" for idx in copied.index]
        decoy_exprs.append(copied)
    decoy_expr = pd.concat(decoy_exprs, axis=0)
    if args.noise_sd > 0:
        true_expr = np.maximum(0, true_expr + rng.normal(0, args.noise_sd, true_expr.shape)).astype("float32")
        decoy_expr = np.maximum(0, decoy_expr + rng.normal(0, args.noise_sd, decoy_expr.shape)).astype("float32")
    sc_expr = pd.concat([true_expr, decoy_expr], axis=0)
    sc_expr.index.name = "cell_id"
    sc_expr.to_csv(export / "sc_expression_normalized.csv")

    sc_meta = pd.DataFrame(
        {
            "cell_id": sc_expr.index.astype(str),
            "cell_type": ["CD4 T cells"] * len(true_expr) + ["Unsupported_CD4_state_decoy"] * len(decoy_expr),
        }
    )
    sc_meta.to_csv(export / "sc_metadata.csv", index=False)
    coords = pd.DataFrame(
        {
            "spot_id": selected,
            "row": spot_meta["center_y"].astype(float).to_numpy(),
            "col": spot_meta["center_x"].astype(float).to_numpy(),
        }
    ).set_index("spot_id")
    coords.to_csv(export / "st_coordinates.csv")
    spot_meta = spot_meta.copy()
    spot_meta["tumor_zone"] = zones
    spot_meta.reset_index(names="spot_id").to_csv(export / "merscope_spot_metadata.csv", index=False)

    pd.Series(genes).to_csv(out_root / "stage1_preprocess" / "hvg_genes.txt", index=False, header=False)
    pd.Series(genes).to_csv(out_root / "stage1_preprocess" / "common_genes.txt", index=False, header=False)
    pd.Series(genes).to_csv(export / "selected_genes.txt", index=False, header=False)

    relabel = pd.DataFrame(
        {
            "cell_id": sc_expr.index.astype(str),
            "orig_type": sc_meta["cell_type"],
            "plugin_type": sc_meta["cell_type"],
            "label_status": "kept",
        }
    )
    decoy_mask = relabel["orig_type"].eq("Unsupported_CD4_state_decoy")
    relabel.loc[decoy_mask, "plugin_type"] = "Unknown_sc_only"
    relabel.loc[decoy_mask, "label_status"] = "forced_unsupported_state_decoy"
    relabel.to_csv(stage3 / "cell_type_relabel.csv", index=False)
    pd.DataFrame(
        [
            {"orig_type": "CD4 T cells", "n_cells": int(len(true_expr)), "support_category": "supported"},
            {"orig_type": "Unsupported_CD4_state_decoy", "n_cells": int(len(decoy_expr)), "support_category": "unsupported"},
        ]
    ).to_csv(stage3 / "type_support.csv", index=False)

    summary = {
        "sample": sample,
        "n_spots": int(len(st_expr)),
        "n_true_cd4_cells": int(len(true_expr)),
        "n_state_decoy_cells": int(len(decoy_expr)),
        "decoy_copies_per_spot": int(args.decoy_copies_per_spot),
        "n_genes": int(len(genes)),
        "state_signal_amplitude": float(args.state_signal_amplitude),
        "noise_sd": float(args.noise_sd),
        "forced_unsupported_type": "Unsupported_CD4_state_decoy",
        "unknown_overview": {"n_unknown": int(decoy_mask.sum()), "unknown_label_effective": "Unknown_sc_only"},
    }
    for path in [stage3 / "stage3_summary.json", result_stage3 / "stage3_summary.json"]:
        path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    _write_config(root, sample)
    print(f"[OK] wrote sample: {sample}")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
