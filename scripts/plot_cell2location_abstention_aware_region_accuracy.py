#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse

_HERE = Path(__file__).resolve()
_ROOT = _HERE.parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from scripts.prepare_cell2location_mouse_brain_stage3b_case import (  # noqa: E402
    read_10x_h5,
    target_markers_from_counts,
    unique_gene_positions,
)


CATEGORY_ORDER = [
    "Astrocyte",
    "Excitatory neuron",
    "Inhibitory neuron",
    "Neuroblast",
    "Microglia",
    "Endothelial",
    "Oligodendrocyte/OPC",
    "Other",
]

SUPPORTED_REGIONS = {
    "Astrocyte": [
        "Astro_AMY",
        "Astro_AMY_CTX",
        "Astro_CTX",
        "Astro_HPC",
        "Astro_HYPO",
        "Astro_STR",
        "Astro_THAL_hab",
        "Astro_THAL_lat",
        "Astro_THAL_med",
        "Astro_WM",
    ],
    "Inhibitory neuron": [
        "Inh_1",
        "Inh_2",
        "Inh_3",
        "Inh_4",
        "Inh_5",
        "Inh_6",
        "Inh_Lamp5",
        "Inh_Meis2_1",
        "Inh_Meis2_2",
        "Inh_Meis2_3",
        "Inh_Meis2_4",
        "Inh_Pvalb",
        "Inh_Sst",
        "Inh_Vip",
    ],
    "Oligodendrocyte/OPC": ["Oligo_1", "Oligo_2", "OPC_1", "OPC_2"],
    "Microglia": ["Micro"],
}

PALETTE = {
    "unsupported": "#1b9e77",
    "supported": "#4c78a8",
    "incorrect": "#d9d9d9",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Plot one-panel abstention-aware region accuracy for cell2location Stage3B case."
    )
    p.add_argument("--project_root", default=".")
    p.add_argument("--sample", default="cell2loc_scan_st8059051_sc_missing_thalamic_excitatory")
    p.add_argument("--section", default="ST8059051")
    p.add_argument("--target", default="thalamic_excitatory")
    p.add_argument("--target_label", default="Thalamic excitatory")
    p.add_argument(
        "--target_region_csv",
        default=(
            "visualizations/cell2location_stage3b_case/thalamic_top15/"
            "cell2location_ST8059051_thalamic_excitatory_marker_region.csv"
        ),
    )
    p.add_argument("--supported_quantile", type=float, default=0.85)
    p.add_argument("--n_markers", type=int, default=30)
    p.add_argument("--out_dir", default="visualizations/cell2location_stage3b_case/thalamic_top15")
    return p.parse_args()


def broad_cell_type(cell_type: str) -> str:
    cell_type = str(cell_type)
    if cell_type.startswith("Astro"):
        return "Astrocyte"
    if cell_type.startswith("Ext_"):
        return "Excitatory neuron"
    if cell_type.startswith("Inh_"):
        return "Inhibitory neuron"
    if cell_type.startswith("Nb_"):
        return "Neuroblast"
    if cell_type.startswith("Micro"):
        return "Microglia"
    if cell_type.startswith("Endo"):
        return "Endothelial"
    if cell_type.startswith("Oligo") or cell_type.startswith("OPC"):
        return "Oligodendrocyte/OPC"
    return "Other"


def aggregate_to_broad_categories(frac: pd.DataFrame) -> pd.DataFrame:
    groups: dict[str, list[str]] = {}
    for col in frac.columns.astype(str):
        groups.setdefault(broad_cell_type(col), []).append(col)
    out = pd.DataFrame(index=frac.index)
    for category, cols in groups.items():
        out[category] = frac[cols].sum(axis=1)
    ordered = [category for category in CATEGORY_ORDER if category in out.columns]
    return out.loc[:, ordered]


def read_fractional(root: Path, sample: str, suffix: str) -> pd.DataFrame:
    path = (
        root
        / "result"
        / sample
        / f"stage4_cytospace{suffix}"
        / "cytospace_output"
        / "fractional_abundances_by_spot.csv"
    )
    if not path.exists():
        raise FileNotFoundError(f"Missing fractional abundance file: {path}")
    frame = pd.read_csv(path, index_col=0)
    frame.index = frame.index.astype(str)
    return frame.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def log_normalize_rows(values: sparse.spmatrix) -> sparse.csr_matrix:
    values = values.astype(np.float32).tocsr()
    totals = np.asarray(values.sum(axis=1)).ravel()
    totals[totals <= 0] = 1.0
    return values.multiply(10000.0 / totals[:, None]).log1p().tocsr()


def load_st_norm(root: Path, section: str) -> tuple[sparse.csr_matrix, list[str], pd.DataFrame]:
    source = root / "data" / "raw" / "cell2location_mouse_brain"
    section_dir = source / "mouse_brain_visium_wo_cloupe_data" / "rawdata" / section
    counts, barcodes, st_genes = read_10x_h5(section_dir / "filtered_feature_bc_matrix.h5")
    positions = pd.read_csv(
        section_dir / "spatial" / "tissue_positions_list.csv",
        header=None,
        names=["spot_id", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"],
    )
    positions["spot_id"] = positions["spot_id"].astype(str)
    barcode_pos = {barcode: i for i, barcode in enumerate(barcodes.astype(str))}
    tissue_spots = [
        spot for spot in positions.loc[positions["in_tissue"] == 1, "spot_id"].tolist() if spot in barcode_pos
    ]
    counts = counts[np.array([barcode_pos[spot] for spot in tissue_spots], dtype=int), :].tocsr()
    return log_normalize_rows(counts), tissue_spots, pd.DataFrame({"gene": st_genes.astype(str)})


def supported_region_table(
    root: Path,
    section: str,
    unsupported_core: pd.Series,
    supported_quantile: float,
    n_markers: int,
) -> pd.DataFrame:
    source = root / "data" / "raw" / "cell2location_mouse_brain"
    sc = ad.read_h5ad(source / "sc.h5ad")
    if sc.raw is None:
        raise ValueError("sc.h5ad lacks raw counts")
    labels = sc.obs["annotation_1"].astype(str)
    raw_counts = sc.raw.X
    if not sparse.issparse(raw_counts):
        raw_counts = sparse.csr_matrix(raw_counts)

    st_norm, tissue_spots, st_gene_frame = load_st_norm(root, section)
    st_genes = st_gene_frame["gene"].to_numpy(dtype=str)
    st_pos = unique_gene_positions(st_genes)
    out = pd.DataFrame(index=pd.Index(tissue_spots, name="spot_id"))
    out["unsupported_core"] = out.index.to_series().isin(unsupported_core.index[unsupported_core])
    out["region_kind"] = "unused"
    out["expected_behavior"] = "unused"
    out["expected_type"] = ""
    out["region_score"] = np.nan
    out.loc[out["unsupported_core"], "region_kind"] = "unsupported_core"
    out.loc[out["unsupported_core"], "expected_behavior"] = "blank"
    out.loc[out["unsupported_core"], "expected_type"] = "Blank"

    candidate_scores: dict[str, np.ndarray] = {}
    candidate_masks: dict[str, np.ndarray] = {}
    for category, drop_types in SUPPORTED_REGIONS.items():
        markers = target_markers_from_counts(raw_counts, sc.raw.var, labels, drop_types, st_pos)
        markers = [gene for gene in markers[:n_markers] if gene in st_pos]
        if len(markers) < 5:
            continue
        idx = np.array([st_pos[gene] for gene in markers], dtype=int)
        score = np.asarray(st_norm[:, idx].mean(axis=1)).ravel()
        cutoff = float(np.quantile(score, supported_quantile))
        candidate_scores[category] = score
        candidate_masks[category] = score >= cutoff

    unsupported = out["unsupported_core"].to_numpy(dtype=bool)
    available = ~unsupported
    score_mat = pd.DataFrame(candidate_scores, index=out.index)
    mask_mat = pd.DataFrame(candidate_masks, index=out.index)
    any_supported = mask_mat.any(axis=1).to_numpy(dtype=bool) & available
    if any_supported.any():
        masked_scores = score_mat.where(mask_mat, -np.inf)
        winner = masked_scores.idxmax(axis=1)
        out.loc[any_supported, "region_kind"] = "supported_marker_region"
        out.loc[any_supported, "expected_behavior"] = "map"
        out.loc[any_supported, "expected_type"] = winner.loc[any_supported]
        out.loc[any_supported, "region_score"] = masked_scores.max(axis=1).loc[any_supported]
    return out


def method_region_scores(region: pd.DataFrame, frac: pd.DataFrame, method: str) -> pd.DataFrame:
    common = region.index.intersection(frac.index)
    region = region.loc[common].copy()
    frac = frac.loc[common].copy()
    broad = aggregate_to_broad_categories(frac)
    row_sum = broad.sum(axis=1)
    mapped = row_sum > 0
    dominant = broad.idxmax(axis=1).where(mapped, "Blank")
    region["method"] = method
    region["mapped"] = mapped
    region["predicted_type"] = dominant
    region["correct"] = False
    unsupported = region["region_kind"].eq("unsupported_core")
    supported = region["region_kind"].eq("supported_marker_region")
    region.loc[unsupported, "correct"] = ~mapped.loc[unsupported]
    region.loc[supported, "correct"] = (
        mapped.loc[supported]
        & dominant.loc[supported].eq(region.loc[supported, "expected_type"])
    )
    return region


def main() -> int:
    args = parse_args()
    root = Path(args.project_root).resolve()
    out_dir = root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    target_path = root / args.target_region_csv
    target = pd.read_csv(target_path, index_col=0)
    target.index = target.index.astype(str)
    target["target_region"] = target["target_region"].astype(bool)
    unsupported_core = target["target_region"]

    region = supported_region_table(
        root,
        args.section,
        unsupported_core,
        args.supported_quantile,
        args.n_markers,
    )
    eval_region = region[region["region_kind"].isin(["unsupported_core", "supported_marker_region"])].copy()
    if eval_region.empty:
        raise ValueError("No evaluation regions were defined")

    baseline = read_fractional(root, args.sample, "_baseline")
    svtuner = read_fractional(root, args.sample, "_stage3b_blank")
    scored = pd.concat(
        [
            method_region_scores(eval_region, baseline, "CytoSPACE"),
            method_region_scores(eval_region, svtuner, "SVTuner"),
        ],
        axis=0,
    )

    total_spots = int(eval_region.shape[0])
    rows = []
    for method, sub in scored.groupby("method"):
        unsupported = sub["region_kind"].eq("unsupported_core")
        supported = sub["region_kind"].eq("supported_marker_region")
        rows.append(
            {
                "method": method,
                "evaluated_spots": total_spots,
                "unsupported_core_spots": int(unsupported.sum()),
                "supported_region_spots": int(supported.sum()),
                "unsupported_correct": int((unsupported & sub["correct"]).sum()),
                "supported_correct": int((supported & sub["correct"]).sum()),
                "unsupported_contribution": float((unsupported & sub["correct"]).sum() / total_spots),
                "supported_contribution": float((supported & sub["correct"]).sum() / total_spots),
                "abstention_aware_accuracy": float(sub["correct"].sum() / total_spots),
                "incorrect_fraction": float((~sub["correct"]).sum() / total_spots),
            }
        )
    metrics = pd.DataFrame(rows)
    method_order = ["CytoSPACE", "SVTuner"]
    metrics["method"] = pd.Categorical(metrics["method"], categories=method_order, ordered=True)
    metrics = metrics.sort_values("method").reset_index(drop=True)

    prefix = f"cell2location_{args.section}_{args.target}_abstention_aware_region_accuracy"
    scored.to_csv(out_dir / f"{prefix}_spots.csv", index_label="spot_id")
    eval_region.to_csv(out_dir / f"{prefix}_regions.csv", index_label="spot_id")
    metrics.to_csv(out_dir / f"{prefix}.csv", index=False)

    fig, ax = plt.subplots(figsize=(4.9, 4.6), dpi=300)
    x = np.arange(len(metrics))
    unsupported_vals = metrics["unsupported_contribution"].to_numpy()
    supported_vals = metrics["supported_contribution"].to_numpy()
    incorrect_vals = metrics["incorrect_fraction"].to_numpy()
    ax.bar(x, unsupported_vals, color=PALETTE["unsupported"], width=0.58, label="Correct blanking in unsupported core")
    ax.bar(
        x,
        supported_vals,
        bottom=unsupported_vals,
        color=PALETTE["supported"],
        width=0.58,
        label="Correct mapping in supported regions",
    )
    for i, row in metrics.iterrows():
        ax.text(
            i,
            min(1.02, row["abstention_aware_accuracy"] + 0.025),
            f"{row['abstention_aware_accuracy']:.2f}",
            ha="center",
            va="bottom",
            fontsize=10,
            weight="bold",
        )
    ax.set_ylim(0, 1.08)
    ax.set_xticks(x)
    ax.set_xticklabels(metrics["method"].astype(str), fontsize=9)
    ax.set_ylabel("Abstention-aware accuracy", fontsize=9)
    ax.set_title("Mapping accuracy with correct blanking credited", fontsize=10.5, weight="bold", pad=8)
    ax.grid(axis="y", color="#d9d9d9", linewidth=0.7)
    ax.set_axisbelow(True)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.legend(frameon=False, fontsize=7.2, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=1)
    png = out_dir / f"{prefix}.png"
    pdf = out_dir / f"{prefix}.pdf"
    fig.savefig(png, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)

    supported_counts = (
        eval_region.loc[eval_region["region_kind"].eq("supported_marker_region"), "expected_type"]
        .value_counts()
        .rename_axis("expected_type")
        .reset_index(name="spots")
    )
    summary = {
        "sample": args.sample,
        "section": args.section,
        "target_label": args.target_label,
        "target_region_definition": str(target_path),
        "unsupported_core_spots": int(eval_region["region_kind"].eq("unsupported_core").sum()),
        "supported_region_spots": int(eval_region["region_kind"].eq("supported_marker_region").sum()),
        "supported_quantile": args.supported_quantile,
        "supported_region_counts": supported_counts.to_dict(orient="records"),
        "metrics": metrics.to_dict(orient="records"),
        "png": str(png),
        "pdf": str(pdf),
    }
    (out_dir / f"{prefix}_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
