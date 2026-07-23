#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse


CATEGORY_ORDER = [
    "Astrocyte",
    "Excitatory neuron",
    "Inhibitory neuron",
    "Neuroblast",
    "Microglia",
    "Endothelial",
    "Oligodendrocyte/OPC",
]

TARGET_TYPES = ["Ext_Thal_1", "Ext_Thal_2"]
TARGET_LABEL = "Thalamic excitatory"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Render single-target Panel B/C for the Stage3B reference-missing stress experiment using the thalamic_top15 case."
    )
    p.add_argument("--project_root", default=".")
    p.add_argument(
        "--case_dir",
        default="visualizations/cell2location_stage3b_case/thalamic_top15",
    )
    p.add_argument(
        "--out_dir",
        default="visualizations/stage3b_reference_missing_stress",
    )
    p.add_argument(
        "--output_prefix",
        default="cell2location_thalamic_top15_reference_missing_panels_bc",
    )
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


def load_forced_assignment(case_dir: Path) -> pd.DataFrame:
    path = case_dir / "cell2location_ST8059051_thalamic_excitatory_forced_assignment_composition.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    composition = pd.read_csv(path)
    required = {"assigned_type", "spots", "fraction"}
    if not required.issubset(composition.columns):
        raise ValueError(f"Unexpected forced-assignment columns: {composition.columns.tolist()}")
    composition = composition.copy()
    composition["assigned_type"] = composition["assigned_type"].astype(str)
    return composition.sort_values("fraction", ascending=False).reset_index(drop=True)


def mean_profile(x: sparse.spmatrix, mask: np.ndarray) -> np.ndarray:
    if int(mask.sum()) == 0:
        return np.full(x.shape[1], np.nan, dtype=float)
    return np.asarray(x[mask, :].mean(axis=0)).ravel().astype(float)


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    ok = np.isfinite(a) & np.isfinite(b)
    if int(ok.sum()) == 0:
        return np.nan
    aa = a[ok]
    bb = b[ok]
    denom = np.linalg.norm(aa) * np.linalg.norm(bb)
    return float(aa @ bb / denom) if denom > 0 else np.nan


def compute_expression_similarity(root: Path, assigned_types: list[str]) -> pd.DataFrame:
    sc_path = root / "data" / "raw" / "cell2location_mouse_brain" / "sc.h5ad"
    sc = ad.read_h5ad(sc_path)
    if sc.raw is None:
        raise ValueError("sc.h5ad lacks raw counts")

    labels = sc.obs["annotation_1"].astype(str)
    raw = sc.raw.X
    if not sparse.issparse(raw):
        raw = sparse.csr_matrix(raw)
    raw = raw.astype(np.float32).tocsr()

    totals = np.asarray(raw.sum(axis=1)).ravel()
    totals[totals <= 0] = 1.0
    norm = raw.multiply(10000.0 / totals[:, None]).log1p().tocsr()

    target_mask = labels.isin(TARGET_TYPES).to_numpy()
    target_profile = mean_profile(norm, target_mask)

    rows = []
    for category in assigned_types:
        category_mask = labels.map(broad_cell_type).eq(category).to_numpy()
        category_mask &= ~target_mask
        rows.append(
            {
                "removed_st_only_type": TARGET_LABEL,
                "remaining_sc_type": category,
                "cosine_similarity": cosine(target_profile, mean_profile(norm, category_mask)),
                "reference_cells": int(category_mask.sum()),
            }
        )
    long = pd.DataFrame(rows)
    return long.pivot(index="removed_st_only_type", columns="remaining_sc_type", values="cosine_similarity")


def plot_panels(misassignment: pd.DataFrame, similarity: pd.DataFrame, out_png: Path, out_pdf: Path) -> None:
    sns.set_theme(style="white")
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.edgecolor": "#333333",
            "axes.linewidth": 1.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig = plt.figure(figsize=(11.0, 3.4), dpi=260, constrained_layout=True)
    grid = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.0])
    ax_b = fig.add_subplot(grid[0, 0])
    ax_c = fig.add_subplot(grid[0, 1])

    sns.heatmap(
        misassignment,
        ax=ax_b,
        cmap="YlOrRd",
        vmin=0,
        vmax=max(float(np.nanmax(misassignment.to_numpy())), 0.01),
        mask=misassignment.isna(),
        linewidths=0.5,
        linecolor="#f1f1f1",
        annot=True,
        fmt=".2f",
        cbar_kws={"label": "Assignment fraction"},
    )
    ax_b.set_title("B. CytoSPACE forced assignment of ST-only region", weight="bold", fontsize=10.5)
    ax_b.set_xlabel("Assigned SC type")
    ax_b.set_ylabel("Removed ST-only type")
    ax_b.set_xticklabels(ax_b.get_xticklabels(), rotation=45, ha="right")
    ax_b.set_yticklabels(ax_b.get_yticklabels(), rotation=0)

    sns.heatmap(
        similarity,
        ax=ax_c,
        cmap="viridis",
        vmin=0,
        vmax=1,
        mask=similarity.isna(),
        linewidths=0.5,
        linecolor="#f1f1f1",
        annot=True,
        fmt=".2f",
        cbar_kws={"label": "Cosine similarity"},
    )
    ax_c.set_title("C. Expression similarity to remaining SC types", weight="bold", fontsize=10.5)
    ax_c.set_xlabel("Remaining SC type")
    ax_c.set_ylabel("Removed ST-only type")
    ax_c.set_xticklabels(ax_c.get_xticklabels(), rotation=45, ha="right")
    ax_c.set_yticklabels(ax_c.get_yticklabels(), rotation=0)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, bbox_inches="tight", facecolor="white")
    fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    root = Path(args.project_root).resolve()
    case_dir = root / args.case_dir
    out_dir = root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    composition = load_forced_assignment(case_dir)
    assigned_order = [category for category in CATEGORY_ORDER if category in set(composition["assigned_type"])]
    extra = [category for category in composition["assigned_type"] if category not in assigned_order]
    assigned_order.extend(extra)

    misassignment = composition.set_index("assigned_type")["fraction"].reindex(assigned_order).to_frame().T
    misassignment.index = [TARGET_LABEL]
    similarity = compute_expression_similarity(root, assigned_order).reindex(columns=assigned_order)

    prefix = args.output_prefix
    misassignment.to_csv(out_dir / f"{prefix}_cytospace_forced_assignment.csv")
    similarity.to_csv(out_dir / f"{prefix}_expression_similarity.csv")
    composition.to_csv(out_dir / f"{prefix}_forced_assignment_composition_source.csv", index=False)

    plot_panels(
        misassignment,
        similarity,
        out_dir / f"{prefix}.png",
        out_dir / f"{prefix}.pdf",
    )
    print(f"[done] {out_dir / f'{prefix}.png'}")
    print(f"[done] {out_dir / f'{prefix}.pdf'}")
    print(misassignment.to_string())
    print(similarity.to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
