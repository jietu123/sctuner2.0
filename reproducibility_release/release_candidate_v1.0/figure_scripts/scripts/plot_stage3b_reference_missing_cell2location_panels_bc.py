#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
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
]

SHORT_LABELS = {
    "Astrocyte": "Astro.",
    "Excitatory neuron": "Excit.",
    "Inhibitory neuron": "Inhib.",
    "Neuroblast": "Neurob.",
    "Microglia": "Micro.",
    "Endothelial": "Endoth.",
    "Oligodendrocyte/OPC": "Oligo/OPC",
}

CASES = [
    {
        "label": "Astrocyte",
        "section": "ST8059051",
        "sample": "cell2loc_scan_st8059051_sc_missing_astrocyte",
        "drop_types": [
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
    },
    {
        "label": "Excitatory neuron",
        "section": "ST8059051",
        "sample": "cell2loc_scan_st8059051_sc_missing_excitatory_neuron",
        "drop_types": [
            "Ext_Amy_1",
            "Ext_Amy_2",
            "Ext_ClauPyr",
            "Ext_Hpc_CA1",
            "Ext_Hpc_CA2",
            "Ext_Hpc_CA3",
            "Ext_Hpc_DG1",
            "Ext_Hpc_DG2",
            "Ext_L23",
            "Ext_L25",
            "Ext_L56",
            "Ext_L5_1",
            "Ext_L5_2",
            "Ext_L5_3",
            "Ext_L6",
            "Ext_L6B",
            "Ext_Med",
            "Ext_Pir",
            "Ext_Thal_1",
            "Ext_Thal_2",
            "Ext_Unk_1",
            "Ext_Unk_2",
            "Ext_Unk_3",
        ],
    },
    {
        "label": "Inhibitory neuron",
        "section": "ST8059051",
        "sample": "cell2loc_scan_st8059051_sc_missing_inhibitory_neuron",
        "drop_types": [
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
    },
    {
        "label": "Neuroblast",
        "section": "ST8059051",
        "sample": "cell2loc_scan_st8059051_sc_missing_neuroblast",
        "drop_types": ["Nb_1", "Nb_2"],
    },
    {
        "label": "Microglia",
        "section": "ST8059051",
        "sample": "cell2loc_scan_st8059051_sc_missing_microglia",
        "drop_types": ["Micro"],
    },
    {
        "label": "Endothelial",
        "section": "ST8059051",
        "sample": "cell2loc_scan_st8059051_sc_missing_endothelial",
        "drop_types": ["Endo"],
    },
    {
        "label": "Oligodendrocyte/OPC",
        "section": "ST8059051",
        "sample": "cell2loc_scan_st8059051_sc_missing_oligodendrocyte_opc",
        "drop_types": ["Oligo_1", "Oligo_2", "OPC_1", "OPC_2"],
    },
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Render updated multi-target Panel B/C for the Stage3B reference-missing stress experiment using available cell2location cases."
    )
    p.add_argument("--project_root", default=".")
    p.add_argument("--target_quantile", type=float, default=0.85)
    p.add_argument("--n_markers", type=int, default=30)
    p.add_argument(
        "--out_dir",
        default="visualizations/stage3b_reference_missing_stress",
    )
    p.add_argument(
        "--output_prefix",
        default="cell2location_reference_missing_multitarget_panels_bc",
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


def aggregate_to_broad_categories(frac: pd.DataFrame) -> pd.DataFrame:
    groups: dict[str, list[str]] = {}
    for col in frac.columns.astype(str):
        category = broad_cell_type(col)
        if category in CATEGORY_ORDER:
            groups.setdefault(category, []).append(col)
    out = pd.DataFrame(index=frac.index)
    for category, cols in groups.items():
        out[category] = frac[cols].sum(axis=1)
    return out.reindex(columns=CATEGORY_ORDER, fill_value=0.0)


def mapping_path(root: Path, sample: str) -> Path:
    candidates = [
        root / "result" / sample / "stage4_cytospace_baseline" / "cytospace_output" / "fractional_abundances_by_spot.csv",
        root / "result" / "cell2location_mouse_brain" / sample / "stage4_cytospace_baseline" / "cytospace_output" / "fractional_abundances_by_spot.csv",
    ]
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError(f"No baseline fractional abundance file found for {sample}")


def log_normalize_rows(values: sparse.spmatrix) -> sparse.csr_matrix:
    values = values.astype(np.float32).tocsr()
    totals = np.asarray(values.sum(axis=1)).ravel()
    totals[totals <= 0] = 1.0
    return values.multiply(10000.0 / totals[:, None]).log1p().tocsr()


def load_st_norm(root: Path, section: str) -> tuple[sparse.csr_matrix, list[str], np.ndarray]:
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
    row_idx = np.array([barcode_pos[spot] for spot in tissue_spots], dtype=int)
    return log_normalize_rows(counts[row_idx, :].tocsr()), tissue_spots, st_genes.astype(str)


def target_region_from_marker(
    root: Path,
    sc: ad.AnnData,
    raw_counts: sparse.spmatrix,
    labels: pd.Series,
    case: dict[str, object],
    target_quantile: float,
    n_markers: int,
) -> pd.Series:
    st_norm, tissue_spots, st_genes = load_st_norm(root, str(case["section"]))
    st_pos = unique_gene_positions(st_genes)
    markers = target_markers_from_counts(
        raw_counts,
        sc.raw.var,
        labels,
        list(case["drop_types"]),
        st_pos,
    )
    markers = [gene for gene in markers[:n_markers] if gene in st_pos]
    if len(markers) < 5:
        raise ValueError(f"Too few marker genes for {case['label']}: {markers}")
    marker_idx = np.array([st_pos[gene] for gene in markers], dtype=int)
    marker_score = np.asarray(st_norm[:, marker_idx].mean(axis=1)).ravel()
    cutoff = float(np.quantile(marker_score, target_quantile))
    return pd.Series(marker_score >= cutoff, index=pd.Index(tissue_spots, name="spot_id"))


def forced_assignment_row(root: Path, case: dict[str, object], target_mask: pd.Series) -> pd.Series:
    frac = pd.read_csv(mapping_path(root, str(case["sample"])), index_col=0)
    frac.index = frac.index.astype(str)
    frac = frac.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    broad = aggregate_to_broad_categories(frac)
    common = target_mask.index.intersection(broad.index)
    region = target_mask.loc[common].astype(bool)
    if int(region.sum()) == 0:
        raise ValueError(f"Empty target marker region for {case['label']}")
    dominant = broad.loc[common].idxmax(axis=1)
    counts = dominant.loc[region].value_counts()
    row = pd.Series(0.0, index=CATEGORY_ORDER, name=str(case["label"]))
    total = float(region.sum())
    for category, count in counts.items():
        if category in row.index:
            row.loc[category] = float(count) / total
    return row


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


def similarity_row(
    norm: sparse.spmatrix,
    labels: pd.Series,
    case: dict[str, object],
) -> pd.Series:
    target_mask = labels.isin(list(case["drop_types"])).to_numpy()
    target_profile = mean_profile(norm, target_mask)
    row = pd.Series(np.nan, index=CATEGORY_ORDER, name=str(case["label"]))
    for category in CATEGORY_ORDER:
        category_mask = labels.map(broad_cell_type).eq(category).to_numpy()
        category_mask &= ~target_mask
        if int(category_mask.sum()) == 0:
            continue
        row.loc[category] = cosine(target_profile, mean_profile(norm, category_mask))
    return row


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

    plot_misassignment = misassignment.copy()
    for label in plot_misassignment.index.intersection(plot_misassignment.columns):
        plot_misassignment.loc[label, label] = np.nan

    short_rows = [SHORT_LABELS.get(x, x) for x in misassignment.index]
    short_cols = [SHORT_LABELS.get(x, x) for x in misassignment.columns]
    plot_misassignment.index = short_rows
    plot_misassignment.columns = short_cols

    plot_similarity = similarity.copy()
    plot_similarity.index = [SHORT_LABELS.get(x, x) for x in similarity.index]
    plot_similarity.columns = [SHORT_LABELS.get(x, x) for x in similarity.columns]

    b_values = plot_misassignment.to_numpy(dtype=float)
    b_vmax = float(np.nanquantile(b_values, 0.95))
    b_vmax = max(b_vmax, 0.05)

    c_values = plot_similarity.to_numpy(dtype=float)
    c_vmin = max(0.45, float(np.nanmin(c_values)) - 0.02)

    fig = plt.figure(figsize=(11.6, 5.0), dpi=300, constrained_layout=True)
    grid = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.16)
    ax_b = fig.add_subplot(grid[0, 0])
    ax_c = fig.add_subplot(grid[0, 1])

    sns.heatmap(
        plot_misassignment,
        ax=ax_b,
        cmap="YlOrRd",
        vmin=0,
        vmax=b_vmax,
        mask=plot_misassignment.isna(),
        linewidths=0.75,
        linecolor="white",
        square=True,
        cbar_kws={"label": "Assignment fraction", "shrink": 0.72, "pad": 0.035, "fraction": 0.045},
    )
    ax_b.set_title("B. Forced assignment of ST-only types", weight="bold", fontsize=11.0, pad=9)
    ax_b.set_xlabel("Assigned SC type")
    ax_b.set_ylabel("Removed ST-only type")
    ax_b.set_xticklabels(ax_b.get_xticklabels(), rotation=45, ha="right", fontsize=9.2)
    ax_b.set_yticklabels(ax_b.get_yticklabels(), rotation=0, fontsize=9.2)

    sns.heatmap(
        plot_similarity,
        ax=ax_c,
        cmap="viridis",
        vmin=c_vmin,
        vmax=1,
        mask=plot_similarity.isna(),
        linewidths=0.75,
        linecolor="white",
        square=True,
        cbar_kws={"label": "Cosine similarity", "shrink": 0.72, "pad": 0.035, "fraction": 0.045},
    )
    ax_c.set_title("C. Similarity to remaining SC types", weight="bold", fontsize=11.0, pad=9)
    ax_c.set_xlabel("Remaining SC type")
    ax_c.set_ylabel("Removed ST-only type")
    ax_c.set_xticklabels(ax_c.get_xticklabels(), rotation=45, ha="right", fontsize=9.2)
    ax_c.set_yticklabels(ax_c.get_yticklabels(), rotation=0, fontsize=9.2)

    for ax in (ax_b, ax_c):
        ax.tick_params(axis="both", length=0)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, bbox_inches="tight", facecolor="white")
    fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    root = Path(args.project_root).resolve()
    out_dir = root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    sc_path = root / "data" / "raw" / "cell2location_mouse_brain" / "sc.h5ad"
    sc = ad.read_h5ad(sc_path)
    if sc.raw is None:
        raise ValueError("sc.h5ad lacks raw counts")
    labels = sc.obs["annotation_1"].astype(str)
    raw_counts = sc.raw.X
    if not sparse.issparse(raw_counts):
        raw_counts = sparse.csr_matrix(raw_counts)
    raw_counts = raw_counts.astype(np.float32).tocsr()
    norm_counts = log_normalize_rows(raw_counts)

    mis_rows = []
    sim_rows = []
    region_rows = []
    for case in CASES:
        target_mask = target_region_from_marker(
            root,
            sc,
            raw_counts,
            labels,
            case,
            args.target_quantile,
            args.n_markers,
        )
        mis_rows.append(forced_assignment_row(root, case, target_mask))
        sim_rows.append(similarity_row(norm_counts, labels, case))
        region_rows.append(
            {
                "label": case["label"],
                "section": case["section"],
                "sample": case["sample"],
                "drop_types": ";".join(case["drop_types"]),
                "target_region_spots": int(target_mask.sum()),
                "st_spots": int(len(target_mask)),
            }
        )

    misassignment = pd.DataFrame(mis_rows).reindex(columns=CATEGORY_ORDER)
    similarity = pd.DataFrame(sim_rows).reindex(columns=CATEGORY_ORDER)
    region_summary = pd.DataFrame(region_rows)

    prefix = args.output_prefix
    misassignment.to_csv(out_dir / f"{prefix}_cytospace_forced_assignment.csv")
    similarity.to_csv(out_dir / f"{prefix}_expression_similarity.csv")
    region_summary.to_csv(out_dir / f"{prefix}_target_region_summary.csv", index=False)
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
