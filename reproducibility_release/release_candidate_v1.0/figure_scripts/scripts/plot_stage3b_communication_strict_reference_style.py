#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import fisher_exact


CLASS_COLORS = {
    "observed_ST_supported": "#2CA25F",
    "high_confidence_false_communication": "#D7301F",
    "insufficient_observed_ST_evidence": "#756BB1",
    "insufficient_ligand_target_prior": "#969696",
}
CLASS_LABELS = {
    "observed_ST_supported": "Observed ST supported",
    "high_confidence_false_communication": "High-confidence false",
    "insufficient_observed_ST_evidence": "Insufficient observed ST",
    "insufficient_ligand_target_prior": "Insufficient prior",
}
PATHWAY_COLORS = {
    "NFkB": "#8C6BB1",
    "TGFb": "#E6550D",
    "TNFa": "#31A354",
    "Hypoxia": "#3182BD",
    "Estrogen": "#DD1C77",
    "No prior": "#BDBDBD",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Strict paper-reference-style plots for Stage3B communication validation."
    )
    parser.add_argument("--project_root", default=".")
    parser.add_argument(
        "--communication_dir",
        default="data/processed/stage3b_spatial_communication/brca_her2_ffpe_plasma",
    )
    parser.add_argument(
        "--validation_dir",
        default=(
            "data/processed/stage3b_communication_downstream_validation/"
            "brca_her2_ffpe_plasma"
        ),
    )
    parser.add_argument(
        "--out_dir",
        default=(
            "visualizations/stage3b_communication_validation/"
            "brca_her2_ffpe_plasma/strict_reference_style"
        ),
    )
    return parser.parse_args()


def pair_label(ligand: str, receptor: str) -> str:
    value = f"{ligand}->{receptor}"
    return (
        value.replace("TGFBR1_TGFBR2", "TGFBR1/2")
        .replace("BMPR1A_BMPR2", "BMPR1A/2")
        .replace("ITGA10_ITGB1", "ITGA10/ITGB1")
    )


def read_coords(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    id_col = "spot_id" if "spot_id" in frame.columns else frame.columns[0]
    frame[id_col] = frame[id_col].astype(str)
    frame = frame.set_index(id_col)
    lower = {column.lower(): column for column in frame.columns}
    if "pixel_col" in lower and "pixel_row" in lower:
        out = frame[[lower["pixel_col"], lower["pixel_row"]]].copy()
    elif "array_col" in lower and "array_row" in lower:
        out = frame[[lower["array_col"], lower["array_row"]]].copy()
    else:
        numeric = [
            column
            for column in frame.columns
            if pd.api.types.is_numeric_dtype(frame[column])
        ]
        if len(numeric) < 2:
            raise ValueError(f"Cannot infer coordinates from {path}")
        out = frame[[numeric[1], numeric[0]]].copy()
    out.columns = ["x", "y"]
    return out.apply(pd.to_numeric, errors="coerce")


def pathway_dotplot_data(validation: pd.DataFrame) -> pd.DataFrame:
    frame = validation.copy()
    frame["pathway"] = frame["assigned_progeny_pathway"].fillna("No prior")
    total_false = int(frame["classification"].eq("high_confidence_false_communication").sum())
    total_nonfalse = len(frame) - total_false
    rows: list[dict[str, object]] = []
    for pathway, group in frame.groupby("pathway", sort=False):
        false_count = int(group["classification"].eq("high_confidence_false_communication").sum())
        pair_count = len(group)
        nonfalse_count = pair_count - false_count
        other_false = total_false - false_count
        other_nonfalse = total_nonfalse - nonfalse_count
        _, pvalue = fisher_exact(
            [[false_count, nonfalse_count], [other_false, other_nonfalse]],
            alternative="greater",
        )
        rows.append(
            {
                "pathway": pathway,
                "selected_pairs": pair_count,
                "false_fraction": false_count / max(pair_count, 1),
                "fisher_p": pvalue,
                "false_pairs": false_count,
            }
        )
    out = pd.DataFrame(rows)
    return out.sort_values(["false_pairs", "selected_pairs"], ascending=False)


def draw_panel_a(axis: plt.Axes, validation: pd.DataFrame) -> pd.DataFrame:
    # Matches SpatialDM Fig. 2c: pathway rows, x = number of pairs,
    # dot colour = fraction, dot size = Fisher exact p-value.
    data = pathway_dotplot_data(validation)
    y = np.arange(len(data))
    scores = -np.log10(data["fisher_p"].clip(lower=1e-6))
    sizes = 18 + 95 * np.clip(scores, 0, 6) / 6
    scatter = axis.scatter(
        data["selected_pairs"],
        y,
        c=data["false_fraction"],
        s=sizes,
        cmap="Reds",
        vmin=0,
        vmax=1,
        edgecolors="#7A1F1F",
        linewidths=0.35,
    )
    axis.set_yticks(y)
    axis.set_yticklabels(data["pathway"], fontsize=7.5)
    axis.invert_yaxis()
    axis.set_xlabel("Number of dropout-induced LR pairs", fontsize=7.5)
    axis.set_ylabel("")
    xmax = max(3, int(data["selected_pairs"].max()) + 1)
    axis.set_xlim(0, xmax)
    axis.set_xticks(np.arange(0, xmax + 1, 1))
    axis.tick_params(axis="x", labelsize=7, length=2)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.grid(False)
    colorbar = plt.colorbar(scatter, ax=axis, fraction=0.048, pad=0.035)
    colorbar.set_label("False-communication fraction", fontsize=7)
    colorbar.ax.tick_params(labelsize=6.5, length=2)

    legend_values = np.array([0.05, 0.01, 0.001])
    legend_scores = -np.log10(legend_values)
    handles = [
        axis.scatter(
            [],
            [],
            s=18 + 95 * min(score, 6) / 6,
            c="#CFCFCF",
            edgecolors="#555555",
            linewidths=0.35,
            label=f"{p:g}",
        )
        for p, score in zip(legend_values, legend_scores)
    ]
    axis.legend(
        handles=handles,
        title="Fisher exact\np value",
        frameon=False,
        fontsize=6.2,
        title_fontsize=6.4,
        loc="lower right",
        bbox_to_anchor=(0.99, 0.01),
        handletextpad=0.7,
        borderaxespad=0,
    )
    return data


def load_local(path: Path) -> dict[str, np.ndarray]:
    data = np.load(path, allow_pickle=True)
    return {key: data[key] for key in data.files}


def pick_false_pair(
    validation: pd.DataFrame,
    dropout: dict[str, np.ndarray],
    svtuner: dict[str, np.ndarray],
) -> tuple[str, str]:
    frame = validation[
        validation["classification"].eq("high_confidence_false_communication")
    ].copy()
    rows: list[dict[str, object]] = []
    ligands = dropout["ligand"].astype(str)
    receptors = dropout["receptor"].astype(str)
    dropout_supported = dropout["supported"].astype(bool)
    svtuner_supported = svtuner["supported"].astype(bool)
    for row in frame.itertuples(index=False):
        mask = (ligands == str(row.ligand)) & (receptors == str(row.receptor))
        if not mask.any():
            continue
        pair_index = int(np.flatnonzero(mask)[0])
        dropout_values = dropout["local_moran"][:, pair_index].astype(float)
        svtuner_values = svtuner["local_moran"][:, pair_index].astype(float)
        dropout_valid = np.isfinite(dropout_values) & dropout_supported
        svtuner_valid = np.isfinite(svtuner_values) & svtuner_supported
        if not dropout_valid.any() or not svtuner_valid.any():
            continue
        threshold = float(np.nanquantile(dropout_values[dropout_valid], 0.95))
        cyto_count = int((dropout_values[dropout_valid] >= threshold).sum())
        svtuner_count = int((svtuner_values[svtuner_valid] >= threshold).sum())
        cyto_mean = float(np.nanmean(dropout_values[dropout_valid]))
        svtuner_mean = float(np.nanmean(svtuner_values[svtuner_valid]))
        rows.append(
            {
                "ligand": str(row.ligand),
                "receptor": str(row.receptor),
                "fixed_hotspot_reduction": cyto_count - svtuner_count,
                "mean_reduction": cyto_mean - svtuner_mean,
                "dropout_induced_target_change": float(
                    row.dropout_induced_target_change
                ),
            }
        )
    ranked = pd.DataFrame(rows)
    ranked = ranked[ranked["mean_reduction"] > 0].sort_values(
        ["fixed_hotspot_reduction", "mean_reduction", "dropout_induced_target_change"],
        ascending=False,
    )
    if ranked.empty:
        ranked = pd.DataFrame(rows).sort_values(
            ["fixed_hotspot_reduction", "dropout_induced_target_change"],
            ascending=False,
        )
    row = ranked.iloc[0]
    return str(row["ligand"]), str(row["receptor"])


def draw_spatial_map(
    axis: plt.Axes,
    coords: pd.DataFrame,
    local: dict[str, np.ndarray],
    pair_index: int,
    title: str,
    vmax: float,
    show_blank: bool,
) -> int:
    spots = local["spots"].astype(str)
    mapped = coords.reindex(pd.Index(spots))
    x = mapped["x"].to_numpy(float)
    y = mapped["y"].to_numpy(float)
    values = local["local_moran"][:, pair_index].astype(float)
    if show_blank:
        supported = local["supported"].astype(bool)
    else:
        supported = np.ones(len(spots), dtype=bool)
    valid = np.isfinite(x) & np.isfinite(y)
    blank = valid & (~supported)
    signal = valid & supported & np.isfinite(values)
    selected = signal & (values >= np.nanquantile(values[signal], 0.90))
    axis.scatter(
        x[valid],
        y[valid],
        c="#FFF7BC",
        s=7,
        marker="h",
        linewidths=0,
        alpha=1.0,
    )
    axis.scatter(
        x[signal],
        y[signal],
        c=np.clip(values[signal], 0, None),
        s=7,
        marker="h",
        cmap="YlOrRd",
        vmin=0,
        vmax=vmax,
        linewidths=0,
    )
    axis.scatter(
        x[blank],
        y[blank],
        c="#D9D9D9",
        s=7,
        marker="h",
        linewidths=0,
        alpha=0.9,
    )
    axis.set_title(title, fontsize=7.3, fontweight="bold", pad=3)
    axis.set_aspect("equal")
    axis.invert_yaxis()
    axis.set_xticks([])
    axis.set_yticks([])
    for spine in axis.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.6)
        spine.set_color("#555555")
    return int(selected.sum())


def draw_panel_b(
    fig: plt.Figure,
    spec,
    communication_dir: Path,
    validation: pd.DataFrame,
) -> dict[str, object]:
    # Matches SpatialDM Fig. 3f: selected local spot maps for one LR pair.
    coords = read_coords(communication_dir / "visium_spot_positions.csv")
    dropout = load_local(communication_dir / "cytospace_dropout_local_lr_statistics.npz")
    svtuner = load_local(communication_dir / "svtuner_stage3b_local_lr_statistics.npz")
    ligand, receptor = pick_false_pair(validation, dropout, svtuner)
    pair_mask = (dropout["ligand"].astype(str) == ligand) & (
        dropout["receptor"].astype(str) == receptor
    )
    if not pair_mask.any():
        raise KeyError(f"Pair not found in local statistics: {ligand}->{receptor}")
    pair_index = int(np.flatnonzero(pair_mask)[0])
    values = np.concatenate(
        [
            dropout["local_moran"][:, pair_index],
            svtuner["local_moran"][:, pair_index],
        ]
    ).astype(float)
    vmax = max(0.01, float(np.nanquantile(np.clip(values, 0, None), 0.985)))
    inner = spec.subgridspec(1, 3, width_ratios=[1, 1, 0.045], wspace=0.08)
    ax1 = fig.add_subplot(inner[0, 0])
    ax2 = fig.add_subplot(inner[0, 1])
    label = pair_label(ligand, receptor)
    n1 = draw_spatial_map(
        ax1,
        coords,
        dropout,
        pair_index,
        f"{label}\nCytoSPACE dropout",
        vmax,
        show_blank=False,
    )
    n2 = draw_spatial_map(
        ax2,
        coords,
        svtuner,
        pair_index,
        f"{label}\nSVTuner Stage3B",
        vmax,
        show_blank=True,
    )
    dropout_values = dropout["local_moran"][:, pair_index].astype(float)
    svtuner_values = svtuner["local_moran"][:, pair_index].astype(float)
    dropout_valid = np.isfinite(dropout_values) & dropout["supported"].astype(bool)
    svtuner_valid = np.isfinite(svtuner_values) & svtuner["supported"].astype(bool)
    hotspot_threshold = float(np.nanquantile(dropout_values[dropout_valid], 0.95))
    dropout_mean = float(np.nanmean(dropout_values[dropout_valid]))
    svtuner_mean = float(np.nanmean(svtuner_values[svtuner_valid]))
    dropout_hotspots = int((dropout_values[dropout_valid] >= hotspot_threshold).sum())
    svtuner_hotspots = int((svtuner_values[svtuner_valid] >= hotspot_threshold).sum())
    mean_reduction = (
        100.0 * (dropout_mean - svtuner_mean) / dropout_mean
        if dropout_mean > 0
        else np.nan
    )
    hotspot_reduction = (
        100.0 * (dropout_hotspots - svtuner_hotspots) / dropout_hotspots
        if dropout_hotspots > 0
        else np.nan
    )
    cax = fig.add_subplot(inner[0, 2])
    sm = plt.cm.ScalarMappable(norm=Normalize(vmin=0, vmax=vmax), cmap="YlOrRd")
    colorbar = fig.colorbar(sm, cax=cax)
    colorbar.set_label("Local Moran", fontsize=7)
    colorbar.ax.tick_params(labelsize=6.5, length=2)
    return {
        "pair": label,
        "cytospace_top10_spots": n1,
        "svtuner_top10_spots": n2,
        "fixed_hotspot_threshold": hotspot_threshold,
        "cytospace_mean_local_moran": dropout_mean,
        "svtuner_mean_local_moran": svtuner_mean,
        "mean_local_moran_reduction_percent": mean_reduction,
        "cytospace_fixed_threshold_hotspots": dropout_hotspots,
        "svtuner_fixed_threshold_hotspots": svtuner_hotspots,
        "fixed_threshold_hotspot_reduction_percent": hotspot_reduction,
    }


def ligand_target_matrix(
    validation: pd.DataFrame,
    targets: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.Series, pd.Series, pd.Series]:
    candidates = validation[
        validation["nichenet_prior_available"].fillna(False).astype(bool)
    ].copy()
    candidates["pair"] = [
        pair_label(l, r) for l, r in zip(candidates["ligand"], candidates["receptor"])
    ]
    class_rank = {
        "high_confidence_false_communication": 0,
        "insufficient_observed_ST_evidence": 1,
        "observed_ST_supported": 2,
        "insufficient_ligand_target_prior": 3,
    }
    candidates["class_rank"] = candidates["classification"].map(class_rank).fillna(9)
    candidates["pathway_rank"] = candidates["assigned_progeny_pathway"].fillna("No prior")
    candidates = candidates.sort_values(
        ["class_rank", "pathway_rank", "dropout_induced_target_change"],
        ascending=[True, True, False],
    )
    ligands = candidates["ligand"].drop_duplicates().tolist()
    selected_targets: list[str] = []
    for ligand in ligands:
        ligand_targets = (
            targets[targets["ligand"].eq(ligand)]
            .sort_values("regulatory_potential", ascending=False)
            .head(4)["target"]
            .astype(str)
            .tolist()
        )
        for target in ligand_targets:
            if target not in selected_targets:
                selected_targets.append(target)
    selected_targets = selected_targets[:42]
    matrix = pd.DataFrame(0.0, index=ligands, columns=selected_targets)
    for row in targets.itertuples(index=False):
        if row.ligand in matrix.index and row.target in matrix.columns:
            matrix.loc[row.ligand, row.target] = max(
                matrix.loc[row.ligand, row.target], float(row.regulatory_potential)
            )
    max_value = matrix.to_numpy().max()
    if max_value > 0:
        matrix = matrix / max_value
    ligand_meta = candidates.drop_duplicates("ligand").set_index("ligand")
    class_by_ligand = ligand_meta["classification"]
    pathway_by_ligand = ligand_meta["assigned_progeny_pathway"].fillna("No prior")
    ligand_score = ligand_meta["dropout_induced_target_change"].astype(float)
    score_max = float(ligand_score.max()) if len(ligand_score) else 0.0
    if score_max > 0:
        ligand_score = ligand_score / score_max
    matrix = matrix.loc[class_by_ligand.index]
    ligand_score = ligand_score.reindex(matrix.index).fillna(0.0)
    # Renoir Fig. 4g/h obtains diagonal-like modules after ranking ligand-target
    # activities by communication domains. Here each ligand contributes its own
    # top targets in row order, producing the same module-preserving visual logic
    # without changing the underlying regulatory-potential values.
    return matrix, ligand_score, class_by_ligand, pathway_by_ligand


def draw_color_strip(
    axis: plt.Axes,
    colors: list[str],
    orientation: str,
) -> None:
    axis.set_axis_off()
    if orientation == "vertical":
        for index, color in enumerate(colors):
            axis.add_patch(Rectangle((0, index), 1, 1, facecolor=color, edgecolor="none"))
        axis.set_xlim(0, 1)
        axis.set_ylim(len(colors), 0)
    else:
        for index, color in enumerate(colors):
            axis.add_patch(Rectangle((index, 0), 1, 1, facecolor=color, edgecolor="none"))
        axis.set_xlim(0, len(colors))
        axis.set_ylim(0, 1)


def draw_panel_c(
    fig: plt.Figure,
    spec,
    validation: pd.DataFrame,
    validation_dir: Path,
) -> tuple[pd.DataFrame, pd.Series, pd.Series, pd.Series]:
    # Matches Renoir Fig. 4g/h: sparse ligand-target heatmap with annotation bars.
    targets = pd.read_csv(validation_dir / "nichenet_top50_targets_used.csv")
    matrix, ligand_score, class_by_ligand, pathway_by_ligand = ligand_target_matrix(
        validation, targets
    )
    inner = spec.subgridspec(
        4,
        6,
        width_ratios=[0.092, 0.018, 0.018, 0.038, 1.0, 0.004],
        height_ratios=[1.0, 0.040, 0.30, 0.16],
        wspace=0.018,
        hspace=0.045,
    )
    label_ax = fig.add_subplot(inner[0, 0])
    class_ax = fig.add_subplot(inner[0, 1])
    pathway_ax = fig.add_subplot(inner[0, 2])
    ligand_ax = fig.add_subplot(inner[0, 3])
    heat_ax = fig.add_subplot(inner[0, 4])
    spacer_ax = fig.add_subplot(inner[0, 5])
    target_ax = fig.add_subplot(inner[1, 4])
    target_label_ax = fig.add_subplot(inner[2, 4])
    cbar_grid = inner[3, 4].subgridspec(2, 1, hspace=1.15)
    target_cbar_ax = fig.add_subplot(cbar_grid[0, 0])
    ligand_cbar_ax = fig.add_subplot(cbar_grid[1, 0])
    spacer_ax.set_axis_off()
    label_ax.set_xlim(0, 1)
    label_ax.set_ylim(matrix.shape[0] - 0.5, -0.5)
    label_ax.set_xticks([])
    label_ax.set_yticks([])
    for index, ligand in enumerate(matrix.index):
        label_ax.text(
            0.985,
            index,
            ligand,
            ha="right",
            va="center",
            fontsize=6.5,
            clip_on=False,
        )
    label_ax.set_ylabel("Ligand", fontsize=7, labelpad=12)
    for spine in label_ax.spines.values():
        spine.set_visible(False)
    draw_color_strip(
        class_ax,
        [CLASS_COLORS.get(value, "#999999") for value in class_by_ligand],
        "vertical",
    )
    draw_color_strip(
        pathway_ax,
        [PATHWAY_COLORS.get(value, "#BDBDBD") for value in pathway_by_ligand],
        "vertical",
    )
    target_colors = sns.color_palette("husl", n_colors=matrix.shape[1]).as_hex()
    draw_color_strip(target_ax, target_colors, "horizontal")

    ligand_image = ligand_ax.imshow(
        ligand_score.to_numpy()[:, None],
        aspect="auto",
        cmap="RdYlBu_r",
        vmin=0,
        vmax=1,
    )
    ligand_ax.set_yticks([])
    ligand_ax.set_xticks([])
    ligand_ax.tick_params(length=0, pad=1.5)
    for spine in ligand_ax.spines.values():
        spine.set_linewidth(0.5)
        spine.set_color("#555555")

    image = heat_ax.imshow(matrix.to_numpy(), aspect="auto", cmap="BuPu", vmin=0, vmax=1)
    heat_ax.set_yticks(np.arange(matrix.shape[0]))
    heat_ax.set_yticklabels([])
    heat_ax.set_xticks([])
    heat_ax.tick_params(length=0)
    heat_ax.set_ylabel("")
    for spine in heat_ax.spines.values():
        spine.set_linewidth(0.5)
        spine.set_color("#555555")

    target_ax.set_axis_off()
    target_label_ax.set_xlim(-0.5, matrix.shape[1] - 0.5)
    target_label_ax.set_ylim(0, 1)
    target_label_ax.set_axis_off()
    for index, target in enumerate(matrix.columns):
        target_label_ax.text(
            index,
            0.98,
            target,
            rotation=90,
            ha="center",
            va="top",
            fontsize=5.4,
            clip_on=False,
        )
    target_label_ax.text(
        (matrix.shape[1] - 1) / 2,
        -0.10,
        "Target",
        ha="center",
        va="top",
        fontsize=7,
        clip_on=False,
    )
    ligand_ax.text(
        0.50,
        -0.085,
        "Ligand\nscore",
        rotation=90,
        ha="center",
        va="top",
        fontsize=6.0,
        transform=ligand_ax.transAxes,
        clip_on=False,
    )

    colorbar = fig.colorbar(image, cax=target_cbar_ax, orientation="horizontal")
    colorbar.set_label("")
    colorbar.ax.tick_params(labelsize=5.6, length=2, pad=0.5)
    target_cbar_ax.text(
        -0.012,
        0.50,
        "Average neighborhood score",
        ha="right",
        va="center",
        fontsize=6.0,
        transform=target_cbar_ax.transAxes,
        clip_on=False,
    )
    ligand_cbar = fig.colorbar(ligand_image, cax=ligand_cbar_ax, orientation="horizontal")
    ligand_cbar.set_label("")
    ligand_cbar.ax.tick_params(labelsize=5.6, length=2, pad=0.5)
    ligand_cbar_ax.text(
        -0.012,
        0.50,
        "Ligand score",
        ha="right",
        va="center",
        fontsize=6.0,
        transform=ligand_cbar_ax.transAxes,
        clip_on=False,
    )
    return matrix, ligand_score, class_by_ligand, pathway_by_ligand


def save_tables(
    out_dir: Path,
    panel_a: pd.DataFrame,
    panel_b: dict[str, object],
    matrix: pd.DataFrame,
    ligand_score: pd.Series,
    class_by_ligand: pd.Series,
    pathway_by_ligand: pd.Series,
) -> None:
    panel_a.to_csv(out_dir / "panel_a_spatialdm_fig2c_pathway_dotplot_values.csv", index=False)
    pd.DataFrame([panel_b]).to_csv(out_dir / "panel_b_spatialdm_fig3f_hotspot_summary.csv", index=False)
    matrix.to_csv(out_dir / "panel_c_renoir_fig4gh_ligand_target_matrix.csv")
    ligand_score.rename("ligand_score").to_csv(
        out_dir / "panel_c_renoir_fig4gh_ligand_score.csv",
        header=True,
    )
    pd.DataFrame(
        {
            "ligand": matrix.index,
            "ligand_score": ligand_score.reindex(matrix.index).values,
            "classification": class_by_ligand.reindex(matrix.index).values,
            "pathway": pathway_by_ligand.reindex(matrix.index).values,
        }
    ).to_csv(out_dir / "panel_c_ligand_annotations.csv", index=False)


def save_panel_b_only(
    communication_dir: Path,
    validation: pd.DataFrame,
    out_dir: Path,
) -> None:
    """Export the existing Panel B visual grammar as a standalone figure."""
    fig = plt.figure(figsize=(7.4, 3.75), dpi=420, facecolor="white")
    grid = fig.add_gridspec(
        1,
        1,
        left=0.075,
        right=0.94,
        top=0.82,
        bottom=0.25,
    )
    panel_b = draw_panel_b(fig, grid[0, 0], communication_dir, validation)
    fig.text(0.025, 0.935, "B", fontsize=12, fontweight="bold", va="top")
    fig.text(
        0.075,
        0.935,
        "Local LR hotspots after reference dropout",
        fontsize=10.5,
        fontweight="bold",
        va="top",
    )
    fig.text(
        0.50,
        0.155,
        (
            "Mean local Moran I: "
            f"{panel_b['cytospace_mean_local_moran']:.4f} "
            "\N{RIGHTWARDS ARROW} "
            f"{panel_b['svtuner_mean_local_moran']:.4f}    |    "
            "Absolute change: "
            f"\N{GREEK CAPITAL LETTER DELTA}I = "
            f"{panel_b['svtuner_mean_local_moran'] - panel_b['cytospace_mean_local_moran']:.4f}"
        ),
        ha="center",
        va="center",
        fontsize=7.5,
        fontweight="bold",
    )
    fig.text(
        0.50,
        0.105,
        (
            "Fixed-threshold hotspots: "
            f"{panel_b['cytospace_fixed_threshold_hotspots']} "
            "\N{RIGHTWARDS ARROW} "
            f"{panel_b['svtuner_fixed_threshold_hotspots']}"
        ),
        ha="center",
        va="center",
        fontsize=7.5,
        fontweight="bold",
    )
    fig.text(
        0.075,
        0.025,
        (
            "Metrics use the common supported-spot set; hotspot threshold is fixed "
            f"at the CytoSPACE-dropout P95 (Local Moran I >= "
            f"{panel_b['fixed_hotspot_threshold']:.3f}). Gray spots are Stage3B-withheld."
        ),
        fontsize=6.3,
        color="#444444",
    )
    png = out_dir / "stage3b_communication_panel_b_local_lr_hotspots.png"
    pdf = out_dir / "stage3b_communication_panel_b_local_lr_hotspots.pdf"
    fig.savefig(png, dpi=420, facecolor="white")
    fig.savefig(pdf, facecolor="white")
    plt.close(fig)
    print(f"[done] {png}")


def plot_figure(root: Path, communication_dir: Path, validation_dir: Path, out_dir: Path) -> None:
    validation = pd.read_csv(validation_dir / "observed_st_downstream_validation.csv")
    out_dir.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="white", context="paper")
    fig = plt.figure(figsize=(11.4, 8.2), dpi=420, facecolor="white")
    grid = fig.add_gridspec(
        2,
        2,
        width_ratios=[0.44, 0.56],
        height_ratios=[0.46, 0.54],
        left=0.065,
        right=0.965,
        top=0.905,
        bottom=0.105,
        wspace=0.28,
        hspace=0.50,
    )
    ax_a = fig.add_subplot(grid[0, 0])
    panel_a = draw_panel_a(ax_a, validation)
    panel_b = draw_panel_b(fig, grid[0, 1], communication_dir, validation)
    matrix, ligand_score, class_by_ligand, pathway_by_ligand = draw_panel_c(
        fig, grid[1, :], validation, validation_dir
    )

    fig.text(0.037, 0.944, "A", fontsize=10, fontweight="bold")
    fig.text(0.065, 0.944, "Pathway enrichment of dropout-induced false communications", fontsize=8.5, fontweight="bold")
    fig.text(0.515, 0.944, "B", fontsize=10, fontweight="bold")
    fig.text(0.543, 0.944, "Local LR hotspots after reference dropout", fontsize=8.5, fontweight="bold")
    fig.text(0.037, 0.470, "C", fontsize=10, fontweight="bold")
    fig.text(0.065, 0.470, "Ligand-target evidence underlying downstream validation", fontsize=8.5, fontweight="bold")

    class_handles = [
        Patch(facecolor=color, edgecolor="none", label=CLASS_LABELS[key])
        for key, color in CLASS_COLORS.items()
    ]
    pathway_handles = [
        Patch(facecolor=color, edgecolor="none", label=key)
        for key, color in PATHWAY_COLORS.items()
    ]
    fig.legend(
        handles=class_handles,
        loc="lower left",
        bbox_to_anchor=(0.065, 0.025),
        frameon=False,
        fontsize=6.2,
        ncol=4,
        columnspacing=0.8,
        handlelength=0.9,
    )
    fig.legend(
        handles=pathway_handles,
        loc="lower right",
        bbox_to_anchor=(0.965, 0.025),
        frameon=False,
        fontsize=6.2,
        ncol=6,
        columnspacing=0.65,
        handlelength=0.9,
    )
    fig.text(
        0.065,
        0.006,
        (
            "BRCA HER2 FFPE Plasma-cell reference dropout. The panels intentionally reuse the "
            "visual grammar of SpatialDM Fig. 2c, SpatialDM Fig. 3f, and Renoir Fig. 4g/h; "
            "statistics come from SVTuner Stage3B communication inference and observed-ST validation."
        ),
        fontsize=5.8,
        color="#444444",
    )
    png = out_dir / "stage3b_communication_strict_reference_style.png"
    pdf = out_dir / "stage3b_communication_strict_reference_style.pdf"
    fig.savefig(png, dpi=420, facecolor="white")
    fig.savefig(pdf, facecolor="white")
    plt.close(fig)
    save_tables(out_dir, panel_a, panel_b, matrix, ligand_score, class_by_ligand, pathway_by_ligand)
    save_panel_b_only(communication_dir, validation, out_dir)
    print(f"[done] {png}")


def main() -> int:
    args = parse_args()
    root = Path(args.project_root).resolve()
    plot_figure(
        root=root,
        communication_dir=root / args.communication_dir,
        validation_dir=root / args.validation_dir,
        out_dir=root / args.out_dir,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
