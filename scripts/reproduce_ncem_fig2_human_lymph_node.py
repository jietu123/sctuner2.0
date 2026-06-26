#!/usr/bin/env python
"""Reproduce NCEM human lymph-node Figure 2A-F from prepared cell2location data."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import anndata as ad
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from scipy import sparse
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests


PANEL_A_RECEIVERS = ("B cells", "FDC", "Mast")
FOCAL_RECEIVER = "B cells"
FOCAL_SENDER = "FDC"
SIGNIFICANCE_THRESHOLD = 0.05
FOLD_CHANGE_THRESHOLD = 0.021671495152134755
VOLCANO_Y_CAP = 14.5
PAPER_STATE_COLORS = ("#2b8cbe", "#41b6c4", "#fdae61", "#4daf4a", "#bdbdbd")
PAPER_CELL_COLORS = (
    "#377eb8",
    "#ff7f00",
    "#4daf4a",
    "#f781bf",
    "#a65628",
    "#984ea3",
    "#999999",
    "#e41a1c",
    "#dede00",
    "#a6cee3",
    "#b2df8a",
    "#fb9a99",
    "#fdbf6f",
    "#cab2d6",
    "#ffff99",
    "#1b9e77",
)


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=root
        / "data"
        / "processed"
        / "ncem_fig2_human_lymph_node"
        / "cell2location_hvg"
        / "cell2location_lymphnode.h5ad",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=root / "visualizations" / "ncem_fig2_human_lymph_node_reproduction",
    )
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--force-statistics", action="store_true")
    parser.add_argument("--force-panel-a", action="store_true")
    return parser.parse_args()


def cell_type_names(adata: ad.AnnData) -> list[str]:
    names = adata.uns["node_type_names"]
    if isinstance(names, dict):
        return list(names.values())
    return list(names)


def add_target_cell(adata: ad.AnnData, names: list[str]) -> None:
    target = np.asarray(adata.obsm["node_types"])
    adata.obs["target_cell"] = pd.Categorical(
        np.asarray(names, dtype=object)[target.argmax(axis=1)],
        categories=names,
    )


def compute_panel_a(
    adata: ad.AnnData,
    names: list[str],
    cache_dir: Path,
    seed: int,
    force: bool,
) -> dict[str, dict[str, np.ndarray | pd.DataFrame]]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    results: dict[str, dict[str, np.ndarray | pd.DataFrame]] = {}

    for receiver in PANEL_A_RECEIVERS:
        slug = receiver.lower().replace(" ", "_")
        embedding_path = cache_dir / f"panel_a_{slug}_embedding.csv"
        composition_path = cache_dir / f"panel_a_{slug}_composition.csv"

        if embedding_path.exists() and composition_path.exists() and not force:
            embedding_df = pd.read_csv(embedding_path)
            composition_df = pd.read_csv(composition_path, index_col=0)
        else:
            subset = adata[adata.obs["target_cell"] == receiver].copy()
            sc.pp.neighbors(subset, n_neighbors=500, n_pcs=50, random_state=seed)
            sc.tl.louvain(subset, random_state=seed, key_added="subcluster")
            sc.tl.umap(subset, random_state=seed)

            embedding_df = pd.DataFrame(
                {
                    "umap_1": subset.obsm["X_umap"][:, 0],
                    "umap_2": subset.obsm["X_umap"][:, 1],
                    "subcluster": subset.obs["subcluster"].astype(str).to_numpy(),
                }
            )
            proportions = pd.DataFrame(
                np.asarray(subset.obsm["proportions"]),
                columns=names,
            )
            proportions["subcluster"] = embedding_df["subcluster"].to_numpy()
            composition_df = proportions.groupby("subcluster", observed=True)[names].mean().T
            embedding_df.to_csv(embedding_path, index=False)
            composition_df.to_csv(composition_path)

        results[receiver] = {
            "embedding": embedding_df,
            "composition": composition_df,
        }

    return results


def sparse_column_variance(matrix: sparse.spmatrix) -> np.ndarray:
    mean = np.asarray(matrix.mean(axis=0)).ravel()
    mean_sq = np.asarray(matrix.power(2).mean(axis=0)).ravel()
    return np.maximum(mean_sq - mean**2, np.finfo(float).tiny)


def compute_sender_receiver_statistics(
    adata: ad.AnnData,
    names: list[str],
    cache_path: Path,
    force: bool,
) -> dict[str, np.ndarray]:
    if cache_path.exists() and not force:
        cached = np.load(cache_path)
        return {key: cached[key] for key in cached.files}

    target = np.asarray(adata.obsm["node_types"], dtype=np.float64)
    proportions = np.asarray(adata.obsm["proportions"], dtype=np.float64)
    # Patsy expands ``target:proportions`` with proportions as the outer
    # dimension and target as the inner dimension. Retaining this order is
    # necessary to reproduce the tensor indexing used by the published code.
    interactions = (proportions[:, :, None] * target[:, None, :]).reshape(target.shape[0], -1)
    design = np.concatenate([target, interactions], axis=1)

    xtx_pinv = np.linalg.pinv(design.T @ design)
    expression = adata.X
    if sparse.issparse(expression):
        expression = expression.tocsr()
        xty = np.asarray(design.T @ expression)
        gene_variance = sparse_column_variance(expression)
    else:
        expression = np.asarray(expression, dtype=np.float64)
        xty = design.T @ expression
        gene_variance = np.maximum(np.var(expression, axis=0), np.finfo(float).tiny)

    parameters = xtx_pinv @ xty
    coefficient_variance = np.diag(xtx_pinv)[:, None] * gene_variance[None, :]
    standard_error = np.sqrt(np.maximum(coefficient_variance, np.finfo(float).tiny))
    z_score = np.abs(parameters / standard_error)
    pvalues = 2.0 * norm.sf(z_score)
    qvalues = multipletests(pvalues.ravel(), method="fdr_bh")[1].reshape(pvalues.shape)

    n_types = len(names)
    interaction_parameters = parameters[n_types:, :].reshape(n_types, n_types, adata.n_vars)
    interaction_pvalues = pvalues[n_types:, :].reshape(n_types, n_types, adata.n_vars)
    interaction_qvalues = qvalues[n_types:, :].reshape(n_types, n_types, adata.n_vars)
    is_significant = interaction_qvalues < SIGNIFICANCE_THRESHOLD

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        cache_path,
        fold_change=interaction_parameters,
        pvalues=interaction_pvalues,
        qvalues=interaction_qvalues,
        is_significant=is_significant,
    )
    return {
        "fold_change": interaction_parameters,
        "pvalues": interaction_pvalues,
        "qvalues": interaction_qvalues,
        "is_significant": is_significant,
    }


def select_focal_genes(
    fold_change: np.ndarray,
    qvalues: np.ndarray,
    names: list[str],
    genes: pd.Index,
) -> list[str]:
    receiver_idx = names.index(FOCAL_RECEIVER)
    sender_idx = names.index(FOCAL_SENDER)
    focal_fc = np.abs(fold_change[receiver_idx, sender_idx, :])
    focal_q = qvalues[receiver_idx, sender_idx, :]
    significant = focal_q <= SIGNIFICANCE_THRESHOLD
    if not significant.any():
        return list(genes[np.argsort(focal_q)[:20]])
    cutoff = focal_fc[significant].max() * 0.2
    selected = significant & (focal_fc >= cutoff)
    selected_idx = np.flatnonzero(selected)
    selected_idx = selected_idx[np.argsort(focal_q[selected_idx])]
    return list(genes[selected_idx])


def plot_panel_a(
    ax_umap: plt.Axes,
    ax_heatmap: plt.Axes,
    receiver: str,
    result: dict[str, np.ndarray | pd.DataFrame],
    show_y_labels: bool,
) -> None:
    embedding = result["embedding"]
    composition = result["composition"]
    assert isinstance(embedding, pd.DataFrame)
    assert isinstance(composition, pd.DataFrame)

    labels = pd.Categorical(embedding["subcluster"])
    palette = PAPER_STATE_COLORS
    for idx, cluster in enumerate(labels.categories):
        mask = labels == cluster
        ax_umap.scatter(
            embedding.loc[mask, "umap_1"],
            embedding.loc[mask, "umap_2"],
            s=4.2,
            color=palette[idx % len(palette)],
            linewidths=0,
            label=str(cluster),
            rasterized=True,
        )
    ax_umap.set_title(receiver, fontsize=7.2, fontweight="normal", pad=1.5)
    ax_umap.set_xticks([])
    ax_umap.set_yticks([])
    ax_umap.set_frame_on(False)
    ax_umap.legend(
        title="state",
        frameon=False,
        fontsize=4.8,
        title_fontsize=5.2,
        markerscale=1.5,
        ncol=1,
        loc="upper right",
        handletextpad=0.2,
        labelspacing=0.15,
        borderaxespad=0.1,
    )

    sns.heatmap(
        composition,
        ax=ax_heatmap,
        cmap="Reds",
        vmin=0,
        cbar=True,
        cbar_kws={"shrink": 0.7, "pad": 0.015},
        xticklabels=True,
        yticklabels=show_y_labels,
        linewidths=0,
    )
    ax_heatmap.set_xlabel(f"{receiver} subclusters", fontsize=5.7, labelpad=1.5)
    ax_heatmap.set_ylabel("Average neighborhood\ncomposition" if show_y_labels else "", fontsize=5.7)
    ax_heatmap.tick_params(axis="x", labelrotation=0, labelsize=5.2, length=0)
    ax_heatmap.tick_params(axis="y", labelsize=4.7, length=0)
    if ax_heatmap.collections and ax_heatmap.collections[0].colorbar is not None:
        colorbar = ax_heatmap.collections[0].colorbar
        colorbar.ax.tick_params(labelsize=4.3, length=1.5, width=0.4)
        colorbar.outline.set_linewidth(0.4)


def coupling_tables(
    fold_change: np.ndarray,
    is_significant: np.ndarray,
    names: list[str],
) -> pd.DataFrame:
    significant_count = is_significant.sum(axis=2)
    l1_norm = (np.abs(fold_change) * is_significant).sum(axis=2)
    rows = []
    for receiver_idx, receiver in enumerate(names):
        for sender_idx, sender in enumerate(names):
            if receiver == sender:
                continue
            rows.append(
                {
                    "receiver": receiver,
                    "sender": sender,
                    "significant_genes": int(significant_count[receiver_idx, sender_idx]),
                    "l1_norm": float(l1_norm[receiver_idx, sender_idx]),
                }
            )
    return pd.DataFrame(rows)


def plot_coupling_network(ax: plt.Axes, coupling: pd.DataFrame, names: list[str]) -> None:
    selected = coupling[coupling["significant_genes"] >= 200].copy()
    graph = nx.DiGraph()
    graph.add_nodes_from(names)
    for row in selected.itertuples(index=False):
        graph.add_edge(row.sender, row.receiver, weight=row.l1_norm)

    angles = np.linspace(0, 2 * np.pi, len(names), endpoint=False)
    angles = angles - angles[names.index("B cells")]
    position = {
        name: np.array([np.cos(angle), np.sin(angle)])
        for name, angle in zip(names, angles)
    }
    node_colors = [PAPER_CELL_COLORS[names.index(node)] for node in graph.nodes]
    weights = np.asarray([data["weight"] for _, _, data in graph.edges(data=True)], dtype=float)
    if weights.size:
        widths = 0.35 + 2.6 * weights / weights.max()
    else:
        widths = []

    nx.draw_networkx_nodes(
        graph,
        position,
        node_color=node_colors,
        node_size=120,
        edgecolors="#444444",
        linewidths=0.35,
        ax=ax,
    )
    nx.draw_networkx_edges(
        graph,
        position,
        width=widths,
        alpha=0.82,
        arrows=True,
        arrowsize=6,
        arrowstyle="-|>",
        connectionstyle="arc3,rad=0.08",
        edge_color="black",
        ax=ax,
    )
    for name, xy in position.items():
        radius = 1.18
        tx, ty = radius * xy
        horizontal = "left" if tx > 0.08 else "right" if tx < -0.08 else "center"
        vertical = "bottom" if ty > 0.08 else "top" if ty < -0.08 else "center"
        ax.text(tx, ty, name, fontsize=5.1, ha=horizontal, va=vertical)
    ax.set_title("b", loc="left", fontweight="bold", fontsize=8, pad=0)
    ax.set_xlim(-1.38, 1.38)
    ax.set_ylim(-1.35, 1.35)
    ax.axis("off")


def heatmap_values(
    fold_change: np.ndarray,
    qvalues: np.ndarray,
    names: list[str],
    genes: pd.Index,
    focal_genes: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    gene_idx = genes.get_indexer(focal_genes)
    receiver_idx = names.index(FOCAL_RECEIVER)
    sender_idx = names.index(FOCAL_SENDER)

    sender_values = fold_change[receiver_idx, :, :][:, gene_idx].copy()
    sender_qvalues = qvalues[receiver_idx, :, :][:, gene_idx]
    sender_values[sender_qvalues > SIGNIFICANCE_THRESHOLD] = 0.0
    sender_df = pd.DataFrame(sender_values, index=names, columns=focal_genes).drop(index=FOCAL_RECEIVER).T

    # Match the published receiver_effect implementation exactly. The
    # fold-change slice and significance mask use different tensor axes in
    # NCEM 0.1.5; changing this would alter the published panel.
    receiver_values = fold_change[sender_idx, :, :][:, gene_idx].copy()
    receiver_qvalues = qvalues[:, sender_idx, :][:, gene_idx]
    receiver_values[receiver_qvalues > SIGNIFICANCE_THRESHOLD] = 0.0
    receiver_df = pd.DataFrame(receiver_values, index=names, columns=focal_genes).drop(index=FOCAL_SENDER).T
    return sender_df, receiver_df


def plot_effect_heatmap(
    ax: plt.Axes,
    values: pd.DataFrame,
    title: str,
    xlabel: str,
    show_y_labels: bool,
    panel_label: str,
) -> None:
    vmax = max(float(np.nanmax(np.abs(values.to_numpy()))), 1e-8)
    sns.heatmap(
        values,
        ax=ax,
        cmap="seismic",
        center=0,
        vmin=-vmax,
        vmax=vmax,
        cbar_kws={
            "label": "Fold change",
            "shrink": 0.62,
            "orientation": "horizontal",
            "location": "top",
            "pad": 0.03,
        },
        yticklabels=show_y_labels,
        xticklabels=True,
    )
    ax.set_title(panel_label, loc="left", fontweight="bold", fontsize=8, pad=2)
    ax.set_xlabel(xlabel, fontsize=5.5, labelpad=1)
    ax.set_ylabel("gene" if show_y_labels else "", fontsize=5.5)
    ax.tick_params(axis="x", rotation=90, labelsize=4.3, length=0)
    ax.tick_params(axis="y", labelsize=4.3, length=0)
    if ax.collections and ax.collections[0].colorbar is not None:
        colorbar = ax.collections[0].colorbar
        colorbar.ax.tick_params(labelsize=4.2, length=1.5, width=0.4)
        colorbar.ax.xaxis.label.set_size(4.8)
        colorbar.outline.set_linewidth(0.4)


def plot_volcano(
    ax: plt.Axes,
    fold_change: np.ndarray,
    qvalues: np.ndarray,
    names: list[str],
    genes: pd.Index,
) -> pd.DataFrame:
    receiver_idx = names.index(FOCAL_RECEIVER)
    sender_idx = names.index(FOCAL_SENDER)
    fc = fold_change[receiver_idx, sender_idx, :]
    q = qvalues[receiver_idx, sender_idx, :]
    neglog_raw = -np.log10(np.clip(q, 1e-300, None))
    neglog = np.minimum(neglog_raw, VOLCANO_Y_CAP)

    significant = q < SIGNIFICANCE_THRESHOLD
    up = significant & (fc >= FOLD_CHANGE_THRESHOLD)
    down = significant & (fc <= -FOLD_CHANGE_THRESHOLD)
    nonsignificant = ~significant
    significant_neutral = significant & ~up & ~down
    ax.scatter(
        fc[nonsignificant],
        neglog[nonsignificant],
        facecolors="white",
        edgecolors="#8a8a8a",
        s=9,
        linewidths=0.45,
        alpha=0.95,
        rasterized=True,
    )
    ax.scatter(
        fc[significant_neutral],
        neglog[significant_neutral],
        facecolors="#a6a6a6",
        edgecolors="#555555",
        s=9,
        linewidths=0.4,
        alpha=0.95,
        rasterized=True,
    )
    ax.scatter(
        fc[up],
        neglog[up],
        facecolors="#e31a1c",
        edgecolors="#8b0000",
        s=11,
        linewidths=0.4,
        alpha=0.95,
        rasterized=True,
    )
    ax.scatter(
        fc[down],
        neglog[down],
        facecolors="#1f4eeb",
        edgecolors="#102878",
        s=11,
        linewidths=0.4,
        alpha=0.95,
        rasterized=True,
    )
    ax.axvline(FOLD_CHANGE_THRESHOLD, color="black", linestyle="--", linewidth=0.45)
    ax.axvline(-FOLD_CHANGE_THRESHOLD, color="black", linestyle="--", linewidth=0.45)
    ax.axhline(-np.log10(SIGNIFICANCE_THRESHOLD), color="black", linestyle="--", linewidth=0.45)
    ax.set_xlabel("log fold change", fontsize=5.5, labelpad=1)
    ax.set_ylabel("-log FDR-corrected p values", fontsize=5.5, labelpad=1)
    ax.set_title("e", loc="left", fontweight="bold", fontsize=8, pad=2)
    ax.tick_params(labelsize=4.7, width=0.4, length=2)
    ax.set_xlim(-0.11, 0.11)
    ax.set_xticks([-0.10, -0.05, 0.00, 0.05, 0.10])
    ax.set_ylim(-0.4, 15.2)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)

    label_candidates = np.flatnonzero(up | down)
    if label_candidates.size:
        rank = label_candidates[np.argsort(q[label_candidates])[:8]]
        for idx in rank:
            ax.annotate(
                str(genes[idx]),
                (fc[idx], neglog[idx]),
                xytext=(2, 2),
                textcoords="offset points",
                fontsize=4.2,
            )

    return pd.DataFrame(
        {
            "gene": genes,
            "fold_change": fc,
            "pvalue_adjusted": q,
            "negative_log10_qvalue": neglog_raw,
            "significant": significant,
            "direction": np.where(up, "up", np.where(down, "down", "other")),
        }
    )


def plot_sender_similarity(
    ax: plt.Axes,
    fold_change: np.ndarray,
    names: list[str],
) -> pd.DataFrame:
    receiver_idx = names.index(FOCAL_RECEIVER)
    correlation = np.corrcoef(fold_change[receiver_idx, :, :])
    linkage_matrix = linkage(correlation, method="average", metric="euclidean")
    order = leaves_list(linkage_matrix)
    ordered = correlation[np.ix_(order, order)]
    ordered_names = np.asarray(names, dtype=object)[order]
    sns.heatmap(
        ordered,
        ax=ax,
        cmap="Purples",
        vmin=0,
        vmax=1,
        square=True,
        xticklabels=ordered_names,
        yticklabels=ordered_names,
        cbar_kws={
            "label": "Correlation",
            "shrink": 0.58,
            "orientation": "horizontal",
            "location": "top",
            "pad": 0.03,
        },
    )
    ax.set_title("f", loc="left", fontweight="bold", fontsize=8, pad=2)
    ax.set_xlabel("Sender", fontsize=5.5, labelpad=1)
    ax.set_ylabel("Sender", fontsize=5.5, labelpad=1)
    ax.tick_params(axis="x", rotation=90, labelsize=4.2, length=0)
    ax.tick_params(axis="y", rotation=0, labelsize=4.2, length=0)
    if ax.collections and ax.collections[0].colorbar is not None:
        colorbar = ax.collections[0].colorbar
        colorbar.ax.tick_params(labelsize=4.2, length=1.5, width=0.4)
        colorbar.ax.xaxis.label.set_size(4.8)
        colorbar.outline.set_linewidth(0.4)
    return pd.DataFrame(correlation, index=names, columns=names)


def save_panel_a_figure(
    panel_a: dict[str, dict[str, np.ndarray | pd.DataFrame]],
    output_path: Path,
) -> None:
    fig = plt.figure(figsize=(8.1, 4.4))
    grid = GridSpec(2, 3, figure=fig, height_ratios=[0.86, 1.14], hspace=0.13, wspace=0.27)
    for column, receiver in enumerate(PANEL_A_RECEIVERS):
        plot_panel_a(
            fig.add_subplot(grid[0, column]),
            fig.add_subplot(grid[1, column]),
            receiver,
            panel_a[receiver],
            show_y_labels=column == 0,
        )
    fig.text(0.01, 0.985, "a", ha="left", va="top", fontsize=9, fontweight="bold")
    fig.savefig(output_path, dpi=400, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = args.output_dir / "intermediate"

    np.random.seed(args.seed)
    sc.settings.verbosity = 2
    sns.set_theme(style="white", context="paper")

    adata = ad.read_h5ad(args.input)
    names = cell_type_names(adata)
    add_target_cell(adata, names)

    panel_a = compute_panel_a(
        adata,
        names,
        cache_dir,
        args.seed,
        args.force_panel_a,
    )
    save_panel_a_figure(panel_a, args.output_dir / "figure2a_cell_states_and_neighborhoods.png")

    statistics = compute_sender_receiver_statistics(
        adata,
        names,
        cache_dir / "sender_receiver_statistics.npz",
        args.force_statistics,
    )
    fold_change = statistics["fold_change"]
    qvalues = statistics["qvalues"]
    is_significant = statistics["is_significant"].astype(bool)

    focal_genes = select_focal_genes(fold_change, qvalues, names, adata.var_names)
    pd.Series(focal_genes, name="gene").to_csv(
        args.output_dir / "figure2_focal_bcell_fdc_genes.csv",
        index=False,
    )
    sender_df, receiver_df = heatmap_values(
        fold_change,
        qvalues,
        names,
        adata.var_names,
        focal_genes,
    )
    sender_df.to_csv(args.output_dir / "figure2c_sender_effect_values.csv")
    receiver_df.to_csv(args.output_dir / "figure2d_receiver_effect_values.csv")

    coupling = coupling_tables(fold_change, is_significant, names)
    coupling.to_csv(args.output_dir / "figure2b_type_coupling.csv", index=False)

    fig_b, ax_b = plt.subplots(figsize=(4.2, 4.0))
    plot_coupling_network(ax_b, coupling, names)
    fig_b.savefig(
        args.output_dir / "figure2b_type_coupling_network.png",
        dpi=400,
        bbox_inches="tight",
        pad_inches=0.03,
    )
    plt.close(fig_b)

    fig_c, ax_c = plt.subplots(figsize=(4.0, 3.4))
    plot_effect_heatmap(
        ax_c,
        sender_df,
        f"Sender effect on {FOCAL_RECEIVER}",
        "Sender cell type",
        True,
        "c",
    )
    fig_c.savefig(
        args.output_dir / "figure2c_sender_effect_b_cells.png",
        dpi=400,
        bbox_inches="tight",
        pad_inches=0.03,
    )
    plt.close(fig_c)

    fig_d, ax_d = plt.subplots(figsize=(4.0, 3.4))
    plot_effect_heatmap(
        ax_d,
        receiver_df,
        f"Receiver effect of {FOCAL_SENDER}",
        "Receiver cell type",
        True,
        "d",
    )
    fig_d.savefig(
        args.output_dir / "figure2d_receiver_effect_fdc.png",
        dpi=400,
        bbox_inches="tight",
        pad_inches=0.03,
    )
    plt.close(fig_d)

    fig_e, ax_e = plt.subplots(figsize=(3.2, 3.4))
    volcano = plot_volcano(ax_e, fold_change, qvalues, names, adata.var_names)
    volcano.to_csv(args.output_dir / "figure2e_bcell_fdc_volcano_values.csv", index=False)
    fig_e.savefig(
        args.output_dir / "figure2e_bcell_fdc_volcano.png",
        dpi=400,
        bbox_inches="tight",
        pad_inches=0.03,
    )
    plt.close(fig_e)

    fig_f, ax_f = plt.subplots(figsize=(3.8, 3.4))
    sender_similarity = plot_sender_similarity(ax_f, fold_change, names)
    sender_similarity.to_csv(args.output_dir / "figure2f_bcell_sender_similarity.csv")
    fig_f.savefig(
        args.output_dir / "figure2f_bcell_sender_similarity.png",
        dpi=400,
        bbox_inches="tight",
        pad_inches=0.03,
    )
    plt.close(fig_f)

    fig = plt.figure(figsize=(12.2, 6.35))
    outer = GridSpec(
        2,
        12,
        figure=fig,
        height_ratios=[1.05, 1.0],
        hspace=0.18,
        wspace=1.15,
        left=0.045,
        right=0.985,
        bottom=0.08,
        top=0.975,
    )
    panel_a_grid = GridSpecFromSubplotSpec(
        2,
        3,
        subplot_spec=outer[0, :8],
        height_ratios=[0.86, 1.14],
        hspace=0.13,
        wspace=0.28,
    )
    for column, receiver in enumerate(PANEL_A_RECEIVERS):
        plot_panel_a(
            fig.add_subplot(panel_a_grid[0, column]),
            fig.add_subplot(panel_a_grid[1, column]),
            receiver,
            panel_a[receiver],
            show_y_labels=column == 0,
        )
    fig.text(0.015, 0.985, "a", ha="left", va="top", fontsize=9, fontweight="bold")
    plot_coupling_network(fig.add_subplot(outer[0, 8:12]), coupling, names)
    plot_effect_heatmap(
        fig.add_subplot(outer[1, 0:3]),
        sender_df,
        f"Sender effect on {FOCAL_RECEIVER}",
        "Sender cell type",
        True,
        "c",
    )
    plot_effect_heatmap(
        fig.add_subplot(outer[1, 3:6]),
        receiver_df,
        f"Receiver effect of {FOCAL_SENDER}",
        "Receiver cell type",
        True,
        "d",
    )
    plot_volcano(
        fig.add_subplot(outer[1, 6:9]),
        fold_change,
        qvalues,
        names,
        adata.var_names,
    )
    plot_sender_similarity(fig.add_subplot(outer[1, 9:12]), fold_change, names)
    fig.savefig(
        args.output_dir / "ncem_figure2a_f_reproduction.png",
        dpi=400,
        bbox_inches="tight",
        pad_inches=0.04,
    )
    plt.close(fig)

    summary = {
        "input": str(args.input),
        "observations": int(adata.n_obs),
        "genes": int(adata.n_vars),
        "cell_types": names,
        "panel_a_receivers": list(PANEL_A_RECEIVERS),
        "focal_receiver": FOCAL_RECEIVER,
        "focal_sender": FOCAL_SENDER,
        "focal_genes": focal_genes,
        "coupling_edges_at_least_200_genes": int((coupling["significant_genes"] >= 200).sum()),
        "significant_bcell_fdc_genes": int(
            is_significant[names.index(FOCAL_RECEIVER), names.index(FOCAL_SENDER), :].sum()
        ),
    }
    (args.output_dir / "reproduction_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
