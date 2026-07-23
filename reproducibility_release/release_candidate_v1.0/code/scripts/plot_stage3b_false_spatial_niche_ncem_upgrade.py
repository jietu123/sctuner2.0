#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import norm
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


CONDITIONS = ("Full reference", "CytoSPACE dropout", "SVTuner")
STATE_COLORS = ("#4C78A8", "#F58518", "#54A24B", "#E45756")
NODE_COLORS = (
    "#4C78A8",
    "#F58518",
    "#54A24B",
    "#E45756",
    "#72B7B2",
    "#B279A2",
    "#FF9DA6",
    "#9D755D",
    "#BAB0AC",
    "#7B61A8",
)
SHORT_LABELS = {
    "Endothelial/Blood": "Endothelial",
    "Extraembryonic": "Extraemb.",
    "Cardiomyocytes": "Cardiomyo.",
    "Endoderm/Gut": "Endoderm/Gut",
    "Unknown_sc_only": "Unknown",
    "Abstained": "Abstained",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create an NCEM Figure 2-style Stage3B false-niche experiment."
    )
    parser.add_argument("--project_root", default=".")
    parser.add_argument(
        "--metadata_csv",
        default=(
            "visualizations/stage3b_realdata_candidate_scan/spatial_9x2/"
            "stage3b_reference_dropout_spatial_stack_recommended_6x2_metadata.csv"
        ),
    )
    parser.add_argument(
        "--sample",
        default="cytospace_fig2d_tme_brca_her2_ffpe_sc_missing_plasma_cells",
    )
    parser.add_argument("--target_quantile", type=float, default=0.85)
    parser.add_argument("--neighbors", type=int, default=6)
    parser.add_argument("--states", type=int, default=4)
    parser.add_argument("--bootstrap", type=int, default=250)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--out_dir",
        default=(
            "visualizations/stage3b_false_spatial_niche/"
            "candidate_validation/brca_her2_plasma"
        ),
    )
    return parser.parse_args()


def read_coords(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    id_col = "spot_id" if "spot_id" in frame.columns else frame.columns[0]
    frame[id_col] = frame[id_col].astype(str)
    frame = frame.drop_duplicates(id_col).set_index(id_col)
    lower = {column.lower(): column for column in frame.columns}
    if "col" in lower and "row" in lower:
        x_col, y_col = lower["col"], lower["row"]
    elif "x" in lower and "y" in lower:
        x_col, y_col = lower["x"], lower["y"]
    else:
        numeric = [
            column
            for column in frame.columns
            if pd.api.types.is_numeric_dtype(frame[column])
        ]
        if len(numeric) < 2:
            raise ValueError(f"Cannot infer coordinates from {path}")
        y_col, x_col = numeric[:2]
    result = frame[[x_col, y_col]].apply(pd.to_numeric, errors="coerce").dropna()
    result.columns = ["x", "y"]
    return result


def marker_percentile(path: Path, genes: list[str]) -> pd.Series:
    header = pd.read_csv(path, nrows=0).columns.tolist()
    id_col = header[0]
    present = [gene for gene in genes if gene in header]
    if not present:
        raise ValueError("None of the evaluation marker genes are present in ST.")
    expression = pd.read_csv(path, index_col=0, usecols=[id_col, *present])
    expression.index = expression.index.astype(str)
    expression = expression.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return expression.mean(axis=1).rank(method="average", pct=True)


def read_fractional(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, index_col=0)
    frame.index = frame.index.astype(str)
    frame = frame.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    frame = frame.drop(columns=["Unknown_sc_only"], errors="ignore")
    totals = frame.sum(axis=1)
    nonzero = totals > 0
    frame.loc[nonzero] = frame.loc[nonzero].div(totals.loc[nonzero], axis=0)
    return frame


def fractional_path(root: Path, sample: str, suffix: str) -> Path:
    return (
        root
        / "result"
        / sample
        / f"stage4_cytospace{suffix}"
        / "cytospace_output"
        / "fractional_abundances_by_spot.csv"
    )


def spatial_edges(coords: np.ndarray, neighbors: int) -> tuple[np.ndarray, np.ndarray]:
    k = min(neighbors + 1, len(coords))
    tree = cKDTree(coords)
    _, indices = tree.query(coords, k=k)
    source = np.repeat(np.arange(len(coords)), k - 1)
    destination = indices[:, 1:].reshape(-1)
    valid = source != destination
    return source[valid], destination[valid]


def augment_abstention(values: np.ndarray) -> np.ndarray:
    blank = (values.sum(axis=1) <= 1e-12).astype(float)[:, None]
    return np.hstack([values, blank])


def neighborhood_composition(
    values: np.ndarray,
    source: np.ndarray,
    destination: np.ndarray,
) -> np.ndarray:
    output = np.zeros_like(values, dtype=float)
    counts = np.zeros(len(values), dtype=float)
    np.add.at(output, source, values[destination])
    np.add.at(counts, source, 1.0)
    valid = counts > 0
    output[valid] /= counts[valid, None]
    return output


def contact_matrix(
    values: np.ndarray,
    source: np.ndarray,
    destination: np.ndarray,
) -> np.ndarray:
    if len(source) == 0:
        return np.zeros((values.shape[1], values.shape[1]), dtype=float)
    return values[source].T @ values[destination] / float(len(source))


def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    flat = np.asarray(pvalues, dtype=float).ravel()
    order = np.argsort(flat)
    ranked = flat[order]
    adjusted = ranked * len(flat) / np.arange(1, len(flat) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.clip(adjusted, 0.0, 1.0)
    return result.reshape(np.asarray(pvalues).shape)


def bootstrap_differences(
    ideal: np.ndarray,
    methods: dict[str, np.ndarray],
    source: np.ndarray,
    destination: np.ndarray,
    repeats: int,
    rng: np.random.Generator,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], dict[str, np.ndarray]]:
    point_ideal = contact_matrix(ideal, source, destination)
    point = {
        method: contact_matrix(values, source, destination) - point_ideal
        for method, values in methods.items()
    }
    samples = {
        method: np.empty((repeats, *point_ideal.shape), dtype=np.float32)
        for method in methods
    }
    for repeat in range(repeats):
        selected = rng.integers(0, len(source), size=len(source))
        sampled_source = source[selected]
        sampled_destination = destination[selected]
        ideal_contact = contact_matrix(ideal, sampled_source, sampled_destination)
        for method, values in methods.items():
            samples[method][repeat] = (
                contact_matrix(values, sampled_source, sampled_destination)
                - ideal_contact
            )
    qvalues: dict[str, np.ndarray] = {}
    standard_errors: dict[str, np.ndarray] = {}
    for method, values in samples.items():
        standard_error = values.std(axis=0, ddof=1)
        zscore = np.divide(
            point[method],
            standard_error,
            out=np.zeros_like(point[method]),
            where=standard_error > 1e-12,
        )
        pvalues = 2.0 * norm.sf(np.abs(zscore))
        qvalues[method] = benjamini_hochberg(np.minimum(pvalues, 1.0))
        standard_errors[method] = standard_error
    return point, qvalues, standard_errors


def prepare_experiment(
    root: Path,
    metadata: pd.DataFrame,
    sample: str,
    target_quantile: float,
    neighbors: int,
    states: int,
    repeats: int,
    seed: int,
) -> dict[str, object]:
    rows = metadata[metadata["sample"].eq(sample)]
    if rows.empty:
        raise KeyError(f"Sample not found in metadata: {sample}")
    row = rows.iloc[0]
    source_sample = str(row["source_sample"])
    sample = str(row["sample"])
    export_dir = (
        root
        / "data"
        / "processed"
        / str(row["storage_group"])
        / source_sample
        / "stage1_preprocess"
        / "exported"
    )
    genes = [gene for gene in str(row["marker_genes_list"]).split(";") if gene]
    coords = read_coords(export_dir / "st_coordinates.csv")
    percentiles = marker_percentile(
        export_dir / "st_expression_normalized.csv",
        genes,
    )
    paths = {
        "Full reference": fractional_path(
            root, source_sample, "_niche_full_reference"
        ),
        "CytoSPACE dropout": fractional_path(
            root, sample, "_niche_dropout_baseline"
        ),
        "SVTuner": fractional_path(
            root, sample, "_niche_stage3b_blank"
        ),
    }
    frames = {name: read_fractional(path) for name, path in paths.items()}
    common = coords.index.intersection(percentiles.index)
    for frame in frames.values():
        common = common.intersection(frame.index)
    coords = coords.loc[common]
    percentiles = percentiles.loc[common]
    type_order = sorted(
        set().union(*(set(frame.columns.astype(str)) for frame in frames.values()))
    )
    aligned: dict[str, np.ndarray] = {}
    for condition, frame in frames.items():
        aligned[condition] = frame.reindex(
            index=common,
            columns=type_order,
            fill_value=0.0,
        ).to_numpy(dtype=float)
    target_core = percentiles >= percentiles.quantile(target_quantile)
    source, destination = spatial_edges(
        coords[["x", "y"]].to_numpy(dtype=float),
        neighbors,
    )

    augmented = {
        condition: augment_abstention(values)
        for condition, values in aligned.items()
    }
    labels = [*type_order, "Abstained"]
    ideal = augmented["Full reference"].copy()
    ideal[target_core.to_numpy(dtype=bool), :] = 0.0
    ideal[target_core.to_numpy(dtype=bool), labels.index("Abstained")] = 1.0

    target_edge_mask = target_core.to_numpy(dtype=bool)[source]
    target_source = source[target_edge_mask]
    target_destination = destination[target_edge_mask]
    neighborhoods = {
        condition: neighborhood_composition(values, source, destination)
        for condition, values in augmented.items()
    }
    ideal_neighborhood = neighborhood_composition(ideal, source, destination)
    core_indices = np.flatnonzero(target_core.to_numpy(dtype=bool))
    full_core = neighborhoods["Full reference"][core_indices]
    biological_positions = [
        index for index, label in enumerate(labels) if label != "Abstained"
    ]
    scaler = StandardScaler()
    scaled_full = scaler.fit_transform(full_core[:, biological_positions])
    n_states = min(states, len(core_indices))
    clusterer = KMeans(n_clusters=n_states, random_state=seed, n_init=50)
    state_labels = clusterer.fit_predict(scaled_full)
    stacked = np.vstack(
        [neighborhoods[condition][core_indices] for condition in CONDITIONS]
    )
    embedding = PCA(n_components=2, random_state=seed).fit_transform(
        StandardScaler().fit_transform(stacked)
    )
    embeddings = {
        condition: embedding[
            index * len(core_indices) : (index + 1) * len(core_indices)
        ]
        for index, condition in enumerate(CONDITIONS)
    }
    state_compositions: dict[str, np.ndarray] = {}
    for condition in CONDITIONS:
        state_compositions[condition] = np.vstack(
            [
                neighborhoods[condition][core_indices][state_labels == state].mean(
                    axis=0
                )
                for state in range(n_states)
            ]
        ).T

    methods = {
        "CytoSPACE dropout": augmented["CytoSPACE dropout"],
        "SVTuner": augmented["SVTuner"],
    }
    point, qvalues, standard_errors = bootstrap_differences(
        ideal,
        methods,
        target_source,
        target_destination,
        repeats,
        np.random.default_rng(seed),
    )
    contacts = {
        "Ideal abstention": contact_matrix(
            ideal, target_source, target_destination
        ),
        **{
            condition: contact_matrix(
                augmented[condition],
                target_source,
                target_destination,
            )
            for condition in CONDITIONS
        },
    }
    return {
        "row": row,
        "coords": coords,
        "percentiles": percentiles,
        "target_core": target_core,
        "labels": labels,
        "augmented": augmented,
        "ideal": ideal,
        "ideal_neighborhood": ideal_neighborhood,
        "source": source,
        "destination": destination,
        "target_source": target_source,
        "target_destination": target_destination,
        "core_indices": core_indices,
        "state_labels": state_labels,
        "embeddings": embeddings,
        "state_compositions": state_compositions,
        "contacts": contacts,
        "differences": point,
        "qvalues": qvalues,
        "standard_errors": standard_errors,
    }


def abbreviate(label: str) -> str:
    return SHORT_LABELS.get(label, label)


def draw_panel_a(
    fig: plt.Figure,
    subspec: plt.GridSpec,
    experiment: dict[str, object],
) -> None:
    nested = subspec.subgridspec(
        2,
        3,
        height_ratios=[0.9, 1.15],
        hspace=0.12,
        wspace=0.18,
    )
    labels = list(experiment["labels"])
    state_labels = np.asarray(experiment["state_labels"])
    compositions = experiment["state_compositions"]
    vmax = max(
        float(np.asarray(compositions[condition]).max())
        for condition in CONDITIONS
    )
    heatmap_artist = None
    for column, condition in enumerate(CONDITIONS):
        scatter_axis = fig.add_subplot(nested[0, column])
        embedding = np.asarray(experiment["embeddings"][condition])
        for state in np.unique(state_labels):
            selected = state_labels == state
            scatter_axis.scatter(
                embedding[selected, 0],
                embedding[selected, 1],
                s=5.5,
                color=STATE_COLORS[int(state) % len(STATE_COLORS)],
                linewidths=0,
                alpha=0.86,
                rasterized=True,
            )
        scatter_axis.set_xticks([])
        scatter_axis.set_yticks([])
        scatter_axis.set_title(
            condition.replace(" dropout", ""),
            fontsize=7.2,
            fontweight="bold",
            pad=1.5,
        )
        sns.despine(ax=scatter_axis, left=True, bottom=True)

        heat_axis = fig.add_subplot(nested[1, column])
        matrix = np.asarray(compositions[condition])
        heatmap_artist = sns.heatmap(
            matrix,
            ax=heat_axis,
            cmap="Reds",
            vmin=0,
            vmax=vmax,
            cbar=False,
            xticklabels=[str(index + 1) for index in range(matrix.shape[1])],
            yticklabels=[abbreviate(label) for label in labels],
            linewidths=0.15,
            linecolor="white",
        )
        heat_axis.tick_params(axis="x", labelrotation=0, labelsize=5.1, length=0)
        heat_axis.tick_params(axis="y", labelrotation=0, labelsize=4.7, length=0)
        heat_axis.set_xlabel("neighborhood state", fontsize=5.4, labelpad=1)
        heat_axis.set_ylabel(
            "neighbor type" if column == 0 else "",
            fontsize=5.4,
            labelpad=2,
        )
        if column > 0:
            heat_axis.set_yticklabels([])
    color_axis = fig.add_axes([0.402, 0.652, 0.006, 0.095])
    fig.colorbar(
        heatmap_artist.collections[0],
        cax=color_axis,
        label="mean neighborhood fraction",
    )
    color_axis.tick_params(labelsize=4.7, length=2)
    color_axis.yaxis.label.set_size(5.2)


def draw_panel_b(
    axis: plt.Axes,
    labels: list[str],
    difference: np.ndarray,
    qvalues: np.ndarray,
) -> None:
    biological = [index for index, label in enumerate(labels) if label != "Abstained"]
    matrix = np.clip(difference[np.ix_(biological, biological)], 0.0, None)
    significance = qvalues[np.ix_(biological, biological)] < 0.05
    weighted = matrix * significance
    number = len(biological)
    angles = np.linspace(np.pi / 2, np.pi / 2 + 2 * np.pi, number, endpoint=False)
    positions = np.c_[np.cos(angles), np.sin(angles)]
    upper = np.triu(weighted + weighted.T, k=1)
    edge_values = upper[upper > 0]
    cutoff = np.quantile(edge_values, 0.25) if len(edge_values) else np.inf
    maximum = float(edge_values.max()) if len(edge_values) else 1.0
    for first in range(number):
        for second in range(first + 1, number):
            value = upper[first, second]
            if value <= cutoff:
                continue
            width = 0.35 + 3.2 * value / maximum
            axis.plot(
                [positions[first, 0], positions[second, 0]],
                [positions[first, 1], positions[second, 1]],
                color="#222222",
                linewidth=width,
                alpha=0.60,
                zorder=1,
            )
    for index, position in enumerate(positions):
        axis.scatter(
            position[0],
            position[1],
            s=125,
            color=NODE_COLORS[index % len(NODE_COLORS)],
            edgecolor="white",
            linewidth=1.0,
            zorder=3,
        )
        horizontal = "left" if position[0] >= 0 else "right"
        axis.text(
            position[0] * 1.18,
            position[1] * 1.18,
            abbreviate(labels[biological[index]]),
            ha=horizontal,
            va="center",
            fontsize=5.1,
        )
    axis.set_xlim(-1.45, 1.45)
    axis.set_ylim(-1.35, 1.35)
    axis.set_aspect("equal")
    axis.axis("off")


def draw_effect_heatmap(
    axis: plt.Axes,
    matrix: np.ndarray,
    labels: list[str],
    title: str,
    limit: float,
) -> None:
    shown = [index for index, label in enumerate(labels) if label != "Abstained"]
    values = matrix[np.ix_(shown, shown)]
    sns.heatmap(
        values,
        ax=axis,
        cmap="bwr",
        center=0,
        vmin=-limit,
        vmax=limit,
        cbar=False,
        xticklabels=[abbreviate(labels[index]) for index in shown],
        yticklabels=[abbreviate(labels[index]) for index in shown],
        linewidths=0.18,
        linecolor="white",
    )
    axis.set_title(title, fontsize=6.4, pad=2)
    axis.tick_params(axis="x", labelrotation=65, labelsize=4.4, length=0)
    axis.tick_params(axis="y", labelrotation=0, labelsize=4.4, length=0)
    axis.set_xlabel("receiver type", fontsize=5.2, labelpad=1)
    axis.set_ylabel("sender type", fontsize=5.2, labelpad=1)


def draw_volcano(
    axis: plt.Axes,
    labels: list[str],
    difference: np.ndarray,
    qvalues: np.ndarray,
) -> None:
    shown = [index for index, label in enumerate(labels) if label != "Abstained"]
    values = difference[np.ix_(shown, shown)].ravel()
    q = qvalues[np.ix_(shown, shown)].ravel()
    y = -np.log10(np.clip(q, 1e-6, 1.0))
    significant = q < 0.05
    colors = np.full(len(values), "#9A9A9A", dtype=object)
    colors[significant & (values > 0)] = "#E41A1C"
    colors[significant & (values < 0)] = "#377EB8"
    axis.scatter(
        values,
        y,
        s=14,
        c=colors,
        edgecolor="white",
        linewidth=0.25,
        alpha=0.88,
    )
    axis.axhline(-np.log10(0.05), color="#777777", linestyle="--", linewidth=0.7)
    axis.axvline(0, color="#777777", linestyle="--", linewidth=0.7)
    pair_labels = [
        f"{abbreviate(labels[source])}\u2192{abbreviate(labels[target])}"
        for source in shown
        for target in shown
    ]
    rank = np.argsort(np.where(significant, np.abs(values) * y, -1.0))[-2:]
    for offset, index in enumerate(rank):
        if not significant[index]:
            continue
        axis.annotate(
            pair_labels[index],
            (values[index], y[index]),
            xytext=(5, 5 + offset * 8),
            textcoords="offset points",
            fontsize=4.6,
            arrowprops={
                "arrowstyle": "-",
                "color": "#666666",
                "linewidth": 0.45,
            },
        )
    axis.set_xlabel("coupling change vs ideal abstention", fontsize=5.3)
    axis.set_ylabel(r"$-\log_{10}$ FDR", fontsize=5.3)
    axis.tick_params(labelsize=4.7)
    sns.despine(ax=axis)


def recovery_matrix(
    ideal: np.ndarray,
    baseline: np.ndarray,
    svtuner: np.ndarray,
) -> np.ndarray:
    baseline_error = np.abs(baseline - ideal)
    svtuner_error = np.abs(svtuner - ideal)
    pair_baseline_error = baseline_error + baseline_error.T
    pair_svtuner_error = svtuner_error + svtuner_error.T
    active = (
        np.abs(ideal)
        + np.abs(ideal.T)
        + np.abs(baseline)
        + np.abs(baseline.T)
        + np.abs(svtuner)
        + np.abs(svtuner.T)
    ) > 1e-7
    nonzero = np.concatenate(
        [
            pair_baseline_error[pair_baseline_error > 1e-10],
            pair_svtuner_error[pair_svtuner_error > 1e-10],
        ]
    )
    scale = float(np.quantile(nonzero, 0.55)) if len(nonzero) else 1.0
    scale = max(scale, 1e-6)
    baseline_score = np.exp(-pair_baseline_error / scale)
    svtuner_score = np.exp(-pair_svtuner_error / scale)
    result = np.full_like(ideal, np.nan, dtype=float)
    upper = np.triu_indices_from(result, k=1)
    lower = np.tril_indices_from(result, k=-1)
    result[upper] = baseline_score[upper]
    result[lower] = svtuner_score[lower]
    result[~active] = np.nan
    return result


def draw_recovery_heatmap(
    axis: plt.Axes,
    labels: list[str],
    ideal: np.ndarray,
    baseline: np.ndarray,
    svtuner: np.ndarray,
) -> None:
    shown = [index for index, label in enumerate(labels) if label != "Abstained"]
    matrix = recovery_matrix(
        ideal[np.ix_(shown, shown)],
        baseline[np.ix_(shown, shown)],
        svtuner[np.ix_(shown, shown)],
    )
    sns.heatmap(
        matrix,
        ax=axis,
        cmap="Purples",
        vmin=0,
        vmax=1,
        cbar=True,
        cbar_kws={"label": "recovery", "shrink": 0.62, "pad": 0.02},
        xticklabels=[abbreviate(labels[index]) for index in shown],
        yticklabels=[abbreviate(labels[index]) for index in shown],
        square=True,
        linewidths=0.15,
        linecolor="white",
        mask=np.isnan(matrix),
    )
    axis.tick_params(axis="x", labelrotation=65, labelsize=4.3, length=0)
    axis.tick_params(axis="y", labelrotation=0, labelsize=4.3, length=0)
    axis.set_xlabel("upper: CytoSPACE", fontsize=5.0, labelpad=1)
    axis.set_ylabel("lower: SVTuner", fontsize=5.0, labelpad=1)
    colorbar = axis.collections[0].colorbar
    colorbar.ax.tick_params(labelsize=4.3, length=2)
    colorbar.ax.yaxis.label.set_size(4.8)


def save_tables(
    experiment: dict[str, object],
    out_dir: Path,
    file_stem: str,
) -> None:
    labels = list(experiment["labels"])
    contacts = experiment["contacts"]
    rows: list[dict[str, object]] = []
    for condition, matrix in contacts.items():
        for source_index, source_type in enumerate(labels):
            for receiver_index, receiver_type in enumerate(labels):
                rows.append(
                    {
                        "condition": condition,
                        "source_type": source_type,
                        "receiver_type": receiver_type,
                        "coupling": matrix[source_index, receiver_index],
                    }
                )
    pd.DataFrame(rows).to_csv(
        out_dir / f"{file_stem}_coupling_matrices.csv",
        index=False,
    )

    significance_rows: list[dict[str, object]] = []
    for method in ("CytoSPACE dropout", "SVTuner"):
        difference = experiment["differences"][method]
        qvalues = experiment["qvalues"][method]
        standard_errors = experiment["standard_errors"][method]
        for source_index, source_type in enumerate(labels):
            for receiver_index, receiver_type in enumerate(labels):
                significance_rows.append(
                    {
                        "method": method,
                        "source_type": source_type,
                        "receiver_type": receiver_type,
                        "coupling_change": difference[source_index, receiver_index],
                        "bootstrap_se": standard_errors[
                            source_index, receiver_index
                        ],
                        "fdr": qvalues[source_index, receiver_index],
                    }
                )
    pd.DataFrame(significance_rows).to_csv(
        out_dir / f"{file_stem}_pair_statistics.csv",
        index=False,
    )


def plot_figure(
    experiment: dict[str, object],
    out_dir: Path,
    file_stem: str,
) -> None:
    sns.set_theme(style="white", context="paper")
    fig = plt.figure(figsize=(11.4, 7.8), dpi=350, facecolor="white")
    outer = fig.add_gridspec(
        2,
        12,
        height_ratios=[1.12, 0.88],
        hspace=0.34,
        wspace=0.85,
        left=0.055,
        right=0.975,
        top=0.95,
        bottom=0.11,
    )
    panel_a_spec = outer[0, :8]
    panel_b_axis = fig.add_subplot(outer[0, 8:])
    panel_c_axis = fig.add_subplot(outer[1, 0:3])
    panel_d_axis = fig.add_subplot(outer[1, 3:6])
    panel_e_axis = fig.add_subplot(outer[1, 6:9])
    panel_f_axis = fig.add_subplot(outer[1, 9:12])

    draw_panel_a(fig, panel_a_spec, experiment)
    labels = list(experiment["labels"])
    differences = experiment["differences"]
    qvalues = experiment["qvalues"]
    draw_panel_b(
        panel_b_axis,
        labels,
        differences["CytoSPACE dropout"],
        qvalues["CytoSPACE dropout"],
    )
    biological = [index for index, label in enumerate(labels) if label != "Abstained"]
    maximum = max(
        float(
            np.abs(differences[method][np.ix_(biological, biological)]).max()
        )
        for method in ("CytoSPACE dropout", "SVTuner")
    )
    draw_effect_heatmap(
        panel_c_axis,
        differences["CytoSPACE dropout"],
        labels,
        "",
        maximum,
    )
    draw_effect_heatmap(
        panel_d_axis,
        differences["SVTuner"],
        labels,
        "",
        maximum,
    )
    draw_volcano(
        panel_e_axis,
        labels,
        differences["CytoSPACE dropout"],
        qvalues["CytoSPACE dropout"],
    )
    contacts = experiment["contacts"]
    draw_recovery_heatmap(
        panel_f_axis,
        labels,
        contacts["Ideal abstention"],
        contacts["CytoSPACE dropout"],
        contacts["SVTuner"],
    )

    fig.text(
        0.052,
        0.978,
        "A",
        fontsize=9,
        fontweight="bold",
        va="top",
    )
    fig.text(
        0.071,
        0.978,
        "Neighborhood states and spatial composition",
        fontsize=7.2,
        fontweight="bold",
        va="top",
    )
    panel_titles = (
        (panel_b_axis, "B", "False cell-type coupling"),
        (panel_c_axis, "C", "Forced sender-receiver effects"),
        (panel_d_axis, "D", "Effects after Stage3B abstention"),
        (panel_e_axis, "E", "Spurious coupling significance"),
        (panel_f_axis, "F", "Coupling recovery"),
    )
    for axis, letter, title in panel_titles:
        axis.text(
            -0.12,
            1.055,
            letter,
            transform=axis.transAxes,
            fontsize=8,
            fontweight="bold",
            va="top",
        )
        axis.text(
            -0.01,
            1.055,
            title,
            transform=axis.transAxes,
            fontsize=6.0,
            fontweight="bold",
            va="top",
        )

    state_legend = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=STATE_COLORS[index],
            markeredgecolor="none",
            markersize=4,
            label=f"state {index + 1}",
        )
        for index in range(len(np.unique(experiment["state_labels"])))
    ]
    fig.legend(
        handles=state_legend,
        loc="upper left",
        bbox_to_anchor=(0.058, 0.905),
        frameon=False,
        ncol=len(state_legend),
        fontsize=4.8,
        handletextpad=0.2,
        columnspacing=0.65,
    )
    source_label = str(experiment["row"]["pair_id"]).replace("_", " ")
    target_label = str(experiment["row"]["target_type"])
    fig.text(
        0.055,
        0.032,
        (
            f"{source_label}, {target_label} reference dropout. Evaluation markers define "
            "the top-15% ST target region only; they are not used by Stage3B. "
            "Ideal abstention credits blanking inside that region."
        ),
        fontsize=5.3,
        color="#444444",
    )
    png = out_dir / f"{file_stem}_ncem_style.png"
    pdf = out_dir / f"{file_stem}_ncem_style.pdf"
    fig.savefig(png, dpi=350, facecolor="white")
    fig.savefig(pdf, facecolor="white")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    root = Path(args.project_root).resolve()
    out_dir = root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    metadata = pd.read_csv(root / args.metadata_csv)
    experiment = prepare_experiment(
        root=root,
        metadata=metadata,
        sample=args.sample,
        target_quantile=args.target_quantile,
        neighbors=args.neighbors,
        states=args.states,
        repeats=args.bootstrap,
        seed=args.seed,
    )
    file_stem = str(experiment["row"]["sample"])
    save_tables(experiment, out_dir, file_stem)
    plot_figure(experiment, out_dir, file_stem)
    target_core = np.asarray(experiment["target_core"], dtype=bool)
    svtuner = np.asarray(experiment["augmented"]["SVTuner"])
    blank_index = list(experiment["labels"]).index("Abstained")
    blank = svtuner[:, blank_index] > 0
    baseline_difference = np.asarray(
        experiment["differences"]["CytoSPACE dropout"]
    )
    svtuner_difference = np.asarray(experiment["differences"]["SVTuner"])
    baseline_qvalue = np.asarray(
        experiment["qvalues"]["CytoSPACE dropout"]
    )
    svtuner_qvalue = np.asarray(experiment["qvalues"]["SVTuner"])
    biological = np.arange(blank_index)
    biological_slice = np.ix_(biological, biological)
    baseline_positive = baseline_difference[biological_slice]
    baseline_positive = baseline_positive[baseline_positive > 0]
    material_threshold = float(
        np.quantile(baseline_positive, 0.90)
    ) if len(baseline_positive) else 0.0
    summary = pd.DataFrame(
        [
            {
                "sample": experiment["row"]["sample"],
                "target_type": experiment["row"]["target_type"],
                "spots": len(target_core),
                "marker_target_spots": int(target_core.sum()),
                "stage3b_blank_spots": int(blank.sum()),
                "target_blank_overlap": int((target_core & blank).sum()),
                "precision": float((target_core & blank).sum() / max(blank.sum(), 1)),
                "recall": float(
                    (target_core & blank).sum() / max(target_core.sum(), 1)
                ),
                "cytospace_absolute_coupling_error": float(
                    np.abs(baseline_difference[biological_slice]).sum()
                ),
                "svtuner_absolute_coupling_error": float(
                    np.abs(svtuner_difference[biological_slice]).sum()
                ),
                "material_effect_threshold": material_threshold,
                "cytospace_material_spurious_pairs": int(
                    (
                        (baseline_qvalue[biological_slice] < 0.05)
                        & (
                            baseline_difference[biological_slice]
                            >= material_threshold
                        )
                    ).sum()
                ),
                "svtuner_material_spurious_pairs": int(
                    (
                        (svtuner_qvalue[biological_slice] < 0.05)
                        & (
                            svtuner_difference[biological_slice]
                            >= material_threshold
                        )
                    ).sum()
                ),
                "uses_target_markers_in_stage3b": False,
            }
        ]
    )
    summary.to_csv(
        out_dir / f"{file_stem}_summary.csv",
        index=False,
    )
    print(summary.to_string(index=False))
    print(
        "[done]",
        out_dir / f"{file_stem}_ncem_style.png",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
