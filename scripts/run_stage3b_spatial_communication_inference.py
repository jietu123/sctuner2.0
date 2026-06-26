from __future__ import annotations

import json
import math
import tarfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats


ROOT = Path(__file__).resolve().parents[1]
SAMPLE = "cytospace_fig2d_tme_brca_her2_ffpe_sc_missing_plasma_cells"
OUT = ROOT / "data/processed/stage3b_spatial_communication/brca_her2_ffpe_plasma"
CACHE = OUT / "cache"

RAW_SC = (
    ROOT
    / "data/raw/cytospace_fig2d_tme/brca_scRNA_GSE176078"
    / "Wu_etal_2021_BRCA_scRNASeq"
)
RAW_ST = (
    ROOT
    / "data/raw/cytospace_fig2d_tme/brca_ST_10x_ffpe"
    / "Visium_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix (1).h5"
)
RAW_SPATIAL = (
    ROOT
    / "data/raw/cytospace_fig2d_tme/brca_ST_10x_ffpe"
    / "Visium_FFPE_Human_Breast_Cancer_spatial.tar.gz"
)
STAGE1 = (
    ROOT
    / "data/processed/cytospace_fig2d_tme"
    / "cytospace_fig2d_tme_brca_her2_ffpe/stage1_preprocess/exported"
)
LR_PATH = (
    ROOT
    / "data/raw/spatial_communication_reference_resources"
    / "liana_1.7.3_omni_resource.csv"
)

CONDITIONS = {
    "full_reference": (
        ROOT
        / "result/cytospace_fig2d_tme_brca_her2_ffpe"
        / "stage4_cytospace_niche_full_reference/cytospace_output"
    ),
    "cytospace_dropout": (
        ROOT
        / f"result/{SAMPLE}"
        / "stage4_cytospace_niche_dropout_baseline/cytospace_output"
    ),
    "svtuner_stage3b": (
        ROOT
        / f"result/{SAMPLE}"
        / "stage4_cytospace_niche_stage3b_blank/cytospace_output"
    ),
}


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    values = np.asarray(pvalues, dtype=float)
    result = np.full(values.shape, np.nan, dtype=float)
    finite = np.isfinite(values)
    if not finite.any():
        return result
    p = values[finite]
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    inverse = np.empty_like(order)
    inverse[order] = np.arange(len(order))
    result[finite] = np.minimum(adjusted[inverse], 1.0)
    return result


def components(entity: str) -> list[str]:
    return [part.strip().upper() for part in str(entity).split("_") if part.strip()]


def load_lr(common_genes: set[str]) -> pd.DataFrame:
    resource = pd.read_csv(LR_PATH)
    resource = resource.loc[
        resource["resource"].eq("consensus"),
        ["source_genesymbol", "target_genesymbol"],
    ].drop_duplicates()
    resource.columns = ["ligand", "receptor"]
    keep = resource["ligand"].map(
        lambda value: all(gene in common_genes for gene in components(value))
    ) & resource["receptor"].map(
        lambda value: all(gene in common_genes for gene in components(value))
    )
    return resource.loc[keep].reset_index(drop=True)


def matrix_market_header_lines(path: Path) -> int:
    lines = 0
    with path.open("rt", encoding="utf-8") as handle:
        for line in handle:
            lines += 1
            if not line.startswith("%"):
                break
    return lines


def extract_selected_sc_counts(selected_cells: list[str]) -> sparse.csr_matrix:
    matrix_cache = CACHE / "selected_sc_raw_counts.npz"
    cells_cache = CACHE / "selected_sc_cells.txt"
    genes_cache = CACHE / "selected_sc_genes.txt"
    if matrix_cache.exists() and cells_cache.exists() and genes_cache.exists():
        cached_cells = cells_cache.read_text(encoding="utf-8").splitlines()
        if cached_cells == selected_cells:
            return sparse.load_npz(matrix_cache).tocsr()

    print("[communication] extracting selected cells from raw SC matrix", flush=True)
    all_cells = (
        RAW_SC / "count_matrix_barcodes.tsv"
    ).read_text(encoding="utf-8").splitlines()
    all_genes = (
        RAW_SC / "count_matrix_genes.tsv"
    ).read_text(encoding="utf-8").splitlines()
    cell_index = {cell: index for index, cell in enumerate(all_cells)}
    missing = [cell for cell in selected_cells if cell not in cell_index]
    if missing:
        raise ValueError(f"{len(missing)} selected cells are absent from raw SC matrix")

    selected_raw_columns = np.array([cell_index[cell] for cell in selected_cells])
    raw_to_selected = np.full(len(all_cells), -1, dtype=np.int32)
    raw_to_selected[selected_raw_columns] = np.arange(
        len(selected_cells), dtype=np.int32
    )
    selected_lookup = raw_to_selected >= 0

    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    data: list[np.ndarray] = []
    path = RAW_SC / "count_matrix_sparse.mtx"
    skiprows = matrix_market_header_lines(path)
    reader = pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        names=["gene", "cell", "count"],
        skiprows=skiprows,
        dtype={"gene": np.int32, "cell": np.int32, "count": np.int32},
        chunksize=2_000_000,
        engine="c",
    )
    processed = 0
    retained = 0
    for chunk in reader:
        gene = chunk["gene"].to_numpy() - 1
        cell = chunk["cell"].to_numpy() - 1
        mask = selected_lookup[cell]
        if mask.any():
            rows.append(gene[mask].astype(np.int32, copy=False))
            cols.append(raw_to_selected[cell[mask]])
            data.append(chunk["count"].to_numpy()[mask].astype(np.float32))
            retained += int(mask.sum())
        processed += len(chunk)
        if processed % 20_000_000 < len(chunk):
            print(
                f"[communication] scanned {processed:,} entries; "
                f"retained {retained:,}",
                flush=True,
            )

    matrix = sparse.coo_matrix(
        (np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
        shape=(len(all_genes), len(selected_cells)),
        dtype=np.float32,
    ).tocsr()
    matrix.sum_duplicates()
    sparse.save_npz(matrix_cache, matrix, compressed=True)
    cells_cache.write_text("\n".join(selected_cells) + "\n", encoding="utf-8")
    genes_cache.write_text(
        "\n".join(gene.upper() for gene in all_genes) + "\n", encoding="utf-8"
    )
    print(
        f"[communication] cached selected SC counts: {matrix.shape}, "
        f"{matrix.nnz:,} nonzero entries",
        flush=True,
    )
    return matrix


def log_normalize_columns(matrix: sparse.csr_matrix) -> sparse.csc_matrix:
    result = matrix.astype(np.float32).tocsc(copy=True)
    totals = np.asarray(result.sum(axis=0)).ravel()
    totals[totals <= 0] = 1.0
    result = result @ sparse.diags(
        (10000.0 / totals).astype(np.float32), format="csc"
    )
    result.data = np.log1p(result.data)
    return result.tocsc()


def read_st() -> tuple[sparse.csr_matrix, list[str], list[str]]:
    with h5py.File(RAW_ST, "r") as handle:
        group = handle["matrix"]
        matrix = sparse.csc_matrix(
            (group["data"][:], group["indices"][:], group["indptr"][:]),
            shape=tuple(int(value) for value in group["shape"][:]),
        ).tocsr()
        genes = [
            value.decode() if isinstance(value, bytes) else str(value)
            for value in group["features"]["name"][:]
        ]
        spots = [
            value.decode() if isinstance(value, bytes) else str(value)
            for value in group["barcodes"][:]
        ]
    return matrix, [gene.upper() for gene in genes], spots


def read_visium_positions() -> pd.DataFrame:
    with tarfile.open(RAW_SPATIAL, "r:gz") as archive:
        member = next(
            member
            for member in archive.getmembers()
            if member.name.endswith("tissue_positions_list.csv")
        )
        handle = archive.extractfile(member)
        if handle is None:
            raise FileNotFoundError(member.name)
        positions = pd.read_csv(handle, header=None)
    positions = positions.iloc[:, :6]
    positions.columns = [
        "spot_id",
        "in_tissue",
        "array_row",
        "array_col",
        "pixel_row",
        "pixel_col",
    ]
    positions["spot_id"] = positions["spot_id"].astype(str)
    return positions.loc[positions["in_tissue"].eq(1)].copy()


def build_visium_graph(spots: list[str]) -> tuple[sparse.csr_matrix, pd.DataFrame]:
    positions = read_visium_positions().set_index("spot_id").loc[spots].reset_index()
    location_to_index = {
        (int(row.array_row), int(row.array_col)): index
        for index, row in positions.iterrows()
    }
    undirected_offsets = [(0, 2), (1, 1), (1, -1)]
    rows = []
    cols = []
    for index, row in positions.iterrows():
        location = (int(row.array_row), int(row.array_col))
        for dr, dc in undirected_offsets:
            neighbor = location_to_index.get((location[0] + dr, location[1] + dc))
            if neighbor is not None:
                rows.extend([index, neighbor])
                cols.extend([neighbor, index])
    adjacency = sparse.coo_matrix(
        (np.ones(len(rows), dtype=np.float32), (rows, cols)),
        shape=(len(spots), len(spots)),
    ).tocsr()
    adjacency.sum_duplicates()
    adjacency.data[:] = 1.0
    return adjacency, positions


def normalized_weight(
    adjacency: sparse.csr_matrix, valid: np.ndarray
) -> sparse.csr_matrix:
    mask = valid.astype(np.float32)
    weight = adjacency.multiply(mask[:, None]).multiply(mask[None, :]).tocsr()
    total = float(weight.sum())
    n_valid = int(valid.sum())
    if total <= 0 or n_valid <= 1:
        raise ValueError("Spatial graph has no valid edges")
    return weight * (n_valid / total)


def entity_matrix(
    spot_gene_expression: np.ndarray,
    gene_to_index: dict[str, int],
    entities: list[str],
) -> np.ndarray:
    output = np.zeros((spot_gene_expression.shape[0], len(entities)), dtype=np.float32)
    for index, entity in enumerate(entities):
        columns = [gene_to_index[gene] for gene in components(entity)]
        if len(columns) == 1:
            output[:, index] = spot_gene_expression[:, columns[0]]
        else:
            # A receptor complex is available only to the extent that all
            # required subunits are present. This matches the conservative
            # complex handling used by communication frameworks.
            output[:, index] = np.min(spot_gene_expression[:, columns], axis=1)
    return output


def standardize(matrix: np.ndarray, valid: np.ndarray) -> np.ndarray:
    output = np.zeros_like(matrix, dtype=np.float32)
    values = matrix[valid].astype(np.float64)
    mean = values.mean(axis=0, keepdims=True)
    std = values.std(axis=0, keepdims=True)
    nonconstant = std.ravel() > 0
    if nonconstant.any():
        output[:, nonconstant] = (
            (matrix[:, nonconstant] - mean[:, nonconstant])
            / std[:, nonconstant]
        ).astype(np.float32)
    output[~valid] = 0.0
    return output


def moran_std(weight: sparse.csr_matrix) -> float:
    n = weight.shape[0]
    numerator = (
        n**2 * weight.multiply(weight.T).sum()
        - 2 * n * (weight.sum(axis=0) @ weight.sum(axis=1)).item()
        + float(weight.sum()) ** 2
    )
    denominator = n**2 * (n - 1) ** 2
    return math.sqrt(max(float(numerator / denominator), 0.0))


def global_moran_statistics(
    ligand_values: np.ndarray,
    receptor_values: np.ndarray,
    pair_ligand: np.ndarray,
    pair_receptor: np.ndarray,
    adjacency: sparse.csr_matrix,
    valid: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, sparse.csr_matrix]:
    weight = normalized_weight(adjacency, valid)
    ligand_z = standardize(ligand_values, valid)
    receptor_z = standardize(receptor_values, valid)
    weighted_ligand = weight @ ligand_z
    weighted_receptor = weight @ receptor_z
    global_moran = (
        weighted_receptor[:, pair_receptor] * ligand_z[:, pair_ligand]
    ).sum(axis=0) / float(weight.sum())
    std = moran_std(weight)
    z_score = global_moran / std if std > 0 else np.zeros(len(pair_ligand))
    p_value = stats.norm.sf(z_score)
    local_moran = 0.5 * (
        weighted_receptor[:, pair_receptor] * ligand_z[:, pair_ligand]
        + weighted_ligand[:, pair_ligand] * receptor_z[:, pair_receptor]
    )
    local_moran[~valid] = np.nan
    return global_moran, z_score, p_value, local_moran, weight


def infer_condition(
    condition: str,
    output: Path,
    selected_cells: list[str],
    selected_sc_norm: sparse.csc_matrix,
    sc_gene_to_row: dict[str, int],
    lr: pd.DataFrame,
    spots: list[str],
    adjacency: sparse.csr_matrix,
    blank_spots: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame | None, np.ndarray]:
    print(f"[communication] reconstructing {condition}", flush=True)
    assignment = pd.read_csv(output / "cell_assignment.csv")
    cell_to_index = {cell: index for index, cell in enumerate(selected_cells)}
    spot_to_index = {spot: index for index, spot in enumerate(spots)}
    missing_cells = set(assignment["cell_id"].astype(str)) - set(cell_to_index)
    if missing_cells:
        raise ValueError(f"{condition}: {len(missing_cells)} assigned cells missing")

    row = assignment["assigned_spot"].map(spot_to_index).to_numpy()
    col = assignment["cell_id"].map(cell_to_index).to_numpy()
    mapping = sparse.coo_matrix(
        (np.ones(len(assignment), dtype=np.float32), (row, col)),
        shape=(len(spots), len(selected_cells)),
    ).tocsr()
    counts = np.asarray(mapping.sum(axis=1)).ravel()

    lr_genes = sorted(
        {
            gene
            for entity in pd.concat([lr["ligand"], lr["receptor"]])
            for gene in components(entity)
        }
    )
    sc_rows = [sc_gene_to_row[gene] for gene in lr_genes]
    spot_gene = (mapping @ selected_sc_norm[sc_rows, :].T).toarray()
    nonempty = counts > 0
    spot_gene[nonempty] /= counts[nonempty, None]
    spot_gene[~nonempty] = np.nan
    np.savez_compressed(
        OUT / f"{condition}_spot_lr_expression.npz",
        expression=spot_gene.astype(np.float32),
        genes=np.asarray(lr_genes),
        spots=np.asarray(spots),
        mapped=np.asarray(nonempty),
    )

    gene_to_index = {gene: index for index, gene in enumerate(lr_genes)}
    ligands = lr["ligand"].drop_duplicates().tolist()
    receptors = lr["receptor"].drop_duplicates().tolist()
    ligand_values = entity_matrix(spot_gene, gene_to_index, ligands)
    receptor_values = entity_matrix(spot_gene, gene_to_index, receptors)
    ligand_index = {entity: index for index, entity in enumerate(ligands)}
    receptor_index = {entity: index for index, entity in enumerate(receptors)}
    pair_ligand = lr["ligand"].map(ligand_index).to_numpy()
    pair_receptor = lr["receptor"].map(receptor_index).to_numpy()

    blank_mask = np.array([spot in blank_spots for spot in spots])
    graph_connected = np.asarray(adjacency.sum(axis=1)).ravel() > 0
    supported = (~blank_mask) & graph_connected
    if condition == "svtuner_stage3b":
        supported &= nonempty
    global_moran, z_score, p_value, local_moran, weight = global_moran_statistics(
        ligand_values,
        receptor_values,
        pair_ligand,
        pair_receptor,
        adjacency,
        supported,
    )

    supported_adjacency = adjacency.multiply(supported[None, :]).tocsr()
    degree = np.asarray(supported_adjacency.sum(axis=1)).ravel()
    degree_safe = degree.copy()
    degree_safe[degree_safe <= 0] = 1.0
    neighbor_receptor = (
        supported_adjacency @ np.nan_to_num(receptor_values)
    ) / degree_safe[:, None]
    intensity = np.sqrt(
        np.maximum(
            ligand_values[:, pair_ligand]
            * neighbor_receptor[:, pair_receptor],
            0.0,
        )
    )
    intensity[~supported] = np.nan
    mean_intensity = np.nanmean(intensity, axis=0)

    np.savez_compressed(
        OUT / f"{condition}_local_lr_statistics.npz",
        local_moran=local_moran.astype(np.float32),
        spots=np.asarray(spots),
        ligand=np.asarray(lr["ligand"]),
        receptor=np.asarray(lr["receptor"]),
        supported=supported,
    )

    result = lr.copy()
    result["condition"] = condition
    result["moran_r"] = global_moran
    result["z_score"] = z_score
    result["p_value"] = p_value
    result["q_value"] = bh_fdr(p_value)
    result["mean_communication_intensity"] = mean_intensity
    result["n_supported_spots"] = int(supported.sum())
    result["n_supported_directed_edges"] = int(weight.nnz)
    result["ligand_nonzero_spots"] = (
        ligand_values[:, pair_ligand][supported] > 0
    ).sum(axis=0)
    result["receptor_nonzero_spots"] = (
        receptor_values[:, pair_receptor][supported] > 0
    ).sum(axis=0)

    whole_result = None
    full_degree = np.asarray(adjacency.sum(axis=1)).ravel()
    full_degree[full_degree <= 0] = 1.0
    neighbor_ligand_full = (
        adjacency @ np.nan_to_num(ligand_values)
    ) / full_degree[:, None]
    neighbor_receptor_full = (
        adjacency @ np.nan_to_num(receptor_values)
    ) / full_degree[:, None]
    target_intensity = 0.5 * (
        np.sqrt(
            np.maximum(
                receptor_values[:, pair_receptor]
                * neighbor_ligand_full[:, pair_ligand],
                0.0,
            )
        )
        + np.sqrt(
            np.maximum(
                ligand_values[:, pair_ligand]
                * neighbor_receptor_full[:, pair_receptor],
                0.0,
            )
        )
    )
    if condition == "svtuner_stage3b":
        target_intensity[:] = np.nan
    else:
        whole_valid = nonempty & graph_connected
        whole_moran, whole_z, whole_p, whole_local, whole_weight = (
            global_moran_statistics(
                ligand_values,
                receptor_values,
                pair_ligand,
                pair_receptor,
                adjacency,
                whole_valid,
            )
        )
        np.savez_compressed(
            OUT / f"{condition}_whole_tissue_local_lr_statistics.npz",
            local_moran=whole_local.astype(np.float32),
            spots=np.asarray(spots),
            ligand=np.asarray(lr["ligand"]),
            receptor=np.asarray(lr["receptor"]),
            supported=whole_valid,
        )
        whole_result = lr.copy()
        whole_result["condition"] = condition
        whole_result["moran_r"] = whole_moran
        whole_result["z_score"] = whole_z
        whole_result["p_value"] = whole_p
        whole_result["q_value"] = bh_fdr(whole_p)
        whole_result["n_spots"] = int(whole_valid.sum())
        whole_result["n_directed_edges"] = int(whole_weight.nnz)

    np.savez_compressed(
        OUT / f"{condition}_target_region_lr_intensity.npz",
        intensity=target_intensity[blank_mask].astype(np.float32),
        spots=np.asarray(spots)[blank_mask],
        ligand=np.asarray(lr["ligand"]),
        receptor=np.asarray(lr["receptor"]),
    )
    return result, whole_result, target_intensity[blank_mask]


def exact_spatial_block_permutation(
    full_target: np.ndarray,
    dropout_target: np.ndarray,
    target_spots: list[str],
    blank_manifest: pd.DataFrame,
    lr: pd.DataFrame,
) -> pd.DataFrame:
    manifest = blank_manifest.set_index("spot_id").loc[target_spots]
    regions = sorted(manifest["region_id"].unique())
    difference = dropout_target - full_target
    region_sums = np.vstack(
        [
            np.nansum(
                difference[manifest["region_id"].to_numpy() == region],
                axis=0,
            )
            for region in regions
        ]
    )
    n_regions = len(regions)
    signs = np.array(
        [
            [1.0 if (mask >> bit) & 1 else -1.0 for bit in range(n_regions)]
            for mask in range(2**n_regions)
        ],
        dtype=np.float32,
    )
    null_means = signs @ region_sums / len(target_spots)
    observed = np.nanmean(difference, axis=0)
    exceedances = np.sum(null_means >= observed[None, :], axis=0)
    p_value = (exceedances + 1.0) / (len(null_means) + 1.0)

    rng = np.random.default_rng(42)
    bootstrap = np.empty((1000, difference.shape[1]), dtype=np.float32)
    for index in range(1000):
        sampled_regions = rng.choice(regions, size=n_regions, replace=True)
        sampled_rows = np.concatenate(
            [
                np.flatnonzero(manifest["region_id"].to_numpy() == region)
                for region in sampled_regions
            ]
        )
        bootstrap[index] = np.nanmean(difference[sampled_rows], axis=0)

    result = lr.copy()
    result["full_reference_target_intensity"] = np.nanmean(full_target, axis=0)
    result["cytospace_dropout_target_intensity"] = np.nanmean(
        dropout_target, axis=0
    )
    result["dropout_induced_target_change"] = observed
    result["spatial_block_permutation_pvalue"] = p_value
    result["spatial_block_permutation_qvalue"] = bh_fdr(p_value)
    result["region_bootstrap_ci_low"] = np.quantile(bootstrap, 0.025, axis=0)
    result["region_bootstrap_ci_high"] = np.quantile(bootstrap, 0.975, axis=0)
    result["n_spatial_regions"] = n_regions
    result["n_target_spots"] = len(target_spots)
    return result


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)

    sc_meta = pd.read_csv(STAGE1 / "sc_metadata.csv")
    selected_cells = sc_meta["cell_id"].astype(str).tolist()
    sc_counts = extract_selected_sc_counts(selected_cells)
    sc_genes = (
        CACHE / "selected_sc_genes.txt"
    ).read_text(encoding="utf-8").splitlines()
    sc_gene_to_row = {gene: index for index, gene in enumerate(sc_genes)}
    sc_norm = log_normalize_columns(sc_counts)

    st_counts, st_genes, st_spots = read_st()
    stage1_coordinates = pd.read_csv(STAGE1 / "st_coordinates.csv")
    spots = stage1_coordinates["spot_id"].astype(str).tolist()
    st_spot_to_col = {spot: index for index, spot in enumerate(st_spots)}
    st_columns = [st_spot_to_col[spot] for spot in spots]
    st_counts = st_counts[:, st_columns]
    common_genes = set(sc_genes) & set(st_genes)
    lr = load_lr(common_genes)
    lr.to_csv(OUT / "liana_consensus_lr_pairs_used.csv", index=False)

    adjacency, positions = build_visium_graph(spots)
    sparse.save_npz(OUT / "visium_six_neighbor_adjacency.npz", adjacency)
    positions.to_csv(OUT / "visium_spot_positions.csv", index=False)

    blank_manifest = pd.read_csv(
        CONDITIONS["svtuner_stage3b"] / "stage3b_blank_spots.csv"
    )
    blank_spots = set(blank_manifest["spot_id"].astype(str))

    results = []
    whole_results = []
    target_results = {}
    for condition, output in CONDITIONS.items():
        supported_result, whole_result, target_intensity = infer_condition(
            condition=condition,
            output=output,
            selected_cells=selected_cells,
            selected_sc_norm=sc_norm,
            sc_gene_to_row=sc_gene_to_row,
            lr=lr,
            spots=spots,
            adjacency=adjacency,
            blank_spots=blank_spots,
        )
        results.append(supported_result)
        target_results[condition] = target_intensity
        if whole_result is not None:
            whole_results.append(whole_result)
    combined = pd.concat(results, ignore_index=True)
    combined.to_csv(OUT / "global_lr_communication_statistics.csv", index=False)
    pd.concat(whole_results, ignore_index=True).to_csv(
        OUT / "whole_tissue_lr_communication_statistics.csv", index=False
    )

    target_spots = [spot for spot in spots if spot in blank_spots]
    target_test = exact_spatial_block_permutation(
        target_results["full_reference"],
        target_results["cytospace_dropout"],
        target_spots,
        blank_manifest,
        lr,
    )
    target_test.to_csv(
        OUT / "target_region_dropout_induced_lr_statistics.csv", index=False
    )

    wide = combined.pivot(
        index=["ligand", "receptor"],
        columns="condition",
        values=["moran_r", "mean_communication_intensity", "q_value"],
    )
    wide.columns = ["__".join(column) for column in wide.columns]
    wide = wide.reset_index()
    wide["dropout_induced_moran_change"] = (
        wide["moran_r__cytospace_dropout"] - wide["moran_r__full_reference"]
    )
    wide["stage3b_supported_moran_change_vs_dropout"] = (
        wide["moran_r__svtuner_stage3b"] - wide["moran_r__cytospace_dropout"]
    )
    wide["dropout_induced_intensity_change"] = (
        wide["mean_communication_intensity__cytospace_dropout"]
        - wide["mean_communication_intensity__full_reference"]
    )
    wide.to_csv(OUT / "three_condition_lr_comparison.csv", index=False)

    whole_wide = pd.concat(whole_results, ignore_index=True).pivot(
        index=["ligand", "receptor"],
        columns="condition",
        values=["moran_r", "q_value"],
    )
    whole_wide.columns = ["__".join(column) for column in whole_wide.columns]
    evidence = (
        target_test.merge(
            whole_wide.reset_index(), on=["ligand", "receptor"], how="left"
        )
        .merge(
            wide[
                [
                    "ligand",
                    "receptor",
                    "moran_r__full_reference",
                    "moran_r__cytospace_dropout",
                    "moran_r__svtuner_stage3b",
                ]
            ],
            on=["ligand", "receptor"],
            how="left",
            suffixes=("_whole", "_supported"),
        )
    )
    evidence["dropout_only_whole_tissue_spatial_signal"] = (
        evidence["q_value__cytospace_dropout"].lt(0.05)
        & evidence["q_value__full_reference"].ge(0.05)
        & evidence["moran_r__cytospace_dropout_whole"].gt(
            evidence["moran_r__full_reference_whole"]
        )
    )
    evidence["positive_target_region_effect"] = (
        evidence["dropout_induced_target_change"].gt(0)
        & evidence["region_bootstrap_ci_low"].gt(0)
    )
    evidence["spatial_blocks_consistent"] = evidence[
        "spatial_block_permutation_pvalue"
    ].le(0.05)
    evidence["supported_region_recovery"] = (
        (
            evidence["moran_r__svtuner_stage3b"]
            - evidence["moran_r__full_reference_supported"]
        ).abs()
        <= (
            evidence["moran_r__cytospace_dropout_supported"]
            - evidence["moran_r__full_reference_supported"]
        ).abs()
    )
    evidence["high_confidence_dropout_induced_false_communication"] = (
        evidence["dropout_only_whole_tissue_spatial_signal"]
        & evidence["positive_target_region_effect"]
        & evidence["spatial_blocks_consistent"]
    )
    evidence.to_csv(OUT / "false_communication_evidence_table.csv", index=False)

    graph_degree = np.asarray(adjacency.sum(axis=1)).ravel()
    summary = {
        "sample": SAMPLE,
        "selected_sc_cells": len(selected_cells),
        "raw_sc_genes": len(sc_genes),
        "raw_st_genes": len(st_genes),
        "full_common_genes": len(common_genes),
        "liana_consensus_lr_pairs": len(lr),
        "spots": len(spots),
        "stage3b_blank_spots": len(blank_spots),
        "stage3b_blank_regions": int(blank_manifest["region_id"].nunique()),
        "visium_undirected_edges": int(adjacency.nnz // 2),
        "visium_degree_min": int(graph_degree.min()),
        "visium_degree_median": float(np.median(graph_degree)),
        "visium_degree_max": int(graph_degree.max()),
        "visium_isolated_spots_excluded_from_spatial_statistics": int(
            (graph_degree == 0).sum()
        ),
        "normalization": "per-cell library size 10000 followed by log1p",
        "complex_rule": "minimum expression across required subunits",
        "spatial_statistic": "bivariate Moran R with analytical one-sided p-value",
        "multiple_testing": "Benjamini-Hochberg within each condition",
        "target_region_test": (
            "exact paired sign permutation across Stage3B spatial regions "
            "plus 1000-region bootstrap"
        ),
        "dropout_only_whole_tissue_pairs": int(
            evidence["dropout_only_whole_tissue_spatial_signal"].sum()
        ),
        "positive_target_region_pairs": int(
            evidence["positive_target_region_effect"].sum()
        ),
        "spatial_block_consistent_pairs": int(
            evidence["spatial_blocks_consistent"].sum()
        ),
        "high_confidence_false_communication_pairs": int(
            evidence[
                "high_confidence_dropout_induced_false_communication"
            ].sum()
        ),
        "stage3b_blank_handling": "excluded before mapping and treated as abstention/NaN",
    }
    (OUT / "communication_inference_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    report = f"""# Stage3B spatial communication inference

## Completed inputs

- Selected HER2+ SC reference cells: {len(selected_cells):,}
- Raw SC genes retained: {len(sc_genes):,}
- Raw ST genes: {len(st_genes):,}
- Full SC-ST common genes: {len(common_genes):,}
- LIANA 1.7.3 consensus LR pairs analyzed: {len(lr):,}
- ST spots: {len(spots):,}
- Stage3B abstained spots: {len(blank_spots):,}

## Fixed analysis rules

1. Every row of `cell_assignment.csv` is treated as one mapped cell instance.
2. SC counts are normalized per cell to 10,000 counts and transformed by
   `log1p` before spot reconstruction.
3. Receptor complexes use the minimum expression across required subunits.
4. Spatial neighbors are the exact six-neighbor Visium lattice; no
   data-dependent distance threshold is used.
5. Global LR spatial association uses bivariate Moran R with an analytical
   one-sided p-value and Benjamini-Hochberg correction.
6. The same LR resource, normalization, graph topology, and statistic are
   used in all three conditions.
7. Stage3B blank spots remain `NaN`/abstained and are excluded from the
   SVTuner communication graph. They are never encoded as zero expression.

## Outputs

- `global_lr_communication_statistics.csv`: all LR results by condition.
- `three_condition_lr_comparison.csv`: direct three-condition comparison.
- `whole_tissue_lr_communication_statistics.csv`: Full-reference and
  CytoSPACE-dropout whole-tissue spatial association.
- `target_region_dropout_induced_lr_statistics.csv`: target-region effect
  sizes, exact spatial-block permutation, and region-bootstrap intervals.
- `false_communication_evidence_table.csv`: preregistered evidence flags for
  every LR pair; no pair is manually selected.
- `*_local_lr_statistics.npz`: spot-level LR spatial statistics for later
  hotspot and direction visualizations.
- `*_spot_lr_expression.npz`: reconstructed LR-gene expression by spot.
- `visium_six_neighbor_adjacency.npz`: shared spatial graph.

## Preregistered screening result

- Dropout-only whole-tissue spatial LR signals:
  {summary['dropout_only_whole_tissue_pairs']:,}
- LR pairs with a positive target-region effect and region-bootstrap
  interval above zero: {summary['positive_target_region_pairs']:,}
- LR pairs consistent across Stage3B spatial blocks:
  {summary['spatial_block_consistent_pairs']:,}
- LR pairs satisfying all three criteria:
  {summary['high_confidence_false_communication_pairs']:,}

These are candidates for false communication, not final biological claims.
They must still fail independent validation against observed ST receptor,
target-gene, and pathway responses before being called unsupported.
No LR pair has been manually selected.
"""
    (OUT / "INFERENCE_REPORT.md").write_text(report, encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(f"[communication] outputs written to {OUT}")


if __name__ == "__main__":
    main()
