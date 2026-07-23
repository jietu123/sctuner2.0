from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pyreadr
from scipy import sparse, stats


ROOT = Path(__file__).resolve().parents[1]
COMM = ROOT / "data/processed/stage3b_spatial_communication/brca_her2_ffpe_plasma"
OUT = ROOT / "data/processed/stage3b_communication_downstream_validation/brca_her2_ffpe_plasma"
ST_H5 = (
    ROOT
    / "data/raw/cytospace_fig2d_tme/brca_ST_10x_ffpe"
    / "Visium_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix (1).h5"
)
NICHENET = (
    ROOT
    / "data/raw/spatial_communication_reference_resources"
    / "nichenet_ligand_target_matrix_v2.rds"
)
PROGENY = (
    ROOT
    / "data/raw/spatial_communication_reference_resources"
    / "omnipath_progeny_annotations_2026-06-25.csv"
)
BLANK_MANIFEST = (
    ROOT
    / "result/cytospace_fig2d_tme_brca_her2_ffpe_sc_missing_plasma_cells"
    / "stage4_cytospace_niche_stage3b_blank/cytospace_output/stage3b_blank_spots.csv"
)


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    values = np.asarray(pvalues, dtype=float)
    result = np.full(values.shape, np.nan)
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


def read_st(spots: list[str]) -> tuple[np.ndarray, list[str]]:
    with h5py.File(ST_H5, "r") as handle:
        group = handle["matrix"]
        matrix = sparse.csc_matrix(
            (group["data"][:], group["indices"][:], group["indptr"][:]),
            shape=tuple(int(value) for value in group["shape"][:]),
        )
        genes = [
            value.decode() if isinstance(value, bytes) else str(value)
            for value in group["features"]["name"][:]
        ]
        raw_spots = [
            value.decode() if isinstance(value, bytes) else str(value)
            for value in group["barcodes"][:]
        ]
    spot_index = {spot: index for index, spot in enumerate(raw_spots)}
    matrix = matrix[:, [spot_index[spot] for spot in spots]]
    totals = np.asarray(matrix.sum(axis=0)).ravel()
    totals[totals <= 0] = 1.0
    matrix = matrix @ sparse.diags(10000.0 / totals, format="csc")
    matrix.data = np.log1p(matrix.data)
    return matrix.T.toarray().astype(np.float32), [gene.upper() for gene in genes]


def zscore_columns(matrix: np.ndarray) -> np.ndarray:
    mean = matrix.mean(axis=0, keepdims=True)
    std = matrix.std(axis=0, keepdims=True)
    output = np.zeros_like(matrix, dtype=np.float32)
    valid = std.ravel() > 0
    output[:, valid] = ((matrix[:, valid] - mean[:, valid]) / std[:, valid]).astype(
        np.float32
    )
    return output


def expand_ring(
    adjacency: sparse.csr_matrix,
    core: np.ndarray,
    all_blank: np.ndarray,
    hops: int = 2,
) -> np.ndarray:
    frontier = core.copy()
    reached = core.copy()
    for _ in range(hops):
        frontier = np.asarray(adjacency @ frontier.astype(np.int8)).ravel() > 0
        frontier &= ~reached
        reached |= frontier
    ring = reached & ~core & ~all_blank
    return ring


def receiver_and_control_masks(
    adjacency: sparse.csr_matrix,
    region_core: np.ndarray,
    all_blank: np.ndarray,
    direction: str,
) -> tuple[np.ndarray, np.ndarray]:
    first_ring = expand_ring(adjacency, region_core, all_blank, hops=1)
    second_ring = expand_ring(adjacency, region_core, all_blank, hops=2)
    third_ring = expand_ring(adjacency, region_core, all_blank, hops=3)
    if direction == "incoming_to_blank_core":
        receiver = region_core
        control = second_ring
    elif direction == "outgoing_to_supported_ring":
        receiver = first_ring
        control = third_ring & ~first_ring
    else:
        raise ValueError(direction)
    control &= ~receiver & ~all_blank
    return receiver, control


def exact_region_contrast(
    score: np.ndarray,
    region_ids: np.ndarray,
    blank_mask: np.ndarray,
    adjacency: sparse.csr_matrix,
    direction: str,
) -> dict[str, float]:
    regions = sorted(np.unique(region_ids[blank_mask]))
    contrasts = []
    core_means = []
    ring_means = []
    ring_sizes = []
    for region in regions:
        region_core = blank_mask & (region_ids == region)
        receiver, control = receiver_and_control_masks(
            adjacency, region_core, blank_mask, direction
        )
        if control.sum() < max(3, receiver.sum() // 2):
            control = expand_ring(adjacency, region_core, blank_mask, hops=4)
            control &= ~receiver & ~blank_mask
        core_mean = float(np.nanmean(score[receiver]))
        ring_mean = float(np.nanmean(score[control]))
        contrasts.append(core_mean - ring_mean)
        core_means.append(core_mean)
        ring_means.append(ring_mean)
        ring_sizes.append(int(control.sum()))
    contrasts_array = np.asarray(contrasts)
    n_regions = len(regions)
    signs = np.array(
        [
            [1.0 if (mask >> bit) & 1 else -1.0 for bit in range(n_regions)]
            for mask in range(2**n_regions)
        ]
    )
    null = signs @ contrasts_array / n_regions
    observed = float(contrasts_array.mean())
    pvalue = float((np.sum(null >= observed) + 1) / (len(null) + 1))

    rng = np.random.default_rng(42)
    bootstrap = np.empty(2000)
    for index in range(len(bootstrap)):
        bootstrap[index] = rng.choice(
            contrasts_array, size=n_regions, replace=True
        ).mean()
    return {
        "core_mean": float(np.mean(core_means)),
        "local_ring_mean": float(np.mean(ring_means)),
        "core_minus_ring": observed,
        "permutation_pvalue": pvalue,
        "bootstrap_ci_low": float(np.quantile(bootstrap, 0.025)),
        "bootstrap_ci_high": float(np.quantile(bootstrap, 0.975)),
        "n_regions": n_regions,
        "median_ring_spots": float(np.median(ring_sizes)),
    }


def mapped_entity_values(
    mapped_file: np.lib.npyio.NpzFile, entities: list[str]
) -> dict[str, np.ndarray]:
    expression = mapped_file["expression"]
    genes = mapped_file["genes"].astype(str)
    gene_index = {gene: index for index, gene in enumerate(genes)}
    output = {}
    for entity in entities:
        columns = [gene_index[gene] for gene in components(entity)]
        if len(columns) == 1:
            output[entity] = expression[:, columns[0]]
        else:
            output[entity] = np.nanmin(expression[:, columns], axis=1)
    return output


def directional_received_signal(
    ligand: np.ndarray,
    receptor: np.ndarray,
    adjacency: sparse.csr_matrix,
    source_mask: np.ndarray | None = None,
) -> np.ndarray:
    source = np.nan_to_num(ligand).copy()
    if source_mask is not None:
        source[~source_mask] = 0.0
        source_neighbors = np.asarray(
            adjacency @ source_mask.astype(np.float32)
        ).ravel()
    else:
        source_neighbors = np.asarray(adjacency.sum(axis=1)).ravel()
    source_neighbors[source_neighbors <= 0] = 1.0
    neighbor_ligand = np.asarray(adjacency @ source).ravel() / source_neighbors
    return np.sqrt(np.maximum(np.nan_to_num(receptor) * neighbor_ligand, 0.0))


def infer_candidate_direction(
    ligand: str,
    receptor: str,
    full_entities: dict[str, np.ndarray],
    dropout_entities: dict[str, np.ndarray],
    adjacency: sparse.csr_matrix,
    blank_mask: np.ndarray,
) -> tuple[str, np.ndarray, np.ndarray, float, float]:
    nonblank = ~blank_mask
    full_incoming = directional_received_signal(
        full_entities[ligand],
        full_entities[receptor],
        adjacency,
        source_mask=nonblank,
    )
    dropout_incoming = directional_received_signal(
        dropout_entities[ligand],
        dropout_entities[receptor],
        adjacency,
        source_mask=nonblank,
    )
    full_outgoing = directional_received_signal(
        full_entities[ligand],
        full_entities[receptor],
        adjacency,
        source_mask=blank_mask,
    )
    dropout_outgoing = directional_received_signal(
        dropout_entities[ligand],
        dropout_entities[receptor],
        adjacency,
        source_mask=blank_mask,
    )
    first_ring = expand_ring(adjacency, blank_mask, blank_mask, hops=1)
    incoming_change = float(
        np.nanmean(dropout_incoming[blank_mask] - full_incoming[blank_mask])
    )
    outgoing_change = float(
        np.nanmean(dropout_outgoing[first_ring] - full_outgoing[first_ring])
    )
    if incoming_change >= outgoing_change:
        return (
            "incoming_to_blank_core",
            full_incoming,
            dropout_incoming,
            incoming_change,
            outgoing_change,
        )
    return (
        "outgoing_to_supported_ring",
        full_outgoing,
        dropout_outgoing,
        incoming_change,
        outgoing_change,
    )


def collect_receiver_values(
    values: np.ndarray,
    region_ids: np.ndarray,
    blank_mask: np.ndarray,
    adjacency: sparse.csr_matrix,
    direction: str,
) -> tuple[np.ndarray, np.ndarray]:
    selected_values = []
    selected_regions = []
    for region in sorted(np.unique(region_ids[blank_mask])):
        region_core = blank_mask & (region_ids == region)
        receiver, _ = receiver_and_control_masks(
            adjacency, region_core, blank_mask, direction
        )
        selected_values.append(values[receiver])
        selected_regions.append(np.full(receiver.sum(), region, dtype=int))
    return np.concatenate(selected_values), np.concatenate(selected_regions)


def receptor_score(
    expression: np.ndarray, gene_index: dict[str, int], receptor: str
) -> np.ndarray:
    columns = [gene_index[gene] for gene in components(receptor)]
    if len(columns) == 1:
        return expression[:, columns[0]]
    return np.min(expression[:, columns], axis=1)


def load_nichenet_targets(
    ligands: list[str], st_genes: set[str], top_n: int = 50
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    matrix = next(iter(pyreadr.read_r(str(NICHENET)).values()))
    matrix.index = matrix.index.astype(str).str.upper()
    matrix.columns = matrix.columns.astype(str).str.upper()
    targets = {}
    rows = []
    for ligand in ligands:
        if ligand not in matrix.columns:
            continue
        series = pd.to_numeric(matrix[ligand], errors="coerce").dropna()
        series = series.loc[series.index.isin(st_genes)]
        series = series.loc[series > 0].nlargest(top_n)
        if series.empty:
            continue
        frame = series.rename("regulatory_potential").reset_index()
        frame.columns = ["target", "regulatory_potential"]
        frame.insert(0, "ligand", ligand)
        targets[ligand] = frame
        rows.append(frame)
    combined = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    return targets, combined


def load_progeny_signatures(st_genes: set[str], top_n: int = 500) -> dict[str, pd.DataFrame]:
    raw = pd.read_csv(PROGENY)
    wide = (
        raw.pivot_table(
            index=["genesymbol", "record_id"],
            columns="label",
            values="value",
            aggfunc="first",
        )
        .reset_index()
        .rename_axis(columns=None)
    )
    wide["genesymbol"] = wide["genesymbol"].astype(str).str.upper()
    wide["weight"] = pd.to_numeric(wide["weight"], errors="coerce")
    wide = wide.dropna(subset=["pathway", "weight"])
    signatures = {}
    for pathway, frame in wide.groupby("pathway"):
        frame = frame.loc[frame["genesymbol"].isin(st_genes)].copy()
        frame = frame.iloc[np.argsort(-frame["weight"].abs().to_numpy())[:top_n]]
        signatures[str(pathway)] = frame[["genesymbol", "weight"]].drop_duplicates(
            "genesymbol"
        )
    return signatures


def weighted_activity(
    z_expression: np.ndarray,
    gene_index: dict[str, int],
    genes: list[str],
    weights: np.ndarray,
) -> np.ndarray:
    columns = [gene_index[gene] for gene in genes]
    normalized = weights / np.sum(np.abs(weights))
    return z_expression[:, columns] @ normalized


def map_ligand_pathway(
    targets: pd.DataFrame, signatures: dict[str, pd.DataFrame]
) -> tuple[str | None, float]:
    if targets is None or targets.empty:
        return None, np.nan
    target_weights = targets.set_index("target")["regulatory_potential"]
    best_pathway = None
    best_score = 0.0
    for pathway, signature in signatures.items():
        signature_weights = signature.set_index("genesymbol")["weight"].abs()
        overlap = target_weights.index.intersection(signature_weights.index)
        if len(overlap) < 2:
            continue
        score = float(
            np.sum(
                target_weights.loc[overlap].to_numpy()
                * signature_weights.loc[overlap].to_numpy()
            )
            / np.sum(target_weights.to_numpy())
        )
        if score > best_score:
            best_score = score
            best_pathway = pathway
    return best_pathway, best_score


def region_bootstrap_correlation_difference(
    full_signal: np.ndarray,
    dropout_signal: np.ndarray,
    response: np.ndarray,
    region: np.ndarray,
) -> tuple[float, float, float, float]:
    valid = np.isfinite(full_signal) & np.isfinite(dropout_signal) & np.isfinite(response)
    if (
        valid.sum() < 10
        or np.std(response[valid]) == 0
        or np.std(full_signal[valid]) == 0
        or np.std(dropout_signal[valid]) == 0
    ):
        return np.nan, np.nan, np.nan, np.nan
    full_corr = float(stats.spearmanr(full_signal[valid], response[valid]).statistic)
    dropout_corr = float(
        stats.spearmanr(dropout_signal[valid], response[valid]).statistic
    )
    regions = np.unique(region[valid])
    rng = np.random.default_rng(42)
    deltas = []
    valid_indices = np.flatnonzero(valid)
    valid_region = region[valid]
    for _ in range(1000):
        sampled = rng.choice(regions, size=len(regions), replace=True)
        positions = np.concatenate(
            [valid_indices[valid_region == item] for item in sampled]
        )
        if (
            np.std(response[positions]) == 0
            or np.std(full_signal[positions]) == 0
            or np.std(dropout_signal[positions]) == 0
        ):
            continue
        a = stats.spearmanr(full_signal[positions], response[positions]).statistic
        b = stats.spearmanr(dropout_signal[positions], response[positions]).statistic
        if np.isfinite(a) and np.isfinite(b):
            deltas.append(b - a)
    if not deltas:
        return full_corr, dropout_corr, np.nan, np.nan
    return (
        full_corr,
        dropout_corr,
        float(np.quantile(deltas, 0.025)),
        float(np.quantile(deltas, 0.975)),
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    evidence = pd.read_csv(COMM / "false_communication_evidence_table.csv")
    candidates = evidence.loc[
        evidence["high_confidence_dropout_induced_false_communication"]
    ].copy()
    pairs = pd.read_csv(COMM / "liana_consensus_lr_pairs_used.csv")
    positions = pd.read_csv(COMM / "visium_spot_positions.csv")
    spots = positions["spot_id"].astype(str).tolist()
    adjacency = sparse.load_npz(COMM / "visium_six_neighbor_adjacency.npz").tocsr()
    blank = pd.read_csv(BLANK_MANIFEST)
    blank_by_spot = blank.set_index("spot_id")
    blank_mask = np.array([spot in blank_by_spot.index for spot in spots])
    region_ids = np.full(len(spots), -1, dtype=int)
    region_ids[blank_mask] = blank_by_spot.loc[
        np.asarray(spots)[blank_mask], "region_id"
    ].to_numpy(dtype=int)

    expression, genes = read_st(spots)
    gene_index = {gene: index for index, gene in enumerate(genes)}
    z_expression = zscore_columns(expression)
    st_gene_set = set(genes)

    ligand_targets, target_table = load_nichenet_targets(
        sorted(candidates["ligand"].unique()), st_gene_set, top_n=50
    )
    target_table.to_csv(OUT / "nichenet_top50_targets_used.csv", index=False)
    signatures = load_progeny_signatures(st_gene_set)

    pathway_activity = {}
    for pathway, signature in signatures.items():
        pathway_activity[pathway] = weighted_activity(
            z_expression,
            gene_index,
            signature["genesymbol"].tolist(),
            signature["weight"].to_numpy(dtype=float),
        )

    full_mapped_file = np.load(
        COMM / "full_reference_spot_lr_expression.npz",
        allow_pickle=True,
    )
    dropout_mapped_file = np.load(
        COMM / "cytospace_dropout_spot_lr_expression.npz",
        allow_pickle=True,
    )
    candidate_entities = sorted(
        set(candidates["ligand"].astype(str))
        | set(candidates["receptor"].astype(str))
    )
    full_entities = mapped_entity_values(full_mapped_file, candidate_entities)
    dropout_entities = mapped_entity_values(
        dropout_mapped_file, candidate_entities
    )

    rows = []
    for _, candidate in candidates.iterrows():
        ligand = str(candidate["ligand"])
        receptor = str(candidate["receptor"])
        (
            communication_direction,
            full_received_signal,
            dropout_received_signal,
            incoming_induced_change,
            outgoing_induced_change,
        ) = infer_candidate_direction(
            ligand,
            receptor,
            full_entities,
            dropout_entities,
            adjacency,
            blank_mask,
        )
        receptor_values = receptor_score(expression, gene_index, receptor)
        receptor_test = exact_region_contrast(
            receptor_values,
            region_ids,
            blank_mask,
            adjacency,
            communication_direction,
        )

        targets = ligand_targets.get(ligand)
        if targets is not None:
            target_values = weighted_activity(
                z_expression,
                gene_index,
                targets["target"].tolist(),
                targets["regulatory_potential"].to_numpy(dtype=float),
            )
            target_test = exact_region_contrast(
                target_values,
                region_ids,
                blank_mask,
                adjacency,
                communication_direction,
            )
            pathway, pathway_overlap = map_ligand_pathway(targets, signatures)
        else:
            target_values = None
            target_test = {
                key: np.nan
                for key in [
                    "core_mean",
                    "local_ring_mean",
                    "core_minus_ring",
                    "permutation_pvalue",
                    "bootstrap_ci_low",
                    "bootstrap_ci_high",
                    "n_regions",
                    "median_ring_spots",
                ]
            }
            pathway, pathway_overlap = None, np.nan

        if pathway is not None:
            pathway_test = exact_region_contrast(
                pathway_activity[pathway],
                region_ids,
                blank_mask,
                adjacency,
                communication_direction,
            )
        else:
            pathway_test = {
                key: np.nan
                for key in [
                    "core_mean",
                    "local_ring_mean",
                    "core_minus_ring",
                    "permutation_pvalue",
                    "bootstrap_ci_low",
                    "bootstrap_ci_high",
                    "n_regions",
                    "median_ring_spots",
                ]
            }

        if target_values is not None:
            full_receiver_values, receiver_regions = collect_receiver_values(
                full_received_signal,
                region_ids,
                blank_mask,
                adjacency,
                communication_direction,
            )
            dropout_receiver_values, _ = collect_receiver_values(
                dropout_received_signal,
                region_ids,
                blank_mask,
                adjacency,
                communication_direction,
            )
            target_receiver_values, _ = collect_receiver_values(
                target_values,
                region_ids,
                blank_mask,
                adjacency,
                communication_direction,
            )
            correlations = region_bootstrap_correlation_difference(
                full_receiver_values,
                dropout_receiver_values,
                target_receiver_values,
                receiver_regions,
            )
        else:
            correlations = (np.nan, np.nan, np.nan, np.nan)

        row = {
            "ligand": ligand,
            "receptor": receptor,
            "inferred_communication_direction": communication_direction,
            "dropout_induced_incoming_change": incoming_induced_change,
            "dropout_induced_outgoing_change": outgoing_induced_change,
            "nichenet_prior_available": targets is not None,
            "nichenet_targets_used": 0 if targets is None else len(targets),
            "assigned_progeny_pathway": pathway,
            "ligand_pathway_overlap_score": pathway_overlap,
            "full_reference_target_response_spearman": correlations[0],
            "cytospace_dropout_target_response_spearman": correlations[1],
            "dropout_minus_full_correlation_ci_low": correlations[2],
            "dropout_minus_full_correlation_ci_high": correlations[3],
        }
        for prefix, test in [
            ("receptor", receptor_test),
            ("target", target_test),
            ("pathway", pathway_test),
        ]:
            for key, value in test.items():
                row[f"{prefix}_{key}"] = value
        rows.append(row)

    result = pd.DataFrame(rows)
    for prefix in ["receptor", "target", "pathway"]:
        result[f"{prefix}_qvalue"] = bh_fdr(
            result[f"{prefix}_permutation_pvalue"].to_numpy()
        )
        result[f"{prefix}_positive_support"] = (
            result[f"{prefix}_core_minus_ring"].gt(0)
            & result[f"{prefix}_permutation_pvalue"].le(0.05)
            & result[f"{prefix}_bootstrap_ci_low"].gt(0)
        )
        result[f"{prefix}_negative_evidence"] = (
            result[f"{prefix}_core_minus_ring"].lt(0)
            & result[f"{prefix}_bootstrap_ci_high"].lt(0)
        )

    result["downstream_positive_support"] = (
        result["target_positive_support"] | result["pathway_positive_support"]
    )
    result["combined_observed_st_pvalue"] = np.nan
    prior_mask = result["nichenet_prior_available"].astype(bool)
    result.loc[prior_mask, "combined_observed_st_pvalue"] = [
        stats.combine_pvalues(
            [
                row["receptor_permutation_pvalue"],
                row["target_permutation_pvalue"],
                row["pathway_permutation_pvalue"],
            ],
            method="fisher",
        ).pvalue
        for _, row in result.loc[prior_mask].iterrows()
    ]
    result["combined_observed_st_qvalue"] = np.nan
    result.loc[prior_mask, "combined_observed_st_qvalue"] = bh_fdr(
        result.loc[prior_mask, "combined_observed_st_pvalue"].to_numpy()
    )
    result["three_layer_positive_concordance"] = (
        result["receptor_positive_support"]
        & result["target_positive_support"]
        & result["pathway_positive_support"]
    )
    result["dropout_concordance_not_improved"] = (
        result["dropout_minus_full_correlation_ci_high"].le(0)
    )
    result["biologically_supported"] = (
        result["three_layer_positive_concordance"]
        & result["combined_observed_st_qvalue"].le(0.05)
    )
    result["high_confidence_false_communication"] = (
        ~result["receptor_positive_support"]
        & ~result["downstream_positive_support"]
        & (
            result["dropout_concordance_not_improved"]
            | result["receptor_negative_evidence"]
            | result["target_negative_evidence"]
            | result["pathway_negative_evidence"]
        )
        & result["nichenet_prior_available"]
    )
    result["classification"] = np.select(
        [
            result["biologically_supported"],
            result["high_confidence_false_communication"],
            ~result["nichenet_prior_available"],
        ],
        [
            "observed_ST_supported",
            "high_confidence_false_communication",
            "insufficient_ligand_target_prior",
        ],
        default="insufficient_observed_ST_evidence",
    )
    result = candidates.merge(result, on=["ligand", "receptor"], how="left")
    result.to_csv(OUT / "observed_st_downstream_validation.csv", index=False)
    result.loc[
        result["classification"].eq("observed_ST_supported")
    ].to_csv(OUT / "observed_st_supported_communications.csv", index=False)
    result.loc[
        result["classification"].eq("high_confidence_false_communication")
    ].to_csv(OUT / "high_confidence_false_communications.csv", index=False)
    (
        result.groupby("classification", as_index=False)
        .size()
        .rename(columns={"size": "n_pairs"})
        .to_csv(OUT / "communication_classification_summary.csv", index=False)
    )

    summary = {
        "input_candidates": int(len(result)),
        "nichenet_prior_available": int(result["nichenet_prior_available"].sum()),
        "receptor_supported": int(result["receptor_positive_support"].sum()),
        "target_program_supported": int(result["target_positive_support"].sum()),
        "pathway_supported": int(result["pathway_positive_support"].sum()),
        "observed_st_supported": int(
            result["classification"].eq("observed_ST_supported").sum()
        ),
        "high_confidence_false_communication": int(
            result["classification"].eq(
                "high_confidence_false_communication"
            ).sum()
        ),
        "insufficient_ligand_target_prior": int(
            result["classification"].eq("insufficient_ligand_target_prior").sum()
        ),
        "insufficient_observed_st_evidence": int(
            result["classification"].eq(
                "insufficient_observed_ST_evidence"
            ).sum()
        ),
        "target_definition": "NicheNet v2 top 50 positive ST-covered targets per ligand",
        "pathway_definition": (
            "PROGENy top 500 absolute-weight ST-covered genes; pathway assigned "
            "by maximum weighted overlap with the fixed NicheNet target set"
        ),
        "local_control": (
            "incoming signals: blank core versus two-hop supported ring; "
            "outgoing signals: first supported receiver ring versus outer ring"
        ),
        "statistical_unit": "Stage3B spatial region",
        "test": "exact paired sign permutation plus 2000-region bootstrap",
        "expression_source": "observed raw ST counts normalized to 10000 and log1p",
    }
    (OUT / "downstream_validation_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    report = f"""# Observed-ST downstream validation

## Independent evidence

All receptor, ligand-target, and pathway-response evidence is calculated from
the observed raw ST count matrix. No mapped-SC expression is used as
downstream validation.

For incoming communication, each Stage3B core is compared with its own
two-hop supported ring. For outgoing communication, the first supported
receiver ring is compared with its outer local ring. Exact paired sign
permutations use the eight regions as the statistical unit, with 2,000
region-bootstrap replicates for uncertainty intervals.

## Fixed priors

- Ligand-target: NicheNet v2, top 50 positive ST-covered targets per ligand.
- Pathway targets: OmniPath PROGENy, top 500 absolute-weight ST-covered genes.
- Pathway assignment: maximum weighted overlap with each ligand's fixed
  NicheNet target set.

## Results

- Input spatial false-communication candidates: {summary['input_candidates']}
- Candidates with NicheNet prior coverage:
  {summary['nichenet_prior_available']}
- Candidates with positive observed-ST receptor support:
  {summary['receptor_supported']}
- Candidates with positive NicheNet target-program support:
  {summary['target_program_supported']}
- Candidates with positive PROGENy pathway support:
  {summary['pathway_supported']}
- Classified as observed-ST supported:
  {summary['observed_st_supported']}
- Classified as high-confidence false communication:
  {summary['high_confidence_false_communication']}
- Insufficient ligand-target prior:
  {summary['insufficient_ligand_target_prior']}
- Insufficient observed-ST evidence:
  {summary['insufficient_observed_st_evidence']}

The final classification is conservative: failure to reach significance alone
does not establish false communication. A high-confidence false call also
requires independent negative evidence or failure of dropout communication to
improve spatial concordance with the observed target response.

Biological support requires concordant positive receptor, NicheNet target, and
PROGENy pathway effects, followed by Fisher evidence combination and
Benjamini-Hochberg correction across candidates with NicheNet coverage.
"""
    (OUT / "VALIDATION_REPORT.md").write_text(report, encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(f"[validation] outputs written to {OUT}")


if __name__ == "__main__":
    main()
