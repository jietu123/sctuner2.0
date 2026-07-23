"""Stage 3B: detect ST regions unsupported by the scRNA-seq reference.

The detector is intentionally independent from Stage3A and Stage4. It does not
consume simulation truth, missing-type labels, or per-dataset whitelists.
Instead, it builds a supported-null distribution from the supplied SC
reference, calibrates every anomaly feature empirically, controls spot-level
false discovery rate, and validates spatial regions by permutation.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

_HERE = Path(__file__).resolve()
_ROOT = _HERE.parents[2]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import numpy as np
import pandas as pd
import yaml
from scipy.optimize import nnls
from scipy.spatial import Delaunay, QhullError
from scipy.stats import norm, rankdata

from src.config import load_project_config_yaml
from src.stages.storage import processed_dir, result_dir


FEATURE_NAMES = (
    "relative_reconstruction_error",
    "cosine_deficit",
    "positive_residual_fraction",
    "residual_concentration",
)


@dataclass(frozen=True)
class Stage3BConfig:
    fdr: float = 0.05
    n_calibration: int = 0
    n_spatial_permutations: int = 200
    random_seed: int = 42
    sc_expr_source: str = "normalized"
    sc_profile_source: str | None = None
    expression_scale: str = "log1p"
    sc_profile_scale: str | None = None
    max_genes: int = 0
    enable_residual_program_branch: bool = True
    residual_program_components: int = 8
    residual_program_candidate_fdr_factor: float = 2.0
    residual_program_min_whole_profile_overlap: float = 0.50
    residual_program_min_reference_orthogonal_score: float = 0.30
    enable_compensatory_diagnostics: bool = True


def _load_yaml(path: Path) -> dict:
    if not path.exists() or path.stat().st_size == 0:
        return {}
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def _resolve_dataset_config(
    project_root: Path,
    sample: str,
    project_cfg: dict,
    explicit: str | None,
) -> tuple[Path, dict]:
    if explicit:
        path = Path(explicit)
        if not path.is_absolute():
            path = project_root / path
    else:
        mapped = (project_cfg.get("dataset_config_map") or {}).get(sample)
        path = project_root / "configs" / "datasets" / (mapped or f"{sample}.yaml")
    return path.resolve(), _load_yaml(path.resolve())


def _expression_path(export_dir: Path, prefix: str, source: str) -> Path:
    candidates: dict[str, tuple[str, ...]] = {
        "normalized": (f"{prefix}_expression_normalized.csv",),
        "data": (f"{prefix}_expression_data.csv", f"{prefix}_expression_normalized.csv"),
        "counts": (f"{prefix}_expression_counts.csv", f"{prefix}_expression_normalized.csv"),
        "auto": (
            f"{prefix}_expression_normalized.csv",
            f"{prefix}_expression_data.csv",
            f"{prefix}_expression_counts.csv",
        ),
    }
    if source not in candidates:
        raise ValueError(f"Unknown expression source: {source}")
    for name in candidates[source]:
        path = export_dir / name
        if path.exists():
            return path
    raise FileNotFoundError(
        f"No {prefix.upper()} expression file found in {export_dir} "
        f"for source={source!r}"
    )


def _read_expression(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, index_col=0)
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)
    if not frame.index.is_unique:
        raise ValueError(f"Expression row identifiers are not unique: {path}")
    numeric = frame.apply(pd.to_numeric, errors="coerce")
    if numeric.isna().any().any():
        raise ValueError(f"Expression matrix contains non-numeric or missing values: {path}")
    if (numeric.to_numpy() < 0).any():
        raise ValueError(f"Stage3B requires nonnegative expression values: {path}")
    return numeric


def _read_coordinates(path: Path, spot_ids: pd.Index) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(index=spot_ids)
    coords = pd.read_csv(path)
    id_candidates = [c for c in coords.columns if str(c).lower() == "spot_id"]
    id_col = id_candidates[0] if id_candidates else coords.columns[0]
    coords[id_col] = coords[id_col].astype(str)
    coords = coords.drop_duplicates(id_col).set_index(id_col)
    return coords.reindex(spot_ids)


def _select_common_genes(
    sc: pd.DataFrame,
    st: pd.DataFrame,
    max_genes: int,
) -> list[str]:
    common = sc.columns.intersection(st.columns, sort=False)
    if common.empty:
        raise ValueError("SC and ST expression matrices have no common genes")
    if max_genes <= 0 or len(common) <= max_genes:
        return common.tolist()

    # This is a compute-budget option, not a biological cutoff. Rank by the
    # combined continuous variance and retain exactly the requested count.
    sc_var = sc.loc[:, common].var(axis=0, ddof=0)
    st_var = st.loc[:, common].var(axis=0, ddof=0)
    rank = (sc_var + st_var).sort_values(ascending=False, kind="stable")
    return rank.index[:max_genes].tolist()


def _compositional_normalize_rows(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    totals = values.sum(axis=1, keepdims=True)
    return np.divide(values, totals, out=np.zeros_like(values), where=totals > 0)


def _linearize_expression(values: np.ndarray, expression_scale: str) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    if expression_scale == "log1p":
        return np.expm1(values)
    if expression_scale == "linear":
        return values
    raise ValueError(f"Unknown expression scale: {expression_scale}")


def build_type_profiles(
    sc_values: np.ndarray,
    labels: Iterable[object],
) -> tuple[np.ndarray, list[str], dict[str, np.ndarray]]:
    labels_array = np.asarray([str(x) for x in labels], dtype=object)
    type_names = sorted(pd.unique(labels_array).tolist())
    if len(type_names) < 2:
        raise ValueError("Stage3B requires at least two SC reference cell types")

    normalized_cells = _compositional_normalize_rows(sc_values)
    cells_by_type: dict[str, np.ndarray] = {}
    profiles: list[np.ndarray] = []
    for cell_type in type_names:
        raw_type_cells = sc_values[labels_array == cell_type]
        type_cells = normalized_cells[labels_array == cell_type]
        if len(type_cells) == 0:
            continue
        cells_by_type[cell_type] = type_cells
        profiles.append(raw_type_cells.sum(axis=0))
    return _compositional_normalize_rows(np.vstack(profiles)), type_names, cells_by_type


def fit_nonnegative_mixtures(
    observations: np.ndarray,
    profiles: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    design = profiles.T
    weights = np.zeros((observations.shape[0], profiles.shape[0]), dtype=np.float64)
    reconstruction = np.zeros_like(observations, dtype=np.float64)
    for i, row in enumerate(observations):
        coef, _ = nnls(design, row)
        total = coef.sum()
        if total > 0:
            coef /= total
        weights[i] = coef
        reconstruction[i] = coef @ profiles
    return weights, reconstruction


def anomaly_features(
    observations: np.ndarray,
    reconstruction: np.ndarray,
) -> np.ndarray:
    residual = observations - reconstruction
    obs_norm = np.linalg.norm(observations, axis=1)
    residual_norm = np.linalg.norm(residual, axis=1)
    rel_error = np.divide(
        residual_norm,
        obs_norm,
        out=np.zeros_like(residual_norm),
        where=obs_norm > 0,
    )

    rec_norm = np.linalg.norm(reconstruction, axis=1)
    cosine = np.divide(
        np.sum(observations * reconstruction, axis=1),
        obs_norm * rec_norm,
        out=np.zeros_like(obs_norm),
        where=(obs_norm * rec_norm) > 0,
    )
    cosine_deficit = 1.0 - np.clip(cosine, 0.0, 1.0)

    positive = np.maximum(residual, 0.0)
    positive_mass = positive.sum(axis=1)
    observed_mass = observations.sum(axis=1)
    positive_fraction = np.divide(
        positive_mass,
        observed_mass,
        out=np.zeros_like(positive_mass),
        where=observed_mass > 0,
    )

    positive_energy = positive * positive
    energy_total = positive_energy.sum(axis=1)
    shares = np.divide(
        positive_energy,
        energy_total[:, None],
        out=np.zeros_like(positive_energy),
        where=energy_total[:, None] > 0,
    )
    concentration = np.sum(shares * shares, axis=1)
    return np.column_stack(
        [rel_error, cosine_deficit, positive_fraction, concentration]
    )


def generate_supported_calibration(
    fitted_weights: np.ndarray,
    cells_by_type: dict[str, np.ndarray],
    type_names: list[str],
    n_samples: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Generate supported pseudo-ST spots using only observed SC cells.

    Mixture compositions are sampled from fitted ST compositions, while each
    type contribution is a randomly drawn SC cell. This preserves the observed
    mixture-complexity distribution and gives a conservative null that includes
    within-type biological heterogeneity.
    """
    weight_rows = fitted_weights[
        rng.integers(0, fitted_weights.shape[0], size=n_samples)
    ]
    n_genes = next(iter(cells_by_type.values())).shape[1]
    pseudo = np.zeros((n_samples, n_genes), dtype=np.float64)
    for type_index, cell_type in enumerate(type_names):
        pool = cells_by_type[cell_type]
        draws = pool[rng.integers(0, len(pool), size=n_samples)]
        pseudo += weight_rows[:, type_index, None] * draws
    return _compositional_normalize_rows(pseudo)


def upper_tail_pvalues(
    observed: np.ndarray,
    calibration: np.ndarray,
) -> np.ndarray:
    calibration = np.sort(np.asarray(calibration, dtype=np.float64))
    observed = np.asarray(observed, dtype=np.float64)
    counts_ge = len(calibration) - np.searchsorted(calibration, observed, side="left")
    return (counts_ge + 1.0) / (len(calibration) + 1.0)


def combine_empirical_features(
    observed_features: np.ndarray,
    calibration_features: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Combine evidence after both reference-null and ST-domain calibration.

    The per-feature empirical probabilities expose disagreement with supported
    pseudo-ST. The decision score is then recentered by the robust location and
    scale of the observed ST feature distribution. This removes global
    cross-platform shifts while retaining upper-tail regional deviations.
    """
    observed_feature_p = np.column_stack(
        [
            upper_tail_pvalues(observed_features[:, j], calibration_features[:, j])
            for j in range(observed_features.shape[1])
        ]
    )
    calibration_feature_p = np.column_stack(
        [
            upper_tail_pvalues(calibration_features[:, j], calibration_features[:, j])
            for j in range(calibration_features.shape[1])
        ]
    )
    decision_feature_count = min(3, observed_features.shape[1])
    robust_z = np.zeros((len(observed_features), decision_feature_count))
    for j in range(decision_feature_count):
        values = observed_features[:, j]
        center = np.median(values)
        mad = np.median(np.abs(values - center))
        scale = 1.4826 * mad
        if scale <= np.finfo(float).eps:
            scale = np.std(values)
        if scale <= np.finfo(float).eps:
            continue
        robust_z[:, j] = (values - center) / scale

    # One-sided Stouffer aggregation: only evidence in the unsupported
    # direction contributes. The resulting normal-tail probability allows BH
    # to resolve tails smaller than 1 / n_calibration.
    observed_score = np.maximum(robust_z, 0.0).sum(axis=1) / np.sqrt(
        decision_feature_count
    )
    combined_p = norm.sf(observed_score)
    return combined_p, observed_score, observed_feature_p


def reference_orthogonal_pvalues(
    observations: np.ndarray,
    reconstruction: np.ndarray,
    profiles: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Score positive residual concentrated in genes absent from SC profiles."""
    positive = np.maximum(observations - reconstruction, 0.0)
    reference_envelope = profiles.max(axis=0)
    if len(reference_envelope) <= 1:
        novelty_weight = np.ones_like(reference_envelope)
    else:
        ranks = rankdata(reference_envelope, method="average")
        novelty_weight = 1.0 - (ranks - 1.0) / (len(reference_envelope) - 1.0)
    score = (positive * novelty_weight[None, :]).sum(axis=1)
    center = np.median(score)
    mad = np.median(np.abs(score - center))
    scale = 1.4826 * mad
    if scale <= np.finfo(float).eps:
        scale = np.std(score)
    if scale <= np.finfo(float).eps:
        pvalues = np.ones(len(score), dtype=np.float64)
    else:
        pvalues = norm.sf((score - center) / scale)
    return pvalues, score, novelty_weight


def robust_upper_tail_pvalues(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    values = np.asarray(values, dtype=np.float64)
    center = np.median(values)
    mad = np.median(np.abs(values - center))
    scale = 1.4826 * mad
    if scale <= np.finfo(float).eps:
        scale = np.std(values)
    if scale <= np.finfo(float).eps:
        z = np.zeros(len(values), dtype=np.float64)
        return np.ones(len(values), dtype=np.float64), z
    z = (values - center) / scale
    return norm.sf(z), z


def residual_program_evidence_models(
    observations: np.ndarray,
    reconstruction: np.ndarray,
    profiles: np.ndarray,
    n_components: int,
) -> tuple[dict[str, np.ndarray], pd.DataFrame, pd.DataFrame]:
    """Discover coherent positive-residual gene programs without target labels.

    Whole-profile reconstruction can miss a partially unsupported cell-type
    program when similar reference types absorb most of the spot. This branch
    factorizes novelty-weighted positive residuals and tests both signed tails
    of the leading residual components. It remains unsupervised: no marker
    whitelist, missing-type name, or simulation truth is consumed.
    """
    positive = np.maximum(observations - reconstruction, 0.0)
    reference_envelope = profiles.max(axis=0)
    if len(reference_envelope) <= 1:
        novelty_weight = np.ones_like(reference_envelope)
    else:
        ranks = rankdata(reference_envelope, method="average")
        novelty_weight = 1.0 - (ranks - 1.0) / (len(reference_envelope) - 1.0)
    weighted = positive * novelty_weight[None, :]
    positive_values = weighted[weighted > 0]
    scale = np.median(positive_values) if len(positive_values) else 1.0
    if scale <= np.finfo(float).eps:
        scale = 1.0
    residual_matrix = np.log1p(weighted / scale)
    residual_matrix = residual_matrix - residual_matrix.mean(axis=0, keepdims=True)
    max_components = max(0, min(n_components, residual_matrix.shape[0] - 1, residual_matrix.shape[1]))
    if max_components == 0:
        return {}, pd.DataFrame(index=np.arange(len(observations))), pd.DataFrame()

    u, singular_values, vt = np.linalg.svd(residual_matrix, full_matrices=False)
    evidence: dict[str, np.ndarray] = {}
    score_columns: dict[str, np.ndarray] = {}
    loading_rows: list[dict[str, object]] = []
    for component in range(max_components):
        raw_score = u[:, component] * singular_values[component]
        for sign_name, sign in (("pos", 1.0), ("neg", -1.0)):
            model_name = f"residual_program_c{component + 1}_{sign_name}"
            signed_score = raw_score * sign
            pvalues, z = robust_upper_tail_pvalues(signed_score)
            evidence[model_name] = pvalues
            score_columns[f"{model_name}_score"] = signed_score
            score_columns[f"{model_name}_z"] = z
        signed_loading = vt[component]
        for gene_index in np.argsort(np.abs(signed_loading))[::-1][:50]:
            loading_rows.append(
                {
                    "component": component + 1,
                    "gene_index": int(gene_index),
                    "loading": float(signed_loading[gene_index]),
                    "abs_loading": float(abs(signed_loading[gene_index])),
                    "novelty_weight": float(novelty_weight[gene_index]),
                }
            )
    scores = pd.DataFrame(score_columns)
    loadings = pd.DataFrame(loading_rows)
    return evidence, scores, loadings


def compensatory_substitution_diagnostics(
    observations: np.ndarray,
    profiles: np.ndarray,
    weights: np.ndarray,
    type_names: list[str],
) -> pd.DataFrame:
    """Score when fitted type usage exceeds that type's own expression support."""
    obs_norm = np.linalg.norm(observations, axis=1, keepdims=True)
    profile_norm = np.linalg.norm(profiles, axis=1)
    support = np.divide(
        observations @ profiles.T,
        obs_norm * profile_norm[None, :],
        out=np.zeros((observations.shape[0], profiles.shape[0]), dtype=np.float64),
        where=(obs_norm * profile_norm[None, :]) > 0,
    )

    def robust_z_columns(values: np.ndarray) -> np.ndarray:
        center = np.median(values, axis=0)
        mad = np.median(np.abs(values - center[None, :]), axis=0)
        scale = 1.4826 * mad
        std = values.std(axis=0)
        scale = np.where(scale <= np.finfo(float).eps, std, scale)
        scale = np.where(scale <= np.finfo(float).eps, 1.0, scale)
        return (values - center[None, :]) / scale[None, :]

    weight_z = robust_z_columns(weights)
    support_z = robust_z_columns(support)
    gap = np.maximum(weight_z, 0.0) - np.maximum(support_z, 0.0)
    gap = np.maximum(gap, 0.0)
    best_index = gap.argmax(axis=1)
    return pd.DataFrame(
        {
            "compensatory_substitution_score": gap[np.arange(len(gap)), best_index],
            "compensatory_substitute_type": [type_names[i] for i in best_index],
            "substitute_weight": weights[np.arange(len(gap)), best_index],
            "substitute_weight_z": weight_z[np.arange(len(gap)), best_index],
            "substitute_self_support_z": support_z[np.arange(len(gap)), best_index],
        }
    )


def largest_component_size(
    candidate_mask: np.ndarray,
    edges: set[tuple[int, int]],
) -> int:
    return max(
        (len(component) for component in connected_components(candidate_mask, edges)),
        default=0,
    )


def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    pvalues = np.asarray(pvalues, dtype=np.float64)
    if pvalues.ndim != 1:
        raise ValueError("pvalues must be one-dimensional")
    n = len(pvalues)
    if n == 0:
        return pvalues.copy()
    order = np.argsort(pvalues, kind="stable")
    ranked = pvalues[order]
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty(n, dtype=np.float64)
    result[order] = np.clip(adjusted, 0.0, 1.0)
    return result


def _coordinate_columns(coords: pd.DataFrame) -> tuple[str, str] | None:
    lower = {str(c).lower(): c for c in coords.columns}
    preferred = (("row", "col"), ("x", "y"), ("array_row", "array_col"))
    for left, right in preferred:
        if left in lower and right in lower:
            return str(lower[left]), str(lower[right])
    numeric = [
        c
        for c in coords.columns
        if pd.api.types.is_numeric_dtype(coords[c]) and coords[c].notna().all()
    ]
    if len(numeric) >= 2:
        return str(numeric[0]), str(numeric[1])
    return None


def build_spatial_edges(coords: pd.DataFrame) -> set[tuple[int, int]]:
    columns = _coordinate_columns(coords)
    if columns is None or len(coords) < 3:
        return set()
    points = coords.loc[:, list(columns)].to_numpy(dtype=np.float64)
    if not np.isfinite(points).all():
        return set()
    try:
        triangles = Delaunay(points).simplices
    except QhullError:
        return set()
    edges: set[tuple[int, int]] = set()
    for triangle in triangles:
        for i, j in ((0, 1), (0, 2), (1, 2)):
            left, right = sorted((int(triangle[i]), int(triangle[j])))
            edges.add((left, right))
    return edges


def connected_components(mask: np.ndarray, edges: set[tuple[int, int]]) -> list[list[int]]:
    selected = set(np.flatnonzero(mask).tolist())
    neighbors: dict[int, set[int]] = {i: set() for i in selected}
    for left, right in edges:
        if left in selected and right in selected:
            neighbors[left].add(right)
            neighbors[right].add(left)
    components: list[list[int]] = []
    while selected:
        start = selected.pop()
        component = [start]
        stack = [start]
        while stack:
            current = stack.pop()
            for neighbor in neighbors[current]:
                if neighbor in selected:
                    selected.remove(neighbor)
                    component.append(neighbor)
                    stack.append(neighbor)
        components.append(sorted(component))
    return components


def validate_spatial_regions(
    candidate_mask: np.ndarray,
    spot_pvalues: np.ndarray,
    edges: set[tuple[int, int]],
    n_permutations: int,
    fdr: float,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, np.ndarray]:
    components = connected_components(candidate_mask, edges)
    assignments = np.full(len(candidate_mask), -1, dtype=int)
    columns = [
        "region_id",
        "n_spots",
        "region_mass",
        "region_pvalue",
        "region_qvalue",
        "is_unsupported_region",
    ]
    if not components:
        return pd.DataFrame(columns=columns), assignments

    masses = np.asarray(
        [
            np.sum(-np.log10(np.clip(spot_pvalues[component], np.finfo(float).tiny, 1.0)))
            for component in components
        ]
    )
    null_max = np.zeros(n_permutations, dtype=np.float64)
    for permutation in range(n_permutations):
        # Permute the candidate decision and its anomaly strength together.
        # Independent randomization would break their defining relationship
        # and inflate null component masses, especially on dense spot graphs.
        order = rng.permutation(len(candidate_mask))
        permuted_mask = candidate_mask[order]
        permuted_p = spot_pvalues[order]
        perm_components = connected_components(permuted_mask, edges)
        if perm_components:
            null_max[permutation] = max(
                np.sum(
                    -np.log10(
                        np.clip(permuted_p[component], np.finfo(float).tiny, 1.0)
                    )
                )
                for component in perm_components
            )
    region_p = np.asarray(
        [(1.0 + np.sum(null_max >= mass)) / (n_permutations + 1.0) for mass in masses]
    )
    # Each component is compared with the maximum component mass from every
    # spatial permutation. These are already max-T family-wise adjusted
    # probabilities; applying BH again would double-correct the same family.
    region_q = region_p.copy()
    rows = []
    for region_index, component in enumerate(components, start=1):
        assignments[component] = region_index
        rows.append(
            {
                "region_id": region_index,
                "n_spots": len(component),
                "region_mass": masses[region_index - 1],
                "region_pvalue": region_p[region_index - 1],
                "region_qvalue": region_q[region_index - 1],
                "is_unsupported_region": bool(region_q[region_index - 1] <= fdr),
            }
        )
    return pd.DataFrame(rows, columns=columns), assignments


def _region_residual_markers(
    spot_ids: pd.Index,
    genes: list[str],
    observations: np.ndarray,
    reconstruction: np.ndarray,
    assignments: np.ndarray,
    regions: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    positive = np.maximum(observations - reconstruction, 0.0)
    supported_ids = set(
        regions.loc[regions["is_unsupported_region"], "region_id"].astype(int).tolist()
    )
    for region_id in sorted(supported_ids):
        member_index = np.flatnonzero(assignments == region_id)
        mean_residual = positive[member_index].mean(axis=0)
        order = np.argsort(mean_residual)[::-1]
        for rank, gene_index in enumerate(order[: min(50, len(order))], start=1):
            rows.append(
                {
                    "region_id": region_id,
                    "rank": rank,
                    "gene": genes[gene_index],
                    "mean_positive_residual": mean_residual[gene_index],
                    "n_spots": len(member_index),
                    "spot_ids": ";".join(spot_ids[member_index].astype(str)),
                }
            )
    return pd.DataFrame(rows)


def run_stage3b(
    project_root: Path,
    sample: str,
    dataset_cfg: dict,
    config: Stage3BConfig,
    output_suffix: str = "",
) -> dict:
    export_dir = processed_dir(project_root, sample, dataset_cfg) / "stage1_preprocess" / "exported"
    profile_source = config.sc_profile_source or config.sc_expr_source
    profile_scale = config.sc_profile_scale or config.expression_scale
    sc_path = _expression_path(export_dir, "sc", profile_source)
    st_path = _expression_path(export_dir, "st", config.sc_expr_source)
    metadata_path = export_dir / "sc_metadata.csv"
    coordinate_path = export_dir / "st_coordinates.csv"
    if not metadata_path.exists():
        raise FileNotFoundError(f"SC metadata not found: {metadata_path}")

    print(f"[Stage3B] Loading SC expression: {sc_path}")
    sc = _read_expression(sc_path)
    print(f"[Stage3B] Loading ST expression: {st_path}")
    st = _read_expression(st_path)
    metadata = pd.read_csv(metadata_path, index_col=0)
    metadata.index = metadata.index.astype(str)
    if "cell_type" not in metadata.columns:
        raise ValueError(f"SC metadata lacks required 'cell_type' column: {metadata_path}")

    shared_cells = sc.index.intersection(metadata.index, sort=False)
    if shared_cells.empty:
        raise ValueError("SC expression and metadata have no matching cell identifiers")
    sc = sc.loc[shared_cells]
    labels = metadata.loc[shared_cells, "cell_type"].astype(str)
    genes = _select_common_genes(sc, st, config.max_genes)
    sc_values = _compositional_normalize_rows(
        _linearize_expression(
            sc.loc[:, genes].to_numpy(dtype=np.float64),
            profile_scale,
        )
    )
    st_values = _compositional_normalize_rows(
        _linearize_expression(
            st.loc[:, genes].to_numpy(dtype=np.float64),
            config.expression_scale,
        )
    )

    profiles, type_names, cells_by_type = build_type_profiles(sc_values, labels)
    fitted_weights, reconstruction = fit_nonnegative_mixtures(st_values, profiles)
    observed_features = anomaly_features(st_values, reconstruction)

    n_calibration = config.n_calibration if config.n_calibration > 0 else len(st)
    rng = np.random.default_rng(config.random_seed)
    pseudo = generate_supported_calibration(
        fitted_weights,
        cells_by_type,
        type_names,
        n_calibration,
        rng,
    )
    _, pseudo_reconstruction = fit_nonnegative_mixtures(pseudo, profiles)
    calibration_features = anomaly_features(pseudo, pseudo_reconstruction)
    spot_p, anomaly_score, feature_p = combine_empirical_features(
        observed_features,
        calibration_features,
    )
    coords = _read_coordinates(coordinate_path, st.index)
    edges = build_spatial_edges(coords)
    orthogonal_p, orthogonal_score, _novelty_weight = reference_orthogonal_pvalues(
        st_values,
        reconstruction,
        profiles,
    )
    whole_profile_evidence = {
        "reconstruction_composite": spot_p,
        "reference_orthogonal_residual": orthogonal_p,
    }
    whole_profile_model_masks: dict[str, np.ndarray] = {}
    whole_profile_component_sizes: dict[str, int] = {}
    for model_name, model_p in whole_profile_evidence.items():
        adjusted_for_search = np.minimum(1.0, model_p * len(whole_profile_evidence))
        model_q = benjamini_hochberg(adjusted_for_search)
        model_mask = model_q <= config.fdr
        whole_profile_model_masks[model_name] = model_mask
        whole_profile_component_sizes[model_name] = largest_component_size(
            model_mask,
            edges,
        )
    selected_whole_profile_model = max(
        whole_profile_evidence,
        key=lambda name: (
            whole_profile_component_sizes[name],
            name == "reconstruction_composite",
        ),
    )
    whole_profile_spot_p = np.minimum(
        1.0,
        whole_profile_evidence[selected_whole_profile_model] * len(whole_profile_evidence),
    )
    whole_profile_spot_q = benjamini_hochberg(whole_profile_spot_p)
    whole_profile_candidate_mask = whole_profile_spot_q <= config.fdr

    whole_profile_regions, whole_profile_assignments = validate_spatial_regions(
        whole_profile_candidate_mask,
        whole_profile_spot_p,
        edges,
        config.n_spatial_permutations,
        config.fdr,
        rng,
    )
    whole_profile_region_ids = set(
        whole_profile_regions.loc[
            whole_profile_regions["is_unsupported_region"], "region_id"
        ]
        .astype(int)
        .tolist()
    )
    whole_profile_region_mask = np.asarray(
        [x in whole_profile_region_ids for x in whole_profile_assignments]
    )

    residual_program_evidence: dict[str, np.ndarray] = {}
    residual_program_scores = pd.DataFrame(index=st.index)
    residual_program_loadings = pd.DataFrame()
    selected_residual_program_model = "disabled"
    residual_program_spot_p = np.ones(len(st), dtype=np.float64)
    residual_program_spot_q = np.ones(len(st), dtype=np.float64)
    residual_program_candidate_mask = np.zeros(len(st), dtype=bool)
    region_columns = [
        "region_id",
        "n_spots",
        "region_mass",
        "region_pvalue",
        "region_qvalue",
        "is_unsupported_region",
    ]
    residual_program_regions = pd.DataFrame(columns=["program_model", *region_columns])
    residual_program_assignments = np.full(len(st), -1, dtype=int)
    residual_program_region_mask = np.zeros(len(st), dtype=bool)
    residual_program_candidate_mask_before_region_gate = np.zeros(len(st), dtype=bool)
    residual_program_region_mask_before_region_gate = np.zeros(len(st), dtype=bool)
    residual_program_region_gate_pass = np.zeros(len(st), dtype=bool)
    residual_program_component_sizes: dict[str, int] = {}
    residual_program_model_count = 0
    if config.enable_residual_program_branch and config.residual_program_components > 0:
        (
            residual_program_evidence,
            residual_program_scores_raw,
            residual_program_loadings,
        ) = residual_program_evidence_models(
            st_values,
            reconstruction,
            profiles,
            config.residual_program_components,
        )
        if residual_program_evidence:
            residual_program_model_count = len(residual_program_evidence)
            residual_program_model_masks: dict[str, np.ndarray] = {}
            for model_name, model_p in residual_program_evidence.items():
                model_q = benjamini_hochberg(model_p)
                model_mask = model_q <= config.fdr
                residual_program_model_masks[model_name] = model_mask
                residual_program_component_sizes[model_name] = largest_component_size(
                    model_mask,
                    edges,
                )
            selected_residual_program_model = max(
                residual_program_evidence,
                key=lambda name: residual_program_component_sizes[name],
            )
            residual_program_spot_p = residual_program_evidence[
                selected_residual_program_model
            ]
            residual_program_spot_q = benjamini_hochberg(residual_program_spot_p)
            residual_candidate_fdr = min(
                1.0,
                config.fdr * max(config.residual_program_candidate_fdr_factor, 1.0),
            )
            residual_program_candidate_mask = residual_program_spot_q <= residual_candidate_fdr
            residual_program_candidate_mask_before_region_gate = (
                residual_program_candidate_mask.copy()
            )
            (
                residual_program_regions,
                residual_program_assignments,
            ) = validate_spatial_regions(
                residual_program_candidate_mask,
                residual_program_spot_p,
                edges,
                config.n_spatial_permutations,
                config.fdr,
                rng,
            )
            if not residual_program_regions.empty:
                residual_program_regions = residual_program_regions.copy()
                residual_program_regions.insert(
                    0,
                    "program_model",
                    selected_residual_program_model,
                )
            residual_program_region_ids = set(
                residual_program_regions.loc[
                    residual_program_regions["is_unsupported_region"], "region_id"
                ]
                .astype(int)
                .tolist()
            )
            residual_program_region_mask = np.asarray(
                [x in residual_program_region_ids for x in residual_program_assignments]
            )
            residual_program_region_mask_before_region_gate = (
                residual_program_region_mask.copy()
            )
            if residual_program_region_ids:
                residual_program_regions = residual_program_regions.copy()
                residual_program_regions["whole_profile_overlap_fraction"] = 0.0
                residual_program_regions["median_reference_orthogonal_score"] = 0.0
                residual_program_regions["passes_residual_region_gate"] = False
                accepted_residual_program_region_ids: set[int] = set()
                for idx, row in residual_program_regions.iterrows():
                    if not bool(row.get("is_unsupported_region", False)):
                        continue
                    region_id = int(row["region_id"])
                    spot_mask = residual_program_assignments == region_id
                    if not np.any(spot_mask):
                        continue
                    overlap_fraction = float(whole_profile_region_mask[spot_mask].mean())
                    median_orthogonal = float(np.median(orthogonal_score[spot_mask]))
                    passes_gate = (
                        overlap_fraction
                        >= config.residual_program_min_whole_profile_overlap
                        or median_orthogonal
                        >= config.residual_program_min_reference_orthogonal_score
                    )
                    residual_program_regions.loc[
                        idx, "whole_profile_overlap_fraction"
                    ] = overlap_fraction
                    residual_program_regions.loc[
                        idx, "median_reference_orthogonal_score"
                    ] = median_orthogonal
                    residual_program_regions.loc[
                        idx, "passes_residual_region_gate"
                    ] = bool(passes_gate)
                    if passes_gate:
                        accepted_residual_program_region_ids.add(region_id)
                residual_program_region_gate_pass = np.asarray(
                    [
                        x in accepted_residual_program_region_ids
                        for x in residual_program_assignments
                    ]
                )
                residual_program_region_mask = (
                    residual_program_region_mask_before_region_gate
                    & residual_program_region_gate_pass
                )
                residual_program_candidate_mask = (
                    residual_program_candidate_mask_before_region_gate
                    & residual_program_region_gate_pass
                )
            residual_program_scores = residual_program_scores_raw.set_index(st.index)

    candidate_mask = whole_profile_candidate_mask | residual_program_candidate_mask
    spot_p = np.minimum(whole_profile_spot_p, residual_program_spot_p)
    spot_q = benjamini_hochberg(spot_p)
    regions, assignments = validate_spatial_regions(
        candidate_mask,
        spot_p,
        edges,
        config.n_spatial_permutations,
        config.fdr,
        rng,
    )
    significant_region_ids = set(
        regions.loc[regions["is_unsupported_region"], "region_id"].astype(int).tolist()
    )
    region_mask = np.asarray([x in significant_region_ids for x in assignments])
    region_mask = region_mask | whole_profile_region_mask | residual_program_region_mask

    calibration_median = np.median(
        calibration_features[:, FEATURE_NAMES.index("positive_residual_fraction")]
    )
    observed_positive = observed_features[
        :, FEATURE_NAMES.index("positive_residual_fraction")
    ]
    unsupported_fraction = np.maximum(
        0.0,
        (observed_positive - calibration_median) / max(1.0 - calibration_median, 1e-12),
    )

    stage3b_dir_name = "stage3b_st_unsupported" + output_suffix
    data_out = processed_dir(project_root, sample, dataset_cfg) / stage3b_dir_name
    result_out = result_dir(project_root, sample, dataset_cfg) / stage3b_dir_name
    data_out.mkdir(parents=True, exist_ok=True)
    result_out.mkdir(parents=True, exist_ok=True)

    profile_frame = pd.DataFrame(profiles, index=type_names, columns=genes)
    profile_frame.index.name = "cell_type"
    profile_frame.to_csv(data_out / "sc_type_profiles.csv")

    weight_frame = pd.DataFrame(fitted_weights, index=st.index, columns=type_names)
    weight_frame.index.name = "spot_id"
    weight_frame.to_csv(data_out / "supported_mixture_weights.csv")

    if config.enable_compensatory_diagnostics:
        compensatory = compensatory_substitution_diagnostics(
            st_values,
            profiles,
            fitted_weights,
            type_names,
        )
        compensatory.index = st.index
        compensatory.index.name = "spot_id"
    else:
        compensatory = pd.DataFrame(index=st.index)

    scores = pd.DataFrame(index=st.index)
    scores.index.name = "spot_id"
    for j, name in enumerate(FEATURE_NAMES):
        scores[name] = observed_features[:, j]
        scores[f"{name}_pvalue"] = feature_p[:, j]
    scores["combined_anomaly_score"] = anomaly_score
    scores["reference_orthogonal_score"] = orthogonal_score
    scores["reference_orthogonal_pvalue"] = orthogonal_p
    scores["whole_profile_selected_model"] = selected_whole_profile_model
    scores["whole_profile_pvalue"] = whole_profile_spot_p
    scores["whole_profile_qvalue"] = whole_profile_spot_q
    scores["is_whole_profile_candidate"] = whole_profile_candidate_mask
    scores["whole_profile_region_id"] = whole_profile_assignments
    scores["is_whole_profile_unsupported_region"] = whole_profile_region_mask
    scores["residual_program_selected_model"] = selected_residual_program_model
    scores["residual_program_pvalue"] = residual_program_spot_p
    scores["residual_program_qvalue"] = residual_program_spot_q
    scores["is_residual_program_candidate_before_region_gate"] = (
        residual_program_candidate_mask_before_region_gate
    )
    scores["is_residual_program_candidate"] = residual_program_candidate_mask
    scores["residual_program_region_id"] = residual_program_assignments
    scores["is_residual_program_unsupported_region_before_region_gate"] = (
        residual_program_region_mask_before_region_gate
    )
    scores["passes_residual_program_region_gate"] = residual_program_region_gate_pass
    scores["is_residual_program_unsupported_region"] = residual_program_region_mask
    if (
        selected_residual_program_model not in {"disabled", "none_significant"}
        and ";" not in selected_residual_program_model
        and not residual_program_scores.empty
    ):
        score_col = f"{selected_residual_program_model}_score"
        z_col = f"{selected_residual_program_model}_z"
        if score_col in residual_program_scores.columns:
            scores["residual_program_score"] = residual_program_scores[score_col]
        if z_col in residual_program_scores.columns:
            scores["residual_program_z"] = residual_program_scores[z_col]
    for column in compensatory.columns:
        scores[column] = compensatory[column]
    scores["selected_evidence_model"] = (
        "whole_profile="
        + selected_whole_profile_model
        + ";residual_program="
        + selected_residual_program_model
    )
    scores["spot_pvalue"] = spot_p
    scores["spot_qvalue"] = spot_q
    scores["is_spot_candidate"] = candidate_mask
    scores["region_id"] = assignments
    scores["is_unsupported_region"] = region_mask
    scores["unsupported_fraction_estimate"] = unsupported_fraction
    for column in coords.columns:
        if column not in scores.columns:
            scores[column] = coords[column]
    scores.to_csv(data_out / "spot_unsupported_scores.csv")
    regions.to_csv(data_out / "unsupported_regions.csv", index=False)
    whole_profile_regions.to_csv(
        data_out / "whole_profile_unsupported_regions.csv",
        index=False,
    )
    residual_program_regions.to_csv(
        data_out / "residual_program_regions.csv",
        index=False,
    )
    if not residual_program_scores.empty:
        residual_program_scores.index.name = "spot_id"
        residual_program_scores.to_csv(data_out / "residual_program_scores.csv")
    if not residual_program_loadings.empty:
        residual_program_loadings["gene"] = residual_program_loadings["gene_index"].map(
            lambda i: genes[int(i)]
        )
        residual_program_loadings.to_csv(
            data_out / "residual_program_gene_loadings.csv",
            index=False,
        )
    if not compensatory.empty:
        compensatory.to_csv(data_out / "compensatory_substitution_scores.csv")

    markers = _region_residual_markers(
        st.index,
        genes,
        st_values,
        reconstruction,
        assignments,
        regions,
    )
    markers.to_csv(data_out / "unsupported_region_residual_genes.csv", index=False)

    summary = {
        "stage": "stage3b_st_unsupported",
        "sample": sample,
        "status": "ok",
        "config": asdict(config),
        "inputs": {
            "sc_expression": str(sc_path),
            "st_expression": str(st_path),
            "sc_metadata": str(metadata_path),
            "st_coordinates": str(coordinate_path),
        },
        "dimensions": {
            "sc_cells": int(len(sc)),
            "st_spots": int(len(st)),
            "genes": int(len(genes)),
            "sc_types": int(len(type_names)),
            "calibration_spots": int(n_calibration),
        },
        "detection": {
            "spot_candidates": int(candidate_mask.sum()),
            "spatial_components": int(len(regions)),
            "unsupported_regions": int(regions["is_unsupported_region"].sum())
            if not regions.empty
            else 0,
            "unsupported_region_spots": int(region_mask.sum()),
            "spatial_edges": int(len(edges)),
            "selected_evidence_model": (
                "whole_profile="
                + selected_whole_profile_model
                + ";residual_program="
                + selected_residual_program_model
            ),
            "whole_profile": {
                "spot_candidates": int(whole_profile_candidate_mask.sum()),
                "unsupported_regions": int(
                    whole_profile_regions["is_unsupported_region"].sum()
                )
                if not whole_profile_regions.empty
                else 0,
                "unsupported_region_spots": int(whole_profile_region_mask.sum()),
                "selected_model": selected_whole_profile_model,
                "model_largest_components": whole_profile_component_sizes,
                "model_count": int(len(whole_profile_evidence)),
            },
            "residual_program": {
                "enabled": bool(config.enable_residual_program_branch),
                "region_gate": {
                    "min_whole_profile_overlap": float(
                        config.residual_program_min_whole_profile_overlap
                    ),
                    "min_reference_orthogonal_score": float(
                        config.residual_program_min_reference_orthogonal_score
                    ),
                    "candidate_spots_before_gate": int(
                        residual_program_candidate_mask_before_region_gate.sum()
                    ),
                    "unsupported_region_spots_before_gate": int(
                        residual_program_region_mask_before_region_gate.sum()
                    ),
                    "accepted_region_spots": int(residual_program_region_mask.sum()),
                },
                "spot_candidates": int(residual_program_candidate_mask.sum()),
                "unsupported_regions": int(
                    residual_program_regions["is_unsupported_region"].sum()
                )
                if not residual_program_regions.empty
                else 0,
                "unsupported_region_spots": int(residual_program_region_mask.sum()),
                "selected_model": selected_residual_program_model,
                "model_largest_components": residual_program_component_sizes,
                "model_count": int(len(residual_program_evidence)),
            },
        },
        "guardrails": {
            "uses_simulation_truth": False,
            "uses_missing_type_labels": False,
            "uses_type_whitelist": False,
            "modifies_stage3a_or_stage4": False,
            "decision_rule": (
                "whole-profile unsupported detection OR unsupervised residual-program "
                "detection; both branches use self-calibrated spot evidence and "
                "spatial permutation, without missing-type labels or cell-type whitelists"
            ),
        },
        "outputs": {
            "spot_scores": str(data_out / "spot_unsupported_scores.csv"),
            "regions": str(data_out / "unsupported_regions.csv"),
            "whole_profile_regions": str(data_out / "whole_profile_unsupported_regions.csv"),
            "residual_program_regions": str(data_out / "residual_program_regions.csv"),
            "residual_genes": str(data_out / "unsupported_region_residual_genes.csv"),
            "residual_program_scores": str(data_out / "residual_program_scores.csv"),
            "residual_program_gene_loadings": str(data_out / "residual_program_gene_loadings.csv"),
            "compensatory_substitution_scores": str(
                data_out / "compensatory_substitution_scores.csv"
            ),
            "mixture_weights": str(data_out / "supported_mixture_weights.csv"),
            "type_profiles": str(data_out / "sc_type_profiles.csv"),
        },
    }
    summary_path = result_out / "stage3b_summary.json"
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(
        "[Stage3B] Done: "
        f"{candidate_mask.sum()} spot candidates, "
        f"{summary['detection']['unsupported_regions']} spatial regions"
    )
    print(f"[Stage3B] Summary: {summary_path}")
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Stage3B: self-calibrated ST-only unsupported-region detection"
    )
    parser.add_argument("--sample", required=True, help="dataset id")
    parser.add_argument("--project_root", default=None, help="project root")
    parser.add_argument(
        "--config",
        default="configs/project_config.yaml",
        help="project config path",
    )
    parser.add_argument(
        "--dataset_config",
        default=None,
        help="dataset config path (overrides auto-detection)",
    )
    parser.add_argument(
        "--output_suffix",
        default="",
        help="suffix appended to the stage3b_st_unsupported output directory",
    )
    parser.add_argument("--fdr", type=float, default=None, help="BH FDR level")
    parser.add_argument(
        "--n_calibration",
        type=int,
        default=None,
        help="pseudo-ST calibration count; <=0 uses the observed spot count",
    )
    parser.add_argument(
        "--n_spatial_permutations",
        type=int,
        default=None,
        help="spatial region permutation count",
    )
    parser.add_argument("--random_seed", type=int, default=None)
    parser.add_argument(
        "--sc_expr_source",
        choices=["normalized", "data", "counts", "auto"],
        default=None,
    )
    parser.add_argument(
        "--sc_profile_source",
        choices=["normalized", "data", "counts", "auto"],
        default=None,
        help="optional SC-only profile source; ST continues to use sc_expr_source",
    )
    parser.add_argument(
        "--expression_scale",
        choices=["log1p", "linear"],
        default=None,
        help="scale of Stage1 expression values before compositional normalization",
    )
    parser.add_argument(
        "--sc_profile_scale",
        choices=["log1p", "linear"],
        default=None,
        help="optional scale override for the SC profile source",
    )
    parser.add_argument(
        "--max_genes",
        type=int,
        default=None,
        help="compute-budget cap; <=0 uses every common gene",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    project_root = (
        Path(args.project_root).absolute() if args.project_root else _ROOT.resolve()
    )
    project_cfg = load_project_config_yaml(project_root, args.config)
    dataset_path, dataset_cfg = _resolve_dataset_config(
        project_root,
        args.sample,
        project_cfg,
        args.dataset_config,
    )
    stage_cfg = dataset_cfg.get("stage3b") or {}
    defaults = asdict(Stage3BConfig())
    merged = {**defaults, **stage_cfg}
    for key in defaults:
        value = getattr(args, key, None)
        if value is not None:
            merged[key] = value
    config = Stage3BConfig(**merged)
    if not 0 < config.fdr < 1:
        raise ValueError("fdr must be between 0 and 1")
    if config.n_spatial_permutations < 1:
        raise ValueError("n_spatial_permutations must be positive")
    print(f"[Stage3B] Dataset config: {dataset_path}")
    run_stage3b(
        project_root,
        args.sample,
        dataset_cfg,
        config,
        output_suffix=args.output_suffix,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
