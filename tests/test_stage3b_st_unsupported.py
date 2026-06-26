import numpy as np
import pandas as pd

from src.stages.stage3b_st_unsupported import (
    anomaly_features,
    build_spatial_edges,
    build_type_profiles,
    combine_empirical_features,
    connected_components,
    fit_nonnegative_mixtures,
    generate_supported_calibration,
    benjamini_hochberg,
    largest_component_size,
    reference_orthogonal_pvalues,
    validate_spatial_regions,
)


def test_benjamini_hochberg_preserves_order_and_monotonicity():
    pvalues = np.array([0.04, 0.001, 0.03, 0.2])
    adjusted = benjamini_hochberg(pvalues)
    assert np.allclose(adjusted, [0.0533333333, 0.004, 0.0533333333, 0.2])
    assert np.all((adjusted >= 0) & (adjusted <= 1))


def test_supported_calibration_and_novel_profile_separate():
    rng = np.random.default_rng(7)
    type_a = rng.normal([8, 7, 0, 0, 0, 0], 0.25, size=(80, 6))
    type_b = rng.normal([0, 0, 8, 7, 0, 0], 0.25, size=(80, 6))
    sc = np.maximum(np.vstack([type_a, type_b]), 0)
    labels = ["A"] * len(type_a) + ["B"] * len(type_b)
    profiles, names, pools = build_type_profiles(sc, labels)

    supported = np.vstack(
        [
            0.7 * profiles[0] + 0.3 * profiles[1],
            0.2 * profiles[0] + 0.8 * profiles[1],
        ]
        * 20
    )
    novel = np.zeros((8, 6))
    novel[:, 4:] = 1.0 / np.sqrt(2.0)
    observed = np.vstack([supported, novel])
    weights, reconstruction = fit_nonnegative_mixtures(observed, profiles)
    features = anomaly_features(observed, reconstruction)

    calibration = generate_supported_calibration(
        weights[: len(supported)],
        pools,
        names,
        400,
        rng,
    )
    _, calibration_reconstruction = fit_nonnegative_mixtures(calibration, profiles)
    calibration_features = anomaly_features(calibration, calibration_reconstruction)
    pvalues, _, _ = combine_empirical_features(features, calibration_features)
    qvalues = benjamini_hochberg(pvalues)

    assert np.median(pvalues[-len(novel) :]) < np.median(pvalues[: len(supported)])
    assert not np.any(qvalues[: len(supported)] <= 0.05)
    assert np.all(qvalues[-len(novel) :] <= 0.05)


def test_delaunay_components_join_adjacent_candidates():
    coords = pd.DataFrame(
        {"row": [0, 0, 1, 1], "col": [0, 1, 0, 1]},
        index=["s1", "s2", "s3", "s4"],
    )
    edges = build_spatial_edges(coords)
    components = connected_components(np.array([True, True, False, False]), edges)
    assert components == [[0, 1]]


def test_spatial_permutation_rejects_isolated_null_components():
    coords = pd.DataFrame(
        [(row, col) for row in range(10) for col in range(10)],
        columns=["row", "col"],
    )
    edges = build_spatial_edges(coords)
    candidate = np.zeros(100, dtype=bool)
    cluster = [
        row * 10 + col
        for row in range(3, 7)
        for col in range(3, 7)
    ]
    candidate[cluster] = True
    pvalues = np.full(100, 0.5)
    pvalues[cluster] = 1e-8

    regions, assignments = validate_spatial_regions(
        candidate,
        pvalues,
        edges,
        n_permutations=200,
        fdr=0.05,
        rng=np.random.default_rng(11),
    )

    significant = regions.loc[regions["is_unsupported_region"]]
    assert len(significant) == 1
    region_id = int(significant.iloc[0]["region_id"])
    assert np.sum(assignments == region_id) == len(cluster)


def test_spatial_permutation_keeps_candidate_strength_coupled():
    coords = pd.DataFrame(
        [(row, col) for row in range(20) for col in range(20)],
        columns=["row", "col"],
    )
    edges = build_spatial_edges(coords)
    candidate = np.zeros(400, dtype=bool)
    cluster = [
        row * 20 + col
        for row in range(5, 15)
        for col in range(5, 15)
    ]
    candidate[cluster] = True
    pvalues = np.full(400, 0.5)
    pvalues[cluster] = np.linspace(1e-10, 1e-4, len(cluster))

    regions, _ = validate_spatial_regions(
        candidate,
        pvalues,
        edges,
        n_permutations=200,
        fdr=0.05,
        rng=np.random.default_rng(17),
    )

    assert regions["is_unsupported_region"].sum() == 1


def test_reference_orthogonal_score_prefers_unrepresented_genes():
    profiles = np.array([[0.5, 0.5, 0.0], [0.4, 0.6, 0.0]])
    observations = np.array([[0.5, 0.5, 0.0], [0.2, 0.2, 0.6]])
    _, reconstruction = fit_nonnegative_mixtures(observations, profiles)

    pvalues, score, weights = reference_orthogonal_pvalues(
        observations,
        reconstruction,
        profiles,
    )

    assert weights[2] > weights[0]
    assert score[1] > score[0]
    assert pvalues[1] < pvalues[0]


def test_largest_component_size_uses_absolute_region_size():
    edges = {(0, 1), (1, 2), (3, 4)}
    assert largest_component_size(
        np.array([True, True, True, False, False]),
        edges,
    ) == 3
    assert largest_component_size(
        np.array([False, False, False, True, True]),
        edges,
    ) == 2
