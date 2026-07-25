import numpy as np
import pytest

from gsua_csb import UserFunctionModel, identifiability_analysis


def _model():
    return UserFunctionModel(
        func=lambda p: p,
        names=["a", "b"],
        range=np.array([[0.0, 10.0], [0.0, 10.0]]),
    )


def test_well_identified_parameter_gets_low_index():
    rng = np.random.default_rng(0)
    # 'a' is tightly clustered around 5 (well identified); 'b' is scattered across the whole range
    # (poorly identified, and correlated with nothing).
    a = 5.0 + rng.normal(0, 0.01, size=200)
    b = rng.uniform(0, 10, size=200)
    estimates = np.column_stack([a, b])
    model = _model()
    result = identifiability_analysis(model, estimates)

    assert result.names == ["a", "b"]
    assert result.range.shape == (2, 2)
    assert result.index.shape == (2,)
    # 'a''s CI should be tight, 'b''s should span nearly the whole prior range.
    assert (result.range[0, 1] - result.range[0, 0]) < (result.range[1, 1] - result.range[1, 0])
    assert result.index[0] < result.index[1]


def test_correction_clips_to_original_range():
    # Estimates concentrated near the boundary can push the CI beyond model.range without clipping.
    rng = np.random.default_rng(1)
    a = rng.normal(0.1, 0.5, size=100)  # some estimates land below 0
    b = rng.uniform(0, 10, size=100)
    estimates = np.column_stack([a, b])
    model = _model()

    uncorrected = identifiability_analysis(model, estimates, correction=False)
    corrected = identifiability_analysis(model, estimates, correction=True)

    assert corrected.range[0, 0] >= model.range[0, 0]
    assert corrected.range[0, 1] <= model.range[0, 1]
    if uncorrected.range[0, 0] < model.range[0, 0]:
        assert corrected.range[0, 0] == model.range[0, 0]


def test_correlated_parameters_get_high_correlation_index():
    rng = np.random.default_rng(2)
    a = rng.uniform(2, 8, size=150)
    b = a + rng.normal(0, 0.05, size=150)  # near-perfectly correlated with 'a'
    estimates = np.column_stack([a, b])
    model = _model()
    result = identifiability_analysis(model, estimates)
    assert abs(result.correlation[0, 1]) > 0.9


def test_outlier_removal_drops_rows():
    rng = np.random.default_rng(3)
    a = rng.normal(5, 0.2, size=100)
    b = rng.normal(5, 0.2, size=100)
    # Inject a handful of extreme outliers.
    a = np.concatenate([a, [50.0, -50.0]])
    b = np.concatenate([b, [50.0, -50.0]])
    estimates = np.column_stack([a, b])
    model = _model()
    result = identifiability_analysis(model, estimates, outlier=True)
    assert result.n_outliers_removed >= 2


def test_no_clustering_by_default():
    rng = np.random.default_rng(4)
    estimates = rng.uniform(0, 10, size=(50, 2))
    model = _model()
    result = identifiability_analysis(model, estimates)
    assert result.cluster.num_clusters == 1
    assert result.cluster.silhouette is None
    assert result.cluster.centers.shape == (1, 2)


def test_clustering_detects_two_well_separated_basins():
    rng = np.random.default_rng(5)
    basin1 = rng.normal([2.0, 2.0], 0.1, size=(60, 2))
    basin2 = rng.normal([8.0, 8.0], 0.1, size=(60, 2))
    estimates = np.vstack([basin1, basin2])
    model = _model()
    result = identifiability_analysis(model, estimates, cluster=True, seed=0)

    assert result.cluster.num_clusters == 2
    assert result.cluster.silhouette is not None
    assert result.cluster.silhouette >= 0.5
    assert result.cluster.centers.shape == (2, 2)
    assert np.sum(result.cluster.cluster_sizes) == 120
    # The two candidate global minima should be close to the two true basin centers.
    dists_to_basin1 = np.linalg.norm(result.cluster.centers - np.array([2.0, 2.0]), axis=1)
    dists_to_basin2 = np.linalg.norm(result.cluster.centers - np.array([8.0, 8.0]), axis=1)
    assert min(dists_to_basin1) < 0.5
    assert min(dists_to_basin2) < 0.5


def test_clustering_rejects_single_unimodal_basin():
    rng = np.random.default_rng(6)
    estimates = rng.normal([5.0, 5.0], 1.0, size=(100, 2))
    model = _model()
    result = identifiability_analysis(model, estimates, cluster=True, seed=0)
    assert result.cluster.num_clusters == 1


def test_fixed_parameter_does_not_produce_nan():
    # A parameter with lb==ub (fixed, e.g. via Model.fix()) has zero range width and, since every
    # repeated estimation holds it at the same value, zero variance -- both the normalization and
    # a correlation against it are mathematically undefined (0/0). This must not poison the rest
    # of the result with NaN cascades (the original motivating case: identifiability_analysis
    # called on repeated estimates right after Model.fix() during the corrective-action step of
    # the semi-automation identification cycle).
    model = _model()
    model.range[1] = [4.0, 4.0]  # fix 'b' at 4.0, matching what Model.fix() would do
    rng = np.random.default_rng(7)
    a = rng.uniform(2, 8, size=100)
    b = np.full(100, 4.0)
    estimates = np.column_stack([a, b])

    result = identifiability_analysis(model, estimates, cluster=True, seed=0)

    assert np.isfinite(result.index[0])  # free parameter 'a' unaffected
    assert result.index[1] == 0.0  # fixed parameter reports index 0, not NaN
    np.testing.assert_allclose(result.range[1], [4.0, 4.0])
    assert np.isnan(result.correlation[0, 1])  # correlation against a zero-variance column is undefined
    assert result.correlation[1, 1] == 1.0  # diagonal stays well-defined
    assert np.all(np.isfinite(result.cluster.centers))
