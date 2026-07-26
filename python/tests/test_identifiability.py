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


def test_cost_filter_drops_bad_fit_runs_before_any_other_statistic():
    # A well-identified parameter should NOT look poorly identified just because a handful of
    # multistart runs failed to converge -- those runs' (effectively arbitrary) parameter values
    # must never reach the correlation/interval/clustering computation at all.
    rng = np.random.default_rng(0)
    good_a = 5.0 + rng.normal(0, 0.01, size=95)
    good_b = 5.0 + rng.normal(0, 0.01, size=95)
    good_cost = rng.uniform(0.0010, 0.00105, size=95)  # all within cost_rtol=0.1 of each other
    bad_a = rng.uniform(0, 10, size=5)
    bad_b = rng.uniform(0, 10, size=5)
    bad_cost = rng.uniform(50, 100, size=5)  # orders of magnitude worse than the good runs

    estimates = np.column_stack([np.concatenate([good_a, bad_a]), np.concatenate([good_b, bad_b])])
    cost = np.concatenate([good_cost, bad_cost])
    model = _model()

    filtered = identifiability_analysis(model, estimates, cost=cost)
    unfiltered = identifiability_analysis(model, estimates)

    assert filtered.n_bad_fit_removed == 5
    # Filtered result should reflect only the tight cluster -> tiny interval, low index.
    assert (filtered.range[0, 1] - filtered.range[0, 0]) < 0.5
    assert filtered.index[0] < 0.3
    # The unfiltered result, contaminated by the 5 bad-fit runs, should look distinctly worse.
    assert unfiltered.index[0] > filtered.index[0]


def test_cost_filter_no_op_when_all_costs_close():
    rng = np.random.default_rng(1)
    estimates = rng.uniform(0, 10, size=(50, 2))
    cost = rng.uniform(0.1, 0.105, size=50)  # all within cost_rtol of each other
    model = _model()
    result = identifiability_analysis(model, estimates, cost=cost)
    assert result.n_bad_fit_removed == 0


def test_cost_filter_length_mismatch_raises():
    model = _model()
    estimates = np.zeros((10, 2))
    with pytest.raises(ValueError, match="cost must have"):
        identifiability_analysis(model, estimates, cost=np.zeros(9))


def test_cost_filter_handles_near_zero_best_cost():
    # A synthetic, near-noiseless fit can have best_cost ~ 0 -- the multiplicative tolerance alone
    # would then be vacuously strict (threshold ~ 0), so cost_atol must keep the filter sane.
    rng = np.random.default_rng(2)
    estimates = rng.normal(5.0, 0.01, size=(20, 2))
    cost = np.full(20, 1e-12)
    model = _model()
    result = identifiability_analysis(model, estimates, cost=cost)
    assert result.n_bad_fit_removed == 0


def test_cost_method_gap_detects_natural_break():
    # No rtol to tune -- the "gap" method should find the same good/bad split on its own by
    # locating the largest jump in sorted (log) cost.
    rng = np.random.default_rng(9)
    good_a = 5.0 + rng.normal(0, 0.01, size=20)
    good_b = 5.0 + rng.normal(0, 0.01, size=20)
    good_cost = rng.uniform(0.0010, 0.0011, size=20)
    bad_a = rng.uniform(0, 10, size=5)
    bad_b = rng.uniform(0, 10, size=5)
    bad_cost = rng.uniform(10, 20, size=5)
    estimates = np.column_stack([np.concatenate([good_a, bad_a]), np.concatenate([good_b, bad_b])])
    cost = np.concatenate([good_cost, bad_cost])
    model = _model()
    result = identifiability_analysis(model, estimates, cost=cost, cost_method="gap")
    assert result.n_bad_fit_removed == 5


def test_cost_method_gap_keeps_everything_without_a_clear_jump():
    rng = np.random.default_rng(11)
    estimates = rng.uniform(0, 10, size=(15, 2))
    cost = np.linspace(1.0, 1.5, 15)  # smooth continuum, no gap
    model = _model()
    result = identifiability_analysis(model, estimates, cost=cost, cost_method="gap")
    assert result.n_bad_fit_removed == 0


def test_unknown_cost_method_raises():
    model = _model()
    estimates = np.zeros((5, 2))
    with pytest.raises(ValueError, match="Unknown cost_method"):
        identifiability_analysis(model, estimates, cost=np.zeros(5), cost_method="bogus")


def test_min_keep_floor_prevents_degenerate_filter():
    # A geometric cost spread means a 10% rtol would keep only the single best-cost run --
    # min_keep should force the filter to leave at least that many runs instead.
    rng = np.random.default_rng(10)
    estimates = rng.uniform(0, 10, size=(10, 2))
    cost = np.array([1.0 * (1.5**i) for i in range(10)])
    model = _model()
    result = identifiability_analysis(model, estimates, cost=cost, cost_rtol=0.1, min_keep=4)
    assert result.n_bad_fit_removed == 6  # kept exactly the min_keep=4 best-cost runs


def test_min_keep_not_triggered_when_filter_already_lenient():
    rng = np.random.default_rng(12)
    estimates = rng.uniform(0, 10, size=(20, 2))
    cost = rng.uniform(1.0, 1.02, size=20)  # all within rtol of each other
    model = _model()
    result = identifiability_analysis(model, estimates, cost=cost, min_keep=3)
    assert result.n_bad_fit_removed == 0


def test_outlier_removal_after_clustering_does_not_destroy_minor_cluster():
    # Regression test for the fix: outlier removal now runs AFTER clustering (and only within the
    # dominant cluster), so a genuine but smaller second basin must survive even with
    # outlier=True -- the old order (global Mahalanobis outlier removal before clustering) risked
    # treating the minor basin as outliers relative to the pooled mean/covariance and deleting it
    # before clustering ever got a chance to find it.
    rng = np.random.default_rng(8)
    dominant = rng.normal([2.0, 2.0], 0.1, size=(90, 2))
    minor = rng.normal([8.0, 8.0], 0.1, size=(10, 2))
    estimates = np.vstack([dominant, minor])
    model = _model()

    result = identifiability_analysis(model, estimates, cluster=True, outlier=True, seed=0)

    assert result.cluster.num_clusters == 2
    dominant_size = int(np.max(result.cluster.cluster_sizes))
    assert result.cluster.dominant_cluster == int(np.argmax(result.cluster.cluster_sizes))
    assert dominant_size >= 80  # the larger basin survived essentially intact
    # Outlier removal was scoped to the dominant cluster only -- it cannot have removed more
    # points than that cluster contained, which would be the signature of the old bug (the whole
    # minor cluster getting swept up as "outliers" relative to the global pooled distribution).
    assert result.n_outliers_removed < dominant_size
