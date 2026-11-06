"""Parity tests for the v0.4.0 cross-module improvements (ported from the MATLAB v1.18 work).

Covers: Model.full_params (#6), the profile_likelihood margin>=1 warning (#3), the clustering
not-evaluated flag (#10), the uncertainty_analysis xdata default (#7), per-row bootstrap noise and
the zero-cost-scale guard (#9), PEResult.initial_points (#5), and per-series sensitivity (#8).
"""
import warnings

import numpy as np
import pytest

import gsua_csb as g


def _two_scale(pars, t):
    t = np.asarray(t, dtype=float)
    a, b, c, d = pars
    return np.vstack([a * np.exp(-b * t), c * (1 - np.exp(-d * t))])


def _two_scale_noarg(pars):
    t = np.arange(1, 41.0)
    a, b, c, d = pars
    return np.vstack([a * np.exp(-b * t), c * (1 - np.exp(-d * t))])


def _model(fix_c=False, domainless=False):
    names = ["a", "b", "c", "d"]
    ranges = np.array([[500.0, 2000.0], [0.01, 0.5], [0.001, 0.1], [0.01, 1.0]])
    if domainless:
        m = g.UserFunctionModel(_two_scale_noarg, names=names, range=ranges,
                                output_names=["count", "frac"])
    else:
        m = g.UserFunctionModel(_two_scale, names=names, range=ranges, domain=np.arange(1, 41.0),
                                output_names=["count", "frac"])
    if fix_c:
        m.fix("c", 0.03)
    return m


# ---- Item 6: Model.full_params ----
def test_full_params_reinserts_fixed_factor():
    m = _model(fix_c=True)
    np.testing.assert_allclose(m.full_params([1000.0, 0.1, 0.3]), [1000.0, 0.1, 0.03, 0.3])
    batch = m.full_params([[1.0, 2, 3], [4.0, 5, 6]])
    assert batch.shape == (2, 4)
    np.testing.assert_allclose(batch[:, 2], [0.03, 0.03])
    with pytest.raises(ValueError):
        m.full_params([1.0, 2.0])


# ---- Item 3: margin>=1 warns ----
def test_profile_likelihood_warns_on_old_margin():
    m = _model()
    xdata = np.arange(1, 41.0)
    with pytest.warns(UserWarning, match="RAW relative noise"):
        g.profile_likelihood(m, xdata, margin=1.08, params=[0], max_iter=2, step=0.2)


# ---- Item 10: clustering flag when too few runs ----
def test_cluster_not_evaluated_flag():
    m = _model()
    estimates = np.array([[1000.0, 0.1, 0.03, 0.3], [1010.0, 0.1, 0.03, 0.3]])  # 2 runs < 3
    with pytest.warns(UserWarning, match="Clustering requested"):
        result = g.identifiability_analysis(m, estimates, cluster=True, seed=0)
    assert result.cluster.clustering_evaluated is False
    assert result.cluster.note != ""
    # cluster=False leaves it not-evaluated but with no note.
    r2 = g.identifiability_analysis(m, estimates, cluster=False)
    assert r2.cluster.clustering_evaluated is False
    assert r2.cluster.note == ""


# ---- Item 7: uncertainty_analysis xdata default for a domain-less model ----
def test_ua_xdata_defaults_to_output_length():
    m = _model(domainless=True)
    M = g.design_matrix(m, 20, seed=0)
    ua = g.uncertainty_analysis(m, M)
    assert ua.xdata is not None
    np.testing.assert_allclose(ua.xdata, np.arange(40.0))


# ---- Item 9: per-row noise + zero-cost-scale guard ----
def test_bootstrap_per_row_noise_runs_and_records():
    m = _model()
    xdata = np.arange(1, 41.0)
    ydata = _two_scale([1000.0, 0.1, 0.03, 0.3], xdata)
    pe = g.parameter_estimation(m, xdata, ydata, n=4, solver="minimize", margin=0.3, seed=2)
    bs = g.parametric_bootstrap(m, xdata, ydata, pe.x, pe.cost, b=4, n=3,
                                noise=["quasipoisson", "none"], seed=1, progress=False)
    assert bs.x.shape == (4, 4)
    assert list(bs.noise) == ["quasipoisson", "none"]


def test_costf_guards_zero_regulator_row():
    ydata = np.array([[1.0, 2.0, 3.0], [0.0, 0.0, 0.0]])
    yf = np.array([[1.0, 2.0, 3.1], [0.0, 0.0, 0.0]])
    assert np.isfinite(g.costf(ydata, yf, np.array([1.0, 0.0]), alpha=2.0))   # row 2 regulator 0
    assert np.isfinite(g.rcostf(ydata, yf, margin=0.1))                        # self-derived reg 0


# ---- Item 5: PEResult.initial_points ----
def test_pe_records_initial_points():
    m = _model()
    xdata = np.arange(1, 41.0)
    ydata = _two_scale([1000.0, 0.1, 0.03, 0.3], xdata)
    pe = g.parameter_estimation(m, xdata, ydata, n=5, solver="least_squares", seed=3)
    assert pe.initial_points.shape == (5, 4)


# ---- Item 8: per-series sensitivity ----
def test_sensitivity_per_series_returns_list():
    m = _model()
    M = g.design_matrix(m, 40, seed=0)
    results = g.sensitivity_analysis(m, M, method="xiao", output_index=None)
    assert isinstance(results, list) and len(results) == 2
    for r in results:
        assert r.Si_vec.shape[0] == m.n_params
