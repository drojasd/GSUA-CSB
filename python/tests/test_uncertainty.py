import numpy as np
import pytest

from gsua_csb import UserFunctionModel, design_matrix, monte_carlo_filter, uncertainty_analysis


def _model(fixed_last=False):
    names = ["a", "b"]
    rng = np.array([[0.0, 1.0], [0.0, 1.0]])
    if fixed_last:
        names = ["a", "b", "c"]
        rng = np.array([[0.0, 1.0], [0.0, 1.0], [7.0, 7.0]])

    def func(p):
        return np.array([p[0] + 2.0 * p[1]])

    return UserFunctionModel(func=func, names=names, range=rng)


def test_uncertainty_analysis_shape_and_default_nominal():
    model = _model()
    M = design_matrix(model, 50, seed=0)
    result = uncertainty_analysis(model, M)
    assert result.Y.shape == (50, 1)
    # Nominal run: params = midpoint of range = [0.5, 0.5] -> 0.5 + 2*0.5 = 1.5
    np.testing.assert_allclose(result.y_nom, [1.5])


def test_uncertainty_analysis_respects_explicit_y_exp():
    model = _model()
    M = design_matrix(model, 20, seed=0)
    result = uncertainty_analysis(model, M, y_exp=[3.0])
    np.testing.assert_allclose(result.y_nom, [3.0])


def test_time_dependent_uncertainty_analysis():
    domain = np.linspace(0, 5, 8)

    def func(p, d):
        return p[0] * np.ones_like(d)

    model = UserFunctionModel(func=func, names=["k"], range=np.array([[1.0, 2.0]]), domain=domain)
    M = design_matrix(model, 10, seed=0)
    result = uncertainty_analysis(model, M)
    assert result.Y.shape == (10, 8)
    assert result.xdata.shape == (8,)


def test_monte_carlo_filter_splits_by_cost_threshold():
    model = _model()
    N = 200
    M = design_matrix(model, N, seed=0)
    ua = uncertainty_analysis(model, M)
    mcf = monte_carlo_filter(model, M, ua.Y, ua.y_nom)

    assert mcf.names == ["a", "b"]
    assert mcf.J.shape == (N,)
    # Every run is strictly classified low, high, or exactly on the threshold.
    assert mcf.low.shape[0] + mcf.high.shape[0] <= N
    assert mcf.low.shape[1] == 2
    assert mcf.high.shape[1] == 2
    # Every row placed in "low" must actually have cost below the reference, and vice versa.
    assert np.all(mcf.J[mcf.J < mcf.j_exp].shape[0] == mcf.low.shape[0])
    assert np.all(mcf.J[mcf.J > mcf.j_exp].shape[0] == mcf.high.shape[0])


def test_monte_carlo_filter_excludes_fixed_parameters():
    model = _model(fixed_last=True)
    M = design_matrix(model, 100, seed=0)
    ua = uncertainty_analysis(model, M)
    mcf = monte_carlo_filter(model, M, ua.Y, ua.y_nom)
    assert mcf.names == ["a", "b"]
    assert mcf.low.shape[1] == 2
    assert mcf.high.shape[1] == 2


def test_monte_carlo_filter_can_be_empty_without_fabricating_data():
    model = _model()
    M = design_matrix(model, 20, seed=0)
    Y = np.ones((20, 1)) * 5.0  # every run identical -> nothing strictly above/below itself
    mcf = monte_carlo_filter(model, M, Y, y_exp=[5.0])
    assert mcf.low.shape == (0, 2)
    assert mcf.high.shape == (0, 2)


def test_monte_carlo_filter_single_time_index():
    domain = np.linspace(0, 5, 6)

    def func(p, d):
        return p[0] * d

    model = UserFunctionModel(func=func, names=["k"], range=np.array([[0.5, 2.0]]), domain=domain)
    M = design_matrix(model, 100, seed=0)
    ua = uncertainty_analysis(model, M)
    mcf = monte_carlo_filter(model, M, ua.Y, ua.y_nom, t_index=3)
    assert mcf.j_exp == ua.y_nom[3]
    assert mcf.names == ["k"]
