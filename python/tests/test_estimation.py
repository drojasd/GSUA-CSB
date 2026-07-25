import numpy as np
import pytest

from gsua_csb import UserFunctionModel, parameter_estimation

XDATA = np.linspace(0, 5, 25)
TRUE_PARAMS = np.array([2.0, 0.5])  # amplitude, decay rate


def _exp_decay(p, d):
    return p[0] * np.exp(-p[1] * d)


def _model(fixed_amplitude=False):
    names = ["amp", "k"]
    rng = np.array([[0.1, 5.0], [0.01, 2.0]])
    if fixed_amplitude:
        rng[0, :] = TRUE_PARAMS[0]
    return UserFunctionModel(func=_exp_decay, names=names, range=rng, domain=XDATA)


@pytest.fixture
def ydata():
    return _exp_decay(TRUE_PARAMS, XDATA)


@pytest.mark.parametrize("solver", ["least_squares", "minimize"])
def test_local_solvers_recover_true_parameters(solver, ydata):
    model = _model()
    result = parameter_estimation(
        model, XDATA, ydata, n=1, solver=solver, initial_points=np.array([[1.5, 0.3]])
    )
    assert result.names == ["amp", "k"]
    assert result.x.shape == (1, 2)
    np.testing.assert_allclose(result.x[0], TRUE_PARAMS, atol=1e-3)
    assert result.cost[0] < 1e-6


@pytest.mark.parametrize("solver", ["differential_evolution", "dual_annealing"])
def test_global_solvers_recover_true_parameters(solver, ydata):
    model = _model()
    result = parameter_estimation(model, XDATA, ydata, n=1, solver=solver, seed=0)
    np.testing.assert_allclose(result.x[0], TRUE_PARAMS, atol=1e-2)


def test_multistart_returns_results_sorted_by_cost(ydata):
    model = _model()
    initial = np.array([[1.5, 0.3], [4.0, 1.5], [0.2, 0.05]])
    result = parameter_estimation(model, XDATA, ydata, n=3, solver="least_squares", initial_points=initial)
    assert result.x.shape == (3, 2)
    assert result.cost.shape == (3,)
    assert np.all(np.diff(result.cost) >= -1e-12)
    # Every attempt should converge to essentially the same true minimum from these start points.
    for row in result.x:
        np.testing.assert_allclose(row, TRUE_PARAMS, atol=1e-3)


def test_fixed_parameter_is_held_at_nominal(ydata):
    model = _model(fixed_amplitude=True)
    result = parameter_estimation(
        model, XDATA, ydata, n=1, solver="least_squares", initial_points=np.array([[2.0, 0.3]])
    )
    assert result.x[0, 0] == TRUE_PARAMS[0]
    np.testing.assert_allclose(result.x[0, 1], TRUE_PARAMS[1], atol=1e-3)


def test_margin_based_cost_used_for_non_least_squares_solvers(ydata):
    model = _model()
    result_mse = parameter_estimation(
        model, XDATA, ydata, n=1, solver="minimize", initial_points=np.array([[1.5, 0.3]]), margin=0.0
    )
    result_rcostf = parameter_estimation(
        model, XDATA, ydata, n=1, solver="minimize", initial_points=np.array([[1.5, 0.3]]), margin=0.1
    )
    # Both should still converge close to the true parameters regardless of which cost is used.
    np.testing.assert_allclose(result_mse.x[0], TRUE_PARAMS, atol=1e-2)
    np.testing.assert_allclose(result_rcostf.x[0], TRUE_PARAMS, atol=1e-2)


def test_default_initial_points_generated_when_not_given(ydata):
    model = _model()
    result = parameter_estimation(model, XDATA, ydata, n=2, solver="least_squares", seed=0)
    assert result.x.shape == (2, 2)


def test_unknown_solver_raises():
    model = _model()
    with pytest.raises(ValueError, match="Unknown solver"):
        parameter_estimation(model, XDATA, np.zeros_like(XDATA), solver="bogus")


def test_initial_points_row_count_mismatch_raises(ydata):
    model = _model()
    with pytest.raises(ValueError, match="initial_points must have"):
        parameter_estimation(model, XDATA, ydata, n=2, initial_points=np.array([[1.5, 0.3]]))
