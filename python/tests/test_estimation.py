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


def test_margin_and_alpha_recorded_on_result(ydata):
    model = _model()
    result = parameter_estimation(model, XDATA, ydata, n=1, solver="minimize", margin=0.15, alpha=1.5)
    assert result.margin == 0.15
    assert result.alpha == 1.5


def test_margin_and_alpha_default_values_recorded(ydata):
    model = _model()
    result = parameter_estimation(model, XDATA, ydata, n=1, solver="least_squares")
    assert result.margin == 0.0
    assert result.alpha == 2.0


def test_log_scale_recovers_true_parameter_and_returns_natural_units():
    # A parameter whose true value sits many orders of magnitude below its bound's upper end --
    # the Okuonghae failure mode (theta: bound [1e-13, 1000], true value 2e-12). log_scale=True
    # searches in log10 space internally, but initial_points and the returned x must stay in
    # natural/linear units -- the transform is an internal search-space detail, not a public
    # contract change.
    true_k = 5e-6
    xdata = np.linspace(0.1, 1.0, 10)

    def f(p, d):
        # linear in k (well-conditioned regardless of k's magnitude) and rescaled by 1e6 so the
        # output -- and hence the residual/cost -- is O(1), not O(k): otherwise least_squares'
        # default absolute tolerances would call an answer 2x off from the truth "converged" simply
        # because the residual is tiny in absolute terms too.
        return p[0] * d * 1e6

    model = UserFunctionModel(
        func=f, names=["k"], range=np.array([[1e-8, 1e3]]), domain=xdata, log_scale=[True]
    )
    ydata = f(np.array([true_k]), xdata)

    # A starting point given in natural units, far from the truth but in the right ballpark on a
    # log scale -- exercises the initial_points -> log10 search-space -> natural-units round trip.
    result = parameter_estimation(
        model, xdata, ydata, n=1, solver="least_squares", initial_points=np.array([[1e-3]])
    )
    np.testing.assert_allclose(result.x[0, 0], true_k, rtol=1e-6)
    assert result.cost[0] < 1e-3


def test_log_scale_nonpositive_lower_bound_raises():
    model = UserFunctionModel(
        func=lambda p, d: p[0] * d, names=["k"], range=np.array([[0.0, 10.0]]), domain=XDATA,
        log_scale=[True],
    )
    with pytest.raises(ValueError, match="strictly positive lower bound"):
        parameter_estimation(model, XDATA, np.zeros_like(XDATA), n=1)


def test_log_scale_defaults_to_all_false_and_matches_prior_behavior(ydata):
    # No log_scale given -- must reproduce the exact pre-log_scale behavior (linear-space search).
    model = _model()
    assert not np.any(model.log_scale)
    result = parameter_estimation(
        model, XDATA, ydata, n=1, solver="least_squares", initial_points=np.array([[1.5, 0.3]])
    )
    np.testing.assert_allclose(result.x[0], TRUE_PARAMS, atol=1e-3)


class TestObjectiveAndRobustness:
    """Parity fixes: negative-margin objective, evaluation-failure handling, opt-in save."""

    @staticmethod
    def _model():
        xd = np.linspace(0.25, 24, 14)
        def pk(p, t):
            ka, ke, V = p
            return ((100 * ka) / (V * (ka - ke)) * (np.exp(-ke * t) - np.exp(-ka * t)))[None, :]
        m = UserFunctionModel(func=pk, names=["ka", "ke", "V"],
                              range=np.array([[0.6, 3.0], [0.05, 0.5], [5.0, 40.0]]),
                              nominal=np.array([1.2, 0.25, 15.0]), domain=xd, output_names=["c"])
        y = m.evaluate(np.array([1.2, 0.25, 15.0])) * (
            1 + 0.05 * np.random.default_rng(0).standard_normal((1, 14)))
        return m, xd, y

    def test_negative_margin_selects_nll_objective(self):
        # margin<0 must run the Gaussian NLL, a different objective than margin>0's regulator
        # cost -- so their reported best costs differ (they are on different scales entirely).
        m, xd, y = self._model()
        neg = parameter_estimation(m, xd, y, n=3, solver="minimize", margin=-0.1, seed=0)
        pos = parameter_estimation(m, xd, y, n=3, solver="minimize", margin=0.1, seed=0)
        assert not np.isclose(neg.cost.min(), pos.cost.min())

    def test_evaluation_failure_does_not_crash_the_run(self):
        xd = np.linspace(0.25, 24, 14)
        def flaky(p, t):
            if p[0] > 2.5:
                raise RuntimeError("unevaluable regime")
            return (p[0] * np.exp(-p[1] * t))[None, :]
        m = UserFunctionModel(func=flaky, names=["a", "r"],
                              range=np.array([[1.0, 4.0], [0.1, 1.0]]),
                              nominal=np.array([2.0, 0.5]), domain=xd, output_names=["o"])
        y = flaky(np.array([2.0, 0.5]), xd)
        for solver in ["least_squares", "minimize"]:
            pe = parameter_estimation(m, xd, y, n=6, solver=solver, seed=0)
            assert np.isfinite(pe.cost.min())          # a good start still wins
            assert pe.x[np.argmin(pe.cost), 0] == pytest.approx(2.0, abs=1e-3)

    def test_save_path_is_opt_in(self, tmp_path):
        m, xd, y = self._model()
        # nothing written without save_path
        before = set(tmp_path.iterdir())
        parameter_estimation(m, xd, y, n=2, seed=0)
        assert set(tmp_path.iterdir()) == before
        # written when asked
        fp = tmp_path / "est.npz"
        parameter_estimation(m, xd, y, n=2, seed=0, save_path=str(fp))
        with np.load(fp, allow_pickle=True) as d:
            assert "x" in d and "cost" in d and d["solver"] == "least_squares"
