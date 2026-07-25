import numpy as np
import pytest

from gsua_csb import UserFunctionModel, profile_likelihood

XDATA = np.linspace(0, 5, 25)
TRUE_PARAMS = np.array([2.0, 0.5])


def _exp_decay(p, d):
    return p[0] * np.exp(-p[1] * d)


def _model(fixed_last=False):
    names = ["amp", "k"]
    rng = np.array([[0.5, 4.0], [0.05, 1.5]])
    if fixed_last:
        names = ["amp", "k", "c"]
        rng = np.array([[0.5, 4.0], [0.05, 1.5], [9.0, 9.0]])
    model = UserFunctionModel(func=_exp_decay, names=names, range=rng, domain=XDATA)
    model.nominal = TRUE_PARAMS if not fixed_last else np.array([2.0, 0.5, 9.0])
    return model


def test_profile_likelihood_bounds_bracket_nominal():
    model = _model()
    result = profile_likelihood(model, XDATA, alpha=0.95, step=0.1, margin=0.1, max_iter=30)
    assert result.names == ["amp", "k"]
    assert result.range.shape == (2, 2)
    # Lower bound below nominal, upper bound above, for every profiled parameter.
    for i in range(2):
        assert result.range[i, 0] <= model.nominal[i]
        assert result.range[i, 1] >= model.nominal[i]
    assert result.threshold > result.nominal_cost


def test_profile_likelihood_nominal_cost_is_near_zero_for_perfect_fit():
    # ydata defaults to the nominal model's own output -> a perfect fit, so nominal_cost should
    # be small (driven only by the log(2*pi*desv) normalization term, not a residual term).
    model = _model()
    result = profile_likelihood(model, XDATA, margin=0.1)
    residual_term_only = np.sum(np.log(2 * np.pi * (model.evaluate(model.nominal, XDATA) * 0.1) ** 2)) / 2
    np.testing.assert_allclose(result.nominal_cost, residual_term_only, rtol=1e-6)


def test_profile_likelihood_respects_fixed_parameters():
    model = _model(fixed_last=True)
    result = profile_likelihood(model, XDATA, margin=0.1, max_iter=20)
    # Fixed parameter is not in the default profiled set -> its range row is untouched.
    np.testing.assert_allclose(result.range[2], [9.0, 9.0])


def test_profile_likelihood_can_restrict_to_specific_params():
    model = _model()
    original_range = model.range.copy()
    result = profile_likelihood(model, XDATA, margin=0.1, params=[0], max_iter=20)
    # Only parameter 0 was profiled -> parameter 1's range row is unchanged.
    np.testing.assert_allclose(result.range[1], original_range[1])
    assert not np.allclose(result.range[0], original_range[0])


def test_profile_likelihood_tighter_alpha_gives_narrower_interval():
    model_loose = _model()
    model_tight = _model()
    loose = profile_likelihood(model_loose, XDATA, alpha=0.99, margin=0.1, max_iter=20)
    tight = profile_likelihood(model_tight, XDATA, alpha=0.5, margin=0.1, max_iter=20)
    width_loose = loose.range[:, 1] - loose.range[:, 0]
    width_tight = tight.range[:, 1] - tight.range[:, 0]
    assert np.all(width_tight <= width_loose)


def test_profile_likelihood_explicit_ydata():
    model = _model()
    ydata = _exp_decay(TRUE_PARAMS, XDATA) + 0.01
    result = profile_likelihood(model, XDATA, ydata=ydata, margin=0.1, max_iter=20)
    assert result.range.shape == (2, 2)
    assert np.all(np.isfinite(result.range))
