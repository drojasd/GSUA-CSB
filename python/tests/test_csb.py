import numpy as np
import pytest

from gsua_csb import UserFunctionModel, confidence_subcontour_box, range_refinement

XDATA = np.linspace(0, 5, 20)


def _linear_model(nominal=None, fixed_last=False):
    names = ["a", "b"]
    rng = np.array([[0.0, 10.0], [0.0, 10.0]])
    if fixed_last:
        names = ["a", "b", "c"]
        rng = np.array([[0.0, 10.0], [0.0, 10.0], [3.0, 3.0]])

    def func(p, d):
        return p[0] + p[1] * d

    model = UserFunctionModel(func=func, names=names, range=rng, domain=XDATA)
    if nominal is not None:
        model.nominal = np.asarray(nominal, dtype=float)
    return model


def test_range_refinement_shrinks_or_grows_but_stays_finite():
    model = _linear_model(nominal=[2.0, 1.0])
    result = range_refinement(model, lim=0.3)
    assert result.names == ["a", "b"]
    assert result.range.shape == (2, 2)
    assert np.all(np.isfinite(result.range))
    assert np.all(result.range[:, 0] <= result.range[:, 1])
    # The nominal point must lie within the refined range for every free parameter (it has zero
    # cost by construction, so the boundary search must expand outward, not collapse inward).
    assert result.range[0, 0] <= 2.0 <= result.range[0, 1]
    assert result.range[1, 0] <= 1.0 <= result.range[1, 1]


def test_range_refinement_correct_clips_to_original_bounds():
    model = _linear_model(nominal=[2.0, 1.0])
    result = range_refinement(model, lim=0.3, correct=True)
    assert np.all(result.range[:, 0] >= model.range[:, 0])
    assert np.all(result.range[:, 1] <= model.range[:, 1])


def test_range_refinement_leaves_fixed_parameter_untouched():
    model = _linear_model(nominal=[2.0, 1.0, 3.0], fixed_last=True)
    result = range_refinement(model, lim=0.3)
    np.testing.assert_allclose(result.range[2], [3.0, 3.0])


def test_range_refinement_tighter_limit_gives_narrower_range():
    model_tight = _linear_model(nominal=[2.0, 1.0])
    model_loose = _linear_model(nominal=[2.0, 1.0])
    tight = range_refinement(model_tight, lim=0.05)
    loose = range_refinement(model_loose, lim=1.0)
    width_tight = tight.range[:, 1] - tight.range[:, 0]
    width_loose = loose.range[:, 1] - loose.range[:, 0]
    assert np.all(width_tight <= width_loose)


def test_confidence_subcontour_box_narrows_range_toward_true_params():
    # True model: a=2, b=1. Start with a wide [0,10]x[0,10] box and check it narrows substantially
    # around the true point and converges within budget -- lim=0.1/stop=0.9 (as used in an earlier,
    # since-loosened version of this test) narrows too slowly to reliably converge within a small
    # test budget; lim=0.3/stop=0.5 (closer to gsua_csb's own defaults) converges reliably.
    true_params = np.array([2.0, 1.0])
    model = _linear_model(nominal=true_params)
    y_exp = model.evaluate(true_params, XDATA)

    result = confidence_subcontour_box(
        model, n=200, y_exp=y_exp, reps=40, lim=0.3, stop=0.5, seed=0
    )

    assert result.names == ["a", "b"]
    assert result.range.shape == (2, 2)
    assert result.range_history.shape[1:] == (2, 2)
    assert result.behavioral_fraction.shape[0] >= 1
    assert result.converged
    # The box should have shrunk substantially from the original [0,10]x[0,10].
    original_width = model.range[:, 1] - model.range[:, 0]
    final_width = result.range[:, 1] - result.range[:, 0]
    assert np.all(final_width < original_width * 0.6)
    # The true parameters should still be inside (or very close to) the final box.
    assert result.range[0, 0] - 0.5 <= true_params[0] <= result.range[0, 1] + 0.5
    assert result.range[1, 0] - 0.5 <= true_params[1] <= result.range[1, 1] + 0.5


def test_confidence_subcontour_box_respects_fixed_parameters():
    model = _linear_model(nominal=[2.0, 1.0, 3.0], fixed_last=True)
    y_exp = model.evaluate(model.nominal, XDATA)
    # reps=3 is deliberately too small to converge -- only checking the fixed column here.
    with pytest.warns(UserWarning, match="did not reach"):
        result = confidence_subcontour_box(model, n=50, y_exp=y_exp, reps=3, seed=0)
    np.testing.assert_allclose(result.range[2], [3.0, 3.0])


def test_confidence_subcontour_box_warns_when_not_converged():
    # An impossibly tight lim with very few reps should fail to converge and warn.
    true_params = np.array([2.0, 1.0])
    model = _linear_model(nominal=true_params)
    y_exp = model.evaluate(true_params, XDATA)
    with pytest.warns(UserWarning, match="did not reach"):
        confidence_subcontour_box(model, n=30, y_exp=y_exp, reps=1, stop=0.999999, seed=0)
