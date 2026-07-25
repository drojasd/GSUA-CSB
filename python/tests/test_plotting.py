import matplotlib

matplotlib.use("Agg")  # headless backend for tests -- must be set before pyplot is imported

import numpy as np
import pytest

from gsua_csb import (
    UserFunctionModel,
    design_matrix,
    identifiability_analysis,
    monte_carlo_filter,
    plot_identifiability_correlation,
    plot_identifiability_index,
    plot_mcf,
    plot_sensitivity_area,
    plot_sensitivity_bar,
    plot_uncertainty,
    sensitivity_analysis,
    uncertainty_analysis,
)

XDATA = np.linspace(0, 5, 20)


def _model():
    def func(p, d):
        return p[0] + p[1] * np.exp(-p[2] * d)

    return UserFunctionModel(
        func=func,
        names=["a", "b", "k"],
        range=np.array([[0.0, 2.0], [0.5, 3.0], [0.1, 1.0]]),
        nominal=np.array([1.0, 2.0, 0.5]),
        domain=XDATA,
    )


def test_plot_uncertainty_draws_ensemble_and_nominal():
    model = _model()
    M = design_matrix(model, 30, seed=0)
    ua = uncertainty_analysis(model, M)
    ax = plot_uncertainty(ua)
    # N ensemble lines + 1 nominal line.
    assert len(ax.lines) == 31
    assert ax.get_legend() is not None


def test_plot_uncertainty_accepts_existing_axes():
    import matplotlib.pyplot as plt

    model = _model()
    M = design_matrix(model, 10, seed=0)
    ua = uncertainty_analysis(model, M)
    _, ax = plt.subplots()
    returned = plot_uncertainty(ua, ax=ax)
    assert returned is ax


def test_plot_sensitivity_bar_one_bar_per_parameter():
    model = _model()
    M = design_matrix(model, 40, seed=0)
    sens = sensitivity_analysis(model, M, method="xiao")
    ax = plot_sensitivity_bar(sens, index="Si")
    assert len(ax.patches) == 3
    assert len(ax.get_yticklabels()) == 3


def test_plot_sensitivity_bar_invalid_index_raises():
    model = _model()
    M = design_matrix(model, 40, seed=0)
    sens = sensitivity_analysis(model, M, method="xiao")
    with pytest.raises(AttributeError):
        plot_sensitivity_bar(sens, index="bogus")


def test_plot_sensitivity_area_stacks_time_dependent_indices():
    model = _model()
    M = design_matrix(model, 40, seed=0)
    sens = sensitivity_analysis(model, M, method="jansen")
    ax = plot_sensitivity_area(sens, XDATA, index="Si_vec")
    assert len(ax.collections) == 3  # one stacked band per parameter


def test_plot_sensitivity_area_raises_for_oat():
    model = _model()
    M = design_matrix(model, 40, seed=0)
    sens = sensitivity_analysis(model, M, method="oat")
    with pytest.raises(ValueError, match="is None"):
        plot_sensitivity_area(sens, XDATA, index="Si_vec")


def test_plot_identifiability_correlation_and_index():
    model = _model()
    rng = np.random.default_rng(0)
    estimates = rng.normal(model.nominal, 0.1, size=(50, 3))
    ia = identifiability_analysis(model, estimates)

    ax_corr = plot_identifiability_correlation(ia)
    assert ax_corr.images  # imshow produced an image artist

    ax_idx = plot_identifiability_index(ia)
    assert len(ax_idx.patches) == 3


def test_plot_identifiability_correlation_handles_fixed_parameter():
    model = _model()
    model.fix(1)
    rng = np.random.default_rng(0)
    estimates = np.column_stack(
        [rng.normal(1.0, 0.1, size=50), np.full(50, model.nominal[1]), rng.normal(0.5, 0.05, size=50)]
    )
    ia = identifiability_analysis(model, estimates)
    ax = plot_identifiability_correlation(ia)
    assert ax.images


def test_plot_mcf_one_panel_per_free_parameter():
    model = _model()
    M = design_matrix(model, 100, seed=0)
    ua = uncertainty_analysis(model, M)
    mcf = monte_carlo_filter(model, M, ua.Y, ua.y_nom)
    axes = plot_mcf(mcf, M)
    assert len(axes) >= 3
    for ax in axes[:3]:
        assert len(ax.lines) == 3  # prior, low, high ECDF steps
